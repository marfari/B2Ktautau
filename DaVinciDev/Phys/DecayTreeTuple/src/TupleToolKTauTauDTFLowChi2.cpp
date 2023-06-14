/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

// local
#include "TupleToolKTauTauDTFLowChi2.h"
#include "TMath.h"
#include <algorithm>

using namespace LHCb;

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolDecayTreeFitter
// Original Authors: Yasmine Amhis, Matt Needham, Patrick Koppenburg
// 30-10-2010, 01-04-2011
// Modified to handle missing particles in B+ -> K+ tau+ tau-, tau -> 3pi nu
// by Aravindhan Venkateswaran.
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolKTauTauDTFLowChi2::TupleToolKTauTauDTFLowChi2( const std::string& type,
                                                    const std::string& name,
                                                    const IInterface* parent )
  : TupleToolBase ( type, name, parent )
{
  declareProperty( "daughtersToConstrain", m_massConstraints,
                   "List of particles to contrain to mass");
  declareProperty( "constrainToOriginVertex", m_constrainToOriginVertex = false,
                   "Do a refit constraining to Origin Vertex (could be PV)");
  declareProperty( "Substitutions", m_map,
                   "PID-substitutions :  { ' decay-component' : 'new-pid' }" ) ;
  declareProperty( "StoreRefittedPVsTwice" ,m_storeAnyway = false,
                   "Store PV even if a refitted version is already the best PV (i.e store twice)." ) ;
  declareProperty( "UpdateDaughters" ,m_updateDaughters = false,
                   "Store updated momenta of tracks in the decay tree." ) ;
  declareProperty( "StateProvider", m_stateprovider ) ;
  declareProperty( "VetoPIDs", m_vetoPIDs,
                   "An optional list of PDG particle codes to skip when filling the tuple." );
  declareProperty( "UseFullTreeInName", m_useFullTreeInName = false,
                   "Use an improved branch naming scheme that includes a full history of the "
                   "parents and grand-parents for each particle. Makes it easier to identify "
                   "cases where the same particle type appears at different levels in the decay tree." );

  //The properties below are "custom", not present in standard TupleToolDecayTreeFitter
  declareProperty( "fillIterInfo", m_fillIterInfo = false, "Fill information from every DTF iteration");
  declareProperty( "maxNumberOfIterations", m_maxNiter = 1000, "Maximum number of iterations during fitting." ); //Standard default value for m_maxNiter is 10
  declareProperty( "maxndiverging", m_maxndiverging = 30, "Maximum number of diverging iterations during fitting." );
  declareProperty( "dChisqQuit", m_dChisqQuit = 2e5, "Maximum number of diverging iterations during fitting." ); //currently not used by code
  declareProperty( "doNuPZFix", m_doNuPZFix = true, "if true, initialize nuPZ to 0 if calculated -ve");
  declareProperty( "tauMass", m_tauMass = 1776.86, "tau mass used in initialization calculation. 1776.8199 MeV is the number used in LHCb MC, 1776.86 is the PDG average");

  declareInterface<IParticleTupleTool>(this);
}

//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::initialize()
{
  StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  // convert the list of names to a list of pids
  m_ppSvc = svc<LHCb::IParticlePropertySvc>("LHCb::ParticlePropertySvc",true) ;
  for ( const auto & S : m_massConstraints )
  {
    const auto prop = m_ppSvc->find( S );
    if (!prop)  Exception("Unknown PID");
    m_massConstraintsPids.push_back(prop->pdgID());
  }

  m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc(), this ) ;
  if ( !m_dva ) return Error("Couldn't get parent IDVAlgorithm", StatusCode::FAILURE);

  m_particleDescendants = tool<IParticleDescendants> ( "ParticleDescendants");

  if ( !m_stateprovider.empty() )
  {
    sc = m_stateprovider.retrieve() ;
    if ( sc.isFailure() ) return sc ;
  }

  if ( m_extraName.empty() )
  {
    const auto en = name() ; // use tool name as prepended name
    const auto d = en.find_last_of(".");
    m_extraName = en.substr(d+1,en.size()-1); // from d to end
    if ( "TupleToolKTauTauDTFLowChi2" == m_extraName )  m_extraName = ""; // user has not chanegd instance name
    info() << "All fields will be prepended with ``" << m_extraName << "''" <<endmsg;
  }

  if ( m_extraName.empty() )
  {
    return Error( "Extraname is empty. Always give an instance name "
                  "to TupleToolKTauTauDTFLowChi2! See doxygen." );
  }

  if ( !m_map.empty() )
  {
    m_substitute = tool<ISubstitutePID>("SubstitutePIDTool",this);
    sc = m_substitute->decodeCode( m_map );
  }

  if ( !m_vetoPIDs.empty() )
  {
    info() << "Will veto PIDs " << m_vetoPIDs << " from filling" << endmsg;
  }

  if ( m_useFullTreeInName )
  {
    info() << "Will use the full decay tree as part of branch names" << endmsg;
  }

  return sc;
}

StatusCode TupleToolKTauTauDTFLowChi2::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  if ( !m_stateprovider.empty() ) { sc = m_stateprovider.release(); }
  return StatusCode{ TupleToolBase::finalize() && sc };
}

//=============================================================================
//  The fill method implementation. This is where our modifications are happening
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fill( const LHCb::Particle* mother, const LHCb::Particle* P, const std::string& head, Tuples::Tuple& tuple ){

//  bool runningOnMC = true;

  if( !P ) return StatusCode::FAILURE;
  if ( P->isBasicParticle() )
  {
    return Error("Do not call TupleToolKTauTauDTFLowChi2 for basic particles. Use Branches. See doxygen.");
  }
  const std::string prefix = fullName(head);
  if (msgLevel(MSG::DEBUG)) debug() << "head ''" << head << "'' prefix ''" << prefix
                                    << "'' extraname ''" << m_extraName << "''" <<endmsg;

  const auto stateprovider = ( m_stateprovider.empty() ? nullptr : &(*m_stateprovider) );

  LHCb::DecayTree tree ( *P ) ;
  LHCb::Particle * treeHead = tree.head();
  
  //Get PV
  std::vector<const VertexBase*> originVtx;

  originVtx = originVertex( mother, P );
  if( originVtx.empty() ){return Error("Can't get an origin vertex");}

  //Get PV
  Gaudi::XYZVector PV( originVtx[0]->position().X(), originVtx[0]->position().Y(), originVtx[0]->position().Z() );//Notice origVtx[0] is being used i.e. only "best PV"

  //Find K+
  const LHCb::Particle * Kplus = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(321) );
  bool kplusFlag = true;

  if ( Kplus == 0 ) 
  {
    Kplus = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-321) );
    kplusFlag = false;
  } 
  if (Kplus == 0)
  {
    return Error("Can't find a K+");
  }
  else 
  {
    debug() << "Kplus :" << *Kplus << endmsg;
    debug() << " end vertex = " << Kplus->endVertex()->position() << endmsg;
    debug() << "   cov = " << Kplus->endVertex()->covMatrix() << endmsg;
  }
  
  // //Get K+ track reference point
  Gaudi::XYZPoint refPoint_Kplus = Kplus->referencePoint(); //ostensibly a point on the K+ trajectory

  // //Other method of getting reference point
  // const LHCb::ProtoParticle *proto_Kplus = Kplus->proto();
  // const LHCb::Track *track_Kplus         = proto_Kplus->track();
  // LHCb::State firstState_Kplus     = track_Kplus->firstState();
  // Gaudi::XYZPoint refPoint_Kplus_alt = firstState_Kplus.position();

  //Get Kplus 3-momentum
  Gaudi::XYZVector p_K = Kplus->momentum().Vect(); //Kplus->momentum() gives 4-momentum

  //Get taus
  //If kplusFlag is true then we want Taup with ID = -15, otherwise Taup with ID = 15

  LHCb::Particle::ConstVector tauList;
  LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList);
  LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList);

  if ( tauList.size() != 2 ) 
  {
    return Error("Can't find two tau decays");
  }

  //Get the energy and momentum of the 3pi

  double E_3pi1 = tauList[0]->momentum().E();
  double E_3pi2 = tauList[1]->momentum().E();

  double E_K    = Kplus->momentum().E();

  Gaudi::XYZVector p_3pi1 = tauList[0]->momentum().Vect();
  Gaudi::XYZVector p_3pi2 = tauList[1]->momentum().Vect();

  Gaudi::LorentzVector p4_3pi1 = tauList[0]->momentum();
  Gaudi::LorentzVector p4_3pi2 = tauList[1]->momentum();

  double pmag_3pi1 = p_3pi1.r();
  double pmag_3pi2 = p_3pi2.r();  

  // Get tau1 and tau2 end vertices
  Gaudi::XYZVector DV1(tauList[0]->endVertex()->position().X(), tauList[0]->endVertex()->position().Y(), tauList[0]->endVertex()->position().Z());
  Gaudi::XYZVector DV2(tauList[1]->endVertex()->position().X(), tauList[1]->endVertex()->position().Y(), tauList[1]->endVertex()->position().Z());

  //Get vertex fitters estimate of B decay vertex
  Gaudi::XYZVector BV(treeHead->endVertex()->position().X(), treeHead->endVertex()->position().Y(), treeHead->endVertex()->position().Z());

  //Get B flight direction
  Gaudi::XYZVector bDir = (BV - PV).Unit();

  //Get K+ momentum perpendicular to the B flight direction
  Gaudi::XYZVector p_K_perp = p_K - (p_K.Dot(bDir))*bDir;

  //Get tau flight directions
  Gaudi::XYZVector tau_dir1_strat0, tau_dir2_strat0; // both taus' directions are initialised based on vertices
  Gaudi::XYZVector tau_dir1_strat1, tau_dir2_strat1; // both taus' directions are initialised based on visible 3pi momenta
  Gaudi::XYZVector tau_dir1_strat2, tau_dir2_strat2; // both taus' directions are initialised based on vertices, using Marseille's solution for the tau momentum
  Gaudi::XYZVector tau_dir1_strat3, tau_dir2_strat3; // tau1 direction is initialised based on vertices and tau2 direction is initialised based on visible 3pi momenta
  Gaudi::XYZVector tau_dir1_strat4, tau_dir2_strat4; // tau1 direction is initialised based on visible 3pi momenta and tau2 direction is initialised based on vertices

  tau_dir1_strat0.SetXYZ(DV1.X() - BV.X(), DV1.Y() - BV.Y(), DV1.Z() - BV.Z());
  tau_dir2_strat0.SetXYZ(DV2.X() - BV.X(), DV2.Y() - BV.Y(), DV2.Z() - BV.Z());
 
  tau_dir1_strat1.SetXYZ( p_3pi1.X(), p_3pi1.Y(), p_3pi1.Z() );
  tau_dir2_strat1.SetXYZ( p_3pi2.X(), p_3pi2.Y(), p_3pi2.Z() );

  tau_dir1_strat2.SetXYZ(DV1.X() - BV.X(), DV1.Y() - BV.Y(), DV1.Z() - BV.Z());
  tau_dir2_strat2.SetXYZ(DV2.X() - BV.X(), DV2.Y() - BV.Y(), DV2.Z() - BV.Z());

  tau_dir1_strat3.SetXYZ(DV1.X() - BV.X(), DV1.Y() - BV.Y(), DV1.Z() - BV.Z());
  tau_dir2_strat3.SetXYZ(p_3pi2.X(), p_3pi2.Y(), p_3pi2.Z());

  tau_dir1_strat4.SetXYZ(p_3pi1.X(), p_3pi1.Y(), p_3pi1.Z());
  tau_dir2_strat4.SetXYZ(DV2.X() - BV.X(), DV2.Y() - BV.Y(), DV2.Z() - BV.Z());

  tau_dir1_strat0 = tau_dir1_strat0.Unit();
  tau_dir2_strat0 = tau_dir2_strat0.Unit();

  tau_dir1_strat1 = tau_dir1_strat1.Unit();
  tau_dir2_strat1 = tau_dir2_strat1.Unit();

  tau_dir1_strat2 = tau_dir1_strat2.Unit();
  tau_dir2_strat2 = tau_dir2_strat2.Unit();

  tau_dir1_strat3 = tau_dir1_strat3.Unit();
  tau_dir2_strat3 = tau_dir2_strat3.Unit();

  tau_dir1_strat4 = tau_dir1_strat4.Unit();
  tau_dir2_strat4 = tau_dir2_strat4.Unit();

  //Get tau direction unit vectors perpendicular to B flight direction
  Gaudi::XYZVector tau_dir1_perp_strat0 = (tau_dir1_strat0 - (tau_dir1_strat0.Dot(bDir))*bDir).Unit();
  Gaudi::XYZVector tau_dir2_perp_strat0 = (tau_dir2_strat0 - (tau_dir2_strat0.Dot(bDir))*bDir).Unit();

  Gaudi::XYZVector tau_dir1_perp_strat1 = (tau_dir1_strat1 - (tau_dir1_strat1.Dot(bDir))*bDir).Unit();
  Gaudi::XYZVector tau_dir2_perp_strat1 = (tau_dir2_strat1 - (tau_dir2_strat1.Dot(bDir))*bDir).Unit();

  Gaudi::XYZVector tau_dir1_perp_strat3 = (tau_dir1_strat3 - (tau_dir1_strat3.Dot(bDir))*bDir).Unit();
  Gaudi::XYZVector tau_dir2_perp_strat3 = (tau_dir2_strat3 - (tau_dir2_strat3.Dot(bDir))*bDir).Unit();

  Gaudi::XYZVector tau_dir1_perp_strat4 = (tau_dir1_strat4 - (tau_dir1_strat4.Dot(bDir))*bDir).Unit();
  Gaudi::XYZVector tau_dir2_perp_strat4 = (tau_dir2_strat4 - (tau_dir2_strat4.Dot(bDir))*bDir).Unit();

  //In plane perpendicular to B flight direction, get angles between tau momenta and K+ momentum
  double cosphi1_strat0 = tau_dir1_perp_strat0.Dot(p_K_perp.Unit());
  double cosphi2_strat0 = tau_dir2_perp_strat0.Dot(p_K_perp.Unit());

  double cosphi1_strat1 = tau_dir1_perp_strat1.Dot(p_K_perp.Unit());
  double cosphi2_strat1 = tau_dir2_perp_strat1.Dot(p_K_perp.Unit());

  double cosphi1_strat3 = tau_dir1_perp_strat3.Dot(p_K_perp.Unit());
  double cosphi2_strat3 = tau_dir2_perp_strat3.Dot(p_K_perp.Unit());

  double cosphi1_strat4 = tau_dir1_perp_strat4.Dot(p_K_perp.Unit());
  double cosphi2_strat4 = tau_dir2_perp_strat4.Dot(p_K_perp.Unit());

  double phi1_strat0 = TMath::ACos(cosphi1_strat0);
  double phi2_strat0 = TMath::ACos(cosphi2_strat0);

  double phi1_strat1 = TMath::ACos(cosphi1_strat1);
  double phi2_strat1 = TMath::ACos(cosphi2_strat1);

  double phi1_strat3 = TMath::ACos(cosphi1_strat3);
  double phi2_strat3 = TMath::ACos(cosphi2_strat3);

  double phi1_strat4 = TMath::ACos(cosphi1_strat4);
  double phi2_strat4 = TMath::ACos(cosphi2_strat4);

  //In this plane, get directions of tau momenta perpendicular to K+ momentum

  Gaudi::XYZVector tau_perp1_perpK_strat0 = tau_dir1_perp_strat0 - (tau_dir1_perp_strat0.Dot(p_K_perp.Unit())*p_K_perp.Unit());
  Gaudi::XYZVector tau_perp2_perpK_strat0 = tau_dir2_perp_strat0 - (tau_dir2_perp_strat0.Dot(p_K_perp.Unit())*p_K_perp.Unit());

  Gaudi::XYZVector tau_perp1_perpK_strat1 = tau_dir1_perp_strat1 - (tau_dir1_perp_strat1.Dot(p_K_perp.Unit())*p_K_perp.Unit());
  Gaudi::XYZVector tau_perp2_perpK_strat1 = tau_dir2_perp_strat1 - (tau_dir2_perp_strat1.Dot(p_K_perp.Unit())*p_K_perp.Unit());

  Gaudi::XYZVector tau_perp1_perpK_strat3 = tau_dir1_perp_strat3 - (tau_dir1_perp_strat3.Dot(p_K_perp.Unit())*p_K_perp.Unit());
  Gaudi::XYZVector tau_perp2_perpK_strat3 = tau_dir2_perp_strat3 - (tau_dir2_perp_strat3.Dot(p_K_perp.Unit())*p_K_perp.Unit());

  Gaudi::XYZVector tau_perp1_perpK_strat4 = tau_dir1_perp_strat4 - (tau_dir1_perp_strat4.Dot(p_K_perp.Unit())*p_K_perp.Unit());
  Gaudi::XYZVector tau_perp2_perpK_strat4 = tau_dir2_perp_strat4 - (tau_dir2_perp_strat4.Dot(p_K_perp.Unit())*p_K_perp.Unit());

  double tau_pPerp_ratio_strat0 = tau_perp1_perpK_strat0.R()/tau_perp2_perpK_strat0.R();
  double tau_pPerp_ratio_strat1 = tau_perp1_perpK_strat1.R()/tau_perp2_perpK_strat1.R();
  double tau_pPerp_ratio_strat3 = tau_perp1_perpK_strat3.R()/tau_perp2_perpK_strat3.R();
  double tau_pPerp_ratio_strat4 = tau_perp1_perpK_strat4.R()/tau_perp2_perpK_strat4.R();

  //Calculate momentum component of taus in this plane
  double pMag_tau1_perp_strat0 = -1*p_K_perp.R()/(cosphi1_strat0 + (cosphi2_strat0*(sin(phi1_strat0)/sin(phi2_strat0))));
  double pMag_tau2_perp_strat0 = pMag_tau1_perp_strat0*tau_pPerp_ratio_strat0;

  double pMag_tau1_perp_strat1 = -1*p_K_perp.R()/(cosphi1_strat1 + (cosphi2_strat1*(sin(phi1_strat1)/sin(phi2_strat1))));
  double pMag_tau2_perp_strat1 = pMag_tau1_perp_strat1*tau_pPerp_ratio_strat1;
  
  double pMag_tau1_perp_strat3 = -1*p_K_perp.R()/(cosphi1_strat3 + (cosphi2_strat3*(sin(phi1_strat3)/sin(phi2_strat3))));
  double pMag_tau2_perp_strat3 = pMag_tau1_perp_strat3*tau_pPerp_ratio_strat3;

  double pMag_tau1_perp_strat4 = -1*p_K_perp.R()/(cosphi1_strat4 + (cosphi2_strat4*(sin(phi1_strat4)/sin(phi2_strat4))));
  double pMag_tau2_perp_strat4 = pMag_tau1_perp_strat4*tau_pPerp_ratio_strat4;
  
  Gaudi::XYZVector p_tau1_perp_strat0 = pMag_tau1_perp_strat0*tau_dir1_perp_strat0;
  Gaudi::XYZVector p_tau2_perp_strat0 = pMag_tau2_perp_strat0*tau_dir2_perp_strat0;

  Gaudi::XYZVector p_tau1_perp_strat1 = pMag_tau1_perp_strat1*tau_dir1_perp_strat1;
  Gaudi::XYZVector p_tau2_perp_strat1 = pMag_tau2_perp_strat1*tau_dir2_perp_strat1;

  Gaudi::XYZVector p_tau1_perp_strat3 = pMag_tau1_perp_strat3*tau_dir1_perp_strat3;
  Gaudi::XYZVector p_tau2_perp_strat3 = pMag_tau2_perp_strat3*tau_dir2_perp_strat3;

  Gaudi::XYZVector p_tau1_perp_strat4 = pMag_tau1_perp_strat4*tau_dir1_perp_strat4;
  Gaudi::XYZVector p_tau2_perp_strat4 = pMag_tau2_perp_strat4*tau_dir2_perp_strat4;

  //Get angles made by tau directions with B flight direction
  double tau_B_cos1_strat0 = (tau_dir1_strat0.Unit()).Dot(bDir.Unit());
  double tau_B_cos2_strat0 = (tau_dir2_strat0.Unit()).Dot(bDir.Unit());

  double tau_B_cos1_strat1 = (tau_dir1_strat1.Unit()).Dot(bDir.Unit());
  double tau_B_cos2_strat1 = (tau_dir2_strat1.Unit()).Dot(bDir.Unit());

  double tau_B_cos1_strat3 = (tau_dir1_strat3.Unit()).Dot(bDir.Unit());
  double tau_B_cos2_strat3 = (tau_dir2_strat3.Unit()).Dot(bDir.Unit());

  double tau_B_cos1_strat4 = (tau_dir1_strat4.Unit()).Dot(bDir.Unit());
  double tau_B_cos2_strat4 = (tau_dir2_strat4.Unit()).Dot(bDir.Unit());

  //Get tau momenta parallel to B flight direction
  double pMag_tau1_long_strat0 = fabs(pMag_tau1_perp_strat0)*tau_B_cos1_strat0/sqrt(1-pow(tau_B_cos1_strat0,2));
  double pMag_tau2_long_strat0 = fabs(pMag_tau2_perp_strat0)*tau_B_cos2_strat0/sqrt(1-pow(tau_B_cos2_strat0,2));

  double pMag_tau1_long_strat1 = fabs(pMag_tau1_perp_strat1)*tau_B_cos1_strat1/sqrt(1-pow(tau_B_cos1_strat1,2));
  double pMag_tau2_long_strat1 = fabs(pMag_tau2_perp_strat1)*tau_B_cos2_strat1/sqrt(1-pow(tau_B_cos2_strat1,2));

  double pMag_tau1_long_strat3 = fabs(pMag_tau1_perp_strat3)*tau_B_cos1_strat3/sqrt(1-pow(tau_B_cos1_strat3,2));
  double pMag_tau2_long_strat3 = fabs(pMag_tau2_perp_strat3)*tau_B_cos2_strat3/sqrt(1-pow(tau_B_cos2_strat3,2));

  double pMag_tau1_long_strat4 = fabs(pMag_tau1_perp_strat4)*tau_B_cos1_strat4/sqrt(1-pow(tau_B_cos1_strat4,2));
  double pMag_tau2_long_strat4 = fabs(pMag_tau2_perp_strat4)*tau_B_cos2_strat4/sqrt(1-pow(tau_B_cos2_strat4,2));

  //Set total tau momentum vector
  Gaudi::XYZVector p_tau1_strat0 = p_tau1_perp_strat0 + (pMag_tau1_long_strat0*bDir);
  Gaudi::XYZVector p_tau2_strat0 = p_tau2_perp_strat0 + (pMag_tau2_long_strat0*bDir);

  Gaudi::XYZVector p_tau1_strat1 = p_tau1_perp_strat1 + (pMag_tau1_long_strat1*bDir);
  Gaudi::XYZVector p_tau2_strat1 = p_tau2_perp_strat1 + (pMag_tau2_long_strat1*bDir);

  Gaudi::XYZVector p_tau1_strat3 = p_tau1_perp_strat3 + (pMag_tau1_long_strat3*bDir);
  Gaudi::XYZVector p_tau2_strat3 = p_tau2_perp_strat3 + (pMag_tau2_long_strat3*bDir);

  Gaudi::XYZVector p_tau1_strat4 = p_tau1_perp_strat4 + (pMag_tau1_long_strat4*bDir);
  Gaudi::XYZVector p_tau2_strat4 = p_tau2_perp_strat4 + (pMag_tau2_long_strat4*bDir);

  //### RD* initialization
  double m_3pi1_sq = pow(E_3pi1, 2) - pow(pmag_3pi1, 2);
  double m_3pi2_sq = pow(E_3pi2, 2) - pow(pmag_3pi2, 2);
  
  double m_3pi1 = 0., m_3pi2 = 0.;

  if(m_3pi1_sq > 0) { m_3pi1 = sqrt(m_3pi1_sq);}
  if(m_3pi2_sq > 0) { m_3pi2 = sqrt(m_3pi2_sq);}

  double thetaMax_tau_3pi1 = TMath::ASin( (pow(m_tauMass,2) - pow(m_3pi1, 2))/(2*m_tauMass*pmag_3pi1) );
  double thetaMax_tau_3pi2 = TMath::ASin( (pow(m_tauMass,2) - pow(m_3pi2, 2))/(2*m_tauMass*pmag_3pi2) );

  double pmag_tau1_rdInit = (pow(m_tauMass,2) + pow(m_3pi1, 2))* pmag_3pi1 * TMath::Cos(thetaMax_tau_3pi1)/(2*(pow(E_3pi1, 2) - (pow(pmag_3pi1,2)*pow(TMath::Cos(thetaMax_tau_3pi1),2))));
  double pmag_tau2_rdInit = (pow(m_tauMass,2) + pow(m_3pi2, 2))* pmag_3pi2 * TMath::Cos(thetaMax_tau_3pi2)/(2*(pow(E_3pi2, 2) - (pow(pmag_3pi2,2)*pow(TMath::Cos(thetaMax_tau_3pi2),2))));

  Gaudi::XYZVector p_tau1_strat2 = pmag_tau1_rdInit * tau_dir1_strat2;
  Gaudi::XYZVector p_tau2_strat2 = pmag_tau2_rdInit * tau_dir2_strat2;
  
  //###

  //###
  if(m_doNuPZFix)
  {
    if(p_tau1_strat0.z() < p_3pi1.z())
    {
      p_tau1_strat0.SetZ(p_3pi1.z());
    }
    if(p_tau2_strat0.z() < p_3pi2.z())
    {
      p_tau2_strat0.SetZ(p_3pi2.z());
    }

    if(p_tau1_strat1.z() < p_3pi1.z())
    {
      p_tau1_strat1.SetZ(p_3pi1.z());
    }
    if(p_tau2_strat1.z() < p_3pi2.z())
    {
      p_tau2_strat1.SetZ(p_3pi2.z());
    }

    if(p_tau1_strat2.z() < p_3pi1.z())
    {
      p_tau1_strat2.SetZ(p_3pi1.z());
    }
    if(p_tau2_strat2.z() < p_3pi2.z())
    {
      p_tau2_strat2.SetZ(p_3pi2.z());
    }

    if(p_tau1_strat3.z() < p_3pi1.z())
    {
      p_tau1_strat3.SetZ(p_3pi1.z());
    }
    if(p_tau2_strat3.z() < p_3pi2.z())
    {
      p_tau2_strat3.SetZ(p_3pi2.z());
    }

    if(p_tau1_strat4.z() < p_3pi1.z())
    {
      p_tau1_strat4.SetZ(p_3pi1.z());
    }
    if(p_tau2_strat4.z() < p_3pi2.z())
    {
      p_tau2_strat4.SetZ(p_3pi2.z());
    }
    
  }

  //###

  //double mTau = m_tauMass;//1776.8199 is the value in truth MC. But PDG value is 1776.86 +- 0.12

  double E_tau1_strat0 = sqrt(pow(m_tauMass, 2) + p_tau1_strat0.Mag2());
  double E_tau2_strat0 = sqrt(pow(m_tauMass, 2) + p_tau2_strat0.Mag2());

  double E_tau1_strat1 = sqrt(pow(m_tauMass, 2) + p_tau1_strat1.Mag2());
  double E_tau2_strat1 = sqrt(pow(m_tauMass, 2) + p_tau2_strat1.Mag2());

  double E_tau1_strat2 = sqrt(pow(m_tauMass, 2) + p_tau1_strat2.Mag2());
  double E_tau2_strat2 = sqrt(pow(m_tauMass, 2) + p_tau2_strat2.Mag2());

  double E_tau1_strat3 = sqrt(pow(m_tauMass, 2) + p_tau1_strat3.Mag2());
  double E_tau2_strat3 = sqrt(pow(m_tauMass, 2) + p_tau2_strat3.Mag2());

  double E_tau1_strat4 = sqrt(pow(m_tauMass, 2) + p_tau1_strat4.Mag2());
  double E_tau2_strat4 = sqrt(pow(m_tauMass, 2) + p_tau2_strat4.Mag2());

  Gaudi::LorentzVector P4_tau1_strat0(p_tau1_strat0.x(), p_tau1_strat0.y(), p_tau1_strat0.z(), E_tau1_strat0);
  Gaudi::LorentzVector P4_tau2_strat0(p_tau2_strat0.x(), p_tau2_strat0.y(), p_tau2_strat0.z(), E_tau2_strat0);

  Gaudi::LorentzVector P4_tau1_strat1(p_tau1_strat1.x(), p_tau1_strat1.y(), p_tau1_strat1.z(), E_tau1_strat1);
  Gaudi::LorentzVector P4_tau2_strat1(p_tau2_strat1.x(), p_tau2_strat1.y(), p_tau2_strat1.z(), E_tau2_strat1);

  Gaudi::LorentzVector P4_tau1_strat2(p_tau1_strat2.x(), p_tau1_strat2.y(), p_tau1_strat2.z(), E_tau1_strat2);
  Gaudi::LorentzVector P4_tau2_strat2(p_tau2_strat2.x(), p_tau2_strat2.y(), p_tau2_strat2.z(), E_tau2_strat2);

  Gaudi::LorentzVector P4_tau1_strat3(p_tau1_strat3.x(), p_tau1_strat3.y(), p_tau1_strat3.z(), E_tau1_strat3);
  Gaudi::LorentzVector P4_tau2_strat3(p_tau2_strat3.x(), p_tau2_strat3.y(), p_tau2_strat3.z(), E_tau2_strat3);

  Gaudi::LorentzVector P4_tau1_strat4(p_tau1_strat4.x(), p_tau1_strat4.y(), p_tau1_strat4.z(), E_tau1_strat4);
  Gaudi::LorentzVector P4_tau2_strat4(p_tau2_strat4.x(), p_tau2_strat4.y(), p_tau2_strat4.z(), E_tau2_strat4);

  LHCb::Particle::Vector nuList;
  LHCb::Particle * tau1 = const_cast<LHCb::Particle*>(tauList[0]);
  LHCb::Particle * tau2 = const_cast<LHCb::Particle*>(tauList[1]);

  // Strategy 0
  nuList.push_back( new LHCb::Particle() );
  nuList.back()->setParticleID( LHCb::ParticleID(16) ); //Looks like the order of neutrino PID's doesn't matter. Check again later to make sure.
  nuList.back()->setMomentum( P4_tau1_strat0 - tauList[0]->momentum() );//NB 4 momenta going in here
  
  tau1->addToDaughters( nuList.back() );
  tau1->setMomentum( P4_tau1_strat0 );
  debug() << "Add nu 1 " << *nuList.back() << endmsg;
  debug() << "After " << *tau1 << endmsg;

  nuList.push_back( new LHCb::Particle() );
  nuList.back()->setParticleID( LHCb::ParticleID(-16) );
  nuList.back()->setMomentum( P4_tau2_strat0 - tauList[1]->momentum() );
  
  tau2->addToDaughters( nuList.back() );
  tau2->setMomentum( P4_tau2_strat0 );
  debug() << "Add nu 2 " << *nuList.back() << endmsg;
  debug() << "After " << *tau2 << endmsg;

  treeHead->setMomentum( P4_tau1_strat0 + P4_tau2_strat0 + Kplus->momentum() );

  TupleMap tMap ; // the temporary data map
  DecayTreeFitter::Fitter fitter_strat0(*treeHead, *originVtx[0], stateprovider );
  float chisq_strat_0 = fit(fitter_strat0, treeHead, originVtx[0], prefix, tMap, tuple, false);

  // Strategy 1
  //remove previously attached neutrino daughters
  tau1->removeFromDaughters(nuList.front());
  tau2->removeFromDaughters(nuList.back());

  //reset 4-momentum of 3pi to the original value
  tau1->setMomentum(p4_3pi1);
  tau2->setMomentum(p4_3pi2);

  nuList.front()->setMomentum( P4_tau1_strat1 - tauList[0]->momentum() );//NB 4 momenta going in here
  nuList.back()->setMomentum( P4_tau2_strat1 - tauList[1]->momentum() );//NB 4 momenta going in here

  tau1->addToDaughters( nuList.front() );
  tau1->setMomentum( P4_tau1_strat1 );

  tau2->addToDaughters( nuList.back() );
  tau2->setMomentum( P4_tau2_strat1 );

  treeHead->setMomentum( P4_tau1_strat1 + P4_tau2_strat1 + Kplus->momentum());

  DecayTreeFitter::Fitter fitter_strat1(*treeHead, *originVtx[0], stateprovider );
  float chisq_strat_1 = fit(fitter_strat1, treeHead, originVtx[0], prefix, tMap, tuple, false);

  // Strategy 2
  //remove previously attached neutrino daughters
  tau1->removeFromDaughters(nuList.front());
  tau2->removeFromDaughters(nuList.back());

  //reset 4-momentum of 3pi to the original value
  tau1->setMomentum(p4_3pi1);
  tau2->setMomentum(p4_3pi2);

  nuList.front()->setMomentum( P4_tau1_strat2 - tauList[0]->momentum() );//NB 4 momenta going in here
  nuList.back()->setMomentum( P4_tau2_strat2 - tauList[1]->momentum() );//NB 4 momenta going in here

  tau1->addToDaughters( nuList.front() );
  tau1->setMomentum( P4_tau1_strat2 );

  tau2->addToDaughters( nuList.back() );
  tau2->setMomentum( P4_tau2_strat2 );

  treeHead->setMomentum( P4_tau1_strat2 + P4_tau2_strat2 + Kplus->momentum() );

  DecayTreeFitter::Fitter fitter_strat2(*treeHead, *originVtx[0], stateprovider) ;
  float chisq_strat_2 = fit(fitter_strat2, treeHead, originVtx[0], prefix, tMap, tuple, false);

  int whichLowestChi2 = -1;

  if((fitter_strat0.status()==0) && (fitter_strat1.status()==0) && (fitter_strat2.status()==0))
  {// 0, 1 and 2 pass
    std::cout<<"All strategies passed "<<std::endl;
    if((chisq_strat_0 < chisq_strat_1) && (chisq_strat_0 < chisq_strat_2))
    {// strategy 0 has lowest chi2
      std::cout<<Form(" Strategy 0 has lowest #chi^{2} of %f ",chisq_strat_0)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat0 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat0 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat0 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat0 );

      treeHead->setMomentum( P4_tau1_strat0 + P4_tau2_strat0 + Kplus->momentum() );

      fit(fitter_strat0, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 0;
    }
    else if((chisq_strat_1 < chisq_strat_0) && (chisq_strat_1 < chisq_strat_2))
    {// strategy 1 has lowest chi2
      std::cout<<Form(" Strategy 1 has lowest #chi^{2} of %f",chisq_strat_1)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat1 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat1 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat1 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat1 );

      treeHead->setMomentum( P4_tau1_strat1 + P4_tau2_strat1 + Kplus->momentum() );

      fit(fitter_strat1, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 1;
    }
    else
    {// strategy 2 has lowest chi2
      std::cout<<Form(" Strategy 2 has lowest #chi^{2} of %f",chisq_strat_2)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat2 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat2 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat2 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat2 );

      treeHead->setMomentum( P4_tau1_strat2 + P4_tau2_strat2 + Kplus->momentum() );

      fit(fitter_strat2, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 2;
    }
  }

  else if((fitter_strat0.status()==0) && (fitter_strat1.status()==0) && (fitter_strat2.status()==1))
  {// 0 and 1 pass
    std::cout<<"Only strategies 0 and 1 pass"<<std::endl;
    if(chisq_strat_0 < chisq_strat_1)
    {// 0 < 1
      std::cout<<Form(" Strategy 0 has lowest #chi^{2} of %f ",chisq_strat_0)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat0 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat0 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat0 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat0 );

      treeHead->setMomentum( P4_tau1_strat0 + P4_tau2_strat0 + Kplus->momentum() );

      fit(fitter_strat0, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 0; 
    }
    else
    {// 1 < 0
      std::cout<<Form(" Strategy 1 has lowest #chi^{2} of %f ",chisq_strat_1)<<std::endl;
    
      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat1 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat1 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat1 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat1 );

      treeHead->setMomentum( P4_tau1_strat1 + P4_tau2_strat1 + Kplus->momentum() );

      fit(fitter_strat1, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 1; 
    }
  }

  else if((fitter_strat0.status()==0) && (fitter_strat1.status()==1) && (fitter_strat2.status()==0))
  {// 0 and 2 pass
    std::cout<<"Only strategies 0 and 2 pass"<<std::endl;
    if(chisq_strat_0 < chisq_strat_2)
    {// 0 < 2
      std::cout<<Form(" Strategy 0 has lowest #chi^{2} of %f ",chisq_strat_0)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat0 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat0 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat0 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat0 );

      treeHead->setMomentum( P4_tau1_strat0 + P4_tau2_strat0 + Kplus->momentum() );

      fit(fitter_strat0, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 0; 
    }
    else
    {// 2 < 0
      std::cout<<Form(" Strategy 2 has lowest #chi^{2} of %f ",chisq_strat_2)<<std::endl;
    
      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat2 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat2 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat2 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat2 );

      treeHead->setMomentum( P4_tau1_strat2 + P4_tau2_strat2 + Kplus->momentum() );

      fit(fitter_strat2, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 2; 
    }
  }

  else if((fitter_strat0.status()==1) && (fitter_strat1.status()==0) && (fitter_strat2.status()==0))
  {// 1 and 2 pass
    std::cout<<"Only strategies 1 and 2 pass"<<std::endl;
    if(chisq_strat_1 < chisq_strat_2)
    {// 1 < 2
      std::cout<<Form(" Strategy 1 has lowest #chi^{2} of %f ",chisq_strat_1)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat1 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat1 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat1 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat1 );

      treeHead->setMomentum( P4_tau1_strat1 + P4_tau2_strat1 + Kplus->momentum() );

      fit(fitter_strat1, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 1; 
    }
    else
    {// 2 < 1
      std::cout<<Form(" Strategy 2 has lowest #chi^{2} of %f ",chisq_strat_2)<<std::endl;

      //remove previously attached neutrino daughters
      tau1->removeFromDaughters(nuList.front());
      tau2->removeFromDaughters(nuList.back());

      //reset 4-momentum of 3pi to the original value
      tau1->setMomentum(p4_3pi1);
      tau2->setMomentum(p4_3pi2);

      nuList.front()->setMomentum( P4_tau1_strat2 - tauList[0]->momentum() );//NB 4 momenta going in here
      nuList.back()->setMomentum( P4_tau2_strat2 - tauList[1]->momentum() );//NB 4 momenta going in here

      tau1->addToDaughters( nuList.front() );
      tau1->setMomentum( P4_tau1_strat2 );

      tau2->addToDaughters( nuList.back() );
      tau2->setMomentum( P4_tau2_strat2 );

      treeHead->setMomentum( P4_tau1_strat2 + P4_tau2_strat2 + Kplus->momentum() );

      fit(fitter_strat2, treeHead, originVtx[0], prefix, tMap, tuple, true);
      whichLowestChi2 = 2; 
    }
  }

  else if((fitter_strat0.status()==0) && (fitter_strat1.status()==1) && (fitter_strat2.status()==1))
  {// only 0 passes
    std::cout<<"Only strategy 0 passes"<<std::endl;
    std::cout<<Form(" Strategy 0 has lowest #chi^{2} of %f ",chisq_strat_0)<<std::endl;

    //remove previously attached neutrino daughters
    tau1->removeFromDaughters(nuList.front());
    tau2->removeFromDaughters(nuList.back());

    //reset 4-momentum of 3pi to the original value
    tau1->setMomentum(p4_3pi1);
    tau2->setMomentum(p4_3pi2);

    nuList.front()->setMomentum( P4_tau1_strat0 - tauList[0]->momentum() );//NB 4 momenta going in here
    nuList.back()->setMomentum( P4_tau2_strat0 - tauList[1]->momentum() );//NB 4 momenta going in here

    tau1->addToDaughters( nuList.front() );
    tau1->setMomentum( P4_tau1_strat0 );

    tau2->addToDaughters( nuList.back() );
    tau2->setMomentum( P4_tau2_strat0 );

    treeHead->setMomentum( P4_tau1_strat0 + P4_tau2_strat0 + Kplus->momentum() );

    fit(fitter_strat0, treeHead, originVtx[0], prefix, tMap, tuple, true);
    whichLowestChi2 = 0; 
  }

  else if((fitter_strat0.status()==1) && (fitter_strat1.status()==0) && (fitter_strat2.status()==1))
  {// only 1 passes
    std::cout<<"Only strategy 1 passes"<<std::endl;
    std::cout<<Form(" Strategy 1 has lowest #chi^{2} of %f ",chisq_strat_1)<<std::endl;
  
    //remove previously attached neutrino daughters
    tau1->removeFromDaughters(nuList.front());
    tau2->removeFromDaughters(nuList.back());

    //reset 4-momentum of 3pi to the original value
    tau1->setMomentum(p4_3pi1);
    tau2->setMomentum(p4_3pi2);

    nuList.front()->setMomentum( P4_tau1_strat1 - tauList[0]->momentum() );//NB 4 momenta going in here
    nuList.back()->setMomentum( P4_tau2_strat1 - tauList[1]->momentum() );//NB 4 momenta going in here

    tau1->addToDaughters( nuList.front() );
    tau1->setMomentum( P4_tau1_strat1 );

    tau2->addToDaughters( nuList.back() );
    tau2->setMomentum( P4_tau2_strat1 );

    treeHead->setMomentum( P4_tau1_strat1 + P4_tau2_strat1 + Kplus->momentum() );

    fit(fitter_strat1, treeHead, originVtx[0], prefix, tMap, tuple, true);
    whichLowestChi2 = 1; 
  }

  else if((fitter_strat0.status()==1) && (fitter_strat1.status()==1) && (fitter_strat2.status()==0))
  {// only 2 passes
    std::cout<<"Only strategy 2 passes"<<std::endl;
    std::cout<<Form(" Strategy 2 has lowest #chi^{2} of %f ",chisq_strat_2)<<std::endl;

    //remove previously attached neutrino daughters
    tau1->removeFromDaughters(nuList.front());
    tau2->removeFromDaughters(nuList.back());

    //reset 4-momentum of 3pi to the original value
    tau1->setMomentum(p4_3pi1);
    tau2->setMomentum(p4_3pi2);

    nuList.front()->setMomentum( P4_tau1_strat2 - tauList[0]->momentum() );//NB 4 momenta going in here
    nuList.back()->setMomentum( P4_tau2_strat2 - tauList[1]->momentum() );//NB 4 momenta going in here

    tau1->addToDaughters( nuList.front() );
    tau1->setMomentum( P4_tau1_strat2 );

    tau2->addToDaughters( nuList.back() );
    tau2->setMomentum( P4_tau2_strat2 );

    treeHead->setMomentum( P4_tau1_strat2 + P4_tau2_strat2 + Kplus->momentum() );

    fit(fitter_strat2, treeHead, originVtx[0], prefix, tMap, tuple, true);
    whichLowestChi2 = 2;  
  }

  else
  {// no strategy passes
    std::cout<<"No strategy passes"<<std::endl;
    std::cout<<Form("Filling failed information of strategy 0 w/ #chi^{2} of %f ",chisq_strat_0)<<std::endl;

    //remove previously attached neutrino daughters
    tau1->removeFromDaughters(nuList.front());
    tau2->removeFromDaughters(nuList.back());

    //reset 4-momentum of 3pi to the original value
    tau1->setMomentum(p4_3pi1);
    tau2->setMomentum(p4_3pi2);

    nuList.front()->setMomentum( P4_tau1_strat0 - tauList[0]->momentum() );//NB 4 momenta going in here
    nuList.back()->setMomentum( P4_tau2_strat0 - tauList[1]->momentum() );//NB 4 momenta going in here

    tau1->addToDaughters( nuList.front() );
    tau1->setMomentum( P4_tau1_strat0 );

    tau2->addToDaughters( nuList.back() );
    tau2->setMomentum( P4_tau2_strat0 );

    treeHead->setMomentum( P4_tau1_strat0 + P4_tau2_strat0 + Kplus->momentum() );

    fit(fitter_strat0, treeHead, originVtx[0], prefix, tMap, tuple, true);
    whichLowestChi2 = -1; 
  }

  //std::cout<<"which fitter has lowest chi2 = "<<whichLowestChi2<<std::endl;
  tuple->column(prefix+"_0_Chi2", chisq_strat_0);
  tuple->column(prefix+"_1_Chi2", chisq_strat_1);
  tuple->column(prefix+"_2_Chi2", chisq_strat_2);
  tuple->column(prefix+"_whichLowChi2", whichLowestChi2);
  tuple->column(prefix+"_refPoint_X", refPoint_Kplus.X());
  tuple->column(prefix+"_refPoint_Y", refPoint_Kplus.Y());
  tuple->column(prefix+"_refPoint_Z", refPoint_Kplus.Z()); 

  return fillTuple(tMap,tuple,prefix); // the actual filling
  //return StatusCode(true);
}

//=============================================================================
// do filling for a given vertex
//=============================================================================
float TupleToolKTauTauDTFLowChi2::fit(DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
                                              const LHCb::VertexBase* pv, const std::string& prefix,
                                              TupleMap& tMap, Tuples::Tuple& tuple,
                                              bool fillMe) const
{
  if (msgLevel(MSG::VERBOSE)) verbose() << "fit " << P << " " << pv << " " << prefix << endmsg ;
  bool test = true ;
  //add mass contraints
  if ( !m_massConstraintsPids.empty() )
  {
    for ( const auto & C : m_massConstraintsPids )
    {
      fitter.setMassConstraint(C);
    }
  }
  // fit
  if (msgLevel(MSG::VERBOSE)) verbose() << "calling Fit" << endmsg ;

  std::vector<float> chisq_iters = {}; //DTF chi-square after every iteration
  std::vector<Gaudi::XYZPoint> BV = {}, DV1 = {}, DV2 = {}; //vertices after every iteration
  std::vector<float> B_M = {}; //B mass after every iteration
  std::vector<Gaudi::XYZVector> p_taup_nu = {}, p_taum_nu = {}; //neutrino kinematics after every iteration
  std::vector<Gaudi::XYZVector> p_taup = {}, p_taum = {}; //tau kinematics after every iteration

  //std::vector<float> BV_X = {}, BV_Y = {}, BV_Z = {};

  // const ParticleBase& pb = P;
  // int posindex = pb.posIndex();

  //calling the actual fit function.
  fitter.fit(m_maxNiter, 0.01, m_maxndiverging, m_dChisqQuit, chisq_iters, 
             BV, DV1, DV2, P, B_M,
             p_taup_nu, p_taum_nu,
             p_taup, p_taum, true);

  // std::cout<<"Initial p_taup_nu: "<<p_taup_nu[0].x()<<" "<<p_taup_nu[0].y()<<" "<<p_taup_nu[0].z()<<" p_taum_nu: "<<p_taum_nu[0].x()<<" "<<p_taum_nu[0].y()<<" "<<p_taum_nu[0].z()<<std::endl;

  // std::cout<<"Initial p_taup: "<<p_taup[0].x()<<" "<<p_taup[0].y()<<" "<<p_taup[0].z()<<" p_taum: "<<p_taum[0].x()<<" "<<p_taum[0].y()<<" "<<p_taum[0].z()<<std::endl;
  if (msgLevel(MSG::VERBOSE)) verbose() << "called Fit" << endmsg ;

  // std::cout<<"Initial BV :"<<BV[0].x()<<" "<<BV[0].y()<<" "<<BV[0].z()<<" "<<" DV1 :"<<DV1[0].x()<<" "<<DV1[0].y()<<" "<<DV1[0].z()<<" DV2 :"<<DV2[0].x()<<" "<<DV2[0].y()<<" "<<DV2[0].z()<<std::endl;

  // std::cout<<"Final BV :"<<BV[fitter.nIter()].x()<<" "<<BV[fitter.nIter()].y()<<" "<<BV[fitter.nIter()].z()<<" "<<" DV1 :"<<DV1[fitter.nIter()].x()<<" "<<DV1[fitter.nIter()].y()<<" "<<DV1[fitter.nIter()].z()<<" DV2 :"<<DV2[fitter.nIter()].x()<<" "<<DV2[fitter.nIter()].y()<<" "<<DV2[fitter.nIter()].z()<<std::endl;

  if(fillMe)
  {
    if(m_fillIterInfo)
    {
        //fill chi2 for each fitter iteration
        fillChi2Iter(fitter, chisq_iters, prefix, tuple); //made by AV 

        //fill decay vertices after each fitter iteration
        fillVtxIter(fitter, prefix, tuple, BV,
                    DV1, DV2, B_M, p_taup_nu, p_taum_nu,
                    p_taup, p_taum); //made by AV

    }

    // fill the fit result
    fillEndVtx(prefix+"_Bp", tMap, BV[fitter.nIter()]);
    fillEndVtx(prefix+"_Taup", tMap, DV1[fitter.nIter()]);
    fillEndVtx(prefix+"_Taum", tMap, DV2[fitter.nIter()]);

    fillDecay(fitter,prefix,tMap );
    fillMomentum(fitter,P,prefix,tMap );
    if (m_constrainToOriginVertex)
    {
      test &= fillPV(pv,prefix,tMap);
      test &= fillLT(fitter,P,prefix,tMap );
    }
    if ( isVerbose() )
    {
      test &= fillDaughters( fitter,P,prefix,tMap );
    }
    if ( m_updateDaughters )
    {
      test &= fillStableDaughters( fitter,P,prefix,tMap );
    }
}
  
  return chisq_iters.back();
}

//=============================================================================
// Fill endvertex information from DTF for intermediates
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillEndVtx(const std::string& prefix, TupleMap& tMap, Gaudi::XYZPoint vtx) const
{
  bool test = true;
  if (msgLevel(MSG::VERBOSE)) verbose() << "fillEndVtx " << prefix << endmsg ;

  test &= insert( prefix+"_ENDVERTEX_X", vtx.X(), tMap );
  test &= insert( prefix+"_ENDVERTEX_Y", vtx.Y(), tMap );
  test &= insert( prefix+"_ENDVERTEX_Z", vtx.Z(), tMap );
  return StatusCode(test);
}
//=============================================================================
// Fill standard stuff
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillPV(const LHCb::VertexBase* pv, const std::string& prefix,
                                       TupleMap& tMap ) const
{
  bool test = true;
  if (msgLevel(MSG::VERBOSE)) verbose() << "FillPV " << prefix << endmsg ;
  if (!pv) Exception("Null PVs cannot happen with ConstrainToOriginVertex!");
  test &= insert( prefix+"_PV_key", pv->key(), tMap );
  if ( isVerbose() )
  {
    test &= insert( prefix+"_PV_X", pv->position().X(), tMap );
    test &= insert( prefix+"_PV_Y", pv->position().Y(), tMap );
    test &= insert( prefix+"_PV_Z", pv->position().Z(), tMap );
  }
  return StatusCode(test);
}
//=============================================================================
// Fill chi2 for all iterations of the fitter
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillChi2Iter(DecayTreeFitter::Fitter& fitter, std::vector<float> chisq_iters,
                                             const std::string& prefix, Tuples::Tuple& tuple ) const
{
  bool test = true;
  if (msgLevel(MSG::VERBOSE)) verbose() << "fillChi2Iter " << prefix << endmsg ;

  test &= tuple->farray(prefix+"_chisq_iters", chisq_iters.begin(), chisq_iters.end(), 
                        prefix+"_nIters", m_maxNiter);
  
  // for(int i=1; i<=fitter.nIter(); i++)
  // {
  //   test &= insert( prefix+"_chisq_iters", chisq_iters[i-1], tMap  );
  // }
  
  return StatusCode(test);
}
//=============================================================================
// Fill standard stuff
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillDecay(const DecayTreeFitter::Fitter& fitter, const std::string& prefix,
                                          TupleMap& tMap ) const
{
  bool test = true;
  if (msgLevel(MSG::VERBOSE)) verbose() << "FillDecay " << prefix << endmsg ;

  test &= insert( prefix+"_status", fitter.status(), tMap );
  test &= insert( prefix+"_nDOF", fitter.nDof(), tMap  );
  test &= insert( prefix+"_chi2", fitter.chiSquare(), tMap  );
  test &= insert( prefix+"_nIter", fitter.nIter(), tMap  );

  return StatusCode(test);
}
//=============================================================================
// Fill B decay vertex after each fitter iteration
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillVtxIter(const DecayTreeFitter::Fitter& fitter, const std::string& prefix,
                                            Tuples::Tuple& tuple, std::vector<Gaudi::XYZPoint> BV,
                                            std::vector<Gaudi::XYZPoint> DV1, std::vector<Gaudi::XYZPoint> DV2,
                                            std::vector<float> B_M,
                                            std::vector<Gaudi::XYZVector> p_taup_nu, std::vector<Gaudi::XYZVector> p_taum_nu,
                                            std::vector<Gaudi::XYZVector> p_taup, std::vector<Gaudi::XYZVector> p_taum) const
{
  bool test = true;

  //if ( isVetoed(P->particleID().pid()) ) { return StatusCode(test); }

  if (msgLevel(MSG::VERBOSE)) verbose() << "fillVtxIter " << prefix << endmsg ;
  //Get the fit parameters
  // const auto    params = fitter.fitParams(P) ;
  // const auto& position = params->position() ;
  std::vector<float> BV_X = {}, BV_Y = {}, BV_Z = {};
  std::vector<float> DV1_X = {}, DV1_Y = {}, DV1_Z = {};
  std::vector<float> DV2_X = {}, DV2_Y = {}, DV2_Z = {};
  std::vector<float> p_taup_nu_X = {}, p_taup_nu_Y = {}, p_taup_nu_Z = {};
  std::vector<float> p_taum_nu_X = {}, p_taum_nu_Y = {}, p_taum_nu_Z = {};
  std::vector<float> p_taup_X = {}, p_taup_Y = {}, p_taup_Z = {};
  std::vector<float> p_taum_X = {}, p_taum_Y = {}, p_taum_Z = {};

  for (Gaudi::XYZPoint iter : BV)
  {
    BV_X.push_back(float(iter.X()));
    BV_Y.push_back(float(iter.Y()));
    BV_Z.push_back(float(iter.Z()));
  }  
  for (Gaudi::XYZPoint iter : DV1)
  {
    DV1_X.push_back(float(iter.X()));
    DV1_Y.push_back(float(iter.Y()));
    DV1_Z.push_back(float(iter.Z()));
  }  
  for (Gaudi::XYZPoint iter : DV2)
  {
    DV2_X.push_back(float(iter.X()));
    DV2_Y.push_back(float(iter.Y()));
    DV2_Z.push_back(float(iter.Z()));
  }
  for (Gaudi::XYZVector iter: p_taup_nu)
  {
    p_taup_nu_X.push_back(float(iter.X()));
    p_taup_nu_Y.push_back(float(iter.Y()));
    p_taup_nu_Z.push_back(float(iter.Z()));
  }  
  for (Gaudi::XYZVector iter: p_taum_nu)
  {
    p_taum_nu_X.push_back(float(iter.X()));
    p_taum_nu_Y.push_back(float(iter.Y()));
    p_taum_nu_Z.push_back(float(iter.Z()));
  }  
  for (Gaudi::XYZVector iter: p_taup)
  {
    p_taup_X.push_back(float(iter.X()));
    p_taup_Y.push_back(float(iter.Y()));
    p_taup_Z.push_back(float(iter.Z()));
  }  
  for (Gaudi::XYZVector iter: p_taum)
  {
    p_taum_X.push_back(float(iter.X()));
    p_taum_Y.push_back(float(iter.Y()));
    p_taum_Z.push_back(float(iter.Z()));
  }  
  // for (float iter : p_taup_nu_X)
  // {
  //   std::cout<<"In fillVtxIter "<<iter<<std::endl;
  // }
  test &= tuple->farray(prefix+"_Bp_ENDVERTEX_X", BV_X.begin(), BV_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Bp_ENDVERTEX_Y", BV_Y.begin(), BV_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Bp_ENDVERTEX_Z", BV_Z.begin(), BV_Z.end(),
                        prefix+"_nIters", m_maxNiter);

  test &= tuple->farray(prefix+"_Taup_ENDVERTEX_X", DV1_X.begin(), DV1_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_ENDVERTEX_Y", DV1_Y.begin(), DV1_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_ENDVERTEX_Z", DV1_Z.begin(), DV1_Z.end(),
                        prefix+"_nIters", m_maxNiter);

  test &= tuple->farray(prefix+"_Taum_ENDVERTEX_X", DV2_X.begin(), DV2_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_ENDVERTEX_Y", DV2_Y.begin(), DV2_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_ENDVERTEX_Z", DV2_Z.begin(), DV2_Z.end(),
                        prefix+"_nIters", m_maxNiter);

  test &= tuple->farray(prefix+"_Taup_nu_PX", p_taup_nu_X.begin(), p_taup_nu_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_nu_PY", p_taup_nu_Y.begin(), p_taup_nu_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_nu_PZ", p_taup_nu_Z.begin(), p_taup_nu_Z.end(),
                        prefix+"_nIters", m_maxNiter);
  
  test &= tuple->farray(prefix+"_Taum_nu_PX", p_taum_nu_X.begin(), p_taum_nu_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_nu_PY", p_taum_nu_Y.begin(), p_taum_nu_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_nu_PZ", p_taum_nu_Z.begin(), p_taum_nu_Z.end(),
                        prefix+"_nIters", m_maxNiter);

  test &= tuple->farray(prefix+"_Taup_PX", p_taup_X.begin(), p_taup_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_PY", p_taup_Y.begin(), p_taup_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_PZ", p_taup_Z.begin(), p_taup_Z.end(),
                        prefix+"_nIters", m_maxNiter);
  
  test &= tuple->farray(prefix+"_Taum_PX", p_taum_X.begin(), p_taum_X.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_PY", p_taum_Y.begin(), p_taum_Y.end(),
                        prefix+"_nIters", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_PZ", p_taum_Z.begin(), p_taum_Z.end(),
                        prefix+"_nIters", m_maxNiter);
  
  // test &= tuple->farray(prefix+"_M_iters", B_M.begin(), B_M.end(),
  //                       prefix+"_nIters", m_maxNiter);
  return StatusCode(test);
}
//=============================================================================
// Fill momentum and mass information
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillMomentum(const DecayTreeFitter::Fitter& fitter, const Particle* P,
                                             const std::string& prefix, TupleMap& tMap ) const
{
  bool test = true;

  if ( isVetoed(P->particleID().pid()) ) { return StatusCode(test); }

  if (msgLevel(MSG::VERBOSE)) verbose() << "FillMomentum " << prefix << endmsg ;
  //Get the fit parameters
  const auto    params = fitter.fitParams(P) ;
  const auto& momentum = params->momentum() ;
  if ( abs( P->particleID().pid()) == 15 ) {//just some checking printouts
    debug() << "Check pars for " << P->particleID().pid() << endmsg;
    debug() << "tau check M = " << momentum.M()  << ", " << momentum.E()*momentum.E() - momentum.Px()*momentum.Px() - momentum.Py()*momentum.Py() - momentum.Pz()*momentum.Pz() << endmsg;
    debug() << "Check p pars = " << momentum.Px() << ", " << momentum.Py() << ", " << momentum.Pz() << ", " << momentum.E() << endmsg;
  }
  test &= insert( prefix+"_M",  momentum.m().value(), tMap  );
  test &= insert( prefix+"_MERR", momentum.m().error(), tMap );
  test &= insert( prefix+"_P", momentum.p().value(), tMap );
  test &= insert( prefix+"_PERR", momentum.p().error(), tMap ) ;//MeV
  test &= insert( prefix+"_PX", momentum.Px(), tMap  );
  test &= insert( prefix+"_PY", momentum.Py(), tMap  );
  test &= insert( prefix+"_PZ", momentum.Pz(), tMap  );
  test &= insert( prefix+"_PE", momentum.E() , tMap  );//MeV
  return StatusCode(test);
}
//=============================================================================
// Fill lifetime information
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillLT(const DecayTreeFitter::Fitter& fitter, const Particle* P,
                                       const std::string& prefix, TupleMap& tMap ) const
{
  bool test = true;

  if ( isVetoed(P->particleID().pid()) ) { return StatusCode(test); }

  if (msgLevel(MSG::VERBOSE)) verbose() << "FillLT " << prefix << endmsg ;
  const auto  tParams     = fitter.fitParams(P);
  const auto& decayLength = tParams->decayLength();
  const auto& ctau        = tParams->ctau();
  test &= insert( prefix+"_ctau", ctau.value(), tMap  );
  test &= insert( prefix+"_ctauErr", ctau.error(), tMap  );
  test &= insert( prefix+"_decayLength", decayLength.value(), tMap  );
  test &= insert( prefix+"_decayLengthErr", decayLength.error(), tMap  );

  return StatusCode(test);
}
//=============================================================================
// Fill lifetime information for non stable daughters
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillDaughters(const DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
                                              const std::string& prefix, TupleMap& tMap ) const
{
  bool test = true;

  if (msgLevel(MSG::VERBOSE)) verbose() << "FillDaughters " << prefix << endmsg ;
  const auto & daughters = ( m_useFullTreeInName ?
                             P->daughtersVector() :
                             m_particleDescendants->descendants(P) );
  if (msgLevel(MSG::DEBUG)) debug() << "for id " << P->particleID().pid()
                                    << " daughter size is " << daughters.size() << endmsg;
  if ( daughters.empty() ) return StatusCode(test);
  std::set<std::string> usedNames;
  unsigned int add = 0;
  for ( const auto& particle : daughters )
  {
    if ( particle->isBasicParticle() ) continue ;
    auto pid = abs(particle->particleID().pid());
    if ( std::find( m_v_ids.begin(), m_v_ids.end(), pid) != m_v_ids.end() ) {
      pid = m_v_ids[0];
    }
    const auto pidName = getName(pid) ;
    auto name = prefix+"_"+pidName ;
    bool renamed = false;
    while ( usedNames.find(name) != usedNames.end() )
    { // fix to bug 88702
      renamed = true;
      if (msgLevel(MSG::VERBOSE)) verbose() << "Found already name " << name
                                            << " trying next " << endmsg;
      name = prefix + "_" + pidName + "_" + boost::lexical_cast<std::string>(add);
      ++add;
    }
    if ( renamed ) Info("Renaming duplicate to "+name,StatusCode::SUCCESS,1);
    usedNames.insert(name);
    test &= fillMomentum( fitter,particle,name,tMap );
    test &= fillLT( fitter,particle,name,tMap);
    if ( m_useFullTreeInName ) { fillDaughters(fitter,particle,name,tMap); }
  }
  return StatusCode(test);
}
//=============================================================================
// Fill lifetime information for stable daughters
//=============================================================================
StatusCode
TupleToolKTauTauDTFLowChi2::fillStableDaughters(const DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
                                         const std::string& prefix, TupleMap& tMap ) const
{
  bool test = true;

  if (msgLevel(MSG::VERBOSE)) verbose() << "FillStableDaughters " << prefix << endmsg ;
  const LHCb::Particle::ConstVector& daughters = P->daughtersVector();
  if (msgLevel(MSG::DEBUG)) debug() << "for id " << P->particleID().pid()
                                    << " daughter size is " << daughters.size() << endmsg;
  if ( daughters.empty() ) return StatusCode(test);
  std::set<std::string> usedNames;
  unsigned int add = 0;
  for ( const auto& particle : daughters )
  {
    if ( !particle->isBasicParticle() )
    {
      const auto pid = abs(particle->particleID().pid());
      const auto pidName = getName(pid) ;
      auto name = prefix+"_"+pidName ;
      bool renamed = false;
      while ( usedNames.find(name) != usedNames.end() )
      { // fix to bug 88702
        if (msgLevel(MSG::VERBOSE)) verbose() << "Found already name " << name
                                              << " trying next " << endmsg ;
        renamed = true;
        name = prefix+"_"+pidName+"_"+boost::lexical_cast<std::string>(add);
        ++add;
      }
      if ( renamed ) Info("Renaming duplicate to "+name,StatusCode::SUCCESS,1);
      usedNames.insert(name);
      test &= fillStableDaughters( fitter, particle, name, tMap );
    }
    else
    {
      // const int pid = particle->particleID().pid();
      auto pid = abs(particle->particleID().pid());
      if ( std::find( m_v_ids.begin(), m_v_ids.end(), pid) != m_v_ids.end() ) {
        pid = m_v_ids[0];
      }
      const auto pidName = getName(pid) ;
      auto name = prefix+"_"+pidName;
      bool renamed = false;
      while ( usedNames.find(name) != usedNames.end() )
      { // fix to bug 88702
        if (msgLevel(MSG::VERBOSE)) verbose() << "Found already name " << name
                                              << " trying next " << endmsg ;
        renamed = true;
        name = prefix+"_"+pidName+"_"+boost::lexical_cast<std::string>(add);
        ++add;
      }
      if ( renamed ) Info("Renaming duplicate to "+name,StatusCode::SUCCESS,1);
      usedNames.insert(name);
      test &= fillTracksMomentum( fitter, particle, name, tMap );
    }
  }
  return StatusCode(test);
}

//=============================================================================
// Fill updated tracks momentum
//=============================================================================
StatusCode
TupleToolKTauTauDTFLowChi2::fillTracksMomentum(const DecayTreeFitter::Fitter& fitter, const Particle* P,
                                        const std::string& prefix, TupleMap& tMap ) const
{
  bool test = true;

  if ( isVetoed(P->particleID().pid()) ) { return StatusCode(test); }

  if (msgLevel(MSG::VERBOSE)) verbose() << "FillTracksMomentum " << prefix << endmsg ;

  // Get the fit parameters
  const auto    params = fitter.fitParams(P) ;
  const auto& momentum = params->momentum() ;

  test &= insert( prefix+"_ID", P->particleID().pid(), tMap  );
  test &= insert( prefix+"_PX", momentum.Px(), tMap  );
  test &= insert( prefix+"_PY", momentum.Py(), tMap  );
  test &= insert( prefix+"_PZ", momentum.Pz(), tMap  );
  test &= insert( prefix+"_PE", momentum.E() , tMap  );//MeV

  return StatusCode(test);
}

//=============================================================================
// append data to TupleMap
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::insert(const std::string& leaf, const double val,
                                       TupleMap& tMap ) const
{
  auto l = tMap.find(leaf);
  if ( l == tMap.end() )
  {  /// first time this is seen. Create
    std::vector<double> vals;
    vals.push_back(val);
    tMap.insert( std::make_pair(leaf,vals) );
  }
  else
  {
    l->second.push_back(val); /// append a to vector
  }
  if (msgLevel(MSG::VERBOSE))
    verbose() << "insert " << leaf << " " << val
              << " size " << l->second.size() << endmsg ;
  return StatusCode::SUCCESS ;
}

//=============================================================================
// actual filling of the Tuple
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::fillTuple(TupleMap& tMap, Tuples::Tuple& tuple,
                                          const std::string& prefix )
{
  bool test = true ;

  if ( UNLIKELY(m_firstTupleFill) )
  {
    // Save the list of keys in the given order for future comparisons
    m_firstTupleKeys.clear();
    m_firstTupleKeys.reserve( tMap.size() );
    for ( const auto& i : tMap ) { m_firstTupleKeys.emplace_back(i.first); }
    // flag having saved the keys
    m_firstTupleFill = false;
  }
  else
  {
    // test against the first set of keys
    test = checkTupleKeys( tMap );
  }

  // if OK, save and continue
  if ( test )
  {
    for ( const auto& t : tMap )
    {
      const auto& leaf = t.first;
      const auto& data = t.second;
      if (msgLevel(MSG::DEBUG))
        debug() << "Filling leaf ``" << leaf << "'' with vector of size "
                << data.size() << endmsg ;
      if ( m_maxPV < data.size() )
        Exception("Seeing data with too many PVs. Have you set MaxPVs?");
      test &= tuple->farray( leaf, data, prefix+"_nPV", m_maxPV);
    }
  }

  return StatusCode(test);
}

//=============================================================================
// Sort Tracks
//=============================================================================
std::set<const LHCb::Track*>
TupleToolKTauTauDTFLowChi2::sortedTracks(const LHCb::VertexBase* vb) const
{
  const LHCb::RecVertex* pv = dynamic_cast<const LHCb::RecVertex*>(vb);
  if (!pv) Exception("Failed to cast PV");
  std::set<const LHCb::Track*> st ;
  for ( const auto& t : pv->tracks() ) { st.insert(t); }
  return st ;
}

//=============================================================================
// Compare PVs, check that one PV's tracks is a subset of the other
//=============================================================================
bool TupleToolKTauTauDTFLowChi2::samePV( const LHCb::VertexBase* vb1,
                                  const LHCb::VertexBase* vb2 ) const
{
  // exception checking. See bug https://savannah.cern.ch/bugs/?100933
  if ( !vb1 && !vb2 )
  {
    Warning("samePV method called with 2 NULL PVs. "
            "The answer is obviously true, but you may want to check the meaning of the question.",
            StatusCode::SUCCESS,1);
    return true ;
  }
  else if ( !vb1 || !vb2 )
  {
    Warning("samePV method called with 1 NULL PV. "
            "The answer is obviously false, but you may want to check the meaning of the question.",
            StatusCode::SUCCESS,1);
    return false ;
  }

  if ( !(vb1->isPrimary()) || !(vb2->isPrimary()) )
  {
    Warning("Non PV VertexBase is being used as PV", StatusCode::SUCCESS, 1).ignore();
    return false ;
  }

  const auto st1 = sortedTracks(vb1);
  const auto st2 = sortedTracks(vb2);

  const bool inc = std::includes(st1.begin(),st1.end(),st2.begin(),st2.end());
  if ( msgLevel(MSG::VERBOSE))
  {
    verbose() << "PV 2 of size " << st2.size() << " is ";
    if (!inc) verbose() << "not ";
    verbose() << "included in PV 1 of size " << st1.size() << endmsg ;
  }
  return inc;
}

//=============================================================================
// get origin vertex
//=============================================================================
std::vector<const VertexBase*>
TupleToolKTauTauDTFLowChi2::originVertex( const Particle* mother, const Particle* P ) const
{
  std::vector<const VertexBase*> oriVx;
  if ( mother == P )
  {// the origin vertex is the primary.
    const auto bpv = m_dva->bestVertex( P );
    if ( bpv )
    {
      oriVx.push_back(bpv);
      if ( msgLevel(MSG::VERBOSE) )
        verbose() << "Pushed back bpv " << bpv << " from "
                  << tesLocation(bpv) << " at "
                  << bpv->position() << endmsg ;
    }
    else if ( m_constrainToOriginVertex)
    {
      Warning( "NULL bestPV while constraining to origin vertex. Fit will be ignored.",
               StatusCode::SUCCESS, 0 ).ignore();
    }

    // Remove code to use other than best PV for now
    // // all the other ones
    // /// @todo : keep only the related ones
    // for ( const auto & pv : m_dva->primaryVertices() )
    // {
    //   if ( m_storeAnyway || !samePV(pv,bpv) )
    //   {
    //     oriVx.push_back(pv);
    //     if ( msgLevel(MSG::VERBOSE) )
    //       verbose() << "Pushed back  pv " << pv << " from "
    //                 << tesLocation(pv) << " at "
    //                 << pv->position() << endmsg ;
    //   }
    //   if ( oriVx.size() >= m_maxPV )
    //   {
    //     Warning("Truncated number of PVs", StatusCode::FAILURE, 0).ignore();
    //     break ;
    //   }
    // }
  }
  else
  {
    const auto & dau = mother->daughters ();
    if ( dau.empty() ) return oriVx ;

    for ( const auto& d : dau )
    {
      if ( P == d )
      {
        oriVx.push_back(mother->endVertex());
        return oriVx ;
      }
    }

    // vertex not yet found, get deeper in the decay:
    for ( const auto& d : dau )
    {
      if ( P != d && !d->isBasicParticle() )
      {
        oriVx = originVertex( d, P );
        if( !oriVx.empty() )
        {
          return oriVx ;  // found
        }
      }
    }
  }
  return oriVx;
}

//=============================================================================
// Convert pid number in names
//=============================================================================
std::string TupleToolKTauTauDTFLowChi2::getName(const int id) const
{
  const auto * prop = m_ppSvc->find( LHCb::ParticleID(id) );
  if (!prop) Exception("Unknown PID");
  //if (msgLevel(MSG::VERBOSE)) verbose() << "ID " << id << " gets name "
  //                                      << Decays::escape(prop->name()) << endmsg ;
  return Decays::escape(prop->name());
}

//=============================================================================
// Substitute
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::substitute(LHCb::DecayTree& tree)
{
  if (msgLevel(MSG::DEBUG)) debug() << "Calling substitute" << endmsg ;
  const auto substituted = m_substitute->substitute ( tree.head() ) ;
  // debugging
  if ( msgLevel(MSG::VERBOSE) || 0 == substituted )
  {
    const auto mp = tree.cloneMap();
    for ( const auto & i : mp )
    {
      if ( i.first->particleID().pid() == i.second->particleID().pid() )
      {
        info() << "A " << getName(i.first->particleID().pid()) << " remains unchanged" << endmsg ;
      }
      else
      {
        info() << "A " << getName(i.first->particleID().pid()) << " is substituted by a "
               << getName(i.second->particleID().pid()) << endmsg ;
      }
    }

  }
  if ( 0 == substituted )
  {
    return Error( "No particles have been substituted. Check your substitution options." );
  }
  return StatusCode::SUCCESS ;
}

//=============================================================================
// Check Mass Constraints
//=============================================================================
StatusCode TupleToolKTauTauDTFLowChi2::checkMassConstraints(const LHCb::DecayTree& tree)
{
  if (!m_first) return StatusCode::SUCCESS ;  // do that only once
  m_first = false ;
  const auto mp = tree.cloneMap();
  for ( const auto & m : m_massConstraintsPids )
  {
    bool found = false ;
    for ( const auto & i : mp )
    {
      if ( m.abspid() == i.second->particleID().abspid() )
      {
        found = true ;
        break ;
      }
    }
    if ( found &&  msgLevel(MSG::VERBOSE) )
      verbose() << "Constraint " << getName(m.pid()) << " was found in tree" << endmsg ;
    if ( !found )
    {
      std::ostringstream mess;
      mess <<  "Constraint " << getName(m.pid())
           << " was not found in tree. Check your options. Maybe also the substitution options.";
      return Error( mess.str() ) ;
    }
  }
  return StatusCode::SUCCESS ;

}

Gaudi::XYZVector TupleToolKTauTauDTFLowChi2::makeTransformation_vec(Gaudi::XYZVector p_K, Gaudi::XYZPoint refPoint, 
                                                             Gaudi::XYZVector theVector, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring p_K into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(p_K.x()/p_K.y());
  //Clockwise rotation about new X axis to align Z axis with p_K
  double theta = -1 * TMath::ATan(p_K.Rho() * (p_K.y()/TMath::Abs(p_K.y())) /p_K.z());

  Gaudi::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  Gaudi::Transform3D myTransformation     = Gaudi::Transform3D(myRotation) * Gaudi::Transform3D(myShift);
  Gaudi::Transform3D myTransformation_inv = myTransformation.Inverse();

  Gaudi::XYZVector theVector_t;
  if(invFlag)
    theVector_t = myTransformation_inv * theVector;
  else 
    theVector_t = myTransformation * theVector;
  
  return theVector_t;
}

Gaudi::XYZPoint TupleToolKTauTauDTFLowChi2::makeTransformation_point(Gaudi::XYZVector p_K, Gaudi::XYZPoint refPoint, 
                                                              Gaudi::XYZPoint thePoint, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring p_K into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(p_K.x()/p_K.y());
  //Clockwise rotation about new X axis to align Z axis with p_K
  double theta = -1 * TMath::ATan(p_K.Rho() * (p_K.y()/TMath::Abs(p_K.y())) /p_K.z());

  Gaudi::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  Gaudi::Transform3D myTransformation     = Gaudi::Transform3D(myRotation) * Gaudi::Transform3D(myShift);
  Gaudi::Transform3D myTransformation_inv = myTransformation.Inverse();

  Gaudi::XYZPoint thePoint_t;
  if(invFlag)
    thePoint_t = myTransformation_inv * thePoint;
  else 
    thePoint_t = myTransformation * thePoint;
  
  return thePoint_t;
}
// Declaration of the Tool Factory
DECLARE_COMPONENT( TupleToolKTauTauDTFLowChi2 )