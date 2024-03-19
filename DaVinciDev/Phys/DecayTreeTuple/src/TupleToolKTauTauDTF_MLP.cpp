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
#include "TupleToolKTauTauDTF_MLP.h"
#include "TMath.h"
#include <algorithm>

#include "Kernel/IVertexFit.h"
#include "TrackInterfaces/ITrackFitter.h"
#include "Event/ODIN.h" // event & run number
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>

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
TupleToolKTauTauDTF_MLP::TupleToolKTauTauDTF_MLP( const std::string& type,
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
  declareProperty( "tauMass", m_tauMass = 1776.86, "tau mass used in initialization calculation. 1776.8199 MeV is the number used in LHCb MC, 1776.86 is the PDG average");
  declareProperty( "kaonMass", m_kaonMass = 493.677, "tau mass used in initialization calculation. 1776.8199 MeV is the number used in LHCb MC, 1776.86 is the PDG average");

  declareInterface<IParticleTupleTool>(this);
}

//=============================================================================
StatusCode TupleToolKTauTauDTF_MLP::initialize()
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
    if ( "TupleToolKTauTauDTF_MLP" == m_extraName )  m_extraName = ""; // user has not chanegd instance name
    info() << "All fields will be prepended with ``" << m_extraName << "''" <<endmsg;
  }

  if ( m_extraName.empty() )
  {
    return Error( "Extraname is empty. Always give an instance name "
                  "to TupleToolKTauTauDTF_MLP! See doxygen." );
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

StatusCode TupleToolKTauTauDTF_MLP::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  if ( !m_stateprovider.empty() ) { sc = m_stateprovider.release(); }
  return StatusCode{ TupleToolBase::finalize() && sc };
}

//=============================================================================
//  The fill method implementation. This is where our modifications are happening
//=============================================================================
StatusCode TupleToolKTauTauDTF_MLP::fill( const LHCb::Particle* mother, const LHCb::Particle* P, const std::string& head, Tuples::Tuple& tuple )
{

//  bool runningOnMC = true;

  if( !P ) return StatusCode::FAILURE;
  if ( P->isBasicParticle() )
  {
    return Error("Do not call TupleToolKTauTauDTF_MLP for basic particles. Use Branches. See doxygen.");
  }
  const std::string prefix = fullName(head);
  if (msgLevel(MSG::DEBUG)) debug() << "head ''" << head << "'' prefix ''" << prefix
                                    << "'' extraname ''" << m_extraName << "''" <<endmsg;

  const auto stateprovider = ( m_stateprovider.empty() ? nullptr : &(*m_stateprovider) );

  LHCb::DecayTree tree ( *P ) ;
  LHCb::Particle * treeHead = tree.head();
  
  // Get PV
  std::vector<const VertexBase*> originVtx;
  originVtx = originVertex(mother, P);
  if( originVtx.empty() ){return Error("Can't get an origin vertex");}
  ROOT::Math::XYZPoint PV( originVtx[0]->position().X(), originVtx[0]->position().Y(), originVtx[0]->position().Z() ); // originVtx[0] is the best PV in the event

  // Get taus and their daughters
  LHCb::Particle::ConstVector tauList;
  if(mother->charge() > 0) // B+ -> tau+ tau- K+
  {
      LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList); // tau+
      LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList); // tau-
  }
  else // B- -> tau- tau+ K-
  {
      LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList); // tau-
      LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList); // tau+
  }
  if ( tauList.size() != 2 ) 
  {
      return Error("Can't find two tau decays");
  }

  LHCb::Particle::ConstVector taup_daus = tauList[0]->daughtersVector();
  LHCb::Particle::ConstVector taum_daus = tauList[1]->daughtersVector();

  // Get K+
  const LHCb::Particle* Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(321) );
  if( Kp == 0 )
  {
      Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-321) );
  }
  if( Kp == 0 )
  {
      return Error("Can't find a K meson in decay tree");
  }

  // Refit tau+ decay vertex
  const IVertexFit* vtxFitter_taup = tool<IVertexFit>( "LoKi::VertexFitter" );
  LHCb::Vertex taup_refittedVertex;
  LHCb::Particle taup_refitted;
  taup_refitted.setParticleID(tauList[0]->particleID());

  vtxFitter_taup->fit(taup_daus, taup_refittedVertex, taup_refitted);

  // Refit tau- decay vertex
  const IVertexFit* vtxFitter_taum = tool<IVertexFit>( "LoKi::VertexFitter" );
  LHCb::Vertex taum_refittedVertex;
  LHCb::Particle taum_refitted;
  taum_refitted.setParticleID(tauList[1]->particleID());

  vtxFitter_taum->fit(taum_daus, taum_refittedVertex, taum_refitted);

  // Get tau+ refitted decay vertex and visible 4-momentum (3pi)+ 
  ROOT::Math::XYZPoint DV1( taup_refitted.endVertex()->position().X(), taup_refitted.endVertex()->position().Y(), taup_refitted.endVertex()->position().Z() );
  Gaudi::LorentzVector P3pi1( taup_refitted.momentum().X(), taup_refitted.momentum().Y(), taup_refitted.momentum().Z(), taup_refitted.momentum().E() );
  ROOT::Math::XYZVector p3pi1( P3pi1.x(), P3pi1.y(), P3pi1.z() );  

  // Get tau- refitted decay vertex and visible 4-momentum (3pi)-
  ROOT::Math::XYZPoint DV2( taum_refitted.endVertex()->position().X(), taum_refitted.endVertex()->position().Y(), taum_refitted.endVertex()->position().Z() );
  Gaudi::LorentzVector P3pi2( taum_refitted.momentum().X(), taum_refitted.momentum().Y(), taum_refitted.momentum().Z(), taum_refitted.momentum().E() );
  ROOT::Math::XYZVector p3pi2( P3pi2.x(), P3pi2.y(), P3pi2.z() );

  // Get K+ reference point and 4-momentum
  ROOT::Math::XYZPoint refPoint_Kp = Kp->referencePoint();
  ROOT::Math::XYZVector pK( Kp->momentum().X(), Kp->momentum().Y(), Kp->momentum().Z() );
  Double_t EK = sqrt( pow(m_kaonMass,2) + pK.Mag2() );

  int dimM = 22;
  CLHEP::HepVector m(dimM);
  m(1) = PV.x(); // the indexing of CLHEP vectors starts with 1
  m(2) = PV.y();
  m(3) = PV.z();
  m(4) = DV1.x();
  m(5) = DV1.y();
  m(6) = DV1.z();
  m(7) = P3pi1.x();
  m(8) = P3pi1.y();
  m(9) = P3pi1.z();
  m(10) = P3pi1.t();
  m(11) = DV2.x();
  m(12) = DV2.y();
  m(13) = DV2.z();
  m(14) = P3pi2.x();
  m(15) = P3pi2.y();
  m(16) = P3pi2.z();
  m(17) = P3pi2.t();
  m(18) = refPoint_Kp.x();
  m(19) = refPoint_Kp.y();
  m(20) = pK.x();
  m(21) = pK.y();
  m(22) = pK.z();

  // Event number
  // Load the ODIN
  const LHCb::ODIN* odin = getIfExists<ODIN>( evtSvc(), LHCb::ODINLocation::Default );
  if ( !odin ) { odin = getIfExists<ODIN>( evtSvc(), LHCb::ODINLocation::Default, false ); }
  if ( !odin ) {
    // should always be available ...
    return Error( "Cannot load the ODIN data object", StatusCode::SUCCESS );
  }
  ULong64_t eventNumber = odin->eventNumber();

  // MLP initialisation
  TMVA::Tools::Instance();
  TMVA::Reader *reader_taup_PX = new TMVA::Reader( "V:Color:Silent" ); 
  TMVA::Reader *reader_taup_PY = new TMVA::Reader( "V:Color:Silent" ); 
  TMVA::Reader *reader_taup_PZ = new TMVA::Reader( "V:Color:Silent" ); 
  TMVA::Reader *reader_taum_PX = new TMVA::Reader( "V:Color:Silent" ); 
  TMVA::Reader *reader_taum_PY = new TMVA::Reader( "V:Color:Silent" ); 
  TMVA::Reader *reader_taum_PZ = new TMVA::Reader( "V:Color:Silent" ); 
  
  TString weightfile_taup_PX;
  TString weightfile_taup_PY;
  TString weightfile_taup_PZ;
  TString weightfile_taum_PX;
  TString weightfile_taum_PY;
  TString weightfile_taum_PZ;

  // Ktautau weights
  // weightfile_taup_PX = "TMVARegression_taup_TRUEP_X_MLP.weights.xml";
  // weightfile_taup_PY = "TMVARegression_taup_TRUEP_Y_MLP.weights.xml";
  // weightfile_taup_PZ = "TMVARegression_taup_TRUEP_Z_MLP.weights.xml";
  // weightfile_taum_PX = "TMVARegression_taum_TRUEP_X_MLP.weights.xml";
  // weightfile_taum_PY = "TMVARegression_taum_TRUEP_Y_MLP.weights.xml";
  // weightfile_taum_PZ = "TMVARegression_taum_TRUEP_Z_MLP.weights.xml";
  weightfile_taup_PX = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP.weights.xml";
  weightfile_taup_PY = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP.weights.xml";
  weightfile_taup_PZ = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP.weights.xml";
  weightfile_taum_PX = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP.weights.xml";
  weightfile_taum_PY = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP.weights.xml";
  weightfile_taum_PZ = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP.weights.xml";

  Float_t MLP_m_vars[dimM];
  Float_t MLP_RPz;
  Float_t MLP_eventNumber;

  reader_taup_PX->AddVariable("df_m_1", &MLP_m_vars[0]);
  reader_taup_PX->AddVariable("df_m_2", &MLP_m_vars[1]);
  reader_taup_PX->AddVariable("df_m_3", &MLP_m_vars[2]);
  reader_taup_PX->AddVariable("df_m_4", &MLP_m_vars[3]);
  reader_taup_PX->AddVariable("df_m_5", &MLP_m_vars[4]);
  reader_taup_PX->AddVariable("df_m_6", &MLP_m_vars[5]);
  reader_taup_PX->AddVariable("df_m_7", &MLP_m_vars[6]);
  reader_taup_PX->AddVariable("df_m_8", &MLP_m_vars[7]);
  reader_taup_PX->AddVariable("df_m_9", &MLP_m_vars[8]);
  reader_taup_PX->AddVariable("df_m_10", &MLP_m_vars[9]);
  reader_taup_PX->AddVariable("df_m_11", &MLP_m_vars[10]);
  reader_taup_PX->AddVariable("df_m_12", &MLP_m_vars[11]);
  reader_taup_PX->AddVariable("df_m_13", &MLP_m_vars[12]);
  reader_taup_PX->AddVariable("df_m_14", &MLP_m_vars[13]);
  reader_taup_PX->AddVariable("df_m_15", &MLP_m_vars[14]);
  reader_taup_PX->AddVariable("df_m_16", &MLP_m_vars[15]);
  reader_taup_PX->AddVariable("df_m_17", &MLP_m_vars[16]);
  reader_taup_PX->AddVariable("df_m_18", &MLP_m_vars[17]);
  reader_taup_PX->AddVariable("df_m_19", &MLP_m_vars[18]);
  reader_taup_PX->AddVariable("df_m_20", &MLP_m_vars[19]);
  reader_taup_PX->AddVariable("df_m_21", &MLP_m_vars[20]);
  reader_taup_PX->AddVariable("df_m_22", &MLP_m_vars[21]);
  reader_taup_PX->AddVariable("Kp_RP_Z", &MLP_RPz);
  reader_taup_PX->AddSpectator("eventNumber", &MLP_eventNumber);

  reader_taup_PY->AddVariable("df_m_1", &MLP_m_vars[0]);
  reader_taup_PY->AddVariable("df_m_2", &MLP_m_vars[1]);
  reader_taup_PY->AddVariable("df_m_3", &MLP_m_vars[2]);
  reader_taup_PY->AddVariable("df_m_4", &MLP_m_vars[3]);
  reader_taup_PY->AddVariable("df_m_5", &MLP_m_vars[4]);
  reader_taup_PY->AddVariable("df_m_6", &MLP_m_vars[5]);
  reader_taup_PY->AddVariable("df_m_7", &MLP_m_vars[6]);
  reader_taup_PY->AddVariable("df_m_8", &MLP_m_vars[7]);
  reader_taup_PY->AddVariable("df_m_9", &MLP_m_vars[8]);
  reader_taup_PY->AddVariable("df_m_10", &MLP_m_vars[9]);
  reader_taup_PY->AddVariable("df_m_11", &MLP_m_vars[10]);
  reader_taup_PY->AddVariable("df_m_12", &MLP_m_vars[11]);
  reader_taup_PY->AddVariable("df_m_13", &MLP_m_vars[12]);
  reader_taup_PY->AddVariable("df_m_14", &MLP_m_vars[13]);
  reader_taup_PY->AddVariable("df_m_15", &MLP_m_vars[14]);
  reader_taup_PY->AddVariable("df_m_16", &MLP_m_vars[15]);
  reader_taup_PY->AddVariable("df_m_17", &MLP_m_vars[16]);
  reader_taup_PY->AddVariable("df_m_18", &MLP_m_vars[17]);
  reader_taup_PY->AddVariable("df_m_19", &MLP_m_vars[18]);
  reader_taup_PY->AddVariable("df_m_20", &MLP_m_vars[19]);
  reader_taup_PY->AddVariable("df_m_21", &MLP_m_vars[20]);
  reader_taup_PY->AddVariable("df_m_22", &MLP_m_vars[21]);
  reader_taup_PY->AddVariable("Kp_RP_Z", &MLP_RPz);
  reader_taup_PY->AddSpectator("eventNumber", &MLP_eventNumber);

  reader_taup_PZ->AddVariable("df_m_1", &MLP_m_vars[0]);
  reader_taup_PZ->AddVariable("df_m_2", &MLP_m_vars[1]);
  reader_taup_PZ->AddVariable("df_m_3", &MLP_m_vars[2]);
  reader_taup_PZ->AddVariable("df_m_4", &MLP_m_vars[3]);
  reader_taup_PZ->AddVariable("df_m_5", &MLP_m_vars[4]);
  reader_taup_PZ->AddVariable("df_m_6", &MLP_m_vars[5]);
  reader_taup_PZ->AddVariable("df_m_7", &MLP_m_vars[6]);
  reader_taup_PZ->AddVariable("df_m_8", &MLP_m_vars[7]);
  reader_taup_PZ->AddVariable("df_m_9", &MLP_m_vars[8]);
  reader_taup_PZ->AddVariable("df_m_10", &MLP_m_vars[9]);
  reader_taup_PZ->AddVariable("df_m_11", &MLP_m_vars[10]);
  reader_taup_PZ->AddVariable("df_m_12", &MLP_m_vars[11]);
  reader_taup_PZ->AddVariable("df_m_13", &MLP_m_vars[12]);
  reader_taup_PZ->AddVariable("df_m_14", &MLP_m_vars[13]);
  reader_taup_PZ->AddVariable("df_m_15", &MLP_m_vars[14]);
  reader_taup_PZ->AddVariable("df_m_16", &MLP_m_vars[15]);
  reader_taup_PZ->AddVariable("df_m_17", &MLP_m_vars[16]);
  reader_taup_PZ->AddVariable("df_m_18", &MLP_m_vars[17]);
  reader_taup_PZ->AddVariable("df_m_19", &MLP_m_vars[18]);
  reader_taup_PZ->AddVariable("df_m_20", &MLP_m_vars[19]);
  reader_taup_PZ->AddVariable("df_m_21", &MLP_m_vars[20]);
  reader_taup_PZ->AddVariable("df_m_22", &MLP_m_vars[21]);
  reader_taup_PZ->AddVariable("Kp_RP_Z", &MLP_RPz);
  reader_taup_PZ->AddSpectator("eventNumber", &MLP_eventNumber);

  reader_taum_PX->AddVariable("df_m_1", &MLP_m_vars[0]);
  reader_taum_PX->AddVariable("df_m_2", &MLP_m_vars[1]);
  reader_taum_PX->AddVariable("df_m_3", &MLP_m_vars[2]);
  reader_taum_PX->AddVariable("df_m_4", &MLP_m_vars[3]);
  reader_taum_PX->AddVariable("df_m_5", &MLP_m_vars[4]);
  reader_taum_PX->AddVariable("df_m_6", &MLP_m_vars[5]);
  reader_taum_PX->AddVariable("df_m_7", &MLP_m_vars[6]);
  reader_taum_PX->AddVariable("df_m_8", &MLP_m_vars[7]);
  reader_taum_PX->AddVariable("df_m_9", &MLP_m_vars[8]);
  reader_taum_PX->AddVariable("df_m_10", &MLP_m_vars[9]);
  reader_taum_PX->AddVariable("df_m_11", &MLP_m_vars[10]);
  reader_taum_PX->AddVariable("df_m_12", &MLP_m_vars[11]);
  reader_taum_PX->AddVariable("df_m_13", &MLP_m_vars[12]);
  reader_taum_PX->AddVariable("df_m_14", &MLP_m_vars[13]);
  reader_taum_PX->AddVariable("df_m_15", &MLP_m_vars[14]);
  reader_taum_PX->AddVariable("df_m_16", &MLP_m_vars[15]);
  reader_taum_PX->AddVariable("df_m_17", &MLP_m_vars[16]);
  reader_taum_PX->AddVariable("df_m_18", &MLP_m_vars[17]);
  reader_taum_PX->AddVariable("df_m_19", &MLP_m_vars[18]);
  reader_taum_PX->AddVariable("df_m_20", &MLP_m_vars[19]);
  reader_taum_PX->AddVariable("df_m_21", &MLP_m_vars[20]);
  reader_taum_PX->AddVariable("df_m_22", &MLP_m_vars[21]);
  reader_taum_PX->AddVariable("Kp_RP_Z", &MLP_RPz);
  reader_taum_PX->AddSpectator("eventNumber", &MLP_eventNumber);

  reader_taum_PY->AddVariable("df_m_1", &MLP_m_vars[0]);
  reader_taum_PY->AddVariable("df_m_2", &MLP_m_vars[1]);
  reader_taum_PY->AddVariable("df_m_3", &MLP_m_vars[2]);
  reader_taum_PY->AddVariable("df_m_4", &MLP_m_vars[3]);
  reader_taum_PY->AddVariable("df_m_5", &MLP_m_vars[4]);
  reader_taum_PY->AddVariable("df_m_6", &MLP_m_vars[5]);
  reader_taum_PY->AddVariable("df_m_7", &MLP_m_vars[6]);
  reader_taum_PY->AddVariable("df_m_8", &MLP_m_vars[7]);
  reader_taum_PY->AddVariable("df_m_9", &MLP_m_vars[8]);
  reader_taum_PY->AddVariable("df_m_10", &MLP_m_vars[9]);
  reader_taum_PY->AddVariable("df_m_11", &MLP_m_vars[10]);
  reader_taum_PY->AddVariable("df_m_12", &MLP_m_vars[11]);
  reader_taum_PY->AddVariable("df_m_13", &MLP_m_vars[12]);
  reader_taum_PY->AddVariable("df_m_14", &MLP_m_vars[13]);
  reader_taum_PY->AddVariable("df_m_15", &MLP_m_vars[14]);
  reader_taum_PY->AddVariable("df_m_16", &MLP_m_vars[15]);
  reader_taum_PY->AddVariable("df_m_17", &MLP_m_vars[16]);
  reader_taum_PY->AddVariable("df_m_18", &MLP_m_vars[17]);
  reader_taum_PY->AddVariable("df_m_19", &MLP_m_vars[18]);
  reader_taum_PY->AddVariable("df_m_20", &MLP_m_vars[19]);
  reader_taum_PY->AddVariable("df_m_21", &MLP_m_vars[20]);
  reader_taum_PY->AddVariable("df_m_22", &MLP_m_vars[21]);
  reader_taum_PY->AddVariable("Kp_RP_Z", &MLP_RPz);
  reader_taum_PY->AddSpectator("eventNumber", &MLP_eventNumber);

  reader_taum_PZ->AddVariable("df_m_1", &MLP_m_vars[0]);
  reader_taum_PZ->AddVariable("df_m_2", &MLP_m_vars[1]);
  reader_taum_PZ->AddVariable("df_m_3", &MLP_m_vars[2]);
  reader_taum_PZ->AddVariable("df_m_4", &MLP_m_vars[3]);
  reader_taum_PZ->AddVariable("df_m_5", &MLP_m_vars[4]);
  reader_taum_PZ->AddVariable("df_m_6", &MLP_m_vars[5]);
  reader_taum_PZ->AddVariable("df_m_7", &MLP_m_vars[6]);
  reader_taum_PZ->AddVariable("df_m_8", &MLP_m_vars[7]);
  reader_taum_PZ->AddVariable("df_m_9", &MLP_m_vars[8]);
  reader_taum_PZ->AddVariable("df_m_10", &MLP_m_vars[9]);
  reader_taum_PZ->AddVariable("df_m_11", &MLP_m_vars[10]);
  reader_taum_PZ->AddVariable("df_m_12", &MLP_m_vars[11]);
  reader_taum_PZ->AddVariable("df_m_13", &MLP_m_vars[12]);
  reader_taum_PZ->AddVariable("df_m_14", &MLP_m_vars[13]);
  reader_taum_PZ->AddVariable("df_m_15", &MLP_m_vars[14]);
  reader_taum_PZ->AddVariable("df_m_16", &MLP_m_vars[15]);
  reader_taum_PZ->AddVariable("df_m_17", &MLP_m_vars[16]);
  reader_taum_PZ->AddVariable("df_m_18", &MLP_m_vars[17]);
  reader_taum_PZ->AddVariable("df_m_19", &MLP_m_vars[18]);
  reader_taum_PZ->AddVariable("df_m_20", &MLP_m_vars[19]);
  reader_taum_PZ->AddVariable("df_m_21", &MLP_m_vars[20]);
  reader_taum_PZ->AddVariable("df_m_22", &MLP_m_vars[21]);
  reader_taum_PZ->AddVariable("Kp_RP_Z", &MLP_RPz);
  reader_taum_PZ->AddSpectator("eventNumber", &MLP_eventNumber);

  reader_taup_PX->BookMVA("MLP", weightfile_taup_PX);
  reader_taup_PY->BookMVA("MLP", weightfile_taup_PY);
  reader_taup_PZ->BookMVA("MLP", weightfile_taup_PZ);
  reader_taum_PX->BookMVA("MLP", weightfile_taum_PX);
  reader_taum_PY->BookMVA("MLP", weightfile_taum_PY);
  reader_taum_PZ->BookMVA("MLP", weightfile_taum_PZ);

  MLP_m_vars[0] = m(1);
  MLP_m_vars[1] = m(2);
  MLP_m_vars[2] = m(3);
  MLP_m_vars[3] = m(4);
  MLP_m_vars[4] = m(5);
  MLP_m_vars[5] = m(6);
  MLP_m_vars[6] = m(7);
  MLP_m_vars[7] = m(8);
  MLP_m_vars[8] = m(9);
  MLP_m_vars[9] = m(10);
  MLP_m_vars[10] = m(11);
  MLP_m_vars[11] = m(12);
  MLP_m_vars[12] = m(13);
  MLP_m_vars[13] = m(14);
  MLP_m_vars[14] = m(15);
  MLP_m_vars[15] = m(16);
  MLP_m_vars[16] = m(17);
  MLP_m_vars[17] = m(18);
  MLP_m_vars[18] = m(19);
  MLP_m_vars[19] = m(20);
  MLP_m_vars[20] = m(21);
  MLP_m_vars[21] = m(22);
  MLP_RPz = refPoint_Kp.z();
  MLP_eventNumber = eventNumber;

  Float_t MLP_taup_PX = reader_taup_PX->EvaluateRegression("MLP")[0];
  Float_t MLP_taup_PY = reader_taup_PY->EvaluateRegression("MLP")[0];
  Float_t MLP_taup_PZ = reader_taup_PZ->EvaluateRegression("MLP")[0];
  Float_t MLP_taum_PX = reader_taum_PX->EvaluateRegression("MLP")[0];
  Float_t MLP_taum_PY = reader_taum_PY->EvaluateRegression("MLP")[0];
  Float_t MLP_taum_PZ = reader_taum_PZ->EvaluateRegression("MLP")[0];

  ROOT::Math::XYZVector ptau1( MLP_taup_PX, MLP_taup_PY, MLP_taup_PZ );
  ROOT::Math::XYZVector ptau2( MLP_taum_PX, MLP_taum_PY, MLP_taum_PZ );
  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
  Double_t Etau1 = sqrt( pow(m_tauMass,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(m_tauMass,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
  Double_t EB = Etau1 + Etau2 + EK;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

  // if(MB_squared > 0)
  // {
  //   std::cout << "MB initial = " << sqrt(MB_squared) << std::endl;
  // }
  // else
  // {
  //   std::cout << "MB initial = " << -sqrt(TMath::Abs(MB_squared)) << std::endl;
  // }

  Gaudi::LorentzVector Pnu1( pnu1.x(), pnu1.y(), pnu1.z(), Enu1 );
  Gaudi::LorentzVector Pnu2( pnu2.x(), pnu2.y(), pnu2.z(), Enu2 );

  Gaudi::LorentzVector Ptau1( ptau1.x(), ptau1.y(), ptau1.z(), Etau1 );
  Gaudi::LorentzVector Ptau2( ptau2.x(), ptau2.y(), ptau2.z(), Etau2 );

  Gaudi::LorentzVector PK( pK.x(), pK.y(), pK.z(), EK );

  Gaudi::LorentzVector PB = Ptau1 + Ptau2 + PK;

  LHCb::Particle * tau1 = const_cast<LHCb::Particle*>(tauList[0]);
  LHCb::Particle * tau2 = const_cast<LHCb::Particle*>(tauList[1]);

  LHCb::Particle::Vector nuList;
  nuList.push_back( new LHCb::Particle() );
  if(mother->charge() > 0)
  {
    nuList.back()->setParticleID( LHCb::ParticleID(-16) ); // antineutrino (partner of a tau+)
  }
  else
  {
    nuList.back()->setParticleID( LHCb::ParticleID(16) ); // neutrino (partner of a tau-)
  }
  nuList.back()->setMomentum( Pnu1 );

  tau1->addToDaughters( nuList.back() );
  tau1->setMomentum( Ptau1 );

  nuList.push_back( new LHCb::Particle() );
  if(mother->charge() > 0)
  {
    nuList.back()->setParticleID( LHCb::ParticleID(16) ); // neutrino (partner of a tau-)
  }
  else
  {
    nuList.back()->setParticleID( LHCb::ParticleID(-16) ); // antineutrino (partner of a tau+)
  }
  nuList.back()->setMomentum( Pnu2 );
  
  tau2->addToDaughters( nuList.back() );
  tau2->setMomentum( Ptau2 );

  treeHead->setMomentum( PB );

  TupleMap tMap ; // the temporary data map
  DecayTreeFitter::Fitter fitter(*treeHead, *originVtx[0], stateprovider ) ;
  if (!fit(fitter, treeHead, originVtx[0], prefix, tMap, tuple, true)) return StatusCode::FAILURE ;

  return fillTuple(tMap,tuple,prefix); // the actual filling
  //return StatusCode(true);
}
//=============================================================================
// do filling for a given vertex
//=============================================================================
StatusCode TupleToolKTauTauDTF_MLP::fit(DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
                                              const LHCb::VertexBase* pv, const std::string& prefix,
                                              TupleMap& tMap, Tuples::Tuple& tuple,
                                              bool fillIfFailed) const
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

  if(fitter.status() == 0 || fillIfFailed)
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
    fillEndVtx(prefix+"_taup", tMap, DV1[fitter.nIter()]);
    fillEndVtx(prefix+"_taum", tMap, DV2[fitter.nIter()]);

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
  

  return StatusCode(test);
}
//=============================================================================
// Fill endvertex information from DTF for intermediates
//=============================================================================
StatusCode TupleToolKTauTauDTF_MLP::fillEndVtx(const std::string& prefix, TupleMap& tMap, Gaudi::XYZPoint vtx) const
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
StatusCode TupleToolKTauTauDTF_MLP::fillPV(const LHCb::VertexBase* pv, const std::string& prefix,
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
StatusCode TupleToolKTauTauDTF_MLP::fillChi2Iter(DecayTreeFitter::Fitter& fitter, std::vector<float> chisq_iters,
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
StatusCode TupleToolKTauTauDTF_MLP::fillDecay(const DecayTreeFitter::Fitter& fitter, const std::string& prefix,
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
StatusCode TupleToolKTauTauDTF_MLP::fillVtxIter(const DecayTreeFitter::Fitter& fitter, const std::string& prefix,
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
StatusCode TupleToolKTauTauDTF_MLP::fillMomentum(const DecayTreeFitter::Fitter& fitter, const Particle* P,
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
StatusCode TupleToolKTauTauDTF_MLP::fillLT(const DecayTreeFitter::Fitter& fitter, const Particle* P,
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
StatusCode TupleToolKTauTauDTF_MLP::fillDaughters(const DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
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
TupleToolKTauTauDTF_MLP::fillStableDaughters(const DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
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
TupleToolKTauTauDTF_MLP::fillTracksMomentum(const DecayTreeFitter::Fitter& fitter, const Particle* P,
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
StatusCode TupleToolKTauTauDTF_MLP::insert(const std::string& leaf, const double val,
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
StatusCode TupleToolKTauTauDTF_MLP::fillTuple(TupleMap& tMap, Tuples::Tuple& tuple,
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
TupleToolKTauTauDTF_MLP::sortedTracks(const LHCb::VertexBase* vb) const
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
bool TupleToolKTauTauDTF_MLP::samePV( const LHCb::VertexBase* vb1,
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
TupleToolKTauTauDTF_MLP::originVertex( const Particle* mother, const Particle* P ) const
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
std::string TupleToolKTauTauDTF_MLP::getName(const int id) const
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
StatusCode TupleToolKTauTauDTF_MLP::substitute(LHCb::DecayTree& tree)
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
StatusCode TupleToolKTauTauDTF_MLP::checkMassConstraints(const LHCb::DecayTree& tree)
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

Gaudi::XYZVector TupleToolKTauTauDTF_MLP::makeTransformation_vec(Gaudi::XYZVector p_K, Gaudi::XYZPoint refPoint, 
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

Gaudi::XYZPoint TupleToolKTauTauDTF_MLP::makeTransformation_point(Gaudi::XYZVector p_K, Gaudi::XYZPoint refPoint, 
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
DECLARE_COMPONENT( TupleToolKTauTauDTF_MLP )