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
#include "B2KtautauDTF.h"
#include "TMath.h"
#include <algorithm>

using namespace LHCb;

//-----------------------------------------------------------------------------
// Implementation file for class : B2KtautauDTF
// Yasmine Amhis, Matt Needham, Patrick Koppenburg
// 30-10-2010, 01-04-2011
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
B2KtautauDTF::B2KtautauDTF( const std::string& type,
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
  declareProperty( "maxNumberOfIterations", m_maxNiter = 1000, "Maximum number of iterations during fitting." );
  declareProperty( "maxndiverging", m_maxndiverging = 15, "Maximum number of diverging iterations during fitting." );
  declareProperty( "dChisqQuit", m_dChisqQuit = 2e5, "Maximum number of diverging iterations during fitting." ); //currently not used by code

  declareInterface<IParticleTupleTool>(this);
}

//=============================================================================
StatusCode B2KtautauDTF::initialize()
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
    if ( "B2KtautauDTF" == m_extraName )  m_extraName = ""; // user has not changed instance name
    info() << "All fields will be prepended with ``" << m_extraName << "''" <<endmsg;
  }

  if ( m_extraName.empty() )
  {
    return Error( "Extraname is empty. Always give an instance name "
                  "to B2KtautauDTF! See doxygen." );
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

StatusCode B2KtautauDTF::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  if ( !m_stateprovider.empty() ) { sc = m_stateprovider.release(); }
  return StatusCode{ TupleToolBase::finalize() && sc };
}

//=============================================================================
//  The fill method implementation. This is where our modifications are happening
//=============================================================================
StatusCode B2KtautauDTF::fill(const LHCb::Particle* mother, const LHCb::Particle* P,
                              const std::string& head, Tuples::Tuple& tuple )
{

//  bool runningOnMC = true;

  if( !P ) return StatusCode::FAILURE;
  if ( P->isBasicParticle() )
  {
    return Error("Do not call B2KtautauDTF for basic particles. Use Branches. See doxygen.");
  }
  const std::string prefix = fullName(head);
  if (msgLevel(MSG::DEBUG)) debug() << "head ''" << head << "'' prefix ''" << prefix
                                    << "'' extraname ''" << m_extraName << "''" <<endmsg;

  const auto stateprovider = ( m_stateprovider.empty() ? nullptr : &(*m_stateprovider) );

  LHCb::DecayTree tree ( *P ) ;
  LHCb::Particle * treeHead = tree.head();

  //#####
  // m_assoc = tool<IParticle2MCAssociator >("MCMatchObjP2MCRelator");
  // const LHCb::MCParticle* treeHead_true = m_assoc->relatedMCP(treeHead, LHCb::MCParticleLocation::Default);
  //#####

  //Arvind's Notes
  //1. Get PV
  //2. Get reference point for K+ track (using LHCb::Particle::referencePoint())
  //3. Get K+ 3-momentum
  //4. With 2 and 3, the coordinate transformation can be defined
  //5. Get 3pi1 and 3pi2 vertices
  //6. With these ingredients, the transformation can be defined. 
  //   Apply the transformation to the necessary points and vectors.
  //7. Do the initialization calculation to get recalculated tau momenta
  //8. Inverse transform back to the LHCb reference frame
  //9. Add a neutrino daughter to each tau, 
  //   and set its momentum to the difference between the 
  //   recalculated tau momentum and original 3pi momentum
  //10. Set tau momenta to recalculated values
  //11. Set B momentum to K momenta + sum of recalculated tau momenta.
  //12. Now run the usual DTF fit.
  //
  //End Arvind Notes
  
  //1. Get PV
  std::vector<const VertexBase*> originVtx;
  originVtx = originVertex( mother, P );
  if( originVtx.empty() ){return Error("Can't get an origin vertex");}
  Gaudi::XYZPoint PV( originVtx[0]->position().X(), originVtx[0]->position().Y(), originVtx[0]->position().Z() ); //Notice origVtx[0] is being used i.e. only "best PV"

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
  
  //2. Get K+ track reference point
  Gaudi::XYZPoint refPoint_Kplus = Kplus->referencePoint(); //ostensibly a point on the K+ trajectory

  //Other method of getting reference point
  const LHCb::ProtoParticle *proto_Kplus = Kplus->proto();
  const LHCb::Track *track_Kplus         = proto_Kplus->track();
  LHCb::State firstState_Kplus     = track_Kplus->firstState();
  Gaudi::XYZPoint refPoint_Kplus_alt = firstState_Kplus.position();

  //Get Kplus 3-momentum + ennergy
  Gaudi::XYZVector Pk = Kplus->momentum().Vect(); //Kplus->momentum() gives 4-momentum
  double Ek = Kplus->momentum().E();

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
  double E3pi1 = tauList[0]->momentum().E();
  double E3pi2 = tauList[1]->momentum().E();
  Gaudi::XYZVector p_3pi1 = tauList[0]->momentum().Vect();
  Gaudi::XYZVector p_3pi2 = tauList[1]->momentum().Vect();

  double pmag_3pi1 = p_3pi1.r();
  double pmag_3pi2 = p_3pi2.r();  

  // Get tau1 and tau2 end vertices
  Gaudi::XYZPoint DV1(tauList[0]->endVertex()->position().X(), tauList[0]->endVertex()->position().Y(), tauList[0]->endVertex()->position().Z());
  Gaudi::XYZPoint DV2(tauList[1]->endVertex()->position().X(), tauList[1]->endVertex()->position().Y(), tauList[1]->endVertex()->position().Z());

  //apply the transformations
  Gaudi::XYZVector Pk_t     = makeTransformation_vec(Pk, refPoint_Kplus, Pk, false);
  Gaudi::XYZVector P3pi1_t  = makeTransformation_vec(Pk, refPoint_Kplus, p_3pi1, false);
  Gaudi::XYZVector P3pi2_t  = makeTransformation_vec(Pk, refPoint_Kplus, p_3pi2, false);
  Gaudi::XYZPoint DV1_t      = makeTransformation_point(Pk, refPoint_Kplus, DV1, false);
  Gaudi::XYZPoint DV2_t      = makeTransformation_point(Pk, refPoint_Kplus, DV2, false);
  Gaudi::XYZPoint refPoint_t = makeTransformation_point(Pk, refPoint_Kplus, refPoint_Kplus, false);
  Gaudi::XYZPoint PV_t       = makeTransformation_point(Pk, refPoint_Kplus, PV, false);

  // tau and 3pi masses
  Double_t mtau = 1776.86;
  Double_t m3pi1 = sqrt( pow(E3pi1,2) -  P3pi1_t.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) -  P3pi2_t.Mag2() );

  //Define the internal variables  
  Double_t a1 = (DV1_t.y())/(DV1_t.x());
  Double_t a2 = (DV2_t.y())/(DV2_t.x());
  Double_t b = (PV_t.y() -a1*PV_t.x())/(a2*PV_t.x() - PV_t.y());
  Double_t c = b*(DV1_t.x())/(DV2_t.x());
  Double_t d = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()));
  Double_t e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
  Double_t f = ( PV_t.x()*sqrt(Pk_t.Mag2()) )/( (1+b)*DV1_t.x() -(1+c)*PV_t.x() );
  Double_t g = c*e + d;
  Double_t h = f*c;
  Double_t i = DV1_t.z() - e*DV1_t.x();
  Double_t j = f*DV1_t.x();

  Double_t x1 = P3pi1_t.x() + a1*P3pi1_t.y() + e*P3pi1_t.z();
  Double_t x2 = b*P3pi2_t.x() + a2*b*P3pi2_t.y() + g*P3pi2_t.z();

  Double_t p1 = 1 + pow(a1,2) + pow(e,2) - pow(x1/E3pi1,2);
  Double_t p2 = 2*e*f - ( pow(mtau,2) + pow(m3pi1,2) + 2*f*P3pi1_t.z() )*(x1/pow(E3pi1,2) );
  Double_t p3 = pow(mtau,2) + pow(f,2) - pow( ( pow(mtau,2) + pow(m3pi1,2) + 2*f*P3pi1_t.z() )/(2*E3pi1), 2);
  Double_t q1 = pow(b,2) + pow(a2*b,2) + pow(g,2) - pow(x2/E3pi2,2);
  Double_t q2 = 2*g*h - ( pow(mtau,2) + pow(m3pi2,2) + 2*h*P3pi2_t.z() )*(x2/pow(E3pi2,2) );
  Double_t q3 = pow(mtau,2) + pow(h,2) - pow( ( pow(mtau,2) + pow(m3pi2,2) + 2*h*P3pi2_t.z() )/(2*E3pi2),2 );

  Double_t Ptau1x_t = (p1*q3 - p3*q1)/(p2*q1 - p1*q2);
  Double_t Ptau1y_t = a1*Ptau1x_t;
  Double_t Ptau1z_t = e*Ptau1x_t + f;
  Double_t Ptau2x_t = b*Ptau1x_t;
  Double_t Ptau2y_t = a2*b*Ptau1x_t;
  Double_t Ptau2z_t = g*Ptau1x_t + h;
  Double_t BVz_t = i - j*(1/Ptau1x_t);

  ROOT::Math::XYZVector Ptau1_t(Ptau1x_t, Ptau1y_t, Ptau1z_t);
  ROOT::Math::XYZVector Ptau2_t(Ptau2x_t, Ptau2y_t, Ptau2z_t); 
  ROOT::Math::XYZPoint BV_t(0., 0., BVz_t);
  ROOT::Math::XYZVector Pnu1_t = Ptau1_t - P3pi1_t;
  ROOT::Math::XYZVector Pnu2_t = Ptau2_t - P3pi2_t;
  ROOT::Math::XYZVector Pb_t = Ptau1_t + Ptau2_t + Pk_t;

  // scale factors
  Double_t L1x, L1y, L1z;
  Double_t L2x, L2y, L2z;
  Double_t Lx, Ly, Lz;

  L1x = Ptau1_t.x()/DV1_t.x();
  L1y = Ptau1_t.y()/DV1_t.y();
  L1z = Ptau1_t.z()/(DV1_t.z() - BV_t.z());
  L2x = Ptau2_t.x()/DV2_t.x();
  L2y = Ptau2_t.y()/DV2_t.y();
  L2z = Ptau2_t.z()/(DV2_t.z() - BV_t.z());
  Lx = -Pb_t.x()/PV_t.x();
  Ly = -Pb_t.y()/PV_t.y();
  Lz = Pb_t.z()/(BV_t.z() - PV_t.z());

  // make transformation back to LHCb frame
  ROOT::Math::XYZVector Ptau1 = makeTransformation_vec(Pk, refPoint_Kplus, Ptau1_t, true);
  ROOT::Math::XYZVector Ptau2 = makeTransformation_vec(Pk, refPoint_Kplus, Ptau2_t, true);
  ROOT::Math::XYZPoint BV = makeTransformation_point(Pk, refPoint_Kplus, BV_t, true);
  ROOT::Math::XYZVector Pnu1 = makeTransformation_vec(Pk, refPoint_Kplus, Pnu1_t, true);
  ROOT::Math::XYZVector Pnu2 = makeTransformation_vec(Pk, refPoint_Kplus, Pnu2_t, true);
  ROOT::Math::XYZVector Pb = makeTransformation_vec(Pk, refPoint_Kplus, Pb_t, true);

  Double_t Etau1 = sqrt( pow(mtau,2) + Ptau1_t.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + Ptau2_t.Mag2() );
  Double_t Enu1 = sqrt(Pnu1_t.Mag2());
  Double_t Enu2 = sqrt(Pnu2_t.Mag2());
  Double_t Eb = Etau1 + Etau2 + Ek;
  Double_t Mb = sqrt( pow(Eb,2) - Pb.Mag2() );
  Double_t tau1_M = sqrt( pow(Etau1,2) - Ptau1_t.Mag2() );
  Double_t tau2_M = sqrt( pow(Etau2,2) - Ptau2_t.Mag2() );
  Double_t nu1_M = sqrt( pow(Enu1,2) - Pnu1_t.Mag2() );
  Double_t nu2_M = sqrt( pow(Enu2,2) - Pnu2_t.Mag2() );

  Gaudi::LorentzVector P4_tau1(Ptau1.x(), Ptau1.y(), Ptau1.z(), Etau1);
  Gaudi::LorentzVector P4_tau2(Ptau2.x(), Ptau2.y(), Ptau2.z(), Etau2);
  Gaudi::LorentzVector P4_nu1(Pnu1.x(), Pnu1.y(), Pnu1.z(), Enu1);
  Gaudi::LorentzVector P4_nu2(Pnu2.x(), Pnu2.y(), Pnu2.z(), Enu2);
  Gaudi::LorentzVector P4_b(Pb.x(), Pb.y(), Pb.z(), Eb);

  // Antineutrino 4-momentum
  LHCb::Particle::Vector nuList;
  nuList.push_back( new LHCb::Particle() );
  nuList.back()->setParticleID( LHCb::ParticleID(16) ); //Looks like the order of neutrino PID's doesn't matter. Check again later to make sure.
  nuList.back()->setMomentum( P4_nu1 );//NB 4 momenta going in here

  // tau+ 4-momentum
  LHCb::Particle * tau1 = const_cast<LHCb::Particle*>(tauList[0]);
  tau1->addToDaughters( nuList.back() );
  tau1->setMomentum( P4_tau1 );

  // Neutrino 4-momentum
  nuList.push_back( new LHCb::Particle() );
  nuList.back()->setParticleID( LHCb::ParticleID(-16) );
  nuList.back()->setMomentum( P4_nu2 );

  // tau- 4-momentum
  LHCb::Particle * tau2 = const_cast<LHCb::Particle*>(tauList[1]);
  tau2->addToDaughters( nuList.back() );
  tau2->setMomentum( P4_tau2 );

  // B+ 4-momentum
  treeHead->setMomentum( P4_b );
  // B+ decay vertex
  const LHCb::Vertex* BV_vertex = new LHCb::Vertex( BV );
  treeHead->setEndVertex( BV_vertex );

  TupleMap tMap ; // the temporary data map
  DecayTreeFitter::Fitter fitter(*treeHead, *originVtx[0], stateprovider ) ;
  if (!fit(fitter, treeHead, originVtx[0], prefix, tMap, tuple, true)) return StatusCode::FAILURE ;

  tuple->column("refPoint_X", refPoint_Kplus.x());
  tuple->column("refPoint_Y", refPoint_Kplus.y());
  tuple->column("refPoint_Z", refPoint_Kplus.z());

  return fillTuple(tMap,tuple,prefix); // the actual filling
  //return StatusCode(true);
}
//=============================================================================
// do filling for a given vertex
//=============================================================================
StatusCode B2KtautauDTF::fit(DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
                            const LHCb::VertexBase* pv, const std::string& prefix,
                            TupleMap& tMap, Tuples::Tuple& tuple, 
                            bool fillIfFailedn) const
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

  //calling the actual fit function.
  fitter.fit(m_maxNiter, 0.01, m_maxndiverging, m_dChisqQuit, chisq_iters, 
             BV, DV1, DV2, P, B_M,
             p_taup_nu, p_taum_nu,
             p_taup, p_taum, true);

  if (msgLevel(MSG::VERBOSE)) verbose() << "called Fit" << endmsg ;

  //fill chi2 for each fitter iteration
  fillChi2Iter(fitter, chisq_iters, prefix, tuple); //made by AV

  //fill decay vertices after each fitter iteration
  fillVtxIter(fitter, prefix, tuple, BV, 
              DV1, DV2, B_M, p_taup_nu, p_taum_nu,
              p_taup, p_taum); //made by AV

  // fill the fit result
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

  return StatusCode(test);
}

//=============================================================================
// Fill standard stuff
//=============================================================================
StatusCode B2KtautauDTF::fillPV(const LHCb::VertexBase* pv, const std::string& prefix,
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
StatusCode B2KtautauDTF::fillChi2Iter(DecayTreeFitter::Fitter& fitter, std::vector<float> chisq_iters,
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
StatusCode B2KtautauDTF::fillDecay(const DecayTreeFitter::Fitter& fitter, const std::string& prefix,
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
StatusCode B2KtautauDTF::fillVtxIter(const DecayTreeFitter::Fitter& fitter, const std::string& prefix,
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
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Bp_ENDVERTEX_Y", BV_Y.begin(), BV_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Bp_ENDVERTEX_Z", BV_Z.begin(), BV_Z.end(),
                        prefix+"_nIters1", m_maxNiter);

  test &= tuple->farray(prefix+"_Taup_ENDVERTEX_X", DV1_X.begin(), DV1_X.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_ENDVERTEX_Y", DV1_Y.begin(), DV1_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_ENDVERTEX_Z", DV1_Z.begin(), DV1_Z.end(),
                        prefix+"_nIters1", m_maxNiter);

  test &= tuple->farray(prefix+"_Taum_ENDVERTEX_X", DV2_X.begin(), DV2_X.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_ENDVERTEX_Y", DV2_Y.begin(), DV2_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_ENDVERTEX_Z", DV2_Z.begin(), DV2_Z.end(),
                        prefix+"_nIters1", m_maxNiter);

  test &= tuple->farray(prefix+"_Taup_nu_PX", p_taup_nu_X.begin(), p_taup_nu_X.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_nu_PY", p_taup_nu_Y.begin(), p_taup_nu_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_nu_PZ", p_taup_nu_Z.begin(), p_taup_nu_Z.end(),
                        prefix+"_nIters1", m_maxNiter);
  
  test &= tuple->farray(prefix+"_Taum_nu_PX", p_taum_nu_X.begin(), p_taum_nu_X.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_nu_PY", p_taum_nu_Y.begin(), p_taum_nu_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_nu_PZ", p_taum_nu_Z.begin(), p_taum_nu_Z.end(),
                        prefix+"_nIters1", m_maxNiter);

  test &= tuple->farray(prefix+"_Taup_PX", p_taup_X.begin(), p_taup_X.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_PY", p_taup_Y.begin(), p_taup_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taup_PZ", p_taup_Z.begin(), p_taup_Z.end(),
                        prefix+"_nIters1", m_maxNiter);
  
  test &= tuple->farray(prefix+"_Taum_PX", p_taum_X.begin(), p_taum_X.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_PY", p_taum_Y.begin(), p_taum_Y.end(),
                        prefix+"_nIters1", m_maxNiter);
  test &= tuple->farray(prefix+"_Taum_PZ", p_taum_Z.begin(), p_taum_Z.end(),
                        prefix+"_nIters1", m_maxNiter);
  
  // test &= tuple->farray(prefix+"_M_iters", B_M.begin(), B_M.end(),
  //                       prefix+"_nIters1", m_maxNiter);
  return StatusCode(test);
}

//=============================================================================
// Fill momentum and mass information
//=============================================================================
StatusCode B2KtautauDTF::fillMomentum(const DecayTreeFitter::Fitter& fitter, const Particle* P,
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
StatusCode B2KtautauDTF::fillLT(const DecayTreeFitter::Fitter& fitter, const Particle* P,
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
StatusCode B2KtautauDTF::fillDaughters(const DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
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
B2KtautauDTF::fillStableDaughters(const DecayTreeFitter::Fitter& fitter, const LHCb::Particle* P,
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
B2KtautauDTF::fillTracksMomentum(const DecayTreeFitter::Fitter& fitter, const Particle* P,
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
StatusCode B2KtautauDTF::insert(const std::string& leaf, const double val,
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
StatusCode B2KtautauDTF::fillTuple(TupleMap& tMap, Tuples::Tuple& tuple,
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
B2KtautauDTF::sortedTracks(const LHCb::VertexBase* vb) const
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
bool B2KtautauDTF::samePV( const LHCb::VertexBase* vb1,
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
B2KtautauDTF::originVertex( const Particle* mother, const Particle* P ) const
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
std::string B2KtautauDTF::getName(const int id) const
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
StatusCode B2KtautauDTF::substitute(LHCb::DecayTree& tree)
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
StatusCode B2KtautauDTF::checkMassConstraints(const LHCb::DecayTree& tree)
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

Gaudi::XYZVector B2KtautauDTF::makeTransformation_vec(Gaudi::XYZVector Pk, Gaudi::XYZPoint refPoint, 
                                                      Gaudi::XYZVector theVector, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  double theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

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

Gaudi::XYZPoint B2KtautauDTF::makeTransformation_point(Gaudi::XYZVector Pk, Gaudi::XYZPoint refPoint, 
                                                              Gaudi::XYZPoint thePoint, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  double theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

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
DECLARE_COMPONENT( B2KtautauDTF )
