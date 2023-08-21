/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

// local
#include "TupleToolKtautauDTF.h"

using namespace LHCb;

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolKtautauDTF
// Yasmine Amhis, Matt Needham, Patrick Koppenburg
// 30-10-2010, 01-04-2011
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolKtautauDTF::TupleToolKtautauDTF( const std::string& type, const std::string& name,
                                                    const IInterface* parent )
    : TupleToolBase( type, name, parent ) {
  declareProperty( "daughtersToConstrain", m_massConstraints, "List of particles to contrain to mass" );
  declareProperty( "constrainToOriginVertex", m_constrainToOriginVertex = false,
                   "Do a refit constraining to Origin Vertex (could be PV)" );
  declareProperty( "Substitutions", m_map, "PID-substitutions :  { ' decay-component' : 'new-pid' }" );
  declareProperty( "StoreRefittedPVsTwice", m_storeAnyway = false,
                   "Store PV even if a refitted version is already the best PV (i.e store twice)." );
  declareProperty( "UpdateDaughters", m_updateDaughters = false, "Store updated momenta of tracks in the decay tree." );
  declareProperty( "StateProvider", m_stateprovider );
  declareProperty( "VetoPIDs", m_vetoPIDs, "An optional list of PDG particle codes to skip when filling the tuple." );
  declareProperty( "UseFullTreeInName", m_useFullTreeInName = false,
                   "Use an improved branch naming scheme that includes a full history of the "
                   "parents and grand-parents for each particle. Makes it easier to identify "
                   "cases where the same particle type appears at different levels in the decay tree." );
  declareInterface<IParticleTupleTool>( this );
}

//=============================================================================
StatusCode TupleToolKtautauDTF::initialize() {
  StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  // convert the list of names to a list of pids
  m_ppSvc = svc<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );
  for ( const auto& S : m_massConstraints ) {
    const auto prop = m_ppSvc->find( S );
    if ( !prop ) Exception( "Unknown PID" );
    m_massConstraintsPids.push_back( prop->pdgID() );
  }

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc(), this );
  if ( !m_dva ) return Error( "Couldn't get parent IDVAlgorithm", StatusCode::FAILURE );

  m_particleDescendants = tool<IParticleDescendants>( "ParticleDescendants" );

  if ( !m_stateprovider.empty() ) {
    sc = m_stateprovider.retrieve();
    if ( sc.isFailure() ) return sc;
  }

  if ( m_extraName.empty() ) {
    const auto en = name(); // use tool name as prepended name
    const auto d  = en.find_last_of( "." );
    m_extraName   = en.substr( d + 1, en.size() - 1 );                 // from d to end
    if ( "TupleToolKtautauDTF" == m_extraName ) m_extraName = ""; // user has not chanegd instance name
    info() << "All fields will be prepended with ``" << m_extraName << "''" << endmsg;
  }

  if ( m_extraName.empty() ) {
    return Error( "Extraname is empty. Always give an instance name "
                  "to TupleToolKtautauDTF! See doxygen." );
  }

  if ( !m_map.empty() ) {
    m_substitute = tool<ISubstitutePID>( "SubstitutePIDTool", this );
    sc           = m_substitute->decodeCode( m_map );
  }

  if ( !m_vetoPIDs.empty() ) { info() << "Will veto PIDs " << m_vetoPIDs << " from filling" << endmsg; }

  if ( m_useFullTreeInName ) { info() << "Will use the full decay tree as part of branch names" << endmsg; }

  return sc;
}

StatusCode TupleToolKtautauDTF::finalize() {
  StatusCode sc = StatusCode::SUCCESS;
  if ( !m_stateprovider.empty() ) { sc = m_stateprovider.release(); }
  return StatusCode{TupleToolBase::finalize() && sc};
}

//=============================================================================
//  The fill method implementation
//=============================================================================
StatusCode TupleToolKtautauDTF::fill( const LHCb::Particle* mother, const LHCb::Particle* P,
                                           const std::string& head, Tuples::Tuple& tuple ) {
  if ( !P ) return StatusCode::FAILURE;
  if ( P->isBasicParticle() ) {
    return Error( "Do not call TupleToolKtautauDTF for basic particles. Use Branches. See doxygen." );
  }
  const std::string prefix = fullName( head );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "head ''" << head << "'' prefix ''" << prefix << "'' extraname ''" << m_extraName << "''" << endmsg;

  const auto stateprovider = ( m_stateprovider.empty() ? nullptr : &( *m_stateprovider ) );

  LHCb::DecayTree tree( *P );
  // substitute
  // if ( m_substitute && !substitute( tree ) ) return StatusCode::FAILURE;
  // if ( !m_massConstraintsPids.empty() && !checkMassConstraints( tree ) ) return StatusCode::FAILURE;

  LHCb::Particle * treeHead = tree.head();

  const LHCb::Particle* Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(321) );
  if( Kp == 0 ){ Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-321)); }

  LHCb::Particle::ConstVector tauList;
  LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList);
  LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList);

  LHCb::Particle* taum = const_cast<LHCb::Particle*>(tauList[0]);
  LHCb::Particle* taup = const_cast<LHCb::Particle*>(tauList[1]);

  LHCb::Particle* antineutrino = new LHCb::Particle();
  antineutrino->setParticleID( LHCb::ParticleID(-16) );
  LHCb::Particle* neutrino = new LHCb::Particle();
  neutrino->setParticleID( LHCb::ParticleID(16) );

  auto P_taup = taup->momentum();
  auto P_taum = taum->momentum();
  auto P_KpTaupTaum = P_taup + P_taum + Kp->momentum();

  auto dir_taup = P_taup.Vect().Unit() ;
  auto dir_taum = P_taum.Vect().Unit() ;

  const double pdgBMass = 5279.34;
  double dM2 = pdgBMass*pdgBMass - P_KpTaupTaum.M2();

  double E_antineutrino = 0.5 * dM2 / (P_KpTaupTaum.E() - P_KpTaupTaum.Vect().Dot(dir_taup));
  double E_neutrino = 0.5 * dM2 / (P_KpTaupTaum.E() - P_KpTaupTaum.Vect().Dot(dir_taum));

  antineutrino->setMomentum( Gaudi::LorentzVector(E_antineutrino*dir_taup.x(), E_antineutrino*dir_taup.y(), E_antineutrino*dir_taup.z(), E_antineutrino) );
  neutrino->setMomentum( Gaudi::LorentzVector(E_neutrino*dir_taum.x(), E_neutrino*dir_taum.y(), E_neutrino*dir_taum.z(), E_neutrino) );

  taup->setMomentum( P_taup + antineutrino->momentum() );
  taup->addToDaughters( antineutrino );

  taum->setMomentum( P_taum + neutrino->momentum() );
  taum->addToDaughters( neutrino );

  // Point of closest approach of kaon, tau+ DV and tau- DV -> BV (B+ endvertex?)

  treeHead->setMomentum( P_taup + P_taum + Kp->momentum() + antineutrino->momentum() + neutrino->momentum() );

  // get origin vertices
  std::vector<const VertexBase*> originVtx;
  TupleMap                       tMap; // the temporary data map

  if ( m_constrainToOriginVertex ) {
    if ( msgLevel( MSG::DEBUG ) ) { debug() << "Constrain the origin vertex" << endmsg; }
    // check for origin vertex
    originVtx = originVertex( mother, P );
    if ( originVtx.empty() ) { return Error( "Can't get an origin vertex" ); }
    if ( msgLevel( MSG::DEBUG ) ) debug() << "PVs: " << originVtx.size() << endmsg;
    for ( const auto& v : originVtx ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Creating DecayTreeFitter on " << treeHead << " " << v << endmsg;

      DecayTreeFitter::Fitter Fitter( *treeHead, *v, stateprovider );
      
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Created DecayTreeFitter" << endmsg;
      if ( !fit( Fitter, treeHead, v, prefix, tMap ) ) return StatusCode::FAILURE;
    }
  } else {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "Do not constrain the origin vertex" << endmsg;
    // Get the Fitter
    DecayTreeFitter::Fitter Fitter( *treeHead, stateprovider );
    if ( !fit( Fitter, treeHead, 0, prefix, tMap ) ) return StatusCode::FAILURE;
  }

  return fillTuple( tMap, tuple, prefix ); // the actual filling
}
//=============================================================================
// do filling for a given vertex
//=============================================================================
StatusCode TupleToolKtautauDTF::fit( DecayTreeFitter::Fitter& Fitter, const LHCb::Particle* P,
                                          const LHCb::VertexBase* pv, const std::string& prefix,
                                          TupleMap& tMap ) const {
  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "fit " << P << " " << pv << " " << prefix << endmsg;
  bool test = true;
  // add mass contraints
  if ( !m_massConstraintsPids.empty() ) {
    for ( const auto& C : m_massConstraintsPids ) { Fitter.setMassConstraint( C ); }
  }
  // fit
  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "calling Fit" << endmsg;
  Fitter.fit();
  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "called Fit" << endmsg;
  // fill the fit result
  fillDecay( Fitter, prefix, tMap ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  fillMomentum( Fitter, P, prefix, tMap ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  if ( m_constrainToOriginVertex ) {
    test &= fillPV( pv, prefix, tMap );
    test &= fillLT( Fitter, P, prefix, tMap );
  }
  if ( isVerbose() ) { test &= fillDaughters( Fitter, P, prefix, tMap ); }
  if ( m_updateDaughters ) { test &= fillStableDaughters( Fitter, P, prefix, tMap ); }

  return StatusCode( test );
}
//=============================================================================
// Fill standard stuff
//=============================================================================
StatusCode TupleToolKtautauDTF::fillPV( const LHCb::VertexBase* pv, const std::string& prefix,
                                             TupleMap& tMap ) const {
  bool test = true;
  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillPV " << prefix << endmsg;
  if ( !pv ) Exception( "Null PVs cannot happen with ConstrainToOriginVertex!" );
  test &= insert( prefix + "_PV_key", pv->key(), tMap );
  if ( isVerbose() ) {
    test &= insert( prefix + "_PV_X", pv->position().X(), tMap );
    test &= insert( prefix + "_PV_Y", pv->position().Y(), tMap );
    test &= insert( prefix + "_PV_Z", pv->position().Z(), tMap );
  }
  return StatusCode( test );
}
//=============================================================================
// Fill standard stuff
//=============================================================================
StatusCode TupleToolKtautauDTF::fillDecay( const DecayTreeFitter::Fitter& Fitter, const std::string& prefix,
                                                TupleMap& tMap ) const {
  bool test = true;
  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillDecay " << prefix << endmsg;

  test &= insert( prefix + "_status", Fitter.status(), tMap );
  test &= insert( prefix + "_nDOF", Fitter.nDof(), tMap );
  test &= insert( prefix + "_chi2", Fitter.chiSquare(), tMap );
  test &= insert( prefix + "_nIter", Fitter.nIter(), tMap );

  return StatusCode( test );
}
//=============================================================================
// Fill momentum and mass information
//=============================================================================
StatusCode TupleToolKtautauDTF::fillMomentum( const DecayTreeFitter::Fitter& Fitter, const Particle* P,
                                                   const std::string& prefix, TupleMap& tMap ) const {
  bool test = true;

  if ( isVetoed( P->particleID().pid() ) ) { return StatusCode( test ); }

  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillMomentum " << prefix << endmsg;
  // Get the fit parameters
  const auto  params   = Fitter.fitParams( P );
  const auto& momentum = params->momentum();

  test &= insert( prefix + "_M", momentum.m().value(), tMap );
  test &= insert( prefix + "_MERR", momentum.m().error(), tMap );
  test &= insert( prefix + "_P", momentum.p().value(), tMap );
  test &= insert( prefix + "_PERR", momentum.p().error(), tMap ); // MeV

  return StatusCode( test );
}
//=============================================================================
// Fill lifetime information
//=============================================================================
StatusCode TupleToolKtautauDTF::fillLT( const DecayTreeFitter::Fitter& Fitter, const Particle* P,
                                             const std::string& prefix, TupleMap& tMap ) const {
  bool test = true;

  if ( isVetoed( P->particleID().pid() ) ) { return StatusCode( test ); }

  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillLT " << prefix << endmsg;
  const auto  tParams     = Fitter.fitParams( P );
  const auto& decayLength = tParams->decayLength();
  const auto& ctau        = tParams->ctau();
  test &= insert( prefix + "_ctau", ctau.value(), tMap );
  test &= insert( prefix + "_ctauErr", ctau.error(), tMap );
  test &= insert( prefix + "_decayLength", decayLength.value(), tMap );
  test &= insert( prefix + "_decayLengthErr", decayLength.error(), tMap );

  return StatusCode( test );
}
//=============================================================================
// Fill lifetime information for non stable daughters
//=============================================================================
StatusCode TupleToolKtautauDTF::fillDaughters( const DecayTreeFitter::Fitter& Fitter, const LHCb::Particle* P,
                                                    const std::string& prefix, TupleMap& tMap ) const {
  bool test = true;

  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillDaughters " << prefix << endmsg;
  const auto& daughters = ( m_useFullTreeInName ? P->daughtersVector() : m_particleDescendants->descendants( P ) );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "for id " << P->particleID().pid() << " daughter size is " << daughters.size() << endmsg;
  if ( daughters.empty() ) return StatusCode( test );
  std::set<std::string> usedNames;
  unsigned int          add = 0;
  for ( const auto& particle : daughters ) {
    if ( particle->isBasicParticle() ) continue;
    const auto pid     = abs( particle->particleID().pid() );
    const auto pidName = getName( pid );
    auto       name    = prefix + "_" + pidName;
    bool       renamed = false;
    while ( usedNames.find( name ) != usedNames.end() ) { // fix to bug 88702
      renamed = true;
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Found already name " << name << " trying next " << endmsg;
      name = prefix + "_" + pidName + "_" + boost::lexical_cast<std::string>( add );
      ++add;
    }
    if ( renamed )
      Info( "Renaming duplicate to " + name, StatusCode::SUCCESS, 1 )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    usedNames.insert( name );
    test &= fillMomentum( Fitter, particle, name, tMap );
    test &= fillLT( Fitter, particle, name, tMap );
    if ( m_useFullTreeInName ) {
      fillDaughters( Fitter, particle, name, tMap ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    }
  }
  return StatusCode( test );
}
//=============================================================================
// Fill lifetime information for stable daughters
//=============================================================================
StatusCode TupleToolKtautauDTF::fillStableDaughters( const DecayTreeFitter::Fitter& Fitter,
                                                          const LHCb::Particle* P, const std::string& prefix,
                                                          TupleMap& tMap ) const {
  bool test = true;

  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillStableDaughters " << prefix << endmsg;
  const LHCb::Particle::ConstVector& daughters = P->daughtersVector();
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "for id " << P->particleID().pid() << " daughter size is " << daughters.size() << endmsg;
  if ( daughters.empty() ) return StatusCode( test );
  std::set<std::string> usedNames;
  unsigned int          add = 0;
  for ( const auto& particle : daughters ) {
    if ( !particle->isBasicParticle() ) {
      const auto pid     = abs( particle->particleID().pid() );
      const auto pidName = getName( pid );
      auto       name    = prefix + "_" + pidName;
      bool       renamed = false;
      while ( usedNames.find( name ) != usedNames.end() ) { // fix to bug 88702
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Found already name " << name << " trying next " << endmsg;
        renamed = true;
        name    = prefix + "_" + pidName + boost::lexical_cast<std::string>( add );
        ++add;
      }
      if ( renamed )
        Info( "Renaming duplicate to " + name, StatusCode::SUCCESS, 1 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      usedNames.insert( name );
      test &= fillStableDaughters( Fitter, particle, name, tMap );
    } else {
      // const int pid = particle->particleID().pid();
      const auto pid     = abs( particle->particleID().pid() );
      const auto pidName = getName( pid );
      auto       name    = prefix + "_" + pidName;
      bool       renamed = false;
      while ( usedNames.find( name ) != usedNames.end() ) { // fix to bug 88702
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Found already name " << name << " trying next " << endmsg;
        renamed = true;
        name    = prefix + "_" + pidName + "_" + boost::lexical_cast<std::string>( add );
        ++add;
      }
      if ( renamed )
        Info( "Renaming duplicate to " + name, StatusCode::SUCCESS, 1 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      usedNames.insert( name );
      test &= fillTracksMomentum( Fitter, particle, name, tMap );
    }
  }
  return StatusCode( test );
}

//=============================================================================
// Fill updated tracks momentum
//=============================================================================
StatusCode TupleToolKtautauDTF::fillTracksMomentum( const DecayTreeFitter::Fitter& Fitter, const Particle* P,
                                                         const std::string& prefix, TupleMap& tMap ) const {
  bool test = true;

  if ( isVetoed( P->particleID().pid() ) ) { return StatusCode( test ); }

  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "FillTracksMomentum " << prefix << endmsg;

  // Get the fit parameters
  const auto  params   = Fitter.fitParams( P );
  const auto& momentum = params->momentum();

  test &= insert( prefix + "_ID", P->particleID().pid(), tMap );
  test &= insert( prefix + "_PX", momentum.Px(), tMap );
  test &= insert( prefix + "_PY", momentum.Py(), tMap );
  test &= insert( prefix + "_PZ", momentum.Pz(), tMap );
  test &= insert( prefix + "_PE", momentum.E(), tMap ); // MeV

  return StatusCode( test );
}

//=============================================================================
// append data to TupleMap
//=============================================================================
StatusCode TupleToolKtautauDTF::insert( const std::string& leaf, const double val, TupleMap& tMap ) const {
  auto l = tMap.find( leaf );
  if ( l == tMap.end() ) { /// first time this is seen. Create
    std::vector<double> vals;
    vals.push_back( val );
    tMap.insert( std::make_pair( leaf, vals ) );
  } else {
    l->second.push_back( val ); /// append a to vector
  }
  if ( msgLevel( MSG::VERBOSE ) )
    verbose() << "insert " << leaf << " " << val << " size " << l->second.size() << endmsg;
  return StatusCode::SUCCESS;
}

//=============================================================================
// actual filling of the Tuple
//=============================================================================
StatusCode TupleToolKtautauDTF::fillTuple( TupleMap& tMap, Tuples::Tuple& tuple, const std::string& prefix ) {
  bool test = true;

  if ( UNLIKELY( m_firstTupleFill ) ) {
    // Save the list of keys in the given order for future comparisons
    m_firstTupleKeys.clear();
    m_firstTupleKeys.reserve( tMap.size() );
    for ( const auto& i : tMap ) { m_firstTupleKeys.emplace_back( i.first ); }
    // flag having saved the keys
    m_firstTupleFill = false;
  } else {
    // test against the first set of keys
    test = checkTupleKeys( tMap );
  }

  // if OK, save and continue
  if ( test ) {
    for ( const auto& t : tMap ) {
      const auto& leaf = t.first;
      const auto& data = t.second;
      if ( msgLevel( MSG::DEBUG ) )
        debug() << "Filling leaf ``" << leaf << "'' with vector of size " << data.size() << endmsg;
      if ( m_maxPV < data.size() ) Exception( "Seeing data with too many PVs. Have you set MaxPVs?" );
      test &= tuple->farray( leaf, data, prefix + "_nPV", m_maxPV );
    }
  }

  return StatusCode( test );
}

//=============================================================================
// Sort Tracks
//=============================================================================
std::set<const LHCb::Track*> TupleToolKtautauDTF::sortedTracks( const LHCb::VertexBase* vb ) const {
  const LHCb::RecVertex* pv = dynamic_cast<const LHCb::RecVertex*>( vb );
  if ( !pv ) Exception( "Failed to cast PV" );
  std::set<const LHCb::Track*> st;
  for ( const auto& t : pv->tracks() ) { st.insert( t ); }
  return st;
}

//=============================================================================
// Compare PVs, check that one PV's tracks is a subset of the other
//=============================================================================
bool TupleToolKtautauDTF::samePV( const LHCb::VertexBase* vb1, const LHCb::VertexBase* vb2 ) const {
  // exception checking. See bug https://savannah.cern.ch/bugs/?100933
  if ( !vb1 && !vb2 ) {
    Warning( "samePV method called with 2 NULL PVs. "
             "The answer is obviously true, but you may want to check the meaning of the question.",
             StatusCode::SUCCESS, 1 )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    return true;
  } else if ( !vb1 || !vb2 ) {
    Warning( "samePV method called with 1 NULL PV. "
             "The answer is obviously false, but you may want to check the meaning of the question.",
             StatusCode::SUCCESS, 1 )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    return false;
  }

  if ( !( vb1->isPrimary() ) || !( vb2->isPrimary() ) ) {
    Warning( "Non PV VertexBase is being used as PV", StatusCode::SUCCESS, 1 ).ignore();
    return false;
  }

  const auto st1 = sortedTracks( vb1 );
  const auto st2 = sortedTracks( vb2 );

  const bool inc = std::includes( st1.begin(), st1.end(), st2.begin(), st2.end() );
  if ( msgLevel( MSG::VERBOSE ) ) {
    verbose() << "PV 2 of size " << st2.size() << " is ";
    if ( !inc ) verbose() << "not ";
    verbose() << "included in PV 1 of size " << st1.size() << endmsg;
  }
  return inc;
}

//=============================================================================
// get origin vertex
//=============================================================================
std::vector<const VertexBase*> TupleToolKtautauDTF::originVertex( const Particle* mother,
                                                                       const Particle* P ) const {
  std::vector<const VertexBase*> oriVx;
  if ( mother == P ) { // the origin vertex is the primary.
    const auto treeHeadv = m_dva->bestVertex( P );
    if ( treeHeadv ) {
      oriVx.push_back( treeHeadv );
      if ( msgLevel( MSG::VERBOSE ) )
        verbose() << "Pushed back treeHeadv " << treeHeadv << " from " << tesLocation( treeHeadv ) << " at " << treeHeadv->position() << endmsg;
    } else if ( m_constrainToOriginVertex ) {
      Warning( "NULL bestPV while constraining to origin vertex. Fit will be ignored.", StatusCode::SUCCESS, 0 )
          .ignore();
    }
    // all the other ones
    /// @todo : keep only the related ones
    for ( const auto& pv : m_dva->primaryVertices() ) {
      if ( m_storeAnyway || !samePV( pv, treeHeadv ) ) {
        oriVx.push_back( pv );
        if ( msgLevel( MSG::VERBOSE ) )
          verbose() << "Pushed back  pv " << pv << " from " << tesLocation( pv ) << " at " << pv->position() << endmsg;
      }
      if ( oriVx.size() >= m_maxPV ) {
        Warning( "Truncated number of PVs", StatusCode::FAILURE, 0 ).ignore();
        break;
      }
    }
  } else {
    const auto& dau = mother->daughters();
    if ( dau.empty() ) return oriVx;

    for ( const auto& d : dau ) {
      if ( P == d ) {
        oriVx.push_back( mother->endVertex() );
        return oriVx;
      }
    }

    // vertex not yet found, get deeper in the decay:
    for ( const auto& d : dau ) {
      if ( P != d && !d->isBasicParticle() ) {
        oriVx = originVertex( d, P );
        if ( !oriVx.empty() ) {
          return oriVx; // found
        }
      }
    }
  }
  return oriVx;
}

//=============================================================================
// Convert pid number in names
//=============================================================================
std::string TupleToolKtautauDTF::getName( const int id ) const {
  const auto* prop = m_ppSvc->find( LHCb::ParticleID( id ) );
  if ( !prop ) Exception( "Unknown PID" );
  // if (msgLevel(MSG::VERBOSE)) verbose() << "ID " << id << " gets name "
  //                                      << Decays::escape(prop->name()) << endmsg ;
  return Decays::escape( prop->name() );
}

//=============================================================================
// Substitute
//=============================================================================
StatusCode TupleToolKtautauDTF::substitute( LHCb::DecayTree& tree ) {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "Calling substitute" << endmsg;
  const auto substituted = m_substitute->substitute( tree.head() );
  // debugging
  if ( msgLevel( MSG::VERBOSE ) || 0 == substituted ) {
    const auto mp = tree.cloneMap();
    for ( const auto& i : mp ) {
      if ( i.first->particleID().pid() == i.second->particleID().pid() ) {
        info() << "A " << getName( i.first->particleID().pid() ) << " remains unchanged" << endmsg;
      } else {
        info() << "A " << getName( i.first->particleID().pid() ) << " is substituted by a "
               << getName( i.second->particleID().pid() ) << endmsg;
      }
    }
  }
  if ( 0 == substituted ) { return Error( "No particles have been substituted. Check your substitution options." ); }
  return StatusCode::SUCCESS;
}

//=============================================================================
// Check Mass Constraints
//=============================================================================
StatusCode TupleToolKtautauDTF::checkMassConstraints( const LHCb::DecayTree& tree ) {
  if ( !m_first ) return StatusCode::SUCCESS; // do that only once
  m_first       = false;
  const auto mp = tree.cloneMap();
  for ( const auto& m : m_massConstraintsPids ) {
    bool found = false;
    for ( const auto& i : mp ) {
      if ( m.abspid() == i.second->particleID().abspid() ) {
        found = true;
        break;
      }
    }
    if ( found && msgLevel( MSG::VERBOSE ) )
      verbose() << "Constraint " << getName( m.pid() ) << " was found in tree" << endmsg;
    if ( !found ) {
      std::ostringstream mess;
      mess << "Constraint " << getName( m.pid() )
           << " was not found in tree. Check your options. Maybe also the substitution options.";
      return Error( mess.str() );
    }
  }
  return StatusCode::SUCCESS;
}

// Declaration of the Tool Factory
DECLARE_COMPONENT( TupleToolKtautauDTF )
