/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "TupleToolBs2MuMuIso2.h"
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

//-----------------------------------------------------------------------------
// Implementation file for class : RelInfoTrackIsolationBDT2
// Based on TupleToolMuonVariables
//
// 2019-01-23 : Maik BECKER (maik.becker@cern.ch)
//              Matteo RAMA (matteo.rama@cern.ch)
// 2020-07-12 : converted to TupleTool format for usage on DSTs by Biplab Dey
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TupleToolBs2MuMuIso2 )

//=============================================================================
// Constructor
//=============================================================================
TupleToolBs2MuMuIso2::TupleToolBs2MuMuIso2( const std::string& type, const std::string& name, const IInterface* parent )
    : TupleToolBase( type, name, parent ) {
  declareInterface<IParticleTupleTool>( this );

  declareProperty( "Variables", m_variables, "List of variables to store (store all if empty)" );
  declareProperty( "Particles", m_particles = {0, 1},
                   "Positions of the two particles in question within the decay descriptor (default is [0,1])" );
  declareProperty( "PVLocation", m_PVLocation = LHCb::RecVertexLocation::Primary,
                   "Input location of primary vertices" );
  declareProperty( "ParticleLocation", m_particleLocation = "/Event/Phys/StdAllNoPIDsPions/Particles",
                   "Input location of particles in TES" );
  declareProperty( "TrackLocation", m_trackLocation = LHCb::TrackLocation::Default,
                   "Input location of tracks in TES (default is Rec/Track/Best)" );
  declareProperty( "BDTLocation", m_BDTLocation = System::getEnv( "TMVAWEIGHTSROOT" ) + "/data/Bs2MuMu_2017/",
                   // m_BDTLocation = "/eos/lhcb/user/m/maik/public/BDTweights/", // Replace the line above with this
                   // line, as long as a new version of TMVAWeights is not released
                   "Input location of BDT classifier files" );
  declareProperty( "isoBDT_long_xmlFile", m_isoBDT_long_xmlFile = "600_0.5_4_30_0.75_0_8_BDT.weights.xml",
                   "Trained classifier for long track isolation" );
  declareProperty( "isoBDT_velo_xmlFile", m_isoBDT_velo_xmlFile = "600_0.5_4_30_0.75_3_1_BDT.weights.xml",
                   "Trained classifier for VELO track isolation" );
  declareProperty( "newBDT_xmlFile", m_newBDT_xmlFile = "1000_1.00_3_30_0.50_c22_BDT.weights.xml",
                   "Trained B->mumu signal classifier" );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode TupleToolBs2MuMuIso2::initialize() {
  if ( msgLevel( MSG::DEBUG ) ) { debug() << "Starting TupleToolBs2MuMuIso2 tool" << endmsg; }

  StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc() );
  if ( !m_dva ) { return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE ); }

  m_dist = m_dva->distanceCalculator();
  if ( !m_dist ) { return Error( "Couldn't get DistanceCalculator", StatusCode::FAILURE ); }

  if ( m_particles.size() != 2 ) {
    return Error( "Cannot calculate isolations for more than two particle", StatusCode::FAILURE );
  }

  m_variables_nominal.push_back( "TRKISOBDTLONG" );
  m_variables_nominal.push_back( "TRKISOBDTLONG_D1MAX1" );
  m_variables_nominal.push_back( "TRKISOBDTLONG_D1MAX2" );
  m_variables_nominal.push_back( "TRKISOBDTLONG_D1MAX3" );
  m_variables_nominal.push_back( "TRKISOBDTLONG_D2MAX1" );
  m_variables_nominal.push_back( "TRKISOBDTLONG_D2MAX2" );
  m_variables_nominal.push_back( "TRKISOBDTLONG_D2MAX3" );
  m_variables_nominal.push_back( "TRKISOBDTVELO" );
  m_variables_nominal.push_back( "TRKISOBDTVELO_D1MAX1" );
  m_variables_nominal.push_back( "TRKISOBDTVELO_D1MAX2" );
  m_variables_nominal.push_back( "TRKISOBDTVELO_D1MAX3" );
  m_variables_nominal.push_back( "TRKISOBDTVELO_D2MAX1" );
  m_variables_nominal.push_back( "TRKISOBDTVELO_D2MAX2" );
  m_variables_nominal.push_back( "TRKISOBDTVELO_D2MAX3" );
  m_variables_nominal.push_back( "BS2MUMUBDT" );
  m_variables_nominal.push_back( "BS2MUMUBDT_INPUT_MINIPS" );
  m_variables_nominal.push_back( "BS2MUMUBDT_INPUT_IPS" );

  if ( m_variables.size() == 0 ) { // add the nominal variables if none set
    for ( int i = 0; i < (int)m_variables_nominal.size(); i++ ) { m_variables.push_back( m_variables_nominal[i] ); }
  }
  // check that the vars in the list
  for ( int i = 0; i < (int)m_variables.size(); i++ ) {
    if ( find( m_variables_nominal.begin(), m_variables_nominal.end(), m_variables[i] ) == m_variables_nominal.end() ) {
      return Error( "Unknown Bs2MuMuIso2 variable name", StatusCode::FAILURE );
    }
  }

  // Initialize isoBDTs
  m_isoBDT_long_reader = new TMVA::Reader( "!Color:!Silent" );
  m_isoBDT_velo_reader = new TMVA::Reader( "!Color:!Silent" );

  TString isoBDT_long_inputVars[m_size_isoBDT_long_inputVars];
  isoBDT_long_inputVars[0] = "log(trk_ips)";
  isoBDT_long_inputVars[1] = "trk_pvdis1";
  isoBDT_long_inputVars[2] = "trk_svdis1";
  isoBDT_long_inputVars[3] = "log(trk_doca1)";
  isoBDT_long_inputVars[4] = "trk_angle1";
  isoBDT_long_inputVars[5] = "trk_fc1";
  isoBDT_long_inputVars[6] = "trk_dphi1";
  isoBDT_long_inputVars[7] = "abs(trk_deta1)";
  isoBDT_long_inputVars[8] = "trk_pt";

  TString isoBDT_velo_inputVars[m_size_isoBDT_velo_inputVars];
  isoBDT_velo_inputVars[0] = "log(trk_ips)";
  isoBDT_velo_inputVars[1] = "trk_pvdis1";
  isoBDT_velo_inputVars[2] = "trk_svdis1";
  isoBDT_velo_inputVars[3] = "log(trk_doca1)";
  isoBDT_velo_inputVars[4] = "trk_angle1";
  isoBDT_velo_inputVars[5] = "trk_fc1";

  for ( int i = 0; i < m_size_isoBDT_long_inputVars; i++ ) {
    m_isoBDT_long_reader->AddVariable( isoBDT_long_inputVars[i], &m_isoBDT_long_inputValues[i] );
  }

  for ( int i = 0; i < m_size_isoBDT_velo_inputVars; i++ ) {
    m_isoBDT_velo_reader->AddVariable( isoBDT_velo_inputVars[i], &m_isoBDT_velo_inputValues[i] );
  }

  m_isoBDT_method = "method";
  m_isoBDT_long_reader->BookMVA( m_isoBDT_method, m_BDTLocation + m_isoBDT_long_xmlFile );
  m_isoBDT_velo_reader->BookMVA( m_isoBDT_method, m_BDTLocation + m_isoBDT_velo_xmlFile );

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Initialized isolation BDT for long tracks" << endmsg;
    debug() << "\tWeightsFile = " << m_isoBDT_long_xmlFile << endmsg;
    debug() << "Initialized isolation BDT for velo tracks" << endmsg;
    debug() << "\tWeightsFile = " << m_isoBDT_velo_xmlFile << endmsg;
  }
  // End of isoBDT initialization

  // Initialize newBDT
  m_newBDT_reader = new TMVA::Reader( "!Color:!Silent" );

  TString newBDT_inputVars[m_size_newBDT_inputVars];
  newBDT_inputVars[0] = "iso_max2_set8_lt";
  newBDT_inputVars[1] = "iso_max2_set1_velnb";
  newBDT_inputVars[2] = "mu_DeltaR";
  newBDT_inputVars[3] = "mu_MINIPS";
  newBDT_inputVars[4] = "B_ACOSDIRA_OWNPV";
  newBDT_inputVars[5] = "B_IPS_OWNPV";
  newBDT_inputVars[6] = "B_ENDVERTEX_CHI2";

  TString newBDT_inputVars_spect[m_size_newBDT_inputVars_spect];
  newBDT_inputVars_spect[0] = "B_MM";

  for ( int i = 0; i < m_size_newBDT_inputVars; i++ ) {
    m_newBDT_reader->AddVariable( newBDT_inputVars[i], &m_newBDT_inputValues[i] );
  }

  m_newBDT_method = "method";
  m_newBDT_reader->AddSpectator( newBDT_inputVars_spect[0], &m_newBDT_inputValues_spect[0] );
  m_newBDT_reader->BookMVA( m_newBDT_method, m_BDTLocation + m_newBDT_xmlFile );

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Initialized BDT of Bs->mumu analysis" << endmsg;
    debug() << "\tWeightsFile = " << m_newBDT_xmlFile << endmsg;
  }
  // End of newBDT initialization

  return StatusCode::SUCCESS;
}

//=============================================================================
// Fill the tuple
//=============================================================================
StatusCode TupleToolBs2MuMuIso2::fill( const LHCb::Particle*, const LHCb::Particle* part, const std::string& head,
                                       Tuples::Tuple& tuple ) {

  if ( !part ) return StatusCode::FAILURE;
  const std::string prefix = fullName( head );

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Inside Bs2MuMu isol fill loop for particle" << part << " with PID " << part->particleID().pid()
            << endmsg;
  }

  m_isoBDT_long        = -9;
  m_isoBDT_long_d1max1 = -9;
  m_isoBDT_long_d1max2 = -9;
  m_isoBDT_long_d1max3 = -9;
  m_isoBDT_long_d2max1 = -9;
  m_isoBDT_long_d2max2 = -9;
  m_isoBDT_long_d2max3 = -9;
  m_isoBDT_velo        = -9;
  m_isoBDT_velo_d1max1 = -9;
  m_isoBDT_velo_d1max2 = -9;
  m_isoBDT_velo_d1max3 = -9;
  m_isoBDT_velo_d2max1 = -9;
  m_isoBDT_velo_d2max2 = -9;
  m_isoBDT_velo_d2max3 = -9;
  m_newBDT             = -9;

  m_daughters.clear();
  m_tracksToIgnore.clear();

  // Get all daughters
  StatusCode scGetAllDaughters = getAllDaughters( part, m_daughters, m_tracksToIgnore );
  if ( !scGetAllDaughters ) return StatusCode( scGetAllDaughters );

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Detecting " << m_daughters.size() << " daughters (basic + intermediate)" << endmsg;
    debug() << "Detecting " << m_tracksToIgnore.size() << " tracks to ignore" << endmsg;
    debug() << "Detecting daughters in question at positions " << m_particles << endmsg;
  }

  // Check whether particles in question are basic
  bool checkBasicness = true;
  if ( !( m_daughters[m_particles[0]]->isBasicParticle() ) ) checkBasicness = false;
  if ( !( m_daughters[m_particles[1]]->isBasicParticle() ) ) checkBasicness = false;

  bool test = true;
  if ( checkBasicness ) {
    // Calculate RelatedInfo
    StatusCode scCalcTrackIsoBDTLong = calcTrackIsoBDTLong( part );
    if ( !scCalcTrackIsoBDTLong ) return StatusCode( scCalcTrackIsoBDTLong );
    StatusCode scCalcTrackIsoBDTVelo = calcTrackIsoBDTVelo( part );
    if ( !scCalcTrackIsoBDTVelo ) return StatusCode( scCalcTrackIsoBDTVelo );
    StatusCode scCalcBDT = calcBDT( part );
    if ( !scCalcBDT ) return StatusCode( scCalcBDT );

    for ( int i = 0; i < (int)m_variables.size(); i++ ) {
      test &= tuple->column( prefix + "_" + m_variables[i] + "_" + to_string( m_particles[0] ) +
                                 to_string( m_particles[1] ),
                             map_name_var( m_variables[i] ) );
    }
  } else {
    return Error( "At least one daughter in question is not a basic particle", StatusCode::FAILURE );
  }
  return StatusCode( test );
}

double TupleToolBs2MuMuIso2::map_name_var( string isolvarname ) {

  if ( isolvarname == "TRKISOBDTLONG" )
    return m_isoBDT_long;
  else if ( isolvarname == "TRKISOBDTLONG_D1MAX1" )
    return m_isoBDT_long_d1max1;
  else if ( isolvarname == "TRKISOBDTLONG_D1MAX2" )
    return m_isoBDT_long_d1max2;
  else if ( isolvarname == "TRKISOBDTLONG_D1MAX3" )
    return m_isoBDT_long_d1max3;
  else if ( isolvarname == "TRKISOBDTLONG_D2MAX1" )
    return m_isoBDT_long_d2max1;
  else if ( isolvarname == "TRKISOBDTLONG_D2MAX2" )
    return m_isoBDT_long_d2max2;
  else if ( isolvarname == "TRKISOBDTLONG_D2MAX3" )
    return m_isoBDT_long_d2max3;
  else if ( isolvarname == "TRKISOBDTVELO" )
    return m_isoBDT_velo;
  else if ( isolvarname == "TRKISOBDTVELO_D1MAX1" )
    return m_isoBDT_velo_d1max1;
  else if ( isolvarname == "TRKISOBDTVELO_D1MAX2" )
    return m_isoBDT_velo_d1max2;
  else if ( isolvarname == "TRKISOBDTVELO_D1MAX3" )
    return m_isoBDT_velo_d1max3;
  else if ( isolvarname == "TRKISOBDTVELO_D2MAX1" )
    return m_isoBDT_velo_d2max1;
  else if ( isolvarname == "TRKISOBDTVELO_D2MAX2" )
    return m_isoBDT_velo_d2max2;
  else if ( isolvarname == "TRKISOBDTVELO_D2MAX3" )
    return m_isoBDT_velo_d2max3;
  else if ( isolvarname == "BS2MUMUBDT" )
    return m_newBDT;
  else if ( isolvarname == "BS2MUMUBDT_INPUT_MINIPS" )
    return m_newBDT_inputValues[3];
  else if ( isolvarname == "BS2MUMUBDT_INPUT_IPS" )
    return m_newBDT_inputValues[5];
  else {
    error() << "Unknown Bs2MuMuIso2 variable name in map_name_var" << endmsg;
    return -999.;
  }
}

//=============================================================================
// Track Isolation Methods
//=============================================================================
StatusCode TupleToolBs2MuMuIso2::calcTrackIsoBDTLong( const Particle* part ) {
  vector<double> isoBDT_long_d1;
  vector<double> isoBDT_long_d2;

  // Declaration of isoBDT input variables
  double trk_ips   = 0;
  double trk_pvdis = 0;
  double trk_svdis = 0;
  double trk_doca  = 0;
  double trk_angle = 0;
  double trk_fc    = 0;
  double trk_dphi  = 0;
  double trk_deta  = 0;
  double trk_pt    = 0;

  // Get PV and SV
  IDVAlgorithm* dva = Gaudi::Utils::getIDVAlgorithm( contextSvc() );
  if ( !dva ) { return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE ); }
  const VertexBase* pv = dva->bestVertex( part );
  const VertexBase* sv = part->endVertex();

  // Get all particles in the event
  Particles* parts = get<Particles>( m_particleLocation );
  if ( !parts ) { return Error( "Couldn't get particle container", StatusCode::FAILURE ); }

  // Start loop over daughters
  int                                   i = 0;
  int                                   j = 0;
  Particle::ConstVector::const_iterator idau;
  for ( idau = m_daughters.begin(); idau != m_daughters.end(); idau++, i++ ) {
    if ( !( i == m_particles[0] || i == m_particles[1] ) ) continue;

    vector<double> isoBDT_long;

    const Track* daughter = ( *idau )->proto()->track();

    // Start loop over all particles in the event
    Particles::const_iterator iparts;
    for ( iparts = parts->begin(); iparts != parts->end(); iparts++ ) {
      const Track* track = ( *iparts )->proto()->track();

      // If track is not a long track, skip the track
      if ( !( track->type() == 3 ) ) continue;
      // If track is equal to one of the daughter tracks, skip the track
      if ( find( m_tracksToIgnore.begin(), m_tracksToIgnore.end(), track ) != m_tracksToIgnore.end() ) continue;

      // Get momenta and positions of the tracks
      Gaudi::XYZVector daughterMomentum = daughter->momentum();
      Gaudi::XYZPoint  daughterPosition = daughter->position();
      Gaudi::XYZVector trackMomentum    = track->momentum();
      Gaudi::XYZPoint  trackPosition    = track->position();

      // Compute isoBDT input variables
      StatusCode scCalcIPToAnyPV = calcIPToAnyPV( *iparts, trk_ips );
      if ( !scCalcIPToAnyPV ) return StatusCode( scCalcIPToAnyPV );

      Gaudi::XYZPoint vertDaughterTrack;
      StatusCode      scGetPerpFeet =
          getPerpFeet( daughterMomentum, daughterPosition, trackMomentum, trackPosition, vertDaughterTrack );
      if ( !scGetPerpFeet ) return StatusCode( scGetPerpFeet );
      trk_pvdis = calcVertexDist( vertDaughterTrack, pv );
      trk_svdis = calcVertexDist( vertDaughterTrack, sv );

      StatusCode scDoca = m_dist->distance( track, daughter, trk_doca );
      if ( !scDoca ) return StatusCode( scDoca );

      StatusCode scInCone = inCone( trackPosition, trackMomentum, daughterPosition, daughterMomentum, trk_doca );
      if ( !scInCone ) return StatusCode( scInCone );

      StatusCode scCalcFC = calcFC( daughterMomentum, trackMomentum, vertDaughterTrack, pv, trk_fc );
      if ( !scCalcFC ) return StatusCode( scCalcFC );

      trk_angle = calcEnclosedAngle( trackMomentum, daughterMomentum );
      trk_dphi  = fabs( trackMomentum.Phi() - daughterMomentum.Phi() );
      trk_deta  = fabs( trackMomentum.Eta() - daughterMomentum.Eta() );
      trk_pt    = track->pt();

      m_isoBDT_long_inputValues[0] = log( trk_ips );
      m_isoBDT_long_inputValues[1] = trk_pvdis;
      m_isoBDT_long_inputValues[2] = trk_svdis;
      m_isoBDT_long_inputValues[3] = log( trk_doca );
      m_isoBDT_long_inputValues[4] = trk_angle;
      m_isoBDT_long_inputValues[5] = trk_fc;
      m_isoBDT_long_inputValues[6] = trk_dphi;
      m_isoBDT_long_inputValues[7] = trk_deta;
      m_isoBDT_long_inputValues[8] = trk_pt;

      // Evaluate the isoBDT
      isoBDT_long.push_back( m_isoBDT_long_reader->EvaluateMVA( m_isoBDT_method ) );
    }

    j++;
    if ( j == 1 ) { isoBDT_long_d1 = isoBDT_long; }
    if ( j == 2 ) { isoBDT_long_d2 = isoBDT_long; }
    if ( j >= 3 ) { return Error( "Cannot calculate isolations for more than two particle", StatusCode::FAILURE ); }
  }

  // Compute the isolation variables
  computeIso( isoBDT_long_d1, m_isoBDT_long_d1max1, m_isoBDT_long_d1max2, m_isoBDT_long_d1max3 );

  computeIso( isoBDT_long_d2, m_isoBDT_long_d2max1, m_isoBDT_long_d2max2, m_isoBDT_long_d2max3 );

  m_isoBDT_long = m_isoBDT_long_d1max1 + m_isoBDT_long_d2max1;

  return StatusCode::SUCCESS;
}

StatusCode TupleToolBs2MuMuIso2::calcTrackIsoBDTVelo( const Particle* part ) {
  vector<double> isoBDT_velo_d1;
  vector<double> isoBDT_velo_d2;

  // Declaration of isoBDT input variables
  double trk_ips   = 0;
  double trk_pvdis = 0;
  double trk_svdis = 0;
  double trk_doca  = 0;
  double trk_angle = 0;
  double trk_fc    = 0;

  // Get all tracks in the event
  Tracks* tracks = get<Tracks>( m_trackLocation );
  if ( !tracks ) { return Error( "Couldn't get track container", StatusCode::FAILURE ); }

  // Start loop over daughters
  int                                   i = 0;
  int                                   j = 0;
  Particle::ConstVector::const_iterator idau;
  for ( idau = m_daughters.begin(); idau != m_daughters.end(); idau++, i++ ) {
    if ( !( i == m_particles[0] || i == m_particles[1] ) ) continue;

    vector<double> isoBDT_velo;

    const Track* daughter = ( *idau )->proto()->track();

    // Start loop over all tracks in the tracks
    Tracks::const_iterator itracks;
    for ( itracks = tracks->begin(); itracks != tracks->end(); itracks++ ) {
      const Track* track = ( *itracks );

      // If track is not a Velo or VeloR track, skip the track
      if ( track->type() != 1 && track->type() != 2 ) continue;
      // If track is backward, skip the track
      if ( track->checkFlag( Track::Backward ) ) continue;

      // Get momenta and positions of the tracks
      Gaudi::XYZVector daughterMomentum = daughter->momentum();
      Gaudi::XYZPoint  daughterPosition = daughter->position();
      Gaudi::XYZVector trackMomentum    = track->momentum();
      Gaudi::XYZPoint  trackPosition    = track->position();

      // Compute isoBDT input variables
      StatusCode scGetTrackMinIP = getTrackMinIP( track, trk_ips );
      if ( !scGetTrackMinIP ) return StatusCode( scGetTrackMinIP );

      StatusCode scDoca = m_dist->distance( track, daughter, trk_doca );
      if ( !scDoca ) return StatusCode( scDoca );

      StatusCode scInCone = inCone( trackPosition, trackMomentum, daughterPosition, daughterMomentum, trk_doca );
      if ( !scInCone ) return StatusCode( scInCone );

      trk_angle = calcEnclosedAngle( trackMomentum, daughterMomentum );

      get_fc_pvdis_svdis( part, daughter, track, trk_fc, trk_pvdis, trk_svdis ); // Cleanup!

      m_isoBDT_velo_inputValues[0] = log( trk_ips );
      m_isoBDT_velo_inputValues[1] = trk_pvdis;
      m_isoBDT_velo_inputValues[2] = trk_svdis;
      m_isoBDT_velo_inputValues[3] = log( trk_doca );
      m_isoBDT_velo_inputValues[4] = trk_angle;
      m_isoBDT_velo_inputValues[5] = trk_fc;

      // Evaluate the isoBDT
      isoBDT_velo.push_back( m_isoBDT_velo_reader->EvaluateMVA( m_isoBDT_method ) );
    }

    j++;
    if ( j == 1 ) { isoBDT_velo_d1 = isoBDT_velo; }
    if ( j == 2 ) { isoBDT_velo_d2 = isoBDT_velo; }
    if ( j >= 3 ) { return Error( "Cannot calculate isolations for more than two particle", StatusCode::FAILURE ); }
  }

  // Compute the isolation variable
  computeIso( isoBDT_velo_d1, m_isoBDT_velo_d1max1, m_isoBDT_velo_d1max2, m_isoBDT_velo_d1max3 );

  computeIso( isoBDT_velo_d2, m_isoBDT_velo_d2max1, m_isoBDT_velo_d2max2, m_isoBDT_velo_d2max3 );

  m_isoBDT_velo = m_isoBDT_velo_d1max1 + m_isoBDT_velo_d2max1;

  return StatusCode::SUCCESS;
}

StatusCode TupleToolBs2MuMuIso2::calcBDT( const Particle* part ) {
  // Declaration of newBDT input variables
  double part_mu_AbsPhi        = 0;
  double part_mu_DeltaPhi      = 0;
  double part_mu_DeltaEta      = 0;
  double part_mu_DeltaR        = 0;
  double part_mu_MINIPS        = 0;
  double part_B_ACOSDIRA_OWNPV = 0;
  double part_B_IPS_OWNPV      = 0;
  double part_B_ENDVERTEX_CHI2 = 0;

  // Get PV
  IDVAlgorithm* dva = Gaudi::Utils::getIDVAlgorithm( contextSvc() );
  if ( !dva ) { return Error( "Couldn't get parent DVAlgorithm", StatusCode::FAILURE ); }
  const VertexBase* pv = dva->bestVertex( part );

  // Get tracks of both daughters in question
  const Track* d1 = m_daughters[m_particles[0]]->proto()->track();
  const Track* d2 = m_daughters[m_particles[1]]->proto()->track();

  // Compute newBDT input variables
  const LoKi::Cuts::BPVIPCHI2 BPVIPCHI2;
  const LoKi::Types::Fun      DIRA = LoKi::Cuts::DIRA( pv );

  part_mu_AbsPhi   = TMath::Abs( d1->momentum().Phi() - d2->momentum().Phi() );
  part_mu_DeltaPhi = TMath::Abs( TMath::Abs( TMath::Abs( part_mu_AbsPhi ) - TMath::ACos( -1 ) ) - TMath::ACos( -1 ) );
  part_mu_DeltaEta = TMath::Abs( d1->momentum().Eta() - d2->momentum().Eta() );
  part_mu_DeltaR   = TMath::Sqrt( part_mu_DeltaPhi * part_mu_DeltaPhi + part_mu_DeltaEta * part_mu_DeltaEta );
  part_mu_MINIPS =
      TMath::Sqrt( TMath::Min( BPVIPCHI2( m_daughters[m_particles[0]] ), BPVIPCHI2( m_daughters[m_particles[1]] ) ) );
  part_B_ACOSDIRA_OWNPV = TMath::ACos( DIRA( part ) );
  part_B_IPS_OWNPV      = TMath::Sqrt( BPVIPCHI2( part ) );
  part_B_ENDVERTEX_CHI2 = part->endVertex()->chi2();

  m_newBDT_inputValues[0] = m_isoBDT_long;
  m_newBDT_inputValues[1] = m_isoBDT_velo;
  m_newBDT_inputValues[2] = part_mu_DeltaR;
  m_newBDT_inputValues[3] = part_mu_MINIPS;
  m_newBDT_inputValues[4] = part_B_ACOSDIRA_OWNPV;
  m_newBDT_inputValues[5] = part_B_IPS_OWNPV;
  m_newBDT_inputValues[6] = part_B_ENDVERTEX_CHI2;

  // Compute newBDT spectator variable
  m_newBDT_inputValues_spect[0] = part->measuredMass();

  // Evaluate newBDT
  m_newBDT = m_newBDT_reader->EvaluateMVA( m_newBDT_method );

  return StatusCode::SUCCESS;
}

//=============================================================================
// Helper Functions
//=============================================================================
StatusCode TupleToolBs2MuMuIso2::getAllDaughters( const Particle* part, Particle::ConstVector& daughters,
                                                  Track::ConstVector& tracksToIgnore ) {
  const Particle::ConstVector daughters_temp = part->daughtersVector();

  bool                                  checkIntegrity = true;
  Particle::ConstVector::const_iterator idau;
  for ( idau = daughters_temp.begin(); idau != daughters_temp.end(); idau++ ) {
    if ( ( *idau )->proto() == NULL && ( *idau )->isBasicParticle() ) {
      checkIntegrity = false;
    } else {
      daughters.push_back( ( *idau ) );

      if ( ( *idau )->isBasicParticle() ) {
        tracksToIgnore.push_back( ( *idau )->proto()->track() );
      } else {
        StatusCode scGetAllDaughters = getAllDaughters( ( *idau ), daughters, tracksToIgnore );
        if ( !scGetAllDaughters ) return StatusCode( scGetAllDaughters );
      }
    }
  }
  if ( !checkIntegrity ) { return Error( "At least one basic daughter particle is NULL", StatusCode::FAILURE ); }

  return StatusCode::SUCCESS;
}

StatusCode TupleToolBs2MuMuIso2::calcIPToAnyPV( const Particle* part, double& ips ) {
  ips           = 100000;
  double ip     = 0.0;
  double ipchi2 = 0.0;

  // Get all PVs
  RecVertex::Container* verts = get<RecVertex::Container>( m_PVLocation );
  if ( !verts ) { return Error( "Couldn't get vertex container", StatusCode::FAILURE ); }

  // Start loop over all PVs in the event
  RecVertex::Container::const_iterator iverts;
  for ( iverts = verts->begin(); iverts != verts->end(); iverts++ ) {
    m_dist->distance( part, ( *iverts ), ip, ipchi2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    if ( ipchi2 < ips ) ips = ipchi2;
  }
  ips = sqrt( ips );

  return StatusCode::SUCCESS;
}

StatusCode TupleToolBs2MuMuIso2::getPerpFeet( Gaudi::XYZVector daughterMom, Gaudi::XYZPoint daughterPos,
                                              Gaudi::XYZVector trackMom, Gaudi::XYZPoint trackPos,
                                              Gaudi::XYZPoint& vert ) {
  const Gaudi::XYZVector diffPos( trackPos - daughterPos );
  const Gaudi::XYZVector unitTrackMom( trackMom.unit() );
  const Gaudi::XYZVector unitDaughterMom( daughterMom.unit() );

  const double d_diffTrack        = diffPos.Dot( unitTrackMom );
  const double d_diffDaughter     = diffPos.Dot( unitDaughterMom );
  const double d_daughterTrack    = unitDaughterMom.Dot( unitTrackMom );
  const double d_trackTrack       = unitTrackMom.Dot( unitTrackMom );
  const double d_daughterDaughter = unitDaughterMom.Dot( unitDaughterMom );
  const double den                = d_daughterTrack * d_daughterTrack - d_daughterDaughter * d_trackTrack;
  if ( fabs( den ) < 1E-27 ) { return Error( "Vertex calculation (perpendicular feet) failed", StatusCode::FAILURE ); }

  const double num = d_diffTrack * d_daughterTrack - d_trackTrack * d_diffDaughter;
  const double mu2 = num / den;
  const double mu1 = ( mu2 * d_daughterTrack - d_diffTrack ) / d_trackTrack;

  const Gaudi::XYZPoint perpFootTrack    = trackPos + unitTrackMom * mu1;
  const Gaudi::XYZPoint perpFootDaughter = daughterPos + unitDaughterMom * mu2;
  vert                                   = ( perpFootTrack + 0.5 * ( perpFootDaughter - perpFootTrack ) );

  return StatusCode::SUCCESS;
}

double TupleToolBs2MuMuIso2::calcVertexDist( Gaudi::XYZPoint track, const VertexBase* vert ) {
  if ( vert == NULL ) return 0;
  const Gaudi::XYZPoint& vtx = vert->position();
  return ( ( track.z() - vtx.z() ) / fabs( track.z() - vtx.z() ) * ( track - vtx ).R() );
}

StatusCode TupleToolBs2MuMuIso2::inCone( Gaudi::XYZPoint o1, Gaudi::XYZVector p1, Gaudi::XYZPoint o2,
                                         Gaudi::XYZVector p2, double& doca ) {
  Gaudi::XYZPoint vert( 0.0, 0.0, 0.0 );
  Gaudi::XYZPoint close1;
  Gaudi::XYZPoint close2;

  StatusCode scClosestPoint = closestPoint( o1, p1, o2, p2, close1, close2, vert );
  if ( !scClosestPoint ) return StatusCode( scClosestPoint );

  doca = ( close1 - close2 ).R();

  return StatusCode::SUCCESS;
}

StatusCode TupleToolBs2MuMuIso2::closestPoint( Gaudi::XYZPoint o1, Gaudi::XYZVector p1, Gaudi::XYZPoint o2,
                                               Gaudi::XYZVector p2, Gaudi::XYZPoint& close1, Gaudi::XYZPoint& close2,
                                               Gaudi::XYZPoint& vert ) {
  const Gaudi::XYZVector v0( o1 - o2 );
  const Gaudi::XYZVector v1( p1.unit() );
  const Gaudi::XYZVector v2( p2.unit() );

  const double d02 = v0.Dot( v2 );
  const double d21 = v2.Dot( v1 );
  const double d01 = v0.Dot( v1 );
  const double d22 = v2.Dot( v2 );
  const double d11 = v1.Dot( v1 );
  const double den = d11 * d22 - d21 * d21;
  if ( fabs( den ) < 1E-27 ) { return Error( "Calculation of closest point failed", StatusCode::FAILURE ); }

  const double num = d02 * d21 - d01 * d22;
  const double mu1 = num / den;
  const double mu2 = ( d02 + d21 * mu1 ) / d22;
  close1           = o1 + v1 * mu1;
  close2           = o2 + v2 * mu2;
  vert             = ( close1 + 0.5 * ( close2 - close1 ) );

  return StatusCode::SUCCESS;
}

double TupleToolBs2MuMuIso2::calcEnclosedAngle( Gaudi::XYZVector p1, Gaudi::XYZVector p2 ) {
  const double den      = p1.R() * p2.R();
  const double cosAngle = p1.Dot( p2 ) / den;
  const double angle    = acos( fabs( cosAngle ) );
  return ( cosAngle < 0 ? ROOT::Math::Pi() - angle : angle );
}

StatusCode TupleToolBs2MuMuIso2::calcFC( Gaudi::XYZVector daughterMom, Gaudi::XYZVector trackMom,
                                         Gaudi::XYZPoint vertDaughterTrack, const VertexBase* PV, double& fc ) {
  const Gaudi::XYZPoint& pv = PV->position();
  const Gaudi::XYZVector p1( trackMom + daughterMom );
  const Gaudi::XYZVector p2( vertDaughterTrack - pv );

  const double angle = sin( calcEnclosedAngle( p1, p2 ) );
  const double num   = p1.R() * angle;
  const double den   = p1.R() * angle + trackMom.Rho() + daughterMom.Rho();
  if ( fabs( den ) < 1E-27 ) { return Error( "Calculation of fc failed", StatusCode::FAILURE ); }

  fc = num / den;

  return StatusCode::SUCCESS;
}

void TupleToolBs2MuMuIso2::computeIso( vector<double> iso_vec, double& iso_max1, double& iso_max2, double& iso_max3 ) {
  iso_max1 = -2;
  iso_max2 = -3;
  iso_max3 = -4;
  double iso_temp;

  for ( int i = 0; i < (int)iso_vec.size(); i++ ) {
    iso_temp = iso_vec[i];

    if ( iso_temp >= iso_max1 ) {
      iso_max3 = iso_max2;
      iso_max2 = iso_max1;
      iso_max1 = iso_temp;
    } else if ( iso_temp >= iso_max2 ) {
      iso_max3 = iso_max2;
      iso_max2 = iso_temp;
    } else if ( iso_temp >= iso_max3 ) {
      iso_max3 = iso_temp;
    }
  }

  return;
}

StatusCode TupleToolBs2MuMuIso2::getTrackMinIP( const Track* track, double& ips ) {
  ips       = 1000000;
  double ip = 1000000;

  Gaudi::XYZPoint     trackPos;
  Gaudi::XYZVector    trackMom;
  Gaudi::SymMatrix6x6 trackCov;
  track->positionAndMomentum( trackPos, trackMom, trackCov );

  Gaudi::XYZVector    trackSlope;
  Gaudi::SymMatrix3x3 trackSlopeCov;
  track->slopes( trackSlope, trackSlopeCov );

  // Get all PVs
  RecVertex::Container* verts = get<RecVertex::Container>( m_PVLocation );
  if ( !verts ) { return Error( "Couldn't get vertex container", StatusCode::FAILURE ); }

  // Start loop over all PVs in the event
  RecVertex::Container::const_iterator iverts;
  for ( iverts = verts->begin(); iverts != verts->end(); iverts++ ) {
    Gaudi::XYZPoint            pvPos        = ( *iverts )->position();
    const Gaudi::SymMatrix3x3& pvCov        = ( *iverts )->covMatrix();
    Gaudi::XYZVector           unitTrackMom = trackMom.Unit();
    Gaudi::XYZVector           docaVec      = unitTrackMom.Cross( pvPos - trackPos );
    double                     temp_ip      = docaVec.R();
    double                     vertChi2     = -1;
    getTrackIPSwrtPV( trackMom, trackPos, trackCov, trackSlopeCov, pvPos, pvCov, vertChi2 );
    double temp_ips = sqrt( vertChi2 );

    if ( temp_ip < ip ) ip = temp_ip;
    if ( temp_ips < ips ) ips = temp_ips;
  }

  return StatusCode::SUCCESS;
}

void TupleToolBs2MuMuIso2::getTrackIPSwrtPV( Gaudi::XYZVector& mom_trk, Gaudi::XYZPoint& pos_trk,
                                             const Gaudi::SymMatrix6x6& cov_trk, Gaudi::SymMatrix3x3 slope_cov_trk,
                                             Gaudi::XYZPoint& pos_pv, const Gaudi::SymMatrix3x3& cov_pv,
                                             double& chi2 ) {
  Gaudi::XYZVector vtx_pos;
  TMatrixD         vtx_Cov( 3, 3 );

  TVector3 x0_trk( pos_trk.x(), pos_trk.y(), pos_trk.z() );
  double   x_pv_arr[3] = {pos_pv.x(), pos_pv.y(), pos_pv.z()};
  TVectorD x_pv( 3, x_pv_arr );

  TMatrixDSym covM_pv( 3 );
  covM_pv[0][0] = cov_pv( 0, 0 );
  covM_pv[0][1] = cov_pv( 0, 1 );
  covM_pv[0][2] = cov_pv( 0, 2 );
  covM_pv[1][0] = cov_pv( 1, 0 );
  covM_pv[1][1] = cov_pv( 1, 1 );
  covM_pv[1][2] = cov_pv( 1, 2 );
  covM_pv[2][0] = cov_pv( 2, 0 );
  covM_pv[2][1] = cov_pv( 2, 1 );
  covM_pv[2][2] = cov_pv( 2, 2 );

  TMatrixDSym covMinv_pv( covM_pv );
  covMinv_pv.Invert();

  TVectorD p_trk( 2 );
  TVectorD qm_trk( 4 );

  p_trk[0]  = mom_trk.x() / mom_trk.z();
  p_trk[1]  = mom_trk.y() / mom_trk.z();
  qm_trk[0] = pos_trk.x();
  qm_trk[1] = pos_trk.y();
  qm_trk[2] = p_trk[0];
  qm_trk[3] = p_trk[1];

  TVector3 v0( x_pv[0], x_pv[1], x_pv[2] );

  TVectorD x_v0( x_pv ); // Store info for the chi2 evaluation below

  // define F(v0,p0_i) as F(v,p_i) evaluated in v=v0, p_i=pm_i // slope
  double F_v0p0_trk_arr[4] = {v0[0] + p_trk[0] * ( x0_trk[2] - v0[2] ), v0[1] + p_trk[1] * ( x0_trk[2] - v0[2] ),
                              p_trk[0], p_trk[1]};

  TVectorD F_v0p0_trk( 4, F_v0p0_trk_arr );

  // define the variation dqm_trk
  TVectorD dqm_trk = qm_trk - F_v0p0_trk;

  // define matrix dF/dv=D. dF/dv is defined in (v0,p0_i)
  TMatrixD D_trk = get_dFdv( p_trk );

  // define matrix dF/dp=E. dF/dp is defined in (v0,p0_i)
  TMatrixD E_trk = get_dFdp( v0, x0_trk );

  // covariance matrix for qm vector
  TMatrixDSym W_trk = get_trkCov( pos_trk, mom_trk, cov_trk, slope_cov_trk );

  W_trk.Invert();

  //
  // Now we have all ingredients to apply Billoir
  //

  TMatrixD D_trk_t( D_trk );
  D_trk_t.T();

  TMatrixD A = D_trk_t * W_trk * D_trk + covMinv_pv;

  TMatrixD B_trk = D_trk_t * W_trk * E_trk;
  TMatrixD B_trk_t( B_trk );
  B_trk_t.T();
  TMatrixD E_trk_t( E_trk );
  E_trk_t.T();
  TMatrixD C_trk = E_trk_t * W_trk * E_trk;

  TVectorD T     = D_trk_t * W_trk * dqm_trk + covMinv_pv * ( x_pv - x_v0 );
  TVectorD U_trk = E_trk_t * W_trk * dqm_trk;

  TMatrixD Cinv_trk( C_trk );
  Cinv_trk.Invert();

  TMatrixD Cov_V = A - B_trk * Cinv_trk * B_trk_t;
  Cov_V.Invert();

  TVectorD vec = T - B_trk * Cinv_trk * U_trk;
  TVectorD dV  = Cov_V * vec;
  vtx_pos.SetXYZ( v0[0] + dV[0], v0[1] + dV[1], v0[2] + dV[2] );
  vtx_Cov = Cov_V;

  TVectorD dp_trk = Cinv_trk * U_trk - Cinv_trk * B_trk_t * dV;

  TVectorD deltaq_trk = qm_trk - ( F_v0p0_trk + D_trk * dV + E_trk * dp_trk );
  chi2                = deltaq_trk * ( W_trk * deltaq_trk ) + dV * ( covMinv_pv * dV );

  return;
}

void TupleToolBs2MuMuIso2::get_fc_pvdis_svdis( const LHCb::Particle* B, const LHCb::Track* mutrack,
                                               const LHCb::Track* track, double& fc, double& pvdis, double& svdis ) {
  Gaudi::XYZPoint  _x_mu = mutrack->position();
  Gaudi::XYZVector _p_mu = mutrack->momentum();

  TVector3 x_mu( _x_mu.x(), _x_mu.y(), _x_mu.z() );
  TVector3 p_mu( _p_mu.x(), _p_mu.y(), _p_mu.z() );

  Gaudi::XYZPoint  _x_trk = track->position();
  Gaudi::XYZVector _p_trk = track->momentum();

  TVector3 x_trk( _x_trk.x(), _x_trk.y(), _x_trk.z() );
  TVector3 p_trk( _p_trk.x(), _p_trk.y(), _p_trk.z() );

  TVector3 u1 = p_mu.Unit();
  TVector3 u2 = p_trk.Unit();
  double   a  = u1.Dot( u2 );

  TVector3 r0 = x_mu - x_trk;
  double   b  = r0.Dot( u1 );
  double   c  = r0.Dot( u2 );
  double   t1 = ( a * c - b ) / ( 1 - a * a );
  double   t2 = ( -a * b + c ) / ( 1 - a * a );

  TVector3 d1  = x_mu + u1 * t1;    // POCA 1
  TVector3 d2  = x_trk + u2 * t2;   // POCA 2
  TVector3 vtx = 0.5 * ( d1 + d2 ); // NB: this approximate estimate of the vertex is used to be consistent
                                    // with the def of the 'classic' track iso (used in 2013)

  const LHCb::VertexBase* bestpv = get_PV_BminIPS( B );

  if ( bestpv == NULL ) { // no proper PV was found. Happens very rarely, but happens
    error() << "TupleToolMuonVariables::get_fc_pvdis_svdis(): Setting pvdis, svdis and fc to nan" << endmsg;
    fc    = std::nan( "" );
    pvdis = std::nan( "" );
    svdis = std::nan( "" );
    return;
  } else {
    TVector3 pv_pos( bestpv->position().x(), bestpv->position().y(), bestpv->position().z() );

    // compute fc
    TVector3 rv   = vtx - pv_pos;
    double   pt   = p_mu.Pt() + p_trk.Pt();
    TVector3 ptot = p_mu + p_trk;
    double   temp = rv.Angle( ptot );
    double   num  = ptot.Mag() * sin( temp );
    double   den  = ( num + pt );
    fc            = num / den;

    // compute psvid
    int sign_pv = vtx.z() - pv_pos.z() >= 0 ? 1 : -1;
    pvdis       = sign_pv * ( vtx - pv_pos ).Mag();

    // compute svdis
    TVector3 b_pos( B->endVertex()->position().x(), B->endVertex()->position().y(), B->endVertex()->position().z() );
    int      sign_sv = vtx.z() - b_pos.z() >= 0 ? 1 : -1;
    svdis            = sign_sv * ( vtx - b_pos ).Mag();

    return;
  }
}

const LHCb::VertexBase* TupleToolBs2MuMuIso2::get_PV_BminIPS( const LHCb::Particle* B ) {

  double b_ips = 1e6;

  double                  ip, ips, fl, fl0, ctau;
  const LHCb::VertexBase* bestpv = NULL;

  // Get the list of primary vertices and store their vtx info in ntuple
  LHCb::RecVertex::Container* verts = NULL;
  if ( exist<LHCb::RecVertex::Container>( m_PVLocation ) ) verts = get<LHCb::RecVertex::Container>( m_PVLocation );
  LHCb::RecVertex::Container::const_iterator iv;

  for ( iv = verts->begin(); iv != verts->end(); iv++ ) {
    const LHCb::VertexBase* aPV = ( *iv );

    get_BIPwrtPV( B, aPV, ip, ips, fl, fl0, ctau );

    if ( ips < b_ips ) {
      bestpv = aPV;
      b_ips  = ips;
    }
  } // end i<PV_num
  if ( bestpv == NULL ) { error() << "TupleToolMuonVariables::get_PV_BminIPS(): No proper PV found!" << endmsg; }
  return bestpv;
}

// Make fit of intersersection between the B direction and PV, minimizing the chi2. Extract the new vtx(B) position,
// slope(B) and the B f.l.
// Compute the IP between the B direction and the PV (before refit).Compute IPS as the chi2 value.
// Compute ctau from the f.l. (see notes below)
//
// Note: B_FD_OWNPV is defined as |x_b-x_pv|, i.e. the equivalent of q0[5] in the code below.
// It would be more correct to define it as
//       (q0+dq)[5], which corresponds to the fl extracted from the 'fit'. However, to be consistent with what DaVinci
//       does, get_BIPwrtPV returns also fl0=|x_b-x_pv|. In the standard ntuples B_CTAU_ps is defined exaclty as
//       B_FD_OWNPV/(c*B_P/B_M)*1e+12. On the other hand, B_TAU is similar to fl/(c*B_P/B_M)*1e+12, even though it is
//       not exactly so (it differs by about 0.1%).
void TupleToolBs2MuMuIso2::get_BIPwrtPV( const LHCb::Particle* B, const LHCb::VertexBase* PV, double& ip, double& ips,
                                         double& fl, double& fl0, double& ctau ) {

  // The vector q of the params to fit is given by x,y,z,slx,sly and f.l.of B
  TVector3 xb_m( B->endVertex()->position().x(), B->endVertex()->position().y(), B->endVertex()->position().z() );
  TVector3 pb_m( B->momentum().x(), B->momentum().y(), B->momentum().z() );

  TVector3 xpv_m( PV->position().x(), PV->position().y(), PV->position().z() );

  fl0                = ( xb_m - xpv_m ).Mag();
  double   q0_arr[6] = {xb_m[0], xb_m[1], xb_m[2], pb_m[0] / pb_m[2], pb_m[1] / pb_m[2], fl0};
  TVectorD q0( 6, q0_arr );

  // take the cov mat of PV and its inverse
  TMatrixDSym                covM_pv( 3 );
  const Gaudi::SymMatrix3x3& cov_pv = PV->covMatrix();
  covM_pv[0][0]                     = cov_pv( 0, 0 );
  covM_pv[0][1]                     = cov_pv( 0, 1 );
  covM_pv[0][2]                     = cov_pv( 0, 2 );
  covM_pv[1][0]                     = cov_pv( 1, 0 );
  covM_pv[1][1]                     = cov_pv( 1, 1 );
  covM_pv[1][2]                     = cov_pv( 1, 2 );
  covM_pv[2][0]                     = cov_pv( 2, 0 );
  covM_pv[2][1]                     = cov_pv( 2, 1 );
  covM_pv[2][2]                     = cov_pv( 2, 2 );
  TMatrixDSym W_pv( covM_pv );
  W_pv.Invert();

  // introduce qm_pv, q0_pv and dq_pv=qm_pv-q0_pv
  double   qm_pv_arr[3] = {PV->position().x(), PV->position().y(), PV->position().z()};
  TVectorD qm_pv( 3, qm_pv_arr );
  double   root         = sqrt( 1 + q0[3] * q0[3] + q0[4] * q0[4] );
  double   q0_pv_arr[3] = {q0[0] - q0[3] / root * q0[5], q0[1] - q0[4] / root * q0[5], q0[2] - 1 / root * q0[5]};
  TVectorD q0_pv( 3, q0_pv_arr );
  TVectorD dq_pv = qm_pv - q0_pv;

  // take the matrix D such that q_pv = q0_pv+D*dq
  TMatrixD D = get_dFdq( q0 );
  TMatrixD D_t( D );
  D_t.T();

  // take the cov mat x,y,z,slx,sly of the B
  TMatrixDSym W_B = get_BCov( B );
  W_B.Invert();

  // define W_B_prime as W_B with the addition of one row and one col (used later)
  TMatrixDSym W_B_prime( 6 );
  for ( int i = 0; i < 6; i++ )
    for ( int j = 0; j < 6; j++ ) {
      if ( i < 5 && j < 5 )
        W_B_prime[i][j] = W_B[i][j];
      else
        W_B_prime[i][j] = 0;
    }

  TMatrixD Cov_dq = D_t * W_pv * D + W_B_prime;
  Cov_dq.Invert();
  TVectorD dq = Cov_dq * ( D_t * W_pv * dq_pv );

  // compute chi2 relative to PV
  TVectorD deltaq_pv = qm_pv - q0_pv - D * dq;
  double   chi2_pv   = deltaq_pv * ( W_pv * deltaq_pv );

  // compute chi2 relative to B (it exploits the fact that q0_b was chosen so that qm_b-q0_b=0)
  double chi2_b = dq * ( W_B_prime * dq );
  double chi2   = chi2_b + chi2_pv;
  ips           = sqrt( chi2 );

  // compute the IP. NB: Do NOT use the _refitted_ vertex and slopes of B to be consistent with DaVinci
  TVector3 v        = pb_m.Unit();
  TVector3 doca_vec = v.Cross( xpv_m - xb_m );
  ip                = doca_vec.Mag();

  // compute the flight length
  fl = q0[5] + dq[5];

  //  double B_P = sqrt(B->pt()*B->pt() + B->momentum().z()*B->momentum().z()); // boh
  double cbg = 2.99792458e+11 * B->p() / B->measuredMass();

  ctau = fl / cbg * 1e+12;

  return;
} // end get_BIPwrtPV

// return dF/dq (3x6 matrix used in get_BIPwrtPV)
TMatrixD TupleToolBs2MuMuIso2::get_dFdq( TVectorD q ) {

  // q={xb,yb,zb,slx_b,sly_b,fl_b}, where xb,yb,zb, slx_b and sly_b are the measured values. fl_b=|xvec_b-xvec_pv|
  double   slx = q[3];
  double   sly = q[4];
  double   fl  = q[5];
  TMatrixD out( 3, 6 );
  out[0][0] = out[1][1] = out[2][2] = 1;
  out[0][1] = out[0][2] = out[1][0] = out[1][2] = out[2][0] = out[2][1] = 0;
  double root1                                                          = sqrt( 1 + q[3] * q[3] + q[4] * q[4] );
  double root2                                                          = pow( root1, 3 );
  out[0][3]                                                             = -( 1 + sly * sly ) / root2 * fl;
  out[0][4]                                                             = slx * sly / root2 * fl;
  out[0][5]                                                             = -slx / root1;
  out[1][3]                                                             = slx * sly / root2 * fl;
  out[1][4]                                                             = -( 1 + slx * slx ) / root2 * fl;
  out[1][5]                                                             = -sly / root1;
  out[2][3]                                                             = slx / root2;
  out[2][4]                                                             = sly / root2;
  out[2][5]                                                             = -1 / root1;

  return out;
} // end get_dFdq

// As get_trkCov but applied to the B cand as opposed to the track
// Trasform the x-p cov mat into x-slopes cov mat
TMatrixDSym TupleToolBs2MuMuIso2::get_BCov( const LHCb::Particle* B ) {

  double x_arr[6] = {B->endVertex()->position().x(),
                     B->endVertex()->position().y(),
                     B->endVertex()->position().z(),
                     B->momentum().x(),
                     B->momentum().y(),
                     B->momentum().z()};

  TVectorD                   x( 6, x_arr );
  const Gaudi::SymMatrix7x7& covMat = B->covMatrix();

  TMatrixDSym covM( 6 );
  covM[0][0] = covMat( 0, 0 );
  covM[0][1] = covMat( 0, 1 );
  covM[0][2] = covMat( 0, 2 );
  covM[0][3] = covMat( 0, 3 );
  covM[0][4] = covMat( 0, 4 );
  covM[0][5] = covMat( 0, 5 );
  covM[1][1] = covMat( 1, 1 );
  covM[1][2] = covMat( 1, 2 );
  covM[1][3] = covMat( 1, 3 );
  covM[1][4] = covMat( 1, 4 );
  covM[1][5] = covMat( 1, 5 );
  covM[2][2] = covMat( 2, 2 );
  covM[2][3] = covMat( 2, 3 );
  covM[2][4] = covMat( 2, 4 );
  covM[2][5] = covMat( 2, 5 );
  covM[3][3] = covMat( 3, 3 );
  covM[3][4] = covMat( 3, 4 );
  covM[3][5] = covMat( 3, 5 );
  covM[4][4] = covMat( 4, 4 );
  covM[4][5] = covMat( 4, 5 );
  covM[5][5] = covMat( 5, 5 );
  for ( int i = 0; i < 6; i++ )
    for ( int j = 0; j < i; j++ ) covM[i][j] = covM[j][i];

  TMatrixD J( 5, 6 ); // J_ij = df_i/d x_j
  for ( int i = 0; i < 5; i++ )
    for ( int j = 0; j < 6; j++ ) J[i][j] = 0;

  J[0][0] = J[1][1] = J[2][2] = 1;
  J[3][3]                     = 1 / x[5];
  J[3][5]                     = -x[3] / ( x[5] * x[5] );
  J[4][4]                     = 1 / x[5];
  J[4][5]                     = -x[4] / ( x[5] * x[5] );

  // cov(x,y,slx,sly) = J*covM*J^T
  TMatrixD    covM_temp  = J * covM;
  TMatrixD    covM_temp2 = covM_temp * J.T();
  TMatrixDSym covM_new( 5 );
  for ( int i = 0; i < 5; i++ )
    for ( int j = 0; j < 5; j++ ) covM_new[i][j] = covM_temp2[i][j];

  return covM_new;
} // end get_BCov

TMatrixD TupleToolBs2MuMuIso2::get_dFdp( TVector3 v0, TVector3 x0 ) {

  // x0 is the vector with the position of the first measurement layer
  TMatrixD out( 4, 2 );

  out[0][0] = out[1][1] = x0[2] - v0[2];
  out[2][0] = out[3][1] = 1;
  out[0][1] = out[1][0] = out[2][1] = out[3][0] = 0;

  return out;
} // end get_dFdp

TMatrixD TupleToolBs2MuMuIso2::get_dFdv( TVectorD p ) {

  TMatrixD out( 4, 3 );

  for ( int i = 0; i < 4; i++ )
    for ( int j = 0; j < 3; j++ ) out[i][j] = 0;

  out[0][0] = 1;
  out[1][1] = 1;
  out[0][2] = -p[0];
  out[1][2] = -p[1];

  return out;
} // end get_dFdv

TMatrixDSym TupleToolBs2MuMuIso2::get_trkCov( Gaudi::XYZPoint& pos_trk, Gaudi::XYZVector& mom_trk,
                                              const Gaudi::SymMatrix6x6& cov_trk,
                                              const Gaudi::SymMatrix3x3& slope_cov_trk ) {

  TVectorD x( 5 );
  x[0] = pos_trk.x();
  x[1] = pos_trk.y();
  x[2] = mom_trk.x();
  x[3] = mom_trk.y();
  x[4] = mom_trk.z();

  TMatrixDSym covM( 5 );

  covM[0][0] = cov_trk( 0, 0 );
  covM[0][1] = cov_trk( 0, 1 );
  covM[0][2] = cov_trk( 0, 2 );
  covM[0][3] = cov_trk( 0, 3 );
  covM[0][4] = cov_trk( 0, 4 );
  covM[1][1] = cov_trk( 1, 1 );
  covM[1][2] = cov_trk( 1, 2 );
  covM[1][3] = cov_trk( 1, 3 );
  covM[1][4] = cov_trk( 1, 4 );
  covM[2][2] = cov_trk( 2, 2 );
  covM[2][3] = cov_trk( 2, 3 );
  covM[2][4] = cov_trk( 2, 4 );
  covM[3][3] = cov_trk( 3, 3 );
  covM[3][4] = cov_trk( 3, 4 );
  covM[4][4] = cov_trk( 4, 4 );

  // Symmetrize the matrix
  for ( int i = 0; i < 5; i++ )
    for ( int j = 0; j < i; j++ ) covM[i][j] = covM[j][i];

  // Define the Jacobian of the {x,y,px,py,pz} -> {x,y,px/pz,py/pz} transformation:
  //{x,y,px/pz,py/pz}={f_i} {x,y,px,py,pz}={x_j}
  TMatrixD J( 4, 5 ); // J_ij = df_i/d x_j

  for ( int i = 0; i < 4; i++ )
    for ( int j = 0; j < 5; j++ ) J[i][j] = 0;

  J[0][0] = J[1][1] = 1;
  J[2][2]           = 1 / x[4];
  J[2][4]           = -x[2] / ( x[4] * x[4] );
  J[3][3]           = 1 / x[4];
  J[3][4]           = -x[3] / ( x[4] * x[4] );

  // cov(x,y,slx,sly) = J*covM*J^T
  TMatrixD    covM_temp  = J * covM;
  TMatrixD    covM_temp2 = covM_temp * J.T();
  TMatrixDSym covM_new( 4 );

  for ( int i = 0; i < 4; i++ )
    for ( int j = 0; j < 4; j++ ) covM_new[i][j] = covM_temp2[i][j];

  for ( int l = 0; l < 4; l++ )
    for ( int m = 0; m < 4; m++ ) covM_new[l][m] = 0;

  covM_new[0][0] = cov_trk( 0, 0 );
  covM_new[0][1] = cov_trk( 0, 1 );
  covM_new[1][1] = cov_trk( 1, 1 );
  covM_new[1][0] = covM_new[0][1];
  covM_new[2][2] = slope_cov_trk( 0, 0 );
  covM_new[2][3] = slope_cov_trk( 0, 1 );
  covM_new[3][3] = slope_cov_trk( 1, 1 );
  covM_new[3][2] = covM_new[2][3];

  return covM_new;
} // end get_trkCov
