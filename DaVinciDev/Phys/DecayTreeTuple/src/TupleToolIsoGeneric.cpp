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
// $Id: TupleToolIsoGeneric.cpp
// Include files

#include "TupleToolIsoGeneric.h"
#include "Event/Particle.h"
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IPVReFitter.h"
#include "TMath.h"
#include <Kernel/DaVinciFun.h>
#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
#include <algorithm>
#include <iostream>
#include <utility>

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolIsoGeneric
// Documentation:
// You can find extensive explanations in this talk:
// https://indico.cern.ch/event/533133/contributions/2172015/attachments/1277028/1895185/C_Isolation_BK.pdf
// TupleToolIsoGeneric              =  tuple.addTupleTool("TupleToolIsoGeneric")
// // you can run the tool on a custom container of particles:
// // TupleToolIsoGeneric.ParticlePath =  AllPions_seq.outputLocation()
// TupleToolIsoGeneric.VerboseMode = True
// // If you are running over MC sample you can add MC truth info. to the tracks
// MCTruth_noniso          = TupleToolMCTruth('MCTruth_noniso')
// MCTruth_noniso.ToolList = ["MCTupleToolHierarchy"]
// LoKiTool_noniso = LoKi__Hybrid__TupleTool("LoKiTool_noniso")
// LoKiTool_noniso.Variables     = { "ETA" : "ETA" , "PHI" : "PHI" , "TRTYPE" : "TRTYPE"}
// // You can add any TupleTool you want as long as it takes as input charged tracks
// TupleToolIsoGeneric.ToolList += ["TupleToolMCTruth/MCTruth_noniso" , "LoKi::Hybrid::TupleTool/LoKiTool_noniso"]
// TupleToolIsoGeneric.addTool(LoKiTool_noniso)
// TupleToolIsoGeneric.addTool(MCTruth_noniso)
// TupleToolIsoGeneric.OutputLevel = 6
//
// Author: Basem KHANJI
//-----------------------------------------------------------------------------
using namespace LHCb;
using namespace DaVinci;

using namespace Gaudi::Units;
using namespace ROOT::Math;
DECLARE_COMPONENT( TupleToolIsoGeneric )
//============================================================================
TupleToolIsoGeneric::TupleToolIsoGeneric( const std::string& type, const std::string& name, const IInterface* parent )
    : TupleToolBase( type, name, parent ) {
  declareInterface<IParticleTupleTool>( this );
}
// =============================================================================
StatusCode TupleToolIsoGeneric::initialize() {

  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc(), this );
  if ( !m_dva ) return Error( "Couldn't get parent DVAlgorithm" );

  m_dist = m_dva->distanceCalculator();

  m_descend = tool<IParticleDescendants>( "ParticleDescendants", this );
  if ( !m_descend ) {
    fatal() << "Unable to retrieve ParticleDescendants tool " << endmsg;
    return StatusCode::FAILURE;
  }

  m_pvReFitter = tool<IPVReFitter>( m_pvReFitterName, this );
  if ( !m_pvReFitter ) { return Error( "Unable to retrieve IPVReFitter instance" ); }

  m_vtxfitter = tool<IVertexFit>( "LoKi::VertexFitter" );
  if ( !m_vtxfitter ) {
    fatal() << "Unable to retrieve Loki VertexFitter" << endmsg;
    return StatusCode::FAILURE;
  }

  m_Reader = new TMVA::Reader( "!Silent" );
  m_Reader->AddVariable( "track_minIPchi2", &track_minIPchi2 );
  m_Reader->AddVariable( "track_eta", &track_eta );
  m_Reader->AddVariable( "track_phi", &track_phi );
  m_Reader->AddVariable( "track_pt", &track_pt );
  m_Reader->AddVariable( "track_MatchChi2", &track_Matchchi2 );
  m_Reader->AddVariable( "track_pvdis_mu", &PVdis_SignalTrack_track );
  m_Reader->AddVariable( "tracksvdis_mu", &SVdis_SignalTrack_track );
  m_Reader->AddVariable( "track_doca_mu", &doca_SignalTrack );
  m_Reader->AddVariable( "track_angle_mu", &angle_SignalTrack );
  m_Reader->AddVariable( "track_fc_mu", &fc_mu );
  m_Reader->BookMVA( "BDTG method", m_BDTLocation.value() + m_weightsName.value() );
  if ( !m_Reader ) { return Error( "Unable to retrieve the IVertexFit tool", StatusCode::FAILURE ); }

  for ( std::vector<std::string>::iterator iName = m_toolNames.begin(); iName != m_toolNames.end(); ++iName ) {
    m_tuple.push_back( tool<IParticleTupleTool>( *iName, this ) );
  }

  return sc;
}
//=============================================================================
StatusCode TupleToolIsoGeneric::fill( const LHCb::Particle* mother, const LHCb::Particle* P, const std::string& head,
                                      Tuples::Tuple& tuple ) {
  Assert( P && m_dva, "No mother or particle, or tools misconfigured." );

  const std::string prefix = fullName( head );
  bool              test   = true;
  if ( !P ) { return StatusCode( test ); }

  if ( !( P->isBasicParticle() ) ) return StatusCode( test );
  if ( mother->isBasicParticle() ) return StatusCode( test );

  int    iso_SignalTrack = 0;
  double minMinBDT_New   = 1.5;

  double Best_PAIR_M_fromiso = -1;
  double Best_PAIR_VCHI2     = -1;
  double Best_PAIR_D         = -1;
  double Best_PAIR_DCHI2     = -1;
  double Best_Vtx_X          = -999;
  double Best_Vtx_Y          = -999;
  double Best_Vtx_Z          = -999;

  double Best_track_minIPchi2         = -999;
  double Best_PVdis_SignalTrack_track = -999;
  double Best_SVdis_SignalTrack_track = -999;
  double Best_doca_SignalTrack        = -999;
  double Best_angle_SignalTrack       = -999;
  double Best_fc_mu                   = -999;

  XYZTVector        NonIsolatedTrack_4Vec;
  const VertexBase* aPV = NULL;
  aPV                   = m_dva->bestVertex( mother ); // bestVertex
  Gaudi::XYZPoint PosPV( 0, 0, 200000 );
  Gaudi::XYZPoint PosSV( 0, 0, 200000 );
  if ( !aPV ) { return StatusCode::SUCCESS; }
  if ( aPV ) { PosPV = aPV->position(); }

  const VertexBase* endv = mother->endVertex();
  if ( !endv ) { return StatusCode::FAILURE; }

  PosSV                 = mother->endVertex()->position();
  double Bd_eta         = mother->momentum().Eta();
  double Bd_phi         = mother->momentum().Phi();
  double Bd_pT          = mother->pt();
  double Sum_of_trackpt = 0;

  XYZTVector SignalParticle4vect = P->momentum();
  double     PAIR_VtxChi2        = -1;
  double     PAIR_Vtx_X          = -999;
  double     PAIR_Vtx_Y          = -999;
  double     PAIR_Vtx_Z          = -999;

  Particle::ConstVector parts_cand = m_descend->descendants( mother );
  if ( parts_cand.empty() ) return StatusCode( test );
  ;
  LHCb::Particle::Range allparts;
  if ( exist<LHCb::Particle::Range>( m_ParticlePath ) ) {
    allparts = get<LHCb::Particle::Range>( m_ParticlePath );
  } else {
    return Error( "Nothing found at " + m_ParticlePath, StatusCode::SUCCESS, 1 );
  }

  const Particle* p_bestNonisotr = NULL;

  std::vector<double> V_track_minIPchi2, V_track_eta, V_track_phi, V_track_pt, V_track_Matchchi2,
      V_PVdis_SignalTrack_track, V_SVdis_SignalTrack_track, V_doca_SignalTrack, V_angle_SignalTrack, V_fc_mu,
      V_PAIR_VtxChi2, V_PAIR_Vtx_X, V_PAIR_Vtx_Y, V_PAIR_Vtx_Z, V_PAIR_D, V_PAIR_DCHI2, V_MinIsoBDT, V_trackxfirst,
      V_trackyfirst, V_trackzfirst;

  LHCb::Particle::Range::const_iterator im;
  for ( im = allparts.begin(); im != allparts.end(); ++im ) {
    debug() << "I am inside the allparts loop now ..." << endmsg;
    const Particle* axp = ( *im );
    debug() << "after defining the particle .." << endmsg;
    const LHCb::ProtoParticle* proto = axp->proto();
    if ( !proto ) continue;
    const LHCb::Track* track = proto->track();
    if ( !track ) continue;
    bool                                  isInDecay = false;
    Particle::ConstVector::const_iterator it_p;
    for ( it_p = parts_cand.begin(); it_p != parts_cand.end(); it_p++ ) {
      if ( !( *it_p )->isBasicParticle() ) continue;
      const LHCb::Track* Track_decay = ( *it_p )->proto()->track();
      // if(!track || !Track_decay) continue;
      if ( Track_decay && track ) {
        debug() << "track of decay particle exist!" << endmsg;
        if ( Track_decay == track ) {
          debug() << " track PID is  ================>>> " << endmsg;
          isInDecay = true;
        }
      }
    }

    if ( isInDecay ) continue;

    track_minIPchi2 = get_MINIPCHI2( axp );
    track_pt        = axp->pt();

    track_eta             = track->momentum().eta();
    track_phi             = track->momentum().phi();
    track_Matchchi2       = track->info( LHCb::Track::FitMatchChi2, -1 );
    XYZTVector track4vect = axp->momentum();

    XYZTVector PAIR_mom  = SignalParticle4vect + track4vect;
    double     PAIR_mass = std::abs( std::sqrt( PAIR_mom.M2() ) );

    Gaudi::XYZPoint  pos_track( track->position() );
    Gaudi::XYZVector mom_track( track->momentum() );
    Gaudi::XYZPoint  pos_SignalTrack( P->proto()->track()->position() );
    Gaudi::XYZVector mom_SignalTrack( P->proto()->track()->momentum() );

    // Calculate the input of ISO variable :
    Gaudi::XYZPoint vtx_SignalTrack( 0., 0., 0. );

    double doca_SignalTrack_d( -1. ), angle_SignalTrack_d( -1. );

    InCone( pos_SignalTrack, mom_SignalTrack, pos_track, mom_track, vtx_SignalTrack, doca_SignalTrack_d,
            angle_SignalTrack_d ); //
    doca_SignalTrack  = (Float_t)doca_SignalTrack_d;
    angle_SignalTrack = (Float_t)angle_SignalTrack_d;

    PVdis_SignalTrack_track =
        ( vtx_SignalTrack.z() - PosPV.z() ) / fabs( vtx_SignalTrack.z() - PosPV.z() ) * ( vtx_SignalTrack - PosPV ).R();
    SVdis_SignalTrack_track =
        ( vtx_SignalTrack.z() - PosSV.z() ) / fabs( vtx_SignalTrack.z() - PosSV.z() ) * ( vtx_SignalTrack - PosSV ).R();
    fc_mu = pointer( vtx_SignalTrack - PosPV, mom_track, mom_SignalTrack );
    // calculate the track-signaltrack vertex quantities:
    LHCb::Vertex vtx_PAIR;
    double       PAIR_distance       = 0;
    double       PAIR_chi2           = 0;
    StatusCode   ok_TrackSignalTrack = m_dist->distance( P, axp, PAIR_distance, PAIR_chi2 );
    StatusCode   PAIR_Vtx_fit        = m_vtxfitter->fit( vtx_PAIR, *P, *axp );

    double trackxfirst = track->position().x();
    double trackyfirst = track->position().y();
    double trackzfirst = track->position().z();
    V_trackxfirst.push_back( trackxfirst );
    V_trackyfirst.push_back( trackyfirst );
    V_trackzfirst.push_back( trackzfirst );

    V_track_minIPchi2.push_back( track_minIPchi2 );
    V_track_eta.push_back( track_eta );
    V_track_phi.push_back( track_phi );
    V_track_pt.push_back( track_pt );
    V_track_Matchchi2.push_back( track_Matchchi2 );
    V_PVdis_SignalTrack_track.push_back( PVdis_SignalTrack_track );
    V_SVdis_SignalTrack_track.push_back( SVdis_SignalTrack_track );
    V_doca_SignalTrack.push_back( doca_SignalTrack );
    V_angle_SignalTrack.push_back( angle_SignalTrack );
    V_fc_mu.push_back( fc_mu );
    double BDT_val_SignalTrack_New = 0;
    BDT_val_SignalTrack_New        = m_Reader->EvaluateMVA( "BDTG method" );

    if ( !PAIR_Vtx_fit ) { PAIR_VtxChi2 = -1; }
    PAIR_VtxChi2 = vtx_PAIR.chi2() / (double)vtx_PAIR.nDoF();
    PAIR_Vtx_X   = vtx_PAIR.position().x();
    PAIR_Vtx_Y   = vtx_PAIR.position().y();
    PAIR_Vtx_Z   = vtx_PAIR.position().z();
    V_PAIR_VtxChi2.push_back( PAIR_VtxChi2 );
    V_PAIR_Vtx_X.push_back( PAIR_Vtx_X );
    V_PAIR_Vtx_Y.push_back( PAIR_Vtx_Y );
    V_PAIR_Vtx_Z.push_back( PAIR_Vtx_Z );
    V_PAIR_D.push_back( PAIR_distance );
    V_PAIR_DCHI2.push_back( PAIR_chi2 );
    V_MinIsoBDT.push_back( BDT_val_SignalTrack_New );
    if ( minMinBDT_New > BDT_val_SignalTrack_New ) {
      Best_PAIR_VCHI2              = PAIR_VtxChi2;
      Best_Vtx_X                   = PAIR_Vtx_X;
      Best_Vtx_Y                   = PAIR_Vtx_Y;
      Best_Vtx_Z                   = PAIR_Vtx_Z;
      Best_PAIR_D                  = PAIR_distance;
      Best_PAIR_DCHI2              = PAIR_chi2;
      minMinBDT_New                = BDT_val_SignalTrack_New;
      Best_PAIR_M_fromiso          = PAIR_mass;
      p_bestNonisotr               = axp;
      Best_track_minIPchi2         = track_minIPchi2;
      Best_PVdis_SignalTrack_track = PVdis_SignalTrack_track;
      Best_SVdis_SignalTrack_track = SVdis_SignalTrack_track;
      Best_doca_SignalTrack        = doca_SignalTrack;
      Best_angle_SignalTrack       = angle_SignalTrack;
      Best_fc_mu                   = fc_mu;
    }
    // Apply giampie's iso cut:
    if ( angle_SignalTrack < 0.27 && fc_mu < 0.60 && doca_SignalTrack < 0.13 && track_minIPchi2 > 3.0 &&
         SVdis_SignalTrack_track > -0.15 && SVdis_SignalTrack_track < 30. && PVdis_SignalTrack_track > 0.5 &&
         PVdis_SignalTrack_track < 40. ) {
      iso_SignalTrack += 1;
    }

    double cone_radius = std::sqrt( std::pow( ( Bd_eta - track_eta ), 2 ) + std::pow( ( Bd_phi - track_phi ), 2 ) );
    if ( cone_radius < 1 ) { Sum_of_trackpt += axp->pt(); }
  }

  double Bd_CDF3 = Bd_pT / ( Bd_pT + Sum_of_trackpt );

  test &= tuple->column( prefix + "_Nb_NIsoTr", iso_SignalTrack );
  test &= tuple->column( prefix + "_IsoMinBDT", minMinBDT_New );
  test &= tuple->column( prefix + "_ConeIso", Bd_CDF3 );
  test &= tuple->column( prefix + "_PAIR_M", Best_PAIR_M_fromiso );
  test &= tuple->column( prefix + "_NIso_PAIR_VtxChi2", Best_PAIR_VCHI2 );
  test &= tuple->column( prefix + "_NIso_DCHI2", Best_PAIR_DCHI2 );
  test &= tuple->column( prefix + "_NIso_DIS", Best_PAIR_D );
  test &= tuple->column( prefix + "_NIsoTr_VTXX", Best_Vtx_X );
  test &= tuple->column( prefix + "_NIsoTr_VTXY", Best_Vtx_Y );
  test &= tuple->column( prefix + "_NIsoTr_VTXZ", Best_Vtx_Z );
  test &= tuple->column( prefix + "_PAIR_IPchi2", Best_track_minIPchi2 );
  test &= tuple->column( prefix + "_NIsoTr_PVDIS", Best_PVdis_SignalTrack_track );
  test &= tuple->column( prefix + "_NIsoTr_SVDIS", Best_SVdis_SignalTrack_track );
  test &= tuple->column( prefix + "_NIsoTr_DOCA", Best_doca_SignalTrack );
  test &= tuple->column( prefix + "_NIsoTr_ANGLE", Best_angle_SignalTrack );
  test &= tuple->column( prefix + "_NIsoTr_FC", Best_fc_mu );

  if ( m_verbose ) {
    test &= tuple->farray( prefix + "_Tr_MinIsoBDT", V_MinIsoBDT, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_PAIR_VTXX", V_PAIR_Vtx_X, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_PAIR_VTXY", V_PAIR_Vtx_Y, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_PAIR_VTXZ", V_PAIR_Vtx_Z, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_PAIR_VTCHI2NDOF", V_PAIR_VtxChi2, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_FirstHitX", V_trackxfirst, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_FirstHitY", V_trackyfirst, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_FirstHitZ", V_trackzfirst, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_DCHI2", V_PAIR_DCHI2, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_DIS", V_PAIR_D, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_MATCHCHI2", V_track_Matchchi2, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_PVDIS", V_PVdis_SignalTrack_track, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_SVDIS", V_SVdis_SignalTrack_track, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_DOCA", V_doca_SignalTrack, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_ANGLE", V_angle_SignalTrack, "NTr", 2000 );
    test &= tuple->farray( prefix + "_Tr_FC_MU", V_fc_mu, "NTr", 2000 );
  }

  if ( p_bestNonisotr == nullptr || !( p_bestNonisotr->proto() ) || !( p_bestNonisotr->proto()->track() ) )
    return StatusCode( test );
  else {
    for ( std::vector<IParticleTupleTool*>::const_iterator iTool = m_tuple.begin(); iTool != m_tuple.end(); ++iTool ) {
      test &= ( *iTool )->fill( NULL, p_bestNonisotr, prefix + "_NIsoTr", tuple );
    }
  }

  return StatusCode( test );
}
//=============================================================================
double TupleToolIsoGeneric::get_MINIPCHI2( const Particle* p ) {
  double minchi2 = -1;

  const RecVertex::Range PV = m_dva->primaryVertices();
  if ( !PV.empty() ) {
    for ( RecVertex::Range::const_iterator pv = PV.begin(); pv != PV.end(); ++pv ) {
      double ip, chi2;

      RecVertex  newPV( **pv );
      StatusCode scfit = m_pvReFitter->remove( p, &newPV );
      if ( !scfit ) {
        err() << "ReFitter fails!" << endmsg;
        continue;
      }

      LHCb::VertexBase* newPVPtr = (LHCb::VertexBase*)&newPV;
      const auto        test2    = m_dist->distance( p, newPVPtr, ip, chi2 );
      if ( test2 ) {
        if ( ( chi2 < minchi2 ) || ( minchi2 < 0. ) ) minchi2 = chi2;
      }
    }
  }
  return minchi2;
}
//=============================================================================
void TupleToolIsoGeneric::InCone( Gaudi::XYZPoint o1, Gaudi::XYZVector p1, Gaudi::XYZPoint o2, Gaudi::XYZVector p2,
                                  Gaudi::XYZPoint& vtx, double& doca, double& angle ) const {

  Gaudi::XYZPoint rv;
  Gaudi::XYZPoint close;
  Gaudi::XYZPoint close_mu;
  bool            fail( false );
  closest_point( o1, p1, o2, p2, close, close_mu, vtx, fail );
  if ( fail ) {
    doca  = -1.;
    angle = -1.;
  } else {
    doca  = ( close - close_mu ).R();
    angle = arcosine( p1, p2 );
  }
}

//============================================================================
double TupleToolIsoGeneric::pointer( Gaudi::XYZVector vertex, Gaudi::XYZVector p, Gaudi::XYZVector p_mu ) const {
  double           pt = p.Rho() + p_mu.Rho();
  Gaudi::XYZVector ptot( p + p_mu );
  double           num = ptot.R() * sin( arcosine( vertex, ptot ) );
  double           den = ( num + pt );
  double           fc  = num / den;
  return fc;
}
//============================================================================
void TupleToolIsoGeneric::closest_point( Gaudi::XYZPoint o, Gaudi::XYZVector p, Gaudi::XYZPoint o_mu,
                                         Gaudi::XYZVector p_mu, Gaudi::XYZPoint& close1, Gaudi::XYZPoint& close2,
                                         Gaudi::XYZPoint& vertex, bool& fail ) const {
  Gaudi::XYZVector v0( o - o_mu );
  Gaudi::XYZVector v1( p.unit() );
  Gaudi::XYZVector v2( p_mu.unit() );
  Gaudi::XYZPoint  temp1( 0., 0., 0. );
  Gaudi::XYZPoint  temp2( 0., 0., 0. );
  fail = false;

  double d02   = v0.Dot( v2 );
  double d21   = v2.Dot( v1 );
  double d01   = v0.Dot( v1 );
  double d22   = v2.Dot( v2 );
  double d11   = v1.Dot( v1 );
  double denom = d11 * d22 - d21 * d21;
  if ( fabs( denom ) <= 0. ) {
    close1 = temp1;
    close2 = temp2;
    fail   = true;
  } else {
    double numer = d02 * d21 - d01 * d22;
    double mu1   = numer / denom;
    double mu2   = ( d02 + d21 * mu1 ) / d22;
    close1       = o + v1 * mu1;
    close2       = o_mu + v2 * mu2;
  }
  vertex = ( close1 + ( close2 - close1 ) * 0.5 );
}
//============================================================================
double TupleToolIsoGeneric::arcosine( Gaudi::XYZVector p1, Gaudi::XYZVector p2 ) const {

  double num    = ( p1.Cross( p2 ) ).R();
  double den    = p1.R() * p2.R();
  double seno   = num / den;
  double coseno = p1.Dot( p2 ) / den;
  double alpha  = asin( fabs( seno ) );
  if ( coseno < 0 ) { alpha = ROOT::Math::Pi() - alpha; }
  return alpha;
}
