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

// Include files

// local
#include "TupleToolVtxDaughters.h"

#include "DecayTreeFitter/Fitter.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IPVReFitter.h"
#include "Kernel/IParticleCombiner.h"
#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
#include <Kernel/IVertexFit.h>

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrackInterfaces/ITrackStateProvider.h"

#include "Event/Particle.h"

using namespace LHCb;
using namespace DecayTreeFitter;

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolVtxDaughters
//
// 2014-10-08 : Adlene Hicheur, first hand
// 2019-12: Adlene Hicheur, upgrade to Triplet vertices, time to commit
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_COMPONENT( TupleToolVtxDaughters )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolVtxDaughters::TupleToolVtxDaughters( const std::string& type, const std::string& name,
                                              const IInterface* parent )
    : TupleToolBase( type, name, parent ), m_dist( 0 ), m_vtx( 0 ), m_dva( 0 ) {
  declareInterface<IParticleTupleTool>( this );
  declareProperty( "FillSubVtxInfo", m_fillSubVtx = false, "Fill complete pair vertex info" );
  declareProperty( "DoThreeProng", m_doThreeProng = false, "Compute three prong vertices" );
  // replaced by Verbose
}

//=============================================================================

StatusCode TupleToolVtxDaughters::initialize() {
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  m_dva = Gaudi::Utils::getIDVAlgorithm( contextSvc(), this );
  if ( !m_dva ) return Error( "Couldn't get parent DVAlgorithm" );

  m_dist = m_dva->distanceCalculator();
  if ( !m_dist ) { return Error( "Unable to retrieve the IDistanceCalculator tool" ); }

  m_vtx = m_dva->vertexFitter();
  if ( !m_vtx ) { return Error( "Unable to retrieve the IVertexFit tool" ); }

  return sc;
}

//=============================================================================
StatusCode TupleToolVtxDaughters::fill( const Particle* mother, const Particle* P, const std::string& head,
                                        Tuples::Tuple& tuple ) {
  const std::string prefix = fullName( head );

  StatusCode sc = StatusCode::SUCCESS;

  Assert( P && m_dist && m_dva && m_vtx, "No mother or particle, or tools misconfigured." );

  if ( isPureNeutralCalo( P ) )
    return Warning( "Will not fill geometry tuple for neutral Calo particles. No worry", StatusCode::SUCCESS, 10 );
  if ( msgLevel( MSG::VERBOSE ) )
    verbose() << "TupleToolVtxDaughters::fill " << mother << " " << P << " " << prefix << endmsg;

  //=========================================================================
  // fill end vertex info
  //=========================================================================
  if ( P->isBasicParticle() || P->daughters().size() < 3 ) {
    if ( msgLevel( MSG::DEBUG ) )
      debug()
          << "TupleToolVtxDaughters does not fill info for basic particles or decays with less than three daughters."
          << endmsg;
  } else {
    //=========================================================================
    if ( msgLevel( MSG::VERBOSE ) && P->endVertex() ) { // https://its.cern.ch/jira/browse/LHCBPS-598
      verbose() << "Before cast : " << P->endVertex() << endmsg;
      verbose() << "Container " << P->endVertex()->parent()->registry()->identifier() << " key "
                << P->endVertex()->key() << endmsg;
    }

    const VertexBase* evtx = P->endVertex();
    if ( !evtx ) {
      fatal() << "Can't retrieve the end vertex for " << prefix
              << ". TupleToolVtxDaughters cannot act on such a particle. Define a branch" << endmsg;
      fatal() << "See https://twiki.cern.ch/twiki/bin/view/LHCb/DaVinciTutorial7#Branches" << endmsg;
      return StatusCode::FAILURE;
    }

    if ( msgLevel( MSG::VERBOSE ) ) { // https://its.cern.ch/jira/browse/LHCBPS-598
      verbose() << "End Vertex : " << *evtx << endmsg;
      verbose() << "Container " << evtx->parent()->registry()->identifier() << " key " << evtx->key() << endmsg;
    }
    // Pairing daughters for subvertices:

    const SmartRefVector<LHCb::Particle>& dau = P->daughters();

    int sizeDau = P->daughters().size();

    // LHCb::DecayTree tree ( *P ) ;

    // Two prongs

    for ( int i = 0; i < sizeDau; ++i ) {
      for ( int j = i + 1; j < sizeDau; ++j ) {
        char buffer[40];
        sprintf( buffer, "_%d%d", i + 1, j + 1 );
        std::string                 trail( buffer );
        const VertexBase*           vtxpair = NULL;
        double                      myMass  = -999.;
        LHCb::Particle              tempmother( ParticleID( 531 ) );
        LHCb::Particle::ConstVector tmp( 2 );
        tmp[0] = dau[i];
        tmp[1] = dau[j];

        vtxpair = ProngVertex( tmp, myMass, tempmother );

        LHCb::Particle::ConstVector newTree;
        newTree.push_back( &tempmother );
        LHCb::Particle newHead( ParticleID( 541 ) );
        newHead.addToDaughters( &tempmother );
        LHCb::Vertex tmpvertex;

        for ( int k = 0; k < sizeDau; k++ ) {
          if ( k != i && k != j ) {
            newTree.push_back( dau[k] );
            newHead.addToDaughters( dau[k] );
          }
        }
        // create the fitter object
        ToolHandle<ITrackStateProvider> m_stateprovider( "TrackStateProvider" );
        StatusCode                      sc            = m_stateprovider.retrieve();
        const ITrackStateProvider*      stateprovider = ( m_stateprovider.empty() ? NULL : &( *m_stateprovider ) );
        DecayTreeFitter::Fitter         myfitter( newHead, stateprovider );
        // fit and get the fit parameters
        myfitter.fit();
        const LHCb::Particle*              pairmother   = &tempmother;
        const Gaudi::Math::ParticleParams* myPairParams = myfitter.fitParams( pairmother );
        // get the decay length with err
        const Gaudi::Math::ValueWithError& declen = myPairParams->decayLength();
        double                             flight = declen.value();
        double                             errlen = declen.error();
        double                             fdchi2 = -999;
        if ( errlen > 1e-8 ) fdchi2 = ( flight / errlen ) * ( flight / errlen );
        const Gaudi::LorentzVector&          p4           = myPairParams->momentum();
        double                               measPairMass = p4.M();
        const Gaudi::Math::Point3DWithError& pairposerr   = myPairParams->position();
        const Gaudi::XYZPoint&               pairpos      = pairposerr.point();
        const Gaudi::SymMatrix3x3&           paircov      = pairposerr.covariance();

        if ( vtxpair ) {

          bool   test = true;
          double dist = 0;
          double chi2 = 0;
          sc          = m_dist->distance( P->endVertex(), vtxpair, dist, chi2 );
          if ( !sc.isFailure() ) {

            test &= tuple->column( prefix + "_SubVtx" + trail + "_FD", flight );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_FDCHI2", fdchi2 );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_Mass", measPairMass );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_DIRA", dira( m_dva->bestVertex( P ), tempmother ) );
          }

          if ( !test ) Warning( "Error in fillFlight " + prefix, StatusCode::SUCCESS, 1 ).ignore();

          if ( m_fillSubVtx ) {

            if ( sc.isFailure() ) { return Warning( "Could not fill Endvertex " + prefix, StatusCode::SUCCESS, 1 ); }
            test = true;
            test &= tuple->column( prefix + "_SubVtx" + trail + "_", pairpos );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_XERR", std::sqrt( paircov( 0, 0 ) ) );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_YERR", std::sqrt( paircov( 1, 1 ) ) );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_ZERR", std::sqrt( paircov( 2, 2 ) ) );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_CHI2", vtxpair->chi2() );
            test &= tuple->column( prefix + "_SubVtx" + trail + "_NDOF", vtxpair->nDoF() );
            test &= tuple->matrix( prefix + "_SubVtx" + trail + "_COV_", paircov );

            if ( !test ) Warning( "Error in filling Vtx info " + prefix + trail, StatusCode::SUCCESS, 1 ).ignore();
          }
        }
      }
    }

    // Three prongs:

    if ( m_doThreeProng && sizeDau > 3 ) {

      for ( int i = 0; i < sizeDau; ++i ) {
        for ( int j = i + 1; j < sizeDau; ++j ) {
          for ( int k = j + 1; k < sizeDau; ++k ) {

            char buffer[40];
            sprintf( buffer, "_%d%d%d", i + 1, j + 1, k + 1 );
            std::string                 trail( buffer );
            const VertexBase*           vtx3prong = NULL;
            double                      myMass    = -999.;
            LHCb::Particle              tempmother( ParticleID( 531 ) );
            LHCb::Particle::ConstVector tmp( 3 );
            tmp[0]    = dau[i];
            tmp[1]    = dau[j];
            tmp[2]    = dau[k];
            vtx3prong = ProngVertex( tmp, myMass, tempmother );
            // create the new head:
            LHCb::Particle::ConstVector newTree;
            newTree.push_back( &tempmother );
            LHCb::Particle newHead( ParticleID( 541 ) );
            newHead.addToDaughters( &tempmother );
            LHCb::Vertex tmpvertex;

            for ( int ibach = 0; ibach < sizeDau; ibach++ ) {
              if ( ibach != i && ibach != j && ibach != k ) {
                newTree.push_back( dau[ibach] );
                newHead.addToDaughters( dau[ibach] );
              }
            }
            // create the fitter object
            ToolHandle<ITrackStateProvider> m_stateprovider( "TrackStateProvider" );
            StatusCode                      sc            = m_stateprovider.retrieve();
            const ITrackStateProvider*      stateprovider = ( m_stateprovider.empty() ? NULL : &( *m_stateprovider ) );
            DecayTreeFitter::Fitter         myfitter( newHead, stateprovider );
            // fit and get the results
            myfitter.fit();
            const LHCb::Particle*              pairmother   = &tempmother;
            const Gaudi::Math::ParticleParams* myPairParams = myfitter.fitParams( pairmother );
            // get the decay length with err
            const Gaudi::Math::ValueWithError& declen = myPairParams->decayLength();
            double                             flight = declen.value();
            double                             errlen = declen.error();
            double                             fdchi2 = -999;
            if ( errlen > 1e-8 ) fdchi2 = ( flight / errlen ) * ( flight / errlen );
            const Gaudi::LorentzVector&          p4           = myPairParams->momentum();
            double                               measPairMass = p4.M();
            const Gaudi::Math::Point3DWithError& pairposerr   = myPairParams->position();
            const Gaudi::XYZPoint&               pairpos      = pairposerr.point();
            const Gaudi::SymMatrix3x3&           paircov      = pairposerr.covariance();
            if ( vtx3prong ) {
              bool   test = true;
              double dist = 0;
              double chi2 = 0;
              sc          = m_dist->distance( P->endVertex(), vtx3prong, dist, chi2 );
              if ( !sc.isFailure() ) {
                test &= tuple->column( prefix + "_SubVtx" + trail + "_FD", flight );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_FDCHI2", fdchi2 );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_Mass", measPairMass );
              }

              if ( !test ) Warning( "Error in fillFlight " + prefix, StatusCode::SUCCESS, 1 ).ignore();

              if ( m_fillSubVtx ) {

                if ( sc.isFailure() ) {
                  return Warning( "Could not fill Endvertex " + prefix, StatusCode::SUCCESS, 1 );
                }
                test = true;
                test &= tuple->column( prefix + "_SubVtx" + trail + "_", pairpos );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_XERR", std::sqrt( paircov( 0, 0 ) ) );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_YERR", std::sqrt( paircov( 1, 1 ) ) );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_ZERR", std::sqrt( paircov( 2, 2 ) ) );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_CHI2", vtx3prong->chi2() );
                test &= tuple->column( prefix + "_SubVtx" + trail + "_NDOF", vtx3prong->nDoF() );
                test &= tuple->matrix( prefix + "_SubVtx" + trail + "_COV_", paircov );

                if ( !test ) Warning( "Error in filling Vtx info " + prefix + trail, StatusCode::SUCCESS, 1 ).ignore();
              }
            }
          }
        }
      }

    } // end doThreeProng
  }

  //=========================================================================

  return sc;
}

//=========================================================================
// fill vertex stuff
//=========================================================================
StatusCode TupleToolVtxDaughters::fillVertex( const LHCb::VertexBase* vtx, const std::string& vtx_name,
                                              Tuples::Tuple& tuple ) const {
  bool test = true;

  // decay vertex information:
  if ( !vtx ) {
    Gaudi::XYZPoint pt( -999., -999., -999. ); // arbitrary point
    test &= tuple->column( vtx_name + "_", pt );
    test &= tuple->column( vtx_name + "_XERR", -999. );
    test &= tuple->column( vtx_name + "_YERR", -999. );
    test &= tuple->column( vtx_name + "_ZERR", -999. );
    test &= tuple->column( vtx_name + "_CHI2", -999. );
    test &= tuple->column( vtx_name + "_NDOF", -1 );
    test &= tuple->matrix( vtx_name + "_COV_", Gaudi::SymMatrix3x3() );
  } else {
    //    std::cout<<"position of passed VTX - fillVertex ++++ "<<vtx->position()<<std::endl;
    test &= tuple->column( vtx_name + "_", vtx->position() );
    const Gaudi::SymMatrix3x3& m = vtx->covMatrix();
    test &= tuple->column( vtx_name + "_XERR", std::sqrt( m( 0, 0 ) ) );
    test &= tuple->column( vtx_name + "_YERR", std::sqrt( m( 1, 1 ) ) );
    test &= tuple->column( vtx_name + "_ZERR", std::sqrt( m( 2, 2 ) ) );
    test &= tuple->column( vtx_name + "_CHI2", vtx->chi2() );
    test &= tuple->column( vtx_name + "_NDOF", vtx->nDoF() );
    test &= tuple->matrix( vtx_name + "_COV_", m );
  }

  // --------------------------------------------------
  if ( !test ) Warning( "Error in fillVertex " + vtx_name, StatusCode::SUCCESS, 1 ).ignore();
  return StatusCode( test );
}
//=========================================================================
// fill flight distance, angle...
//=========================================================================
StatusCode TupleToolVtxDaughters::fillFlight( const LHCb::VertexBase* PairVtx, const Particle* P,
                                              const std::string& prefix, Tuples::Tuple& tuple,
                                              const std::string& trail ) const {
  bool test = true;
  // --------------------------------------------------
  if ( !PairVtx ) {
    test &= tuple->column( prefix + "_FD" + trail, -999. );
    test &= tuple->column( prefix + "_FDCHI2" + trail, -999. );
    test &= tuple->column( prefix + "_DIRA" + trail, -999. );
  } else {

    // flight distance
    double     dist = 0;
    double     chi2 = 0;
    StatusCode sc   = m_dist->distance( P->endVertex(), PairVtx, dist, chi2 );
    if ( sc.isFailure() ) return sc;

    test &= tuple->column( prefix + "_FD" + trail, dist );
    test &= tuple->column( prefix + "_FDCHI2" + trail, chi2 );
    // --------------------------------------------------
    // cosine of (flight distance) dot (momentum):
    // find the origin vertex. Either the primary or the origin in the
    // decay
    //   test &= tuple->column( prefix + "_DIRA"+trail, dira(oriVtx,P) );
  }

  if ( !test ) Warning( "Error in fillFlight " + prefix, StatusCode::SUCCESS, 1 ).ignore();
  return StatusCode( test );
}
// =====================================================
// find pair vertex in the decay chain
// =====================================================

const VertexBase* TupleToolVtxDaughters::ProngVertex( LHCb::Particle::ConstVector daus, double& Mass,
                                                      LHCb::Particle& PairMom ) const {

  LHCb::Vertex   DauVtx;
  LHCb::Particle DummyPart( ParticleID( 531 ) );

  StatusCode scFit = m_vtx->fit( daus, DauVtx, DummyPart );
  // StatusCode scFit = m_vtx->fit(DauPairVtx,tmp);
  if ( !scFit ) {
    Warning( "Pair Vtx Fit error" ).ignore();
    return 0;
  }

  const VertexBase* myVtx = DummyPart.endVertex();
  Mass                    = DummyPart.measuredMass();
  PairMom                 = DummyPart;
  return myVtx;
}
