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
// Include files

// local
#include "TupleToolKinematic.h"

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/Particle.h"

#include "Kernel/IParticleTransporter.h"

//-----------------------------------------------------------------------------
// Implementation file for class : EventInfoTupleTool
//
// 2007-11-07 : Jeremie Borel
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_COMPONENT( TupleToolKinematic )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolKinematic::TupleToolKinematic( const std::string& type, const std::string& name, const IInterface* parent )
    : TupleToolBase( type, name, parent ), m_transporter(), m_transporterName( "ParticleTransporter:PUBLIC" ) {
  declareInterface<IParticleTupleTool>( this );

  declareProperty( "Transporter", m_transporterName );
}

//=============================================================================

StatusCode TupleToolKinematic::initialize() {
  const StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;

  m_transporter = tool<IParticleTransporter>( m_transporterName, this );

  return sc;
}

//=============================================================================

StatusCode TupleToolKinematic::fill( const LHCb::Particle* mother, const LHCb::Particle* P, const std::string& head,
                                     Tuples::Tuple& tuple ) {
  const std::string prefix = fullName( head );

  bool test = true;
  if ( P ) {
    test &= tuple->column( prefix + "_P", P->p() );
    test &= tuple->column( prefix + "_PT", P->pt() );

    // momentum components
    test &= tuple->column( prefix + "_P", P->momentum() );

    // reference point:
    if ( isVerbose() ) test &= tuple->column( prefix + "_REFP", P->referencePoint() );
    // mass before fit (what CombinationCut cuts on)
    // if(isVerbose() && !(P->isBasicParticle()) ) test &= tuple->column( prefix+"_PreFitMass", preFitMass(P) );

    if ( !P->isBasicParticle() || P->particleID().pid() == 111 )
      test &= tuple->column( prefix + "_MM", P->measuredMass() );
    if ( !P->isBasicParticle() || P->particleID().pid() == 111 )
      test &= tuple->column( prefix + "_MMERR", P->measuredMassErr() );
    test &= tuple->column( prefix + "_M", P->momentum().M() );

    // Compute momentum at origin vertex
    if ( isVerbose() && mother && P->isBasicParticle() && P->charge() != 0 ) {
      const LHCb::Vertex* originvtx = originVertex( mother, P );
      if ( originvtx ) {
        double zvtx = originvtx->position().Z();

        LHCb::Particle transParticle;
        if ( !m_transporter )
          Error( "null pointer m_transporter !!!!" ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        test &= m_transporter->transport( P, zvtx, transParticle );

        test &= tuple->column( prefix + "_AtVtx_P", transParticle.momentum() );
        test &= tuple->column( prefix + "_AtVtx_P", transParticle.p() );
        test &= tuple->column( prefix + "_AtVtx_PT", transParticle.pt() );
      }

      // Compute momentum at Secondary Vertex
      const LHCb::Vertex* SV = mother->endVertex();
      if ( SV ) {
        double zSV = SV->position().Z();

        LHCb::Particle transParticle_SV;
        if ( !m_transporter )
          Error( "null pointer m_transporter !!!!" ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        test &= m_transporter->transport( P, zSV, transParticle_SV );

        test &= tuple->column( prefix + "_AtSV_P", transParticle_SV.momentum() );
        test &= tuple->column( prefix + "_AtSV_P", transParticle_SV.p() );
        test &= tuple->column( prefix + "_AtSV_PT", transParticle_SV.pt() );
      }
    }
  } else {
    return StatusCode::FAILURE;
  }
  return StatusCode( test );
}

//=============================================================================

double TupleToolKinematic::preFitMass( const LHCb::Particle* p ) const {
  Gaudi::LorentzVector Mom;
  for ( SmartRefVector<LHCb::Particle>::const_iterator d = p->daughters().begin(); d != p->daughters().end(); ++d ) {
    Mom += ( *d )->momentum();
  }
  return Mom.M();
}

// =====================================================
// find origin vertex in the decay chain
// =====================================================
const LHCb::Vertex* TupleToolKinematic::originVertex( const LHCb::Particle* top, const LHCb::Particle* P ) const {
  // this used to pass back zero if P was a basic particle.
  if ( top == P || top->isBasicParticle() ) return nullptr;

  const auto& dau = top->daughters();
  if ( dau.empty() ) return nullptr;

  for ( const auto& d : dau ) {
    if ( P == d ) { // I found the daughter
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << "It's a daughter, returning mother's endvertex : " << endmsg;
      return top->endVertex();
    }
  }

  // vertex not yet found, get deeper in the decay:
  for ( const auto& d : dau ) {
    const auto vv = originVertex( d, P );
    if ( !vv ) continue;
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Went up : " << vv << endmsg;
    return vv;
  }
  return nullptr;
}
