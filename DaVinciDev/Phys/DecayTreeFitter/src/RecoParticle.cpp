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
#include <stdio.h>

#include "Event/Particle.h"
#include "FitParams.h"
#include "RecoParticle.h"

namespace DecayTreeFitter {
  extern int vtxverbose;

  RecoParticle::RecoParticle( const LHCb::Particle& bc, const ParticleBase* mother ) : ParticleBase( bc, mother ) {}

  RecoParticle::~RecoParticle() {}

  std::string RecoParticle::parname( int index ) const { return ParticleBase::parname( index + 4 ); }

  ErrCode RecoParticle::projectConstraint( Constraint::Type type, const FitParams& fitparams, Projection& p ) const {
    ErrCode status;
    switch ( type ) {
    case Constraint::track:
    case Constraint::photon:
      status |= projectRecoConstraint( fitparams, p );
      break;
    default:
      status |= ParticleBase::projectConstraint( type, fitparams, p );
    }
    return status;
  }

  double RecoParticle::chiSquare( const FitParams* fitparams ) const {
    // project
    Projection p( fitparams->dim(), dimM() );
    projectRecoConstraint( *fitparams, p );
    return p.chiSquare();
  }
} // namespace DecayTreeFitter
