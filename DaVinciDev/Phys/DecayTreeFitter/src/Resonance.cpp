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
#include <algorithm>
#include <iomanip>

#include "Resonance.h"

namespace DecayTreeFitter {

  extern int vtxverbose;

  Resonance::Resonance( const LHCb::Particle& bc, const ParticleBase* mother, const Configuration& config )
      : InternalParticle( bc, mother, config ) {}

  Resonance::~Resonance() {}

  ErrCode Resonance::initPar1( FitParams* fitparams ) {
    ErrCode status;
    for ( daucontainer::const_iterator it = begin(); it != end(); ++it ) status |= ( *it )->initPar1( fitparams );
    return status;
  }

  ErrCode Resonance::initPar2( FitParams* fitparams ) {
    ErrCode status;
    for ( daucontainer::const_iterator it = begin(); it != end(); ++it ) status |= ( *it )->initPar2( fitparams );
    initMom( fitparams );
    return status;
  }

  std::string Resonance::parname( int index ) const { return ParticleBase::parname( index + 4 ); }

} // namespace DecayTreeFitter
