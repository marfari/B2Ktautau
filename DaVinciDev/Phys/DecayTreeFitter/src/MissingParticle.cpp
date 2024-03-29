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
/**
 To finish the implementation, this is what we need to do: We would like to make the dimension of the missing particle
 depend on the mass constraint. That means that any call to 'setMassConstraint' should invalidate the 'indices'. The
 creation of indices is done with 'DecayChain::updateIndex'. We could call this from 'initPar' (or from any other place
 that initialized the constraints. That also means that 'setMassConstraint' should change Fitter.m_fitStatus to
 Unfitted'. Finally, to make sure that we don't break anything, we need a test case.

*/

#include "MissingParticle.h"
#include "Event/Particle.h"
#include "FitParams.h"

namespace DecayTreeFitter {

  extern int vtxverbose;

  MissingParticle::MissingParticle( const LHCb::Particle& bc, const ParticleBase* mother )
      : ParticleBase( bc, mother ) {
    // this will be one of the very few particles for which we adjust
    // the dimension if there is a constraint
  }

  MissingParticle::~MissingParticle() {}

  ErrCode MissingParticle::initPar1( FitParams* fitpar ) {
    // take them from the bc
    Gaudi::LorentzVector p4       = particle().momentum();
    int                  momindex = momIndex();

    fitpar->par()( momindex + 1 ) = p4.x();
    fitpar->par()( momindex + 2 ) = p4.y();
    fitpar->par()( momindex + 3 ) = p4.z();
    if ( hasEnergy() ) fitpar->par()( momindex + 4 ) = p4.t();
    return ErrCode();
  }

  std::string MissingParticle::parname( int index ) const { return ParticleBase::parname( index + 4 ); }
} // namespace DecayTreeFitter
