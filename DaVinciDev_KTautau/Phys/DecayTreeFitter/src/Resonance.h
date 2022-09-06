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
#ifndef __VTK_INTERNALRESONANCE_HH__
#define __VTK_INTERNALRESONANCE_HH__

#include "InternalParticle.h"

namespace DecayTreeFitter {
  class FitParams;

  class Resonance : public InternalParticle {
  public:
    Resonance( const LHCb::Particle& bc, const ParticleBase* mother, const Configuration& config );
    virtual ~Resonance();

    int         dim() const override { return 4; }
    int         type() const override { return kResonance; }
    std::string parname( int index ) const override;

    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;

    int  posIndex() const override { return mother()->posIndex(); }
    int  momIndex() const override { return index(); }
    int  lenIndex() const override { return -1; }
    bool hasPosition() const override { return false; }

  private:
  };

} // namespace DecayTreeFitter

#endif
