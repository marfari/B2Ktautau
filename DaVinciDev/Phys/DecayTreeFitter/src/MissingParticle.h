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
#ifndef VTK_MISSINGPARTICLE_HH
#define VTK_MISSINGPARTICLE_HH

#include "ParticleBase.h"

namespace DecayTreeFitter {

  class MissingParticle : public ParticleBase {
  public:
    MissingParticle( const LHCb::Particle& bc, const ParticleBase* mother );
    virtual ~MissingParticle();

    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override { return ErrCode::success; }

    std::string parname( int index ) const override;
    int         dim() const override { return hasMassConstraint() ? 3 : 4; }
    int         momIndex() const override { return index(); }
    bool        hasEnergy() const override { return hasMassConstraint() ? false : true; }
    int         type() const override { return kMissingParticle; }
    void        addToConstraintList( constraintlist& /*alist*/, int /*depth*/ ) const override {}
  };

} // namespace DecayTreeFitter
#endif
