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
#ifndef RECOPARTICLE_HH
#define RECOPARTICLE_HH

#include "ParticleBase.h"

namespace DecayTreeFitter {

  class RecoParticle : public ParticleBase {
  public:
    RecoParticle( const LHCb::Particle& bc, const ParticleBase* mother );
    virtual ~RecoParticle();

    virtual int dimM() const = 0; // dimension of the measurement
    ErrCode     initPar1( FitParams* ) override { return ErrCode::success; }
    // virtual ErrCode initCov(FitParams*) const ;
    std::string parname( int index ) const override;
    int         dim() const override { return 3; } //(px,py,pz)

    int  momIndex() const override { return index(); }
    bool hasEnergy() const override { return false; }

    virtual ErrCode projectRecoConstraint( const FitParams& fitparams, Projection& p ) const = 0;
    ErrCode         projectConstraint( Constraint::Type, const FitParams&, Projection& ) const override;
    double          chiSquare( const FitParams* fitparams ) const override;
  };

} // namespace DecayTreeFitter
#endif
