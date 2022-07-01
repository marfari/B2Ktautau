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
#ifndef RECOPHOTON_H
#define RECOPHOTON_H

#include "RecoParticle.h"

namespace DecayTreeFitter {

  class RecoPhoton : public RecoParticle {
  public:
    RecoPhoton( const LHCb::Particle& bc, const ParticleBase* mother );
    virtual ~RecoPhoton();

    int     dimM() const override { return 3; }
    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;

    ErrCode initCov( FitParams* ) const override;
    int     type() const override { return kRecoPhoton; }
    ErrCode projectRecoConstraint( const FitParams&, Projection& ) const override;
    ErrCode updCache();

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      alist.push_back( Constraint( this, Constraint::photon, depth, dimM() ) );
    }

  private:
    virtual ErrCode     initParPhoton( FitParams*, const Gaudi::XYZPoint& motherpos ) const;
    double              m_z;
    Gaudi::Vector3      m_m;
    Gaudi::SymMatrix3x3 m_V;
  };

} // namespace DecayTreeFitter
#endif
