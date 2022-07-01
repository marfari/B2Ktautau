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
#ifndef __VTK_CONVERTEDPHOTON_HH__
#define __VTK_CONVERTEDPHOTON_HH__

#include "InternalParticle.h"

namespace DecayTreeFitter {
  class ConvertedPhoton : public InternalParticle {
  public:
    ConvertedPhoton( const LHCb::Particle& bc, const ParticleBase* mother, const Configuration& config );
    virtual ~ConvertedPhoton() {}

    int type() const override { return kConvertedPhoton; }

    ErrCode projectConversionPositionConstraint( const FitParams&, Projection& ) const;
    ErrCode projectConversionMassConstraint( const FitParams& fitparams, Projection& p ) const;
    ErrCode projectConstraint( Constraint::Type type, const FitParams& fitparams, Projection& p ) const override;
    void    addToConstraintList( constraintlist& alist, int depth ) const override;

    ErrCode initCov( FitParams* fitpars ) const override;

  private:
    void updateCache( const FitParams& pars );

  private:
    double m_conversionZ;
    double m_conversionZCov;
  };

} // namespace DecayTreeFitter

#endif
