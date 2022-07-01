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
#ifndef _VTK_JETMOMENTUM_HH_
#define _VTK_JETMOMENTUM_HH_

#include "ParticleBase.h"
#include <GaudiKernel/GenericVectorTypes.h>
#include <GaudiKernel/SymmetricMatrixTypes.h>

namespace DecayTreeFitter {
  class JetMomentum : public ParticleBase {
  public:
    JetMomentum( const LHCb::Particle& bc, const ParticleBase* mother );
    virtual ~JetMomentum();

    // the number of parameters
    int dim() const override { return 4; } // px,py,pz,E)

    // the number of 'measurements'
    int     dimM() const { return 4; }
    ErrCode projectJetMomentum( const FitParams&, Projection& ) const;
    ErrCode projectConstraint( Constraint::Type, const FitParams&, Projection& ) const override;

    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override { return ErrCode::success; }
    int     type() const override { return kJetMomentum; }

    int  momIndex() const override { return index(); }
    bool hasEnergy() const override { return true; }

    void   updCache();
    double chiSquare( const FitParams* fitparams ) const override;

    std::string parname( int index ) const override { return ParticleBase::parname( index + 4 ); }

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      alist.push_back( Constraint( this, Constraint::externalmomentum, depth, dimM() ) );
    }

  protected: // I hate this, so we need to change the design ...
    // cache
    Gaudi::Vector4      m_m; // 'measurement' (px,py,pz,E)
    Gaudi::SymMatrix4x4 m_V; // covariance in measurement
  };

} // namespace DecayTreeFitter

#endif
