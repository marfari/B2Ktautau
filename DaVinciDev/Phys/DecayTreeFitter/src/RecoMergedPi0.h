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
#ifndef _VTK_RECOMERGEDPI0_HH_
#define _VTK_RECOMERGEDPI0_HH_

#include "ParticleBase.h"

namespace LHCb {
  class CaloMomentum;
}

namespace DecayTreeFitter {

  class RecoMergedPi0 : public ParticleBase {
  public:
    RecoMergedPi0( const LHCb::Particle& bc, const ParticleBase* mother );
    ~RecoMergedPi0();

    // the number of parameters
    int dim() const override { return hasMassConstraint() ? 3 : 4; }

    // the number of 'measurements'
    int dimM() const { return dim(); }

    // does it have an energy component?
    bool hasEnergy() const override { return !hasMassConstraint(); }

    // project the constraint
    ErrCode projectPi0Constraint( const FitParams&, Projection& ) const;
    ErrCode projectConstraint( Constraint::Type type, const FitParams& fitparams, Projection& p ) const override {
      ErrCode status;
      switch ( type ) {
      case Constraint::btacomposite:
        status |= projectPi0Constraint( fitparams, p );
        break;
      default:
        status |= ParticleBase::projectConstraint( type, fitparams, p );
      }
      return status;
    }

    std::string parname( int index ) const override { return ParticleBase::parname( index + 4 ); }

    ErrCode initParPi0( FitParams* );
    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;
    ErrCode initCov( FitParams* ) const override;
    int     type() const override { return kRecoMergedPi0; }
    double  chiSquare( const FitParams* fitparams ) const override;

    int momIndex() const override { return index(); }

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      alist.push_back( Constraint( this, Constraint::btacomposite, depth, dimM() ) );
    }

  protected:
    LHCb::CaloMomentum* m_calomom;
  };

} // namespace DecayTreeFitter

#endif
