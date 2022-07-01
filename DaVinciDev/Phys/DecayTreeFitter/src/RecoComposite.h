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
#ifndef _VTK_EXTERNALBTAPARTICLE_HH_
#define _VTK_EXTERNALBTAPARTICLE_HH_

#include "ParticleBase.h"
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Vector.h>

namespace DecayTreeFitter {

  class RecoComposite : public ParticleBase {
  public:
    RecoComposite( const LHCb::Particle& bc, const ParticleBase* mother );
    virtual ~RecoComposite();

    // the number of parameters
    int dim() const override { return m_hasEnergy ? 8 : 7; } // (x,y,z,t,px,py,pz,(E))

    // the number of 'measurements'
    int     dimM() const { return m_hasEnergy ? 7 : 6; }
    ErrCode projectRecoComposite( const FitParams&, Projection& ) const;
    ErrCode projectConstraint( Constraint::Type, const FitParams&, Projection& ) const override;

    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;
    int     type() const override { return kRecoComposite; }

    int posIndex() const override { return index(); }
    int lenIndex() const override { return index() + 3; }
    int momIndex() const override { return index() + 4; }

    bool hasEnergy() const override { return m_hasEnergy; }
    bool hasPosition() const override { return true; }

    virtual void updCache();
    double       chiSquare( const FitParams* fitparams ) const override;

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      alist.push_back( Constraint( this, Constraint::btacomposite, depth, dimM() ) );
      alist.push_back( Constraint( this, Constraint::geometric, depth, 3 ) );
    }

  protected: // I hate this, so we need to change the design ...
    // cache
    CLHEP::HepVector    m_m;       // 'measurement' (x,y,zpx,py,pz,E)
    CLHEP::HepSymMatrix m_matrixV; // covariance in measurement
    bool                m_hasEnergy;
  };

} // namespace DecayTreeFitter

#endif
