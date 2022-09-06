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
#ifndef RECORESONANCE_H
#define RECORESONANCE_H

#include "RecoComposite.h"

namespace DecayTreeFitter {

  class RecoResonance : public RecoComposite {
  public:
    RecoResonance( const LHCb::Particle& bc, const ParticleBase* mother );
    virtual ~RecoResonance();

    int dim() const override { return hasEnergy() ? 4 : 3; } // (px,py,pz,(E))

    ErrCode projectConstraint( Constraint::Type, const FitParams&, Projection& ) const override;
    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;
    int     type() const override { return kRecoResonance; }

    int posIndex() const override { return mother()->posIndex(); }
    int momIndex() const override { return index(); }
    int lenIndex() const override { return -1; }

    std::string parname( int index ) const override;

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      alist.push_back( Constraint( this, Constraint::btaresonance, depth, dimM() ) );
    }

  private:
  };

} // namespace DecayTreeFitter

#endif
