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
#ifndef INTERNALPARTICLE_H
#define INTERNALPARTICLE_H

#include "ParticleBase.h"
#include <vector>

namespace DecayTreeFitter {

  class InternalParticle : public ParticleBase {
  public:
    InternalParticle( const LHCb::Particle& bc, const ParticleBase* mother, const Configuration& config );

    int dim() const override { return mother() ? 8 : 7; }

    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;
    int     type() const override { return kInternalParticle; }

    // parameter definition
    int         posIndex() const override { return index(); }
    int         lenIndex() const override { return mother() ? index() + 3 : -1; }
    int         momIndex() const override { return mother() ? index() + 4 : index() + 3; }
    bool        hasEnergy() const override { return true; }
    bool        hasPosition() const override { return true; }
    std::string parname( int index ) const override;

    // constraints
    ErrCode projectKineConstraint( const FitParams&, Projection& ) const;
    ErrCode projectLifeTimeConstraint( const FitParams&, Projection& ) const;
    ErrCode projectMassConstraintTwoBody( const FitParams& fitparams, Projection& p ) const;
    ErrCode projectConstraint( Constraint::Type type, const FitParams& fitparams, Projection& p ) const override;

    // some of that other stuff
    void addToConstraintList( constraintlist& alist, int depth ) const override;

    // bool swapMotherDaughter(FitParams* fitparams, const ParticleBase* newmother) ;

  protected:
    ErrCode initMom( FitParams* fitparams ) const;

  private:
    bool m_lifetimeconstraint;
  };

} // namespace DecayTreeFitter

#endif
