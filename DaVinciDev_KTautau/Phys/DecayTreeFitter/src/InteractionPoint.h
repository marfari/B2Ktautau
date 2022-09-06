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
#ifndef DECAYTREEFITTER_INTERACTIONPOINT_H
#define DECAYTREEFITTER_INTERACTIONPOINT_H

#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "ParticleBase.h"

namespace LHCb {
  class VertexBase;
}

namespace DecayTreeFitter {

  class InteractionPoint : public ParticleBase {
  public:
    InteractionPoint( const LHCb::VertexBase& ipvertex, const LHCb::Particle& daughter, const Configuration& config );

    int     dim() const override { return 3; } // (x,y,z)
    ErrCode initPar1( FitParams* ) override;
    ErrCode initPar2( FitParams* ) override;
    ErrCode initCov( FitParams* ) const override;

    int type() const override { return kInteractionPoint; }

    double chiSquare( const FitParams* par ) const override;

    ErrCode projectIPConstraint( const FitParams& fitpar, Projection& ) const;
    ErrCode projectConstraint( Constraint::Type, const FitParams&, Projection& ) const override;

    void addToConstraintList( constraintlist& alist, int depth ) const override;

    int posIndex() const override { return index(); }

  private:
    Gaudi::Vector3      m_ipPos;    // interaction point position
    Gaudi::SymMatrix3x3 m_ipCov;    // cov matrix
    Gaudi::SymMatrix3x3 m_ipCovInv; // inverse of cov matrix
  };

} // namespace DecayTreeFitter

#endif
