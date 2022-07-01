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
#ifndef VTK_MERGEDCONSTRAINT_HH
#define VTK_MERGEDCONSTRAINT_HH

#include "Constraint.h"
#include <vector>

namespace DecayTreeFitter {
  class MergedConstraint : public Constraint {
  public:
    typedef std::vector<Constraint*> constraintlist;

    MergedConstraint() : Constraint( Constraint::merged ) {}
    virtual ~MergedConstraint() {}

    MergedConstraint( const constraintlist& list ) : Constraint( Constraint::merged ), m_list( list ) {
      int d( 0 );
      for ( constraintlist::iterator it = m_list.begin(); it != m_list.end(); ++it ) d += ( *it )->dim();
      setDim( d );
    }

    ErrCode project( const FitParams& fitpar, Projection& p ) const override;

    void push_back( Constraint* c ) {
      m_list.push_back( c );
      setDim( dim() + c->dim() );
      setNIter( std::max( nIter(), c->nIter() ) );
    }

    void print( std::ostream& os = std::cout ) const override;

  private:
    constraintlist m_list;
  };

} // namespace DecayTreeFitter

#endif
