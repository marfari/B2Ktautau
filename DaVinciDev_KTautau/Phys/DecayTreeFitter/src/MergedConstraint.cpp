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
#include "MergedConstraint.h"
#include "FitParams.h"
#include "ParticleBase.h"

namespace DecayTreeFitter {

  ErrCode MergedConstraint::project( const FitParams& fitpar, Projection& p ) const {
    ErrCode status;
    for ( constraintlist::const_iterator it = m_list.begin(); it != m_list.end(); ++it ) {
      status |= ( *it )->project( fitpar, p );
      p.incrementOffset( ( *it )->dim() );
    }

    return status;
  }

  void MergedConstraint::print( std::ostream& os ) const {
    os << "Merged constraint: " << std::endl;
    for ( constraintlist::const_iterator it = m_list.begin(); it != m_list.end(); ++it ) {
      os << "          ";
      ( *it )->print( os );
    }
  }

} // namespace DecayTreeFitter
