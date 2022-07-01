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
#include "VtxErrCode.h"

namespace DecayTreeFitter {
  std::ostream& operator<<( std::ostream& os, const ErrCode& ec ) {
    unsigned int flag = ec.flag();
    os << "flag=" << flag << " ";
    if ( flag ) {
      if ( flag & ErrCode::pocafailure ) os << "pocafailure ";
      if ( flag & ErrCode::baddistance ) os << "baddistance ";
      if ( flag & ErrCode::inversionerror ) os << "inversionerror ";
      if ( flag & ErrCode::badsetup ) os << "badsetup ";
      if ( flag & ErrCode::divergingconstraint ) os << "divergingconstraint ";
      if ( flag & ErrCode::slowdivergingfit ) os << "slowdivergingfit ";
      if ( flag & ErrCode::fastdivergingfit ) os << "fastdivergingfit ";
      if ( flag & ErrCode::fastdivergingfit ) os << "filtererror ";
    } else {
      os << "success ";
    }
    return os;
  }
} // namespace DecayTreeFitter
