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
#ifndef __VTK_HELIXUTILS_HH__
#define __VTK_HELIXUTILS_HH__

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include <string>

typedef Gaudi::XYZVector BField;

namespace DecayTreeFitter {

  class VtkHelixUtils {
  public:
    enum VertexCoor { in_x = 0, in_y, in_z, in_px, in_py, in_pz };
    enum HelixCoor { ex_d0 = 0, ex_phi0 = 1, ex_omega = 2, ex_z0 = 3, ex_tanDip = 4, ex_flt = 5 };

    static void vertexFromHelix( const CLHEP::HepVector& helixpar, const BField& fieldmap, CLHEP::HepVector& vertexpar,
                                 int& charge );

    static void helixFromVertex( const CLHEP::HepVector& vertexpar, int charge, const BField& fieldmap,
                                 CLHEP::HepVector& helixpar, CLHEP::HepMatrix& jacobian );

    static void helixFromVertexNumerical( const CLHEP::HepVector& vertexpar, int charge, const BField& fieldmap,
                                          CLHEP::HepVector& helixpar, CLHEP::HepMatrix& jacobian );

    static std::string helixParName( int i );
    static std::string vertexParName( int i );
    static void        printHelixPar( const CLHEP::HepVector& helixpar );
    static void        printVertexPar( const CLHEP::HepVector& vertexpar, int charge );

    static double helixPoca( const CLHEP::HepVector& helixpar1, const CLHEP::HepVector& helixpar2, double& flt1,
                             double& flt2, Gaudi::XYZPoint& v, bool parallel = false );
    static double helixPoca( const CLHEP::HepVector& helixpar, const Gaudi::XYZPoint& point, double& flt );
    static double phidomain( const double phi );
  };

} // namespace DecayTreeFitter

#endif
