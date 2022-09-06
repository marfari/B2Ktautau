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

#ifndef TUPLETOOLRADUPDOWNASYM_H
#define TUPLETOOLRADUPDOWNASYM_H 1

// Include files
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "GaudiKernel/Vector4DTypes.h"
#include "Kernel/IParticleTupleTool.h" // Interface
#include "Math/Boost.h"

/** @class TupleToolRadUpDownAsym TupleToolRadUpDownAsym.h
 * @Original author: C. Potterat
 * @Code Mantainer: A. Puig
 * @Date: 05-02-2019
 *
 * \brief Add relative angles of a photon with respect to different combination of daughters to DecayTreeTuple.
 *	  For 3 body decays plus photon
 *
 * Tuple columns:
 * - head_UpDownCosTheta(XY): Angle between gamma and combination X-Y of particles (x3 combinations)
 * - head_UpDownID(X): ID of each particle (x3 particles)
 */

class TupleToolRadUpDownAsym : public TupleToolBase, virtual public IParticleTupleTool {

public:
  /// Standard constructor
  TupleToolRadUpDownAsym( const std::string& type, const std::string& name, const IInterface* parent );

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  bool UpDownTheta( const LHCb::Particle* top, const LHCb::Particle* part, std::vector<double>& angleValue,
                    std::vector<double>& pidValue ) const;

  bool                               isEndTree( const LHCb::Particle* p ) const;
  std::vector<const LHCb::Particle*> getTree( const LHCb::Particle* P ) const;
};

#endif // TUPLETOOLRADUPDOWNASYM_H
