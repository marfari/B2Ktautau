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
// Include files

// from LHCb
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IEventTupleTool.h"

#include "LoKi/BeamSpot.h"

/** @class TupleToolBeamSpot DecayTreeTuple/TupleToolBeamSpot.cpp
 *
 * Add beamspot position (from the Velo resolver position)
 *
 * Tuple columns:
 *  - BeamX
 *  - BeamY
 *
 * @see LoKi::BeamSpot
 *
 * @author Pieter David <pieter.david@nikhef.nl>
 * @date   2015-07-10
 */
class TupleToolBeamSpot : public TupleToolBase, virtual public IEventTupleTool {
public:
  /// Standard constructor
  TupleToolBeamSpot( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~TupleToolBeamSpot(); ///< Destructor

  StatusCode initialize() override;
  StatusCode finalize() override;

  // Interface methods
  StatusCode fill( Tuples::Tuple& tuple ) override;

private:
  double                          m_bound;
  std::unique_ptr<LoKi::BeamSpot> m_beamspot;
};
