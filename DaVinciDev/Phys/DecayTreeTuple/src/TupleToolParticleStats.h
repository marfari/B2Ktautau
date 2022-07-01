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
#ifndef TUPLETOOLRECOSTATS_H
#define TUPLETOOLRECOSTATS_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IEventTupleTool.h" // Interface

#include "Event/RecSummary.h"
#include "Event/Track.h"

/** @class TupleToolParticleStats TupleToolParticleStats.h
 *
 *  Fills Reco stats, from RecSummary
 *
 *  @author Patrick Koppenburg, Fatima Soomro, Jibo He
 *  @date   2009-02-11
 */
class TupleToolParticleStats : public TupleToolBase, virtual public IEventTupleTool {

public:
  /// Standard constructor
  TupleToolParticleStats( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~TupleToolParticleStats();          ///< Destructor
  StatusCode fill( Tuples::Tuple& ) override; ///< Fill tuple
  //  int unusedVelo();
private:
  std::vector<std::string> m_locations; ///< Hlt locations to look at
};

#endif // TUPLETOOLRECOSTATS_H
