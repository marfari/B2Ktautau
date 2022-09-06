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

// from Gaudi
#include "TupleToolParticleStats.h"
#include "Event/ODIN.h"
#include "Event/Track.h"
#include "Event/VeloCluster.h"

// local
#include "TupleToolParticleStats.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolParticleStats
//
// 2009-02-11 : Patrick Koppenburg
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolParticleStats::TupleToolParticleStats( const std::string& type, const std::string& name,
                                                const IInterface* parent )
    : TupleToolBase( type, name, parent ), m_locations() {
  declareInterface<IEventTupleTool>( this );
  declareProperty( "InputLocations", m_locations, "Locations to look at" );
  m_locations.push_back( "StdLooseAllPhotons" );
}

//=============================================================================
// Destructor
//=============================================================================
TupleToolParticleStats::~TupleToolParticleStats() {}
//=============================================================================

//=============================================================================
StatusCode TupleToolParticleStats::fill( Tuples::Tuple& tup ) {
  const std::string prefix = fullName();
  bool              test   = true;
  for ( std::vector<std::string>::const_iterator l = m_locations.begin(); l != m_locations.end(); ++l ) {
    if ( !( exist<LHCb::Particles>( "/Event/Phys/" + *l + "/Particles" ) ) ) {
      error() << " your configuration is problematic, tthe location you requested: " << endmsg;
      error() << " >>>  " << ( "/Event/Phys/" + *l + "/Particles" ) << " <<<  doesnt seem to have been filled."
              << endmsg;
      continue;
    }
    test &= tup->column( prefix + "NumberOf" + *l, number<LHCb::Particles>( "/Event/Phys/" + *l + "/Particles" ) );
  }
  return StatusCode( test );
}

// Declaration of the Tool Factory
DECLARE_COMPONENT( TupleToolParticleStats )
