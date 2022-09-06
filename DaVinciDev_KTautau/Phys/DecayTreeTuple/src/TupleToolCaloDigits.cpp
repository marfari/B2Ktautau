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
// $Id: TupleToolCaloDigits.cpp,v 0.1 2015-08-05 13:50:45 bursche $
// Include files

// local
#include "TupleToolCaloDigits.h"

#include "GaudiAlg/ITupleTool.h"
#include "GaudiAlg/Tuple.h"

#include "Event/CaloDigit.h"

#include "GaudiKernel/IRegistry.h" // IOpaqueAddress

#include "CaloDet/DeCalorimeter.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolCaloDigits
//
// 2015-08-05 Albert Bursche
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_COMPONENT( TupleToolCaloDigits )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolCaloDigits::TupleToolCaloDigits( const std::string& type, const std::string& name, const IInterface* parent )
    : TupleToolBase( type, name, parent ) {
  declareInterface<IEventTupleTool>( this );
  declareProperty( "DigitLocation", m_DigitLocation );
  declareProperty( "CaloName", m_CaloName );
  declareProperty( "CaloLocation", m_CaloLocation );
  declareProperty( "MaxDigits", m_maxSize );
  declareProperty( "AutoConfigure", m_auto_configure );
}
//=============================================================================

StatusCode TupleToolCaloDigits::initialize() {
  if ( m_auto_configure ) {
    if ( m_CaloName.compare( "Spd" ) == 0 ) {
      m_CaloLocation  = DeCalorimeterLocation::Spd;
      m_DigitLocation = LHCb::CaloDigitLocation::Spd;
      return StatusCode::SUCCESS;
    }
    if ( m_CaloName.compare( "Prs" ) == 0 ) {
      m_CaloLocation  = DeCalorimeterLocation::Prs;
      m_DigitLocation = LHCb::CaloDigitLocation::Prs;
      return StatusCode::SUCCESS;
    }
    if ( m_CaloName.compare( "Ecal" ) == 0 ) {
      m_CaloLocation  = DeCalorimeterLocation::Ecal;
      m_DigitLocation = LHCb::CaloDigitLocation::Ecal;
      return StatusCode::SUCCESS;
    }
    if ( m_CaloName.compare( "Hcal" ) == 0 ) {
      m_CaloLocation  = DeCalorimeterLocation::Hcal;
      m_DigitLocation = LHCb::CaloDigitLocation::Hcal;
      return StatusCode::SUCCESS;
    }
    // Feel free to add more locations from Event/CaloDigit.h
    fatal() << "Invalid CaloName. Mus be one of  Spd,Prs,Ecal,Hcal if AutoConfigure is used (i.e. true)" << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}
StatusCode TupleToolCaloDigits::finalize() { return StatusCode::SUCCESS; }
StatusCode TupleToolCaloDigits::fill( Tuples::Tuple& tuple ) {
  const std::string prefix = fullName();
  LHCb::CaloDigits* digits = getIfExists<LHCb::CaloDigits>( m_DigitLocation );
  if ( !digits ) return StatusCode::SUCCESS;
  DeCalorimeter* calo = getDet<DeCalorimeter>( m_CaloLocation );
  for ( const auto& digit : *digits ) {
    const LHCb::CaloCellID cellID = digit->cellID();
    m_index.push_back( cellID.index() );
    m_calo.push_back( cellID.calo() );
    m_area.push_back( cellID.area() );
    m_row.push_back( cellID.row() );
    m_column.push_back( cellID.col() );

    auto p = calo->cellCenter( cellID );
    m_xs.push_back( p.X() );
    m_ys.push_back( p.Y() );
    m_zs.push_back( p.Z() );
    m_es.push_back( digit->e() );
    double sinTheta = sqrt( p.X() * p.X() + p.Y() * p.Y() ) / sqrt( p.X() * p.X() + p.Y() * p.Y() + p.Z() * p.Z() );
    m_ets.push_back( digit->e() * sinTheta );
    if ( msgLevel( MSG::DEBUG ) ) debug() << cellID.toString() << " has an energy of " << m_es.back() << " \n";
  }
  if ( msgLevel( MSG::DEBUG ) )
    debug() << " saved " << m_index.size() << " digits to n tuple " + m_CaloName + "Digit." << endmsg;
  if ( m_verbose ) {
    tuple
        ->farray( m_extraName + m_CaloName + "DigitIndex", std::begin( m_index ), std::end( m_index ),
                  m_extraName + m_CaloName + "DigitNi", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "DigitCalo", std::begin( m_calo ), std::end( m_calo ),
                  m_extraName + m_CaloName + "DigitNc", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "DigitArea", std::begin( m_area ), std::end( m_area ),
                  m_extraName + m_CaloName + "DigitNa", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "DigitRow", std::begin( m_row ), std::end( m_row ),
                  m_extraName + m_CaloName + "DigitNr", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "DigitColumn", std::begin( m_column ), std::end( m_column ),
                  m_extraName + m_CaloName + "DigitNco", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "Digit_X", std::begin( m_xs ), std::end( m_xs ),
                  m_extraName + m_CaloName + "DigitNx", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "Digit_Y", std::begin( m_ys ), std::end( m_ys ),
                  m_extraName + m_CaloName + "DigitNy", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "Digit_Z", std::begin( m_zs ), std::end( m_zs ),
                  m_extraName + m_CaloName + "DigitNz", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "DigitEnergy", std::begin( m_es ), std::end( m_es ),
                  m_extraName + m_CaloName + "DigitNe", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tuple
        ->farray( m_extraName + m_CaloName + "DigitTransverseEnergy", std::begin( m_ets ), std::end( m_ets ),
                  m_extraName + m_CaloName + "DigitNe", m_maxSize )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  }
  tuple
      ->column( m_extraName + m_CaloName + "DigitEnergySum",
                std::accumulate( std::begin( m_es ), std::end( m_es ), 0.0 ) )
      .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  tuple
      ->column( m_extraName + m_CaloName + "DigitTransverseEnergySum",
                std::accumulate( std::begin( m_ets ), std::end( m_ets ), 0.0 ) )
      .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  m_index.clear();
  m_calo.clear();
  m_area.clear();
  m_row.clear();
  m_column.clear();
  m_xs.clear();
  m_ys.clear();
  m_zs.clear();
  m_es.clear();
  m_ets.clear();

  return StatusCode::SUCCESS;
}
