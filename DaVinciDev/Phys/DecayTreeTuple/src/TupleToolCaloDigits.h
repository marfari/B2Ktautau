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
// $Id: TupleToolCaloDigits.h,v 0.1 2015-08-05 bursche $
#ifndef TUPLETOOLCALODIGITS_H
#define TUPLETOOLCALODIGITS_H 1

// Include files
// from Gaudi
#include "CaloDet/DeCalorimeter.h"
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IEventTupleTool.h" // Interface
#include "TupleToolCaloDigits.h"
class ITupleTool;

/** @class TupleToolCaloDigits TupleToolCaloDigits.h
 *
 * \brief Add Calorimter digits to DecayTreeTuple
 *
 * Tuple columns:
 * - Digit energy
 * - Cell ID
 * - Cell area, column and row
 * \sa DecayTreeTuple
 *  @author Albert Bursche
 *  @date   2015-08-05
 */
class TupleToolCaloDigits : public TupleToolBase, virtual public IEventTupleTool {
public:
  /// Standard constructor
  TupleToolCaloDigits( const std::string& type, const std::string& name, const IInterface* parent );

  ~TupleToolCaloDigits() override = default; ///< Destructor

  StatusCode fill( Tuples::Tuple& ) override;
  StatusCode initialize() override;
  StatusCode finalize() override;

private:
  std::string m_DigitLocation = LHCb::CaloDigitLocation::Spd;
  std::string m_CaloName      = "Spd";
  std::string m_CaloLocation  = DeCalorimeterLocation::Spd;

  std::vector<int>           m_index;
  std::vector<char>          m_calo;
  std::vector<unsigned char> m_area;
  std::vector<unsigned char> m_row;
  std::vector<unsigned char> m_column;
  std::vector<float>         m_xs;
  std::vector<float>         m_ys;
  std::vector<float>         m_zs;
  std::vector<float>         m_es;
  std::vector<float>         m_ets;

  unsigned long m_maxSize        = 4096;
  bool          m_auto_configure = false;
};
#endif //
