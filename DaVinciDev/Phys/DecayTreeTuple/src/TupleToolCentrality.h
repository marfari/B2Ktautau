/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// $Id: TupleToolCentrality.h,v 1.0 2020-11-12 garcia $
#ifndef TUPLETOOLCENTRALITY_H
#define TUPLETOOLCENTRALITY_H 1

// Include files
// from Gaudi
#include "CaloDet/DeCalorimeter.h"
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IEventTupleTool.h" // Interface
#include "TupleToolCentrality.h"
class ITupleTool;

/** @class TupleToolCentrality TupleToolCentrality.h
 *
 * \brief Add Centrality info to DecayTreeTuple
 *
 * Tuple columns:
 * - Digit energy (verbose mode)
 * - Centrality percentile
 * - Number of participant nucleons
 * - Number of binary nucleon-nucleon collisions
 * - Impact parameter
 * - Error in number of participant nucleons
 * - Error in number of binary nucleon-nucleon collisions
 * - Error in the impact parameter
 * \sa DecayTreeTuple
 *  @author Felipe Garcia
 *  @date   2020-11-12
 */
class TupleToolCentrality : public TupleToolBase, virtual public IEventTupleTool {
public:
  /// Standard constructor
  TupleToolCentrality( const std::string& type, const std::string& name, const IInterface* parent );

  ~TupleToolCentrality() override = default; ///< Destructor

  StatusCode fill( Tuples::Tuple& ) override;
  StatusCode initialize() override;
  StatusCode finalize() override;

private:
  std::string                  m_DigitLocation = LHCb::CaloDigitLocation::Ecal;
  std::string                  m_CaloName      = "Ecal";
  Gaudi::Property<std::string> m_CollSyst{this, "CollSyst", "PbPb",
                                          "Collision system. By default 'PbPb' (case insensitive)."};
  std::string                  m_CaloLocation = DeCalorimeterLocation::Ecal;

  std::vector<int>           m_index;
  std::vector<char>          m_calo;
  std::vector<unsigned char> m_area;
  std::vector<unsigned char> m_row;
  std::vector<unsigned char> m_column;
  std::vector<float>         m_es;

  double m_TES;
  double m_centrality;
  double m_npart_final;
  double m_ncoll_final;
  double m_b_final;
  int    m_count;
  int m_i = 0; // counter for specific centrality percentile in 1% intervals. The higher m_i the more central the event
  int m_j = 0; // counter for specific centrality class asked by the user. The higher m_j the more central the event;
  int m_switch = 1;
  Gaudi::Property<int> m_nclasses{this, "NClasses", 10, "Requested number of centrality classes. By default 10."};
  std::vector<float>   m_ecalcuts;
  std::vector<float>   m_npart;
  std::vector<float>   m_npart_sys;
  std::vector<float>   m_npart_stat;
  double               m_npart_error;
  std::vector<float>   m_ncoll;
  std::vector<float>   m_ncoll_sys;
  std::vector<float>   m_ncoll_stat;
  double               m_ncoll_error;
  std::vector<float>   m_b;
  std::vector<float>   m_b_sys;
  std::vector<float>   m_b_stat;
  double               m_b_error;
  std::vector<float>   m_npart_temp = {};
  std::vector<float>   m_ncoll_temp = {};
  std::vector<float>   m_b_temp     = {};

  Gaudi::Property<unsigned long> m_maxSize{this, "MaxDigits", 8192,
                                           "Vector size to read the ECAL cells. By default 8192."};

  Gaudi::Property<bool> m_debug{
      this, "Debugging", false,
      "Boolean to print some debugging variables (counters for the classification m_i and m_j). By default false."};
};
#endif //
