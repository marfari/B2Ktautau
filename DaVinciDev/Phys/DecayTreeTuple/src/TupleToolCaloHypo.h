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
// $Id: TupleToolCalohypo.h,v 1.1 2010-05-21 17:03:11 odescham Exp $
#ifndef _TUPLETOOLECALOHYPO_H
#define _TUPLETOOLECALOHYPO_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
// Interface
#include "CaloDAQ/ICaloDigitFilterTool.h"
#include "CaloInterfaces/ICaloRelationsGetter.h"
#include "CaloUtils/ICaloElectron.h"
#include "Event/Track.h"
#include "Kernel/IParticleTupleTool.h"
#include "Relations/IRelationWeighted2D.h"

#include "CaloInterfaces/ICaloHypoEstimator.h"

//============================================================================
class TupleToolCaloHypo : public TupleToolBase, virtual public IParticleTupleTool {
  //==========================================================================
public:
  // Standard constructor
  TupleToolCaloHypo( const std::string& type, const std::string& name, const IInterface* parent );

  StatusCode initialize() override;
  virtual ~TupleToolCaloHypo(){};
  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  double                   incidence( const LHCb::CaloHypo* hypo, bool straight = false ) const;
  std::vector<std::string> m_dataList;
  bool                     useData( int id, int mask = 0x3 );
  ICaloHypoEstimator*      m_estimator;
  bool                     m_brem;
  bool                     m_charged;
  bool                     m_daugh;
  ICaloRelationsGetter*    m_tables;
  ICaloElectron*           m_caloElectron;
  ICaloDigitFilterTool*    m_pileup;
};
#endif // _TUPLETOOLECALOHYPO_H
