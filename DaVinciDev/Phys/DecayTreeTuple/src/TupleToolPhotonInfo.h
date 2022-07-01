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
// $Id: TupleToolPhotonInfo.h,v 1.1 2010-05-21 17:03:11 odescham Exp $
#ifndef _TUPLETOOLEPHOTONINFO_H
#define _TUPLETOOLEPHOTONINFO_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
// Interface
#include "Kernel/IParticleTupleTool.h"

#include "CaloUtils/CaloParticle.h"

//============================================================================
class TupleToolPhotonInfo : public TupleToolBase, virtual public IParticleTupleTool {
  //==========================================================================
public:
  // Standard constructor
  TupleToolPhotonInfo( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~TupleToolPhotonInfo(){}; ///< Destructor

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  int m_PhotonID;
};
#endif // _TUPLETOOLEPHOTONINFO_H
