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
// $Id: TupleToolBremInfo.h
#ifndef _TUPLETOOLBREMINFO_H
#define _TUPLETOOLBREMINFO_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
// Interfaces
#include "Kernel/IBremAdder.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/IParticleTupleTool.h"
#include "Kernel/ParticleProperty.h"
#include <memory>

//============================================================================
class TupleToolBremInfo : public TupleToolBase, virtual public IParticleTupleTool {
  //==========================================================================

public:
  // Standard constructor
  TupleToolBremInfo( const std::string& type, const std::string& name, const IInterface* parent );

  StatusCode initialize() override;
  virtual ~TupleToolBremInfo(){}; ///< Destructor

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  const LHCb::ParticleProperty* property( const LHCb::ParticleID pid ) {
    return ( m_ppsvc ) ? m_ppsvc->find( pid ) : NULL;
  };
  LHCb::IParticlePropertySvc* m_ppsvc;
  IBremAdder*                 m_adder;
  std::vector<std::string>    m_parts;
  std::vector<unsigned int>   m_pids;
  bool                        m_dst;
};

#endif // _TUPLETOOLBREMINFO_H
