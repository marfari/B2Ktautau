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
// $Id: TupleToolVeto.h
#ifndef _TUPLETOOLVETO_H
#define _TUPLETOOLVETO_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
// Interfaces
#include "Kernel/IParticleTupleTool.h"
#include "Kernel/IParticleVeto.h"
#include <memory>

//============================================================================

class TupleToolVeto : public TupleToolBase, virtual public IParticleTupleTool {

public:
  /// Standard constructor
  TupleToolVeto( const std::string& type, const std::string& name, const IInterface* parent );

  /// Initialise
  StatusCode initialize() override;

  /// Destructor
  virtual ~TupleToolVeto() {}

public:
  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  IParticleVeto*                                  m_check;
  std::string                                     m_part;
  unsigned int                                    m_pid;
  std::map<std::string, std::vector<std::string>> m_veto;
  std::map<std::string, std::vector<std::string>> m_vetoOther;
};

#endif // _TUPLETOOLVETO_H
