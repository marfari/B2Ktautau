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
#ifndef FSOOMRO_TUPLETOOLDIRA_H
#define FSOOMRO_TUPLETOOLDIRA_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h" // Interface
#include "Math/SMatrix.h"

#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>

#include "Event/Particle.h"
#include "Event/RecVertex.h"
#include "Event/Vertex.h"

class TupleToolDira : public TupleToolBase, virtual public IParticleTupleTool {

private:
  typedef ROOT::Math::SMatrix<double, 10, 10, ROOT::Math::MatRepSym<double, 10>> SymMatrix10x10;

public:
  /// Standard constructor
  TupleToolDira( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~TupleToolDira(){}; ///< Destructor

  StatusCode initialize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  StatusCode fillDIRAError( const LHCb::VertexBase* primVtx, const LHCb::Particle* P, const std::string head,
                            Tuples::Tuple& tuple ) const;
  void       calculateDIRAError( const LHCb::Particle*, const LHCb::VertexBase*, double&, double& ) const;

  //  IContextTool* m_context;
  IDVAlgorithm* m_context;
};

#endif // FSOOMRO_TUPLETOOLDIRA_H
