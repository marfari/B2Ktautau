/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#ifndef TUPLETOOLVTXDAUGHTERS_H
#define TUPLETOOLVTXDAUGHTERS_H 1

// Include files
// from Gaudi
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/IParticleTupleTool.h" // Interface

struct IVertexFit;
struct IDistanceCalculator;
struct IDVAlgorithm;

namespace LHCb {
  class Particle;
  class VertexBase;
  class Vertex;
} // namespace LHCb

namespace DecayTreeFitter {
  class Fitter;
}

/** @class TupleToolVtxDaughters TupleToolVtxDaughters.h
 *
 * \brief Compute vertices of pairs and triplet of daughters + FD info between these and the mother vertex for
DecayTreeTuple.
 *

 *
 * \sa DecayTreeTuple
 *
 *  @author Adlene Hicheur
 *  @date   2014-10-07, first hand
//  2019-12: upgrade to Triplet vertices, time to commit
 */
class TupleToolVtxDaughters : public TupleToolBase, virtual public IParticleTupleTool {

public:
  /// Standard constructor
  TupleToolVtxDaughters( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~TupleToolVtxDaughters(){}; ///< Destructor

  StatusCode initialize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  /// prong vertex
  const LHCb::VertexBase* ProngVertex( LHCb::Particle::ConstVector daus, double& mass, LHCb::Particle& PairMom ) const;

  /// fill end vertex stuff
  StatusCode fillVertex( const LHCb::VertexBase* vtx, const std::string& vtx_name, Tuples::Tuple& ) const;

  /// fill flight
  StatusCode fillFlight( const LHCb::VertexBase* pairVtx, const LHCb::Particle* P, const std::string& prefix,
                         Tuples::Tuple& tuple, const std::string& trail = "" ) const;

  /*   StatusCode fillFlight(const Gaudi::XYZPoint& pairVtx,
                        const LHCb::Particle* P,
                        const std::string& prefix,
                        Tuples::Tuple& tuple,
                        const std::string& trail = "" ) const ;*/

  double dira( const LHCb::VertexBase* oriVtx, const LHCb::Particle P ) const {
    if ( !oriVtx ) {
      Exception( "Wrong use of dira" );
      return -1501.;
    }
    const LHCb::Vertex* evtx = P.endVertex();
    if ( !evtx ) {
      Warning( "Cannot find end vertex", StatusCode::SUCCESS, 1 ).ignore();
      return -999.;
    }
    const Gaudi::XYZVector& A = P.momentum().Vect();
    const Gaudi::XYZVector  B = evtx->position() - oriVtx->position();
    return A.Dot( B ) / std::sqrt( A.Mag2() * B.Mag2() );
  }

private:
  const IDistanceCalculator* m_dist = nullptr;
  const IVertexFit*          m_vtx  = nullptr;

  bool m_fillSubVtx;
  bool m_doThreeProng;

  // bool m_fillMother;

  IDVAlgorithm* m_dva = nullptr;
};
#endif // TUPLETOOLVTXDAUGHTERS_H
