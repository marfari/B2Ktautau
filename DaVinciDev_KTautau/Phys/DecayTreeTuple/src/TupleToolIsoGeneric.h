/*****************************************************************************\
* (c) Copyright 2020 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// $Id: TupleToolIsoGeneric.h
#ifndef TUPLETOOLISOGENERIC_H
#define TUPLETOOLISOGENERIC_H 1

// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"

#include "Kernel/IParticleDescendants.h"
#include "Kernel/IVertexFit.h"
#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDVAlgorithm.h>

#include "TMVA/Reader.h"

struct IDVAlgorithm;
struct IDistanceCalculator;
struct IPVReFitter;
struct IVertexFit;
// from local

/**
 *
 *
 *  Tuple designed for isolation study
 *  @ Basem KHANJI
 *  @ date   01/09/2020
 */

class TupleToolIsoGeneric : public TupleToolBase, virtual public IParticleTupleTool {

public:
  /// Standard constructor
  TupleToolIsoGeneric( const std::string& type, const std::string& name, const IInterface* parent );

  // virtual StatusCode initialize();
  StatusCode initialize() override;

  // virtual ~TupleToolIsoGeneric( ){}; ///< Destructor not needed
  ~TupleToolIsoGeneric() {}

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

  double calculateThetaL( const LHCb::Particle* particle );
  double calculateThetaK( const LHCb::Particle* particle );
  double calculatePhi( const LHCb::Particle* particle );
  double calculateTransThetaTr( const LHCb::Particle* particle );
  double calculateTransPhiTr( const LHCb::Particle* particle );
  double calculateTransThetaV( const LHCb::Particle* particle );

private:
  Float_t track_minIPchi2, track_eta, track_phi, track_pt, track_Matchchi2, PVdis_SignalTrack_track,
      SVdis_SignalTrack_track, doca_SignalTrack, angle_SignalTrack, fc_mu;
  std::string                  aux_name;
  IDVAlgorithm*                m_dva            = nullptr;
  const IDistanceCalculator*   m_dist           = nullptr;
  IParticleDescendants*        m_descend        = nullptr;
  std::string                  m_pvReFitterName = "LoKi::PVReFitter:PUBLIC";
  IPVReFitter*                 m_pvReFitter     = nullptr;
  IVertexFit*                  m_vtxfitter      = nullptr;
  Gaudi::Property<std::string> m_ParticlePath{this, "ParticlePath", "/Event/Phys/StdAllNoPIDsPions/Particles",
                                              "Underlying particles"};
  std::vector<std::string>     toolNames_default;
  Gaudi::Property<std::vector<std::string>> m_toolNames{this, "ToolNames", toolNames_default,
                                                        "Tools to apply on the none isolated particle"};
  Gaudi::Property<std::string>              m_weightsName{this, "WeightsFile", "IsolationTools_KMuNu.xml",
                                             "Name of the weight file for Isotool"};
  Gaudi::Property<std::string> m_BDTLocation{this, "BDTLocation", System::getEnv( "TMVAWEIGHTSROOT" ) + "/data/",
                                             "getting weight file for Isotool"};
  Gaudi::Property<bool>        m_verbose{this, "VerboseMode", false, "Fill verbose Isolation information"};

  std::vector<IParticleTupleTool*> m_tuple;
  TMVA::Reader*                    m_Reader = nullptr;
  double                           get_MINIPCHI2( const LHCb::Particle* p );
  bool                             isTrackInDecay( const LHCb::Track* track );

  void InCone( Gaudi::XYZPoint o1, Gaudi::XYZVector p1, Gaudi::XYZPoint o2, Gaudi::XYZVector p2, Gaudi::XYZPoint& vtx,
               double& doca, double& angle ) const;

  double pointer( Gaudi::XYZVector vtx, Gaudi::XYZVector P_tr, Gaudi::XYZVector P_mu ) const;

  double arcosine( Gaudi::XYZVector p1, Gaudi::XYZVector p2 ) const;

  void closest_point( Gaudi::XYZPoint o, Gaudi::XYZVector p, Gaudi::XYZPoint o_mu, Gaudi::XYZVector p_mu,
                      Gaudi::XYZPoint& close1, Gaudi::XYZPoint& close2, Gaudi::XYZPoint& vertex, bool& fail ) const;
};
//===========================================================================//
#endif // USER_ISOBDT_SMALL_H
