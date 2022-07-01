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
#ifndef TUPLETOOLBS2MUMUISO2_H
#define TUPLETOOLBS2MUMUISO2_H 1

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Event/Particle.h"
#include "Kernel/DaVinciAlgorithm.h"
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IDistanceCalculator.h"
#include "Kernel/IParticleDictTool.h"
#include "Kernel/IParticleTupleTool.h"
#include "Kernel/IRelatedPVFinder.h"
#include "LoKi/AParticleCuts.h"
#include "LoKi/ParticleContextCuts.h"
#include "LoKi/ParticleCuts.h"
#include "TMVA/Reader.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVector3.h"
#include "TVectorD.h"
#include <algorithm>
#include <iostream>

using namespace std;
using namespace LHCb;

/** @class RelInfoTrackIsolationBDT2 RelInfoTrackIsolationBDT2.h
 *
 *  Converted from TupleToolMuonVariables
 *
 *  @author Maik BECKER (maik.becker@cern.ch), Matteo RAMA (matteo.rama@cern.ch)
 *  @date   2019-01-23
 *  2020-07-12 : converted to TupleTool format for usage on DSTs by Biplab Dey
 */

struct IDVAlgorithm;
struct IDistanceCalculator;
class IParticleDictTool;

class TupleToolBs2MuMuIso2 : public TupleToolBase, virtual public IParticleTupleTool {

public:
  // Constructor
  TupleToolBs2MuMuIso2( const string& type, const string& name, const IInterface* parent );

  // Initialization
  StatusCode initialize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  // Declare RelatedInfoTool properties
  vector<int>    m_particles;
  vector<string> m_variables_nominal;
  vector<string> m_variables;
  string         m_PVLocation;
  string         m_particleLocation;
  string         m_trackLocation;
  string         m_BDTLocation;
  string         m_isoBDT_long_xmlFile;
  string         m_isoBDT_velo_xmlFile;
  string         m_newBDT_xmlFile;

  IDVAlgorithm*              m_dva;
  const IDistanceCalculator* m_dist;

  TMVA::Reader*    m_isoBDT_long_reader;
  TMVA::Reader*    m_isoBDT_velo_reader;
  static const int m_size_isoBDT_long_inputVars = 9;
  static const int m_size_isoBDT_velo_inputVars = 6;
  float            m_isoBDT_long_inputValues[m_size_isoBDT_long_inputVars];
  float            m_isoBDT_velo_inputValues[m_size_isoBDT_velo_inputVars];
  string           m_isoBDT_method;

  TMVA::Reader*    m_newBDT_reader;
  static const int m_size_newBDT_inputVars       = 7;
  static const int m_size_newBDT_inputVars_spect = 1;
  float            m_newBDT_inputValues[m_size_newBDT_inputVars];
  float            m_newBDT_inputValues_spect[m_size_newBDT_inputVars_spect];
  string           m_newBDT_method;

  Particle::ConstVector m_daughters;
  Track::ConstVector    m_tracksToIgnore;

  // Declare RelatedInfo variables
  double m_isoBDT_long;
  double m_isoBDT_long_d1max1;
  double m_isoBDT_long_d1max2;
  double m_isoBDT_long_d1max3;
  double m_isoBDT_long_d2max1;
  double m_isoBDT_long_d2max2;
  double m_isoBDT_long_d2max3;
  double m_isoBDT_velo;
  double m_isoBDT_velo_d1max1;
  double m_isoBDT_velo_d1max2;
  double m_isoBDT_velo_d1max3;
  double m_isoBDT_velo_d2max1;
  double m_isoBDT_velo_d2max2;
  double m_isoBDT_velo_d2max3;
  double m_newBDT;

  double map_name_var( string isolvarname ); // maps the isol name to the variable value

  //===========================================================================
  // Track Isolation Methods
  //===========================================================================
  StatusCode calcTrackIsoBDTLong( const Particle* part );

  StatusCode calcTrackIsoBDTVelo( const Particle* part );

  StatusCode calcBDT( const Particle* part );

  //===========================================================================
  // Helper Functions
  //===========================================================================
  StatusCode getAllDaughters( const Particle* part, Particle::ConstVector& daughters,
                              Track::ConstVector& tracksToIgnore );

  StatusCode calcIPToAnyPV( const Particle* part, double& ips );

  //===========================================================================
  // Calculates the perpendicular feet on the track's and lepton's trajectory
  //
  // The formulae are derived from:
  // Let p1 and p2 be the points where the intersect (p1 - p2) is equal to the
  // doca-vector of the track and the lepton satisfying
  //    ( p1 - p2 ) * track_mom == 0
  // && ( p1 - p2 ) * lepton_mom == 0
  // Solve equations for scalars lepton1 and lepton2 which are the translations
  // along the unit momentum vectors
  // => p1 = track_pos + lepton1 * Unit(track_mom) and
  //    p2 = lepton_pos + lepton2 * Unit(lepton_mom)
  //
  // Additionally, a "vertex" is calculated as the mid-point
  // on the doca line segment
  //===========================================================================
  StatusCode getPerpFeet( Gaudi::XYZVector daughterMom, Gaudi::XYZPoint daughterPos, Gaudi::XYZVector trackMom,
                          Gaudi::XYZPoint trackPos, Gaudi::XYZPoint& vert );

  double calcVertexDist( Gaudi::XYZPoint track, const VertexBase* vert );

  StatusCode inCone( Gaudi::XYZPoint o1, Gaudi::XYZVector p1, Gaudi::XYZPoint o2, Gaudi::XYZVector p2, double& doca );

  StatusCode closestPoint( Gaudi::XYZPoint o1, Gaudi::XYZVector p1, Gaudi::XYZPoint o2, Gaudi::XYZVector p2,
                           Gaudi::XYZPoint& close1, Gaudi::XYZPoint& close2, Gaudi::XYZPoint& vert );

  double calcEnclosedAngle( Gaudi::XYZVector p1, Gaudi::XYZVector p2 );

  //===========================================================================
  // Calculates fc value which is a combination of momenta, transverse
  // momenta and angles
  //
  // For more information see Bs2mumu roadmap or
  // G. Mancinelli & J. Serrano LHCb-INT-2010-011
  // XYZVector::R() equals length of the vector,
  // XYZVector::Rho() equals length projection on X-Y plane
  // => like transverse momentum
  //===========================================================================
  StatusCode calcFC( Gaudi::XYZVector daughterMom, Gaudi::XYZVector trackMom, Gaudi::XYZPoint vertDaughterTrack,
                     const VertexBase* PV, double& fc );

  void computeIso( vector<double> iso_vec, double& iso_max1, double& iso_max2, double& iso_max3 );

  StatusCode getTrackMinIP( const Track* track, double& ips );

  void get_fc_pvdis_svdis( const Particle* B, const Track* mutrack, const Track* track, double& fc, double& pvdis,
                           double& svdis );

  const VertexBase* get_PV_BminIPS( const Particle* B );

  void get_BIPwrtPV( const Particle* B, const VertexBase* PV, double& ip, double& ips, double& fl, double& fl0,
                     double& ctau );

  TMatrixD get_dFdq( TVectorD q );

  TMatrixDSym get_BCov( const Particle* B );

  void getTrackIPSwrtPV( Gaudi::XYZVector& mom_trk, Gaudi::XYZPoint& pos_trk, const Gaudi::SymMatrix6x6& cov_trk,
                         Gaudi::SymMatrix3x3 slope_cov_trk, Gaudi::XYZPoint& pos_pv, const Gaudi::SymMatrix3x3& cov_pv,
                         double& chi2 );

  TMatrixDSym get_trkCov( Gaudi::XYZPoint& pos_trk, Gaudi::XYZVector& mom_trk, const Gaudi::SymMatrix6x6& cov_trk,
                          const Gaudi::SymMatrix3x3& slope_cov_trk );

  TMatrixD get_dFdp( TVector3 v0, TVector3 x0 );

  TMatrixD get_dFdv( TVectorD p );
};
#endif // TUPLETOOLBS2MUMUISO2_H
