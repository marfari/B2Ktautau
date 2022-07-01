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
#ifndef RECOTRACK_H
#define RECOTRACK_H

#include "Configuration.h"
#include "Event/State.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "RecoParticle.h"

struct ITrackStateProvider;
namespace LHCb {
  class TrackTraj;
  class Track;
} // namespace LHCb

namespace DecayTreeFitter {

  class RecoTrack : public RecoParticle {
  public:
    RecoTrack( const LHCb::Particle& bc, const ParticleBase* mother, const Configuration& config );
    virtual ~RecoTrack();

    ErrCode initPar2( FitParams* ) override;
    ErrCode initCov( FitParams* ) const override;
    int     dimM() const override { return 5; }
    int     type() const override { return kRecoTrack; }

    ErrCode projectRecoConstraint( const FitParams&, Projection& ) const override;
    ErrCode updCache( const FitParams& fitparams );
    // tatic void setApplyCovCorrection(bool b=true) { gApplyCovCorrection = b ; }
    // static void correctCov(HepSymMatrix& V) ;

    int nFinalChargedCandidates() const override { return 1; }

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      alist.push_back( Constraint( this, Constraint::track, depth, dimM() ) );
    }
    // ErrCode updFltToMother(const FitParams& fitparams) ;
    void               setFlightLength( double flt ) { m_flt = flt; }
    const LHCb::Track& track() const { return *m_track; }
    const LHCb::State& state() const { return m_state; }

    // return a trajectory (declared by base class)
    const LHCb::Trajectory* trajectory() const override;

    // return a tracktraj
    const LHCb::TrackTraj* tracktraj() const;

  private:
    const Gaudi::XYZVector     m_bfield;
    const LHCb::Track*         m_track;
    const ITrackStateProvider* m_stateprovider;
    bool                       m_useTrackTraj;
    const LHCb::TrackTraj*     m_tracktraj;
    bool                       m_ownstracktraj;
    bool                       m_cached;
    double                     m_flt;
    LHCb::State                m_state;
    double                     m_bremEnergy;
    double                     m_bremEnergyCov;
    double                     m_maxCovVeloTrack =
        75000. * 75000; // Maximum value for the initial covariance matrix' elements, only for VELO tracks
  };

} // namespace DecayTreeFitter
#endif
