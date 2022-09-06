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
#ifndef INTERNALRECOTRACK_H
#define INTERNALRECOTRACK_H

#include "Event/State.h"
#include "InternalParticle.h"

struct ITrackStateProvider;
namespace LHCb {
  class TrackTraj;
  class Track;
} // namespace LHCb

namespace DecayTreeFitter {

  class InternalRecoTrack : public InternalParticle {
  public:
    InternalRecoTrack( const LHCb::Particle& bc, const ParticleBase* mother, const Configuration& config );
    ~InternalRecoTrack();
    int         type() const override { return kInternalRecoTrack; }
    virtual int dimM() const { return 5; }

    virtual ErrCode projectRecoConstraint( const FitParams&, Projection& ) const;
    ErrCode         updCache( const FitParams& fitparams );

    void addToConstraintList( constraintlist& alist, int depth ) const override {
      InternalParticle::addToConstraintList( alist, depth );
      alist.push_back( Constraint( this, Constraint::track, depth, dimM() ) );
    }

    ErrCode projectConstraint( Constraint::Type, const FitParams&, Projection& ) const override;

  private:
    // ErrCode updFltToMother(const FitParams& fitparams) ;
    void               setFlightLength( double flt ) { m_flt = flt; }
    const LHCb::Track& track() const { return *m_track; }
    const LHCb::State& state() const { return m_state; }

  private:
    const LHCb::Track*         m_track;
    const ITrackStateProvider* m_stateprovider;
    bool                       m_useTrackTraj;
    const LHCb::TrackTraj*     m_tracktraj;
    bool                       m_cached;
    double                     m_flt;
    LHCb::State                m_state;
  };
} // namespace DecayTreeFitter

#endif
