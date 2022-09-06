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
#ifndef DECAYTREEFITTER_CONFIG_H
#define DECAYTREEFITTER_CONFIG_H

struct ITrackStateProvider;

namespace DecayTreeFitter {
  class Configuration {
  public:
    Configuration( bool forceFitAll = true, const ITrackStateProvider* stateprovider = 0, bool useTrackTraj = true )
        : m_forceFitAll( forceFitAll ), m_stateprovider( stateprovider ), m_useTrackTraj( useTrackTraj ) {}

    const ITrackStateProvider* stateProvider() const { return m_stateprovider; }

    bool forceFitAll() const { return m_forceFitAll; }
    bool useTrackTraj() const { return m_stateprovider && m_useTrackTraj; }

  private:
    bool                       m_forceFitAll;
    const ITrackStateProvider* m_stateprovider;
    bool                       m_useTrackTraj;
  };
} // namespace DecayTreeFitter

#endif
