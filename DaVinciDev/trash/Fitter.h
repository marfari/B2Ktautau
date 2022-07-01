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
// $Id: Fitter.h,v 1.7 2010-05-29 13:39:56 ibelyaev Exp $
// ============================================================================
#ifndef DECAYTREEFITTER_FITTER_HH
#define DECAYTREEFITTER_FITTER_HH
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
// ============================================================================
// GaudiKernel
// ============================================================================
#include "GaudiKernel/SmartIF.h"
// ============================================================================
// Track interfaces
// ============================================================================
#include "TrackInterfaces/ITrackStateProvider.h"
// ============================================================================
// LHCbMath
// ============================================================================
#include "LHCbMath/ParticleParams.h"
#include "LHCbMath/ValueWithError.h"
// ============================================================================
// DaVinciTypes
// ============================================================================
#include "Kernel/DecayTree.h"
// ============================================================================
namespace DecayTreeFitter {
  class DecayChain;
  class FitParams;
  class ParticleBase;
  class ChiSquare;
} // namespace DecayTreeFitter
class IAlgorithm;
class IAlgContextSvc;
// ============================================================================
namespace DecayTreeFitter {
  // ==========================================================================
  /** @class Fitter DecayTreeFitter/Fitter.h
   *  ``Decay-Tree-Fitter''`
   *  @author Wouter Hulsbergen  Wouter.Hulsbergen@nikhef.nl
   */
  class GAUDI_API Fitter {
  public:
    // ========================================================================
    enum FitStatus { UnFitted = -1, Success = 0, Failed, BadInput, NonConverged };
    // ========================================================================
  public:
    // ========================================================================
    /// constructor from the particle (decay head)
    Fitter( const LHCb::Particle& bc, const bool forceFitAll = true, const ITrackStateProvider* extrapolator = 0, const bool useTrackTraj = true );
    /// constructor from the particle (decay head)
    Fitter( const LHCb::Particle& bc, const ITrackStateProvider* extrapolator, const bool forceFitAll = true, const bool useTrackTraj = true );
    /// constructor from the particle (decay head) and primary vertex
    Fitter( const LHCb::Particle& bc, const LHCb::VertexBase& pv, const bool forceFitAll = true,
            const ITrackStateProvider* extrapolator = 0, const bool useTrackTraj = true );
    /// constructor from the particle (decay head) and primary vertex
    Fitter( const LHCb::Particle& bc, const LHCb::VertexBase& pv, const ITrackStateProvider* extrapolator,
            const bool forceFitAll = true, const bool useTrackTraj = true );
    /// destructor
    ~Fitter(); // destructor
    // ========================================================================
  public:
    // ========================================================================
    /// Add or remove a mass constraint
    void setMassConstraint( const LHCb::Particle& cand, bool add = true );
    /// Add a constraint to a mass different from the property table mass
    void setMassConstraint( const LHCb::Particle& cand, double mass );
    /// Add or remove a mass constraintfor a certain ParticleID
    void setMassConstraint( const LHCb::ParticleID& pid, bool add = true );
    /// Add a constraint to a mass different from the property table mass
    void setMassConstraint( const LHCb::ParticleID& pid, double mass );
    /// Fit the decay tree
    void fit( int maxNumberOfIterations = 10, double relDeltaChisquareConverged = 0.01, int maxndiverging = 3, double dChisqQuit = 1.e6);
    void fit( int maxNumberOfIterations, double relDeltaChisquareConverged, int maxndiverging, double dChisqQuit, 
              std::vector<float>& chisq_iters, std::vector<Gaudi::XYZPoint>& BV, std::vector<Gaudi::XYZPoint>& DV1, std::vector<Gaudi::XYZPoint>& DV2,
              const LHCb::Particle* P, std::vector<float> B_M);
 
    /// Fit just one step
    void fitOneStep();
    /// Print the result of the fit
    void print() const { fillStream( std::cout ); }
    /// Print the result of the fit
    std::ostream& fillStream( std::ostream& s ) const;
    /// The top level particle that is fitted
    const LHCb::Particle* particle() const { return m_particle; }
    /** Currently the only accessor to the actual fitted data
     *  @param p (INPUT) the particle
     *  @retrn the fitted parameters ( 0 for invaild parameters/fits)
     */
    const Gaudi::Math::ParticleParams* fitParams( const LHCb::Particle* p = 0 ) const;
    /// Total chisquare
    double chiSquare() const { return m_chiSquare; }
    /// Total number of DOFs
    int nDof() const;
    /// Status of fit
    FitStatus status() const { return m_status; }
    /// Number of iterations used by vertex fit
    int nIter() const { return m_niter; }
    /// get the chisquare of everything 'downstream' of a particle
    ChiSquare chiSquare( const LHCb::Particle& p ) const;
    /**  Compute the decay length sum of two particles in
     *   the decay tree (useful for e.g. B->DD)
     */
    Gaudi::Math::ValueWithError decayLengthSum( const LHCb::Particle&, const LHCb::Particle& ) const;
    /** return an updated decay tree.
     *  this is not a final solution. will
     *  try to move more info to Particle
     */
    LHCb::DecayTree getFittedTree() const;
    /** methods to retrieve the result in terms of LHCb::Particles
     * (note: mother vertex is not updated, and decay length cannot be
     * stored anywhere. Use fitParams instead
     */
    LHCb::Particle getFitted() const;
    /** methods to retrieve the result in terms of LHCb::Particles
     * (note: mother vertex is not updated, and decay length cannot be
     * stored anywhere. Use fitParams instead
     */
    LHCb::Particle getFitted( const LHCb::Particle& cand ) const;
    /// update a particlular candidate in the tree
    bool updateCand( LHCb::Particle& cand ) const;
    /// update a particlular candidate in the tree
    bool updateTree( LHCb::Particle& cand ) const;
    /// error code
    int errCode() { return m_errCode; }
    /// set the verbosity level (for debugging only)
    static void setVerbose( int i );

    // ========================================================================
  public:
    // ========================================================================
    /// get the extrapolator
    const ITrackStateProvider* extrapolator() const { return m_extrapolator; }
    /// set the track extrapolator
    void setStateProvider( const ITrackStateProvider* extrapolator );
    // ========================================================================
  protected:
    // ========================================================================
    // expert interface. not yet for real consumption
    Gaudi::Math::ParticleParams fitParams( const ParticleBase& pb ) const;
    /// Name of a particle in the decay tree
    std::string name( const LHCb::Particle& cand ) const;
    // ========================================================================
    Gaudi::Math::ValueWithError decayLengthSum( const ParticleBase&, const ParticleBase& ) const;
    // ========================================================================
    DecayChain*       decaychain() { return m_decaychain; }
    FitParams*        fitparams() { return m_fitparams; }
    const DecayChain* decaychain() const { return m_decaychain; }
    const FitParams*  fitparams() const { return m_fitparams; }
    // ========================================================================
    double globalChiSquare() const;
    // ========================================================================
    // must be moved to derived class or so ...
    double add( const LHCb::Particle& cand );
    double remove( const LHCb::Particle& cand );
    void   updateIndex();
    // ========================================================================
    LHCb::Particle* fittedCand( const LHCb::Particle& cand, LHCb::Particle* headoftree ) const;
    void            updateCand( const ParticleBase& pb, LHCb::Particle& cand ) const;
    // ========================================================================
  private:
    // ========================================================================
    /// default constructor is disabled
    Fitter(); //  default constructor is disabled
    /// copy constructor is disabled
    Fitter( const Fitter& ); //     copy constructor is disabled
    /// assignement operator is disabled
    Fitter& operator=( const Fitter& ); // assignement operator is disabled
    // ========================================================================
  private:
    // ========================================================================
    /// Get context service
    const IAlgContextSvc* algContextSvc() const;
    /// Get message service
    IMessageSvc* msgService() const;
    /// Print a message
    void print( const std::string& msg, const MSG::Level level = MSG::INFO ) const;
    /// Get current active algorithm
    const IAlgorithm* getAlg() const;
    // ========================================================================
  private:
    // ========================================================================
    const LHCb::Particle* m_particle;
    DecayChain*           m_decaychain;
    FitParams*            m_fitparams;
    FitStatus             m_status;
    double                m_chiSquare;
    int                   m_niter;
    int                   m_errCode;
    // ========================================================================
    typedef std::map<const LHCb::Particle*, Gaudi::Math::ParticleParams> Map;
    mutable Map                                                          m_map;
    // ========================================================================
    /// track extrapolator (if needed)
    SmartIF<ITrackStateProvider> m_extrapolator; // track extrapolator
    // ========================================================================
  };
  // ==========================================================================
} // namespace DecayTreeFitter
// ============================================================================
// The END
// ============================================================================
#endif
