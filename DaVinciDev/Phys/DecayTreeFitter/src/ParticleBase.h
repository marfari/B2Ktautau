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
#ifndef __VTX_PARTICLEBASE_HH__
#define __VTX_PARTICLEBASE_HH__

#include "Constraint.h"
#include "DecayTreeFitter/ChiSquare.h"
#include "Projection.h"
#include <string>
#include <vector>
#include "Event/Particle.h"

namespace LHCb {
  class Particle;
  class ParticleID;
  class ParticleProperty;
  class Trajectory;
} // namespace LHCb

namespace DecayTreeFitter {
  class FitParams;
  class Configuration;

  class ParticleBase {
  public:
    enum ParticleType {
      kInteractionPoint,
      kRecoComposite,
      kRecoResonance,
      kInternalParticle,
      kRecoTrack,
      kResonance,
      kRecoPhoton,
      kMissingParticle,
      kJetMomentum,
      kInternalRecoTrack,
      kRecoMergedPi0,
      kConvertedPhoton
    };
    typedef std::vector<ParticleBase*> ParticleContainer;

    // 'default' constructor
    ParticleBase( const LHCb::Particle& bc, const ParticleBase* mother );

    // constructor used for InteractionPoint
    ParticleBase( const std::string& name );

    virtual ~ParticleBase();

    static ParticleBase* createParticle( const LHCb::Particle& bc, const ParticleBase* mother,
                                         const Configuration& config );

    virtual int           dim() const = 0;
    virtual void          updateIndex( int& offset );
    virtual ErrCode       initPar1( FitParams* ) = 0; // init everything that does not need mother vtx
    virtual ErrCode       initPar2( FitParams* ) = 0; // everything else
    virtual ErrCode       initCov( FitParams* ) const;
    virtual std::string   parname( int index ) const;
    virtual std::ostream& fillStream( std::ostream&, const FitParams* ) const;
    void                  print( const FitParams* fitpar ) const { fillStream( std::cout, fitpar ); }

    const ParticleBase*   locate( const LHCb::Particle& bc ) const;
    void                  locate( const LHCb::ParticleID& pid, ParticleContainer& result );
    const LHCb::Particle& particle() const { return *m_particle; }

    int                 index() const { return m_index; }
    const ParticleBase* mother() const { return m_mother; }
    const std::string&  name() const { return m_name; }

    virtual ErrCode projectGeoConstraint( const FitParams&, Projection& ) const;
    virtual ErrCode projectMassConstraint( const FitParams&, Projection& ) const;
    virtual ErrCode projectConstraint( Constraint::Type, const FitParams&, Projection& ) const;

    // indices to fit parameters
    virtual int type() const = 0;
    virtual int posIndex() const { return -1; }
    virtual int lenIndex() const { return -1; }
    virtual int momIndex() const { return -1; }

    // does the particle have a 3-momentum or a 4-momentum ?
    // virtual bool hasEnergy() const { return true; } // default

    virtual bool hasEnergy() const  
    { 
      if( (particle().particleID().pid() == 15) || (particle().particleID().pid() == -15) )
      {
        return false;
      }
      else 
      {
        return true;
      }
    }

    // does the particle have is own decay vertex ? (resonances and
    // recoparticles do not)
    virtual bool hasPosition() const { return false; }

    int eneIndex() const { return hasEnergy() ? momIndex() + 3 : -1; }

    // calculates the global chisquare (pretty useless)
    virtual double chiSquare( const FitParams* ) const;

    // access to particle PDT parameters
    double pdtMass() const { return m_pdtMass; }
    double pdtWidth() const { return m_pdtWidth; }
    double pdtCLifeTime() const { return m_pdtCLifeTime; }
    double pdtTau() const { return m_pdtMass > 0 ? m_pdtCLifeTime / m_pdtMass : 0; }
    int    charge() const { return m_charge; }

    // return a trajectory
    virtual const LHCb::Trajectory* trajectory() const { return 0; }

    // access to daughters
    typedef std::vector<ParticleBase*>   daucontainer;
    typedef daucontainer::const_iterator const_iterator;

    const daucontainer& daughters() const { return m_daughters; }
    const_iterator      begin() const { return m_daughters.begin(); }
    const_iterator      end() const { return m_daughters.end(); }
    ParticleBase*       addDaughter( const LHCb::Particle&, const Configuration& config );
    void                removeDaughter( const ParticleBase* pb );

    typedef std::vector<std::pair<const ParticleBase*, int>> indexmap;
    virtual void                                             retrieveIndexMap( indexmap& anindexmap ) const;

    void setMother( const ParticleBase* m ) { m_mother = m; }

    typedef std::vector<DecayTreeFitter::Constraint> constraintlist;
    virtual void                                     addToConstraintList( constraintlist& alist, int depth ) const = 0;
    virtual int                                      nFinalChargedCandidates() const;
    void                                             setParticle( const LHCb::Particle* bc ) { m_particle = bc; }

    // collect all particles emitted from vertex with position posindex
    void collectVertexDaughters( daucontainer& particles, int posindex );
    // set the mass constraint for this particle. return true if value changed
    bool setMassConstraint( bool add ) {
      std::swap( add, m_hasMassConstraint );
      return add != m_hasMassConstraint;
    }
    // set the mass of the mass constraint (use with care!)
    void setMassConstraint( double mass ) {
      m_hasMassConstraint = true;
      m_pdtMass           = mass;
    }

    ChiSquare chiSquare( const FitParams& params ) const;

    bool hasMassConstraint() const { return m_hasMassConstraint; }

  protected:
    static double pdtCLifeTime( const LHCb::Particle& bc );
    static bool   isAResonance( const LHCb::ParticleProperty& bc );
    static double bFieldOverC() { return 0; } // Bz/c
    ErrCode       initTau( FitParams* par ) const;
    void          makeName( const LHCb::Particle& bc );
    daucontainer& daughters() { return m_daughters; }

  protected:
    void setIndex( int i ) { m_index = i; }
    void setName( const std::string& n ) { m_name = n; }

  private:
    const LHCb::Particle*         m_particle;
    const ParticleBase*           m_mother;
    ParticleContainer             m_daughters;
    const LHCb::ParticleProperty* m_prop;
    int                           m_index;
    double                        m_pdtMass;      // cached mass
    double                        m_pdtWidth;     // particle width (for mass constraints)
    double                        m_pdtCLifeTime; // cached lifetime
    int                           m_charge;       // charge
    std::string                   m_name;
    bool                          m_hasMassConstraint;
  };

} // namespace DecayTreeFitter

#endif
