
#pragma once

// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/Transform3DTypes.h" //Added by Arvind, to get physics vector transformations in Gaudi working
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "DecayTreeTupleBase/TupleToolBase.h"

#include "Kernel/IParticleTupleTool.h"
#include "Kernel/ISubstitutePID.h"
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IDVAlgorithm.h"
#include "Kernel/ParticleProperty.h"
#include "Kernel/IParticleDescendants.h"
#include "Kernel/Escape.h"
#include "Kernel/IParticlePropertySvc.h"

#include "Kernel/IParticle2MCAssociator.h"

#include "TrackInterfaces/ITrackStateProvider.h"

#include "DecayTreeFitter/Fitter.h"

#include "Event/RecVertex.h"
#include "Event/Particle.h"

#include "LoKi/ParticleProperties.h"

// boost
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

// STL
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

// CLHEP
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Vector.h>

struct IParticleDescendants;

namespace LHCb
{
  class ParticleID;
  class IParticlePropertySvc;
  class VertexBase;
  class Particle;
  class RecVertex ;
  class DecayTree ;
}

struct IDVAlgorithm;

class TupleToolDecayFit : public TupleToolBase, virtual public IParticleTupleTool
{
    private:
        typedef std::map<std::string,std::string> SubstitutionMap;
        typedef std::map< std::string, std::vector<double> > TupleMap;
    public:
        // Standard constructor
        TupleToolDecayFit( const std::string& type,
                           const std::string& name,
                           const IInterface* parent );

        // Destructor
        ~TupleToolDecayFit() = default;

        StatusCode initialize() override;

        StatusCode finalize() override;

        StatusCode fill( const LHCb::Particle*, const LHCb::Particle*,
                         const std::string&, Tuples::Tuple& ) override;


    private:
        std::vector<const LHCb::VertexBase*> originVertex( const LHCb::Particle*, const LHCb::Particle* ) const;

        /// Get the TES location for a data object
        template<class TYPE>
        inline std::string tesLocation( const TYPE * obj ) const
        {
            return ( obj && obj->parent() && obj->parent()->registry() ?
                    obj->parent()->registry()->identifier() : "NotInTES" );
        }

    private:
        IDVAlgorithm* m_dva = nullptr;

        ToolHandle<ITrackStateProvider> m_stateprovider{ "TrackStateProvider" };

        bool m_constrainToOriginVertex;

};
