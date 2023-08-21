
#include "TupleToolDecayFit2.h"
#include "TMath.h"
#include <algorithm>
#include "TMinuit.h"

using namespace LHCb;
using namespace std;

TupleToolDecayFit2::TupleToolDecayFit2( const std::string& type, const std::string& name, const IInterface* parent ):TupleToolBase( type, name, parent )
{
    declareProperty( "StateProvider", m_stateprovider ) ;
    declareProperty( "constrainToOriginVertex", m_constrainToOriginVertex = false,
                    "Do a refit constraining to Origin Vertex (could be PV)");

    declareInterface<IParticleTupleTool>(this);
}

StatusCode TupleToolDecayFit2::initialize()
{
    StatusCode sc = TupleToolBase::initialize();
    if( sc.isFailure() ) return sc;

    m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc(), this );
    if ( !m_dva ) return Error("Couldn't get parent IDVAlgorithm", StatusCode::FAILURE);

    return sc;
}

StatusCode TupleToolDecayFit2::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
    return StatusCode{ TupleToolBase::finalize() && sc };
}

StatusCode TupleToolDecayFit2::fill( const LHCb::Particle* mother, const LHCb::Particle* P,
                                    const std::string& head, Tuples::Tuple& tuple )
{
    LHCb::DecayTree tree(*P);
    LHCb::Particle* treeHead = tree.head();

    // 1) Initialise vector of measured parameters m in LHCb frame
    int dimM = 32;
    // Get primary vertex
    std::vector<const VertexBase*> originVtx;
    originVtx = originVertex(mother, P);
    if( originVtx.empty() ){return Error("Can't get an origin vertex");}
    ROOT::Math::XYZPoint PV( originVtx[0]->position().X(), originVtx[0]->position().Y(), originVtx[0]->position().Z() ); // originVtx[0] is the best PV in the event

    // Get tau+ decay vertex + tau+ pions momenta
    const LHCb::Particle* taup = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-15) );
    ROOT::Math::XYZPoint DV1( taup->endVertex()->position().X(), taup->endVertex()->position().Y(), taup->endVertex()->position().Z() );
    LHCb::Particle::ConstVector taup_pions = taup->daughtersVector();

    const LHCb::ProtoParticle* taup_pi1_proto = taup_pions[0]->proto();
    const LHCb::ProtoParticle* taup_pi2_proto = taup_pions[1]->proto();
    const LHCb::ProtoParticle* taup_pi3_proto = taup_pions[2]->proto();

    const LHCb::Track* track1 = taup_pi1_proto->track();
    const LHCb::Track* track2 = taup_pi2_proto->track();
    const LHCb::Track* track3 = taup_pi3_proto->track();

    LHCb::State state1;
    LHCb::State state2;
    LHCb::State state3;

    Double_t ztolerance = 0.001 * Gaudi::Units::mm;

    StatusCode sc1 = m_stateprovider->state( state1, *track1, DV1.z(), ztolerance );
    StatusCode sc2 = m_stateprovider->state( state2, *track2, DV1.z(), ztolerance );
    StatusCode sc3 = m_stateprovider->state( state3, *track3, DV1.z(), ztolerance );

    ROOT::Math::XYZVector P1 = state1.momentum();
    ROOT::Math::XYZVector P2 = state2.momentum();
    ROOT::Math::XYZVector P3 = state3.momentum();

    Gaudi::LorentzVector mom1( taup_pions[0]->momentum().X(), taup_pions[0]->momentum().Y(), taup_pions[0]->momentum().Z(), taup_pions[0]->momentum().E() );

    std::cout << "p1x state = " << P1.x() << endl;
    std::cout << "p1y state = " << P1.y() << endl;
    std::cout << "p1z state = " << P1.z() << endl;

    std::cout << "p1x part = " << mom1.x() << endl;
    std::cout << "p1y part = " << mom1.y() << endl;
    std::cout << "p1z part = " << mom1.z() << endl;

    Gaudi::SymMatrix6x6 taup_pi1_cov = state1.posMomCovariance();
    Gaudi::SymMatrix6x6 taup_pi2_cov = state2.posMomCovariance();
    Gaudi::SymMatrix6x6 taup_pi3_cov = state3.posMomCovariance();

    // Get tau- decay vertex + tau- pions momenta
    // const LHCb::Particle* taum = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(15) );
    // ROOT::Math::XYZPoint DV2( taum->endVertex()->position().X(), taum->endVertex()->position().Y(), taum->endVertex()->position().Z() );
    // LHCb::Particle::ConstVector taup_pions = taup->daughtersVector();

    // const LHCb::ProtoParticle* taup_pi1_proto = taup_pions[0]->proto();
    // const LHCb::ProtoParticle* taup_pi2_proto = taup_pions[1]->proto();
    // const LHCb::ProtoParticle* taup_pi3_proto = taup_pions[2]->proto();

    // const LHCb::Track* track1 = taup_pi1_proto->track();
    // const LHCb::Track* track2 = taup_pi2_proto->track();
    // const LHCb::Track* track3 = taup_pi3_proto->track();

    // LHCb::State state1;
    // LHCb::State state2;
    // LHCb::State state3;

    // Double_t ztolerance = 0.001 * Gaudi::Units::mm;

    // StatusCode sc1 = m_stateprovider->state( state1, *track1, DV1.z(), ztolerance );
    // StatusCode sc2 = m_stateprovider->state( state2, *track2, DV1.z(), ztolerance );
    // StatusCode sc3 = m_stateprovider->state( state3, *track3, DV1.z(), ztolerance );

    // ROOT::Math::XYZVector P1 = state1.momentum();
    // ROOT::Math::XYZVector P2 = state2.momentum();
    // ROOT::Math::XYZVector P3 = state3.momentum();

    // Gaudi::SymMatrix6x6 taup_pi1_cov = state1.posMomCovariance();
    // Gaudi::SymMatrix6x6 taup_pi2_cov = state2.posMomCovariance();
    // Gaudi::SymMatrix6x6 taup_pi3_cov = state3.posMomCovariance();


    return StatusCode(true);
}

std::vector<const VertexBase*> TupleToolDecayFit2::originVertex( const Particle* mother, const Particle* P ) const
{
    std::vector<const VertexBase*> oriVx;

    if ( mother == P )
    {// the origin vertex is the primary.
        const auto bpv = m_dva->bestVertex( P );
        if ( bpv )
        {
            oriVx.push_back(bpv);
            if ( msgLevel(MSG::VERBOSE) )
                verbose() << "Pushed back bpv " << bpv << " from "
                        << tesLocation(bpv) << " at "
                        << bpv->position() << endmsg ;
        }
        else if ( m_constrainToOriginVertex)
        {
            Warning( "NULL bestPV while constraining to origin vertex. Fit will be ignored.",
                    StatusCode::SUCCESS, 0 ).ignore();
        }
    }
    else
    {
        const auto & dau = mother->daughters ();
        if ( dau.empty() ) return oriVx ;

        for ( const auto& d : dau )
        {
            if ( P == d )
            {
                oriVx.push_back(mother->endVertex());
                return oriVx ;
            }
        }

        // vertex not yet found, get deeper in the decay:
        for ( const auto& d : dau )
        {
            if ( P != d && !d->isBasicParticle() )
            {
                oriVx = originVertex( d, P );
                if( !oriVx.empty() )
                {
                    return oriVx ;  // found
                }
            }
        }
    }
    return oriVx;
}

// Declaration of the Tool Factory
DECLARE_COMPONENT( TupleToolDecayFit2 )