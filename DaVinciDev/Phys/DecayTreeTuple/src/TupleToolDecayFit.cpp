
#include "TupleToolDecayFit.h"
#include "TMath.h"
#include <algorithm>

#include "Kernel/IVertexFit.h"
#include "TrackInterfaces/ITrackFitter.h"

using namespace LHCb;
using namespace std;

TupleToolDecayFit::TupleToolDecayFit( const std::string& type, const std::string& name, const IInterface* parent ):TupleToolBase( type, name, parent )
{
    declareProperty( "StateProvider", m_stateprovider ) ;
    declareProperty( "constrainToOriginVertex", m_constrainToOriginVertex = false,
                    "Do a refit constraining to Origin Vertex (could be PV)");

    declareInterface<IParticleTupleTool>(this);
}

StatusCode TupleToolDecayFit::initialize()
{
    StatusCode sc = TupleToolBase::initialize();
    if( sc.isFailure() ) return sc;

    m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc(), this );
    if ( !m_dva ) return Error("Couldn't get parent IDVAlgorithm", StatusCode::FAILURE);

    return sc;
}

StatusCode TupleToolDecayFit::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
    return StatusCode{ TupleToolBase::finalize() && sc };
}

StatusCode TupleToolDecayFit::fill( const LHCb::Particle* mother, const LHCb::Particle* P,
                                    const std::string& head, Tuples::Tuple& tuple )
{
    LHCb::DecayTree tree(*P);
    LHCb::Particle* treeHead = tree.head();

    // 1) Initialise vector of measured parameters m in LHCb frame
    int dimM = 22;

    // Get primary vertex
    std::vector<const VertexBase*> originVtx;
    originVtx = originVertex(mother, P);
    if( originVtx.empty() ){return Error("Can't get an origin vertex");}
    ROOT::Math::XYZPoint PV( originVtx[0]->position().X(), originVtx[0]->position().Y(), originVtx[0]->position().Z() ); // originVtx[0] is the best PV in the event
 
    // Get taus and their daughters
    LHCb::Particle::ConstVector tauList;
    if(mother->charge() > 0) // B+ -> tau+ tau- K+
    {
        LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList); // tau+
        LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList); // tau-
    }
    else // B- -> tau- tau+ K-
    {
        LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList); // tau-
        LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList); // tau+
    }
    if ( tauList.size() != 2 ) 
    {
        return Error("Can't find two tau decays");
    }

    LHCb::Particle::ConstVector taup_daus = tauList[0]->daughtersVector();
    LHCb::Particle::ConstVector taum_daus = tauList[1]->daughtersVector();

    // Get K+
    const LHCb::Particle* Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(321) );
    if( Kp == 0 )
    {
        Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-321) );
    }
    if( Kp == 0 )
    {
        return Error("Can't find a K meson in decay tree");
    }

    // Refit tau+ decay vertex
    const IVertexFit* vtxFitter_taup = tool<IVertexFit>( "LoKi::VertexFitter" );
    LHCb::Vertex taup_refittedVertex;
    LHCb::Particle taup_refitted;
    taup_refitted.setParticleID(tauList[0]->particleID()); // tau+ (B+) or tau- (B-)

    vtxFitter_taup->fit(taup_daus, taup_refittedVertex, taup_refitted);

    // Refit tau- decay vertex
    const IVertexFit* vtxFitter_taum = tool<IVertexFit>( "LoKi::VertexFitter" );
    LHCb::Vertex taum_refittedVertex;
    LHCb::Particle taum_refitted;
    taum_refitted.setParticleID(tauList[1]->particleID()); // tau- (B+) pr tau+ (B-)

    vtxFitter_taum->fit(taum_daus, taum_refittedVertex, taum_refitted);

    // Get tau+ refitted decay vertex and visible 4-momentum (3pi)+ 
    ROOT::Math::XYZPoint DV1( taup_refitted.endVertex()->position().X(), taup_refitted.endVertex()->position().Y(), taup_refitted.endVertex()->position().Z() );
    Gaudi::LorentzVector P3pi1( taup_refitted.momentum().X(), taup_refitted.momentum().Y(), taup_refitted.momentum().Z(), taup_refitted.momentum().E() );

    // Get tau- refitted decay vertex and visible 4-momentum (3pi)-
    ROOT::Math::XYZPoint DV2( taum_refitted.endVertex()->position().X(), taum_refitted.endVertex()->position().Y(), taum_refitted.endVertex()->position().Z() );
    Gaudi::LorentzVector P3pi2( taum_refitted.momentum().X(), taum_refitted.momentum().Y(), taum_refitted.momentum().Z(), taum_refitted.momentum().E() );

    // Get K+ reference point and 4-momentum
    ROOT::Math::XYZPoint refPoint_Kp = Kp->referencePoint();
    ROOT::Math::XYZVector pK( Kp->momentum().X(), Kp->momentum().Y(), Kp->momentum().Z() );

    // const LHCb::ProtoParticle *Kp_proto = Kp->proto();
    // const LHCb::Track *Kp_track = Kp_proto->track();
    // const LHCb::State *Kp_state = Kp_track->stateAt(LHCb::State::Location(1));

    // ROOT::Math::XYZPoint refPoint_Kp( Kp_state->position().X(), Kp_state->position().Y(), Kp_state->position().Z() );
    // ROOT::Math::XYZVector pK( Kp_state->momentum().X(), Kp_state->momentum().Y(), Kp_state->momentum().Z() );

    // There are 23 known parameters
    CLHEP::HepVector m(dimM);
    m(1) = PV.x(); // the indexing of CLHEP vectors starts with 1
    m(2) = PV.y();
    m(3) = PV.z();
    m(4) = DV1.x();
    m(5) = DV1.y();
    m(6) = DV1.z();
    m(7) = P3pi1.x();
    m(8) = P3pi1.y();
    m(9) = P3pi1.z();
    m(10) = P3pi1.t();
    m(11) = DV2.x();
    m(12) = DV2.y();
    m(13) = DV2.z();
    m(14) = P3pi2.x();
    m(15) = P3pi2.y();
    m(16) = P3pi2.z();
    m(17) = P3pi2.t();
    m(18) = refPoint_Kp.x();
    m(19) = refPoint_Kp.y();
    m(20) = pK.x();
    m(21) = pK.y();
    m(22) = pK.z();

    // 2) Initialise covariance matrix of measured quantities in LHCb frame 
    // Get PV covariance matrix
    // auto pv = m_dva->bestVertex( P );
    Gaudi::SymMatrix3x3 PV_cov = originVtx[0]->covMatrix(); // pv[0] is the best PV
    // Get DV1 + P3pi1 covariance matrix
    Gaudi::SymMatrix7x7 taup_cov = taup_refitted.covMatrix();
    // Get DV2 + P3pi2 covariance matrix
    Gaudi::SymMatrix7x7 taum_cov = taum_refitted.covMatrix();
    // Get RP + PK covariance matrix
    Gaudi::SymMatrix7x7 Kp_cov = Kp->covMatrix();
    // Gaudi::SymMatrix6x6 Kp_cov = Kp_state->posMomCovariance();

    CLHEP::HepSymMatrix V(dimM); // covariance matrix of measured quantities (PV,tau+,tau-,K+ RP_T, B+ P4)
    for(int i = 1; i <= dimM; i++)
    {
        for(int j = 1; j <= dimM; j++)
        {
            V(i,j) = 0.;
        }
    }
    // Primary vertex (1,2,3)
    for(int i = 1; i <= 3; i++)
    {
        for(int j = 1; j <= 3; j++)
        {
            V(i, j) = PV_cov(i-1,j-1);
        }
    }
    // DV1 and P3pi1 (4,5,6, 7,8,9,10)
    for(int i = 4; i <= 10; i++)
    {
        for(int j = 4; j <= 10; j++)
        {
            V(i, j) = taup_cov(i-4,j-4);
        }
    }
    // DV2 and p3pi2 (11,12,13 ,14,15,16,17)
    for(int i = 11; i <= 17; i++)
    {
        for(int j = 11; j <= 17; j++)
        {
            V(i, j) = taum_cov(i-11,j-11);
        }
    }
    // RP and pK (18,19, 20,21,22)
    for(int i = 18; i <= 22; i++ )
    {
        for(int j = 18; j <= 22; j++)
        {
            if( (i < 20) && (j < 20) )
            {
                V(i,j) = Kp_cov(i-18,j-18);
            }
            else if( (i < 20) && (j >= 20) )
            {
                V(i,j) = Kp_cov(i-18,j-17);
            }
            else if( (i >= 20) && (j < 20) )
            {
                V(i,j) = Kp_cov(i-17,j-18);
            }
            else if( (i >= 20) && (j >= 20) )
            {
                V(i,j) = Kp_cov(i-17,j-17);
            }
        }
    }

    // Measurement vector m
    for(int i = 1; i <= dimM; i++)
    {
        tuple->column(Form("df_m_%i",i), m(i));
    }

    // std::vector<Double_t> m_data;
    // for(int i = 0; i < dimM; i++)
    // {
    //     m_data.push_back(m(i+1));
    // }
    // const std::vector<Double_t>& m_data_const = m_data; 
    // tuple->farray( "df_m", m_data_const.begin(), m_data_const.end(), "df_dimM", dimM+1 );

    // Covariance matrix V
    for(int i = 1; i <= dimM; i++)
    {
        for(int j = 1; j <= dimM; j++)
        {
            tuple->column(Form("df_V_%i_%i",i,j), V(i,j));
        }
    }

    // std::vector< std::vector<Double_t> > V_data;
    // for(int i = 0; i < dimM; i++)
    // {
    //     std::vector<Double_t> temp_vec;
    //     for(int j = 0; j < dimM; j++)
    //     {
    //         temp_vec.push_back(V(i+1,j+1));
    //     }
    //     V_data.push_back(temp_vec);
    // }
    // const std::vector< std::vector<Double_t> >& V_data_const = V_data;
    // tuple->fmatrix( "df_V", V_data_const, dimM, dimM, "df_dimM", dimM+1 );

    tuple->column("Kp_RP_X", refPoint_Kp.x());
    tuple->column("Kp_RP_Y", refPoint_Kp.y());
    tuple->column("Kp_RP_Z", refPoint_Kp.z());

    return StatusCode(true);
}

std::vector<const VertexBase*> TupleToolDecayFit::originVertex( const Particle* mother, const Particle* P ) const
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
DECLARE_COMPONENT( TupleToolDecayFit )