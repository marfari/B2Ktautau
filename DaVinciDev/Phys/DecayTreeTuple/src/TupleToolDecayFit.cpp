
#include "TupleToolDecayFit.h"
#include "TMath.h"
#include <algorithm>
#include "TMinuit.h"

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
    int dim = 24;
    // Get primary vertex
    std::vector<const VertexBase*> originVtx;
    originVtx = originVertex(mother, P);
    if( originVtx.empty() ){return Error("Can't get an origin vertex");}
    ROOT::Math::XYZPoint PV( originVtx[0]->position().X(), originVtx[0]->position().Y(), originVtx[0]->position().Z() ); // originVtx[0] is the best PV in the event

    // Get tau+ decay vertex and visible 4-momentum (3pi)
    const LHCb::Particle* taup = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-15) );
    ROOT::Math::XYZPoint DV1( taup->endVertex()->position().X(), taup->endVertex()->position().Y(), taup->endVertex()->position().Z() );
    Gaudi::LorentzVector P3pi1( taup->momentum().X(), taup->momentum().Y(), taup->momentum().Z(), taup->momentum().E() );

    // Get tau- decay vertex and visible 4-momentum (3pi)
    const LHCb::Particle* taum = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(15) );
    ROOT::Math::XYZPoint DV2 ( taum->endVertex()->position().X(), taum->endVertex()->position().Y(), taum->endVertex()->position().Z() );
    Gaudi::LorentzVector P3pi2( taum->momentum().X(), taum->momentum().Y(), taum->momentum().Z(), taum->momentum().E() );

    // Get B+ visible 4-momentum (6piK)
    const LHCb::Particle* Bp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(521) );
    if( Bp == 0 )
    {
        Bp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-521) );
    }
    if( Bp == 0 )
    {
        return Error("Can't find a B meson in decay tree");
    }
    Gaudi::LorentzVector P6piK( Bp->momentum().X(), Bp->momentum().Y(), Bp->momentum().Z(), Bp->momentum().E() );

    // Get reference point in K+ trajectory and K+ 4-momentum
    const LHCb::Particle* Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(321) );
    if( Kp == 0 )
    {
        Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-321) );
    }
    if( Kp == 0 )
    {
        return Error("Can't find a K meson in decay tree");
    }
   
    ROOT::Math::XYZPoint refPoint_Kp = Kp->referencePoint();
    Gaudi::LorentzVector PK( Kp->momentum().X(), Kp->momentum().Y(), Kp->momentum().Z(), Kp->momentum().E() );

    // There are 24 known parameters
    CLHEP::HepVector m(dim);
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
    m(20) = refPoint_Kp.z();
    m(21) = PK.x();
    m(22) = PK.y();
    m(23) = PK.z();
    m(24) = PK.t();

    // 2) Initialise covariance matrix of measured quantities in LHCb frame 

    // Get PV covariance matrix
    auto pv = m_dva->bestVertex( P );
    Gaudi::SymMatrix3x3 PV_cov = pv[0].covMatrix(); // pv[0] is the best PV
    // Get DV1 + P3pi1 covariance matrix
    Gaudi::SymMatrix7x7 taup_cov = taup->covMatrix();
    // Get DV2 + P3pi2 covariance matrix
    Gaudi::SymMatrix7x7 taum_cov = taum->covMatrix();
    // Get covariance matrix of reference point in K+ trajectory + PK
    Gaudi::SymMatrix3x3 RP_cov = Kp->posCovMatrix();
    Gaudi::SymMatrix4x4 Kp_mom_cov = Kp->momCovMatrix();

    CLHEP::HepSymMatrix V(dim, dim); // covariance matrix of measured quantities
    // Primary vertex
    for(int i = 1; i <= 3; i++)
    {
        for(int j = 1; j <= 3; j++)
        {
            V(i, j) = PV_cov(i-1,j-1);
        }
    }
    // DV1 and P3pi1
    for(int i = 4; i <= 10; i++)
    {
        for(int j = 4; j <= 10; j++)
        {
            V(i, j) = taup_cov(i-4,j-4);
        }
    }
    // DV2 and P3pi2
    for(int i = 11; i <= 17; i++)
    {
        for(int j = 11; j <= 17; j++)
        {
            V(i, j) = taum_cov(i-11,j-11);
        }
    }
    // RP 
    for(int i = 18; i <= 20; i++)
    {
        for(int j = 18; j <= 20; j++)
        {
            V(i, j) = RP_cov(i-18,j-18);
        }
    }
    // PK
    for(int i = 21; i <= 24; i++)
    {
        for(int j = 21; j <= 24; j++ )
        {
            V(i, j) = Kp_mom_cov(i-21,j-21);
        }
    }

    // (?) Change of variables E3pi1 -> m3pi1^2 , E3pi2->m3pi2^2 and E6piK -> m6piK^2 // check if this is necessary in the long run, will use energy for now

    // Measurement vector
    for(int i = 1; i <= dim; i++)
    {
        tuple->column(Form("df_m_%i",i), m(i));
    }

    // Covariance matrix
    for(int i = 1; i <= dim; i++)
    {
        for(int j = 1; j <= dim; j++)
        {
            tuple->column(Form("df_V_%i%i",i,j), V(i,j));
        }
    }

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