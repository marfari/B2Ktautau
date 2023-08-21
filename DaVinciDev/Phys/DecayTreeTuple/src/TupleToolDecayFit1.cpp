
#include "TupleToolDecayFit1.h"
#include "TMath.h"
#include <algorithm>
#include "TMinuit.h"

using namespace LHCb;
using namespace std;

TupleToolDecayFit1::TupleToolDecayFit1( const std::string& type, const std::string& name, const IInterface* parent ):TupleToolBase( type, name, parent )
{
    declareProperty( "StateProvider", m_stateprovider ) ;
    declareProperty( "constrainToOriginVertex", m_constrainToOriginVertex = false,
                    "Do a refit constraining to Origin Vertex (could be PV)");

    declareInterface<IParticleTupleTool>(this);
}

StatusCode TupleToolDecayFit1::initialize()
{
    StatusCode sc = TupleToolBase::initialize();
    if( sc.isFailure() ) return sc;

    m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc(), this );
    if ( !m_dva ) return Error("Couldn't get parent IDVAlgorithm", StatusCode::FAILURE);

    return sc;
}

StatusCode TupleToolDecayFit1::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
    return StatusCode{ TupleToolBase::finalize() && sc };
}

StatusCode TupleToolDecayFit1::fill( const LHCb::Particle* mother, const LHCb::Particle* P,
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
    LHCb::Particle::ConstVector tauList;
    LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(-15), tauList); // tau+
    LHCb::DecayTree::findInTree( treeHead, LHCb::ParticleID(15), tauList); // tau-
    if ( tauList.size() != 2 ) 
    {
        return Error("Can't find two tau decays");
    }

    const LHCb::Particle* taup = tauList[0];    
    ROOT::Math::XYZPoint DV1( taup->endVertex()->position().X(), taup->endVertex()->position().Y(), taup->endVertex()->position().Z() );
    LHCb::Particle::ConstVector taup_pions = taup->daughtersVector();
    Gaudi::LorentzVector P1( taup_pions[0]->momentum().X(), taup_pions[0]->momentum().Y(), taup_pions[0]->momentum().Z(), taup_pions[0]->momentum().E() );
    Gaudi::LorentzVector P2( taup_pions[1]->momentum().X(), taup_pions[1]->momentum().Y(), taup_pions[1]->momentum().Z(), taup_pions[1]->momentum().E() );
    Gaudi::LorentzVector P3( taup_pions[2]->momentum().X(), taup_pions[2]->momentum().Y(), taup_pions[2]->momentum().Z(), taup_pions[2]->momentum().E() );

    // Get tau- decay vertex + tau- pions momenta
    // const LHCb::Particle* taum = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(15) );
    const LHCb::Particle* taum = tauList[1];
    ROOT::Math::XYZPoint DV2( taum->endVertex()->position().X(), taum->endVertex()->position().Y(), taum->endVertex()->position().Z() );
    LHCb::Particle::ConstVector taum_pions = taum->daughtersVector();
    Gaudi::LorentzVector P4( taum_pions[0]->momentum().X(), taum_pions[0]->momentum().Y(), taum_pions[0]->momentum().Z(), taum_pions[0]->momentum().E() );
    Gaudi::LorentzVector P5( taum_pions[1]->momentum().X(), taum_pions[1]->momentum().Y(), taum_pions[1]->momentum().Z(), taum_pions[1]->momentum().E() );
    Gaudi::LorentzVector P6( taum_pions[2]->momentum().X(), taum_pions[2]->momentum().Y(), taum_pions[2]->momentum().Z(), taum_pions[2]->momentum().E() );

    // Get K+ reference point + 4-momentum
    const LHCb::Particle* Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(321) );
    if( Kp == 0 )
    {
        Kp = LHCb::DecayTree::findFirstInTree( treeHead, LHCb::ParticleID(-321) );
    }
    if( Kp == 0 )
    {
        return Error("Can't find a K meson in decay tree");
    }
    ROOT::Math::XYZPoint RPK = Kp->referencePoint();
    Gaudi::LorentzVector PK( Kp->momentum().X(), Kp->momentum().Y(), Kp->momentum().Z(), Kp->momentum().E() );

    // There are 24 known parameters
    CLHEP::HepVector m(dimM);
    m(1) = PV.x(); // the indexing of CLHEP vectors starts with 1
    m(2) = PV.y();
    m(3) = PV.z();
    m(4) = DV1.x();
    m(5) = DV1.y();
    m(6) = DV1.z();
    m(7) = P1.x();
    m(8) = P1.y();
    m(9) = P1.z();
    m(10) = P2.x();
    m(11) = P2.y();
    m(12) = P2.z();
    m(13) = P3.x();
    m(14) = P3.y();
    m(15) = P3.z();
    m(16) = DV2.x();
    m(17) = DV2.y();
    m(18) = DV2.z();
    m(19) = P4.x();
    m(20) = P4.y();
    m(21) = P4.z();
    m(22) = P5.x();
    m(23) = P5.y();
    m(24) = P5.z();
    m(25) = P6.x();
    m(26) = P6.y();
    m(27) = P6.z();
    m(28) = RPK.x();
    m(29) = RPK.y();
    m(30) = PK.x();
    m(31) = PK.y();
    m(32) = PK.z();

    // 2) Initialise covariance matrix of measured quantities in LHCb frame 
    // Get PV covariance matrix
    Gaudi::SymMatrix3x3 PV_cov = originVtx[0]->covMatrix(); // pv[0] is the best PV
    // Get DV1 covariance matrix
    Gaudi::SymMatrix3x3 DV1_cov = taup->posCovMatrix();
    // Get tau+ pions momentum covariance matrix
    Gaudi::SymMatrix4x4 pi1_mom_cov = taup_pions[0]->momCovMatrix();
    Gaudi::SymMatrix4x4 pi2_mom_cov = taup_pions[1]->momCovMatrix();
    Gaudi::SymMatrix4x4 pi3_mom_cov = taup_pions[2]->momCovMatrix();
    // Get DV2 covariance matrix
    Gaudi::SymMatrix3x3 DV2_cov = taum->posCovMatrix();
    // Get tau- pions momentum covariance matrix
    Gaudi::SymMatrix4x4 pi4_mom_cov = taum_pions[0]->momCovMatrix();
    Gaudi::SymMatrix4x4 pi5_mom_cov = taum_pions[1]->momCovMatrix();
    Gaudi::SymMatrix4x4 pi6_mom_cov = taum_pions[2]->momCovMatrix();
    // Get K+ reference point + momentum covariance matrix
    Gaudi::SymMatrix7x7 Kp_cov = Kp->covMatrix();

    Gaudi::SymMatrix7x7 taup_cov = taup->covMatrix();
    Gaudi::SymMatrix7x7 taum_cov = taum->covMatrix();
    Gaudi::SymMatrix7x7 Bp_cov = mother->covMatrix();

    CLHEP::HepSymMatrix V(dimM); 
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
            V(i,j) = PV_cov(i-1,j-1);
        }
    }
    // DV1 (4,5,6)
    for(int i = 4; i <= 6; i++)
    {
        for(int j = 4; j <= 6; j++)
        {
            V(i,j) = DV1_cov(i-4,j-4);
        }
    }
    // p1 (7,8,9)
    for(int i = 7; i <= 9; i++)
    {
        for(int j = 7; j <= 9; j++)
        {
            V(i,j) = pi1_mom_cov(i-7,j-7);
        }
    }
    // p2 (10,11,12)
    for(int i = 10; i <= 12; i++)
    {
        for(int j = 10; j <= 12; j++)
        {
            V(i,j) = pi2_mom_cov(i-10,j-10);
        }
    }
    // p3 (13,14,15)
    for(int i = 13; i <= 15; i++)
    {
        for(int j = 13; j <= 15; j++)
        {
            V(i,j) = pi3_mom_cov(i-13,j-13);
        }
    }
    // DV2 (16,17,18)
    for(int i = 16; i <= 18; i++)
    {
        for(int j = 16; j <= 18; j++)
        {
            V(i,j) = DV2_cov(i-16,j-16);
        }
    }
    // p4 (19,20,21)
    for(int i = 19; i <= 21; i++)
    {
        for(int j = 19; j <= 21; j++)
        {
            V(i,j) = pi4_mom_cov(i-19,j-19);
        }
    }
    // p5 (22,23,24)
    for(int i = 22; i <= 24; i++)
    {
        for(int j = 22; j <= 24; j++)
        {
            V(i,j) = pi5_mom_cov(i-22,j-22);
        }
    }
    // p6 (25,26,27)
    for(int i = 25; i <= 27; i++)
    {
        for(int j = 25; j <= 27; j++)
        {
            V(i,j) = pi6_mom_cov(i-25,j-25);
        }
    }
    // RP + pK (28,29, 30,31,32)
    for(int i = 28; i <= 32; i++ )
    {
        for(int j = 28; j <= 32; j++)
        {
            if( (i < 30) && (j < 30) )
            {
                V(i,j) = Kp_cov(i-28,j-28);
            }
            else if( (i < 30) && (j >= 30) )
            {
                V(i,j) = Kp_cov(i-28,j-27);
            }
            else if( (i >= 30) && (j < 30) )
            {
                V(i,j) = Kp_cov(i-27,j-28);
            }
            else if( (i >= 30) && (j >= 30) )
            {
                V(i,j) = Kp_cov(i-27,j-27);
            }
        }
    }

    // Measurement vector m
    for(int i = 1; i <= dimM; i++)
    {
        tuple->column(Form("df_m_%i",i), m(i));
    }

    // Covariance matrix V
    for(int i = 1; i <= dimM; i++)
    {
        for(int j = 1; j <= dimM; j++)
        {
            tuple->column(Form("df_V_%i_%i",i,j), V(i,j));
        }
    }
    
    tuple->column("Kp_RP_X", RPK.x());
    tuple->column("Kp_RP_Y", RPK.y());
    tuple->column("Kp_RP_Z", RPK.z());

    for(int i = 0; i < 7; i++)
    {
        for(int j = 0; j < 7; j++)
        {
            tuple->column(Form("taup_cov_%i_%i",i,j), taup_cov(i,j));
            tuple->column(Form("taum_cov_%i_%i",i,j), taum_cov(i,j));
            tuple->column(Form("Bp_cov_%i_%i",i,j), Bp_cov(i,j));
        }
    }

    return StatusCode(true);
}

std::vector<const VertexBase*> TupleToolDecayFit1::originVertex( const Particle* mother, const Particle* P ) const
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
DECLARE_COMPONENT( TupleToolDecayFit1 )