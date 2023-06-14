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
#include <algorithm>
#include <boost/foreach.hpp>
#include <iomanip>

#include "Event/Particle.h"
#include "FitParams.h"
#include "InternalParticle.h"
#include "MissingParticle.h"
#include "RecoTrack.h"

#include "LHCbMath/GeomFun.h"
#include "LHCbMath/Line.h"

namespace DecayTreeFitter {

  extern int vtxverbose;

  inline bool sortByType( const ParticleBase* lhs, const ParticleBase* rhs ) {
    int  lhstype = lhs->type();
    int  rhstype = rhs->type();
    bool rc      = false;
    if ( lhstype == rhstype && lhstype == ParticleBase::kRecoTrack )
      rc = lhs->particle().pt() > rhs->particle().pt();
    else if ( lhs->particle().daughters().size() > 0 && rhs->particle().daughters().size() > 0 )
      rc = lhs->nFinalChargedCandidates() > rhs->nFinalChargedCandidates();
    else
      rc = lhstype < rhstype;
    return rc;
  }

  InternalParticle::InternalParticle( const LHCb::Particle& bc, const ParticleBase* mother,
                                      const Configuration& config )
      : ParticleBase( bc, mother ), m_lifetimeconstraint( false ) {
    BOOST_FOREACH ( const LHCb::Particle* daughter, bc.daughters() )
      addDaughter( *daughter, config );
    // copy constraints
    m_lifetimeconstraint = false; // bc && bc->constraint(BtaConstraint::Life) ;
  }

  bool compTrkTransverseMomentum( const RecoTrack* lhs, const RecoTrack* rhs ) {
    return lhs->particle().pt() > rhs->particle().pt();
  }

  ErrCode InternalParticle::initPar1( FitParams* fitparams ) {
    // This is the most complicated part of the vertexer: an
    // initialization that always works.

    // There are two ways out: If the LHCb::Particle was vertexed
    // before, we can rely on the existing vertex (case A). If not, we
    // need to estimate the vertex position from the daughters; that
    // is very complicated (case B). The momentum is always
    // initialized from the sum of the daughter four-vectors. In the
    // end, it doesn't really matter.

    // FIX ME: Currently, this scheme does not work for B->K*D0, with
    // D0->pi0Ks, because the D0 is initialized before there is a B
    // vertex.

    if ( vtxverbose >= 3 )
      std::cout << "InternalParticle::initPar: " << particle().particleID().pid() << " " << daughters().size() << " "
                << hasPosition() << " " << posIndex() << std::endl;

    ErrCode status;
    int     posindex = posIndex();

    // logic check: we do not want to call this routine for resonances.
    assert( hasPosition() );

    // Start with origin
    for ( int row = 1; row <= 3; ++row ) fitparams->par()( row + posindex ) = 0;

    // Step 1: pre-initialization of all daughters
    for ( daucontainer::const_iterator it = begin(); it != end(); ++it ) status |= ( *it )->initPar1( fitparams );

    // Step 2: initialize the vertex. if we are lucky, we had a
    // 'btaresonant daughter, and we are already done.
    if ( fitparams->par()( posindex + 1 ) == 0 && fitparams->par()( posindex + 2 ) == 0 &&
         fitparams->par()( posindex + 3 ) == 0 ) {

      const LHCb::Vertex* vtx = particle().endVertex();
      if ( vtx && vtx->position() != Gaudi::XYZPoint( 0, 0, 0 )
           // && vtx->technique() ==
      ) {
        // we found an existing valid vertex. that's fine as well ...
        Gaudi::XYZPoint point          = vtx->position();
        fitparams->par( posindex + 1 ) = point.x();
        fitparams->par( posindex + 2 ) = point.y();
        fitparams->par( posindex + 3 ) = point.z();
        if ( vtxverbose >= 2 ) std::cout << "using existing vertex: " << point << std::endl;

      } else {
        // Case B: the hard way ... use the daughters to estimate the
        // vertex. First we check if there are sufficient tracks
        // attached to this vertex. If so, estimate the poca of the
        // two tracks with the highest momentum. This will work for
        // the majority of the cases. If there are not sufficient
        // tracks, add the composites and take the two with the best
        // doca.

        // create a vector with all daughters that constitute a
        // 'trajectory' (ie tracks, composites and daughters of
        // resonances.)
        daucontainer alldaughters;
        collectVertexDaughters( alldaughters, posindex );
        if ( vtxverbose >= 2 ) {
          std::cout << "number of daughters for initializing vertex: " << name() << " " << alldaughters.size()
                    << std::endl;
          if ( vtx ) std::cout << "Original vertex pos: " << vtx->position() << std::endl;
        }

        // select daughters that are either charged, or have an initialized vertex
        daucontainer            vtxdaughters;
        std::vector<RecoTrack*> trkdaughters;
        for ( daucontainer::const_iterator it = alldaughters.begin(); it != alldaughters.end(); ++it ) {
          if ( ( *it )->type() == ParticleBase::kRecoTrack ) {
            trkdaughters.push_back( static_cast<RecoTrack*>( *it ) );
          } else if ( ( *it )->hasPosition() && fitparams->par( ( *it )->posIndex() + 1 ) != 0 ) {
            vtxdaughters.push_back( *it );
          }
        }

        if ( trkdaughters.size() >= 2 ) {
          // sort in pT. not very efficient, but it works.
          if ( trkdaughters.size() > 2 )
            std::sort( trkdaughters.begin(), trkdaughters.end(), compTrkTransverseMomentum );
          // now, just take the first two ...
          RecoTrack* dau1 = trkdaughters[0];
          RecoTrack* dau2 = trkdaughters[1];

          // get the poca of the two statevectors
          const LHCb::State&                                   state1 = dau1->state();
          const LHCb::State&                                   state2 = dau2->state();
          Gaudi::Math::Line<Gaudi::XYZPoint, Gaudi::XYZVector> line1( state1.position(), state1.slopes() );
          Gaudi::Math::Line<Gaudi::XYZPoint, Gaudi::XYZVector> line2( state2.position(), state2.slopes() );
          double                                               mu1( 0 ), mu2( 0 );
          Gaudi::Math::closestPointParams( line1, line2, mu1, mu2 );
          Gaudi::XYZPoint p1               = line1.position( mu1 );
          Gaudi::XYZPoint p2               = line2.position( mu2 );
          fitparams->par()( posindex + 1 ) = 0.5 * ( p1.x() + p2.x() );
          fitparams->par()( posindex + 2 ) = 0.5 * ( p1.y() + p2.y() );
          fitparams->par()( posindex + 3 ) = 0.5 * ( p1.z() + p2.z() );
          dau1->setFlightLength( mu1 );
          dau2->setFlightLength( mu2 );

        } else if ( trkdaughters.size() + vtxdaughters.size() >= 2 ) {
          // FIXME: this still needs a proper implementation.
          //
          // we get here if there is at most one track from the vertex
          // and if the composite daughter is long-lived. we could
          // 'vertex' the daughter and the track, but if the daughter
          // happens to be a V0, then that will now work very well
          // either. for now, if there is a track, take the first
          // point on the track. then if there is a daughter that is
          // more upstream, take that.
          //
          bool success = false;
          if ( !trkdaughters.empty() ) {
            const LHCb::State& state         = trkdaughters.front()->state();
            fitparams->par()( posindex + 1 ) = state.position().x();
            fitparams->par()( posindex + 2 ) = state.position().y();
            fitparams->par()( posindex + 3 ) = state.position().z();
            success                          = true;
          }
          for ( const auto& dau : vtxdaughters ) {
            const int dauposindex = dau->posIndex();
            if ( !success || ( fitparams->par()( posindex + 3 ) > fitparams->par()( dauposindex + 3 ) ) ) {
              success = true;
              for ( int i = 1; i <= 3; ++i ) fitparams->par()( posindex + i ) = fitparams->par()( dauposindex + i );
            }
          }
          /*

          // that's unfortunate: no enough charged tracks from this
          // vertex. need all daughters. create trajectories and use
          // normal TrkPoca.

          std::vector<LHCb::StateZTraj> trajectories ;
          for(vector<RecoTrack*>::const_iterator it = trkdaughters.begin() ;
          it != trkdaughters.end() ; ++it)
          trajectories.push_back( (*it)->traj() ) ;

          trajectories.push_back(&((*it)->particle().trkAbsFit()->traj())) ;

          std::vector<TrkLineTraj> linetrajectories ; // store trajectories of composites
          linetrajectories.reserve(  vtxdaughters.size() ) ;
          for(daucontainer::const_iterator it = vtxdaughters.begin() ;
          it != vtxdaughters.end() ; ++it) {
          //std::cout << (*it)->particle().pdtEntry()->name() << std::endl ;
          int dauposindex = (*it)->posIndex() ;
          int daumomindex = (*it)->momIndex() ;
          Gaudi::XYZPoint point(fitparams->par()(dauposindex+1),
          fitparams->par()(dauposindex+2),
          fitparams->par()(dauposindex+3)) ;
          Hep3Vector direction(fitparams->par()(daumomindex+1),
          fitparams->par()(daumomindex+2),
          fitparams->par()(daumomindex+3)) ;
          linetrajectories.push_back(TrkLineTraj(point,direction,1) ) ;
          trajectories.push_back(&(linetrajectories.back())) ;
          //daupoint = point ;
          }

          // we select the two trajectories with the best poca
          double docabest(99999);
          TrkErrCode pocastatus ;
          for( std::vector<const Trajectory*>::iterator it1 = trajectories.begin() ;
          it1 != trajectories.end(); ++it1 )
          for( std::vector<const Trajectory*>::iterator it2 = trajectories.begin() ;
          it2 != it1; ++it2 ) {
          TrkPoca poca(**it1,0.,**it2, 0.);
          Hep3Vector dir1 = (*it1)->direction(poca.flt1());
          Hep3Vector dir2 = (*it2)->direction(poca.flt2());
          double doca = poca.doca() ;
          if(fabs(doca)<fabs(docabest)) {
          Gaudi::XYZPoint pnt1 = (*it1)->position(poca.flt1());
          Gaudi::XYZPoint pnt2 = (*it2)->position(poca.flt2());
          fitparams->par()(posindex+1) = 0.5*(pnt1.x()+pnt2.x()) ;
          fitparams->par()(posindex+2) = 0.5*(pnt1.y()+pnt2.y()) ;
          fitparams->par()(posindex+3) = 0.5*(pnt1.z()+pnt2.z()) ;
          docabest = doca ;
          pocastatus = poca.status() ;
          }
          }
          */
        } else if ( mother() && mother()->posIndex() >= 0 ) {

          // let's hope the mother was initialized
          int posindexmother = mother()->posIndex();
          for ( int ipos = 1; ipos <= 3; ++ipos ) {
            fitparams->par()( posindex + ipos ) = fitparams->par()( posindexmother + ipos );
          }
        } else {
          // something is wrong!
          std::cout << "There are not sufficient geometric constraints to fit "
                    << "this decay tree. Perhaps you should add a beam constraint. " << std::endl;
          //<< particle().constraint(BtaConstraint::Beam)
          //   << std::endl
          //   << treeprinter.print(*bc()) << endmsg ;
          status |= ErrCode::badsetup;
        }
      }
    }

    // step 3: do the post initialization step of all daughters
    for ( daucontainer::const_iterator it = daughters().begin(); it != daughters().end(); ++it )
      ( *it )->initPar2( fitparams );

    // step 4: initialize the momentum by adding up the daughter 4-vectors
    initMom( fitparams );

    if ( vtxverbose >= 3 )
      std::cout << "End of initpar: " << name() << " (" << fitparams->par()( posindex + 1 ) << ","
                << fitparams->par()( posindex + 2 ) << "," << fitparams->par()( posindex + 3 ) << ")" << std::endl;

    return status;
  }

  ErrCode InternalParticle::initPar2( FitParams* fitparams ) {
    // FIX ME: in the unfortunate case (the B-->D0K*- above) that our
    // vertex is still the origin, we copy the mother vertex.
    int posindex = posIndex();
    if ( hasPosition() && mother() && fitparams->par( posindex + 1 ) == 0 && fitparams->par( posindex + 2 ) == 0 &&
         fitparams->par( posindex + 3 ) == 0 ) {
      int posindexmom = mother()->posIndex();
      for ( int irow = 1; irow <= 3; ++irow ) fitparams->par( posindex + irow ) = fitparams->par( posindexmom + irow );
    }
    // step 5: initialize the lifetime
    return initTau( fitparams );
  }

  ErrCode InternalParticle::initMom( FitParams* fitparams ) const {
    int momindex = momIndex();
    // reset
    for ( int irow = 1; irow <= 4; ++irow ) fitparams->par( momindex + irow ) = 0;

    // now add daughter momenta
    for ( daucontainer::const_iterator it = begin(); it != end(); ++it ) {
      int    daumomindex = ( *it )->momIndex();
      double e2( 0 );
      int maxrow = ( *it )->hasEnergy() ? 4 : 3; 
      for ( int irow = 1; irow <= maxrow; ++irow ) {
        double px = fitparams->par()( daumomindex + irow );
        e2 += px * px;
        fitparams->par( momindex + irow ) += px;
      }
      if ( maxrow == 3 ) { // if daughter does not have an energy
        double mass = ( *it )->pdtMass();
        fitparams->par( momindex + 4 ) += std::sqrt( e2 + mass * mass ); 
      }
    }
    
    // if there is a mass constraint, ignore what we have just
    // computed for the energy, and just insert the mass
    if ( hasMassConstraint() ) { 
      double mass = pdtMass();
      double p2   = 0;
      for ( int irow = 1; irow <= 3; ++irow ) {
        double px = fitparams->par( momindex + irow );
        p2 += px * px;
      }
      fitparams->par( momindex + 4 ) = std::sqrt( p2 + mass * mass );
    }

    return ErrCode::success;
  }

  ErrCode InternalParticle::projectKineConstraint( const FitParams& fitparams, Projection& p ) const {
    // these are in fact four independent constraints. i'll filter
    // them as one, making the code simpler at the expense of a bit of
    // CPU.

    // first add the mother
    int momindex = momIndex();
    for ( int imom = 1; imom <= 4; ++imom ) {
      p.r( imom )                  = fitparams.par()( momindex + imom );
      p.H( imom, momindex + imom ) = 1;
    }

    // now add the daughters
    for ( daucontainer::const_iterator it = daughters().begin(); it != daughters().end(); ++it ) {
      int    daulenindex = ( *it )->lenIndex();
      int    daumomindex = ( *it )->momIndex();
      double mass        = ( *it )->pdtMass();
      double e2          = mass * mass;
      int    maxrow      = ( *it )->hasEnergy() ? 4 : 3;
      for ( int imom = 1; imom <= maxrow; ++imom ) {
        double px = fitparams.par()( daumomindex + imom );
        e2 += px * px;
        p.r( imom ) += -px;
        p.H( imom, daumomindex + imom ) = -1;
      }
      if ( maxrow == 3 ) { // if daughter does not have an energy
        // treat the energy for particles that are parameterized with p3
        double energy = sqrt( e2 );
        p.r( 4 ) += -energy;
        for ( int jmom = 1; jmom <= 3; ++jmom ) {
          double px                    = fitparams.par()( daumomindex + jmom );
          p.H( 4, daumomindex + jmom ) = -px / energy;
        }
      } else if ( false && daulenindex >= 0 && ( *it )->charge() != 0 ) {

        double       tau          = fitparams.par()( daulenindex + 1 );
        double       lambda       = bFieldOverC() * ( *it )->charge();
        double       px0          = fitparams.par()( daumomindex + 1 );
        double       py0          = fitparams.par()( daumomindex + 2 );
        double       pt0          = sqrt( px0 * px0 + py0 * py0 );
        const double posprecision = 1e-4; // 1mu
        if ( fabs( pt0 * lambda * tau * tau ) > posprecision ) {
          double sinlt = sin( lambda * tau );
          double coslt = cos( lambda * tau );
          double px    = px0 * coslt - py0 * sinlt;
          double py    = py0 * coslt + px0 * sinlt;
          p.r( 1 ) += px0 - px;
          p.r( 2 ) += py0 - py;
          p.H( 1, daumomindex + 1 ) += 1 - coslt;
          p.H( 1, daumomindex + 2 ) += sinlt;
          p.H( 1, daulenindex + 1 ) += lambda * py;
          p.H( 2, daumomindex + 1 ) += -sinlt;
          p.H( 2, daumomindex + 2 ) += 1 - coslt;
          p.H( 2, daulenindex + 1 ) += -lambda * px;
        }
      }
    }
    
    return ErrCode::success;
  }

  ErrCode InternalParticle::projectLifeTimeConstraint( const FitParams&, Projection& ) const {
    std::cout << "Not yet implemented lifetime constraint!" << std::endl;
    // int lenindex = lenIndex() ;
    //     assert(lenindex>=0) ;
    //     double tau = pdtTau() ;
    //     p.r(1)            = fitparams.par()(lenindex+1) - tau ;
    //     p.Vfast(1,1)      = tau*tau ;
    //     p.H(1,lenindex+1) = 1 ;
    return ErrCode::success;
  }

  //   ErrCode InternalParticle::projectScaleFactorConstraint( const FitParams& fitparams, Projection& p) const { // Added by Maria
  //   // This is necessary for neutrino initialisatin calculation to work
  //   // L1x = Ptau1x / (DV1x - BVx) ; L1y = Ptau1y / (DV1y - BVy) ; L1z = Ptau1z / (Dv1z - BVz)
  //   // L2x = Ptau2x / (DV2x - BVx) ; L2y = Ptau2y / (DV2y - BVy) ; L2z = Ptau2z / (Dv2z - BVz)
  //   // Lx = Pbx / (BVx - PVx) ; Ly = Pby / (BVy - PVy) ; Lz = Pbz / (BVz - PVz)
  //   // we want the scale factors to be the same for the 3 momentum components:
  //   // L1x - L1y = 0 and L1x - L1z = 0
  //   // L2x - L2y = 0 and L2x - L2z = 0
  //   // Lx - Ly = 0 and Lx - Lz = 0

  //   // iterate over daughters
  //   // primary vertex see DTF tuple; what are the PV indices?

  //   // B+ momentum and position indices
  //   int B_momindex = momIndex();
  //   int B_posindex = posIndex();

  //   //std::cout << "Particle ID = " << particle().particleID().pid() << std::endl;

  //   // tau+ and tau- momentum and position indices
  //   // check that tau+ = daughters()[1] and tau- = daughters()[2]
  //   int tau1_momindex = daughters()[1]->momIndex();
  //   int tau2_momindex = daughters()[2]->momIndex();
  //   int tau1_posindex = daughters()[1]->posIndex();
  //   int tau2_posindex = daughters()[2]->posIndex();

  //   // Primary vertex
  //   // the coordinates of the PV are the last 3 parameters in the model 
  //   int PV_posindex = 51; // only true for our decay model; there should be a way to get this info generally but idk how
  //   double PVx = fitparams.par()( PV_posindex + 1 );
  //   double PVy = fitparams.par()( PV_posindex + 2 );
  //   double PVz = fitparams.par()( PV_posindex + 3 );
    
  //   // Momentum
  //   double Ptau1x = fitparams.par()( tau1_momindex + 1 );
  //   double Ptau1y = fitparams.par()( tau1_momindex + 2 );
  //   double Ptau1z = fitparams.par()( tau1_momindex + 3 );

  //   double Ptau2x = fitparams.par()( tau2_momindex + 1 );
  //   double Ptau2y = fitparams.par()( tau2_momindex + 2 );
  //   double Ptau2z = fitparams.par()( tau2_momindex + 3 );

  //   double Pbx = fitparams.par()( B_momindex + 1 );
  //   double Pby = fitparams.par()( B_momindex + 2 );
  //   double Pbz = fitparams.par()( B_momindex + 3 );

  //   // Decay vertices
  //   double DV1x = fitparams.par()( tau1_posindex + 1 );
  //   double DV1y = fitparams.par()( tau1_posindex + 2 );
  //   double DV1z = fitparams.par()( tau1_posindex + 3 );

  //   double DV2x = fitparams.par()( tau2_posindex + 1 );
  //   double DV2y = fitparams.par()( tau2_posindex + 2 );
  //   double DV2z = fitparams.par()( tau2_posindex + 3 );

  //   double BVx = fitparams.par()( B_posindex + 1 );
  //   double BVy = fitparams.par()( B_posindex + 2 );
  //   double BVz = fitparams.par()( B_posindex + 3 );

  //   // Scale factors 
  //   // check if SF denominator is zero. These are cases in which the 3pi combination comes from the BV (thrown away after truth-matc)
  //   double d1x = DV1x - BVx;
  //   double d1y = DV1y - BVy;
  //   double d1z = DV1z - BVz;
  //   double d2x = DV2x - BVx;
  //   double d2y = DV2y - BVy;
  //   double d2z = DV2z - BVz;
  //   double dx = BVx - PVx;
  //   double dy = BVy - PVy;
  //   double dz = BVz - PVz;

  //   double eps = pow(10,-6);

  //   if( ( abs(d1x) < eps ) || ( abs(d1y) < eps ) || ( abs(d1z) < eps ) || ( abs(d2x) < eps ) || ( abs(d2y) < eps ) || ( abs(d2z) < eps ) || ( abs(dx) < eps ) || ( abs(dy < eps) ) || ( abs(dz) < eps ) )
  //   {
  //     std::cout << "Denominator of SF is zero. Setting it to 1" << std::endl;
  //     if( abs(d1x) < eps ){d1x = 1.;}
  //     if( abs(d1y) < eps ){d1y = 1.;}
  //     if( abs(d1z) < eps ){d1z = 1.;}
  //     if( abs(d2x) < eps ){d2x = 1.;}
  //     if( abs(d2y) < eps ){d2y = 1.;}
  //     if( abs(d2z) < eps ){d2z = 1.;}
  //     if( abs(dx) < eps ){dx = 1.;}
  //     if( abs(dy) < eps ){dy = 1.;}
  //     if( abs(dz) < eps ){dz = 1.;}
  //   }

  //   double L1x = Ptau1x/d1x;
  //   double L1y = Ptau1y/d1y;
  //   double L1z = Ptau1z/d1z;

  //   double L2x = Ptau2x/d2x;
  //   double L2y = Ptau2y/d2y;
  //   double L2z = Ptau2z/d2z;

  //   double Lx = Pbx/dx;
  //   double Ly = Pby/dy;
  //   double Lz = Pbz/dz;

  //   // L1x - L1y
  //   p.r( 1 ) = L1x - L1y;

  //   p.H( 1, tau1_momindex + 1 ) = 1/d1x; // Ptau1x
  //   p.H( 1, tau1_posindex + 1 ) = -Ptau1x/pow( d1x, 2 ); // DV1x
  //   p.H( 1, B_posindex + 1 ) = Ptau1x/pow( d1x, 2 ); // BVx
  //   p.H( 1, tau1_momindex + 2 ) = -1/d1y; // Ptau1y
  //   p.H( 1, tau1_posindex + 2 ) = Ptau1y/pow( d1y, 2 ); // DV1y
  //   p.H( 1, B_posindex + 2 ) = -Ptau1y/pow( d1y, 2 ); // BVy

  //   // L1x - L1z
  //   p.r( 2 ) = L1x - L1z;

  //   p.H( 2, tau1_momindex + 1 ) = 1/d1x; // Ptau1x
  //   p.H( 2, tau1_posindex + 1 ) = -Ptau1x/pow( d1x, 2 ); // DV1x
  //   p.H( 2, B_posindex + 1 ) = Ptau1x/pow( d1x, 2 ); // BVx
  //   p.H( 2, tau1_momindex + 3 ) = -1/d1z; // Ptau1z
  //   p.H( 2, tau1_posindex + 3 ) = Ptau1z/pow( d1z, 2 ); // DV1z
  //   p.H( 2, B_posindex + 3 ) = -Ptau1z/pow( d1z, 2 ); // BVz

  //   // L2x - L2y
  //   p.r( 3 ) = L2x - L2y;

  //   p.H( 3, tau2_momindex + 1 ) = 1/d2x; // Ptau2x
  //   p.H( 3, tau2_posindex + 1 ) = -Ptau2x/pow( d2x, 2 ); // DV2x
  //   p.H( 3, B_posindex + 1 ) = Ptau2x/pow( d2x, 2 ); // BVx
  //   p.H( 3, tau2_momindex + 2 ) = -1/d2y; // Ptau2y
  //   p.H( 3, tau2_posindex + 2 ) = Ptau2y/pow( d2y, 2 ); // DV2y
  //   p.H( 3, B_posindex + 2 ) = -Ptau2y/pow( d2y, 2 ); // BVy
 
  //   // L2x - L2z
  //   p.r( 4 ) = L2x - L2z;

  //   p.H( 4, tau2_momindex + 1 ) = 1/d2x; // Ptau2x
  //   p.H( 4, tau2_posindex + 1 ) = -Ptau2x/pow( d2x, 2 ); // DV2x
  //   p.H( 4, B_posindex + 1 ) = Ptau2x/pow( d2x, 2 ); // BVx
  //   p.H( 4, tau2_momindex + 3 ) = -1/d2z; // Ptau2z
  //   p.H( 4, tau2_posindex + 3 ) = Ptau2z/pow( d2z, 2 ); // DV2z
  //   p.H( 4, B_posindex + 3 ) = -Ptau2z/pow( d2z, 2 ); // BVz

  //   // Lx - Ly
  //   p.r( 5 ) = Lx - Ly;

  //   p.H( 5, B_momindex + 1 ) = 1/dx; // Pbx
  //   p.H( 5, B_posindex + 1 ) = -Pbx/pow( dx, 2 ); // BVx
  //   p.H( 5, PV_posindex + 1 ) = Pbx/pow( dx, 2 ); // PVx
  //   p.H( 5, B_momindex + 2 ) = -1/dy; // Pby
  //   p.H( 5, B_posindex + 2 ) = Pby/pow( dy, 2 ); // BVy
  //   p.H( 5, PV_posindex + 2 ) = -Pby/pow( dy, 2 ); // PVy

  //   // Lx - Lz
  //   p.r( 6 ) = Lx - Lz;

  //   p.H( 6, B_momindex + 1 ) = 1/dx; // Pbx
  //   p.H( 6, B_posindex + 1 ) = -Pbx/pow( dx, 2 ); // BVx
  //   p.H( 6, PV_posindex + 1 ) = Pbx/pow( dx, 2 ); // PVx
  //   p.H( 6, B_momindex + 3 ) = -1/dz; // Pbz
  //   p.H( 6, B_posindex + 3 ) = Pbz/pow( dz, 2 ); // BVz
  //   p.H( 6, PV_posindex + 3 ) = -Pbz/pow( dz, 2 ); // PVz
    
  //   return ErrCode::success;
  // }

  ErrCode InternalParticle::projectConstraint( Constraint::Type type, const FitParams& fitparams,
                                               Projection& p ) const {
    ErrCode status;
    switch ( type ) {
    case Constraint::mass:
    case Constraint::massEnergy:
      //       if( m_daughters.size()==2 &&
      //    !m_daughters.front()->hasEnergy() &&
      //    !m_daughters.back()->hasEnergy() )
      //  status |= projectMassConstraintTwoBody(fitparams,p) ;
      ///      else
      status |= projectMassConstraint( fitparams, p );
      // chisq = filterMassConstraintOnDaughters(fitpar) ;
      break;
    case Constraint::geometric:
      status |= projectGeoConstraint( fitparams, p );
      break;
    case Constraint::kinematic:
      status |= projectKineConstraint( fitparams, p );
      break;
    case Constraint::lifetime:
      status |= projectLifeTimeConstraint( fitparams, p );
      break;
    default:
      status |= ParticleBase::projectConstraint( type, fitparams, p );
    }
    return status;
  }

  ErrCode InternalParticle::projectMassConstraintTwoBody( const FitParams& fitparams, Projection& p ) const {
    // we can also apply the constraint to the daughters. that may
    // work better if the opening angle is small.

    // m^2 = ma^1 + mb^2 + 2 * (Ea*Eb - pxa*pxb - pya*pyb - pza*pzb )

    ParticleBase* d1 = daughters()[0];
    ParticleBase* d2 = daughters()[1];

    assert( d1->hasEnergy() == false && d2->hasEnergy() == false );

    double mass = pdtMass();
    double m1   = d1->pdtMass();
    double m2   = d2->pdtMass();

    int momindex1 = d1->momIndex();
    int momindex2 = d2->momIndex();

    // initialize the value
    double px1 = fitparams.par()( momindex1 + 1 );
    double py1 = fitparams.par()( momindex1 + 2 );
    double pz1 = fitparams.par()( momindex1 + 3 );

    double px2 = fitparams.par()( momindex2 + 1 );
    double py2 = fitparams.par()( momindex2 + 2 );
    double pz2 = fitparams.par()( momindex2 + 3 );

    double E1 = std::sqrt( m1 * m1 + px1 * px1 + py1 * py1 + pz1 * pz1 );
    double E2 = std::sqrt( m2 * m2 + px2 * px2 + py2 * py2 + pz2 * pz2 );

    p.r( 1 ) = m1 * m1 + m2 * m2 + 2 * ( E1 * E2 - px1 * px2 - py1 * py2 - pz1 * pz2 ) - mass * mass;

    // calculate the projection matrix
    p.H( 1, momindex1 + 1 ) = 2 * ( E2 * px1 / E1 - px2 );
    p.H( 1, momindex1 + 2 ) = 2 * ( E2 * py1 / E1 - py2 );
    p.H( 1, momindex1 + 3 ) = 2 * ( E2 * pz1 / E1 - pz2 );
    p.H( 1, momindex2 + 1 ) = 2 * ( E1 * px2 / E2 - px1 );
    p.H( 1, momindex2 + 2 ) = 2 * ( E1 * py2 / E2 - py1 );
    p.H( 1, momindex2 + 3 ) = 2 * ( E1 * pz2 / E2 - pz1 );

    // set the variance in the residual
    double width    = pdtWidth();
    p.Vfast( 1, 1 ) = 4 * mass * mass * width * width;

    return ErrCode::success;
  }

  void InternalParticle::addToConstraintList( constraintlist& alist, int depth ) const {
    // first the daughters
    for ( daucontainer::const_iterator it = daughters().begin(); it != daughters().end(); ++it )
      ( *it )->addToConstraintList( alist, depth - 1 );

    // double geoprecision  = 1e-5 ; // 1mu
    // double massprecision = 4*pdtMass()*pdtMass()*1e-5 ; // 0.01 MeV

    // the lifetime constraint
    if ( lenIndex() >= 0 && m_lifetimeconstraint )
      alist.push_back( Constraint( this, Constraint::lifetime, depth, 1 ) );
    // the kinematic constraint
    if ( momIndex() >= 0 ) alist.push_back( Constraint( this, Constraint::kinematic, depth, 4 ) );
    // the geometric constraint
    if ( mother() && lenIndex() >= 0 ) alist.push_back( Constraint( this, Constraint::geometric, depth, 3, 3 ) );
    // the scalefactor constraint (apply it only if it is a B+ meson)
    //if ( (particle().particleID().pid() == 521) || (particle().particleID().pid() == -521) ) alist.push_back( Constraint( this, Constraint::scalefactor, depth, 6 ) );
    // the mass constraint. FIXME: move to ParticleBase
    if ( hasMassConstraint() ) alist.push_back( Constraint( this, Constraint::mass, depth, 1, 3 ) );
  }

  std::string InternalParticle::parname( int thisindex ) const {
    int id = thisindex;
    // skip the lifetime parameter if there is no mother
    if ( !mother() && id >= 3 ) ++id;
    return ParticleBase::parname( id );
  }

  //   struct printfunctor : public unary_function<ParticleBase*,void>
  //   {
  //     printfunctor(const FitParams* fitpar) : _arg(fitpar)  {}
  //     void operator() (const ParticleBase* x) { x->print(_arg) ; }
  //     const FitParams* _arg ;
  //   };

  //   bool InternalParticle::swapMotherDaughter(FitParams* fitparams,
  //          const ParticleBase* newmother)
  //   {
  //     // routine that switches momentum vectors in a vertex, used for
  //     // reconstructing the tagging vertex.
  //     assert((newmother->type()==kBtaComposite||newmother->type()==kBtaResonance)) ;
  //     daucontainer::iterator it = std::find(m_daughters.begin(),m_daughters.end(),newmother) ;
  //     assert( it != m_daughters.end() ) ;

  //     // now start substituting
  //     // 1. assign the missing particle
  //     // 2.
  //     // 3. swap the momenta around

  //     int dummy = newmother->index() ;
  //     ParticleBase* missingparticle = new MissingParticle(0,this) ;
  //     missingparticle->updateIndex(dummy) ;

  //     // invert tau
  //     if( newmother->lenIndex()>=0 && tauIndex()>=0) {
  //       double tau = fitparams->par()(newmother->lenIndex()+1) ;
  //       fitparams->par()(lenIndex()+1) = -tau ;
  //     }

  //     // set the btacandidate
  //     setCandidate( newmother->bc() ) ;

  //     // do the momentum
  //     int momindex = momIndex() ;
  //     int momindexmother = newmother->momIndex() ;
  //     int momindexmissing = missingparticle->momIndex() ;

  //     int maxrow = newmother->hasEnergy() && hasEnergy() ? 4 : 3 ;
  //     double energy2(0) ;
  //     for( int row=1; row<=maxrow; ++row) {
  //       double pxin  = fitparams->par()(momindexmother+row) ;
  //       double pxout = fitparams->par()(momindex      +row) ;
  //       // the new missing particle
  //       fitparams->par()(momindexmissing+row) = 2*pxin - pxout ;
  //       fitparams->par()(momindex+row) = pxin ;
  //       energy2 += pxin*pxin ;
  //     }

  //     if( newmother->hasEnergy() && hasEnergy() ) {
  //       double mass = newmother->pdtMass() ;
  //       double Ein  = sqrt(energy2+mass*mass) ;
  //       double Eout = fitparams->par()(momindex+4) ;
  //       fitparams->par()(momindexmissing+4) = 2*Ein - Eout ;
  //       fitparams->par()(momindex+4)        = Ein ;
  //     }

  //     ParticleBase* newmothercopy = *it ;
  //     *it = missingparticle ;
  //     delete newmothercopy ;

  //     return true ;
  //   }

} // namespace DecayTreeFitter
