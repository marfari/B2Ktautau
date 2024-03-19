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
// $Id: Fitter.cpp,v 1.16 2010/06/09 17:55:46 ibelyaev Exp $
// ============================================================================

// STL
#include <iomanip>
#include <sstream>
#include <stdio.h>

// Boost
#include <boost/foreach.hpp>

// gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IAlgContextSvc.h"
#include "GaudiKernel/IAlgorithm.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PhysicalConstants.h"

// Event model
#include "Event/Particle.h"

#include "DecayTreeFitter/Fitter.h"

#include "Configuration.h"
#include "DecayChain.h"
#include "FitParams.h"
#include "ParticleBase.h"

namespace DecayTreeFitter {

  extern int vtxverbose;
  // ==========================================================================
  // constructor from the particle (decay head)
  // ==========================================================================
  Fitter::Fitter( const LHCb::Particle& bc, const bool forceFitAll, const ITrackStateProvider* extrapolator,
                  const bool useTrackTraj )
      : m_particle( &bc )
      , m_decaychain( 0 )
      , m_fitparams( 0 )
      , m_status( UnFitted )
      , m_chiSquare( -1 )
      , m_niter( -1 )
      , m_errCode( 0 )
      , m_extrapolator( const_cast<ITrackStateProvider*>( extrapolator ) ) {
    Configuration config( forceFitAll, extrapolator, useTrackTraj );
    // build the tree
    m_decaychain = new DecayChain( bc, config );
    // allocate the fit parameters
    m_fitparams = new FitParams( m_decaychain->dim() );
  }
  // ==========================================================================
  // constructor from the particle (decay head)
  // ==========================================================================
  Fitter::Fitter( const LHCb::Particle& bc, const ITrackStateProvider* extrapolator, const bool forceFitAll,
                  const bool useTrackTraj )
      : m_particle( &bc )
      , m_decaychain( 0 )
      , m_fitparams( 0 )
      , m_status( UnFitted )
      , m_chiSquare( -1 )
      , m_niter( -1 )
      , m_errCode( 0 )
      , m_extrapolator( const_cast<ITrackStateProvider*>( extrapolator ) ) {
    Configuration config( forceFitAll, extrapolator, useTrackTraj );
    // build the tree
    m_decaychain = new DecayChain( bc, config );
    // allocate the fit parameters
    m_fitparams = new FitParams( m_decaychain->dim() );
  }
  // ==========================================================================
  // constructor from the particle (decay head) and primary vertex
  // ==========================================================================
  Fitter::Fitter( const LHCb::Particle& bc, const LHCb::VertexBase& pv, const bool forceFitAll,
                  const ITrackStateProvider* extrapolator, const bool useTrackTraj )
      : m_particle( &bc )
      , m_decaychain( 0 )
      , m_fitparams( 0 )
      , m_status( UnFitted )
      , m_chiSquare( -1 )
      , m_niter( -1 )
      , m_errCode( 0 )
      //
      , m_extrapolator( const_cast<ITrackStateProvider*>( extrapolator ) )
  //
  {
    Configuration config( forceFitAll, extrapolator, useTrackTraj );
    m_decaychain = new DecayChain( bc, pv, config );
    m_fitparams  = new FitParams( m_decaychain->dim() );
  }
  // ==========================================================================
  // constructor from the particle (decay head) and primary vertex
  // ==========================================================================
  Fitter::Fitter( const LHCb::Particle& bc, const LHCb::VertexBase& pv, const ITrackStateProvider* extrapolator,
                  const bool forceFitAll, const bool useTrackTraj )
      : m_particle( &bc )
      , m_decaychain( 0 )
      , m_fitparams( 0 )
      , m_status( UnFitted )
      , m_chiSquare( -1 )
      , m_niter( -1 )
      , m_errCode( 0 )
      //
      , m_extrapolator( const_cast<ITrackStateProvider*>( extrapolator ) )
  //
  {
    Configuration config( forceFitAll, extrapolator, useTrackTraj );
    m_decaychain = new DecayChain( bc, pv, config );
    m_fitparams  = new FitParams( m_decaychain->dim() );
  }
  // ==========================================================================
  // set the track extrapolator
  // ==========================================================================
  void Fitter::setStateProvider( const ITrackStateProvider* extrapolator ) {
    m_extrapolator = const_cast<ITrackStateProvider*>( extrapolator );
  }

  void Fitter::setVerbose( int i ) { vtxverbose = i; }

  Fitter::~Fitter() {
    delete m_decaychain;
    delete m_fitparams;
  }
  void Fitter::fit( int nitermax, double reldChisqConv, int maxndiverging, double dChisqQuit) {
    std::vector<float> chisq_iters = {};
    std::vector<Gaudi::XYZPoint> BV = {};
    std::vector<Gaudi::XYZPoint> DV1 = {};
    std::vector<Gaudi::XYZPoint> DV2 = {};
    std::vector<float> B_M = {};
    std::vector<Gaudi::XYZVector> p_taup_nu = {}, p_taum_nu = {};
    std::vector<Gaudi::XYZVector> p_taup = {}, p_taum = {};
    const LHCb::Particle* P = 0;

    fit(nitermax, reldChisqConv, maxndiverging, dChisqQuit, chisq_iters,
        BV, DV1, DV2, P, B_M,
        p_taup_nu, p_taum_nu,
        p_taup, p_taum, false);
  }
  void Fitter::fit( int nitermax, double reldChisqConv, int maxndiverging, double dChisqQuit, std::vector<float>& chisq_iters,
                    std::vector<Gaudi::XYZPoint>& BV, std::vector<Gaudi::XYZPoint>& DV1, std::vector<Gaudi::XYZPoint>& DV2,
                    const LHCb::Particle* P, std::vector<float> B_M,
                    std::vector<Gaudi::XYZVector>& p_taup_nu, std::vector<Gaudi::XYZVector>& p_taum_nu,
                    std::vector<Gaudi::XYZVector>& p_taup, std::vector<Gaudi::XYZVector>& p_taum,
                    bool fillArrays) 
  {
    m_map.clear();
    const ParticleBase* pb = 0 != P ? m_decaychain->locate( *P ) : m_decaychain->cand();

    int posindex = pb->posIndex() - 2;
    //const int maxndiverging = 3;
    // const double dChisqQuit = nDof() ; // if chi2 increases by more than this --> fit failed

    // initialize
    m_chiSquare = -1;
    // m_errCode.reset() ;
    m_errCode = 0;

    if ( m_status == UnFitted ) { m_errCode = m_decaychain->init( *m_fitparams ).flag(); }

    // std::cout<<"After initialization printing params "<<std::endl;
    // m_fitparams->print();
    
    int posindex_Taup = 0;
    int posindex_Taum = 0;
    if(fillArrays)
    {
      //Get Kaon
      const LHCb::Particle * Kplus = LHCb::DecayTree::findFirstInTree( P, LHCb::ParticleID(321) );
      if ( Kplus == 0 ) 
      {
        Kplus = LHCb::DecayTree::findFirstInTree( P, LHCb::ParticleID(-321) );
      }
      //Get taus
      LHCb::Particle::ConstVector tauList;
      //If Kaon is K+, then we want Taup to have ID = -15, Taum to be +15
      //Otherwise, Taup should be the +15 and Taum should be -15

      if((Kplus->particleID()).pid() == 321)
      {
        LHCb::DecayTree::findInTree( P, LHCb::ParticleID(-15), tauList); // -15, 411
        LHCb::DecayTree::findInTree( P, LHCb::ParticleID(15), tauList); // 15, -411
      }
      else if((Kplus->particleID()).pid() == -321)
      {
        LHCb::DecayTree::findInTree( P, LHCb::ParticleID(15), tauList); // 15, -411
        LHCb::DecayTree::findInTree( P, LHCb::ParticleID(-15), tauList); // -15, 411
      }

      const ParticleBase* pb_Taup = m_decaychain->locate( *(tauList[0]) );
      const ParticleBase* pb_Taum = m_decaychain->locate( *(tauList[1]) );

      posindex_Taup = pb_Taup->posIndex(); 
      posindex_Taum = pb_Taum->posIndex();
      //NB: The numbering in the fitparams array starts from 1
      BV.push_back(Gaudi::XYZPoint( m_fitparams->par()(posindex+1), m_fitparams->par()(posindex+2), m_fitparams->par()(posindex+3) ));
      DV1.push_back(Gaudi::XYZPoint( m_fitparams->par()(posindex_Taup+1), m_fitparams->par()(posindex_Taup+2), m_fitparams->par()(posindex_Taup+3) ));
      DV2.push_back(Gaudi::XYZPoint( m_fitparams->par()(posindex_Taum+1), m_fitparams->par()(posindex_Taum+2), m_fitparams->par()(posindex_Taum+3) ));
      p_taup_nu.push_back(Gaudi::XYZVector(m_fitparams->par()(12+1), m_fitparams->par()(13+1), m_fitparams->par()(14+1)));
      p_taum_nu.push_back(Gaudi::XYZVector(m_fitparams->par()(32+1), m_fitparams->par()(33+1), m_fitparams->par()(34+1)));
      p_taup.push_back(Gaudi::XYZVector(m_fitparams->par()(19+1), m_fitparams->par()(20+1), m_fitparams->par()(21+1)));
      p_taum.push_back(Gaudi::XYZVector(m_fitparams->par()(39+1), m_fitparams->par()(40+1), m_fitparams->par()(41+1)));
      // std::cout<<"Pushing back Taup_nu: "<<m_fitparams->par()(12+1)<<" "<<m_fitparams->par()(13+1)<<" "<<m_fitparams->par()(14+1)<<std::endl;
      // std::cout<<"Pushing back Taup_ENDVTX: "<<m_fitparams->par()(15+1)<<" "<<m_fitparams->par()(16+1)<<" "<<m_fitparams->par()(17+1)<<std::endl;
      // std::cout<<"Pushing back Taum_nu: "<<m_fitparams->par()(32+1)<<" "<<m_fitparams->par()(33+1)<<" "<<m_fitparams->par()(34+1)<<std::endl;
      // std::cout<<"Pushing back Taum_ENDVTX: "<<m_fitparams->par()(35+1)<<" "<<m_fitparams->par()(36+1)<<" "<<m_fitparams->par()(37+1)<<std::endl;
      //B_M.push_back(fitParams(*pb).momentum().m().value());
      // BV_X.push_back(m_fitparams->par()(posindex+1));
      // BV_Y.push_back(m_fitparams->par()(posindex+2));
      // BV_Z.push_back(m_fitparams->par()(posindex+3));
      //std::cout<<"After initialization pushing back BV_X = "<<BV_X.back()<<" BV_Y = "<<BV_Y.back()<<" BV_Z = "<<BV_Z.back()<<std::endl;

      //m_fitparams->chiSquare() here is 0
      chisq_iters.push_back(m_fitparams->chiSquare()); //Pushing back anyway to have a same length array as the vertices
      //std::cout<<"After initialization Pushing back "<<m_fitparams->chiSquare()<<std::endl;
    }

    
    // if(m_errCode.failure()) {
    if ( 0 != m_errCode ) {
      // the input tracks are too far apart
      m_status = BadInput;

    } else {
      // reset the status flag
      m_status = UnFitted;

      int    ndiverging        = 0;
      bool   finished          = false;
      double prevabsdeltachisq = -1;

      for ( m_niter = 0; m_niter < nitermax && !finished; ++m_niter ) {

        CLHEP::HepVector prevpar   = m_fitparams->par();
        bool             firstpass = m_niter == 0;
        m_errCode                  = m_decaychain->filter( *m_fitparams, firstpass ).flag();
        double chisq               = m_fitparams->chiSquare();
        double deltachisq          = chisq - m_chiSquare;
        // if chi2 increases by more than this --> fit failed
        int          ndof       = nDof();
        // std::cout << "ndof = " << ndof << std::endl;
        const double dChisqQuit = std::max( double( 2 * ndof ), 2 * m_chiSquare );
        const double dChisqConv = reldChisqConv * std::max( double( ndof ), std::min( chisq, m_chiSquare ) );

        if ( 0 != m_errCode ) {
          finished = true;
          m_status = Failed;
        } else {
          if ( m_niter > 0 ) {
            if ( std::abs( deltachisq ) < dChisqConv ) {
              m_chiSquare = chisq;
              m_status    = Success;
              finished    = true;
            } else if ( m_niter > 1 && deltachisq > dChisqQuit ) {
              //m_fitparams->par() = prevpar;
              m_status           = Failed;
              m_errCode          = ErrCode::fastdivergingfit;
              finished           = true;
            } else if ( m_niter > 2 && std::abs( deltachisq ) > prevabsdeltachisq && ++ndiverging >= maxndiverging ) {
              //m_fitparams->par() = prevpar;
              m_status           = NonConverged;
              m_errCode          = ErrCode::slowdivergingfit;
              finished           = true;
            } else if ( deltachisq > 0 ) {
              // make a half step and reduce stepsize to prevent oscillations
              m_fitparams->par() += prevpar;
              m_fitparams->par() *= 0.5;
              // if( m_niter>10) m_fitparams->scale() *= 0.5 ;
            }
          }
          if ( deltachisq < 0 )
          {
            ndiverging = 0; // start over with counting
          } 
          if ( !finished ) {
            m_chiSquare       = chisq;
            prevabsdeltachisq = std::abs( deltachisq );
          }
        }

        if ( vtxverbose >= 1 ) {
          std::ostringstream mess;
          mess << "DecayTreeFitter::Fitter: step, chiSquare: " << std::setw( 3 ) << m_niter << std::setw( 3 )
               << m_status << std::setw( 3 ) << nDof() << std::setprecision( 6 ) << std::setw( 12 ) << chisq
               << std::setw( 12 ) << deltachisq;
          print( mess.str(), MSG::VERBOSE );
        }

        if ( vtxverbose >= 4 ) {
          fillStream( std::cout );
          std::cout << "press a key ...." << std::endl;
          getchar();
        }
        
        if(fillArrays)
        {
          chisq_iters.push_back(chisq);
          BV.push_back(Gaudi::XYZPoint( m_fitparams->par()(posindex+1), m_fitparams->par()(posindex+2), m_fitparams->par()(posindex+3) ));
          DV1.push_back(Gaudi::XYZPoint( m_fitparams->par()(posindex_Taup+1), m_fitparams->par()(posindex_Taup+2), m_fitparams->par()(posindex_Taup+3) ));
          DV2.push_back(Gaudi::XYZPoint( m_fitparams->par()(posindex_Taum+1), m_fitparams->par()(posindex_Taum+2), m_fitparams->par()(posindex_Taum+3) ));
          p_taup_nu.push_back(Gaudi::XYZVector(m_fitparams->par()(12+1), m_fitparams->par()(13+1), m_fitparams->par()(14+1)));
          p_taum_nu.push_back(Gaudi::XYZVector(m_fitparams->par()(32+1), m_fitparams->par()(33+1), m_fitparams->par()(34+1)));
          p_taup.push_back(Gaudi::XYZVector(m_fitparams->par()(19+1), m_fitparams->par()(20+1), m_fitparams->par()(21+1)));
          p_taum.push_back(Gaudi::XYZVector(m_fitparams->par()(39+1), m_fitparams->par()(40+1), m_fitparams->par()(41+1)));
          //std::cout<<"Pushing back Taup_nu: "<<m_fitparams->par()(12)<<" "<<m_fitparams->par()(13)<<" "<<m_fitparams->par()(14)<<std::endl;
          //std::cout<<"Iteration = "<<m_niter<<" ChiSquare = "<<chisq<<std::endl; -> the values of chisq for each iteration
        }

        // std::cout<<"Iteration "<<m_niter+1<<" printing params "<<std::endl;
        // m_fitparams->print();
        //B_M.push_back(fitParams(*pb).momentum().m().value());
        // BV_X.push_back(m_fitparams->par()(posindex+1));
        // BV_Y.push_back(m_fitparams->par()(posindex+2));
        // BV_Z.push_back(m_fitparams->par()(posindex+3));
        // std::cout<<"Iteration "<<m_niter<<" Pushing back "<<chisq<<std::endl;
        // std::cout<<"Iteration "<<m_niter<<" pushing back BV_X = "<<BV_X.back()<<" BV_Y = "<<BV_Y.back()<<" BV_Z = "<<BV_Z.back()<<std::endl;

      }
    
      if ( m_niter == nitermax && m_status != Success ) { 
        m_status = NonConverged; 
      }

      // m_decaychain->mother()->forceP4Sum(*m_fitparams) ;

      if ( !( m_fitparams->testCov() ) ) {
        static int counter( 10 );
        if ( --counter >= 0 ) {
          print( "DecayTreeFitter::Fitter: Warning: Error matrix not positive definite."
                 " Changing status to failed.",
                 MSG::DEBUG );
        }
        m_status = Failed;
      }

    }
  }

  void Fitter::fitOneStep() {
    bool firstpass = m_status == UnFitted;
    if ( firstpass ) m_decaychain->init( *m_fitparams );
    m_decaychain->filter( *m_fitparams, firstpass );
    m_chiSquare = m_fitparams->chiSquare();
    if ( vtxverbose >= 1 ) {
      std::ostringstream mess;
      mess << "DecayTreeFitter::Fitter::fitOneStep(): " << m_status << " " << firstpass << " " << m_chiSquare;
      print( mess.str(), MSG::VERBOSE );
    }
    m_status = Success;
  }

  std::ostream& Fitter::fillStream( std::ostream& os ) const {
    // std::ostringstream s ;
    os << "****** chisq,ndof,ncontraints,niter,status,errcode: " << chiSquare() << " " << nDof() << " "
       << m_fitparams->nConstraints() << " " << nIter() << " " << status() << " " << m_errCode << std::endl;
    m_decaychain->mother()->fillStream( os, m_fitparams );
    os << "***************************************" << std::endl;
    return os;
  }

  int Fitter::nDof() const { return m_fitparams->nDof(); }

  double Fitter::add( const LHCb::Particle& cand ) {
    // first obtain a map
    // ParticleBase::indexmap indexmap ;
    // m_decaychain->mother()->addToMap( indexmap ) ;
    // add this track

    ParticleBase* bp     = m_decaychain->mother()->addDaughter( cand, Configuration() );
    int           offset = m_fitparams->dim();
    double        deltachisq( 999 );
    if ( bp ) {
      bp->updateIndex( offset );
      // reassign indices
      // int offset(0) ;
      // m_decaychain->mother()->updIndex(offset) ;
      // resize the fitparams
      m_fitparams->resize( offset );
      // initialize the added track, filter and return the delta chisq
      bp->initPar1( m_fitparams );
      bp->initPar2( m_fitparams );
      bp->initCov( m_fitparams );

      ParticleBase::constraintlist constraints;
      bp->addToConstraintList( constraints, 0 );
      double  chisq = m_fitparams->chiSquare();
      ErrCode status;
      for ( ParticleBase::constraintlist::const_iterator it = constraints.begin(); it != constraints.end(); ++it )
        status |= it->filter( *m_fitparams );

      deltachisq  = m_fitparams->chiSquare() - chisq;
      m_chiSquare = m_fitparams->chiSquare();

      // we want this somewhere else, but too much work now
      decaychain()->initConstraintList();

    } else {
      std::ostringstream mess;
      mess << "DecayTreeFitter::Fitter: WARNING: Cannot add track to this vertex ..." << m_decaychain->mother()->type();
      print( mess.str(), MSG::WARNING );
    }
    return deltachisq;
  }

  double Fitter::remove( const LHCb::Particle& cand ) {
    ParticleBase* pb = const_cast<ParticleBase*>( m_decaychain->locate( cand ) );
    ErrCode       status;
    double        dchisq( 0 );
    if ( pb ) {
      // filter with negative weight
      ParticleBase::constraintlist constraints;
      pb->addToConstraintList( constraints, 0 );
      double chisq = m_fitparams->chiSquare();
      for ( ParticleBase::constraintlist::iterator it = constraints.begin(); it != constraints.end(); ++it ) {
        it->setWeight( -1 );
        status |= it->filter( *m_fitparams );
      }
      dchisq = chisq - m_fitparams->chiSquare();
      // remove
      ParticleBase* themother = const_cast<ParticleBase*>( pb->mother() );
      if ( themother ) themother->removeDaughter( pb );
    }
    return dchisq;
  }

  void Fitter::setMassConstraint( const LHCb::Particle& bc, bool add ) {
    m_decaychain->setMassConstraint( bc, add );
    m_status = UnFitted;
  }

  void Fitter::setMassConstraint( const LHCb::Particle& bc, double mass ) {
    m_decaychain->setMassConstraint( bc, mass );
    m_status = UnFitted;
  }

  void Fitter::setMassConstraint( const LHCb::ParticleID& pid, bool add ) {
    m_decaychain->setMassConstraint( pid, add );
    m_status = UnFitted;
  }
  void Fitter::setMassConstraint( const LHCb::ParticleID& pid, double mass ) {
    m_decaychain->setMassConstraint( pid, mass );
    m_status = UnFitted;
  }

  void Fitter::updateIndex() {
    int offset = 0;
    m_decaychain->mother()->updateIndex( offset );
    m_fitparams->resize( offset );
  }

  double Fitter::globalChiSquare() const { return m_decaychain->chiSquare( m_fitparams ); }

  Gaudi::Math::ParticleParams Fitter::fitParams( const ParticleBase& pb ) const {
    int posindex = pb.posIndex();
    // hack: for tracks and photons, use the production vertex
    if ( posindex < 0 && pb.mother() ) posindex = pb.mother()->posIndex();
    int                  momindex = pb.momIndex();
    int                  lenindex = pb.lenIndex();
    Gaudi::XYZPoint      pos( m_fitparams->par()( posindex + 1 ), m_fitparams->par()( posindex + 2 ),
                         m_fitparams->par()( posindex + 3 ) );
    Gaudi::LorentzVector p4;
    p4.SetPx( m_fitparams->par()( momindex + 1 ) );
    p4.SetPy( m_fitparams->par()( momindex + 2 ) );
    p4.SetPz( m_fitparams->par()( momindex + 3 ) );

    Gaudi::SymMatrix8x8 cov8;
    double              decaylength = lenindex >= 0 ? m_fitparams->par()( lenindex + 1 ) : 0;

    if ( pb.hasEnergy() ) {
      // if particle has energy, get full p4 from fitparams
      p4.SetE( m_fitparams->par()( momindex + 4 ) );
      int parmap[8];
      for ( int i = 0; i < 3; ++i ) parmap[i] = posindex + i;
      for ( int i = 0; i < 4; ++i ) parmap[i + 3] = momindex + i;
      parmap[7]  = lenindex;
      int maxrow = lenindex >= 0 ? 8 : 7;
      for ( int row = 0; row < maxrow; ++row )
        for ( int col = 0; col <= row; ++col )
          cov8( row, col ) = m_fitparams->cov()( parmap[row] + 1, parmap[col] + 1 );
    } else {
      // if not, use the pdttable mass
      Gaudi::SymMatrix7x7 cov7;
      int                 parmap[8];
      for ( int i = 0; i < 3; ++i ) parmap[i] = posindex + i;
      for ( int i = 0; i < 3; ++i ) parmap[i + 3] = momindex + i;
      parmap[6]  = lenindex;
      int maxrow = lenindex >= 0 ? 7 : 6;
      for ( int row = 0; row < maxrow; ++row )
        for ( int col = 0; col <= row; ++col )
          cov7( row, col ) = m_fitparams->cov()( parmap[row] + 1, parmap[col] + 1 );

      // now fill the jacobian
      double mass    = pb.pdtMass();
      double energy2 = mass * mass;
      for ( int row = 0; row < 3; ++row ) {
        double px = m_fitparams->par()( momindex + row + 1 );
        energy2 += px * px;
      }
      double energy = std::sqrt( energy2 );

      ROOT::Math::SMatrix<double, 8, 7> jacobian;
      for ( int col = 0; col < 7; ++col ) jacobian( col, col ) = 1;
      for ( int col = 3; col < 6; ++col ) jacobian( 6, col ) = m_fitparams->par()( parmap[col] + 1 ) / energy;

      p4.SetE( energy );
      cov8 = ROOT::Math::Similarity( jacobian, cov7 );
    }
    
    //
    // VtxFitParams vtxfitparams(pb.charge(),pos,p4,decaylength,cov8) ;
    // return vtxfitparams ;
    return Gaudi::Math::ParticleParams( pos, p4, decaylength, cov8 );
  }

  const Gaudi::Math::ParticleParams* Fitter::fitParams( const LHCb::Particle* cand ) const {
    Map::const_iterator ifind = m_map.find( cand );
    if ( m_map.end() != ifind ) { return &ifind->second; } // RETURN
    //
    const ParticleBase* pb = 0 != cand ? m_decaychain->locate( *cand ) : m_decaychain->cand();
    //
    if ( 0 == pb ) { return NULL; } // RETURN
    //
    return &( m_map.insert( Map::value_type( cand, fitParams( *pb ) ) ).first->second );
  }

  Gaudi::Math::ValueWithError Fitter::decayLengthSum( const LHCb::Particle& candA, const LHCb::Particle& candB ) const {
    // returns the decaylengthsum of two particles (use ful for e.g. B->DD)
    Gaudi::Math::ValueWithError rc( 0, 0 );
    const ParticleBase*         pbA = m_decaychain->locate( candA );
    const ParticleBase*         pbB = m_decaychain->locate( candB );
    if ( pbA && pbB && pbA->mother() && pbB->mother() ) {
      int lenindexA = pbA->lenIndex();
      int lenindexB = pbB->lenIndex();
      if ( lenindexA >= 0 && lenindexB >= 0 ) {
        double lenA = m_fitparams->par()( lenindexA + 1 );
        double lenB = m_fitparams->par()( lenindexB + 1 );
        double cov  = m_fitparams->cov().fast( lenindexA + 1, lenindexA + 1 ) +
                     m_fitparams->cov().fast( lenindexB + 1, lenindexB + 1 ) +
                     2 * m_fitparams->cov().fast( lenindexA + 1, lenindexB + 1 );

        // rc = VtxDoubleErr(lenA+lenB,std::sqrt(cov)) ;
        rc = Gaudi::Math::ValueWithError( lenA + lenB, cov );
      }
    }
    return rc;
  }

  std::string Fitter::name( const LHCb::Particle& cand ) const {
    const ParticleBase* pb = m_decaychain->locate( cand );
    return pb ? pb->name() : "Not found";
  }

  LHCb::Particle Fitter::getFitted() const {
    LHCb::Particle thecand = *particle();
    // const ParticleBase* pb = m_decaychain->locate(thecand) ;
    const ParticleBase* pb = m_decaychain->locate( *particle() );
    if ( pb )
      updateCand( *pb, thecand );
    else {
      print( "DecayTreeFitter::Fitter: ERROR: Cannot find particle in tree", MSG::ERROR );
    }
    return thecand;
  }

  LHCb::Particle Fitter::getFitted( const LHCb::Particle& cand ) const {
    // clone the particle
    LHCb::Particle thecand = cand;
    // find the ParticleBase corresponding to the originaERROR: l particle.
    const ParticleBase* pb = m_decaychain->locate( cand );
    if ( pb )
      updateCand( *pb, thecand );
    else {
      print( "DecayTreeFitter::Fitter: ERROR: Cannot find particle in tree", MSG::ERROR );
    }
    return thecand;
  }

  //   LHCb::Particle
  //   Fitter::getFitted(const LHCb::Particle& cand) const
  //   {
  //     LHCb::Particle thecand = cand ;
  //     updateCand(thecand) ;
  //     return thecand ;
  //   }

  LHCb::Particle* Fitter::fittedCand( const LHCb::Particle& /*cand*/, LHCb::Particle* /*headOfTree*/ ) const {
    print( "DecayTreeFitter::Fitter::fittedCand: WARNING: Not yet implemented", MSG::WARNING );
    return 0;
    // assigns fitted parameters to candidate in tree
    // LHCb::Particle* acand = const_cast<LHCb::Particle*>(headOfTree->cloneInTree(cand)) ;
    // updateCand(*acand) ;
    // return acand ;
  }

  LHCb::DecayTree Fitter::getFittedTree() const {
    // clone the decay tree
    LHCb::DecayTree tree( *m_particle );
    // update the tree.
    // the easy version will only work once we have a proper 'isCloneOf'.
    // updateTree( *tree.head() ) ;
    for ( LHCb::DecayTree::CloneMap::const_iterator it = tree.cloneMap().begin(); it != tree.cloneMap().end(); ++it ) {
      const ParticleBase* pb = m_decaychain->locate( *( it->first ) );
      if ( pb != 0 ) updateCand( *pb, *( it->second ) );
    }
    return tree;
  }

  /*
    std::string myvtxprint(const LHCb::Particle & cand,
    const ComIOManip & format) {
    std::ostringstream stream;
    const PdtEntry * pdtEntry = cand.pdtEntry();
    if (0 != pdtEntry){
    stream << pdtEntry->name();
    if(!cand.isComposite())
    stream << "(" << cand.uid() << ")" ;
    stream << std::ends;
    }
    HepString result(stream.str());
    //delete [] stream.str();           // must take ownership here
    return result;
    }
  */

  void Fitter::updateCand( const ParticleBase& pb, LHCb::Particle& particle ) const {
    // assigns fitted parameters to a candidate
    // VtxFitParams vtxpar = fitParams(pb) ;
    Gaudi::Math::ParticleParams vtxpar = fitParams( pb );

    // update everything inside the particle. don't update the vertex, for now.
    particle.setMomentum( vtxpar.momentum() ); // default
    particle.setReferencePoint( vtxpar.position() );
    particle.setMomCovMatrix( vtxpar.momCovMatrix() );
    particle.setPosCovMatrix( vtxpar.posCovMatrix() );
    particle.setPosMomCovMatrix( vtxpar.momPosCov() );
    //
    // update the chi2 of global fit if this is the head of the tree
    if ( &pb == m_decaychain->cand() ) {
      particle.eraseInfo( LHCb::Particle::Chi2OfParticleReFitter );
      particle.addInfo( LHCb::Particle::Chi2OfParticleReFitter, chiSquare() );
    }
    // update the vertex as well, if this is the head of the tree
    if ( &pb == m_decaychain->cand() ) {
      LHCb::Vertex* vertex = particle.endVertex();
      if ( !vertex ) {
        // create a new vertex
        vertex = new LHCb::Vertex();
        // don't add any daughters. I don't understand why we need them anyway.
        particle.setEndVertex( vertex );
      }
      vertex->setTechnique( LHCb::Vertex::LastGlobal );
      vertex->setChi2AndDoF( chiSquare(), nDof() );
      vertex->setPosition( vtxpar.position() );
      vertex->setCovMatrix( vtxpar.posCovMatrix() );
    }

    //  if( pb ==m_decaychain->cand() ) {
    //    vtx = new VtxVertex(chiSquare(),nDof(),vtxpar.pos(),vtxpar.posCov(),vtxpar.xp4Cov()) ;
    //    vtx->setStatus(FitStatus::VtxStatus(status())) ;
    //    vtx->setType(FitStatus::Geometric) ;
    //  } else {
    //    // all updated daughters are reset to unfitted, but leave the chisquare
    //    double chisq = cand.decayVtx() ? cand.decayVtx()->chiSquared() : 0 ;
    //    int ndof     = cand.decayVtx() ? cand.decayVtx()->nDof() : 0 ;
    //    vtx = new VtxVertex(chisq,ndof,vtxpar.pos(),vtxpar.posCov(),vtxpar.xp4Cov()) ;
    //    vtx->setStatus(FitStatus::UnFitted) ;
    //  }
    //       }
    //      cand.setTrajectory(vtxpar,vtx) ;
  }

  bool Fitter::updateCand( LHCb::Particle& particle ) const {
    // assigns fitted parameters to a candidate
    const ParticleBase* pb = m_decaychain->locate( particle );
    if ( pb ) updateCand( *pb, particle );
    // else {
    // this error message does not make sense, because it is also
    // triggered for daughters that were simply not refitted. we
    // have to do something about that.
    //       VtxPrintTree printer(&myvtxprint) ;
    //       ErrMsg(error) << "cann't find candidate " << std::endl
    //       << printer.print(cand)
    //       << "in tree " << std::endl
    //       << printer.print(*_bc)
    //       << endmsg;
    return pb != 0;
  }

  bool Fitter::updateTree( LHCb::Particle& p ) const {
    bool rc;
    if ( ( rc = updateCand( p ) ) ) {
      BOOST_FOREACH ( const LHCb::Particle* daughter, p.daughters() ) {
        updateTree( const_cast<LHCb::Particle&>( *daughter ) );
      }
    }
    return rc;
  }

  // Get the chisquare of a particular particle in the decay chain
  ChiSquare Fitter::chiSquare( const LHCb::Particle& particle ) const {
    return m_decaychain->chiSquare( particle, *m_fitparams );
  }

  IMessageSvc* Fitter::msgService() const {
    static IMessageSvc* msgSvc = NULL;
    if ( !msgSvc ) {
      ISvcLocator* svcLocator = Gaudi::svcLocator();
      if ( svcLocator ) {
        svcLocator->service( "MessageSvc", msgSvc ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      }
    }
    return msgSvc;
  }

  const IAlgContextSvc* Fitter::algContextSvc() const {
    static const IAlgContextSvc* asvc = NULL;
    if ( !asvc ) {
      ISvcLocator* svcLocator = Gaudi::svcLocator();
      if ( svcLocator ) {
        svcLocator->service( "AlgContextSvc", asvc ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      }
    }
    return asvc;
  }

  const IAlgorithm* Fitter::getAlg() const {
    // Locate the context service
    const IAlgContextSvc* asvc = algContextSvc();

    // return the current alg
    return ( asvc ? asvc->currentAlg() : NULL );
  }

  void Fitter::print( const std::string& msg, const MSG::Level level ) const {
    // Try and get current alg
    const IAlgorithm* alg = getAlg();
    // get the message svc
    IMessageSvc* msgSvc = msgService();
    if ( msgSvc && alg ) {
      // Use gaudi messaging
      MsgStream log( msgSvc, alg->name() );
      log << level << msg << endmsg;
    } else {
      // Gaudi not available, so use basic cout...
      std::cout << msg << std::endl;
    }
  }

} // namespace DecayTreeFitter

// ============================================================================
// The END
// ============================================================================
