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

// Include files
#include "gsl/gsl_sys.h"

// from Gaudi
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/Vector3DTypes.h"

// local
#include "TupleToolRadUpDownAsym.h"

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/Particle.h"

using namespace LHCb;
//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolRadUpDownAsym
//
// Original author: Cedric Potterat
// Code Maintainer: Albert Puig
// Date: 05-02-2019
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_COMPONENT( TupleToolRadUpDownAsym )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolRadUpDownAsym::TupleToolRadUpDownAsym( const std::string& type, const std::string& name,
                                                const IInterface* parent )
    : TupleToolBase( type, name, parent ) {
  declareInterface<IParticleTupleTool>( this );
}

//=============================================================================
// Fill
//=============================================================================
StatusCode TupleToolRadUpDownAsym::fill( const LHCb::Particle* top, const LHCb::Particle* part, const std::string& head,
                                         Tuples::Tuple& tuple ) {
  const std::string prefix = fullName( head );

  std::vector<double> angleValue;
  std::vector<double> pidValue;

  angleValue.push_back( -999 );
  angleValue.push_back( -999 );
  angleValue.push_back( -999 );
  pidValue.push_back( -999 );
  pidValue.push_back( -999 );
  pidValue.push_back( -999 );

  bool test = true;
  test      = UpDownTheta( top, part, angleValue, pidValue );
  if ( !test ) return StatusCode::SUCCESS;
  // fill the tuple:
  test &= tuple->column( prefix + "_UpDownCosTheta12", angleValue[0] );
  test &= tuple->column( prefix + "_UpDownCosTheta13", angleValue[1] );
  test &= tuple->column( prefix + "_UpDownCosTheta23", angleValue[2] );
  test &= tuple->column( prefix + "_UpDownID1", pidValue[0] );
  test &= tuple->column( prefix + "_UpDownID2", pidValue[1] );
  test &= tuple->column( prefix + "_UpDownID3", pidValue[2] );

  return StatusCode( test );
}
//=============================================================================
// Fill
//=============================================================================

bool TupleToolRadUpDownAsym::isEndTree( const LHCb::Particle* p ) const {
  if ( p->isBasicParticle() ) return true;
  // converted photons, pi0 and V0s are considered as basic particles
  if ( (int)p->particleID().abspid() == 22 ) return true;   // gamma
  if ( (int)p->particleID().abspid() == 111 ) return true;  // pi0
  if ( (int)p->particleID().abspid() == 310 ) return true;  // KS0
  if ( (int)p->particleID().abspid() == 130 ) return true;  // KL0
  if ( (int)p->particleID().abspid() == 3122 ) return true; // Lambda
  return false;
}

std::vector<const LHCb::Particle*> TupleToolRadUpDownAsym::getTree( const LHCb::Particle* P ) const {
  std::vector<const LHCb::Particle*> tree;
  if ( isEndTree( P ) ) {
    tree.push_back( P );
  } else {
    SmartRefVector<LHCb::Particle> daughters = P->daughters(); // local copy to sort
    for ( SmartRefVector<LHCb::Particle>::const_iterator idau = daughters.begin(); idau != daughters.end(); ++idau ) {
      std::vector<const LHCb::Particle*> temp = getTree( *idau );
      for ( std::vector<const LHCb::Particle*>::iterator dd = temp.begin(); temp.end() != dd; ++dd ) {
        tree.push_back( *dd );
      }
    }
  }
  return tree;
}

bool TupleToolRadUpDownAsym::UpDownTheta( const LHCb::Particle* top, const LHCb::Particle* part,
                                          std::vector<double>& angleValue, std::vector<double>& pidValue ) const {

  if ( top != part ) return false;

  SmartRefVector<LHCb::Particle> dau = part->daughters();

  if ( dau.empty() ) return false;
  if ( dau.size() != 2 ) return false;

  Gaudi::LorentzVector gamma;
  Gaudi::LorentzVector res;
  Gaudi::LorentzVector h1;
  Gaudi::LorentzVector h2;
  Gaudi::LorentzVector h3;
  double               pidh1 = -999;
  double               pidh2 = -999;
  double               pidh3 = -999;

  for ( SmartRefVector<LHCb::Particle>::const_iterator it = dau.begin(); dau.end() != it; ++it ) {
    if ( ( *it )->particleID().abspid() == 22 ) // found the gamma
      gamma = Gaudi::LorentzVector( ( *it )->momentum() );
    else { // found the resonance
      res                                              = Gaudi::LorentzVector( ( *it )->momentum() );
      const std::vector<const LHCb::Particle*> dau_res = getTree( *it );
      if ( dau_res.empty() ) return false;
      if ( dau_res.size() != 3 ) return false;
      double fac = ( ( *it )->particleID().pid() > 0 ) ? 1. : -1.;

      for ( std::vector<const LHCb::Particle*>::const_iterator it2 = dau_res.begin(); dau_res.end() != it2; ++it2 ) {
        if ( ( *it2 )->charge() * fac < 0 ) {
          h3    = Gaudi::LorentzVector( ( *it2 )->momentum() );
          pidh3 = ( *it2 )->particleID().pid();
        } else if ( pidh1 == -999 ) {
          h1    = Gaudi::LorentzVector( ( *it2 )->momentum() );
          pidh1 = ( *it2 )->particleID().pid();
        } else {
          h2    = Gaudi::LorentzVector( ( *it2 )->momentum() );
          pidh2 = ( *it2 )->particleID().pid();
        }
      }
    }
  }

  ROOT::Math::Boost boost( res.BoostToCM() );

  const Gaudi::LorentzVector boostedh1    = boost( h1 );
  const Gaudi::LorentzVector boostedh2    = boost( h2 );
  const Gaudi::LorentzVector boostedh3    = boost( h3 );
  const Gaudi::LorentzVector boostedgamma = boost( gamma );

  const Gaudi::LorentzVector res12 = Gaudi::LorentzVector( boostedh1 + boostedh2 );
  const Gaudi::LorentzVector res13 = Gaudi::LorentzVector( boostedh1 + boostedh3 );
  const Gaudi::LorentzVector res23 = Gaudi::LorentzVector( boostedh2 + boostedh3 );

  double s12 = res12.M2();
  double s13 = res13.M2();
  double s23 = res23.M2();

  // Calculate cos theta

  const Gaudi::XYZVector h1V    = Gaudi::XYZVector( boostedh1 ).unit();
  const Gaudi::XYZVector h2V    = Gaudi::XYZVector( boostedh2 ).unit();
  const Gaudi::XYZVector h3V    = Gaudi::XYZVector( boostedh3 ).unit();
  const Gaudi::XYZVector gammaV = Gaudi::XYZVector( boostedgamma ).unit();

  const Gaudi::XYZVector n12 = ( h1V.Cross( h2V ) ).unit();
  const Gaudi::XYZVector n13 = ( h1V.Cross( h3V ) ).unit();
  const Gaudi::XYZVector n23 = ( h2V.Cross( h3V ) ).unit();

  double costheta12 = n12.Dot( -gammaV );
  double costheta13 = n13.Dot( -gammaV );
  double costheta23 = n23.Dot( -gammaV );
  // Correct sign and return
  angleValue[0] = ( s13 > s23 ) ? costheta12 : -costheta12;
  angleValue[1] = ( s12 > s23 ) ? costheta13 : -costheta13;
  angleValue[2] = ( s12 > s13 ) ? costheta23 : -costheta23;

  pidValue[0] = pidh1;
  pidValue[1] = pidh2;
  pidValue[2] = pidh3;

  return true;
}
