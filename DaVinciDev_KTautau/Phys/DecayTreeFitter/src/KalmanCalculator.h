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
#ifndef VTK_KALMANCALCULATOR_HH
#define VTK_KALMANCALCULATOR_HH

#include "FitParams.h"
#include "VtxErrCode.h"
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace DecayTreeFitter {

  class KalmanCalculator {
  public:
    ErrCode init( const CLHEP::HepVector& value, const CLHEP::HepMatrix& G, const FitParams& fitparams,
                  const CLHEP::HepSymMatrix* V = 0, int weight = 1 );
    void    updatePar( FitParams& fitparams );
    void    updatePar( const CLHEP::HepVector& prediction, FitParams& fitparams );
    void    updateCov( FitParams& fitparams );
    double  chisq() const { return m_chisq; }

  private:
    int                     m_nconstraints; // dimension of the constraint
    int                     m_nparameters;  // dimension of the state
    const CLHEP::HepVector* m_value;
    const CLHEP::HepMatrix* m_matrixG;
    CLHEP::HepSymMatrix     m_matrixR;    // cov of residual
    CLHEP::HepSymMatrix     m_matrixRinv; // inverse of cov of residual
    CLHEP::HepMatrix        m_matrixK;    // kalman gain matrix
    double                  m_chisq;
    int                     m_ierr;
    // some temporary results
    CLHEP::HepMatrix m_matrixCGT;
  };
} // namespace DecayTreeFitter

#endif
