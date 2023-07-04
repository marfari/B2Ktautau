
using namespace std;

// Global variables
int dimM = 22; // number of known parameters
int dimX = 22; // number of unkown parameters

TVectorD m( dimM );
TMatrixDSym V( dimM );
TMatrixD W( dimM, dimM );
TMatrixD Vprime( dimM, dimM );
TVectorD x0( dimX );
Double_t RPz;
Double_t chi2, MB, MB_error;
Int_t status;

double mtau = 1776.86; // MeV (PDG 2023)
double mK = 493.677; // MeV (PDG 2023) 

// Functions
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);

TVectorD x_initial_estimate( TVectorD m );
TVectorD x_initial_estimate2( TVectorD m, ROOT::Math::XYZPoint BV );

TVectorD h( TVectorD x );
TMatrixD dh_dx( TVectorD x );
double chisquare( const double* x_values );

void minimize( Double_t X[dimX], Double_t XERR[dimX], int init, ROOT::Math::XYZPoint BV );

void decay_fit()
{
    TFile* f = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2018/DVntuple_MC_2018_MagUp_visible.root");
    TTree* t = (TTree*)f->Get("ntuple/DecayTree");

    UInt_t num_entries = t->GetEntries(); 
    Double_t m_vars[dimM];
    Double_t V_vars[dimM][dimM];
    Double_t kp_pz;
    Double_t BVx, BVy, BVz;

    for(int i = 1; i <= dimM; i++)
    {
        t->SetBranchAddress(Form("df_m_%i",i), &m_vars[i-1] );
        for(int j = 1; j <= dimM; j++)
        {
            t->SetBranchAddress(Form("df_V_%i_%i",i,j), &V_vars[i-1][j-1]);
        }
    }
    t->SetBranchAddress("df_RPz", &RPz);
    t->SetBranchAddress("Bp_ENDVERTEX_X", &BVx);
    t->SetBranchAddress("Bp_ENDVERTEX_Y", &BVy);
    t->SetBranchAddress("Bp_ENDVERTEX_Z", &BVz);

    TFile* fout = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2018/decay_fit_tuple_2018.root", "RECREATE");
    TTree* tout = new TTree("DecayTree", "DecayTree");

    Double_t X[dimX];
    Double_t XERR[dimX];
  
    TString name_x[] = {
      "df_BVx",
      "df_BVy",
      "df_BVz",
      "df_Bp_PX",
      "df_Bp_PY",
      "df_Bp_PZ",
      "df_taup_PX",
      "df_taup_PY",
      "df_taup_PZ",
      "df_antinutau_PX",
      "df_antinutau_PY",
      "df_antinutau_PZ",
      "df_taum_PX",
      "df_taum_PY",
      "df_taum_PZ",
      "df_nutau_PX",
      "df_nutau_PY",
      "df_nutau_PZ",
      "df_L1",
      "df_L2",
      "df_L",
      "df_LK"
    };

    TString name_xerr[] = {
      "df_BVx_err",
      "df_BVy_err",
      "df_BVz_err",
      "df_Bp_PX_err",
      "df_Bp_PY_err",
      "df_Bp_PZ_err",
      "df_taup_PX_err",
      "df_taup_PY_err",
      "df_taup_PZ_err",
      "df_antinutau_PX_err",
      "df_antinutau_PY_err",
      "df_antinutau_PZ_err",
      "df_taum_PX_err",
      "df_taum_PY_err",
      "df_taum_PZ_err",
      "df_nutau_PX_err",
      "df_nutau_PY_err",
      "df_nutau_PZ_err",
      "df_L1_err",
      "df_L2_err",
      "df_L_err",
      "df_LK_err"
    };

    for(int i = 0; i < dimX; i++)
    {
      tout->Branch(name_x[i], &X[i]);
      tout->Branch(name_xerr[i], &XERR[i]);
    }
    tout->Branch("df_chi2", &chi2);
    tout->Branch("df_Bp_M", &MB);
    tout->Branch("df_Bp_M_err", &MB_error);
    tout->Branch("df_status", &status);

    TMatrixDSym taup_cov(7);

    // Loop over events
    for(int evt = 0; evt < num_entries; evt++)
    {
      t->GetEntry(evt);
      
      // 1)  Measured parameters, m, and their covariance matrix V
      // m = {PV, DV1, P3pi1, DV2, P3pi1, RP, PK}
      for(int i = 0; i < dimM; i++)
      {
          m(i) = m_vars[i];

          for(int j = 0; j < dimM; j++)
          {
            V(i,j) = V_vars[i][j];
          }
      }  

      // TMatrixD corr_E(dimM,dimM);
      // for(int i = 0; i < dimM; i++)
      // {
      //   for(int j = 0; j < dimM; j++)
      //   {
      //     corr_E(i,j) = V(i,j)/(sqrt(V(i,i))*sqrt(V(j,j)));
      //   }
      // }

      // Make change of variables E -> m^2 in m
      Double_t E3pi1 = m(9);
      Double_t E3pi2 = m(16);
      ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
      ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
      Double_t m3pi1_squared = pow(E3pi1,2) - p3pi1.Mag2();
      Double_t m3pi2_squared = pow(E3pi2,2) - p3pi2.Mag2();

      m(9) = m3pi1_squared; // E3pi1 -> m3pi1^2
      m(16) = m3pi2_squared; // E3pi2 -> m3pi2^2

      TMatrixDSym Jacob( dimM ); // Jacobian matrix of the transformation
      for(int i = 0; i < dimM; i++) // initiate matrix as the identity
      {
        for(int j = 0; j < dimM; j++)
        {
          if(i == j)
          {
            Jacob(i,j) = 1;
          }
          else
          {
            Jacob(i,j) = 0; 
          }
        }
      }
      Jacob(9,9) = 2*E3pi1;
      Jacob(9,6) = -2*p3pi1.x();
      Jacob(9,7) = -2*p3pi1.y();
      Jacob(9,8) = -2*p3pi1.z();

      Jacob(16,16) = 2*E3pi2;
      Jacob(16,13) = -2*p3pi2.x();
      Jacob(16,14) = -2*p3pi2.y();
      Jacob(16,15) = -2*p3pi2.z();
 
      TMatrixD Jacob_transpose = Jacob;
      Jacob_transpose.T();
      Vprime = Jacob*(V*Jacob_transpose); // transform covariance matrix V' = J V J^T
      Vprime(9,9) = abs(Vprime(9,9));
      Vprime(16,16) = abs(Vprime(16,16));

      // TMatrixD corr_m2(dimM,dimM);
      // for(int i = 0; i < dimM; i++)
      // {
      //   for(int j = 0; j < dimM; j++)
      //   {
      //     corr_m2(i,j) = Vprime(i,j)/(sqrt(Vprime(i,i))*sqrt(Vprime(j,j)));
      //   }
      // }

      // weights matrix
      W = Vprime;
      W.Invert();

      // BV (necessary for Marseille initialisation)
      ROOT::Math::XYZPoint BV( BVx, BVy, BVz );
      minimize( X, XERR, 2, BV );

      tout->Fill();
    }
    fout->cd();
    tout->Write();
    fout->Close();
}

void minimize( Double_t X[dimX], Double_t XERR[dimX], int init, ROOT::Math::XYZPoint BV )
{
  // Initial values for the unkown parameters x0
  if(init == 2) 
  {
     x0 = x_initial_estimate2( m, BV );
  }
  else if(init == 0)
  {
    x0 = x_initial_estimate( m );
  }

  Double_t x0_vars[dimX], x0_err[dimX];
  for(int i = 0; i < dimX; i++)
  {
    x0_vars[i] = x0(i); // x0_vars is a double, x0 is a TVectorD; the function to minimise must receive a double as input
    x0_err[i] = 0.1*abs(x0_vars[i]); // initial errors on x0; they are used as the first step size in the minimisation; considering a 10% error now
  }

  // x0.Print();
  // TVectorD hx = h(x0);
  // hx.Print();
  // TVectorD r = m-hx;
  // r.Print();
  // cout << chisquare(x0_vars) << endl;

  // Boundaries for x0 
  Double_t x0_min[dimX];
  x0_min[0] = -100;
  x0_min[1] = -100;
  x0_min[2] = -1000;
  x0_min[3] = -500000;
  x0_min[4] = -500000;
  x0_min[5] = -5000;
  x0_min[6] = -100000;
  x0_min[7] = -100000;
  x0_min[8] = -1000;
  x0_min[9] = -100000;
  x0_min[10] = -100000;
  x0_min[11] = -1000;
  x0_min[12] = -100000;
  x0_min[13] = -100000;
  x0_min[14] = -1000;
  x0_min[15] = -100000;
  x0_min[16] = -100000;
  x0_min[17] = -1000;
  x0_min[18] = -10;
  x0_min[19] = -10;
  x0_min[20] = -10;
  x0_min[21] = -10;

  Double_t x0_max[dimX];
  x0_max[0] = 100;
  x0_max[1] = 100;
  x0_max[2] = 1000;
  x0_max[3] = 500000;
  x0_max[4] = 500000;
  x0_max[5] = 5000000;
  x0_max[6] = 100000;
  x0_max[7] = 100000;
  x0_max[8] = 1000000;
  x0_max[9] = 100000;
  x0_max[10] = 100000;
  x0_max[11] = 1000000;
  x0_max[12] = 100000;
  x0_max[13] = 100000;
  x0_max[14] = 1000000;
  x0_max[15] = 100000;
  x0_max[16] = 100000;
  x0_max[17] = 1000000;
  x0_max[18] = 10;
  x0_max[19] = 10;
  x0_max[20] = 10;
  x0_max[21] = 10;

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");    
  ROOT::Math::Functor func(&chisquare, dimX);

  min->SetMaxIterations(1000000);  
  min->SetPrintLevel(1);
  min->SetPrecision(pow(10,-15));
  // min->SetTolerance(0.001);
  min->SetMaxFunctionCalls(1000000);

  min->SetFunction(func);

  string name[] = {
    "BVx",
    "BVy",
    "BVz",
    "Bp_PX",
    "Bp_PY",
    "Bp_PZ",
    "taup_PX",
    "taup_PY",
    "taup_PZ",
    "antinutau_PX",
    "antinutau_PY",
    "antinutau_PZ",
    "taum_PX",
    "taum_PY",
    "taum_PZ",
    "nutau_PX",
    "nutau_PY",
    "nutau_PZ",
    "L1",
    "L2",
    "L",
    "LK"
  };

  for(int i = 0; i < dimM; i++)
  {
    min->SetVariable(i, name[i], x0[i], x0_err[i]);
    min->SetVariableLimits(i, x0_min[i], x0_max[i]);
  }

  min->Minimize();

  min->PrintResults();

  // Return / save results from the fit
  chi2 = min->MinValue();
  status = min->Status();

  const double *xMin = min->X();
  const double *xMin_errors = min->Errors();
  for(int i = 0; i < dimX; i++)
  {
    X[i] = xMin[i];
    XERR[i] = xMin_errors[i];
  }

  ROOT::Math::XYZVector pB_min( X[3], X[4], X[5] );
  // double EB_min = X[6];
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  double EK = sqrt( pow(mK,2) + pK.Mag2() );
  ROOT::Math::XYZVector ptau1_min( X[6], X[7], X[8] );
  ROOT::Math::XYZVector ptau2_min( X[12], X[13], X[14] );
  double Etau1_min = sqrt( pow(mtau,2) + ptau1_min.Mag2() );
  double Etau2_min = sqrt( pow(mtau,2) + ptau2_min.Mag2() );

  // Estimate error on MB using error propagation:
  double pKx_err = sqrt( Vprime(19,19) );
  double pKy_err = sqrt( Vprime(20,20) );
  double pKz_err = sqrt( Vprime(21,21) );
  double EK_err = sqrt( pow(pK.x()/EK,2)*pow(pKx_err,2) + pow(pK.y()/EK,2)*pow(pKy_err,2) + pow(pK.z()/EK,2)*pow(pKz_err,2) );
  // cout << "EK = " << EK << " +/- " << EK_err << endl;

  double ptau1x_min_err = XERR[6];
  double ptau1y_min_err = XERR[7];
  double ptau1z_min_err = XERR[8];
  double Etau1_min_err = sqrt( pow(ptau1_min.x()/Etau1_min,2)*pow(ptau1x_min_err,2) + pow(ptau1_min.y()/Etau1_min,2)*pow(ptau1y_min_err,2) + pow(ptau1_min.z()/Etau1_min,2)*pow(ptau1z_min_err,2) );
  // cout << "Etau1_min = " << Etau1_min << " +/- " << Etau1_min_err << endl;

  double ptau2x_min_err = XERR[12];
  double ptau2y_min_err = XERR[13];
  double ptau2z_min_err = XERR[14];
  double Etau2_min_err = sqrt( pow(ptau2_min.x()/Etau2_min,2)*pow(ptau2x_min_err,2) + pow(ptau2_min.y()/Etau2_min,2)*pow(ptau2y_min_err,2) + pow(ptau2_min.z()/Etau2_min,2)*pow(ptau2z_min_err,2) );
  // cout << "Etau2_min = " << Etau2_min << " +/- " << Etau2_min_err << endl;

  double EB_min = EK + Etau1_min + Etau2_min;
  double EB_min_err = sqrt( pow(EK_err,2) + pow(Etau1_min_err,2) + pow(Etau2_min_err,2) );
  // cout << "EB_min = " << EB_min << " +/- " << EB_min_err << endl;

  if( pow(EB_min,2) > pB_min.Mag2() )
  {
    MB = sqrt( pow(EB_min,2) - pB_min.Mag2() );
  }
  else
  {
    MB = -sqrt( pB_min.Mag2() - pow(EB_min,2) );
  }

  double pBx_min_err = XERR[3]; 
  double pBy_min_err = XERR[4]; 
  double pBz_min_err = XERR[5]; 
  MB_error = sqrt( pow(EB_min/MB,2)*pow(EB_min_err,2) + pow(pB_min.x()/MB,2)*pow(pBx_min_err,2) + pow(pB_min.y()/MB,2)*pow(pBy_min_err,2) + pow(pB_min.z()/MB,2)*pow(pBz_min_err,2) );

  // cout << "initial chi2 = " << chisquare(x0_vars) << endl;
  // cout << "final chi2 = " << min->MinValue() << endl;
  cout << "MB = " << MB << " +/- " << MB_error << endl;
  cout << "Status = " << status << endl;

  // TVectorD x_min(dimX);
  // for(int i = 0; i < dimX; i++)
  // {
  //   x_min(i) = xMin[i];
  // }

  // x_min.Print();

  // TVectorD hx_min = h(x_min);
  // hx_min.Print();

  // TVectorD r_min = m-hx_min;
  // r_min.Print();
}

double chisquare( const double* x_values )
{
    TVectorD x(dimX);
    for(int i = 0; i < dimX; i++)
    {
      x(i) = x_values[i];
    }

    TVectorD hx = h( x );

    TVectorD r = m - h(x);

    double chi2 = r*(W*r);

    // if(chi2 < 0)
    // {
    //   return 100000000;
    // }

    return chi2;
}

TVectorD h( TVectorD x )
{
  // h(x)
  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  // double EB = x(6);
  ROOT::Math::XYZVector ptau1( x(6), x(7), x(8) );
  ROOT::Math::XYZVector pnu1( x(9), x(10), x(11) );
  ROOT::Math::XYZVector ptau2( x(12), x(13), x(14) );
  ROOT::Math::XYZVector pnu2( x(15), x(16), x(17) );
  double L1 = x(18);
  double L2 = x(19);
  double L = x(20);
  double LK = x(21);

  double Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  double Enu1 = sqrt( pnu1.Mag2() );
  double Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  double Enu2 = sqrt( pnu2.Mag2() );

  ROOT::Math::PxPyPzEVector Ptau1( ptau1.x(), ptau1.y(), ptau1.z(), Etau1 );
  ROOT::Math::PxPyPzEVector Ptau2( ptau2.x(), ptau2.y(), ptau2.z(), Etau2 );
  ROOT::Math::PxPyPzEVector Pnu1( pnu1.x(), pnu1.y(), pnu1.z(), Enu1 );
  ROOT::Math::PxPyPzEVector Pnu2( pnu2.x(), pnu2.y(), pnu2.z(), Enu2 );

  ROOT::Math::XYZPoint PV( BV.x() - L*pB.x(), BV.y() - L*pB.y(), BV.z() - L*pB.z() );
  ROOT::Math::XYZPoint DV1( BV.x() + L1*ptau1.x(), BV.y() + L1*ptau1.y(), BV.z() + L1*ptau1.z() );
  ROOT::Math::XYZVector p3pi1 = ptau1 - pnu1;
  // double E3pi1 = Etau1 - Enu1;
  double m3pi1_squared = pow(mtau,2) -2*Ptau1.Dot(Pnu1);
  ROOT::Math::XYZPoint DV2( BV.x() + L2*ptau2.x(), BV.y() + L2*ptau2.y(), BV.z() + L2*ptau2.z() );
  ROOT::Math::XYZVector p3pi2 = ptau2 - pnu2;
  // double E3pi2 = Etau2 - Enu2;
  double m3pi2_squared = pow(mtau,2) -2*Ptau2.Dot(Pnu2);
  ROOT::Math::XYZPoint RP( BV.x() + LK*( pB.x() - ptau1.x() - ptau2.x() ), BV.y() + LK*( pB.y() - ptau1.y() - ptau2.y() ), RPz );
  ROOT::Math::XYZVector pK = pB - ptau1 - ptau2;
  // double EK = EB - Etau1 - Etau2;

  TVectorD h( dimM );
  h(0) = PV.x();
  h(1) = PV.y();
  h(2) = PV.z();
  h(3) = DV1.x();
  h(4) = DV1.y();
  h(5) = DV1.z();
  h(6) = p3pi1.x();
  h(7) = p3pi1.y();
  h(8) = p3pi1.z();
  h(9) = m3pi1_squared;
  // h(9) = E3pi1;
  h(10) = DV2.x();
  h(11) = DV2.y();
  h(12) = DV2.z();
  h(13) = p3pi2.x();
  h(14) = p3pi2.y();
  h(15) = p3pi2.z();
  h(16) = m3pi2_squared;
  // h(16) = E3pi2;
  h(17) = RP.x();
  h(18) = RP.y();
  // h(19) = RP.z();
  h(19) = pK.x();
  h(20) = pK.y();
  h(21) = pK.z();
  // h(22) = EK;

  return h;
}

TVectorD x_initial_estimate2( TVectorD m, ROOT::Math::XYZPoint BV ) // Marseille initialisation for x
{
  // Builds an initial estimate for x based on the analytical calculations and on the known parameters in m

  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  double m3pi1_squared = m(9);
  // double E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  // double E3pi2 = m(16);
  double m3pi2_squared = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  // double EK = m(22);

  double EK = sqrt( pow(mK,2) + pK.Mag2() );

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  double E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  double E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );

  double m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  double m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  double theta1 = asin( ( pow(mtau,2) - pow(m3pi1,2) )/( 2*mtau*sqrt( p3pi1.Mag2() ) ) );
  double theta2 = asin( ( pow(mtau,2) - pow(m3pi2,2) )/( 2*mtau*sqrt( p3pi2.Mag2() ) ) );
  
  double ptau1_mag = ( (pow(mtau,2) + pow(m3pi1,2))*sqrt(p3pi1.Mag2())*cos(theta1) )/( 2*( pow(E3pi1,2) - p3pi1.Mag2()*pow(cos(theta1),2) ) );
  double ptau2_mag = ( (pow(mtau,2) + pow(m3pi2,2))*sqrt(p3pi2.Mag2())*cos(theta2) )/( 2*( pow(E3pi2,2) - p3pi2.Mag2()*pow(cos(theta2),2) ) );

  ROOT::Math::XYZVector ptau1 = ptau1_mag*u1;
  ROOT::Math::XYZVector ptau2 = ptau1_mag*u2;

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  // ROOT::Math::XYZVector u = (BV - PV).Unit();
  // double pB_mag = sqrt( (pK + ptau1 + ptau2).Mag2() );
  // ROOT::Math::XYZVector pB = pB_mag*u;

  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;

  double Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  double Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  // double EB = Etau1 + Etau2 + EK;

  double Enu1 = sqrt( pnu1.Mag2() );
  double Enu2 = sqrt( pnu2.Mag2() );
  double EB = E3pi1 + Enu1 + E3pi2 + Enu2 + EK;

  double L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  double L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  double LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );
  double L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  // x0(6) = EB;
  x0(6) = ptau1.x();
  x0(7) = ptau1.y();
  x0(8) = ptau1.z();
  x0(9) = pnu1.x();
  x0(10) = pnu1.y();
  x0(11) = pnu1.z();
  x0(12) = ptau2.x();
  x0(13) = ptau2.y();
  x0(14) = ptau2.z();
  x0(15) = pnu2.x();
  x0(16) = pnu2.y();
  x0(17) = pnu2.z();
  x0(18) = L1;
  x0(19) = L2;
  x0(20) = L;
  x0(21) = LK;

  return x0;
}

TVectorD x_initial_estimate( TVectorD m ) // Original initialisation for x (based on Anne Keune's thesis)
{
  // Builds an initial estimate for x based on the analytical calculations and on the known parameters in m

  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  double E3pi1 = m(6);
  ROOT::Math::XYZVector p3pi1( m(7), m(8), m(9) );
  // double E3pi1 = m(9);
  // double m3pi1_squared = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  double E3pi2 = m(13);
  ROOT::Math::XYZVector p3pi2( m(14), m(15), m(16) );
  // double E3pi2 = m(16);
  // double m3pi2_squared = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  // double EK = m(22);

  // double E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  // double E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );
  double EK = sqrt( pow(mK,2) + pK.Mag2() );

  double m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  double m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  ROOT::Math::XYZPoint PV_t = makeTransformation_point( pK, RP, PV, false );
  ROOT::Math::XYZPoint DV1_t = makeTransformation_point( pK, RP, DV1, false );
  ROOT::Math::XYZVector p3pi1_t = makeTransformation_vec( pK, RP, p3pi1, false );
  ROOT::Math::XYZPoint DV2_t = makeTransformation_point( pK, RP, DV2, false );
  ROOT::Math::XYZVector p3pi2_t = makeTransformation_vec( pK, RP, p3pi2, false );
  ROOT::Math::XYZVector pK_t = makeTransformation_vec( pK, RP, pK, false );

  double a1 = (DV1_t.y())/(DV1_t.x());
  double a2 = (DV2_t.y())/(DV2_t.x());
  double b = (PV_t.y() -a1*PV_t.x())/(a2*PV_t.x() - PV_t.y());
  double c = b*(DV1_t.x())/(DV2_t.x());
  double d = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()));
  double e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
  double f = ( PV_t.x()*sqrt(pK_t.Mag2()) )/( (1+b)*DV1_t.x() -(1+c)*PV_t.x() );
  double g = c*e + d;
  double h = f*c;
  double i = DV1_t.z() - e*DV1_t.x();
  double j = f*DV1_t.x();

  double x1 = p3pi1_t.x() + a1*p3pi1_t.y() + e*p3pi1_t.z();
  double x2 = b*p3pi2_t.x() + a2*b*p3pi2_t.y() + g*p3pi2_t.z();

  double p1 = 1 + pow(a1,2) + pow(e,2) - pow(x1/E3pi1,2);
  double p2 = 2*e*f - ( pow(mtau,2) + pow(m3pi1,2) + 2*f*p3pi1_t.z() )*(x1/pow(E3pi1,2) );
  double p3 = pow(mtau,2) + pow(f,2) - pow( ( pow(mtau,2) + pow(m3pi1,2) + 2*f*p3pi1_t.z() )/(2*E3pi1), 2);
  double q1 = pow(b,2) + pow(a2*b,2) + pow(g,2) - pow(x2/E3pi2,2);
  double q2 = 2*g*h - ( pow(mtau,2) + pow(m3pi2,2) + 2*h*p3pi2_t.z() )*(x2/pow(E3pi2,2) );
  double q3 =  pow(mtau,2) + pow(h,2) - pow( ( pow(mtau,2) + pow(m3pi2,2) + 2*h*p3pi2_t.z() )/(2*E3pi2),2 );

  double Ptau1x_t = (p1*q3 - p3*q1)/(p2*q1 - p1*q2);
  double Ptau1y_t = a1*Ptau1x_t;
  double Ptau1z_t = e*Ptau1x_t + f;
  double Ptau2x_t = b*Ptau1x_t;
  double Ptau2y_t = a2*b*Ptau1x_t;
  double Ptau2z_t = g*Ptau1x_t + h;
  double BVz_t = i - j*(1/Ptau1x_t);

  ROOT::Math::XYZVector ptau1_t( Ptau1x_t, Ptau1y_t, Ptau1z_t );
  ROOT::Math::XYZVector ptau2_t( Ptau2x_t, Ptau2y_t, Ptau2z_t );
  ROOT::Math::XYZPoint BV_t(0, 0, BVz_t);

  // transform back to LHCb frame
  ROOT::Math::XYZVector ptau1 = makeTransformation_vec( pK, RP, ptau1_t, true );
  ROOT::Math::XYZVector ptau2 = makeTransformation_vec( pK, RP, ptau2_t, true );
  ROOT::Math::XYZPoint BV = makeTransformation_point( pK, RP, BV_t, true );

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;

  double Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  double Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  double Enu1 = sqrt( pnu1.Mag2() );
  double Enu2 = sqrt( pnu2.Mag2() );
  double EB = Etau1 + Etau2 + EK;

  double L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  double L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  double L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );
  double LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  // x0(6) = EB;
  x0(6) = ptau1.x();
  x0(7) = ptau1.y();
  x0(8) = ptau1.z();
  x0(9) = pnu1.x();
  x0(10) = pnu1.y();
  x0(11) = pnu1.z();
  x0(12) = ptau2.x();
  x0(13) = ptau2.y();
  x0(14) = ptau2.z();
  x0(15) = pnu2.x();
  x0(16) = pnu2.y();
  x0(17) = pnu2.z();
  x0(18) = L1;
  x0(19) = L2;
  x0(20) = L;
  x0(21) = LK;

  return x0;
}

ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  double theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

  ROOT::Math::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  ROOT::Math::Transform3D myTransformation     = ROOT::Math::Transform3D(myRotation) * ROOT::Math::Transform3D(myShift);
  ROOT::Math::Transform3D myTransformation_inv = myTransformation.Inverse();

  ROOT::Math::XYZVector theVector_t;
  if(invFlag)
    theVector_t = myTransformation_inv * theVector;
  else 
    theVector_t = myTransformation * theVector;
  
  return theVector_t;
}

ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  double theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

  ROOT::Math::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  ROOT::Math::Transform3D myTransformation     = ROOT::Math::Transform3D(myRotation) * ROOT::Math::Transform3D(myShift);
  ROOT::Math::Transform3D myTransformation_inv = myTransformation.Inverse();

  ROOT::Math::XYZPoint thePoint_t;
  if(invFlag)
    thePoint_t = myTransformation_inv * thePoint;
  else 
    thePoint_t = myTransformation * thePoint;
  
  return thePoint_t;
}