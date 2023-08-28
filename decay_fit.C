
using namespace std;

////////////////////////////////////////             m6piK^2 parametrisation from track parametrisation           ////////////////////////////////////////////////////////

// Global variables
int dimM = 23; // number of known parameters
int dimX = 23; // number of unkown parameters

int dimM_before = 32;

TVectorD m( dimM_before ); // vector of measured parameters
TMatrixDSym V( dimM_before ); // covariance matrix of m

TVectorD mprime( dimM );
TMatrixDSym Vprime( dimM );
TMatrixDSym W( dimM ); // weights matrix, W=V^-1
TVectorD x0( dimX ); // initial estimate for the vector of unkown parameters

TMatrixDSym taup_cov(7);
TMatrixDSym taum_cov(7);
TMatrixDSym Bp_cov(7);

Double_t RPz = 0; // z-component of the reference point in the K+ trajectory (fixed)

// Things saved from the minimisation:
Double_t chi2, MB, MB_err, taup_PE, taum_PE, taup_M, taum_M; 
TVectorD x_m(dimX); // vector of unkown parameters after the minimisation
TMatrixD C_m(dimX,dimX); // covariance matrix of the fit parameters
Int_t status; // status of the fit

Double_t mtau = 1776.86; // MeV (PDG 2023)
Double_t mpion = 139.57039; // MeV
Double_t mkaon = 493.677; // MeV

// Functions
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);
TVectorD x_initial_estimate( TVectorD m ); // Original initialisation
TVectorD x_initial_estimate0( TVectorD mprime, ROOT::Math::XYZPoint BV ); // B->K*tautau initialisation; taus direction from vertices
TVectorD x_initial_estimate1( TVectorD mprime, ROOT::Math::XYZPoint BV ); // Using B->K* tautau initialisation; taus direction based visible 3pi momenta
TVectorD x_initial_estimate2( TVectorD m, ROOT::Math::XYZPoint BV ); // Marseille's initialisation
TVectorD h( TVectorD x );
TMatrixD dh_dx( TVectorD x );
Double_t chisquare( const Double_t* x_values );
void minimize( ROOT::Math::XYZPoint BV, int init );
TVectorD transform_m( TVectorD m );
TMatrixDSym transform_V( TVectorD m, TMatrixDSym V, TMatrixDSym taup_cov, TMatrixDSym taum_cov );

void decay_fit(int year, TString MC_files, TString RS_DATA_files, TString WS_DATA_files, int species, int component, int line)
{
  TFileCollection* fc;
  if(species == 0)
  { 
    fc = new TFileCollection("MC", "MC", MC_files, 1, line);
  }
  else if(species == 1)
  {
    fc = new TFileCollection("MC", "MC", RS_DATA_files, 1, line);
  }
  else if(species == 2)
  {
    fc = new TFileCollection("MC", "MC", WS_DATA_files, 1, line);
  }

  TChain* t = new TChain("DecayTree");
  t->AddFileInfoList((TCollection*)fc->GetList());

  UInt_t num_entries = t->GetEntries();
  Double_t m_vars[dimM_before];
  Double_t V_vars[dimM_before][dimM_before];
  Double_t BVx, BVy, BVz;

  Double_t taup_cov_vars[7][7];
  Double_t taum_cov_vars[7][7];
  Double_t Bp_cov_vars[7][7];

  for(int i = 1; i <= dimM_before; i++)
  {
      t->SetBranchAddress(Form("df_m_%i",i), &m_vars[i-1] );
      for(int j = 1; j <= dimM_before; j++)
      {
          t->SetBranchAddress(Form("df_V_%i_%i",i,j), &V_vars[i-1][j-1]);
      }
  }
  t->SetBranchAddress("Kp_RP_Z", &RPz);
  t->SetBranchAddress("Bp_ENDVERTEX_X", &BVx);
  t->SetBranchAddress("Bp_ENDVERTEX_Y", &BVy);
  t->SetBranchAddress("Bp_ENDVERTEX_Z", &BVz);

  for(int i = 0; i < 7; i++)
  {
    for(int j = 0; j < 7; j++)
    {
      t->SetBranchAddress(Form("taup_cov_%i_%i",i,j), &taup_cov_vars[i][j]);
      t->SetBranchAddress(Form("taum_cov_%i_%i",i,j), &taum_cov_vars[i][j]);
      t->SetBranchAddress(Form("Bp_cov_%i_%i",i,j), &Bp_cov_vars[i][j]);
    }
  }

  TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/Component_%i/fit_result_%i.root",year,species,component,line), "RECREATE");
  TTree* tout = new TTree("DecayTree", "DecayTree");

  TString name_x[] = {
    "df_BVx",
    "df_BVy",
    "df_BVz",
    "df_Bp_PX",
    "df_Bp_PY",
    "df_Bp_PZ",
    "df_Bp_M2",
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

  for(int i = 0; i < dimX; i++)
  {
    tout->Branch(name_x[i], &x_m(i));
    for(int j = 0; j < dimX; j++)
    {
      tout->Branch(Form("df_C_%i_%i",i,j), &C_m(i,j));
    }
  }
  tout->Branch("df_chi2", &chi2);
  tout->Branch("df_Bp_M", &MB);
  tout->Branch("df_Bp_MERR", &MB_err);
  tout->Branch("df_status", &status);
  tout->Branch("df_taup_PE", &taup_PE);
  tout->Branch("df_taum_PE", &taum_PE);
  tout->Branch("df_taup_M", &taup_M);
  tout->Branch("df_taum_M", &taum_M);

  // Loop over events
  for(int evt = 0; evt < num_entries; evt++)
  {
    t->GetEntry(evt);
    
    // 1)  Measured parameters, m, and their covariance matrix V in the track parametrisation
    // m = (PV,DV1,p1,p2,p3,DV2,p4,p5,p6,RP_T,pK) -> 32 parameters
    for(int i = 0; i < dimM_before; i++)
    {
        m(i) = m_vars[i];
        for(int j = 0; j < dimM_before; j++)
        {
          V(i,j) = V_vars[i][j];
        }
    }  
    // m.Print();
    // V.Print();

    // Eigenvalues of V in track parametrisation
    // TVectorD V_eigen_values(dimM_before);  
    // TMatrixD V_eigen_vectors = V.EigenVectors(V_eigen_values);
    // V_eigen_values.Print();

    // for(int i = 0; i < 7; i++)
    // {
    //   for(int j = 0; j < 7; j++)
    //   {
    //     taup_cov(i,j) = taup_cov_vars[i][j];
    //     taum_cov(i,j) = taum_cov_vars[i][j];
    //     Bp_cov(i,j) = Bp_cov_vars[i][j];
    //   }
    // }
    // taup_cov.Print();
    // taum_cov.Print();
    // Bp_cov.Print();

    // Vprime.SetTol(pow(10,-23));
    // W.SetTol(pow(10,-23));
    // mprime = transform_m(m);
    // Vprime = transform_V(m,V,taup_cov,taum_cov);
    // mprime.Print();
    // Vprime.Print();

    ////////////////////////////////////   Tikhonov regularisation   /////////////////////////////////////////////////
    // 1) Get eigenvalues of V' (before regularisation)
    // TVectorD V_eigen_values(dimM);  
    // TMatrixD V_eigen_vectors = Vprime.EigenVectors(V_eigen_values);
    // V_eigen_values.Print();

    // Correlation matrix of m' (before regularisation)
    // TMatrixD corr(dimM,dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //   for(int j = 0; j < dimM; j++)
    //   {
    //     corr(i,j) = Vprime(i,j)/(sqrt(Vprime(i,i))*sqrt(Vprime(j,j)));
    //   }
    // }
    // corr.Print();

    // // 2) Check if these sub-matrices have negative eigenvalues and if so apply Tikhonov regularisation to the sub-matrix: add a constant lambda=1.1*abs(q) to diagonal elements, where q is the value of the negative eigenvalue
    // // note: the eigenvalues of the sub-matrices are a sub-group of the eigenvalues of the large matric V
    // Double_t lambda = 10*V_eigen_values(dimM-1);
    // for(int i = 0; i < dimM; i++)
    // {
    //   if( (lambda < abs(V_eigen_values(i))) )
    //   {
    //     lambda = 10*abs( V_eigen_values(i) );
    //   }
    // }
    // cout << lambda << endl;

    // // 3) Regularise V'
    // for(int i = 0; i < dimM; i++)
    // {
    //   Vprime(i,i) += lambda; 
    // }
    // Vprime.Print();
    
    // Eigenvalues of V' (after regularisation)
    // TVectorD V_eigen_values_a(dimM);  
    // TMatrixD V_eigen_vectors_a = Vprime.EigenVectors(V_eigen_values_a);
    // V_eigen_values_a.Print();

    // Correlation matrix of m' (after regularisation)
    // TMatrixD corr_p(dimM,dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //   for(int j = 0; j < dimM; j++)
    //   {
    //     corr_p(i,j) = Vprime(i,j)/(sqrt(Vprime(i,i))*sqrt(Vprime(j,j)));
    //   }
    // }
    // corr_p.Print();

    // ////////////////////////////////////////////// CHANGE OF VARIABLES E -> m^2 ///////////////////////////////////////////
    // // Make change of variables E -> m^2 in m
    // Double_t E3pi1 = m(9);
    // Double_t E3pi2 = m(16);
    // Double_t E6piK = m(22);
    // ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
    // ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
    // ROOT::Math::XYZVector p6piK( m(19), m(20), m(21) );
    // Double_t m3pi1_squared = pow(E3pi1,2) - p3pi1.Mag2();
    // Double_t m3pi2_squared = pow(E3pi2,2) - p3pi2.Mag2();
    // Double_t m6piK_squared = pow(E6piK,2) - p6piK.Mag2();

    // m(9) = m3pi1_squared; // E3pi1 -> m3pi1^2
    // m(16) = m3pi2_squared; // E3pi2 -> m3pi2^2
    // m(22) = m6piK_squared; // E6piK -> m6piK^2

    // TMatrixD Jacob( dimM, dimM ); // Jacobian matrix of the transformation
    // for(int i = 0; i < dimM; i++) // initiate matrix as the identity
    // {
    //   for(int j = 0; j < dimM; j++)
    //   {
    //     if(i == j)
    //     {
    //       Jacob(i,j) = 1;
    //     }
    //     else
    //     {
    //       Jacob(i,j) = 0; 
    //     }
    //   }
    // }
    // Jacob(9,9) = 2*E3pi1;
    // Jacob(9,6) = -2*p3pi1.x();
    // Jacob(9,7) = -2*p3pi1.y();
    // Jacob(9,8) = -2*p3pi1.z();

    // Jacob(16,16) = 2*E3pi2;
    // Jacob(16,13) = -2*p3pi2.x();
    // Jacob(16,14) = -2*p3pi2.y();
    // Jacob(16,15) = -2*p3pi2.z();

    // Jacob(22,22) = 2*E6piK;
    // Jacob(22,19) = -2*p6piK.x();
    // Jacob(22,20) = -2*p6piK.y();
    // Jacob(22,21) = -2*p6piK.z();

    // TMatrixD Jacob_transpose = Jacob;
    // // Jacob.Print();
    // Jacob_transpose.T();
    // Vprime = Jacob*(V*Jacob_transpose); // transform covariance matrix V' = J V J^T (approximation to 1st order in Taylor expansion)
    // // Vprime(9,9) = abs(Vprime(9,9)); // force this to be positive for now; why is it negative in the 1st place?
    // // Vprime(16,16) = abs(Vprime(16,16));
    // // Vprime(22,22) = abs(Vprime(22,22));
    // // m.Print();
    // // Vprime.Print();

    // // Eigenvalues of V' (after regularisation)
    // // TVectorD Vp_eigen_values(dimM);  
    // // TMatrixD Vp_eigen_vectors = Vprime.EigenVectors(Vp_eigen_values);
    // // Vp_eigen_values.Print();

    // // // Correlation matrix of m (after change of variables)
    // // TMatrixD corr_m2(dimM,dimM);
    // // for(int i = 0; i < dimM; i++)
    // // {
    // //   for(int j = 0; j < dimM; j++)
    // //   {
    // //     corr_m2(i,j) = Vprime(i,j)/(sqrt(Vprime(i,i))*sqrt(Vprime(j,j)));
    // //   }
    // // }
    // // corr_m2.Print();
    // // return;
    // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // // Weights matrix W = V'^-1
    W = Vprime;
    W.Invert();
    // W.Print();

    ROOT::Math::XYZPoint BV( BVx, BVy, BVz ); // BV (necessary for Marseille initialisation)

    // 2130
    // best sequence so far
    int init = 2; // Marseille
    minimize( BV, init );
    if( (status != 0) && (init == 2) )
    {
      init = 1; // K*tautau w/ visible 3pi momenta
      cout << "status before init = 1 : " << status << endl;
      minimize( BV, init );
    }
    if( (status != 0) && (init == 1) )
    {
      init = -1; // Anne Keune's
      minimize( BV, init );
    }
    if( (status != 0) && (init == -1) )
    {
      init = 0; // K*tautau from vertices
      minimize( BV, init );
    }
  
    // 012
    // int init = 0; // K*tautau from vertices
    // minimize( BV, init );
    // if( (status != 0) && (init == 0) )
    // {
    //   init = 1; // K*tautau w/ visible 3pi momenta
    //   minimize( BV, init );
    // }
    // if( (status != 0) && (init == 1) )
    // {
    //   init = 2; // Marseille's
    //   minimize( BV, init );
    // }
  
    // 210
    // int init = 2; 
    // minimize( BV, init );
    // if( (status != 0) && (init == 2) )
    // {
    //   init = 1; 
    //   minimize( BV, init );
    // }
    // if( (status != 0) && (init == 1) )
    // {
    //   init = 0; // Marseille's
    //   minimize( BV, init );
    // }

    tout->Fill();
  }
  cout << "FINISHED" << endl;

  fout->cd();
  tout->Write();
  fout->Close();
}

void minimize( ROOT::Math::XYZPoint BV, int init )
{
  // 1) This function initialises the vector of unkown parameters x with the result of analytical calculations; by specifying the value of init we choose which analytical calculation we want to use: 
  //   init = 0 : original calculations based on Anne Keune thesis; does not use the BV offline estimate 
  //   init = 2 : uses the Marseille initialisation; it uses the BV offline estimate and assumes that the neutrino takes as much momentum as possible from the tau 
  // 2) It sets up the minimizer: defining the initial values and bounds on the x parameters  
  // 3) It builds the chi^2 function with x0 by calling the function chisquare 
  // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

  // Initial values for the unkown parameters x0
  if(init == 0) // B->K*tautau initialisation; taus direction from vertices
  {
    x0 = x_initial_estimate0( mprime, BV );
  }
  else if(init == 1)
  {
    x0 = x_initial_estimate1( mprime, BV );
  }
  else if(init == 2) 
  {
     x0 = x_initial_estimate2( mprime, BV );
  }
  else if(init == -1)
  {
    x0 = x_initial_estimate( mprime ); // Original initialisation
  }

  Double_t x0_vars[dimX], x0_err[dimX];
  for(int i = 0; i < dimX; i++)
  {
    x0_vars[i] = x0(i); // x0_vars is a Double_t, x0 is a TVectorD; the function to minimise must receive a Double_t as input
    x0_err[i] = 0.1*abs(x0_vars[i]); // initial errors on x0; they are used as the first step size in the minimisation; considering a 10% error for now
  }

  // x0.Print();
  // mprime.Print();
  // TVectorD hx = h(x0);
  // hx.Print();
  // TVectorD r = mprime-hx;
  // r.Print();
  // cout << chisquare(x0_vars) << endl;

  // Boundaries for x0 
  Double_t x0_min[] = {
    -100, // BVx
    -100, // BVy
    -1000, // BVz
    -500000, // pBx
    -500000, // pBy
    -1000., // pBz
    0., // MB^2
    -100000, // ptau1x
    -100000, // ptau1y
    -1000, // ptau1z
    -100000, // pnu1x
    -100000, // pnu1y
    -1000, // pnu1z
    -100000, // ptau2x
    -100000, // ptau2y
    -1000, // ptau2z
    -100000, // pnu2x
    -100000, // pnu2y
    -1000, // pnu2z
    -10, // L1
    -10, // L2
    -10, // L
    -10 // LK
  };
  if( sizeof(x0_min)/sizeof(x0_min[0]) != dimX ){ 
    cout << "ERROR: x0_min does not have the proper dimension" << endl;
    return; 
  }

  Double_t x0_max[] = {
    100, // BVx
    100, // BVy
    100, // BVz
    500000, // pBx
    500000, // pBy
    5000000, // pBz
    100000000, // MB^2
    100000, // ptau1x
    100000, // ptau1y
    1000000, // ptau1z
    100000, // pnu1x
    100000, // pnu1y
    1000000, // pnu1z
    100000, // ptau2x
    100000, // ptau2y
    1000000, // ptau2z
    100000, // pnu2x
    100000, // pnu2y
    1000000, // pnu2z
    10, // L1
    10, // L2
    10, // L
    10 // LK
  };
  if( sizeof(x0_min)/sizeof(x0_min[0]) != dimX ){ 
    cout << "ERROR: x0_max does not have the proper dimension" << endl;
    return; 
  }

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");    
  ROOT::Math::Functor func(&chisquare, dimX);

  min->SetMaxIterations(100000000);  
  min->SetMaxFunctionCalls(100000000);
  min->SetPrintLevel(1);
  // min->SetPrecision(pow(10,-15));
  min->SetTolerance(0.000000001);
  min->SetFunction(func);

  string name[] = {
    "BVx",
    "BVy",
    "BVz",
    "Bp_PX",
    "Bp_PY",
    "Bp_PZ",
    "Bp_M2",
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
  if( sizeof(name)/sizeof(name[0]) != dimX ){ 
    cout << "ERROR: name does not have the proper dimension" << endl;
    return; 
  }

  for(int i = 0; i < dimM; i++)
  {
    min->SetVariable(i, name[i], x0[i], x0_err[i]);
    // min->SetVariableLimits(i, x0_min[i], x0_max[i]);
  }
  min->Minimize();
  // min->PrintResults();

  // Return / save results from the fit
  chi2 = min->MinValue();
  status = min->Status();

  const Double_t *xMin = min->X();
  for(int i = 0; i < dimX; i++)
  {
    x_m(i) = xMin[i];
    for(int j = 0; j < dimX; j++)
    {
      C_m(i,j) = min->CovMatrix(i,j);
    }
  }
  // x_m.Print();
  // TVectorD hx_m = h(x_m);
  // hx_m.Print();
  // TVectorD r_m = mprime-hx_m;
  // r_m.Print();
  // C_m.Print();

  // Calculate MB and its error
  Double_t MB_squared_m = x_m(6);
  if( MB_squared_m > 0 )
  {
    MB = sqrt( MB_squared_m );
  }
  else
  {
    MB = -sqrt( abs(MB_squared_m) );
  }
  Double_t dMB_squared = sqrt(C_m(6,6));
  MB_err = dMB_squared/(2*MB);

  ROOT::Math::XYZVector ptau1_m( x_m(7), x_m(8), x_m(9) );
  ROOT::Math::XYZVector pnu1_m( x_m(10), x_m(11), x_m(12) );
  ROOT::Math::XYZVector ptau2_m( x_m(13), x_m(14), x_m(15) );
  ROOT::Math::XYZVector pnu2_m( x_m(16), x_m(17), x_m(18) );

  Double_t m3pi1_squared = mprime(9);
  Double_t m3pi2_squared = mprime(16);
  ROOT::Math::XYZVector p3pi1( mprime(6), mprime(7), mprime(8) );
  ROOT::Math::XYZVector p3pi2( mprime(13), mprime(14), mprime(15) );
  Double_t E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  Double_t E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );

  Double_t Enu1_m = sqrt( pnu1_m.Mag2() ); 
  Double_t Enu2_m = sqrt( pnu2_m.Mag2() ); 

  taup_PE = E3pi1 + Enu1_m;
  taum_PE = E3pi2 + Enu2_m;

  taup_M = sqrt( pow(taup_PE,2) - ptau1_m.Mag2() );
  taum_M = sqrt( pow(taum_PE,2) - ptau2_m.Mag2() );

  // cout << "initial chi2 = " << chisquare(x0_vars) << endl;
  cout << "final chi2 = " << min->MinValue() << endl;
  cout << "MB = " << MB << " +/- " << MB_err << endl;
  cout << "Status = " << status << endl;

  // Correlation matrix of fit parameters
  // TMatrixD Corr_x(dimX,dimX);
  // for(int i = 0; i < dimX; i++)
  // {
  //   for(int j = 0; j < dimX; j++)
  //   {
  //     Corr_x(i,j) = C_m(i,j)/( sqrt(C_m(i,i))*sqrt(C_m(j,j)) );
  //   }
  // }
  // Corr_x.Print();

}

Double_t chisquare( const Double_t* x_values )
{
  // Computes the chi^2
  TVectorD x(dimX);
  for(int i = 0; i < dimX; i++)
  {
    x(i) = x_values[i];
  }

  TVectorD hx = h( x );

  TVectorD r = mprime - h(x);

  Double_t chi2 = r*(W*r);

  return chi2;
}

TVectorD h( TVectorD x )
{
  // The model; writes the known parameters in terms of the unkown parameters using the model constraints
  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  Double_t MB_squared = x(6);
  ROOT::Math::XYZVector ptau1( x(7), x(8), x(9) );
  ROOT::Math::XYZVector pnu1( x(10), x(11), x(12) );
  ROOT::Math::XYZVector ptau2( x(13), x(14), x(15) );
  ROOT::Math::XYZVector pnu2( x(16), x(17), x(18) );
  Double_t L1 = x(19);
  Double_t L2 = x(20);
  Double_t L = x(21);
  Double_t LK = x(22);

  Double_t EB = sqrt( MB_squared + pB.Mag2() );

  // Tau and neutrino mass constraints:
  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  ROOT::Math::PxPyPzEVector Ptau1( ptau1.x(), ptau1.y(), ptau1.z(), Etau1 );
  ROOT::Math::PxPyPzEVector Ptau2( ptau2.x(), ptau2.y(), ptau2.z(), Etau2 );
  ROOT::Math::PxPyPzEVector Pnu1( pnu1.x(), pnu1.y(), pnu1.z(), Enu1 );
  ROOT::Math::PxPyPzEVector Pnu2( pnu2.x(), pnu2.y(), pnu2.z(), Enu2 );
  ROOT::Math::PxPyPzEVector PB( pB.x(), pB.y(), pB.z(), EB);

  // pB must point back to the PV:
  ROOT::Math::XYZPoint PV( BV.x() - L*pB.x(), BV.y() - L*pB.y(), BV.z() - L*pB.z() );
  // ptau1 must point back to BV:
  ROOT::Math::XYZPoint DV1( BV.x() + L1*ptau1.x(), BV.y() + L1*ptau1.y(), BV.z() + L1*ptau1.z() );
  // 4-momentum conservation in DV1:
  ROOT::Math::XYZVector p3pi1 = ptau1 - pnu1;
  // Double_t E3pi1 = Etau1 - Enu1; 
  Double_t m3pi1_squared = pow(mtau,2) -2*Ptau1.Dot(Pnu1); 
  // ptau2 must point back to DV2:
  ROOT::Math::XYZPoint DV2( BV.x() + L2*ptau2.x(), BV.y() + L2*ptau2.y(), BV.z() + L2*ptau2.z() );
  // 4-momentum conservation in DV2:
  ROOT::Math::XYZVector p3pi2 = ptau2 - pnu2;
  // Double_t E3pi2 = Etau2 - Enu2; 
  Double_t m3pi2_squared = pow(mtau,2) -2*Ptau2.Dot(Pnu2);  
  // B+ must lie in the K+ trajectory:
  ROOT::Math::XYZPoint RP( BV.x() + LK*( pB.x() - ptau1.x() - ptau2.x() ), BV.y() + LK*( pB.y() - ptau1.y() - ptau2.y() ), RPz );
  // 4-momentum conservation in BV:
  ROOT::Math::XYZVector p6piK = pB - pnu1 - pnu2;
  // Double_t E6piK = EB - Enu1 - Enu2;
  Double_t m6piK_squared = MB_squared - 2*PB.Dot(Pnu1) - 2*PB.Dot(Pnu2) + 2*Pnu1.Dot(Pnu2);

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
  // h(9) = E3pi1; 
  h(9) = m3pi1_squared;
  h(10) = DV2.x();
  h(11) = DV2.y();
  h(12) = DV2.z();
  h(13) = p3pi2.x();
  h(14) = p3pi2.y();
  h(15) = p3pi2.z(); 
  // h(16) = E3pi2;
  h(16) = m3pi2_squared;
  h(17) = RP.x();
  h(18) = RP.y(); // h(19) = RP.z();
  h(19) = p6piK.x();
  h(20) = p6piK.y();
  h(21) = p6piK.z();
  // h(22) = E6piK;
  h(22) = m6piK_squared;

  return h;
}

TVectorD x_initial_estimate( TVectorD m ) // Original initialisation for x (based on Anne Keune's thesis)
{
  // Builds an initial estimate for x based on the analytical calculations and on the known parameters in m

  ROOT::Math::XYZPoint PV( mprime(0), mprime(1), mprime(2) );
  ROOT::Math::XYZPoint DV1( mprime(3), mprime(4), mprime(5) );
  ROOT::Math::XYZVector p3pi1( mprime(6), mprime(7), mprime(8) );
  Double_t m3pi1_squared = mprime(9);
  ROOT::Math::XYZPoint DV2( mprime(10), mprime(11), mprime(12) );
  ROOT::Math::XYZVector p3pi2( mprime(13), mprime(14), mprime(15) );
  Double_t m3pi2_squared = mprime(16);
  ROOT::Math::XYZPoint RP( mprime(17), mprime(18), RPz );
  ROOT::Math::XYZVector p6piK( mprime(19), mprime(20), mprime(21) );
  Double_t m6piK_squared = mprime(22);

  Double_t E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  Double_t E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );
  Double_t E6piK = sqrt( m6piK_squared + p6piK.Mag2() );

  Double_t EK = E6piK - E3pi1 - E3pi2;
  ROOT::Math::XYZVector pK = p6piK - p3pi1 - p3pi2;

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  // Make transformation from LHCb to reference frame where the z-axis points along the K+ trajectory
  ROOT::Math::XYZPoint PV_t = makeTransformation_point( pK, RP, PV, false );
  ROOT::Math::XYZPoint DV1_t = makeTransformation_point( pK, RP, DV1, false );
  ROOT::Math::XYZVector p3pi1_t = makeTransformation_vec( pK, RP, p3pi1, false );
  ROOT::Math::XYZPoint DV2_t = makeTransformation_point( pK, RP, DV2, false );
  ROOT::Math::XYZVector p3pi2_t = makeTransformation_vec( pK, RP, p3pi2, false );
  ROOT::Math::XYZVector pK_t = makeTransformation_vec( pK, RP, pK, false );

  Double_t a1 = (DV1_t.y())/(DV1_t.x());
  Double_t a2 = (DV2_t.y())/(DV2_t.x());
  Double_t b = (PV_t.y() - a1*PV_t.x())/(a2*PV_t.x() - PV_t.y());
  Double_t c = b*(DV1_t.x())/(DV2_t.x());
  Double_t d = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()));
  Double_t e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
  Double_t f = ( PV_t.x()*sqrt(pK_t.Mag2()) )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
  Double_t g = c*e + d;
  Double_t h = f*c;
  Double_t i = DV1_t.z() - e*DV1_t.x();
  Double_t j = f*DV1_t.x();

  Double_t x1 = p3pi1_t.x() + a1*p3pi1_t.y() + e*p3pi1_t.z();
  Double_t x2 = b*p3pi2_t.x() + a2*b*p3pi2_t.y() + g*p3pi2_t.z();

  Double_t p1 = 1 + pow(a1,2) + pow(e,2) - pow(x1/E3pi1,2);
  Double_t p2 = 2*e*f - ( pow(mtau,2) + pow(m3pi1,2) + 2*f*p3pi1_t.z() )*(x1/pow(E3pi1,2) );
  Double_t p3 = pow(mtau,2) + pow(f,2) - pow( ( pow(mtau,2) + pow(m3pi1,2) + 2*f*p3pi1_t.z() )/(2*E3pi1), 2);
  Double_t q1 = pow(b,2) + pow(a2*b,2) + pow(g,2) - pow(x2/E3pi2,2);
  Double_t q2 = 2*g*h - ( pow(mtau,2) + pow(m3pi2,2) + 2*h*p3pi2_t.z() )*(x2/pow(E3pi2,2) );
  Double_t q3 =  pow(mtau,2) + pow(h,2) - pow( ( pow(mtau,2) + pow(m3pi2,2) + 2*h*p3pi2_t.z() )/(2*E3pi2),2 );

  Double_t Ptau1x_t = (p1*q3 - p3*q1)/(p2*q1 - p1*q2);
  Double_t Ptau1y_t = a1*Ptau1x_t;
  Double_t Ptau1z_t = e*Ptau1x_t + f;
  Double_t Ptau2x_t = b*Ptau1x_t;
  Double_t Ptau2y_t = a2*b*Ptau1x_t;
  Double_t Ptau2z_t = g*Ptau1x_t + h;
  Double_t BVz_t = i - j*(1/Ptau1x_t);

  ROOT::Math::XYZVector ptau1_t( Ptau1x_t, Ptau1y_t, Ptau1z_t );
  ROOT::Math::XYZVector ptau2_t( Ptau2x_t, Ptau2y_t, Ptau2z_t );
  ROOT::Math::XYZPoint BV_t(0, 0, BVz_t);

  // transform back to LHCb frame
  ROOT::Math::XYZVector ptau1 = makeTransformation_vec( pK, RP, ptau1_t, true );
  ROOT::Math::XYZVector ptau2 = makeTransformation_vec( pK, RP, ptau2_t, true );
  ROOT::Math::XYZPoint BV = makeTransformation_point( pK, RP, BV_t, true );

  // 3-momentum conservation in DV1, DV2 and BV
  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;

  // Neutrino and tau mass constraints 
  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  // Energy conservation in BV
  Double_t EB = Etau1 + Etau2 + EK;
  Double_t mB_squared = pow(EB,2) - pB.Mag2();

  Double_t L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  Double_t L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  Double_t L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );
  Double_t LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  x0(6) = mB_squared;
  x0(7) = ptau1.x();
  x0(8) = ptau1.y();
  x0(9) = ptau1.z();
  x0(10) = pnu1.x();
  x0(11) = pnu1.y();
  x0(12) = pnu1.z();
  x0(13) = ptau2.x();
  x0(14) = ptau2.y();
  x0(15) = ptau2.z();
  x0(16) = pnu2.x();
  x0(17) = pnu2.y();
  x0(18) = pnu2.z();
  x0(19) = L1;
  x0(20) = L2;
  x0(21) = L;
  x0(22) = LK;

  return x0;
}

TVectorD x_initial_estimate0( TVectorD mprime, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; taus direction based on vertices
{
  // Builds an initial estimate for x based on the Mariseille analytical calculations; it uses the offline estimate for BV  
  // (this initialisation does not apply the constraint: pB must point back to the PV) 

  ROOT::Math::XYZPoint PV( mprime(0), mprime(1), mprime(2) );
  ROOT::Math::XYZPoint DV1( mprime(3), mprime(4), mprime(5) );
  ROOT::Math::XYZVector p3pi1( mprime(6), mprime(7), mprime(8) ); 
  // Double_t E3pi1 = mprime(9); 
  Double_t m3pi1_squared = mprime(9); 
  ROOT::Math::XYZPoint DV2( mprime(10), mprime(11), mprime(12) );
  ROOT::Math::XYZVector p3pi2( mprime(13), mprime(14), mprime(15) );
  // Double_t E3pi2 = mprime(16); 
  Double_t m3pi2_squared = mprime(16); 
  ROOT::Math::XYZPoint RP( mprime(17), mprime(18), RPz );
  ROOT::Math::XYZVector p6piK( mprime(19), mprime(20), mprime(21) );
  // Double_t E6piK = mprime(22);
  Double_t m6piK_squared = mprime(22);

  Double_t E6piK = sqrt( m6piK_squared + p6piK.Mag2() );
  Double_t E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  Double_t E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );

  // Get kaon momentum
  ROOT::Math::XYZVector pK = p6piK - p3pi1 - p3pi2;
  Double_t EK = E6piK - E3pi1 - E3pi2;

  // Magnitude of 3pi momenta
  Double_t p3pi1_mag = p3pi1.r();
  Double_t p3pi2_mag = p3pi2.r();

  // B+ flight direction
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  // Get K+ momentum perpendicular to the B+ flight direction
  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  // Get tau flight directions (from vertices)
  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() );
  tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() );

  tau_dir1 = tau_dir1.Unit();
  tau_dir2 = tau_dir2.Unit();

  // Get tau direction unit vectors perpendicular to B+ flight direction
  ROOT::Math::XYZVector tau_dir1_perp = ( tau_dir1 - (tau_dir1.Dot(bDir))*bDir ).Unit();
  ROOT::Math::XYZVector tau_dir2_perp = ( tau_dir2 - (tau_dir2.Dot(bDir))*bDir ).Unit();

  // In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
  Double_t cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit());
  Double_t cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit());

  Double_t phi1 = TMath::ACos(cosphi1);
  Double_t phi2 = TMath::ACos(cosphi2);

  // In this place, get directions of tau momenta perpendicular to K+ momentum
  ROOT::Math::XYZVector tau_perp1_perpK = tau_dir1_perp - ( tau_dir1_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );
  ROOT::Math::XYZVector tau_perp2_perpK = tau_dir2_perp - ( tau_dir2_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );

  Double_t tau_perp_ratio = (tau_perp1_perpK.R())/(tau_perp2_perpK.R());

  // Calculate momentum component of taus in this plane
  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi1 + (cosphi2*(sin(phi1)/sin(phi2)))) );
  Double_t pMag_tau2_perp = pMag_tau1_perp*tau_perp_ratio;

  ROOT::Math::XYZVector p_tau1_perp = pMag_tau1_perp*tau_dir1_perp;
  ROOT::Math::XYZVector p_tau2_perp = pMag_tau2_perp*tau_dir2_perp;

  // Get angles made by tau directions with B+ flight direction
  Double_t tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit());
  Double_t tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit());

  // Get tau momenta parallel to B+ flight direction
  Double_t pMag_tau1_long = fabs(pMag_tau1_perp)*tau_B_cos1/sqrt( 1 - pow(tau_B_cos1,2) );
  Double_t pMag_tau2_long = fabs(pMag_tau2_perp)*tau_B_cos2/sqrt( 1 - pow(tau_B_cos2,2) );

  // Total tau momentum vector
  ROOT::Math::XYZVector ptau1 = p_tau1_perp + (pMag_tau1_long*bDir);
  ROOT::Math::XYZVector ptau2 = p_tau2_perp + (pMag_tau2_long*bDir);

  // Get the rest of x
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  Double_t EB = Etau1 + Etau2 + EK;
  Double_t mB_squared = pow(EB,2) - pB.Mag2();

  Double_t L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  Double_t L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  Double_t LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );
  Double_t L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  x0(6) = mB_squared;
  x0(7) = ptau1.x();
  x0(8) = ptau1.y();
  x0(9) = ptau1.z();
  x0(10) = pnu1.x();
  x0(11) = pnu1.y();
  x0(12) = pnu1.z();
  x0(13) = ptau2.x();
  x0(14) = ptau2.y();
  x0(15) = ptau2.z();
  x0(16) = pnu2.x();
  x0(17) = pnu2.y();
  x0(18) = pnu2.z();
  x0(19) = L1;
  x0(20) = L2;
  x0(21) = L;
  x0(22) = LK;

  return x0;
}

TVectorD x_initial_estimate1( TVectorD mprime, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; taus direction based visible 3pi momenta
{
  // Builds an initial estimate for x based on the Mariseille analytical calculations; it uses the offline estimate for BV  
  // (this initialisation does not apply the constraint: pB must point back to the PV) 

  ROOT::Math::XYZPoint PV( mprime(0), mprime(1), mprime(2) );
  ROOT::Math::XYZPoint DV1( mprime(3), mprime(4), mprime(5) );
  ROOT::Math::XYZVector p3pi1( mprime(6), mprime(7), mprime(8) ); 
  // Double_t E3pi1 = mprime(9); 
  Double_t m3pi1_squared = mprime(9); 
  ROOT::Math::XYZPoint DV2( mprime(10), mprime(11), mprime(12) );
  ROOT::Math::XYZVector p3pi2( mprime(13), mprime(14), mprime(15) );
  // Double_t E3pi2 = mprime(16); 
  Double_t m3pi2_squared = mprime(16); 
  ROOT::Math::XYZPoint RP( mprime(17), mprime(18), RPz );
  ROOT::Math::XYZVector p6piK( mprime(19), mprime(20), mprime(21) );
  // Double_t E6piK = mprime(22);
  Double_t m6piK_squared = mprime(22);

  Double_t E6piK = sqrt( m6piK_squared + p6piK.Mag2() );
  Double_t E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  Double_t E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );

  // Get kaon momentum
  ROOT::Math::XYZVector pK = p6piK - p3pi1 - p3pi2;
  Double_t EK = E6piK - E3pi1 - E3pi2;

  // Magnitude of 3pi momenta
  Double_t p3pi1_mag = p3pi1.r();
  Double_t p3pi2_mag = p3pi2.r();

  // B+ flight direction
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  // Get K+ momentum perpendicular to the B+ flight direction
  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  // Get tau flight directions (from vertices)
  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() );
  tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() );

  tau_dir1 = tau_dir1.Unit();
  tau_dir2 = tau_dir2.Unit();

  // Get tau direction unit vectors perpendicular to B+ flight direction
  ROOT::Math::XYZVector tau_dir1_perp = ( tau_dir1 - (tau_dir1.Dot(bDir))*bDir ).Unit();
  ROOT::Math::XYZVector tau_dir2_perp = ( tau_dir2 - (tau_dir2.Dot(bDir))*bDir ).Unit();

  // In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
  Double_t cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit());
  Double_t cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit());

  Double_t phi1 = TMath::ACos(cosphi1);
  Double_t phi2 = TMath::ACos(cosphi2);

  // In this place, get directions of tau momenta perpendicular to K+ momentum
  ROOT::Math::XYZVector tau_perp1_perpK = tau_dir1_perp - ( tau_dir1_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );
  ROOT::Math::XYZVector tau_perp2_perpK = tau_dir2_perp - ( tau_dir2_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );

  Double_t tau_perp_ratio = (tau_perp1_perpK.R())/(tau_perp2_perpK.R());

  // Calculate momentum component of taus in this plane
  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi1 + (cosphi2*(sin(phi1)/sin(phi2)))) );
  Double_t pMag_tau2_perp = pMag_tau1_perp*tau_perp_ratio;

  ROOT::Math::XYZVector p_tau1_perp = pMag_tau1_perp*tau_dir1_perp;
  ROOT::Math::XYZVector p_tau2_perp = pMag_tau2_perp*tau_dir2_perp;

  // Get angles made by tau directions with B+ flight direction
  Double_t tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit());
  Double_t tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit());

  // Get tau momenta parallel to B+ flight direction
  Double_t pMag_tau1_long = fabs(pMag_tau1_perp)*tau_B_cos1/sqrt( 1 - pow(tau_B_cos1,2) );
  Double_t pMag_tau2_long = fabs(pMag_tau2_perp)*tau_B_cos2/sqrt( 1 - pow(tau_B_cos2,2) );

  // Total tau momentum vector
  ROOT::Math::XYZVector ptau1 = p_tau1_perp + (pMag_tau1_long*bDir);
  ROOT::Math::XYZVector ptau2 = p_tau2_perp + (pMag_tau2_long*bDir);

  // Get the rest of x
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  Double_t EB = Etau1 + Etau2 + EK;
  Double_t mB_squared = pow(EB,2) - pB.Mag2();

  Double_t L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  Double_t L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  Double_t LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );
  Double_t L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  x0(6) = mB_squared;
  x0(7) = ptau1.x();
  x0(8) = ptau1.y();
  x0(9) = ptau1.z();
  x0(10) = pnu1.x();
  x0(11) = pnu1.y();
  x0(12) = pnu1.z();
  x0(13) = ptau2.x();
  x0(14) = ptau2.y();
  x0(15) = ptau2.z();
  x0(16) = pnu2.x();
  x0(17) = pnu2.y();
  x0(18) = pnu2.z();
  x0(19) = L1;
  x0(20) = L2;
  x0(21) = L;
  x0(22) = LK;

  return x0;
}

TVectorD x_initial_estimate2( TVectorD mprime, ROOT::Math::XYZPoint BV ) 
{
  // Builds an initial estimate for x based on the Mariseille analytical calculations; it uses the offline estimate for BV  
  // (this initialisation does not apply the constraint: pB must point back to the PV) 

  ROOT::Math::XYZPoint PV( mprime(0), mprime(1), mprime(2) );
  ROOT::Math::XYZPoint DV1( mprime(3), mprime(4), mprime(5) );
  ROOT::Math::XYZVector p3pi1( mprime(6), mprime(7), mprime(8) ); 
  // Double_t E3pi1 = mprime(9); 
  Double_t m3pi1_squared = mprime(9); 
  ROOT::Math::XYZPoint DV2( mprime(10), mprime(11), mprime(12) );
  ROOT::Math::XYZVector p3pi2( mprime(13), mprime(14), mprime(15) );
  // Double_t E3pi2 = mprime(16); 
  Double_t m3pi2_squared = mprime(16); 
  ROOT::Math::XYZPoint RP( mprime(17), mprime(18), RPz );
  ROOT::Math::XYZVector p6piK( mprime(19), mprime(20), mprime(21) );
  // Double_t E6piK = mprime(22);
  Double_t m6piK_squared = mprime(22);

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  Double_t E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  Double_t E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );
  Double_t E6piK = sqrt( m6piK_squared + p6piK.Mag2() );

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  // Use the maximum value of theta: neutrino takes the maximum portion of momentum from the tau
  Double_t theta1 = asin( ( pow(mtau,2) - pow(m3pi1,2) )/( 2*mtau*sqrt( p3pi1.Mag2() ) ) );
  Double_t theta2 = asin( ( pow(mtau,2) - pow(m3pi2,2) )/( 2*mtau*sqrt( p3pi2.Mag2() ) ) );
  
  Double_t ptau1_mag = ( (pow(mtau,2) + pow(m3pi1,2))*sqrt(p3pi1.Mag2())*cos(theta1) )/( 2*( pow(E3pi1,2) - p3pi1.Mag2()*pow(cos(theta1),2) ) );
  Double_t ptau2_mag = ( (pow(mtau,2) + pow(m3pi2,2))*sqrt(p3pi2.Mag2())*cos(theta2) )/( 2*( pow(E3pi2,2) - p3pi2.Mag2()*pow(cos(theta2),2) ) );

  ROOT::Math::XYZVector ptau1 = ptau1_mag*u1;
  ROOT::Math::XYZVector ptau2 = ptau1_mag*u2;

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  ROOT::Math::XYZVector pB = p6piK + pnu1 + pnu2;
  Double_t EB = E6piK + Enu1 + Enu2;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

  ROOT::Math::XYZVector pK = p6piK - p3pi1 - p3pi2;

  // ROOT::Math::XYZVector pK = p6piK - p3pi1 - p3pi2;
  Double_t L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  Double_t L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  Double_t LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );
  Double_t L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  x0(6) = MB_squared;
  x0(7) = ptau1.x();
  x0(8) = ptau1.y();
  x0(9) = ptau1.z();
  x0(10) = pnu1.x();
  x0(11) = pnu1.y();
  x0(12) = pnu1.z();
  x0(13) = ptau2.x();
  x0(14) = ptau2.y();
  x0(15) = ptau2.z();
  x0(16) = pnu2.x();
  x0(17) = pnu2.y();
  x0(18) = pnu2.z();
  x0(19) = L1;
  x0(20) = L2;
  x0(21) = L;
  x0(22) = LK;

  return x0;
}

ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag)
{
  Double_t deltaX = refPoint.x();
  Double_t deltaY = refPoint.y();
  Double_t deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  Double_t phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  Double_t theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

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
  Double_t deltaX = refPoint.x();
  Double_t deltaY = refPoint.y();
  Double_t deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  Double_t phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  Double_t theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

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

TVectorD transform_m( TVectorD m )
{
  TVectorD mp( dimM );

  mp(0) = m(0); // PVx
  mp(1) = m(1); // PVy
  mp(2) = m(2); // PVz
  mp(3) = m(3); // DV1x
  mp(4) = m(4); // DV1y
  mp(5) = m(5); // DV1z

  ROOT::Math::XYZVector p1( m(6), m(7), m(8) );
  ROOT::Math::XYZVector p2( m(9), m(10), m(11) );
  ROOT::Math::XYZVector p3( m(12), m(13), m(14) );

  mp(6) = p1.x() + p2.x() + p3.x(); // p3pi1x
  mp(7) = p1.y() + p2.y() + p3.y(); // p3pi1y
  mp(8) = p1.z() + p2.z() + p3.z(); // p3pi1z

  Double_t E1 = sqrt( pow(mpion,2) + p1.Mag2() );
  Double_t E2 = sqrt( pow(mpion,2) + p2.Mag2() );
  Double_t E3 = sqrt( pow(mpion,2) + p3.Mag2() );
  
  ROOT::Math::PxPyPzEVector P1( p1.x(), p1.y(), p1.z(), E1 );
  ROOT::Math::PxPyPzEVector P2( p2.x(), p2.y(), p2.z(), E2 );
  ROOT::Math::PxPyPzEVector P3( p3.x(), p3.y(), p3.z(), E3 );

  mp(9) = 3*pow(mpion,2) + 2*P1.Dot(P2) + 2*P1.Dot(P3) + 2*P2.Dot(P3) ; // m3pi1^2

  mp(10) = m(15); // DV2x
  mp(11) = m(16); // DV2y
  mp(12) = m(17); // DV2z

  ROOT::Math::XYZVector p4( m(18), m(19), m(20) );
  ROOT::Math::XYZVector p5( m(21), m(22), m(23) );
  ROOT::Math::XYZVector p6( m(24), m(25), m(26) );

  mp(13) = p4.x() + p5.x() + p6.x(); // p3pi2x
  mp(14) = p4.y() + p5.y() + p6.y(); // p3pi2y
  mp(15) = p4.z() + p5.z() + p6.z(); // p3pi2z

  Double_t E4 = sqrt( pow(mpion,2) + p4.Mag2() );
  Double_t E5 = sqrt( pow(mpion,2) + p5.Mag2() );
  Double_t E6 = sqrt( pow(mpion,2) + p6.Mag2() );

  ROOT::Math::PxPyPzEVector P4( p4.x(), p4.y(), p4.z(), E4 );
  ROOT::Math::PxPyPzEVector P5( p5.x(), p5.y(), p5.z(), E5 );
  ROOT::Math::PxPyPzEVector P6( p6.x(), p6.y(), p6.z(), E6 );

  mp(16) = 3*pow(mpion,2) + 2*P4.Dot(P5) + 2*P4.Dot(P6) + 2*P5.Dot(P6); // m3pi2^2

  mp(17) = m(27); // RPx
  mp(18) = m(28); // RPy

  ROOT::Math::XYZVector pK( m(29), m(30), m(31) );

  mp(19) = p1.x() + p2.x() + p3.x() + p4.x() + p5.x() + p6.x() + pK.x(); // p6piKx
  mp(20) = p1.y() + p2.y() + p3.y() + p4.y() + p5.y() + p6.y() + pK.y(); // p6piKy
  mp(21) = p1.z() + p2.z() + p3.z() + p4.z() + p5.z() + p6.z() + pK.z(); // p6piKz

  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );
  ROOT::Math::PxPyPzEVector PK( pK.x(), pK.y(), pK.z(), EK );

  mp(22) = pow(mkaon,2) + 6*pow(mpion,2) + 2*( P1.Dot(P2) + P1.Dot(P3) + P2.Dot(P3) + P4.Dot(P5) + P4.Dot(P6) + P5.Dot(P6) + P1.Dot(P4) + P1.Dot(P5) + P1.Dot(P6) + P2.Dot(P4) + P2.Dot(P5) + P2.Dot(P6) + P3.Dot(P4) + P3.Dot(P5) + P3.Dot(P6) + P1.Dot(PK) + P2.Dot(PK) + P3.Dot(PK) + P4.Dot(PK) + P5.Dot(PK) + P6.Dot(PK) )  ; // m6piK^2

  return mp;

}

TMatrixDSym transform_V( TVectorD m, TMatrixDSym V, TMatrixDSym taup_cov, TMatrixDSym taum_cov )
{
  TMatrixD Vp( dimM, dimM );

  ROOT::Math::XYZVector p1( m(6), m(7), m(8) );
  ROOT::Math::XYZVector p2( m(9), m(10), m(11) );
  ROOT::Math::XYZVector p3( m(12), m(13), m(14) );

  ROOT::Math::XYZVector p4( m(18), m(19), m(20) );
  ROOT::Math::XYZVector p5( m(21), m(22), m(23) );
  ROOT::Math::XYZVector p6( m(24), m(25), m(26) );

  ROOT::Math::XYZVector pK( m(29), m(30), m(31) );

  Double_t E1 = sqrt( pow(mpion,2) + p1.Mag2() );
  Double_t E2 = sqrt( pow(mpion,2) + p2.Mag2() );
  Double_t E3 = sqrt( pow(mpion,2) + p3.Mag2() );
  Double_t E4 = sqrt( pow(mpion,2) + p4.Mag2() );
  Double_t E5 = sqrt( pow(mpion,2) + p5.Mag2() );
  Double_t E6 = sqrt( pow(mpion,2) + p6.Mag2() );
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  // Jacobian matrix
  TMatrixD J( dimM, dimM_before );
  for(int i = 0; i < dimM; i++)
  {
    for(int j = 0; j < dimM_before; j++)
    {
      J(i,j) = 0.;
    }
  }

  // PV
  J(0,0) = 1.; // PVx
  J(1,1) = 1.; // PVy
  J(2,2) = 1.; // PVz

  // DV1
  J(3,3) = 1.; // DV1x
  J(4,4) = 1.; // DV1y
  J(5,5) = 1.; // DV1z

  // DV2
  J(10,15) = 1.; // DV2x
  J(11,16) = 1.; // DV2y
  J(12,17) = 1.; // DV2z

  // RP
  J(17,27) = 1.; // RPx
  J(18,28) = 1.; // RPy

  // p3pi1x
  J(6,6) = 1.; // p1x
  J(6,9) = 1.; // p2x
  J(6,12) = 1.; // p3x

  // p3pi1y
  J(7,7) = 1.; // p1y
  J(7,10) = 1.; // p2y
  J(7,13) = 1.; // p3y

  // p3pi1z
  J(8,8) = 1.; // p1z
  J(8,11) = 1.; // p3y
  J(8,14) = 1.; // p3z

  // m3pi1^2
  J(9,6) = 2*( ((E2+E3)/E1)*p1.x() - p2.x() - p3.x() ); // p1x
  J(9,7) = 2*( ((E2+E3)/E1)*p1.y() - p2.y() - p3.y() ); // p1y
  J(9,8) = 2*( ((E2+E3)/E1)*p1.z() - p2.z() - p3.z() ); // p1z
  J(9,9) = 2*( ((E1+E3)/E2)*p2.x() - p1.x() - p3.x() ); // p2x
  J(9,10) = 2*( ((E1+E3)/E2)*p2.y() - p1.y() - p3.y() ); // p2y
  J(9,11) = 2*( ((E1+E3)/E2)*p2.z() - p1.z() - p3.z() ); // p2z
  J(9,12) = 2*( ((E1+E2)/E3)*p3.x() - p1.x() - p2.x() ); // p3x
  J(9,13) = 2*( ((E1+E2)/E3)*p3.y() - p1.y() - p2.y() ); // p3y
  J(9,14) = 2*( ((E1+E2)/E3)*p3.z() - p1.z() - p2.z() ); // p3z

  // p3pi2x
  J(13,18) = 1.; // p4x
  J(13,21) = 1.; // p5x
  J(13,24) = 1.; // p6x

  // p3pi1y
  J(14,19) = 1.; // p4y
  J(14,22) = 1.; // p5y
  J(14,25) = 1.; // p6y

  // p3pi1z
  J(15,20) = 1.; // p4z
  J(15,23) = 1.; // p5z
  J(15,26) = 1.; // p6z

  // m3pi2^2
  J(16,18) = 2*( ((E5+E6)/E4)*p4.x() - p5.x() - p6.x() ); // p4x
  J(16,19) = 2*( ((E5+E6)/E4)*p4.y() - p5.y() - p6.y() ); // p4y
  J(16,20) = 2*( ((E5+E6)/E4)*p4.z() - p5.z() - p6.z() ); // p4z
  J(16,21) = 2*( ((E4+E6)/E5)*p5.x() - p4.x() - p6.x() ); // p5x
  J(16,22) = 2*( ((E4+E6)/E5)*p5.y() - p4.y() - p6.y() ); // p5y
  J(16,23) = 2*( ((E4+E6)/E5)*p5.z() - p4.z() - p6.z() ); // p5z
  J(16,24) = 2*( ((E4+E5)/E6)*p6.x() - p4.x() - p5.x() ); // p6x
  J(16,25) = 2*( ((E4+E5)/E6)*p6.y() - p4.y() - p5.y() ); // p6y
  J(16,26) = 2*( ((E4+E5)/E6)*p6.z() - p4.z() - p5.z() ); // p6z

  // p6piKx
  J(19,6) = 1.; // p1x
  J(19,9) = 1.; // p2x
  J(19,12) = 1.; // p3x
  J(19,18) = 1.; // p4x 
  J(19,21) = 1.; // p5x
  J(19,24) = 1.; // p6x
  J(19,29) = 1.; // pKx

  // p6piKy
  J(20,7) = 1.; // p1y
  J(20,10) = 1.; // p2y
  J(20,13) = 1.; // p3y
  J(20,19) = 1.; // p4y
  J(20,22) = 1.; // p5y
  J(20,25) = 1.; // p6y
  J(20,30) = 1.; // pKy

  // p6piKz
  J(21,8) = 1.; // p1z
  J(21,11) = 1.; // p2z
  J(21,14) = 1.; // p3z
  J(21,20) = 1.; // p4z
  J(21,23) = 1.; // p5z
  J(21,26) = 1.; // p6z
  J(21,31) = 1.; // pKz

  // m6piK^2
  J(22,6) = 2*( ((E2+E3+E4+E5+E6+EK)/E1)*p1.x() - p2.x() - p3.x() - p4.x() - p5.x() - p6.x() - pK.x() ); // p1x
  J(22,7) = 2*( ((E2+E3+E4+E5+E6+EK)/E1)*p1.y() - p2.y() - p3.y() - p4.y() - p5.y() - p6.y() - pK.y() ); // p1y
  J(22,8) = 2*( ((E2+E3+E4+E5+E6+EK)/E1)*p1.z() - p2.z() - p3.z() - p4.z() - p5.z() - p6.z() - pK.z() ); // p1z
  J(22,9) = 2*( ((E1+E3+E4+E5+E6+EK)/E2)*p2.x() - p1.x() - p3.x() - p4.x() - p5.x() - p6.x() - pK.x() ); // p2x
  J(22,10) = 2*( ((E1+E3+E4+E5+E6+EK)/E2)*p2.y() - p1.y() - p3.y() - p4.y() - p5.y() - p6.y() - pK.y() ); // p2y
  J(22,11) = 2*( ((E1+E3+E4+E5+E6+EK)/E2)*p2.z() - p1.z() - p3.z() - p4.z() - p5.z() - p6.z() - pK.z() ); // p2z
  J(22,12) = 2*( ((E1+E2+E4+E5+E6+EK)/E3)*p3.x() - p1.x() - p2.x() - p4.x() - p5.x() - p6.x() - pK.x() ); // p3x
  J(22,13) = 2*( ((E1+E2+E4+E5+E6+EK)/E3)*p3.y() - p1.y() - p2.y() - p4.y() - p5.y() - p6.y() - pK.y() ); // p3y
  J(22,14) = 2*( ((E1+E2+E4+E5+E6+EK)/E3)*p3.z() - p1.z() - p2.z() - p4.z() - p5.z() - p6.z() - pK.z() ); // p3z
  J(22,18) = 2*( ((E1+E2+E3+E5+E6+EK)/E4)*p4.x() - p1.x() - p2.x() - p3.x() - p5.x() - p6.x() - pK.x() ); // p4x
  J(22,19) = 2*( ((E1+E2+E3+E5+E6+EK)/E4)*p4.x() - p1.y() - p2.y() - p3.y() - p5.y() - p6.y() - pK.y() ); // p4y
  J(22,20) = 2*( ((E1+E2+E3+E5+E6+EK)/E4)*p4.z() - p1.z() - p2.z() - p3.z() - p5.z() - p6.z() - pK.z() ); // p4z
  J(22,21) = 2*( ((E1+E2+E3+E4+E6+EK)/E5)*p5.x() - p1.x() - p2.x() - p3.x() - p4.x() - p6.x() - pK.x() ); // p5x
  J(22,22) = 2*( ((E1+E2+E3+E4+E6+EK)/E5)*p5.y() - p1.y() - p2.y() - p3.y() - p4.y() - p6.y() - pK.y() ); // p5y
  J(22,23) = 2*( ((E1+E2+E3+E4+E6+EK)/E5)*p5.z() - p1.z() - p2.z() - p3.z() - p4.z() - p6.z() - pK.z() ); // p5z
  J(22,24) = 2*( ((E1+E2+E3+E4+E5+EK)/E6)*p6.x() - p1.x() - p2.x() - p3.x() - p4.x() - p5.x() - pK.x() ); // p6x
  J(22,25) = 2*( ((E1+E2+E3+E4+E5+EK)/E6)*p6.y() - p1.y() - p2.y() - p3.y() - p4.y() - p5.y() - pK.y() ); // p6y
  J(22,26) = 2*( ((E1+E2+E3+E4+E5+EK)/E6)*p6.z() - p1.z() - p2.z() - p3.z() - p4.z() - p5.z() - pK.z() ); // p6z
  J(22,29) = 2*( ((E1+E2+E3+E4+E5+E6)/EK)*pK.x() - p1.x() - p2.x() - p3.x() - p4.x() - p5.x() - p6.x() ); // pKx
  J(22,30) = 2*( ((E1+E2+E3+E4+E5+E6)/EK)*pK.y() - p1.y() - p2.y() - p3.y() - p4.y() - p5.y() - p6.y() ); // pKy
  J(22,31) = 2*( ((E1+E2+E3+E4+E5+E6)/EK)*pK.z() - p1.z() - p2.z() - p3.z() - p4.z() - p5.z() - p6.z() ); // pKz

  // J.Print();

  TMatrixD J_transpose = J;
  J_transpose.T();

  Vp = J*(V*J_transpose);

  // DV vs p3pi correlations
  // for(int i = 0; i < 3; i++)
  // {
  //   for(int j = 0; j < 4; j++)
  //   {
  //     Vp(i+3,j+6) = taup_cov(i,j+3);
  //     Vp(i+10,j+13) = taum_cov(i,j+3);
  //   }
  // }

  // for(int i = 0; i < 4; i++)
  // {
  //   for(int j = 0; j < 3; j++)
  //   {
  //     Vp(i+6,j+3) = taup_cov(i+3,j);
  //     Vp(i+13,j+10) = taum_cov(i+3,j);
  //   }
  // }

  TMatrixDSym Vp_sym(dimM);
  Vp_sym.SetTol(pow(10,-23));
  for(int i = 0; i < dimM; i++)
  {
    for(int j = 0; j < dimM; j++)
    {
      Vp_sym(i,j) = Vp(i,j);
    }
  }

  return Vp_sym;
}