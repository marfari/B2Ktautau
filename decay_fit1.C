
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
Double_t chi2, MB, MB_err;
TVectorD x_m(dimX);
TMatrixD C_m(dimX,dimX);
Int_t status;
Int_t tikho;

Double_t mtau = 1776.86; // MeV (PDG 2023)
Double_t mK = 493.677; // MeV (PDG 2023) 

// Functions
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);

TVectorD x_initial_estimate( TVectorD m );
TVectorD x_initial_estimate2( TVectorD m, ROOT::Math::XYZPoint BV );

TVectorD h( TVectorD x );
TMatrixD dh_dx( TVectorD x );
Double_t chisquare( const Double_t* x_values );

void minimize( ROOT::Math::XYZPoint BV, int init );

void decay_fit1(int year, TString MC_files, int component, int line)
{
  TFileCollection* fc = new TFileCollection("MC", "MC", MC_files, 1, line);
  TChain* t = new TChain("DecayTree");
  t->AddFileInfoList((TCollection*)fc->GetList());

  // TFile* f = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2018/DVntuple_MC_2018_MagUp_visible.root");
  // TTree* t = (TTree*)f->Get("ntuple/DecayTree");
  // t->AddFriend("DecayTree", "/panfs/felician/B2Ktautau/ROOT_Sim/2018/transformed_m_V.root");

  UInt_t num_entries = t->GetEntries(); 
  Double_t m_vars[dimM];
  Double_t V_vars[dimM][dimM];
  Double_t BVx, BVy, BVz;

  for(int i = 1; i <= dimM; i++)
  {
      t->SetBranchAddress(Form("df_m_%i",i), &m_vars[i-1] );
      for(int j = 1; j <= dimM; j++)
      {
          t->SetBranchAddress(Form("df_V_%i_%i",i,j), &V_vars[i-1][j-1]);
      }
  }

  // int i_list[] = {4,5,6, 7,8,9,10, 11,12,13, 14,15,16,17, 20,21,22,23};
  // int N = sizeof(i_list)/sizeof(i_list[0]);
  // Double_t m_prime_vars[N];
  // Double_t V_prime_vars[N][N];
  // for(int i = 0; i < N; i++){
  //   t->SetBranchAddress(Form("df_m_prime_%i",i_list[i]), &m_prime_vars[i]);
  //   for(int j = 0; j < N; j++)
  //   {
  //     t->SetBranchAddress(Form("df_V_prime_%i_%i",i_list[i],i_list[j]), &V_prime_vars[i][j]);
  //   }
  // }

  t->SetBranchAddress("df_RPz", &RPz);
  t->SetBranchAddress("Bp_ENDVERTEX_X", &BVx);
  t->SetBranchAddress("Bp_ENDVERTEX_Y", &BVy);
  t->SetBranchAddress("Bp_ENDVERTEX_Z", &BVz);

  TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Component_%i/fit_result_%i.root",year,component,line), "RECREATE");
  TTree* tout = new TTree("DecayTree", "DecayTree");

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
  tout->Branch("df_tikho", &tikho);

  // Loop over events
  for(int evt = 0; evt < num_entries; evt++)
  {
    t->GetEntry(evt);
    
    // 1)  Measured parameters, m, and their covariance matrix V
    // m = {PV, DV1, p3pi1, m3pi1^2, DV2, p3pi2, m3pi2^2, RP_T, pK}
    for(int i = 0; i < dimM; i++)
    {
        m(i) = m_vars[i];
        for(int j = 0; j < dimM; j++)
        {
          V(i,j) = V_vars[i][j];
        }
    }  
    // m.Print();
    // V.Print();

    // Compute eigenvalues of the 3pi+ and 3pi- momentum sub-matrices:
    TMatrixD V_3pi1(7,7);
    for(int i = 0; i < 7; i++)
    {
      for(int j = 0; j < 7; j++)
      {
        V_3pi1(i,j) = V(i+3,j+3);
      }
    }
    // V_3pi1.Print();
    TVectorD V_eigen_values_3pi1(dimM);  
    TMatrixD V_eigen_vectors_3pi1 = V_3pi1.EigenVectors(V_eigen_values_3pi1);
    // V_eigen_values_3pi1.Print();

    TMatrixD V_3pi2(7,7);
    for(int i = 0; i < 7; i++)
    {
      for(int j = 0; j < 7; j++)
      {
        V_3pi2(i,j) = V(i+10,j+10);
      }
    }
    // V_3pi2.Print();
    TVectorD V_eigen_values_3pi2(dimM);  
    TMatrixD V_eigen_vectors_3pi2 = V_3pi2.EigenVectors(V_eigen_values_3pi2);
    // V_eigen_values_3pi2.Print();

    // Check if there are negative eigenvalues and if so apply Tikhonov regularisation to the sub-matrix
    Double_t lambda1 = 0.;
    Double_t lambda2 = 0.;
    for(int i = 0; i < 7; i++)
    {
      if( V_eigen_values_3pi1(i) < 0 )
      {
        lambda1 = 1.5*abs( V_eigen_values_3pi1(i) );
      }
    }

    for(int i = 0; i < 7; i++)
    {
      if( V_eigen_values_3pi2(i) < 0 )
      {
        lambda2 = 1.5*abs( V_eigen_values_3pi2(i) );
      }
    }

    for(int i = 0; i < 7; i++)
    {
      // V(i+3,i+3) += lambda1;
      // V(i+10,i+10) += lambda2;
      V(i+6,i+6) += lambda1;
      V(i+13,i+13) += lambda2;

      V_3pi1(i,i) += lambda1;
      V_3pi2(i,i) += lambda2;
    }
    // cout << "l1 = " << lambda1 << endl;
    // cout << "l2 = " << lambda2 << endl;
    // V_3pi1.Print();
    // V_3pi2.Print();

    if( (lambda1 == 0.) && (lambda2 == 0.) )
    {
      tikho = 0;
    }
    else if( (lambda1 != 0.) && (lambda2 != 0.) )
    {
      tikho = 2;
    }
    else
    {
      tikho = 1;
    }
    // V.Print();

    // TVectorD V_eigen_values(dimM);  
    // TMatrixD V_eigen_vectors = V.EigenVectors(V_eigen_values);
    // V_eigen_values.Print();

    // Correlation matrix of m (before change of variables)
    // TMatrixD corr_E(dimM,dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //   for(int j = 0; j < dimM; j++)
    //   {
    //     corr_E(i,j) = V(i,j)/(sqrt(V(i,i))*sqrt(V(j,j)));
    //   }
    // }
    // corr_E.Print();

    // Make change of variables E -> m^2 in m
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
    // Jacob.Print();
    // Jacob_transpose.T();
    // Vprime = Jacob*(V*Jacob_transpose); // transform covariance matrix V' = J V J^T (approximation to 1st order in Taylor expansion)
    // Vprime(9,9) = abs(Vprime(9,9)); // force this to be positive for now; why is it negative in the 1st place?
    // Vprime(16,16) = abs(Vprime(16,16));
    // Vprime(22,22) = abs(Vprime(22,22));
    // m.Print();
    // Vprime.Print();

    // Update m and V with change of variables transformation
    // m(9) = m_prime_vars[6];
    // m(16) = m_prime_vars[13];
    // m(22) = m_prime_vars[17];

    // for(int i = 0; i < N; i++)
    // {
    //   for(int j = 0; j < N; j++)
    //   {
    //     V(i_list[i]-1,i_list[j]-1) = V_prime_vars[i][j];
    //   }
    // }

    // Correlation matrix of m (after change of variables)
    // TMatrixD corr_m2(dimM,dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //   for(int j = 0; j < dimM; j++)
    //   {
    //     corr_m2(i,j) = Vprime(i,j)/(sqrt(Vprime(i,i))*sqrt(Vprime(j,j)));
    //   }
    // }
    // corr_m2.Print();
    // return;

    // Weights matrix W = V'^-1
    W = V;
    W.Invert();
    // W.Print();

    ROOT::Math::XYZPoint BV( BVx, BVy, BVz ); // BV (necessary for Marseille initialisation)
    minimize( BV, 2 );

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
  //   init = 2 : uses the Marseille initialisation; it uses the BV offline estimate and assumes that the neutrino takes as much momentum as possible form the tau 
  // 2) It sets up the minimizer: defining the initial values and bounds of the x parameters  
  // 3) It builds the chi^2 function with x0 by calling the function chisquare 
  // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

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
    x0_vars[i] = x0(i); // x0_vars is a Double_t, x0 is a TVectorD; the function to minimise must receive a Double_t as input
    x0_err[i] = 0.1*abs(x0_vars[i]); // initial errors on x0; they are used as the first step size in the minimisation; considering a 10% error for now
  }

  // x0.Print();
  // TVectorD hx = h(x0);
  // hx.Print();
  // TVectorD r = m-hx;
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

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");    
  ROOT::Math::Functor func(&chisquare, dimX);

  min->SetMaxIterations(1000000);  
  min->SetPrintLevel(1);
  min->SetPrecision(pow(10,-15));
  // min->SetTolerance(0.001);
  min->SetMaxFunctionCalls(100000000);
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
  if( sizeof(name)/sizeof(name[0]) != dimX ){ 
    cout << "ERROR: name does not have the proper dimension" << endl;
    return; 
  }

  for(int i = 0; i < dimM; i++)
  {
    min->SetVariable(i, name[i], x0[i], x0_err[i]);
    min->SetVariableLimits(i, x0_min[i], x0_max[i]);
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
  // TVectorD r_m = m-hx_m;
  // r_m.Print();
  // C_m.Print();

  // Calculate MB and its error
  // Double_t MB_squared_m = x_m(6);
  ROOT::Math::XYZVector ptau1_m( x_m(6), x_m(7), x_m(8) );
  ROOT::Math::XYZVector pnu1_m( x_m(9), x_m(10), x_m(11) );
  ROOT::Math::XYZVector ptau2_m( x_m(12), x_m(13), x_m(14) );
  ROOT::Math::XYZVector pnu2_m( x_m(15), x_m(16), x_m(17) );
  Double_t Etau1_m = sqrt( pow(mtau,2) + ptau1_m.Mag2() );
  Double_t Enu1_m = sqrt( pnu1_m.Mag2() );
  Double_t Etau2_m = sqrt( pow(mtau,2) + ptau2_m.Mag2() );
  Double_t Enu2_m = sqrt( pnu2_m.Mag2() );

  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  Double_t EK = sqrt( pow(mK,2) + pK.Mag2() );

  // // Energy convervation in BV
  Double_t EB = Etau1_m + Etau2_m + EK;
  ROOT::Math::XYZVector pB_m( x_m(3), x_m(4), x_m(5) );
  Double_t MB_squared_m = pow(EB,2) - pB_m.Mag2();

  if( MB_squared_m > 0 )
  {
    MB = sqrt( MB_squared_m );
  }
  else
  {
    MB = -sqrt( abs(MB_squared_m) );
  }
  // Double_t dMB_squared = sqrt(C_m(6,6));
  // MB_err = dMB_squared/(2*MB);

  Double_t dptau1x = C_m(6,6);
  Double_t dptau1y = C_m(7,7);
  Double_t dptau1z = C_m(8,8);
  Double_t dptau1xy = C_m(6,7);
  Double_t dptau1xz = C_m(6,8);
  Double_t dptau1yz = C_m(7,8);
  Double_t dEtau1 = sqrt( pow(ptau1_m.x()/Etau1_m,2)*dptau1x + pow(ptau1_m.y()/Etau1_m,2)*dptau1y + pow(ptau1_m.z()/Etau1_m,2)*dptau1z + 2*( (ptau1_m.x()*ptau1_m.y())/pow(Etau1_m,2) )*dptau1xy + 2*( (ptau1_m.x()*ptau1_m.z())/pow(Etau1_m,2) )*dptau1xz + 2*( (ptau1_m.y()*ptau1_m.z())/pow(Etau1_m,2) )*dptau1yz );
  cout << "Etau1 = " << Etau1_m << " +/- " << dEtau1 << endl;

  Double_t dptau2x = C_m(12,12);
  Double_t dptau2y = C_m(13,13);
  Double_t dptau2z = C_m(14,14);
  Double_t dptau2xy = C_m(12,13);
  Double_t dptau2xz = C_m(12,14);
  Double_t dptau2yz = C_m(13,14);
  Double_t dEtau2 = sqrt( pow(ptau2_m.x()/Etau2_m,2)*dptau2x + pow(ptau2_m.y()/Etau2_m,2)*dptau2y + pow(ptau2_m.z()/Etau2_m,2)*dptau2z + 2*( (ptau2_m.x()*ptau2_m.y())/pow(Etau2_m,2) )*dptau2xy + 2*( (ptau2_m.x()*ptau2_m.z())/pow(Etau2_m,2) )*dptau2xz + 2*( (ptau2_m.y()*ptau2_m.z())/pow(Etau2_m,2) )*dptau2yz );
  cout << "Etau2 = " << Etau2_m << " +/- " << dEtau2 << endl;

  Double_t dpKx = V(19,19);
  Double_t dpKy = V(20,20);
  Double_t dpKz = V(21,21);
  Double_t dpKxy = V(19,20);
  Double_t dpKxz = V(19,21);
  Double_t dpKyz = V(20,21);
  Double_t dEK = sqrt( pow(pK.x()/EK,2)*dpKx + pow(pK.y()/EK,2)*dpKy + pow(pK.z()/EK,2)*dpKz + 2*( (pK.x()*pK.y())/pow(EK,2) )*dpKxy + 2*( (pK.x()*pK.z())/pow(EK,2) )*dpKxz + 2*( (pK.y()*pK.z())/pow(EK,2) )*dpKyz );
  cout << "EK = " << EK << " +/- " << dEK << endl;

  Double_t dEB = sqrt( pow(dEtau1,2) + pow(dEtau2,2) + pow(dEK,2) );
  cout << "EB = " << EB << " +/- " << dEB << endl;

  // EB + pB,pB
  Double_t dpBx = C_m(3,3);
  Double_t dpBy = C_m(4,4);
  Double_t dpBz = C_m(5,5);
  Double_t dpBxy = C_m(3,4);
  Double_t dpBxz = C_m(3,5);
  Double_t dpByz = C_m(4,5);

  Double_t EB_pBpB = pow(2*EB,2)*dEB + pow(2*pB_m.x(),2)*dpBx + pow(2*pB_m.y(),2)*dpBy + pow(2*pB_m.z(),2)*dpBz + 8*pB_m.x()*pB_m.y()*dpBxy + 8*pB_m.x()*pB_m.z()*dpBxz + 8*pB_m.y()*pB_m.z()*dpByz;

  // ptau1,pB
  Double_t dpt1xBx = C_m(3,6);
  Double_t dpt1xBy = C_m(4,6);
  Double_t dpt1xBz = C_m(5,6);
  Double_t dpt1yBx = C_m(3,7);
  Double_t dpt1yBy = C_m(4,7);
  Double_t dpt1yBz = C_m(5,7);
  Double_t dpt1zBx = C_m(3,8);
  Double_t dpt1zBy = C_m(4,8);
  Double_t dpt1zBz = C_m(5,8);

  Double_t ptau1pB = -8*EB*pB_m.x()*(ptau1_m.x()/Etau1_m)*dpt1xBx -8*EB*pB_m.y()*(ptau1_m.x()/Etau1_m)*dpt1xBy -8*EB*pB_m.z()*(ptau1_m.x()/Etau1_m)*dpt1xBz -8*EB*pB_m.x()*(ptau1_m.y()/Etau1_m)*dpt1yBx -8*EB*pB_m.y()*(ptau1_m.y()/Etau1_m)*dpt1yBy -8*EB*pB_m.z()*(ptau1_m.y()/Etau1_m)*dpt1yBz -8*EB*pB_m.x()*(ptau1_m.z()/Etau1_m)*dpt1zBx -8*EB*pB_m.y()*(ptau1_m.z()/Etau1_m)*dpt1zBy -8*EB*pB_m.z()*(ptau1_m.z()/Etau1_m)*dpt1zBz;

  // ptau2,pB
  Double_t dpt2xBx = C_m(3,12);
  Double_t dpt2xBy = C_m(4,12);
  Double_t dpt2xBz = C_m(5,12);
  Double_t dpt2yBx = C_m(3,13);
  Double_t dpt2yBy = C_m(4,13);
  Double_t dpt2yBz = C_m(5,13);
  Double_t dpt2zBx = C_m(3,14);
  Double_t dpt2zBy = C_m(4,14);
  Double_t dpt2zBz = C_m(5,14);

  Double_t ptau2pB = -8*EB*pB_m.x()*(ptau2_m.x()/Etau2_m)*dpt2xBx -8*EB*pB_m.y()*(ptau2_m.x()/Etau2_m)*dpt2xBy -8*EB*pB_m.z()*(ptau2_m.x()/Etau2_m)*dpt2xBz -8*EB*pB_m.x()*(ptau2_m.y()/Etau2_m)*dpt2yBx -8*EB*pB_m.y()*(ptau2_m.y()/Etau2_m)*dpt2yBy -8*EB*pB_m.z()*(ptau2_m.y()/Etau2_m)*dpt2yBz -8*EB*pB_m.x()*(ptau2_m.z()/Etau2_m)*dpt2zBx -8*EB*pB_m.y()*(ptau2_m.z()/Etau2_m)*dpt2zBy -8*EB*pB_m.z()*(ptau2_m.z()/Etau2_m)*dpt2zBz;

  // ptau1,ptau1
  Double_t dpt1x = C_m(6,6);
  Double_t dpt1y = C_m(7,7);
  Double_t dpt1z = C_m(8,8);
  Double_t dpt1xy = C_m(6,7);
  Double_t dpt1xz = C_m(6,8);
  Double_t dpt1yz = C_m(7,8);

  Double_t ptau1ptau1 = pow(2*EB*(ptau1_m.x()/Etau1_m),2)*dpt1x + pow(2*EB*(ptau1_m.y()/Etau1_m),2)*dpt1y + pow(2*EB*(ptau1_m.z()/Etau1_m),2)*dpt1z + 8*pow(EB/Etau1_m,2)*ptau1_m.x()*ptau1_m.y()*dpt1xy + 8*pow(EB/Etau1_m,2)*ptau1_m.x()*ptau1_m.z()*dpt1xz + 8*pow(EB/Etau1_m,2)*ptau1_m.y()*ptau1_m.z()*dpt1yz;

  // ptau2,ptau2
  Double_t dpt2x = C_m(12,12);
  Double_t dpt2y = C_m(13,13);
  Double_t dpt2z = C_m(14,14);
  Double_t dpt2xy = C_m(12,13);
  Double_t dpt2xz = C_m(12,14);
  Double_t dpt2yz = C_m(13,14);

  Double_t ptau2ptau2 = pow(2*EB*(ptau2_m.x()/Etau2_m),2)*dpt2x + pow(2*EB*(ptau2_m.y()/Etau2_m),2)*dpt2y + pow(2*EB*(ptau2_m.z()/Etau2_m),2)*dpt2z + 8*pow(EB/Etau2_m,2)*ptau2_m.x()*ptau2_m.y()*dpt2xy + 8*pow(EB/Etau2_m,2)*ptau2_m.x()*ptau2_m.z()*dpt2xz + 8*pow(EB/Etau2_m,2)*ptau2_m.y()*ptau2_m.z()*dpt2yz;

  // pK,pK
  Double_t dpKx = V(19,19);
  Double_t dpKy = V(20,20);
  Double_t dpKz = V(21,21);
  Double_t dpKxy = V(19,20);
  Double_t dpKxz = V(19,21);
  Double_t dpKyz = V(20,21);

  Double_t pKpK = pow(2*EB*(pK.x()/EK),2)*dpKx + pow(2*EB*(pK.y()/EK),2)*dpKy + pow(2*EB*(pK.z()/EK),2)*dpKz + 8*pow(EB/EK,2)*pK.x()*pK.y()*dpKxy + 8*pow(EB/EK,2)*pK.x()*pK.z()*dpKxz + 8*pow(EB/EK,2)*pK.y()*pK.z()*dpKyz;

  // ptau1,prau2
  Double_t dpt1xt2x = C_m(6,12);
  Double_t dpt1xt2y = C_m(7,12);
  Double_t dpt1xt2z = C_m(8,12);
  Double_t dpt1yt2x = C_m(6,13);
  Double_t dpt1yt2y = C_m(7,13);
  Double_t dpt1yt2z = C_m(8,13);
  Double_t dpt1zt2x = C_m(6,14);
  Double_t dpt1zt2y = C_m(7,14);
  Double_t dpt1zt2z = C_m(8,14);

  Double_t ptau1ptau2 = 8*pow(EB,2)*(ptau1_m.x()/Etau1_m)*(ptau2_m.x()*Etau2_m)*dpt1xt2x + 8*pow(EB,2)*(ptau1_m.x()/Etau1_m)*(ptau2_m.y()*Etau2_m)*dpt1xt2y + 8*pow(EB,2)*(ptau1_m.x()/Etau1_m)*(ptau2_m.z()*Etau2_m)*dpt1xt2z + 8*pow(EB,2)*(ptau1_m.y()/Etau1_m)*(ptau2_m.x()*Etau2_m)*dpt1yt2x + 8*pow(EB,2)*(ptau1_m.y()/Etau1_m)*(ptau2_m.y()*Etau2_m)*dpt1yt2y + 8*pow(EB,2)*(ptau1_m.y()/Etau1_m)*(ptau2_m.z()*Etau2_m)*dpt1yt2z + 8*pow(EB,2)*(ptau1_m.z()/Etau1_m)*(ptau2_m.x()*Etau2_m)*dpt1zt2x + 8*pow(EB,2)*(ptau1_m.z()/Etau1_m)*(ptau2_m.y()*Etau2_m)*dpt1zt2y + 8*pow(EB,2)*(ptau1_m.z()/Etau1_m)*(ptau2_m.z()*Etau2_m)*dpt1zt2z ;

  Double_t MB_err_squared = sqrt( EB_pBpB + ptau1pB + ptau2pB + ptau1ptau1 + ptau2ptau2 + pKpK + ptau1ptau2 );

  // MB_err = sqrt( pow(EB/MB,2)*pow(dEB,2) + pow(pB_m.x()/MB,2)*dpBx + pow(pB_m.y()/MB,2)*dpBy + pow(pB_m.z()/MB,2)*dpBz + 2*( (pB_m.x()*pB_m.y())/pow(MB,2) )*dpBxy + 2*( (pB_m.x()*pB_m.z())/pow(MB,2) )*dpBxz + 2*( (pB_m.y()*pB_m.z())/pow(MB,2) )*dpByz + 2*(EB/MB)*(ptau1_m.x()/Etau1_m)*(-pB_m.x()/MB)*dpt1xBx + 2*(EB/MB)*(ptau1_m.x()/Etau1_m)*(-pB_m.y()/MB)*dpt1xBy + 2*(EB/MB)*(ptau1_m.x()/Etau1_m)*(-pB_m.z()/MB)*dpt1xBz + 2*(EB/MB)*(ptau1_m.y()/Etau1_m)*(-pB_m.x()/MB)*dpt1yBx + 2*(EB/MB)*(ptau1_m.y()/Etau1_m)*(-pB_m.y()/MB)*dpt1yBy + 2*(EB/MB)*(ptau1_m.y()/Etau1_m)*(-pB_m.z()/MB)*dpt1yBz + 2*(EB/MB)*(ptau1_m.z()/Etau1_m)*(-pB_m.x()/MB)*dpt1zBx + 2*(EB/MB)*(ptau1_m.z()/Etau1_m)*(-pB_m.y()/MB)*dpt1zBy + 2*(EB/MB)*(ptau1_m.z()/Etau1_m)*(-pB_m.z()/MB)*dpt1zBz + 2*(EB/MB)*(ptau2_m.x()/Etau2_m)*(-pB_m.x()/MB)*dpt2xBx + 2*(EB/MB)*(ptau2_m.x()/Etau2_m)*(-pB_m.y()/MB)*dpt2xBy + 2*(EB/MB)*(ptau2_m.x()/Etau2_m)*(-pB_m.z()/MB)*dpt2xBz + 2*(EB/MB)*(ptau2_m.y()/Etau2_m)*(-pB_m.x()/MB)*dpt2yBx + 2*(EB/MB)*(ptau2_m.y()/Etau2_m)*(-pB_m.y()/MB)*dpt2yBy + 2*(EB/MB)*(ptau2_m.y()/Etau2_m)*(-pB_m.z()/MB)*dpt2yBz + 2*(EB/MB)*(ptau2_m.z()/Etau2_m)*(-pB_m.x()/MB)*dpt2zBx + 2*(EB/MB)*(ptau2_m.z()/Etau2_m)*(-pB_m.y()/MB)*dpt2zBy + 2*(EB/MB)*(ptau2_m.z()/Etau2_m)*(-pB_m.z()/MB)*dpt2zBz + 2*(EB/MB)*(ptau1_m.x()/Etau1_m)*(EB/MB)*(ptau2_m.x()/Etau2_m)*dpt1xt2x + 2*(EB/MB)*(ptau1_m.x()/Etau1_m)*(EB/MB)*(ptau2_m.y()/Etau2_m)*dpt1xt2y + 2*(EB/MB)*(ptau1_m.x()/Etau1_m)*(EB/MB)*(ptau2_m.z()/Etau2_m)*dpt1xt2z + 2*(EB/MB)*(ptau1_m.y()/Etau1_m)*(EB/MB)*(ptau2_m.x()/Etau2_m)*dpt1yt2x + 2*(EB/MB)*(ptau1_m.y()/Etau1_m)*(EB/MB)*(ptau2_m.y()/Etau2_m)*dpt1yt2y + 2*(EB/MB)*(ptau1_m.y()/Etau1_m)*(EB/MB)*(ptau2_m.z()/Etau2_m)*dpt1yt2z + 2*(EB/MB)*(ptau1_m.z()/Etau1_m)*(EB/MB)*(ptau2_m.x()/Etau2_m)*dpt1zt2x + 2*(EB/MB)*(ptau1_m.z()/Etau1_m)*(EB/MB)*(ptau2_m.y()/Etau2_m)*dpt1zt2y + 2*(EB/MB)*(ptau1_m.z()/Etau1_m)*(EB/MB)*(ptau2_m.z()/Etau2_m)*dpt1zt2z);

  // cout << "initial chi2 = " << chisquare(x0_vars) << endl;
  // cout << "final chi2 = " << min->MinValue() << endl;
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

  TVectorD r = m - h(x);

  Double_t chi2 = r*(W*r);

  return chi2;
}

TVectorD h( TVectorD x )
{
  // The model; writes the known parameters in terms of the unkown parameters using the constraints
  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  // Double_t MB_squared = x(6);
  ROOT::Math::XYZVector ptau1( x(6), x(7), x(8) );
  ROOT::Math::XYZVector pnu1( x(9), x(10), x(11) );
  ROOT::Math::XYZVector ptau2( x(12), x(13), x(14) );
  ROOT::Math::XYZVector pnu2( x(15), x(16), x(17) );
  Double_t L1 = x(18);
  Double_t L2 = x(19);
  Double_t L = x(20);
  Double_t LK = x(21);

  // Tau and neutrino mass constraints:
  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );
  // Double_t EB = sqrt( MB_squared + pB.Mag2() );

  ROOT::Math::PxPyPzEVector Ptau1( ptau1.x(), ptau1.y(), ptau1.z(), Etau1 );
  ROOT::Math::PxPyPzEVector Ptau2( ptau2.x(), ptau2.y(), ptau2.z(), Etau2 );
  ROOT::Math::PxPyPzEVector Pnu1( pnu1.x(), pnu1.y(), pnu1.z(), Enu1 );
  ROOT::Math::PxPyPzEVector Pnu2( pnu2.x(), pnu2.y(), pnu2.z(), Enu2 );
  // ROOT::Math::PxPyPzEVector PB( pB.x(), pB.y(), pB.z(), EB);

  // pB must point back to the PV:
  ROOT::Math::XYZPoint PV( BV.x() - L*pB.x(), BV.y() - L*pB.y(), BV.z() - L*pB.z() );
  // ptau1 must point back to DV1:
  ROOT::Math::XYZPoint DV1( BV.x() + L1*ptau1.x(), BV.y() + L1*ptau1.y(), BV.z() + L1*ptau1.z() );
  // 4-momentum conservation in DV1:
  ROOT::Math::XYZVector p3pi1 = ptau1 - pnu1;
  Double_t E3pi1 = Etau1 - Enu1;
  // Double_t m3pi1_squared = pow(mtau,2) -2*Ptau1.Dot(Pnu1); 
  // ptau2 must point back to DV2:
  ROOT::Math::XYZPoint DV2( BV.x() + L2*ptau2.x(), BV.y() + L2*ptau2.y(), BV.z() + L2*ptau2.z() );
  // 4-momentum conservation in DV2:
  ROOT::Math::XYZVector p3pi2 = ptau2 - pnu2;
  // Double_t m3pi2_squared = pow(mtau,2) -2*Ptau2.Dot(Pnu2);  
  Double_t E3pi2 = Etau2 - Enu2;
  // B+ must lie in the K+ trajectory:
  ROOT::Math::XYZPoint RP( BV.x() + LK*( pB.x() - ptau1.x() - ptau2.x() ), BV.y() + LK*( pB.y() - ptau1.y() - ptau2.y() ), RPz );
  // 4-momentum conservation in BV:
  // ROOT::Math::XYZVector p6piK = pB - pnu1 - pnu2;
  ROOT::Math::XYZVector pK = pB - ptau1 - ptau2;
  // Double_t m6piK_squared = MB_squared - 2*PB.Dot(Pnu1) - 2*PB.Dot(Pnu2) + 2*Pnu1.Dot(Pnu2);
  // Double_t E6piK = EB - Enu1 - Enu2;

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
  // h(9) = m3pi1_squared;
  h(9) = E3pi1;
  h(10) = DV2.x();
  h(11) = DV2.y();
  h(12) = DV2.z();
  h(13) = p3pi2.x();
  h(14) = p3pi2.y();
  h(15) = p3pi2.z();
  // h(16) = m3pi2_squared;
  h(16) = E3pi2;
  h(17) = RP.x();
  h(18) = RP.y();
  // h(19) = RP.z();
  h(19) = pK.x();
  h(20) = pK.y();
  h(21) = pK.z();
  // h(22) = m6piK_squared;
  // h(22) = E6piK;

  return h;
}

TVectorD x_initial_estimate2( TVectorD m, ROOT::Math::XYZPoint BV ) // Marseille initialisation for x
{
  // Builds an initial estimate for x based on the Mariseille analytical calculations; it uses the offline estimate for BV  
  // (this initialisation does not apply the constraint: pB must point back to the PV) 

  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  // Double_t m3pi1_squared = m(9); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  // Double_t m3pi2_squared = m(16); 
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  // ROOT::Math::XYZVector p6piK( m(19), m(20), m(21) );
  // Double_t m6piK_squared = m(22);
  // Double_t E6piK = m(22);

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  // Double_t E3pi1 = sqrt( m3pi1_squared + p3pi1.Mag2() );
  // Double_t E3pi2 = sqrt( m3pi2_squared + p3pi2.Mag2() );
  // Double_t E6piK = sqrt( m6piK_squared + p6piK.Mag2() );

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

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

  // ROOT::Math::XYZVector pB = p6piK + pnu1 + pnu2;
  // Double_t EB = E6piK + Enu1 + Enu2;

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t EK = sqrt( pow(mK,2) + pK.Mag2() );

  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
  Double_t EB = Etau1 + Etau2 + EK;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector p6piK( m(19), m(20), m(21) );
  Double_t E6piK = m(22);

  Double_t EK = E6piK - E3pi1 - E3pi2;
  ROOT::Math::XYZVector pK = p6piK - p3pi1 - p3pi2;

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  ROOT::Math::XYZPoint PV_t = makeTransformation_point( pK, RP, PV, false );
  ROOT::Math::XYZPoint DV1_t = makeTransformation_point( pK, RP, DV1, false );
  ROOT::Math::XYZVector p3pi1_t = makeTransformation_vec( pK, RP, p3pi1, false );
  ROOT::Math::XYZPoint DV2_t = makeTransformation_point( pK, RP, DV2, false );
  ROOT::Math::XYZVector p3pi2_t = makeTransformation_vec( pK, RP, p3pi2, false );
  ROOT::Math::XYZVector pK_t = makeTransformation_vec( pK, RP, pK, false );

  Double_t a1 = (DV1_t.y())/(DV1_t.x());
  Double_t a2 = (DV2_t.y())/(DV2_t.x());
  Double_t b = (PV_t.y() -a1*PV_t.x())/(a2*PV_t.x() - PV_t.y());
  Double_t c = b*(DV1_t.x())/(DV2_t.x());
  Double_t d = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()));
  Double_t e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
  Double_t f = ( PV_t.x()*sqrt(pK_t.Mag2()) )/( (1+b)*DV1_t.x() -(1+c)*PV_t.x() );
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

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );
  Double_t EB = Etau1 + Etau2 + EK;

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
  x0(6) = EB;
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