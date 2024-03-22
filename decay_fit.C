using namespace std;

////////////////////////////////////////           default TMINUIT           ////////////////////////////////////////////////////////

// Global variables
int dimM = 22; // number of known parameters
int dimX = 22; // number of unkown parameters

TVectorD m( dimM ); // vector of measured parameters
TMatrixDSym V( dimM ); // covariance matrix of m

TMatrixDSym W( dimM ); // weights matrix, W=V^-1
TVectorD x0( dimX ); // initial estimate for the vector of unkown parameters

Double_t RPz = 0; // z-component of the reference point in the K+ trajectory (fixed)
// Double_t RPz_t = 0; // truth-match

// Things saved from the minimisation:
Double_t chi2, MB, MB_err, taup_PE, taum_PE, taup_M, taum_M, deltaT; 
TVectorD x_m(dimX); // vector of unkown parameters after the minimisation
TVectorD x_err_m(dimX); // vector of x errors from minimisation
TMatrixDSym C_m(dimX); // covariance matrix of the fit parameters
Int_t status, init, nIter; // status of the fit
ULong64_t eventNumber;

Double_t mtau = 1776.86; // MeV (PDG 2023)
Double_t mpion = 139.57039; // MeV
Double_t mkaon = 493.677; // MeV

// Functions
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);

TVectorD x_initial_estimate_0( TVectorD m );
TVectorD x_initial_estimate_1( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_2( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_3( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_4( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_5( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_6( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_7( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_8( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_9( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_MLP( TVectorD m, ROOT::Math::XYZPoint BV, Int_t species );

void lowest_chi2(ROOT::Math::XYZPoint BV, Int_t species, TVectorD m);

TVectorD model( TVectorD x );
TVectorD linear_model( TVectorD x );
TMatrixD dh_dx( TVectorD x );
Double_t chisquare( const Double_t* x_values );
void minimize( ROOT::Math::XYZPoint BV, int init, int species, TVectorD m );
TVectorD transform_m( TVectorD m );
TMatrixDSym transform_V( TVectorD m, TMatrixDSym V );
TMatrixDSym scale_tau_DVs_uncertainties(TMatrixDSym V, Double_t factor);
TMatrixDSym scale_uncertainties(TMatrixDSym V, Double_t factor);
TMatrixDSym invert_sum(TMatrixDSym V1, TMatrixDSym V2);
Double_t B_mass();
Double_t B_mass_error();
TMatrixDSym bootstrap_errors(Int_t n_rep, Bool_t init, ROOT::Math::XYZPoint BV, Int_t species, Double_t m_vars[dimM], Double_t V_vars[dimM*dimM]);
TMatrixD H_matrix( TVectorD x );

Bool_t useLinear = false;

void decay_fit(int year, TString RECO_files, int species, int line)
{
  TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files, 1, line);
  TChain* t = new TChain("DecayTree");
  t->AddFileInfoList((TCollection*)fc->GetList());

  UInt_t num_entries = t->GetEntries();
  Double_t m_vars[dimM];
  Double_t V_vars[dimM][dimM];
  Double_t BVx, BVy, BVz;

  // t->SetBranchAddress("df_m", m_vars);
  // t->SetBranchAddress("df_V", V_vars);

  for(int i = 0; i < dimM; i++)
  {
    t->SetBranchAddress(Form("df_m_%i",i+1), &m_vars[i]);
    for(int j = 0; j < dimM; j++)
    {
      t->SetBranchAddress(Form("df_V_%i_%i",i+1,j+1), &V_vars[i][j]);
    }
  }

  t->SetBranchAddress("Kp_RP_Z", &RPz);
  t->SetBranchAddress("Bp_ENDVERTEX_X", &BVx);
  t->SetBranchAddress("Bp_ENDVERTEX_Y", &BVy);
  t->SetBranchAddress("Bp_ENDVERTEX_Z", &BVz);
  t->SetBranchAddress("eventNumber", &eventNumber);

  TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/%i.root",year,species,line), "RECREATE");
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
    "df_LK",
  };

  TString name_x_err[] = {
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
    tout->Branch(name_x[i], &x_m(i));
    tout->Branch(name_x_err[i], &x_err_m(i));
    for(int j = 0; j < dimX; j++)
    {
      tout->Branch(Form("df_C_%i_%i",i,j), &C_m(i,j));
    }
  }
  tout->Branch("df_chi2", &chi2);
  tout->Branch("df_init", &init);
  tout->Branch("df_Bp_M", &MB);
  tout->Branch("df_Bp_MERR", &MB_err);
  tout->Branch("df_status", &status);
  tout->Branch("df_time", &deltaT);
  tout->Branch("df_nIter", &nIter);

  TStopwatch watch_total;
  
  // Loop over events
  for(int evt = 0; evt < 100; evt++)
  {
    TStopwatch watch;
    t->GetEntry(evt);
    
    // 1)  Measured parameters, m, and their covariance matrix V in the track parametrisation
    // m = (PV,DV1,p1,p2,p3,DV2,p4,p5,p6,RP_T,pK) -> 32 parameters
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

    // Eigenvalues of V
    // TVectorD V_eigen_values(dimM);  
    // TMatrixD V_eigen_vectors = V.EigenVectors(V_eigen_values);
    // V_eigen_values.Print();

    // Correlation matrix
    // TMatrixDSym V_corr(dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //   for(int j = 0; j < dimM; j++)
    //   {
    //     V_corr(i,j) = V(i,j)/sqrt( V(i,i)*V(j,j) );
    //   }
    // }
    // V_corr.Print();

    // scale uncertainties on tau DVs
    // V = scale_uncertainties(V,-0.999999999);

    // 3) Weights matrix W = V'^-1
    W = V;
    W.Invert();
    // W.Print();

    // cout << "Determinant of V is " << V.Determinant() << endl;
    // cout << "Determinant of W is " << W.Determinant() << endl;

    // Eigenvalues of W
    // TVectorD W_eigen_values(dimM);  
    // TMatrixD W_eigen_vectors = W.EigenVectors(W_eigen_values);
    // W_eigen_values.Print();

    // WV ~ I ?
    // TMatrixD identity = W*V;
    // identity.Print();
    
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz ); // offline estimate of BV (necessary for some initialisations)
    // ROOT::Math::XYZPoint BV( RPx_t, RPy_t, RPz_t ); // truth-match

    // Estimate errors on the fit parameters with bootstrap method
    // Double_t* m_vars_double[] = {};
    // Double_t* V_vars_double[] = {};
    // for(int i = 0; i < dimM; i++)
    // {
    //   m_vars_double[i] = &m_vars[i];
    //   for(int j = 0; j < dimM; j++)
    //   { 
    //     V_vars_double[i][j] = V_vars[i][j];
    //   }
    // }
    // Double_t V_vars_double[dimM*dimM];
    // for(int i = 0; i < dimM*dimM; i++)
    // {
    //   V_vars_double[i] = V_vars[i/dimM][i%dimM];   
    // }
    // C_m = bootstrap_errors(100,-2,BV,species, m_vars, V_vars_double);
    // C_m.Print();

    // TVectorD C_eigen_values(dimX);  
    // TMatrixD C_eigen_vectors = C_m.EigenVectors(C_eigen_values);
    // C_eigen_values.Print();

    // MLP initialisation
    minimize(BV, -2, species, m);

    // Lowest chi2 initialisation
    // lowest_chi2(BV, species, m);

    // DTF sequence 012 
    // minimize(BV, 1, species, m);
    // if(status > 1)
    // {
    //   minimize(BV, 2, species, m);
    // }
    // if(status > 1)
    // {
    //   minimize(BV, 3, species, m);
    // }

    // x_m.Print();
    // x_err_m.Print();

    // MB_err = B_mass_error();

    cout << "init = " << init << endl;
    cout << "final chi2 = " << chi2 << endl;
    cout << "MB = " << MB << " +/- " << MB_err << endl;
    cout << "Status = " << status << endl;

    watch.Stop();
    deltaT = watch.RealTime();
    tout->Fill();
  }
  watch_total.Stop();

  cout << "FINISHED in " << watch_total.RealTime() << " seconds." << endl;

  fout->cd();
  tout->Write();
  fout->Close();
}

void minimize( ROOT::Math::XYZPoint BV, int init, int species, TVectorD m )
{
  // 1) This function initialises the vector of unkown parameters x with the result of analytical calculations; by specifying the value of init we choose which analytical calculation we want to use: 
  //   init = 0 : original calculations based on Anne Keune thesis; does not use the BV offline estimate 
  //   init = 2 : uses the Marseille initialisation; it uses the BV offline estimate and assumes that the neutrino takes as much momentum as possible from the tau 
  // 2) It sets up the minimizer: defining the initial values and bounds on the x parameters  
  // 3) It builds the chi^2 function with x0 by calling the function chisquare 
  // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

  // Initial values for the unkown parameters x0
  if(init == 0) // original calculations = Anne Keune
  {
    x0 = x_initial_estimate_0( m ); 
  }
  else if(init == 1) // K*tautau calculations; tau momentum direction initialised based on vertices 
  {
    x0 = x_initial_estimate_1( m, BV );
  }
  else if(init == 2) // K*tautau calculations; tau momentum direction initialised based on 3pi visible momentum 
  {
    x0 = x_initial_estimate_2( m, BV );
  }
  else if(init == 3) // RD calculations 
  {
    x0 = x_initial_estimate_3( m, BV ); 
  }
  else if(init == 4) // K*tautau; tau+ from vertices, tau- from pions
  {
    x0 = x_initial_estimate_4( m, BV ); 
  }
  else if(init == 5) // K*tautau; tau- from vertices, tau+ from pions
  {
    x0 = x_initial_estimate_5( m, BV ); 
  }
  else if(init == 6) // Mixes Marseille (tau+) and K*tautau vertices (tau-)
  {
    x0 = x_initial_estimate_6( m, BV ); 
  }
  else if(init == 7)  // Mixes Marseille (tau-) and K*tautau vertices (tau+)
  {
    x0 = x_initial_estimate_7( m, BV ); 
  }
  else if(init == 8) // Mixes Marseille (tau+) and K*tautau pions (tau-)
  {
    x0 = x_initial_estimate_8( m, BV ); 
  }
  else if(init == 9) // Mixes Marseille (tau-) and K*tautau pions (tau+)
  {
    x0 = x_initial_estimate_9( m, BV ); 
  }
  else if(init == -2)
  {
    x0 = x_initial_estimate_MLP( m, BV, species );
  }

  Double_t x0_vars[dimX], x0_err[dimX];
  for(int i = 0; i < dimX; i++)
  {
    x0_vars[i] = x0(i); // x0_vars is a Double_t, x0 is a TVectorD; the function to minimise must receive a Double_t as input
    x0_err[i] = 0.1*abs(x0_vars[i]); // initial errors on x0; they are used as the first step size in the minimisation; considering a 10% error for now
  }

  // Initial step size
  // x0_err[0] = 0.05; // BVx
  // x0_err[1] = 0.05; // BVy
  // x0_err[2] = 0.5; // BVz
  // x0_err[4] = 1000; // pbx
  // x0_err[5] = 1000; // pby
  // x0_err[6] = 10000; // pbz
  // x0_err[7] = 100; // ptau1x
  // x0_err[8] = 100; // ptau1y
  // x0_err[9] = 1000; // ptau1z
  // x0_err[10] = 100; // pnu1x
  // x0_err[11] = 100; // pnu1y
  // x0_err[12] = 1000; // pnu1z
  // x0_err[13] = 100; // ptau2x
  // x0_err[14] = 100; // ptau2y
  // x0_err[15] = 1000; // ptau2z
  // x0_err[16] = 100; // pnu2x
  // x0_err[17] = 100; // pnu2y
  // x0_err[18] = 1000; // pnu2z
  // x0_err[19] = 2.5*pow(10,-7); // L1
  // x0_err[20] = 2.5*pow(10,-7); // L2
  // x0_err[21] = 2.5*pow(10,-5); // L
  // x0_err[22] = 1*pow(10,-4); // LK

  // x0.Print();
  // m.Print();
  // hx.Print();
  // TVectorD r = m-hx;
  // r.Print();

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "MigradImproved");    
  ROOT::Math::Functor func(&chisquare, dimX);

  min->SetMaxIterations(10000000);  
  min->SetMaxFunctionCalls(10000000);
  min->SetPrintLevel(0);
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

  for(int i = 0; i < dimX; i++)
  {
    min->SetVariable(i, name[i], x0[i], x0_err[i]);  
  }

  min->Minimize();
  // min->Hesse();

  // Return / save results from the fit
  chi2 = min->MinValue();
  status = min->Status();
  nIter = min->NIterations();

  const Double_t *xMin = min->X();
  for(int i = 0; i < dimX; i++)
  {
    x_m(i) = xMin[i];
  }
  MB = B_mass();

}

TMatrixD H_matrix( TVectorD x )
{
  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  ROOT::Math::XYZVector ptau1( x(6), x(7), x(8) );
  ROOT::Math::XYZVector pnu1( x(9), x(10), x(11) );
  ROOT::Math::XYZVector ptau2( x(12), x(13), x(14) );
  ROOT::Math::XYZVector pnu2( x(15), x(16), x(17) );
  Double_t L1 = x(18);
  Double_t L2 = x(19);
  Double_t L = x(20);
  Double_t LK = x(21);

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );

  TMatrixD H(dimM,dimX);
  // BVx
  H(0,0) = 1;
  H(0,3) = 1;
  H(0,10) = 1;
  H(0,17) = 1;

  // BVy
  H(1,1) = 1;
  H(1,4) = 1;
  H(1,11) = 1;
  H(1,18) = 1;

  // BVz
  H(2,2) = 1;
  H(2,5) = 1;
  H(2,12) = 1;

  // pBx
  H(3,0) = -L;
  H(3,19) = 1;

  // pBy
  H(4,1) = -L;
  H(4,20) = 1;

  // pBz
  H(5,2) = -L;
  H(5,21) = 1;

  // ptau1x
  H(6,3) = L1;
  H(6,6) = 1;
  H(6,9) = ptau1.x()/Etau1;
  H(6,19) = -1;

  // ptau1y
  H(7,4) = L1;
  H(7,7) = 1;
  H(7,9) = ptau1.y()/Etau1;
  H(7,20) = -1;

  // ptau1z
  H(8,5) = L1;
  H(8,8) = 1;
  H(8,9) = ptau1.z()/Etau1;
  H(8,21) = -1;

  // pnu1x
  H(9,6) = -1;
  H(9,9) = -pnu1.x()/Enu1;

  // pnu1y
  H(10,7) = -1;
  H(10,9) = -pnu1.y()/Enu1;

  // pnu1z
  H(11,8) = -1;
  H(11,9) = -pnu1.z()/Enu1;

  // ptau2x
  H(12,10) = L2;
  H(12,13) = 1;
  H(12,16) = ptau2.x()/Etau2;
  H(12,19) = -1;

  // ptau2y
  H(13,11) = L2;
  H(13,14) = 1;
  H(13,16) = ptau2.y()/Etau2;
  H(13,20) = -1;

  // ptau2z
  H(14,12) = L2;
  H(14,15) = 1;
  H(14,16) = ptau2.z()/Etau2;
  H(14,21) = -1;

  // pnu2x
  H(15,13) = -1;
  H(15,16) = -pnu2.x()/Enu2;

  // pnu2y
  H(16,14) = -1;
  H(16,16) = -pnu2.y()/Enu2;

  // pnu2z
  H(17,15) = -1;
  H(17,16) = -pnu2.z()/Enu2;

  // L1
  H(18,3) = ptau1.x();
  H(18,4) = ptau1.y();
  H(18,5) = ptau1.z();

  // L2
  H(19,10) = ptau2.x();
  H(19,11) = ptau2.y();
  H(19,12) = ptau2.z();

  // L
  H(20,0) = -pB.x();
  H(20,1) = -pB.y();
  H(20,2) = -pB.z();

  // LK
  H(21,17) = -pK.x();
  H(21,18) = -pK.y();

  return H;
}


TMatrixDSym bootstrap_errors(Int_t n_rep, Bool_t init, ROOT::Math::XYZPoint BV, Int_t species, Double_t m_vars[dimM], Double_t V_vars[dimM*dimM])
{
  // TFile* fout = new TFile("/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_10/x_gaussians.root", "RECREATE");
  TTree* t = new TTree("t", "t");

  TString x_name[] = {
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

  Double_t x_vars[dimX];
  for(int i = 0; i < dimX; i++)
  {
    t->Branch(x_name[i], &x_vars[i]);
  }

  TString m_name[] = {
    "PVx",
    "PVy",
    "PVz",
    "DV1x",
    "DV1y",
    "DV1z",
    "p3pi1x",
    "p3pi1y",
    "p3pi1z",
    "E3pi1",
    "DV2x",
    "DV2y",
    "DV2z",
    "p3pi2x",
    "p3pi2y",
    "p3pi2z",
    "E3pi2",
    "RPx",
    "RPy",
    "pKx",
    "pKy",
    "pKz"
  };

  TVectorD m_random(dimM);
  for(int i = 0; i < dimM; i++)
  {
    t->Branch(m_name[i], &m_random[i]);
  }

  ROOT::Math::GSLRandomEngine rnd;
  rnd.Initialize();

  // Double_t m_vars[dimM];
  // Double_t V_vars[dimM*dimM];
  Double_t m_random_vars[dimM];

  for(int n = 0; n < n_rep; n++)
  {
    // Sample m from a random multivariate Gaussian
    rnd.GaussianND(dimM, m_vars, V_vars, m_random_vars);

    for(int i = 0; i < dimM; i++)
    {
      m_random(i) = m_random_vars[i];
    }
    // m.Print();
    // m_random.Print();

    minimize(BV, -2, species, m_random);

    for(int k = 0; k < dimX; k++)
    {
      x_vars[k] = x_m(k);
    }

    t->Fill();
  }

  TMatrixDSym V_errors(dimX);
  for(int i = 0; i < dimX; i++)
  {
    t->Draw(x_name[i]+Form(" >> h_%i",i));
    TH1D* h = (TH1D*)gPad->GetPrimitive(Form("h_%i",i));
    x_err_m[i] = h->GetStdDev();
    for(int j = 0; j < dimX; j++)
    {
      t->Draw(x_name[i]+" : "+x_name[j]+Form(" >> h_%i_%i",i,j));
      TH2D* h = (TH2D*)gDirectory->Get(Form("h_%i_%i",i,j));
      Double_t C_ij = h->GetCorrelationFactor();
      Double_t V_ii = pow(x_err_m[i],2);
      Double_t V_jj = pow(x_err_m[j],2);
      V_errors(i,j) = C_ij*sqrt( V_ii*V_jj );
    }
  }

  // It's because it hasnt computed all x_err_i yet
  for(int i = 0; i < dimX; i++)
  {
    for(int j = 0; j < i; j++)
    {
      V_errors(j,i) = V_errors(i,j);
    }
  }

  // TVectorD m_err(dimM);
  // for(int i = 0; i < dimM; i++)
  // {
  //   m_err(i) = sqrt(V(i,i));
  // }
  // m_err.Print();

  // fout->cd();
  // t->Write();
  // fout->Close();

  return V_errors;
}

Double_t B_mass()
{
  ROOT::Math::XYZVector pB_m( x_m(3), x_m(4), x_m(5) );
  ROOT::Math::XYZVector ptau1_m( x_m(6), x_m(7), x_m(8) );
  ROOT::Math::XYZVector ptau2_m( x_m(12), x_m(13), x_m(14) );
  Double_t Etau1_m = sqrt( pow(mtau,2) + ptau1_m.Mag2() );
  Double_t Etau2_m = sqrt( pow(mtau,2) + ptau2_m.Mag2() );
  ROOT::Math::XYZVector pK_m = pB_m - ptau1_m - ptau2_m;
  // ROOT::Math::XYZVector pK_m( m(19), m(20), m(21) );
  Double_t EK_m = sqrt( pow(mkaon,2) + pK_m.Mag2() );
  Double_t EB_m = Etau1_m + Etau2_m + EK_m; // Energy convervation in BV
  Double_t MB_squared_m = pow(EB_m,2) - pB_m.Mag2();

  if( MB_squared_m > 0 )
  {
    MB = sqrt( MB_squared_m );
  }
  else
  {
    MB = -sqrt( abs(MB_squared_m) );
  }

  return MB;
}

Double_t B_mass_error()
{
  ROOT::Math::XYZVector pB_m( x_m(3), x_m(4), x_m(5) );
  ROOT::Math::XYZVector ptau1_m( x_m(6), x_m(7), x_m(8) );
  ROOT::Math::XYZVector ptau2_m( x_m(12), x_m(13), x_m(14) );
  Double_t Etau1_m = sqrt( pow(mtau,2) + ptau1_m.Mag2() );
  Double_t Etau2_m = sqrt( pow(mtau,2) + ptau2_m.Mag2() );
  ROOT::Math::XYZVector pK_m = pB_m - ptau1_m - ptau2_m;
  // ROOT::Math::XYZVector pK_m( m(19), m(20), m(21) );
  Double_t EK_m = sqrt( pow(mkaon,2) + pK_m.Mag2() );
  Double_t EB_m = Etau1_m + Etau2_m + EK_m; // Energy convervation in BV
  
  TMatrixD Cprime(26,26);
  for(int i = 0; i < dimX; i++)
  {
    for(int j = 0; j < dimX; j++)
    {
      Cprime(i,j) = C_m(i,j);
    }
  }

  // Add pK
  TMatrixD J1(26,26);
  for(int i = 0; i < dimX; i++)
  {
    for(int j = 0; j < dimX; j++)
    {
      if(i == j)
      {
        J1(i,j) = 1.;
      }
      else
      {
        J1(i,j) = 0.;
      }
    }
  }
  J1(22,3) = 1;
  J1(23,4) = 1;
  J1(24,5) = 1;

  J1(22,6) = -1;
  J1(23,7) = -1;
  J1(24,8) = -1;

  J1(22,12) = -1;
  J1(23,13) = -1;
  J1(24,14) = -1;

  J1(25,25) = 1;

  TMatrixD J1_transpose = J1;
  J1_transpose.T();

  Cprime = J1*(Cprime*J1_transpose);
  // Cprime.Print();

  // add mB^2 part
  TMatrixD J2(26,26);

  for(int i = 0; i < 25; i++)
  {
    for(int j = 0; j < 25; j++)
    {
      if(i == j)
      {
        J2(i,j) = 1.;
      }
      else
      {
        J2(i,j) = 0.;
      }
    }
  }

  J2(25,3) = -2*pB_m.x();
  J2(25,4) = -2*pB_m.y();
  J2(25,5) = -2*pB_m.z();

  J2(25,6) = 2*EB_m*(ptau1_m.x()/Etau1_m);
  J2(25,7) = 2*EB_m*(ptau1_m.y()/Etau1_m);
  J2(25,8) = 2*EB_m*(ptau1_m.z()/Etau1_m);

  J2(25,12) = 2*EB_m*(ptau2_m.x()/Etau2_m);
  J2(25,13) = 2*EB_m*(ptau2_m.y()/Etau2_m);
  J2(25,14) = 2*EB_m*(ptau2_m.z()/Etau2_m);

  J2(25,22) = 2*EB_m*(pK_m.x()/EK_m);
  J2(25,23) = 2*EB_m*(pK_m.y()/EK_m);
  J2(25,24) = 2*EB_m*(pK_m.z()/EK_m);

  TMatrixD J2_transpose = J2;
  J2_transpose.T();

  Cprime = J2*(Cprime*J2_transpose);
  // Cprime.Print();

  Double_t dMB_squared = sqrt(Cprime(25,25));
  Double_t dMB = dMB_squared/(2*abs(MB));
  
  return dMB;
}

Double_t chisquare( const Double_t* x_values )
{
  // Computes the chi^2
  TVectorD x(dimX);
  for(int i = 0; i < dimX; i++)
  {
    x(i) = x_values[i];
  }

  TVectorD hx = model(x);

  TVectorD r = m - hx;
  
  Double_t chi2 = r*(W*r);

  return chi2;
}

TVectorD model( TVectorD x )
{
  // The model; writes the known parameters in terms of the unkown parameters using the model constraints
  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  ROOT::Math::XYZVector ptau1( x(6), x(7), x(8) );
  ROOT::Math::XYZVector pnu1( x(9), x(10), x(11) );
  ROOT::Math::XYZVector ptau2( x(12), x(13), x(14) );
  ROOT::Math::XYZVector pnu2( x(15), x(16), x(17) );
  Double_t L1 = x(18);
  Double_t L2 = x(19);
  Double_t L = x(20);
  Double_t LK = x(21);

  // Tau and neutrino mass constraints (4):
  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  ROOT::Math::PxPyPzEVector Ptau1( ptau1.x(), ptau1.y(), ptau1.z(), Etau1 );
  ROOT::Math::PxPyPzEVector Ptau2( ptau2.x(), ptau2.y(), ptau2.z(), Etau2 );
  ROOT::Math::PxPyPzEVector Pnu1( pnu1.x(), pnu1.y(), pnu1.z(), Enu1 );
  ROOT::Math::PxPyPzEVector Pnu2( pnu2.x(), pnu2.y(), pnu2.z(), Enu2 );

  // // 3-momentum conservation in BV (3):
  // ROOT::Math::XYZVector pK = pB - ptau1 - ptau2;

  // // pB must point back to the PV (2):
  // ROOT::Math::XYZPoint PV( BV.x() - L*pB.x(), BV.y() - L*pB.y(), BV.z() - L*pB.z() );

  // // ptau1 must point back to the BV (2):
  // ROOT::Math::XYZPoint DV1( BV.x() + L1*ptau1.x(), BV.y() + L1*ptau1.y(), BV.z() + L1*ptau1.z() );

  // // 4-momentum conservation in DV1 (4):
  // ROOT::Math::PxPyPzEVector P3pi1 = Ptau1 - Pnu1;
  // ROOT::Math::XYZVector p3pi1( P3pi1.X(), P3pi1.Y(), P3pi1.Z() );
  // Double_t E3pi1 = P3pi1.E();

  // // ptau2 must point back to DV2 (2):
  // ROOT::Math::XYZPoint DV2( BV.x() + L2*ptau2.x(), BV.y() + L2*ptau2.y(), BV.z() + L2*ptau2.z() );

  // // 4-momentum conservation in DV2 (4):
  // ROOT::Math::PxPyPzEVector P3pi2 = Ptau2 - Pnu2;
  // ROOT::Math::XYZVector p3pi2( P3pi2.X(), P3pi2.Y(), P3pi2.Z() );
  // Double_t E3pi2 = P3pi2.E();

  // // B+ must lie in the K+ trajectory (2):
  // ROOT::Math::XYZPoint RP( BV.x() - LK*pK.x(), BV.y() - LK*pK.y(), RPz );

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
  h(9) = E3pi1;
  h(10) = DV2.x();
  h(11) = DV2.y();
  h(12) = DV2.z();
  h(13) = p3pi2.x();
  h(14) = p3pi2.y();
  h(15) = p3pi2.z(); 
  h(16) = E3pi2;
  h(17) = RP.x();
  h(18) = RP.y(); 
  h(19) = pK.x();
  h(20) = pK.y();
  h(21) = pK.z();

  return h;
}

TVectorD linear_model( TVectorD x )
{
  return model(x0) + H_matrix(x0)*(x - x0);
}

void lowest_chi2(ROOT::Math::XYZPoint BV, Int_t species, TVectorD m)
{
  // Init 0
  minimize(BV, 0, species, m);
  Int_t status0 = status;
  Double_t MB_0 = MB;
  Double_t MB_err0 = MB_err;
  TVectorD X0 = x_m;
  TVectorD X0_ERR = x_err_m;
  Double_t chi2_0 = chi2;
  Int_t nIter_0 = nIter;

  // cout << "status0 = " << status0 << endl;
  // cout << "MB0 = " << MB_0 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 1
  minimize(BV, 1, species, m);
  Int_t status1 = status;
  Double_t MB_1 = MB;
  Double_t MB_err1 = MB_err;
  TVectorD X1 = x_m;
  TVectorD X1_ERR = x_err_m;
  Double_t chi2_1 = chi2;
  Int_t nIter_1 = nIter;

  // cout << "status1 = " << status1 << endl;
  // cout << "MB1 = " << MB_1 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 2
  minimize(BV, 2, species, m);
  Int_t status2 = status;
  Double_t MB_2 = MB;
  Double_t MB_err2 = MB_err;
  TVectorD X2 = x_m;
  TVectorD X2_ERR = x_err_m;
  Double_t chi2_2 = chi2;
  Int_t nIter_2 = nIter;

  // cout << "status2 = " << status2 << endl;
  // cout << "MB2 = " << MB_2 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 3
  minimize(BV, 3, species, m);
  Int_t status3 = status;
  Double_t MB_3 = MB;
  Double_t MB_err3 = MB_err;
  TVectorD X3 = x_m;
  TVectorD X3_ERR = x_err_m;
  Double_t chi2_3 = chi2;
  Int_t nIter_3 = nIter;

  // cout << "status3 = " << status3 << endl;
  // cout << "MB3 = " << MB_3 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 4
  minimize(BV, 4, species, m);
  Int_t status4 = status;
  Double_t MB_4 = MB;
  Double_t MB_err4 = MB_err;
  TVectorD X4 = x_m;
  TVectorD X4_ERR = x_err_m;
  Double_t chi2_4 = chi2;
  Int_t nIter_4 = nIter;

  // cout << "status4 = " << status4 << endl;
  // cout << "MB4 = " << MB_4 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 5
  minimize(BV, 5, species, m);
  Int_t status5 = status;
  Double_t MB_5 = MB;
  Double_t MB_err5 = MB_err;
  TVectorD X5 = x_m;
  TVectorD X5_ERR = x_err_m;
  Double_t chi2_5 = chi2;
  Int_t nIter_5 = nIter;

  // cout << "status5 = " << status5 << endl;
  // cout << "MB5 = " << MB_5 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 6
  minimize(BV, 6, species, m);
  Int_t status6 = status;
  Double_t MB_6 = MB;
  Double_t MB_err6 = MB_err;
  TVectorD X6 = x_m;
  TVectorD X6_ERR = x_err_m;
  Double_t chi2_6 = chi2;
  Int_t nIter_6 = nIter;

  // cout << "status6 = " << status6 << endl;
  // cout << "MB6 = " << MB_6 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 7
  minimize(BV, 7, species, m);
  Int_t status7 = status;
  Double_t MB_7 = MB;
  Double_t MB_err7 = MB_err;
  TVectorD X7 = x_m;
  TVectorD X7_ERR = x_err_m;
  Double_t chi2_7 = chi2;
  Int_t nIter_7 = nIter;

  // cout << "status7 = " << status7 << endl;
  // cout << "MB7 = " << MB_7 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 8
  minimize(BV, 8, species, m);
  Int_t status8 = status;
  Double_t MB_8 = MB;
  Double_t MB_err8 = MB_err;
  TVectorD X8 = x_m;
  TVectorD X8_ERR = x_err_m;
  Double_t chi2_8 = chi2;
  Int_t nIter_8 = nIter;

  // cout << "status8 = " << status8 << endl;
  // cout << "MB8 = " << MB_8 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  // Init 9
  minimize(BV, 9, species, m);
  Int_t status9 = status;
  Double_t MB_9 = MB;
  Double_t MB_err9 = MB_err;
  TVectorD X9 = x_m;
  TVectorD X9_ERR = x_err_m;
  Double_t chi2_9 = chi2;
  Int_t nIter_9 = nIter;

  // cout << "status9 = " << status9 << endl;
  // cout << "MB9 = " << MB_9 << endl;
  // cout << " <<<<<<<<<<<<<<<<<< " << endl;

  Int_t status_vec[10] = { status0, status1, status2, status3, status4, status5, status6, status7, status8, status9 };
  Double_t MB_vec[10] = { MB_0, MB_1, MB_2, MB_3, MB_4, MB_5, MB_6, MB_7, MB_8, MB_9 }; 
  Double_t MB_err_vec[10] = { MB_err0, MB_err1, MB_err2, MB_err3, MB_err4, MB_err5, MB_err6, MB_err7, MB_err8, MB_err9 };
  Double_t chi2_vec[10] = { chi2_0, chi2_1, chi2_2, chi2_3, chi2_4, chi2_5, chi2_6, chi2_7, chi2_8, chi2_9 };
  Int_t nIter_vec[10] = { nIter_0, nIter_1, nIter_2, nIter_3, nIter_4, nIter_5, nIter_6, nIter_7, nIter_8, nIter_9 };

  std::vector<TVectorD> X_vec;
  X_vec.push_back(X0);
  X_vec.push_back(X1);
  X_vec.push_back(X2);
  X_vec.push_back(X3);
  X_vec.push_back(X4);
  X_vec.push_back(X5);
  X_vec.push_back(X6);
  X_vec.push_back(X7);
  X_vec.push_back(X8);
  X_vec.push_back(X9);

  std::vector<TVectorD> X_ERR_vec;
  X_ERR_vec.push_back(X0_ERR);
  X_ERR_vec.push_back(X1_ERR);
  X_ERR_vec.push_back(X2_ERR);
  X_ERR_vec.push_back(X3_ERR);
  X_ERR_vec.push_back(X4_ERR);
  X_ERR_vec.push_back(X5_ERR);
  X_ERR_vec.push_back(X6_ERR);
  X_ERR_vec.push_back(X7_ERR);
  X_ERR_vec.push_back(X8_ERR);
  X_ERR_vec.push_back(X9_ERR);

  bool all_fail = false;
  if( (status0 > 1) && (status1 > 1) && (status2 > 1) && (status3 > 1) && (status4 > 1) && (status5 > 1) && (status6 > 1) && (status7 > 1) && (status8 > 1) && (status9 > 1) )
  {
      all_fail = true;
  }
  
  Double_t chi2_min = 100000000000000000;
  Int_t i_min = 0;

  if(all_fail) // If all fail, returns the state of init 0
  {
    i_min = 0;
  }
  else // if one or more passes, from the ones that pass return the one that gives the lowest value for the chi^2
  {
      for(int i = 0; i < 10; i++)
      {
          if( ( (status_vec[i] == 0) || (status_vec[i] == 1) ) && (chi2_vec[i] < chi2_min) )
          {
              chi2_min = chi2_vec[i];
              i_min = i;
          }
      }   
  }
  
  init = i_min;
  status = status_vec[i_min];
  MB = MB_vec[i_min];
  MB_err = MB_err_vec[i_min];
  x_m = X_vec[i_min];
  x_err_m = X_ERR_vec[i_min];
  nIter = nIter_vec[i_min];
  chi2 = chi2_vec[i_min];

}


TVectorD x_initial_estimate_0( TVectorD m ) // Original initialisation for x (based on Anne Keune's thesis)
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );

  Double_t EK  = sqrt( pow(mkaon,2) + pK.Mag2() );
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

TVectorD x_initial_estimate_1( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; taus direction based on vertices
{

  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9); 
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16); 
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

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
  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
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

TVectorD x_initial_estimate_2( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; taus direction based on visible 3pi momenta
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9); 
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16); 
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );

  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

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
  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
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

TVectorD x_initial_estimate_3( TVectorD m, ROOT::Math::XYZPoint BV ) // RD initialisation
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  // Use the maximum value of theta: neutrino takes the maximum portion of momentum from the tau
  Double_t theta1 = asin( ( pow(mtau,2) - pow(m3pi1,2) )/( 2*mtau*sqrt( p3pi1.Mag2() ) ) );
  Double_t theta2 = asin( ( pow(mtau,2) - pow(m3pi2,2) )/( 2*mtau*sqrt( p3pi2.Mag2() ) ) );
  
  Double_t ptau1_mag = ( (pow(mtau,2) + pow(m3pi1,2))*sqrt(p3pi1.Mag2())*cos(theta1) )/( 2*( pow(E3pi1,2) - p3pi1.Mag2()*pow(cos(theta1),2) ) );
  Double_t ptau2_mag = ( (pow(mtau,2) + pow(m3pi2,2))*sqrt(p3pi2.Mag2())*cos(theta2) )/( 2*( pow(E3pi2,2) - p3pi2.Mag2()*pow(cos(theta2),2) ) );

  ROOT::Math::XYZVector ptau1 = ptau1_mag*u1;
  ROOT::Math::XYZVector ptau2 = ptau2_mag*u2;

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
  Double_t EB = EK + Etau1 + Etau2;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_4( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; tau+ direction from vertices ; tau- direction from pions
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  // Magnitude of 3pi momenta
  Double_t p3pi1_mag = p3pi1.r();
  Double_t p3pi2_mag = p3pi2.r();

  // B+ flight direction
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  // Get K+ momentum perpendicular to the B+ flight direction
  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  // Get tau flight directions (tau+ from vertices, tau- from pions)
  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() );
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
  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
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
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t EB = Etau1 + Etau2 + EK;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_5( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; tau- direction from vertices ; tau+ direction from pions
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  // Magnitude of 3pi momenta
  Double_t p3pi1_mag = p3pi1.r();
  Double_t p3pi2_mag = p3pi2.r();

  // B+ flight direction
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  // Get K+ momentum perpendicular to the B+ flight direction
  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  // Get tau flight directions (tau+ from vertices, tau- from pions)
  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() );
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
  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
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
  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t EB = Etau1 + Etau2 + EK;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_6( TVectorD m, ROOT::Math::XYZPoint BV ) // Mixes Marseille (tau+) and K* tau tau vertices (tau-)
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  //////////////////////////////////////////////////////////////////////////////////////////// tau+ is Marseille
  Double_t theta1 = asin( ( pow(mtau,2) - pow(m3pi1,2) )/( 2*mtau*sqrt( p3pi1.Mag2() ) ) );

  Double_t ptau1_mag = ( (pow(mtau,2) + pow(m3pi1,2))*sqrt(p3pi1.Mag2())*cos(theta1) )/( 2*( pow(E3pi1,2) - p3pi1.Mag2()*pow(cos(theta1),2) ) );

  ROOT::Math::XYZVector ptau1 = ptau1_mag*u1;

  ////////////////////////////////////////////////////////////////////////////////////////// tau- is K*tautau vertices
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() );
  tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() );

  tau_dir1 = tau_dir1.Unit();
  tau_dir2 = tau_dir2.Unit();

  ROOT::Math::XYZVector tau_dir1_perp = ( tau_dir1 - (tau_dir1.Dot(bDir))*bDir ).Unit();
  ROOT::Math::XYZVector tau_dir2_perp = ( tau_dir2 - (tau_dir2.Dot(bDir))*bDir ).Unit();

  Double_t cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit());
  Double_t cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit());

  Double_t phi1 = TMath::ACos(cosphi1);
  Double_t phi2 = TMath::ACos(cosphi2);

  ROOT::Math::XYZVector tau_perp1_perpK = tau_dir1_perp - ( tau_dir1_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );
  ROOT::Math::XYZVector tau_perp2_perpK = tau_dir2_perp - ( tau_dir2_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );

  Double_t tau_perp_ratio = (tau_perp1_perpK.R())/(tau_perp2_perpK.R());

  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
  Double_t pMag_tau2_perp = pMag_tau1_perp*tau_perp_ratio;

  ROOT::Math::XYZVector p_tau1_perp = pMag_tau1_perp*tau_dir1_perp;
  ROOT::Math::XYZVector p_tau2_perp = pMag_tau2_perp*tau_dir2_perp;

  Double_t tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit());
  Double_t tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit());

  Double_t pMag_tau1_long = fabs(pMag_tau1_perp)*tau_B_cos1/sqrt( 1 - pow(tau_B_cos1,2) );
  Double_t pMag_tau2_long = fabs(pMag_tau2_perp)*tau_B_cos2/sqrt( 1 - pow(tau_B_cos2,2) );

  ROOT::Math::XYZVector ptau2 = p_tau2_perp + (pMag_tau2_long*bDir);

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
  Double_t EB = EK + Etau1 + Etau2;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_7( TVectorD m, ROOT::Math::XYZPoint BV ) // Mixes Marseille (tau-) and K* tau tau vertices (tau+)
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  //////////////////////////////////////////////////////////////////////////////////////////// tau+ is Marseille
  Double_t theta2 = asin( ( pow(mtau,2) - pow(m3pi2,2) )/( 2*mtau*sqrt( p3pi2.Mag2() ) ) );

  Double_t ptau2_mag = ( (pow(mtau,2) + pow(m3pi2,2))*sqrt(p3pi2.Mag2())*cos(theta2) )/( 2*( pow(E3pi2,2) - p3pi2.Mag2()*pow(cos(theta2),2) ) );

  ROOT::Math::XYZVector ptau2 = ptau2_mag*u2;

  ////////////////////////////////////////////////////////////////////////////////////////// tau- is K*tautau vertices
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() );
  tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() );

  tau_dir1 = tau_dir1.Unit();
  tau_dir2 = tau_dir2.Unit();

  ROOT::Math::XYZVector tau_dir1_perp = ( tau_dir1 - (tau_dir1.Dot(bDir))*bDir ).Unit();
  ROOT::Math::XYZVector tau_dir2_perp = ( tau_dir2 - (tau_dir2.Dot(bDir))*bDir ).Unit();

  Double_t cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit());
  Double_t cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit());

  Double_t phi1 = TMath::ACos(cosphi1);
  Double_t phi2 = TMath::ACos(cosphi2);

  ROOT::Math::XYZVector tau_perp1_perpK = tau_dir1_perp - ( tau_dir1_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );
  ROOT::Math::XYZVector tau_perp2_perpK = tau_dir2_perp - ( tau_dir2_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );

  Double_t tau_perp_ratio = (tau_perp1_perpK.R())/(tau_perp2_perpK.R());

  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
  Double_t pMag_tau2_perp = pMag_tau1_perp*tau_perp_ratio;

  ROOT::Math::XYZVector p_tau1_perp = pMag_tau1_perp*tau_dir1_perp;
  ROOT::Math::XYZVector p_tau2_perp = pMag_tau2_perp*tau_dir2_perp;

  Double_t tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit());
  Double_t tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit());

  Double_t pMag_tau1_long = fabs(pMag_tau1_perp)*tau_B_cos1/sqrt( 1 - pow(tau_B_cos1,2) );
  Double_t pMag_tau2_long = fabs(pMag_tau2_perp)*tau_B_cos2/sqrt( 1 - pow(tau_B_cos2,2) );

  ROOT::Math::XYZVector ptau1 = p_tau1_perp + (pMag_tau1_long*bDir);

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
  Double_t EB = EK + Etau1 + Etau2;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_8( TVectorD m, ROOT::Math::XYZPoint BV ) // Mixes Marseille (tau+) and K* tau tau pions (tau-)
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  //////////////////////////////////////////////////////////////////////////////////////////// tau+ is Marseille
  Double_t theta1 = asin( ( pow(mtau,2) - pow(m3pi1,2) )/( 2*mtau*sqrt( p3pi1.Mag2() ) ) );

  Double_t ptau1_mag = ( (pow(mtau,2) + pow(m3pi1,2))*sqrt(p3pi1.Mag2())*cos(theta1) )/( 2*( pow(E3pi1,2) - p3pi1.Mag2()*pow(cos(theta1),2) ) );

  ROOT::Math::XYZVector ptau1 = ptau1_mag*u1;

  ////////////////////////////////////////////////////////////////////////////////////////// tau- is K*tautau vertices
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() );
  tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() );

  tau_dir1 = tau_dir1.Unit();
  tau_dir2 = tau_dir2.Unit();

  ROOT::Math::XYZVector tau_dir1_perp = ( tau_dir1 - (tau_dir1.Dot(bDir))*bDir ).Unit();
  ROOT::Math::XYZVector tau_dir2_perp = ( tau_dir2 - (tau_dir2.Dot(bDir))*bDir ).Unit();

  Double_t cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit());
  Double_t cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit());

  Double_t phi1 = TMath::ACos(cosphi1);
  Double_t phi2 = TMath::ACos(cosphi2);

  ROOT::Math::XYZVector tau_perp1_perpK = tau_dir1_perp - ( tau_dir1_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );
  ROOT::Math::XYZVector tau_perp2_perpK = tau_dir2_perp - ( tau_dir2_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );

  Double_t tau_perp_ratio = (tau_perp1_perpK.R())/(tau_perp2_perpK.R());

  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
  Double_t pMag_tau2_perp = pMag_tau1_perp*tau_perp_ratio;

  ROOT::Math::XYZVector p_tau1_perp = pMag_tau1_perp*tau_dir1_perp;
  ROOT::Math::XYZVector p_tau2_perp = pMag_tau2_perp*tau_dir2_perp;

  Double_t tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit());
  Double_t tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit());

  Double_t pMag_tau1_long = fabs(pMag_tau1_perp)*tau_B_cos1/sqrt( 1 - pow(tau_B_cos1,2) );
  Double_t pMag_tau2_long = fabs(pMag_tau2_perp)*tau_B_cos2/sqrt( 1 - pow(tau_B_cos2,2) );

  ROOT::Math::XYZVector ptau2 = p_tau2_perp + (pMag_tau2_long*bDir);

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
  Double_t EB = EK + Etau1 + Etau2;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_9( TVectorD m, ROOT::Math::XYZPoint BV ) // Mixes Marseille (tau-) and K* tau tau pions (tau+)
{
  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
  Double_t E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );    
  Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

  ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
  ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

  Double_t m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  Double_t m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  //////////////////////////////////////////////////////////////////////////////////////////// tau+ is Marseille
  Double_t theta2 = asin( ( pow(mtau,2) - pow(m3pi2,2) )/( 2*mtau*sqrt( p3pi2.Mag2() ) ) );

  Double_t ptau2_mag = ( (pow(mtau,2) + pow(m3pi2,2))*sqrt(p3pi2.Mag2())*cos(theta2) )/( 2*( pow(E3pi2,2) - p3pi2.Mag2()*pow(cos(theta2),2) ) );

  ROOT::Math::XYZVector ptau2 = ptau2_mag*u2;

  ////////////////////////////////////////////////////////////////////////////////////////// tau- is K*tautau vertices
  ROOT::Math::XYZVector bDir = (BV - PV).Unit();

  ROOT::Math::XYZVector pK_perp = pK - (pK.Dot(bDir))*bDir;

  ROOT::Math::XYZVector tau_dir1, tau_dir2;
  tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() );
  tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() );

  tau_dir1 = tau_dir1.Unit();
  tau_dir2 = tau_dir2.Unit();

  ROOT::Math::XYZVector tau_dir1_perp = ( tau_dir1 - (tau_dir1.Dot(bDir))*bDir ).Unit();
  ROOT::Math::XYZVector tau_dir2_perp = ( tau_dir2 - (tau_dir2.Dot(bDir))*bDir ).Unit();

  Double_t cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit());
  Double_t cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit());

  Double_t phi1 = TMath::ACos(cosphi1);
  Double_t phi2 = TMath::ACos(cosphi2);

  ROOT::Math::XYZVector tau_perp1_perpK = tau_dir1_perp - ( tau_dir1_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );
  ROOT::Math::XYZVector tau_perp2_perpK = tau_dir2_perp - ( tau_dir2_perp.Dot(pK_perp.Unit())*pK_perp.Unit() );

  Double_t tau_perp_ratio = (tau_perp1_perpK.R())/(tau_perp2_perpK.R());

  Double_t pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(sin(phi1)/sin(phi2))) );
  Double_t pMag_tau2_perp = pMag_tau1_perp*tau_perp_ratio;

  ROOT::Math::XYZVector p_tau1_perp = pMag_tau1_perp*tau_dir1_perp;
  ROOT::Math::XYZVector p_tau2_perp = pMag_tau2_perp*tau_dir2_perp;

  Double_t tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit());
  Double_t tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit());

  Double_t pMag_tau1_long = fabs(pMag_tau1_perp)*tau_B_cos1/sqrt( 1 - pow(tau_B_cos1,2) );
  Double_t pMag_tau2_long = fabs(pMag_tau2_perp)*tau_B_cos2/sqrt( 1 - pow(tau_B_cos2,2) );

  ROOT::Math::XYZVector ptau1 = p_tau1_perp + (pMag_tau1_long*bDir);

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;

  Double_t Enu1 = sqrt( pnu1.Mag2() );
  Double_t Enu2 = sqrt( pnu2.Mag2() );

  Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

  ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
  Double_t EB = EK + Etau1 + Etau2;
  Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

TVectorD x_initial_estimate_MLP( TVectorD m, ROOT::Math::XYZPoint BV, Int_t species )
{
    TMVA::Tools::Instance();
    TMVA::Reader *reader_taup_PX = new TMVA::Reader( "V:Color:Silent" ); 
    TMVA::Reader *reader_taup_PY = new TMVA::Reader( "V:Color:Silent" ); 
    TMVA::Reader *reader_taup_PZ = new TMVA::Reader( "V:Color:Silent" ); 
    TMVA::Reader *reader_taum_PX = new TMVA::Reader( "V:Color:Silent" ); 
    TMVA::Reader *reader_taum_PY = new TMVA::Reader( "V:Color:Silent" ); 
    TMVA::Reader *reader_taum_PZ = new TMVA::Reader( "V:Color:Silent" ); 
    
    TString weightfile_taup_PX;
    TString weightfile_taup_PY;
    TString weightfile_taup_PZ;
    TString weightfile_taum_PX;
    TString weightfile_taum_PY;
    TString weightfile_taum_PZ;

    if((species == 10) || (species == 11) || (species == 12) || (species == 1) || (species == 2) || (species ==3))
    {
      weightfile_taup_PX = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP.weights.xml";
      weightfile_taup_PY = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP.weights.xml";
      weightfile_taup_PZ = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP.weights.xml";
      weightfile_taum_PX = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP.weights.xml";
      weightfile_taum_PY = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP.weights.xml";
      weightfile_taum_PZ = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP.weights.xml";
      // weightfile_taup_PX = "/home/felician/B2Ktautau/MLP_init_weights/Ktautau/TMVARegression_taup_TRUEP_X_MLP.weights.xml";
      // weightfile_taup_PY = "/home/felician/B2Ktautau/MLP_init_weights/Ktautau/TMVARegression_taup_TRUEP_Y_MLP.weights.xml";
      // weightfile_taup_PZ = "/home/felician/B2Ktautau/MLP_init_weights/Ktautau/TMVARegression_taup_TRUEP_Z_MLP.weights.xml";
      // weightfile_taum_PX = "/home/felician/B2Ktautau/MLP_init_weights/Ktautau/TMVARegression_taum_TRUEP_X_MLP.weights.xml";
      // weightfile_taum_PY = "/home/felician/B2Ktautau/MLP_init_weights/Ktautau/TMVARegression_taum_TRUEP_Y_MLP.weights.xml";
      // weightfile_taum_PZ = "/home/felician/B2Ktautau/MLP_init_weights/Ktautau/TMVARegression_taum_TRUEP_Z_MLP.weights.xml";
    }
    else if(species == 4)
    {
      weightfile_taup_PX = "/home/felician/B2Ktautau/MLP_init_weights/DDK/TMVARegression_Dp_TRUEP_X_MLP.weights.xml";
      weightfile_taup_PY = "/home/felician/B2Ktautau/MLP_init_weights/DDK/TMVARegression_Dp_TRUEP_Y_MLP.weights.xml";
      weightfile_taup_PZ = "/home/felician/B2Ktautau/MLP_init_weights/DDK/TMVARegression_Dp_TRUEP_Z_MLP.weights.xml";
      weightfile_taum_PX = "/home/felician/B2Ktautau/MLP_init_weights/DDK/TMVARegression_Dm_TRUEP_X_MLP.weights.xml";
      weightfile_taum_PY = "/home/felician/B2Ktautau/MLP_init_weights/DDK/TMVARegression_Dm_TRUEP_Y_MLP.weights.xml";
      weightfile_taum_PZ = "/home/felician/B2Ktautau/MLP_init_weights/DDK/TMVARegression_Dm_TRUEP_Z_MLP.weights.xml";
    }

    Float_t MLP_m_vars[dimM];
    Float_t MLP_RPz;
    Float_t MLP_eventNumber;

    reader_taup_PX->AddVariable("df_m_1", &MLP_m_vars[0]);
    reader_taup_PX->AddVariable("df_m_2", &MLP_m_vars[1]);
    reader_taup_PX->AddVariable("df_m_3", &MLP_m_vars[2]);
    reader_taup_PX->AddVariable("df_m_4", &MLP_m_vars[3]);
    reader_taup_PX->AddVariable("df_m_5", &MLP_m_vars[4]);
    reader_taup_PX->AddVariable("df_m_6", &MLP_m_vars[5]);
    reader_taup_PX->AddVariable("df_m_7", &MLP_m_vars[6]);
    reader_taup_PX->AddVariable("df_m_8", &MLP_m_vars[7]);
    reader_taup_PX->AddVariable("df_m_9", &MLP_m_vars[8]);
    reader_taup_PX->AddVariable("df_m_10", &MLP_m_vars[9]);
    reader_taup_PX->AddVariable("df_m_11", &MLP_m_vars[10]);
    reader_taup_PX->AddVariable("df_m_12", &MLP_m_vars[11]);
    reader_taup_PX->AddVariable("df_m_13", &MLP_m_vars[12]);
    reader_taup_PX->AddVariable("df_m_14", &MLP_m_vars[13]);
    reader_taup_PX->AddVariable("df_m_15", &MLP_m_vars[14]);
    reader_taup_PX->AddVariable("df_m_16", &MLP_m_vars[15]);
    reader_taup_PX->AddVariable("df_m_17", &MLP_m_vars[16]);
    reader_taup_PX->AddVariable("df_m_18", &MLP_m_vars[17]);
    reader_taup_PX->AddVariable("df_m_19", &MLP_m_vars[18]);
    reader_taup_PX->AddVariable("df_m_20", &MLP_m_vars[19]);
    reader_taup_PX->AddVariable("df_m_21", &MLP_m_vars[20]);
    reader_taup_PX->AddVariable("df_m_22", &MLP_m_vars[21]);
    reader_taup_PX->AddVariable("Kp_RP_Z", &MLP_RPz);
    reader_taup_PX->AddSpectator("eventNumber", &MLP_eventNumber);

    reader_taup_PY->AddVariable("df_m_1", &MLP_m_vars[0]);
    reader_taup_PY->AddVariable("df_m_2", &MLP_m_vars[1]);
    reader_taup_PY->AddVariable("df_m_3", &MLP_m_vars[2]);
    reader_taup_PY->AddVariable("df_m_4", &MLP_m_vars[3]);
    reader_taup_PY->AddVariable("df_m_5", &MLP_m_vars[4]);
    reader_taup_PY->AddVariable("df_m_6", &MLP_m_vars[5]);
    reader_taup_PY->AddVariable("df_m_7", &MLP_m_vars[6]);
    reader_taup_PY->AddVariable("df_m_8", &MLP_m_vars[7]);
    reader_taup_PY->AddVariable("df_m_9", &MLP_m_vars[8]);
    reader_taup_PY->AddVariable("df_m_10", &MLP_m_vars[9]);
    reader_taup_PY->AddVariable("df_m_11", &MLP_m_vars[10]);
    reader_taup_PY->AddVariable("df_m_12", &MLP_m_vars[11]);
    reader_taup_PY->AddVariable("df_m_13", &MLP_m_vars[12]);
    reader_taup_PY->AddVariable("df_m_14", &MLP_m_vars[13]);
    reader_taup_PY->AddVariable("df_m_15", &MLP_m_vars[14]);
    reader_taup_PY->AddVariable("df_m_16", &MLP_m_vars[15]);
    reader_taup_PY->AddVariable("df_m_17", &MLP_m_vars[16]);
    reader_taup_PY->AddVariable("df_m_18", &MLP_m_vars[17]);
    reader_taup_PY->AddVariable("df_m_19", &MLP_m_vars[18]);
    reader_taup_PY->AddVariable("df_m_20", &MLP_m_vars[19]);
    reader_taup_PY->AddVariable("df_m_21", &MLP_m_vars[20]);
    reader_taup_PY->AddVariable("df_m_22", &MLP_m_vars[21]);
    reader_taup_PY->AddVariable("Kp_RP_Z", &MLP_RPz);
    reader_taup_PY->AddSpectator("eventNumber", &MLP_eventNumber);

    reader_taup_PZ->AddVariable("df_m_1", &MLP_m_vars[0]);
    reader_taup_PZ->AddVariable("df_m_2", &MLP_m_vars[1]);
    reader_taup_PZ->AddVariable("df_m_3", &MLP_m_vars[2]);
    reader_taup_PZ->AddVariable("df_m_4", &MLP_m_vars[3]);
    reader_taup_PZ->AddVariable("df_m_5", &MLP_m_vars[4]);
    reader_taup_PZ->AddVariable("df_m_6", &MLP_m_vars[5]);
    reader_taup_PZ->AddVariable("df_m_7", &MLP_m_vars[6]);
    reader_taup_PZ->AddVariable("df_m_8", &MLP_m_vars[7]);
    reader_taup_PZ->AddVariable("df_m_9", &MLP_m_vars[8]);
    reader_taup_PZ->AddVariable("df_m_10", &MLP_m_vars[9]);
    reader_taup_PZ->AddVariable("df_m_11", &MLP_m_vars[10]);
    reader_taup_PZ->AddVariable("df_m_12", &MLP_m_vars[11]);
    reader_taup_PZ->AddVariable("df_m_13", &MLP_m_vars[12]);
    reader_taup_PZ->AddVariable("df_m_14", &MLP_m_vars[13]);
    reader_taup_PZ->AddVariable("df_m_15", &MLP_m_vars[14]);
    reader_taup_PZ->AddVariable("df_m_16", &MLP_m_vars[15]);
    reader_taup_PZ->AddVariable("df_m_17", &MLP_m_vars[16]);
    reader_taup_PZ->AddVariable("df_m_18", &MLP_m_vars[17]);
    reader_taup_PZ->AddVariable("df_m_19", &MLP_m_vars[18]);
    reader_taup_PZ->AddVariable("df_m_20", &MLP_m_vars[19]);
    reader_taup_PZ->AddVariable("df_m_21", &MLP_m_vars[20]);
    reader_taup_PZ->AddVariable("df_m_22", &MLP_m_vars[21]);
    reader_taup_PZ->AddVariable("Kp_RP_Z", &MLP_RPz);
    reader_taup_PZ->AddSpectator("eventNumber", &MLP_eventNumber);

    reader_taum_PX->AddVariable("df_m_1", &MLP_m_vars[0]);
    reader_taum_PX->AddVariable("df_m_2", &MLP_m_vars[1]);
    reader_taum_PX->AddVariable("df_m_3", &MLP_m_vars[2]);
    reader_taum_PX->AddVariable("df_m_4", &MLP_m_vars[3]);
    reader_taum_PX->AddVariable("df_m_5", &MLP_m_vars[4]);
    reader_taum_PX->AddVariable("df_m_6", &MLP_m_vars[5]);
    reader_taum_PX->AddVariable("df_m_7", &MLP_m_vars[6]);
    reader_taum_PX->AddVariable("df_m_8", &MLP_m_vars[7]);
    reader_taum_PX->AddVariable("df_m_9", &MLP_m_vars[8]);
    reader_taum_PX->AddVariable("df_m_10", &MLP_m_vars[9]);
    reader_taum_PX->AddVariable("df_m_11", &MLP_m_vars[10]);
    reader_taum_PX->AddVariable("df_m_12", &MLP_m_vars[11]);
    reader_taum_PX->AddVariable("df_m_13", &MLP_m_vars[12]);
    reader_taum_PX->AddVariable("df_m_14", &MLP_m_vars[13]);
    reader_taum_PX->AddVariable("df_m_15", &MLP_m_vars[14]);
    reader_taum_PX->AddVariable("df_m_16", &MLP_m_vars[15]);
    reader_taum_PX->AddVariable("df_m_17", &MLP_m_vars[16]);
    reader_taum_PX->AddVariable("df_m_18", &MLP_m_vars[17]);
    reader_taum_PX->AddVariable("df_m_19", &MLP_m_vars[18]);
    reader_taum_PX->AddVariable("df_m_20", &MLP_m_vars[19]);
    reader_taum_PX->AddVariable("df_m_21", &MLP_m_vars[20]);
    reader_taum_PX->AddVariable("df_m_22", &MLP_m_vars[21]);
    reader_taum_PX->AddVariable("Kp_RP_Z", &MLP_RPz);
    reader_taum_PX->AddSpectator("eventNumber", &MLP_eventNumber);

    reader_taum_PY->AddVariable("df_m_1", &MLP_m_vars[0]);
    reader_taum_PY->AddVariable("df_m_2", &MLP_m_vars[1]);
    reader_taum_PY->AddVariable("df_m_3", &MLP_m_vars[2]);
    reader_taum_PY->AddVariable("df_m_4", &MLP_m_vars[3]);
    reader_taum_PY->AddVariable("df_m_5", &MLP_m_vars[4]);
    reader_taum_PY->AddVariable("df_m_6", &MLP_m_vars[5]);
    reader_taum_PY->AddVariable("df_m_7", &MLP_m_vars[6]);
    reader_taum_PY->AddVariable("df_m_8", &MLP_m_vars[7]);
    reader_taum_PY->AddVariable("df_m_9", &MLP_m_vars[8]);
    reader_taum_PY->AddVariable("df_m_10", &MLP_m_vars[9]);
    reader_taum_PY->AddVariable("df_m_11", &MLP_m_vars[10]);
    reader_taum_PY->AddVariable("df_m_12", &MLP_m_vars[11]);
    reader_taum_PY->AddVariable("df_m_13", &MLP_m_vars[12]);
    reader_taum_PY->AddVariable("df_m_14", &MLP_m_vars[13]);
    reader_taum_PY->AddVariable("df_m_15", &MLP_m_vars[14]);
    reader_taum_PY->AddVariable("df_m_16", &MLP_m_vars[15]);
    reader_taum_PY->AddVariable("df_m_17", &MLP_m_vars[16]);
    reader_taum_PY->AddVariable("df_m_18", &MLP_m_vars[17]);
    reader_taum_PY->AddVariable("df_m_19", &MLP_m_vars[18]);
    reader_taum_PY->AddVariable("df_m_20", &MLP_m_vars[19]);
    reader_taum_PY->AddVariable("df_m_21", &MLP_m_vars[20]);
    reader_taum_PY->AddVariable("df_m_22", &MLP_m_vars[21]);
    reader_taum_PY->AddVariable("Kp_RP_Z", &MLP_RPz);
    reader_taum_PY->AddSpectator("eventNumber", &MLP_eventNumber);

    reader_taum_PZ->AddVariable("df_m_1", &MLP_m_vars[0]);
    reader_taum_PZ->AddVariable("df_m_2", &MLP_m_vars[1]);
    reader_taum_PZ->AddVariable("df_m_3", &MLP_m_vars[2]);
    reader_taum_PZ->AddVariable("df_m_4", &MLP_m_vars[3]);
    reader_taum_PZ->AddVariable("df_m_5", &MLP_m_vars[4]);
    reader_taum_PZ->AddVariable("df_m_6", &MLP_m_vars[5]);
    reader_taum_PZ->AddVariable("df_m_7", &MLP_m_vars[6]);
    reader_taum_PZ->AddVariable("df_m_8", &MLP_m_vars[7]);
    reader_taum_PZ->AddVariable("df_m_9", &MLP_m_vars[8]);
    reader_taum_PZ->AddVariable("df_m_10", &MLP_m_vars[9]);
    reader_taum_PZ->AddVariable("df_m_11", &MLP_m_vars[10]);
    reader_taum_PZ->AddVariable("df_m_12", &MLP_m_vars[11]);
    reader_taum_PZ->AddVariable("df_m_13", &MLP_m_vars[12]);
    reader_taum_PZ->AddVariable("df_m_14", &MLP_m_vars[13]);
    reader_taum_PZ->AddVariable("df_m_15", &MLP_m_vars[14]);
    reader_taum_PZ->AddVariable("df_m_16", &MLP_m_vars[15]);
    reader_taum_PZ->AddVariable("df_m_17", &MLP_m_vars[16]);
    reader_taum_PZ->AddVariable("df_m_18", &MLP_m_vars[17]);
    reader_taum_PZ->AddVariable("df_m_19", &MLP_m_vars[18]);
    reader_taum_PZ->AddVariable("df_m_20", &MLP_m_vars[19]);
    reader_taum_PZ->AddVariable("df_m_21", &MLP_m_vars[20]);
    reader_taum_PZ->AddVariable("df_m_22", &MLP_m_vars[21]);
    reader_taum_PZ->AddVariable("Kp_RP_Z", &MLP_RPz);
    reader_taum_PZ->AddSpectator("eventNumber", &MLP_eventNumber);

    reader_taup_PX->BookMVA("MLP", weightfile_taup_PX);
    reader_taup_PY->BookMVA("MLP", weightfile_taup_PY);
    reader_taup_PZ->BookMVA("MLP", weightfile_taup_PZ);
    reader_taum_PX->BookMVA("MLP", weightfile_taum_PX);
    reader_taum_PY->BookMVA("MLP", weightfile_taum_PY);
    reader_taum_PZ->BookMVA("MLP", weightfile_taum_PZ);

    MLP_m_vars[0] = m(0);
    MLP_m_vars[1] = m(1);
    MLP_m_vars[2] = m(2);
    MLP_m_vars[3] = m(3);
    MLP_m_vars[4] = m(4);
    MLP_m_vars[5] = m(5);
    MLP_m_vars[6] = m(6);
    MLP_m_vars[7] = m(7);
    MLP_m_vars[8] = m(8);
    MLP_m_vars[9] = m(9);
    MLP_m_vars[10] = m(10);
    MLP_m_vars[11] = m(11);
    MLP_m_vars[12] = m(12);
    MLP_m_vars[13] = m(13);
    MLP_m_vars[14] = m(14);
    MLP_m_vars[15] = m(15);
    MLP_m_vars[16] = m(16);
    MLP_m_vars[17] = m(17);
    MLP_m_vars[18] = m(18);
    MLP_m_vars[19] = m(19);
    MLP_m_vars[20] = m(20);
    MLP_m_vars[21] = m(21);
    MLP_RPz = RPz;
    MLP_eventNumber = eventNumber;

    Float_t MLP_taup_PX = reader_taup_PX->EvaluateRegression("MLP")[0];
    Float_t MLP_taup_PY = reader_taup_PY->EvaluateRegression("MLP")[0];
    Float_t MLP_taup_PZ = reader_taup_PZ->EvaluateRegression("MLP")[0];
    Float_t MLP_taum_PX = reader_taum_PX->EvaluateRegression("MLP")[0];
    Float_t MLP_taum_PY = reader_taum_PY->EvaluateRegression("MLP")[0];
    Float_t MLP_taum_PZ = reader_taum_PZ->EvaluateRegression("MLP")[0];

    ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
    ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
    ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) ); 
    Double_t E3pi1 = m(9);
    ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
    ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
    Double_t E3pi2 = m(16);
    ROOT::Math::XYZPoint RP( m(17), m(18), RPz );
    ROOT::Math::XYZVector pK( m(19), m(20), m(21) );  
    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

    ROOT::Math::XYZVector ptau1( MLP_taup_PX, MLP_taup_PY, MLP_taup_PZ );
    ROOT::Math::XYZVector ptau2( MLP_taum_PX, MLP_taum_PY, MLP_taum_PZ );
    ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
    ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );
    ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
    Double_t EB = Etau1 + Etau2 + EK;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

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

  for(int i = 0; i < dimM; i++)
  {
    mp(i) = m(i);
  }

  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi1 = m(9);
  Double_t E3pi2 = m(16);

  mp(9) = pow(E3pi1,2) - p3pi1.Mag2();
  mp(16) = pow(E3pi2,2) - p3pi2.Mag2();

  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );
  mp(21) = pow(pK.x(), 2) + pow( pK.y(), 2) + pow( pK.z(), 2);

  return mp;

}

TMatrixDSym transform_V( TVectorD m, TMatrixDSym V )
{
  TMatrixD Vp( dimM, dimM );


  // Jacobian matrix
  TMatrixD J( dimM, dimM );
  for(int i = 0; i < dimM; i++)
  {
    for(int j = 0; j < dimM; j++)
    {
      if(i == j)
      {
        J(i,j) = 1.;
      }
      else
      {
        J(i,j) = 0.;
      }
    }
  }

  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  Double_t E3pi1 = m(9);
  Double_t E3pi2 = m(16);

  // m3pi1^2
  J(9,6) = -2*p3pi1.x();
  J(9,7) = -2*p3pi1.y();
  J(9,8) = -2*p3pi1.z();
  J(9,9) = 2*E3pi1;

  // m3pi2^2
  J(16,13) = -2*p3pi2.x();
  J(16,14) = -2*p3pi2.y();
  J(16,15) = -2*p3pi2.z();
  J(16,16) = 2*E3pi2;

  ROOT::Math::XYZVector pK( m(19), m(20), m(21) );

  // pK^2
  J(21,19) = 2*pK.x();
  J(21,20) = 2*pK.y();
  J(21,21) = 2*pK.z();

  // J.Print();

  TMatrixD J_transpose = J;
  J_transpose.T();

  Vp = J*(V*J_transpose);

  TMatrixDSym Vp_sym(dimM);
  // Vp_sym.SetTol(pow(10,-23));
  for(int i = 0; i < dimM; i++)
  {
    for(int j = 0; j < dimM; j++)
    {
      Vp_sym(i,j) = Vp(i,j);
    }
  }

  // Vp_sym(22,22) = pow(10,-12);

  return Vp_sym;
}

TMatrixDSym scale_tau_DVs_uncertainties(TMatrixDSym V, Double_t factor)
{
  // scale uncertainties on the tau DVs

  // 1) Find correlation matrix
  TMatrixD V_corr(dimM,dimM);
  for(int i = 0; i < dimM; i++)
  {
    for(int j = 0; j < dimM; j++)
    {
      V_corr(i,j) = V(i,j)/sqrt( V(i,i)*V(j,j) );
    }
  }

  // tau+ DV (3,4,5)
  for(int i = 3; i < 6; i++)
  {
    Double_t sigma_i = sqrt(V(i,i));
    sigma_i += factor*sigma_i;

    for(int j = 3; j < 6; j++)
    {
      Double_t sigma_j = sqrt(V(j,j));
      sigma_j += factor*sigma_j;
      

      V(i,j) = V_corr(i,j)*sigma_i*sigma_j;
    }
  }

  // tau- DV (10,11,12)
  for(int i = 10; i < 13; i++)
  {
    Double_t sigma_i = sqrt(V(i,i));
    sigma_i += factor*sigma_i;

    for(int j = 10; j < 13; j++)
    {
      Double_t sigma_j = sqrt(V(j,j));
      sigma_j += factor*sigma_j;
      

      V(i,j) = V_corr(i,j)*sigma_i*sigma_j;
    }
  }

  return V;
}

TMatrixDSym scale_uncertainties(TMatrixDSym V, Double_t factor)
{
  // scale uncertainties on the tau DVs

  // 1) Find correlation matrix
  TMatrixD V_corr(dimM,dimM);
  for(int i = 0; i < dimM; i++)
  {
    for(int j = 0; j < dimM; j++)
    {
      V_corr(i,j) = V(i,j)/sqrt( V(i,i)*V(j,j) );
    }
  }

  for(int i = 0; i < dimM; i++)
  {
    Double_t sigma_i = sqrt(V(i,i));
    sigma_i += factor*sigma_i;

    for(int j = 0; j < dimM; j++)
    {
      Double_t sigma_j = sqrt(V(j,j));
      sigma_j += factor*sigma_j;

      V(i,j) = V_corr(i,j)*sigma_i*sigma_j;
    }
  }

  return V;
}

TMatrixDSym invert_sum(TMatrixDSym V1, TMatrixDSym V2)
{
  // V1 has diagonal blocks
  // V2 has off-diagonal blocks

  TMatrixDSym V_inverse(dimM);

  // 1) Write V2 has the sum of rank 1 matrices
  // For this we will use the SVD decomposition
  // SVD decomposition: V2 = U S V^T = sum_i s_i u_i v_i^T (= sum of rank 1 matrices)
  TDecompSVD *svd = new TDecompSVD(V2);  
  TMatrixD U = svd->GetU(); // orthogonal matrix
  TVectorD S = svd->GetSig(); // diagonal matrix (vector w/ diagonal entries)
  TMatrixD V = svd->GetV(); // orthogonal matrix

  // 2) Invert simpler matrix V1; here inversion works well, V1*V1^-1 ~ I
  TMatrixDSym V1_inverse = V1;
  V1_inverse.Invert();

  TMatrixD Bk(dimM,dimM);
  TMatrixD trace_matrix(dimM,dimM);
  Double_t gk = 0;
  TMatrixD Wk(dimM,dimM);

  for(int k = 0; k < dimM; k++)
  {
    TMatrixD u = U.GetSub(0,dimM-1,k,k); // columns of U
    TMatrixD v = V.GetSub(0,dimM-1,k,k); // columns of V = rows of V^T

    if(k == 0)
    {
      Wk = V1_inverse;
    }
    else
    {
      // TMatrixDSym bk = S(k)*(u*v);
      Bk += S(k)*(u*v.T());

      trace_matrix = Bk*Wk;
      for(int i = 0; i < dimM; i++)
      {
        gk += trace_matrix(i,i);
      }

      Wk -= (1/(1+gk))*(Wk*(Bk*Wk));
    }
  }

  return V_inverse;
}