// Decay fit with 24 constraints to 23 unknowns

using namespace std;

// Global variables 
int dimM = 22; // number of measured parameters
int dimX = 23; // number of unknown parameters
int dimC = 24; // number of exact constraints
int eq_flag = 0; // index of the equation
Double_t initial_tolerance = pow(10,-10);

TVectorD m(dimM);
TMatrixDSym V(dimM);
TMatrixDSym W(dimM);
TVectorD xm(dimM);
TVectorD xu(dimX);
TVectorD lambda(dimC);

TVectorD x0(dimM+dimX+dimC);
Double_t RPz = 0.; // z-component of the RP on the K+ trajectory (fixed)
TVectorD x0_current(dimM+dimX+dimC), x1_current(dimM+dimX+dimC), x2_current(dimM+dimX+dimC), x3_current(dimM+dimX+dimC), x4_current(dimM+dimX+dimC), x5_current(dimM+dimX+dimC), x6_current(dimM+dimX+dimC), x7_current(dimM+dimX+dimC), x8_current(dimM+dimX+dimC), x9_current(dimM+dimX+dimC), x_true_current(dimM+dimX+dimC), x_MLP_current(dimM+dimX+dimC);

// Fit results saved
Double_t MB, MB_err, F_tolerance, chi2, deltaT;
TVectorD X(dimM+dimX+dimC);
TVectorD X_ERR(dimM+dimX+dimC);
TVectorD F(dimM+dimX+dimC);
TMatrixDSym U(dimM+dimX+dimC);
TMatrixDSym U_corr(dimM+dimX+dimC);
Int_t status, init, nIter;
ULong64_t eventNumber;
// Bool_t saddle_point = false;
TVectorD meas_ratio(dimM);

// Constants
Double_t mtau;
Double_t mkaon = 493.677;

// Functions
void solve( ROOT::Math::XYZPoint BV, int init,  ROOT::Math::GSLMultiRootFinder* solver, Double_t tolerance, Int_t species );
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);

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
TVectorD x_initial_estimate_true( TVectorD m );
TVectorD x_initial_estimate_MLP( TVectorD m, ROOT::Math::XYZPoint BV, Int_t species );

void sequence(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder * solver, Double_t tolerance, Int_t species);
void lowest_chi2(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder *solver, Double_t tolerance, Int_t species);

Double_t lagrangian( TVectorD x );
Double_t chisquare( TVectorD xm );
TVectorD exact_constraints( TVectorD x );
TMatrixDSym U_matrix();
vector<double> range(double min, double max, size_t N);
void scan_lagrangian( TVectorD X, Int_t index, Int_t npoints, TString x_name, Int_t year, Int_t species, Int_t line );
void scan_chisquare( TVectorD X, Int_t index, Int_t npoints, TString x_name, Int_t year, Int_t species, Int_t line);

double equations( const double* x_vars );
double eq1( const double* x );
double eq2( const double* x );
double eq3( const double* x );
double eq4( const double* x );
double eq5( const double* x );
double eq6( const double* x );
double eq7( const double* x );
double eq8( const double* x );
double eq9( const double* x );
double eq10( const double* x );
double eq11( const double* x );
double eq12( const double* x );
double eq13( const double* x );
double eq14( const double* x );
double eq15( const double* x );
double eq16( const double* x );
double eq17( const double* x );
double eq18( const double* x );
double eq19( const double* x );
double eq20( const double* x );
double eq21( const double* x );
double eq22( const double* x );
double eq23( const double* x );
double eq24( const double* x );
double eq25( const double* x );
double eq26( const double* x );
double eq27( const double* x );
double eq28( const double* x );
double eq29( const double* x );
double eq30( const double* x );
double eq31( const double* x );
double eq32( const double* x );
double eq33( const double* x );
double eq34( const double* x );
double eq35( const double* x );
double eq36( const double* x );
double eq37( const double* x );
double eq38( const double* x );
double eq39( const double* x );
double eq40( const double* x );
double eq41( const double* x );
double eq42( const double* x );
double eq43( const double* x );
double eq44( const double* x );
double eq45( const double* x );
double eq46( const double* x );
double eq47( const double* x );
double eq48( const double* x );
double eq49( const double* x );
double eq50( const double* x );
double eq51( const double* x );
double eq52( const double* x );
double eq53( const double* x );
double eq54( const double* x );
double eq55( const double* x );
double eq56( const double* x );
double eq57( const double* x );
double eq58( const double* x );
double eq59( const double* x );
double eq60( const double* x );
double eq61( const double* x );
double eq62( const double* x );
double eq63( const double* x );
double eq64( const double* x );
double eq65( const double* x );
double eq66( const double* x );
double eq67( const double* x );
double eq68( const double* x );
double eq69( const double* x );

void DECAY_FIT(int year, TString RECO_files, int species, int line)
{   
    gErrorIgnoreLevel = 2000; // Used to remove "The iteration has converged" comments from ROOT finder

    if( (species == 4) || (species == 5) )
    {
        mtau = 1869.66; // D meson 
    }
    else
    {
        mtau = 1776.86; // tau lepton
    }

    // Retrieve m and V from ntuple
    TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files, 1, line);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    Double_t m_vars[dimM];
    Double_t V_vars[dimM][dimM];
    Double_t BVx, BVy, BVz; // offline estimate for the BV is needed to have a first estimate for the unknown parameters using some of the initialisations

    for(int i = 0; i < dimM; i++)
    {
        t->SetBranchAddress(Form("df_m_%i",i+1), &m_vars[i]);
        for(int j = 0; j < dimM; j++)
        {
        t->SetBranchAddress(Form("df_V_%i_%i",i+1,j+1), &V_vars[i][j]);
        }
    }

    t->SetBranchAddress("Kp_RP_Z", &RPz); 
    t->SetBranchAddress("eventNumber", &eventNumber);
    t->SetBranchAddress("Bp_ENDVERTEX_X", &BVx);
    t->SetBranchAddress("Bp_ENDVERTEX_Y", &BVy);
    t->SetBranchAddress("Bp_ENDVERTEX_Z", &BVz);    

    // Save result of the fit (xm,xu) in a .root file (x = (xm,xu))
    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/%i.root",year,species,line), "RECREATE");
    TTree* tout = new TTree("DecayTree", "DecayTree");

    TString name_x[] = {
        "df_PVx",
        "df_PVy",
        "df_PVz",
        "df_DV1x",
        "df_DV1y",
        "df_DV1z",
        "df_3pi1_PX",
        "df_3pi1_PY",
        "df_3pi1_PZ",
        "df_3pi1_PE",
        "df_DV2x",
        "df_DV2y",
        "df_DV2z",
        "df_3pi2_PX",
        "df_3pi2_PY",
        "df_3pi2_PZ",
        "df_3pi2_PE",
        "df_RPx",
        "df_RPy",
        "df_Kp_PX",
        "df_Kp_PY",
        "df_Kp_PZ",
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
        "df_taup_PE",
        "df_antinutau_PX",
        "df_antinutau_PY",
        "df_antinutau_PZ",
        "df_antinutau_PE",
        "df_taum_PX",
        "df_taum_PY",
        "df_taum_PZ",
        "df_taum_PE",
        "df_nutau_PX",
        "df_nutau_PY",
        "df_nutau_PZ",
        "df_nutau_PE"
    };

    TString name_x_err[] = {
        "df_PVx_err",
        "df_PVy_err",
        "df_PVz_err",
        "df_DV1x_err",
        "df_DV1y_err",
        "df_DV1z_err",
        "df_3pi1_PX_err",
        "df_3pi1_PY_err",
        "df_3pi1_PZ_err",
        "df_3pi1_PE_err",
        "df_DV2x_err",
        "df_DV2y_err",
        "df_DV2z_err",
        "df_3pi2_PX_err",
        "df_3pi2_PY_err",
        "df_3pi2_PZ_err",
        "df_3pi2_PE_err",
        "df_RPx_err",
        "df_RPy_err",
        "df_Kp_PX_err",
        "df_Kp_PY_err",
        "df_Kp_PZ_err",
        "df_BVx_err",
        "df_BVy_err",
        "df_BVz_err",
        "df_Bp_PX_err",
        "df_Bp_PY_err",
        "df_Bp_PZ_err",
        "df_Bp_M2_err",
        "df_taup_PX_err",
        "df_taup_PY_err",
        "df_taup_PZ_err",
        "df_taup_PE_err",
        "df_antinutau_PX_err",
        "df_antinutau_PY_err",
        "df_antinutau_PZ_err",
        "df_antinutau_PE_err",
        "df_taum_PX_err",
        "df_taum_PY_err",
        "df_taum_PZ_err",
        "df_taum_PE_err",
        "df_nutau_PX_err",
        "df_nutau_PY_err",
        "df_nutau_PZ_err",
        "df_nutau_PE_err"
    };

    for(int i = 0; i < dimM+dimX; i++)
    {
        tout->Branch(name_x[i], &X(i));
        tout->Branch(name_x_err[i], &X_ERR(i));
    }

    TString name_F[] = {
        "df_dL_dPVx",
        "df_dL_dPVy",
        "df_dL_dPVz",
        "df_dL_dDV1x",
        "df_dL_dDV1y",
        "df_dL_DV1z",
        "df_dL_dp3pi1x",
        "df_dL_dp3pi1y",
        "df_dL_dp3pi1z",
        "df_dL_dE3pi1",
        "df_dL_dDV2x",
        "df_dL_dDV2y",
        "df_dL_dDV2z",
        "df_dL_dp3pi2x",
        "df_dL_dp3pi2y",
        "df_dL_dp3pi2z",
        "df_dL_dE3pi2",
        "df_dL_dRPx",
        "df_dL_dRPy",
        "df_dL_dpKx",
        "df_dL_dpKy",
        "df_dL_dpKz",
        "df_dL_dBVx",
        "df_dL_dBVy",
        "df_dL_dBVz",
        "df_dL_dpBx",
        "df_dL_dpBy",
        "df_dL_dpBz",
        "df_dL_dMB2",
        "df_dL_dptau1x",
        "df_dL_dptau1y",
        "df_dL_dptau1z",
        "df_dL_dEtau1",
        "df_dL_dpnu1x",
        "df_dL_dpnu1y",
        "df_dL_dpnu1z",
        "df_dL_dEnu1",
        "df_dL_dptau2x",
        "df_dL_dptau2y",
        "df_dL_dptau2z",
        "df_dL_dEtau2",
        "df_dL_dpnu2x",
        "df_dL_dpnu2y",
        "df_dL_dpnu2z",
        "df_dL_dEnu2",
        "df_pB_pointing_to_PV_XZ",
        "df_pB_pointing_to_PV_YZ",
        "df_ptau1_pointing_to_BV_XZ",
        "df_ptau1_pointing_to_BV_YZ",
        "df_px_conservation_in_DV1",
        "df_py_conservation_in_DV1",
        "df_pz_conservation_in_DV1",
        "df_E_conservation_in_DV1",
        "df_taup_mass_constraint",
        "df_antinu_mass_constraint",
        "df_ptau2_pointing_to_BV_XZ",
        "df_ptau2_pointing_to_BV_YZ",
        "df_px_conservation_in_DV2",
        "df_py_conservation_in_DV2",
        "df_pz_conservation_in_DV2",
        "df_E_conservation_in_DV2",
        "df_taum_mass_constraint",
        "df_nu_mass_constraint",
        "df_BV_in_K_trajectory_XZ",
        "df_BV_in_K_trajectory_YZ",
        "df_px_conservation_in_BV",
        "df_py_conservation_in_BV",
        "df_pz_conservation_in_BV",
        "df_E_conservation_in_BV"
    };

    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        tout->Branch(name_F[i], &F(i));
    }

    tout->Branch("df_init", &init);
    tout->Branch("df_Bp_M", &MB);
    tout->Branch("df_Bp_MERR", &MB_err);
    tout->Branch("df_status", &status);
    tout->Branch("df_F_tolerance", &F_tolerance);
    tout->Branch("df_chi2", &chi2);
    tout->Branch("df_nIter", &nIter);
    tout->Branch("df_time", &deltaT);
    // tout->Branch("df_isSaddlePoint", &saddle_point);

    UInt_t num_entries = t->GetEntries();
    
    TStopwatch watch_total;

    // Loop over events
    for(int evt = 0; evt < num_entries; evt++)
    {
        TStopwatch watch;
        t->GetEntry(evt);

        // 1) Fill m and V
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

        // 2) Build weights matrix W = V^-1
        W = V;
        W.Invert();
        // W.Print();

        // TMatrixD identity = W*V;
        // identity.Print();

        // 3) Define the system of equations
        ROOT::Math::GSLMultiRootFinder * solver = new ROOT::Math::GSLMultiRootFinder("HybridS");
        solver->AddFunction(eq1,dimM+dimX+dimC);
        solver->AddFunction(eq2,dimM+dimX+dimC);
        solver->AddFunction(eq3,dimM+dimX+dimC);
        solver->AddFunction(eq4,dimM+dimX+dimC);
        solver->AddFunction(eq5,dimM+dimX+dimC);
        solver->AddFunction(eq6,dimM+dimX+dimC);
        solver->AddFunction(eq7,dimM+dimX+dimC);
        solver->AddFunction(eq8,dimM+dimX+dimC);
        solver->AddFunction(eq9,dimM+dimX+dimC);
        solver->AddFunction(eq10,dimM+dimX+dimC);
        solver->AddFunction(eq11,dimM+dimX+dimC);
        solver->AddFunction(eq12,dimM+dimX+dimC);
        solver->AddFunction(eq13,dimM+dimX+dimC);
        solver->AddFunction(eq14,dimM+dimX+dimC);
        solver->AddFunction(eq15,dimM+dimX+dimC);
        solver->AddFunction(eq16,dimM+dimX+dimC);
        solver->AddFunction(eq17,dimM+dimX+dimC);
        solver->AddFunction(eq18,dimM+dimX+dimC);
        solver->AddFunction(eq19,dimM+dimX+dimC);
        solver->AddFunction(eq20,dimM+dimX+dimC);
        solver->AddFunction(eq21,dimM+dimX+dimC);
        solver->AddFunction(eq22,dimM+dimX+dimC);
        solver->AddFunction(eq23,dimM+dimX+dimC);
        solver->AddFunction(eq24,dimM+dimX+dimC);
        solver->AddFunction(eq25,dimM+dimX+dimC);
        solver->AddFunction(eq26,dimM+dimX+dimC);
        solver->AddFunction(eq27,dimM+dimX+dimC);
        solver->AddFunction(eq28,dimM+dimX+dimC);
        solver->AddFunction(eq29,dimM+dimX+dimC);
        solver->AddFunction(eq30,dimM+dimX+dimC);
        solver->AddFunction(eq31,dimM+dimX+dimC);
        solver->AddFunction(eq32,dimM+dimX+dimC);
        solver->AddFunction(eq33,dimM+dimX+dimC);
        solver->AddFunction(eq34,dimM+dimX+dimC);
        solver->AddFunction(eq35,dimM+dimX+dimC);
        solver->AddFunction(eq36,dimM+dimX+dimC);
        solver->AddFunction(eq37,dimM+dimX+dimC);
        solver->AddFunction(eq38,dimM+dimX+dimC);
        solver->AddFunction(eq39,dimM+dimX+dimC);
        solver->AddFunction(eq40,dimM+dimX+dimC);
        solver->AddFunction(eq41,dimM+dimX+dimC);
        solver->AddFunction(eq42,dimM+dimX+dimC);
        solver->AddFunction(eq43,dimM+dimX+dimC);
        solver->AddFunction(eq44,dimM+dimX+dimC);
        solver->AddFunction(eq45,dimM+dimX+dimC);
        solver->AddFunction(eq46,dimM+dimX+dimC);
        solver->AddFunction(eq47,dimM+dimX+dimC);
        solver->AddFunction(eq48,dimM+dimX+dimC);
        solver->AddFunction(eq49,dimM+dimX+dimC);
        solver->AddFunction(eq50,dimM+dimX+dimC);
        solver->AddFunction(eq51,dimM+dimX+dimC);
        solver->AddFunction(eq52,dimM+dimX+dimC);
        solver->AddFunction(eq53,dimM+dimX+dimC);
        solver->AddFunction(eq54,dimM+dimX+dimC);
        solver->AddFunction(eq55,dimM+dimX+dimC);
        solver->AddFunction(eq56,dimM+dimX+dimC);
        solver->AddFunction(eq57,dimM+dimX+dimC);
        solver->AddFunction(eq58,dimM+dimX+dimC);
        solver->AddFunction(eq59,dimM+dimX+dimC);
        solver->AddFunction(eq60,dimM+dimX+dimC);
        solver->AddFunction(eq61,dimM+dimX+dimC);  
        solver->AddFunction(eq62,dimM+dimX+dimC);  
        solver->AddFunction(eq63,dimM+dimX+dimC);  
        solver->AddFunction(eq64,dimM+dimX+dimC);  
        solver->AddFunction(eq65,dimM+dimX+dimC);  
        solver->AddFunction(eq66,dimM+dimX+dimC);  
        solver->AddFunction(eq67,dimM+dimX+dimC);  
        solver->AddFunction(eq68,dimM+dimX+dimC);  
        solver->AddFunction(eq69,dimM+dimX+dimC);  
        solver->SetPrintLevel(0);

        // 4) Solve the system of equations (69)
        ROOT::Math::XYZPoint BV( BVx, BVy, BVz );

        // Lowest chi^2 approach (with 10 different initialisations)
        // lowest_chi2(BV, solver, initial_tolerance, species);
        // if(status != 0)
        // {
        //     lowest_chi2(BV, solver, pow(10,-6), species);       
        // }
        // if(status != 0)
        // {
        //     lowest_chi2(BV, solver, pow(10,-3), species);
        // }
        // if(status != 0)
        // {
        //     lowest_chi2(BV, solver, 1, species);
        // }
        // if(status != 0)
        // {
        //     lowest_chi2(BV, solver, 100, species);
        // }

        // MLP approach
        // solve( BV, -2, solver, initial_tolerance, species );
        // if(status != 0)
        // {
        //     x_MLP_current = X;
        //     solve( BV, -2, solver, pow(10,-6), species );
        // }
        // if(status != 0)
        // {
        //     x_MLP_current = X;
        //     solve( BV, -2, solver, pow(10,-3), species );
        // }
        // if(status != 0)
        // {
        //     x_MLP_current = X;
        //     solve( BV, -2, solver, 1, species );
        // }
        // if(status != 0)
        // {
        //     x_MLP_current = X;
        //     solve( BV, -2, solver, 100, species );
        // }

        // DTF sequence approach
        sequence(BV, solver, initial_tolerance, species);
        if(status != 0)
        {
            sequence(BV, solver, pow(10,-6), species);       
        }
        if(status != 0)
        {
            sequence(BV, solver, pow(10,-3), species);
        }
        if(status != 0)
        {
            sequence(BV, solver, 1, species);
        }
        if(status != 0)
        {
            sequence(BV, solver, 100, species);
        }

        // 5) Error estimation
        U = U_matrix();
        for(int i = 0; i < dimM+dimX+dimC; i++)
        {
            X_ERR(i) = sqrt(U(i,i));
        }

        Double_t dMB_squared = sqrt(U(28,28));
        MB_err = dMB_squared/(2*abs(MB));    

        cout << "FINAL" << endl;
        cout << "init = " << init << endl;
        cout << "status == " << status << endl;
        cout << "sum_Fi = " << F_tolerance << endl;
        cout << "MB = " << MB << " +/- " << MB_err << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "#iterations = " << nIter << endl;  

        // X.Print();
        // X_ERR.Print();
    
        // TVectorD V_errors(dimM);
        // for(int i = 0; i < dimM; i++)
        // {
        //     V_errors(i) = sqrt(V(i,i));
        // }
        // m.Print();
        // V_errors.Print();

        // for(int i = 0; i < dimM; i++)
        // {
        //     meas_ratio(i) = X_ERR(i)/V_errors(i);
        // }
        // meas_ratio.Print();

        // TMatrixDSym U_corr(dimM+dimX+dimC);
        // for(int i = 0; i < dimM+dimX+dimC; i++)
        // {
        //     for(int j = 0; j < dimM+dimX+dimC; j++)
        //     {
        //         U_corr(i,j) = U(i,j)/sqrt( U(i,i)*U(j,j) );
        //     }
        // }
        // U_corr.Print();

        // TVectorD U_eigen_values(dimM+dimX+dimC);  
        // TMatrixD U_eigen_vectors = U.EigenVectors(U_eigen_values);
        // U_eigen_values.Print();

        /*
        TString NAME[] = {
            "PVx (mm)",
            "PVy (mm)",
            "PVz (mm)",
            "DV1x (mm)",
            "DV1y (mm)",
            "DV1z (mm)",
            "p3pi1x (MeV)",
            "p3pi1y (MeV)",
            "p3pi1z (MeV)",
            "E3pi1 (MeV)",
            "DV2x (mm)",
            "DV2y (mm)",
            "DV2z (mm)",
            "p3pi2x (MeV)",
            "p3pi2y (MeV)",
            "p3pi2z (MeV)",
            "E3pi2 (MeV)",
            "RPx (mm)",
            "RPy (mm)",
            "pKx (MeV)",
            "pKy (MeV)",
            "pKz (MeV)",
            "BVx (mm)",
            "BVy (mm)",
            "BVz (mm)",
            "pBx (MeV)",
            "pBy (MeV)",
            "pBz (MeV)",
            "MB^{2} (MeV^{2})",
            "ptau1x (MeV)",
            "ptau1y (MeV)",
            "ptau1z (MeV)",
            "pnu1x (MeV)",
            "pnu1y (MeV)",
            "pnu1z (MeV)",
            "ptau2x (MeV)",
            "ptau2y (MeV)",
            "ptau2z (MeV)",
            "pnu2x (MeV)",
            "pnu2y (MeV)",
            "pnu2z (MeV)",
            "lambda1",
            "lambda2",
            "lambda3",
            "lambda4",
            "lambda5",
            "lambda6",
            "lambda7",
            "lambda8",
            "lambda9",
            "lambda10",
            "lambda11",
            "lambda12",
            "lambda13",
            "lambda14",
            "lambda15",
            "lambda16",
            "lambda17",
            "lambda18",
            "lambda19",
            "lambda20"
        };
        if( sizeof(NAME)/sizeof(NAME[0]) != dimM+dimX+dimC )
        {
            cout << "NAME size is not correct" << endl;
            return;
        }

        for(int i = 0; i < dimM+dimX+dimC; i++)
        {
            scan_lagrangian(X, i, 100, NAME[i], year, species, line);
        }
        for(int i = 0; i < dimM; i++)
        {
            scan_chisquare(X, i, 100, NAME[i], year, species, line);
        }
        */

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

void solve( ROOT::Math::XYZPoint BV, int init,  ROOT::Math::GSLMultiRootFinder* solver, Double_t tolerance, Int_t species )
{
    // 1) This function initialises the vector of unkown parameters x with the result of analytical calculations
    // 2) It sets up the minimizer: defining the initial values and bounds on the x parameters  
    // 3) It builds the chi^2 function with x0 by calling the function chisquare 
    // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

    // 1) Initial values for x=(xm,xu,lambda)
    if(init == 0)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_0( m ); // original calculations = Anne Keune
        }
        else
        {       
            x0 = x0_current; // where it left 
        }
    }
    else if(init == 1)
    {   
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_1( m, BV ); // K*tautau calculations; tau momentum direction initialised based on vertices (uses offline estimate for BV)
        }
        else
        {
            x0 = x1_current; // where it left 
        }
    }
    else if(init == 2)
    {   
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_2( m, BV ); // K*tautau calculations; tau momentum direction initialised based on 3pi visible momentum (uses offline estimate for BV) 
        }
        else
        {
            x0 = x2_current; // where it left 
        }
    }
    else if(init == 3)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_3( m, BV ); // RD calculations (uses offline estimate for BV)
        }
        else
        {
            x0 = x3_current; // where it left 
        }
    }
    else if(init == 4)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_4( m, BV ); // K*tautau; tau+ from vertices, tau- from pions
        }
        else
        {
            x0 = x4_current; // where it left 
        }
    }
    else if(init == 5)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_5( m, BV ); // K*tautau; tau- from vertices, tau+ from pions
        }
        else
        {
            x0 = x5_current; // where it left
        }
    }
    else if(init == 6)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_6( m, BV ); // Mixes Marseille (tau+) and K*tautau vertices (tau-)
        }
        else
        {
            x0 = x6_current; // where it left 
        }
    }
    else if(init == 7)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_7( m, BV ); // Mixes Marseille (tau-) and K*tautau vertices (tau+)
        }
        else
        {
            x0 = x7_current; // where it left 
        }
    }
    else if(init == 8)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_8( m, BV ); // Mixes Marseille (tau+) and K*tautau pions (tau-)
        }
        else
        {
            x0 = x8_current; // where it left 
        }
    }
    else if(init == 9)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_9( m, BV ); // Mixes Marseille (tau-) and K*tautau pions (tau+)
        }
        else
        {
            x0 = x9_current; // where it left 
        }
    }
    else if(init == -1)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_true( m );
        }
        else
        {
            x0 = x_true_current;
        }
        
    }
    else if(init == -2)
    {
        if(tolerance == initial_tolerance)
        {
            x0 = x_initial_estimate_MLP( m, BV, species );
        }
        else
        {
            x0 = x_MLP_current;
        }
    }
    // x0.Print();

    Double_t x0_vars[dimM+dimX+dimC];
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        x0_vars[i] = x0(i); 
    }

    // 2) Solve system of equations
    solver->Solve(x0_vars, 10000, tolerance);

    // 3) Retrieve results
    const double* x_results = solver->X();
    const double* f_vals = solver->FVal();
    status = solver->Status();
    if(tolerance == initial_tolerance)
    { 
        nIter = solver->Iterations();
    }
    else
    {
        nIter += solver->Iterations();
    }

    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        X(i) = x_results[i];
        F(i) = f_vals[i];
    }
    // X.Print();
    // F.Print();

    Double_t MB_squared = X(dimM+6);
    if( MB_squared > 0 )
    {
        MB = sqrt( MB_squared );
    }
    else
    {
        MB = -sqrt( abs(MB_squared) );
    }

    F_tolerance = 0;
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        F_tolerance += abs(F(i));
    }

    TVectorD xm(dimM);
    for(int i = 0; i < dimM; i++)
    {
        xm(i) = X(i);
    }
    chi2 = chisquare(xm);

}

void lowest_chi2(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder *solver, Double_t tolerance, Int_t species)
{
        // Lowest sum
        // Init 0
        solve( BV, 0, solver, tolerance, species );
        x0_current = X;
        Int_t status0 = status;
        Double_t MB_0 = MB;
        Double_t MB_err0 = MB_err;
        TVectorD X0 = X;
        TVectorD X0_ERR = X_ERR;
        TVectorD F0 = F;
        Double_t F_tol_0 = F_tolerance;
        Double_t chi2_0 = chi2;
        Int_t nIter_0 = nIter;

        // cout << "MB_0 = " << MB_0 << endl;
        // cout << "sum_0 = " << F_tol_0 << endl;
        // cout << "chi2_0 = " << chi2_0 << endl;
        // cout << "status_0 = " << status0 << endl;
        // cout << "nIter_0 = " << nIter_0 << endl;

        // Init 1
        solve( BV, 1, solver, tolerance, species );
        x1_current = X;
        Int_t status1 = status;
        Double_t MB_1 = MB;
        Double_t MB_err1 = MB_err;
        TVectorD X1 = X;
        TVectorD X1_ERR = X_ERR;
        TVectorD F1 = F;
        Double_t F_tol_1 = F_tolerance;
        Double_t chi2_1 = chi2;
        Int_t nIter_1 = nIter;
    
        // cout << "MB_1 = " << MB_1 << endl;
        // cout << "sum_1 = " << F_tol_1 << endl;
        // cout << "chi2_1 = " << chi2_1 << endl;
        // cout << "status_1 = " << status1 << endl;
        // cout << "nIter_1 = " << nIter_1 << endl;

        // Init 2
        solve( BV, 2, solver, tolerance, species );
        x2_current = X;
        Int_t status2 = status;
        Double_t MB_2 = MB;
        Double_t MB_err2 = MB_err;
        TVectorD X2 = X;
        TVectorD X2_ERR = X_ERR;
        TVectorD F2 = F;
        Double_t F_tol_2 = F_tolerance;
        Double_t chi2_2 = chi2;
        Int_t nIter_2 = nIter;

        // cout << "MB_2 = " << MB_2 << endl;
        // cout << "sum_2 = " << F_tol_2 << endl;
        // cout << "chi2_2 = " << chi2_2 << endl;
        // cout << "status_2 = " << status2 << endl;
        // cout << "nIter_2 = " << nIter_2 << endl;

        // Init 3
        solve( BV, 3, solver, tolerance, species );
        x3_current = X;
        Int_t status3 = status;
        Double_t MB_3 = MB;
        Double_t MB_err3 = MB_err;
        TVectorD X3 = X;
        TVectorD X3_ERR = X_ERR;
        TVectorD F3 = F; 
        Double_t F_tol_3 = F_tolerance; 
        Double_t chi2_3 = chi2; 
        Int_t nIter_3 = nIter;

        // cout << "MB_3 = " << MB_3 << endl;
        // cout << "sum_3 = " << F_tol_3 << endl;
        // cout << "chi2_3 = " << chi2_3 << endl;
        // cout << "status_3 = " << status3 << endl;
        // cout << "nIter_3 = " << nIter_3 << endl;

        // Init 4
        solve( BV, 4, solver, tolerance, species );
        x4_current = X;
        Int_t status4 = status;
        Double_t MB_4 = MB;
        Double_t MB_err4 = MB_err;
        TVectorD X4 = X;
        TVectorD X4_ERR = X_ERR;
        TVectorD F4 = F;
        Double_t F_tol_4 = F_tolerance;
        Double_t chi2_4 = chi2;
        Int_t nIter_4 = nIter;

        // cout << "MB_4 = " << MB_4 << endl;
        // cout << "sum_4 = " << F_tol_4 << endl;
        // cout << "chi2_4 = " << chi2_4 << endl;
        // cout << "status_4 = " << status4 << endl;
        // cout << "nIter_4 = " << nIter_4 << endl;

        // Init 5
        solve( BV, 5, solver, tolerance, species );
        x5_current = X;
        Int_t status5 = status;
        Double_t MB_5 = MB;
        Double_t MB_err5 = MB_err;
        TVectorD X5 = X;
        TVectorD X5_ERR = X_ERR;
        TVectorD F5 = F;
        Double_t F_tol_5 = F_tolerance;
        Double_t chi2_5 = chi2;
        Int_t nIter_5 = nIter;

        // cout << "MB_5 = " << MB_5 << endl;
        // cout << "sum_5 = " << F_tol_5 << endl;
        // cout << "chi2_5 = " << chi2_5 << endl;
        // cout << "status_5 = " << status5 << endl;
        // cout << "nIter_5 = " << nIter_5 << endl;

        // Init 6
        solve( BV, 6, solver, tolerance, species );
        x6_current = X;
        Int_t status6 = status;
        Double_t MB_6 = MB;
        Double_t MB_err6 = MB_err;
        TVectorD X6 = X;
        TVectorD X6_ERR = X_ERR;
        TVectorD F6 = F;
        Double_t F_tol_6 = F_tolerance;
        Double_t chi2_6 = chi2;
        Int_t nIter_6 = nIter;

        // cout << "MB_6 = " << MB_6 << endl;
        // cout << "sum_6 = " << F_tol_6 << endl;
        // cout << "chi2_6 = " << chi2_6 << endl;
        // cout << "status_6 = " << status6 << endl;
        // cout << "nIter_6 = " << nIter_6 << endl;

        // Init 7
        solve( BV, 7, solver, tolerance, species );
        x7_current = X;
        Int_t status7 = status;
        Double_t MB_7 = MB;
        Double_t MB_err7 = MB_err;
        TVectorD X7 = X;
        TVectorD X7_ERR = X_ERR;
        TVectorD F7 = F;
        Double_t F_tol_7 = F_tolerance;
        Double_t chi2_7 = chi2;
        Int_t nIter_7 = nIter;

        // cout << "MB_7 = " << MB_7 << endl;
        // cout << "sum_7 = " << F_tol_7 << endl;
        // cout << "chi2_7 = " << chi2_7 << endl;
        // cout << "status_7 = " << status7 << endl;
        // cout << "nIter_7 = " << nIter_7 << endl;

        // Init 8
        solve( BV, 8, solver, tolerance, species );
        x8_current = X;
        Int_t status8 = status;
        Double_t MB_8 = MB;
        Double_t MB_err8 = MB_err;
        TVectorD X8 = X;
        TVectorD X8_ERR = X_ERR;
        TVectorD F8 = F;
        Double_t F_tol_8 = F_tolerance;
        Double_t chi2_8 = chi2;
        Int_t nIter_8 = nIter;

        // cout << "MB_8 = " << MB_8 << endl;
        // cout << "sum_8 = " << F_tol_8 << endl;
        // cout << "chi2_8 = " << chi2_8 << endl;
        // cout << "status_8 = " << status8 << endl;
        // cout << "nIter_8 = " << nIter_8 << endl;

        // Init 9
        solve( BV, 9, solver, tolerance, species );
        x9_current = X;
        Int_t status9 = status;
        Double_t MB_9 = MB;
        Double_t MB_err9 = MB_err;
        TVectorD X9 = X;
        TVectorD X9_ERR = X_ERR;
        TVectorD F9 = F;
        Double_t F_tol_9 = F_tolerance;
        Double_t chi2_9 = chi2;
        Int_t nIter_9 = nIter;

        // cout << "MB_9 = " << MB_9 << endl;
        // cout << "sum_9 = " << F_tol_9 << endl;
        // cout << "chi2_9 = " << chi2_9 << endl;
        // cout << "status_9 = " << status9 << endl;
        // cout << "nIter_9 = " << nIter_9 << endl;

        Int_t status_vec[10] = { status0, status1, status2, status3, status4, status5, status6, status7, status8, status9 };
        Double_t MB_vec[10] = { MB_0, MB_1, MB_2, MB_3, MB_4, MB_5, MB_6, MB_7, MB_8, MB_9 }; 
        Double_t MB_err_vec[10] = { MB_err0, MB_err1, MB_err2, MB_err3, MB_err4, MB_err5, MB_err6, MB_err7, MB_err8, MB_err9 };
        Double_t F_tol_vec[10] = { F_tol_0, F_tol_1, F_tol_2, F_tol_3, F_tol_4, F_tol_5, F_tol_6, F_tol_7, F_tol_8, F_tol_9 };
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

        std::vector<TVectorD> F_vec;
        F_vec.push_back(F0);
        F_vec.push_back(F1);
        F_vec.push_back(F2);
        F_vec.push_back(F3);
        F_vec.push_back(F4);
        F_vec.push_back(F5);
        F_vec.push_back(F6);
        F_vec.push_back(F7);
        F_vec.push_back(F8);
        F_vec.push_back(F9);

/*         Double_t F_min = 100;
        Int_t i_min = 0;
        for(int i = 0; i < 10; i++)
        {
            if( (F_tol_vec[i] < F_min) )
            {
                F_min = F_tol_vec[i];
                i_min = i;
            }
        }   */ 

        bool all_fail = false;
        if( (status0 != 0) && (status1 != 0) && (status2 != 0) && (status3 != 0) && (status4 != 0) && (status5 != 0) && (status6 != 0) && (status7 != 0) && (status8 != 0) && (status9 != 0) )
        {
            all_fail = true;
        }

        
        Double_t chi2_min = 100000000000000000;
        Double_t F_min = 100;
        Int_t i_min = 0;

        if(all_fail) // If all fail, returns the one with the lowest value for the sum
        {
            for(int i = 0; i < 10; i++)
            {
                if( (F_tol_vec[i] < F_min) )
                {
                    F_min = F_tol_vec[i];
                    i_min = i;
                }
            }   
        }
        else // if one or more passes, from the ones that pass return the one that gives the lowest value for the chi^2
        {
            for(int i = 0; i < 10; i++)
            {
                if( (status_vec[i] == 0) && (chi2_vec[i] < chi2_min) )
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
        X = X_vec[i_min];
        X_ERR = X_ERR_vec[i_min];
        F = F_vec[i_min];
        F_tolerance = F_tol_vec[i_min];
        nIter = nIter_vec[i_min];
        chi2 = chi2_vec[i_min];

}

void sequence(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder *solver, Double_t tolerance, Int_t species)
{
    init = 1;  // K*tautau vertex
    solve( BV, init, solver, tolerance, species );
    x1_current = X;
    if(status != 0)
    {
        init = 2; // K*tautau pions
        solve( BV, init, solver, tolerance, species );
        x2_current = X;
    }
    if(status != 0)
    {
        init = 3; // Marseille
        solve( BV, init, solver, tolerance, species );
        x3_current = X;
    }

    // init = 3; // Marseille
    // solve( BV, init, solver, tolerance, species );
    // if(status != 0)
    // {
    //     init = 0; // Original
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 9; // Marseille (tau-) + K*tautau pions (tau+)
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 8; // Marseille (tau+) + K*tautau pions (tau-)
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 7; // Marseille (tau-) + K*tautau vertex (tau+)
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 6; // Marseille (tau+) + K*tautau vertex (tau-)
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 2; // K*tautau pions
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 5; // K*tautau vertex (tau-) + K*tautau pions (tau+)
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 4; // K*tautau vertex (tau+) + K*tautau pions (tau-)
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = 1; // K*tautau vertex
    //     solve( BV, init, solver, tolerance, species );
    // }
    // if(status != 0)
    // {
    //     init = -1;
    // }
}

void scan_chisquare( TVectorD X, Int_t index, Int_t npoints, TString x_name, Int_t year, Int_t species, Int_t line)
{
    if(index > dimM)
    {
        cout << "Wrong index" << endl;
        return;
    }

    TVectorD xm(dimM);
    for(int i = 0; i < dimM; i++)
    {
        xm(i) = X(i);
    }

    Double_t x_val = xm(index);
    Double_t x_min = x_val - X_ERR(index);
    Double_t x_max = x_val + X_ERR(index);

    vector<Double_t> scan_range = range(x_min, x_max, npoints);

    Double_t x[npoints];
    Double_t y[npoints];

    for(int i = 0; i < npoints; i++)
    {
        x[i] = scan_range[i];
        xm(index) = x[i];
        y[i] = chisquare(xm);
    }

    TGraph* g = new TGraph(npoints, x, y);

    TCanvas c;
    c.cd();
    gPad->SetGrid(1,1);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kBlue);
    g->Draw("AP");
    g->GetXaxis()->SetTitle(x_name);
    g->GetYaxis()->SetTitle("#chi^{2}");
    g->SetTitle("");

    TLine* line1 = new TLine(x_val, g->GetYaxis()->GetXmin(), x_val, g->GetYaxis()->GetXmax());
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->Draw("same");

    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/Plots_chisquare/%i_%i.gif",year,species,index,line));
    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/Plots_chisquare/%i_%i.pdf",year,species,index,line));

    X(index) = x_val;
}

void scan_lagrangian( TVectorD X, Int_t index, Int_t npoints, TString x_name, Int_t year, Int_t species, Int_t line)
{
    Double_t x_val = X(index);
    Double_t x_min = x_val - X_ERR(index);
    Double_t x_max = x_val + X_ERR(index);

    vector<Double_t> scan_range = range(x_min, x_max, npoints);

    Double_t x[npoints];
    Double_t y[npoints];

    for(int i = 0; i < npoints; i++)
    {
        x[i] = scan_range[i];
        X(index) = x[i];
        y[i] = log(lagrangian(X));
    }

    TGraph* g = new TGraph(npoints, x, y);

    TCanvas c;
    c.cd();
    gPad->SetGrid(1,1);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kBlue);
    g->Draw("AP");
    g->GetXaxis()->SetTitle(x_name);
    g->GetYaxis()->SetTitle("ln(L)");
    g->SetTitle("");

    TLine* line2 = new TLine(x_val, g->GetYaxis()->GetXmin(), x_val, g->GetYaxis()->GetXmax());
    line2->SetLineColor(kRed);
    line2->SetLineWidth(2);
    line2->Draw("same");

    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/Plots_lagrangian/%i_%i.gif",year,species,index,line));
    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/Plots_lagrangian/%i_%i.pdf",year,species,index,line));

    X(index) = x_val;
}

vector<double> range(double min, double max, size_t N) {
    vector<double> range;
    double delta = (max-min)/double(N-1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}

TMatrixDSym U_matrix()
{
    TMatrixDSym U_sym(dimM+dimX+dimC);

    // Matrix B
    TMatrixD B(dimM+dimX+dimC,dimM);
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        for(int j = 0; j < dimM; j++)
        {
            if(i < dimM)
            {
                B(i,j) = -2*W(i,j);
            }
            else
            {
                B(i,j) = 0.;
            }
        }
    }
    // B.Print();

    TMatrixD A(dimM+dimX+dimC,dimM+dimX+dimC);

    // Measured parameters (22)
    ROOT::Math::XYZPoint PV( X(0), X(1), X(2) );
    ROOT::Math::XYZPoint DV1( X(3), X(4), X(5) );
    ROOT::Math::XYZVector p3pi1( X(6), X(7), X(8) );
    Double_t E3pi1 = X(9);
    ROOT::Math::XYZPoint DV2( X(10), X(11), X(12) );
    ROOT::Math::XYZVector p3pi2( X(13), X(14), X(15) );
    Double_t E3pi2 = X(16);
    ROOT::Math::XYZPoint RP( X(17), X(18), RPz );
    ROOT::Math::XYZVector pK( X(19), X(20), X(21) );
    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

    // Unknown parameters (19)
    ROOT::Math::XYZPoint BV( X(dimM), X(dimM+1), X(dimM+2) );
    ROOT::Math::XYZVector pB( X(dimM+3), X(dimM+4), X(dimM+5) );
    Double_t MB_squared = X(dimM+6);
    ROOT::Math::XYZVector ptau1( X(dimM+7), X(dimM+8), X(dimM+9) );
    Double_t Etau1 = X(dimM+10);
    ROOT::Math::XYZVector pnu1( X(dimM+11), X(dimM+12), X(dimM+13) );
    Double_t Enu1 = X(dimM+14);
    ROOT::Math::XYZVector ptau2( X(dimM+15), X(dimM+16), X(dimM+17) );
    Double_t Etau2 = X(dimM+18);
    ROOT::Math::XYZVector pnu2( X(dimM+19), X(dimM+20), X(dimM+21) );
    Double_t Enu2 = X(dimM+22);
    Double_t EB = sqrt( MB_squared + pB.Mag2() );

    // Lagrange multipliers (20)
    TVectorD l(dimC);
    for(int i = 0; i < dimC; i++)
    {
        l(i) = X(dimM+dimX+i);
    }
    int L = dimM+dimX;

    ///////////////////////////// xm ///////////////////////////
    // F0 (PVx)
    for(int i = 0; i < dimM; i++)
    {
        A(0,i) = 2*W(0,i);
    }
    A(0,dimM+5) = l(0);
    A(0,L) = pB.z();

    // F1 (PVy)
    for(int i = 0; i < dimM; i++)
    {
        A(1,i) = 2*W(1,i);
    }
    A(1,dimM+5) = l(1);
    A(1,L+1) = pB.z();

    // F2 (PVz)
    for(int i = 0; i < dimM; i++)
    {
        A(2,i) = 2*W(2,i);
    }
    A(2,dimM+3) = -l(0);
    A(2,dimM+4) = -l(1);
    A(2,L) = -pB.x();
    A(2,L+1) = -pB.y();

    // F3 (DV1x)
    for(int i = 0; i < dimM; i++)
    {
        A(3,i) = 2*W(3,i);
    }
    A(3,dimM+9) = -l(2);
    A(3,L+2) = -ptau1.z();

    // F4 (DV1y)
    for(int i = 0; i < dimM; i++)
    {
        A(4,i) = 2*W(4,i);
    }
    A(4,dimM+9) = -l(3);
    A(4,L+3) = -ptau1.z();

    // F5 (DV1z)
    for(int i = 0; i < dimM; i++)
    {
            A(5,i) = 2*W(5,i);
    }
    A(5,dimM+7) = l(2);
    A(5,dimM+8) = l(3);
    A(5,L+2) = ptau1.x();
    A(5,L+3) = ptau1.y();

    // F6,7,8,9 (P3pi1)
    for(int i = 0; i < dimM; i++)
    {
        A(6,i) = 2*W(6,i);
        A(7,i) = 2*W(7,i);
        A(8,i) = 2*W(8,i);
        A(9,i) = 2*W(9,i);
    }
    A(6,L+4) = -1;
    A(7,L+5) = -1;
    A(8,L+6) = -1;
    A(9,L+7) = -1;

    // F10 (DV2x)
    for(int i = 0; i < dimM; i++)
    {

        A(10,i) = 2*W(10,i);
    }
    A(10,dimM+17) = -l(10);
    A(10,L+10) = -ptau2.z();

    // F11 (DV2y)
    for(int i = 0; i < dimM; i++)
    {
        A(11,i) = 2*W(11,i);
    }
    A(11,dimM+17) = -l(11);
    A(11,L+11) = -ptau2.z();

    // F12 (DV2z)
    for(int i = 0; i < dimM; i++)
    {
        A(12,i) = 2*W(12,i);
    }
    A(12,dimM+15) = l(10);
    A(12,dimM+16) = l(11);
    A(12,L+10) = ptau2.x();
    A(12,L+11) = ptau2.y();

    // F13,14,15,16 (P3pi2)
    for(int i = 0; i < dimM; i++)
    {
        A(13,i) = 2*W(13,i);
        A(14,i) = 2*W(14,i);
        A(15,i) = 2*W(15,i);
        A(16,i) = 2*W(16,i);
    }
    A(13,L+12) = -1;
    A(14,L+13) = -1;
    A(15,L+14) = -1;
    A(16,L+15) = -1;

    // F17 (RPx)
    A(17,21) = 2*W(17,21) + l(18);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 21)
        {
            A(17,i) = 2*W(17,i);
        }
    }
    A(17,L+18) = pK.z();

    // F18 (RPy)
    A(18,21) = 2*W(18,21) + l(19);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 21)
        {
            A(18,i) = 2*W(18,i);
        }
    }
    A(18,L+19) = pK.z();

    // F19 (pKx)
    A(19,19) = 2*W(19,19) - ((pow(EK,2) - pow(pK.x(),2))/pow(EK,3))*l(23);
    A(19,20) = 2*W(19,20) + ((pK.x()*pK.y()/pow(EK,3)))*l(23);
    A(19,21) = 2*W(19,21) + ((pK.x()*pK.z()/pow(EK,3)))*l(23);
    for(int i = 0; i < dimM; i++)
    {
        if((i!=19) && (i!=20) && (i!=21))
        {
            A(19,i) = 2*W(19,i);
        }
    }
    A(19,dimM+2) = l(18);
    A(19,L+18) = BV.z()-RP.z();
    A(19,L+20) = -1;
    A(19,L+23) = -pK.x()/EK;

    // F20 (pKy)
    A(20,19) = 2*W(20,19) + ((pK.x()*pK.y())/pow(EK,3))*l(23);
    A(20,20) = 2*W(20,20) - ((pow(EK,2) - pow(pK.y(),2))/pow(EK,3))*l(23);
    A(20,21) = 2*W(20,21) + ((pK.y()*pK.z())/pow(EK,3))*l(23);
    for(int i = 0; i < dimM; i++)
    {
        if((i!=19) && (i!=20) && (i!=21))
        {
            A(20,i) = 2*W(20,i);
        }
    }
    A(20,dimM+2) = l(19);
    A(20,L+19) = BV.z()-RP.z();
    A(20,L+21) = -1;
    A(20,L+23) = -pK.y()/EK;

    // F21 (pKz)
    A(21,17) = 2*W(21,17) + l(18);
    A(21,18) = 2*W(21,18) + l(19);
    A(21,19) = 2*W(21,19) + ((pK.x()*pK.z())/pow(EK,3))*l(23);
    A(21,20) = 2*W(21,20) + ((pK.y()*pK.z())/pow(EK,3))*l(23);
    A(21,21) = 2*W(21,21) - ((pow(EK,2) - pow(pK.z(),2))/pow(EK,3))*l(23);
    for(int i = 0; i < dimM; i++)
    {
        if((i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(21,i) = 2*W(21,i);
        }
    }
    A(21,dimM) = -l(18);
    A(21,dimM+1) = -l(19);
    A(21,L+18) = -(BV.x()-RP.x());
    A(21,L+19) = -(BV.y()-RP.y());
    A(21,L+22) = -1;
    A(21,L+23) = -pK.z()/EK;

    //////////////////////////////// xu /////////////////////////
    // F22 (BVx)
    A(22,21) = -l(18);
    A(22,dimM+5) = -l(0);
    A(22,dimM+9) = l(2);
    A(22,dimM+17) = l(10);
    A(22,L) = -pB.z();
    A(22,L+2) = ptau1.z();
    A(22,L+10) = ptau2.z();
    A(22,L+18) = -pK.z();

    // F23 (BVy)
    A(23,21) = -l(19);
    A(23,dimM+5) = -l(1);
    A(23,dimM+9) = l(3);
    A(23,dimM+17) = l(11);
    A(23,L+1) = -pB.z();
    A(23,L+3) = ptau1.z();
    A(23,L+11) = ptau2.z();
    A(23,L+19) = -pK.z();

    // F24 (BVz)
    A(24,19) = l(18);
    A(24,20) = l(19);
    A(24,dimM+3) = l(0);
    A(24,dimM+4) = l(1);
    A(24,dimM+7) = -l(2);
    A(24,dimM+8) = -l(3);
    A(24,dimM+15) = -l(10);
    A(24,dimM+16) = -l(11);
    A(24,L) = pB.x();
    A(24,L+1) = pB.y();
    A(24,L+2) = -ptau1.x();
    A(24,L+3) = -ptau1.y();
    A(24,L+10) = -ptau2.x();
    A(24,L+11) = -ptau2.y();
    A(24,L+18) = pK.x();
    A(24,L+19) = pK.y();

    // F25 (pBx)
    A(25,2) = -l(0);
    A(25,dimM+2) = l(0);
    A(25,dimM+3) = ((pow(EB,2) - pow(pB.x(),2))/pow(EB,3))*l(23);
    A(25,dimM+4) = -((pB.x()*pB.y())/pow(EB,3))*l(23);
    A(25,dimM+5) = -((pB.x()*pB.z())/pow(EB,3))*l(23);
    A(25,dimM+6) = -(pB.x()/(2*pow(EB,3)))*l(23);
    A(25,L) = BV.z()-PV.z();
    A(25,L+20) = 1;
    A(25,L+23) = pB.x()/EB;

    // F26 (pBy)
    A(26,2) = -l(1);
    A(26,dimM+2) = l(1);
    A(26,dimM+3) = -((pB.x()*pB.y())/pow(EB,3))*l(23);
    A(26,dimM+4) = ((pow(EB,2) - pow(pB.y(),2))/pow(EB,3))*l(23);
    A(26,dimM+5) = -((pB.y()*pB.z())/pow(EB,3))*l(23);
    A(26,dimM+6) = -(pB.y()/(2*pow(EB,3)))*l(23);
    A(26,L+1) = BV.z()-PV.z();
    A(26,L+21) = 1;
    A(26,L+23) = pB.y()/EB;

    // F27 (pBz)
    A(27,0) = l(0);
    A(27,1) = l(1);
    A(27,dimM) = -l(0);
    A(27,dimM+1) = -l(1);
    A(27,dimM+3) = -((pB.x()*pB.z())/pow(EB,3))*l(23);
    A(27,dimM+4) = -((pB.y()*pB.z())/pow(EB,3))*l(23);
    A(27,dimM+5) = ((pow(EB,2) - pow(pB.z(),2))/pow(EB,3))*l(23);
    A(27,dimM+6) = -(pB.z()/(2*pow(EB,3)))*l(23);
    A(27,L) = -(BV.x()-PV.x());
    A(27,L+1) = -(BV.y()-PV.y());
    A(27,L+22) = 1;
    A(27,L+23) = pB.z()/EB;

    // F28 (MB_squared)
    A(28,dimM+3) = -(pB.x()/(2*pow(EB,3)))*l(23);
    A(28,dimM+4) = -(pB.y()/(2*pow(EB,3)))*l(23);
    A(28,dimM+5) = -(pB.z()/(2*pow(EB,3)))*l(23);
    A(28,dimM+6) = -(1./(4*pow(EB,3)))*l(23);
    A(28,L+23) = 1./(2*EB);

    // F29 (ptau1x)
    A(29,5) = l(2);
    A(29,dimM+2) = -l(2);
    A(29,dimM+7) = -(( pow(mtau,2) + pow(ptau1.y(),2) + pow(ptau1.z(),2) )/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(29,dimM+8) = ((ptau1.x()*ptau1.y())/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(29,dimM+9) = ((ptau1.x()*ptau1.z())/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(29,L+2) = DV1.z()-BV.z();
    A(29,L+4) = 1;
    A(29,L+8) = - ptau1.x()/sqrt(pow(mtau,2) + ptau1.Mag2());
    A(29,L+20) = -1;

    // F30 (ptau1y)
    A(30,5) = l(3);
    A(30,dimM+2) = -l(3);
    A(30,dimM+7) = ((ptau1.x()*ptau1.y())/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(30,dimM+8) = -(( pow(mtau,2) + pow(ptau1.x(),2) + pow(ptau1.z(),2) )/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(30,dimM+9) = ((ptau1.y()*ptau1.z())/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(30,L+3) = DV1.z()-BV.z();
    A(30,L+5) = 1;
    A(30,L+8) = -ptau1.y()/sqrt(pow(mtau,2) + ptau1.Mag2());
    A(30,L+21) = -1;

    // F31 (ptau1z)
    A(31,3) = -l(2);
    A(31,4) = -l(3);
    A(31,dimM) = l(2);
    A(31,dimM+1) = l(3);
    A(31,dimM+7) = ((ptau1.x()*ptau1.z())/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(31,dimM+8) = ((ptau1.y()*ptau1.z())/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(31,dimM+9) = -(( pow(mtau,2) + pow(ptau1.x(),2) + pow(ptau1.y(),2) )/pow( pow(mtau,2) + ptau1.Mag2(), 1.5))*l(8);
    A(31,L+2) = -(DV1.x()-BV.x());
    A(31,L+3) = -(DV1.y()-BV.y());
    A(31,L+6) = 1;
    A(31,L+8) = -ptau1.z()/sqrt(pow(mtau,2) + ptau1.Mag2());
    A(31,L+22) = -1;

    // F32 (Etau1)
    A(32,L+7) = 1.;
    A(32,L+8) = 1.;
    A(32,L+23) = -1;

    // F33 (pnu1x)
    A(33,dimM+11) = -(( pow(pnu1.y(),2) + pow(pnu1.z(),2) )/pow( pnu1.Mag2(), 1.5))*l(9);
    A(33,dimM+12) = ((pnu1.x()*pnu1.y())/pow( pnu1.Mag2(), 1.5))*l(9);
    A(33,dimM+13) = ((pnu1.x()*pnu1.z())/pow( pnu1.Mag2(), 1.5))*l(9);
    A(33,L+4) = -1;
    A(33,L+9) = -pnu1.x()/sqrt(pnu1.Mag2());

    // F34 (pnu1y)
    A(34,dimM+11) = ((pnu1.x()*pnu1.y())/pow( pnu1.Mag2(), 1.5))*l(9);
    A(34,dimM+12) = -(( pow(pnu1.x(),2) + pow(pnu1.z(),2) )/pow( pnu1.Mag2(), 1.5))*l(9);
    A(34,dimM+13) = ((pnu1.y()*pnu1.z())/pow( pnu1.Mag2(), 1.5))*l(9);
    A(34,L+5) = -1;
    A(34,L+9) = -pnu1.y()/sqrt(pnu1.Mag2());

    // F35 (pnu1z)
    A(35,dimM+11) = ((pnu1.x()*pnu1.z())/pow( pnu1.Mag2(), 1.5))*l(9);
    A(35,dimM+12) = ((pnu1.y()*pnu1.z())/pow( pnu1.Mag2(), 1.5))*l(9);
    A(35,dimM+13) = -((pow(pnu1.x(),2) + pow(pnu1.y(),2))/pow( pnu1.Mag2(), 1.5))*l(9);
    A(35,L+6) = -1;
    A(35,L+9) = -pnu1.z()/sqrt(pnu1.Mag2());

    // F36 (Enu1)
    A(36,L+7) = -1;
    A(36,L+9) = 1;

    // F37 (ptau2x)
    A(37,12) = l(10);
    A(37,dimM+2) = -l(10);
    A(37,dimM+15) = -((pow(mtau,2) + pow(ptau2.y(),2) + pow(ptau2.z(),2))/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(37,dimM+16) = ((ptau2.x()*ptau2.y())/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(37,dimM+17) = ((ptau2.x()*ptau2.z())/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(37,L+10) = DV2.z()-BV.z();
    A(37,L+12) = 1;
    A(37,L+16) = -ptau2.x()/sqrt( pow(mtau,2) + ptau2.Mag2() );
    A(37,L+20) = -1;

    // F38 (ptau2y)
    A(38,12) = l(11);
    A(38,dimM+2) = -l(11);
    A(38,dimM+15) = ((ptau2.x()*ptau2.y())/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(38,dimM+16) = -(( pow(mtau,2) + pow(ptau2.x(),2) + pow(ptau2.z(),2) )/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(38,dimM+17) = ((ptau2.y()*ptau2.z())/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(38,L+11) = DV2.z()-BV.z();
    A(38,L+13) = 1;
    A(38,L+16) = -ptau2.y()/sqrt( pow(mtau,2) + ptau2.Mag2() );
    A(38,L+21) = -1;

    // F39 (ptau2z)
    A(39,10) = -l(10);
    A(39,11) = -l(11);
    A(39,dimM) = l(10);
    A(39,dimM+1) = l(11);
    A(39,dimM+15) = ((ptau2.x()*ptau2.z())/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(39,dimM+16) = ((ptau2.y()*ptau2.z())/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(39,dimM+17) = -((pow(mtau,2) + pow(ptau2.x(),2) + pow(ptau2.y(),2))/pow( pow(mtau,2) + ptau2.Mag2(), 1.5))*l(16);
    A(39,L+10) = -(DV2.x()-BV.x());
    A(39,L+11) = -(DV2.y()-BV.y());
    A(39,L+14) = 1;
    A(39,L+16) = -ptau2.z()/sqrt( pow(mtau,2) + ptau2.Mag2() );
    A(39,L+22) = -1;

    // F40 (Etau2)
    A(40,L+15) = 1;
    A(40,L+16) = 1;
    A(40,L+23) = -1;

    // F41 (pnu2x)
    A(41,dimM+19) = -((pow(pnu2.y(),2) + pow(pnu2.z(),2))/pow( pnu2.Mag2(), 1.5))*l(17);
    A(41,dimM+20) = ((pnu2.x()*pnu2.y())/pow( pnu2.Mag2(), 1.5))*l(17);
    A(41,dimM+21) = ((pnu2.x()*pnu2.z())/pow( pnu2.Mag2(), 1.5))*l(17);
    A(41,L+12) = -1;
    A(41,L+17) = -pnu2.x()/sqrt( pnu2.Mag2() );

    // F42 (pnu2y)
    A(42,dimM+19) = ((pnu2.x()*pnu2.y())/pow( pnu2.Mag2(), 1.5))*l(17);
    A(42,dimM+20) = -(( pow(pnu2.x(),2) + pow(pnu2.z(),2) )/pow( pnu2.Mag2(), 1.5))*l(17);
    A(42,dimM+21) = ((pnu2.y()*pnu2.z())/pow( pnu2.Mag2(), 1.5))*l(17);
    A(42,L+13) = -1;
    A(42,L+17) = -pnu2.y()/sqrt( pnu2.Mag2() );

    // F43 (pnu2z)
    A(43,dimM+19) = ((pnu2.x()*pnu2.z())/pow( pnu2.Mag2(), 1.5))*l(17);
    A(43,dimM+20) = ((pnu2.y()*pnu2.z())/pow( pnu2.Mag2(), 1.5))*l(17);
    A(43,dimM+21) = -(( pow(pnu2.x(),2) + pow(pnu2.y(),2) )/pow( pnu2.Mag2(), 1.5))*l(17);
    A(43,L+14) = -1;
    A(43,L+17) = -pnu2.z()/sqrt( pnu2.Mag2() );

    // F44 (Enu2)
    A(44,L+15) = -1;
    A(44,L+17) = 1;

    ///////////////////////////// lambda /////////////////////////
    // F45 (l0)
    A(45,0) = pB.z();
    A(45,2) = -pB.x();
    A(45,dimM) = -pB.z();
    A(45,dimM+2) = pB.x();
    A(45,dimM+3) = BV.z()-PV.z();
    A(45,dimM+5) = -(BV.x()-PV.x());

    // F46 (l1)
    A(46,1) = pB.z();
    A(46,2) = -pB.y();
    A(46,dimM+1) = -pB.z();
    A(46,dimM+2) = pB.y();
    A(46,dimM+4) = BV.z()-PV.z();
    A(46,dimM+5) = -(BV.y()-PV.y());

    // F47 (l2)
    A(47,3) = -ptau1.z();
    A(47,5) = ptau1.x();
    A(47,dimM) = ptau1.z();
    A(47,dimM+2) = -ptau1.x();
    A(47,dimM+7) = DV1.z()-BV.z();
    A(47,dimM+9) = -(DV1.x()-BV.x());

    // F48 (l3)
    A(48,4) = -ptau1.z();
    A(48,5) = ptau1.y();
    A(48,dimM+1) = ptau1.z();
    A(48,dimM+2) = -ptau1.y();
    A(48,dimM+8) = DV1.z()-BV.z();
    A(48,dimM+9) = -(DV1.y()-BV.y());

    // F49 (l4)
    A(49,dimM+7) = 1;
    A(49,6) = -1;
    A(49,dimM+11) = -1;

    // F50 (l5)
    A(50,dimM+8) = 1;
    A(50,7) = -1;
    A(50,dimM+12) = -1;

    // F51 (l6)
    A(51,dimM+9) = 1;
    A(51,8) = -1;
    A(51,dimM+13) = -1;

    // F52 (l7)
    A(52,dimM+10) = 1;
    A(52,9) = -1;
    A(52,dimM+14) = -1;

    // F53 (l8)
    A(53,dimM+7) = -ptau1.x()/Etau1; 
    A(53,dimM+8) = -ptau1.y()/Etau1; 
    A(53,dimM+9) = -ptau1.z()/Etau1; 
    A(53,dimM+10) = 1;

    // F54 (l9)
    A(54,dimM+11) = -pnu1.x()/Enu1;
    A(54,dimM+12) = -pnu1.y()/Enu1;
    A(54,dimM+13) = -pnu1.z()/Enu1;
    A(54,dimM+14) = 1;

    // F55 (l10)
    A(55,10) = -ptau2.z();
    A(55,12) = ptau2.x();
    A(55,dimM) = ptau2.z();
    A(55,dimM+2) = -ptau2.x();
    A(55,dimM+15) = DV2.z()-BV.z();
    A(55,dimM+17) = -(DV2.x()-BV.x());

    // F56 (l11)
    A(56,11) = -ptau2.z();
    A(56,12) = ptau2.y();
    A(56,dimM+1) = ptau2.z();
    A(56,dimM+2) = -ptau2.y();
    A(56,dimM+16) = DV2.z()-BV.z();
    A(56,dimM+17) = -(DV2.y()-BV.y());

    // F57 (l12)
    A(57,dimM+15) = 1;
    A(57,13) = -1;
    A(57,dimM+19) = -1;

    // F58 (l13)
    A(58,dimM+16) = 1;
    A(58,14) = -1;
    A(58,dimM+20) = -1;

    // F59 (l14)
    A(59,dimM+17) = 1;
    A(59,15) = -1;
    A(59,dimM+21) = -1;

    // F60 (l15)
    A(60,dimM+18) = 1;
    A(60,16) = -1;
    A(60,dimM+22) = -1;

    // F61 (l16)
    A(61,dimM+15) = -ptau2.x()/Etau2;
    A(61,dimM+16) = -ptau2.y()/Etau2;
    A(61,dimM+17) = -ptau2.z()/Etau2;
    A(61,dimM+18) = 1;

    // F62 (l17)
    A(62,dimM+19) = -pnu2.x()/Enu2;
    A(62,dimM+20) = -pnu2.y()/Enu2;
    A(62,dimM+21) = -pnu2.z()/Enu2;
    A(62,dimM+22) = 1;

    // F63 (l18)
    A(63,17) = pK.z();
    A(63,19) = BV.z()-RP.z();
    A(63,21) = -(BV.x()-RP.x());
    A(63,dimM) = -pK.z();
    A(63,dimM+2) = pK.x();

    // F64 (l19)
    A(64,18) = pK.z();
    A(64,20) = BV.z()-RP.z();
    A(64,21) = -(BV.y()-RP.y());
    A(64,dimM+1) = -pK.z();
    A(64,dimM+2) = pK.y(); 

    // F65 (l20)
    A(65,dimM+3) = 1;
    A(65,dimM+7) = -1;
    A(65,dimM+15) = -1;
    A(65,19) = -1;

    // F66 (l21)
    A(66,dimM+4) = 1;
    A(66,dimM+8) = -1;
    A(66,dimM+16) = -1;
    A(66,20) = -1;

    // F67 (l22)
    A(67,dimM+5) = 1;
    A(67,dimM+9) = -1;
    A(67,dimM+17) = -1;
    A(67,21) = -1;

    // F68 (l23)
    A(68,19) = -pK.x()/EK;
    A(68,20) = -pK.y()/EK;
    A(68,21) = -pK.z()/EK;
    A(68,dimM+3) = pB.x()/EB;
    A(68,dimM+4) = pB.y()/EB;
    A(68,dimM+5) = pB.z()/EB;
    A(68,dimM+6) = 1./(2*EB);
    A(68,dimM+10) = -1;
    A(68,dimM+18) = -1;
    // A.Print();

    // TMatrixDSym A_sym(dimM+dimX+dimC);
    // for(int i = 0; i < dimM+dimX+dimC; i++)
    // {
    //     for(int j = 0; j < dimM+dimX+dimC; j++)
    //     {
    //         A_sym(i,j) = A(i,j);
    //     }
    // }

    // cout << "|A| = " << A_sym.Determinant() << endl;

    // TVectorD A_eigen_values(dimM+dimX+dimC);  
    // TMatrixD A_eigen_vectors = A.EigenVectors(A_eigen_values);
    // A_eigen_values.Print();

    // Double_t value = 1;
    // Double_t value_before = 1;
    // for(int i = 0; i < dimM+dimX+dimC; i++)
    // {
    //     value_before = value;
    //     value *= A_eigen_values(i);

    //     // if there are at least 2 consecutive entries with different sign, there are positive and negative eigenvalues -> it's a saddle point
    //     if((i!=0) && (TMath::Sign(1,value_before) != TMath::Sign(1,value)))
    //     {
    //         saddle_point = true;
    //     }
    // }

    TMatrixD A_inverse = A;
    A_inverse.Invert();

    // TMatrixD identity = A_inverse*A_sym;
    // identity.Print();

    // C = - A^-1 B
    TMatrixD C = (A_inverse*B);
    C *= -1;
    // C.Print();

    // Covariance matrix of fit parameters U = C V C^T
    TMatrixD C_transpose = C;
    C_transpose.T();

    TMatrixD U = C*(V*C_transpose);

    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        for(int j = 0; j < dimM+dimX+dimC; j++)
        {
            U_sym(i,j) = U(i,j);
        }
    }
    // U_sym.Print();

    return U_sym;
}

double eq1( const double* x_vars )
{
    eq_flag = 0;
    return equations(x_vars);
}
double eq2( const double* x_vars )
{
    eq_flag = 1;
    return equations(x_vars);
}
double eq3( const double* x_vars )
{
    eq_flag = 2;
    return equations(x_vars);
}
double eq4( const double* x_vars )
{
    eq_flag = 3;
    return equations(x_vars);
}
double eq5( const double* x_vars )
{
    eq_flag = 4;
    return equations(x_vars);
}
double eq6( const double* x_vars )
{
    eq_flag = 5;
    return equations(x_vars);
}
double eq7( const double* x_vars )
{
    eq_flag = 6;
    return equations(x_vars);
}
double eq8( const double* x_vars )
{
    eq_flag = 7;
    return equations(x_vars);
}
double eq9( const double* x_vars )
{
    eq_flag = 8;
    return equations(x_vars);
}
double eq10( const double* x_vars )
{
    eq_flag = 9;
    return equations(x_vars);
}
double eq11( const double* x_vars )
{
    eq_flag = 10;
    return equations(x_vars);
}
double eq12( const double* x_vars )
{
    eq_flag = 11;
    return equations(x_vars);
}
double eq13( const double* x_vars )
{
    eq_flag = 12;
    return equations(x_vars);
}
double eq14( const double* x_vars )
{
    eq_flag = 13;
    return equations(x_vars);
}
double eq15( const double* x_vars )
{
    eq_flag = 14;
    return equations(x_vars);
}
double eq16( const double* x_vars )
{
    eq_flag = 15;
    return equations(x_vars);
}
double eq17( const double* x_vars )
{
    eq_flag = 16;
    return equations(x_vars);
}
double eq18( const double* x_vars )
{
    eq_flag = 17;
    return equations(x_vars);
}
double eq19( const double* x_vars )
{
    eq_flag = 18;
    return equations(x_vars);
}
double eq20( const double* x_vars )
{
    eq_flag = 19;
    return equations(x_vars);
}
double eq21( const double* x_vars )
{
    eq_flag = 20;
    return equations(x_vars);
}
double eq22( const double* x_vars )
{
    eq_flag = 21;
    return equations(x_vars);
}
double eq23( const double* x_vars )
{
    eq_flag = 22;
    return equations(x_vars);
}
double eq24( const double* x_vars )
{
    eq_flag = 23;
    return equations(x_vars);
}
double eq25( const double* x_vars )
{
    eq_flag = 24;
    return equations(x_vars);
}
double eq26( const double* x_vars )
{
    eq_flag = 25;
    return equations(x_vars);
}
double eq27( const double* x_vars )
{
    eq_flag = 26;
    return equations(x_vars);
}
double eq28( const double* x_vars )
{
    eq_flag = 27;
    return equations(x_vars);
}
double eq29( const double* x_vars )
{
    eq_flag = 28;
    return equations(x_vars);
}
double eq30( const double* x_vars )
{
    eq_flag = 29;
    return equations(x_vars);
}
double eq31( const double* x_vars )
{
    eq_flag = 30;
    return equations(x_vars);
}
double eq32( const double* x_vars )
{
    eq_flag = 31;
    return equations(x_vars);
}
double eq33( const double* x_vars )
{
    eq_flag = 32;
    return equations(x_vars);
}
double eq34( const double* x_vars )
{
    eq_flag = 33;
    return equations(x_vars);
}
double eq35( const double* x_vars )
{
    eq_flag = 34;
    return equations(x_vars);
}
double eq36( const double* x_vars )
{
    eq_flag = 35;
    return equations(x_vars);
}
double eq37( const double* x_vars )
{
    eq_flag = 36;
    return equations(x_vars);
}
double eq38( const double* x_vars )
{
    eq_flag = 37;
    return equations(x_vars);
}
double eq39( const double* x_vars )
{
    eq_flag = 38;
    return equations(x_vars);
}
double eq40( const double* x_vars )
{
    eq_flag = 39;
    return equations(x_vars);
}
double eq41( const double* x_vars )
{
    eq_flag = 40;
    return equations(x_vars);
}
double eq42( const double* x_vars )
{
    eq_flag = 41;
    return equations(x_vars);
}
double eq43( const double* x_vars )
{
    eq_flag = 42;
    return equations(x_vars);
}
double eq44( const double* x_vars )
{
    eq_flag = 43;
    return equations(x_vars);
}
double eq45( const double* x_vars )
{
    eq_flag = 44;
    return equations(x_vars);
}
double eq46( const double* x_vars )
{
    eq_flag = 45;
    return equations(x_vars);
}
double eq47( const double* x_vars )
{
    eq_flag = 46;
    return equations(x_vars);
}
double eq48( const double* x_vars )
{
    eq_flag = 47;
    return equations(x_vars);
}
double eq49( const double* x_vars )
{
    eq_flag = 48;
    return equations(x_vars);
}
double eq50( const double* x_vars )
{
    eq_flag = 49;
    return equations(x_vars);
}
double eq51( const double* x_vars )
{
    eq_flag = 50;
    return equations(x_vars);
}
double eq52( const double* x_vars )
{
    eq_flag = 51;
    return equations(x_vars);
}
double eq53( const double* x_vars )
{
    eq_flag = 52;
    return equations(x_vars);
}
double eq54( const double* x_vars )
{
    eq_flag = 53;
    return equations(x_vars);
}
double eq55( const double* x_vars )
{
    eq_flag = 54;
    return equations(x_vars);
}
double eq56( const double* x_vars )
{
    eq_flag = 55;
    return equations(x_vars);
}
double eq57( const double* x_vars )
{
    eq_flag = 56;
    return equations(x_vars);
}
double eq58( const double* x_vars )
{
    eq_flag = 57;
    return equations(x_vars);
}
double eq59( const double* x_vars )
{
    eq_flag = 58;
    return equations(x_vars);
}
double eq60( const double* x_vars )
{
    eq_flag = 59;
    return equations(x_vars);
}
double eq61( const double* x_vars )
{
    eq_flag = 60;
    return equations(x_vars);
}
double eq62( const double* x_vars )
{
    eq_flag = 61;
    return equations(x_vars);
}
double eq63( const double* x_vars )
{
    eq_flag = 62;
    return equations(x_vars);
}
double eq64( const double* x_vars )
{
    eq_flag = 63;
    return equations(x_vars);
}
double eq65( const double* x_vars )
{
    eq_flag = 64;
    return equations(x_vars);
}
double eq66( const double* x_vars )
{
    eq_flag = 65;
    return equations(x_vars);
}
double eq67( const double* x_vars )
{
    eq_flag = 66;
    return equations(x_vars);
}
double eq68( const double* x_vars )
{
    eq_flag = 67;
    return equations(x_vars);
}
double eq69( const double* x_vars )
{
    eq_flag = 68;
    return equations(x_vars);
}

double equations( const double* x_vars )
{
    // 69 equations to be solved
    TVectorD x(dimM+dimX+dimC);
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        x(i) = x_vars[i]; // I need to unnormalise x because the equations are written with the dimensional quantities
    }

    // Measured parameters (22)
    ROOT::Math::XYZPoint PV( x(0), x(1), x(2) );
    ROOT::Math::XYZPoint DV1( x(3), x(4), x(5) );
    ROOT::Math::XYZVector p3pi1( x(6), x(7), x(8) );
    Double_t E3pi1 = x(9);
    ROOT::Math::XYZPoint DV2( x(10), x(11), x(12) );
    ROOT::Math::XYZVector p3pi2( x(13), x(14), x(15) );
    Double_t E3pi2 = x(16);
    ROOT::Math::XYZPoint RP( x(17), x(18), RPz );
    ROOT::Math::XYZVector pK( x(19), x(20), x(21) );
    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

    // Unknown parameters (19)
    ROOT::Math::XYZPoint BV( x(dimM), x(dimM+1), x(dimM+2) );
    ROOT::Math::XYZVector pB( x(dimM+3), x(dimM+4), x(dimM+5) );
    Double_t MB_squared = x(dimM+6);
    ROOT::Math::XYZVector ptau1( x(dimM+7), x(dimM+8), x(dimM+9) );
    Double_t Etau1 = x(dimM+10);
    ROOT::Math::XYZVector pnu1( x(dimM+11), x(dimM+12), x(dimM+13) );
    Double_t Enu1 = x(dimM+14);
    ROOT::Math::XYZVector ptau2( x(dimM+15), x(dimM+16), x(dimM+17) );
    Double_t Etau2 = x(dimM+18);
    ROOT::Math::XYZVector pnu2( x(dimM+19), x(dimM+20), x(dimM+21) );
    Double_t Enu2 = x(dimM+22);
    Double_t EB = sqrt( MB_squared + pB.Mag2() );

    // Lagrange multipliers (20)
    TVectorD l(dimC);
    for(int i = 0; i < dimC;i++)
    {
        l(i) = x(dimM+dimX+i);
    }

    // Get terms sum_i Wki (mi - xm^i) that are common for derivatives of the chi2 wrt xm^k 
    TVectorD chi2_sum(dimM);
    for(int k = 0; k < dimM; k++)
    {
        for(int i = 0; i < dimM; i++)
        {
            chi2_sum(k) += W(k,i)*(m(i) - x(i));
        }
    }

    // Write 61 equations:
    TVectorD eqs(dimM+dimX+dimC);
    // PV
    eqs(0) = -2.*chi2_sum(0) + pB.z()*l(0);
    eqs(1) = -2.*chi2_sum(1) + pB.z()*l(1);
    eqs(2) = -2.*chi2_sum(2) - pB.x()*l(0) - pB.y()*l(1);
    // DV1
    eqs(3) = -2.*chi2_sum(3) - ptau1.z()*l(2);
    eqs(4) = -2.*chi2_sum(4) - ptau1.z()*l(3);
    eqs(5) = -2.*chi2_sum(5) + ptau1.x()*l(2) + ptau1.y()*l(3);
    // P3pi1
    eqs(6) = -2.*chi2_sum(6) - l(4);
    eqs(7) = -2.*chi2_sum(7) - l(5);
    eqs(8) = -2.*chi2_sum(8) - l(6);
    eqs(9) = -2.*chi2_sum(9) - l(7);
    // DV2
    eqs(10) = -2.*chi2_sum(10) - ptau2.z()*l(10);
    eqs(11) = -2.*chi2_sum(11) - ptau2.z()*l(11);
    eqs(12) = -2.*chi2_sum(12) + ptau2.x()*l(10) + ptau2.y()*l(11);
    // P3pi2
    eqs(13) = -2.*chi2_sum(13) - l(12);
    eqs(14) = -2.*chi2_sum(14) - l(13);
    eqs(15) = -2.*chi2_sum(15) - l(14);
    eqs(16) = -2.*chi2_sum(16) - l(15);
    // RP_T
    eqs(17) = -2.*chi2_sum(17) + pK.z()*l(18);
    eqs(18) = -2.*chi2_sum(18) + pK.z()*l(19);
    // pK
    eqs(19) = -2.*chi2_sum(19) - l(20) - (pK.x()/EK)*l(23) + (BV.z()-RP.z())*l(18);
    eqs(20) = -2.*chi2_sum(20) - l(21) - (pK.y()/EK)*l(23) + (BV.z()-RP.z())*l(19);
    eqs(21) = -2.*chi2_sum(21) - l(22) - (pK.z()/EK)*l(23) - (BV.x()-RP.x())*l(18) - (BV.y()-RP.y())*l(19);
    // BV
    eqs(dimM) = -pB.z()*l(0) + ptau1.z()*l(2) + ptau2.z()*l(10) - pK.z()*l(18);
    eqs(dimM+1) = -pB.z()*l(1) + ptau1.z()*l(3) + ptau2.z()*l(11) -pK.z()*l(19);
    eqs(dimM+2) = pB.x()*l(0) + pB.y()*l(1) - ptau1.x()*l(2) - ptau1.y()*l(3) - ptau2.x()*l(10) - ptau2.y()*l(11) + pK.x()*l(18) + pK.y()*l(19);
    // pB, mB^2
    eqs(dimM+3) = (BV.z()-PV.z())*l(0) + l(20) + (pB.x()/EB)*l(23);
    eqs(dimM+4) = (BV.z()-PV.z())*l(1) + l(21) + (pB.y()/EB)*l(23);
    eqs(dimM+5) = -(BV.x()-PV.x())*l(0) - (BV.y()-PV.y())*l(1) + l(22) + (pB.z()/EB)*l(23);
    eqs(dimM+6) = (1./(2.*EB))*l(23);
    // ptau1
    eqs(dimM+7) = (DV1.z()-BV.z())*l(2) + l(4) - (ptau1.x()/sqrt(pow(mtau,2) + ptau1.Mag2()))*l(8) - l(20);
    eqs(dimM+8) = (DV1.z()-BV.z())*l(3) + l(5) - (ptau1.y()/sqrt(pow(mtau,2) + ptau1.Mag2()))*l(8) - l(21);
    eqs(dimM+9) = -(DV1.x()-BV.x())*l(2) - (DV1.y()-BV.y())*l(3) + l(6) - (ptau1.z()/sqrt(pow(mtau,2) + ptau1.Mag2()))*l(8) - l(22);
    eqs(dimM+10) = l(7) + l(8) - l(23);
    // pnu1
    eqs(dimM+11) = -l(4) - (pnu1.x()/sqrt(pnu1.Mag2()))*l(9);
    eqs(dimM+12) = -l(5) - (pnu1.y()/sqrt(pnu1.Mag2()))*l(9);
    eqs(dimM+13) = -l(6) - (pnu1.z()/sqrt(pnu1.Mag2()))*l(9);
    eqs(dimM+14) = -l(7) + l(9);
    // ptau2
    eqs(dimM+15) = (DV2.z()-BV.z())*l(10) + l(12) - (ptau2.x()/sqrt(pow(mtau,2) + ptau2.Mag2()))*l(16) - l(20);
    eqs(dimM+16) = (DV2.z()-BV.z())*l(11) + l(13) - (ptau2.y()/sqrt(pow(mtau,2) + ptau2.Mag2()))*l(16) - l(21);
    eqs(dimM+17) = -(DV2.x()-BV.x())*l(10) - (DV2.y()-BV.y())*l(11) + l(14) - (ptau2.z()/sqrt(pow(mtau,2) + ptau2.Mag2()))*l(16) - l(22);
    eqs(dimM+18) = l(15) + l(16) - l(23);
    // pnu2
    eqs(dimM+19) = -l(12) - (pnu2.x()/sqrt(pnu2.Mag2()))*l(17);
    eqs(dimM+20) = -l(13) - (pnu2.y()/sqrt(pnu2.Mag2()))*l(17);
    eqs(dimM+21) = -l(14) - (pnu2.z()/sqrt(pnu2.Mag2()))*l(17);
    eqs(dimM+22) = -l(15) + l(17);
    // pB must point back to the PV (2)
    eqs(dimM+dimX) = pB.x()*(BV.z() - PV.z()) - pB.z()*(BV.x() - PV.x());
    eqs(dimM+dimX+1) = pB.y()*(BV.z() - PV.z()) - pB.z()*(BV.y() - PV.y());
    // ptau1 must point back to the BV (2)
    eqs(dimM+dimX+2) = ptau1.x()*(DV1.z() - BV.z()) - ptau1.z()*(DV1.x() - BV.x());
    eqs(dimM+dimX+3) = ptau1.y()*(DV1.z() - BV.z()) - ptau1.z()*(DV1.y() - BV.y());
    // 4-momentum conservation in DV1 (4)
    eqs(dimM+dimX+4) = ptau1.x() - p3pi1.x() - pnu1.x();
    eqs(dimM+dimX+5) = ptau1.y() - p3pi1.y() - pnu1.y();
    eqs(dimM+dimX+6) = ptau1.z() - p3pi1.z() - pnu1.z();
    eqs(dimM+dimX+7) = Etau1 - E3pi1 - Enu1;
    // tau+ and anti-nu mass constraints (2)
    eqs(dimM+dimX+8) = Etau1 - sqrt(pow(mtau,2) + ptau1.Mag2());
    eqs(dimM+dimX+9) = Enu1 - sqrt(pnu1.Mag2());
    // ptau2 must point back to the BV (2)
    eqs(dimM+dimX+10) = ptau2.x()*(DV2.z() - BV.z()) - ptau2.z()*(DV2.x() - BV.x());
    eqs(dimM+dimX+11) = ptau2.y()*(DV2.z() - BV.z()) - ptau2.z()*(DV2.y() - BV.y());
    // 4-momentum conservation in DV2 (4)
    eqs(dimM+dimX+12) = ptau2.x() - p3pi2.x() - pnu2.x();
    eqs(dimM+dimX+13) = ptau2.y() - p3pi2.y() - pnu2.y();
    eqs(dimM+dimX+14) = ptau2.z() - p3pi2.z() - pnu2.z();
    eqs(dimM+dimX+15) = Etau2 - E3pi2 - Enu2;
    // tau- and nu mass constraints (2)
    eqs(dimM+dimX+16) = Etau2 - sqrt(pow(mtau,2) + ptau2.Mag2());
    eqs(dimM+dimX+17) = Enu2 - sqrt(pnu2.Mag2());
    // BV must lie in K+ trajectory (2)
    eqs(dimM+dimX+18) = pK.x()*(BV.z() - RP.z()) - pK.z()*(BV.x() - RP.x());
    eqs(dimM+dimX+19) = pK.y()*(BV.z() - RP.z()) - pK.z()*(BV.y() - RP.y());
    // 4-momentum conservation in BV (4)
    eqs(dimM+dimX+20) = pB.x() - ptau1.x() - ptau2.x() - pK.x();
    eqs(dimM+dimX+21) = pB.y() - ptau1.y() - ptau2.y() - pK.y();
    eqs(dimM+dimX+22) = pB.z() - ptau1.z() - ptau2.z() - pK.z();
    eqs(dimM+dimX+23) = EB - Etau1 - Etau2 - EK;
 
    if( (eq_flag < 0) || (eq_flag > dimM+dimX+dimC))
    {
        cout << "Not a valid flag; should be between 0 - 69" << endl;
        return 0;
    }
    else
    {
        return eqs(eq_flag);
    }
}

Double_t lagrangian( TVectorD x )
{
    // 1) Separate x into xm, xu and lambda
    TVectorD xm(dimM);
    for(int i = 0; i < dimM; i++)
    {
        xm(i) = x(i);
    }

    TVectorD xu(dimX);
    for(int i = 0; i < dimX; i++)
    {
        xu(i) = x(dimM+i);
    }

    TVectorD lambda(dimC);
    for(int i = 0; i < dimC; i++)
    {
        lambda(i) = x(dimM+dimX+i);
    }

    // 2) Build the chi^2
    Double_t chi2 = chisquare(xm);

    // 3) Build the vector of exact constraints g(xm,xu)
    TVectorD g = exact_constraints( x );

    // 4) Build the Lagrangian
    Double_t L = chi2 + lambda*g;
    
    return L;
}

TVectorD exact_constraints( TVectorD x )
{
    ROOT::Math::XYZPoint PV( x(0), x(1), x(2) );
    ROOT::Math::XYZPoint DV1( x(3), x(4), x(5) );
    ROOT::Math::XYZVector p3pi1( x(6), x(7), x(8) );
    Double_t E3pi1 = x(9);
    ROOT::Math::XYZPoint DV2( x(10), x(11), x(12) );
    ROOT::Math::XYZVector p3pi2( x(13), x(14), x(15) );
    Double_t E3pi2 = x(16);
    ROOT::Math::XYZPoint RP( x(17), x(18), RPz );
    ROOT::Math::XYZVector pK( x(19), x(20), x(21) );
    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

    ROOT::Math::XYZPoint BV( x(dimM), x(dimM+1), x(dimM+2) );
    ROOT::Math::XYZVector pB( x(dimM+3), x(dimM+4), x(dimM+5) );
    Double_t MB_squared = x(dimM+6);
    ROOT::Math::XYZVector ptau1( x(dimM+7), x(dimM+8), x(dimM+9) );
    ROOT::Math::XYZVector pnu1( x(dimM+10), x(dimM+11), x(dimM+12) );
    ROOT::Math::XYZVector ptau2( x(dimM+13), x(dimM+14), x(dimM+15) );
    ROOT::Math::XYZVector pnu2( x(dimM+16), x(dimM+17), x(dimM+18) );
    Double_t EB = sqrt( MB_squared + pB.Mag2() );
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );

    TVectorD g(dimC);
    // pB must point back to the PV (2)
    g(0) = pB.x()*( BV.z() - PV.z() ) - pB.z()*( BV.x() - PV.x() );
    g(1) = pB.y()*( BV.z() - PV.z() ) - pB.z()*( BV.y() - PV.y() );
    // ptau1 must point back to the BV (2)
    g(2) = ptau1.x()*( DV1.z() - BV.z() ) - ptau1.z()*( DV1.x() - BV.x() );
    g(3) = ptau1.y()*( DV1.z() - BV.z() ) - ptau1.z()*( DV1.y() - BV.y() );
    // 4-momentum conservation in DV1 (4)
    g(4) = ptau1.x() - p3pi1.x() - pnu1.x();
    g(5) = ptau1.y() - p3pi1.y() - pnu1.y();
    g(6) = ptau1.z() - p3pi1.z() - pnu1.z();
    g(7) = Etau1 - E3pi1 - Enu1;
    // tau+ and anti-nu mass constraints (2)
    g(8) = Etau1 - sqrt( pow(mtau,2) + ptau1.Mag2() );
    g(9) = Enu1 - sqrt( pnu1.Mag2() );
    // ptau2 must point back to the BV (2)
    g(10) = ptau2.x()*( DV2.z() - BV.z() ) - ptau2.z()*( DV2.x() - BV.x() );
    g(11) = ptau2.y()*( DV2.z() - BV.z() ) - ptau2.z()*( DV2.y() - BV.y() );
    // 4-momentum conservation in DV2 (4)
    g(12) = ptau2.x() - p3pi2.x() - pnu2.x();
    g(13) = ptau2.y() - p3pi2.y() - pnu2.y();
    g(14) = ptau2.z() - p3pi2.z() - pnu2.z();
    g(15) = Etau2 - E3pi2 - Enu2;
    // tau- and nu mass constraints (2)
    g(16) = Etau2 - sqrt( pow(mtau,2) + ptau2.Mag2() );
    g(17) = Enu2 - sqrt( pnu2.Mag2() );
    // BV must lie in K+ trajectory (2)
    g(18) = pK.x()*( BV.z() - RP.z() ) - pK.z()*( BV.x() - RP.x() );
    g(19) = pK.y()*( BV.z() - RP.z() ) - pK.z()*( BV.y() - RP.y() );
    // 4-momentum conservation in BV (4)
    g(20) = pB.x() - ptau1.x() - ptau2.x() - pK.x();
    g(21) = pB.y() - ptau1.y() - ptau2.y() - pK.y();
    g(22) = pB.z() - ptau1.z() - ptau2.z() - pK.z();
    g(23) = EB - Etau1 - Etau2 - EK;

    return g;
}

Double_t chisquare( TVectorD xm )
{
    // Computes the chi^2 = (m - xm)^T W (m - xm)
    TVectorD r = m - xm;

    Double_t chi2 = r*(W*r);

    return chi2;
}


TVectorD x_initial_estimate_0( TVectorD m ) // Original initialisation for x (based on Anne Keune's thesis)
{
    // Builds an initial estimate for x based on the analytical calculations and on the known parameters in m
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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
    Double_t BVz_t = i - j*(1./Ptau1x_t);

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
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_1( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; taus direction based on vertices
{
    // Builds an initial estimate for x using B->K* tautau initialisation; taus direction based on vertices
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );

    Double_t EB = Etau1 + Etau2 + EK;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_2( TVectorD m, ROOT::Math::XYZPoint BV )  // B->K* tautau initialisation; taus direction based on visible 3pi momenta
{
    // Builds an initial estimate for x using B->K* tautau initialisation; taus direction based on visible 3pi momenta
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    // Get tau flight directions (from pions visible momenta)
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
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );

    Double_t EB = Etau1 + Etau2 + EK;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_3( TVectorD m, ROOT::Math::XYZPoint BV ) // RD initialisation
{
    // Builds an initial estimate for x based on the Marseille analytical calculations; it uses the offline estimate for BV  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

    Double_t Enu1 = Etau1 - E3pi1;
    Double_t Enu2 = Etau2 - E3pi2;

    // Double_t Enu1 = sqrt( pnu1.Mag2() );
    // Double_t Enu2 = sqrt( pnu2.Mag2() );

    ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
    Double_t EB = EK + Etau1 + Etau2;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_4( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; tau+ direction from vertices ; tau- direction from pions
{
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_5( TVectorD m, ROOT::Math::XYZPoint BV ) // B->K* tautau initialisation; tau- direction from vertices ; tau+ direction from pions
{
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_6( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau+) and K* tau tau vertices (tau-)
{  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_7( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau-) and K* tau tau vertices (tau+)
{  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_8( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau+) and K* tau tau pions (tau-)
{  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_9( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau-) and K* tau tau pions (tau+)
{  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
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

    x0(dimM) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

    return x0;
}

TVectorD x_initial_estimate_true( TVectorD m )
{  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
    x0(dimM) = 0;
    x0(dimM+1) = 0;
    x0(dimM+2) = 0;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    return x0;
}

TVectorD x_initial_estimate_MLP( TVectorD m, ROOT::Math::XYZPoint BV, Int_t species )
{
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu 
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

    if((species == 10) || (species == 11) || (species == 12) || (species == 2) || (species == 3))
    {
        weightfile_taup_PX = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP.weights.xml";
        weightfile_taup_PY = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP.weights.xml";
        weightfile_taup_PZ = "/home/avenkate/jobs/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP.weights.xml";
        weightfile_taum_PX = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP.weights.xml";
        weightfile_taum_PY = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP.weights.xml";
        weightfile_taum_PZ = "/home/avenkate/jobs/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP.weights.xml";
    }
    else if(species == 4)
    {
        weightfile_taup_PX = "/home/avenkate/jobs/DDK_MLP_Train_taup_PX/dataset/weights/TMVARegression_Dp_TRUEP_X_MLP.weights.xml";
        weightfile_taup_PY = "/home/avenkate/jobs/DDK_MLP_Train_taup_PY/dataset/weights/TMVARegression_Dp_TRUEP_Y_MLP.weights.xml";
        weightfile_taup_PZ = "/home/avenkate/jobs/DDK_MLP_Train_taup_PZ/dataset/weights/TMVARegression_Dp_TRUEP_Z_MLP.weights.xml";
        weightfile_taum_PX = "/home/avenkate/jobs/DDK_MLP_Train_taum_PX/dataset/weights/TMVARegression_Dm_TRUEP_X_MLP.weights.xml";
        weightfile_taum_PY = "/home/avenkate/jobs/DDK_MLP_Train_taum_PY/dataset/weights/TMVARegression_Dm_TRUEP_Y_MLP.weights.xml";
        weightfile_taum_PZ = "/home/avenkate/jobs/DDK_MLP_Train_taum_PZ/dataset/weights/TMVARegression_Dm_TRUEP_Z_MLP.weights.xml";
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

    x0(dimM) = BV.x(); 
    x0(dimM+1) = BV.y(); 
    x0(dimM+2) = BV.z() ; // -0.26
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = Etau1;
    x0(dimM+11) = pnu1.x();
    x0(dimM+12) = pnu1.y();
    x0(dimM+13) = pnu1.z();
    x0(dimM+14) = Enu1;
    x0(dimM+15) = ptau2.x();
    x0(dimM+16) = ptau2.y();
    x0(dimM+17) = ptau2.z();
    x0(dimM+18) = Etau2;
    x0(dimM+19) = pnu2.x();
    x0(dimM+20) = pnu2.y();
    x0(dimM+21) = pnu2.z();
    x0(dimM+22) = Enu2;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1./10000.; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1./1000.; // px conservation in DV1
    x0(dimM+dimX+5) = 1./1000.; // py conservation in DV1
    x0(dimM+dimX+6) = 1./10000.; // pz conservation in DV1
    x0(dimM+dimX+7) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+8) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+9) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+10) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+11) = 1./10000.; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+12) = 1./1000.; // px conservation in DV2
    x0(dimM+dimX+13) = 1./1000.; // py conservation in DV2
    x0(dimM+dimX+14) = 1./10000.; // pz conservation in DV2
    x0(dimM+dimX+15) = 1./10000.; // E conservation in DV2
    x0(dimM+dimX+16) = 1./10000.; // Tau mass constraint
    x0(dimM+dimX+17) = 1./10000.; // Nu mass constraint
    x0(dimM+dimX+18) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+19) = 1./1000.; // pKx/(BVx - RPx)
    x0(dimM+dimX+20) = 1./10000.; // px vonservation in BV
    x0(dimM+dimX+21) = 1./10000.; // py conservation in BV
    x0(dimM+dimX+22) = 1./100000.; // pz conservation in BV
    x0(dimM+dimX+23) = 1./100000.; // E conservation in BV

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