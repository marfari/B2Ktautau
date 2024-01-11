// Decay fit with 4 constraints to 3 unknowns
// Every constraint that is easy to apply by parameter substitution is applied in this way
// Remaining unknowns: (BVz,ptau1z,ptau2z)
// Remaining constraints:
// - pB must point back to the PV (2)
// - Energy conservation in DV1 and in DV2 (2)

using namespace std;

// Global variables 
int dimM = 22; // number of measured parameters
int dimX = 3; // number of unknown parameters
int dimC = dimX+1; // number of exact constraints
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
Double_t BVz_true, ptau1z_true, ptau2z_true;
TVectorD x0_current(dimM+dimX+dimC), x1_current(dimM+dimX+dimC), x2_current(dimM+dimX+dimC), x3_current(dimM+dimX+dimC), x4_current(dimM+dimX+dimC), x5_current(dimM+dimX+dimC), x6_current(dimM+dimX+dimC), x7_current(dimM+dimX+dimC), x8_current(dimM+dimX+dimC), x9_current(dimM+dimX+dimC), x_true_current;

// Fit results saved
Double_t MB, MB_err, F_tolerance, chi2, L;
TVectorD X(dimM+dimX+dimC);
TVectorD X_ERR(dimM+dimX+dimC);
TVectorD F(dimM+dimX+dimC);
TMatrixDSym U(dimM+dimX+dimC);
TMatrixDSym U_corr(dimM+dimX+dimC);
TMatrixDSym D(dimM+dimX+dimC);
TMatrixDSym D_inverse(dimM+dimX+dimC);
Int_t status, init, nIter;
TVectorD X_total(22+23);
TVectorD X_ERR_total(22+23);

// Constants
Double_t mtau = 1776.86;
Double_t mkaon = 493.677;

// Functions
void solve( ROOT::Math::XYZPoint BV, int init,  ROOT::Math::GSLMultiRootFinder* solver, Double_t tolerance );
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
void sequence(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder * solver, Double_t tolerance);
void lowest_chi2(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder *solver, Double_t tolerance);
Double_t lagrangian( TVectorD x );
Double_t chisquare( TVectorD xm );
TVectorD exact_constraints( TVectorD x );
TMatrixDSym U_matrix();
TMatrixDSym D_matrix();
vector<double> range(double min, double max, size_t N);
void scan_lagrangian( TVectorD X, Int_t index, Int_t npoints, TString x_name, Int_t year, Int_t species, Int_t line );
void scan_chisquare( TVectorD X, Int_t index, Int_t npoints, TString x_name, Int_t year, Int_t species, Int_t line);
TVectorD normalise( TVectorD a, bool norm );

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

void DECAY_FIT3(int year, TString RECO_files, int species, int line)
{   
    // Retrieve m and V from ntuple
    TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files, 1, line);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    UInt_t num_entries = t->GetEntries();
    Double_t m_vars[dimM];
    Double_t V_vars[dimM][dimM];
    Double_t BVx, BVy, BVz; // offline estimate for the BV is needed to have a first estimate for the unknown parameters using some of the initialisations
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

    // truth-match quantities
    t->SetBranchAddress("Bp_TRUEENDVERTEX_Z", &BVz_true);
    t->SetBranchAddress("taup_TRUEP_Z", &ptau1z_true);
    t->SetBranchAddress("taum_TRUEP_Z", &ptau2z_true);

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
        "df_3p1_PX_err",
        "df_3p1_PY_err",
        "df_3p1_PZ_err",
        "df_3p1_PE_err",
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

    for(int i = 0; i < 22+23; i++)
    {
        tout->Branch(name_x[i], &X_total(i));
        tout->Branch(name_x_err[i], &X_ERR_total(i));
    }

    for(int i = 0; i < dimM+dimX+dimC; i++)
    { 
        tout->Branch(Form("df_F_%i",i), &F(i));
    }

    tout->Branch("df_init", &init);
    tout->Branch("df_Bp_M", &MB);
    tout->Branch("df_Bp_MERR", &MB_err);
    tout->Branch("df_status", &status);
    tout->Branch("df_F_tolerance", &F_tolerance);
    tout->Branch("df_chi2", &chi2);
    tout->Branch("df_L", &L);
    tout->Branch("df_nIter", &nIter);

    // Loop over events
    for(int evt = 0; evt <  1; evt++)
    {
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

        // Matrix D for normalisation
        // D = D_matrix();
        // D_inverse = D;
        // D_inverse.Invert();

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
        solver->SetPrintLevel(0);

        // 4) Solve the system of equations (29)
        ROOT::Math::XYZPoint BV( BVx, BVy, BVz );
        // init = 0;
        // solve( BV, init, solver, pow(10,-6) );

        lowest_chi2(BV, solver, initial_tolerance); // initial_tolerance = pow(10,-10)
        if(status != 0)
        {
            lowest_chi2(BV, solver, pow(10,-6));   
        }
        if(status != 0)
        {
            lowest_chi2(BV, solver, pow(10,-3));
        }
        if(status != 0)
        {
            lowest_chi2(BV, solver, 1);
        }
        if(status != 0)
        {
            lowest_chi2(BV, solver, 100);
        }

        // truth initialisation
        // solve( BV, -1, solver, initial_tolerance ); // initial_tolerance = pow(10,-10)
        // if(status != 0)
        // {
        //     solve( BV, -1, solver, pow(10,-6) );   
        // }
        // if(status != 0)
        // {
        //     solve( BV, -1, solver, pow(10,-3) );
        // }
        // if(status != 0)
        // {
        //     solve( BV, -1, solver, 1 );
        // }
        // if(status != 0)
        // {
        //     solve( BV, -1, solver, 100 );
        // }

        cout << "FINAL" << endl;
        cout << "init = " << init << endl;
        cout << "status == " << status << endl;
        cout << "sum_Fi = " << F_tolerance << endl;
        cout << "MB = " << MB << " +/- " << MB_err << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "#iterations = " << nIter << endl;        

        // F.Print();
        // X.Print();
        // X_total.Print();
        // X_ERR.Print();

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
            "BVz (mm)",
            "ptau1z (MeV)",
            "ptau2z (MeV)",
            "lambda1",
            "lambda2",
            "lambda3",
            "lambda4"
        };
        if( sizeof(NAME)/sizeof(NAME[0]) != dimM+dimX+dimC )
        {
            cout << "NAME size is not correct" << endl;
            return;
        }

        // for(int i = 0; i < dimM+dimX+dimC; i++)
        // {
        //     scan_lagrangian(X, i, 100, NAME[i], year, species, line);
        // }
        // for(int i = 0; i < dimM; i++)
        // {
        //     scan_chisquare(X, i, 100, NAME[i], year, species, line);
        // }

        tout->Fill();
    }
    cout << "FINISHED" << endl;

    fout->cd();
    tout->Write();
    fout->Close();
}

void solve( ROOT::Math::XYZPoint BV, int init,  ROOT::Math::GSLMultiRootFinder* solver, Double_t tolerance )
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
    // x0.Print();

    Double_t x0_vars[dimM+dimX+dimC];
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        x0_vars[i] = x0(i); 
    }

    // 2) Solve system of equations
    solver->Solve(x0_vars, 1000, tolerance);

    // 3) Retrieve results
    const double* x_results = solver->X();
    const double* f_vals = solver->FVal();
    status = solver->Status();
    if(tolerance == pow(10,-10))
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

    if((init == -1) && (tolerance == initial_tolerance))
    {
        x_true_current = X;
    }

    // Measured parameters (22)
    ROOT::Math::XYZPoint PV_m( X(0), X(1), X(2) );
    ROOT::Math::XYZPoint DV1_m( X(3), X(4), X(5) );
    ROOT::Math::XYZVector p3pi1_m( X(6), X(7), X(8) );
    Double_t E3pi1_m = X(9);
    ROOT::Math::XYZPoint DV2_m( X(10), X(11), X(12) );
    ROOT::Math::XYZVector p3pi2_m( X(13), X(14), X(15) );
    Double_t E3pi2_m = X(16);
    ROOT::Math::XYZPoint RP_m( X(17), X(18), RPz );
    ROOT::Math::XYZVector pK_m( X(19), X(20), X(21) );
    Double_t EK_m = sqrt( pow(mkaon,2) + pK_m.Mag2() );

    // Unknown parameters
    Double_t BVz_m = X(dimM);
    Double_t ptau1z_m = X(dimM+1);
    Double_t ptau2z_m = X(dimM+2);

    // BV must lie in K+ trajectory
    Double_t BVx_m = RP_m.x() + (pK_m.x()/pK_m.z())*(BVz_m - RP_m.z());
    Double_t BVy_m = RP_m.y() + (pK_m.y()/pK_m.z())*(BVz_m - RP_m.z()); 
    ROOT::Math::XYZPoint BV_m( BVx_m, BVy_m, BVz_m );
    // ptau1 must point back to the BV
    Double_t ptau1x_m = ((DV1_m.x() - BV_m.x())/(DV1_m.z() - BV_m.z()))*ptau1z_m;
    Double_t ptau1y_m = ((DV1_m.y() - BV_m.y())/(DV1_m.z() - BV_m.z()))*ptau1z_m;
    ROOT::Math::XYZVector ptau1_m( ptau1x_m, ptau1y_m, ptau1z_m );
    // ptau2 must point back to the BV
    Double_t ptau2x_m = ((DV2_m.x() - BV_m.x())/(DV2_m.z() - BV_m.z()))*ptau2z_m;
    Double_t ptau2y_m = ((DV2_m.y() - BV_m.y())/(DV2_m.z() - BV_m.z()))*ptau2z_m;
    ROOT::Math::XYZVector ptau2_m( ptau2x_m, ptau2y_m, ptau2z_m );
    // Tau mass constraints
    Double_t Etau1_m = sqrt( pow(mtau,2) + ptau1_m.Mag2() );
    Double_t Etau2_m = sqrt( pow(mtau,2) + ptau2_m.Mag2() );
    // 3-momentum conservation in DV1 and in DV2
    ROOT::Math::XYZVector pnu1_m = ptau1_m - p3pi1_m;
    ROOT::Math::XYZVector pnu2_m = ptau2_m - p3pi2_m;
    // Neutrino mass constraints
    Double_t Enu1_m = sqrt( pnu1_m.Mag2() );
    Double_t Enu2_m = sqrt( pnu2_m.Mag2() );
    // 4-momentum conservation in BV
    ROOT::Math::XYZVector pB_m = ptau1_m + ptau2_m + pK_m;
    Double_t EB_m = Etau1_m + Etau2_m + EK_m;

    U = U_matrix();
    // for(int i = 0; i < dimM+dimX+dimC; i++)
    // {
    //     X_ERR(i) = sqrt(U(i,i));
    // }
    // X_ERR.Print();

    Double_t MB_squared = pow(EB_m,2) - pB_m.Mag2();
    if( MB_squared > 0 )
    {
        MB = sqrt( MB_squared );
    }
    else
    {
        MB = -sqrt( abs(MB_squared) );
    }

    for(int i = 0; i < dimM; i++)
    {
        X_total(i) = X(i);
    }
    X_total(dimM) = BV_m.x();
    X_total(dimM+1) = BV_m.y();
    X_total(dimM+2) = BV_m.z();
    X_total(dimM+3) = pB_m.x();
    X_total(dimM+4) = pB_m.y();
    X_total(dimM+5) = pB_m.z();
    X_total(dimM+6) = MB_squared;
    X_total(dimM+7) = ptau1_m.x();
    X_total(dimM+8) = ptau1_m.y();
    X_total(dimM+9) = ptau1_m.z();
    X_total(dimM+10) = Etau1_m;
    X_total(dimM+11) = pnu1_m.x();
    X_total(dimM+12) = pnu1_m.y();
    X_total(dimM+13) = pnu1_m.z();
    X_total(dimM+14) = Enu1_m;
    X_total(dimM+15) = ptau2_m.x();
    X_total(dimM+16) = ptau2_m.y();
    X_total(dimM+17) = ptau2_m.z();
    X_total(dimM+18) = Etau2_m;
    X_total(dimM+19) = pnu2_m.x();
    X_total(dimM+20) = pnu2_m.y();
    X_total(dimM+21) = pnu2_m.z();
    X_total(dimM+22) = Enu2_m;

    // Double_t dMB_squared = sqrt(U(28,28));
    // MB_err = dMB_squared/(2*abs(MB));

    // cout << "MB = " << MB << " +/- " << MB_err << endl;

    // TVectorD V_errors(dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //     V_errors(i) = sqrt(V(i,i));
    // }
    // m.Print();
    // V_errors.Print();

    // TVectorD diff(dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //     diff(i) = X(i)-m(i);
    // }
    // diff.Print();

    // TVectorD ratio(dimM);
    // for(int i = 0; i < dimM; i++)
    // {
    //     ratio(i) = X_ERR(i)/V_errors(i);
    // }
    // ratio.Print();

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
    L = lagrangian(X);

}

void lowest_chi2(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder *solver, Double_t tolerance)
{
        // If one or more initialisations pass with the current tolerance, it chooses the ones with the lowest chi^2
        // It saves the state from each initialisation in a TVectorD

        // Init 0
        solve( BV, 0, solver, tolerance );
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

        // Init 1
        solve( BV, 1, solver, tolerance );
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

        // Init 2
        solve( BV, 2, solver, tolerance );
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

        // Init 3
        solve( BV, 3, solver, tolerance );
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

        // Init 4
        solve( BV, 4, solver, tolerance );
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

        // Init 5
        solve( BV, 5, solver, tolerance );
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

        // Init 6
        solve( BV, 6, solver, tolerance );
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

        // Init 7
        solve( BV, 7, solver, tolerance );
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

        // Init 8
        solve( BV, 8, solver, tolerance );
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

        // Init 9
        solve( BV, 9, solver, tolerance );
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

        bool all_fail = false;
        if( (status0 != 0) && (status1 != 0) && (status2 != 0) && (status3 != 0) && (status4 != 0) && (status5 != 0) && (status6 != 0) && (status7 != 0) && (status8 != 0) && (status9 != 0) )
        {
            all_fail = true;
        }

        Double_t chi2_min = 100;
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

        // Double_t F_min = 100;
        // for(int i = 0; i < 10; i++)
        // {
        //     if( (F_tol_vec[i] < F_min) )
        //     {
        //         F_min = F_tol_vec[i];
        //         i_min = i;
        //     }
        // }   

        init = i_min;
        status = status_vec[i_min];
        MB = MB_vec[i_min];
        MB_err = MB_err_vec[i_min];
        X = X_vec[i_min];
        X_ERR = X_ERR_vec[i_min];
        F = F_vec[i_min];
        F_tolerance = F_tol_vec[i_min];
        chi2 = chi2_vec[i_min];
        nIter = nIter_vec[i_min];
}

void sequence(ROOT::Math::XYZPoint BV, ROOT::Math::GSLMultiRootFinder *solver, Double_t tolerance)
{
    init = 3; // Marseille
    solve( BV, init, solver, tolerance );
    if(status != 0)
    {
        init = 0; // Original
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 9; // Marseille (tau-) + K*tautau pions (tau+)
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 8; // Marseille (tau+) + K*tautau pions (tau-)
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 7; // Marseille (tau-) + K*tautau vertex (tau+)
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 6; // Marseille (tau+) + K*tautau vertex (tau-)
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 2; // K*tautau pions
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 5; // K*tautau vertex (tau-) + K*tautau pions (tau+)
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 4; // K*tautau vertex (tau+) + K*tautau pions (tau-)
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = 1; // K*tautau vertex
        solve( BV, init, solver, tolerance );
    }
    if(status != 0)
    {
        init = -1;
    }
}

TMatrixDSym D_matrix()
{
    TVectorD d(dimM+dimX+dimC);
    for(int i = 0; i < dimM; i++)
    {
        d(i) = sqrt(V(i,i));
    }

    d(dimM) = 0.5*pow(10,-3); // BVx
    d(dimM+1) = 0.5*pow(10,-3); // BVy
    d(dimM+2) = 50*pow(10,-3); // BVz
    d(dimM+3) = 0.005*10000; // pBx
    d(dimM+4) = 0.005*10000; // pBy
    d(dimM+5) = 0.005*100000; // pBz
    d(dimM+6) = 2*1000*100; // mB^2
    d(dimM+7) = 0.005*1000; // ptau1x
    d(dimM+8) = 0.005*1000; // ptau1y
    d(dimM+9) = 0.005*10000; // ptau1z
    d(dimM+10) = 0.005*1000; // pnu1x
    d(dimM+11) = 0.005*1000; // pnu1y
    d(dimM+12) = 0.005*10000; // pnu1z
    d(dimM+13) = 0.005*1000; // ptau2x
    d(dimM+14) = 0.005*1000; // ptau2y
    d(dimM+15) = 0.005*10000; // ptau2z
    d(dimM+16) = 0.005*1000; // pnu2x
    d(dimM+17) = 0.005*1000; // pnu2y
    d(dimM+18) = 0.005*10000; // pnu2z

    d(dimM+dimX) = 1./1000; // pBx/(BVx - PVx)
    d(dimM+dimX+1) = 1./1000.; // pBx/(BVx - PVx)
    d(dimM+dimX+2) = 1./1000; // ptau1x/(DV1x - BVx)
    d(dimM+dimX+3) = 1./1000; // ptau1x/(DV1x - BVx)
    d(dimM+dimX+4) = 1./100; // px conservation in DV1
    d(dimM+dimX+5) = 1./100; // py conservation in DV1
    d(dimM+dimX+6) = 1./1000; // pz conservation in DV1
    d(dimM+dimX+7) = 1./1000; // E conservation in DV1
    d(dimM+dimX+8) = 1./1000; // ptau2x/(DV2x - BVx)
    d(dimM+dimX+9) = 1./1000; // ptau2x/(DV2x - BVx)
    d(dimM+dimX+10) = 1./100; // px conservation in DV2
    d(dimM+dimX+11) = 1./100; // py conservation in DV2
    d(dimM+dimX+12) = 1./1000; // pz conservation in DV2
    d(dimM+dimX+13) = 1./1000; // E conservation in DV2
    d(dimM+dimX+14) = 1./100; // pKx/(BVx - RPx)
    d(dimM+dimX+15) = 1./100; // pKx/(BVx - RPx)
    d(dimM+dimX+16) = 1./1000; // px vonservation in BV
    d(dimM+dimX+17) = 1./1000; // py conservation in BV
    d(dimM+dimX+18) = 1./10000; // pz conservation in BV
    d(dimM+dimX+19) = 1./10000; // E conservation in BV
    // d.Print();

    TMatrixDSym Dmatrix(dimM+dimX+dimC);
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        Dmatrix(i,i) = d(i);
    }
    return Dmatrix;
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

    TMatrixDSym A(dimM+dimX+dimC);

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

    // Unknown parameters (3)
    Double_t BVz = X(dimM);
    Double_t ptau1z = X(dimM+1);
    Double_t ptau2z = X(dimM+2);

    // BV must lie in K+ trajectory
    Double_t BVx = RP.x() + (pK.x()/pK.z())*(BVz - RP.z());
    Double_t BVy = RP.y() + (pK.y()/pK.z())*(BVz - RP.z()); 
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz );
    // ptau1 must point back to the BV
    Double_t ptau1x = ((DV1.x() - BV.x())/(DV1.z() - BV.z()))*ptau1z;
    Double_t ptau1y = ((DV1.y() - BV.y())/(DV1.z() - BV.z()))*ptau1z;
    ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );
    // ptau2 must point back to the BV
    Double_t ptau2x = ((DV2.x() - BV.x())/(DV2.z() - BV.z()))*ptau2z;
    Double_t ptau2y = ((DV2.y() - BV.y())/(DV2.z() - BV.z()))*ptau2z;
    ROOT::Math::XYZVector ptau2( ptau2x, ptau2y, ptau2z );
    // Tau mass constraints
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    // 3-momentum conservation in DV1 and in DV2
    ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
    ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
    // Neutrino mass constraints
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );
    // 4-momentum conservation in BV
    ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
    Double_t EB = Etau1 + Etau2 + EK;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    // Lagrange multipliers (20)
    TVectorD l(dimC);
    for(int i = 0; i < dimC;i++)
    {
        l(i) = X(dimM+dimX+i);
    }
    
    ///////////////////////////// xm ///////////////////////////
    // PVx, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if(i != 21)
        {
            A(0,i) = 2*W(0,i);
        }
    }
    A(0,21) = 2*W(0,21) + l(0);
    A(0,23) = l(0);
    A(0,24) = l(0);
    A(0,25) = pB.z();

    // PVy, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if(i != 21)
        {
            A(1,i) = 2*W(1,i);
        }
    }
    A(1,21) = 2*W(1,21) + l(1);
    A(1,23) = l(1);
    A(1,24) = l(1);
    A(1,26) = pB.z();

    // PVz, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if( (i!=3) && (i!=4) && (i!=5) && (i!=10) && (i!=11) && (i!=12) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(2,i) = 2*W(2,i);
        }
    }
    A(2,3) = 2*W(2,3) + (l(0)*ptau1.z())/(BV.z() - DV1.z());
    A(2,4) = 2*W(2,4) + (l(1)*ptau1.z())/(BV.z() - DV1.z());
    A(2,5) = 2*W(2,5) - ((l(0)*ptau1.x() + l(1)*ptau1.y())/(BV.z() - DV1.z()));
    A(2,10) = 2*W(2,10) + (l(0)*ptau2.z())/(BV.z() - DV2.z());
    A(2,11) = 2*W(2,11) + (l(1)*ptau2.z())/(BV.z() - DV2.z());
    A(2,12) = 2*W(2,12) - ((l(0)*ptau2.x() + l(1)*ptau2.y())/(BV.z() - DV2.z()));
    A(2,17) = 2*W(2,17) + l(0)*(ptau1.z()/(-BV.z() + DV1.z()) + ptau2.z()/(-BV.z() + DV2.z()));
    A(2,18) = 2*W(2,18) + l(1)*(ptau1.z()/(-BV.z() + DV1.z()) + ptau2.z()/(-BV.z() + DV2.z()));
    A(2,19) = 2*W(2,19) + l(0)*(-1. + (ptau1.z()*(BV.z() - RP.z()))/((-BV.z() + DV1.z())*pK.z()) + (ptau2.z()*(BV.z() - RP.z()))/((-BV.z() + DV2.z())*pK.z()));
    A(2,20) = 2*W(2,20) + l(1)*(-1. + (ptau1.z()*(BV.z() - RP.z()))/((-BV.z() + DV1.z())*pK.z()) + (ptau2.z()*(BV.z() - RP.z()))/((-BV.z() + DV2.z())*pK.z()));
    A(2,21) = 2*W(2,21) + ((l(0)*pK.x() + l(1)*pK.y())*(-(DV2.z()*ptau1.z()) - DV1.z()*ptau2.z() + BV.z()*(ptau1.z() + ptau2.z()))*(BV.z() - RP.z()))/((BV.z() - DV1.z())*(BV.z() - DV2.z())*pow(pK.z(),2));
    A(2,22) = (DV2.z()*(-(l(0)*pK.z()*ptau1.x()) - l(1)*pK.z()*ptau1.y() + l(0)*pK.x()*ptau1.z() + l(1)*pK.y()*ptau1.z()) + DV1.z()*(-(l(0)*pK.z()*ptau2.x()) - l(1)*pK.z()*ptau2.y() + l(0)*pK.x()*ptau2.z() + l(1)*pK.y()*ptau2.z()) + BV.z()*(l(0)*pK.z()*(ptau1.x() + ptau2.x()) + l(1)*pK.z()*(ptau1.y() + ptau2.y()) - l(0)*pK.x()*(ptau1.z() + ptau2.z()) - l(1)*pK.y()*(ptau1.z() + ptau2.z())))/((BV.z() - DV1.z())*(BV.z() - DV2.z())*pK.z());
    A(2,23) = ( (-BV.x() + DV1.x())*l(0) + (-BV.y() + DV1.y())*l(1))/(BV.z() - DV1.z());
    A(2,24) = ( (-BV.x() + DV2.x())*l(0) + (-BV.y() + DV2.y())*l(1))/(BV.z() - DV2.z());
    A(2,25) = -pB.x();
    A(2,26) = -pB.y();

    // DV1x, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(3,i) = 2*W(3,i);
        }
    }
    A(3,3) = 2*W(3,3) + (l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.x(),2) - pow(Enu1,3)*pow(ptau1.x(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(3,4) = 2*W(3,4) + (l(2)*((pnu1.x()*pnu1.y())/pow(Enu1,3) - (ptau1.x()*ptau1.y())/pow(Etau1,3))*pow(ptau1.z(),2))/pow(BV.z() - DV1.z(),2);
    A(3,5) = 2*W(3,5) + (ptau1.z()*(-(BV.z()*pow(Enu1,3)*pow(Etau1,3)*l(0)) - 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*ptau1.x() - pow(Etau1,3)*l(2)*pow(pnu1.x(),2)*ptau1.x() + pow(Enu1,3)*l(2)*pow(ptau1.x(),3) + pow(Enu1,2)*pow(Etau1,3)*l(2)*(pnu1.x() + ptau1.x()) - pow(Etau1,3)*l(2)*pnu1.x()*pnu1.y()*ptau1.y() + pow(Enu1,3)*l(2)*ptau1.x()*pow(ptau1.y(),2) + pow(Enu1,3)*pow(Etau1,3)*l(0)*PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(3,6) = 2*W(3,6) - ((l(2)*(pow(Enu1,2) - pow(pnu1.x(),2))*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3)));
    A(3,7) = 2*W(3,7) + (l(2)*pnu1.x()*pnu1.y()*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3));
    A(3,8) = 2*W(3,8) + (l(2)*pnu1.x()*pnu1.z()*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3));
    A(3,17) = 2*W(3,17) + (l(2)*(-(pow(Enu1,3)*pow(Etau1,2)) + pow(Enu1,2)*pow(Etau1,3) - pow(Etau1,3)*pow(pnu1.x(),2) + pow(Enu1,3)*pow(ptau1.x(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(3,18) = 2*W(3,18) + (l(2)*(-((pnu1.x()*pnu1.y())/pow(Enu1,3)) + (ptau1.x()*ptau1.y())/pow(Etau1,3))*pow(ptau1.z(),2))/pow(BV.z() - DV1.z(),2);   
    A(3,19) = 2*W(3,19) + (l(2)*(-(pow(Enu1,3)*pow(Etau1,2)) + pow(Enu1,2)*pow(Etau1,3) - pow(Etau1,3)pow(pnu1.x(),2) + pow(Enu1,3)*pow(ptau1.x(),2))*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z());
    A(3,20) = 2*W(3,20) + (l(2)*(-(pow(Etau1,3)*pnu1.x()*pnu1.y()) + pow(Enu1,3)*ptau1.x()*ptau1.y())*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z());
    A(3,21) = 2*W(3,21) + (l(2)*(pow(Enu1,3)*pow(Etau1,2)*pK.x() - pow(Enu1,2)*pow(Etau1,3)*pK.x() + pow(Etau1,3)*pK.x()*pow(pnu1.x(),2) + pow(Etau1,3)*pK.y()*pnu1.x()*pnu1.y() - pow(Enu1,3)*pK.x()*pow(ptau1.x(),2) - pow(Enu1,3)*pK.y()*ptau1.x()*ptau1.y())*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pow(pK.z(),2));
    A(3,22) = -((ptau1.z()*(-(DV1.z()*pow(Enu1,3)*pow(Etau1,3)*l(0)*pK.z()) - 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*pK.z()*ptau1.x() - pow(Etau1,3)*l(2)*pK.z()*pow(pnu1.x(),2)*ptau1.x() + pow(Enu1,3)*l(2)*pK.z()*pow(ptau1.x(),3) - pow(Etau1,3)*l(2)*pK.z()*pnu1.x()*pnu1.y()*ptau1.y() + pow(Enu1,3)*l(2)*pK.z()*ptau1.x()*pow(ptau1.y(),2) + pow(Enu1,3)*pow(Etau1,2)*l(2)*pK.x()*ptau1.z() + pow(Etau1,3)*l(2)*pK.x()*pow(pnu1.x(),2)*ptau1.z() + pow(Etau1,3)*l(2)*pK.y()*pnu1.x()*pnu1.y()*ptau1.z() - pow(Enu1,3)*l(2)*pK.x()*pow(ptau1.x(),2)*ptau1.z() - pow(Enu1,3)*l(2)*pK.y()*ptau1.x()*ptau1.y()*ptau1.z() + pow(Enu1,2)*pow(Etau1,3)*l(2)*(pK.z()*(pnu1.x() + ptau1.x()) - pK.x()*ptau1.z()) + pow(Enu1,3)*pow(Etau1,3)*l(0)*pK.z()*PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z()));
    A(3,23) = (((BV.z() - DV1.z())*l(2)*pnu1.x())/sqrt(pow(Enu1,2)) + ((BV.z() - DV1.z())*l(2)*ptau1.x())/Enu1 - (2*(BV.z() - DV1.z())*l(2)*ptau1.x())/Etau1 - (l(2)*pnu1.x()*(BV.x()*pnu1.x() - DV1.x()*pnu1.x() + BV.y()*pnu1.y() - DV1.y()*pnu1.y() + BV.z()*pnu1.z() - DV1.z()*pnu1.z())*ptau1.z())/pow(Enu1,3) + (l(2)*ptau1.x()*ptau1.z()*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + BV.y()*ptau1.y() - DV1.y()*ptau1.y() + BV.z()*ptau1.z() - DV1.z()*ptau1.z()))/pow(Etau1,3) - (BV.z() - DV1.z())*l(0)*(BV.z() - PV.z()))/pow(BV.z() - DV1.z(),2); 
    A(3,25) =(ptau1.z()*(BV.z() - PV.z()))/(-BV.z() + DV1.z());
    A(3,27) = ((-(pnu1.x()/Enu1) + ptau1.x()/Etau1)*ptau1.z())/(-BV.z() + DV1.z());

    // DV1y, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(4,i) = 2*W(4,i);
        }
    }
    A(4,4) = 2*W(4,4) + (l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + po)w(Etau1,3)*pow(pnu1.y(),2) - pow(Enu1,3)*pow(ptau1.y(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(4,5) = 2*W(4,5) + (ptau1.z()*(-(BV.z()*pow(Enu1,3)*pow(Etau1,3)*l(1)) - pow(Etau1,3)*l(2)*pnu1.x()*pnu1.y()*ptau1.x() - 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*ptau1.y() - pow(Etau1,3)*l(2)*pow(pnu1.y(),2)*ptau1.y() + pow(Enu1,3)*l(2)*pow(ptau1.x(),2)*ptau1.y() + pow(Enu1,3)*l(2)*pow(ptau1.y(),3) + pow(Enu1,2)*pow(Etau1,3)*l(2)*(pnu1.y() + ptau1.y()) + pow(Enu1,3)*pow(Etau1,3)*l(1)*PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(4,6) = 2*W(4,6) + (l(2)*pnu1.x()*pnu1.y()*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3));
    A(4,7) = 2*W(4,7) - ((l(2)*(pow(Enu1,2) - pow(pnu1.y(),2))*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3)));
    A(4,8) = 2*W(4,8) + (l(2)*pnu1.y()*pnu1.z()*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3));
    A(4,17) = 2*W(4,17) + (l(2)*(-((pnu1.x()*pnu1.y())/pow(Enu1,3)) + (ptau1.x()*ptau1.y())/pow(Etau1,3))*pow(ptau1.z(),2))/pow(BV.z() - DV1.z(),2);
    A(4,18) = 2*W(4,18) + (l(2)*(-(pow(Enu1,3)*pow(Etau1,2)) + pow(Enu1,2)*pow(Etau1,3) - pow(Etau1,3)*pow(pnu1.y(),2) + pow(Enu1,3)*pow(ptau1.y(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(4,19) = 2*W(4,19) + (l(2)*(-(pow(Etau1,3)*pnu1.x()*pnu1.y()) + pow(Enu1,3)*ptau1.x()*ptau1.y())*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z());
    A(4,20) = 2*W(4,20) + (l(2)*(-(pow(Enu1,3)*pow(Etau1,2)) + pow(Enu1,2)*pow(Etau1,3) - pow(Etau1,3)*pow(pnu1.y(),2) + pow(Enu1,3)*pow(ptau1.y(),2))*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z());
    A(4,21) = 2*W(4,21) + (l(2)*(pow(Enu1,3)*pow(Etau1,2)*pK.y() - pow(Enu1,2)*pow(Etau1,3)*pK.y() + pow(Etau1,3)*pK.x()*pnu1.x()*pnu1.y() + pow(Etau1,3)*pK.y()*pow(pnu1.y(),2) - pow(Enu1,3)*pK.x()*ptau1.x()*ptau1.y() - pow(Enu1,3)*pK.y()*pow(ptau1.y(),2))*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pow(pK.z(),2));
    A(4,22) =  -((ptau1.z()*(-(DV1.z()*pow(Enu1,3)*pow(Etau1,3)*l(1)*pK.z()) - pow(Etau1,3)*l(2)*pK.z()*pnu1.x()*pnu1.y()*ptau1.x() - 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*pK.z()*ptau1.y() - pow(Etau1,3)*l(2)*pK.z()*pow(pnu1.y(),2)*ptau1.y() + pow(Enu1,3)*l(2)*pK.z()*pow(ptau1.x(),2)*ptau1.y() + pow(Enu1,3)*l(2)*pK.z()*pow(ptau1.y(),3) + pow(Enu1,3)*pow(Etau1,2)*l(2)*pK.y()*ptau1.z() + pow(Etau1,3)*l(2)*pK.x()*pnu1.x()*pnu1.y()*ptau1.z() + pow(Etau1,3)*l(2)*pK.y()*pow(pnu1.y(),2)*ptau1.z() - pow(Enu1,3)*l(2)*pK.x()*ptau1.x()*ptau1.y()*ptau1.z() - pow(Enu1,3)*l(2)*pK.y()*pow(ptau1.y(),2)*ptau1.z() + pow(Enu1,2)*pow(Etau1,3)*l(2)*(pK.z()*(pnu1.y() + ptau1.y()) - pK.y()*ptau1.z()) + pow(Enu1,3)*pow(Etau1,3)*l(1)*pK.z()*PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z()));
    A(4,23) = (((BV.z() - DV1.z())*l(2)*pnu1.y())/Enu1 + ((BV.z() - DV1.z())*l(2)*ptau1.y())/Enu1 - (2*(BV.z() - DV1.z())*l(2)*ptau1.y())/Etau1 - (l(2)*pnu1.y()*(BV.x()*pnu1.x() - DV1.x()*pnu1.x() + BV.y()*pnu1.y() - DV1.y()*pnu1.y() + BV.z()*pnu1.z() - DV1.z()*pnu1.z())*ptau1.z())/pow(Enu1,3) + (l(2)*ptau1.y()*ptau1.z()*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + BV.y()*ptau1.y() - DV1.y()*ptau1.y() + BV.z()*ptau1.z() - DV1.z()*ptau1.z()))/pow(Etau1,3) - (BV.z() - DV1.z())*l(1)*(BV.z() - PV.z()))/pow(BV.z() - DV1.z(),2);
    A(4,26) = (ptau1.z()*(BV.z() - PV.z()))/(-BV.z() + DV1.z());
    A(4,27) = ((-(pnu1.y()/pow(Enu1,3)) + ptau1.y()/pow(Etau1,3))*ptau1.z())/(-BV.z() + DV1.z());

    // DV1z, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(5,i) = 2*W(5,i);
        }
    }
    A(5,5) = 2*W(5,5) + ((l(2)*pow(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y(),2))/pow(Enu1,3) + (3*l(2)*(pow(ptau1.x(),2) + pow(ptau1.y(),2)))/Etau1 - (l(2)*pow(pow(ptau1.x(),2) + pow(ptau1.y(),2),2))/pow(Etau1,3) - (l(2)*(2*pnu1.x()*ptau1.x() + pow(ptau1.x(),2) + ptau1.y()*(2*pnu1.y() + ptau1.y())))/Enu1 + 2*l(0)*ptau1.x()*(BV.z() - PV.z()) + 2*l(1)*ptau1.y()*(BV.z() - PV.z()))/pow(BV.z() - DV1.z(),2);
    A(5,6) = 2*W(5,6) + (l(2)*(pow(Enu1,2)*ptau1.x() - pnu1.x()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y())))/((BV.z() - DV1.z())*pow(Enu1,3));
    A(5,7) = 2*W(5,7) + (l(2)*(-(pnu1.x()*pnu1.y()*ptau1.x()) + (pow(Enu1,2) - pow(pnu1.y(),2))*ptau1.y()))/((BV.z() - DV1.z())*pow(Enu1,3));
    A(5,8) = 2*W(5,8) - ((l(2)*pnu1.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()))/((BV.z() - DV1.z())*pow(Enu1,3)));
    A(5,17) = 2*W(5,17) + (ptau1.z()*(BV.z()*pow(Enu1,3)*pow(Etau1,3)*l(0) + 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*ptau1.x() + pow(Etau1,3)*l(2)*pow(pnu1.x(),2)*ptau1.x() - pow(Enu1,3)*l(2)*pow(ptau1.x(),3) - pow(Enu1,2)*pow(Etau1,3)*l(2)*(pnu1.x() + ptau1.x()) + pow(Etau1,3)*l(2)*pnu1.x()*pnu1.y()*ptau1.y() - pow(Enu1,3)*l(2)*ptau1.x()*pow(ptau1.y(),2) - pow(Enu1,3)*pow(Etau1,3)*l(0)*PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(5,18) = 2*W(5,18) + (ptau1.z()*(BV.z()*pow(Enu1,3)*pow(Etau1,3)*l(1) + pow(Etau1,3)*l(2)*pnu1.x()*pnu1.y()*ptau1.x() + 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*ptau1.y() + pow(Etau1,3)*l(2)*pow(pnu1.y(),2)*ptau1.y() - pow(Enu1,3)*l(2)*pow(ptau1.x(),2)*ptau1.y() - pow(Enu1,3)*l(2)*pow(ptau1.y(),3) - pow(Enu1,2)*pow(Etau1,3)*l(2)*(pnu1.y() + ptau1.y()) - pow(Enu1,3)*pow(Etau1,3)*l(1)*PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3));
    A(5,19) = 2*W(5,19) + (ptau1.z()*(BV.z()*pow(Enu1,3)*pow(Etau1,3)*l(0) + 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*ptau1.x() + pow(Etau1,3)*l(2)*pow(pnu1x,2)*ptau1.x() - pow(Enu1,3)*l(2)*pow(ptau1.x(),3) - pow(Enu1,2)*pow(Etau1,3)*l(2)*(pnu1.x() + ptau1.x()) + pow(Etau1,3)*l(2)*pnu1.x()*pnu1.y()*ptau1.y() - pow(Enu1,3)*l(2)*ptau1.x()*pow(ptau1.y(),2) - pow(Enu1,3)*pow(Etau1,3)*l(0)*PV.z())*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z());
    A(5,20) = 2*W(5,20) + (ptau1.z()*(BV.z()*pow(Enu1,3)*pow(Etau1,3)*l(1) + pow(Etau1,3)*l(2)*pnu1.x()*pnu1.y()*ptau1.x() + 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*ptau1.y() + pow(Etau1,3)*l(2)*pow(pnu1.y(),2)*ptau1.y() - pow(Enu1,3)*l(2)*pow(ptau1.x(),2)*ptau1.y() - pow(Enu1,3)*l(2)*pow(ptau1.y(),3) - pow(Enu1,2)*pow(Etau1,3)*l(2)*(pnu1.y() + ptau1.y()) - pow(Enu1,3)*pow(Etau1,3)*l(1)*PV.z())*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z());
    A(5,21) = 2*W(5,21) - ((ptau1.z()*(BV.z()*pow(Enu1,3)*pow(Etau1,3)*(l(0)*pK.x() + l(1)*pK.y()) + 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*pK.x()*ptau1.x() + pow(Etau1,3)*l(2)*pK.x()*pow(pnu1.x(),2)*ptau1.x() + pow(Etau1,3)*l(2)*pK.y()*pnu1.x()*pnu1.y()*ptau1.x() - pow(Enu1,3)*l(2)*pK.x()*pow(ptau1.x(),3) + 2*pow(Enu1,3)*pow(Etau1,2)*l(2)*pK.y()*ptau1.y() + pow(Etau1,3)*l(2)*pK.x()*pnu1.x()*pnu1.y()*ptau1.y() + pow(Etau1,3)*l(2)*pK.y()*pow(pnu1.y(),2)*ptau1.y() - pow(Enu1,3)*l(2)*pK.y()*pow(ptau1.x(),2)*ptau1.y() - pow(Enu1,3)*l(2)*pK.x()*ptau1.x()*pow(ptau1.y(),2) - pow(Enu1,3)*l(2)*pK.y()*pow(ptau1.y(),3) - pow(Enu1,2)*pow(Etau1,3)*l(2)*(pK.x()*(pnu1.x() + ptau1.x()) + pK.y()*(pnu1.y() + ptau1.y())) - pow(Enu1,3)*pow(Etau1,3)*l(0)*pK.x()*PV.z() - pow(Enu1,3)*pow(Etau1,3)*l(1)*pK.y()*PV.z())*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pow(pK.z(),2)));
    A(5,22) = ((BV.z() - DV1.z())*l(0)*ptau1.x() + (BV.z() - DV1.z())*l(2)*ptau1.y() - (l(2)*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y())*(pK.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()) - (pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()))/(pow(Enu1,3)*pK.z()) + (l(2)*(pow(ptau1.x(),2) + pow(ptau1.y(),2))*(pK.z()*(pow(ptau1.x(),2) + pow(ptau1.y(),2)) - (pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z()))/(pow(Etau1,3)*pK.z()) + (l(2)*(-3*pK.z()*(pow(ptau1.x(),2) + pow(ptau1.y(),2)) + 2*(pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z()))/(Etau1*pK.z()) + (l(2)*(pK.z()*(2*pnu1.x()*ptau1.x() + pow(ptau1.x(),2) + 2*pnu1.y()*ptau1.y() + pow(ptau1.y(),2)) - (pK.x()*(pnu1.x() + ptau1.x()) + pK.y()*(pnu1.y() + ptau1.y()))*ptau1.z()))/(Enu1*pK.z()) - 2*l(0)*ptau1.x()*(BV.z() - PV.z()) - 2*l(1)*ptau1.y()*(BV.z() - PV.z()) + (l(0)*pK.x()*ptau1.z()*(BV.z() - PV.z()))/pK.z() + (l(1)*pK.y()*ptau1.z()*(BV.z() - PV.z()))/pK.z())/pow(BV.z() - DV1.z(),2);
    A(5,23) = ((2*l(2)*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + (BV.y() - DV1.y())*ptau1.y()))/Etau1 + (l(2)*(BV.x()*pnu1.x() - DV1.x()*pnu1.x() + BV.y()*pnu1.y() - DV1.y()*pnu1.y() + BV.z()*pnu1.z() - DV1.z()*pnu1.z())*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()))/pow(Enu1,3) - (l(2)*(BV.x()*(pnu1.x() + ptau1.x()) - DV1.x()*(pnu1.x() + ptau1.x()) + (BV.y() - DV1.y())*(pnu1.y() + ptau1.y())))/Enu1 - (l(2)*(pow(ptau1.x(),2) + pow(ptau1.y(),2))*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + BV.y()*ptau1.y() - DV1.y()*ptau1.y() + BV.z()*ptau1.z() - DV1.z()*ptau1.z()))/pow(Etau1,3) + (BV.x() - DV1.x())*l(0)*(BV.z() - PV.z()) + (BV.y() - DV1.y())*l(1)*(BV.z() - PV.z()))/pow(BV.z() - DV1.z(),2);
    A(5,25) = (ptau1.x()*(BV.z() - PV.z()))/(BV.z() - DV1.z());
    A(5,26) = (ptau1.y()*(BV.z() - PV.z()))/(BV.z() - DV1.z());
    A(5,27) = -(((pnu1.x()*ptau1.x())/Enu1 - pow(ptau1.x(),2)/Etau1 + (pnu1.y()*ptau1.y())/Enu1 - pow(ptau1.y(),2)/Etau1)/(BV.z() - DV1.z()));

    // p3pi1x, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(6,i) = 2*W(6,i);
        }
    }
    A(6,6) = 2*W(6,6) + (l(2)*(-pow(Enu1,2) + pow(pnu1.x(),2)))/pow(Enu1,3);
    A(6,7) = 2*W(6,7) + (l(2)*pnu1.x()*pnu1.y())/pow(Enu1,3);
    A(6,8) = 2*W(6,8) + (l(2)*pnu1.x()*pnu1.z())/pow(Enu1,3);
    A(6,17) = 2*W(6,17) + (l(2)*(pow(Enu1,2) - pow(pnu1.x(),2))*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3));
    A(6,18) = 2*W(6,18) + (l(2)*pnu1.x()*pnu1.y()*ptau1.z())/((-BV.z() + DV1.z())*pow(Enu1,3));
    A(6,19) = 2*W(6,19) + (l(2)*(pow(Enu1,2) - pow(pnu1.x(),2))*ptau1.z()*(BV.z() - RP.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pK.z());
    A(6,20) = 2*W(6,20) + (l(2)*pnu1.x()*pnu1.y()*ptau1.z()*(BV.z() - RP.z()))/((-BV.z() + DV1.z())*pow(Enu1,3)*pK.z());
    A(6,21) = 2*W(6,21) - ((l(2)*(pow(Enu1,2)*pK.x() - pnu1.x()*(pK.x()*pnu1.x() + pK.y()*pnu1.y()))*ptau1.z()*(BV.z() - RP.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pow(pK.z(),2)));
    A(6,22) = -((l(2)*(-(pK.z()*pnu1.x()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y())) + pnu1.x()*(pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z() + pow(Enu1,2)*(pK.z()*ptau1.x() - pK.x()*ptau1.z())))/((BV.z() - DV1.z())*pow(Enu1,3)*pK.z()));
    A(6,23) = (l(2)*(BV.x()*(pow(Enu1,2) - pow(pnu1.x(),2)) + DV1.x()*(-pow(Enu1,2) + pow(pnu1.x(),2)) + pnu1.x()*(-(BV.y()*pnu1.y()) + DV1.y()*pnu1.y() - BV.z()*pnu1.z() + DV1.z()*pnu1.z())))/((BV.z() - DV1.z())*pow(Enu1,3));
    A(6,27) = pnu1.x()/Enu1;

    // p3pi1y, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(7,i) = 2*W(7,i);
        }
    }
    A(7,7) = 2*W(7,7) + (l(2)*(-pow(Enu1,2) + pow(pnu1.y(),2)))/pow(Enu1,3);
    A(7,8) = 2*W(7,8) + (l(2)*pnu1.y()*pnu1.z())/pow(Enu1,3);
    A(7,17) = 2*W(7,17) + (l(2)*pnu1.x()*pnu1.y()*ptau1.z())/((-BV.z() + DV1.z())*pow(Enu1,3));
    A(7,18) = 2*W(7,18) + (l(2)*(pow(Enu1,2) - pow(pnu1.y(),2))*ptau1.z())/((BV.z() - DV1.z())*pow(Enu1,3));
    A(7,19) = 2*W(7,19) + (l(2)*pnu1.x()*pnu1.y()*ptau1.z()*(BV.z() - RP.z()))/((-BV.z() + DV1.z())*pow(Enu1,3)*pK.z());
    A(7,20) = 2*W(7,20) + (l(2)*(pow(Enu1,2) - pow(pnu1.y(),2))*ptau1.z()*(BV.z() - RP.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pK.z());
    A(7,21) = 2*W(7,21) - ((l(2)*(pow(Enu1,2)*pK.y() - pnu1.y()*(pK.x()*pnu1.x() + pK.y()*pnu1.y()))*ptau1.z()*(BV.z(9) - RP.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pow(pK.z(),2)));
    A(7,22) = (l(2)*(pK.z()*pnu1.x()*pnu1.y()*ptau1.x() + pK.z*()(-pow(Enu1,2) + pow(pnu1.y(),2))*ptau1.y() + pow(Enu1,2)*pK.y()*ptau1.z() - pnu1.y()*(pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pK.z());
    A(7,23) = (l(2)*(BV.y()*(pow(Enu1,2) - pow(pnu1.y(),2)) + DV1.y()*(-pow(Enu1,2) + pow(pnu1.y(),2)) + pnu1.y()*(-(BV.x()*pnu1.x()) + DV1.x()*pnu1.x() - BV.z()*pnu1.z() + DV1.z()*pnu1.z())))/((BV.z() - DV1.z())*pow(Enu1,3));
    A(7,27) = pnu1.y()/Enu1;

    // p3pi1z, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(8,i) = 2*W(8,i);
        }
    }   
    A(8,8) = 2*W(8,8) + (l(2)*(-pow(Enu1,2) + pow(pnu1.z(),2)))/pow(Enu1,3);
    A(8,17) = 2*W(8,17) + (l(2)*pnu1.x()*pnu1.z()*ptau1.z())/((-BV.z() + DV1.z())*pow(Enu1,3));
    A(8,18) = 2*W(8,18) + (l(2)*pnu1.y()*pnu1.z()*ptau1.z())/((-BV.z() + DV1.z())*pow(Enu1,3));
    A(8,19) = 2*W(8,19) + (l(2)*pnu1.x()*pnu1.z()*ptau1.z()*(BV.z() - RP.z()))/((-BV.z() + DV1.z())*pow(Enu1,3)*pK.z());
    A(8,20) = 2*W(8,20) + (l(2)*pnu1.y()*pnu1.z()*ptau1.z()*(BV.z() - RP.z()))/((-BV.z() + DV1.z())*pow(Enu1,3)*pK.z());
    A(8,21) = 2*W(8,21) + (l(2)*(pK.x()*pnu1.x() + pK.y()*pnu1.y())*pnu1.z()*ptau1.z()*(BV.z() - RP.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pow(pK.z(),2));
    A(8,22) = (l(2)*pnu1.z()*(pK.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()) - (pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()))/((BV.z() - DV1.z())*pow(Enu1,3)*pK.z());
    A(8,23) = (l(2)*((-(BV.x()*pnu1.x()) + DV1.x()*pnu1.x() - BV.y()*pnu1.y() + DV1.y()*pnu1.y())*pnu1.z() + BV.z()*(pow(Enu1,2) - pow(pnu1.z(),2)) + DV1.z()*(-pow(Enu1,2) + pow(pnu1.z(),2))))/((BV.z() - DV1.z())*pow(Enu1,3));
    A(8,27) = pnu1.z()/Enu1;

    // E3pi1, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        A(9,i) = 2*W(9,i);
    } 
    A(9,27) = -1;

    // DV2x, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(10,i) = 2*W(10,i);
        }
    }
    A(10,10) = 2*W(10,10) + (l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.x(),2) - pow(Enu2,3)*pow(ptau2.x(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(10,11) = 2*W(10,11) + (l(3)*((pnu2.x()*pnu2.y())/pow(Enu2,3) - (ptau2.x()*ptau2.y())/pow(Etau2,3))*pow(ptau2.z(),2))/pow(BV.z() - DV2.z(),2);
    A(10,12) = 2*W(10,12) + (ptau2.z()*(-(BV.z()*pow(Enu2,3)*pow(Etau2,3)*l(0)) - 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*ptau2.x() - pow(Etau2,3)*l(3)*pow(pnu2.x(),2)*ptau2.x() + pow(Enu2,3)*l(3)*pow(ptau2.x(),3) + pow(Enu2,2)*pow(Etau2,3)*l(3)*(pnu2.x() + ptau2.x()) - pow(Etau2,3)*l(3)*pnu2.x()*pnu2.y()*ptau2.y() + pow(Enu2,3)*l(3)*ptau2.x()*pow(ptau2.y(),2) + pow(Enu2,3)*pow(Etau2,3)*l(0)*PV.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(10,13) = 2*W(10,13) - ((l(3)*(pow(Enu2,2) - pow(pnu2.x(),2))*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3)));
    A(10,14) = 2*W(10,14) + (l(3)*pnu2.x()*pnu2.y()*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3));
    A(10,15) = 2*W(10,15) + (l(3)*pnu2.x()*pnu2.z()*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3));
    A(10,17) = 2*W(10,17) + (l(3)*(-(pow(Enu2,3)*pow(Etau2,2)) + pow(Enu2,2)*pow(Etau2,3) - pow(Etau2,3)*pow(pnu2.x(),2) + pow(Enu2,3)*pow(ptau2.x(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(10,18) = 2*W(10,18) + (l(3)*(-((pnu2.x()*pnu2.y())/pow(Enu2,3)) + (ptau2.x()*ptau2.y())/pow(Etau2,3))*pow(ptau2.z(),2))/pow(BV.z() - DV2.z(),2);
    A(10,19) = 2*W(10,19) + (l(3)*(-(pow(Enu2,3)*pow(Etau2,2)) + pow(Enu2,2)*pow(Etau2,3) - pow(Etau2,3)*pow(pnu2.x(),2) + pow(Enu2,3)*pow(ptau2.x(),2))*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z());
    A(10,20) = 2*W(10,20) + (l(3)*(-(pow(Etau2,3)*pnu2.x()*pnu2.y()) + pow(Enu2,3)*ptau2.x()*ptau2.y())*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z());
    A(10,21) = 2*W(10,21) + (l(3)*(pow(Enu2,3)*pow(Etau2,2)*pK.x() - pow(Enu2,2)*pow(Etau2,3)*pK.x() + pow(Etau2,3)*pK.x()*pow(pnu2.x(),2) + pow(Etau2,3)*pK.y()*pnu2.x()*pnu2.y() - pow(Enu2,3)*pK.x()*pow(ptau2.x(),2) - pow(Enu2,3)*pK.y()*ptau2.x()*ptau2.y())*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pow(pK.z(),2));
    A(10,22) = -((ptau2.z()*(-(DV2.z()*pow(Enu2,3)*pow(Etau2,3)*l(0)*pK.z()) - 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*pK.z()*ptau2.x() - pow(Etau2,3)*l(3)*pK.z()*pow(pnu2.x(),2)*ptau2.x() + pow(Enu2,3)*l(3)*pK.z()*pow(ptau2.x(),3) - pow(Etau2,3)*l(3)*pK.z()*pnu2.x()*pnu2.y()*ptau2.y() + pow(Enu2,3)*l(3)*pK.z()*ptau2.x()*pow(ptau2.y(),2) + pow(Enu2,3)*pow(Etau2,2)*l(3)*pK.x()*ptau2.z() + pow(Etau2,3)*l(3)*pK.x()*pow(pnu2.x(),2)*ptau2.z() + pow(Etau2,3)*l(3)*pK.y()*pnu2.x()*pnu2.y()*ptau2.z() - pow(Enu2,3)*l(3)*pK.x()*pow(ptau2.x(),2)*ptau2.z() - pow(Enu2,3)*l(3)*pK.y()*ptau2.x()*ptau2.y()*ptau2.z() + pow(Enu2,2)*pow(Etau2,3)*l(3)*(pK.z()*(pnu2.x() + ptau2.x()) - pK.x()*ptau2.z()) + pow(Enu2,3)*pow(Etau2,3)*l(0)*pK.z()*PV.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z()));
    A(10,24) = (((BV.z() - DV2.z())*l(3)*pnu2.x())/Enu2 + ((BV.z() - DV2.z())*l(3)*ptau2.x())/Enu2 - (2*(BV.z() - DV2.z())*l(3)*ptau2.x())/Etau2 - (l(3)*pnu2.x()*(BV.x()*pnu2.x() - DV2.x()*pnu2.x() + BV.y()*pnu2.y() - DV2.y()*pnu2.y() + BV.z()*pnu2.z() - DV2.z()*pnu2.z())*ptau2.z())/pow(Enu2,3) + (l(3)*ptau2.x()*ptau2.z()*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + BV.y()*ptau2.y() - DV2.y()*ptau2.y() + BV.z()*ptau2.z() - DV2.z()*ptau2.z()))/pow(Etau2,3) - (BV.z() - DV2.z())*l(0)*(BV.z() - PV.z()))/pow(BV.z() - DV2.z(),2);
    A(10,25) = (ptau2.z()*(BV.z() - PV.z()))/(-BV.z() + DV2.z());
    A(10,28) = ((-(pnu2.x()/Enu2) + ptau2.x()/Etau2)*ptau2.z())/(-BV.z() + DV2.z());

    // DV2y, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=10) && (i!=11) && (i!=12) & (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(11,i) = 2*W(11,i);
        }
    }
    A(11,11) = 2*W(11,11) + (l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.y(),2) - pow(Enu2,3)*pow(ptau2.y(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(11,12) = 2*W(11,12) + (ptau2.z()*(-(BV.z()*pow(Enu2,3)*pow(Etau2,3)*l(1)) - pow(Etau2,3)*l(3)*pnu2.x()*pnu2.y()*ptau2.x() - 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*ptau2.y() - pow(Etau2,3)*l(3)*pow(pnu2y,2)*ptau2.y() + pow(Enu2,3)*l(3)*pow(ptau2.x(),2)*ptau2.y() + pow(Enu2,3)*l(3)*pow(ptau2.y(),3) + pow(Enu2,2)*pow(Etau2,3)*l(3)*(pnu2.y() + ptau2.y()) + pow(Enu2,3)*pow(Etau2,3)*l(1)*PV.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(11,13) = 2*W(11,13) + (l(3)*pnu2.x()*pnu2.y()*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3));
    A(11,14) = 2*W(11,14) - ((l(3)*(pow(Enu2,2) - pow(pnu2.y(),2))*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3)));
    A(11,15) = 2*W(11,15) + (l(3)*pnu2.y()*pnu2.z()*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3));
    A(11,17) = 2*W(11,17) + (l(3)*(-((pnu2.x()*pnu2.y())/pow(Enu2,3)) + (ptau2.x()*ptau2.y())/pow(Etau2,3))*pow(ptau2.z(),2))/pow(BV.z() - DV2.z(),2);
    A(11,18) = 2*W(11,18) + (l(3)*(-(pow(Enu2,3)*pow(Etau2,2)) + pow(Enu2,2)*pow(Etau2,3) - pow(Etau2,3)*pow(pnu2.y(),2) + pow(Enu2,3)*pow(ptau2.y(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(11,19) = 2*W(11,19) + (l(3)*(-(pow(Etau2,3)*pnu2.x()*pnu2.y()) + pow(Enu2,3)*ptau2.x()*ptau2.y())*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z());
    A(11,20) = 2*W(11,20) + (l(3)*(-(pow(Enu2,3)*pow(Etau2,2)) + pow(Enu2,2)*pow(Etau2,3) - pow(Etau2,3)*pow(pnu2.y(),2) + pow(Enu2,3)*pow(ptau2.y(),2))*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z());
    A(11,21) = 2*W(11,21) + (l(3)*(pow(Enu2,3)*pow(Etau2,2)*pK.y() - pow(Enu2,2)*pow(Etau2,3)*pK.y() + pow(Etau2,3)*pK.x()*pnu2.x()*pnu2.y() + pow(Etau2,3)*pK.y()*pow(pnu2.y(),2) - pow(Enu2,3)*pK.x()*ptau2.x()*ptau2.y() - pow(Enu2,3)*pK.y()*pow(ptau2.y(),2))*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pow(pK.z(),2));
    A(11,22) = -((ptau2.z()*(-(DV2.z()*pow(Enu2,3)*pow(Etau2,3)*l(1)*pK.z()) - pow(Etau2,3)*l(3)*pK.z()*pnu2.x()*pnu2.y()*ptau2.x() - 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*pK.z()*ptau2.y() - pow(Etau2,3)*l(3)*pK.z()*pow(pnu2.y(),2)*ptau2.y() + pow(Enu2,3)*l(3)*pK.z()*pow(ptau2.x(),2)*ptau2.y() + pow(Enu2,3)*l(3)*pK.z()*pow(ptau2.y(),3) + pow(Enu2,3)*pow(Etau2,2)*l(3)*pK.y()*ptau2.z() + pow(Etau2,3)*l(3)*pK.x()*pnu2.x()*pnu2.y()*ptau2.z() + pow(Etau2,3)*l(3)*pK.y()*pow(pnu2.y(),2)*ptau2.z() - pow(Enu2,3)*l(3)*pK.x()*ptau2.x()*ptau2.y()*ptau2.z() - pow(Enu2,3)*l(3)*pK.y()*pow(ptau2.y(),2)*ptau2.z() + pow(Enu2,2)*pow(Etau2,3)*l(3)*(pK.z()*(pnu2.y() + ptau2.y()) - pK.y()*ptau2.z()) + pow(Enu2,3)*pow(Etau2,3)*l(1)*pK.z()*PV.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z()));
    A(11,24) = (((BV.z() - DV2.z())*l(3)*pnu2.y())/Enu2 + ((BV.z() - DV2.z())*l(3)*ptau2.y())/Enu2 - (2*(BV.z() - DV2.z())*l(3)*ptau2.y())/Etau2 - (l(3)*pnu2.y()*(BV.x()*pnu2.x() - DV2.x()*pnu2.x() + BV.y()*pnu2.y() - DV2.y()*pnu2.y() + BV.z()*pnu2.z() - DV2.z()*pnu2.z())*ptau2.z())/pow(Enu2,3) + (l(3)*ptau2.y()*ptau2.z()*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + BV.y()*ptau2.y() - DV2.y()*ptau2.y() + BV.z()*ptau2.z() - DV2.z()*ptau2.z()))/pow(Etau2,3) - (BV.z() - DV2.z())*l(1)*(BV.z() - PV.z()))/pow(BV.z() - DV2.z(),2);
    A(11,26) = (ptau2.z()*(BV.z() - PV.z()))/(-BV.z() + DV2.z());
    A(11,28) = ((-(pnu2.y()/Enu2) + ptau2.y()/Etau2)*ptau2.z())/(-BV.z() + DV2.z());

    // DV2z, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21))
        {
            A(12,i) = 2*W(12,i);
        }
    }
    A(12,12) = 2*W(12,12) + ((l(3)*pow(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y(),2))/pow(Enu2,3) + (3*l(3)*(pow(ptau2.x(),2) + pow(ptau2.y(),2)))/Etau2 - (l(3)*pow(pow(ptau2.x(),2) + pow(ptau2.y(),2),2))/pow(Etau2,3) - (l(3)*(2*pnu2.x()*ptau2.x() + pow(ptau2.x(),2) + ptau2.y()*(2*pnu2.y() + ptau2.y())))/Enu2 + 2*l(0)*ptau2.x()*(BV.z() - PV.z()) + 2*l(1)*ptau2.y()*(BV.z() - PV.z()))/pow(BV.z() - DV2.z(),2);
    A(12,13) = 2*W(12,13) + (l(3)*(pow(Enu2,2)*ptau2.x() - pnu2.x()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y())))/((BV.z() - DV2.z())*pow(Enu2,3));
    A(12,14) = 2*W(12,14) + (l(3)*(-(pnu2.x()*pnu2.y()*ptau2.x()) + (pow(Enu2,2) - pow(pnu2.y(),2))*ptau2.y()))/((BV.z() - DV2.z())*pow(Enu2,3));
    A(12,15) = 2*W(12,15) - ((l(3)*pnu2.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()))/((BV.z() - DV2.z())*pow(Enu2,3)));
    A(12,17) = 2*W(12,17) + (ptau2.z()*(BV.z()*pow(Enu2,3)*pow(Etau2,3)*l(0) + 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*ptau2.x() + pow(Etau2,3)*l(3)*pow(pnu2.x(),2)*ptau2.x() - pow(Enu2,3)*l(3)*pow(ptau2.x(),3) - pow(Enu2,2)*pow(Etau2,3)*l(3)*(pnu2.x() + ptau2.x()) + pow(Etau2,3)*l(3)*pnu2.x()*pnu2.y()*ptau2.y() - pow(Enu2,3)*l(3)*ptau2.x()*pow(ptau2.y(),2) - pow(Enu2,3)*pow(Etau2,3)*l(0)*PV.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(12,18) = 2*W(12,18) + (ptau2.z()*(BV.z()*pow(Enu2,3)*pow(Etau2,3)*l(1) + pow(Etau2,3)*l(3)*pnu2.x()*pnu2.y()*ptau2.x() + 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*ptau2.y() + pow(Etau2,3)*l(3)*pow(pnu2.y(),2)*ptau2.y() - pow(Enu2,3)*l(3)*pow(ptau2.x(),2)*ptau2.y() - pow(Enu2,3)*l(3)*pow(ptau2.y(),3) - pow(Enu2,2)*pow(Etau2,3)*l(3)*(pnu2.y() + ptau2.y()) - pow(Enu2,3)*pow(Etau2,3)*l(1)*PV.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(12,19) = 2*W(12,19) + (ptau2.z()*(BV.z()*pow(Enu2,3)*pow(Etau2,3)*l(0) + 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*ptau2.x() + pow(Etau2,3)*l(3)*pow(pnu2.x(),2)*ptau2.x() - pow(Enu2,3)*l(3)*pow(ptau2.x(),3) - pow(Enu2,2)*pow(Etau2,3)*l(3)*(pnu2.x() + ptau2.x()) + pow(Etau2,3)*l(3)*pnu2.x()*pnu2.y()*ptau2.y() - pow(Enu2,3)*l(3)*ptau2.x()*pow(ptau2.y(),2) - pow(Enu2,3)*pow(Etau2,3)*l(0)*PV.z())*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z());
    A(12,20) = 2*W(12,20) + (ptau2.z()*(BV.z()*pow(Enu2,3)*pow(Etau2,3)*l(1) + pow(Etau2,3)*l(3)*pnu2.x()*pnu2.y()*ptau2.x() + 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*ptau2.y() + pow(Etau2,3)*l(3)*pow(pnu2.y(),2)*ptau2.y() - pow(Enu2,3)*l(3)*pow(ptau2.(),2)*ptau2.y() - pow(Enu2,3)*l(3)*pow(ptau2.y(),3) - pow(Enu2,2)*pow(Etau2,3)*l(3)*(pnu2.y() + ptau2.y()) - pow(Enu2,3)*pow(Etau2,3)*l(1)*PV.z())*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z());
    A(12,21) = 2*W(12,21) - ((ptau2.z()*(BV.z()*pow(Enu2,3)*pow(Etau2,3)*(l(0)*pK.x() + l(1)*pK.y()) + 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*pK.x()*ptau2.x() + pow(Etau2,3)*l(3)*pK.x()*pow(pnu2.x(),2)*ptau2.x() + pow(Etau2,3)*l(3)*pK.y()*pnu2.x()*pnu2.y()*ptau2.x() - pow(Enu2,3)*l(3)*pK.x()*pow(ptau2.x(),3) + 2*pow(Enu2,3)*pow(Etau2,2)*l(3)*pK.y()*ptau2.y() + pow(Etau2,3)*l(3)*pK.x()*pnu2.x()*pnu2.y()*ptau2.y() + pow(Etau2,3)*l(3)*pK.y()*pow(pnu2.y(),2)*ptau2.y() - pow(Enu2,3)*l(3)*pK.y()*pow(ptau2.x(),2)*ptau2.y() - pow(Enu2,3)*l(3)*pK.x()*ptau2.x()*pow(ptau2.y(),2) - pow(Enu2,3)*l(3)*pK.y()*pow(ptau2.y(),3) - pow(Enu2,2)*pow(Etau2,3)*l(3)*(pK.x()*(pnu2.x() + ptau2.x()) + pK.y()*(pnu2.y() + ptau2.y())) - pow(Enu2,3)*pow(Etau2,3)*l(0)*pK.x()*PV.z() - pow(Enu2,3)*pow(Etau2,3)*l(1)*pK.y()*PV.z())*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pow(pK.z(),2)));
    A(12,22) = ((BV.z() - DV2.z())*l(0)*ptau2.x() + (BV.z() - DV2.z())*l(1)*ptau2.y() - (l(3)*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y())*(pK.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()) - (pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()))/(pow(Enu2,3)*pK.z()) + (l(3)*(pow(ptau2.x(),2) + pow(ptau2.y(),2))*(pK.z()*(pow(ptau2.x(),2) + pow(ptau2.y(),2)) - (pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z()))/(pow(Etau2,3)*pK.z()) + (l(3)*(-3*pK.z()*(pow(ptau2.x(),2) + pow(ptau2.y(),2)) + 2*(pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z()))/(Etau2*pK.z()) + (l(3)*(pK.z()*(2*pnu2.x()*ptau2.x() + pow(ptau2.x(),2) + 2*pnu2.y()*ptau2.y() + pow(ptau2.y(),2)) - (pK.x()*(pnu2.x() + ptau2.x()) + pK.y()*(pnu2.y() + ptau2.y()))*ptau2.z()))/(Enu2*pK.z()) - 2*l(0)*ptau2.x()*(BV.z() - PV.z()) - 2*l(1)*ptau2.y()*(BV.z() - PV.z()) + (l(0)*pK.x()*ptau2.z()*(BV.z() - PV.z()))/pK.z() + (l(1)*pK.y()*ptau2.z()*(BV.z() - PV.z()))/pK.z())/pow(BV.z() - DV2.z(),2);
    A(12,24) = ((2*l(3)*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + (BV.y() - DV2.y())*ptau2.y()))/Etau2 + (l(3)*(BV.x()*pnu2.x() - DV2.x()*pnu2.x() + BV.y()*pnu2.y() - DV2.y()*pnu2.y() + BV.z()*pnu2.z() - DV2.z()*pnu2.z())*(pnu2x*ptau2x + pnu2y*ptau2y))/pow(Enu2,3) - (l(3)*(BV.x()*(pnu2.x() + ptau2.x()) - DV2.x()*(pnu2.x() + ptau2.x()) + (BV.y() - DV2.y())*(pnu2.y() + ptau2.y())))/Enu2 - (l(3)*(pow(ptau2.x(),2) + pow(ptau2.y(),2))*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + BV.y()*ptau2.y() - DV2.y()*ptau2.y() + BV.z()*ptau2.z() - DV2.z()*ptau2.z()))/pow(Etau2,3) + (BV.x() - DV2.x())*l(0)*(BV.z() - PV.z()) + (BV.y() - DV2.y())*l(1)*(BV.z() - PV.z()))/pow(BV.z() - DV2.z(),2);
    A(12,25) = (ptau2.x()*(BV.z() - PV.z()))/(BV.z() - DV2.z());
    A(12,26) = (ptau2.y()*(BV.z() - PV.z()))/(BV.z() - DV2.z());
    A(12,28) = -(((pnu2.x()*ptau2.x())/Enu2 - pow(ptau2.x(),2)/Etau2 + (pnu2.y()*ptau2.y())/Enu2 - pow(ptau2.y(),2)/Etau2)/(BV.z() - DV2.z()));

    // p3pi2x, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(13,i) = 2*W(13,i);
        }
    }
    A(13,13) = 2*W(13,13) + (l(3)*(-pow(Enu2,2) + pow(pnu2.x(),2)))/pow(Enu2,3);
    A(13,14) = 2*W(13,14) + (l(3)*pnu2.x()*pnu2.y())/pow(Enu2,3);
    A(13,15) = 2*W(13,15) + (l(3)*pnu2.x()*pnu2.z())/pow(Enu2,3);
    A(13,17) = 2*W(13,17) + (l(3)*(pow(Enu2,2) - pow(pnu2.x(),2))*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3));
    A(13,18) = 2*W(13,18) + (l(3)*pnu2.x()*pnu2.y()*ptau2.z())/((-BV.z() + DV2.z())*pow(Enu2,3));
    A(13,19) = 2*W(13,19) + (l(3)*(pow(Enu2,2) - pow(pnu2.x(),2))*ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pK.z());
    A(13,20) = 2*W(13,20) + (l(3)*pnu2.x()*pnu2.y()*ptau2.z()*(BV.z() - RP.z()))/((-BV.z() + DV2.z())*pow(Enu2,3)*pK.z());
    A(13,21) = 2*W(13,21) + -((l(3)*(pow(Enu2,2)*pK.x() - pnu2.x()*(pK.x()*pnu2.x() + pK.y()*pnu2.y()))*ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pow(pK.z(),2)));
    A(13,22) = -((l(3)*(-(pK.z()*pnu2.x()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y())) + pnu2.x()*(pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z() + pow(Enu2,2)*(pK.z()*ptau2.x() - pK.x()*ptau2.z())))/((BV.z() - DV2.z())*pow(Enu2,3)*pK.z()));
    A(13,24) = (l(3)*(BV.x()*(pow(Enu2,2) - pow(pnu2.x(),2)) + DV2.x()*(-pow(Enu2,2) + pow(pnu2.x(),2)) + pnu2.x()*(-(BV.y()*pnu2.y()) + DV2.y()*pnu2.y() - BV.z()*pnu2.z() + DV2.z()*pnu2.z())))/((BV.z() - DV2.z())*pow(Enu2,3));
    A(13,28) = pnu2.x()/Enu2;

    // p3pi2y, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(14,i) = 2*W(14,i);
        }
    }
    A(14,14) = 2*W(14,14) + (l(3)*(-pow(Enu2,2) + pow(pnu2.y(),2)))/pow(Enu2,3);
    A(14,15) = 2*W(14,15) + (l(3)*pnu2.y()*pnu2.z())/pow(Enu2,3);
    A(14,17) = 2*W(14,17) + (l(3)*pnu2.x()*pnu2.y()*ptau2.z())/((-BV.z() + DV2.z())*pow(Enu2,3));
    A(14,18) = 2*W(14,18) + (l(3)*(pow(Enu2,2) - pow(pnu2.y(),2))*ptau2.z())/((BV.z() - DV2.z())*pow(Enu2,3));
    A(14,19) = 2*W(14,19) + (l(3)*pnu2.x()*pnu2.y()*ptau2.z()*(BV.z() - RP.z()))/((-BV.z() + DV2.z())*pow(Enu2,3)*pK.z());
    A(14,20) = 2*W(14,20) + (l(3)*(pow(Enu2,2) - pow(pnu2.y(),2))*ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pK.z());
    A(14,21) = 2*W(14,21) + -((l(3)*(pow(Enu2,2)*pK.y() - pnu2.y()*(pK.x()*pnu2.x() + pK.y()*pnu2.y()))*ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pow(pK.z(),2)));
    A(14,22) = (l(3)*(pK.z()*pnu2.x()*pnu2.y()*ptau2.x() + pK.z()*(-pow(Enu2,2) + pow(pnu2.y(),2))*ptau2.y() + pow(Enu2,2)*pK.y()*ptau2.z() - pnu2.y()*(pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pK.z());
    A(14,24) = (l(3)*(BV.y()*(pow(Enu2,2) - pow(pnu2.y(),2)) + DV2.y()*(-pow(Enu2,2) + pow(pnu2.y(),2)) + pnu2.y()*(-(BV.x()*pnu2.x()) + DV2.x()*pnu2.x() - BV.z()*pnu2.z() + DV2.z()*pnu2.z())))/((BV.z() - DV2.z())*pow(Enu2,3));
    A(14,28) = pnu2.y()/Enu2;

    // p3pi2z, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(15,i) = 2*W(15,i);
        }
    }
    A(15,15) = 2*W(15,15) + (l(3)*(-pow(Enu2,2) + pow(pnu2.z(),2)))/pow(Enu2,3);
    A(15,17) = 2*W(15,17) + (l(3)*pnu2.x()*pnu2.z()*ptau2.z())/((-BV.z() + DV2.z())*pow(Enu2,3));
    A(15,18) = 2*W(15,18) + (l(3)*pnu2.y()*pnu2.z()*ptau2.z())/((-BV.z() + DV2.z())*pow(Enu2,3));
    A(15,19) = 2*W(15,19) + (l(3)*pnu2.x()*pnu2.z()*ptau2.z()*(BV.z() - RP.z()))/((-BV.z() + DV2.z())*pow(Enu2,3)*pK.z()); 
    A(15,20) = 2*W(15,20) + (l(3)*pnu2.y()*pnu2.z()*ptau2.()z*(BV.z() - RP.z()))/((-BV.z() + DV2.z())*pow(Enu2,3)*pK.z());
    A(15,21) = 2*W(15,21) + (l(3)*(pK.x()*pnu2.x() + pK.y()*pnu2.y())*pnu2.z()*ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pow(pK.z(),2));
    A(15,22) = (l(3)*pnu2.z()*(pK.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()) - (pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()))/((BV.z() - DV2.z())*pow(Enu2,3)*pK.z());
    A(15,24) = (l(3)*((-(BV.x()*pnu2.x()) + DV2.x()*pnu2.x() - BV.y()*pnu2.y() + DV2.y()*pnu2.y())*pnu2.z() + BV.z()*(pow(Enu2,2) - pow(pnu2.z(),2)) + DV2.z()*(-pow(Enu2,2) + pow(pnu2.z(),2))))/((BV.z() - DV2.z())*pow(Enu2,3));
    A(15,28) = pnu2.z()/Enu2;

    // E3pi2, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        A(16,i) = 2*W(16,i);
    }
    A(16,28) = -1;


    // RPx, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(17,i) = 2*W(17,i);
        }
    }
    A(17,17) = 2*W(17,17) + (l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.x(),2) - pow(Enu1,3)*pow(ptau1.x(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.x(),2) - pow(Enu2,3)*pow(ptau2.x(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(17,18) = 2*W(17,18) + (l(2)*((pnu1.x()*pnu1.y())/pow(Enu1,3) - (ptau1.x()*ptau1.y())/pow(Etau1,3))*pow(ptau1.z(),2))/pow(BV.z() - DV1.z(),2) + (l(3)*((pnu2.x()*pnu2.y())/pow(Enu2,3) - (ptau2.x()*ptau2.y())/pow(Etau2,3))*pow(ptau2.z(),2))/pow(BV.z() - DV2.z(),2);
    A(17,19) = 2*W(17,19) + ((pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.x(),2) - pow(Enu1,3)*pow(ptau1.x(),2))*pow(ptau1.z(),2) + pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.x(),2) - pow(Enu2,3)*pow(ptau2.x(),2))*pow(ptau2.z(),2))*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(BV.z() - DV2.z(),2)*pow(Enu1,3)*pow(Enu2,3)*pow(Etau1,3)*pow(Etau2,3)*pK.z());
    A(17,20) = 2*W(17,20) + (((l(3)*pnu1.x()*pnu1.y()*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)) - (l(2)*ptau1.x()*ptau1.y()*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Etau1,3)) + (l(3)*pnu2.x()*pnu2.y()*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)) - (l(3)*ptau2.x()*ptau2.y()*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Etau2,3)))*(BV.z() - RP.z()))/pK.z();
    A(17,21) = 2*W(17,21) + -l(0) + (l(2)*(-(pow(Enu1,3)*pow(Etau1,2)*pK.x()) + pow(Enu1,2)*pow(Etau1,3)*pK.x() - pow(Etau1,3)*pK.x()*pow(pnu1.x(),2) - pow(Etau1,3)*pK.y()*pnu1.x()*pnu1.y() + pow(Enu1,3)*pK.x()*pow(ptau1.x(),2) + pow(Enu1,3)*pK.y()*ptau1.x()*ptau1.y())*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pow(pK.z(),2)) + (l(3)*(-(pow(Enu2,3)*pow(Etau2,2)*pK.x()) + pow(Enu2,2)*pow(Etau2,3)*pK.x() - pow(Etau2,3)*pK.x()*pow(pnu2.x(),2) - pow(Etau2,3)*pK.y()*pnu2.x()*pnu2.y() + pow(Enu2,3)*pK.x()*pow(ptau2.x(),2) + pow(Enu2,3)*pK.y()*ptau2.x()*ptau2.y())*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pow(pK.z(),2));
    A(17,22) = (l(2)*ptau1.z()*(pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.x() - 2*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.x() + pow(Enu1,3)*pow(Etau1,2)*pK.x()*ptau1.z() + pow(Enu1,2)*pow(Etau1,3)*(pK.z()*ptau1.x() - pK.x()*ptau1.z()) - pow(Etau1,3)*pnu1.x()*(pK.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()) - (pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()) + pow(Enu1,3)*ptau1.x()*(pK.z()*(pow(ptau1.x(),2) + pow(ptau1.y(),2)) - (pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z())))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z()) + (l(3)*ptau2.z()*(pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.x() - 2*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.x() + pow(Enu2,3)*pow(Etau2,2)*pK.x()*ptau2.z() + pow(Enu2,2)*pow(Etau2,3)*(pK.z()*ptau2.x() - pK.x()*ptau2.z()) - pow(Etau2,3)*pnu2.x()*(pK.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()) - (pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()) + pow(Enu2,3)*ptau2.x()*(pK.z()*(pow(ptau2.x(),2) + pow(ptau2.y(),2)) - (pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z())))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z()) + l(0)*(ptau1.z()/(BV.z() - DV1.z()) + ptau2.z()/(BV.z() - DV2.z()) - ((pow(BV.z() - DV2.z(),2)*ptau1.z() + pow(BV.z() - DV1.z(),2)*ptau2.z())*(BV.z() - PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(BV.z() - DV2.z(),2)));
    A(17,23) = (-(((BV.z() - DV1.z())*l(2)*pnu1.x())/Enu1) - ((BV.z() - DV1.z())*l(2)*ptau1.x())/Enu1 + (2*(BV.z() - DV1.z())*l(2)*ptau1.x())/Etau1 + (l(2)*pnu1.x()*(BV.x()*pnu1.x() - DV1.x()*pnu1.x() + BV.y()*pnu1.y() - DV1.y()*pnu1.y() + BV.z()*pnu1.z() - DV1.z()*pnu1.z())*ptau1.z())/pow(Enu1,3) - (l(2)*ptau1.x()*ptau1.z()*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + BV.y()*ptau1.y() - DV1.y()*ptau1.y() + BV.z()*ptau1.z() - DV1.z()*ptau1.z()))/pow(Etau1,3) + (BV.z() - DV1.z())*l(0)*(DV1.z() - PV.z()))/pow(BV.z() - DV1.z(),2);
    A(17,24) = (-(((BV.z() - DV2.z())*l(3)*pnu2.x())/Enu2) - ((BV.z() - DV2.z())*l(3)*ptau2.x())/Enu2 + (2*(BV.z() - DV2.z())*l(3)*ptau2.x())/Etau2 + (l(3)*pnu2.x()*(BV.x()*pnu2.x() - DV2.x()*pnu2.x() + BV.y()*pnu2.y() - DV2.y()*pnu2.y() + BV.z()*pnu2.z() - DV2.z()*pnu2.z())*ptau2.z())/pow(Enu2,3) - (l(3)*ptau2.x()*ptau2.z()*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + BV.y()*ptau2.y() - DV2.y()*ptau2.y() + BV.z()*ptau2.z() - DV2.z()*ptau2.z()))/pow(Etau2,3) + (BV.z() - DV2.z())*l(0)*(DV2.z() - PV.z()))/pow(BV.z() - DV2.z(),2);
    A(17,25) = -pB.z() + (pnu2.x()*ptau2.z())/((-BV.z() + DV2.z())*Enu2) + (ptau2.x()*ptau2.z())/((BV.z() - DV2.z())*Etau2) + (ptau1.z()/(BV.z() - DV1.z()) + ptau2.z()/(BV.z() - DV2.z()))*(BV.z() - PV.z())
    A(17,27) = ((pnu1.x()/Enu1 - ptau1.x()/Etau1)*ptau1.z())/(-BV.z() + DV1.z());
    A(17,28) = ((pnu2.x()/Enu2 - ptau2.x()/Etau2)*ptau2.z())/(-BV.z() + DV2.z());

    // RPy, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(18,i) = 2*W(18,i);
        }
    }
    A(18,18) = 2*W(18,18) + (l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.y(),2) - pow(Enu1,3)*pow(ptau1.y(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.y(),2) - pow(Enu2,3)*pow(ptau2.y(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3));
    A(18,19) = 2*W(18,19) + (((l(2)*pnu1.x()*pnu1.y()*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)) - (l(2)*ptau1.x()*ptau1.y()*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Etau1,3)) + (l(3)*pnu2.x()*pnu2.y()*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)) - (l(3)*ptau2.x()*ptau2.y()*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Etau2,3)))*(BV.z() - RP.z()))/pK.z();
    A(18,20) = 2*W(19,20) + ((pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.y(),2) - pow(Enu1,3)*pow(ptau1.y(),2))*pow(ptau1.z(),2) + pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.y(),2) - pow(Enu2,3)*pow(ptau2.y(),2))*pow(ptau2.z(),2))*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(BV.z() - DV2.z(),2)*pow(Enu1,3)*pow(Enu2,3)*pow(Etau1,3)*pow(Etau2,3)*pK.z());
    A(18,21) = 2*W(18,21) + -l(1) + (l(2)*(-(pow(Enu1,3)*pow(Etau1,2)*pK.y()) + pow(Enu1,2)*pow(Etau1,3)*pK.y() - pow(Etau1,3)*pK.x()*pnu1.x()*pnu1.y() - pow(Etau1,3)*pK.y()*pow(pnu1.y(),2) + pow(Enu1,3)*pK.x()*ptau1.x()*ptau1.y() + pow(Enu1,3)*pK.y()*pow(ptau1.y(),2))*pow(ptau1.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pow(pK.z(),2)) + (l(3)*(-(pow(Enu2,3)*pow(Etau2,2)*pK.y()) + pow(Enu2,2)*pow(Etau2,3)*pK.y() - pow(Etau2,3)*pK.x()*pnu2.x()*pnu2.y() - pow(Etau2,3)*pK.y()*pow(pnu2.y(),2) + pow(Enu2,3)*pK.x()*ptau2.x()*ptau2.y() + pow(Enu2,3)*pK.y()*pow(ptau2.y(),2))*pow(ptau2.z(),2)*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pow(pK.z(),2));
    A(18,22) = (l(2)*ptau1.z()*(pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.y() - 2*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.y() + pow(Enu1,3)*pow(Etau1,2)*pK.y()*ptau1.z() + pow(Enu1,2)*pow(Etau1,3)*(pK.z()*ptau1.y() - pK.y()*ptau1.z()) - pow(Etau1,3)*pnu1.y()*(pK.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()) - (pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()) + pow(Enu1,3)*ptau1.y()*(pK.z()*(pow(ptau1.x(),2) + pow(ptau1.y(),2)) - (pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z())))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)*pK.z()) + (l(3)*ptau2.z()*(pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.y() - 2*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.y() + pow(Enu2,3)*pow(Etau2,2)*pK.y()*ptau2.z() + pow(Enu2,2)*pow(Etau2,3)*(pK.z()*ptau2.y() - pK.y()*ptau2.z()) - pow(Etau2,3)*pnu2.y()*(pK.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()) - (pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()) + pow(Enu2,3)*ptau2.y()*(pK.z()*(pow(ptau2.x(),2) + pow(ptau2.y(),2)) - (pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z())))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)*pK.z()) + l(1)*(ptau1.z()/(BV.z() - DV1.z()) + ptau2.z()/(BV.z() - DV2.z()) - ((pow(BV.z() - DV2.z(),2)*ptau1.z() + pow(BV.z() - DV1.z(),2)*ptau2.z())*(BV.z() - PV.z()))/(pow(BV.z() - DV1.z(),2)*pow(BV.z() - DV2.z(),2)));
    A(18,23) = (-(((BV.z() - DV1.z())*l(2)*pnu1.y())/Enu1) - ((BV.z() - DV1.z())*l(2)*ptau1.y())/Enu1 + (2*(BV.z() - DV1.z())*l(2)*ptau1.y())/Etau1 + (l(2)*pnu1.y()*(BV.x()*pnu1.x() - DV1.x()*pnu1.x() + BV.y()*pnu1.y() - DV1.y()*pnu1.y() + BV.z()*pnu1.z() - DV1.z()*pnu1.z())*ptau1.z())/pow(Enu1,3) - (l(2)*ptau1.y()*ptau1.z()*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + BV.y()*ptau1.y() - DV1.y()*ptau1.y() + BV.z()*ptau1.z() - DV1.z()*ptau1.z()))/pow(Etau1,3) + (BV.z() - DV1.z())*l(1)*(DV1.z() - PV.z()))/pow(BV.z() - DV1.z(),2);
    A(18,24) = (-(((BV.z() - DV2.z())*l(3)*pnu2.y())/Enu2) - ((BV.z() - DV2.z())*l(3)*ptau2.y())/Enu2 + (2*(BV.z() - DV2.z())*l(3)*ptau2.y())/Etau2 + (l(3)*pnu2.y()*(BV.x()*pnu2.x() - DV2.x()*pnu2.x() + BV.y()*pnu2.y() - DV2.y()*pnu2.y() + BV.z()*pnu2.z() - DV2.z()*pnu2.z())*ptau2.z())/pow(Enu2,3) - (l(3)*ptau2.y()*ptau2.z()*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + BV.y()*ptau2.y() - DV2.y()*ptau2.y() + BV.z()*ptau2.z() - DV2.z()*ptau2.z()))/pow(Etau2,3) + (BV.z() - DV2.z())*l(1)*(DV2.z() - PV.z()))/pow(BV.z() - DV2.z(),2);
    A(18,26) = -pB.z() + (ptau1.z()/(BV.z() - DV1.z()) + ptau2.z()/(BV.z() - DV2.z()))*(BV.z() - PV.z());
    A(18,27) = ((pnu1.y()/Enu1 - ptau1.y()/Etau1)*ptau1.z())/(-BV.z() + DV1.z());
    A(18,28) = ((pnu2.y()/Enu2 - ptau2.y()/Etau2)*ptau2.z())/(-BV.z() + DV2.z());

    // pKx,gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(19,i) = 2*W(19,i);
        }
    }
    A(19,19) = 2*W(19,19) + (((l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.x(),2) - pow(Enu1,3)*pow(ptau1.x(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.x(),2) - pow(Enu2,3)*pow(ptau2.x(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)))*pow(BV.z() - RP.z(),2))/pow(pK.z(),2);
    A(19,20) = 2*W(19,20) + (((l(2)*(pow(Etau1,3)*pnu1.x()*pnu1.y() - pow(Enu1,3)*ptau1.x()*ptau1.y())*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*(pow(Etau2,3)*pnu2.x()*pnu2.y() - pow(Enu2,3)*ptau2.x()*ptau2.y())*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)))*pow(BV.z() - RP.z(),2))/pow(pK.z(),2);
    A(19,21) = 2*W(19,21) + ((l(0)*pK.z()*(pB.z() - pK.z() + ((DV2.z()*ptau1.z() + DV1.z()*ptau2.z() - BV.z()*(ptau1.z() + ptau2.z()))*(BV.z() - PV.z()))/((BV.z() - DV1.z())*(BV.z() - DV2.z()))) + (l(2)*ptau1.z()*((BV.z() - DV1.z())*pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.x() - (BV.z() - DV1.z())*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.x() - pow(Enu1,3)*pow(Etau1,2)*pK.x()*ptau1.z()*(BV.z() - RP.z()) + pow(Enu1,2)*pow(Etau1,3)*pK.x()*ptau1.z()*(BV.z() - RP.z()) - pow(Etau1,3)*pnu1.x()*(pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()*(BV.z() - RP.z()) + pow(Enu1,3)*ptau1.x()*(pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z()*(BV.z() - RP.z())))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*ptau2.z()*((BV.z() - DV2.z())*pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.x() - (BV.z() - DV2.z())*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.x() - pow(Enu2,3)*pow(Etau2,2)*pK.x()*ptau2.z()*(BV.z() - RP.z()) + pow(Enu2,2)*pow(Etau2,3)*pK.x()*ptau2.z()*(BV.z() - RP.z()) - pow(Etau2,3)*pnu2.x()*(pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()*(BV.z() - RP.z()) + pow(Enu2,3)*ptau2.x()*(pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z()*(BV.z() - RP.z())))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)))*(BV.z() - RP.z()))/pow(pK.z(),3);
    A(19,22) = (-((l(2)*ptau1.z()*((BV.z() - DV1.z())*pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.x() - (BV.z() - DV1.z())*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.x() - pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.x()*(BV.z() - RP.z()) + 2*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.x()*(BV.z() - RP.z()) - pow(Enu1,3)*pow(Etau1,2)*pK.x()*ptau1.z()*(BV.z() - RP.z()) - pow(Enu1,2)*pow(Etau1,3)*(pK.z()*ptau1.x() - pK.x()*ptau1.z())*(BV.z() - RP.z()) + pow(Etau1,3)*pnu1.x()*(pK.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()) - (pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z())*(BV.z() - RP.z()) - pow(Enu1,3)*ptau1.x()*(pK.z()*(pow(ptau1.x(),2) + pow(ptau1.y(),2)) - (pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z())*(BV.z() - RP.z())))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3))) - (l(3)*ptau2.z()*((BV.z() - DV2.z())*pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.x() - (BV.z() - DV2.z())*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.x() - pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.x()*(BV.z() - RP.z()) + 2*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.x()*(BV.z() - RP.z()) - pow(Enu2,3)*pow(Etau2,2)*pK.x()*ptau2.z()*(BV.z() - RP.z()) - pow(Enu2,2)*pow(Etau2,3)*(pK.z()*ptau2.x() - pK.x()*ptau2.z())*(BV.z() - RP.z()) + pow(Etau2,3)*pnu2.x()*(pK.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()) - (pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z())*(BV.z() - RP.z()) - pow(Enu2,3)*ptau2.x()*(pK.z()*(pow(ptau2.x(),2) + pow(ptau2.y(),2)) - (pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z())*(BV.z() - RP.z())))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)) + l(0)*pK.z()*(-pB.z() + pK.z() + (ptau1.z()*(BV.z() - RP.z()))/(BV.z() - DV1.z()) + (ptau2.z()*(BV.z() - RP.z()))/(BV.z() - DV2.z()) - ((BV.z() - PV.z())*(DV1.z()*pow(DV2.z(),2)*ptau1.z() + pow(DV1.z(),2)*ptau2.z()*(DV2.z() - RP.z()) + 2*BV.z()*DV2.z()*ptau1.z()*RP.z() - pow(DV2.z(),2)*ptau1.z()*RP.z() - 2*BV.z()*DV1.z()*(DV2.z()*(ptau1.z() + ptau2.z()) - ptau2.z()*RP.z()) + pow(BV.z(),2)*(DV1.z()*ptau1.z() + DV2.z()*ptau2.z() - (ptau1.z() + ptau2.z())*RP.z())))/(pow(BV.z() - DV1.z(),2)*pow(BV.z() - DV2.z(),2))))/pow(pK.z(),2);
    A(19,23) = ((-(((BV.z() - DV1.z())*l(2)*pnu1.x())/Enu1) - ((BV.z() - DV1.z())*l(2)*ptau1.x())/Enu1 + (2*(BV.z() - DV1.z())*l(2)*ptau1.x())/Etau1 + (l(2)*pnu1.x()*(BV.x()*pnu1.x() - DV1.x()*pnu1.x() + BV.y()*pnu1.y() - DV1.y()*pnu1.y() + BV.z()*pnu1.z() - DV1.z()*pnu1.z())*ptau1.z())/pow(Enu1,3) - (l(2)*ptau1.x()*ptau1.z()*(BV.x()*ptau1.x() - DV1.x()*ptau1.x() + BV.y()*ptau1.y() - DV1.y()*ptau1.y() + BV.z()*ptau1.z() - DV1.z()*ptau1.z()))/pow(Etau1,3) + (BV.z() - DV1.z())*l(0)*(DV1.z() - PV.z()))*(BV.z() - RP.z()))/(pow(BV.z() - DV1.z(),2)*pK.z());
    A(19,24) = ((-(((BV.z() - DV2.z())*l(3)*pnu2.x())/Enu2) - ((BV.z() - DV2.z())*l(3)*ptau2.x())/Enu2 + (2*(BV.z() - DV2.z())*l(3)*ptau2.x())/Etau2 + (l(3)*pnu2.x()*(BV.x()*pnu2.x() - DV2.x()*pnu2.x() + BV.y()*pnu2.y() - DV2.y()*pnu2.y() + BV.z()*pnu2.z() - DV2.z()*pnu2.z())*ptau2.z())/pow(Enu2,3) - (l(3)*ptau2.x()*ptau2.z()*(BV.x()*ptau2.x() - DV2.x()*ptau2.x() + BV.y()*ptau2.y() - DV2.y()*ptau2.y() + BV.z()*ptau2.z() - DV2.z()*ptau2.z()))/pow(Etau2,3) + (BV.z() - DV2.z())*l(0)*(DV2.z() - PV.z()))*(BV.z() - RP.z()))/(pow(BV.z() - DV2.z(),2)*pK.z());
    A(19,25) = (BV.z() - PV.z())*(1 + (ptau1.z()*(BV.z() - RP.z()))/((BV.z() - DV1.z())*pK.z()) + (ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*pK.z())) + (pB.z()*(-BV.z() + RP.z()))/pK.z();
    A(19,27) = ((-(Etau1*pnu1.x()) + Enu1*ptau1.x())*ptau1.z()*(BV.z() - RP.z()))/((BV.z() - DV1.z())*Enu1*Etau1*pK.z());
    A(19,28) = ((-(Etau2*pnu2.x()) + Enu2*ptau2.x())*ptau2.z()*(BV.z() - RP.z()))/((BV.z() - DV2.z())*Enu2*Etau2*pK.z());

    // pKy, gamma_j
    for(int i = 0; i < dimM; i++)
    {
        if((i!=2) && (i!=3) && (i!=4) && (i!=5) && (i!=6) && (i!=7) && (i!=8) && (i!=10) && (i!=11) && (i!=12) && (i!=13) && (i!=14) && (i!=15) && (i!=17) && (i!=18) && (i!=19) && (i!=20) && (i!=21) )
        {
            A(20,i) = 2*W(20,i);
        }
    }
    A(20,20) = 2*W(20,20) + (((l(2)*(pow(Enu1,3)*pow(Etau1,2) - pow(Enu1,2)*pow(Etau1,3) + pow(Etau1,3)*pow(pnu1.y(),2) - pow(Enu1,3)*pow(ptau1.y(),2))*pow(ptau1.z(),2))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*(pow(Enu2,3)*pow(Etau2,2) - pow(Enu2,2)*pow(Etau2,3) + pow(Etau2,3)*pow(pnu2.y(),2) - pow(Enu2,3)*pow(ptau2.y(),2))*pow(ptau2.z(),2))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)))*pow(BV.z() - RP.z(),2))/pow(pK.z(),2);
    A(20,21) = 2*W(20,21) + ((l(1)*pK.z()*(pB.z() - pK.z() + ((DV2.z()*ptau1.z() + DV1.z()*ptau2.z() - BV.z()*(ptau1.z() + ptau2.z()))*(BV.z() - PV.z()))/((BV.z() - DV1.z())*(BV.z() - DV2.z()))) + (l(2)*ptau1.z()*((BV.z() - DV1.z())*pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.y() - (BV.z() - DV1.z())*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.y() - pow(Enu1,3)*pow(Etau1,2)*pK.y()*ptau1.z()*(BV.z() - RP.z()) + pow(Enu1,2)*pow(Etau1,3)*pK.y()*ptau1.z()*(BV.z() - RP.z()) - pow(Etau1,3)*pnu1.y()*(pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z()*(BV.z() - RP.z()) + pow(Enu1,3)*ptau1.y()*(pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z()*(BV.z() - RP.z())))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3)) + (l(3)*ptau2.z()*((BV.z() - DV2.z())*pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.y() - (BV.z() - DV2.z())*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.y() - pow(Enu2,3)*pow(Etau2,2)*pK.y()*ptau2.z()*(BV.z() - RP.z()) + pow(Enu2,2)*pow(Etau2,3)*pK.y()*ptau2.z()*(BV.z() - RP.z()) - pow(Etau2,3)*pnu2.y()*(pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z()*(BV.z() - RP.z()) + pow(Enu2,3)*ptau2.y()*(pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z()*(BV.z() - RP.z())))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)))*(BV.z() - RP.z()))/pow(pK.z(),3);
    A(20,22) = (-((l(2)*ptau1.z()*((BV.z() - DV1.z())*pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.y() - (BV.z() - DV1.z())*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.y() - pow(Enu1,2)*pow(Etau1,3)*pK.z()*pnu1.y()*(BV.z() - RP.z()) + 2*pow(Enu1,3)*pow(Etau1,2)*pK.z()*ptau1.y()*(BV.z() - RP.z()) - pow(Enu1,3)*pow(Etau1,2)*pK.y()*ptau1.z()*(BV.z() - RP.z()) - pow(Enu1,2)*pow(Etau1,3)*(pK.z()*ptau1.y() - pK.y()*ptau1.z())*(BV.z() - RP.z()) + pow(Etau1,3)*pnu1.y()*(pK.z()*(pnu1.x()*ptau1.x() + pnu1.y()*ptau1.y()) - (pK.x()*pnu1.x() + pK.y()*pnu1.y())*ptau1.z())*(BV.z() - RP.z()) - pow(Enu1,3)*ptau1.y()*(pK.z()*(pow(ptau1.x(),2) + pow(ptau1.y(),2)) - (pK.x()*ptau1.x() + pK.y()*ptau1.y())*ptau1.z())*(BV.z() - RP.z())))/(pow(BV.z() - DV1.z(),2)*pow(Enu1,3)*pow(Etau1,3))) - (l(3)*ptau2.z()*((BV.z() - DV2.z())*pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.y() - (BV.z() - DV2.z())*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.y() - pow(Enu2,2)*pow(Etau2,3)*pK.z()*pnu2.y()*(BV.z() - RP.z()) + 2*pow(Enu2,3)*pow(Etau2,2)*pK.z()*ptau2.y()*(BV.z() - RP.z()) - pow(Enu2,3)*pow(Etau2,2)*pK.y()*ptau2.z()*(BV.z() - RP.z()) - pow(Enu2,2)*pow(Etau2,3)*(pK.z()*ptau2.y() - pK.y()*ptau2.z())*(BV.z() - RP.z()) + pow(Etau2,3)*pnu2.y()*(pK.z()*(pnu2.x()*ptau2.x() + pnu2.y()*ptau2.y()) - (pK.x()*pnu2.x() + pK.y()*pnu2.y())*ptau2.z())*(BV.z() - RP.z()) - pow(Enu2,3)*ptau2.y()*(pK.z()*(pow(ptau2.x(),2) + pow(ptau2.y(),2)) - (pK.x()*ptau2.x() + pK.y()*ptau2.y())*ptau2.z())*(BV.z() - RP.z())))/(pow(BV.z() - DV2.z(),2)*pow(Enu2,3)*pow(Etau2,3)) + l(1)*pK.z()*(-pB.z() + pK.z() + (ptau1.z()*(BV.z() - RP.z()))/(BV.z() - DV1.z()) + (ptau2.z()*(BV.z() - RP.z()))/(BV.z() - DV2.z()) - ((BV.z() - PV.z())*(DV1.z()*pow(DV2.z(),2)*ptau1.z() + pow(DV1.z(),2)*ptau2.z()*(DV2.z() - RP.z()) + 2*BV.z()*DV2.z()*ptau1.z()*RP.z() - pow(DV2.z(),2)*ptau1.z()*RP.z() - 2*BV.z()*DV1.z()*(DV2.z()*(ptau1.z() + ptau2.z()) - ptau2.z()*RP.z()) + pow(BV.z(),2)*(DV1.z()*ptau1.z() + DV2.z()*ptau2.z() - (ptau1.z() + ptau2.z())*RP.z())))/(pow(BV.z() - DV1.z(),2)*pow(BV.z() - DV2.z(),2))))/pow(pK.z(),2);




    A(21,21) = 2*W(21,21) - ((pow(EK,2) - pow(pK.z(),2))/pow(EK,3))*l(19);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 21)
        {
            A(21,i) = 2*W(21,i);
        }
    }
    A(21,24) = l(15)/pow(BV.z()-RP.z(),2);
    A(21,L+15) = -1./(BV.z()-RP.z());
    A(21,L+18) = -1;
    A(21,L+19) = -pK.z()/EK;

    //////////////////////////////// xu /////////////////////////
    A(22,0) = -(2*pB.x()/pow(BV.x()-PV.x(),3))*(l(0)+l(1));
    A(22,3) = -(2*ptau1.x()/pow(DV1.x()-BV.x(),3))*(l(2)+l(3));
    A(22,10) = -(2*ptau2.x()/pow(DV2.x()-BV.x(),3))*(l(8)+l(9));
    A(22,17) = -(2*pK.x()/pow(BV.x()-RP.x(),3))*(l(14)+l(15));
    A(22,19) = -(l(14)+l(15))/pow(BV.x()-RP.x(),2);
    A(22,22) = (2*pB.x()/pow(BV.x()-PV.x(),3))*(l(0)+l(1)) + (2*ptau1.x()/pow(DV1.x()-BV.x(),3))*(l(2)+l(3)) + (2*ptau2.x()/pow(DV2.x()-BV.x(),3))*(l(8)+l(9)) + (2*pK.x()/pow(BV.x()-RP.x(),3))*(l(14)+l(15));
    A(22,25) = -(l(0)+l(1))/pow(BV.x()-PV.x(),2);
    A(22,29) = (l(2)+l(3))/pow(DV1.x()-BV.x(),2);
    A(22,35) = (l(8)+l(9))/pow(DV2.x()-BV.x(),2);
    A(22,L) = -pB.x()/pow(BV.x()-PV.x(),2);
    A(22,L+1) = -pB.x()/pow(BV.x()-PV.x(),2);
    A(22,L+2) = ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(22,L+3) = ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(22,L+8) = ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(22,L+9) = ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(22,L+14) = -pK.x()/pow(BV.x()-RP.x(),2);
    A(22,L+15) = -pK.x()/pow(BV.x()-RP.x(),2);

    A(23,1) = (2*pB.y()/pow(BV.y()-PV.y(),3))*l(0);
    A(23,4) = (2*ptau1.y()/pow(DV1.y()-BV.y(),3))*l(2);
    A(23,11) = (2*ptau2.y()/pow(DV2.y()-BV.y(),3))*l(8);
    A(23,18) = (2*pK.y()/pow(BV.y()-RP.y(),3))*l(14);
    A(23,20) = l(14)/pow(BV.y()-RP.y(),2);
    A(23,23) = -(2*pB.y()/pow(BV.y()-PV.y(),3))*l(0) - (2*ptau1.y()/pow(DV1.y()-BV.y(),3))*l(2) - (2*ptau2.y()/pow(DV2.y()-BV.y(),3))*l(8) - (2*pK.y()/pow(BV.y()-RP.y(),3))*l(14);
    A(23,26) = l(0)/pow(BV.y()-PV.y(),2);
    A(23,30) = -l(2)/pow(DV1.y()-BV.y(),2);
    A(23,36) = -l(8)/pow(DV2.y()-BV.y(),2);
    A(23,L) = pB.y()/pow(BV.y()-PV.y(),2);
    A(23,L+2) = -ptau1.y()/pow(DV1.y()-BV.y(),2);
    A(23,L+8) = -ptau2.y()/pow(DV2.y()-BV.y(),2);
    A(23,L+14) = pK.y()/pow(BV.y()-RP.y(),2);

    A(24,2) = (2*pB.z()/pow(BV.z()-PV.z(),3))*l(1);
    A(24,5) = (2*ptau1.z()/pow(DV1.z()-BV.z(),3))*l(3);
    A(24,12) = (2*ptau2.z()/pow(DV2.z()-BV.z(),3))*l(9);
    A(24,21) = l(15)/pow(BV.z()-RP.z(),2);
    A(24,24) = -(2*pB.z()/pow(BV.z()-PV.z(),3))*l(1) - (2*ptau1.z()/pow(DV1.z()-BV.z(),3))*l(3) - (2*ptau2.z()/pow(DV2.z()-BV.z(),3))*l(9) - (2*pK.z()/pow(BV.z()-RP.z(),3))*l(15);
    A(24,27) = l(1)/pow(BV.z()-PV.z(),2);
    A(24,31) = -l(3)/pow(DV1.z()-BV.z(),2);
    A(24,37) = -l(9)/pow(DV2.z()-BV.z(),2);
    A(24,L+1) = pB.z()/pow(BV.z()-PV.z(),2);
    A(24,L+3) = -ptau1.z()/pow(DV1.z()-BV.z(),2);
    A(24,L+9) = -ptau2.z()/pow(DV2.z()-BV.z(),2);
    A(24,L+15) = pK.z()/pow(BV.z()-RP.z(),2);

    A(25,0) = (l(0)+l(1))/pow(BV.x()-PV.x(),2);
    A(25,22) = -(l(0)+l(1))/pow(BV.x()-PV.x(),2);
    A(25,L) = 1./(BV.x()-PV.x());
    A(25,L+1) = 1./(BV.x()-PV.x());
    A(25,L+16) = 1;

    A(26,1) = -l(0)/pow(BV.y()-PV.y(),2);
    A(26,23) = l(0)/pow(BV.y()-PV.y(),2);
    A(26,L) = -1./(BV.y()-PV.y());
    A(26,L+17) = 1;

    A(27,2) = -l(1)/pow(BV.z()-PV.z(),2);
    A(27,24) = l(1)/pow(BV.z()-PV.z(),2);
    A(27,L+1) = -1./(BV.z()-PV.z());
    A(27,L+18) = 1;

    A(28,25) = (pB.x()/EB)*l(19);
    A(28,26) = (pB.y()/EB)*l(19);
    A(28,27) = (pB.z()/EB)*l(19);
    A(28,28) = -l(19)/(4*pow(EB,3));
    A(28,L+19) = 1./(2*EB);

    A(29,3) = -(l(2)+l(3))/pow(DV1.x()-BV.x(),2);
    A(29,22) = (l(2)+l(3))/pow(DV1.x()-BV.x(),2);
    A(29,29) = ((pow(Etau1,2) - pow(ptau1.x(),2))/pow(Etau1,3))*(l(7)-l(19));
    A(29,L+2) = 1./(DV1.x()-BV.x());
    A(29,L+3) = 1./(DV1.x()-BV.x());
    A(29,L+4) = 1;
    A(29,L+7) = ptau1.x()/Etau1;
    A(29,L+16) = -1;
    A(29,L+19) = -ptau1.x()/Etau1;

    A(30,4) = l(2)/pow(DV1.y()-BV.y(),2);
    A(30,23) = -l(2)/pow(DV1.y()-BV.y(),2);
    A(30,30) = ((pow(Etau1,2) - pow(ptau1.y(),2))/pow(Etau1,3))*(l(7)-l(19));
    A(30,L+2) = -1./(DV1.y()-BV.y());
    A(30,L+5) = 1;
    A(30,L+7) = ptau1.y()/Etau1;
    A(30,L+17) = -1;
    A(30,L+19) = -ptau1.y()/Etau1;

    A(31,5) = l(3)/pow(DV1.z()-BV.z(),2);
    A(31,24) = -l(3)/pow(DV1.z()-BV.z(),2);
    A(31,31) = ((pow(Etau1,2) - pow(ptau1.z(),2))/pow(Etau1,3))*(l(7)-l(19));
    A(31,L+3) = -1./(DV1.z()-BV.z());
    A(31,L+6) = 1;
    A(31,L+7) = ptau1.z()/Etau1;
    A(31,L+18) = -1;
    A(31,L+19) = -ptau1.z()/Etau1;

    A(32,32) = -((pow(Enu1,2) - pow(pnu1.x(),2))/pow(Enu1,3))*l(7);
    A(32,L+4) = -1;
    A(32,L+7) = -pnu1.x()/Enu1;

    A(33,33) = -((pow(Enu1,2) - pow(pnu1.y(),2))/pow(Enu1,3))*l(7);
    A(33,L+5) = -1;
    A(33,L+7) = -pnu1.y()/Enu1;

    A(34,34) = -((pow(Enu1,2) - pow(pnu1.z(),2))/pow(Enu1,3))*l(7);
    A(34,L+6) = -1;
    A(34,L+7) = -pnu1.z()/Enu1;

    A(35,10) = -(l(8)+l(9))/pow(DV2.x()-BV.x(),2);
    A(35,22) = (l(8)+l(9))/pow(DV2.x()-BV.x(),2);
    A(35,35) = ((pow(Etau2,2) - pow(ptau2.x(),2))/pow(Etau2,3))*(l(13)-l(19));
    A(35,L+8) = 1./(DV2.x()-BV.x());
    A(35,L+9) = 1./(DV2.x()-BV.x());
    A(35,L+10) = 1;
    A(35,L+13) = ptau2.x()/Etau2;
    A(35,L+16) = -1;
    A(35,L+19) = -ptau2.x()/Etau2;

    A(36,11) = l(8)/pow(DV2.y()-BV.y(),2);
    A(36,23) = -l(8)/pow(DV2.y()-BV.y(),2);
    A(36,36) = ((pow(Etau2,2) - pow(ptau2.y(),2))/pow(Etau2,3))*(l(13)-l(19));
    A(36,L+8) = -1./(DV2.y()-BV.y());
    A(36,L+11) = 1;
    A(36,L+13) = ptau2.y()/Etau2;
    A(36,L+17) = -1;
    A(36,L+19) = -ptau2.y()/Etau2;

    A(37,12) = l(9)/pow(DV2.z()-BV.z(),2);
    A(37,24) = -l(9)/pow(DV2.z()-BV.z(),2);
    A(37,37) = ((pow(Etau2,2) - pow(ptau2.z(),2))/pow(Etau2,3))*(l(13)-l(19));
    A(37,L+9) = -1./(DV2.z()-BV.z());
    A(37,L+12) = 1;
    A(37,L+13) = ptau2.z()/Etau2;
    A(37,L+18) = -1;
    A(37,L+19) = -ptau2.z()/Etau2;

    A(38,38) = -((pow(Enu2,2) - pow(pnu2.x(),2))/pow(Enu2,3))*l(13);
    A(38,L+10) = -1;
    A(38,L+13) = -pnu2.x()/Enu2;

    A(39,39) = -((pow(Enu2,2) - pow(pnu2.y(),2))/pow(Enu2,3))*l(13);
    A(39,L+11) = -1;
    A(39,L+13) = -pnu2.y()/Enu2;

    A(40,40) = -((pow(Enu2,2) - pow(pnu2.z(),2))/pow(Enu2,3))*l(13);
    A(40,L+12) = -1;
    A(40,L+13) = -pnu2.z()/Enu2;

    ///////////////////////////// lambda /////////////////////////
    A(41,0) = pB.x()/pow(BV.x()-PV.x(),2);
    A(41,1) = -pB.y()/pow(BV.y()-PV.y(),2);
    A(41,22) = -pB.x()/pow(BV.x()-PV.x(),2);
    A(41,23) = pB.y()/pow(BV.y()-PV.y(),2);
    A(41,25) = 1./(BV.x()-PV.x());
    A(41,26) = -1./(BV.y()-PV.y());

    A(42,0) = pB.x()/pow(BV.x()-PV.x(),2);
    A(42,2) = -pB.z()/pow(BV.z()-PV.z(),2);
    A(42,22) = -pB.x()/pow(BV.x()-PV.x(),2);
    A(42,24) = pB.z()/pow(BV.z()-PV.z(),2);
    A(42,25) = 1./(BV.x()-PV.x());
    A(42,27) = -1./(BV.z()-PV.z());

    A(43,3) = -ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(43,4) = ptau1.y()/pow(DV1.y()-BV.y(),2);
    A(43,22) = ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(43,23) = -ptau1.y()/pow(DV1.y()-BV.y(),2);
    A(43,29) = 1./(DV1.x()-BV.x());
    A(43,30) = -1./(DV1.y()-BV.y());

    A(44,3) = -ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(44,5) = ptau1.z()/pow(DV1.z()-BV.z(),2);
    A(44,22) = ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(44,24) = -ptau1.z()/pow(DV1.z()-BV.z(),2);
    A(44,29) = 1./(DV1.x()-BV.x());
    A(44,31) = -1./(DV1.z()-BV.z());

    A(45,29) = 1;
    A(45,6) = -1;
    A(45,32) = -1;

    A(46,30) = 1;
    A(46,7) = -1;
    A(46,33) = -1;

    A(47,31) = 1;
    A(47,8) = -1;
    A(47,34) = -1;

    A(48,29) = ptau1.x()/Etau1;
    A(48,30) = ptau1.y()/Etau1;
    A(48,31) = ptau1.z()/Etau1;
    A(48,9) = -1;
    A(48,32) = -pnu1.x()/Enu1;
    A(48,33) = -pnu1.y()/Enu1;
    A(48,34) = -pnu1.z()/Enu1;

    A(49,10) = -ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(49,11) = ptau2.y()/pow(DV2.y()-BV.y(),2);
    A(49,22) = ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(49,23) = -ptau2.y()/pow(DV2.y()-BV.y(),2);
    A(49,35) = 1./(DV2.x()-BV.x());
    A(49,36) = -1./(DV2.y()-BV.y());

    A(50,10) = -ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(50,12) = ptau2.z()/pow(DV2.z()-BV.z(),2);
    A(50,22) = ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(50,24) = -ptau2.z()/pow(DV2.z()-BV.z(),2);
    A(50,35) = 1./(DV2.x()-BV.x());
    A(50,37) = -1./(DV2.z()-BV.z());

    A(51,35) = 1;
    A(51,13) = -1;
    A(51,38) = -1;

    A(52,36) = 1;
    A(52,14) = -1;
    A(52,39) = -1;

    A(53,37) = 1;
    A(53,15) = -1;
    A(53,40) = -1;

    A(54,35) = ptau2.x()/Etau2;
    A(54,36) = ptau2.y()/Etau2;
    A(54,37) = ptau2.z()/Etau2;
    A(54,16) = -1;
    A(54,38) = -pnu2.x()/Enu2;
    A(54,39) = -pnu2.y()/Enu2;
    A(54,40) = -pnu2.z()/Enu2;

    A(55,17) = pK.x()/pow(BV.x()-RP.x(),2);
    A(55,18) = -pK.y()/pow(BV.y()-RP.y(),2);
    A(55,22) = -pK.x()/pow(BV.x()-RP.x(),2);
    A(55,23) = pK.y()/pow(BV.y()-RP.y(),2);
    A(55,19) = 1./(BV.x()-RP.x());
    A(55,20) = -1./(BV.y()-RP.y());

    A(56,17) = pK.x()/pow(BV.x()-RP.x(),2);
    A(56,22) = -pK.x()/pow(BV.x()-RP.x(),2);
    A(56,24) = pK.z()/pow(BV.z()-RP.z(),2);
    A(56,19) = 1./(BV.x()-RP.x());
    A(56,21) = -1./(BV.z()-RP.z());

    A(57,25) = 1;
    A(57,29) = -1;
    A(57,35) = -1;
    A(57,19) = -1;

    A(58,26) = 1;
    A(58,30) = -1;
    A(58,36) = -1;
    A(58,20) = -1;

    A(59,27) = 1;
    A(59,31) = -1;
    A(59,37) = -1;
    A(59,21) = -1;

    A(60,25) = pB.x()/EB;
    A(60,26) = pB.y()/EB;
    A(60,27) = pB.z()/EB;
    A(60,28) = 1./(2*EB);
    A(60,29) = -ptau1.x()/Etau1;
    A(60,30) = -ptau1.y()/Etau1;
    A(60,31) = -ptau1.z()/Etau1;
    A(60,35) = -ptau2.x()/Etau2;
    A(60,36) = -ptau2.y()/Etau2;
    A(60,37) = -ptau2.z()/Etau2;
    A(60,19) = -pK.x()/EK;
    A(60,20) = -pK.y()/EK;
    A(60,21) = -pK.z()/EK;
    */

    // A.Print();

    TMatrixDSym A_sym(dimM+dimX+dimC);
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        for(int j = 0; j < dimM+dimX+dimC; j++)
        {
            A_sym(i,j) = A(i,j);
        }
    }

    // cout << "|A| = " << A_sym.Determinant() << endl;

    // TVectorD A_eigen_values(dimM+dimX+dimC);  
    // TMatrixD A_eigen_vectors = A_sym.EigenVectors(A_eigen_values);
    // A_eigen_values.Print();

    TMatrixD A_inverse = A_sym;
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

double equations( const double* x_vars )
{
    TVectorD x(dimM+dimX+dimC);
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        x(i) = x_vars[i];
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

    // Unknown parameters (3)
    Double_t BVz = x(dimM);
    Double_t ptau1z = x(dimM+1);
    Double_t ptau2z = x(dimM+2);

    // BV must lie in K+ trajectory
    Double_t BVx = RP.x() + (pK.x()/pK.z())*(BVz - RP.z());
    Double_t BVy = RP.y() + (pK.y()/pK.z())*(BVz - RP.z()); 
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz );
    // ptau1 must point back to the BV
    Double_t ptau1x = ((DV1.x() - BV.x())/(DV1.z() - BV.z()))*ptau1z;
    Double_t ptau1y = ((DV1.y() - BV.y())/(DV1.z() - BV.z()))*ptau1z;
    ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );
    // ptau2 must point back to the BV
    Double_t ptau2x = ((DV2.x() - BV.x())/(DV2.z() - BV.z()))*ptau2z;
    Double_t ptau2y = ((DV2.y() - BV.y())/(DV2.z() - BV.z()))*ptau2z;
    ROOT::Math::XYZVector ptau2( ptau2x, ptau2y, ptau2z );
    // Tau mass constraints
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    // 3-momentum conservation in DV1 and in DV2
    ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
    ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
    // Neutrino mass constraints
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );
    // 4-momentum conservation in BV
    ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
    Double_t EB = Etau1 + Etau2 + EK;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    // Lagrange multipliers 
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
    eqs(0) = -2.*chi2_sum(0) + l(0)*pB.z();
    eqs(1) = -2.*chi2_sum(1) + l(1)*pB.z();
    eqs(2) = -2.*chi2_sum(2) - l(0)*pB.x() - l(1)*pB.y();
    // DV1
    eqs(3) = -2.*chi2_sum(3) + l(0)*((ptau1.z()*(BV.z()-PV.z()))/(DV1.z()-BV.z())) + l(2)*(ptau1.z()/(DV1.z()-BV.z()))*(-pnu1.x()/Enu1 + ptau1.x()/Etau1);
    eqs(4) = -2.*chi2_sum(4) + l(1)*((ptau1.z()*(BV.z()-PV.z()))/(DV1.z()-BV.z())) + l(2)*(ptau1.z()/(DV1.z()-BV.z()))*(-pnu1.y()/Enu1 + ptau1.y()/Etau1);
    eqs(5) = -2.*chi2_sum(5) - l(0)*((ptau1.x()*(BV.z()-PV.z()))/(DV1.z()-BV.z())) - l(1)*((ptau1.y()*(BV.z()-PV.z()))/(DV1.z()-BV.z())) + l(2)*(1./(DV1.z()-BV.z()))*((ptau1.x()*pnu1.x() + ptau1.y()*pnu1.y())/Enu1 - (ptau1.x()*ptau1.x() + ptau1.y()*ptau1.y())/Etau1);
    // P3pi1
    eqs(6) = -2.*chi2_sum(6) + l(2)*(pnu1.x()/Enu1);
    eqs(7) = -2.*chi2_sum(7) + l(2)*(pnu1.y()/Enu1);
    eqs(8) = -2.*chi2_sum(8) + l(2)*(pnu1.z()/Enu1);
    eqs(9) = -2.*chi2_sum(9) - l(2);
    // DV2
    eqs(10) = -2.*chi2_sum(10) + l(0)*((ptau2.z()*(BV.z()-PV.z()))/(DV2.z()-BV.z())) + l(3)*(ptau2.z()/(DV2.z()-BV.z()))*(-pnu2.x()/Enu2 + ptau2.x()/Etau2);
    eqs(11) = -2.*chi2_sum(11) + l(1)*((ptau2.z()*(BV.z()-PV.z()))/(DV2.z()-BV.z())) + l(3)*(ptau2.z()/(DV2.z()-BV.z()))*(-pnu2.y()/Enu2 + ptau2.y()/Etau2);
    eqs(12) = -2.*chi2_sum(12) - l(0)*((ptau2.x()*(BV.z()-PV.z()))/(DV2.z()-BV.z())) - l(1)*((ptau2.y()*(BV.z()-PV.z()))/(DV2.z()-BV.z())) + l(3)*(1./(DV2.z()-BV.z()))*((ptau2.x()*pnu2.x() + ptau2.y()*pnu2.y())/Enu2 - (ptau2.x()*ptau2.x() + ptau2.y()*ptau2.y())/Etau2);
    // P3pi2
    eqs(13) = -2.*chi2_sum(13) + l(3)*(pnu2.x()/Enu2);
    eqs(14) = -2.*chi2_sum(14) + l(3)*(pnu2.y()/Enu2);
    eqs(15) = -2.*chi2_sum(15) + l(3)*(pnu2.z()/Enu2);
    eqs(16) = -2.*chi2_sum(16) - l(3);
    // RP_T
    eqs(17) = -2.*chi2_sum(17) + l(0)*(-pB.z() - (BV.z()-PV.z())*(ptau1.z()/(DV1.z()-BV.z()) + ptau2.z()/(DV2.z()-BV.z()))) + l(2)*(ptau1.z()/(DV1.z()-BV.z()))*(pnu1.x()/Enu1 - ptau1.x()/Etau1) + l(3)*(ptau2.z()/(DV2.z()-BV.z()))*(pnu2.x()/Enu2 - ptau2.x()/Etau2);
    eqs(18) = -2.*chi2_sum(18) + l(1)*(-pB.z() - (BV.z()-PV.z())*(ptau1.z()/(DV1.z()-BV.z()) + ptau2.z()/(DV2.z()-BV.z()))) + l(2)*(ptau1.z()/(DV1.z()-BV.z()))*(pnu1.y()/Enu1 - ptau1.y()/Etau1) + l(3)*(ptau2.z()/(DV2.z()-BV.z()))*(pnu2.y()/Enu2 - ptau2.y()/Etau2);
    // pK
    eqs(19) = -2.*chi2_sum(19) + l(0)*( (BV.z()-PV.z())*(1. - (ptau1.z()/pK.z())*((BV.z()-RP.z())/(DV1.z()-BV.z())) - (ptau2.z()/pK.z())*((BV.z()-RP.z())/(DV2.z()-BV.z()))) - (BV.z()-RP.z())*(pB.z()/pK.z()) ) + l(2)*(ptau1.z()/pK.z())*((BV.z()-RP.z())/(DV1.z()-BV.z()))*(pnu1.x()/Enu1 - ptau1.x()/Etau1) + l(3)*(ptau2.z()/pK.z())*((BV.z()-RP.z())/(DV2.z()-BV.z()))*(pnu2.x()/Enu2 - ptau2.x()/Etau2);
    eqs(20) = -2.*chi2_sum(20) + l(1)*( (BV.z()-PV.z())*(1. - (ptau1.z()/pK.z())*((BV.z()-RP.z())/(DV1.z()-BV.z())) - (ptau2.z()/pK.z())*((BV.z()-RP.z())/(DV2.z()-BV.z()))) - (BV.z()-RP.z())*(pB.z()/pK.z()) ) + l(2)*(ptau1.z()/pK.z())*((BV.z()-RP.z())/(DV1.z()-BV.z()))*(pnu1.y()/Enu1 - ptau1.y()/Etau1) + l(3)*(ptau2.z()/pK.z())*((BV.z()-RP.z())/(DV2.z()-BV.z()))*(pnu2.y()/Enu2 - ptau2.y()/Etau2);
    eqs(21) = -2.*chi2_sum(21) + l(0)*( PV.x()-RP.x() + (BV.z()-PV.z())*(BV.z()-RP.z())*(pK.x()/pow(pK.z(),2))*(ptau1.z()/(DV1.z()-BV.z()) + ptau2.z()/(DV2.z()-BV.z())) - (pK.x()/pK.z())*(BV.z()-RP.z()) + (pK.x()*pB.z()/pow(pK.z(),2))*(BV.z()-RP.z()) ) + l(1)*( PV.y()-RP.y() + (BV.z()-PV.z())*(BV.z()-RP.z())*(pK.y()/pow(pK.z(),2))*(ptau1.z()/(DV1.z()-BV.z()) + ptau2.z()/(DV2.z()-BV.z())) - (pK.y()/pK.z())*(BV.z()-RP.z()) + (pK.y()*pB.z()/pow(pK.z(),2))*(BV.z()-RP.z())  ) + l(2)*(ptau1.z()/pow(pK.z(),2))*((BV.z()-RP.z())/(DV1.z()-BV.z()))*( -(pK.x()*pnu1.x() + pK.y()*pnu1.y())/Enu1 + (pK.x()*ptau1.x() + pK.y()*ptau1.y())/Etau1 ) + l(3)*(ptau2.z()/pow(pK.z(),2))*((BV.z()-RP.z())/(DV2.z()-BV.z()))*( -(pK.x()*pnu2.x() + pK.y()*pnu2.y())/Enu2 + (pK.x()*ptau2.x() + pK.y()*ptau2.y())/Etau2 );                                                                     
    // BVz
    eqs(dimM) = l(0)*( pB.x() - (pK.x()*pB.z()/pK.z()) + (BV.z()-PV.z())*( -(pK.x()*ptau1.z()/(pK.z()*(DV1.z()-BV.z()))) - (pK.x()*ptau2.z()/(pK.z()*(DV2.z()-BV.z()))) + (ptau1.x()/(DV1.z()-BV.z())) + (ptau2.x()/(DV2.z()-BV.z())) ) ) + l(1)*( pB.y() - (pK.y()*pB.z()/pK.z()) + (BV.z()-PV.z())*( -(pK.y()*ptau1.z()/(pK.z()*(DV1.z()-BV.z()))) - (pK.y()*ptau2.z()/(pK.z()*(DV2.z()-BV.z()))) + (ptau1.y()/(DV1.z()-BV.z())) + (ptau2.y()/(DV2.z()-BV.z())) ) ) + l(2)*( (-(pK.x()*ptau1.z())/(pK.z()*(DV1.z()-BV.z())) + ptau1.x()/(DV1.z()-BV.z()) )*(ptau1.x()/Etau1 - pnu1.x()/Enu1) + (-(pK.y()*ptau1.z())/(pK.z()*(DV1.z()-BV.z())) + ptau1.y()/(DV1.z()-BV.z()))*(ptau1.y()/Etau1 - pnu1.y()/Enu1) ) + l(3)*( (-(pK.x()*ptau2.z())/(pK.z()*(DV2.z()-BV.z())) + ptau2.x()/(DV2.z()-BV.z()))*(ptau2.x()/Etau2 - pnu2.x()/Enu2) + (-(pK.y()*ptau2.z())/(pK.z()*(DV2.z()-BV.z())) + ptau2.y()/(DV2.z()-BV.z()))*(ptau2.y()/Etau2 - pnu2.y()/Enu2) );
    // ptau1z
    eqs(dimM+1) = l(0)*( PV.x()-RP.x() + (BV.z()-PV.z())*(DV1.x()-BV.x())/(DV1.z()-BV.z()) - (pK.x()/pK.z())*(BV.z()-RP.z()) ) + l(1)*( PV.y()-RP.y() + (BV.z()-PV.z())*(DV1.y()-BV.y())/(DV1.z()-BV.z()) - (pK.y()/pK.z())*(BV.z()-RP.z()) ) + l(2)*( -(1./Enu1)*( pnu1.z() + pnu1.x()*((DV1.x()-BV.x())/(DV1.z()-BV.z())) + pnu1.y()*((DV1.y()-BV.y())/(DV1.z()-BV.z())) ) + (1./Etau1)*( ptau1.z() + ptau1.x()*((DV1.x()-BV.x())/(DV1.z()-BV.z())) + ptau1.y()*((DV1.y()-BV.y())/(DV1.z()-BV.z())) ) );
    // ptau2z
    eqs(dimM+2) = l(0)*( PV.x()-RP.x() + (BV.z()-PV.z())*(DV2.x()-BV.x())/(DV2.z()-BV.z()) - (pK.x()/pK.z())*(BV.z()-RP.z()) ) + l(1)*( PV.y()-RP.y() + (BV.z()-PV.z())*(DV2.y()-BV.y())/(DV2.z()-BV.z()) - (pK.y()/pK.z())*(BV.z()-RP.z()) ) + l(3)*( -(1./Enu2)*( pnu2.z() + pnu2.x()*((DV2.x()-BV.x())/(DV2.z()-BV.z())) + pnu2.y()*((DV2.y()-BV.y())/(DV2.z()-BV.z())) ) + (1./Etau2)*( ptau2.z() + ptau2.x()*((DV2.x()-BV.x())/(DV2.z()-BV.z())) + ptau2.y()*((DV2.y()-BV.y())/(DV2.z()-BV.z())) ) );
    // pB must point back to the PV (2)
    eqs(dimM+dimX) = pB.x()*( BV.z() - PV.z() ) - pB.z()*( BV.x() - PV.x() );
    eqs(dimM+dimX+1) = pB.y()*( BV.z() - PV.z() ) - pB.z()*( BV.y() - PV.y() );
    // Energy conservation in DV1 and in DV2
    eqs(dimM+dimX+2) = Etau1 - E3pi1 - Enu1;
    eqs(dimM+dimX+3) = Etau2 - E3pi2 - Enu2;

    // normalise equations
    // eqs = normalise(eqs, true);
    // eqs = D*eqs;
    // eqs.Print();
 
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

TVectorD normalise( TVectorD a, bool norm )
{
    if(norm) // normalise
    {
        for(int i = 0; i < dimM+dimX+dimC; i++)
        {
            a(i) = a(i)/pow(10,int(log10(abs(a(i)))));
        }
    }
    else // un-normalise
    {
        for(int i = 0; i < dimM+dimX+dimC; i++)
        {
            a(i) = a(i)*pow(10,int(log10(abs(a(i)))));
        }
    }
    return a;
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

    // Unknown parameters (3)
    Double_t BVz = x(dimM);
    Double_t ptau1z = x(dimM+1);
    Double_t ptau2z = x(dimM+2);

    // BV must lie in K+ trajectory
    Double_t BVx = RP.x() + (pK.x()/pK.z())*(BVz - RP.z());
    Double_t BVy = RP.y() + (pK.y()/pK.z())*(BVz - RP.z()); 
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz );
    // ptau1 must point back to the BV
    Double_t ptau1x = ((DV1.x() - BV.x())/(DV1.z() - BV.z()))*ptau1z;
    Double_t ptau1y = ((DV1.y() - BV.y())/(DV1.z() - BV.z()))*ptau1z;
    ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );
    // ptau2 must point back to the BV
    Double_t ptau2x = ((DV2.x() - BV.x())/(DV2.z() - BV.z()))*ptau2z;
    Double_t ptau2y = ((DV2.y() - BV.y())/(DV2.z() - BV.z()))*ptau2z;
    ROOT::Math::XYZVector ptau2( ptau2x, ptau2y, ptau2z );
    // Tau mass constraints
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    // 3-momentum conservation in DV1 and in DV2
    ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
    ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
    // Neutrino mass constraints
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );
    // 4-momentum conservation in BV
    ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
    Double_t EB = Etau1 + Etau2 + EK;

    TVectorD g(dimC);
    // pB must point back to the PV
    g(0) = pB.x()*( BV.z() - PV.z() ) - pB.z()*( BV.x() - PV.x() );
    g(1) = pB.y()*( BV.z() - PV.z() ) - pB.z()*( BV.y() - PV.y() );
    // Energy conservation in DV1 and in DV2
    g(2) = Etau1 - E3pi1 - Enu1;
    g(3) = Etau2 - E3pi2 - Enu2;

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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    return x0;
}

TVectorD x_initial_estimate_1( TVectorD m, ROOT::Math::XYZPoint BV ) 
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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    return x0;
}

TVectorD x_initial_estimate_2( TVectorD m, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; taus direction based visible 3pi momenta
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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    return x0;
}

TVectorD x_initial_estimate_3( TVectorD m, ROOT::Math::XYZPoint BV ) 
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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    // normalise X
    // x0 = normalise(x0, true);
    // x0 = D_inverse*x0;

    return x0;
}

TVectorD x_initial_estimate_4( TVectorD m, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; tau+ direction from vertices ; tau- direction from pions
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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    return x0;
}

TVectorD x_initial_estimate_5( TVectorD m, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; tau- direction from vertices ; tau+ direction from pions
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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

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

    ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
    ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

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

    ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
    ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

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

    ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
    ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

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

    ROOT::Math::XYZVector u1 = (DV1 - BV).Unit();
    ROOT::Math::XYZVector u2 = (DV2 - BV).Unit();

    Double_t EK = sqrt( pow(mkaon,2) + pK.Mag2() );

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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

    return x0;
}

TVectorD x_initial_estimate_true( TVectorD m )  // Unmeasured quantities are initialised using true values
{  
    TVectorD x0(dimM+dimX+dimC);

    // 1) Initialise xm
    for(int i = 0; i < dimM; i++)
    {
        x0(i) = m(i);
    }

    // 2) Initialise xu
    x0(dimM) = BVz_true;
    x0(dimM+1) = ptau1z_true;
    x0(dimM+2) = ptau2z_true;

    // lambda
    x0(dimM+dimX) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1./10000.; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1./10000.; // E conservation in DV1
    x0(dimM+dimX+3) = 1./10000.; // E conservation in DV2

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