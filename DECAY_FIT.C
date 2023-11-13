using namespace std;

// Global variables 
int dimM = 22; // number of measured parameters
int dimX = 19; // number of unknown parameters
int dimC = 20; // number of exact constraints
int eq_flag = 0; // index of the equation

TVectorD m(dimM);
TMatrixDSym V(dimM);
TMatrixDSym W(dimM);
TVectorD xm(dimM);
TVectorD xu(dimX);
TVectorD lambda(dimC);

TVectorD x0(dimM+dimX+dimC);
Double_t RPz; // z-component of the RP on the K+ trajectory (fixed)

// Fit results saved
Double_t MB, MB_err;
TVectorD X(dimM+dimX+dimC);
TVectorD X_ERR(dimM+dimX+dimC);
TVectorD F(dimM+dimX+dimC);
TMatrixDSym C(dimM+dimX+dimC);
TMatrixDSym C_corr(dimM+dimX+dimC);
Int_t status, init;

// Constants
Double_t mtau = 1776.86;
Double_t mkaon = 493.677;

// Functions
void solve( ROOT::Math::XYZPoint BV, int init,  ROOT::Math::GSLMultiRootFinder* solver );
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);
TVectorD x_initial_estimate_0( TVectorD m );
TVectorD x_initial_estimate_1( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_2( TVectorD m, ROOT::Math::XYZPoint BV );
TVectorD x_initial_estimate_3( TVectorD m, ROOT::Math::XYZPoint BV );
Double_t lagrangian( const Double_t* x_values );
Double_t chisquare( TVectorD xm );
TVectorD exact_constraints( TVectorD x );
TMatrixDSym U_matrix();
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

void DECAY_FIT(int year, TString MC_files, TString RS_DATA_files, TString WS_DATA_files, int species, int component, int line)
{
    // Retrieve m and V from ntuple
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

    // Save result of the fit (xm,xu) in a .root file (x = (xm,xu))
    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/Component_%i/fit_result_%i.root",year,species,component,line), "RECREATE");
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
        "df_antinutau_PX",
        "df_antinutau_PY",
        "df_antinutau_PZ",
        "df_taum_PX",
        "df_taum_PY",
        "df_taum_PZ",
        "df_nutau_PX",
        "df_nutau_PY",
        "df_nutau_PZ"
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
        "df_antinutau_PX_err",
        "df_antinutau_PY_err",
        "df_antinutau_PZ_err",
        "df_taum_PX_err",
        "df_taum_PY_err",
        "df_taum_PZ_err",
        "df_nutau_PX_err",
        "df_nutau_PY_err",
        "df_nutau_PZ_err"
    };

    for(int i = 0; i < dimM+dimX; i++)
    {
        tout->Branch(name_x[i], &X(i));
        // tout->Branch(name_x_err[i], &X_ERR(i));
        // for(int j = 0; j < dimX; j++)
        // {
        //     tout->Branch(Form("df_C_%i_%i",i,j), &C(i,j));
        // }
    }

    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        tout->Branch(Form("df_F_%i",i), &F(i));
    }

    // tout->Branch("df_init", &init);
    tout->Branch("df_Bp_M", &MB);
    // tout->Branch("df_Bp_MERR", &MB_err);
    tout->Branch("df_status", &status);

    // Loop over events
    for(int evt = 0; evt < 1; evt++)
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
        solver->SetPrintLevel(1);

        // 4) Solve the system of euqations (61)
        // init = 0:  original calculations = Anne Keune (29%)
        // init = 1: K*tautau calculations; tau momentum direction initialised based on vertices (uses offline estimate for BV) (13%)
        // init = 2: // K*tautau calculations; tau momentum direction initialised based on 3pi visible momentum (uses offline estimate for BV) (8%)
        // init = 3: RD calculations (uses offline estimate for BV) (58%)

        ROOT::Math::XYZPoint BV( BVx, BVy, BVz );

        // 3102
        solve( BV, 3, solver );
        if(status != 0)
        {
            solve( BV, 1, solver );
        }
        if(status != 0)
        {
            solve( BV, 0, solver );
        }
        if(status != 0)
        {
            solve( BV, 2, solver );
        }

        // 1230 (62%)
        // solve( BV, 1, solver );
        // if(status != 0)
        // {
        //     solve( BV, 2, solver );
        // }
        // if(status != 0)
        // {
        //     solve( BV, 3, solver );
        // }
        // if(status != 0)
        // {
        //     solve( BV, 0, solver );
        // }

        // 3012 (38%)
        // solve( BV, 3, solver );
        // if(status != 0)
        // {
        //     solve( BV, 0, solver );
        // }
        // if(status != 0)
        // {
        //     solve( BV, 1, solver );
        // }
        // if(status != 0)
        // {
        //     solve( BV, 2, solver );
        // }

        tout->Fill();
    }
    cout << "FINISHED" << endl;

    fout->cd();
    tout->Write();
    fout->Close();
}

void solve( ROOT::Math::XYZPoint BV, int init,  ROOT::Math::GSLMultiRootFinder* solver )
{
    // 1) This function initialises the vector of unkown parameters x with the result of analytical calculations
    // 2) It sets up the minimizer: defining the initial values and bounds on the x parameters  
    // 3) It builds the chi^2 function with x0 by calling the function chisquare 
    // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

    // 1) Initial values for x=(xm,xu,lambda)
    if(init == 0)
    {
        x0 = x_initial_estimate_0( m ); // original calculations = Anne Keune
    }
    else if(init == 1)
    {   
        x0 = x_initial_estimate_1( m, BV ); // K*tautau calculations; tau momentum direction initialised based on vertices (uses offline estimate for BV)
    }
    else if(init == 2)
    {
        x0 = x_initial_estimate_2( m, BV ); // K*tautau calculations; tau momentum direction initialised based on 3pi visible momentum (uses offline estimate for BV) 
    }
    else if(init == 3)
    {
        x0 = x_initial_estimate_3( m, BV ); // RD calculations (uses offline estimate for BV)
    }
    // x0.Print();

    Double_t x0_vars[dimM+dimX+dimC], x0_err[dimM+dimX+dimC];
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        x0_vars[i] = x0(i); // x0_vars is a Double_t, x0 is a TVectorD; the function to minimise must receive a Double_t as input
        x0_err[i] = 0.1*abs(x0_vars[i]); // initial errors on x0; they are used as the first step size 
    }

    // 2) Solve system of equations
    solver->Solve(x0_vars, 10000, 1);

    // 3) Retrieve results
    const double* x_results = solver->X();
    const double* f_vals = solver->FVal();
    status = solver->Status();
    cout << "status = " << status << endl;

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
    cout << "MB = " << MB << endl;

    TMatrixDSym U = U_matrix();

}

TMatrixDSym U_matrix()
{
    TMatrixDSym U_cov(dimM+dimX+dimC);

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
    ROOT::Math::XYZVector pnu1( X(dimM+10), X(dimM+11), X(dimM+12) );
    ROOT::Math::XYZVector ptau2( X(dimM+13), X(dimM+14), X(dimM+15) );
    ROOT::Math::XYZVector pnu2( X(dimM+16), X(dimM+17), X(dimM+18) );

    Double_t EB = sqrt( MB_squared + pB.Mag2() );
    // tau and neutrino mass constraints are applied by simple parameter substitution
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );

    // Lagrange multipliers (20)
    TVectorD l(dimC);
    for(int i = 0; i < dimC;i++)
    {
        l(i) = X(dimM+dimX+i);
    }

    int L = dimM+dimX;

    ///////////////////////////// xm ///////////////////////////
    A(0,0) = 2*W(0,0) + (2*pB.x()/pow(BV.x()-PV.x(),3))*(l(0)+l(1));
    for(int i = 0; i < dimM; i++)
    {
        if(i != 0)
        {
            A(0,i) = 2*W(0,i);
        }
    }
    A(0,22) = -(2*pB.x()/pow(BV.x()-PV.x(),3))*(l(0) + l(1));
    A(0,25) = (l(0) + l(1))/pow(BV.x()-PV.x(),2);
    A(0,L) = pB.x()/pow(BV.x()-pB.x(),2);
    A(0,L+1) = pB.x()/pow(BV.x()-pB.x(),2);

    A(1,1) = 2*W(1,1) - (2*pB.y()/pow(BV.y()-PV.y(),3))*l(0);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 1)
        {
            A(1,i) = 2*W(1,i);
        }
    }
    A(1,23) = (2*pB.y()/pow(BV.y()-PV.y(),3))*l(0);
    A(1,L) = -pB.y()/pow(BV.y()-pB.y(),2);

    A(2,2) = 2*W(2,2) - (2*pB.z()/pow(BV.z()-PV.z(),3))*l(1);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 2)
        {
            A(2,i) = 2*W(2,i);
        }
    }
    A(2,24) = (2*pB.z()/pow(BV.z()-PV.z(),3))*l(1);
    A(2,L+1) = -pB.z()/pow(BV.z()-PV.z(),2);

    A(3,3) = 2*W(3,3) + (2*ptau1.x()/pow(DV1.x()-BV.x(),3))*(l(2) + l(3));
    for(int i = 0; i < dimM; i++)
    {
        if(i != 3)
        {
            A(3,i) = 2*W(3,i);
        }
    }
    A(3,22) = -(2*ptau1.x()/pow(DV1.x()-BV.x(),3))*(l(2) + l(3));
    A(3,29) = -(l(2) + l(3))/pow(DV1.x()-BV.x(),2);
    A(3,L+2) = -ptau1.x()/pow(DV1.x()-BV.x(),2);
    A(3,L+3) = -ptau1.x()/pow(DV1.x()-BV.x(),2);

    A(4,4) = 2*W(4,4) - (2*ptau1.y()/pow(DV1.y()-BV.y(),3))*l(2);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 4)
        {
            A(4,i) = 2*W(4,i);
        }
    }
    A(4,23) = (2*ptau1.y()/pow(DV1.y()-BV.y(),3))*l(2);
    A(4,30) = l(2)/pow(DV1.y()-BV.y(),2);
    A(4,L+2) = ptau1.y()/pow(DV1.y()-BV.y(),2);

    A(5,5) = 2*W(5,5) - (2*ptau1.z()/pow(DV1.z()-BV.z(),3))*l(3);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 5)
        {
            A(5,i) = 2*W(5,i);
        }
    }
    A(5,24) = (2*ptau1.z()/pow(DV1.z()-BV.z(),3))*l(3);
    A(5,31) = l(3)/pow(DV1.z()-BV.z(),2);
    A(5,L+3) = ptau1.z()/pow(DV1.z()-BV.z(),2);

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

    A(10,10) = 2*W(10,10) + (2*ptau2.x()/pow(DV2.x()-BV.x(),3))*(l(8) + l(9));
    for(int i = 0; i < dimM; i++)
    {
        if(i != 10)
        {
            A(10,i) = 2*W(10,i);
        }
    }
    A(10,22) = -(2*ptau2.x()/pow(DV2.x()-BV.x(),3))*(l(8) + l(9));
    A(10,35) = -(l(8) + l(9))/pow(DV2.x()-BV.x(),2);
    A(10,L+8) = -ptau2.x()/pow(DV2.x()-BV.x(),2);
    A(10,L+9) = -ptau2.x()/pow(DV2.x()-BV.x(),2);

    A(11,11) = 2*W(11,11) - (2*ptau2.y()/pow(DV2.y()-BV.y(),3))*l(8);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 11)
        {
            A(11,i) = 2*W(11,i);
        }
    }
    A(11,23) = (2*ptau2.y()/pow(DV2.y()-BV.y(),3))*l(8);
    A(11,36) = l(8)/pow(DV2.y()-BV.y(),2);
    A(11,L+8) = ptau2.y()/pow(DV2.y()-BV.y(),2);

    A(12,12) = 2*W(12,12) - (2*ptau2.z()/pow(DV2.z()-BV.z(),3))*l(9);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 12)
        {
            A(12,i) = 2*W(12,i);
        }
    }
    A(12,24) = (2*ptau2.z()/pow(DV2.z()-BV.z(),3))*l(9);
    A(12,37) = l(9)/pow(DV2.z()-BV.z(),2);
    A(12,L+9) = ptau2.z()/pow(DV2.z()-BV.z(),2);

    for(int i = 0; i < dimM; i++)
    {
        A(13,i) = 2*W(13,i);
        A(14,i) = 2*W(14,i);
        A(15,i) = 2*W(15,i);
        A(16,i) = 2*W(16,i);
    }
    A(13,L+10) = -1;
    A(14,L+11) = -1;
    A(15,L+12) = -1;
    A(16,L+13) = -1;

    A(17,17) = 2*W(17,17) + (2*pK.x()/pow(BV.x()-RP.x(),3))*(l(14) + l(15));
    for(int i = 0; i < dimM; i++)
    {
        if((i != 17) && (i != 19))
        {
            A(17,i) = 2*W(17,i);
        }
    }
    A(17,22) = -(2*pK.x()/pow(BV.x()-RP.x(),3))*(l(14) + l(15));
    A(17,19) = 2*W(17,19) + (l(14) + l(15))/pow(BV.x()-RP.x(),2);
    A(17,L+14) = pK.x()/pow(BV.x()-RP.x(),2);
    A(17,L+15) = pK.x()/pow(BV.x()-RP.x(),2);

    A(18,18) = 2*W(18,18) - (2*pK.y()/pow(BV.y()-RP.y(),3))*l(14);
    for(int i = 0; i < dimM; i++)
    {
        if((i != 18) && (i != 20))
        {
            A(18,i) = 2*W(18,i);
        }
    }
    A(18,23) = (2*pK.y()/pow(BV.y()-RP.y(),3))*l(14);
    A(18,20) = 2*W(18,20) - l(14)/pow(BV.y()-RP.y(),2);
    A(18,L+14) = pK.y()/pow(BV.y()-RP.y(),2);

    A(19,19) = 2*W(19,19) - ((pow(EK,2) - pow(pK.x(),2))/pow(EK,3))*l(19);
    for(int i = 0; i < dimM; i++)
    {
        if((i != 19) && (i != 17))
        {
            A(19,i) = 2*W(19,i);
        }
    }
    A(19,17) = 2*W(19,17) + (l(14) + l(15))/pow(BV.x()-RP.x(),2);
    A(19,22) = -(l(14) + l(15))/pow(BV.x()-RP.x(),2);
    A(19,L+14) = 1/(BV.x()-RP.x());
    A(19,L+15) = 1/(BV.x()-RP.x());
    A(19,L+16) = -1;
    A(19,L+19) = -pK.x()/EK;

    A(20,20) = 2*W(20,20) - ((pow(EK,2) - pow(pK.y(),2))/pow(EK,3))*l(19);
    for(int i = 0; i < dimM; i++)
    {
        if((i != 20) && (i != 18))
        {
            A(20,i) = 2*W(20,i);
        }
    }
    A(20,18) = 2*W(20,18) - l(14)/pow(BV.y()-RP.y(),2);
    A(20,23) = l(14)/pow(BV.y()-RP.y(),2);
    A(20,L+14) = -1/(BV.y()-RP.y());
    A(20,L+17) = -1;
    A(20,L+19) = -pK.y()/EK;

    A(21,21) = 2*W(21,21) - ((pow(EK,2) - pow(pK.z(),2))/pow(EK,3))*l(19);
    for(int i = 0; i < dimM; i++)
    {
        if(i != 21)
        {
            A(21,i) = 2*W(21,i);
        }
    }
    A(21,24) = l(15)/pow(BV.z()-RP.z(),2);
    A(21,L+15) = -1/(BV.z()-RP.z());
    A(21,L+18) = -1;
    A(21,L+19) = -pK.z()/EK;

    //////////////////////////////// xu /////////////////////////
    A(22,0) = -(2*pB.x()/pow(BV.x()-PV.x(),3))*(l(0)+l(1));
    A(22,3) = -(2*ptau1.x()/pow(DV1.x()-BV.x(),3))*(L(2)+l(3));
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
    A(23,18) = -(2*pK.y()/pow(BV.y()-RP.y(),3))*l(14);
    A(23,20) = l(14)/pow(BV.y()-RP.y(),2);
    A(23,23) = -(2*pB.y()/pow(BV.y()-RP.y(),3))*l(0) - (2*ptau1.y()/pow(DV1.y()-BV.y(),3))*l(2) - (2*ptau2.y()/pow(DV2.y()-BV.y(),3))*l(8) - (2*pK.y()/pow(BV.y()-RP.y(),3))*l(14);
    A(23,26) = l(0)/pow(BV.y()-PV.y(),2);
    A(23,30) = -l(2)/pow(DV1.y()-BV.y(),2);
    A(23,36) = -l(8)/pow(DV2.y()-BV.y(),2);
    A(23,L) = pB.y()/pow(BV.y()-PV.y(),2);
    A(23,L+2) = -ptau1.y()/pow(DV1.y()-BV.y(),2);
    A(23,L+8) = -ptau2.y()/pow(DV2.y()-BV.y(),2);
    A(23,L+14) = pK.y()/pow(BV.y()-RP.y(),2);

    A(24,2) = (2*pB.z()/pow(BV.z()-RP.z(),3))*l(1);
    A(24,5) = (2*ptau1.z()/pow(DV1.z()-BV.z(),3))*l(3);
    A(24,12) = (2*ptau2.z()/pow(DV2.z()-BV.z(),3))*l(9);
    A(24,19) = l(15)/pow(BV.z()-RP.z(),2);
    A(24,24) = -(2*pB.z()/pow(BV.z()-RP.z(),3))*l(1) - (2*ptau1.z()/pow(DV1.z()-BV.z(),3))*l(3) - (2*ptau2.z()/pow(DV2.z()-BV.z(),3))*l(9) - (2*pK.z()/pow(BV.z()-RP.z(),3))*l(15);
    A(24,27) = l(1)/pow(BV.z()-PV.z(),2);
    A(24,31) = -l(3)/pow(DV1.z()-BV.z(),2);
    A(24,37) = -l(9)/pow(DV2.z()-BV.z(),2);
    A(24,L+1) = pB.z()/pow(BV.z()-PV.z(),2);
    A(24,L+3) = -ptau1.z()/pow(DV1.z()-BV.z(),2);
    A(24,L+9) = -ptau2.z()/pow(DV2.z()-BV.z(),2);
    A(24,L+15) = pK.z()/pow(BV.z()-RP.z(),2);

    A(25,0) = (l(0)+l(1))/pow(BV.x()-PV.x(),2);
    A(25,22) = -(l(0)+l(1))/pow(BV.x()-PV.x(),2);
    A(25,L) = 1/(BV.x()-PV.x());
    A(25,L+1) = 1/(BV.x()-PV.x());
    A(25,L+16) = 1;

    A(26,1) = -l(0)/pow(BV.y()-PV.y(),2);
    A(26,23) = l(0)/pow(BV.y()-PV.y(),2);
    A(26,L) = -1/(BV.y()-PV.y());
    A(26,L+17) = 1;

    A(27,2) = -l(1)/pow(BV.z()-PV.z(),2);
    A(27,24) = l(1)/pow(BV.z()-PV.z(),2);
    A(27,L+1) = -1/(BV.z()-PV.z());
    A(27,L+18) = 1;

    A(28,28) = -l(19)/(2*pow(EB,3));

    A(29,3) = -(l(2)+l(3))/pow(DV1.x()-BV.x(),2);
    A(29,22) = (l(2)+l(3))/pow(DV1.x()-BV.x(),2);
    A(29,29) = ((pow(Etau1,2) - pow(ptau1.x(),2))/pow(Etau1,3))*(l(7)-l(19));
    A(29,L+2) = 1/(DV1.x()-BV.x());
    A(29,L+3) = 1/(DV1.x()-BV.x());
    A(29,L+4) = 1;
    A(29,L+7) = ptau1.x()/Etau1;
    A(29,L+16) = -1;
    A(29,L+19) = -ptau1.x()/Etau1;

    A(30,4) = l(2)/pow(DV1.y()-BV.y(),2);
    A(30,23) = -l(2)/pow(DV1.y()-BV.y(),2);
    A(30,30) = ((pow(Etau1,2) - pow(ptau1.y(),2))/pow(Etau1,3))*(l(7)-l(19));
    A(30,L+2) = -1/(DV1.y()-BV.y());
    A(30,L+5) = 1;
    A(30,L+7) = ptau1.y()/Etau1;
    A(30,L+17) = -1;
    A(30,L+19) = -ptau1.y()/Etau1;

    A(31,5) = l(3)/pow(DV1.z()-BV.z(),2);
    A(31,24) = -l(3)/pow(DV1.z()-BV.z(),2);
    A(31,31) = ((pow(Etau1,2) - pow(ptau1.z(),2))/pow(Etau1,3))*(l(7)-l(19));
    A(31,L+3) = -1/(DV1.z()-BV.z());
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
    A(35,L+8) = 1/(DV2.x()-BV.x());
    A(35,L+9) = 1/(DV2.x()-BV.x());
    A(35,L+10) = 1;
    A(35,L+13) = ptau2.x()/Etau2;
    A(35,L+16) = -1;
    A(35,L+19) = -ptau2.x()/Etau2;

    A(36,11) = l(8)/pow(DV2.y()-BV.y(),2);
    A(36,23) = -l(8)/pow(DV2.y()-BV.y(),2);
    A(36,36) = ((pow(Etau2,2) - pow(ptau2.y(),2))/pow(Etau2,3))*(l(13)-l(19));
    A(36,L+8) = -1/(DV2.y()-BV.y());
    A(36,L+11) = 1;
    A(36,L+13) = ptau2.y()/Etau2;
    A(36,L+17) = -1;
    A(36,L+19) = -ptau2.y()/Etau2;

    A(37,12) = l(9)/pow(DV2.z()-BV.z(),2);
    A(37,24) = -l(9)/pow(DV2.z()-BV.z(),2);
    A(37,37) = ((pow(Etau2,2) - pow(ptau2.z(),2))/pow(Etau2,3))*(l(13)-l(19));
    A(37,L+9) = -1/(DV2.z()-BV.z());
    A(37,L+12) = 1;
    A(37,L+13) = ptau2.z()/Etau2;
    A(37,L+18) = -1;
    A(37,L+19) = -ptau2.z()/Etau2;

    A(38,38) = -((pow(Enu2,2) - pow(pnu2.x(),2))/pow(Enu2,3))*l(13);
    A(38,L+10) = -1;
    A(38,L+13) = pnu2.x()/Enu2;

    A(39,39) = -((pow(Enu2,2) - pow(pnu2.y(),2))/pow(Enu2,3))*l(13);
    A(39,L+11) = -1;
    A(39,L+13) = -pnu2.y()/Enu2;

    A(40,40) = -((pow(Enu2,2) - pow(pnu2.z(),2))/pow(Enu2,3))*l(13);
    A(40,L+12) = -1;
    A(40,L+13) = -pnu2.z()/Enu2;
  
    A.Print();

    return U_cov;
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

double equations( const double* x_vars )
{
    // 69 equations to be solved
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

    // Unknown parameters (19)
    ROOT::Math::XYZPoint BV( x(dimM), x(dimM+1), x(dimM+2) );
    ROOT::Math::XYZVector pB( x(dimM+3), x(dimM+4), x(dimM+5) );
    Double_t MB_squared = x(dimM+6);
    ROOT::Math::XYZVector ptau1( x(dimM+7), x(dimM+8), x(dimM+9) );
    ROOT::Math::XYZVector pnu1( x(dimM+10), x(dimM+11), x(dimM+12) );
    ROOT::Math::XYZVector ptau2( x(dimM+13), x(dimM+14), x(dimM+15) );
    ROOT::Math::XYZVector pnu2( x(dimM+16), x(dimM+17), x(dimM+18) );

    Double_t EB = sqrt( MB_squared + pB.Mag2() );
    // tau and neutrino mass constraints are applied by simple parameter substitution
    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );

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
    eqs(0) = -2*chi2_sum(0) + ( pB.x()/pow(BV.x()-PV.x(),2) )*(l(0) + l(1));
    eqs(1) = -2*chi2_sum(1) - ( pB.y()/pow(BV.y()-PV.y(),2) )*l(0);
    eqs(2) = -2*chi2_sum(2) - ( pB.z()/pow(BV.z()-PV.z(),2) )*l(1);
    // DV1
    eqs(3) = -2*chi2_sum(3) - ( ptau1.x()/pow(DV1.x()-BV.x(),2) )*( l(2) + l(3) );
    eqs(4) = -2*chi2_sum(4) + ( ptau1.y()/pow(DV1.y()-BV.y(),2) )*l(2);
    eqs(5) = -2*chi2_sum(5) + ( ptau1.z()/pow(DV1.z()-BV.z(),2) )*l(3);
    // P3pi1
    eqs(6) = -2*chi2_sum(6) - l(4);
    eqs(7) = -2*chi2_sum(7) - l(5);
    eqs(8) = -2*chi2_sum(8) - l(6);
    eqs(9) = -2*chi2_sum(9) - l(7);
    // DV2
    eqs(10) = -2*chi2_sum(10) - ( ptau2.x()/pow(DV2.x()-BV.x(),2) )*(l(8) + l(9));
    eqs(11) = -2*chi2_sum(11) + ( ptau2.y()/pow(DV2.y()-BV.y(),2) )*l(8);
    eqs(12) = -2*chi2_sum(12) + ( ptau2.z()/pow(DV2.z()-BV.z(),2) )*l(9);
    // P3pi2
    eqs(13) = -2*chi2_sum(13) - l(10);
    eqs(14) = -2*chi2_sum(14) - l(11);
    eqs(15) = -2*chi2_sum(15) - l(12);
    eqs(16) = -2*chi2_sum(16) - l(13);
    // RP_T
    eqs(17) = -2*chi2_sum(17) + ( pK.x()/pow(BV.x()-RP.x(),2) )*(l(14) + l(15));
    eqs(18) = -2*chi2_sum(18) - ( pK.y()/pow(BV.y()-RP.y(),2) )*l(14);
    // pK
    eqs(19) = -2*chi2_sum(19) + (1/(BV.x() - RP.x()))*(l(14) + l(15)) - l(16) - (pK.x()/EK)*l(19);
    eqs(20) = -2*chi2_sum(20) - (1/(BV.y() - RP.y()))*l(14) - l(17) - (pK.y()/EK)*l(19);
    eqs(21) = -2*chi2_sum(21) - (1/(BV.z() - RP.z()))*l(15) - l(18) - (pK.z()/EK)*l(19);
    // BV
    eqs(dimM) = -( pB.x()/pow(BV.x()-PV.x(),2) )*(l(0) + l(1)) + ( ptau1.x()/pow(DV1.x()-BV.x(),2) )*(l(2) + l(3)) + ( ptau2.x()/pow(DV2.x()-BV.x(),2) )*(l(8) + l(9)) - ( pK.x()/pow(BV.x()-RP.x(),2) )*(l(14) + l(15));
    eqs(dimM+1) = ( pB.y()/pow(BV.y()-PV.y(),2) )*l(0) - ( ptau1.y()/pow(DV1.y()-BV.y(),2) )*l(2) - ( ptau2.y()/pow(DV2.y()-BV.y(),2) )*l(8) + ( pK.y()/pow(BV.y()-RP.y(),2) )*l(14);
    eqs(dimM+2) = ( pB.z()/pow(BV.z()-PV.z(),2) )*l(1) - ( ptau1.z()/pow(DV1.z()-BV.z(),2) )*l(3) - ( ptau2.z()/pow(DV2.z()-BV.z(),2) )*l(9) + ( pK.z()/pow(BV.z()-RP.z(),2) )*l(15);
    // pB, mB^2
    eqs(dimM+3) = (1/(BV.x() - PV.x()))*(l(0) + l(1)) + l(16);
    eqs(dimM+4) = -(1/(BV.y() - PV.y()))*l(0) + l(17);
    eqs(dimM+5) = -(1/(BV.z() - PV.z()))*l(1) + l(18);
    eqs(dimM+6) = (1/(2*EB))*l(19);
    // ptau1
    eqs(dimM+7) = (1/(DV1.x() - BV.x()))*(l(2) + l(3)) + l(4) + (ptau1.x()/Etau1)*l(7) - l(16) - (ptau1.x()/Etau1)*l(19);
    eqs(dimM+8) = -(1/(DV1.y() - BV.y()))*l(2) + l(5) + (ptau1.y()/Etau1)*l(7) - l(17) - (ptau1.y()/Etau1)*l(19);
    eqs(dimM+9) = -(1/(DV1.z() - BV.z()))*l(3) + l(6) + (ptau1.z()/Etau1)*l(7) - l(18) - (ptau1.z()/Etau1)*l(19);
    // eqs(dimM+10) = l(7) + l(8) - l(23);
    // pnu1
    eqs(dimM+10) = -l(4) - (pnu1.x()/Enu1)*l(7);
    eqs(dimM+11) = -l(5) - (pnu1.y()/Enu1)*l(7);
    eqs(dimM+12) = -l(6) - (pnu1.z()/Enu1)*l(7);
    // eqs(dimM+14) = -l(7) + l(9);
    // ptau2
    eqs(dimM+13) = (1/(DV2.x() - BV.x()))*(l(8) + l(9)) + l(10) + (ptau2.x()/Etau2)*l(13) - l(16) - (ptau2.x()/Etau2)*l(19);
    eqs(dimM+14) = -(1/(DV2.y() - BV.y()))*l(8) + l(11) + (ptau2.y()/Etau2)*l(13) - l(17) - (ptau2.y()/Etau2)*l(19);
    eqs(dimM+15) = -(1/(DV2.z() - BV.z()))*l(9) + l(12) + (ptau2.z()/Etau2)*l(13) - l(18) - (ptau2.z()/Etau2)*l(19);
    // eqs(dimM+18) = l(15) + l(16) - l(23);
    // pnu2
    eqs(dimM+16) = -l(10) - (pnu2.x()/Enu2)*l(13);
    eqs(dimM+17) = -l(11) - (pnu2.y()/Enu2)*l(13);
    eqs(dimM+18) = -l(12) - (pnu2.z()/Enu2)*l(13);
    // eqs(dimM+22) = -l(15) + l(17);
    // pB must point back to the PV (2)
    eqs(dimM+dimX) = pB.x()/( BV.x() - PV.x() ) - pB.y()/( BV.y() - PV.y() );
    eqs(dimM+dimX+1) = pB.x()/( BV.x() - PV.x() ) - pB.z()/( BV.z() - PV.z() );
    // ptau1 must point back to the BV (2)
    eqs(dimM+dimX+2) = ptau1.x()/( DV1.x() - BV.x() ) - ptau1.y()/( DV1.y() - BV.y() );
    eqs(dimM+dimX+3) = ptau1.x()/( DV1.x() - BV.x() ) - ptau1.z()/( DV1.z() - BV.z() );
    // 4-momentum conservation in DV1 (4)
    eqs(dimM+dimX+4) = ptau1.x() - p3pi1.x() - pnu1.x();
    eqs(dimM+dimX+5) = ptau1.y() - p3pi1.y() - pnu1.y();
    eqs(dimM+dimX+6) = ptau1.z() - p3pi1.z() - pnu1.z();
    eqs(dimM+dimX+7) = Etau1 - E3pi1 - Enu1;
    // tau+ and anti-nu mass constraints (2)
    // eqs(dimM+dimX+8) = Etau1 - sqrt( pow(mtau,2) + ptau1.Mag2() );
    // eqs(dimM+dimX+9) = Enu1 - sqrt( pnu1.Mag2() );
    // ptau2 must point back to the BV (2)
    eqs(dimM+dimX+8) = ptau2.x()/( DV2.x() - BV.x() ) - ptau2.y()/( DV2.y() - BV.y() );
    eqs(dimM+dimX+9) = ptau2.x()/( DV2.x() - BV.x() ) - ptau2.z()/( DV2.z() - BV.z() );
    // 4-momentum conservation in DV2 (4)
    eqs(dimM+dimX+10) = ptau2.x() - p3pi2.x() - pnu2.x();
    eqs(dimM+dimX+11) = ptau2.y() - p3pi2.y() - pnu2.y();
    eqs(dimM+dimX+12) = ptau2.z() - p3pi2.z() - pnu2.z();
    eqs(dimM+dimX+13) = Etau2 - E3pi2 - Enu2;
    // tau- and nu mass constraints (2)
    // eqs(dimM+dimX+16) = Etau2 - sqrt( pow(mtau,2) + ptau2.Mag2() );
    // eqs(dimM+dimX+17) = Enu2 - sqrt( pnu2.Mag2() );
    // BV must lie in K+ trajectory (2)
    eqs(dimM+dimX+14) = pK.x()/( BV.x() - RP.x() ) - pK.y()/( BV.y() - RP.y() );
    eqs(dimM+dimX+15) = pK.x()/( BV.x() - RP.x() ) - pK.z()/( BV.z() - RP.z() );
    // 4-momentum conservation in BV (4)
    eqs(dimM+dimX+16) = pB.x() - ptau1.x() - ptau2.x() - pK.x();
    eqs(dimM+dimX+17) = pB.y() - ptau1.y() - ptau2.y() - pK.y();
    eqs(dimM+dimX+18) = pB.z() - ptau1.z() - ptau2.z() - pK.z();
    eqs(dimM+dimX+19) = EB - Etau1 - Etau2 - EK;
 
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

Double_t lagrangian( const Double_t* x_values )
{
    TVectorD x(dimM+dimX+dimC);
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        x(i) = x_values[i];
    }

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
    Double_t EB = sqrt( MB_squared + pB.Mag2() );
    ROOT::Math::XYZVector ptau1( x(dimM+7), x(dimM+8), x(dimM+9) );
    Double_t Etau1 = x(dimM+10);
    ROOT::Math::XYZVector pnu1( x(dimM+11), x(dimM+12), x(dimM+13) );
    Double_t Enu1 = x(dimM+14);
    ROOT::Math::XYZVector ptau2( x(dimM+15), x(dimM+16), x(dimM+17) );
    Double_t Etau2 = x(dimM+18);
    ROOT::Math::XYZVector pnu2( x(dimM+19), x(dimM+20), x(dimM+21) );
    Double_t Enu2 = x(dimM+22);

    TVectorD g(dimC);
    // pB must point back to the PV (2)
    g(0) = pB.x()*( BV.x() - PV.x() ) - pB.y()*( BV.y() - PV.y() );
    g(1) = pB.x()*( BV.x() - PV.x() ) - pB.z()*( BV.z() - PV.z() );
    // ptau1 must point back to the BV (2)
    g(2) = ptau1.x()*( DV1.x() - BV.x() ) - ptau1.y()*( DV1.y() - BV.y() );
    g(3) = ptau1.x()*( DV1.x() - BV.x() ) - ptau1.z()*( DV1.z() - BV.z() );
    // 4-momentum conservation in DV1 (4)
    g(4) = ptau1.x() - p3pi1.x() - pnu1.x();
    g(5) = ptau1.y() - p3pi1.y() - pnu1.y();
    g(6) = ptau1.z() - p3pi1.z() - pnu1.z();
    g(7) = Etau1 - E3pi1 - Enu1;
    // tau+ and anti-nu mass constraints (2)
    g(8) = Etau1 - sqrt( pow(mtau,2) + ptau1.Mag2() );
    g(9) = Enu1 - sqrt( pnu1.Mag2() );
    // ptau2 must point back to the BV (2)
    g(10) = ptau2.x()*( DV2.x() - BV.x() ) - ptau2.y()*( DV2.y() - BV.y() );
    g(11) = ptau2.x()*( DV2.x() - BV.x() ) - ptau2.z()*( DV2.z() - BV.z() );
    // 4-momentum conservation in DV2 (4)
    g(12) = ptau2.x() - p3pi2.x() - pnu2.x();
    g(13) = ptau2.y() - p3pi2.y() - pnu2.y();
    g(14) = ptau2.z() - p3pi2.z() - pnu2.z();
    g(15) = Etau2 - E3pi2 - Enu2;
    // tau- and nu mass constraints (2)
    g(16) = Etau2 - sqrt( pow(mtau,2) + ptau2.Mag2() );
    g(17) = Enu2 - sqrt( pnu2.Mag2() );
    // BV must lie in K+ trajectory (2)
    g(18) = pK.x()*( BV.x() - RP.x() ) - pK.y()*( BV.y() - RP.y() );
    g(19) = pK.x()*( BV.x() - RP.x() ) - pK.z()*( BV.z() - RP.z() );
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
    x0(dimM+10) = pnu1.x();
    x0(dimM+11) = pnu1.y();
    x0(dimM+12) = pnu1.z();
    x0(dimM+13) = ptau2.x();
    x0(dimM+14) = ptau2.y();
    x0(dimM+15) = ptau2.z();
    x0(dimM+16) = pnu2.x();
    x0(dimM+17) = pnu2.y();
    x0(dimM+18) = pnu2.z();

    // 3) Initialise lambda (to 1)
    // for(int i = 0; i < dimC; i++)
    // {
    //     x0(dimM+dimX+i) = 1.;
    // }

    x0(dimM+dimX) = 1/10000; // pBx/(BVx - PVx)
    x0(dimM+dimX+1) = 1/10000; // pBx/(BVx - PVx)
    x0(dimM+dimX+2) = 1/10000; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+3) = 1/10000; // ptau1x/(DV1x - BVx)
    x0(dimM+dimX+4) = 1/1000;
    x0(dimM+dimX+5) = 1/1000;
    x0(dimM+dimX+6) = 1/10000;
    x0(dimM+dimX+7) = 1/10000;
    x0(dimM+dimX+8) = 1/10000; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+9) = 1/10000; // ptau2x/(DV2x - BVx)
    x0(dimM+dimX+10) = 1/1000;
    x0(dimM+dimX+11) = 1/1000;
    x0(dimM+dimX+12) = 1/10000;
    x0(dimM+dimX+13) = 1/10000;
    x0(dimM+dimX+14) = 1/1000;
    x0(dimM+dimX+15) = 1/1000;
    x0(dimM+dimX+16) = 1/10000;
    x0(dimM+dimX+17) = 1/10000;
    x0(dimM+dimX+18) = 1/100000;
    x0(dimM+dimX+19) = 1/100000;

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
    x0(dimM+10) = pnu1.x();
    x0(dimM+11) = pnu1.y();
    x0(dimM+12) = pnu1.z();
    x0(dimM+13) = ptau2.x();
    x0(dimM+14) = ptau2.y();
    x0(dimM+15) = ptau2.z();
    x0(dimM+16) = pnu2.x();
    x0(dimM+17) = pnu2.y();
    x0(dimM+18) = pnu2.z();

    // 3) Initialise lambda (to 1)
    // for(int i = 0; i < dimC; i++)
    // {
    //     x0(dimM+dimX+i) = 1.;
    // }

    x0(dimM+dimX) = 1/10000; 
    x0(dimM+dimX+1) = 1/10000;
    x0(dimM+dimX+2) = 1/1000;
    x0(dimM+dimX+3) = 1/1000;
    x0(dimM+dimX+4) = 1/1000;
    x0(dimM+dimX+5) = 1/1000;
    x0(dimM+dimX+6) = 1/10000;
    x0(dimM+dimX+7) = 1/10000;
    x0(dimM+dimX+8) = 1/1000;
    x0(dimM+dimX+9) = 1/1000;
    x0(dimM+dimX+10) = 1/1000;
    x0(dimM+dimX+11) = 1/1000;
    x0(dimM+dimX+12) = 1/10000;
    x0(dimM+dimX+13) = 1/10000;
    x0(dimM+dimX+14) = 1/1000;
    x0(dimM+dimX+15) = 1/1000;
    x0(dimM+dimX+16) = 1/10000;
    x0(dimM+dimX+17) = 1/10000;
    x0(dimM+dimX+18) = 1/100000;
    x0(dimM+dimX+19) = 1/100000;

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
    x0(dimM+10) = pnu1.x();
    x0(dimM+11) = pnu1.y();
    x0(dimM+12) = pnu1.z();
    x0(dimM+13) = ptau2.x();
    x0(dimM+14) = ptau2.y();
    x0(dimM+15) = ptau2.z();
    x0(dimM+16) = pnu2.x();
    x0(dimM+17) = pnu2.y();
    x0(dimM+18) = pnu2.z();

    // 3) Initialise lambda (to 1)
    // for(int i = 0; i < dimC; i++)
    // {
    //     x0(dimM+dimX+i) = 1.;
    // }

    x0(dimM+dimX) = 1/10000; 
    x0(dimM+dimX+1) = 1/10000;
    x0(dimM+dimX+2) = 1/1000;
    x0(dimM+dimX+3) = 1/1000;
    x0(dimM+dimX+4) = 1/1000;
    x0(dimM+dimX+5) = 1/1000;
    x0(dimM+dimX+6) = 1/10000;
    x0(dimM+dimX+7) = 1/10000;
    x0(dimM+dimX+8) = 1/1000;
    x0(dimM+dimX+9) = 1/1000;
    x0(dimM+dimX+10) = 1/1000;
    x0(dimM+dimX+11) = 1/1000;
    x0(dimM+dimX+12) = 1/10000;
    x0(dimM+dimX+13) = 1/10000;
    x0(dimM+dimX+14) = 1/1000;
    x0(dimM+dimX+15) = 1/1000;
    x0(dimM+dimX+16) = 1/10000;
    x0(dimM+dimX+17) = 1/10000;
    x0(dimM+dimX+18) = 1/100000;
    x0(dimM+dimX+19) = 1/100000;

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

    Double_t Enu1 = sqrt( pnu1.Mag2() );
    Double_t Enu2 = sqrt( pnu2.Mag2() );

    Double_t Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
    Double_t Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );

    ROOT::Math::XYZVector pB = pK + ptau1 + ptau2;
    Double_t EB = EK + Etau1 + Etau2;
    Double_t MB_squared = pow(EB,2) - pB.Mag2();

    x0(dimM+0) = BV.x();
    x0(dimM+1) = BV.y();
    x0(dimM+2) = BV.z();
    x0(dimM+3) = pB.x();
    x0(dimM+4) = pB.y();
    x0(dimM+5) = pB.z();
    x0(dimM+6) = MB_squared;
    x0(dimM+7) = ptau1.x();
    x0(dimM+8) = ptau1.y();
    x0(dimM+9) = ptau1.z();
    x0(dimM+10) = pnu1.x();
    x0(dimM+11) = pnu1.y();
    x0(dimM+12) = pnu1.z();
    x0(dimM+13) = ptau2.x();
    x0(dimM+14) = ptau2.y();
    x0(dimM+15) = ptau2.z();
    x0(dimM+16) = pnu2.x();
    x0(dimM+17) = pnu2.y();
    x0(dimM+18) = pnu2.z();

    // 3) Initialise lambda (to 1)
    // for(int i = 0; i < dimC; i++)
    // {
    //     x0(dimM+dimX+i) = 1.;
    // }

    x0(dimM+dimX) = 1/10000; 
    x0(dimM+dimX+1) = 1/10000;
    x0(dimM+dimX+2) = 1/1000;
    x0(dimM+dimX+3) = 1/1000;
    x0(dimM+dimX+4) = 1/1000;
    x0(dimM+dimX+5) = 1/1000;
    x0(dimM+dimX+6) = 1/10000;
    x0(dimM+dimX+7) = 1/10000;
    x0(dimM+dimX+8) = 1/1000;
    x0(dimM+dimX+9) = 1/1000;
    x0(dimM+dimX+10) = 1/1000;
    x0(dimM+dimX+11) = 1/1000;
    x0(dimM+dimX+12) = 1/10000;
    x0(dimM+dimX+13) = 1/10000;
    x0(dimM+dimX+14) = 1/1000;
    x0(dimM+dimX+15) = 1/1000;
    x0(dimM+dimX+16) = 1/10000;
    x0(dimM+dimX+17) = 1/10000;
    x0(dimM+dimX+18) = 1/100000;
    x0(dimM+dimX+19) = 1/100000;

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