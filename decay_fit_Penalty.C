
using namespace std;

// Global variables 
int dimM = 22; // number of measured parameters
int dimX = 23; // number of unknown parameters
int dimC = 24; // number of constraints

TVectorD m(dimM);
TMatrixDSym V(dimM);
TMatrixDSym W(dimM);

TVectorD x0(dimM+dimX);
TVectorD x0_current(dimM+dimX), x1_current(dimM+dimX), x2_current(dimM+dimX), x3_current(dimM+dimX), x4_current(dimM+dimX), x5_current(dimM+dimX), x6_current(dimM+dimX), x7_current(dimM+dimX), x8_current(dimM+dimX), x9_current(dimM+dimX), x_true_current(dimM+dimX);
Double_t RPz = 0.; // z-component of the RP on the K+ trajectory (fixed)

// Fit results
Double_t MB, MB_err, chi2, sum_of_constraints, deltaT;
TVectorD X(dimM+dimX);
TVectorD X_ERR(dimM+dimX);
Int_t status, init, cov_status, nIter;
TVectorD g(dimC);
TMatrixD C(dimM+dimX,dimM+dimX);

// Constants
Double_t mtau = 1776.86;
Double_t mkaon = 493.677;

// Functions
void minimize( ROOT::Math::XYZPoint BV, int init );
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

Double_t chisquare( TVectorD x );
TVectorD exact_constraints( TVectorD x );
Double_t function_to_minimize( const Double_t* x_values );
Double_t g_sum(TVectorD g);
void sequence(ROOT::Math::XYZPoint BV);
void lowest_chi2(ROOT::Math::XYZPoint BV);
Bool_t fisrtTrial = true;

void decay_fit_Penalty(int year, TString RECO_files, int species, int line)
{
    // Retrieve m and V from ntuple
    TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files, 1, line);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    Float_t m_vars[dimM];
    Float_t V_vars[dimM][dimM];
    Double_t BVx, BVy, BVz; // offline estimate for the BV is needed to have a first estimate for the unknown parameters using some of the initialisations

    t->SetBranchAddress("df_m", m_vars);
    t->SetBranchAddress("df_V", V_vars);
    t->SetBranchAddress("Kp_RP_Z", &RPz);
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

    TString name_g[] = {
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

    for(int i = 0; i < dimM+dimX; i++)
    {
        tout->Branch(name_x[i], &X(i));
        tout->Branch(name_x_err[i], &X_ERR(i));
    }

    for(int i = 0; i < dimC; i++)
    {
        tout->Branch(name_g[i], &g(i));
    }

    tout->Branch("df_init", &init);
    tout->Branch("df_Bp_M", &MB);
    tout->Branch("df_Bp_MERR", &MB_err);
    tout->Branch("df_status", &status);
    tout->Branch("df_cov_status", &cov_status);
    tout->Branch("df_chi2", &chi2);
    tout->Branch("df_g_sum", &sum_of_constraints);
    tout->Branch("df_nIter", &nIter);
    tout->Branch("df_time", &deltaT);

    for(int i = 0; i < dimC; i++)
    {
        tout->Branch(Form("df_F_%i",i), &g(i));
    }

    TStopwatch watch;

    // Loop over events
    UInt_t num_entries = t->GetEntries();
    for(int evt = 0; evt < 10; evt++)
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

        ROOT::Math::XYZPoint BV( BVx, BVy, BVz ); // offline estimate of BV (necessary for some initialisations)

        // for(int init = 0; init < 10; init++)
        // {
        //     minimize(BV, init);
        
        //     cout << "init = " << init << endl;
        //     cout << "chi2 = " << chi2 << endl;
        //     cout << "MB = " << MB << " +/- " << MB_err << endl;
        //     cout << "Status = " << status << endl;
        //     cout << "Cov matrix status = " << cov_status << endl;
        //     cout << "sum of constraints = " << sum_of_constraints << endl;
        //     cout << "Niter = " << nIter << endl;
        // }
        // sequence(BV);

        lowest_chi2(BV);
        // if( sum_of_constraints > 10 )
        // {
        //     fisrtTrial = false;
        //     lowest_chi2(BV);
        // }
    
        // cout << "init = " << init << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "MB = " << MB << " +/- " << MB_err << endl;
        cout << "Status = " << status << endl;
        // cout << "Cov matrix status = " << cov_status << endl;
        cout << "sum of constraints = " << sum_of_constraints << endl;
        cout << "Niter = " << nIter << endl;

        g.Print();

        cout << " ///////////////////////////////////////////////////////////////////////// " << endl;

        // X.Print();
        // X_ERR.Print();
        // TVectorD g_final = exact_constraints(X);
        // g_final.Print();
    
        // TVectorD diff(dimM);
        // for(int i = 0; i < dimM; i++)
        // {
        //     diff(i) = m(i) - X(i);
        // }
        // diff.Print();
        // g.Print();

        tout->Fill();
    }
    watch.Stop();
    deltaT = watch.RealTime();

    cout << "FINISHED" << endl;

    cout << "FINISHED in " << deltaT << " seconds." << endl;
    cout << "Looped over " << num_entries << " events." << endl;

    fout->cd();
    tout->Write();
    fout->Close();
}

void lowest_chi2(ROOT::Math::XYZPoint BV)
{
        // Lowest sum
        // Init 0
        minimize( BV, 0 );
        x0_current = X;
        Int_t status0 = status;
        Double_t MB_0 = MB;
        Double_t MB_err0 = MB_err;
        TVectorD X0 = X;
        TVectorD X0_ERR = X_ERR;
        Double_t chi2_0 = chi2;
        Int_t nIter_0 = nIter;
        Double_t sum0 = sum_of_constraints;
        TVectorD g0 = g;

        // cout << "MB_0 = " << MB_0 << endl;
        // cout << "sum_0 = " << F_tol_0 << endl;
        // cout << "chi2_0 = " << chi2_0 << endl;
        // cout << "status_0 = " << status0 << endl;
        // cout << "nIter_0 = " << nIter_0 << endl;

        // Init 1
        minimize( BV, 1 );
        x1_current = X;
        Int_t status1 = status;
        Double_t MB_1 = MB;
        Double_t MB_err1 = MB_err;
        TVectorD X1 = X;
        TVectorD X1_ERR = X_ERR;
        Double_t chi2_1 = chi2;
        Int_t nIter_1 = nIter;
        Double_t sum1 = sum_of_constraints;
        TVectorD g1 = g;
    
        // cout << "MB_1 = " << MB_1 << endl;
        // cout << "sum_1 = " << F_tol_1 << endl;
        // cout << "chi2_1 = " << chi2_1 << endl;
        // cout << "status_1 = " << status1 << endl;
        // cout << "nIter_1 = " << nIter_1 << endl;

        // Init 2
        minimize( BV, 2 );
        x2_current = X;
        Int_t status2 = status;
        Double_t MB_2 = MB;
        Double_t MB_err2 = MB_err;
        TVectorD X2 = X;
        TVectorD X2_ERR = X_ERR;
        Double_t chi2_2 = chi2;
        Int_t nIter_2 = nIter;
        Double_t sum2 = sum_of_constraints;
        TVectorD g2 = g;

        // cout << "MB_2 = " << MB_2 << endl;
        // cout << "sum_2 = " << F_tol_2 << endl;
        // cout << "chi2_2 = " << chi2_2 << endl;
        // cout << "status_2 = " << status2 << endl;
        // cout << "nIter_2 = " << nIter_2 << endl;

        // Init 3
        minimize( BV, 3 );
        x3_current = X;
        Int_t status3 = status;
        Double_t MB_3 = MB;
        Double_t MB_err3 = MB_err;
        TVectorD X3 = X;
        TVectorD X3_ERR = X_ERR;
        Double_t chi2_3 = chi2;
        Int_t nIter_3 = nIter;
        Double_t sum3 = sum_of_constraints;
        TVectorD g3 = g;

        // cout << "MB_3 = " << MB_3 << endl;
        // cout << "sum_3 = " << F_tol_3 << endl;
        // cout << "chi2_3 = " << chi2_3 << endl;
        // cout << "status_3 = " << status3 << endl;
        // cout << "nIter_3 = " << nIter_3 << endl;

        // Init 4
        minimize( BV, 4 );
        x4_current = X;
        Int_t status4 = status;
        Double_t MB_4 = MB;
        Double_t MB_err4 = MB_err;
        TVectorD X4 = X;
        TVectorD X4_ERR = X_ERR;
        Double_t chi2_4 = chi2;
        Int_t nIter_4 = nIter;
        Double_t sum4 = sum_of_constraints;
        TVectorD g4 = g;

        // cout << "MB_4 = " << MB_4 << endl;
        // cout << "sum_4 = " << F_tol_4 << endl;
        // cout << "chi2_4 = " << chi2_4 << endl;
        // cout << "status_4 = " << status4 << endl;
        // cout << "nIter_4 = " << nIter_4 << endl;

        // Init 5
        minimize( BV, 5 );
        x5_current = X;
        Int_t status5 = status;
        Double_t MB_5 = MB;
        Double_t MB_err5 = MB_err;
        TVectorD X5 = X;
        TVectorD X5_ERR = X_ERR;
        Double_t chi2_5 = chi2;
        Int_t nIter_5 = nIter;
        Double_t sum5 = sum_of_constraints;
        TVectorD g5 = g;

        // cout << "MB_5 = " << MB_5 << endl;
        // cout << "sum_5 = " << F_tol_5 << endl;
        // cout << "chi2_5 = " << chi2_5 << endl;
        // cout << "status_5 = " << status5 << endl;
        // cout << "nIter_5 = " << nIter_5 << endl;

        // Init 6
        minimize( BV, 6 );
        x6_current = X;
        Int_t status6 = status;
        Double_t MB_6 = MB;
        Double_t MB_err6 = MB_err;
        TVectorD X6 = X;
        TVectorD X6_ERR = X_ERR;
        Double_t chi2_6 = chi2;
        Int_t nIter_6 = nIter;
        Double_t sum6 = sum_of_constraints;
        TVectorD g6 = g;

        // cout << "MB_6 = " << MB_6 << endl;
        // cout << "sum_6 = " << F_tol_6 << endl;
        // cout << "chi2_6 = " << chi2_6 << endl;
        // cout << "status_6 = " << status6 << endl;
        // cout << "nIter_6 = " << nIter_6 << endl;

        // Init 7
        minimize( BV, 7 );
        x7_current = X;
        Int_t status7 = status;
        Double_t MB_7 = MB;
        Double_t MB_err7 = MB_err;
        TVectorD X7 = X;
        TVectorD X7_ERR = X_ERR;
        Double_t chi2_7 = chi2;
        Int_t nIter_7 = nIter;
        Double_t sum7 = sum_of_constraints;
        TVectorD g7 = g;

        // cout << "MB_7 = " << MB_7 << endl;
        // cout << "sum_7 = " << F_tol_7 << endl;
        // cout << "chi2_7 = " << chi2_7 << endl;
        // cout << "status_7 = " << status7 << endl;
        // cout << "nIter_7 = " << nIter_7 << endl;

        // Init 8
        minimize( BV, 8 );
        x8_current = X;
        Int_t status8 = status;
        Double_t MB_8 = MB;
        Double_t MB_err8 = MB_err;
        TVectorD X8 = X;
        TVectorD X8_ERR = X_ERR;
        Double_t chi2_8 = chi2;
        Int_t nIter_8 = nIter;
        Double_t sum8 = sum_of_constraints;
        TVectorD g8 = g;

        // cout << "MB_8 = " << MB_8 << endl;
        // cout << "sum_8 = " << F_tol_8 << endl;
        // cout << "chi2_8 = " << chi2_8 << endl;
        // cout << "status_8 = " << status8 << endl;
        // cout << "nIter_8 = " << nIter_8 << endl;

        // Init 9
        minimize( BV, 9 );
        x9_current = X;
        Int_t status9 = status;
        Double_t MB_9 = MB;
        Double_t MB_err9 = MB_err;
        TVectorD X9 = X;
        TVectorD X9_ERR = X_ERR;
        Double_t chi2_9 = chi2;
        Int_t nIter_9 = nIter;
        Double_t sum9 = sum_of_constraints;
        TVectorD g9 = g;

        // cout << "MB_9 = " << MB_9 << endl;
        // cout << "sum_9 = " << F_tol_9 << endl;
        // cout << "chi2_9 = " << chi2_9 << endl;
        // cout << "status_9 = " << status9 << endl;
        // cout << "nIter_9 = " << nIter_9 << endl;

        Int_t status_vec[10] = { status0, status1, status2, status3, status4, status5, status6, status7, status8, status9 };
        Double_t MB_vec[10] = { MB_0, MB_1, MB_2, MB_3, MB_4, MB_5, MB_6, MB_7, MB_8, MB_9 }; 
        Double_t MB_err_vec[10] = { MB_err0, MB_err1, MB_err2, MB_err3, MB_err4, MB_err5, MB_err6, MB_err7, MB_err8, MB_err9 };
        Double_t chi2_vec[10] = { chi2_0, chi2_1, chi2_2, chi2_3, chi2_4, chi2_5, chi2_6, chi2_7, chi2_8, chi2_9 };
        Int_t nIter_vec[10] = { nIter_0, nIter_1, nIter_2, nIter_3, nIter_4, nIter_5, nIter_6, nIter_7, nIter_8, nIter_9 };
        Double_t sum[10] = { sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9 };

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

        std::vector<TVectorD> g_vec;
        g_vec.push_back(g0);
        g_vec.push_back(g1);
        g_vec.push_back(g2);
        g_vec.push_back(g3);
        g_vec.push_back(g4);
        g_vec.push_back(g5);
        g_vec.push_back(g6);
        g_vec.push_back(g7);
        g_vec.push_back(g8);
        g_vec.push_back(g9);

        bool all_fail = false;
        if( (status0 != 0) && (status1 != 0) && (status2 != 0) && (status3 != 0) && (status4 != 0) && (status5 != 0) && (status6 != 0) && (status7 != 0) && (status8 != 0) && (status9 != 0) )
        {
            all_fail = true;
        }

        Double_t chi2_min = 100000000000000000;
        Double_t sum_min = 100;
        Int_t i_min = 0;

        if(all_fail) // If all fail, returns the one with the lowest value for the sum
        {
            for(int i = 0; i < 10; i++)
            {
                if( (sum[i] < sum_min) )
                {
                    sum_min = sum[i];
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
        g = g_vec[i_min];
        sum_of_constraints = sum[i_min];
        nIter = nIter_vec[i_min];
        chi2 = chi2_vec[i_min];

}

void sequence(ROOT::Math::XYZPoint BV)
{
    init = 3; // Marseille
    minimize( BV, init );
    if(status != 0)
    {
        init = 0; // Original
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 9; // Marseille (tau-) + K*tautau pions (tau+)
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 8; // Marseille (tau+) + K*tautau pions (tau-)
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 7; // Marseille (tau-) + K*tautau vertex (tau+)
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 6; // Marseille (tau+) + K*tautau vertex (tau-)
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 2; // K*tautau pions
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 5; // K*tautau vertex (tau-) + K*tautau pions (tau+)
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 4; // K*tautau vertex (tau+) + K*tautau pions (tau-)
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = 1; // K*tautau vertex
        minimize( BV, init );
    }
    if(status != 0)
    {
        init = -1;
    }
}

void minimize( ROOT::Math::XYZPoint BV, int init )
{
    // 1) This function initialises the vector of unkown parameters x with the result of analytical calculations
    // 2) It sets up the minimizer: defining the initial values and bounds on the x parameters  
    // 3) It builds the chi^2 function with x0 by calling the function chisquare 
    // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

    // 1) Initial values for x=(xm,xu,lambda)
    if(init == 0)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_0( m ); // original calculations = Anne Keune
        }
        else
        {
            x0 = x0_current;
        }
    }
    else if(init == 1)
    {   
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_1( m, BV ); // K*tautau calculations; tau momentum direction initialised based on vertices (uses offline estimate for BV)
        }
        else
        {
            x0 = x1_current;
        }
    }
    else if(init == 2)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_2( m, BV ); // K*tautau calculations; tau momentum direction initialised based on 3pi visible momentum (uses offline estimate for BV) 
        }
        else
        {
            x0 = x2_current;
        }
    }
    else if(init == 3)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_3( m, BV ); // RD calculations (uses offline estimate for BV)
        }
        else
        {
            x0 = x3_current;
        }
    }
    else if(init == 4)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_4( m, BV ); // K*tautau; tau+ from vertices, tau- from pions
        }
        else
        {
            x0 = x4_current;
        }
    }
    else if(init == 5)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_5( m, BV ); // K*tautau; tau- from vertices, tau+ from pions
        }
        else
        {
            x0 = x5_current;
        }
    }
    else if(init == 6)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_6( m, BV ); // Mixes Marseille (tau+) and K*tautau vertices (tau-)            
        }
        else
        {
            x0 = x6_current;
        }
    }
    else if(init == 7)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_7( m, BV ); // Mixes Marseille (tau-) and K*tautau vertices (tau+)
        }
        else
        {
            x0 = x7_current;
        }
    }
    else if(init == 8)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_8( m, BV ); // Mixes Marseille (tau+) and K*tautau pions (tau-)
        }
        else
        {
            x0 = x8_current;
        }
    }
    else if(init == 9)
    {
        if(fisrtTrial)
        {
            x0 = x_initial_estimate_9( m, BV ); // Mixes Marseille (tau-) and K*tautau pions (tau+)
        }
        else
        {
            x0 = x9_current;
        }
    }
    // x0.Print();

    // cout << "initial chi2 = " << chisquare(x0) << endl;
    // cout << "initial g = " << endl;
    // TVectorD g0 = exact_constraints(x0);
    // g0.Print();

    Double_t x0_vars[dimM+dimX], x0_err[dimM+dimX];
    for(int i = 0; i < dimM+dimX; i++)
    {
        x0_vars[i] = x0(i); // x0_vars is a Double_t, x0 is a TVectorD; the function to minimise must receive a Double_t as input
        x0_err[i] = 0.1*abs(x0_vars[i]); // initial errors on x0; they are used as the first step size in the minimisation; considering a 10% error for now
    }

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");    
    ROOT::Math::Functor func(&function_to_minimize, dimM+dimX);

    min->SetMaxIterations(10000000);  
    min->SetMaxFunctionCalls(10000000);
    min->SetPrintLevel(0);
    min->SetFunction(func);

    string name[] = {
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
        "pKz",
        "BVx",
        "BVy",
        "BVz",
        "pBx",
        "pBy",
        "pBz",
        "MB_squared",
        "ptau1x",
        "ptau1y",
        "ptau1z",
        "Etau1",
        "pnu1x",
        "pnu1y",
        "pnu1z",
        "Enu1",
        "ptau2x",
        "ptau2y",
        "ptau2z",
        "Etau2",
        "pnu2x",
        "pnu2y",
        "pnu2z",
        "Enu2"
    };
    if( sizeof(name)/sizeof(name[0]) != dimM+dimX ){ 
        cout << "ERROR: name does not have the proper dimension" << endl;
        return; 
    }

    for(int i = 0; i < dimM+dimX; i++)
    {
        min->SetVariable(i, name[i], x0[i], x0_err[i]);
    }

    min->Minimize();
    // min->Hesse();

    // Return / save results from the fit
    chi2 = min->MinValue();
    status = min->Status();
    cov_status = min->CovMatrixStatus();
    nIter = min->NIterations();

    const Double_t *xMin = min->X();
    for(int i = 0; i < dimM+dimX; i++)
    {
        X(i) = xMin[i];
        for(int j = 0; j < dimM+dimX; j++)
        {
            C(i,j) = min->CovMatrix(i,j);
        }
        X_ERR(i) = sqrt(C(i,i));
    }

    Double_t MB_squared = X(dimM+6);
    if( MB_squared > 0 )
    {
        MB = sqrt( MB_squared );
    }
    else
    {
        MB = -sqrt( abs(MB_squared) );
    }

    Double_t dMB_squared = sqrt(C(dimM+6,dimM+6));
    MB_err = dMB_squared/(2*abs(MB));

    g = exact_constraints(X);
    sum_of_constraints = g_sum(g);

}

Double_t function_to_minimize( const Double_t* x_values )
{
    TVectorD x(dimM+dimX);
    for(int i = 0; i < dimM+dimX; i++)
    {
        x(i) = x_values[i];
    }

    // Measured parameters
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

    // Unknown parameters
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

    Double_t f = chisquare(x);
    TVectorD g = exact_constraints(x);

    TVectorD delta_xu(dimX);
    delta_xu(0) = 0.1;
    delta_xu(1) = 0.1;
    delta_xu(2) = 10;
    delta_xu(3) = 10;
    delta_xu(4) = 10;
    delta_xu(5) = 100;
    delta_xu(6) = 100;
    delta_xu(7) = 1;
    delta_xu(8) = 1;
    delta_xu(9) = 100;
    delta_xu(10) = 100;
    delta_xu(11) = 1;
    delta_xu(12) = 1;
    delta_xu(13) = 100;
    delta_xu(14) = 100;
    delta_xu(15) = 1;
    delta_xu(16) = 1;
    delta_xu(17) = 100;
    delta_xu(18) = 100;
    delta_xu(19) = 1;
    delta_xu(20) = 1;
    delta_xu(21) = 100;
    delta_xu(22) = 100;

    TVectorD delta_g(dimC); // errors squared on the constraints
    delta_g(0) = pB.z()*V(0,0)*pB.z() - pB.z()*V(0,2)*pB.x() - pB.x()*V(2,0)*pB.z() + pB.z()*pow(delta_xu(0),2)*pB.z() + pB.x()*pow(delta_xu(2),2)*pB.x() + (BV.z() - PV.z())*pow(delta_xu(3),2)*(BV.z()-PV.z()) + (BV.x()-PV.x())*pow(delta_xu(5),2)*(BV.x()-PV.x());
    delta_g(1) = pB.z()*V(1,1)*pB.z() - pB.z()*V(1,2)*pB.y() - pB.y()*V(2,1)*pB.z() + pB.z()*pow(delta_xu(1),2)*pB.z() + pB.y()*pow(delta_xu(2),2)*pB.y() + (BV.z() - PV.z())*pow(delta_xu(4),2)*(BV.z()-PV.z()) + (BV.y()-PV.y())*pow(delta_xu(5),2)*(BV.y()-PV.y());
    delta_g(2) = ptau1.z()*V(3,3)*ptau1.z() - ptau1.z()*V(3,5)*ptau1.x() - ptau1.x()*V(5,3)*ptau1.z() + ptau1.z()*pow(delta_xu(0),2)*ptau1.z() + ptau1.x()*pow(delta_xu(2),2)*ptau1.x() + (DV1.z()-BV.z())*pow(delta_xu(7),2)*(DV1.z()-BV.z()) + (DV1.x()-BV.x())*pow(delta_xu(9),2)*(DV1.x()-BV.x());
    delta_g(3) = ptau1.z()*V(4,4)*ptau1.z() - ptau1.z()*V(4,5)*ptau1.y() - ptau1.y()*V(5,4)*ptau1.z() + ptau1.z()*pow(delta_xu(1),2)*ptau1.z() + ptau1.y()*pow(delta_xu(2),2)*ptau1.y() + (DV1.z()-BV.z())*pow(delta_xu(8),2)*(DV1.z()-BV.z()) + (DV1.y()-BV.y())*pow(delta_xu(9),2)*(DV1.y()-BV.y());
    delta_g(4) = V(6,6) + pow(delta_xu(7),2) + pow(delta_xu(11),2);
    delta_g(5) = V(7,7) + pow(delta_xu(8),2) + pow(delta_xu(12),2);
    delta_g(6) = V(8,8) + pow(delta_xu(9),2) + pow(delta_xu(13),2);
    delta_g(7) = V(9,9) + pow(delta_xu(10),2) + pow(delta_xu(14),2);
    delta_g(8) = (ptau1.x()/sqrt( pow(mtau,2) + ptau1.Mag2()  ))*pow(delta_xu(7),2)*(ptau1.x()/sqrt( pow(mtau,2) + ptau1.Mag2()  )) + (ptau1.y()/sqrt( pow(mtau,2) + ptau1.Mag2()  ))*pow(delta_xu(8),2)*(ptau1.y()/sqrt( pow(mtau,2) + ptau1.Mag2()  )) + (ptau1.z()/sqrt( pow(mtau,2) + ptau1.Mag2()  ))*pow(delta_xu(9),2)*(ptau1.z()/sqrt( pow(mtau,2) + ptau1.Mag2()  )) + pow(delta_xu(10),2);
    delta_g(9) = (pnu1.x()/sqrt( pnu1.Mag2() ))*pow(delta_xu(11),2)*(pnu1.x()/sqrt( pnu1.Mag2() )) + (pnu1.y()/sqrt( pnu1.Mag2() ))*pow(delta_xu(12),2)*(pnu1.y()/sqrt( pnu1.Mag2() )) + (pnu1.z()/sqrt( pnu1.Mag2() ))*pow(delta_xu(13),2)*(pnu1.z()/sqrt( pnu1.Mag2() )) + pow(delta_xu(14),2);
    delta_g(10) = ptau2.z()*V(15,15)*ptau2.z() - ptau2.z()*V(15,17)*ptau2.x() - ptau2.x()*V(17,15)*ptau2.z() + ptau2.z()*pow(delta_xu(0),2)*ptau2.z() + ptau2.x()*pow(delta_xu(2),2)*ptau2.z() + (DV2.z()-BV.z())*pow(delta_xu(15),2)*(DV2.z()-BV.z()) + (DV2.x()-BV.x())*pow(delta_xu(17),2)*(DV2.x()-BV.x());
    delta_g(11) = ptau2.z()*V(16,16)*ptau2.z() - ptau2.z()*V(16,17)*ptau2.y() - ptau2.y()*V(17,16)*ptau2.z() + ptau2.z()*pow(delta_xu(1),2)*ptau2.z() + ptau2.y()*pow(delta_xu(2),2)*ptau2.y() + (DV2.z()-BV.z())*pow(delta_xu(16),2)*(DV2.z()-BV.z()) + (DV2.y()-BV.y())*pow(delta_xu(17),2)*(DV2.y()-BV.y());
    delta_g(12) = V(13,13) + pow(delta_xu(15),2) + pow(delta_xu(19),2);
    delta_g(13) = V(14,14) + pow(delta_xu(16),2) + pow(delta_xu(20),2);
    delta_g(14) = V(15,15) + pow(delta_xu(17),2) + pow(delta_xu(21),2);
    delta_g(15) = V(16,16) + pow(delta_xu(18),2) + pow(delta_xu(22),2);
    delta_g(16) = (ptau2.x()/sqrt( pow(mtau,2) + ptau2.Mag2()  ))*pow(delta_xu(15),2)*(ptau2.x()/sqrt( pow(mtau,2) + ptau2.Mag2()  )) + (ptau2.y()/sqrt( pow(mtau,2) + ptau2.Mag2()  ))*pow(delta_xu(16),2)*(ptau2.y()/sqrt( pow(mtau,2) + ptau2.Mag2()  )) + (ptau2.z()/sqrt( pow(mtau,2) + ptau2.Mag2()  ))*pow(delta_xu(17),2)*(ptau2.z()/sqrt( pow(mtau,2) + ptau2.Mag2()  )) + pow(delta_xu(18),2);
    delta_g(17) = (pnu2.x()/sqrt( pnu2.Mag2() ))*pow(delta_xu(19),2)*(pnu2.x()/sqrt( pnu2.Mag2() )) + (pnu2.y()/sqrt( pnu2.Mag2() ))*pow(delta_xu(20),2)*(pnu2.y()/sqrt( pnu2.Mag2() )) + (pnu2.z()/sqrt( pnu2.Mag2() ))*pow(delta_xu(21),2)*(pnu2.z()/sqrt( pnu2.Mag2() )) + pow(delta_xu(22),2);
    delta_g(18) = pK.z()*V(17,17)*pK.z() + pK.z()*V(17,19)*(BV.z()-RP.z()) + (BV.z()-RP.z())*V(19,17)*pK.z() - pK.z()*V(17,21)*(BV.x()-RP.x()) - (BV.x()-RP.x())*V(21,17)*pK.z() + (BV.z()-RP.z())*V(19,19)*(BV.z()-RP.z()) + (BV.x()-RP.x())*V(21,21)*(BV.x()-RP.x()) - (BV.z()-RP.z())*V(19,21)*(BV.x()-RP.x()) - (BV.x()-RP.x())*V(21,19)*(BV.z()-RP.z()) + pK.z()*pow(delta_xu(0),2)*pK.z() + pK.x()*pow(delta_xu(2),2)*pK.x();
    delta_g(19) = pK.z()*V(18,18)*pK.z() + pK.z()*V(18,20)*(BV.z()-RP.z()) + (BV.z()-RP.z())*V(20,18)*pK.z() - pK.z()*V(18,21)*(BV.y()-RP.y()) - (BV.y()-RP.y())*V(21,18)*pK.z() + (BV.z()-RP.z())*V(20,20)*(BV.z()-RP.z()) + (BV.y()-RP.y())*V(21,21)*(BV.y()-RP.y()) - (BV.z()-RP.z())*V(20,21)*(BV.y()-RP.y()) - (BV.y()-RP.y())*V(21,20)*(BV.z()-RP.z()) + pK.z()*pow(delta_xu(1),2)*pK.z() + pK.y()*pow(delta_xu(2),2)*pK.y();
    delta_g(20) = pow(delta_xu(3),2) + pow(delta_xu(7),2) + pow(delta_xu(15),2) + pow(delta_xu(19),2);
    delta_g(21) = pow(delta_xu(4),2) + pow(delta_xu(8),2) + pow(delta_xu(16),2) + pow(delta_xu(20),2);
    delta_g(22) = pow(delta_xu(5),2) + pow(delta_xu(9),2) + pow(delta_xu(17),2) + pow(delta_xu(21),2);
    delta_g(23) = (pB.x()/EB)*pow(delta_xu(3),2)*(pB.x()/EB) + (pB.y()/EB)*pow(delta_xu(4),2)*(pB.y()/EB) + (pB.z()/EB)*pow(delta_xu(5),2)*(pB.z()/EB) + (1./(2*EB))*pow(delta_xu(6),2)*(1./(2*EB)) + pow(delta_xu(10),2) + pow(delta_xu(18),2) + (pK.x()/EK)*V(19,19)*(pK.x()/EK) + (pK.x()/EK)*V(19,20)*(pK.y()/EK) + (pK.y()/EK)*V(20,19)*(pK.x()/EK) + (pK.x()/EK)*V(19,21)*(pK.z()/EK) + (pK.z()/EK)*V(21,19)*(pK.x()/EK) + (pK.y()/EK)*V(20,21)*(pK.z()/EK) + (pK.z()/EK)*V(21,20)*(pK.y()/EK) + (pK.y()/EK)*V(20,20)*(pK.y()/EK) + (pK.z()/EK)*V(21,21)*(pK.z()/EK);

    // delta_g.Print();
    // g.Print();

    Double_t P = 0;
    for(int i = 0; i < dimC; i++)
    {
        P += pow(g(i),2)/delta_g(i);
    }
    f += P;

    // cout << "chi2 = " << chisquare(x) << endl;
    // cout << "sum = " << P << endl;

    return f;
}

Double_t chisquare( TVectorD x )
{
    TVectorD xm(dimM);
    for(int i = 0; i < dimM; i++)
    {
        xm(i) = x(i);
    }

    // Computes the chi^2 = (m - xm)^T W (m - xm)
    TVectorD r = m - xm;

    Double_t chi2 = r*(W*r);

    return abs(chi2);
}

TVectorD exact_constraints( TVectorD x )
{
    // Measured parameters
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

    // Unknown parameters
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

    TVectorD g(dimC);
    // pB must point back to the PV (2)
    g(0) = pB.x()*(BV.z() - PV.z()) - pB.z()*(BV.x() - PV.x());
    g(1) = pB.y()*(BV.z() - PV.z()) - pB.z()*(BV.y() - PV.y());
    // ptau1 must point back to the BV (2)
    g(2) = ptau1.x()*(DV1.z() - BV.z()) - ptau1.z()*(DV1.x() - BV.x());
    g(3) = ptau1.y()*(DV1.z() - BV.z()) - ptau1.z()*(DV1.y() - BV.y());
    // 4-momentum conservation in DV1 (4)
    g(4) = ptau1.x() - p3pi1.x() - pnu1.x();
    g(5) = ptau1.y() - p3pi1.y() - pnu1.y();
    g(6) = ptau1.z() - p3pi1.z() - pnu1.z();
    g(7) = Etau1 - E3pi1 - Enu1;
    // tau+ and anti-nu mass constraints (2)
    g(8) = Etau1 - sqrt(pow(mtau,2) + ptau1.Mag2());
    g(9) = Enu1 - sqrt(pnu1.Mag2());
    // ptau2 must point back to the BV (2)
    g(10) = ptau2.x()*(DV2.z() - BV.z()) - ptau2.z()*(DV2.x() - BV.x());
    g(11) = ptau2.y()*(DV2.z() - BV.z()) - ptau2.z()*(DV2.y() - BV.y());
    // 4-momentum conservation in DV2 (4)
    g(12) = ptau2.x() - p3pi2.x() - pnu2.x();
    g(13) = ptau2.y() - p3pi2.y() - pnu2.y();
    g(14) = ptau2.z() - p3pi2.z() - pnu2.z();
    g(15) = Etau2 - E3pi2 - Enu2;
    // tau- and nu mass constraints (2)
    g(16) = Etau2 - sqrt(pow(mtau,2) + ptau2.Mag2());
    g(17) = Enu2 - sqrt(pnu2.Mag2());
    // BV must lie in K+ trajectory (2)
    g(18) = pK.x()*(BV.z() - RP.z()) - pK.z()*(BV.x() - RP.x());
    g(19) = pK.y()*(BV.z() - RP.z()) - pK.z()*(BV.y() - RP.y());
    // 4-momentum conservation in BV (4)
    g(20) = pB.x() - ptau1.x() - ptau2.x() - pK.x();
    g(21) = pB.y() - ptau1.y() - ptau2.y() - pK.y();
    g(22) = pB.z() - ptau1.z() - ptau2.z() - pK.z();
    g(23) = EB - Etau1 - Etau2 - EK;

    return g;
}

Double_t g_sum(TVectorD g)
{
    Double_t sum = 0;
    for(int i = 0; i < dimC; i++)
    {
        sum += abs(g(i));
    }
    return sum;
}

///////////////////////////////////////////////////////// INITIALISATIONS //////////////////////////////////////////////////////////////////////////////////////////

TVectorD x_initial_estimate_0( TVectorD m ) // Original initialisation for x (based on Anne Keune's thesis)
{
    // Builds an initial estimate for x based on the analytical calculations and on the known parameters in m
    TVectorD x0(dimM+dimX);

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

    return x0;
}

TVectorD x_initial_estimate_1( TVectorD m, ROOT::Math::XYZPoint BV ) 
{
    // Builds an initial estimate for x using B->K* tautau initialisation; taus direction based on vertices
    TVectorD x0(dimM+dimX);

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


    return x0;
}

TVectorD x_initial_estimate_2( TVectorD m, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; taus direction based visible 3pi momenta
{
    // Builds an initial estimate for x using B->K* tautau initialisation; taus direction based on visible 3pi momenta
    TVectorD x0(dimM+dimX);

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


    return x0;
}

TVectorD x_initial_estimate_3( TVectorD m, ROOT::Math::XYZPoint BV ) 
{
    // Builds an initial estimate for x based on the Marseille analytical calculations; it uses the offline estimate for BV  
    TVectorD x0(dimM+dimX);

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

    return x0;
}

TVectorD x_initial_estimate_4( TVectorD m, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; tau+ direction from vertices ; tau- direction from pions
{
    TVectorD x0(dimM+dimX);

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

    return x0;
}

TVectorD x_initial_estimate_5( TVectorD m, ROOT::Math::XYZPoint BV ) // Using B->K* tautau initialisation; tau- direction from vertices ; tau+ direction from pions
{
    TVectorD x0(dimM+dimX);

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


    return x0;
}

TVectorD x_initial_estimate_6( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau+) and K* tau tau vertices (tau-)
{  
    TVectorD x0(dimM+dimX);

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


    return x0;
}

TVectorD x_initial_estimate_7( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau-) and K* tau tau vertices (tau+)
{  
    TVectorD x0(dimM+dimX);

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

    return x0;
}

TVectorD x_initial_estimate_8( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau+) and K* tau tau pions (tau-)
{  
    TVectorD x0(dimM+dimX);

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

    return x0;
}

TVectorD x_initial_estimate_9( TVectorD m, ROOT::Math::XYZPoint BV )  // Mixes Marseille (tau-) and K* tau tau pions (tau+)
{  
    TVectorD x0(dimM+dimX);

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