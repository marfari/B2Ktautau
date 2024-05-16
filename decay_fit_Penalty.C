
using namespace std;

// Global variables 
int dimM = 22; // number of measured parameters
int dimX = 4; // number of unknown parameters
int dimC = 5; // number of constraints

TVectorD m(dimM);
TMatrixDSym V(dimM);
TMatrixDSym W(dimM);

TVectorD x0(dimM+dimX);
TVectorD x0_current(dimM+dimX), x1_current(dimM+dimX), x2_current(dimM+dimX), x3_current(dimM+dimX), x4_current(dimM+dimX), x5_current(dimM+dimX), x6_current(dimM+dimX), x7_current(dimM+dimX), x8_current(dimM+dimX), x9_current(dimM+dimX), x_true_current(dimM+dimX), x_current(dimM+dimX), x_initial_estimate_MLP_init(dimM+dimX);
Double_t status_current, cov_status_current;
Double_t RPz = 0.; // z-component of the RP on the K+ trajectory (fixed)

// Fit results
Double_t MB, MB_err, chi2, sum_of_constraints, deltaT;
TVectorD X(dimM+dimX);
TVectorD X_ERR(dimM+dimX);
TVectorD X_tot(dimM+23);
TVectorD X_ERR_tot(dimM+23);
std::vector<Double_t> chi2_iter, chi2_iter_0, chi2_iter_1, chi2_iter_2, chi2_iter_3, chi2_iter_4, chi2_iter_5, chi2_iter_6, chi2_iter_7, chi2_iter_8, chi2_iter_9, chi2_iter_10;
Int_t status, init, cov_status, nIter, STATUS;
ULong64_t eventNumber;

TVectorD g(dimC);
TMatrixD C(dimM+dimX,dimM+dimX);

// Constants
Double_t mtau = 1776.86;
Double_t mkaon = 493.677;
Double_t factor = 1;
Double_t factor_init;

Double_t initial_tolerance = pow(10,-10);

// Functions
void minimize( ROOT::Math::XYZPoint BV, int init, int species, double factor );
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
TVectorD x_initial_estimate_MLP( TVectorD m, ROOT::Math::XYZPoint BV, Int_t species );

Double_t chisquare( TVectorD x );
TVectorD exact_constraints( TVectorD x );
Double_t function_to_minimize( const Double_t* x_values );
Double_t g_sum(TVectorD g);
void sequence(ROOT::Math::XYZPoint BV);
void lowest_chi2(ROOT::Math::XYZPoint BV);
Bool_t firstTrial = true;
void run_block(Double_t eps, Int_t max_iter, Double_t r, Int_t init, Int_t species, ROOT::Math::XYZPoint BV_0);
void equations(TVectorD X);

void decay_fit_Penalty(int year, TString RECO_files, int species, int line)
{
    // Retrieve m and V from ntuple
    TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files, 1, line);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    Double_t m_vars[dimM];
    Double_t V_vars[dimM][dimM];
    Double_t BVx_0, BVy_0, BVz_0; // offline estimate for the BV is needed to have a first estimate for the unknown parameters using some of the initialisations

    for(int i = 0; i < dimM; i++)
    {
        t->SetBranchAddress(Form("df_m_%i",i+1), &m_vars[i]);
        for(int j = 0; j < dimM; j++)
        {
        t->SetBranchAddress(Form("df_V_%i_%i",i+1,j+1), &V_vars[i][j]);
        }
    }

    t->SetBranchAddress("Kp_RP_Z", &RPz);
    t->SetBranchAddress("Bp_ENDVERTEX_X", &BVx_0);
    t->SetBranchAddress("Bp_ENDVERTEX_Y", &BVy_0);
    t->SetBranchAddress("Bp_ENDVERTEX_Z", &BVz_0);    
    t->SetBranchAddress("eventNumber", &eventNumber);

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
        "df_E_conservation_in_DV1",
        "df_E_conservation_in_DV2",
        "df_E_conservation_in_BV"
    };

    for(int i = 0; i < dimM+23; i++)
    {
        tout->Branch(name_x[i], &X_tot(i));
        // tout->Branch(name_x_err[i], &X_ERR_tot(i));
    }

    for(int i = 0; i < dimC; i++)
    {
        tout->Branch(name_g[i], &g(i));
    }

    tout->Branch("df_Bp_M", &MB);
    tout->Branch("df_Bp_MERR", &MB_err);
    tout->Branch("df_status", &status);
    tout->Branch("df_cov_status", &cov_status);
    tout->Branch("df_chi2", &chi2);
    tout->Branch("df_tolerance", &sum_of_constraints);
    tout->Branch("df_nIter", &nIter);
    tout->Branch("df_time", &deltaT);
    tout->Branch("df_STATUS", &STATUS);
    tout->Branch("df_chi2_iter", &chi2_iter);

    for(int i = 0; i < dimC; i++)
    {
        tout->Branch(Form("df_g_%i",i), &g(i));
    }

    TStopwatch watch_total;

    // Loop over events
    UInt_t num_entries = t->GetEntries();
    for(int evt = 0; evt < 10; evt++)
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

        ROOT::Math::XYZPoint BV_0( BVx_0, BVy_0, BVz_0 ); // offline estimate of BV (necessary for some initialisations)

        Int_t max_iter = 10;
        int entry = 0;
        Double_t r = 10.;

        // 1st trial
        minimize(BV_0, -2, species, factor);
        x_current = X;
        status_current = status;
        factor *= r;
        firstTrial = false;
        Double_t x_vars_init[dimM+dimX];
        for(int i = 0; i < dimM+dimX; i++)
        {
            x_vars_init[i] = X(i); 
        }
        chi2_iter.push_back(function_to_minimize(x_vars_init));

        // Tolerance : 10-4
        run_block(pow(10,-6), 10, 10, -2, species, BV_0);

        if( ((status == 0) || (status == 1)) && (sum_of_constraints < pow(10,-4)) )
        {
            STATUS = 0;
        }
        else
        {
            STATUS = 1;
        }
        equations(X);

        // cout << "init = " << init << endl;
        cout << "status == " << status << endl;
        cout << "cov_status == " << cov_status << endl;
        cout << "STATUS == " << STATUS << endl;
        cout << "sum of constraints = " << sum_of_constraints << endl;
        cout << "MB = " << MB << " +/- " << MB_err << endl;
        cout << " ///////////////////////////////////////////////////////////////////////// " << endl;

        // Fill X_total
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

        Double_t BVz = X(dimM);
        Double_t ptau1z = X(dimM+1);
        Double_t ptau2z = X(dimM+2);
        Double_t MB_squared = X(dimM+3);

        // BV must be in the K+ trajectory
        Double_t BVx = RP.x() + (pK.x()/pK.z())*( BVz - RP.z() );
        Double_t BVy = RP.y() + (pK.y()/pK.z())*( BVz - RP.z() );
        ROOT::Math::XYZPoint BV( BVx, BVy, BVz );

        // ptau1 must point back to the BV
        Double_t ptau1x = ((DV1.x()-BV.x())/(DV1.z()-BV.z()))*ptau1z;
        Double_t ptau1y = ((DV1.y()-BV.y())/(DV1.z()-BV.z()))*ptau1z;
        ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );

        // ptau2 must point back to the BV
        Double_t ptau2x = ((DV2.x()-BV.x())/(DV2.z()-BV.z()))*ptau2z;
        Double_t ptau2y = ((DV2.y()-BV.y())/(DV2.z()-BV.z()))*ptau2z;
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

        // 3-momentum conservation in BV
        ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
        Double_t EB = sqrt( MB_squared + pB.Mag2() );

        for(int i = 0; i < dimM; i++)
        {
            X_tot(i) = X(i);
        }

        X_tot(dimM) = BV.x();
        X_tot(dimM+1) = BV.y();
        X_tot(dimM+2) = BV.z();
        X_tot(dimM+3) = pB.x();
        X_tot(dimM+4) = pB.y();
        X_tot(dimM+5) = pB.z();
        X_tot(dimM+6) = MB_squared;
        X_tot(dimM+7) = ptau1.x();
        X_tot(dimM+8) = ptau1.y();
        X_tot(dimM+9) = ptau1.z();
        X_tot(dimM+10) = Etau1;
        X_tot(dimM+11) = pnu1.x();
        X_tot(dimM+12) = pnu1.y();
        X_tot(dimM+13) = pnu1.z();
        X_tot(dimM+14) = Enu1;
        X_tot(dimM+15) = ptau2.x();
        X_tot(dimM+16) = ptau2.y();
        X_tot(dimM+17) = ptau2.z();
        X_tot(dimM+18) = Etau2;
        X_tot(dimM+19) = pnu2.x();
        X_tot(dimM+20) = pnu2.y();
        X_tot(dimM+21) = pnu2.z();
        X_tot(dimM+22) = Enu2;

        // X_tot.Print();

        watch.Stop();
        deltaT = watch.RealTime();
        tout->Fill();

        firstTrial = true;
        factor = 1;
        chi2_iter.clear();
    }
    watch_total.Stop();

    cout << "FINISHED in " << watch_total.RealTime() << " seconds." << endl;

    fout->cd();
    tout->Write();
    fout->Close();
}

void run_block(Double_t eps, Int_t max_iter, Double_t r, Int_t init, Int_t species, ROOT::Math::XYZPoint BV_0)
{
    Double_t x_vars[dimM+dimX];
    Int_t entry = 0;
    while( (sum_of_constraints > eps) && (entry < max_iter))
    {
        minimize(BV_0, init, species, factor);

        if(status > 1)
        {
            // return value before Minuit failure
            X = x_current;  
            status = status_current;
            cov_status = cov_status_current;
            break;
        }

        for(int i = 0; i < dimM+dimX; i++)
        {
            x_vars[i] = X(i); 
        }
        chi2_iter.push_back(function_to_minimize(x_vars));

        factor *= r;
        entry += 1;
        x_current = X;
        status_current = status;
        cov_status_current = cov_status;
    }
}

void minimize( ROOT::Math::XYZPoint BV, int init, int species, Double_t factor )
{
    // 1) This function initialises the vector of unkown parameters x with the result of analytical calculations
    // 2) It sets up the minimizer: defining the initial values and bounds on the x parameters  
    // 3) It builds the chi^2 function with x0 by calling the function chisquare 
    // 4) It does the minimisation and updates the values of the parameters that will be saved in a TTree 

    // 1) Initial values for x=(xm,xu,lambda)
    if(init == 0)
    {
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
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
        if(firstTrial)
        {
            x0 = x_initial_estimate_9( m, BV ); // Mixes Marseille (tau-) and K*tautau pions (tau+)
        }
        else
        {
            x0 = x9_current;
        }
    }
    else if(init == -2)
    {
        if(firstTrial)
        {   
            x0 = x_initial_estimate_MLP(m, BV, species);
        }
        else
        {
            x0 = x_current;
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
        "BVz",
        "Etau1",
        "Etau2",
        "MB_squared"
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
    min->Hesse();

    // Return / save results from the fit
    chi2 = min->MinValue();
    status = min->Status();
    cov_status = min->CovMatrixStatus();
    if(firstTrial)
    {
        nIter = min->NIterations();
    }
    else
    {
        nIter += min->NIterations();
    }
    
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

    Double_t MB_squared = X(dimM+3);
    if( MB_squared > 0 )
    {
        MB = sqrt( MB_squared );
    }
    else
    {
        MB = -sqrt( abs(MB_squared) );
    }

    Double_t dMB_squared = sqrt(C(dimM+3,dimM+3));
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
    Double_t BVz = x(dimM);
    Double_t ptau1z = x(dimM+1);
    Double_t ptau2z = x(dimM+2);
    Double_t MB_squared = x(dimM+3);

    // BV must be in the K+ trajectory
    Double_t BVx = RP.x() + (pK.x()/pK.z())*( BVz - RP.z() );
    Double_t BVy = RP.y() + (pK.y()/pK.z())*( BVz - RP.z() );
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz );

    // ptau1 must point back to the BV
    Double_t ptau1x = ((DV1.x()-BV.x())/(DV1.z()-BV.z()))*ptau1z;
    Double_t ptau1y = ((DV1.y()-BV.y())/(DV1.z()-BV.z()))*ptau1z;
    ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );

    // ptau2 must point back to the BV
    Double_t ptau2x = ((DV2.x()-BV.x())/(DV2.z()-BV.z()))*ptau2z;
    Double_t ptau2y = ((DV2.y()-BV.y())/(DV2.z()-BV.z()))*ptau2z;
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

    // 3-momentum conservation in BV
    ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
    Double_t EB = sqrt( MB_squared + pB.Mag2() );

    Double_t f = chisquare(x);
    TVectorD g = exact_constraints(x);

    TVectorD delta_g(dimC); // errors squared on the constraints
    delta_g(0) = 0.2890; // mm*MeV
    delta_g(1) = 0.3078; // mm*MeV
    delta_g(2) = 0.2595; // MeV
    delta_g(3) = 0.2567; // MeV
    delta_g(4) = 0.3000; // MeV

    for(int i = 0; i < dimC; i++)
    {
        delta_g(i) = delta_g(i)/factor;
    }

    // delta_g(0) = pB.z()*V(0,0)*pB.z() - pB.z()*V(0,2)*pB.x() - pB.x()*V(2,0)*pB.z() + pB.x()*pow(delta_xu(0),2)*pB.x();
    // delta_g(1) = pB.z()*V(1,1)*pB.z() - pB.z()*V(1,2)*pB.y() - pB.y()*V(2,1)*pB.z() + pB.y()*pow(delta_xu(0),2)*pB.y();
    // delta_g(2) = V(9,9) + 
    // delta_g(3) = V(16,16) + pow(delta_xu(18),2) + pow(delta_xu(22),2);
    // delta_g(4) = (pB.x()/EB)*pow(delta_xu(3),2)*(pB.x()/EB) + (pB.y()/EB)*pow(delta_xu(4),2)*(pB.y()/EB) + (pB.z()/EB)*pow(delta_xu(5),2)*(pB.z()/EB) + (1./(2*EB))*pow(delta_xu(6),2)*(1./(2*EB)) + pow(delta_xu(10),2) + pow(delta_xu(18),2) + (pK.x()/EK)*V(19,19)*(pK.x()/EK) + (pK.x()/EK)*V(19,20)*(pK.y()/EK) + (pK.y()/EK)*V(20,19)*(pK.x()/EK) + (pK.x()/EK)*V(19,21)*(pK.z()/EK) + (pK.z()/EK)*V(21,19)*(pK.x()/EK) + (pK.y()/EK)*V(20,21)*(pK.z()/EK) + (pK.z()/EK)*V(21,20)*(pK.y()/EK) + (pK.y()/EK)*V(20,20)*(pK.y()/EK) + (pK.z()/EK)*V(21,21)*(pK.z()/EK);

    // delta_g.Print();
    // g.Print();

    Double_t P = 0;
    for(int i = 0; i < dimC; i++)
    {
        P += pow(g(i),2)/pow(delta_g(i),2);
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
    Double_t BVz = x(dimM);
    Double_t ptau1z = x(dimM+1);
    Double_t ptau2z = x(dimM+2);
    Double_t MB_squared = x(dimM+3);

    // BV must be in the K+ trajectory
    Double_t BVx = RP.x() + (pK.x()/pK.z())*( BVz - RP.z() );
    Double_t BVy = RP.y() + (pK.y()/pK.z())*( BVz - RP.z() );
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz );

    // ptau1 must point back to the BV
    Double_t ptau1x = ((DV1.x()-BV.x())/(DV1.z()-BV.z()))*ptau1z;
    Double_t ptau1y = ((DV1.y()-BV.y())/(DV1.z()-BV.z()))*ptau1z;
    ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );

    // ptau2 must point back to the BV
    Double_t ptau2x = ((DV2.x()-BV.x())/(DV2.z()-BV.z()))*ptau2z;
    Double_t ptau2y = ((DV2.y()-BV.y())/(DV2.z()-BV.z()))*ptau2z;
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

    // 3-momentum conservation in BV
    ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK;
    Double_t EB = sqrt( MB_squared + pB.Mag2() );

    TVectorD g(dimC);
    // pB must point back to the PV (2)
    g(0) = pB.x()*(BV.z() - PV.z()) - pB.z()*(BV.x() - PV.x());
    g(1) = pB.y()*(BV.z() - PV.z()) - pB.z()*(BV.y() - PV.y());
    // Energy conservation in DV1, DV2 and in BV
    g(2) = Etau1 - E3pi1 - Enu1;
    g(3) = Etau2 - E3pi2 - Enu2;
    g(4) = EB - Etau1 - Etau2 - EK;

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

TVectorD x_initial_estimate_MLP( TVectorD m, ROOT::Math::XYZPoint BV, Int_t species )
{
    TVectorD x0(dimM+dimX);

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

    if((species == 10) || (species == 11) || (species == 12) || (species == 1) || (species == 2) || (species ==3))
    {
      weightfile_taup_PX = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP.weights.xml";
      weightfile_taup_PY = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP.weights.xml";
      weightfile_taup_PZ = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP.weights.xml";
      weightfile_taum_PX = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP.weights.xml";
      weightfile_taum_PY = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP.weights.xml";
      weightfile_taum_PZ = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP.weights.xml";
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

    x0(dimM) = BV.z();
    x0(dimM+1) = ptau1.z();
    x0(dimM+2) = ptau2.z();
    x0(dimM+3) = MB_squared;

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

void equations(TVectorD X)
{
    // Measured parameters
    ROOT::Math::XYZPoint PV( X(0), X(1), X(2) );
    ROOT::Math::XYZPoint DV1( X(3), X(4), X(5) );
    ROOT::Math::XYZVector p3pi1( X(6), X(7), X(8) );
    Double_t E3pi1 = X(9);
    ROOT::Math::XYZPoint DV2( X(10), X(11), X(12) );
    ROOT::Math::XYZVector p3pi2( X(13), X(14), X(15) );
    Double_t E3pi2 = X(16);
    ROOT::Math::XYZPoint RP( X(17), X(18), RPz );
    ROOT::Math::XYZVector pK( X(19), X(20), X(21) );

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
    Double_t c1 = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()));
    Double_t c2 = b*(DV1_t.x())/(DV2_t.x());
    Double_t d1 = ( (1+b)*(DV1_t.z() - PV_t.z()) + c1*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c2)*PV_t.x() );
    Double_t d2 = ( PV_t.x()*sqrt(pK_t.Mag2()) )/( (1+b)*DV1_t.x() - (1+c2)*PV_t.x() );
    Double_t e1 = c1 + c2*d1;
    Double_t e2 = c2*d2;

    Double_t x1 = p3pi1_t.x() + a1*p3pi1_t.y() + d1*p3pi1_t.z();
    Double_t x2 = b*p3pi2_t.x() + a2*b*p3pi2_t.y() + e1*p3pi2_t.z();

    Double_t g = 1 + pow(a1,2) + pow(d1,2) - pow(x1/E3pi1,2);
    Double_t h = 2*d1*d2 - ( pow(mtau,2) + pow(m3pi1,2) + 2*e2*p3pi1_t.z() )*(x1/pow(E3pi1,2) );
    Double_t j = pow(mtau,2) + pow(d2,2) - pow( ( pow(mtau,2) + pow(m3pi1,2) + 2*d2*p3pi1_t.z() )/(2*E3pi1), 2);

    Double_t m = pow(b,2) + pow(a2*b,2) + pow(e1,2) - pow(x2/E3pi2,2);
    Double_t n = 2*e1*e2 - ( pow(mtau,2) + pow(m3pi2,2) + 2*e2*p3pi2_t.z() )*(x2/pow(E3pi2,2) );
    Double_t o =  pow(mtau,2) + pow(e2,2) - pow( ( pow(mtau,2) + pow(m3pi2,2) + 2*e2*p3pi2_t.z() )/(2*E3pi2),2 );

    Double_t BVz = X(dimM);
    Double_t ptau1z = X(dimM+1);

    // BV must be in the K+ trajectory
    Double_t BVx = RP.x() + (pK.x()/pK.z())*( BVz - RP.z() );
    Double_t BVy = RP.y() + (pK.y()/pK.z())*( BVz - RP.z() );
    ROOT::Math::XYZPoint BV( BVx, BVy, BVz );

    // ptau1 must point back to the BV
    Double_t ptau1x = ((DV1.x()-BV.x())/(DV1.z()-BV.z()))*ptau1z;
    Double_t ptau1y = ((DV1.y()-BV.y())/(DV1.z()-BV.z()))*ptau1z;
    ROOT::Math::XYZVector ptau1( ptau1x, ptau1y, ptau1z );

    ROOT::Math::XYZVector ptau1_t = makeTransformation_vec( pK, RP, ptau1, false );

    Double_t x = ptau1_t.x();
    
    cout << "Eq1 value at solution = " << g*pow(x,2) + h*x + j << endl;
    cout << "Eq2 value at solution = " << m*pow(x,2) + n*x + o << endl;

}

/*
void lowest_chi2(Int_t max_iter, Double_t r, Double_t eps,  Int_t species, ROOT::Math::XYZPoint BV)
{
        Int_t entry = 0;
        Double_t factor = 1;

        // Lowest sum
        // Init 0
        minimize(BV, 0, species, factor);
        x_current = X;
        status_current = status;
        factor *= r;
        firstTrial = false;
        Double_t x_vars_init[dimM+dimX];
        for(int i = 0; i < dimM+dimX; i++)
        {
            x_vars_init[i] = X(i); 
        }
        chi2_iter_0.push_back(function_to_minimize(x_vars_init));
        run_block(eps, r, max_iter, 0, species, BV);

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


        // 1st trial
        

        // Tolerance : 10-4
        

        if( ((status == 0) || (status == 1)) && (sum_of_constraints < pow(10,-4)) )
        {
            STATUS = 0;
        }
        else
        {
            STATUS = 1;
        }

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
*/

/*
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
*/