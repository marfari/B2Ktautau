std::vector<Double_t> B_Bp_Ktautau(Double_t BDT1, Double_t BDT2, TString the_case, TTree* t_3pi3pipi0_2016, TTree* t_3pi3pi2pi0_2016, TTree* t_sig, TTree* t_bkg);

vector<double> range(double min, double max, size_t N);

void rough_sensitivity_estimate_2D(TString the_case)
{
    ///////////////////////////////////////////////////////////////// Input files /////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3pi3pi MC
    TFile* f_sig_iso = new TFile("/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/isolation_test_dataset_sig.root");
    TFile* f_sig_kin = new TFile("/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/topology_test_dataset_sig.root");

    TTree* t_sig = (TTree*)f_sig_iso->Get("XGBoost/DecayTree");
    TTree* t_sig1 = (TTree*)f_sig_kin->Get("XGBoost/DecayTree");
    t_sig->AddFriend(t_sig1);

    // 3pi3pi pi0 MC
    TFileCollection *fc_3pi3pipi0_2016 = new TFileCollection("fc_3pi3pipi0_2016", "fc_3pi3pipi0_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_11/fit_results.txt");
    TFileCollection *fc_3pi3pipi0_2017 = new TFileCollection("fc_3pi3pipi0_2017", "fc_3pi3pipi0_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_11/fit_results.txt");
    TFileCollection *fc_3pi3pipi0_2018 = new TFileCollection("fc_3pi3pipi0_2018", "fc_3pi3pipi0_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_11/fit_results.txt");

    TChain *t_3pi3pipi0_2016 = new TChain("DecayTree");
    TChain *t_3pi3pipi0_2017 = new TChain("DecayTree");
    TChain *t_3pi3pipi0_2018 = new TChain("DecayTree");

    t_3pi3pipi0_2016->AddFileInfoList((TCollection*)fc_3pi3pipi0_2016->GetList());
    t_3pi3pipi0_2017->AddFileInfoList((TCollection*)fc_3pi3pipi0_2017->GetList());
    t_3pi3pipi0_2018->AddFileInfoList((TCollection*)fc_3pi3pipi0_2018->GetList());

    t_3pi3pipi0_2016->GetEntries();
    t_3pi3pipi0_2017->GetEntries();
    t_3pi3pipi0_2018->GetEntries();

    t_3pi3pipi0_2016->Add(t_3pi3pipi0_2017);
    t_3pi3pipi0_2016->Add(t_3pi3pipi0_2018);

    TFileCollection *fc1_3pi3pipi0_2016 = new TFileCollection("fc1_3pi3pipi0_2016", "fc1_3pi3pipi0_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_11/bdt_output.txt");
    TFileCollection *fc1_3pi3pipi0_2017 = new TFileCollection("fc1_3pi3pipi0_2017", "fc1_3pi3pipi0_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_11/bdt_output.txt");
    TFileCollection *fc1_3pi3pipi0_2018 = new TFileCollection("fc1_3pi3pipi0_2018", "fc1_3pi3pipi0_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_11/bdt_output.txt");

    TChain *t1_3pi3pipi0_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_3pi3pipi0_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_3pi3pipi0_2018 = new TChain("XGBoost/DecayTree");

    t1_3pi3pipi0_2016->AddFileInfoList((TCollection*)fc1_3pi3pipi0_2016->GetList());
    t1_3pi3pipi0_2017->AddFileInfoList((TCollection*)fc1_3pi3pipi0_2017->GetList());
    t1_3pi3pipi0_2018->AddFileInfoList((TCollection*)fc1_3pi3pipi0_2018->GetList());

    t1_3pi3pipi0_2016->GetEntries();
    t1_3pi3pipi0_2017->GetEntries();
    t1_3pi3pipi0_2018->GetEntries();

    t1_3pi3pipi0_2016->Add(t1_3pi3pipi0_2017);
    t1_3pi3pipi0_2016->Add(t1_3pi3pipi0_2018);

    t_3pi3pipi0_2016->AddFriend(t1_3pi3pipi0_2016);

    // 3pi3pi 2pi0 MC
    TFileCollection *fc_3pi3pi2pi0_2016 = new TFileCollection("fc_3pi3pi2pi0_2016", "fc_3pi3pi2pi0_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_12/fit_results.txt");
    TFileCollection *fc_3pi3pi2pi0_2017 = new TFileCollection("fc_3pi3pi2pi0_2017", "fc_3pi3pi2pi0_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_12/fit_results.txt");
    TFileCollection *fc_3pi3pi2pi0_2018 = new TFileCollection("fc_3pi3pi2pi0_2018", "fc_3pi3pi2pi0_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_12/fit_results.txt");

    TChain *t_3pi3pi2pi0_2016 = new TChain("DecayTree");
    TChain *t_3pi3pi2pi0_2017 = new TChain("DecayTree");
    TChain *t_3pi3pi2pi0_2018 = new TChain("DecayTree");

    t_3pi3pi2pi0_2016->AddFileInfoList((TCollection*)fc_3pi3pi2pi0_2016->GetList());
    t_3pi3pi2pi0_2017->AddFileInfoList((TCollection*)fc_3pi3pi2pi0_2017->GetList());
    t_3pi3pi2pi0_2018->AddFileInfoList((TCollection*)fc_3pi3pi2pi0_2018->GetList());

    t_3pi3pi2pi0_2016->GetEntries();
    t_3pi3pi2pi0_2017->GetEntries();
    t_3pi3pi2pi0_2018->GetEntries();

    t_3pi3pi2pi0_2016->Add(t_3pi3pi2pi0_2017);
    t_3pi3pi2pi0_2016->Add(t_3pi3pi2pi0_2018);

    TFileCollection *fc1_3pi3pi2pi0_2016 = new TFileCollection("fc1_3pi3pi2pi0_2016", "fc1_3pi3pi2pi0_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_12/bdt_output.txt");
    TFileCollection *fc1_3pi3pi2pi0_2017 = new TFileCollection("fc1_3pi3pi2pi0_2017", "fc1_3pi3pi2pi0_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_12/bdt_output.txt");
    TFileCollection *fc1_3pi3pi2pi0_2018 = new TFileCollection("fc1_3pi3pi2pi0_2018", "fc1_3pi3pi2pi0_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_12/bdt_output.txt");

    TChain *t1_3pi3pi2pi0_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_3pi3pi2pi0_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_3pi3pi2pi0_2018 = new TChain("XGBoost/DecayTree");

    t1_3pi3pi2pi0_2016->AddFileInfoList((TCollection*)fc1_3pi3pi2pi0_2016->GetList());
    t1_3pi3pi2pi0_2017->AddFileInfoList((TCollection*)fc1_3pi3pi2pi0_2017->GetList());
    t1_3pi3pi2pi0_2018->AddFileInfoList((TCollection*)fc1_3pi3pi2pi0_2018->GetList());

    t1_3pi3pi2pi0_2016->GetEntries();
    t1_3pi3pi2pi0_2017->GetEntries();
    t1_3pi3pi2pi0_2018->GetEntries();

    t1_3pi3pi2pi0_2016->Add(t1_3pi3pi2pi0_2017);
    t1_3pi3pi2pi0_2016->Add(t1_3pi3pi2pi0_2018);

    t_3pi3pi2pi0_2016->AddFriend(t1_3pi3pi2pi0_2016);

    // Background test sample
    TFile* f_bkg_iso = new TFile("/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/isolation_test_dataset_bkg.root");
    TFile* f_bkg_kin = new TFile("/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/topology_test_dataset_bkg.root");

    TTree* t_bkg = (TTree*)f_bkg_iso->Get("XGBoost/DecayTree");
    TTree* t_bkg1 = (TTree*)f_bkg_kin->Get("XGBoost/DecayTree");
    t_bkg->AddFriend(t_bkg1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Int_t N = 30;
    Double_t BR_values[N][N], eps_sig_values[N][N], N_sig_values[N][N], N_bkg_values[N][N];
    Double_t BDT1_values[N][N], BDT2_values[N][N];

    vector<double> bdts = range(0, 1, N);
    Double_t BDT1, BDT2;
    Double_t BDT1_flat_values[N*N], BDT2_flat_values[N*N], BR_flat_values[N*N], N_sig_flat_values[N*N], N_bkg_flat_values[N*N], eps_sig_flat_values[N*N];

    Int_t a = 0;
    for(int i = 0; i < N; i++)
    {
        BDT1 = bdts[i];
        for(int j = 0; j < N; j++)
        {
            BDT2 = bdts[j];

            BDT1_values[i][j] = BDT1;
            BDT2_values[i][j] = BDT2;

            std::vector<Double_t> results = B_Bp_Ktautau(BDT1, BDT2, the_case, t_3pi3pipi0_2016, t_3pi3pi2pi0_2016, t_sig, t_bkg);
            N_sig_values[i][j]= results[0];
            N_bkg_values[i][j] = results[1];
            eps_sig_values[i][j] = results[2];
            BR_values[i][j] = results[3];

            if(eps_sig_values[i][j] == 0)
            {
                BR_values[i][j] = 0;
            }

            cout << "BDT1 = " << BDT1 << " | BDT2 = " << BDT2 << " | BR = " << BR_values[i][j] << " | N_sig = " << N_sig_values[i][j] << " | N_bkg = " << N_bkg_values[i][j]  << " | eps_sig = " << eps_sig_values[i][j] << endl;
        
            BDT1_flat_values[a] = BDT1;
            BDT2_flat_values[a] = BDT2;
            BR_flat_values[a] = BR_values[i][j];
            N_sig_flat_values[a] = N_sig_values[i][j];
            N_bkg_flat_values[a] = N_bkg_values[i][j];
            eps_sig_flat_values[a] = eps_sig_values[i][j];

            a +=1;
        }
    }

    // Find BDT cut values that minimum the branching fraction
    Int_t index = 0;
    for(int i = 0; i < N*N; i++)
    {
        if( (BR_flat_values[i] < BR_flat_values[index]) && (BR_flat_values[i] != 0) )
        {
            index = i;
        }
    }

    cout << Form("The BDT cuts (BDT1,BDT2) = (%.3lf,%.3lf) minimuse the BR(B->Ktautau) to %.10lf at 95& C.L.", BDT1_flat_values[index], BDT2_flat_values[index], BR_flat_values[index]) << endl;
    cout << Form("N_sig = %.1f", N_sig_flat_values[index]) << endl;
    cout << Form("N_bkg = %.1f", N_bkg_flat_values[index]) << endl;
    cout << Form("eps_sig = %.4f percent", eps_sig_flat_values[index]*100) << endl;

    TString name;
    if(the_case == "zero_bkg_3pi3pi")
    {
        name = "Case 1A: 3#pi3#pi only, zero background (B=0)";
    }
    else if(the_case == "zero_bkg_all_mc")
    {
        name = "Case 2A: All MC, zero background (B=0)";
    }
    else if(the_case == "comb_bkg_3pi3pi")
    {
        name = "Case 1B : 3#pi3#pi only, combinatorial background";
    }
    else if(the_case == "comb_bkg_all_mc")
    {
        name = "Case 2B: All MC, combinatorial background";
    }

    TCanvas c;
    c.cd();
    if( !((the_case == "zero_bkg_3pi3pi") || (the_case == "zero_bkg_all_mc")) )
    {
        c.SetTheta(45.);
    }
    TGraph2D *g = new TGraph2D("g", "g", N*N, BDT1_flat_values, BDT2_flat_values, BR_flat_values);
    g->SetTitle(name+"; Isolation: BDT1 > x; Topology: BDT2 > x; B(B^{+} #rightarrow K^{+} #tau^{+} #tau^{-})");
    gStyle->SetPalette(1);
    g->SetMarkerStyle(20);
    g->Draw("pcol");
    // auto mark = new TMarker(BDT1_flat_values[index],BDT2_flat_values[index],BR_flat_values[index]);
    // mark->SetMarkerStyle(29);
    // mark->Draw("same");
    c.SaveAs("/panfs/felician/B2Ktautau/workflow/rough_sensitivity_2D/BR_vs_bdt1_bdt2_"+the_case+".pdf");

}


std::vector<Double_t> B_Bp_Ktautau(Double_t BDT1, Double_t BDT2, TString the_case, TTree* t_3pi3pipi0_2016, TTree* t_3pi3pi2pi0_2016, TTree* t_sig, TTree* t_bkg)
{
    // External inputs:
    // 1. Branching fractions (from PDG)
    Double_t B_tau_3pi_nu = 9.31/100; // +/- 0.05/100
    Double_t B_tau_3pi_pi0_nu = 4.62/100; // +/- 0.05/100
    Double_t B_Bp_D0bar_Dsp = 9.0/1000; // +/- 0.9/1000
    Double_t B_D0bar_Kp_pi = 3.947/100; // +/- 0.030/100
    Double_t B_Dsp_Kp_Km_pi = 5.37/100; // +/- 0.10/100
    // 2. Signal yield of normalisation channel (from mass fit)
    Double_t N_sig_norm_2016 = 9737.;
    Double_t N_sig_norm_2017 = 10630.;
    Double_t N_sig_norm_2018 = 12273.;
    Double_t N_sig_norm = N_sig_norm_2016+N_sig_norm_2017+N_sig_norm_2018;
    // 3. Efficiency of normalisation channel
    Double_t eps_sig_norm_2016 = (14.50/100)*(2.256/100)*(55.81/100)*(64.11/100)*(99.94/100)*(99.6/100);
    Double_t eps_sig_norm_2017 = (14.48/100)*(2.278/100)*(58.54/100)*(62.15/100)*(99.97/100)*(99.6/100);
    Double_t eps_sig_norm_2018 = (14.58/100)*(2.286/100)*(50.76/100)*(61.94/100)*(99.97/100)*(99.6/100);
    Double_t eps_sig_norm = (10.8/34.2)*eps_sig_norm_2016 + (12./34.2)*eps_sig_norm_2017 + (11.4/34.2)*eps_sig_norm_2018;
    // cout << "2016 = " << eps_sig_norm_2016*100 << endl;
    // cout << "2017 = " << eps_sig_norm_2017*100 << endl;
    // cout << "2018 = " << eps_sig_norm_2018*100 << endl;
    // cout << "Norm pre-BDT sig eff = " << eps_sig_norm*100 << endl;
    // 4. Ktautau signal pre-selection efficiency
    Double_t eps_acc_2016 = (5.318/100);
    Double_t eps_acc_2017 = (5.321/100);
    Double_t eps_acc_2018 = (5.330/100);

    Double_t eps_strip_2016 = (1.028/100);
    Double_t eps_strip_2017 = (1.056/100);
    Double_t eps_strip_2018 = (1.050/100);

    Double_t eps_sig_ktautau_3pi3pi_2016 = eps_acc_2016*eps_strip_2016*(62.67/100)*(37.47/100)*(92.60/100)*(96.29/100)*(100./100);
    Double_t eps_sig_ktautau_3pi3pi_2017 = eps_acc_2017*eps_strip_2017*(61.63/100)*(41.17/100)*(89.63/100)*(96.60/100)*(100./100);
    Double_t eps_sig_ktautau_3pi3pi_2018 = eps_acc_2018*eps_strip_2018*(62.00/100)*(34.79/100)*(90.19/100)*(96.57/100)*(100./100);

    Double_t eps_sig_ktautau_all_mc_2016 = eps_acc_2016*eps_strip_2016*(61.12/100)*(36.75/100)*(85.12/100)*(96.57/100)*(100./100);
    Double_t eps_sig_ktautau_all_mc_2017 = eps_acc_2017*eps_strip_2017*(61.49/100)*(40.61/100)*(85.16/100)*(96.78/100)*(99.5/100);
    Double_t eps_sig_ktautau_all_mc_2018 = eps_acc_2018*eps_strip_2018*(61.66/100)*(34.44/100)*(84.76/100)*(96.79/100)*(99.5/100);

    Double_t eps_sig_ktautau_3pi3pi = (28.5/165)*eps_sig_ktautau_3pi3pi_2016 + (60.0/165)*eps_sig_ktautau_3pi3pi_2017 + (76.5/165)*eps_sig_ktautau_3pi3pi_2018;
    Double_t eps_sig_ktautau_all_mc = (28.5/165)*eps_sig_ktautau_all_mc_2016 + (60.0/165)*eps_sig_ktautau_all_mc_2017 + (76.5/165)*eps_sig_ktautau_all_mc_2018;

    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC, 2016) = " << eps_sig_ktautau_3pi3pi_2016*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC, 2017) = " << eps_sig_ktautau_3pi3pi_2017*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC, 2018) = " << eps_sig_ktautau_3pi3pi_2018*100 << endl;

    // cout << "Ktautau pre-BDT sig eff (All MC, 2016) = " << eps_sig_ktautau_all_mc_2016*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (All MC, 2017) = " << eps_sig_ktautau_all_mc_2017*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (All MC, 2018) = " << eps_sig_ktautau_all_mc_2018*100 << endl;

    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC) = " << eps_sig_ktautau_3pi3pi*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (all MC) = " << eps_sig_ktautau_all_mc*100 << endl;
    
    // 5. Estimate for the number of Ktautau signal events in the signal region (pre-BDT)
    Double_t N_sig_Ktautau = 0.;
    Double_t eps_sig_ktautau = 0.; 
    Double_t N_bkg_Ktautau = 0.;

    Bool_t use3pi3piMC = false;
    Bool_t zeroBkg = false;
    Bool_t combBkg = false;
    Bool_t physBkg = false;

    if((the_case == "zero_bkg_3pi3pi") || (the_case == "comb_bkg_3pi3pi") || (the_case == "phys_bkg_3pi3pi"))
    {
        use3pi3piMC = true;
    }
    if((the_case == "zero_bkg_3pi3pi") || (the_case == "zero_bkg_all_mc"))
    {
        zeroBkg = true;
    }
    if((the_case == "comb_bkg_3pi3pi") || (the_case == "comb_bkg_all_mc"))
    {
        combBkg = true;
    }
    if((the_case == "phys_bkg_3pi3pi") || (the_case == "phys_bkg_all_mc"))
    {
        physBkg = true;
    }


    if(use3pi3piMC)
    {
        eps_sig_ktautau = eps_sig_ktautau_3pi3pi;
    }
    else
    {
        eps_sig_ktautau = eps_sig_ktautau_all_mc;
    }

    if(zeroBkg)
    {
        N_sig_Ktautau = 3.;
    }
    else if(combBkg)
    {
        Double_t N_bkg_ktautau_2016 = 25939425.;
        Double_t N_bkg_ktautau_2017 = 27336738.;
        Double_t N_bkg_ktautau_2018 = 30862612.;
        N_bkg_Ktautau = N_bkg_ktautau_2016+N_bkg_ktautau_2017+N_bkg_ktautau_2018;
    }

    Double_t bdt_sig_eff_den, bdt_sig_eff_num;
    if(use3pi3piMC)
    {
        bdt_sig_eff_den = t_sig->GetEntries();
        bdt_sig_eff_num = t_sig->GetEntries(Form("(BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2));
    } 
    else
    {
        bdt_sig_eff_den = (t_3pi3pipi0_2016->GetEntries("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800)")) + (t_3pi3pi2pi0_2016->GetEntries("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800)")) + (t_sig->GetEntries());
        bdt_sig_eff_num = (t_3pi3pipi0_2016->GetEntries(Form("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800) && (BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2))) + (t_3pi3pi2pi0_2016->GetEntries(Form("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800) && (BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2))) + (t_sig->GetEntries(Form("(BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2)));
    }

    Double_t eps_sig_ktautau_bdt = bdt_sig_eff_num/bdt_sig_eff_den;

    eps_sig_ktautau *= eps_sig_ktautau_bdt;

    Double_t eps_bkg_ktautau_bdt;
    if(combBkg)
    {
        Double_t bdt_bkg_eff_den = t_bkg->GetEntries();
        Double_t bdt_bkg_eff_num = t_bkg->GetEntries(Form("(BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2));

        eps_bkg_ktautau_bdt = bdt_bkg_eff_num/bdt_bkg_eff_den;

        N_bkg_Ktautau *= eps_bkg_ktautau_bdt;

        N_sig_Ktautau = 0.5*( pow(1.64,2) + sqrt(pow(1.64,4) + 4*pow(1.64,2)*N_bkg_Ktautau) ); 
    }

    Double_t B_Bp_Ktautau;
    if(use3pi3piMC)
    {
        B_Bp_Ktautau = (N_sig_Ktautau/eps_sig_ktautau)*(eps_sig_norm/N_sig_norm)*((B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)/pow(B_tau_3pi_nu,2));
    }
    else
    {
        B_Bp_Ktautau = (N_sig_Ktautau/eps_sig_ktautau)*(eps_sig_norm/N_sig_norm)*((B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)/pow(B_tau_3pi_nu+B_tau_3pi_pi0_nu,2));
    }
    
    std::vector<Double_t> results;
    results.push_back(N_sig_Ktautau);
    results.push_back(N_bkg_Ktautau);
    results.push_back(eps_sig_ktautau);
    results.push_back(B_Bp_Ktautau);

    return results;
}

vector<double> range(double min, double max, size_t N) {
    vector<double> range;
    double delta = (max-min)/double(N-1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}