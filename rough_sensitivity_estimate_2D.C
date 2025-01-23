std::vector<Double_t> B_Bp_Ktautau(Double_t BDT1, Double_t BDT2, TString the_case);
vector<double> range(double min, double max, size_t N);

void rough_sensitivity_estimate_2D(TString the_case)
{
    Int_t N = 20;
    Double_t BR_values[N][N], eps_sig_values[N][N], N_sig_values[N][N];
    Double_t BDT1_values[N][N], BDT2_values[N][N];

    vector<double> bdts = range(0, 1, N);
    Double_t BDT1, BDT2;
    Double_t BDT1_flat_values[N*N], BDT2_flat_values[N*N], BR_flat_values[N*N];

    Int_t a = 0;
    for(int i = 0; i < N; i++)
    {
        BDT1 = bdts[i];
        for(int j = 0; j < N; j++)
        {
            BDT2 = bdts[j];

            BDT1_values[i][j] = BDT1;
            BDT2_values[i][j] = BDT2;

            std::vector<Double_t> results = B_Bp_Ktautau(BDT1, BDT2, the_case);
            N_sig_values[i][j]= results[0];
            eps_sig_values[i][j] = results[1];
            BR_values[i][j] = results[2];
            if(eps_sig_values[i][j] == 0)
            {
                BR_values[i][j] = 0;
            }

            // cout << "BDT1 = " << BDT1 << " | BDT2 = " << BDT2 << " | BR = " << BR_values[i][j] << " | N_sig = " << N_sig_values[i][j] << " | eps_sig = " << eps_sig_values[i][j] << endl;
        
            BDT1_flat_values[a] = BDT1;
            BDT2_flat_values[a] = BDT2;
            BR_flat_values[a] = BR_values[i][j];

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

    TString name;
    if(the_case == "zero_bkg_3pi3pi")
    {
        name = "Case 1A: 3#pi3#pi only, zero background (B=0)";
    }
    else if(the_case == "zero_bkg_all_mc")
    {
        name = "Case 2A: All MC, zero background (B=0)";
    }

    TCanvas c;
    c.cd();
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


std::vector<Double_t> B_Bp_Ktautau(Double_t BDT1, Double_t BDT2, TString the_case)
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
    Double_t eps_sig_norm_2016 = (14.50/100)*(2.256/100)*(55.8/100)*(99.97/100)*(11.39/100);
    Double_t eps_sig_norm_2017 = (14.48/100)*(2.278/100)*(58.5/100)*(99.98/100)*(11.91/100);
    Double_t eps_sig_norm_2018 = (14.58/100)*(2.286/100)*(50.8/100)*(99.97/100)*(10.97/100);
    Double_t eps_sig_norm = (1.6/5.4)*eps_sig_norm_2016 + (1.7/5.4)*eps_sig_norm_2017 + (2.1/5.4)*eps_sig_norm_2018;
    // cout << "Norm pre-BDT sig eff = " << eps_sig_norm*100 << endl;
    // 4. Ktautau signal pre-selection efficiency
    Double_t eps_sig_ktautau_3pi3pi_2016 = (5.318/100)*(0.645/100)*(37.5/100)*(93.0/100)*(96.3/100)*(79.0/100);
    Double_t eps_sig_ktautau_3pi3pi_2017 = (5.321/100)*(0.654/100)*(40.9/100)*(90.0/100)*(96.3/100)*(79.0/100);
    Double_t eps_sig_ktautau_3pi3pi_2018 = (5.330/100)*(0.651/100)*(34.8/100)*(90.2/100)*(96.5/100)*(79.1/100);

    Double_t eps_sig_ktautau_all_mc_2016 = (5.318/100)*(0.628/100)*(36.8/100)*(85.1/100)*(96.6/100)*(75.8/100);
    Double_t eps_sig_ktautau_all_mc_2017 = (5.321/100)*(0.649/100)*(40.6/100)*(85.2/100)*(96.5/100)*(75.4/100);
    Double_t eps_sig_ktautau_all_mc_2018 = (5.330/100)*(0.648/100)*(34.4/100)*(84.8/100)*(96.8/100)*(75.4/100);

    Double_t eps_sig_ktautau_3pi3pi = (1.6/5.4)*eps_sig_ktautau_3pi3pi_2016 + (1.7/5.4)*eps_sig_ktautau_3pi3pi_2017 + (2.1/5.4)*eps_sig_ktautau_3pi3pi_2018;
    Double_t eps_sig_ktautau_all_mc = (1.6/5.4)*eps_sig_ktautau_all_mc_2016 + (1.7/5.4)*eps_sig_ktautau_all_mc_2017 + (2.1/5.4)*eps_sig_ktautau_all_mc_2018;

    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC, 2016) = " << eps_sig_ktautau_3pi3pi_2016*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC, 2017) = " << eps_sig_ktautau_3pi3pi_2017*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC, 2018) = " << eps_sig_ktautau_3pi3pi_2018*100 << endl;

    // cout << "Ktautau pre-BDT sig eff (All MC, 2016) = " << eps_sig_ktautau_all_mc_2016*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (All MC, 2017) = " << eps_sig_ktautau_all_mc_2017*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (All MC, 2018) = " << eps_sig_ktautau_all_mc_2018*100 << endl;

    // cout << "Ktautau pre-BDT sig eff (3pi3pi MC) = " << eps_sig_ktautau_3pi3pi*100 << endl;
    // cout << "Ktautau pre-BDT sig eff (all MC) = " << eps_sig_ktautau_all_mc*100 << endl;

    // 5. Estimate for the number of Ktautau signal events in the signal region
    Double_t N_sig_Ktautau, eps_sig_ktautau;
    if(the_case == "zero_bkg_3pi3pi")
    {
        N_sig_Ktautau = 3.;
        eps_sig_ktautau = eps_sig_ktautau_3pi3pi;
    }
    else if(the_case == "zero_bkg_all_mc")
    {
        N_sig_Ktautau = 3.;
        eps_sig_ktautau = eps_sig_ktautau_all_mc;
    }

    ///////////////////////////////////////////////////////////////// Input files /////////////////////////////////////////////////////////////////////////////////////////////////////
    TFileCollection *fc_3pi3pi_2016 = new TFileCollection("fc_3pi3pi_2016", "fc_3pi3pi_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_10/fit_results.txt");
    TFileCollection *fc_3pi3pi_2017 = new TFileCollection("fc_3pi3pi_2017", "fc_3pi3pi_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_10/fit_results.txt");
    TFileCollection *fc_3pi3pi_2018 = new TFileCollection("fc_3pi3pi_2018", "fc_3pi3pi_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_10/fit_results.txt");

    TChain *t_3pi3pi_2016 = new TChain("DecayTree");
    TChain *t_3pi3pi_2017 = new TChain("DecayTree");
    TChain *t_3pi3pi_2018 = new TChain("DecayTree");

    t_3pi3pi_2016->AddFileInfoList((TCollection*)fc_3pi3pi_2016->GetList());
    t_3pi3pi_2017->AddFileInfoList((TCollection*)fc_3pi3pi_2017->GetList());
    t_3pi3pi_2018->AddFileInfoList((TCollection*)fc_3pi3pi_2018->GetList());

    t_3pi3pi_2016->GetEntries();
    t_3pi3pi_2017->GetEntries();
    t_3pi3pi_2018->GetEntries();

    t_3pi3pi_2016->Add(t_3pi3pi_2017);
    t_3pi3pi_2016->Add(t_3pi3pi_2018);

    TFileCollection *fc1_3pi3pi_2016 = new TFileCollection("fc1_3pi3pi_2016", "fc1_3pi3pi_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_10/bdt_output.txt");
    TFileCollection *fc1_3pi3pi_2017 = new TFileCollection("fc1_3pi3pi_2017", "fc1_3pi3pi_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_10/bdt_output.txt");
    TFileCollection *fc1_3pi3pi_2018 = new TFileCollection("fc1_3pi3pi_2018", "fc1_3pi3pi_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_10/bdt_output.txt");

    TChain *t1_3pi3pi_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_3pi3pi_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_3pi3pi_2018 = new TChain("XGBoost/DecayTree");

    t1_3pi3pi_2016->AddFileInfoList((TCollection*)fc1_3pi3pi_2016->GetList());
    t1_3pi3pi_2017->AddFileInfoList((TCollection*)fc1_3pi3pi_2017->GetList());
    t1_3pi3pi_2018->AddFileInfoList((TCollection*)fc1_3pi3pi_2018->GetList());

    t1_3pi3pi_2016->GetEntries();
    t1_3pi3pi_2017->GetEntries();
    t1_3pi3pi_2018->GetEntries();

    t1_3pi3pi_2016->Add(t1_3pi3pi_2017);
    t1_3pi3pi_2016->Add(t1_3pi3pi_2018);

    t_3pi3pi_2016->AddFriend(t1_3pi3pi_2016);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Double_t bdt_sig_eff_den = t_3pi3pi_2016->GetEntries("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800)");
    Double_t bdt1_sig_eff_num = t_3pi3pi_2016->GetEntries(Form("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800) && (BDT1 > %f)",BDT1));
    Double_t bdt2_sig_eff_num = t_3pi3pi_2016->GetEntries(Form("(df_status==0) && (df_Bp_M > 4700) && (df_Bp_M < 5800) && (BDT2 > %f)",BDT2));

    Double_t eps_sig_ktautau_bdt1 = bdt1_sig_eff_num/bdt_sig_eff_den;
    Double_t eps_sig_ktautau_bdt2 = bdt2_sig_eff_num/bdt_sig_eff_den;

    eps_sig_ktautau *= eps_sig_ktautau_bdt1*eps_sig_ktautau_bdt2;

    Double_t B_Bp_Ktautau;
    if(the_case == "zero_bkg_3pi3pi")
    {
        B_Bp_Ktautau = (N_sig_Ktautau/eps_sig_ktautau)*(eps_sig_norm/N_sig_norm)*((B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)/pow(B_tau_3pi_nu,2));
    }
    else if(the_case == "zero_bkg_all_mc")
    {
        B_Bp_Ktautau = (N_sig_Ktautau/eps_sig_ktautau)*(eps_sig_norm/N_sig_norm)*((B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)/pow(B_tau_3pi_nu+B_tau_3pi_pi0_nu,2));
    }

    std::vector<Double_t> results;
    results.push_back(N_sig_Ktautau);
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