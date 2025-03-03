std::vector<Double_t> B_Bp_Ktautau(Double_t BDT1, Double_t BDT2, TString the_case, TTree* t_mc, TTree* t_gen_3pi3pi, TTree* t_gen_all, TTree* t_norm_mc, TTree* t_norm_gen, TTree* t_ws);

vector<double> range(double min, double max, size_t N);

void rough_sensitivity_estimate_2D(TString the_case)
{
    ///////////////////////////////////////////////////////////////// Input files /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Ktautau MC after all selections
    TFile* f_mc = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_1/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_mc = (TTree*)f_mc->Get("DecayTree");

    // Ktautau gen MC
    TFileCollection* fc_gen_2016 = new TFileCollection("fc_gen_2016", "fc_gen_2016", "Files_on_grid/MC_2016.txt");
    TFileCollection* fc_gen_2017 = new TFileCollection("fc_gen_2017", "fc_gen_2017", "Files_on_grid/MC_2017.txt");
    TFileCollection* fc_gen_2018 = new TFileCollection("fc_gen_2018", "fc_gen_2018", "Files_on_grid/MC_2018.txt");

    // 2016
    TChain* t_gen_3pi3pi = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
    TChain* t_gen_3pi3pi_pi0_2016 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
    TChain* t_gen_3pipi0_3pi_2016 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
    TChain* t_gen_3pi3pi2pi0_2016 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
    TChain* t_gen_all = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");

    t_gen_3pi3pi->AddFileInfoList((TCollection*)fc_gen_2016->GetList());
    t_gen_3pi3pi_pi0_2016->AddFileInfoList((TCollection*)fc_gen_2016->GetList());
    t_gen_3pipi0_3pi_2016->AddFileInfoList((TCollection*)fc_gen_2016->GetList());
    t_gen_3pi3pi2pi0_2016->AddFileInfoList((TCollection*)fc_gen_2016->GetList());
    t_gen_all->AddFileInfoList((TCollection*)fc_gen_2016->GetList());

    cout << "Ktautau gen 2016" << endl;
    cout << t_gen_3pi3pi->GetEntries() << endl;
    cout << t_gen_3pi3pi_pi0_2016->GetEntries() << endl;
    cout << t_gen_3pipi0_3pi_2016->GetEntries() << endl;
    cout << t_gen_3pi3pi2pi0_2016->GetEntries() << endl;
    cout << t_gen_all->GetEntries() << endl;

    t_gen_all->Add(t_gen_3pi3pi_pi0_2016);
    t_gen_all->Add(t_gen_3pipi0_3pi_2016);
    t_gen_all->Add(t_gen_3pi3pi2pi0_2016);

    // 2017
    TChain* t_gen_3pi3pi_2017 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
    TChain* t_gen_3pi3pi_pi0_2017 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
    TChain* t_gen_3pipi0_3pi_2017 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
    TChain* t_gen_3pi3pi2pi0_2017 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
    TChain* t_gen_all_2017 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");

    t_gen_3pi3pi_2017->AddFileInfoList((TCollection*)fc_gen_2017->GetList());
    t_gen_3pi3pi_pi0_2017->AddFileInfoList((TCollection*)fc_gen_2017->GetList());
    t_gen_3pipi0_3pi_2017->AddFileInfoList((TCollection*)fc_gen_2017->GetList());
    t_gen_3pi3pi2pi0_2017->AddFileInfoList((TCollection*)fc_gen_2017->GetList());
    t_gen_all_2017->AddFileInfoList((TCollection*)fc_gen_2017->GetList());

    cout << "Ktautau gen 2017" << endl;
    cout << t_gen_3pi3pi_2017->GetEntries() << endl;
    cout << t_gen_3pi3pi_pi0_2017->GetEntries() << endl;
    cout << t_gen_3pipi0_3pi_2017->GetEntries() << endl;
    cout << t_gen_3pi3pi2pi0_2017->GetEntries() << endl;
    cout << t_gen_all_2017->GetEntries() << endl;

    t_gen_all_2017->Add(t_gen_3pi3pi_pi0_2017);
    t_gen_all_2017->Add(t_gen_3pipi0_3pi_2017);
    t_gen_all_2017->Add(t_gen_3pi3pi2pi0_2017);

    // 2018
    TChain* t_gen_3pi3pi_2018 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
    TChain* t_gen_3pi3pi_pi0_2018 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
    TChain* t_gen_3pipi0_3pi_2018 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
    TChain* t_gen_3pi3pi2pi0_2018 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
    TChain* t_gen_all_2018 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");

    t_gen_3pi3pi_2018->AddFileInfoList((TCollection*)fc_gen_2018->GetList());
    t_gen_3pi3pi_pi0_2018->AddFileInfoList((TCollection*)fc_gen_2018->GetList());
    t_gen_3pipi0_3pi_2018->AddFileInfoList((TCollection*)fc_gen_2018->GetList());
    t_gen_3pi3pi2pi0_2018->AddFileInfoList((TCollection*)fc_gen_2018->GetList());
    t_gen_all_2018->AddFileInfoList((TCollection*)fc_gen_2018->GetList());

    cout << "Ktautau 2018" << endl;
    cout << t_gen_3pi3pi_2018->GetEntries() << endl;
    cout << t_gen_3pi3pi_pi0_2018->GetEntries() << endl;
    cout << t_gen_3pipi0_3pi_2018->GetEntries() << endl;
    cout << t_gen_3pi3pi2pi0_2018->GetEntries() << endl;
    cout << t_gen_all_2018->GetEntries() << endl;

    t_gen_all_2018->Add(t_gen_3pi3pi_pi0_2018);
    t_gen_all_2018->Add(t_gen_3pipi0_3pi_2018);
    t_gen_all_2018->Add(t_gen_3pi3pi2pi0_2018);

    t_gen_3pi3pi->Add(t_gen_3pi3pi_2017);
    t_gen_3pi3pi->Add(t_gen_3pi3pi_2018);

    t_gen_all->Add(t_gen_all_2017);
    t_gen_all->Add(t_gen_all_2018);

    // WS data after all selections
    TFile* f_ws = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_ws = (TTree*)f_ws->Get("DecayTree");

    // Norm MC after all selections
    TFile* f_norm_mc = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_72/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_norm_mc = (TTree*)f_norm_mc->Get("DecayTree");

    // Gen norm MC
    TFileCollection* fc_gen_norm_2016 = new TFileCollection("fc_gen_norm_2016", "fc_gen_norm_2016", "Files_on_grid/MC_D0Dps_2016.txt");
    TFileCollection* fc_gen_norm_2017 = new TFileCollection("fc_gen_norm_2017", "fc_gen_norm_2017", "Files_on_grid/MC_D0Dps_2017.txt");
    TFileCollection* fc_gen_norm_2018 = new TFileCollection("fc_gen_norm_2018", "fc_gen_norm_2018", "Files_on_grid/MC_D0Dps_2018.txt");

    TChain* t_norm_gen = new TChain("mc_ntuple/MCDecayTree");
    TChain* t_norm_gen_2017 = new TChain("mc_ntuple/MCDecayTree");
    TChain* t_norm_gen_2018 = new TChain("mc_ntuple/MCDecayTree");

    t_norm_gen->AddFileInfoList((TCollection*)fc_gen_norm_2016->GetList());
    t_norm_gen_2017->AddFileInfoList((TCollection*)fc_gen_norm_2017->GetList());
    t_norm_gen_2018->AddFileInfoList((TCollection*)fc_gen_norm_2018->GetList());

    cout << "DDs gen" << endl;
    cout << t_norm_gen->GetEntries() << endl;
    cout << t_norm_gen_2017->GetEntries() << endl;
    cout << t_norm_gen_2018->GetEntries() << endl;

    t_norm_gen->Add(t_norm_gen_2017);
    t_norm_gen->Add(t_norm_gen_2018);

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

            std::vector<Double_t> results = B_Bp_Ktautau(BDT1, BDT2, the_case, t_mc, t_gen_3pi3pi, t_gen_all, t_norm_mc, t_norm_gen, t_ws);
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


std::vector<Double_t> B_Bp_Ktautau(Double_t BDT1, Double_t BDT2, TString the_case, TTree* t_mc, TTree* t_gen_3pi3pi, TTree* t_gen_all, TTree* t_norm_mc, TTree* t_norm_gen, TTree* t_ws)
{
    // External inputs:
    // 1. Branching fractions (from PDG)
    Double_t B_tau_3pi_nu = 9.31/100; // +/- 0.05/100
    Double_t B_tau_3pi_pi0_nu = 4.62/100; // +/- 0.05/100
    Double_t B_Bp_D0bar_Dsp = 9.0/1000; // +/- 0.9/1000
    Double_t B_D0bar_Kp_pi = 3.947/100; // +/- 0.030/100
    Double_t B_Dsp_Kp_Km_pi = 5.37/100; // +/- 0.10/100
    // 2. Signal yield of normalisation channel (from mass fit)
    Double_t N_sig_norm = 32858;
    // 3. Efficiency of normalisation channel
    Double_t N_norm_num = t_norm_mc->GetEntries();
    Double_t N_norm_den = t_norm_gen->GetEntries();
    Double_t eps_sig_norm_post_acc = N_norm_num/N_norm_den;

    Double_t eps_acc_norm_2016 = 14.50/100;
    Double_t eps_acc_norm_2017 = 14.48/100;
    Double_t eps_acc_norm_2018 = 14.58/100;
    Double_t eps_acc_norm = (10.8/34.2)*eps_acc_norm_2016 + (12./34.2)*eps_acc_norm_2017 + (11.4/34.2)*eps_acc_norm_2018;

    Double_t eps_sig_norm = eps_acc_norm*eps_sig_norm_post_acc;
    // cout << "Norm pre-BDT sig eff = " << eps_sig_norm*100 << endl;

    // 4. Ktautau signal pre-selection efficiency
    Double_t eps_acc_2016 = (5.318/100);
    Double_t eps_acc_2017 = (5.321/100);
    Double_t eps_acc_2018 = (5.330/100);
    Double_t eps_acc = (28.5/165)*eps_acc_2016 + (60.0/165)*eps_acc_2017 + (76.5/165)*eps_acc_2018;

    Double_t eps_strip_2016 = (1.028/100);
    Double_t eps_strip_2017 = (1.056/100);
    Double_t eps_strip_2018 = (1.050/100);
    Double_t eps_strip = (28.5/165)*eps_strip_2016 + (60.0/165)*eps_strip_2017 + (76.5/165)*eps_strip_2018;

    Double_t N_3pi3pi_num = t_mc->GetEntries("component==0");
    Double_t N_3pi3pi_den = t_gen_3pi3pi->GetEntries();
    Double_t eps_post_strip_3pi3pi = N_3pi3pi_num/N_3pi3pi_den;

    Double_t N_all_num = t_mc->GetEntries();
    Double_t N_all_den = t_gen_all->GetEntries();
    Double_t eps_post_strip_all = N_all_num/N_all_den;

    Double_t eps_sig_ktautau_3pi3pi = eps_acc*eps_strip*eps_post_strip_3pi3pi;
    Double_t eps_sig_ktautau_all_mc = eps_acc*eps_strip*eps_post_strip_all;

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
        N_bkg_Ktautau = 19006409;
    }

    Double_t bdt_sig_eff_den, bdt_sig_eff_num;
    if(use3pi3piMC)
    {
        bdt_sig_eff_den = t_mc->GetEntries("component==0");
        bdt_sig_eff_num = t_mc->GetEntries(Form("(BDT1 > %f) && (BDT2 > %f) && (component==0)",BDT1,BDT2));
    } 
    else
    {
        bdt_sig_eff_den = t_mc->GetEntries();
        bdt_sig_eff_num = t_mc->GetEntries(Form("(BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2));
    }

    Double_t eps_sig_ktautau_bdt = bdt_sig_eff_num/bdt_sig_eff_den;

    eps_sig_ktautau *= eps_sig_ktautau_bdt;

    Double_t eps_bkg_ktautau_bdt;
    if(combBkg)
    {
        Double_t bdt_bkg_eff_den = t_ws->GetEntries();
        Double_t bdt_bkg_eff_num = t_ws->GetEntries(Form("(BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2));

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