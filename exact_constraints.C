Double_t eps_error(Double_t Num, Double_t Den);
double FWHM(TH1* h);
std::vector<TH1D*> gsl_pull_distributions(TChain* t, TCut pass_GSL);

Int_t dimM = 22;
Int_t dimX = 23;
Int_t dimC = 24;

// Name of the true variable
TString true_name[] = {
    "Bp_TRUEORIGINVERTEX_X",
    "Bp_TRUEORIGINVERTEX_Y",
    "Bp_TRUEORIGINVERTEX_Z",
    "taup_TRUEENDVERTEX_X",
    "taup_TRUEENDVERTEX_Y",
    "taup_TRUEENDVERTEX_Z",
    "(taup_pi1_TRUEP_X + taup_pi2_TRUEP_X + taup_pi3_TRUEP_X)",
    "(taup_pi1_TRUEP_Y + taup_pi2_TRUEP_Y + taup_pi3_TRUEP_Y)",
    "(taup_pi1_TRUEP_Z + taup_pi2_TRUEP_Z + taup_pi3_TRUEP_Z)",
    "(taup_pi1_TRUEP_E + taup_pi2_TRUEP_E + taup_pi3_TRUEP_E)",
    "taum_TRUEENDVERTEX_X",
    "taum_TRUEENDVERTEX_Y",
    "taum_TRUEENDVERTEX_Z",
    "(taum_pi1_TRUEP_X + taum_pi2_TRUEP_X + taum_pi3_TRUEP_X)",
    "(taum_pi1_TRUEP_Y + taum_pi2_TRUEP_Y + taum_pi3_TRUEP_Y)",
    "(taum_pi1_TRUEP_Z + taum_pi2_TRUEP_Z + taum_pi3_TRUEP_Z)",
    "(taum_pi1_TRUEP_E + taum_pi2_TRUEP_E + taum_pi3_TRUEP_E)",
    "Bp_TRUEENDVERTEX_X",
    "Bp_TRUEENDVERTEX_Y",
    "Kp_TRUEP_X",
    "Kp_TRUEP_Y",
    "Kp_TRUEP_Z",
    "Bp_TRUEENDVERTEX_X",
    "Bp_TRUEENDVERTEX_Y",
    "Bp_TRUEENDVERTEX_Z",
    "Bp_TRUEP_X",
    "Bp_TRUEP_Y",
    "Bp_TRUEP_Z",
    "5279.41",
    "taup_TRUEP_X",
    "taup_TRUEP_Y",
    "taup_TRUEP_Z",
    "taup_TRUEP_E",
    "(taup_TRUEP_X - taup_pi1_TRUEP_X - taup_pi2_TRUEP_X - taup_pi3_TRUEP_X)",
    "(taup_TRUEP_Y - taup_pi1_TRUEP_Y - taup_pi2_TRUEP_Y - taup_pi3_TRUEP_Y)",
    "(taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z)",
    "(taup_TRUEP_E - taup_pi1_TRUEP_E - taup_pi2_TRUEP_E - taup_pi3_TRUEP_E)",
    "taum_TRUEP_X",
    "taum_TRUEP_Y",
    "taum_TRUEP_Z",
    "taum_TRUEP_E",
    "(taum_TRUEP_X - taum_pi1_TRUEP_X - taum_pi2_TRUEP_X - taum_pi3_TRUEP_X)",
    "(taum_TRUEP_Y - taum_pi1_TRUEP_Y - taum_pi2_TRUEP_Y - taum_pi3_TRUEP_Y)",
    "(taum_TRUEP_Z - taum_pi1_TRUEP_Z - taum_pi2_TRUEP_Z - taum_pi3_TRUEP_Z)",
    "(taum_TRUEP_E - taum_pi1_TRUEP_E - taum_pi2_TRUEP_E - taum_pi3_TRUEP_E)"
};

// Name of the fitted error variable
TString df_err_name[] = {
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
    "df_Bp_MERR",
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

// name of the fitted variable
TString df_name[] = {
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
    "df_Bp_M",
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

TString label_name[] = {
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
    "MB (MeV)",
    "ptau1x (MeV)",
    "ptau1y (MeV)",
    "ptau1z (MeV)",
    "Etau1 (MeV)",
    "pnu1x (MeV)",
    "pnu1y (MeV)",
    "pnu1z (MeV)",
    "Enu1 (MeV)",
    "ptau2x (MeV)",
    "ptau2y (MeV)",
    "ptau2z (MeV)",
    "Etau2 (MeV)",
    "pnu2x (MeV)",
    "pnu2y (MeV)",
    "pnu2z (MeV)",
    "Enu2 (MeV)"
};

TString title_name[] = {
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
    "MB",
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


void exact_constraints()
{
    gStyle->SetOptStat(0);

    // Ktautau MC
    TFileCollection* fc_mc_2016 = new TFileCollection("fc_mc_2016", "fc_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1/pre_sel_tree.txt");
    TFileCollection* fc_mc_2017 = new TFileCollection("fc_mc_2017", "fc_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1/pre_sel_tree.txt");
    TFileCollection* fc_mc_2018 = new TFileCollection("fc_mc_2018", "fc_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1/pre_sel_tree.txt");

    TFileCollection* fc1_mc_2016 = new TFileCollection("fc1_mc_2016", "fc1_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1/fit_results.txt");
    TFileCollection* fc1_mc_2017 = new TFileCollection("fc1_mc_2017", "fc1_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1/fit_results.txt");
    TFileCollection* fc1_mc_2018 = new TFileCollection("fc1_mc_2018", "fc1_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1/fit_results.txt");

    TChain* t_mc_2016 = new TChain("DecayTree");
    TChain* t_mc_2017 = new TChain("DecayTree");
    TChain* t_mc_2018 = new TChain("DecayTree");

    TChain* t1_mc_2016 = new TChain("DecayTree");
    TChain* t1_mc_2017 = new TChain("DecayTree");
    TChain* t1_mc_2018 = new TChain("DecayTree");

    t_mc_2016->AddFileInfoList((TCollection*)fc_mc_2016->GetList());
    t_mc_2017->AddFileInfoList((TCollection*)fc_mc_2017->GetList());
    t_mc_2018->AddFileInfoList((TCollection*)fc_mc_2018->GetList());

    t1_mc_2016->AddFileInfoList((TCollection*)fc1_mc_2016->GetList());
    t1_mc_2017->AddFileInfoList((TCollection*)fc1_mc_2017->GetList());
    t1_mc_2018->AddFileInfoList((TCollection*)fc1_mc_2018->GetList());

    Int_t N_presel_MC_2016 = t_mc_2016->GetEntries();
    Int_t N_presel_MC_2017 = t_mc_2017->GetEntries();
    Int_t N_presel_MC_2018 = t_mc_2018->GetEntries();

    Int_t N_fit_MC_2016 = t1_mc_2016->GetEntries();
    Int_t N_fit_MC_2017 = t1_mc_2017->GetEntries();
    Int_t N_fit_MC_2018 = t1_mc_2018->GetEntries();

    if( (N_presel_MC_2016 != N_fit_MC_2016) || (N_presel_MC_2017 != N_fit_MC_2017) || (N_presel_MC_2018 != N_fit_MC_2018) )
    {
        cout << "Wrong number of entries" << endl;
        return;
    }
    else
    {
        t_mc_2016->AddFriend(t1_mc_2016);
        t_mc_2017->AddFriend(t1_mc_2017);
        t_mc_2018->AddFriend(t1_mc_2018);
    }

    // RS data
    TFileCollection* fc_rs_2016 = new TFileCollection("fc_rs_2016", "fc_rs_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_2/pre_sel_tree.txt", 100);
    TFileCollection* fc_rs_2017 = new TFileCollection("fc_rs_2017", "fc_rs_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_2/pre_sel_tree.txt", 100);
    TFileCollection* fc_rs_2018 = new TFileCollection("fc_rs_2018", "fc_rs_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_2/pre_sel_tree.txt", 100);

    TFileCollection* fc1_rs_2016 = new TFileCollection("fc1_rs_2016", "fc1_rs_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt", 100);
    TFileCollection* fc1_rs_2017 = new TFileCollection("fc1_rs_2017", "fc1_rs_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt", 100);
    TFileCollection* fc1_rs_2018 = new TFileCollection("fc1_rs_2018", "fc1_rs_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt", 100);

    TChain* t_rs_2016 = new TChain("DecayTree");
    TChain* t_rs_2017 = new TChain("DecayTree");
    TChain* t_rs_2018 = new TChain("DecayTree");

    TChain* t1_rs_2016 = new TChain("DecayTree");
    TChain* t1_rs_2017 = new TChain("DecayTree");
    TChain* t1_rs_2018 = new TChain("DecayTree");

    t_rs_2016->AddFileInfoList((TCollection*)fc_rs_2016->GetList());
    t_rs_2017->AddFileInfoList((TCollection*)fc_rs_2017->GetList());
    t_rs_2018->AddFileInfoList((TCollection*)fc_rs_2018->GetList());

    t1_rs_2016->AddFileInfoList((TCollection*)fc1_rs_2016->GetList());
    t1_rs_2017->AddFileInfoList((TCollection*)fc1_rs_2017->GetList());
    t1_rs_2018->AddFileInfoList((TCollection*)fc1_rs_2018->GetList());

    Int_t N_presel_RS_2016 = t_rs_2016->GetEntries();
    Int_t N_presel_RS_2017 = t_rs_2017->GetEntries();
    Int_t N_presel_RS_2018 = t_rs_2018->GetEntries();

    Int_t N_fit_RS_2016 = t1_rs_2016->GetEntries();
    Int_t N_fit_RS_2017 = t1_rs_2017->GetEntries();
    Int_t N_fit_RS_2018 = t1_rs_2018->GetEntries();

    if( (N_presel_RS_2016 != N_fit_RS_2016) || (N_presel_RS_2017 != N_fit_RS_2017) || (N_presel_RS_2018 != N_fit_RS_2018) )
    {
        cout << "Wrong number of entries" << endl;
        return;
    }
    else
    {
        t_rs_2016->AddFriend(t1_rs_2016);
        t_rs_2017->AddFriend(t1_rs_2017);
        t_rs_2018->AddFriend(t1_rs_2018);
    }

    // WS data
    TFileCollection* fc_ws_2016 = new TFileCollection("fc_ws_2016", "fc_ws_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 100);
    TFileCollection* fc_ws_2017 = new TFileCollection("fc_ws_2017", "fc_ws_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 100);
    TFileCollection* fc_ws_2018 = new TFileCollection("fc_ws_2018", "fc_ws_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 100);

    TFileCollection* fc1_ws_2016 = new TFileCollection("fc1_ws_2016", "fc1_ws_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 100);
    TFileCollection* fc1_ws_2017 = new TFileCollection("fc1_ws_2017", "fc1_ws_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 100);
    TFileCollection* fc1_ws_2018 = new TFileCollection("fc1_ws_2018", "fc1_ws_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 100);

    TChain* t_ws_2016 = new TChain("DecayTree");
    TChain* t_ws_2017 = new TChain("DecayTree");
    TChain* t_ws_2018 = new TChain("DecayTree");

    TChain* t1_ws_2016 = new TChain("DecayTree");
    TChain* t1_ws_2017 = new TChain("DecayTree");
    TChain* t1_ws_2018 = new TChain("DecayTree");

    t_ws_2016->AddFileInfoList((TCollection*)fc_ws_2016->GetList());
    t_ws_2017->AddFileInfoList((TCollection*)fc_ws_2017->GetList());
    t_ws_2018->AddFileInfoList((TCollection*)fc_ws_2018->GetList());

    t1_ws_2016->AddFileInfoList((TCollection*)fc1_ws_2016->GetList());
    t1_ws_2017->AddFileInfoList((TCollection*)fc1_ws_2017->GetList());
    t1_ws_2018->AddFileInfoList((TCollection*)fc1_ws_2018->GetList());

    Int_t N_presel_WS_2016 = t_ws_2016->GetEntries();
    Int_t N_presel_WS_2017 = t_ws_2017->GetEntries();
    Int_t N_presel_WS_2018 = t_ws_2018->GetEntries();

    Int_t N_fit_WS_2016 = t1_ws_2016->GetEntries();
    Int_t N_fit_WS_2017 = t1_ws_2017->GetEntries();
    Int_t N_fit_WS_2018 = t1_ws_2018->GetEntries();

    if( (N_presel_WS_2016 != N_fit_WS_2016) || (N_presel_WS_2017 != N_fit_WS_2017) || (N_presel_WS_2018 != N_fit_WS_2018) )
    {
        cout << "Wrong number of entries" << endl;
        return;
    }
    else
    {
        t_ws_2016->AddFriend(t1_ws_2016);
        t_ws_2017->AddFriend(t1_ws_2017);
        t_ws_2018->AddFriend(t1_ws_2018);
    }

    TCut pass_GSL = "df_status==0";
    TCut pass_DTF = "Bp_dtf_status[0]==0";

    /////////////////////////////////////////////////////////////////////////// 3pi3pi, 3pi3pipi0, 3pi3pi2pi0 mass distributions ///////////////////////////////////////////////////////////////////////////////
    // 2016
    TCanvas c_2016 = new TCanvas();
    c_2016.cd();
    TH1D* h_3pi3pi_2016 = new TH1D("h_3pi3pi_2016", "h_3pi3pi_2016", 100, 4000, 8000);
    TH1D* h_3pi3pipi0_2016 = new TH1D("h_3pi3pipi0_2016", "h_3pi3pipi0_2016", 100, 4000, 8000);
    TH1D* h_3pi3pi2pi0_2016 = new TH1D("h_3pi3pi2pi0_2016", "h_3pi3pi2pi0_2016", 100, 4000, 8000);
    TH1D* h_all_mc_2016 = new TH1D("h_all_mc_2016", "h_all_mc_2016", 100, 4000, 8000);

    t_mc_2016->Draw("df_Bp_M >> h_3pi3pi_2016", pass_GSL+"component==0");
    t_mc_2016->Draw("df_Bp_M >> h_3pi3pipi0_2016", pass_GSL+"component==1");
    t_mc_2016->Draw("df_Bp_M >> h_3pi3pi2pi0_2016", pass_GSL+"component==2");
    t_mc_2016->Draw("df_Bp_M >> h_all_mc_2016", pass_GSL);

    h_3pi3pi_2016->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_2016->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_2016->SetTitle("2016");

    h_3pi3pi_2016->SetLineColor(kBlue);
    h_3pi3pi_2016->SetFillColorAlpha(kBlue, 0.25);

    h_3pi3pipi0_2016->SetLineColor(kCyan+1);
    h_3pi3pipi0_2016->SetFillColorAlpha(kCyan+1, 0.25);

    h_3pi3pi2pi0_2016->SetLineColor(kGreen+1);
    h_3pi3pi2pi0_2016->SetFillColorAlpha(kGreen+1, 0.25);

    h_3pi3pi_2016->DrawNormalized();
    h_3pi3pipi0_2016->DrawNormalized("same");
    h_3pi3pi2pi0_2016->DrawNormalized("same");

    TLegend* leg_2016 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_2016->AddEntry(h_3pi3pi_2016, "3#pi3#pi MC", "lf");
    leg_2016->AddEntry(h_3pi3pipi0_2016, "3#pi3#pi #pi^{0} MC", "lf");
    leg_2016->AddEntry(h_3pi3pi2pi0_2016, "3#pi3#pi 2#pi^{0} MC", "lf");
    leg_2016->Draw("same");
    c_2016.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_components_2016.pdf");

    // 2017
    TCanvas c_2017 = new TCanvas();
    c_2017.cd();
    TH1D* h_3pi3pi_2017 = new TH1D("h_3pi3pi_2017", "h_3pi3pi_2017", 100, 4000, 8000);
    TH1D* h_3pi3pipi0_2017 = new TH1D("h_3pi3pipi0_2017", "h_3pi3pipi0_2017", 100, 4000, 8000);
    TH1D* h_3pi3pi2pi0_2017 = new TH1D("h_3pi3pi2pi0_2017", "h_3pi3pi2pi0_2017", 100, 4000, 8000);
    TH1D* h_all_mc_2017 = new TH1D("h_all_mc_2017", "h_all_mc_2017", 100, 4000, 8000);

    t_mc_2017->Draw("df_Bp_M >> h_3pi3pi_2017", pass_GSL+"component==0");
    t_mc_2017->Draw("df_Bp_M >> h_3pi3pipi0_2017", pass_GSL+"component==1");
    t_mc_2017->Draw("df_Bp_M >> h_3pi3pi2pi0_2017", pass_GSL+"component==2");
    t_mc_2017->Draw("df_Bp_M >> h_all_mc_2017", pass_GSL);

    h_3pi3pi_2017->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_2017->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_2017->SetTitle("2017");

    h_3pi3pi_2017->SetLineColor(kBlue);
    h_3pi3pi_2017->SetFillColorAlpha(kBlue, 0.25);

    h_3pi3pipi0_2017->SetLineColor(kCyan+1);
    h_3pi3pipi0_2017->SetFillColorAlpha(kCyan+1, 0.25);

    h_3pi3pi2pi0_2017->SetLineColor(kGreen+1);
    h_3pi3pi2pi0_2017->SetFillColorAlpha(kGreen+1, 0.25);

    h_3pi3pi_2017->DrawNormalized();
    h_3pi3pipi0_2017->DrawNormalized("same");
    h_3pi3pi2pi0_2017->DrawNormalized("same");

    TLegend* leg_2017 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_2017->AddEntry(h_3pi3pi_2017, "3#pi3#pi MC", "lf");
    leg_2017->AddEntry(h_3pi3pipi0_2017, "3#pi3#pi #pi^{0} MC", "lf");
    leg_2017->AddEntry(h_3pi3pi2pi0_2017, "3#pi3#pi 2#pi^{0} MC", "lf");
    leg_2017->Draw("same");
    c_2017.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_components_2017.pdf");

    // 2018
    TCanvas c_2018 = new TCanvas();
    c_2018.cd();
    TH1D* h_3pi3pi_2018 = new TH1D("h_3pi3pi_2018", "h_3pi3pi_2018", 100, 4000, 8000);
    TH1D* h_3pi3pipi0_2018 = new TH1D("h_3pi3pipi0_2018", "h_3pi3pipi0_2018", 100, 4000, 8000);
    TH1D* h_3pi3pi2pi0_2018 = new TH1D("h_3pi3pi2pi0_2018", "h_3pi3pi2pi0_2018", 100, 4000, 8000);
    TH1D* h_all_mc_2018 = new TH1D("h_all_mc_2018", "h_all_mc_2018", 100, 4000, 8000);

    t_mc_2018->Draw("df_Bp_M >> h_3pi3pi_2018", pass_GSL+"component==0");
    t_mc_2018->Draw("df_Bp_M >> h_3pi3pipi0_2018", pass_GSL+"component==1");
    t_mc_2018->Draw("df_Bp_M >> h_3pi3pi2pi0_2018", pass_GSL+"component==2");
    t_mc_2018->Draw("df_Bp_M >> h_all_mc_2018", pass_GSL);

    h_3pi3pi_2018->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_2018->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_2018->SetTitle("2018");

    h_3pi3pi_2018->SetLineColor(kBlue);
    h_3pi3pi_2018->SetFillColorAlpha(kBlue, 0.25);

    h_3pi3pipi0_2018->SetLineColor(kCyan+1);
    h_3pi3pipi0_2018->SetFillColorAlpha(kCyan+1, 0.25);

    h_3pi3pi2pi0_2018->SetLineColor(kGreen+1);
    h_3pi3pi2pi0_2018->SetFillColorAlpha(kGreen+1, 0.25);

    h_3pi3pi_2018->DrawNormalized();
    h_3pi3pipi0_2018->DrawNormalized("same");
    h_3pi3pi2pi0_2018->DrawNormalized("same");

    TLegend* leg_2018 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_2018->AddEntry(h_3pi3pi_2018, "3#pi3#pi MC", "lf");
    leg_2018->AddEntry(h_3pi3pipi0_2018, "3#pi3#pi #pi^{0} MC", "lf");
    leg_2018->AddEntry(h_3pi3pi2pi0_2018, "3#pi3#pi 2#pi^{0} MC", "lf");
    leg_2018->Draw("same");
    c_2018.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_components_2018.pdf");

    std::ofstream file("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_components_passing_rate.tex");
    file << " \\begin{table}[!htbp]" << std::endl;
    file << " \\centering " << std::endl;
    file << " \\begin{tabular}{|c|c|c|c|c|c|}" << std::endl;
    file << " \\hline" << std::endl;
    file << "  & $3\\pi3\\pi$ MC & $3\\pi3\\pi \\pi^0$ MC & $3\\pi3\\pi 2\\pi^0$ MC & All MC \\\\ " << std::endl;
    file << "\\hline" << std::endl;

    Double_t N_3pi3pi_2016_num = t_mc_2016->GetEntries(pass_GSL+"component==0");
    Double_t N_3pi3pi_2017_num = t_mc_2017->GetEntries(pass_GSL+"component==0");
    Double_t N_3pi3pi_2018_num = t_mc_2018->GetEntries(pass_GSL+"component==0");

    Double_t N_3pi3pipi0_2016_num = t_mc_2016->GetEntries(pass_GSL+"component==1");
    Double_t N_3pi3pipi0_2017_num = t_mc_2017->GetEntries(pass_GSL+"component==1");
    Double_t N_3pi3pipi0_2018_num = t_mc_2018->GetEntries(pass_GSL+"component==1");

    Double_t N_3pi3pi2pi0_2016_num = t_mc_2016->GetEntries(pass_GSL+"component==2");
    Double_t N_3pi3pi2pi0_2017_num = t_mc_2017->GetEntries(pass_GSL+"component==2");
    Double_t N_3pi3pi2pi0_2018_num = t_mc_2018->GetEntries(pass_GSL+"component==2");

    Double_t N_all_mc_2016_num = t_mc_2016->GetEntries(pass_GSL);
    Double_t N_all_mc_2017_num = t_mc_2017->GetEntries(pass_GSL);
    Double_t N_all_mc_2018_num = t_mc_2018->GetEntries(pass_GSL);

    Double_t N_3pi3pi_2016_den = t_mc_2016->GetEntries("component==0");
    Double_t N_3pi3pi_2017_den = t_mc_2017->GetEntries("component==0");
    Double_t N_3pi3pi_2018_den = t_mc_2018->GetEntries("component==0");

    Double_t N_3pi3pipi0_2016_den = t_mc_2016->GetEntries("component==1");
    Double_t N_3pi3pipi0_2017_den = t_mc_2017->GetEntries("component==1");
    Double_t N_3pi3pipi0_2018_den = t_mc_2018->GetEntries("component==1");

    Double_t N_3pi3pi2pi0_2016_den = t_mc_2016->GetEntries("component==2");
    Double_t N_3pi3pi2pi0_2017_den = t_mc_2017->GetEntries("component==2");
    Double_t N_3pi3pi2pi0_2018_den = t_mc_2018->GetEntries("component==2");

    Double_t N_all_mc_2016_den = t_mc_2016->GetEntries();
    Double_t N_all_mc_2017_den = t_mc_2017->GetEntries();
    Double_t N_all_mc_2018_den = t_mc_2018->GetEntries();

    file << Form("2016 & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %1.f$\\,\\% \\\\ ", (N_3pi3pi_2016_num/N_3pi3pi_2016_den)*100, eps_error(N_3pi3pi_2016_num,N_3pi3pi_2016_den)*100, (N_3pi3pipi0_2016_num/N_3pi3pipi0_2016_den)*100, eps_error(N_3pi3pipi0_2016_num,N_3pi3pipi0_2016_den)*100, (N_3pi3pi2pi0_2016_num/N_3pi3pi2pi0_2016_den)*100, eps_error(N_3pi3pi2pi0_2016_num,N_3pi3pi2pi0_2016_den)*100, (N_all_mc_2016_num/N_all_mc_2016_den)*100, eps_error(N_all_mc_2016_num,N_all_mc_2016_den)*100 ) << std::endl;
    file << Form("2017 & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %1.f$\\,\\% \\\\ ", (N_3pi3pi_2017_num/N_3pi3pi_2017_den)*100, eps_error(N_3pi3pi_2017_num,N_3pi3pi_2017_den)*100, (N_3pi3pipi0_2017_num/N_3pi3pipi0_2017_den)*100, eps_error(N_3pi3pipi0_2017_num,N_3pi3pipi0_2017_den)*100, (N_3pi3pi2pi0_2017_num/N_3pi3pi2pi0_2017_den)*100, eps_error(N_3pi3pi2pi0_2017_num,N_3pi3pi2pi0_2017_den)*100, (N_all_mc_2017_num/N_all_mc_2017_den)*100, eps_error(N_all_mc_2017_num,N_all_mc_2017_den)*100 ) << std::endl;
    file << Form("2018 & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %1.f$\\,\\% \\\\ ", (N_3pi3pi_2018_num/N_3pi3pi_2018_den)*100, eps_error(N_3pi3pi_2018_num,N_3pi3pi_2018_den)*100, (N_3pi3pipi0_2018_num/N_3pi3pipi0_2018_den)*100, eps_error(N_3pi3pipi0_2018_num,N_3pi3pipi0_2018_den)*100, (N_3pi3pi2pi0_2018_num/N_3pi3pi2pi0_2018_den)*100, eps_error(N_3pi3pi2pi0_2018_num,N_3pi3pi2pi0_2018_den)*100, (N_all_mc_2018_num/N_all_mc_2018_den)*100, eps_error(N_all_mc_2018_num,N_all_mc_2018_den)*100 ) << std::endl;

    file << "\\hline" << std::endl;
    file << "\\end{tabular}" << std::endl;
    file << "\\end{table}" << std::endl;
    file.close();

    std::ofstream file1("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_components_resolution.tex");
    file1 << " \\begin{table}[!htbp]" << std::endl;
    file1 << " \\centering " << std::endl;
    file1 << " \\begin{tabular}{|c|c|c|c|c|c|}" << std::endl;
    file1 << " \\hline" << std::endl;
    file1 << "  & $3\\pi3\\pi$ MC & $3\\pi3\\pi \\pi^0$ MC & $3\\pi3\\pi 2\\pi^0$ MC & All MC \\\\ " << std::endl;
    file1 << "\\hline" << std::endl;

    Double_t res_3pi3pi_2016 = FWHM(h_3pi3pi_2016)/2.4;
    Double_t res_3pi3pipi0_2016 = FWHM(h_3pi3pipi0_2016)/2.4;
    Double_t res_3pi3pi2pi0_2016 = FWHM(h_3pi3pi2pi0_2016)/2.4;
    Double_t res_all_mc_2016 = FWHM(h_all_mc_2016)/2.4;

    Double_t res_3pi3pi_2017 = FWHM(h_3pi3pi_2017)/2.4;
    Double_t res_3pi3pipi0_2017 = FWHM(h_3pi3pipi0_2017)/2.4;
    Double_t res_3pi3pi2pi0_2017 = FWHM(h_3pi3pi2pi0_2017)/2.4;
    Double_t res_all_mc_2017 = FWHM(h_all_mc_2017)/2.4;

    Double_t res_3pi3pi_2018 = FWHM(h_3pi3pi_2018)/2.4;
    Double_t res_3pi3pipi0_2018 = FWHM(h_3pi3pipi0_2018)/2.4;
    Double_t res_3pi3pi2pi0_2018 = FWHM(h_3pi3pi2pi0_2018)/2.4;
    Double_t res_all_mc_2018 = FWHM(h_all_mc_2018)/2.4;

    file1 << Form("2016 & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV \\\\", res_3pi3pi_2016, res_3pi3pipi0_2016, res_3pi3pi2pi0_2016, res_all_mc_2016);
    file1 << Form("2017 & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV \\\\", res_3pi3pi_2017, res_3pi3pipi0_2017, res_3pi3pi2pi0_2017, res_all_mc_2017);
    file1 << Form("2018 & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV \\\\", res_3pi3pi_2018, res_3pi3pipi0_2018, res_3pi3pi2pi0_2018, res_all_mc_2018);

    file1 << "\\hline" << std::endl;
    file1 << "\\end{tabular}" << std::endl;
    file1 << "\\end{table}" << std::endl;
    file1.close();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // 3pi3pi vs RS vs WS data
    TH1D* h_rs_2016 = new TH1D("h_rs_2016", "h_rs_2016", 100, 4000, 8000);
    TH1D* h_rs_2017 = new TH1D("h_rs_2017", "h_rs_2017", 100, 4000, 8000);
    TH1D* h_rs_2018 = new TH1D("h_rs_2018", "h_rs_2018", 100, 4000, 8000);

    TH1D* h_ws_2016 = new TH1D("h_ws_2016", "h_ws_2016", 100, 4000, 8000);
    TH1D* h_ws_2017 = new TH1D("h_ws_2017", "h_ws_2017", 100, 4000, 8000);
    TH1D* h_ws_2018 = new TH1D("h_ws_2018", "h_ws_2018", 100, 4000, 8000);

    t_rs_2016->Draw("df_Bp_M >> h_rs_2016", pass_GSL);
    t_rs_2017->Draw("df_Bp_M >> h_rs_2017", pass_GSL);
    t_rs_2018->Draw("df_Bp_M >> h_rs_2018", pass_GSL);

    t_ws_2016->Draw("df_Bp_M >> h_ws_2016", pass_GSL);
    t_ws_2017->Draw("df_Bp_M >> h_ws_2017", pass_GSL);
    t_ws_2018->Draw("df_Bp_M >> h_ws_2018", pass_GSL);

    h_rs_2016->SetLineColor(kBlack);
    h_rs_2016->SetFillColorAlpha(kBlack, 0.25);

    h_rs_2017->SetLineColor(kBlack);
    h_rs_2017->SetFillColorAlpha(kBlack, 0.25);

    h_rs_2018->SetLineColor(kBlack);
    h_rs_2018->SetFillColorAlpha(kBlack, 0.25);

    h_ws_2016->SetLineColor(kRed);
    h_ws_2016->SetFillColorAlpha(kRed, 0.25);

    h_ws_2017->SetLineColor(kRed);
    h_ws_2017->SetFillColorAlpha(kRed, 0.25);

    h_ws_2018->SetLineColor(kRed);
    h_ws_2018->SetFillColorAlpha(kRed, 0.25);

    TCanvas c1_2016 = new TCanvas();
    c1_2016.cd();
    h_3pi3pi_2016->DrawNormalized();
    h_rs_2016->DrawNormalized("same");
    h_ws_2016->DrawNormalized("same");

    TLegend* leg1_2016 = new TLegend(0.7,0.7,0.9,0.9);
    leg1_2016->AddEntry(h_3pi3pi_2016, "3#pi3#pi MC", "lf");
    leg1_2016->AddEntry(h_rs_2016, "RS data", "lf");
    leg1_2016->AddEntry(h_ws_2016, "WS data", "lf");
    leg1_2016->Draw("same");
    c1_2016.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_vs_data_2016.pdf");

    TCanvas c1_2017 = new TCanvas();
    c1_2017.cd();
    h_3pi3pi_2017->DrawNormalized();
    h_rs_2017->DrawNormalized("same");
    h_ws_2017->DrawNormalized("same");

    TLegend* leg1_2017 = new TLegend(0.7,0.7,0.9,0.9);
    leg1_2017->AddEntry(h_3pi3pi_2017, "3#pi3#pi MC", "lf");
    leg1_2017->AddEntry(h_rs_2017, "RS data", "lf");
    leg1_2017->AddEntry(h_ws_2017, "WS data", "lf");
    leg1_2017->Draw("same");
    c1_2017.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_vs_data_2017.pdf");

    TCanvas c1_2018 = new TCanvas();
    c1_2018.cd();
    h_3pi3pi_2018->DrawNormalized();
    h_rs_2018->DrawNormalized("same");
    h_ws_2018->DrawNormalized("same");

    TLegend* leg1_2018 = new TLegend(0.7,0.7,0.9,0.9);
    leg1_2018->AddEntry(h_3pi3pi_2018, "3#pi3#pi MC", "lf");
    leg1_2018->AddEntry(h_rs_2018, "RS data", "lf");
    leg1_2018->AddEntry(h_ws_2018, "WS data", "lf");
    leg1_2018->Draw("same");
    c1_2018.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_vs_data_2018.pdf");

    std::ofstream file2("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_vs_data_passing_rate.tex");
    file2 << " \\begin{table}[!htbp]" << std::endl;
    file2 << " \\centering " << std::endl;
    file2 << " \\begin{tabular}{|c|c|c|c|}" << std::endl;
    file2 << " \\hline" << std::endl;
    file2 << "  & $3\\pi3\\pi$ MC & RS data & WS data \\\\ " << std::endl;
    file2 << "\\hline" << std::endl;    

    Double_t N_rs_2016_num = t_rs_2016->GetEntries(pass_GSL);
    Double_t N_rs_2017_num = t_rs_2017->GetEntries(pass_GSL);
    Double_t N_rs_2018_num = t_rs_2018->GetEntries(pass_GSL);

    Double_t N_ws_2016_num = t_ws_2016->GetEntries(pass_GSL);
    Double_t N_ws_2017_num = t_ws_2017->GetEntries(pass_GSL);
    Double_t N_ws_2018_num = t_ws_2018->GetEntries(pass_GSL);

    Double_t N_rs_2016_den = t_rs_2016->GetEntries();
    Double_t N_rs_2017_den = t_rs_2017->GetEntries();
    Double_t N_rs_2018_den = t_rs_2018->GetEntries();

    Double_t N_ws_2016_den = t_ws_2016->GetEntries();
    Double_t N_ws_2017_den = t_ws_2017->GetEntries();
    Double_t N_ws_2018_den = t_ws_2018->GetEntries();

    file2 << Form("2016 & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% \\\\", (N_3pi3pi_2016_num/N_3pi3pi_2016_den)*100, eps_error(N_3pi3pi_2016_num,N_3pi3pi_2016_den)*100, (N_rs_2016_num/N_rs_2016_den)*100, eps_error(N_rs_2016_num,N_rs_2016_den)*100, (N_ws_2016_num/N_ws_2016_den)*100, eps_error(N_ws_2016_num,N_ws_2016_den)*100 );
    file2 << Form("2017 & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% \\\\", (N_3pi3pi_2017_num/N_3pi3pi_2017_den)*100, eps_error(N_3pi3pi_2017_num,N_3pi3pi_2017_den)*100, (N_rs_2017_num/N_rs_2017_den)*100, eps_error(N_rs_2017_num,N_rs_2017_den)*100, (N_ws_2017_num/N_ws_2017_den)*100, eps_error(N_ws_2017_num,N_ws_2017_den)*100 );
    file2 << Form("2018 & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% \\\\", (N_3pi3pi_2018_num/N_3pi3pi_2018_den)*100, eps_error(N_3pi3pi_2018_num,N_3pi3pi_2018_den)*100, (N_rs_2018_num/N_rs_2018_den)*100, eps_error(N_rs_2018_num,N_rs_2018_den)*100, (N_ws_2018_num/N_ws_2018_den)*100, eps_error(N_ws_2018_num,N_ws_2018_den)*100 );

    file2 << "\\hline" << std::endl;
    file2 << "\\end{tabular}" << std::endl;
    file2 << "\\end{table}" << std::endl;
    file2.close();

    TH1D* h_chi2_3pi3pi_2016 = new TH1D("h_chi2_3pi3pi_2016", "h_chi2_3pi3pi_2016", 100, 0, 30);
    TH1D* h_chi2_3pi3pi_2017 = new TH1D("h_chi2_3pi3pi_2017", "h_chi2_3pi3pi_2017", 100, 0, 30);
    TH1D* h_chi2_3pi3pi_2018 = new TH1D("h_chi2_3pi3pi_2018", "h_chi2_3pi3pi_2018", 100, 0, 30);

    t_mc_2016->Draw("df_chi2 >> h_chi2_3pi3pi_2016", pass_GSL+"component==0");
    t_mc_2017->Draw("df_chi2 >> h_chi2_3pi3pi_2017", pass_GSL+"component==0");
    t_mc_2018->Draw("df_chi2 >> h_chi2_3pi3pi_2018", pass_GSL+"component==0");

    TH1D* h_chi2_rs_2016 = new TH1D("h_chi2_rs_2016", "h_chi2_rs_2016", 100, 0, 30);
    TH1D* h_chi2_rs_2017 = new TH1D("h_chi2_rs_2017", "h_chi2_rs_2017", 100, 0, 30);
    TH1D* h_chi2_rs_2018 = new TH1D("h_chi2_rs_2018", "h_chi2_rs_2018", 100, 0, 30);

    t_rs_2016->Draw("df_chi2 >> h_chi2_rs_2016", pass_GSL);
    t_rs_2017->Draw("df_chi2 >> h_chi2_rs_2017", pass_GSL);
    t_rs_2018->Draw("df_chi2 >> h_chi2_rs_2018", pass_GSL);

    TH1D* h_chi2_ws_2016 = new TH1D("h_chi2_ws_2016", "h_chi2_ws_2016", 100, 0, 30);
    TH1D* h_chi2_ws_2017 = new TH1D("h_chi2_ws_2017", "h_chi2_ws_2017", 100, 0, 30);
    TH1D* h_chi2_ws_2018 = new TH1D("h_chi2_ws_2018", "h_chi2_ws_2018", 100, 0, 30);

    t_ws_2016->Draw("df_chi2 >> h_chi2_ws_2016", pass_GSL);
    t_ws_2017->Draw("df_chi2 >> h_chi2_ws_2017", pass_GSL);
    t_ws_2018->Draw("df_chi2 >> h_chi2_ws_2018", pass_GSL);

    h_chi2_3pi3pi_2016->GetXaxis()->SetTitle("GSL #chi^{2}");
    h_chi2_3pi3pi_2016->GetYaxis()->SetTitle("Normalized entries / (0.3)");
    h_chi2_3pi3pi_2016->SetTitle("2016");

    h_chi2_3pi3pi_2017->GetXaxis()->SetTitle("GSL #chi^{2}");
    h_chi2_3pi3pi_2017->GetYaxis()->SetTitle("Normalized entries / (0.3)");
    h_chi2_3pi3pi_2017->SetTitle("2017");

    h_chi2_3pi3pi_2018->GetXaxis()->SetTitle("GSL #chi^{2}");
    h_chi2_3pi3pi_2018->GetYaxis()->SetTitle("Normalized entries / (0.3)");
    h_chi2_3pi3pi_2018->SetTitle("2018");

    h_chi2_3pi3pi_2016->SetLineColor(kBlue);
    h_chi2_3pi3pi_2016->SetFillColorAlpha(kBlue, 0.25);

    h_chi2_3pi3pi_2017->SetLineColor(kBlue);
    h_chi2_3pi3pi_2017->SetFillColorAlpha(kBlue, 0.25);

    h_chi2_3pi3pi_2018->SetLineColor(kBlue);
    h_chi2_3pi3pi_2018->SetFillColorAlpha(kBlue, 0.25);

    h_chi2_rs_2016->SetLineColor(kBlack);
    h_chi2_rs_2016->SetFillColorAlpha(kBlack, 0.25);

    h_chi2_rs_2017->SetLineColor(kBlack);
    h_chi2_rs_2017->SetFillColorAlpha(kBlack, 0.25);

    h_chi2_rs_2018->SetLineColor(kBlack);
    h_chi2_rs_2018->SetFillColorAlpha(kBlack, 0.25);

    h_chi2_ws_2016->SetLineColor(kRed);
    h_chi2_ws_2016->SetFillColorAlpha(kRed, 0.25);

    h_chi2_ws_2017->SetLineColor(kRed);
    h_chi2_ws_2017->SetFillColorAlpha(kRed, 0.25);

    h_chi2_ws_2018->SetLineColor(kRed);
    h_chi2_ws_2018->SetFillColorAlpha(kRed, 0.25);

    TCanvas c2_2016 = new TCanvas();
    c2_2016.cd();
    h_chi2_3pi3pi_2016->DrawNormalized();
    h_chi2_rs_2016->DrawNormalized("same");
    h_chi2_ws_2016->DrawNormalized("same");

    TLegend* leg2_2016 = new TLegend(0.7,0.7,0.9,0.9);
    leg2_2016->AddEntry(h_chi2_3pi3pi_2016, "3#pi3#pi MC", "lf");
    leg2_2016->AddEntry(h_chi2_rs_2016, "RS data", "lf");
    leg2_2016->AddEntry(h_chi2_ws_2016, "WS data", "lf");
    leg2_2016->Draw("same");
    c2_2016.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/chi2_mc_vs_data_2016.pdf");

    TCanvas c2_2017 = new TCanvas();
    c2_2017.cd();
    h_chi2_3pi3pi_2017->DrawNormalized();
    h_chi2_rs_2017->DrawNormalized("same");
    h_chi2_ws_2017->DrawNormalized("same");

    TLegend* leg2_2017 = new TLegend(0.7,0.7,0.9,0.9);
    leg2_2017->AddEntry(h_chi2_3pi3pi_2017, "3#pi3#pi MC", "lf");
    leg2_2017->AddEntry(h_chi2_rs_2017, "RS data", "lf");
    leg2_2017->AddEntry(h_chi2_ws_2017, "WS data", "lf");
    leg2_2017->Draw("same");
    c2_2017.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/chi2_mc_vs_data_2017.pdf");

    TCanvas c2_2018 = new TCanvas();
    c2_2018.cd();
    h_chi2_3pi3pi_2018->DrawNormalized();
    h_chi2_rs_2018->DrawNormalized("same");
    h_chi2_ws_2018->DrawNormalized("same");

    TLegend* leg2_2018 = new TLegend(0.7,0.7,0.9,0.9);
    leg2_2018->AddEntry(h_chi2_3pi3pi_2018, "3#pi3#pi MC", "lf");
    leg2_2018->AddEntry(h_chi2_rs_2018, "RS data", "lf");
    leg2_2018->AddEntry(h_chi2_ws_2018, "WS data", "lf");
    leg2_2018->Draw("same");
    c2_2018.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/chi2_mc_vs_data_2018.pdf");

    TH1D* h_init_3pi3pi_2016 = new TH1D("h_init_3pi3pi_2016", "h_init_3pi3pi_2016", 5, 0, 5);
    TH1D* h_init_3pi3pi_2017 = new TH1D("h_init_3pi3pi_2017", "h_init_3pi3pi_2017", 5, 0, 5);
    TH1D* h_init_3pi3pi_2018 = new TH1D("h_init_3pi3pi_2018", "h_init_3pi3pi_2018", 5, 0, 5);

    t_mc_2016->Draw("df_init >> h_init_3pi3pi_2016", pass_GSL+"component==0");
    t_mc_2017->Draw("df_init >> h_init_3pi3pi_2017", pass_GSL+"component==0");
    t_mc_2018->Draw("df_init >> h_init_3pi3pi_2018", pass_GSL+"component==0");

    TH1D* h_init_rs_2016 = new TH1D("h_init_rs_2016", "h_init_rs_2016", 5, 0, 5);
    TH1D* h_init_rs_2017 = new TH1D("h_init_rs_2017", "h_init_rs_2017", 5, 0, 5);
    TH1D* h_init_rs_2018 = new TH1D("h_init_rs_2018", "h_init_rs_2018", 5, 0, 5);

    t_rs_2016->Draw("df_init >> h_init_rs_2016", pass_GSL);
    t_rs_2017->Draw("df_init >> h_init_rs_2017", pass_GSL);
    t_rs_2018->Draw("df_init >> h_init_rs_2018", pass_GSL);

    TH1D* h_init_ws_2016 = new TH1D("h_init_ws_2016", "h_init_ws_2016", 5, 0, 5);
    TH1D* h_init_ws_2017 = new TH1D("h_init_ws_2017", "h_init_ws_2017", 5, 0, 5);
    TH1D* h_init_ws_2018 = new TH1D("h_init_ws_2018", "h_init_ws_2018", 5, 0, 5);

    t_ws_2016->Draw("df_init >> h_init_ws_2016", pass_GSL);
    t_ws_2017->Draw("df_init >> h_init_ws_2017", pass_GSL);
    t_ws_2018->Draw("df_init >> h_init_ws_2018", pass_GSL);

    h_init_rs_2016->GetXaxis()->SetTitle("Initialisation picked");
    h_init_rs_2016->GetYaxis()->SetTitle("Normalized entries / (1)");
    h_init_rs_2016->SetTitle("2016");

    h_init_rs_2017->GetXaxis()->SetTitle("Initialisation picked");
    h_init_rs_2017->GetYaxis()->SetTitle("Normalized entries / (1)");
    h_init_rs_2017->SetTitle("2017");

    h_init_rs_2018->GetXaxis()->SetTitle("Initialisation picked");
    h_init_rs_2018->GetYaxis()->SetTitle("Normalized entries / (1)");
    h_init_rs_2018->SetTitle("2018");

    h_init_3pi3pi_2016->SetLineColor(kBlue);
    h_init_3pi3pi_2016->SetFillColorAlpha(kBlue, 0.25);

    h_init_3pi3pi_2017->SetLineColor(kBlue);
    h_init_3pi3pi_2017->SetFillColorAlpha(kBlue, 0.25);

    h_init_3pi3pi_2018->SetLineColor(kBlue);
    h_init_3pi3pi_2018->SetFillColorAlpha(kBlue, 0.25);

    h_init_rs_2016->SetLineColor(kBlack);
    h_init_rs_2016->SetFillColorAlpha(kBlack, 0.25);

    h_init_rs_2017->SetLineColor(kBlack);
    h_init_rs_2017->SetFillColorAlpha(kBlack, 0.25);

    h_init_rs_2018->SetLineColor(kBlack);
    h_init_rs_2018->SetFillColorAlpha(kBlack, 0.25);

    h_init_ws_2016->SetLineColor(kRed);
    h_init_ws_2016->SetFillColorAlpha(kRed, 0.25);

    h_init_ws_2017->SetLineColor(kRed);
    h_init_ws_2017->SetFillColorAlpha(kRed, 0.25);

    h_init_ws_2018->SetLineColor(kRed);
    h_init_ws_2018->SetFillColorAlpha(kRed, 0.25);

    TCanvas c3_2016 = new TCanvas();
    c3_2016.cd();
    h_init_rs_2016->DrawNormalized();
    h_init_3pi3pi_2016->DrawNormalized("same");
    h_init_ws_2016->DrawNormalized("same");

    TLegend* leg3_2016 = new TLegend(0.7,0.7,0.9,0.9);
    leg3_2016->AddEntry(h_init_3pi3pi_2016, "3#pi3#pi MC", "lf");
    leg3_2016->AddEntry(h_init_rs_2016, "RS data", "lf");
    leg3_2016->AddEntry(h_init_ws_2016, "WS data", "lf");
    leg3_2016->Draw("same");
    c3_2016.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/init_mc_vs_data_2016.pdf");

    TCanvas c3_2017 = new TCanvas();
    c3_2017.cd();
    h_init_rs_2017->DrawNormalized();
    h_init_3pi3pi_2017->DrawNormalized("same");
    h_init_ws_2017->DrawNormalized("same");

    TLegend* leg3_2017 = new TLegend(0.7,0.7,0.9,0.9);
    leg3_2017->AddEntry(h_init_3pi3pi_2017, "3#pi3#pi MC", "lf");
    leg3_2017->AddEntry(h_init_rs_2017, "RS data", "lf");
    leg3_2017->AddEntry(h_init_ws_2017, "WS data", "lf");
    leg3_2017->Draw("same");
    c3_2017.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/init_mc_vs_data_2017.pdf");

    TCanvas c3_2018 = new TCanvas();
    c3_2018.cd();
    h_init_rs_2018->DrawNormalized();
    h_init_3pi3pi_2018->DrawNormalized("same");
    h_init_ws_2018->DrawNormalized("same");

    TLegend* leg3_2018 = new TLegend(0.7,0.7,0.9,0.9);
    leg3_2018->AddEntry(h_init_3pi3pi_2018, "3#pi3#pi MC", "lf");
    leg3_2018->AddEntry(h_init_rs_2018, "RS data", "lf");
    leg3_2018->AddEntry(h_init_ws_2018, "WS data", "lf");
    leg3_2018->Draw("same");
    c3_2018.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/init_mc_vs_data_2018.pdf");

    TH1D* h_min_3pi3pi_2016 = new TH1D("h_min_3pi3pi_2016", "h_min_3pi3pi_2016", 5, 0, 5);
    TH1D* h_min_3pi3pi_2017 = new TH1D("h_min_3pi3pi_2017", "h_min_3pi3pi_2017", 5, 0, 5);
    TH1D* h_min_3pi3pi_2018 = new TH1D("h_min_3pi3pi_2018", "h_min_3pi3pi_2018", 5, 0, 5);

    t_mc_2016->Draw("df_N_local_min >> h_min_3pi3pi_2016", pass_GSL+"component==0");
    t_mc_2017->Draw("df_N_local_min >> h_min_3pi3pi_2017", pass_GSL+"component==0");
    t_mc_2018->Draw("df_N_local_min >> h_min_3pi3pi_2018", pass_GSL+"component==0");

    TH1D* h_min_rs_2016 = new TH1D("h_min_rs_2016", "h_min_rs_2016", 5, 0, 5);
    TH1D* h_min_rs_2017 = new TH1D("h_min_rs_2017", "h_min_rs_2017", 5, 0, 5);
    TH1D* h_min_rs_2018 = new TH1D("h_min_rs_2018", "h_min_rs_2018", 5, 0, 5);

    t_rs_2016->Draw("df_N_local_min >> h_min_rs_2016", pass_GSL);
    t_rs_2017->Draw("df_N_local_min >> h_min_rs_2017", pass_GSL);
    t_rs_2018->Draw("df_N_local_min >> h_min_rs_2018", pass_GSL);

    TH1D* h_min_ws_2016 = new TH1D("h_min_ws_2016", "h_min_ws_2016", 5, 0, 5);
    TH1D* h_min_ws_2017 = new TH1D("h_min_ws_2017", "h_min_ws_2017", 5, 0, 5);
    TH1D* h_min_ws_2018 = new TH1D("h_min_ws_2018", "h_min_ws_2018", 5, 0, 5);

    t_ws_2016->Draw("df_N_local_min >> h_min_ws_2016", pass_GSL);
    t_ws_2017->Draw("df_N_local_min >> h_min_ws_2017", pass_GSL);
    t_ws_2018->Draw("df_N_local_min >> h_min_ws_2018", pass_GSL);

    h_min_rs_2016->GetXaxis()->SetTitle("Number of local minima");
    h_min_rs_2016->GetYaxis()->SetTitle("Normalized entries / (1)");
    h_min_rs_2016->SetTitle("2016");

    h_min_rs_2017->GetXaxis()->SetTitle("Number of local minima");
    h_min_rs_2017->GetYaxis()->SetTitle("Normalized entries / (1)");
    h_min_rs_2017->SetTitle("2017");

    h_min_rs_2018->GetXaxis()->SetTitle("Number of local minima");
    h_min_rs_2018->GetYaxis()->SetTitle("Normalized entries / (1)");
    h_min_rs_2018->SetTitle("2018");

    h_min_3pi3pi_2016->SetLineColor(kBlue);
    h_min_3pi3pi_2016->SetFillColorAlpha(kBlue, 0.25);

    h_min_3pi3pi_2017->SetLineColor(kBlue);
    h_min_3pi3pi_2017->SetFillColorAlpha(kBlue, 0.25);

    h_min_3pi3pi_2018->SetLineColor(kBlue);
    h_min_3pi3pi_2018->SetFillColorAlpha(kBlue, 0.25);

    h_min_rs_2016->SetLineColor(kBlack);
    h_min_rs_2016->SetFillColorAlpha(kBlack, 0.25);

    h_min_rs_2017->SetLineColor(kBlack);
    h_min_rs_2017->SetFillColorAlpha(kBlack, 0.25);

    h_min_rs_2018->SetLineColor(kBlack);
    h_min_rs_2018->SetFillColorAlpha(kBlack, 0.25);

    h_min_ws_2016->SetLineColor(kRed);
    h_min_ws_2016->SetFillColorAlpha(kRed, 0.25);

    h_min_ws_2017->SetLineColor(kRed);
    h_min_ws_2017->SetFillColorAlpha(kRed, 0.25);

    h_min_ws_2018->SetLineColor(kRed);
    h_min_ws_2018->SetFillColorAlpha(kRed, 0.25);

    TCanvas c4_2016 = new TCanvas();
    c4_2016.cd();
    h_min_rs_2016->DrawNormalized();
    h_min_3pi3pi_2016->DrawNormalized("same");
    h_min_ws_2016->DrawNormalized("same");

    TLegend* leg4_2016 = new TLegend(0.7,0.7,0.9,0.9);
    leg4_2016->AddEntry(h_min_3pi3pi_2016, "3#pi3#pi MC", "lf");
    leg4_2016->AddEntry(h_min_rs_2016, "RS data", "lf");
    leg4_2016->AddEntry(h_min_ws_2016, "WS data", "lf");
    leg4_2016->Draw("same");
    c4_2016.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/Nlocalmin_mc_vs_data_2016.pdf");

    TCanvas c4_2017 = new TCanvas();
    c4_2017.cd();
    h_min_rs_2017->DrawNormalized();
    h_min_3pi3pi_2017->DrawNormalized("same");
    h_min_ws_2017->DrawNormalized("same");

    TLegend* leg4_2017 = new TLegend(0.7,0.7,0.9,0.9);
    leg4_2017->AddEntry(h_min_3pi3pi_2017, "3#pi3#pi MC", "lf");
    leg4_2017->AddEntry(h_min_rs_2017, "RS data", "lf");
    leg4_2017->AddEntry(h_min_ws_2017, "WS data", "lf");
    leg4_2017->Draw("same");
    c4_2017.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/Nlocalmin_mc_vs_data_2017.pdf");

    TCanvas c4_2018 = new TCanvas();
    c4_2018.cd();
    h_min_rs_2018->DrawNormalized();
    h_min_3pi3pi_2018->DrawNormalized("same");
    h_min_ws_2018->DrawNormalized("same");

    TLegend* leg4_2018 = new TLegend(0.7,0.7,0.9,0.9);
    leg4_2018->AddEntry(h_min_3pi3pi_2018, "3#pi3#pi MC", "lf");
    leg4_2018->AddEntry(h_min_rs_2018, "RS data", "lf");
    leg4_2018->AddEntry(h_min_ws_2018, "WS data", "lf");
    leg4_2018->Draw("same");
    c4_2018.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/Nlocalmin_mc_vs_data_2018.pdf");

    ///////////////////////////////////////////////////////////////////////////////////////////////// B+ mass in error categories /////////////////////////////////////////////////////////////////////////////////////////////////////////
    TChain* t_mc = new TChain("DecayTree");
    t_mc->AddFileInfoList((TCollection*)fc_mc_2016->GetList());
    t_mc->Add(t_mc_2017);
    t_mc->Add(t_mc_2018);

    TChain* t1_mc = new TChain("DecayTree");
    t1_mc->AddFileInfoList((TCollection*)fc1_mc_2016->GetList());
    t1_mc->Add(t1_mc_2017);
    t1_mc->Add(t1_mc_2018);
    t_mc->AddFriend(t1_mc);

    TChain* t_rs = new TChain("DecayTree");
    t_rs->AddFileInfoList((TCollection*)fc_rs_2016->GetList());
    t_rs->Add(t_rs_2017);
    t_rs->Add(t_rs_2018);

    TChain* t1_rs = new TChain("DecayTree");
    t1_rs->AddFileInfoList((TCollection*)fc1_rs_2016->GetList());
    t1_rs->GetEntries();
    t1_rs_2017->GetEntries();
    t1_rs_2018->GetEntries();
    t1_rs->Add(t1_rs_2017);
    t1_rs->Add(t1_rs_2018);
    t_rs->AddFriend(t1_rs);

    TChain* t_ws = new TChain("DecayTree");
    t_ws->AddFileInfoList((TCollection*)fc_ws_2016->GetList());
    t_ws->Add(t_ws_2017);
    t_ws->Add(t_ws_2018);

    TChain* t1_ws = new TChain("DecayTree");
    t1_ws->AddFileInfoList((TCollection*)fc1_ws_2016->GetList());    
    t1_ws->Add(t1_ws_2017);
    t1_ws->Add(t1_ws_2018);
    t_ws->AddFriend(t1_ws);

    TH1D* h_merr_mc = new TH1D("h_merr_mc", "h_merr_mc", 100, 0, 1000);
    TH1D* h_merr_rs = new TH1D("h_merr_rs", "h_merr_rs", 100, 0, 1000);
    TH1D* h_merr_ws = new TH1D("h_merr_ws", "h_merr_ws", 100, 0, 1000);

    t_mc->Draw("df_Bp_MERR >> h_merr_mc", pass_GSL+"component==0");
    t_rs->Draw("df_Bp_MERR >> h_merr_rs", pass_GSL);
    t_ws->Draw("df_Bp_MERR >> h_merr_ws", pass_GSL);

    h_merr_ws->GetXaxis()->SetTitle("#Delta m_{B} (MeV)");
    h_merr_ws->GetYaxis()->SetTitle("Normalized entries / (10 MeV)");
    h_merr_ws->SetTitle("2016-2018");

    h_merr_mc->SetLineColor(kBlue);
    h_merr_mc->SetFillColorAlpha(kBlue, 0.25);

    h_merr_rs->SetLineColor(kBlack);
    h_merr_rs->SetFillColorAlpha(kBlack, 0.25);

    h_merr_ws->SetLineColor(kRed);
    h_merr_ws->SetFillColorAlpha(kRed, 0.25);

    TCanvas c5 = new TCanvas();
    c5.cd();
    h_merr_ws->DrawNormalized();
    h_merr_mc->DrawNormalized("same");
    h_merr_rs->DrawNormalized("same");

    TLegend* leg5 = new TLegend(0.7,0.7,0.9,0.9);
    leg5->AddEntry(h_merr_mc, "3#pi3#pi MC", "lf");
    leg5->AddEntry(h_merr_rs, "RS data", "lf");
    leg5->AddEntry(h_merr_ws, "WS data", "lf");
    leg5->Draw("same");
    c5.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_mc_vs_data.pdf");

    TCut error0 = "(df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)";
    TCut error1 = "(df_Bp_MERR > 100) && (df_Bp_MERR <= 250)";
    TCut error2 = "(df_Bp_MERR > 250)";

    TH1D* h_3pi3pi_0 = new TH1D("h_3pi3pi_0", "h_3pi3pi_0", 100, 4000, 8000);
    TH1D* h_3pi3pi_1 = new TH1D("h_3pi3pi_1", "h_3pi3pi_1", 100, 4000, 8000);
    TH1D* h_3pi3pi_2 = new TH1D("h_3pi3pi_2", "h_3pi3pi_2", 100, 4000, 8000);
    TH1D* h_3pi3pi = new TH1D("h_3pi3pi", "h_3pi3pi", 100, 4000, 8000);

    t_mc->Draw("df_Bp_M >> h_3pi3pi_0", pass_GSL+"(component==0)"+error0);
    t_mc->Draw("df_Bp_M >> h_3pi3pi_1", pass_GSL+"(component==0)"+error1);
    t_mc->Draw("df_Bp_M >> h_3pi3pi_2", pass_GSL+"(component==0)"+error2);
    t_mc->Draw("df_Bp_M >> h_3pi3pi", pass_GSL+"(component==0)");

    TH1D* h_3pi3pipi0_0 = new TH1D("h_3pi3pipi0_0", "h_3pi3pipi0_0", 100, 4000, 8000);
    TH1D* h_3pi3pipi0_1 = new TH1D("h_3pi3pipi0_1", "h_3pi3pipi0_1", 100, 4000, 8000);
    TH1D* h_3pi3pipi0_2 = new TH1D("h_3pi3pipi0_2", "h_3pi3pipi0_2", 100, 4000, 8000);

    t_mc->Draw("df_Bp_M >> h_3pi3pipi0_0", pass_GSL+"(component==1)"+error0);
    t_mc->Draw("df_Bp_M >> h_3pi3pipi0_1", pass_GSL+"(component==1)"+error1);
    t_mc->Draw("df_Bp_M >> h_3pi3pipi0_2", pass_GSL+"(component==1)"+error2);

    TH1D* h_3pi3pi2pi0_0 = new TH1D("h_3pi3pi2pi0_0", "h_3pi3pi2pi0_0", 100, 4000, 8000);
    TH1D* h_3pi3pi2pi0_1 = new TH1D("h_3pi3pi2pi0_1", "h_3pi3pi2pi0_1", 100, 4000, 8000);
    TH1D* h_3pi3pi2pi0_2 = new TH1D("h_3pi3pi2pi0_2", "h_3pi3pi2pi0_2", 100, 4000, 8000);

    t_mc->Draw("df_Bp_M >> h_3pi3pi2pi0_0", pass_GSL+"(component==2)"+error0);
    t_mc->Draw("df_Bp_M >> h_3pi3pi2pi0_1", pass_GSL+"(component==2)"+error1);
    t_mc->Draw("df_Bp_M >> h_3pi3pi2pi0_2", pass_GSL+"(component==2)"+error2);

    TH1D* h_all_mc_0 = new TH1D("h_all_mc_0", "h_all_mc_0", 100, 4000, 8000);
    TH1D* h_all_mc_1 = new TH1D("h_all_mc_1", "h_all_mc_1", 100, 4000, 8000);
    TH1D* h_all_mc_2 = new TH1D("h_all_mc_2", "h_all_mc_2", 100, 4000, 8000);

    t_mc->Draw("df_Bp_M >> h_all_mc_0", pass_GSL+error0);
    t_mc->Draw("df_Bp_M >> h_all_mc_1", pass_GSL+error1);
    t_mc->Draw("df_Bp_M >> h_all_mc_2", pass_GSL+error2);

    TH1D* h_rs_0 = new TH1D("h_rs_0", "h_rs_0", 100, 4000, 8000);
    TH1D* h_rs_1 = new TH1D("h_rs_1", "h_rs_1", 100, 4000, 8000);
    TH1D* h_rs_2 = new TH1D("h_rs_2", "h_rs_2", 100, 4000, 8000);
    TH1D* h_rs = new TH1D("h_rs", "h_rs", 100, 4000, 8000);

    t_rs->Draw("df_Bp_M >> h_rs_0", pass_GSL+error0);
    t_rs->Draw("df_Bp_M >> h_rs_1", pass_GSL+error1);
    t_rs->Draw("df_Bp_M >> h_rs_2", pass_GSL+error2);
    t_rs->Draw("df_Bp_M >> h_rs", pass_GSL);

    TH1D* h_ws_0 = new TH1D("h_ws_0", "h_ws_0", 100, 4000, 8000);
    TH1D* h_ws_1 = new TH1D("h_ws_1", "h_ws_1", 100, 4000, 8000);
    TH1D* h_ws_2 = new TH1D("h_ws_2", "h_ws_2", 100, 4000, 8000);
    TH1D* h_ws = new TH1D("h_ws", "h_ws", 100, 4000, 8000);

    t_ws->Draw("df_Bp_M >> h_ws_0", pass_GSL+error0);
    t_ws->Draw("df_Bp_M >> h_ws_1", pass_GSL+error1);
    t_ws->Draw("df_Bp_M >> h_ws_2", pass_GSL+error2);
    t_ws->Draw("df_Bp_M >> h_ws", pass_GSL);
    
    h_3pi3pi_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_0->SetTitle("3#pi3#pi MC: 2016-2018");

    h_3pi3pipi0_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pipi0_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pipi0_0->SetTitle("3#pi3#pi #pi^{0} MC: 2016-2018");

    h_3pi3pi2pi0_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi2pi0_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi2pi0_0->SetTitle("3#pi3#pi 2#pi^{0} MC: 2016-2018");

    h_all_mc_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_all_mc_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_all_mc_0->SetTitle("All MC: 2016-2018");

    h_rs_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_rs_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_rs_0->SetTitle("RS data: 2016-2018");

    h_ws_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ws_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ws_0->SetTitle("WS data: 2016-2018");

    TCanvas c6 = new TCanvas();
    c6.cd();

    h_3pi3pi_0->SetLineColor(kBlue);
    h_3pi3pi_0->SetFillColorAlpha(kBlue, 0.25);

    h_3pi3pi_1->SetLineColor(kGreen+1);
    h_3pi3pi_1->SetFillColorAlpha(kGreen+1, 0.25);

    h_3pi3pi_2->SetLineColor(kRed);
    h_3pi3pi_2->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pi_0->DrawNormalized();
    h_3pi3pi_1->DrawNormalized("same");
    h_3pi3pi_2->DrawNormalized("same");

    TLegend* leg6 = new TLegend(0.7,0.7,0.9,0.9);
    leg6->AddEntry(h_3pi3pi_0, "#Delta m_{B} #in [0,100] MeV");
    leg6->AddEntry(h_3pi3pi_1, "#Delta m_{B} #in ]100,250] MeV");
    leg6->AddEntry(h_3pi3pi_2, "#Delta m_{B} > 250 MeV");
    leg6->Draw("same");
    c6.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_3pi3pi.pdf");

    TCanvas c7 = new TCanvas();
    c7.cd();

    h_3pi3pipi0_0->SetLineColor(kBlue);
    h_3pi3pipi0_0->SetFillColorAlpha(kBlue, 0.25);

    h_3pi3pipi0_1->SetLineColor(kGreen+1);
    h_3pi3pipi0_1->SetFillColorAlpha(kGreen+1, 0.25);

    h_3pi3pipi0_2->SetLineColor(kRed);
    h_3pi3pipi0_2->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pipi0_0->DrawNormalized();
    h_3pi3pipi0_1->DrawNormalized("same");
    h_3pi3pipi0_2->DrawNormalized("same");

    leg6->Draw("same");
    c7.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_3pi3pipi0.pdf");

    TCanvas c8 = new TCanvas();
    c8.cd();

    h_3pi3pi2pi0_0->SetLineColor(kBlue);
    h_3pi3pi2pi0_0->SetFillColorAlpha(kBlue, 0.25);

    h_3pi3pi2pi0_1->SetLineColor(kGreen+1);
    h_3pi3pi2pi0_1->SetFillColorAlpha(kGreen+1, 0.25);

    h_3pi3pi2pi0_2->SetLineColor(kRed);
    h_3pi3pi2pi0_2->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pi2pi0_0->DrawNormalized();
    h_3pi3pi2pi0_1->DrawNormalized("same");
    h_3pi3pi2pi0_2->DrawNormalized("same");

    leg6->Draw("same");
    c8.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_3pi3pi2pi0.pdf");

    TCanvas c9 = new TCanvas();
    c9.cd();

    h_all_mc_0->SetLineColor(kBlue);
    h_all_mc_0->SetFillColorAlpha(kBlue, 0.25);

    h_all_mc_1->SetLineColor(kGreen+1);
    h_all_mc_1->SetFillColorAlpha(kGreen+1, 0.25);

    h_all_mc_2->SetLineColor(kRed);
    h_all_mc_2->SetFillColorAlpha(kRed, 0.25);

    h_all_mc_0->DrawNormalized();
    h_all_mc_1->DrawNormalized("same");
    h_all_mc_2->DrawNormalized("same");

    leg6->Draw("same");
    c9.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_all_mc.pdf");

    TCanvas c10 = new TCanvas();
    c10.cd();

    h_rs_0->SetLineColor(kBlue);
    h_rs_0->SetFillColorAlpha(kBlue, 0.25);

    h_rs_1->SetLineColor(kGreen+1);
    h_rs_1->SetFillColorAlpha(kGreen+1, 0.25);

    h_rs_2->SetLineColor(kRed);
    h_rs_2->SetFillColorAlpha(kRed, 0.25);

    h_rs_0->DrawNormalized();
    h_rs_1->DrawNormalized("same");
    h_rs_2->DrawNormalized("same");

    leg6->Draw("same");
    c10.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_rs.pdf");

    TCanvas c11 = new TCanvas();
    c11.cd();

    h_ws_0->SetLineColor(kBlue);
    h_ws_0->SetFillColorAlpha(kBlue, 0.25);

    h_ws_1->SetLineColor(kGreen+1);
    h_ws_1->SetFillColorAlpha(kGreen+1, 0.25);

    h_ws_2->SetLineColor(kRed);
    h_ws_2->SetFillColorAlpha(kRed, 0.25);

    h_ws_0->DrawNormalized();
    h_ws_1->DrawNormalized("same");
    h_ws_2->DrawNormalized("same");

    leg6->Draw("same");
    c11.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_ws.pdf");

    std::ofstream file3("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_fractions.tex");
    file3 << " \\begin{table}[!htbp]" << std::endl;
    file3 << " \\centering " << std::endl;
    file3 << "\\tiny" << std::endl;    
    file3 << " \\begin{tabular}{|c|c|c|c|c|c|c|}" << std::endl;
    file3 << " \\hline" << std::endl;
    file3 << "  & $3\\pi3\\pi$ MC & $3\\pi3\\pi \\pi^0$ MC & $3\\pi3\\pi 2\\pi^0$ MC & All MC & RS data & WS data \\\\ " << std::endl;
    file3 << "\\hline" << std::endl;    

    Double_t N_3pi3pi_0_num = t_mc->GetEntries(pass_GSL+"component==0"+error0);
    Double_t N_3pi3pi_1_num = t_mc->GetEntries(pass_GSL+"component==0"+error1);
    Double_t N_3pi3pi_2_num = t_mc->GetEntries(pass_GSL+"component==0"+error2);
    Double_t N_3pi3pi_den = t_mc->GetEntries(pass_GSL+"component==0");

    Double_t N_3pi3pipi0_0_num = t_mc->GetEntries(pass_GSL+"component==1"+error0);
    Double_t N_3pi3pipi0_1_num = t_mc->GetEntries(pass_GSL+"component==1"+error1);
    Double_t N_3pi3pipi0_2_num = t_mc->GetEntries(pass_GSL+"component==1"+error2);
    Double_t N_3pi3pipi0_den = t_mc->GetEntries(pass_GSL+"component==1");

    Double_t N_3pi3pi2pi0_0_num = t_mc->GetEntries(pass_GSL+"component==2"+error0);
    Double_t N_3pi3pi2pi0_1_num = t_mc->GetEntries(pass_GSL+"component==2"+error1);
    Double_t N_3pi3pi2pi0_2_num = t_mc->GetEntries(pass_GSL+"component==2"+error2);
    Double_t N_3pi3pi2pi0_den = t_mc->GetEntries(pass_GSL+"component==2");

    Double_t N_all_mc_0_num = t_mc->GetEntries(pass_GSL+error0);
    Double_t N_all_mc_1_num = t_mc->GetEntries(pass_GSL+error1);
    Double_t N_all_mc_2_num = t_mc->GetEntries(pass_GSL+error2);
    Double_t N_all_mc_den = t_mc->GetEntries(pass_GSL);

    Double_t N_rs_0_num = t_rs->GetEntries(pass_GSL+error0);
    Double_t N_rs_1_num = t_rs->GetEntries(pass_GSL+error1);
    Double_t N_rs_2_num = t_rs->GetEntries(pass_GSL+error2);
    Double_t N_rs_den = t_rs->GetEntries(pass_GSL);

    Double_t N_ws_0_num = t_ws->GetEntries(pass_GSL+error0);
    Double_t N_ws_1_num = t_ws->GetEntries(pass_GSL+error1);
    Double_t N_ws_2_num = t_ws->GetEntries(pass_GSL+error2);
    Double_t N_ws_den = t_ws->GetEntries(pass_GSL);

    file3 << Form(" $\\Delta m_B \\in [0,100]$\\,MeV & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.2f \\pm %.2f$\\,\\% & $%.2f \\pm %.2f$\\,\\% \\\\", (N_3pi3pi_0_num/N_3pi3pi_den)*100, eps_error(N_3pi3pi_0_num,N_3pi3pi_den)*100, (N_3pi3pipi0_0_num/N_3pi3pipi0_den)*100, eps_error(N_3pi3pipi0_0_num,N_3pi3pipi0_den)*100, (N_3pi3pi2pi0_0_num/N_3pi3pi2pi0_den)*100, eps_error(N_3pi3pi2pi0_0_num,N_3pi3pi2pi0_den)*100, (N_all_mc_0_num/N_all_mc_den)*100, eps_error(N_all_mc_0_num,N_all_mc_den)*100, (N_rs_0_num/N_rs_den)*100, eps_error(N_rs_0_num,N_rs_den)*100, (N_ws_0_num/N_ws_den)*100, eps_error(N_ws_0_num,N_ws_den)*100 );
    file3 << Form(" $\\Delta m_B \\in ]100,250]$\\,MeV & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% &  $%.2f \\pm %.2f$\\,\\% & $%.2f \\pm %.2f$\\,\\% \\\\", (N_3pi3pi_1_num/N_3pi3pi_den)*100, eps_error(N_3pi3pi_1_num,N_3pi3pi_den)*100, (N_3pi3pipi0_1_num/N_3pi3pipi0_den)*100, eps_error(N_3pi3pipi0_1_num,N_3pi3pipi0_den)*100, (N_3pi3pi2pi0_1_num/N_3pi3pi2pi0_den)*100, eps_error(N_3pi3pi2pi0_1_num,N_3pi3pi2pi0_den)*100, (N_all_mc_1_num/N_all_mc_den)*100, eps_error(N_all_mc_1_num,N_all_mc_den)*100, (N_rs_1_num/N_rs_den)*100, eps_error(N_rs_1_num,N_rs_den)*100, (N_ws_1_num/N_ws_den)*100, eps_error(N_ws_1_num,N_ws_den)*100 );
    file3 << Form(" $\\Delta m_B > 250$\\,MeV & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.1f \\pm %.1f$\\,\\% & $%.2f \\pm %.2f$\\,\\% & $%.2f \\pm %.2f$\\,\\% \\\\", (N_3pi3pi_2_num/N_3pi3pi_den)*100, eps_error(N_3pi3pi_2_num,N_3pi3pi_den)*100, (N_3pi3pipi0_2_num/N_3pi3pipi0_den)*100, eps_error(N_3pi3pipi0_2_num,N_3pi3pipi0_den)*100, (N_3pi3pi2pi0_2_num/N_3pi3pi2pi0_den)*100, eps_error(N_3pi3pi2pi0_2_num,N_3pi3pi2pi0_den)*100, (N_all_mc_2_num/N_all_mc_den)*100, eps_error(N_all_mc_2_num,N_all_mc_den)*100, (N_rs_2_num/N_rs_den)*100, eps_error(N_rs_2_num,N_rs_den)*100, (N_ws_2_num/N_ws_den)*100, eps_error(N_ws_2_num,N_ws_den)*100 );

    file3 << "\\hline" << std::endl;
    file3 << "\\end{tabular}" << std::endl;
    file3 << "\\end{table}" << std::endl;
    file3.close();

    std::ofstream file4("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_error_categories_resolution.tex");
    file4 << " \\begin{table}[!htbp]" << std::endl;
    file4 << " \\centering " << std::endl;
    file4 << " \\begin{tabular}{|c|c|c|c|c|}" << std::endl;
    file4 << " \\hline" << std::endl;
    file4 << "  & $3\\pi3\\pi$ MC & $3\\pi3\\pi \\pi^0$ MC & $3\\pi3\\pi 2\\pi^0$ MC & All MC \\\\ " << std::endl;
    file4 << "\\hline" << std::endl;    

    Double_t res_3pi3pi_0 = FWHM(h_3pi3pi_0)/2.4;
    Double_t res_3pi3pi_1 = FWHM(h_3pi3pi_1)/2.4;
    Double_t res_3pi3pi_2 = FWHM(h_3pi3pi_2)/2.4;

    Double_t res_3pi3pipi0_0 = FWHM(h_3pi3pipi0_0)/2.4;
    Double_t res_3pi3pipi0_1 = FWHM(h_3pi3pipi0_1)/2.4;
    Double_t res_3pi3pipi0_2 = FWHM(h_3pi3pipi0_2)/2.4;

    Double_t res_3pi3pi2pi0_0 = FWHM(h_3pi3pi2pi0_0)/2.4;
    Double_t res_3pi3pi2pi0_1 = FWHM(h_3pi3pi2pi0_1)/2.4;
    Double_t res_3pi3pi2pi0_2 = FWHM(h_3pi3pi2pi0_2)/2.4;

    Double_t res_all_mc_0 = FWHM(h_all_mc_0)/2.4;
    Double_t res_all_mc_1 = FWHM(h_all_mc_1)/2.4;
    Double_t res_all_mc_2 = FWHM(h_all_mc_2)/2.4;

    file4 << Form(" $\\Delta m_B \\in [0,100]$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV \\\\ ", res_3pi3pi_0, res_3pi3pipi0_0, res_3pi3pi2pi0_0, res_all_mc_0);
    file4 << Form(" $\\Delta m_B \\in ]100,150]$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV \\\\ ", res_3pi3pi_1, res_3pi3pipi0_1, res_3pi3pi2pi0_1, res_all_mc_1);
    file4 << Form(" $\\Delta m_B > 250$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV & $%.0f \\pm 40$\\,MeV \\\\ ", res_3pi3pi_2, res_3pi3pipi0_2, res_3pi3pi2pi0_2, res_all_mc_2);

    file4 << "\\hline" << std::endl;
    file4 << "\\end{tabular}" << std::endl;
    file4 << "\\end{table}" << std::endl;
    file4.close();

    auto c1 = new TCanvas(); c1->DrawFrame(0.,0.,1.,1.);
    auto c2 = new TCanvas(); c2->DrawFrame(0.,0.,2.,2.);
    auto c3 = new TCanvas(); c3->DrawFrame(0.,0.,3.,3.);
    auto c4 = new TCanvas(); c4->DrawFrame(0.,0.,4.,4.);

    auto FourPads = new TCanvas("FourPads","FourPads");
    FourPads->Divide(2,2);
    FourPads->cd(1) ; c1->DrawClonePad();

    h_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi->SetTitle("All events: 2016-2018");

    h_3pi3pi->SetLineColor(kBlue);
    h_3pi3pi->SetFillColorAlpha(kBlue, 0.25);

    h_rs->SetLineColor(kBlack);
    h_rs->SetFillColorAlpha(kBlack, 0.25);

    h_ws->SetLineColor(kRed);
    h_ws->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pi->DrawNormalized();
    h_rs->DrawNormalized("same");
    h_ws->DrawNormalized("same");

    TLegend* leg7 = new TLegend(0.7,0.7,0.9,0.9);
    leg7->AddEntry(h_3pi3pi, "3#pi3#pi MC", "lf");
    leg7->AddEntry(h_rs, "RS data", "lf");
    leg7->AddEntry(h_ws, "WS data", "lf");
    leg7->Draw("same");

    FourPads->cd(2) ; c2->DrawClonePad();

    h_3pi3pi_0->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_0->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_0->SetTitle("#Delta m_{B} #in [0,100] MeV: 2016-2018");

    h_3pi3pi_0->SetLineColor(kBlue);
    h_3pi3pi_0->SetFillColorAlpha(kBlue, 0.25);

    h_rs_0->SetLineColor(kBlack);
    h_rs_0->SetFillColorAlpha(kBlack, 0.25);

    h_ws_0->SetLineColor(kRed);
    h_ws_0->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pi_0->DrawNormalized();
    h_rs_0->DrawNormalized("same");
    h_ws_0->DrawNormalized("same");
    leg7->Draw("same");
    
    FourPads->cd(3) ; c3->DrawClonePad();

    h_3pi3pi_1->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_1->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_1->SetTitle("#Delta m_{B} #in ]100,250] MeV: 2016-2018");

    h_3pi3pi_1->SetLineColor(kBlue);
    h_3pi3pi_1->SetFillColorAlpha(kBlue, 0.25);

    h_rs_1->SetLineColor(kBlack);
    h_rs_1->SetFillColorAlpha(kBlack, 0.25);

    h_ws_1->SetLineColor(kRed);
    h_ws_1->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pi_1->DrawNormalized();
    h_rs_1->DrawNormalized("same");
    h_ws_1->DrawNormalized("same");
    leg7->Draw("same");

    FourPads->cd(4) ; c4->DrawClonePad();

    h_3pi3pi_2->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_2->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_2->SetTitle("#Delta m_{B} > 250 MeV: 2016-2018");

    h_3pi3pi_2->SetLineColor(kBlue);
    h_3pi3pi_2->SetFillColorAlpha(kBlue, 0.25);

    h_rs_2->SetLineColor(kBlack);
    h_rs_2->SetFillColorAlpha(kBlack, 0.25);

    h_ws_2->SetLineColor(kRed);
    h_ws_2->SetFillColorAlpha(kRed, 0.25);

    h_3pi3pi_2->DrawNormalized();
    h_rs_2->DrawNormalized("same");
    h_ws_2->DrawNormalized("same");
    leg7->Draw("same");

    FourPads->SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_vs_data_error_categories.pdf");

    gsl_pull_distributions(t_mc, pass_GSL);

    /*

    Bool_t isMC = false;
    if( (species == 10) || (species == 11) || (species == 12) || (species == 1) || (species==4) || (species == 9))
    {
        isMC = true;
    }

    int dimM = 22; // number of measured parameters
    int dimX = 23; // number of unknown parameters
    int dimC = 24; // number of exact constraints

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

    std::vector<TH1D> histos;
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        histos.push_back(TH1D(Form("h_%i",i), Form("h_%i",i), 100, -pow(10,-6), pow(10,-6)) );
        t->Draw(Form("df_F_%i >> h_%i",i,i),df_status);
    }

    TString name[] = {
        "dL/dPVx (1/mm)",
        "dL/dPVy (1/mm)",
        "dL/dPVz (1/mm)",
        "dL/dDV1x (1/mm)",
        "dL/dDV1y (1/mm)",
        "dL/DV1z (1/mm)",
        "dL/dp3pi1x (1/MeV)",
        "dL/dp3pi1y (1/MeV)",
        "dL/dp3pi1z (1/MeV)",
        "dL/dE3pi1 (1/MeV)",
        "dL/dDV2x (1/mm)",
        "dL/dDV2y (1/mm)",
        "dL/dDV2z (1/mm)",
        "dL/dp3pi2x (1/MeV)",
        "dL/dp3pi2y (1/MeV)",
        "dL/dp3pi2z (1/MeV)",
        "dL/dE3pi2 (1/MeV)",
        "dL/dRPx (1/mm)",
        "dL/dRPy (1/mm)",
        "dL/dpKx (1/MeV)",
        "dL/dpKy (1/MeV)",
        "dL/dpKz (1/MeV)",
        "dL/dBVx (1/mm)",
        "dL/dBVy (1/mm)",
        "dL/dBVz (1/mm)",
        "dL/dpBx (1/MeV)",
        "dL/dpBy (1/MeV)",
        "dL/dpBz (1/MeV)",
        "dL/dM_{B}^{2} (1/MeV^{2})",
        "dL/dptau1x (1/MeV)",
        "dL/dptau1y (1/MeV)",
        "dL/dptau1z (1/MeV)",
        "dL/dEtau1 (1/MeV)",
        "dL/dpnu1x (1/MeV)",
        "dL/dpnu1y (1/MeV)",
        "dL/dpnu1z (1/MeV)",
        "dL/dEnu1 (1/MeV)",
        "dL/dptau2x (1/MeV)",
        "dL/dptau2y (1/MeV)",
        "dL/dptau2z (1/MeV)",
        "dL/dEtau2 (1/MeV)",
        "dL/dpnu2x (1/MeV)",
        "dL/dpnu2y (1/MeV)",
        "dL/dpnu2z (1/MeV)",
        "dL/dEnu2 (1/MeV)",
        "BV must lie in K+ trajectory (x,z) (MeV/mm)",
        "BV must lie in K+ trajectory (y,z) (MeV/mm)",
        "ptau1 must point back to BV (x,z) (MeV/mm)",
        "ptau1 must point back to BV (y,z) (MeV/mm)",
        "ptau2 must point back to DV2 (x,z) (MeV/mm)",
        "ptau2 must point back to DV2 (y,z) (MeV/mm)",
        "pB must point back to PV (x,z) (MeV/mm)",
        "pB must point back to PV (y,z) (MeV/mm)",
        "px conservation in DV1 (MeV)",
        "py conservation in DV1 (MeV)",
        "pz conservation in DV1 (MeV)",
        "E conservation in DV1 (MeV)",
        "Tau+ mass constraint (MeV)",
        "Anti-nu mass constraint (MeV)",
        "px conservation in DV2 (MeV)",
        "py conservation in DV2 (MeV)",
        "pz conservation in DV2 (MeV)",
        "E conservation in DV2 (MeV)",
        "Tau- mass constraint (MeV)",
        "Nu mass constraint (MeV)",
        "px conservation in BV (MeV)",
        "py conservation in BV (MeV)",
        "pz conservation in BV (MeV)",
        "E conservation in BV (MeV)"
    };

    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        TCanvas c = new TCanvas();
        histos[i].GetXaxis()->SetTitle(name[i]);
        histos[i].GetYaxis()->SetTitle( TString::Format("Events /(%g)",(histos[i].GetXaxis()->GetXmax() - histos[i].GetXaxis()->GetXmin())/100) );
        histos[i].SetTitle(name[i]);
        histos[i].Draw();
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/Species_%i/F%i.pdf",species,i));
    }

    Double_t var_min[] = {
        0.4, // PVx
        -0.4, // PVy
        -200,  // PVz
        -10, // DV1x
        -10, // DV1y
        -100, // DV1z
        -10000, // p3pi1x
        -10000, // p3pi1y
        0, // p3pi1z
        0, // E3pi1
        -10, // DV2x
        -10, // DV2y
        -100, // DV2z
        -10000, // p3pi2x
        -10000, // p3pi2y
        0, // p3pi2z
        0, // E3pi2
        -10, // RPx
        -10, // RPy
        -50000, // pKx
        -50000, // pKy
        0, // pKz
        -10, // BVx
        -10, // BVy
        -200, // BVz
        -100000, // pBx
        -100000, // pBy
        0, // pBz
        0, // mB^2
        -20000, // ptau1x
        -20000, // ptau1y
        0, // ptau1z
        0, // Etau1
        -10000, // pnu1x
        -10000, // pnu1y
        -100, // pnu1z
        -100, // Enu1
        -20000, // ptau2x
        -20000, // ptau2y
        0, // ptau2z
        0, // Etau2
        -10000, // pnu2x
        -10000, // pnu2y
        -100, // pnu2z
        -100 // Enu2
    };

    Double_t var_max[] = {
        1.4, // PVx
        0.4, // PVy
        200,  // PVz
        10, // DV1x
        10, // DV1y
        200, // DV1z
        20000, // p3pi1x
        20000, // p3pi1y
        200000, // p3pi1z
        200000, // E3pi1
        10, // DV2x
        10, // DV2y
        200, // DV2z
        20000, // p3pi2x
        20000, // p3pi2y
        200000, // p3pi2z
        200000, // E3pi2
        10, // RPx
        10, // RPy
        50000, // pKx
        50000, // pKy
        500000, // pKz
        10, // BVx
        10, // BVy
        300, // BVz
        100000, // pBx
        100000, // pBy
        1000000, // pBz
        pow(10,8), // mB^2
        20000, // ptau1x
        20000, // ptau1y
        200000, // ptau1z
        200000, // Etau1
        10000, // pnu1x
        10000, // pnu1y
        100000, // pnu1z
        100000, // Enu1
        20000, // ptau2x
        20000, // ptau2y
        500000, // ptau2z
        500000, // Etau2
        10000, // pnu2x
        10000, // pnu2y
        100000, // pnu2z
        100000 // Enu2
    };

    // Bias in fit parameter

    TString DTF_name[] = {
        "Bp_dtf_PV_X[0]",
        "Bp_dtf_PV_Y[0]",
        "Bp_dtf_PV_Z[0]",
        "Bp_dtf_taup_ENDVERTEX_X[0]",
        "Bp_dtf_taup_ENDVERTEX_Y[0]",
        "Bp_dtf_taup_ENDVERTEX_Z[0]",
        "(Bp_dtf_tauminus_piplus_PX[0] + Bp_dtf_tauminus_piplus_0_PX[0] + Bp_dtf_tauminus_piplus_1_PX[0])",
        "(Bp_dtf_tauminus_piplus_PY[0] + Bp_dtf_tauminus_piplus_0_PY[0] + Bp_dtf_tauminus_piplus_1_PY[0])",
        "(Bp_dtf_tauminus_piplus_PZ[0] + Bp_dtf_tauminus_piplus_0_PZ[0] + Bp_dtf_tauminus_piplus_1_PZ[0])",
        "(Bp_dtf_tauminus_piplus_PE[0] + Bp_dtf_tauminus_piplus_0_PE[0] + Bp_dtf_tauminus_piplus_1_PE[0])",
        "Bp_dtf_taum_ENDVERTEX_X[0]",
        "Bp_dtf_taum_ENDVERTEX_Y[0]",
        "Bp_dtf_taum_ENDVERTEX_Z[0]",
        "(Bp_dtf_tauminus_0_piplus_PX[0] + Bp_dtf_tauminus_0_piplus_0_PX[0] + Bp_dtf_tauminus_0_piplus_1_PX[0])",
        "(Bp_dtf_tauminus_0_piplus_PY[0] + Bp_dtf_tauminus_0_piplus_0_PY[0] + Bp_dtf_tauminus_0_piplus_1_PY[0])",
        "(Bp_dtf_tauminus_0_piplus_PZ[0] + Bp_dtf_tauminus_0_piplus_0_PZ[0] + Bp_dtf_tauminus_0_piplus_1_PZ[0])",
        "(Bp_dtf_tauminus_0_piplus_PE[0] + Bp_dtf_tauminus_0_piplus_0_PE[0] + Bp_dtf_tauminus_0_piplus_1_PE[0])",
        "Bp_dtf_Bp_ENDVERTEX_X[0]",
        "Bp_dtf_Bp_ENDVERTEX_Y[0]",
        "Bp_dtf_Kplus_PX[0]",
        "Bp_dtf_Kplus_PY[0]",
        "Bp_dtf_Kplus_PZ[0]",
        "Bp_dtf_Bp_ENDVERTEX_X[0]",
        "Bp_dtf_Bp_ENDVERTEX_Y[0]",
        "Bp_dtf_Bp_ENDVERTEX_Z[0]",
        "Bp_dtf_PX[0]",
        "Bp_dtf_PY[0]",
        "Bp_dtf_PZ[0]",
        "pow(Bp_dtf_M[0],2)",
        "Bp_dtf_tauminus_PX[0]",
        "Bp_dtf_tauminus_PY[0]",
        "Bp_dtf_tauminus_PZ[0]",
        "Bp_dtf_tauminus_PE[0]",
        "(Bp_dtf_tauminus_PX[0] - Bp_dtf_tauminus_piplus_PX[0] - Bp_dtf_tauminus_piplus_0_PX[0] - Bp_dtf_tauminus_piplus_1_PX[0])",
        "(Bp_dtf_tauminus_PY[0] - Bp_dtf_tauminus_piplus_PY[0] - Bp_dtf_tauminus_piplus_0_PY[0] - Bp_dtf_tauminus_piplus_1_PY[0])",
        "(Bp_dtf_tauminus_PZ[0] - Bp_dtf_tauminus_piplus_PZ[0] - Bp_dtf_tauminus_piplus_0_PZ[0] - Bp_dtf_tauminus_piplus_1_PZ[0])",
        "(Bp_dtf_tauminus_PE[0] - Bp_dtf_tauminus_piplus_PE[0] - Bp_dtf_tauminus_piplus_0_PE[0] - Bp_dtf_tauminus_piplus_1_PE[0])",
        "Bp_dtf_tauminus_0_PX[0]",
        "Bp_dtf_tauminus_0_PY[0]",
        "Bp_dtf_tauminus_0_PZ[0]",
        "Bp_dtf_tauminus_0_PE[0]",
        "(Bp_dtf_tauminus_0_PX[0] - Bp_dtf_tauminus_0_piplus_PX[0] - Bp_dtf_tauminus_0_piplus_0_PX[0] - Bp_dtf_tauminus_0_piplus_1_PX[0])",
        "(Bp_dtf_tauminus_0_PY[0] - Bp_dtf_tauminus_0_piplus_PY[0] - Bp_dtf_tauminus_0_piplus_0_PY[0] - Bp_dtf_tauminus_0_piplus_1_PY[0])",
        "(Bp_dtf_tauminus_0_PZ[0] - Bp_dtf_tauminus_0_piplus_PZ[0] - Bp_dtf_tauminus_0_piplus_0_PZ[0] - Bp_dtf_tauminus_0_piplus_1_PZ[0])",
        "(Bp_dtf_tauminus_0_PE[0] - Bp_dtf_tauminus_0_piplus_PE[0] - Bp_dtf_tauminus_0_piplus_0_PE[0] - Bp_dtf_tauminus_0_piplus_1_PE[0])"
    };

    std::vector<TH1D> df_histos;
    std::vector<TH1D> DTF_histos;

    for(int i = 0; i < dimM+dimX; i++)
    {
        // Double_t x_min = t->GetMinimum(df_name[i]);
        // Double_t x_max = t->GetMaximum(df_name[i]);

        // if(df_name[i] == "df_BVx")
        // {
        //     x_min = -6;
        // }
        // else if(df_name[i] == "df_BVy")
        // {
        //     x_min = -8;
        // }
        // else if(df_name[i] == "df_BVz")
        // {
        //     x_min = -150;
        // }
        // else if(df_name[i] == "df_Bp_PX")
        // {
        //     x_min = -50000;
        //     x_max = 50000;
        // }
        // else if(df_name[i] == "df_Bp_PY")
        // {
        //     x_min = -50000;
        //     x_max = 50000;
        // }
        // else if(df_name[i] == "df_Bp_PZ")
        // {
        //     x_min = 0;
        //     x_max = 1000000;
        // }
        // else if(df_name[i] == "df_taum_PX")
        // {
        //     x_min = -15000;
        //     x_max = 15000;
        // }
        // else if(df_name[i] == "df_taum_PY")
        // {
        //     x_min = -15000;
        //     x_max = 15000;
        // }
        // else if(df_name[i] == "df_taum_PZ")
        // {
        //     x_min = -1000;
        //     x_max = 1500000;
        // }
        // else if(df_name[i] == "df_taum_PE")
        // {
        //     x_min = -1000;
        //     x_max = 1500000;
        // }
        // else if(df_name[i] == "df_nutau_PX")
        // {
        //     x_min = -15000;
        //     x_max = 15000;
        // }
        // else if(df_name[i] == "df_nutau_PY")
        // {
        //     x_min = -15000;
        //     x_max = 15000;
        // }
        // else if(df_name[i] == "df_antinutau_PZ")
        // {
        //     x_min = -1000;
        //     x_max = 200000;
        // }
        // else if(df_name[i] == "df_antinutau_PE")
        // {
        //     x_min = -1000;
        //     x_max = 200000;
        // }

        // else if(df_name[i] == "df_nutau_PZ")
        // {
        //     x_min = -1000;
        //     x_max = 200000;
        // }
        // else if(df_name[i] == "df_nutau_PE")
        // {
        //     x_min = -1000;
        //     x_max = 200000;
        // }

        df_histos.push_back( TH1D(Form("df_h_%i",i), Form("df_h_%i",i), 100, -range[i], range[i]) );
        DTF_histos.push_back( TH1D(Form("DTF_h_%i",i), Form("DTF_h_%i",i), 100, -range[i], range[i]) );
    }

    for(int i = 0; i < dimM+dimX; i++)
    {
        if(species == 4)
        {
            t->Draw(df_name[i]+" - "+true_name_DDK[i]+Form(" >> df_h_%i",i),df_status);
            t->Draw(DTF_name_DDK[i]+" - "+true_name_DDK[i]+Form(" >> DTF_h_%i",i),passDTF);                
        }
        else if(species == 9)
        {
            t->Draw(df_name[i]+" - "+true_name_D0D0K[i]+Form(" >> df_h_%i",i),df_status);
            //t->Draw(DTF_name_DDK[i]+" - "+true_name_DDK[i]+Form(" >> DTF_h_%i",i),df_status);  
        }
        else
        {
            // t->Draw(df_name[i]+Form(" >> df_h_%i",i),df_status);
            t->Draw(df_name[i]+" - "+true_name[i]+Form(" >> df_h_%i",i),df_status);
            t->Draw(DTF_name[i]+" - "+true_name[i]+Form(" >> DTF_h_%i",i),passDTF);
        }
    }

    for(int i = 0; i < dimM+dimX;i++)
    {
        TCanvas c3 = new TCanvas();
        c3.cd();
        df_histos[i].GetXaxis()->SetTitle(name1[i]);
        df_histos[i].GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(df_histos[i].GetXaxis()->GetXmax() - df_histos[i].GetXaxis()->GetXmin())/100) );
        df_histos[i].SetTitle(" ");
        DTF_histos[i].GetXaxis()->SetTitle(name1[i]);
        DTF_histos[i].GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(DTF_histos[i].GetXaxis()->GetXmax() - DTF_histos[i].GetXaxis()->GetXmin())/100) );

        DTF_histos[i].SetLineColor(kBlack);
        DTF_histos[i].SetFillColorAlpha(kBlack, 0.25);
        df_histos[i].SetLineColor(kBlue);
        df_histos[i].SetFillColorAlpha(kBlue, 0.25);

        df_histos[i].SetTitle("Fit - true");
        df_histos[i].DrawNormalized();
        DTF_histos[i].DrawNormalized("same");

        if( (df_histos[i].GetMaximum())/(df_histos[i].Integral()) > (DTF_histos[i].GetMaximum())/(DTF_histos[i].Integral()) )
        {
            df_histos[i].SetTitle("Fit - true");
            df_histos[i].DrawNormalized();
            DTF_histos[i].DrawNormalized("same");
        }
        else
        {
            DTF_histos[i].SetTitle("Fit - true");
            DTF_histos[i].DrawNormalized();
            df_histos[i].DrawNormalized("same");
        }

        TLegend* leg3 = new TLegend(0.6, 0.6, 0.85, 0.85);
        leg3->SetBorderSize(0);
        leg3->SetTitle(0.03);
        leg3->AddEntry(&DTF_histos[i],Form("DTF: m %0.2f w %0.2f",DTF_histos[i].GetMean(),DTF_histos[i].GetStdDev()),"f");
        leg3->AddEntry(&df_histos[i],Form("Fitter: m %0.2f w %0.2f",df_histos[i].GetMean(),df_histos[i].GetStdDev()),"f");
        leg3->Draw("same");

        c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/Species_%i/var_%i.pdf",species,i));
    }


    gStyle->SetOptStat(1);

    std::vector<TH1D> df_histos1;
    for(int i = 0; i < dimM+dimX; i++)
    {
        df_histos1.push_back( TH1D(Form("df_h1_%i",i), Form("df_h1_%i",i), 100, -10, 10) );
    }
    */
}

Double_t eps_error(Double_t Num, Double_t Den)
{
    return (Num/Den)*sqrt( 1./Num + 1./Den );
}

double FWHM(TH1* h1)
{
   int bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
   int bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
   double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
   return fwhm;
}

std::vector<TH1D*> gsl_pull_distributions(TChain* t, TCut pass_GSL)
{
    gStyle->SetOptStat(1);
    std::vector<TH1D*> pull_histos;
    std::vector<TCanvas> pull_c;
    for(int i = 0; i < dimM+dimX; i++)
    {
        pull_histos.push_back( new TH1D(Form("h_gsl_pull_%i",i), Form("h_gsl_pull_%i",i), 100, -10, 10) );

        t->Draw( "("+df_name[i]+" - "+true_name[i]+")/"+df_err_name[i]+Form(" >> h_gsl_pull_%i",i), pass_GSL+"component==0");

        pull_histos[i]->GetXaxis()->SetTitle(label_name[i]);
        pull_histos[i]->GetYaxis()->SetTitle(Form("Entries / (%.2f)", (pull_histos[i]->GetXaxis()->GetXmax() - pull_histos[i]->GetXaxis()->GetXmin())/pull_histos[i]->GetNbinsX() ));
        pull_histos[i]->SetTitle(title_name[i]+" : 3#pi3#pi MC (2016-2018)");

        pull_histos[i]->SetLineColor(kBlue);
        pull_histos[i]->SetFillColorAlpha(kBlue, 0.25);

        TCanvas c = new TCanvas(Form("c_%i",i), Form("c_%i",i));
        c.cd();
        pull_histos[i]->Draw();
        c.SaveAs("/panfs/felician/B2Ktautau/workflow/exact_constraints/GSL_pulls/"+title_name[i]+".pdf");
    }
    return pull_histos;
}