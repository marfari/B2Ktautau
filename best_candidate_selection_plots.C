
void best_candidate_selection_plots()
{
    gStyle->SetOptStat(0);
    //////////////////////////////////////////////////////////////// Ktautau ///////////////////////////////////////////////////////////////////////////
    // TM MC
    TFile* f_ktautau_mc_tm = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_1/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_ktautau_mc_tm = (TTree*)f_ktautau_mc_tm->Get("DecayTree");

    // Un-TM MC
    TFile* f_ktautau_mc = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_ktautau_mc = (TTree*)f_ktautau_mc->Get("DecayTree");

    // Best MC
    TFile* f_ktautau_mc_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_10/is_best_candidate_branch.root");
    TTree* t_ktautau_mc_best = (TTree*)f_ktautau_mc_best->Get("DecayTree");

    Int_t N_mc = t_ktautau_mc->GetEntries();
    Int_t N_mc_best = t_ktautau_mc_best->GetEntries();
    if(N_mc != N_mc_best)
    {
        cout << "Wrong number of entries in Ktautau MC" << endl;
        return;
    }
    else
    {
        t_ktautau_mc->AddFriend(t_ktautau_mc_best);
    }

    // RS data
    TFile* f_ktautau_rs_data = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_ktautau_rs_data = (TTree*)f_ktautau_rs_data->Get("DecayTree");

    // Best RS data
    // TFile* f_ktautau_rs_data_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_2/is_best_candidate_branch.root");
    // TTree* t_ktautau_rs_data_best = (TTree*)f_ktautau_rs_data_best->Get("DecayTree");

    // Int_t N_rs_data = t_ktautau_rs_data->GetEntries();
    // Int_t N_rs_data_best = t_ktautau_rs_data_best->GetEntries();
    // if(N_rs_data != N_rs_data_best)
    // {
    //     cout << "Wrong number of entries" << endl;
    //     return;
    // }
    // else
    // {
    //     t_ktautau_rs_data->AddFriend(t_ktautau_rs_data_best);
    // }

    // WS data
    TFile* f_ktautau_ws_data = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_ktautau_ws_data = (TTree*)f_ktautau_ws_data->Get("DecayTree");

    // WS data best
    // TFile* f_ktautau_ws_data_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_3/is_best_candidate_branch.root");
    // TTree* t_ktautau_ws_data_best = (TTree*)f_ktautau_ws_data_best->Get("DecayTree");

    // Int_t N_ws = t_ktautau_ws_data->GetEntries();
    // Int_t N_ws_best = t_ktautau_ws_data_best->GetEntries();
    // if(N_ws != N_ws_best)
    // {
    //     cout << "Wrong number of entries" << endl;
    //     return;
    // }
    // else
    // {
    //     t_ktautau_ws_data->AddFriend(t_ktautau_ws_data_best);
    // }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////// DDs //////////////////////////////////////////////////////////////////////////////
    // TM MC
    TFileCollection* fc_dds_mc_tm_2016 = new TFileCollection("fc_dds_mc_tm_2016", "fc_dds_mc_tm_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_7/pre_sel_tree.txt");
    TFileCollection* fc_dds_mc_tm_2017 = new TFileCollection("fc_dds_mc_tm_2017", "fc_dds_mc_tm_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_7/pre_sel_tree.txt");
    TFileCollection* fc_dds_mc_tm_2018 = new TFileCollection("fc_dds_mc_tm_2018", "fc_dds_mc_tm_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_7/pre_sel_tree.txt");

    TChain* t_dds_mc_tm = new TChain("DecayTree");
    TChain* t_dds_mc_tm_2017 = new TChain("DecayTree");
    TChain* t_dds_mc_tm_2018 = new TChain("DecayTree");

    t_dds_mc_tm->AddFileInfoList((TCollection*)fc_dds_mc_tm_2016->GetList());
    t_dds_mc_tm_2017->AddFileInfoList((TCollection*)fc_dds_mc_tm_2017->GetList());
    t_dds_mc_tm_2018->AddFileInfoList((TCollection*)fc_dds_mc_tm_2018->GetList());

    t_dds_mc_tm->GetEntries();
    t_dds_mc_tm_2017->GetEntries();
    t_dds_mc_tm_2018->GetEntries();

    t_dds_mc_tm->Add(t_dds_mc_tm_2017);
    t_dds_mc_tm->Add(t_dds_mc_tm_2018);

    // Un-TM MC
    TFileCollection* fc_dds_mc_2016 = new TFileCollection("fc_dds_mc_2016", "fc_dds_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_72/pre_sel_tree.txt");
    TFileCollection* fc_dds_mc_2017 = new TFileCollection("fc_dds_mc_2017", "fc_dds_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_72/pre_sel_tree.txt");
    TFileCollection* fc_dds_mc_2018 = new TFileCollection("fc_dds_mc_2018", "fc_dds_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_72/pre_sel_tree.txt");

    TChain* t_dds_mc = new TChain("DecayTree");
    TChain* t_dds_mc_2017 = new TChain("DecayTree");
    TChain* t_dds_mc_2018 = new TChain("DecayTree");

    t_dds_mc->AddFileInfoList((TCollection*)fc_dds_mc_2016->GetList());
    t_dds_mc_2017->AddFileInfoList((TCollection*)fc_dds_mc_2017->GetList());
    t_dds_mc_2018->AddFileInfoList((TCollection*)fc_dds_mc_2018->GetList());

    t_dds_mc->GetEntries();
    t_dds_mc_2017->GetEntries();
    t_dds_mc_2018->GetEntries();

    t_dds_mc->Add(t_dds_mc_2017);
    t_dds_mc->Add(t_dds_mc_2018);

    // MC best
    TFile* fc_dds_mc_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_72/is_best_candidate_branch.root");
    TChain* t_dds_mc_best = (TChain*)fc_dds_mc_best->Get("DecayTree");

    Int_t N_dds_mc = t_dds_mc->GetEntries();
    Int_t N_dds_mc_best = t_dds_mc_best->GetEntries();
    if(N_dds_mc != N_dds_mc_best)
    {
        cout << "Wrong number of entries in DDs MC" << endl;
        return;
    }
    else
    {
        t_dds_mc->AddFriend(t_dds_mc_best);
    }

    // Data
    TFileCollection* fc_dds_data_2016 = new TFileCollection("fc_dds_data_2016", "fc_dds_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_8/pre_sel_tree.txt");
    TFileCollection* fc_dds_data_2017 = new TFileCollection("fc_dds_data_2017", "fc_dds_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_8/pre_sel_tree.txt");
    TFileCollection* fc_dds_data_2018 = new TFileCollection("fc_dds_data_2018", "fc_dds_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_8/pre_sel_tree.txt");

    TChain* t_dds_data = new TChain("DecayTree");
    TChain* t_dds_data_2017 = new TChain("DecayTree");
    TChain* t_dds_data_2018 = new TChain("DecayTree");

    t_dds_data->AddFileInfoList((TCollection*)fc_dds_data_2016->GetList());
    t_dds_data_2017->AddFileInfoList((TCollection*)fc_dds_data_2017->GetList());
    t_dds_data_2018->AddFileInfoList((TCollection*)fc_dds_data_2018->GetList());

    t_dds_data->GetEntries();
    t_dds_data_2017->GetEntries();
    t_dds_data_2018->GetEntries();

    t_dds_data->Add(t_dds_data_2017);
    t_dds_data->Add(t_dds_data_2018);

    // Data best
    TFile* f_dds_data_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_8/is_best_candidate_branch.root");
    TChain* t_dds_data_best = (TChain*)f_dds_data_best->Get("DecayTree");

    Int_t N_dds_data = t_dds_data->GetEntries();
    Int_t N_dds_data_best = t_dds_data_best->GetEntries();
    if(N_dds_data != N_dds_data_best)
    {
        cout << "Wrong number of entries in DDs data" << endl;
        return;
    }
    else
    {
        t_dds_data->AddFriend(t_dds_data_best);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////// Cocktail MCs ////////////////////////////////////////////////////////////////////////////////////
    // B+ -> DD K+
    TFile* f_BuDDKp = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_BuDDKp = (TTree*)f_BuDDKp->Get("DecayTree");

    TFile* f_BuDDKp_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_100/is_best_candidate_branch.root");
    TTree* t_BuDDKp_best = (TTree*)f_BuDDKp_best->Get("DecayTree");

    Int_t N_100 = t_BuDDKp->GetEntries();
    Int_t N_100_best = t_BuDDKp_best->GetEntries();
    if(N_100 != N_100_best)
    {
        cout << "Wrong number of entries in BuDDKp MC" << endl;
        return;
    }
    else
    {
        t_BuDDKp->AddFriend(t_BuDDKp_best);
    }

    // B0 -> DD K+
    TFile* f_BdDDKp = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_BdDDKp = (TTree*)f_BdDDKp->Get("DecayTree");

    TFile* f_BdDDKp_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_110/is_best_candidate_branch.root");
    TTree* t_BdDDKp_best = (TTree*)f_BdDDKp_best->Get("DecayTree");

    Int_t N_110 = t_BdDDKp->GetEntries();
    Int_t N_110_best = t_BdDDKp_best->GetEntries();
    if(N_110 != N_110_best)
    {
        cout << "Wrong number of entries in BdDDKp MC" << endl;
        return;
    }
    else
    {
        t_BdDDKp->AddFriend(t_BdDDKp_best);
    }

    // Bs -> DD K+
    TFile* f_BsDDKp = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_BsDDKp = (TTree*)f_BsDDKp->Get("DecayTree");

    TFile* f_BsDDKp_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_120/is_best_candidate_branch.root");
    TTree* t_BsDDKp_best = (TTree*)f_BsDDKp_best->Get("DecayTree");

    Int_t N_120 = t_BsDDKp->GetEntries();
    Int_t N_120_best = t_BsDDKp_best->GetEntries();
    if(N_120 != N_120_best)
    {
        cout << "Wrong number of entries in BsDDKp MC" << endl;
        return;
    }
    else
    {
        t_BsDDKp->AddFriend(t_BsDDKp_best);
    }

    // B+ -> DD K0
    TFile* f_BuDDK0 = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_BuDDK0 = (TTree*)f_BuDDK0->Get("DecayTree");

    TFile* f_BuDDK0_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_130/is_best_candidate_branch.root");
    TTree* t_BuDDK0_best = (TTree*)f_BuDDK0_best->Get("DecayTree");

    Int_t N_130 = t_BuDDK0->GetEntries();
    Int_t N_130_best = t_BuDDK0_best->GetEntries();
    if(N_130 != N_130_best)
    {
        cout << "Wrong number of entries in BuDDK0 MC" << endl;
        return;
    }
    else
    {
        t_BuDDK0->AddFriend(t_BuDDK0_best);
    }

    // B+ -> DD 
    TFile* f_BuDD = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt1_0_bdt2_0.root");
    TTree* t_BuDD = (TTree*)f_BuDD->Get("DecayTree");

    TFile* f_BuDD_best = new TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_150/is_best_candidate_branch.root");
    TTree* t_BuDD_best = (TTree*)f_BuDD_best->Get("DecayTree");

    Int_t N_150 = t_BuDD->GetEntries();
    Int_t N_150_best = t_BuDD_best->GetEntries();
    if(N_150 != N_150_best)
    {
        cout << "Wrong number of entries in BuDD MC" << endl;
        return;
    }
    else
    {
        t_BuDD->AddFriend(t_BuDD_best);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Ktautau MC: TM vs best candidate selection
    TH1D* h_ktautau_mc_3pi3pi_tm = new TH1D("h_ktautau_mc_3pi3pi_tm", "h_ktautau_mc_3pi3pi_tm", 100, 4000, 8000);
    TH1D* h_ktautau_mc_3pi3pipi0_tm = new TH1D("h_ktautau_mc_3pi3pipi0_tm", "h_ktautau_mc_3pi3pipi0_tm", 100, 4000, 8000);
    TH1D* h_ktautau_mc_3pi3pi2pi0_tm = new TH1D("h_ktautau_mc_3pi3pi2pi0_tm", "h_ktautau_mc_3pi3pi2pi0_tm", 100, 4000, 8000);
    TH1D* h_ktautau_mc_all_tm = new TH1D("h_ktautau_mc_all_tm", "h_ktautau_mc_all_tm", 100, 4000, 8000);

    TH1D* h_ktautau_mc_3pi3pi = new TH1D("h_ktautau_mc_3pi3pi", "h_ktautau_mc_3pi3pi", 100, 4000, 8000);
    TH1D* h_ktautau_mc_3pi3pipi0 = new TH1D("h_ktautau_mc_3pi3pipi0", "h_ktautau_mc_3pi3pipi0", 100, 4000, 8000);
    TH1D* h_ktautau_mc_3pi3pi2pi0 = new TH1D("h_ktautau_mc_3pi3pi2pi0", "h_ktautau_mc_3pi3pi2pi0", 100, 4000, 8000);
    TH1D* h_ktautau_mc_all = new TH1D("h_ktautau_mc_all", "h_ktautau_mc_all", 100, 4000, 8000);

    t_ktautau_mc_tm->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi_tm", "component==0");
    t_ktautau_mc_tm->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0_tm", "component==1");
    t_ktautau_mc_tm->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0_tm", "component==2");
    t_ktautau_mc_tm->Draw("df_Bp_M >> h_ktautau_mc_all_tm");

    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi", "(is_best_cand == 1) && (component==0)");
    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0", "(is_best_cand == 1) && (component==1)");
    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0", "(is_best_cand == 1) && (component==2)");
    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all", "(is_best_cand == 1)");

    auto c1 = new TCanvas();
    auto c2 = new TCanvas(); 
    auto c3 = new TCanvas(); 
    auto c4 = new TCanvas(); 

    auto FourPads = new TCanvas("FourPads","FourPads");
    FourPads->Divide(2,2);
    FourPads->cd(1) ; c1->DrawClonePad();

    h_ktautau_mc_3pi3pi_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pi_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pi_tm->SetTitle("3#pi3#pi MC: 2016-2018");

    h_ktautau_mc_3pi3pi_tm->SetLineColor(kBlue);
    h_ktautau_mc_3pi3pi_tm->SetFillColorAlpha(kBlue, 0.25);

    h_ktautau_mc_3pi3pi->SetLineColor(kBlack);
    h_ktautau_mc_3pi3pi->SetFillColorAlpha(kBlack, 0.25);

    h_ktautau_mc_3pi3pi_tm->DrawNormalized();
    h_ktautau_mc_3pi3pi->DrawNormalized("same");

    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(h_ktautau_mc_3pi3pi_tm, "Truth-matched", "lf");
    leg->AddEntry(h_ktautau_mc_3pi3pi, "Best candidate", "lf");
    leg->Draw("same");

    FourPads->cd(2) ; c2->DrawClonePad();

    h_ktautau_mc_3pi3pipi0_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pipi0_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pipi0_tm->SetTitle("3#pi3#pi #pi^{0} MC: 2016-2018");

    h_ktautau_mc_3pi3pipi0_tm->SetLineColor(kBlue);
    h_ktautau_mc_3pi3pipi0_tm->SetFillColorAlpha(kBlue, 0.25);

    h_ktautau_mc_3pi3pipi0->SetLineColor(kBlack);
    h_ktautau_mc_3pi3pipi0->SetFillColorAlpha(kBlack, 0.25);

    h_ktautau_mc_3pi3pipi0_tm->DrawNormalized();
    h_ktautau_mc_3pi3pipi0->DrawNormalized("same");
    leg->Draw("same");
    
    FourPads->cd(3) ; c3->DrawClonePad();

    h_ktautau_mc_3pi3pi2pi0_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pi2pi0_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pi2pi0_tm->SetTitle("3#pi3#pi 2#pi^{0} MC: 2016-2018");

    h_ktautau_mc_3pi3pi2pi0_tm->SetLineColor(kBlue);
    h_ktautau_mc_3pi3pi2pi0_tm->SetFillColorAlpha(kBlue, 0.25);

    h_ktautau_mc_3pi3pi2pi0->SetLineColor(kBlack);
    h_ktautau_mc_3pi3pi2pi0->SetFillColorAlpha(kBlack, 0.25);

    h_ktautau_mc_3pi3pi2pi0_tm->DrawNormalized();
    h_ktautau_mc_3pi3pi2pi0->DrawNormalized("same");
    leg->Draw("same");

    FourPads->cd(4) ; c4->DrawClonePad();

    h_ktautau_mc_all_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_all_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_all_tm->SetTitle("All MC: 2016-2018");

    h_ktautau_mc_all_tm->SetLineColor(kBlue);
    h_ktautau_mc_all_tm->SetFillColorAlpha(kBlue, 0.25);

    h_ktautau_mc_all->SetLineColor(kBlack);
    h_ktautau_mc_all->SetFillColorAlpha(kBlack, 0.25);

    h_ktautau_mc_all_tm->DrawNormalized();
    h_ktautau_mc_all->DrawNormalized("same");
    leg->Draw("same");

    FourPads->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_TM_vs_best.pdf");

    TCut ktautau_truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";

    TH1D* h_ktautau_mc_3pi3pi_not_TM = new TH1D("h_ktautau_mc_3pi3pi_not_TM", "h_ktautau_mc_3pi3pi_not_TM", 100, 4000, 8000);
    TH1D* h_ktautau_mc_3pi3pipi0_not_TM = new TH1D("h_ktautau_mc_3pi3pipi0_not_TM", "h_ktautau_mc_3pi3pipi0_not_TM", 100, 4000, 8000);
    TH1D* h_ktautau_mc_3pi3pi2pi0_not_TM = new TH1D("h_ktautau_mc_3pi3pi2pi0_not_TM", "h_ktautau_mc_3pi3pi2pi0_not_TM", 100, 4000, 8000);
    TH1D* h_ktautau_mc_all_not_TM = new TH1D("h_ktautau_mc_all_not_TM", "h_ktautau_mc_all_not_TM", 100, 4000, 8000);

    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi_not_TM", "(is_best_cand == 1) && (component==0)"+!ktautau_truthMatch);
    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0_not_TM", "(is_best_cand == 1) && (component==1)"+!ktautau_truthMatch);
    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0_not_TM", "(is_best_cand == 1) && (component==2)"+!ktautau_truthMatch);
    t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all_not_TM", "(is_best_cand == 1)"+!ktautau_truthMatch);

    auto d1 = new TCanvas(); 
    auto d2 = new TCanvas(); 
    auto d3 = new TCanvas(); 
    auto d4 = new TCanvas(); 

    auto FourPads1 = new TCanvas("FourPads1","FourPads1");
    FourPads1->Divide(2,2);
    FourPads1->cd(1) ; d1->DrawClonePad();

    h_ktautau_mc_3pi3pi_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pi_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pi_not_TM->SetTitle("3#pi3#pi MC: 2016-2018");

    h_ktautau_mc_3pi3pi_not_TM->SetLineColor(kCyan+1);
    h_ktautau_mc_3pi3pi_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    h_ktautau_mc_3pi3pi_not_TM->DrawNormalized();

    TLegend* leg1 = new TLegend(0.5,0.8,0.9,0.9);
    leg1->AddEntry(h_ktautau_mc_3pi3pi_not_TM, "Best candidate (not pass TM)", "lf");
    leg1->Draw("same");

    FourPads1->cd(2) ; d2->DrawClonePad();

    h_ktautau_mc_3pi3pipi0_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pipi0_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pipi0_not_TM->SetTitle("3#pi3#pi #pi^{0} MC: 2016-2018");

    h_ktautau_mc_3pi3pipi0_not_TM->SetLineColor(kCyan+1);
    h_ktautau_mc_3pi3pipi0_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    h_ktautau_mc_3pi3pipi0_not_TM->DrawNormalized();
    leg1->Draw("same");

    FourPads1->cd(3) ; d3->DrawClonePad();

    h_ktautau_mc_3pi3pi2pi0_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pi2pi0_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pi2pi0_not_TM->SetTitle("3#pi3#pi 2#pi^{0} MC: 2016-2018");

    h_ktautau_mc_3pi3pi2pi0_not_TM->SetLineColor(kCyan+1);
    h_ktautau_mc_3pi3pi2pi0_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    h_ktautau_mc_3pi3pi2pi0_not_TM->DrawNormalized();
    leg1->Draw("same");

    FourPads1->cd(4) ; d4->DrawClonePad();

    h_ktautau_mc_all_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_all_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_all_not_TM->SetTitle("All MC: 2016-2018");

    h_ktautau_mc_all_not_TM->SetLineColor(kCyan+1);
    h_ktautau_mc_all_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    h_ktautau_mc_all_not_TM->DrawNormalized();
    leg1->Draw("same");

    FourPads1->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_best_notTM.pdf");

    // DDs MC: TM vs best candidate
    TH1D* h_dds_mc_tm = new TH1D("h_dds_mc_tm", "h_dds_mc_tm", 100, 5235, 5355);
    TH1D* h_dds_mc = new TH1D("h_dds_mc", "h_dds_mc", 100, 5235, 5355);

    t_dds_mc_tm->Draw("Bp_dtf_M[0] >> h_dds_mc_tm");
    t_dds_mc->Draw("Bp_dtf_M[0] >> h_dds_mc", "is_best_cand == 1");

    h_dds_mc_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_mc_tm->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_mc_tm->SetTitle("DDs MC");

    h_dds_mc_tm->SetLineColor(kBlue);
    h_dds_mc_tm->SetFillColorAlpha(kBlue, 0.25);

    h_dds_mc->SetLineColor(kBlack);
    h_dds_mc->SetFillColorAlpha(kBlack, 0.25);

    TCanvas c5;
    c5.cd();
    h_dds_mc_tm->DrawNormalized();
    h_dds_mc->DrawNormalized("same");
    leg->Draw("same");
    c5.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_mc_TM_vs_best.pdf");

    TCut dds_truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(Dsp_K1_TRUEID) == 321) && (abs(Dsp_K2_TRUEID) == 321) && (abs(D0bar_pi_TRUEID) == 211) && (abs(Dsp_pi_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(Dsp_TRUEID) == 431) && (abs(Bp_TRUEID) == 521)";
    
    TH1D* h_dds_mc_not_tm = new TH1D("h_dds_mc_not_tm", "h_dds_mc_not_tm", 100, 5235, 5355);

    t_dds_mc->Draw("Bp_dtf_M[0] >> h_dds_mc_not_tm", "(is_best_cand == 1)"+!dds_truthMatch);

    h_dds_mc_not_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_mc_not_tm->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_mc_not_tm->SetTitle("DDs MC");

    h_dds_mc_not_tm->SetLineColor(kCyan+1);
    h_dds_mc_not_tm->SetFillColorAlpha(kCyan+1, 0.25);

    TCanvas d5;
    d5.cd();
    h_dds_mc_not_tm->DrawNormalized();
    leg1->Draw("same");
    d5.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_mc_best_notTM.pdf");

    /*
    // Ktautau data
    TH1D* h_rs_data = new TH1D("h_rs_data", "h_rs_data", 100, 4000, 8000);
    TH1D* h_rs_data_best = new TH1D("h_rs_data_best", "h_rs_data_best", 100, 4000, 8000);

    t_ktautau_rs_data->Draw("df_Bp_M >> h_rs_data");
    // t_ktautau_rs_data->Draw("df_Bp_M >> h_rs_data_best", "(is_best_cand == 1)");

    h_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_rs_data->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_rs_data->SetTitle("RS data");

    h_rs_data->SetLineColor(kRed);
    h_rs_data->SetFillColorAlpha(kRed, 0.25);

    // h_rs_data_best->SetLineColor(kBlack);
    // h_rs_data_best->SetFillColorAlpha(kBlack, 0.25);

    TCanvas c6;
    c6.cd();
    h_rs_data->DrawNormalized();
    // h_rs_data_best->DrawNormalized("same");
    TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2->AddEntry(h_rs_data, "All candidates");
    // leg2->AddEntry(h_rs_data_best, "Best candidate");
    leg2->Draw("same");
    c6.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/rs_data_multiple_vs_best.pdf");


    TH1D* h_ws_data = new TH1D("h_ws_data", "h_ws_data", 100, 4000, 8000);
    TH1D* h_ws_data_best = new TH1D("h_ws_data_best", "h_ws_data_best", 100, 4000, 8000);

    t_ktautau_ws_data->Draw("df_Bp_M >> h_ws_data");
    // t_ktautau_ws_data->Draw("df_Bp_M >> h_ws_data_best", "(is_best_cand == 1)");

    h_ws_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ws_data->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ws_data->SetTitle("WS data");

    h_ws_data->SetLineColor(kRed);
    h_ws_data->SetFillColorAlpha(kRed, 0.25);

    // h_ws_data_best->SetLineColor(kBlack);
    // h_ws_data_best->SetFillColorAlpha(kBlack, 0.25);

    TCanvas c7;
    c7.cd();
    h_ws_data->DrawNormalized();
    // h_ws_data_best->DrawNormalized("same");
    leg2->Draw("same");
    c7.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ws_data_multiple_vs_best.pdf");
    */

    // DDs data
    TH1D* h_dds_data = new TH1D("h_dds_data", "h_dds_data", 100, 5235, 5355);
    TH1D* h_dds_data_best = new TH1D("h_dds_data_best", "h_dds_data_best", 100, 5235, 5355);

    t_dds_data->Draw("Bp_dtf_M[0] >> h_dds_data");
    t_dds_data->Draw("Bp_dtf_M[0] >> h_dds_data_best", "(is_best_cand == 1)");

    h_dds_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_data->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_data->SetTitle("DDs data");

    h_dds_data->SetLineColor(kRed);
    h_dds_data->SetFillColorAlpha(kRed, 0.25);

    h_dds_data_best->SetLineColor(kBlack);
    h_dds_data_best->SetFillColorAlpha(kBlack, 0.25);

    TCanvas c8;
    c8.cd();
    h_dds_data->DrawNormalized();
    h_dds_data_best->DrawNormalized("same");
    TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2->AddEntry(h_dds_data, "All candidates");
    leg2->AddEntry(h_dds_data_best, "Best candidate");
    leg2->Draw("same");
    c8.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_data_multiple_vs_best.pdf");

    // Cocktail MCs
    TH1D* h_BuD0D0Kp_norm_noBDT = new TH1D("h_BuD0D0Kp_norm_noBDT", "h_BuD0D0Kp_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BuDpDmKp_norm_noBDT = new TH1D("h_BuDpDmKp_norm_noBDT", "h_BuDpDmKp_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BuDsDsKp_norm_noBDT = new TH1D("h_BuDsDsKp_norm_noBDT", "h_BuDsDsKp_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BdDpD0Kp_norm_noBDT = new TH1D("h_BdDpD0Kp_norm_noBDT", "h_BdDpD0Kp_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BsDsD0Kp_norm_noBDT = new TH1D("h_BsDsD0Kp_norm_noBDT", "h_BsDsD0Kp_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BuD0DpK0_norm_noBDT = new TH1D("h_BuD0DpK0_norm_noBDT", "h_BuD0DpK0_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BuD0Ds_norm_noBDT = new TH1D("h_BuD0Ds_norm_noBDT", "h_BuD0Ds_norm_noBDT", 100, 4000, 8000);
    TH1D* h_BuD0Dp_norm_noBDT = new TH1D("h_BuD0Dp_norm_noBDT", "h_BuD0Dp _norm_noBDT", 100, 4000, 8000);

    t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_norm_noBDT", "(species == 100) && (is_best_cand == 1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_norm_noBDT", "(species == 101) && (is_best_cand == 1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_norm_noBDT", "(species == 102) && (is_best_cand == 1)");
    t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_norm_noBDT", "(species == 110) && (is_best_cand == 1)");  
    t_BsDDKp->Draw("df_Bp_M >> h_BsDsD0Kp_norm_noBDT", "(species == 120) && (is_best_cand == 1)");
    t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_norm_noBDT", "(species == 130) && (is_best_cand == 1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_norm_noBDT", "(species == 150) && (is_best_cand == 1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_norm_noBDT", "(species == 151) && (is_best_cand == 1)");

    h_ktautau_mc_3pi3pi->SetLineColor(kBlue);
    h_ktautau_mc_3pi3pi->SetFillColorAlpha(kBlue, 0.25);

    h_ktautau_mc_all->SetLineColor(kBlack);
    h_ktautau_mc_all->SetFillColorAlpha(kBlack, 0.25);

    h_BuD0D0Kp_norm_noBDT->SetLineColor(kCyan+1);
    h_BuD0D0Kp_norm_noBDT->SetFillColorAlpha(kCyan+1, 0.25);

    h_BuDpDmKp_norm_noBDT->SetLineColor(kGreen+1);
    h_BuDpDmKp_norm_noBDT->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuDsDsKp_norm_noBDT->SetLineColor(kViolet+1);
    h_BuDsDsKp_norm_noBDT->SetFillColorAlpha(kViolet+1, 0.25);

    h_BdDpD0Kp_norm_noBDT->SetLineColor(kViolet+2);
    h_BdDpD0Kp_norm_noBDT->SetFillColorAlpha(kViolet+2, 0.25);

    h_BsDsD0Kp_norm_noBDT->SetLineColor(kPink+7);    
    h_BsDsD0Kp_norm_noBDT->SetFillColorAlpha(kPink+7, 0.25);

    h_BuD0DpK0_norm_noBDT->SetLineColor(kRed);
    h_BuD0DpK0_norm_noBDT->SetFillColorAlpha(kRed, 0.25);

    h_BuD0Ds_norm_noBDT->SetLineColor(kOrange-3);
    h_BuD0Ds_norm_noBDT->SetFillColorAlpha(kOrange-3, 0.25);

    h_BuD0Dp_norm_noBDT->SetLineColor(kOrange+4);
    h_BuD0Dp_norm_noBDT->SetFillColorAlpha(kOrange+4, 0.25);

    TCanvas d6;
    d6.cd();
    h_ktautau_mc_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_ktautau_mc_3pi3pi->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_ktautau_mc_3pi3pi->SetTitle("");
    h_ktautau_mc_3pi3pi->DrawNormalized();
    h_ktautau_mc_all->DrawNormalized("same");
    h_BuD0D0Kp_norm_noBDT->DrawNormalized("same");
    h_BuDpDmKp_norm_noBDT->DrawNormalized("same");
    h_BuDsDsKp_norm_noBDT->DrawNormalized("same");
    h_BdDpD0Kp_norm_noBDT->DrawNormalized("same");
    h_BsDsD0Kp_norm_noBDT->DrawNormalized("same");
    h_BuD0DpK0_norm_noBDT->DrawNormalized("same");
    h_BuD0Ds_norm_noBDT->DrawNormalized("same");
    h_BuD0Dp_norm_noBDT->DrawNormalized("same");
    TLegend* leg3 = new TLegend(0.65,0.3,0.9,0.9);
    leg3->AddEntry(h_ktautau_mc_3pi3pi, "3#pi3#pi MC");
    leg3->AddEntry(h_ktautau_mc_all, "All MC");
    leg3->AddEntry(h_BuD0D0Kp_norm_noBDT, "B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+}");
    leg3->AddEntry(h_BuDpDmKp_norm_noBDT, "B^{+} #rightarrow D^{+} D^{-} K^{+}");
    leg3->AddEntry(h_BuDsDsKp_norm_noBDT, "B^{+} #rightarrow D^{+}_{s} D^{+}_{s} K^{+}");
    leg3->AddEntry(h_BdDpD0Kp_norm_noBDT, "B^{0} #rightarrow D^{-} D^{0} K^{+}");
    leg3->AddEntry(h_BsDsD0Kp_norm_noBDT, "B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+}");
    leg3->AddEntry(h_BuD0DpK0_norm_noBDT, "B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0}");
    leg3->AddEntry(h_BuD0Ds_norm_noBDT, "B^{+} #rightarrow  #bar{D}^{0} D^{+}_{s}");
    leg3->AddEntry(h_BuD0Dp_norm_noBDT, "B^{+} #rightarrow  #bar{D}^{0} D^{+}");
    leg3->Draw("same");
    d6.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_best_bmass.pdf");

    TH1D* h_3pi3pi_norm_BDTcut = new TH1D("h_3pi3pi_norm_BDTcut", "h_3pi3pi_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_all_mc_norm_BDTcut = new TH1D("h_all_mc_norm_BDTcut", "h_all_mc_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BuD0D0Kp_norm_BDTcut = new TH1D("h_BuD0D0Kp_norm_BDTcut", "h_BuD0D0Kp_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BuDpDmKp_norm_BDTcut = new TH1D("h_BuDpDmKp_norm_BDTcut", "h_BuDpDmKp_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BuDsDsKp_norm_BDTcut = new TH1D("h_BuDsDsKp_norm_BDTcut", "h_BuDsDsKp_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BdDpD0Kp_norm_BDTcut = new TH1D("h_BdDpD0Kp_norm_BDTcut", "h_BdDpD0Kp_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BsDsD0Kp_norm_BDTcut = new TH1D("h_BsDsD0Kp_norm_BDTcut", "h_BsDsD0Kp_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BuD0DpK0_norm_BDTcut = new TH1D("h_BuD0DpK0_norm_BDTcut", "h_BuD0DpK0_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BuD0Ds_norm_BDTcut = new TH1D("h_BuD0Ds_norm_BDTcut", "h_BuD0Ds_norm_BDTcut", 100, 4000, 8000);
    TH1D* h_BuD0Dp_norm_BDTcut = new TH1D("h_BuD0Dp_norm_BDTcut", "h_BuD0Dp _norm_BDTcut", 100, 4000, 8000);

    t_ktautau_mc->Draw("df_Bp_M >> h_3pi3pi_norm_BDTcut", "(component==0) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_ktautau_mc->Draw("df_Bp_M >> h_all_mc_norm_BDTcut", "(is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_norm_BDTcut", "(species == 100) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_norm_BDTcut", "(species == 101) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_norm_BDTcut", "(species == 102) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_norm_BDTcut", "(species == 110) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");  
    t_BsDDKp->Draw("df_Bp_M >> h_BsDsD0Kp_norm_BDTcut", "(species == 120) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_norm_BDTcut", "(species == 130) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_norm_BDTcut", "(species == 150) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_norm_BDTcut", "(species == 151) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)");

    h_3pi3pi_norm_BDTcut->SetLineColor(kBlue);
    h_3pi3pi_norm_BDTcut->SetFillColorAlpha(kBlue, 0.25);

    h_all_mc_norm_BDTcut->SetLineColor(kBlack);
    h_all_mc_norm_BDTcut->SetFillColorAlpha(kBlack, 0.25);

    h_BuD0D0Kp_norm_BDTcut->SetLineColor(kCyan+1);
    h_BuD0D0Kp_norm_BDTcut->SetFillColorAlpha(kCyan+1, 0.25);

    h_BuDpDmKp_norm_BDTcut->SetLineColor(kGreen+1);
    h_BuDpDmKp_norm_BDTcut->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuDsDsKp_norm_BDTcut->SetLineColor(kViolet+1);
    h_BuDsDsKp_norm_BDTcut->SetFillColorAlpha(kViolet+1, 0.25);

    h_BdDpD0Kp_norm_BDTcut->SetLineColor(kViolet+2);
    h_BdDpD0Kp_norm_BDTcut->SetFillColorAlpha(kViolet+2, 0.25);

    h_BsDsD0Kp_norm_BDTcut->SetLineColor(kPink+7);    
    h_BsDsD0Kp_norm_BDTcut->SetFillColorAlpha(kPink+7, 0.25);

    h_BuD0DpK0_norm_BDTcut->SetLineColor(kRed);
    h_BuD0DpK0_norm_BDTcut->SetFillColorAlpha(kRed, 0.25);

    h_BuD0Ds_norm_BDTcut->SetLineColor(kOrange-3);
    h_BuD0Ds_norm_BDTcut->SetFillColorAlpha(kOrange-3, 0.25);

    h_BuD0Dp_norm_BDTcut->SetLineColor(kOrange+4);
    h_BuD0Dp_norm_BDTcut->SetFillColorAlpha(kOrange+4, 0.25);

    TCanvas d7;
    d7.cd();
    h_3pi3pi_norm_BDTcut->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_3pi3pi_norm_BDTcut->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_3pi3pi_norm_BDTcut->SetTitle("");
    h_3pi3pi_norm_BDTcut->DrawNormalized();
    h_all_mc_norm_BDTcut->DrawNormalized("same");
    h_BuD0D0Kp_norm_BDTcut->DrawNormalized("same");
    h_BuDpDmKp_norm_BDTcut->DrawNormalized("same");
    h_BuDsDsKp_norm_BDTcut->DrawNormalized("same");
    h_BdDpD0Kp_norm_BDTcut->DrawNormalized("same");
    h_BsDsD0Kp_norm_BDTcut->DrawNormalized("same");
    h_BuD0DpK0_norm_BDTcut->DrawNormalized("same");
    h_BuD0Ds_norm_BDTcut->DrawNormalized("same");
    h_BuD0Dp_norm_BDTcut->DrawNormalized("same");
    leg3->Draw("same");
    d7.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_best_bmass_BDTcut.pdf");

    // Ktautau MC: 3 MC components
    TCanvas d8;
    d8.cd();
    h_ktautau_mc_3pi3pi->SetTitle("B^{+} mass");
    h_ktautau_mc_3pi3pipi0->SetLineColor(kOrange-3);
    h_ktautau_mc_3pi3pipi0->SetFillColorAlpha(kOrange-3, 0.25);

    h_ktautau_mc_3pi3pi2pi0->SetLineColor(kRed);
    h_ktautau_mc_3pi3pi2pi0->SetFillColorAlpha(kRed, 0.25);

    h_ktautau_mc_3pi3pi->DrawNormalized();
    h_ktautau_mc_3pi3pipi0->DrawNormalized("same");
    h_ktautau_mc_3pi3pi2pi0->DrawNormalized("same");

    TLegend* leg4 = new TLegend(0.7,0.7,0.9,0.9);
    leg4->AddEntry(h_ktautau_mc_3pi3pi, "3#pi3#pi MC");
    leg4->AddEntry(h_ktautau_mc_3pi3pipi0, "3#pi3#pi #pi^{0} MC");
    leg4->AddEntry(h_ktautau_mc_3pi3pi2pi0, "3#pi3#pi 2#pi^{0} MC");
    leg4->Draw("same");
    d8.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_comp.pdf");

    TH1D* h_3pi3pi_iso = new TH1D("h_3pi3pi_iso", "h_3pi3pi_iso", 100, 0, 1);
    TH1D* h_3pi3pipi0_iso = new TH1D("h_3pi3pipi0_iso", "h_3pi3pipi0_iso", 100, 0, 1);
    TH1D* h_3pi3pi2pi0_iso = new TH1D("h_3pi3pi2pi0_iso", "h_3pi3pi2pi0_iso", 100, 0, 1);

    t_ktautau_mc->Draw("BDT1 >> h_3pi3pi_iso", "(component==0) && (is_best_cand == 1)");
    t_ktautau_mc->Draw("BDT1 >> h_3pi3pipi0_iso", "(component==1) && (is_best_cand == 1)");
    t_ktautau_mc->Draw("BDT1 >> h_3pi3pi2pi0_iso", "(component==2) && (is_best_cand == 1)");

    h_3pi3pi_iso->GetXaxis()->SetTitle("BDT1");
    h_3pi3pi_iso->GetYaxis()->SetTitle("Normalized entries / (0.01)");
    h_3pi3pi_iso->SetTitle("Isolation MVA");

    h_3pi3pi_iso->SetLineColor(kBlue);
    h_3pi3pi_iso->SetFillColorAlpha(kBlue, 0.25);
    
    h_3pi3pipi0_iso->SetLineColor(kOrange-3);
    h_3pi3pipi0_iso->SetFillColorAlpha(kOrange-3, 0.25);

    h_3pi3pi2pi0_iso->SetLineColor(kRed);
    h_3pi3pi2pi0_iso->SetFillColorAlpha(kRed, 0.25);
    
    TCanvas d9;
    d9.cd();
    h_3pi3pi_iso->DrawNormalized();
    h_3pi3pipi0_iso->DrawNormalized("same");
    h_3pi3pi2pi0_iso->DrawNormalized("same");
    TLegend* leg5 = new TLegend(0.4,0.7,0.6,0.9);
    leg5->AddEntry(h_3pi3pi_iso, "3#pi3#pi MC");
    leg5->AddEntry(h_3pi3pipi0_iso, "3#pi3#pi #pi^{0} MC");
    leg5->AddEntry(h_3pi3pi2pi0_iso, "3#pi3#pi 2#pi^{0} MC");
    leg5->Draw("same");
    d9.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_comp_bdt1.pdf");

    TH1D* h_3pi3pi_kin = new TH1D("h_3pi3pi_kin", "h_3pi3pi_kin", 100, 0, 1);
    TH1D* h_3pi3pipi0_kin = new TH1D("h_3pi3pipi0_kin", "h_3pi3pipi0_kin", 100, 0, 1);
    TH1D* h_3pi3pi2pi0_kin = new TH1D("h_3pi3pi2pi0_kin", "h_3pi3pi2pi0_kin", 100, 0, 1);

    t_ktautau_mc->Draw("BDT2 >> h_3pi3pi_kin", "(component==0) && (is_best_cand == 1)");
    t_ktautau_mc->Draw("BDT2 >> h_3pi3pipi0_kin", "(component==1) && (is_best_cand == 1)");
    t_ktautau_mc->Draw("BDT2 >> h_3pi3pi2pi0_kin", "(component==2) && (is_best_cand == 1)");

    h_3pi3pi_kin->GetXaxis()->SetTitle("BDT2");
    h_3pi3pi_kin->GetYaxis()->SetTitle("Normalized entries / (0.01)");
    h_3pi3pi_kin->SetTitle("Topological MVA");

    h_3pi3pi_kin->SetLineColor(kBlue);
    h_3pi3pi_kin->SetFillColorAlpha(kBlue, 0.25);
    
    h_3pi3pipi0_kin->SetLineColor(kOrange-3);
    h_3pi3pipi0_kin->SetFillColorAlpha(kOrange-3, 0.25);

    h_3pi3pi2pi0_kin->SetLineColor(kRed);
    h_3pi3pi2pi0_kin->SetFillColorAlpha(kRed, 0.25);
    
    TCanvas d10;
    d10.cd();
    h_3pi3pi_kin->DrawNormalized();
    h_3pi3pipi0_kin->DrawNormalized("same");
    h_3pi3pi2pi0_kin->DrawNormalized("same");
    leg5->Draw("same");
    d10.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_comp_bdt2.pdf");

    // Cocktail MCs
    auto e1 = new TCanvas(); 
    auto e2 = new TCanvas(); 
    auto e3 = new TCanvas(); 
    auto e4 = new TCanvas(); 
    auto e5 = new TCanvas(); 
    auto e6 = new TCanvas();
    auto e7 = new TCanvas(); 
    auto e8 = new TCanvas(); 

    auto EightPads = new TCanvas("EightPads","EightPads");
    EightPads->Divide(2,4);
    EightPads->cd(1) ; e1->DrawClonePad();

    TH1D* h_BuD0D0Kp_DD = new TH1D("h_BuD0D0Kp_DD", "h_BuD0D0Kp_DD", 100, 4000, 8000);
    TH1D* h_BuD0D0Kp_DsD = new TH1D("h_BuD0D0Kp_DsD", "h_BuD0D0Kp_DsD", 100, 4000, 8000);
    TH1D* h_BuD0D0Kp_DDs = new TH1D("h_BuD0D0Kp_DDs", "h_BuD0D0Kp_DDs", 100, 4000, 8000);
    TH1D* h_BuD0D0Kp_DsDs = new TH1D("h_BuD0D0Kp_DsDs", "h_BuD0D0Kp_DsDs", 100, 4000, 8000);

    t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DD", "(species==100) && (component==0) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DsD", "(species==100) && (component==1) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DDs", "(species==100) && (component==2) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DsDs", "(species==100) && (component==3) && (is_best_cand==1)");

    h_BuD0D0Kp_DsDs->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BuD0D0Kp_DsDs->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BuD0D0Kp_DsDs->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+} MC");

    h_BuD0D0Kp_DD->SetLineColor(kBlue);
    h_BuD0D0Kp_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BuD0D0Kp_DsD->SetLineColor(kBlack);
    h_BuD0D0Kp_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BuD0D0Kp_DDs->SetLineColor(kRed);
    h_BuD0D0Kp_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BuD0D0Kp_DsDs->SetLineColor(kGreen+1);
    h_BuD0D0Kp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuD0D0Kp_DsDs->DrawNormalized();
    h_BuD0D0Kp_DsD->DrawNormalized("same");
    h_BuD0D0Kp_DDs->DrawNormalized("same");
    h_BuD0D0Kp_DD->DrawNormalized("same");

    TLegend* leg6 = new TLegend(0.8,0.6,0.9,0.9);
    leg6->AddEntry(h_BuD0D0Kp_DD, "DD");
    leg6->AddEntry(h_BuD0D0Kp_DsD, "Dstar D");
    leg6->AddEntry(h_BuD0D0Kp_DDs, "D Dstar");
    leg6->AddEntry(h_BuD0D0Kp_DsDs, "Dstar Dstar");
    leg6->Draw("same");

    EightPads->cd(2) ; e2->DrawClonePad();

    TH1D* h_BuDpDmKp_DD = new TH1D("h_BuDpDmKp_DD", "h_BuDpDmKp_DD", 100, 4000, 8000);
    TH1D* h_BuDpDmKp_DsD = new TH1D("h_BuDpDmKp_DsD", "h_BuDpDmKp_DsD", 100, 4000, 8000);
    TH1D* h_BuDpDmKp_DDs = new TH1D("h_BuDpDmKp_DDs", "h_BuDpDmKp_DDs", 100, 4000, 8000);
    TH1D* h_BuDpDmKp_DsDs = new TH1D("h_BuDpDmKp_DsDs", "h_BuDpDmKp_DsDs", 100, 4000, 8000);

    t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DD", "(species==101) && (component==0) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DsD", "(species==101) && (component==1) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DDs", "(species==101) && (component==2) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DsDs", "(species==101) && (component==3) && (is_best_cand==1)");

    h_BuDpDmKp_DsD->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BuDpDmKp_DsD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BuDpDmKp_DsD->SetTitle("B^{+} #rightarrow D^{+} D^{-} K^{+} MC");

    h_BuDpDmKp_DD->SetLineColor(kBlue);
    h_BuDpDmKp_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BuDpDmKp_DsD->SetLineColor(kBlack);
    h_BuDpDmKp_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BuDpDmKp_DDs->SetLineColor(kRed);
    h_BuDpDmKp_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BuDpDmKp_DsDs->SetLineColor(kGreen+1);
    h_BuDpDmKp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuDpDmKp_DsD->DrawNormalized();
    h_BuDpDmKp_DD->DrawNormalized("same");
    h_BuDpDmKp_DDs->DrawNormalized("same");
    h_BuDpDmKp_DsDs->DrawNormalized("same");
    leg6->Draw("same");
    
    EightPads->cd(3) ; e3->DrawClonePad();

    TH1D* h_BuDsDsKp_DD = new TH1D("h_BuDsDsKp_DD", "h_BuDsDsKp_DD", 100, 4000, 8000);
    TH1D* h_BuDsDsKp_DsD = new TH1D("h_BuDsDsKp_DsD", "h_BuDsDsKp_DsD", 100, 4000, 8000);
    TH1D* h_BuDsDsKp_DDs = new TH1D("h_BuDsDsKp_DDs", "h_BuDsDsKp_DDs", 100, 4000, 8000);
    TH1D* h_BuDsDsKp_DsDs = new TH1D("h_BuDsDsKp_DsDs", "h_BuDsDsKp_DsDs", 100, 4000, 8000);

    t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DD", "(species==102) && (component==0) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DsD", "(species==102) && (component==1) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DDs", "(species==102) && (component==2) && (is_best_cand==1)");
    t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DsDs", "(species==102) && (component==3) && (is_best_cand==1)");

    h_BuDsDsKp_DDs->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BuDsDsKp_DDs->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BuDsDsKp_DDs->SetTitle("B^{+} #rightarrow D^{+}_{s} D^{-}_{s} K^{+} MC");

    h_BuDsDsKp_DD->SetLineColor(kBlue);
    h_BuDsDsKp_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BuDsDsKp_DsD->SetLineColor(kBlack);
    h_BuDsDsKp_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BuDsDsKp_DDs->SetLineColor(kRed);
    h_BuDsDsKp_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BuDsDsKp_DsDs->SetLineColor(kGreen+1);
    h_BuDsDsKp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuDsDsKp_DDs->DrawNormalized();
    h_BuDsDsKp_DsD->DrawNormalized("same");
    h_BuDsDsKp_DD->DrawNormalized("same");
    h_BuDsDsKp_DsDs->DrawNormalized("same");
    leg6->Draw("same");

    EightPads->cd(4) ; e4->DrawClonePad();

    TH1D* h_BdDpD0Kp_DD = new TH1D("h_BdDpD0Kp_DD", "h_BdDpD0Kp_DD", 100, 4000, 8000);
    TH1D* h_BdDpD0Kp_DsD = new TH1D("h_BdDpD0Kp_DsD", "h_BdDpD0Kp_DsD", 100, 4000, 8000);
    TH1D* h_BdDpD0Kp_DDs = new TH1D("h_BdDpD0Kp_DDs", "h_BdDpD0Kp_DDs", 100, 4000, 8000);
    TH1D* h_BdDpD0Kp_DsDs = new TH1D("h_BdDpD0Kp_DsDs", "h_BdDpD0Kp_DsDs", 100, 4000, 8000);

    t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DD", "(species==110) && (component==0) && (is_best_cand==1)");
    t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DsD", "(species==110) && (component==1) && (is_best_cand==1)");
    t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DDs", "(species==110) && (component==2) && (is_best_cand==1)");
    t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DsDs", "(species==110) && (component==3) && (is_best_cand==1)");

    h_BdDpD0Kp_DD->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BdDpD0Kp_DD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BdDpD0Kp_DD->SetTitle("B^{0} #rightarrow D^{-} D^{0} K^{+} MC");

    h_BdDpD0Kp_DD->SetLineColor(kBlue);
    h_BdDpD0Kp_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BdDpD0Kp_DsD->SetLineColor(kBlack);
    h_BdDpD0Kp_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BdDpD0Kp_DDs->SetLineColor(kRed);
    h_BdDpD0Kp_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BdDpD0Kp_DsDs->SetLineColor(kGreen+1);
    h_BdDpD0Kp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BdDpD0Kp_DD->DrawNormalized();
    h_BdDpD0Kp_DsD->DrawNormalized("same");
    h_BdDpD0Kp_DDs->DrawNormalized("same");
    h_BdDpD0Kp_DsDs->DrawNormalized("same");
    leg6->Draw("same");

    EightPads->cd(5) ; e5->DrawClonePad();

    TH1D* h_BsDpD0Kp_DD = new TH1D("h_BsDpD0Kp_DD", "h_BsDpD0Kp_DD", 100, 4000, 8000);
    TH1D* h_BsDpD0Kp_DsD = new TH1D("h_BsDpD0Kp_DsD", "h_BsDpD0Kp_DsD", 100, 4000, 8000);
    TH1D* h_BsDpD0Kp_DDs = new TH1D("h_BsDpD0Kp_DDs", "h_BsDpD0Kp_DDs", 100, 4000, 8000);
    TH1D* h_BsDpD0Kp_DsDs = new TH1D("h_BsDpD0Kp_DsDs", "h_BsDpD0Kp_DsDs", 100, 4000, 8000);

    t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DD", "(species==120) && (component==0) && (is_best_cand==1)");
    t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DsD", "(species==120) && (component==1) && (is_best_cand==1)");
    t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DDs", "(species==120) && (component==2) && (is_best_cand==1)");
    t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DsDs", "(species==120) && (component==3) && (is_best_cand==1)");

    h_BsDpD0Kp_DsDs->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BsDpD0Kp_DsDs->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BsDpD0Kp_DsDs->SetTitle("B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+} MC");

    h_BsDpD0Kp_DD->SetLineColor(kBlue);
    h_BsDpD0Kp_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BsDpD0Kp_DsD->SetLineColor(kBlack);
    h_BsDpD0Kp_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BsDpD0Kp_DDs->SetLineColor(kRed);
    h_BsDpD0Kp_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BsDpD0Kp_DsDs->SetLineColor(kGreen+1);
    h_BsDpD0Kp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BsDpD0Kp_DsDs->DrawNormalized();
    h_BsDpD0Kp_DsD->DrawNormalized("same");
    h_BsDpD0Kp_DDs->DrawNormalized("same");
    h_BsDpD0Kp_DD->DrawNormalized("same");
    leg6->Draw("same");

    EightPads->cd(6) ; e6->DrawClonePad();

    TH1D* h_BuD0DpK0_DD = new TH1D("h_BuD0DpK0_DD", "h_BuD0DpK0_DD", 100, 4000, 8000);
    TH1D* h_BuD0DpK0_DsD = new TH1D("h_BuD0DpK0_DsD", "h_BuD0DpK0_DsD", 100, 4000, 8000);
    TH1D* h_BuD0DpK0_DDs = new TH1D("h_BuD0DpK0_DDs", "h_BuD0DpK0_DDs", 100, 4000, 8000);
    TH1D* h_BuD0DpK0_DsDs = new TH1D("h_BuD0DpK0_DsDs", "h_BuD0DpK0_DsDs", 100, 4000, 8000);

    t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DD", "(species==130) && (component==0) && (is_best_cand==1)");
    t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DsD", "(species==130) && (component==1) && (is_best_cand==1)");
    t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DDs", "(species==130) && (component==2) && (is_best_cand==1)");
    t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DsDs", "(species==130) && (component==3) && (is_best_cand==1)");

    h_BuD0DpK0_DD->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BuD0DpK0_DD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BuD0DpK0_DD->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0} MC");

    h_BuD0DpK0_DD->SetLineColor(kBlue);
    h_BuD0DpK0_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BuD0DpK0_DsD->SetLineColor(kBlack);
    h_BuD0DpK0_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BuD0DpK0_DDs->SetLineColor(kRed);
    h_BuD0DpK0_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BuD0DpK0_DsDs->SetLineColor(kGreen+1);
    h_BuD0DpK0_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuD0DpK0_DD->DrawNormalized();
    h_BuD0DpK0_DsD->DrawNormalized("same");
    h_BuD0DpK0_DDs->DrawNormalized("same");
    h_BuD0DpK0_DsDs->DrawNormalized("same");
    leg6->Draw("same");

    EightPads->cd(7) ; e7->DrawClonePad();

    TH1D* h_BuD0Ds_DD = new TH1D("h_BuD0Ds_DD", "h_BuD0Ds_DD", 100, 4000, 8000);
    TH1D* h_BuD0Ds_DsD = new TH1D("h_BuD0Ds_DsD", "h_BuD0Ds_DsD", 100, 4000, 8000);
    TH1D* h_BuD0Ds_DDs = new TH1D("h_BuD0Ds_DDs", "h_BuD0Ds_DDs", 100, 4000, 8000);
    TH1D* h_BuD0Ds_DsDs = new TH1D("h_BuD0Ds_DsDs", "h_BuD0Ds_DsDs", 100, 4000, 8000);

    t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DD", "(species==150) && (component==0) && (is_best_cand==1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DsD", "(species==150) && (component==1) && (is_best_cand==1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DDs", "(species==150) && (component==2) && (is_best_cand==1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DsDs", "(species==150) && (component==3) && (is_best_cand==1)");

    h_BuD0Ds_DsD->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BuD0Ds_DsD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BuD0Ds_DsD->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{+}_{s} MC");

    h_BuD0Ds_DD->SetLineColor(kBlue);
    h_BuD0Ds_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BuD0Ds_DsD->SetLineColor(kBlack);
    h_BuD0Ds_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BuD0Ds_DDs->SetLineColor(kRed);
    h_BuD0Ds_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BuD0Ds_DsDs->SetLineColor(kGreen+1);
    h_BuD0Ds_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuD0Ds_DsD->DrawNormalized();
    h_BuD0Ds_DD->DrawNormalized("same");
    h_BuD0Ds_DDs->DrawNormalized("same");
    h_BuD0Ds_DsDs->DrawNormalized("same");
    leg6->Draw("same");

    EightPads->cd(8) ; e8->DrawClonePad();

    TH1D* h_BuD0Dp_DD = new TH1D("h_BuD0Dp_DD", "h_BuD0Dp_DD", 100, 4000, 8000);
    TH1D* h_BuD0Dp_DsD = new TH1D("h_BuD0Dp_DsD", "h_BuD0Dp_DsD", 100, 4000, 8000);
    TH1D* h_BuD0Dp_DDs = new TH1D("h_BuD0Dp_DDs", "h_BuD0Dp_DDs", 100, 4000, 8000);
    TH1D* h_BuD0Dp_DsDs = new TH1D("h_BuD0Dp_DsDs", "h_BuD0Dp_DsDs", 100, 4000, 8000);

    t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DD", "(species==151) && (component==0) && (is_best_cand==1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DsD", "(species==151) && (component==1) && (is_best_cand==1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DDs", "(species==151) && (component==2) && (is_best_cand==1)");
    t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DsDs", "(species==151) && (component==3) && (is_best_cand==1)");

    h_BuD0Dp_DD->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_BuD0Dp_DD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    h_BuD0Dp_DD->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{+} MC");

    h_BuD0Dp_DD->SetLineColor(kBlue);
    h_BuD0Dp_DD->SetFillColorAlpha(kBlue, 0.25);

    h_BuD0Dp_DsD->SetLineColor(kBlack);
    h_BuD0Dp_DsD->SetFillColorAlpha(kBlack, 0.25);

    h_BuD0Dp_DDs->SetLineColor(kRed);
    h_BuD0Dp_DDs->SetFillColorAlpha(kRed, 0.25);

    h_BuD0Dp_DsDs->SetLineColor(kGreen+1);
    h_BuD0Dp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuD0Dp_DD->DrawNormalized();
    h_BuD0Dp_DsD->DrawNormalized("same");
    h_BuD0Dp_DDs->DrawNormalized("same");
    h_BuD0Dp_DsDs->DrawNormalized("same");
    leg6->Draw("same");

    EightPads->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_DD_components.pdf");

    auto b1 = new TCanvas(); 
    auto b2 = new TCanvas(); 
    auto b3 = new TCanvas(); 
    auto b4 = new TCanvas(); 
    auto b5 = new TCanvas(); 

    auto SixPads = new TCanvas("SixPads","SixPads");
    SixPads->Divide(2,3);
    SixPads->cd(1) ; b1->DrawClonePad();

    TH1D* h_BuD0D0Kp_iso = new TH1D("h_BuD0D0Kp_iso", "h_BuD0D0Kp_iso", 100, 0, 1);
    TH1D* h_BuDpDmKp_iso = new TH1D("h_BuDpDmKp_iso", "h_BuDpDmKp_iso", 100, 0, 1);
    TH1D* h_BuDsDsKp_iso = new TH1D("h_BuDsDsKp_iso", "h_BuDsDsKp_iso", 100, 0, 1);

    t_BuDDKp->Draw("BDT1 >> h_BuD0D0Kp_iso", "(species==100) && (is_best_cand==1)");
    t_BuDDKp->Draw("BDT1 >> h_BuDpDmKp_iso", "(species==101) && (is_best_cand==1)");
    t_BuDDKp->Draw("BDT1 >> h_BuDsDsKp_iso", "(species==102) && (is_best_cand==1)");

    h_BuDpDmKp_iso->GetXaxis()->SetTitle("BDT1");
    h_BuDpDmKp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BuDpDmKp_iso->SetTitle("B^{+} #rightarrow D D K^{+} MC");

    h_BuD0D0Kp_iso->SetLineColor(kCyan+1);
    h_BuD0D0Kp_iso->SetFillColorAlpha(kCyan+1, 0.25);

    h_BuDpDmKp_iso->SetLineColor(kGreen+1);
    h_BuDpDmKp_iso->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuDsDsKp_iso->SetLineColor(kViolet+1);
    h_BuDsDsKp_iso->SetFillColorAlpha(kViolet+1, 0.25);

    h_BuDpDmKp_iso->DrawNormalized();
    h_BuD0D0Kp_iso->DrawNormalized("same");
    h_BuDsDsKp_iso->DrawNormalized("same");

    TLegend* leg7 = new TLegend(0.35,0.6,0.65,0.9);
    leg7->AddEntry(h_BuD0D0Kp_iso, "B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+}");
    leg7->AddEntry(h_BuDpDmKp_iso, "B^{+} #rightarrow D^{+} D^{-} K^{+}");
    leg7->AddEntry(h_BuDsDsKp_iso, "B^{+} #rightarrow D^{+}_{s} D^{-}_{s} K^{+}");
    leg7->Draw("same");

    SixPads->cd(2) ; b2->DrawClonePad();

    TH1D* h_BdDpD0Kp_iso = new TH1D("h_BdDpD0Kp_iso", "h_BdDpD0Kp_iso", 100, 0, 1);

    t_BdDDKp->Draw("BDT1 >> h_BdDpD0Kp_iso", "(species==110) && (is_best_cand==1)");

    h_BdDpD0Kp_iso->GetXaxis()->SetTitle("BDT1");
    h_BdDpD0Kp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BdDpD0Kp_iso->SetTitle("B^{0} #rightarrow D D K^{+} MC");

    h_BdDpD0Kp_iso->SetLineColor(kViolet+2);
    h_BdDpD0Kp_iso->SetFillColorAlpha(kViolet+2, 0.25);

    h_BdDpD0Kp_iso->DrawNormalized();

    TLegend* leg8 = new TLegend(0.35,0.8,0.65,0.9);
    leg8->AddEntry(h_BdDpD0Kp_iso, "B^{0} #rightarrow D^{-} D^{0} K^{+}");
    leg8->Draw("same");

    SixPads->cd(3) ; b3->DrawClonePad();

    TH1D* h_BsDsD0Kp_iso = new TH1D("h_BsDsD0Kp_iso", "h_BsDsD0Kp_iso", 100, 0, 1);

    t_BsDDKp->Draw("BDT1 >> h_BsDsD0Kp_iso", "(species==120) && (is_best_cand==1)");

    h_BsDsD0Kp_iso->GetXaxis()->SetTitle("BDT1");
    h_BsDsD0Kp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BsDsD0Kp_iso->SetTitle("B^{0}_{s} #rightarrow D D K^{+} MC");

    h_BsDsD0Kp_iso->SetLineColor(kPink+7);
    h_BsDsD0Kp_iso->SetFillColorAlpha(kPink+7, 0.25);

    h_BsDsD0Kp_iso->DrawNormalized();

    TLegend* leg9 = new TLegend(0.35,0.8,0.65,0.9);
    leg9->AddEntry(h_BsDsD0Kp_iso, "B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+}");
    leg9->Draw("same");

    SixPads->cd(4) ; b4->DrawClonePad();

    TH1D* h_BuD0DpK0_iso = new TH1D("h_BuD0DpK0_iso", "h_BuD0DpK0_iso", 100, 0, 1);

    t_BuDDK0->Draw("BDT1 >> h_BuD0DpK0_iso", "(species==130) && (is_best_cand==1)");

    h_BuD0DpK0_iso->GetXaxis()->SetTitle("BDT1");
    h_BuD0DpK0_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BuD0DpK0_iso->SetTitle("B^{+} #rightarrow D D K^{0} MC");

    h_BuD0DpK0_iso->SetLineColor(kRed);
    h_BuD0DpK0_iso->SetFillColorAlpha(kRed, 0.25);

    h_BuD0DpK0_iso->DrawNormalized();

    TLegend* leg10 = new TLegend(0.35,0.8,0.65,0.9);
    leg10->AddEntry(h_BuD0DpK0_iso, "B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0}");
    leg10->Draw("same");

    SixPads->cd(5) ; b5->DrawClonePad();

    TH1D* h_BuD0Ds_iso = new TH1D("h_BuD0Ds_iso", "h_BuD0Ds_iso", 100, 0, 1);
    TH1D* h_BuD0Dp_iso = new TH1D("h_BuD0Dp_iso", "h_BuD0Dp_iso", 100, 0, 1);

    t_BuDD->Draw("BDT1 >> h_BuD0Ds_iso", "(species==150) && (is_best_cand==1)");
    t_BuDD->Draw("BDT1 >> h_BuD0Dp_iso", "(species==151) && (is_best_cand==1)");

    h_BuD0Dp_iso->GetXaxis()->SetTitle("BDT1");
    h_BuD0Dp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BuD0Dp_iso->SetTitle("B^{+} #rightarrow D D MC");

    h_BuD0Ds_iso->SetLineColor(kOrange-3);
    h_BuD0Ds_iso->SetFillColorAlpha(kOrange-3, 0.25);

    h_BuD0Dp_iso->SetLineColor(kOrange+3);
    h_BuD0Dp_iso->SetFillColorAlpha(kOrange+3, 0.25);

    h_BuD0Dp_iso->DrawNormalized();
    h_BuD0Ds_iso->DrawNormalized("same");

    TLegend* leg11 = new TLegend(0.35,0.7,0.65,0.9);
    leg11->AddEntry(h_BuD0Ds_iso, "B^{+} #rightarrow #bar{D}^{0} D^{+}_{s}");
    leg11->AddEntry(h_BuD0Dp_iso, "B^{+} #rightarrow #bar{D}^{0} D^{+}");
    leg11->Draw("same");
    SixPads->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_BDT1.pdf");






































    auto a1 = new TCanvas(); 
    auto a2 = new TCanvas(); 
    auto a3 = new TCanvas(); 
    auto a4 = new TCanvas(); 
    auto a5 = new TCanvas(); 

    auto SixPads1 = new TCanvas("SixPads1","SixPads1");
    SixPads1->Divide(2,3);
    SixPads1->cd(1) ; a1->DrawClonePad();

    TH1D* h_BuD0D0Kp_kin = new TH1D("h_BuD0D0Kp_kin", "h_BuD0D0Kp_kin", 100, 0, 1);
    TH1D* h_BuDpDmKp_kin = new TH1D("h_BuDpDmKp_kin", "h_BuDpDmKp_kin", 100, 0, 1);
    TH1D* h_BuDsDsKp_kin = new TH1D("h_BuDsDsKp_kin", "h_BuDsDsKp_kin", 100, 0, 1);

    t_BuDDKp->Draw("BDT2 >> h_BuD0D0Kp_kin", "(species==100) && (is_best_cand==1)");
    t_BuDDKp->Draw("BDT2 >> h_BuDpDmKp_kin", "(species==101) && (is_best_cand==1)");
    t_BuDDKp->Draw("BDT2 >> h_BuDsDsKp_kin", "(species==102) && (is_best_cand==1)");

    h_BuDsDsKp_kin->GetXaxis()->SetTitle("BDT2");
    h_BuDsDsKp_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BuDsDsKp_kin->SetTitle("B^{+} #rightarrow D D K^{+} MC");

    h_BuD0D0Kp_kin->SetLineColor(kCyan+1);
    h_BuD0D0Kp_kin->SetFillColorAlpha(kCyan+1, 0.25);

    h_BuDpDmKp_kin->SetLineColor(kGreen+1);
    h_BuDpDmKp_kin->SetFillColorAlpha(kGreen+1, 0.25);

    h_BuDsDsKp_kin->SetLineColor(kViolet+1);
    h_BuDsDsKp_kin->SetFillColorAlpha(kViolet+1, 0.25);

    h_BuDsDsKp_kin->DrawNormalized();
    h_BuDpDmKp_kin->DrawNormalized("same");
    h_BuD0D0Kp_kin->DrawNormalized("same");
    leg7->Draw("same");

    SixPads1->cd(2) ; a2->DrawClonePad();

    TH1D* h_BdDpD0Kp_kin = new TH1D("h_BdDpD0Kp_kin", "h_BdDpD0Kp_kin", 100, 0, 1);

    t_BdDDKp->Draw("BDT2 >> h_BdDpD0Kp_kin", "(species==110) && (is_best_cand==1)");

    h_BdDpD0Kp_kin->GetXaxis()->SetTitle("BDT2");
    h_BdDpD0Kp_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BdDpD0Kp_kin->SetTitle("B^{0} #rightarrow D D K^{+} MC");

    h_BdDpD0Kp_kin->SetLineColor(kViolet+2);
    h_BdDpD0Kp_kin->SetFillColorAlpha(kViolet+2, 0.25);

    h_BdDpD0Kp_kin->DrawNormalized();
    leg8->Draw("same");

    SixPads1->cd(3) ; a3->DrawClonePad();

    TH1D* h_BsDsD0Kp_kin = new TH1D("h_BsDsD0Kp_kin", "h_BsDsD0Kp_kin", 100, 0, 1);

    t_BsDDKp->Draw("BDT2 >> h_BsDsD0Kp_kin", "(species==120) && (is_best_cand==1)");

    h_BsDsD0Kp_kin->GetXaxis()->SetTitle("BDT2");
    h_BsDsD0Kp_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BsDsD0Kp_kin->SetTitle("B^{0}_{s} #rightarrow D D K^{+} MC");

    h_BsDsD0Kp_kin->SetLineColor(kPink+7);
    h_BsDsD0Kp_kin->SetFillColorAlpha(kPink+7, 0.25);

    h_BsDsD0Kp_kin->DrawNormalized();
    leg9->Draw("same");

    SixPads1->cd(4) ; a4->DrawClonePad();

    TH1D* h_BuD0DpK0_kin = new TH1D("h_BuD0DpK0_kin", "h_BuD0DpK0_kin", 100, 0, 1);

    t_BuDDK0->Draw("BDT2 >> h_BuD0DpK0_kin", "(species==130) && (is_best_cand==1)");

    h_BuD0DpK0_kin->GetXaxis()->SetTitle("BDT2");
    h_BuD0DpK0_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BuD0DpK0_kin->SetTitle("B^{+} #rightarrow D D K^{0} MC");

    h_BuD0DpK0_kin->SetLineColor(kRed);
    h_BuD0DpK0_kin->SetFillColorAlpha(kRed, 0.25);

    h_BuD0DpK0_kin->DrawNormalized();
    leg10->Draw("same");

    SixPads1->cd(5) ; a5->DrawClonePad();

    TH1D* h_BuD0Ds_kin = new TH1D("h_BuD0Ds_kin", "h_BuD0Ds_kin", 100, 0, 1);
    TH1D* h_BuD0Dp_kin = new TH1D("h_BuD0Dp_kin", "h_BuD0Dp_kin", 100, 0, 1);

    t_BuDD->Draw("BDT2 >> h_BuD0Ds_kin", "(species==150) && (is_best_cand==1)");
    t_BuDD->Draw("BDT2 >> h_BuD0Dp_kin", "(species==151) && (is_best_cand==1)");

    h_BuD0Ds_kin->GetXaxis()->SetTitle("BDT2");
    h_BuD0Ds_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    h_BuD0Ds_kin->SetTitle("B^{+} #rightarrow D D MC");

    h_BuD0Ds_kin->SetLineColor(kOrange-3);
    h_BuD0Ds_kin->SetFillColorAlpha(kOrange-3, 0.25);

    h_BuD0Dp_kin->SetLineColor(kOrange+3);
    h_BuD0Dp_kin->SetFillColorAlpha(kOrange+3, 0.25);

    h_BuD0Ds_kin->DrawNormalized();
    h_BuD0Dp_kin->DrawNormalized("same");
    leg11->Draw("same");

    SixPads1->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_BDT2.pdf");
}


