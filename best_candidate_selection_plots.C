Double_t eps_error(Double_t Num, Double_t Den);

void best_candidate_selection_plots()
{
    gStyle->SetOptStat(0);
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
    TFileCollection* fc_dds_mc_best_2016 = new TFileCollection("fc_dds_mc_best_2016", "fc_dds_mc_best_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_72/multiple_events.txt");
    TFileCollection* fc_dds_mc_best_2017 = new TFileCollection("fc_dds_mc_best_2017", "fc_dds_mc_best_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_72/multiple_events.txt");
    TFileCollection* fc_dds_mc_best_2018 = new TFileCollection("fc_dds_mc_best_2018", "fc_dds_mc_best_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_72/multiple_events.txt");

    TChain* t_dds_mc_best = new TChain("DecayTree");
    TChain* t_dds_mc_best_2017 = new TChain("DecayTree");
    TChain* t_dds_mc_best_2018 = new TChain("DecayTree");

    t_dds_mc_best->AddFileInfoList((TCollection*)fc_dds_mc_best_2016->GetList());
    t_dds_mc_best_2017->AddFileInfoList((TCollection*)fc_dds_mc_best_2017->GetList());
    t_dds_mc_best_2018->AddFileInfoList((TCollection*)fc_dds_mc_best_2018->GetList());

    t_dds_mc_best->GetEntries();
    t_dds_mc_best_2017->GetEntries();
    t_dds_mc_best_2018->GetEntries();

    t_dds_mc_best->Add(t_dds_mc_best_2017);
    t_dds_mc_best->Add(t_dds_mc_best_2018);

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
    TFileCollection* fc_dds_data_best_2016 = new TFileCollection("fc_dds_data_best_2016", "fc_dds_data_best_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_8/multiple_events.txt");
    TFileCollection* fc_dds_data_best_2017 = new TFileCollection("fc_dds_data_best_2017", "fc_dds_data_best_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_8/multiple_events.txt");
    TFileCollection* fc_dds_data_best_2018 = new TFileCollection("fc_dds_data_best_2018", "fc_dds_data_best_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_8/multiple_events.txt");

    TChain* t_dds_data_best = new TChain("DecayTree");
    TChain* t_dds_data_best_2017 = new TChain("DecayTree");
    TChain* t_dds_data_best_2018 = new TChain("DecayTree");

    t_dds_data_best->AddFileInfoList((TCollection*)fc_dds_data_best_2016->GetList());
    t_dds_data_best_2017->AddFileInfoList((TCollection*)fc_dds_data_best_2017->GetList());
    t_dds_data_best_2018->AddFileInfoList((TCollection*)fc_dds_data_best_2018->GetList());

    t_dds_data_best->GetEntries();
    t_dds_data_best_2017->GetEntries();
    t_dds_data_best_2018->GetEntries();

    t_dds_data_best->Add(t_dds_data_best_2017);
    t_dds_data_best->Add(t_dds_data_best_2018);

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

    ////////////////////// Normalisation mode selections
    TFile* f_dds_mc = new TFile("/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_7/selections.root");
    TObjString* dds_mc_truthMatch_exp = (TObjString*)f_dds_mc->Get("truthMatch");
    TObjString* dds_trigger_exp = (TObjString*)f_dds_mc->Get("trigger");
    TObjString* dds_mc_rectangular_cuts_exp = (TObjString*)f_dds_mc->Get("rectangular_cuts");
    TObjString* dds_pass_mass_fit_exp = (TObjString*)f_dds_mc->Get("pass_mass_fit");
    TObjString* dds_fit_range_exp = (TObjString*)f_dds_mc->Get("fit_range");

    TCut dds_mc_truthMatch(dds_mc_truthMatch_exp->GetString().Data());
    TCut dds_trigger(dds_trigger_exp->GetString().Data());
    TCut dds_mc_rectangular_cuts(dds_mc_rectangular_cuts_exp->GetString().Data());
    TCut dds_pass_mass_fit(dds_pass_mass_fit_exp->GetString().Data());
    TCut dds_fit_range(dds_fit_range_exp->GetString().Data());
    TCut dds_best_cand("(is_best_cand == 1)");

    TFile* f_dds_data = new TFile("/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_8/selections.root");
    TObjString* dds_data_rectangular_cuts_exp = (TObjString*)f_dds_data->Get("rectangular_cuts");
    TCut dds_data_rectangular_cuts(dds_data_rectangular_cuts_exp->GetString().Data());

    TCut dds_mc_selections = dds_trigger+dds_mc_rectangular_cuts+dds_pass_mass_fit+dds_fit_range;
    TCut dds_data_selections = dds_trigger+dds_data_rectangular_cuts+dds_pass_mass_fit+dds_fit_range;

    //////////////////////////////// 
    // Average number of candidates per event
    Int_t N = 5;
    TH1D* h_nCand_dds_mc = new TH1D("h_nCand_dds_mc", "h_nCand_dds_mc", N, 0, N);
    TH1D* h_nCand_dds_data = new TH1D("h_nCand_dds_data", "h_nCand_dds_data", N, 0, N);

    t_dds_mc->Draw("n_candidates >> h_nCand_dds_mc", dds_mc_selections);
    t_dds_data->Draw("n_candidates >> h_nCand_dds_data", dds_data_selections);

    TCanvas B;
    B.cd();
    B.SetLeftMargin(0.15);
    h_nCand_dds_mc->GetXaxis()->SetTitle("Number of candidates per event");
    h_nCand_dds_mc->GetYaxis()->SetTitle("Entries / (1)");
    h_nCand_dds_mc->SetTitle(Form("Average number of candidates per event (B^{+} #rightarrow #bar{D}^{0} D^{+}_{s} MC): %.4f",h_nCand_dds_mc->GetMean()));
    h_nCand_dds_mc->Draw();
    B.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/n_candidates_dds_mc.pdf");

    TCanvas B1;
    B1.cd();
    B1.SetLeftMargin(0.15);
    h_nCand_dds_data->GetXaxis()->SetTitle("Number of candidates per event");
    h_nCand_dds_data->GetYaxis()->SetTitle("Entries / (1)");
    h_nCand_dds_data->SetTitle(Form("Average number of candidates per event (B^{+} #rightarrow #bar{D}^{0} D^{+}_{s} data): %.4f",h_nCand_dds_data->GetMean()));
    h_nCand_dds_data->Draw();
    B1.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/n_candidates_dds_data.pdf");

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
    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(h_dds_mc_tm, "Truth-matched", "lf");
    leg->AddEntry(h_dds_mc, "Randomly picked candidates", "lf");
    leg->Draw("same");
    c5.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_mc_TM_vs_best.pdf");
    
    TH1D* h_dds_mc_not_tm = new TH1D("h_dds_mc_not_tm", "h_dds_mc_not_tm", 100, 5235, 5355);

    t_dds_mc->Draw("Bp_dtf_M[0] >> h_dds_mc_not_tm", dds_best_cand+!dds_mc_truthMatch);

    h_dds_mc_not_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_mc_not_tm->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_mc_not_tm->SetTitle("DDs MC");

    h_dds_mc_not_tm->SetLineColor(kCyan+1);
    h_dds_mc_not_tm->SetFillColorAlpha(kCyan+1, 0.25);

    TCanvas d5;
    d5.cd();
    h_dds_mc_not_tm->DrawNormalized();
    TLegend* leg1 = new TLegend(0.5,0.8,0.9,0.9);
    leg1->AddEntry(h_dds_mc_not_tm, "Randomly picked candidates (not pass TM)", "lf");
    leg1->Draw("same");
    d5.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_mc_best_notTM.pdf");

    Double_t metric_num_norm = t_dds_mc->GetEntries(dds_mc_truthMatch+dds_best_cand);
    Double_t metric_den_norm = t_dds_mc->GetEntries(dds_mc_truthMatch);
    Double_t metric_norm = metric_num_norm/metric_den_norm;
    cout << "Metric (norm) = " << metric_norm*100 << " +/- " << eps_error(metric_num_norm,metric_den_norm) << " \\%" << endl;

    TH1D* h_dds_mc_not_best = new TH1D("h_dds_mc_not_best", "h_dds_mc_not_best", 100, 5235, 5355);
    t_dds_mc->Draw("Bp_dtf_M[0] >> h_dds_mc_not_best", !dds_best_cand);

    h_dds_mc_not_best->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_mc_not_best->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_mc_not_best->SetTitle("DDs MC");

    h_dds_mc_not_best->SetLineColor(kRed);
    h_dds_mc_not_best->SetFillColorAlpha(kRed, 0.25);

    TLegend* leg13 = new TLegend(0.5,0.8,0.9,0.9);
    leg13->AddEntry(h_dds_mc_not_best, "Candidates not picked", "lf");

    TCanvas d11;
    d11.cd();
    h_dds_mc_not_best->DrawNormalized();
    leg13->Draw("same");
    d11.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_mc_best_not_best.pdf");

    // DDs data
    TH1D* h_dds_data = new TH1D("h_dds_data", "h_dds_data", 100, 5235, 5355);
    TH1D* h_dds_data_best = new TH1D("h_dds_data_best", "h_dds_data_best", 100, 5235, 5355);

    t_dds_data->Draw("Bp_dtf_M[0] >> h_dds_data", !dds_best_cand);
    t_dds_data->Draw("Bp_dtf_M[0] >> h_dds_data_best", dds_best_cand);

    h_dds_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_data->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_data->SetTitle("DDs data");

    h_dds_data_best->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_dds_data_best->GetYaxis()->SetTitle("Normalized entries / (1.2 MeV)");
    h_dds_data_best->SetTitle("DDs data");

    h_dds_data->SetLineColor(kRed);
    h_dds_data->SetFillColorAlpha(kRed, 0.25);

    h_dds_data_best->SetLineColor(kBlack);
    h_dds_data_best->SetFillColorAlpha(kBlack, 0.25);

    TCanvas c8;
    c8.cd();
    h_dds_data_best->DrawNormalized();
    TLegend* leg22 = new TLegend(0.7,0.7,0.9,0.9);
    leg22->AddEntry(h_dds_data_best, "Randomly picked candidates");
    leg22->Draw("same");
    c8.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_data_best.pdf");

    TCanvas c9;
    c9.cd();
    h_dds_data->DrawNormalized();
    TLegend* leg23 = new TLegend(0.7,0.7,0.9,0.9);
    leg23->AddEntry(h_dds_data, "Candidates not picked");
    leg23->Draw("same");
    c9.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/dds_data_not_best.pdf");

}

Double_t eps_error(Double_t Num, Double_t Den)
{
    return (Num/Den)*sqrt( 1./Num + 1./Den );
}

    // //////////////////////////////////////////////////////////////// Ktautau ///////////////////////////////////////////////////////////////////////////
    // // TM MC
    // TFileCollection *fc_tm_ktautau_mc_2016 = new TFileCollection("fc_tm_ktautau_mc_2016", "fc_tm_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1/pre_sel_tree.txt");
    // TFileCollection *fc_tm_ktautau_mc_2017 = new TFileCollection("fc_tm_ktautau_mc_2017", "fc_tm_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1/pre_sel_tree.txt");
    // TFileCollection *fc_tm_ktautau_mc_2018 = new TFileCollection("fc_tm_ktautau_mc_2018", "fc_tm_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1/pre_sel_tree.txt");

    // TChain* t_tm_ktautau_mc = new TChain("DecayTree");
    // TChain* t_tm_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t_tm_ktautau_mc_2018 = new TChain("DecayTree");

    // t_tm_ktautau_mc->AddFileInfoList((TCollection*)fc_tm_ktautau_mc_2016->GetList());
    // t_tm_ktautau_mc_2017->AddFileInfoList((TCollection*)fc_tm_ktautau_mc_2017->GetList());
    // t_tm_ktautau_mc_2018->AddFileInfoList((TCollection*)fc_tm_ktautau_mc_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_tm_ktautau_mc->GetEntries() << endl;
    // cout << t_tm_ktautau_mc_2017->GetEntries() << endl;
    // cout << t_tm_ktautau_mc_2018->GetEntries() << endl;

    // t_tm_ktautau_mc->Add(t_tm_ktautau_mc_2017);
    // t_tm_ktautau_mc->Add(t_tm_ktautau_mc_2018);

    // TFileCollection *fc1_tm_ktautau_mc_2016 = new TFileCollection("fc1_tm_ktautau_mc_2016", "fc1_tm_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1/fit_results.txt");
    // TFileCollection *fc1_tm_ktautau_mc_2017 = new TFileCollection("fc1_tm_ktautau_mc_2017", "fc1_tm_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1/fit_results.txt");
    // TFileCollection *fc1_tm_ktautau_mc_2018 = new TFileCollection("fc1_tm_ktautau_mc_2018", "fc1_tm_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1/fit_results.txt");

    // TChain* t1_tm_ktautau_mc = new TChain("DecayTree");
    // TChain* t1_tm_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t1_tm_ktautau_mc_2018 = new TChain("DecayTree");

    // t1_tm_ktautau_mc->AddFileInfoList((TCollection*)fc1_tm_ktautau_mc_2016->GetList());
    // t1_tm_ktautau_mc_2017->AddFileInfoList((TCollection*)fc1_tm_ktautau_mc_2017->GetList());
    // t1_tm_ktautau_mc_2018->AddFileInfoList((TCollection*)fc1_tm_ktautau_mc_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_tm_ktautau_mc->GetEntries() << endl;
    // cout << t1_tm_ktautau_mc_2017->GetEntries() << endl;
    // cout << t1_tm_ktautau_mc_2018->GetEntries() << endl;

    // t1_tm_ktautau_mc->Add(t1_tm_ktautau_mc_2017);
    // t1_tm_ktautau_mc->Add(t1_tm_ktautau_mc_2018);

    // TFileCollection *fc3_tm_ktautau_mc_2016 = new TFileCollection("fc3_tm_ktautau_mc_2016", "fc3_tm_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_1/bdt_output.txt");
    // TFileCollection *fc3_tm_ktautau_mc_2017 = new TFileCollection("fc3_tm_ktautau_mc_2017", "fc3_tm_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_1/bdt_output.txt");
    // TFileCollection *fc3_tm_ktautau_mc_2018 = new TFileCollection("fc3_tm_ktautau_mc_2018", "fc3_tm_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_1/bdt_output.txt");

    // TChain* t3_tm_ktautau_mc = new TChain("DecayTree");
    // TChain* t3_tm_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t3_tm_ktautau_mc_2018 = new TChain("DecayTree");

    // t3_tm_ktautau_mc->AddFileInfoList((TCollection*)fc3_tm_ktautau_mc_2016->GetList());
    // t3_tm_ktautau_mc_2017->AddFileInfoList((TCollection*)fc3_tm_ktautau_mc_2017->GetList());
    // t3_tm_ktautau_mc_2018->AddFileInfoList((TCollection*)fc3_tm_ktautau_mc_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_tm_ktautau_mc->GetEntries() << endl;
    // cout << t3_tm_ktautau_mc_2017->GetEntries() << endl;
    // cout << t3_tm_ktautau_mc_2018->GetEntries() << endl;

    // t3_tm_ktautau_mc->Add(t3_tm_ktautau_mc_2017);
    // t3_tm_ktautau_mc->Add(t3_tm_ktautau_mc_2018);

    // TFileCollection *fc4_tm_ktautau_mc_2016 = new TFileCollection("fc4_tm_ktautau_mc_2016", "fc4_tm_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1/invariant_mass_tree.txt");
    // TFileCollection *fc4_tm_ktautau_mc_2017 = new TFileCollection("fc4_tm_ktautau_mc_2017", "fc4_tm_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1/invariant_mass_tree.txt");
    // TFileCollection *fc4_tm_ktautau_mc_2018 = new TFileCollection("fc4_tm_ktautau_mc_2018", "fc4_tm_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1/invariant_mass_tree.txt");

    // TChain* t4_tm_ktautau_mc = new TChain("DecayTree");
    // TChain* t4_tm_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t4_tm_ktautau_mc_2018 = new TChain("DecayTree");

    // t4_tm_ktautau_mc->AddFileInfoList((TCollection*)fc4_tm_ktautau_mc_2016->GetList());
    // t4_tm_ktautau_mc_2017->AddFileInfoList((TCollection*)fc4_tm_ktautau_mc_2017->GetList());
    // t4_tm_ktautau_mc_2018->AddFileInfoList((TCollection*)fc4_tm_ktautau_mc_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_tm_ktautau_mc->GetEntries() << endl;
    // cout << t4_tm_ktautau_mc_2017->GetEntries() << endl;
    // cout << t4_tm_ktautau_mc_2018->GetEntries() << endl;

    // t4_tm_ktautau_mc->Add(t4_tm_ktautau_mc_2017);
    // t4_tm_ktautau_mc->Add(t4_tm_ktautau_mc_2018);

    // Int_t N_tm_mc = t_tm_ktautau_mc->GetEntries();
    // Int_t N1_tm_mc = t1_tm_ktautau_mc->GetEntries();
    // Int_t N3_tm_mc = t3_tm_ktautau_mc->GetEntries();
    // Int_t N4_tm_mc = t4_tm_ktautau_mc->GetEntries();
    // if((N_tm_mc == N1_tm_mc) && (N_tm_mc == N3_tm_mc) && (N_tm_mc == N4_tm_mc))
    // {
    //     t_tm_ktautau_mc->AddFriend(t1_tm_ktautau_mc, "gsl");
    //     t_tm_ktautau_mc->AddFriend(t3_tm_ktautau_mc, "sklearn");
    //     t_tm_ktautau_mc->AddFriend(t4_tm_ktautau_mc, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // Un-TM MC
    // TFileCollection *fc_ktautau_mc_2016 = new TFileCollection("fc_ktautau_mc_2016", "fc_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_10/pre_sel_tree.txt");
    // TFileCollection *fc_ktautau_mc_2017 = new TFileCollection("fc_ktautau_mc_2017", "fc_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_10/pre_sel_tree.txt");
    // TFileCollection *fc_ktautau_mc_2018 = new TFileCollection("fc_ktautau_mc_2018", "fc_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_10/pre_sel_tree.txt");

    // TChain* t_ktautau_mc = new TChain("DecayTree");
    // TChain* t_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t_ktautau_mc_2018 = new TChain("DecayTree");

    // t_ktautau_mc->AddFileInfoList((TCollection*)fc_ktautau_mc_2016->GetList());
    // t_ktautau_mc_2017->AddFileInfoList((TCollection*)fc_ktautau_mc_2017->GetList());
    // t_ktautau_mc_2018->AddFileInfoList((TCollection*)fc_ktautau_mc_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_ktautau_mc->GetEntries() << endl;
    // cout << t_ktautau_mc_2017->GetEntries() << endl;
    // cout << t_ktautau_mc_2018->GetEntries() << endl;

    // t_ktautau_mc->Add(t_ktautau_mc_2017);
    // t_ktautau_mc->Add(t_ktautau_mc_2018);

    // TFileCollection *fc1_ktautau_mc_2016 = new TFileCollection("fc1_ktautau_mc_2016", "fc1_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_10/fit_results.txt");
    // TFileCollection *fc1_ktautau_mc_2017 = new TFileCollection("fc1_ktautau_mc_2017", "fc1_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_10/fit_results.txt");
    // TFileCollection *fc1_ktautau_mc_2018 = new TFileCollection("fc1_ktautau_mc_2018", "fc1_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_10/fit_results.txt");

    // TChain* t1_ktautau_mc = new TChain("DecayTree");
    // TChain* t1_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t1_ktautau_mc_2018 = new TChain("DecayTree");

    // t1_ktautau_mc->AddFileInfoList((TCollection*)fc1_ktautau_mc_2016->GetList());
    // t1_ktautau_mc_2017->AddFileInfoList((TCollection*)fc1_ktautau_mc_2017->GetList());
    // t1_ktautau_mc_2018->AddFileInfoList((TCollection*)fc1_ktautau_mc_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_ktautau_mc->GetEntries() << endl;
    // cout << t1_ktautau_mc_2017->GetEntries() << endl;
    // cout << t1_ktautau_mc_2018->GetEntries() << endl;

    // t1_ktautau_mc->Add(t1_ktautau_mc_2017);
    // t1_ktautau_mc->Add(t1_ktautau_mc_2018);

    // TFileCollection *fc2_ktautau_mc_2016 = new TFileCollection("fc2_ktautau_mc_2016", "fc2_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_10/multiple_events.txt");
    // TFileCollection *fc2_ktautau_mc_2017 = new TFileCollection("fc2_ktautau_mc_2017", "fc2_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_10/multiple_events.txt");
    // TFileCollection *fc2_ktautau_mc_2018 = new TFileCollection("fc2_ktautau_mc_2018", "fc2_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_10/multiple_events.txt");

    // TChain* t2_ktautau_mc = new TChain("DecayTree");
    // TChain* t2_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t2_ktautau_mc_2018 = new TChain("DecayTree");

    // t2_ktautau_mc->AddFileInfoList((TCollection*)fc2_ktautau_mc_2016->GetList());
    // t2_ktautau_mc_2017->AddFileInfoList((TCollection*)fc2_ktautau_mc_2017->GetList());
    // t2_ktautau_mc_2018->AddFileInfoList((TCollection*)fc2_ktautau_mc_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t2_ktautau_mc->GetEntries() << endl;
    // cout << t2_ktautau_mc_2017->GetEntries() << endl;
    // cout << t2_ktautau_mc_2018->GetEntries() << endl;

    // t2_ktautau_mc->Add(t2_ktautau_mc_2017);
    // t2_ktautau_mc->Add(t2_ktautau_mc_2018);

    // TFileCollection *fc3_ktautau_mc_2016 = new TFileCollection("fc3_ktautau_mc_2016", "fc3_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_10/bdt_output.txt");
    // TFileCollection *fc3_ktautau_mc_2017 = new TFileCollection("fc3_ktautau_mc_2017", "fc3_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_10/bdt_output.txt");
    // TFileCollection *fc3_ktautau_mc_2018 = new TFileCollection("fc3_ktautau_mc_2018", "fc3_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_10/bdt_output.txt");

    // TChain* t3_ktautau_mc = new TChain("DecayTree");
    // TChain* t3_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t3_ktautau_mc_2018 = new TChain("DecayTree");

    // t3_ktautau_mc->AddFileInfoList((TCollection*)fc3_ktautau_mc_2016->GetList());
    // t3_ktautau_mc_2017->AddFileInfoList((TCollection*)fc3_ktautau_mc_2017->GetList());
    // t3_ktautau_mc_2018->AddFileInfoList((TCollection*)fc3_ktautau_mc_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_ktautau_mc->GetEntries() << endl;
    // cout << t3_ktautau_mc_2017->GetEntries() << endl;
    // cout << t3_ktautau_mc_2018->GetEntries() << endl;

    // t3_ktautau_mc->Add(t3_ktautau_mc_2017);
    // t3_ktautau_mc->Add(t3_ktautau_mc_2018);

    // TFileCollection *fc4_ktautau_mc_2016 = new TFileCollection("fc4_ktautau_mc_2016", "fc4_ktautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_10/invariant_mass_tree.txt");
    // TFileCollection *fc4_ktautau_mc_2017 = new TFileCollection("fc4_ktautau_mc_2017", "fc4_ktautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_10/invariant_mass_tree.txt");
    // TFileCollection *fc4_ktautau_mc_2018 = new TFileCollection("fc4_ktautau_mc_2018", "fc4_ktautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_10/invariant_mass_tree.txt");

    // TChain* t4_ktautau_mc = new TChain("DecayTree");
    // TChain* t4_ktautau_mc_2017 = new TChain("DecayTree");
    // TChain* t4_ktautau_mc_2018 = new TChain("DecayTree");

    // t4_ktautau_mc->AddFileInfoList((TCollection*)fc4_ktautau_mc_2016->GetList());
    // t4_ktautau_mc_2017->AddFileInfoList((TCollection*)fc4_ktautau_mc_2017->GetList());
    // t4_ktautau_mc_2018->AddFileInfoList((TCollection*)fc4_ktautau_mc_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_ktautau_mc->GetEntries() << endl;
    // cout << t4_ktautau_mc_2017->GetEntries() << endl;
    // cout << t4_ktautau_mc_2018->GetEntries() << endl;

    // t4_ktautau_mc->Add(t4_ktautau_mc_2017);
    // t4_ktautau_mc->Add(t4_ktautau_mc_2018);

    // Int_t N_mc = t_ktautau_mc->GetEntries();
    // Int_t N1_mc = t1_ktautau_mc->GetEntries();
    // Int_t N2_mc = t2_ktautau_mc->GetEntries();
    // Int_t N3_mc = t3_ktautau_mc->GetEntries();
    // Int_t N4_mc = t4_ktautau_mc->GetEntries();
    // if((N_mc == N1_mc) && (N_mc == N2_mc) && (N_mc == N3_mc) && (N_mc == N4_mc))
    // {
    //     t_ktautau_mc->AddFriend(t1_ktautau_mc, "gsl");
    //     t_ktautau_mc->AddFriend(t2_ktautau_mc, "best_cand");
    //     t_ktautau_mc->AddFriend(t3_ktautau_mc, "sklearn");
    //     t_ktautau_mc->AddFriend(t4_ktautau_mc, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // RS data
    // TFileCollection *fc0_rs_2016 = new TFileCollection("fc0_rs_2016", "fc0_rs_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_2/pre_sel_tree.txt", 50);
    // TFileCollection *fc0_rs_2017 = new TFileCollection("fc0_rs_2017", "fc0_rs_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_2/pre_sel_tree.txt", 50);
    // TFileCollection *fc0_rs_2018 = new TFileCollection("fc0_rs_2018", "fc0_rs_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_2/pre_sel_tree.txt", 50);

    // TChain* t0_ktautau_rs = new TChain("DecayTree");
    // TChain* t0_ktautau_rs_2017 = new TChain("DecayTree");
    // TChain* t0_ktautau_rs_2018 = new TChain("DecayTree");

    // t0_ktautau_rs->AddFileInfoList((TCollection*)fc0_rs_2016->GetList());
    // t0_ktautau_rs_2017->AddFileInfoList((TCollection*)fc0_rs_2017->GetList());
    // t0_ktautau_rs_2018->AddFileInfoList((TCollection*)fc0_rs_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t0_ktautau_rs->GetEntries() << endl;
    // cout << t0_ktautau_rs_2017->GetEntries() << endl;
    // cout << t0_ktautau_rs_2018->GetEntries() << endl;

    // t0_ktautau_rs->Add(t0_ktautau_rs_2017);
    // t0_ktautau_rs->Add(t0_ktautau_rs_2018);

    // TFileCollection *fc_rs_2016 = new TFileCollection("fc_rs_2016", "fc_rs_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt", 50);
    // TFileCollection *fc_rs_2017 = new TFileCollection("fc_rs_2017", "fc_rs_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt", 50);
    // TFileCollection *fc_rs_2018 = new TFileCollection("fc_rs_2018", "fc_rs_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt", 50);

    // TChain* t_ktautau_rs_data = new TChain("DecayTree");
    // TChain* t_ktautau_rs_data_2017 = new TChain("DecayTree");
    // TChain* t_ktautau_rs_data_2018 = new TChain("DecayTree");

    // t_ktautau_rs_data->AddFileInfoList((TCollection*)fc_rs_2016->GetList());
    // t_ktautau_rs_data_2017->AddFileInfoList((TCollection*)fc_rs_2017->GetList());
    // t_ktautau_rs_data_2018->AddFileInfoList((TCollection*)fc_rs_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t_ktautau_rs_data->GetEntries() << endl;
    // cout << t_ktautau_rs_data_2017->GetEntries() << endl;
    // cout << t_ktautau_rs_data_2018->GetEntries() << endl;

    // t_ktautau_rs_data->Add(t_ktautau_rs_data_2017);
    // t_ktautau_rs_data->Add(t_ktautau_rs_data_2018);

    // TFileCollection *fc1_rs_2016 = new TFileCollection("fc1_rs_2016", "fc1_rs_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_2/multiple_events.txt", 50);
    // TFileCollection *fc1_rs_2017 = new TFileCollection("fc1_rs_2017", "fc1_rs_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_2/multiple_events.txt", 50);
    // TFileCollection *fc1_rs_2018 = new TFileCollection("fc1_rs_2018", "fc1_rs_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_2/multiple_events.txt", 50);

    // TChain* t1_ktautau_rs = new TChain("DecayTree");
    // TChain* t1_ktautau_rs_2017 = new TChain("DecayTree");
    // TChain* t1_ktautau_rs_2018 = new TChain("DecayTree");

    // t1_ktautau_rs->AddFileInfoList((TCollection*)fc1_rs_2016->GetList());
    // t1_ktautau_rs_2017->AddFileInfoList((TCollection*)fc1_rs_2017->GetList());
    // t1_ktautau_rs_2018->AddFileInfoList((TCollection*)fc1_rs_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t1_ktautau_rs->GetEntries() << endl;
    // cout << t1_ktautau_rs_2017->GetEntries() << endl;
    // cout << t1_ktautau_rs_2018->GetEntries() << endl;

    // t1_ktautau_rs->Add(t1_ktautau_rs_2017);
    // t1_ktautau_rs->Add(t1_ktautau_rs_2018);

    // TFileCollection *fc2_rs_2016 = new TFileCollection("fc2_rs_2016", "fc2_rs_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_2/invariant_mass_tree.txt", 50);
    // TFileCollection *fc2_rs_2017 = new TFileCollection("fc2_rs_2017", "fc2_rs_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_2/invariant_mass_tree.txt", 50);
    // TFileCollection *fc2_rs_2018 = new TFileCollection("fc2_rs_2018", "fc2_rs_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_2/invariant_mass_tree.txt", 50);

    // TChain* t2_ktautau_rs = new TChain("DecayTree");
    // TChain* t2_ktautau_rs_2017 = new TChain("DecayTree");
    // TChain* t2_ktautau_rs_2018 = new TChain("DecayTree");

    // t2_ktautau_rs->AddFileInfoList((TCollection*)fc2_rs_2016->GetList());
    // t2_ktautau_rs_2017->AddFileInfoList((TCollection*)fc2_rs_2017->GetList());
    // t2_ktautau_rs_2018->AddFileInfoList((TCollection*)fc2_rs_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t2_ktautau_rs->GetEntries() << endl;
    // cout << t2_ktautau_rs_2017->GetEntries() << endl;
    // cout << t2_ktautau_rs_2018->GetEntries() << endl;

    // t2_ktautau_rs->Add(t2_ktautau_rs_2017);
    // t2_ktautau_rs->Add(t2_ktautau_rs_2018);

    // Int_t N0_rs = t0_ktautau_rs->GetEntries();
    // Int_t N_rs = t_ktautau_rs_data->GetEntries();
    // Int_t N1_rs = t1_ktautau_rs->GetEntries();
    // Int_t N2_rs = t2_ktautau_rs->GetEntries();
    // if((N_rs == N0_rs) && (N_rs == N1_rs) && (N_rs == N2_rs))
    // {
    //     t_ktautau_rs_data->AddFriend(t1_ktautau_rs, "best_cand");
    //     t_ktautau_rs_data->AddFriend(t2_ktautau_rs, "mass");
    //     t_ktautau_rs_data->AddFriend(t0_ktautau_rs, "pre_sel");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // WS data
    // TFileCollection *fc0_ws_2016 = new TFileCollection("fc0_ws_2016", "fc0_ws_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 50);
    // TFileCollection *fc0_ws_2017 = new TFileCollection("fc0_ws_2017", "fc0_ws_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 50);
    // TFileCollection *fc0_ws_2018 = new TFileCollection("fc0_ws_2018", "fc0_ws_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 50);

    // TChain* t0_ktautau_ws = new TChain("DecayTree");
    // TChain* t0_ktautau_ws_2017 = new TChain("DecayTree");
    // TChain* t0_ktautau_ws_2018 = new TChain("DecayTree");

    // t0_ktautau_ws->AddFileInfoList((TCollection*)fc0_ws_2016->GetList());
    // t0_ktautau_ws_2017->AddFileInfoList((TCollection*)fc0_ws_2017->GetList());
    // t0_ktautau_ws_2018->AddFileInfoList((TCollection*)fc0_ws_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t0_ktautau_ws->GetEntries() << endl;
    // cout << t0_ktautau_ws_2017->GetEntries() << endl;
    // cout << t0_ktautau_ws_2018->GetEntries() << endl;

    // t0_ktautau_ws->Add(t0_ktautau_ws_2017);
    // t0_ktautau_ws->Add(t0_ktautau_ws_2018);

    // TFileCollection *fc_ws_2016 = new TFileCollection("fc_ws_2016", "fc_ws_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 50);
    // TFileCollection *fc_ws_2017 = new TFileCollection("fc_ws_2017", "fc_ws_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 50);
    // TFileCollection *fc_ws_2018 = new TFileCollection("fc_ws_2018", "fc_ws_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 50);

    // TChain* t_ktautau_ws_data = new TChain("DecayTree");
    // TChain* t_ktautau_ws_data_2017 = new TChain("DecayTree");
    // TChain* t_ktautau_ws_data_2018 = new TChain("DecayTree");

    // t_ktautau_ws_data->AddFileInfoList((TCollection*)fc_ws_2016->GetList());
    // t_ktautau_ws_data_2017->AddFileInfoList((TCollection*)fc_ws_2017->GetList());
    // t_ktautau_ws_data_2018->AddFileInfoList((TCollection*)fc_ws_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t_ktautau_ws_data->GetEntries() << endl;
    // cout << t_ktautau_ws_data_2017->GetEntries() << endl;
    // cout << t_ktautau_ws_data_2018->GetEntries() << endl;

    // t_ktautau_ws_data->Add(t_ktautau_ws_data_2017);
    // t_ktautau_ws_data->Add(t_ktautau_ws_data_2018);

    // TFileCollection *fc1_ws_2016 = new TFileCollection("fc1_ws_2016", "fc1_ws_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_3/multiple_events.txt", 50);
    // TFileCollection *fc1_ws_2017 = new TFileCollection("fc1_ws_2017", "fc1_ws_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_3/multiple_events.txt", 50);
    // TFileCollection *fc1_ws_2018 = new TFileCollection("fc1_ws_2018", "fc1_ws_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_3/multiple_events.txt", 50);

    // TChain* t1_ktautau_ws = new TChain("DecayTree");
    // TChain* t1_ktautau_ws_2017 = new TChain("DecayTree");
    // TChain* t1_ktautau_ws_2018 = new TChain("DecayTree");

    // t1_ktautau_ws->AddFileInfoList((TCollection*)fc1_ws_2016->GetList());
    // t1_ktautau_ws_2017->AddFileInfoList((TCollection*)fc1_ws_2017->GetList());
    // t1_ktautau_ws_2018->AddFileInfoList((TCollection*)fc1_ws_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t1_ktautau_ws->GetEntries() << endl;
    // cout << t1_ktautau_ws_2017->GetEntries() << endl;
    // cout << t1_ktautau_ws_2018->GetEntries() << endl;

    // t1_ktautau_ws->Add(t1_ktautau_ws_2017);
    // t1_ktautau_ws->Add(t1_ktautau_ws_2018);

    // TFileCollection *fc2_ws_2016 = new TFileCollection("fc2_ws_2016", "fc2_ws_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_3/invariant_mass_tree.txt", 50);
    // TFileCollection *fc2_ws_2017 = new TFileCollection("fc2_ws_2017", "fc2_ws_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_3/invariant_mass_tree.txt", 50);
    // TFileCollection *fc2_ws_2018 = new TFileCollection("fc2_ws_2018", "fc2_ws_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_3/invariant_mass_tree.txt", 50);

    // TChain* t2_ktautau_ws = new TChain("DecayTree");
    // TChain* t2_ktautau_ws_2017 = new TChain("DecayTree");
    // TChain* t2_ktautau_ws_2018 = new TChain("DecayTree");

    // t2_ktautau_ws->AddFileInfoList((TCollection*)fc2_ws_2016->GetList());
    // t2_ktautau_ws_2017->AddFileInfoList((TCollection*)fc2_ws_2017->GetList());
    // t2_ktautau_ws_2018->AddFileInfoList((TCollection*)fc2_ws_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t2_ktautau_ws->GetEntries() << endl;
    // cout << t2_ktautau_ws_2017->GetEntries() << endl;
    // cout << t2_ktautau_ws_2018->GetEntries() << endl;

    // t2_ktautau_ws->Add(t1_ktautau_ws_2017);
    // t2_ktautau_ws->Add(t1_ktautau_ws_2018);

    // Int_t N0_ws = t_ktautau_ws_data->GetEntries();
    // Int_t N_ws = t_ktautau_ws_data->GetEntries();
    // Int_t N1_ws = t1_ktautau_ws->GetEntries();
    // Int_t N2_ws = t2_ktautau_ws->GetEntries();
    // if((N_ws == N0_ws) && (N_ws == N1_ws) && (N_ws == N2_ws))
    // {
    //     t_ktautau_ws_data->AddFriend(t1_ktautau_ws, "best_cand");
    //     t_ktautau_ws_data->AddFriend(t2_ktautau_ws, "mass");
    //     t_ktautau_ws_data->AddFriend(t0_ktautau_ws, "pre_sel");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
        // ///////////////////////////////////////////////////////// Cocktail MCs ////////////////////////////////////////////////////////////////////////////////////
    // // B+ -> DD K+
    // TFileCollection *fc_BuDDKp_2016 = new TFileCollection("fc_BuDDKp_2016", "fc_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_100/pre_sel_tree.txt");
    // TFileCollection *fc_BuDDKp_2017 = new TFileCollection("fc_BuDDKp_2017", "fc_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_100/pre_sel_tree.txt");
    // TFileCollection *fc_BuDDKp_2018 = new TFileCollection("fc_BuDDKp_2018", "fc_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_100/pre_sel_tree.txt");

    // TChain* t_BuDDKp = new TChain("DecayTree");
    // TChain* t_BuDDKp_2017 = new TChain("DecayTree");
    // TChain* t_BuDDKp_2018 = new TChain("DecayTree");

    // t_BuDDKp->AddFileInfoList((TCollection*)fc_BuDDKp_2016->GetList());
    // t_BuDDKp_2017->AddFileInfoList((TCollection*)fc_BuDDKp_2017->GetList());
    // t_BuDDKp_2018->AddFileInfoList((TCollection*)fc_BuDDKp_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_BuDDKp->GetEntries() << endl;
    // cout << t_BuDDKp_2017->GetEntries() << endl;
    // cout << t_BuDDKp_2018->GetEntries() << endl;

    // t_BuDDKp->Add(t_BuDDKp_2017);
    // t_BuDDKp->Add(t_BuDDKp_2018);

    // TFileCollection *fc1_BuDDKp_2016 = new TFileCollection("fc1_BuDDKp_2016", "fc1_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_100/fit_results.txt");
    // TFileCollection *fc1_BuDDKp_2017 = new TFileCollection("fc1_BuDDKp_2017", "fc1_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_100/fit_results.txt");
    // TFileCollection *fc1_BuDDKp_2018 = new TFileCollection("fc1_BuDDKp_2018", "fc1_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_100/fit_results.txt");

    // TChain* t1_BuDDKp = new TChain("DecayTree");
    // TChain* t1_BuDDKp_2017 = new TChain("DecayTree");
    // TChain* t1_BuDDKp_2018 = new TChain("DecayTree");

    // t1_BuDDKp->AddFileInfoList((TCollection*)fc1_BuDDKp_2016->GetList());
    // t1_BuDDKp_2017->AddFileInfoList((TCollection*)fc1_BuDDKp_2017->GetList());
    // t1_BuDDKp_2018->AddFileInfoList((TCollection*)fc1_BuDDKp_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_BuDDKp->GetEntries() << endl;
    // cout << t1_BuDDKp_2017->GetEntries() << endl;
    // cout << t1_BuDDKp_2018->GetEntries() << endl;

    // t1_BuDDKp->Add(t1_BuDDKp_2017);
    // t1_BuDDKp->Add(t1_BuDDKp_2018);

    // TFileCollection *fc2_BuDDKp_2016 = new TFileCollection("fc2_BuDDKp_2016", "fc2_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_100/multiple_events.txt");
    // TFileCollection *fc2_BuDDKp_2017 = new TFileCollection("fc2_BuDDKp_2017", "fc2_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_100/multiple_events.txt");
    // TFileCollection *fc2_BuDDKp_2018 = new TFileCollection("fc2_BuDDKp_2018", "fc2_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_100/multiple_events.txt");

    // TChain* t2_BuDDKp = new TChain("DecayTree");
    // TChain* t2_BuDDKp_2017 = new TChain("DecayTree");
    // TChain* t2_BuDDKp_2018 = new TChain("DecayTree");

    // t2_BuDDKp->AddFileInfoList((TCollection*)fc2_BuDDKp_2016->GetList());
    // t2_BuDDKp_2017->AddFileInfoList((TCollection*)fc2_BuDDKp_2017->GetList());
    // t2_BuDDKp_2018->AddFileInfoList((TCollection*)fc2_BuDDKp_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t2_BuDDKp->GetEntries() << endl;
    // cout << t2_BuDDKp_2017->GetEntries() << endl;
    // cout << t2_BuDDKp_2018->GetEntries() << endl;

    // t2_BuDDKp->Add(t2_BuDDKp_2017);
    // t2_BuDDKp->Add(t2_BuDDKp_2018);

    // TFileCollection *fc3_BuDDKp_2016 = new TFileCollection("fc3_BuDDKp_2016", "fc3_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_100/bdt_output.txt");
    // TFileCollection *fc3_BuDDKp_2017 = new TFileCollection("fc3_BuDDKp_2017", "fc3_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_100/bdt_output.txt");
    // TFileCollection *fc3_BuDDKp_2018 = new TFileCollection("fc3_BuDDKp_2018", "fc3_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_100/bdt_output.txt");

    // TChain* t3_BuDDKp = new TChain("DecayTree");
    // TChain* t3_BuDDKp_2017 = new TChain("DecayTree");
    // TChain* t3_BuDDKp_2018 = new TChain("DecayTree");

    // t3_BuDDKp->AddFileInfoList((TCollection*)fc3_BuDDKp_2016->GetList());
    // t3_BuDDKp_2017->AddFileInfoList((TCollection*)fc3_BuDDKp_2017->GetList());
    // t3_BuDDKp_2018->AddFileInfoList((TCollection*)fc3_BuDDKp_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_BuDDKp->GetEntries() << endl;
    // cout << t3_BuDDKp_2017->GetEntries() << endl;
    // cout << t3_BuDDKp_2018->GetEntries() << endl;

    // t3_BuDDKp->Add(t3_BuDDKp_2017);
    // t3_BuDDKp->Add(t3_BuDDKp_2018);

    // TFileCollection *fc4_BuDDKp_2016 = new TFileCollection("fc4_BuDDKp_2016", "fc4_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_100/invariant_mass_tree.txt");
    // TFileCollection *fc4_BuDDKp_2017 = new TFileCollection("fc4_BuDDKp_2017", "fc4_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_100/invariant_mass_tree.txt");
    // TFileCollection *fc4_BuDDKp_2018 = new TFileCollection("fc4_BuDDKp_2018", "fc4_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_100/invariant_mass_tree.txt");

    // TChain* t4_BuDDKp = new TChain("DecayTree");
    // TChain* t4_BuDDKp_2017 = new TChain("DecayTree");
    // TChain* t4_BuDDKp_2018 = new TChain("DecayTree");

    // t4_BuDDKp->AddFileInfoList((TCollection*)fc4_BuDDKp_2016->GetList());
    // t4_BuDDKp_2017->AddFileInfoList((TCollection*)fc4_BuDDKp_2017->GetList());
    // t4_BuDDKp_2018->AddFileInfoList((TCollection*)fc4_BuDDKp_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_BuDDKp->GetEntries() << endl;
    // cout << t4_BuDDKp_2017->GetEntries() << endl;
    // cout << t4_BuDDKp_2018->GetEntries() << endl;

    // t4_BuDDKp->Add(t4_BuDDKp_2017);
    // t4_BuDDKp->Add(t4_BuDDKp_2018);

    // Int_t N_BuDDKp = t_BuDDKp->GetEntries();
    // Int_t N1_BuDDKp = t1_BuDDKp->GetEntries();
    // Int_t N2_BuDDKp = t2_BuDDKp->GetEntries();
    // Int_t N3_BuDDKp = t3_BuDDKp->GetEntries();
    // Int_t N4_BuDDKp = t4_BuDDKp->GetEntries();
    // if((N_BuDDKp == N1_BuDDKp) && (N_BuDDKp == N2_BuDDKp) && (N_BuDDKp == N3_BuDDKp) && (N_BuDDKp == N4_BuDDKp))
    // {
    //     t_BuDDKp->AddFriend(t1_BuDDKp, "gsl");
    //     t_BuDDKp->AddFriend(t2_BuDDKp, "best_cand");
    //     t_BuDDKp->AddFriend(t3_BuDDKp, "sklearn");
    //     t_BuDDKp->AddFriend(t4_BuDDKp, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // B0 -> DD K+
    // TFileCollection *fc_BdDDKp_2016 = new TFileCollection("fc_BdDDKp_2016", "fc_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_110/pre_sel_tree.txt");
    // TFileCollection *fc_BdDDKp_2017 = new TFileCollection("fc_BdDDKp_2017", "fc_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_110/pre_sel_tree.txt");
    // TFileCollection *fc_BdDDKp_2018 = new TFileCollection("fc_BdDDKp_2018", "fc_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_110/pre_sel_tree.txt");

    // TChain* t_BdDDKp = new TChain("DecayTree");
    // TChain* t_BdDDKp_2017 = new TChain("DecayTree");
    // TChain* t_BdDDKp_2018 = new TChain("DecayTree");

    // t_BdDDKp->AddFileInfoList((TCollection*)fc_BdDDKp_2016->GetList());
    // t_BdDDKp_2017->AddFileInfoList((TCollection*)fc_BdDDKp_2017->GetList());
    // t_BdDDKp_2018->AddFileInfoList((TCollection*)fc_BdDDKp_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_BdDDKp->GetEntries() << endl;
    // cout << t_BdDDKp_2017->GetEntries() << endl;
    // cout << t_BdDDKp_2018->GetEntries() << endl;

    // t_BdDDKp->Add(t_BdDDKp_2017);
    // t_BdDDKp->Add(t_BdDDKp_2018);

    // TFileCollection *fc1_BdDDKp_2016 = new TFileCollection("fc1_BdDDKp_2016", "fc1_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_110/fit_results.txt");
    // TFileCollection *fc1_BdDDKp_2017 = new TFileCollection("fc1_BdDDKp_2017", "fc1_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_110/fit_results.txt");
    // TFileCollection *fc1_BdDDKp_2018 = new TFileCollection("fc1_BdDDKp_2018", "fc1_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_110/fit_results.txt");

    // TChain* t1_BdDDKp = new TChain("DecayTree");
    // TChain* t1_BdDDKp_2017 = new TChain("DecayTree");
    // TChain* t1_BdDDKp_2018 = new TChain("DecayTree");

    // t1_BdDDKp->AddFileInfoList((TCollection*)fc1_BdDDKp_2016->GetList());
    // t1_BdDDKp_2017->AddFileInfoList((TCollection*)fc1_BdDDKp_2017->GetList());
    // t1_BdDDKp_2018->AddFileInfoList((TCollection*)fc1_BdDDKp_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_BdDDKp->GetEntries() << endl;
    // cout << t1_BdDDKp_2017->GetEntries() << endl;
    // cout << t1_BdDDKp_2018->GetEntries() << endl;

    // t1_BdDDKp->Add(t1_BdDDKp_2017);
    // t1_BdDDKp->Add(t1_BdDDKp_2018);

    // TFileCollection *fc2_BdDDKp_2016 = new TFileCollection("fc2_BdDDKp_2016", "fc2_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_110/multiple_events.txt");
    // TFileCollection *fc2_BdDDKp_2017 = new TFileCollection("fc2_BdDDKp_2017", "fc2_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_110/multiple_events.txt");
    // TFileCollection *fc2_BdDDKp_2018 = new TFileCollection("fc2_BdDDKp_2018", "fc2_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_110/multiple_events.txt");

    // TChain* t2_BdDDKp = new TChain("DecayTree");
    // TChain* t2_BdDDKp_2017 = new TChain("DecayTree");
    // TChain* t2_BdDDKp_2018 = new TChain("DecayTree");

    // t2_BdDDKp->AddFileInfoList((TCollection*)fc2_BdDDKp_2016->GetList());
    // t2_BdDDKp_2017->AddFileInfoList((TCollection*)fc2_BdDDKp_2017->GetList());
    // t2_BdDDKp_2018->AddFileInfoList((TCollection*)fc2_BdDDKp_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t2_BdDDKp->GetEntries() << endl;
    // cout << t2_BdDDKp_2017->GetEntries() << endl;
    // cout << t2_BdDDKp_2018->GetEntries() << endl;

    // t2_BdDDKp->Add(t2_BdDDKp_2017);
    // t2_BdDDKp->Add(t2_BdDDKp_2018);

    // TFileCollection *fc3_BdDDKp_2016 = new TFileCollection("fc3_BdDDKp_2016", "fc3_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_110/bdt_output.txt");
    // TFileCollection *fc3_BdDDKp_2017 = new TFileCollection("fc3_BdDDKp_2017", "fc3_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_110/bdt_output.txt");
    // TFileCollection *fc3_BdDDKp_2018 = new TFileCollection("fc3_BdDDKp_2018", "fc3_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_110/bdt_output.txt");

    // TChain* t3_BdDDKp = new TChain("DecayTree");
    // TChain* t3_BdDDKp_2017 = new TChain("DecayTree");
    // TChain* t3_BdDDKp_2018 = new TChain("DecayTree");

    // t3_BdDDKp->AddFileInfoList((TCollection*)fc3_BdDDKp_2016->GetList());
    // t3_BdDDKp_2017->AddFileInfoList((TCollection*)fc3_BdDDKp_2017->GetList());
    // t3_BdDDKp_2018->AddFileInfoList((TCollection*)fc3_BdDDKp_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_BdDDKp->GetEntries() << endl;
    // cout << t3_BdDDKp_2017->GetEntries() << endl;
    // cout << t3_BdDDKp_2018->GetEntries() << endl;

    // t3_BdDDKp->Add(t3_BdDDKp_2017);
    // t3_BdDDKp->Add(t3_BdDDKp_2018);

    // TFileCollection *fc4_BdDDKp_2016 = new TFileCollection("fc4_BdDDKp_2016", "fc4_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_110/invariant_mass_tree.txt");
    // TFileCollection *fc4_BdDDKp_2017 = new TFileCollection("fc4_BdDDKp_2017", "fc4_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_110/invariant_mass_tree.txt");
    // TFileCollection *fc4_BdDDKp_2018 = new TFileCollection("fc4_BdDDKp_2018", "fc4_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_110/invariant_mass_tree.txt");

    // TChain* t4_BdDDKp = new TChain("DecayTree");
    // TChain* t4_BdDDKp_2017 = new TChain("DecayTree");
    // TChain* t4_BdDDKp_2018 = new TChain("DecayTree");

    // t4_BdDDKp->AddFileInfoList((TCollection*)fc4_BdDDKp_2016->GetList());
    // t4_BdDDKp_2017->AddFileInfoList((TCollection*)fc4_BdDDKp_2017->GetList());
    // t4_BdDDKp_2018->AddFileInfoList((TCollection*)fc4_BdDDKp_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_BdDDKp->GetEntries() << endl;
    // cout << t4_BdDDKp_2017->GetEntries() << endl;
    // cout << t4_BdDDKp_2018->GetEntries() << endl;

    // t4_BdDDKp->Add(t4_BdDDKp_2017);
    // t4_BdDDKp->Add(t4_BdDDKp_2018);

    // Int_t N_BdDDKp = t_BdDDKp->GetEntries();
    // Int_t N1_BdDDKp = t1_BdDDKp->GetEntries();
    // Int_t N2_BdDDKp = t2_BdDDKp->GetEntries();
    // Int_t N3_BdDDKp = t3_BdDDKp->GetEntries();
    // Int_t N4_BdDDKp = t4_BdDDKp->GetEntries();
    // if((N_BdDDKp == N1_BdDDKp) && (N_BdDDKp == N2_BdDDKp) && (N_BdDDKp == N3_BdDDKp) && (N_BdDDKp == N4_BdDDKp))
    // {
    //     t_BdDDKp->AddFriend(t1_BdDDKp, "gsl");
    //     t_BdDDKp->AddFriend(t2_BdDDKp, "best_cand");
    //     t_BdDDKp->AddFriend(t3_BdDDKp, "sklearn");
    //     t_BdDDKp->AddFriend(t4_BdDDKp, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // Bs -> DD K+
    // TFileCollection *fc_BsDDKp_2016 = new TFileCollection("fc_BsDDKp_2016", "fc_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_120/pre_sel_tree.txt");
    // TFileCollection *fc_BsDDKp_2017 = new TFileCollection("fc_BsDDKp_2017", "fc_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_120/pre_sel_tree.txt");
    // TFileCollection *fc_BsDDKp_2018 = new TFileCollection("fc_BsDDKp_2018", "fc_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_120/pre_sel_tree.txt");

    // TChain* t_BsDDKp = new TChain("DecayTree");
    // TChain* t_BsDDKp_2017 = new TChain("DecayTree");
    // TChain* t_BsDDKp_2018 = new TChain("DecayTree");

    // t_BsDDKp->AddFileInfoList((TCollection*)fc_BsDDKp_2016->GetList());
    // t_BsDDKp_2017->AddFileInfoList((TCollection*)fc_BsDDKp_2017->GetList());
    // t_BsDDKp_2018->AddFileInfoList((TCollection*)fc_BsDDKp_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_BsDDKp->GetEntries() << endl;
    // cout << t_BsDDKp_2017->GetEntries() << endl;
    // cout << t_BsDDKp_2018->GetEntries() << endl;

    // t_BsDDKp->Add(t_BsDDKp_2017);
    // t_BsDDKp->Add(t_BsDDKp_2018);

    // TFileCollection *fc1_BsDDKp_2016 = new TFileCollection("fc1_BsDDKp_2016", "fc1_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_120/fit_results.txt");
    // TFileCollection *fc1_BsDDKp_2017 = new TFileCollection("fc1_BsDDKp_2017", "fc1_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_120/fit_results.txt");
    // TFileCollection *fc1_BsDDKp_2018 = new TFileCollection("fc1_BsDDKp_2018", "fc1_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_120/fit_results.txt");

    // TChain* t1_BsDDKp = new TChain("DecayTree");
    // TChain* t1_BsDDKp_2017 = new TChain("DecayTree");
    // TChain* t1_BsDDKp_2018 = new TChain("DecayTree");

    // t1_BsDDKp->AddFileInfoList((TCollection*)fc1_BsDDKp_2016->GetList());
    // t1_BsDDKp_2017->AddFileInfoList((TCollection*)fc1_BsDDKp_2017->GetList());
    // t1_BsDDKp_2018->AddFileInfoList((TCollection*)fc1_BsDDKp_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_BsDDKp->GetEntries() << endl;
    // cout << t1_BsDDKp_2017->GetEntries() << endl;
    // cout << t1_BsDDKp_2018->GetEntries() << endl;

    // t1_BsDDKp->Add(t1_BsDDKp_2017);
    // t1_BsDDKp->Add(t1_BsDDKp_2018);

    // TFileCollection *fc2_BsDDKp_2016 = new TFileCollection("fc2_BsDDKp_2016", "fc2_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_120/multiple_events.txt");
    // TFileCollection *fc2_BsDDKp_2017 = new TFileCollection("fc2_BsDDKp_2017", "fc2_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_120/multiple_events.txt");
    // TFileCollection *fc2_BsDDKp_2018 = new TFileCollection("fc2_BsDDKp_2018", "fc2_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_120/multiple_events.txt");

    // TChain* t2_BsDDKp = new TChain("DecayTree");
    // TChain* t2_BsDDKp_2017 = new TChain("DecayTree");
    // TChain* t2_BsDDKp_2018 = new TChain("DecayTree");

    // t2_BsDDKp->AddFileInfoList((TCollection*)fc2_BsDDKp_2016->GetList());
    // t2_BsDDKp_2017->AddFileInfoList((TCollection*)fc2_BsDDKp_2017->GetList());
    // t2_BsDDKp_2018->AddFileInfoList((TCollection*)fc2_BsDDKp_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t2_BsDDKp->GetEntries() << endl;
    // cout << t2_BsDDKp_2017->GetEntries() << endl;
    // cout << t2_BsDDKp_2018->GetEntries() << endl;

    // t2_BsDDKp->Add(t2_BsDDKp_2017);
    // t2_BsDDKp->Add(t2_BsDDKp_2018);

    // TFileCollection *fc3_BsDDKp_2016 = new TFileCollection("fc3_BsDDKp_2016", "fc3_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_120/bdt_output.txt");
    // TFileCollection *fc3_BsDDKp_2017 = new TFileCollection("fc3_BsDDKp_2017", "fc3_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_120/bdt_output.txt");
    // TFileCollection *fc3_BsDDKp_2018 = new TFileCollection("fc3_BsDDKp_2018", "fc3_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_120/bdt_output.txt");

    // TChain* t3_BsDDKp = new TChain("DecayTree");
    // TChain* t3_BsDDKp_2017 = new TChain("DecayTree");
    // TChain* t3_BsDDKp_2018 = new TChain("DecayTree");

    // t3_BsDDKp->AddFileInfoList((TCollection*)fc3_BsDDKp_2016->GetList());
    // t3_BsDDKp_2017->AddFileInfoList((TCollection*)fc3_BsDDKp_2017->GetList());
    // t3_BsDDKp_2018->AddFileInfoList((TCollection*)fc3_BsDDKp_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_BsDDKp->GetEntries() << endl;
    // cout << t3_BsDDKp_2017->GetEntries() << endl;
    // cout << t3_BsDDKp_2018->GetEntries() << endl;

    // t3_BsDDKp->Add(t3_BsDDKp_2017);
    // t3_BsDDKp->Add(t3_BsDDKp_2018);

    // TFileCollection *fc4_BsDDKp_2016 = new TFileCollection("fc4_BsDDKp_2016", "fc4_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_120/invariant_mass_tree.txt");
    // TFileCollection *fc4_BsDDKp_2017 = new TFileCollection("fc4_BsDDKp_2017", "fc4_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_120/invariant_mass_tree.txt");
    // TFileCollection *fc4_BsDDKp_2018 = new TFileCollection("fc4_BsDDKp_2018", "fc4_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_120/invariant_mass_tree.txt");

    // TChain* t4_BsDDKp = new TChain("DecayTree");
    // TChain* t4_BsDDKp_2017 = new TChain("DecayTree");
    // TChain* t4_BsDDKp_2018 = new TChain("DecayTree");

    // t4_BsDDKp->AddFileInfoList((TCollection*)fc4_BsDDKp_2016->GetList());
    // t4_BsDDKp_2017->AddFileInfoList((TCollection*)fc4_BsDDKp_2017->GetList());
    // t4_BsDDKp_2018->AddFileInfoList((TCollection*)fc4_BsDDKp_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_BsDDKp->GetEntries() << endl;
    // cout << t4_BsDDKp_2017->GetEntries() << endl;
    // cout << t4_BsDDKp_2018->GetEntries() << endl;

    // t4_BsDDKp->Add(t4_BsDDKp_2017);
    // t4_BsDDKp->Add(t4_BsDDKp_2018);

    // Int_t N_BsDDKp = t_BsDDKp->GetEntries();
    // Int_t N1_BsDDKp = t1_BsDDKp->GetEntries();
    // Int_t N2_BsDDKp = t2_BsDDKp->GetEntries();
    // Int_t N3_BsDDKp = t3_BsDDKp->GetEntries();
    // Int_t N4_BsDDKp = t4_BsDDKp->GetEntries();
    // if((N_BsDDKp == N1_BsDDKp) && (N_BsDDKp == N2_BsDDKp) && (N_BsDDKp == N3_BsDDKp) && (N_BsDDKp == N4_BsDDKp))
    // {
    //     t_BsDDKp->AddFriend(t1_BsDDKp, "gsl");
    //     t_BsDDKp->AddFriend(t2_BsDDKp, "best_cand");
    //     t_BsDDKp->AddFriend(t3_BsDDKp, "sklearn");
    //     t_BsDDKp->AddFriend(t4_BsDDKp, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // B+ -> DD K0
    // TFileCollection *fc_BuDDK0_2016 = new TFileCollection("fc_BuDDK0_2016", "fc_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_130/pre_sel_tree.txt");
    // TFileCollection *fc_BuDDK0_2017 = new TFileCollection("fc_BuDDK0_2017", "fc_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_130/pre_sel_tree.txt");
    // TFileCollection *fc_BuDDK0_2018 = new TFileCollection("fc_BuDDK0_2018", "fc_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_130/pre_sel_tree.txt");

    // TChain* t_BuDDK0 = new TChain("DecayTree");
    // TChain* t_BuDDK0_2017 = new TChain("DecayTree");
    // TChain* t_BuDDK0_2018 = new TChain("DecayTree");

    // t_BuDDK0->AddFileInfoList((TCollection*)fc_BuDDK0_2016->GetList());
    // t_BuDDK0_2017->AddFileInfoList((TCollection*)fc_BuDDK0_2017->GetList());
    // t_BuDDK0_2018->AddFileInfoList((TCollection*)fc_BuDDK0_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_BuDDK0->GetEntries() << endl;
    // cout << t_BuDDK0_2017->GetEntries() << endl;
    // cout << t_BuDDK0_2018->GetEntries() << endl;

    // t_BuDDK0->Add(t_BuDDK0_2017);
    // t_BuDDK0->Add(t_BuDDK0_2018);

    // TFileCollection *fc1_BuDDK0_2016 = new TFileCollection("fc1_BuDDK0_2016", "fc1_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_130/fit_results.txt");
    // TFileCollection *fc1_BuDDK0_2017 = new TFileCollection("fc1_BuDDK0_2017", "fc1_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_130/fit_results.txt");
    // TFileCollection *fc1_BuDDK0_2018 = new TFileCollection("fc1_BuDDK0_2018", "fc1_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_130/fit_results.txt");

    // TChain* t1_BuDDK0 = new TChain("DecayTree");
    // TChain* t1_BuDDK0_2017 = new TChain("DecayTree");
    // TChain* t1_BuDDK0_2018 = new TChain("DecayTree");

    // t1_BuDDK0->AddFileInfoList((TCollection*)fc1_BuDDK0_2016->GetList());
    // t1_BuDDK0_2017->AddFileInfoList((TCollection*)fc1_BuDDK0_2017->GetList());
    // t1_BuDDK0_2018->AddFileInfoList((TCollection*)fc1_BuDDK0_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_BuDDK0->GetEntries() << endl;
    // cout << t1_BuDDK0_2017->GetEntries() << endl;
    // cout << t1_BuDDK0_2018->GetEntries() << endl;

    // t1_BuDDK0->Add(t1_BuDDK0_2017);
    // t1_BuDDK0->Add(t1_BuDDK0_2018);

    // TFileCollection *fc2_BuDDK0_2016 = new TFileCollection("fc2_BuDDK0_2016", "fc2_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_130/multiple_events.txt");
    // TFileCollection *fc2_BuDDK0_2017 = new TFileCollection("fc2_BuDDK0_2017", "fc2_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_130/multiple_events.txt");
    // TFileCollection *fc2_BuDDK0_2018 = new TFileCollection("fc2_BuDDK0_2018", "fc2_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_130/multiple_events.txt");

    // TChain* t2_BuDDK0 = new TChain("DecayTree");
    // TChain* t2_BuDDK0_2017 = new TChain("DecayTree");
    // TChain* t2_BuDDK0_2018 = new TChain("DecayTree");

    // t2_BuDDK0->AddFileInfoList((TCollection*)fc2_BuDDK0_2016->GetList());
    // t2_BuDDK0_2017->AddFileInfoList((TCollection*)fc2_BuDDK0_2017->GetList());
    // t2_BuDDK0_2018->AddFileInfoList((TCollection*)fc2_BuDDK0_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t2_BuDDK0->GetEntries() << endl;
    // cout << t2_BuDDK0_2017->GetEntries() << endl;
    // cout << t2_BuDDK0_2018->GetEntries() << endl;

    // t2_BuDDK0->Add(t2_BuDDK0_2017);
    // t2_BuDDK0->Add(t2_BuDDK0_2018);

    // TFileCollection *fc3_BuDDK0_2016 = new TFileCollection("fc3_BuDDK0_2016", "fc3_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_130/bdt_output.txt");
    // TFileCollection *fc3_BuDDK0_2017 = new TFileCollection("fc3_BuDDK0_2017", "fc3_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_130/bdt_output.txt");
    // TFileCollection *fc3_BuDDK0_2018 = new TFileCollection("fc3_BuDDK0_2018", "fc3_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_130/bdt_output.txt");

    // TChain* t3_BuDDK0 = new TChain("DecayTree");
    // TChain* t3_BuDDK0_2017 = new TChain("DecayTree");
    // TChain* t3_BuDDK0_2018 = new TChain("DecayTree");

    // t3_BuDDK0->AddFileInfoList((TCollection*)fc3_BuDDK0_2016->GetList());
    // t3_BuDDK0_2017->AddFileInfoList((TCollection*)fc3_BuDDK0_2017->GetList());
    // t3_BuDDK0_2018->AddFileInfoList((TCollection*)fc3_BuDDK0_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_BuDDK0->GetEntries() << endl;
    // cout << t3_BuDDK0_2017->GetEntries() << endl;
    // cout << t3_BuDDK0_2018->GetEntries() << endl;

    // t3_BuDDK0->Add(t3_BuDDK0_2017);
    // t3_BuDDK0->Add(t3_BuDDK0_2018);

    // TFileCollection *fc4_BuDDK0_2016 = new TFileCollection("fc4_BuDDK0_2016", "fc4_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_130/invariant_mass_tree.txt");
    // TFileCollection *fc4_BuDDK0_2017 = new TFileCollection("fc4_BuDDK0_2017", "fc4_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_130/invariant_mass_tree.txt");
    // TFileCollection *fc4_BuDDK0_2018 = new TFileCollection("fc4_BuDDK0_2018", "fc4_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_130/invariant_mass_tree.txt");

    // TChain* t4_BuDDK0 = new TChain("DecayTree");
    // TChain* t4_BuDDK0_2017 = new TChain("DecayTree");
    // TChain* t4_BuDDK0_2018 = new TChain("DecayTree");

    // t4_BuDDK0->AddFileInfoList((TCollection*)fc4_BuDDK0_2016->GetList());
    // t4_BuDDK0_2017->AddFileInfoList((TCollection*)fc4_BuDDK0_2017->GetList());
    // t4_BuDDK0_2018->AddFileInfoList((TCollection*)fc4_BuDDK0_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_BuDDK0->GetEntries() << endl;
    // cout << t4_BuDDK0_2017->GetEntries() << endl;
    // cout << t4_BuDDK0_2018->GetEntries() << endl;

    // t4_BuDDK0->Add(t4_BuDDK0_2017);
    // t4_BuDDK0->Add(t4_BuDDK0_2018);

    // Int_t N_BuDDK0 = t_BuDDK0->GetEntries();
    // Int_t N1_BuDDK0 = t1_BuDDK0->GetEntries();
    // Int_t N2_BuDDK0 = t2_BuDDK0->GetEntries();
    // Int_t N3_BuDDK0 = t3_BuDDK0->GetEntries();
    // Int_t N4_BuDDK0 = t4_BuDDK0->GetEntries();
    // if((N_BuDDK0 == N1_BuDDK0) && (N_BuDDK0 == N2_BuDDK0) && (N_BuDDK0 == N3_BuDDK0) && (N_BuDDK0 == N4_BuDDK0))
    // {
    //     t_BuDDK0->AddFriend(t1_BuDDK0, "gsl");
    //     t_BuDDK0->AddFriend(t2_BuDDK0, "best_cand");
    //     t_BuDDK0->AddFriend(t3_BuDDK0, "sklearn");
    //     t_BuDDK0->AddFriend(t4_BuDDK0, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }

    // // B+ -> DD 
    // TFileCollection *fc_BuDD_2016 = new TFileCollection("fc_BuDD_2016", "fc_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_150/pre_sel_tree.txt");
    // TFileCollection *fc_BuDD_2017 = new TFileCollection("fc_BuDD_2017", "fc_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_150/pre_sel_tree.txt");
    // TFileCollection *fc_BuDD_2018 = new TFileCollection("fc_BuDD_2018", "fc_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_150/pre_sel_tree.txt");

    // TChain* t_BuDD = new TChain("DecayTree");
    // TChain* t_BuDD_2017 = new TChain("DecayTree");
    // TChain* t_BuDD_2018 = new TChain("DecayTree");

    // t_BuDD->AddFileInfoList((TCollection*)fc_BuDD_2016->GetList());
    // t_BuDD_2017->AddFileInfoList((TCollection*)fc_BuDD_2017->GetList());
    // t_BuDD_2018->AddFileInfoList((TCollection*)fc_BuDD_2018->GetList());

    // cout << "Pre-selection files" << endl;
    // cout << t_BuDD->GetEntries() << endl;
    // cout << t_BuDD_2017->GetEntries() << endl;
    // cout << t_BuDD_2018->GetEntries() << endl;

    // t_BuDD->Add(t_BuDD_2017);
    // t_BuDD->Add(t_BuDD_2018);

    // TFileCollection *fc1_BuDD_2016 = new TFileCollection("fc1_BuDD_2016", "fc1_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_150/fit_results.txt");
    // TFileCollection *fc1_BuDD_2017 = new TFileCollection("fc1_BuDD_2017", "fc1_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_150/fit_results.txt");
    // TFileCollection *fc1_BuDD_2018 = new TFileCollection("fc1_BuDD_2018", "fc1_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_150/fit_results.txt");

    // TChain* t1_BuDD = new TChain("DecayTree");
    // TChain* t1_BuDD_2017 = new TChain("DecayTree");
    // TChain* t1_BuDD_2018 = new TChain("DecayTree");

    // t1_BuDD->AddFileInfoList((TCollection*)fc1_BuDD_2016->GetList());
    // t1_BuDD_2017->AddFileInfoList((TCollection*)fc1_BuDD_2017->GetList());
    // t1_BuDD_2018->AddFileInfoList((TCollection*)fc1_BuDD_2018->GetList());

    // cout << "Fit results" << endl;
    // cout << t1_BuDD->GetEntries() << endl;
    // cout << t1_BuDD_2017->GetEntries() << endl;
    // cout << t1_BuDD_2018->GetEntries() << endl;

    // t1_BuDD->Add(t1_BuDD_2017);
    // t1_BuDD->Add(t1_BuDD_2018);

    // TFileCollection *fc2_BuDD_2016 = new TFileCollection("fc2_BuDD_2016", "fc2_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_150/multiple_events.txt");
    // TFileCollection *fc2_BuDD_2017 = new TFileCollection("fc2_BuDD_2017", "fc2_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_150/multiple_events.txt");
    // TFileCollection *fc2_BuDD_2018 = new TFileCollection("fc2_BuDD_2018", "fc2_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_150/multiple_events.txt");

    // TChain* t2_BuDD = new TChain("DecayTree");
    // TChain* t2_BuDD_2017 = new TChain("DecayTree");
    // TChain* t2_BuDD_2018 = new TChain("DecayTree");

    // t2_BuDD->AddFileInfoList((TCollection*)fc2_BuDD_2016->GetList());
    // t2_BuDD_2017->AddFileInfoList((TCollection*)fc2_BuDD_2017->GetList());
    // t2_BuDD_2018->AddFileInfoList((TCollection*)fc2_BuDD_2018->GetList());

    // cout << "Best candidate" << endl;
    // cout << t2_BuDD->GetEntries() << endl;
    // cout << t2_BuDD_2017->GetEntries() << endl;
    // cout << t2_BuDD_2018->GetEntries() << endl;

    // t2_BuDD->Add(t2_BuDD_2017);
    // t2_BuDD->Add(t2_BuDD_2018);

    // TFileCollection *fc3_BuDD_2016 = new TFileCollection("fc3_BuDD_2016", "fc3_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_150/bdt_output.txt");
    // TFileCollection *fc3_BuDD_2017 = new TFileCollection("fc3_BuDD_2017", "fc3_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_150/bdt_output.txt");
    // TFileCollection *fc3_BuDD_2018 = new TFileCollection("fc3_BuDD_2018", "fc3_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_150/bdt_output.txt");

    // TChain* t3_BuDD = new TChain("DecayTree");
    // TChain* t3_BuDD_2017 = new TChain("DecayTree");
    // TChain* t3_BuDD_2018 = new TChain("DecayTree");

    // t3_BuDD->AddFileInfoList((TCollection*)fc3_BuDD_2016->GetList());
    // t3_BuDD_2017->AddFileInfoList((TCollection*)fc3_BuDD_2017->GetList());
    // t3_BuDD_2018->AddFileInfoList((TCollection*)fc3_BuDD_2018->GetList());

    // cout << "Sklearn" << endl;
    // cout << t3_BuDD->GetEntries() << endl;
    // cout << t3_BuDD_2017->GetEntries() << endl;
    // cout << t3_BuDD_2018->GetEntries() << endl;

    // t3_BuDD->Add(t3_BuDD_2017);
    // t3_BuDD->Add(t3_BuDD_2018);

    // TFileCollection *fc4_BuDD_2016 = new TFileCollection("fc4_BuDD_2016", "fc4_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_150/invariant_mass_tree.txt");
    // TFileCollection *fc4_BuDD_2017 = new TFileCollection("fc4_BuDD_2017", "fc4_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_150/invariant_mass_tree.txt");
    // TFileCollection *fc4_BuDD_2018 = new TFileCollection("fc4_BuDD_2018", "fc4_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_150/invariant_mass_tree.txt");

    // TChain* t4_BuDD = new TChain("DecayTree");
    // TChain* t4_BuDD_2017 = new TChain("DecayTree");
    // TChain* t4_BuDD_2018 = new TChain("DecayTree");

    // t4_BuDD->AddFileInfoList((TCollection*)fc4_BuDD_2016->GetList());
    // t4_BuDD_2017->AddFileInfoList((TCollection*)fc4_BuDD_2017->GetList());
    // t4_BuDD_2018->AddFileInfoList((TCollection*)fc4_BuDD_2018->GetList());

    // cout << "Masses" << endl;
    // cout << t4_BuDD->GetEntries() << endl;
    // cout << t4_BuDD_2017->GetEntries() << endl;
    // cout << t4_BuDD_2018->GetEntries() << endl;

    // t4_BuDD->Add(t4_BuDD_2017);
    // t4_BuDD->Add(t4_BuDD_2018);

    // Int_t N_BuDD = t_BuDD->GetEntries();
    // Int_t N1_BuDD = t1_BuDD->GetEntries();
    // Int_t N2_BuDD = t2_BuDD->GetEntries();
    // Int_t N3_BuDD = t3_BuDD->GetEntries();
    // Int_t N4_BuDD = t4_BuDD->GetEntries();
    // if((N_BuDD == N1_BuDD) && (N_BuDD == N2_BuDD) && (N_BuDD == N3_BuDD) && (N_BuDD == N4_BuDD))
    // {
    //     t_BuDD->AddFriend(t1_BuDD, "gsl");
    //     t_BuDD->AddFriend(t2_BuDD, "best_cand");
    //     t_BuDD->AddFriend(t3_BuDD, "sklearn");
    //     t_BuDD->AddFriend(t4_BuDD, "mass");
    // }
    // else
    // {
    //     cout << "wrong number of entries" << endl;
    //     return;
    // }
    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // // Ktautau MC: TM vs best candidate selection vs 2nd best candidate
    // TH1D* h_ktautau_mc_3pi3pi_tm = new TH1D("h_ktautau_mc_3pi3pi_tm", "h_ktautau_mc_3pi3pi_tm", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pipi0_tm = new TH1D("h_ktautau_mc_3pi3pipi0_tm", "h_ktautau_mc_3pi3pipi0_tm", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pi2pi0_tm = new TH1D("h_ktautau_mc_3pi3pi2pi0_tm", "h_ktautau_mc_3pi3pi2pi0_tm", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_all_tm = new TH1D("h_ktautau_mc_all_tm", "h_ktautau_mc_all_tm", 100, 4000, 8000);

    // TH1D* h_ktautau_mc_3pi3pi = new TH1D("h_ktautau_mc_3pi3pi", "h_ktautau_mc_3pi3pi", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pipi0 = new TH1D("h_ktautau_mc_3pi3pipi0", "h_ktautau_mc_3pi3pipi0", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pi2pi0 = new TH1D("h_ktautau_mc_3pi3pi2pi0", "h_ktautau_mc_3pi3pi2pi0", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_all = new TH1D("h_ktautau_mc_all", "h_ktautau_mc_all", 100, 4000, 8000);

    // TH1D* h_ktautau_mc_3pi3pi_2nd_best = new TH1D("h_ktautau_mc_3pi3pi_2nd_best", "h_ktautau_mc_3pi3pi_2nd_best", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pipi0_2nd_best = new TH1D("h_ktautau_mc_3pi3pipi0_2nd_best", "h_ktautau_mc_3pi3pipi0_2nd_best", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pi2pi0_2nd_best = new TH1D("h_ktautau_mc_3pi3pi2pi0_2nd_best", "h_ktautau_mc_3pi3pi2pi0_2nd_best", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_all_2nd_best = new TH1D("h_ktautau_mc_all_2nd_best", "h_ktautau_mc_all_2nd_best", 100, 4000, 8000);

    // t_tm_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi_tm", "component==0"+ktautau_mc_post_sel_cuts);
    // t_tm_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0_tm", "component==1"+ktautau_mc_post_sel_cuts);
    // t_tm_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0_tm", "component==2"+ktautau_mc_post_sel_cuts);
    // t_tm_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all_tm", ktautau_mc_post_sel_cuts);

    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi", "(is_best_cand == 1) && (component==0)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0", "(is_best_cand == 1) && (component==1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0", "(is_best_cand == 1) && (component==2)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all", "(is_best_cand == 1)"+ktautau_mc_post_sel_cuts);

    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi_2nd_best", "(is_second_best_cand == 1) && (component==0)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0_2nd_best", "(is_second_best_cand == 1) && (component==1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0_2nd_best", "(is_second_best_cand == 1) && (component==2)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all_2nd_best", "(is_second_best_cand == 1)"+ktautau_mc_post_sel_cuts);

    // auto c1 = new TCanvas();
    // auto c2 = new TCanvas(); 
    // auto c3 = new TCanvas(); 
    // auto c4 = new TCanvas(); 

    // auto FourPads = new TCanvas("FourPads","FourPads");
    // FourPads->Divide(2,2);
    // FourPads->cd(1) ; c1->DrawClonePad();

    // h_ktautau_mc_3pi3pi_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi_tm->SetTitle("3#pi3#pi MC: 2016-2018");

    // h_ktautau_mc_3pi3pi_tm->SetLineColor(kBlue);
    // h_ktautau_mc_3pi3pi_tm->SetFillColorAlpha(kBlue, 0.25);

    // h_ktautau_mc_3pi3pi->SetLineColor(kBlack);
    // h_ktautau_mc_3pi3pi->SetFillColorAlpha(kBlack, 0.25);

    // h_ktautau_mc_3pi3pi_2nd_best->SetLineColor(kGreen+1);
    // h_ktautau_mc_3pi3pi_2nd_best->SetFillColorAlpha(kGreen+1, 0.25);

    // h_ktautau_mc_3pi3pi_tm->DrawNormalized();
    // h_ktautau_mc_3pi3pi->DrawNormalized("same");
    // h_ktautau_mc_3pi3pi_2nd_best->DrawNormalized("same");

    // TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    // leg->AddEntry(h_ktautau_mc_3pi3pi_tm, "Truth-matched", "lf");
    // leg->AddEntry(h_ktautau_mc_3pi3pi, "Best candidate", "lf");
    // leg->AddEntry(h_ktautau_mc_3pi3pi_2nd_best, "2nd best candidate");
    // leg->Draw("same");

    // FourPads->cd(2) ; c2->DrawClonePad();

    // h_ktautau_mc_3pi3pipi0_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pipi0_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pipi0_tm->SetTitle("3#pi3#pi #pi^{0} MC: 2016-2018");

    // h_ktautau_mc_3pi3pipi0_tm->SetLineColor(kBlue);
    // h_ktautau_mc_3pi3pipi0_tm->SetFillColorAlpha(kBlue, 0.25);

    // h_ktautau_mc_3pi3pipi0->SetLineColor(kBlack);
    // h_ktautau_mc_3pi3pipi0->SetFillColorAlpha(kBlack, 0.25);

    // h_ktautau_mc_3pi3pipi0_2nd_best->SetLineColor(kGreen+1);
    // h_ktautau_mc_3pi3pipi0_2nd_best->SetFillColorAlpha(kGreen+1, 0.25);

    // h_ktautau_mc_3pi3pipi0_tm->DrawNormalized();
    // h_ktautau_mc_3pi3pipi0->DrawNormalized("same");
    // h_ktautau_mc_3pi3pipi0_2nd_best->DrawNormalized("same");
    // leg->Draw("same");
    
    // FourPads->cd(3) ; c3->DrawClonePad();

    // h_ktautau_mc_3pi3pi2pi0_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi2pi0_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi2pi0_tm->SetTitle("3#pi3#pi 2#pi^{0} MC: 2016-2018");

    // h_ktautau_mc_3pi3pi2pi0_tm->SetLineColor(kBlue);
    // h_ktautau_mc_3pi3pi2pi0_tm->SetFillColorAlpha(kBlue, 0.25);

    // h_ktautau_mc_3pi3pi2pi0->SetLineColor(kBlack);
    // h_ktautau_mc_3pi3pi2pi0->SetFillColorAlpha(kBlack, 0.25);

    // h_ktautau_mc_3pi3pi2pi0_2nd_best->SetLineColor(kGreen+1);
    // h_ktautau_mc_3pi3pi2pi0_2nd_best->SetFillColorAlpha(kGreen+1, 0.25);

    // h_ktautau_mc_3pi3pi2pi0_tm->DrawNormalized();
    // h_ktautau_mc_3pi3pi2pi0->DrawNormalized("same");
    // h_ktautau_mc_3pi3pi2pi0_2nd_best->DrawNormalized("same");
    // leg->Draw("same");

    // FourPads->cd(4) ; c4->DrawClonePad();

    // h_ktautau_mc_all_tm->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_all_tm->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_all_tm->SetTitle("All MC: 2016-2018");

    // h_ktautau_mc_all_tm->SetLineColor(kBlue);
    // h_ktautau_mc_all_tm->SetFillColorAlpha(kBlue, 0.25);

    // h_ktautau_mc_all->SetLineColor(kBlack);
    // h_ktautau_mc_all->SetFillColorAlpha(kBlack, 0.25);

    // h_ktautau_mc_all_2nd_best->SetLineColor(kGreen+1);
    // h_ktautau_mc_all_2nd_best->SetFillColorAlpha(kGreen+1, 0.25);

    // h_ktautau_mc_all_tm->DrawNormalized();
    // h_ktautau_mc_all->DrawNormalized("same");
    // h_ktautau_mc_all_2nd_best->DrawNormalized("same");
    // leg->Draw("same");

    // FourPads->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_TM_vs_best_vs_2nd_best.pdf");

    // // Ktautau MC: not best candidate
    // TH1D* h_ktautau_mc_3pi3pi_not_best = new TH1D("h_ktautau_mc_3pi3pi_not_best", "h_ktautau_mc_3pi3pi_not_best", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pipi0_not_best = new TH1D("h_ktautau_mc_3pi3pipi0_not_best", "h_ktautau_mc_3pi3pipi0_not_best", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pi2pi0_not_best = new TH1D("h_ktautau_mc_3pi3pi2pi0_not_best", "h_ktautau_mc_3pi3pi2pi0_not_best", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_all_not_best = new TH1D("h_ktautau_mc_all_not_best", "h_ktautau_mc_all_not_best", 100, 4000, 8000);

    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi_not_best", "(is_best_cand == 0) && (component==0)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0_not_best", "(is_best_cand == 0) && (component==1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0_not_best", "(is_best_cand == 0) && (component==2)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all_not_best", "(is_best_cand == 0)"+ktautau_mc_post_sel_cuts);

    // auto g1 = new TCanvas(); 
    // auto g2 = new TCanvas(); 
    // auto g3 = new TCanvas(); 
    // auto g4 = new TCanvas(); 

    // auto FourPads2 = new TCanvas("FourPads2","FourPads2");
    // FourPads2->Divide(2,2);
    // FourPads2->cd(1) ; g1->DrawClonePad();

    // h_ktautau_mc_3pi3pi_not_best->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi_not_best->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi_not_best->SetTitle("3#pi3#pi MC: 2016-2018");

    // h_ktautau_mc_3pi3pi_not_best->SetLineColor(kRed);
    // h_ktautau_mc_3pi3pi_not_best->SetFillColorAlpha(kRed, 0.25);

    // h_ktautau_mc_3pi3pi_not_best->DrawNormalized();

    // TLegend* leg12 = new TLegend(0.5,0.8,0.9,0.9);
    // leg12->AddEntry(h_ktautau_mc_3pi3pi_not_best, "Not best candidate", "lf");
    // leg12->Draw("same");

    // FourPads2->cd(2) ; g2->DrawClonePad();

    // h_ktautau_mc_3pi3pipi0_not_best->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pipi0_not_best->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pipi0_not_best->SetTitle("3#pi3#pi #pi^{0} MC: 2016-2018");

    // h_ktautau_mc_3pi3pipi0_not_best->SetLineColor(kRed);
    // h_ktautau_mc_3pi3pipi0_not_best->SetFillColorAlpha(kRed, 0.25);

    // h_ktautau_mc_3pi3pipi0_not_best->DrawNormalized();
    // leg12->Draw("same");

    // FourPads2->cd(3) ; g3->DrawClonePad();

    // h_ktautau_mc_3pi3pi2pi0_not_best->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi2pi0_not_best->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi2pi0_not_best->SetTitle("3#pi3#pi 2#pi^{0} MC: 2016-2018");

    // h_ktautau_mc_3pi3pi2pi0_not_best->SetLineColor(kRed);
    // h_ktautau_mc_3pi3pi2pi0_not_best->SetFillColorAlpha(kRed, 0.25);

    // h_ktautau_mc_3pi3pi2pi0_not_best->DrawNormalized();
    // leg12->Draw("same");

    // FourPads2->cd(4) ; g4->DrawClonePad();

    // h_ktautau_mc_all_not_best->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_all_not_best->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_all_not_best->SetTitle("All MC: 2016-2018");

    // h_ktautau_mc_all_not_best->SetLineColor(kRed);
    // h_ktautau_mc_all_not_best->SetFillColorAlpha(kRed, 0.25);

    // h_ktautau_mc_all_not_best->DrawNormalized();
    // leg12->Draw("same");

    // FourPads2->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_not_best.pdf");


    // // Ktautau MC: mis-reconstructed i.e. pass best candidate selection but no truth match cuts
    // TCut ktautau_truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";

    // TH1D* h_ktautau_mc_3pi3pi_not_TM = new TH1D("h_ktautau_mc_3pi3pi_not_TM", "h_ktautau_mc_3pi3pi_not_TM", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pipi0_not_TM = new TH1D("h_ktautau_mc_3pi3pipi0_not_TM", "h_ktautau_mc_3pi3pipi0_not_TM", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_3pi3pi2pi0_not_TM = new TH1D("h_ktautau_mc_3pi3pi2pi0_not_TM", "h_ktautau_mc_3pi3pi2pi0_not_TM", 100, 4000, 8000);
    // TH1D* h_ktautau_mc_all_not_TM = new TH1D("h_ktautau_mc_all_not_TM", "h_ktautau_mc_all_not_TM", 100, 4000, 8000);

    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi_not_TM", "(is_best_cand == 1) && (component==0)"+!ktautau_truthMatch+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pipi0_not_TM", "(is_best_cand == 1) && (component==1)"+!ktautau_truthMatch+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_3pi3pi2pi0_not_TM", "(is_best_cand == 1) && (component==2)"+!ktautau_truthMatch+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_ktautau_mc_all_not_TM", "(is_best_cand == 1)"+!ktautau_truthMatch+ktautau_mc_post_sel_cuts);

    // Double_t metric_num = t_ktautau_mc->GetEntries(ktautau_truthMatch+"is_best_cand == 1");
    // Double_t metric_den = t_ktautau_mc->GetEntries(ktautau_truthMatch);
    // Double_t metric = metric_num/metric_den;
    // cout << "Metric (sig) = " << metric*100  << " +/- " << eps_error(metric_num, metric_den)*100 << " \\%" << endl;

    // // Double_t metric_num_1 = t_BuDDKp->GetEntries(ktautau_truthMatch+"is_best_cand == 1");
    // // Double_t metric_den_1 = t_BuDDKp->GetEntries(ktautau_truthMatch);
    // // Double_t metric_1 = metric_num_1/metric_den_1;
    // // cout << "Metric (BuDDKp) = " << metric_1*100  << " +/- " << eps_error(metric_num_1, metric_den_1)*100 << " \\%" << endl;

    // // Double_t metric_num_2 = t_BdDDKp->GetEntries(ktautau_truthMatch+"is_best_cand == 1");
    // // Double_t metric_den_2 = t_BdDDKp->GetEntries(ktautau_truthMatch);
    // // Double_t metric_2 = metric_num_2/metric_den_2;
    // // cout << "Metric (BdDDKp) = " << metric_2*100  << " +/- " << eps_error(metric_num_2, metric_den_2)*100 << " \\%" << endl;

    // // Double_t metric_num_3 = t_BsDDKp->GetEntries(ktautau_truthMatch+"is_best_cand == 1");
    // // Double_t metric_den_3 = t_BsDDKp->GetEntries(ktautau_truthMatch);
    // // Double_t metric_3 = metric_num_3/metric_den_3;
    // // cout << "Metric (BsDDKp) = " << metric_3*100  << " +/- " << eps_error(metric_num_3, metric_den_3)*100 << " \\%" << endl;

    // // Double_t metric_num_4 = t_BuDDK0->GetEntries(ktautau_truthMatch+"is_best_cand == 1");
    // // Double_t metric_den_4 = t_BuDDK0->GetEntries(ktautau_truthMatch);
    // // Double_t metric_4 = metric_num_4/metric_den_4;
    // // cout << "Metric (BuDDK0) = " << metric_4*100  << " +/- " << eps_error(metric_num_4, metric_den_4)*100 << " \\%" << endl;

    // // Double_t metric_num_5 = t_BuDD->GetEntries(ktautau_truthMatch+"is_best_cand == 1");
    // // Double_t metric_den_5 = t_BuDD->GetEntries(ktautau_truthMatch);
    // // Double_t metric_5 = metric_num_5/metric_den_5;
    // // cout << "Metric (BuDD) = " << metric_5*100  << " +/- " << eps_error(metric_num_5, metric_den_5)*100 << " \\%" << endl;

    // auto d1 = new TCanvas(); 
    // auto d2 = new TCanvas(); 
    // auto d3 = new TCanvas(); 
    // auto d4 = new TCanvas(); 

    // auto FourPads1 = new TCanvas("FourPads1","FourPads1");
    // FourPads1->Divide(2,2);
    // FourPads1->cd(1) ; d1->DrawClonePad();

    // h_ktautau_mc_3pi3pi_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi_not_TM->SetTitle("3#pi3#pi MC: 2016-2018");

    // h_ktautau_mc_3pi3pi_not_TM->SetLineColor(kCyan+1);
    // h_ktautau_mc_3pi3pi_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    // h_ktautau_mc_3pi3pi_not_TM->DrawNormalized();

    // TLegend* leg1 = new TLegend(0.5,0.8,0.9,0.9);
    // leg1->AddEntry(h_ktautau_mc_3pi3pi_not_TM, "Best candidate (not pass TM)", "lf");
    // leg1->Draw("same");

    // FourPads1->cd(2) ; d2->DrawClonePad();

    // h_ktautau_mc_3pi3pipi0_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pipi0_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pipi0_not_TM->SetTitle("3#pi3#pi #pi^{0} MC: 2016-2018");

    // h_ktautau_mc_3pi3pipi0_not_TM->SetLineColor(kCyan+1);
    // h_ktautau_mc_3pi3pipi0_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    // h_ktautau_mc_3pi3pipi0_not_TM->DrawNormalized();
    // leg1->Draw("same");

    // FourPads1->cd(3) ; d3->DrawClonePad();

    // h_ktautau_mc_3pi3pi2pi0_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi2pi0_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi2pi0_not_TM->SetTitle("3#pi3#pi 2#pi^{0} MC: 2016-2018");

    // h_ktautau_mc_3pi3pi2pi0_not_TM->SetLineColor(kCyan+1);
    // h_ktautau_mc_3pi3pi2pi0_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    // h_ktautau_mc_3pi3pi2pi0_not_TM->DrawNormalized();
    // leg1->Draw("same");

    // FourPads1->cd(4) ; d4->DrawClonePad();

    // h_ktautau_mc_all_not_TM->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_all_not_TM->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_all_not_TM->SetTitle("All MC: 2016-2018");

    // h_ktautau_mc_all_not_TM->SetLineColor(kCyan+1);
    // h_ktautau_mc_all_not_TM->SetFillColorAlpha(kCyan+1, 0.25);

    // h_ktautau_mc_all_not_TM->DrawNormalized();
    // leg1->Draw("same");

    // FourPads1->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_best_notTM.pdf");

    // // Ktautau data
    // TH1D* h_rs_data = new TH1D("h_rs_data", "h_rs_data", 100, 4000, 8000);
    // TH1D* h_rs_data_best = new TH1D("h_rs_data_best", "h_rs_data_best", 100, 4000, 8000);

    // t_ktautau_rs_data->Draw("df_Bp_M >> h_rs_data", "(is_best_cand == 0)"+ktautau_data_post_sel_cuts);
    // t_ktautau_rs_data->Draw("df_Bp_M >> h_rs_data_best", "(is_best_cand == 1)"+ktautau_data_post_sel_cuts);

    // h_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_rs_data->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_rs_data->SetTitle("RS data");

    // h_rs_data->SetLineColor(kRed);
    // h_rs_data->SetFillColorAlpha(kRed, 0.25);

    // h_rs_data_best->SetLineColor(kBlack);
    // h_rs_data_best->SetFillColorAlpha(kBlack, 0.25);

    // TCanvas c6;
    // c6.cd();
    // h_rs_data->DrawNormalized();
    // h_rs_data_best->DrawNormalized("same");
    // TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    // leg2->AddEntry(h_rs_data_best, "Best candidate");
    // leg2->AddEntry(h_rs_data, "Not best candidate");
    // leg2->Draw("same");
    // c6.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/rs_data_multiple_vs_best.pdf");


    // TH1D* h_ws_data = new TH1D("h_ws_data", "h_ws_data", 100, 4000, 8000);
    // TH1D* h_ws_data_best = new TH1D("h_ws_data_best", "h_ws_data_best", 100, 4000, 8000);

    // t_ktautau_ws_data->Draw("df_Bp_M >> h_ws_data", "(is_best_cand == 0)"+ktautau_data_post_sel_cuts);
    // t_ktautau_ws_data->Draw("df_Bp_M >> h_ws_data_best", "(is_best_cand == 1)"+ktautau_data_post_sel_cuts);

    // h_ws_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ws_data->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ws_data->SetTitle("WS data");

    // h_ws_data->SetLineColor(kRed);
    // h_ws_data->SetFillColorAlpha(kRed, 0.25);

    // h_ws_data_best->SetLineColor(kBlack);
    // h_ws_data_best->SetFillColorAlpha(kBlack, 0.25);

    // TCanvas c7;
    // c7.cd();
    // h_ws_data->DrawNormalized();
    // h_ws_data_best->DrawNormalized("same");
    // leg2->Draw("same");
    // c7.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ws_data_multiple_vs_best.pdf");

    // // Cocktail MCs
    // TH1D* h_BuD0D0Kp_norm_noBDT = new TH1D("h_BuD0D0Kp_norm_noBDT", "h_BuD0D0Kp_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BuDpDmKp_norm_noBDT = new TH1D("h_BuDpDmKp_norm_noBDT", "h_BuDpDmKp_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BuDsDsKp_norm_noBDT = new TH1D("h_BuDsDsKp_norm_noBDT", "h_BuDsDsKp_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BdDpD0Kp_norm_noBDT = new TH1D("h_BdDpD0Kp_norm_noBDT", "h_BdDpD0Kp_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BsDsD0Kp_norm_noBDT = new TH1D("h_BsDsD0Kp_norm_noBDT", "h_BsDsD0Kp_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BuD0DpK0_norm_noBDT = new TH1D("h_BuD0DpK0_norm_noBDT", "h_BuD0DpK0_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BuD0Ds_norm_noBDT = new TH1D("h_BuD0Ds_norm_noBDT", "h_BuD0Ds_norm_noBDT", 100, 4000, 8000);
    // TH1D* h_BuD0Dp_norm_noBDT = new TH1D("h_BuD0Dp_norm_noBDT", "h_BuD0Dp _norm_noBDT", 100, 4000, 8000);

    // t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_norm_noBDT", "(species == 100) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_norm_noBDT", "(species == 101) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_norm_noBDT", "(species == 102) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_norm_noBDT", "(species == 110) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);  
    // t_BsDDKp->Draw("df_Bp_M >> h_BsDsD0Kp_norm_noBDT", "(species == 120) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_norm_noBDT", "(species == 130) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_norm_noBDT", "(species == 150) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_norm_noBDT", "(species == 151) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);

    // h_ktautau_mc_3pi3pi->SetLineColor(kBlue);
    // h_ktautau_mc_3pi3pi->SetFillColorAlpha(kBlue, 0.25);

    // h_ktautau_mc_all->SetLineColor(kBlack);
    // h_ktautau_mc_all->SetFillColorAlpha(kBlack, 0.25);

    // h_BuD0D0Kp_norm_noBDT->SetLineColor(kCyan+1);
    // h_BuD0D0Kp_norm_noBDT->SetFillColorAlpha(kCyan+1, 0.25);

    // h_BuDpDmKp_norm_noBDT->SetLineColor(kGreen+1);
    // h_BuDpDmKp_norm_noBDT->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuDsDsKp_norm_noBDT->SetLineColor(kViolet+1);
    // h_BuDsDsKp_norm_noBDT->SetFillColorAlpha(kViolet+1, 0.25);

    // h_BdDpD0Kp_norm_noBDT->SetLineColor(kViolet+2);
    // h_BdDpD0Kp_norm_noBDT->SetFillColorAlpha(kViolet+2, 0.25);

    // h_BsDsD0Kp_norm_noBDT->SetLineColor(kPink+7);    
    // h_BsDsD0Kp_norm_noBDT->SetFillColorAlpha(kPink+7, 0.25);

    // h_BuD0DpK0_norm_noBDT->SetLineColor(kRed);
    // h_BuD0DpK0_norm_noBDT->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0Ds_norm_noBDT->SetLineColor(kOrange-3);
    // h_BuD0Ds_norm_noBDT->SetFillColorAlpha(kOrange-3, 0.25);

    // h_BuD0Dp_norm_noBDT->SetLineColor(kOrange+4);
    // h_BuD0Dp_norm_noBDT->SetFillColorAlpha(kOrange+4, 0.25);

    // TCanvas d6;
    // d6.cd();
    // h_ktautau_mc_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_ktautau_mc_3pi3pi->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_ktautau_mc_3pi3pi->SetTitle("");
    // h_ktautau_mc_3pi3pi->DrawNormalized();
    // h_ktautau_mc_all->DrawNormalized("same");
    // h_BuD0D0Kp_norm_noBDT->DrawNormalized("same");
    // h_BuDpDmKp_norm_noBDT->DrawNormalized("same");
    // h_BuDsDsKp_norm_noBDT->DrawNormalized("same");
    // h_BdDpD0Kp_norm_noBDT->DrawNormalized("same");
    // h_BsDsD0Kp_norm_noBDT->DrawNormalized("same");
    // h_BuD0DpK0_norm_noBDT->DrawNormalized("same");
    // h_BuD0Ds_norm_noBDT->DrawNormalized("same");
    // h_BuD0Dp_norm_noBDT->DrawNormalized("same");
    // TLegend* leg3 = new TLegend(0.65,0.3,0.9,0.9);
    // leg3->AddEntry(h_ktautau_mc_3pi3pi, "3#pi3#pi MC");
    // leg3->AddEntry(h_ktautau_mc_all, "All MC");
    // leg3->AddEntry(h_BuD0D0Kp_norm_noBDT, "B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+}");
    // leg3->AddEntry(h_BuDpDmKp_norm_noBDT, "B^{+} #rightarrow D^{+} D^{-} K^{+}");
    // leg3->AddEntry(h_BuDsDsKp_norm_noBDT, "B^{+} #rightarrow D^{+}_{s} D^{+}_{s} K^{+}");
    // leg3->AddEntry(h_BdDpD0Kp_norm_noBDT, "B^{0} #rightarrow D^{-} D^{0} K^{+}");
    // leg3->AddEntry(h_BsDsD0Kp_norm_noBDT, "B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+}");
    // leg3->AddEntry(h_BuD0DpK0_norm_noBDT, "B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0}");
    // leg3->AddEntry(h_BuD0Ds_norm_noBDT, "B^{+} #rightarrow  #bar{D}^{0} D^{+}_{s}");
    // leg3->AddEntry(h_BuD0Dp_norm_noBDT, "B^{+} #rightarrow  #bar{D}^{0} D^{+}");
    // leg3->Draw("same");
    // d6.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_best_bmass.pdf");

    // TH1D* h_3pi3pi_norm_BDTcut = new TH1D("h_3pi3pi_norm_BDTcut", "h_3pi3pi_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_all_mc_norm_BDTcut = new TH1D("h_all_mc_norm_BDTcut", "h_all_mc_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BuD0D0Kp_norm_BDTcut = new TH1D("h_BuD0D0Kp_norm_BDTcut", "h_BuD0D0Kp_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BuDpDmKp_norm_BDTcut = new TH1D("h_BuDpDmKp_norm_BDTcut", "h_BuDpDmKp_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BuDsDsKp_norm_BDTcut = new TH1D("h_BuDsDsKp_norm_BDTcut", "h_BuDsDsKp_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BdDpD0Kp_norm_BDTcut = new TH1D("h_BdDpD0Kp_norm_BDTcut", "h_BdDpD0Kp_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BsDsD0Kp_norm_BDTcut = new TH1D("h_BsDsD0Kp_norm_BDTcut", "h_BsDsD0Kp_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BuD0DpK0_norm_BDTcut = new TH1D("h_BuD0DpK0_norm_BDTcut", "h_BuD0DpK0_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BuD0Ds_norm_BDTcut = new TH1D("h_BuD0Ds_norm_BDTcut", "h_BuD0Ds_norm_BDTcut", 100, 4000, 8000);
    // TH1D* h_BuD0Dp_norm_BDTcut = new TH1D("h_BuD0Dp_norm_BDTcut", "h_BuD0Dp _norm_BDTcut", 100, 4000, 8000);

    // t_ktautau_mc->Draw("df_Bp_M >> h_3pi3pi_norm_BDTcut", "(component==0) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("df_Bp_M >> h_all_mc_norm_BDTcut", "(is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_norm_BDTcut", "(species == 100) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_norm_BDTcut", "(species == 101) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_norm_BDTcut", "(species == 102) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_norm_BDTcut", "(species == 110) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);  
    // t_BsDDKp->Draw("df_Bp_M >> h_BsDsD0Kp_norm_BDTcut", "(species == 120) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_norm_BDTcut", "(species == 130) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_norm_BDTcut", "(species == 150) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_norm_BDTcut", "(species == 151) && (is_best_cand == 1) && (BDT1 > 0.9) && (BDT2 > 0.9)"+ktautau_mc_post_sel_cuts);

    // h_3pi3pi_norm_BDTcut->SetLineColor(kBlue);
    // h_3pi3pi_norm_BDTcut->SetFillColorAlpha(kBlue, 0.25);

    // h_all_mc_norm_BDTcut->SetLineColor(kBlack);
    // h_all_mc_norm_BDTcut->SetFillColorAlpha(kBlack, 0.25);

    // h_BuD0D0Kp_norm_BDTcut->SetLineColor(kCyan+1);
    // h_BuD0D0Kp_norm_BDTcut->SetFillColorAlpha(kCyan+1, 0.25);

    // h_BuDpDmKp_norm_BDTcut->SetLineColor(kGreen+1);
    // h_BuDpDmKp_norm_BDTcut->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuDsDsKp_norm_BDTcut->SetLineColor(kViolet+1);
    // h_BuDsDsKp_norm_BDTcut->SetFillColorAlpha(kViolet+1, 0.25);

    // h_BdDpD0Kp_norm_BDTcut->SetLineColor(kViolet+2);
    // h_BdDpD0Kp_norm_BDTcut->SetFillColorAlpha(kViolet+2, 0.25);

    // h_BsDsD0Kp_norm_BDTcut->SetLineColor(kPink+7);    
    // h_BsDsD0Kp_norm_BDTcut->SetFillColorAlpha(kPink+7, 0.25);

    // h_BuD0DpK0_norm_BDTcut->SetLineColor(kRed);
    // h_BuD0DpK0_norm_BDTcut->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0Ds_norm_BDTcut->SetLineColor(kOrange-3);
    // h_BuD0Ds_norm_BDTcut->SetFillColorAlpha(kOrange-3, 0.25);

    // h_BuD0Dp_norm_BDTcut->SetLineColor(kOrange+4);
    // h_BuD0Dp_norm_BDTcut->SetFillColorAlpha(kOrange+4, 0.25);

    // TCanvas d7;
    // d7.cd();
    // h_3pi3pi_norm_BDTcut->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_3pi3pi_norm_BDTcut->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_3pi3pi_norm_BDTcut->SetTitle("");
    // h_3pi3pi_norm_BDTcut->DrawNormalized();
    // h_all_mc_norm_BDTcut->DrawNormalized("same");
    // h_BuD0D0Kp_norm_BDTcut->DrawNormalized("same");
    // h_BuDpDmKp_norm_BDTcut->DrawNormalized("same");
    // h_BuDsDsKp_norm_BDTcut->DrawNormalized("same");
    // h_BdDpD0Kp_norm_BDTcut->DrawNormalized("same");
    // h_BsDsD0Kp_norm_BDTcut->DrawNormalized("same");
    // h_BuD0DpK0_norm_BDTcut->DrawNormalized("same");
    // h_BuD0Ds_norm_BDTcut->DrawNormalized("same");
    // h_BuD0Dp_norm_BDTcut->DrawNormalized("same");
    // leg3->Draw("same");
    // d7.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_best_bmass_BDTcut.pdf");

    // // Ktautau MC: 3 MC components
    // TCanvas d8;
    // d8.cd();
    // h_ktautau_mc_3pi3pi->SetTitle("B^{+} mass");
    // h_ktautau_mc_3pi3pipi0->SetLineColor(kOrange-3);
    // h_ktautau_mc_3pi3pipi0->SetFillColorAlpha(kOrange-3, 0.25);

    // h_ktautau_mc_3pi3pi2pi0->SetLineColor(kRed);
    // h_ktautau_mc_3pi3pi2pi0->SetFillColorAlpha(kRed, 0.25);

    // h_ktautau_mc_3pi3pi->DrawNormalized();
    // h_ktautau_mc_3pi3pipi0->DrawNormalized("same");
    // h_ktautau_mc_3pi3pi2pi0->DrawNormalized("same");

    // TLegend* leg4 = new TLegend(0.7,0.7,0.9,0.9);
    // leg4->AddEntry(h_ktautau_mc_3pi3pi, "3#pi3#pi MC");
    // leg4->AddEntry(h_ktautau_mc_3pi3pipi0, "3#pi3#pi #pi^{0} MC");
    // leg4->AddEntry(h_ktautau_mc_3pi3pi2pi0, "3#pi3#pi 2#pi^{0} MC");
    // leg4->Draw("same");
    // d8.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_comp.pdf");

    // TH1D* h_3pi3pi_iso = new TH1D("h_3pi3pi_iso", "h_3pi3pi_iso", 100, 0, 1);
    // TH1D* h_3pi3pipi0_iso = new TH1D("h_3pi3pipi0_iso", "h_3pi3pipi0_iso", 100, 0, 1);
    // TH1D* h_3pi3pi2pi0_iso = new TH1D("h_3pi3pi2pi0_iso", "h_3pi3pi2pi0_iso", 100, 0, 1);

    // t_ktautau_mc->Draw("BDT1 >> h_3pi3pi_iso", "(component==0) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("BDT1 >> h_3pi3pipi0_iso", "(component==1) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("BDT1 >> h_3pi3pi2pi0_iso", "(component==2) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);

    // h_3pi3pi_iso->GetXaxis()->SetTitle("BDT1");
    // h_3pi3pi_iso->GetYaxis()->SetTitle("Normalized entries / (0.01)");
    // h_3pi3pi_iso->SetTitle("Isolation MVA");

    // h_3pi3pi_iso->SetLineColor(kBlue);
    // h_3pi3pi_iso->SetFillColorAlpha(kBlue, 0.25);
    
    // h_3pi3pipi0_iso->SetLineColor(kOrange-3);
    // h_3pi3pipi0_iso->SetFillColorAlpha(kOrange-3, 0.25);

    // h_3pi3pi2pi0_iso->SetLineColor(kRed);
    // h_3pi3pi2pi0_iso->SetFillColorAlpha(kRed, 0.25);
    
    // TCanvas d9;
    // d9.cd();
    // h_3pi3pi_iso->DrawNormalized();
    // h_3pi3pipi0_iso->DrawNormalized("same");
    // h_3pi3pi2pi0_iso->DrawNormalized("same");
    // TLegend* leg5 = new TLegend(0.4,0.7,0.6,0.9);
    // leg5->AddEntry(h_3pi3pi_iso, "3#pi3#pi MC");
    // leg5->AddEntry(h_3pi3pipi0_iso, "3#pi3#pi #pi^{0} MC");
    // leg5->AddEntry(h_3pi3pi2pi0_iso, "3#pi3#pi 2#pi^{0} MC");
    // leg5->Draw("same");
    // d9.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_comp_bdt1.pdf");

    // TH1D* h_3pi3pi_kin = new TH1D("h_3pi3pi_kin", "h_3pi3pi_kin", 100, 0, 1);
    // TH1D* h_3pi3pipi0_kin = new TH1D("h_3pi3pipi0_kin", "h_3pi3pipi0_kin", 100, 0, 1);
    // TH1D* h_3pi3pi2pi0_kin = new TH1D("h_3pi3pi2pi0_kin", "h_3pi3pi2pi0_kin", 100, 0, 1);

    // t_ktautau_mc->Draw("BDT2 >> h_3pi3pi_kin", "(component==0) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("BDT2 >> h_3pi3pipi0_kin", "(component==1) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);
    // t_ktautau_mc->Draw("BDT2 >> h_3pi3pi2pi0_kin", "(component==2) && (is_best_cand == 1)"+ktautau_mc_post_sel_cuts);

    // h_3pi3pi_kin->GetXaxis()->SetTitle("BDT2");
    // h_3pi3pi_kin->GetYaxis()->SetTitle("Normalized entries / (0.01)");
    // h_3pi3pi_kin->SetTitle("Topological MVA");

    // h_3pi3pi_kin->SetLineColor(kBlue);
    // h_3pi3pi_kin->SetFillColorAlpha(kBlue, 0.25);
    
    // h_3pi3pipi0_kin->SetLineColor(kOrange-3);
    // h_3pi3pipi0_kin->SetFillColorAlpha(kOrange-3, 0.25);

    // h_3pi3pi2pi0_kin->SetLineColor(kRed);
    // h_3pi3pi2pi0_kin->SetFillColorAlpha(kRed, 0.25);
    
    // TCanvas d10;
    // d10.cd();
    // h_3pi3pi_kin->DrawNormalized();
    // h_3pi3pipi0_kin->DrawNormalized("same");
    // h_3pi3pi2pi0_kin->DrawNormalized("same");
    // leg5->Draw("same");
    // d10.SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/ktautau_mc_comp_bdt2.pdf");

    // // Cocktail MCs
    // auto e1 = new TCanvas(); 
    // auto e2 = new TCanvas(); 
    // auto e3 = new TCanvas(); 
    // auto e4 = new TCanvas(); 
    // auto e5 = new TCanvas(); 
    // auto e6 = new TCanvas();
    // auto e7 = new TCanvas(); 
    // auto e8 = new TCanvas(); 

    // auto EightPads = new TCanvas("EightPads","EightPads");
    // EightPads->Divide(2,4);
    // EightPads->cd(1) ; e1->DrawClonePad();

    // TH1D* h_BuD0D0Kp_DD = new TH1D("h_BuD0D0Kp_DD", "h_BuD0D0Kp_DD", 100, 4000, 8000);
    // TH1D* h_BuD0D0Kp_DsD = new TH1D("h_BuD0D0Kp_DsD", "h_BuD0D0Kp_DsD", 100, 4000, 8000);
    // TH1D* h_BuD0D0Kp_DDs = new TH1D("h_BuD0D0Kp_DDs", "h_BuD0D0Kp_DDs", 100, 4000, 8000);
    // TH1D* h_BuD0D0Kp_DsDs = new TH1D("h_BuD0D0Kp_DsDs", "h_BuD0D0Kp_DsDs", 100, 4000, 8000);

    // t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DD", "(species==100) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DsD", "(species==100) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DDs", "(species==100) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuD0D0Kp_DsDs", "(species==100) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0D0Kp_DsDs->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BuD0D0Kp_DsDs->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BuD0D0Kp_DsDs->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+} MC");

    // h_BuD0D0Kp_DD->SetLineColor(kBlue);
    // h_BuD0D0Kp_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BuD0D0Kp_DsD->SetLineColor(kBlack);
    // h_BuD0D0Kp_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BuD0D0Kp_DDs->SetLineColor(kRed);
    // h_BuD0D0Kp_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0D0Kp_DsDs->SetLineColor(kGreen+1);
    // h_BuD0D0Kp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuD0D0Kp_DsDs->DrawNormalized();
    // h_BuD0D0Kp_DsD->DrawNormalized("same");
    // h_BuD0D0Kp_DDs->DrawNormalized("same");
    // h_BuD0D0Kp_DD->DrawNormalized("same");

    // TLegend* leg6 = new TLegend(0.8,0.6,0.9,0.9);
    // leg6->AddEntry(h_BuD0D0Kp_DD, "DD");
    // leg6->AddEntry(h_BuD0D0Kp_DsD, "Dstar D");
    // leg6->AddEntry(h_BuD0D0Kp_DDs, "D Dstar");
    // leg6->AddEntry(h_BuD0D0Kp_DsDs, "Dstar Dstar");
    // leg6->Draw("same");

    // EightPads->cd(2) ; e2->DrawClonePad();

    // TH1D* h_BuDpDmKp_DD = new TH1D("h_BuDpDmKp_DD", "h_BuDpDmKp_DD", 100, 4000, 8000);
    // TH1D* h_BuDpDmKp_DsD = new TH1D("h_BuDpDmKp_DsD", "h_BuDpDmKp_DsD", 100, 4000, 8000);
    // TH1D* h_BuDpDmKp_DDs = new TH1D("h_BuDpDmKp_DDs", "h_BuDpDmKp_DDs", 100, 4000, 8000);
    // TH1D* h_BuDpDmKp_DsDs = new TH1D("h_BuDpDmKp_DsDs", "h_BuDpDmKp_DsDs", 100, 4000, 8000);

    // t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DD", "(species==101) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DsD", "(species==101) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DDs", "(species==101) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDpDmKp_DsDs", "(species==101) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuDpDmKp_DsD->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BuDpDmKp_DsD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BuDpDmKp_DsD->SetTitle("B^{+} #rightarrow D^{+} D^{-} K^{+} MC");

    // h_BuDpDmKp_DD->SetLineColor(kBlue);
    // h_BuDpDmKp_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BuDpDmKp_DsD->SetLineColor(kBlack);
    // h_BuDpDmKp_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BuDpDmKp_DDs->SetLineColor(kRed);
    // h_BuDpDmKp_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BuDpDmKp_DsDs->SetLineColor(kGreen+1);
    // h_BuDpDmKp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuDpDmKp_DsD->DrawNormalized();
    // h_BuDpDmKp_DD->DrawNormalized("same");
    // h_BuDpDmKp_DDs->DrawNormalized("same");
    // h_BuDpDmKp_DsDs->DrawNormalized("same");
    // leg6->Draw("same");
    
    // EightPads->cd(3) ; e3->DrawClonePad();

    // TH1D* h_BuDsDsKp_DD = new TH1D("h_BuDsDsKp_DD", "h_BuDsDsKp_DD", 100, 4000, 8000);
    // TH1D* h_BuDsDsKp_DsD = new TH1D("h_BuDsDsKp_DsD", "h_BuDsDsKp_DsD", 100, 4000, 8000);
    // TH1D* h_BuDsDsKp_DDs = new TH1D("h_BuDsDsKp_DDs", "h_BuDsDsKp_DDs", 100, 4000, 8000);
    // TH1D* h_BuDsDsKp_DsDs = new TH1D("h_BuDsDsKp_DsDs", "h_BuDsDsKp_DsDs", 100, 4000, 8000);

    // t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DD", "(species==102) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DsD", "(species==102) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DDs", "(species==102) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("df_Bp_M >> h_BuDsDsKp_DsDs", "(species==102) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuDsDsKp_DDs->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BuDsDsKp_DDs->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BuDsDsKp_DDs->SetTitle("B^{+} #rightarrow D^{+}_{s} D^{-}_{s} K^{+} MC");

    // h_BuDsDsKp_DD->SetLineColor(kBlue);
    // h_BuDsDsKp_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BuDsDsKp_DsD->SetLineColor(kBlack);
    // h_BuDsDsKp_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BuDsDsKp_DDs->SetLineColor(kRed);
    // h_BuDsDsKp_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BuDsDsKp_DsDs->SetLineColor(kGreen+1);
    // h_BuDsDsKp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuDsDsKp_DDs->DrawNormalized();
    // h_BuDsDsKp_DsD->DrawNormalized("same");
    // h_BuDsDsKp_DD->DrawNormalized("same");
    // h_BuDsDsKp_DsDs->DrawNormalized("same");
    // leg6->Draw("same");

    // EightPads->cd(4) ; e4->DrawClonePad();

    // TH1D* h_BdDpD0Kp_DD = new TH1D("h_BdDpD0Kp_DD", "h_BdDpD0Kp_DD", 100, 4000, 8000);
    // TH1D* h_BdDpD0Kp_DsD = new TH1D("h_BdDpD0Kp_DsD", "h_BdDpD0Kp_DsD", 100, 4000, 8000);
    // TH1D* h_BdDpD0Kp_DDs = new TH1D("h_BdDpD0Kp_DDs", "h_BdDpD0Kp_DDs", 100, 4000, 8000);
    // TH1D* h_BdDpD0Kp_DsDs = new TH1D("h_BdDpD0Kp_DsDs", "h_BdDpD0Kp_DsDs", 100, 4000, 8000);

    // t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DD", "(species==110) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DsD", "(species==110) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DDs", "(species==110) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BdDDKp->Draw("df_Bp_M >> h_BdDpD0Kp_DsDs", "(species==110) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BdDpD0Kp_DD->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BdDpD0Kp_DD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BdDpD0Kp_DD->SetTitle("B^{0} #rightarrow D^{-} D^{0} K^{+} MC");

    // h_BdDpD0Kp_DD->SetLineColor(kBlue);
    // h_BdDpD0Kp_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BdDpD0Kp_DsD->SetLineColor(kBlack);
    // h_BdDpD0Kp_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BdDpD0Kp_DDs->SetLineColor(kRed);
    // h_BdDpD0Kp_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BdDpD0Kp_DsDs->SetLineColor(kGreen+1);
    // h_BdDpD0Kp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BdDpD0Kp_DD->DrawNormalized();
    // h_BdDpD0Kp_DsD->DrawNormalized("same");
    // h_BdDpD0Kp_DDs->DrawNormalized("same");
    // h_BdDpD0Kp_DsDs->DrawNormalized("same");
    // leg6->Draw("same");

    // EightPads->cd(5) ; e5->DrawClonePad();

    // TH1D* h_BsDpD0Kp_DD = new TH1D("h_BsDpD0Kp_DD", "h_BsDpD0Kp_DD", 100, 4000, 8000);
    // TH1D* h_BsDpD0Kp_DsD = new TH1D("h_BsDpD0Kp_DsD", "h_BsDpD0Kp_DsD", 100, 4000, 8000);
    // TH1D* h_BsDpD0Kp_DDs = new TH1D("h_BsDpD0Kp_DDs", "h_BsDpD0Kp_DDs", 100, 4000, 8000);
    // TH1D* h_BsDpD0Kp_DsDs = new TH1D("h_BsDpD0Kp_DsDs", "h_BsDpD0Kp_DsDs", 100, 4000, 8000);

    // t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DD", "(species==120) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DsD", "(species==120) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DDs", "(species==120) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BsDDKp->Draw("df_Bp_M >> h_BsDpD0Kp_DsDs", "(species==120) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BsDpD0Kp_DsDs->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BsDpD0Kp_DsDs->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BsDpD0Kp_DsDs->SetTitle("B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+} MC");

    // h_BsDpD0Kp_DD->SetLineColor(kBlue);
    // h_BsDpD0Kp_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BsDpD0Kp_DsD->SetLineColor(kBlack);
    // h_BsDpD0Kp_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BsDpD0Kp_DDs->SetLineColor(kRed);
    // h_BsDpD0Kp_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BsDpD0Kp_DsDs->SetLineColor(kGreen+1);
    // h_BsDpD0Kp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BsDpD0Kp_DsDs->DrawNormalized();
    // h_BsDpD0Kp_DsD->DrawNormalized("same");
    // h_BsDpD0Kp_DDs->DrawNormalized("same");
    // h_BsDpD0Kp_DD->DrawNormalized("same");
    // leg6->Draw("same");

    // EightPads->cd(6) ; e6->DrawClonePad();

    // TH1D* h_BuD0DpK0_DD = new TH1D("h_BuD0DpK0_DD", "h_BuD0DpK0_DD", 100, 4000, 8000);
    // TH1D* h_BuD0DpK0_DsD = new TH1D("h_BuD0DpK0_DsD", "h_BuD0DpK0_DsD", 100, 4000, 8000);
    // TH1D* h_BuD0DpK0_DDs = new TH1D("h_BuD0DpK0_DDs", "h_BuD0DpK0_DDs", 100, 4000, 8000);
    // TH1D* h_BuD0DpK0_DsDs = new TH1D("h_BuD0DpK0_DsDs", "h_BuD0DpK0_DsDs", 100, 4000, 8000);

    // t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DD", "(species==130) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DsD", "(species==130) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DDs", "(species==130) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDK0->Draw("df_Bp_M >> h_BuD0DpK0_DsDs", "(species==130) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0DpK0_DD->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BuD0DpK0_DD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BuD0DpK0_DD->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0} MC");

    // h_BuD0DpK0_DD->SetLineColor(kBlue);
    // h_BuD0DpK0_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BuD0DpK0_DsD->SetLineColor(kBlack);
    // h_BuD0DpK0_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BuD0DpK0_DDs->SetLineColor(kRed);
    // h_BuD0DpK0_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0DpK0_DsDs->SetLineColor(kGreen+1);
    // h_BuD0DpK0_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuD0DpK0_DD->DrawNormalized();
    // h_BuD0DpK0_DsD->DrawNormalized("same");
    // h_BuD0DpK0_DDs->DrawNormalized("same");
    // h_BuD0DpK0_DsDs->DrawNormalized("same");
    // leg6->Draw("same");

    // EightPads->cd(7) ; e7->DrawClonePad();

    // TH1D* h_BuD0Ds_DD = new TH1D("h_BuD0Ds_DD", "h_BuD0Ds_DD", 100, 4000, 8000);
    // TH1D* h_BuD0Ds_DsD = new TH1D("h_BuD0Ds_DsD", "h_BuD0Ds_DsD", 100, 4000, 8000);
    // TH1D* h_BuD0Ds_DDs = new TH1D("h_BuD0Ds_DDs", "h_BuD0Ds_DDs", 100, 4000, 8000);
    // TH1D* h_BuD0Ds_DsDs = new TH1D("h_BuD0Ds_DsDs", "h_BuD0Ds_DsDs", 100, 4000, 8000);

    // t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DD", "(species==150) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DsD", "(species==150) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DDs", "(species==150) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Ds_DsDs", "(species==150) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0Ds_DsD->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BuD0Ds_DsD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BuD0Ds_DsD->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{+}_{s} MC");

    // h_BuD0Ds_DD->SetLineColor(kBlue);
    // h_BuD0Ds_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BuD0Ds_DsD->SetLineColor(kBlack);
    // h_BuD0Ds_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BuD0Ds_DDs->SetLineColor(kRed);
    // h_BuD0Ds_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0Ds_DsDs->SetLineColor(kGreen+1);
    // h_BuD0Ds_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuD0Ds_DsD->DrawNormalized();
    // h_BuD0Ds_DD->DrawNormalized("same");
    // h_BuD0Ds_DDs->DrawNormalized("same");
    // h_BuD0Ds_DsDs->DrawNormalized("same");
    // leg6->Draw("same");

    // EightPads->cd(8) ; e8->DrawClonePad();

    // TH1D* h_BuD0Dp_DD = new TH1D("h_BuD0Dp_DD", "h_BuD0Dp_DD", 100, 4000, 8000);
    // TH1D* h_BuD0Dp_DsD = new TH1D("h_BuD0Dp_DsD", "h_BuD0Dp_DsD", 100, 4000, 8000);
    // TH1D* h_BuD0Dp_DDs = new TH1D("h_BuD0Dp_DDs", "h_BuD0Dp_DDs", 100, 4000, 8000);
    // TH1D* h_BuD0Dp_DsDs = new TH1D("h_BuD0Dp_DsDs", "h_BuD0Dp_DsDs", 100, 4000, 8000);

    // t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DD", "(species==151) && (component==0) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DsD", "(species==151) && (component==1) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DDs", "(species==151) && (component==2) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("df_Bp_M >> h_BuD0Dp_DsDs", "(species==151) && (component==3) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0Dp_DD->GetXaxis()->SetTitle("m_{B} (MeV)");
    // h_BuD0Dp_DD->GetYaxis()->SetTitle("Normalized entries / (40 MeV)");
    // h_BuD0Dp_DD->SetTitle("B^{+} #rightarrow #bar{D}^{0} D^{+} MC");

    // h_BuD0Dp_DD->SetLineColor(kBlue);
    // h_BuD0Dp_DD->SetFillColorAlpha(kBlue, 0.25);

    // h_BuD0Dp_DsD->SetLineColor(kBlack);
    // h_BuD0Dp_DsD->SetFillColorAlpha(kBlack, 0.25);

    // h_BuD0Dp_DDs->SetLineColor(kRed);
    // h_BuD0Dp_DDs->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0Dp_DsDs->SetLineColor(kGreen+1);
    // h_BuD0Dp_DsDs->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuD0Dp_DD->DrawNormalized();
    // h_BuD0Dp_DsD->DrawNormalized("same");
    // h_BuD0Dp_DDs->DrawNormalized("same");
    // h_BuD0Dp_DsDs->DrawNormalized("same");
    // leg6->Draw("same");

    // EightPads->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_DD_components.pdf");

    // auto b1 = new TCanvas(); 
    // auto b2 = new TCanvas(); 
    // auto b3 = new TCanvas(); 
    // auto b4 = new TCanvas(); 
    // auto b5 = new TCanvas(); 

    // auto SixPads = new TCanvas("SixPads","SixPads");
    // SixPads->Divide(2,3);
    // SixPads->cd(1) ; b1->DrawClonePad();

    // TH1D* h_BuD0D0Kp_iso = new TH1D("h_BuD0D0Kp_iso", "h_BuD0D0Kp_iso", 100, 0, 1);
    // TH1D* h_BuDpDmKp_iso = new TH1D("h_BuDpDmKp_iso", "h_BuDpDmKp_iso", 100, 0, 1);
    // TH1D* h_BuDsDsKp_iso = new TH1D("h_BuDsDsKp_iso", "h_BuDsDsKp_iso", 100, 0, 1);

    // t_BuDDKp->Draw("BDT1 >> h_BuD0D0Kp_iso", "(species==100) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("BDT1 >> h_BuDpDmKp_iso", "(species==101) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("BDT1 >> h_BuDsDsKp_iso", "(species==102) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuDpDmKp_iso->GetXaxis()->SetTitle("BDT1");
    // h_BuDpDmKp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BuDpDmKp_iso->SetTitle("B^{+} #rightarrow D D K^{+} MC");

    // h_BuD0D0Kp_iso->SetLineColor(kCyan+1);
    // h_BuD0D0Kp_iso->SetFillColorAlpha(kCyan+1, 0.25);

    // h_BuDpDmKp_iso->SetLineColor(kGreen+1);
    // h_BuDpDmKp_iso->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuDsDsKp_iso->SetLineColor(kViolet+1);
    // h_BuDsDsKp_iso->SetFillColorAlpha(kViolet+1, 0.25);

    // h_BuDpDmKp_iso->DrawNormalized();
    // h_BuD0D0Kp_iso->DrawNormalized("same");
    // h_BuDsDsKp_iso->DrawNormalized("same");

    // TLegend* leg7 = new TLegend(0.35,0.6,0.65,0.9);
    // leg7->AddEntry(h_BuD0D0Kp_iso, "B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+}");
    // leg7->AddEntry(h_BuDpDmKp_iso, "B^{+} #rightarrow D^{+} D^{-} K^{+}");
    // leg7->AddEntry(h_BuDsDsKp_iso, "B^{+} #rightarrow D^{+}_{s} D^{-}_{s} K^{+}");
    // leg7->Draw("same");

    // SixPads->cd(2) ; b2->DrawClonePad();

    // TH1D* h_BdDpD0Kp_iso = new TH1D("h_BdDpD0Kp_iso", "h_BdDpD0Kp_iso", 100, 0, 1);

    // t_BdDDKp->Draw("BDT1 >> h_BdDpD0Kp_iso", "(species==110) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BdDpD0Kp_iso->GetXaxis()->SetTitle("BDT1");
    // h_BdDpD0Kp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BdDpD0Kp_iso->SetTitle("B^{0} #rightarrow D D K^{+} MC");

    // h_BdDpD0Kp_iso->SetLineColor(kViolet+2);
    // h_BdDpD0Kp_iso->SetFillColorAlpha(kViolet+2, 0.25);

    // h_BdDpD0Kp_iso->DrawNormalized();

    // TLegend* leg8 = new TLegend(0.35,0.8,0.65,0.9);
    // leg8->AddEntry(h_BdDpD0Kp_iso, "B^{0} #rightarrow D^{-} D^{0} K^{+}");
    // leg8->Draw("same");

    // SixPads->cd(3) ; b3->DrawClonePad();

    // TH1D* h_BsDsD0Kp_iso = new TH1D("h_BsDsD0Kp_iso", "h_BsDsD0Kp_iso", 100, 0, 1);

    // t_BsDDKp->Draw("BDT1 >> h_BsDsD0Kp_iso", "(species==120) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BsDsD0Kp_iso->GetXaxis()->SetTitle("BDT1");
    // h_BsDsD0Kp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BsDsD0Kp_iso->SetTitle("B^{0}_{s} #rightarrow D D K^{+} MC");

    // h_BsDsD0Kp_iso->SetLineColor(kPink+7);
    // h_BsDsD0Kp_iso->SetFillColorAlpha(kPink+7, 0.25);

    // h_BsDsD0Kp_iso->DrawNormalized();

    // TLegend* leg9 = new TLegend(0.35,0.8,0.65,0.9);
    // leg9->AddEntry(h_BsDsD0Kp_iso, "B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+}");
    // leg9->Draw("same");

    // SixPads->cd(4) ; b4->DrawClonePad();

    // TH1D* h_BuD0DpK0_iso = new TH1D("h_BuD0DpK0_iso", "h_BuD0DpK0_iso", 100, 0, 1);

    // t_BuDDK0->Draw("BDT1 >> h_BuD0DpK0_iso", "(species==130) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0DpK0_iso->GetXaxis()->SetTitle("BDT1");
    // h_BuD0DpK0_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BuD0DpK0_iso->SetTitle("B^{+} #rightarrow D D K^{0} MC");

    // h_BuD0DpK0_iso->SetLineColor(kRed);
    // h_BuD0DpK0_iso->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0DpK0_iso->DrawNormalized();

    // TLegend* leg10 = new TLegend(0.35,0.8,0.65,0.9);
    // leg10->AddEntry(h_BuD0DpK0_iso, "B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0}");
    // leg10->Draw("same");

    // SixPads->cd(5) ; b5->DrawClonePad();

    // TH1D* h_BuD0Ds_iso = new TH1D("h_BuD0Ds_iso", "h_BuD0Ds_iso", 100, 0, 1);
    // TH1D* h_BuD0Dp_iso = new TH1D("h_BuD0Dp_iso", "h_BuD0Dp_iso", 100, 0, 1);

    // t_BuDD->Draw("BDT1 >> h_BuD0Ds_iso", "(species==150) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("BDT1 >> h_BuD0Dp_iso", "(species==151) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0Dp_iso->GetXaxis()->SetTitle("BDT1");
    // h_BuD0Dp_iso->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BuD0Dp_iso->SetTitle("B^{+} #rightarrow D D MC");

    // h_BuD0Ds_iso->SetLineColor(kOrange-3);
    // h_BuD0Ds_iso->SetFillColorAlpha(kOrange-3, 0.25);

    // h_BuD0Dp_iso->SetLineColor(kOrange+3);
    // h_BuD0Dp_iso->SetFillColorAlpha(kOrange+3, 0.25);

    // h_BuD0Dp_iso->DrawNormalized();
    // h_BuD0Ds_iso->DrawNormalized("same");

    // TLegend* leg11 = new TLegend(0.35,0.7,0.65,0.9);
    // leg11->AddEntry(h_BuD0Ds_iso, "B^{+} #rightarrow #bar{D}^{0} D^{+}_{s}");
    // leg11->AddEntry(h_BuD0Dp_iso, "B^{+} #rightarrow #bar{D}^{0} D^{+}");
    // leg11->Draw("same");
    // SixPads->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_BDT1.pdf");

    // auto a1 = new TCanvas(); 
    // auto a2 = new TCanvas(); 
    // auto a3 = new TCanvas(); 
    // auto a4 = new TCanvas(); 
    // auto a5 = new TCanvas(); 

    // auto SixPads1 = new TCanvas("SixPads1","SixPads1");
    // SixPads1->Divide(2,3);
    // SixPads1->cd(1) ; a1->DrawClonePad();

    // TH1D* h_BuD0D0Kp_kin = new TH1D("h_BuD0D0Kp_kin", "h_BuD0D0Kp_kin", 100, 0, 1);
    // TH1D* h_BuDpDmKp_kin = new TH1D("h_BuDpDmKp_kin", "h_BuDpDmKp_kin", 100, 0, 1);
    // TH1D* h_BuDsDsKp_kin = new TH1D("h_BuDsDsKp_kin", "h_BuDsDsKp_kin", 100, 0, 1);

    // t_BuDDKp->Draw("BDT2 >> h_BuD0D0Kp_kin", "(species==100) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("BDT2 >> h_BuDpDmKp_kin", "(species==101) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDDKp->Draw("BDT2 >> h_BuDsDsKp_kin", "(species==102) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuDsDsKp_kin->GetXaxis()->SetTitle("BDT2");
    // h_BuDsDsKp_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BuDsDsKp_kin->SetTitle("B^{+} #rightarrow D D K^{+} MC");

    // h_BuD0D0Kp_kin->SetLineColor(kCyan+1);
    // h_BuD0D0Kp_kin->SetFillColorAlpha(kCyan+1, 0.25);

    // h_BuDpDmKp_kin->SetLineColor(kGreen+1);
    // h_BuDpDmKp_kin->SetFillColorAlpha(kGreen+1, 0.25);

    // h_BuDsDsKp_kin->SetLineColor(kViolet+1);
    // h_BuDsDsKp_kin->SetFillColorAlpha(kViolet+1, 0.25);

    // h_BuDsDsKp_kin->DrawNormalized();
    // h_BuDpDmKp_kin->DrawNormalized("same");
    // h_BuD0D0Kp_kin->DrawNormalized("same");
    // leg7->Draw("same");

    // SixPads1->cd(2) ; a2->DrawClonePad();

    // TH1D* h_BdDpD0Kp_kin = new TH1D("h_BdDpD0Kp_kin", "h_BdDpD0Kp_kin", 100, 0, 1);

    // t_BdDDKp->Draw("BDT2 >> h_BdDpD0Kp_kin", "(species==110) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BdDpD0Kp_kin->GetXaxis()->SetTitle("BDT2");
    // h_BdDpD0Kp_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BdDpD0Kp_kin->SetTitle("B^{0} #rightarrow D D K^{+} MC");

    // h_BdDpD0Kp_kin->SetLineColor(kViolet+2);
    // h_BdDpD0Kp_kin->SetFillColorAlpha(kViolet+2, 0.25);

    // h_BdDpD0Kp_kin->DrawNormalized();
    // leg8->Draw("same");

    // SixPads1->cd(3) ; a3->DrawClonePad();

    // TH1D* h_BsDsD0Kp_kin = new TH1D("h_BsDsD0Kp_kin", "h_BsDsD0Kp_kin", 100, 0, 1);

    // t_BsDDKp->Draw("BDT2 >> h_BsDsD0Kp_kin", "(species==120) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BsDsD0Kp_kin->GetXaxis()->SetTitle("BDT2");
    // h_BsDsD0Kp_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BsDsD0Kp_kin->SetTitle("B^{0}_{s} #rightarrow D D K^{+} MC");

    // h_BsDsD0Kp_kin->SetLineColor(kPink+7);
    // h_BsDsD0Kp_kin->SetFillColorAlpha(kPink+7, 0.25);

    // h_BsDsD0Kp_kin->DrawNormalized();
    // leg9->Draw("same");

    // SixPads1->cd(4) ; a4->DrawClonePad();

    // TH1D* h_BuD0DpK0_kin = new TH1D("h_BuD0DpK0_kin", "h_BuD0DpK0_kin", 100, 0, 1);

    // t_BuDDK0->Draw("BDT2 >> h_BuD0DpK0_kin", "(species==130) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0DpK0_kin->GetXaxis()->SetTitle("BDT2");
    // h_BuD0DpK0_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BuD0DpK0_kin->SetTitle("B^{+} #rightarrow D D K^{0} MC");

    // h_BuD0DpK0_kin->SetLineColor(kRed);
    // h_BuD0DpK0_kin->SetFillColorAlpha(kRed, 0.25);

    // h_BuD0DpK0_kin->DrawNormalized();
    // leg10->Draw("same");

    // SixPads1->cd(5) ; a5->DrawClonePad();

    // TH1D* h_BuD0Ds_kin = new TH1D("h_BuD0Ds_kin", "h_BuD0Ds_kin", 100, 0, 1);
    // TH1D* h_BuD0Dp_kin = new TH1D("h_BuD0Dp_kin", "h_BuD0Dp_kin", 100, 0, 1);

    // t_BuDD->Draw("BDT2 >> h_BuD0Ds_kin", "(species==150) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);
    // t_BuDD->Draw("BDT2 >> h_BuD0Dp_kin", "(species==151) && (is_best_cand==1)"+ktautau_mc_post_sel_cuts);

    // h_BuD0Ds_kin->GetXaxis()->SetTitle("BDT2");
    // h_BuD0Ds_kin->GetYaxis()->SetTitle("Normalized entries / (0.01 MeV)");
    // h_BuD0Ds_kin->SetTitle("B^{+} #rightarrow D D MC");

    // h_BuD0Ds_kin->SetLineColor(kOrange-3);
    // h_BuD0Ds_kin->SetFillColorAlpha(kOrange-3, 0.25);

    // h_BuD0Dp_kin->SetLineColor(kOrange+3);
    // h_BuD0Dp_kin->SetFillColorAlpha(kOrange+3, 0.25);

    // h_BuD0Ds_kin->DrawNormalized();
    // h_BuD0Dp_kin->DrawNormalized("same");
    // leg11->Draw("same");

    // SixPads1->SaveAs("/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/cocktail_mcs_BDT2.pdf");