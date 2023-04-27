
void compare_BDTs(int year, int tmva_type, int background_proxy, int num_data_dst_files, TString MC_files, TString MC_0_files, TString MC_1_files, TString MC_2_files, TString RS_DATA_files, TString WS_DATA_files, TString MC_BDT_files, TString MC_0_BDT_files, TString MC_1_BDT_files, TString MC_2_BDT_files, TString RS_DATA_BDT_files, TString WS_DATA_BDT_files)
{
    TFileCollection* fc_mc = new TFileCollection("MC", "MC", MC_files);
    TFileCollection* fc_mc_0 = new TFileCollection("MC_0", "MC_0", MC_0_files);
    TFileCollection* fc_mc_1 = new TFileCollection("MC_1", "MC_1", MC_1_files);
    TFileCollection* fc_mc_2 = new TFileCollection("MC_2", "MC_2", MC_2_files);
    TFileCollection* fc_rs_data = new TFileCollection("rs_data", "rs_data", RS_DATA_files, num_data_dst_files);
    TFileCollection* fc_ws_data = new TFileCollection("ws_data", "ws_data", WS_DATA_files, num_data_dst_files);

    TChain* t_mc = new TChain("DecayTree");
    TChain* t_mc_0 = new TChain("DecayTree");
    TChain* t_mc_1 = new TChain("DecayTree");
    TChain* t_mc_2 = new TChain("DecayTree");
    TChain* t_RS_data = new TChain("DecayTree");
    TChain* t_WS_data = new TChain("DecayTree");

    t_mc->AddFileInfoList((TCollection*)fc_mc->GetList());
    t_mc_0->AddFileInfoList((TCollection*)fc_mc_0->GetList());
    t_mc_1->AddFileInfoList((TCollection*)fc_mc_1->GetList());
    t_mc_2->AddFileInfoList((TCollection*)fc_mc_2->GetList());
    t_RS_data->AddFileInfoList((TCollection*)fc_rs_data->GetList());
    t_WS_data->AddFileInfoList((TCollection*)fc_ws_data->GetList());

    TFileCollection* fc_mc_bdt = new TFileCollection("MC_BDT", "MC_BDT", MC_BDT_files);
    TFileCollection* fc_mc_0_bdt = new TFileCollection("MC_0_BDT", "MC_0_BDT", MC_0_BDT_files);
    TFileCollection* fc_mc_1_bdt = new TFileCollection("MC_1_BDT", "MC_1_BDT", MC_1_BDT_files);
    TFileCollection* fc_mc_2_bdt = new TFileCollection("MC_2_BDT", "MC_2_BDT", MC_2_BDT_files);
    TFileCollection* fc_rsdata_bdt = new TFileCollection("RS_DATA_BDT", "RS_DATA_BDT", RS_DATA_BDT_files, num_data_dst_files);
    TFileCollection* fc_wsdata_bdt = new TFileCollection("WS_DATA_BDT", "WS_DATA_BDT", WS_DATA_BDT_files, num_data_dst_files);

    TChain* t_mc_bdt = new TChain("DecayTree");
    TChain* t_mc_0_bdt = new TChain("DecayTree");
    TChain* t_mc_1_bdt = new TChain("DecayTree");
    TChain* t_mc_2_bdt = new TChain("DecayTree");
    TChain* t_rsdata_bdt = new TChain("DecayTree");
    TChain* t_wsdata_bdt = new TChain("DecayTree");

    t_mc_bdt->AddFileInfoList((TCollection*)fc_mc_bdt->GetList());
    t_mc_0_bdt->AddFileInfoList((TCollection*)fc_mc_0_bdt->GetList());
    t_mc_1_bdt->AddFileInfoList((TCollection*)fc_mc_1_bdt->GetList());
    t_mc_2_bdt->AddFileInfoList((TCollection*)fc_mc_2_bdt->GetList());
    t_rsdata_bdt->AddFileInfoList((TCollection*)fc_rsdata_bdt->GetList());
    t_wsdata_bdt->AddFileInfoList((TCollection*)fc_wsdata_bdt->GetList());

    t_mc->AddFriend(t_mc_bdt);
    t_mc_0->AddFriend(t_mc_0_bdt);
    t_mc_1->AddFriend(t_mc_1_bdt);
    t_mc_2->AddFriend(t_mc_2_bdt);
    t_RS_data->AddFriend(t_rsdata_bdt);
    t_WS_data->AddFriend(t_wsdata_bdt);

    TCut sigMass = "(Bp_ConsBp_seq_12_M[0]>4700) && (Bp_ConsBp_seq_12_M[0]<5800)";

    // 3pi3pi MC vs WS data (signal region)
    TH1D* h_signal = new TH1D("h_signal", "h_signal", 100, -1, 1);
    TH1D* h_background = new TH1D("h_background", "h_background", 100, -1, 1);

    TH1D* h_3pi_3pi = new TH1D("h_3pi_3pi", "h_3pi_3pi", 100, -1, 1);
    TH1D* h_3pi_3pi_pi0 = new TH1D("h_3pi_3pi_pi0", "h_3pi_3pi_pi0", 100, -1, 1);
    TH1D* h_3pi_3pi_2pi0 = new TH1D("h_3pi_3pi_2pi0", "h_3pi_3pi_2pi0", 100, -1, 1);

    TH1D* h_RS = new TH1D("h_RS", "h_RS", 100, -1, 1);
    TH1D* h_WS = new TH1D("h_WS", "h_WS", 100, -1, 1);

    cout << "signal vs background" << endl;
    cout << "signal" << endl;
    t_mc_0->Draw("BDT_response >> h_signal", sigMass);
    cout << "background" << endl;
    t_WS_data->Draw("BDT_response >> h_background", sigMass);

    cout << "MC components" << endl;
    cout << "3pi 3pi" << endl;
    t_mc_0->Draw("BDT_response >> h_3pi_3pi");
    cout << "3pi 3pi pi0" << endl;
    t_mc_1->Draw("BDT_response >> h_3pi_3pi_pi0");
    cout << "3pi 3pi 2pi0" << endl;
    t_mc_2->Draw("BDT_response >> h_3pi_3pi_2pi0");

    cout << "RS vs WS data" << endl;
    cout << "RS data" << endl;
    t_RS_data->Draw("BDT_response >> h_RS");
    cout << "WS data" << endl;
    t_WS_data->Draw("BDT_response >> h_WS");

    gStyle->SetOptStat(0);
    TCanvas c;
    c.cd();
    h_background->GetXaxis()->SetTitle("BDT response");
    h_background->GetYaxis()->SetTitle("Normalised entries (100 bins)");
    h_background->SetTitle("(truthMatch)+passDTF+sigMass+trigger+others");
    h_signal->SetLineColor(kBlue);
    h_signal->SetFillColorAlpha(kBlue, 0.25);
    h_background->SetLineColor(kRed);
    h_background->SetFillColorAlpha(kRed, 0.25);
    h_background->DrawNormalized();
    h_signal->DrawNormalized("same");
    TLegend* leg;
    leg = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg->AddEntry(h_signal, "Signal", "f");
    leg->AddEntry(h_background, "Background", "f");
    leg->SetTextSize(0.03);
    leg->Draw("same");
    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201%i/sig_bkg_tmva_type%i_bkgProxy%i_numDST%i_comparison.gif",year,tmva_type,background_proxy,num_data_dst_files));
    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201%i/sig_bkg_tmva_type%i_bkgProxy%i_numDST%i_comparison.pdf",year,tmva_type,background_proxy,num_data_dst_files));

    TCanvas c1;
    c1.cd();
    h_3pi_3pi_2pi0->GetXaxis()->SetTitle("BDT response");
    h_3pi_3pi_2pi0->GetYaxis()->SetTitle("Normalised entries (100 bins)");
    h_3pi_3pi_2pi0->SetTitle("truthMatch+passDTF+trigger+others");
    h_3pi_3pi_2pi0->SetLineColor(kRed);
    h_3pi_3pi_2pi0->SetFillColorAlpha(kRed, 0.25);
    h_3pi_3pi_pi0->SetLineColor(kBlue);
    h_3pi_3pi_pi0->SetFillColorAlpha(kBlue, 0.25);
    h_3pi_3pi->SetLineColor(kBlack);
    h_3pi_3pi->SetFillColorAlpha(kBlack, 0.25);
    h_3pi_3pi_2pi0->DrawNormalized();
    h_3pi_3pi_pi0->DrawNormalized("same");
    h_3pi_3pi->DrawNormalized("same");
    cout << "3pi3pi " << h_3pi_3pi->GetEntries() << endl;
    cout << "3pi3pi pi0 " << h_3pi_3pi_pi0->GetEntries() << endl;
    cout << "3pi3pi 2pi0 " << h_3pi_3pi_2pi0->GetEntries() << endl;

    TLegend* leg1;
    leg1 = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg1->AddEntry(h_3pi_3pi, "3#pi 3#pi", "f");
    leg1->AddEntry(h_3pi_3pi_pi0, "3#pi 3#pi #pi^{0}", "f");
    leg1->AddEntry(h_3pi_3pi_2pi0, "3#pi 3#pi 2#pi^{0}", "f");
    leg1->SetTextSize(0.03);
    leg1->Draw("same");
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201%i/mc_components_tmva_type%i_bkgProxy%i_numDST%i_comparison.gif",year,tmva_type,background_proxy,num_data_dst_files));
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201%i/mc_components_tmva_type%i_bkgProxy%i_numDST%i_comparison.pdf",year,tmva_type,background_proxy,num_data_dst_files));

    TCanvas c2;
    c2.cd();
    h_RS->GetXaxis()->SetTitle("BDT response");
    h_RS->GetYaxis()->SetTitle("Normalised entries (100 bins)");
    h_RS->SetTitle("passDTF+trigger+others");
    h_RS->SetLineColor(kBlue);
    h_RS->SetFillColorAlpha(kBlue, 0.25);
    h_WS->SetLineColor(kRed);
    h_WS->SetFillColorAlpha(kRed, 0.25);
    h_RS->DrawNormalized();
    h_WS->DrawNormalized("same");
    TLegend* leg2;
    leg2 = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg2->AddEntry(h_RS, "RS data", "f");
    leg2->AddEntry(h_WS, "WS data", "f");
    leg2->SetTextSize(0.03);
    leg2->Draw("same");
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201%i/RS_WS_tmva_type%i_bkgProxy%i_numDST%i_comparison.gif",year,tmva_type,background_proxy,num_data_dst_files));
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201%i/RS_WS_tmva_type%i_bkgProxy%i_numDST%i_comparison.pdf",year,tmva_type,background_proxy,num_data_dst_files));

    return;
}