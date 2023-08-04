
void compare_fits()
{
    TCut df_status = "df_status==0";
    TCut passDTF = "Bp_dtf_12_status==0";

    // TTree from DaVinci
    TFileCollection* fc = new TFileCollection("MC", "MC", "/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/2018/Component_-1/pre_sel_tree.txt");
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    // TTree with results from the standalone fitter
    TFileCollection* fc1 = new TFileCollection("MC", "MC", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Component_-1/fit_results.txt");
    TChain* t1 = new TChain("DecayTree");
    t1->AddFileInfoList((TCollection*)fc1->GetList());
    t->AddFriend(t1);

    TH1D* h1_mb = new TH1D("h1_mb", "h1_mb", 100, 4000, 8000);
    TH1D* h2_mb = new TH1D("h2_mb", "h2_mb", 100, 4000, 8000);

    t->Draw("Bp_dtf_12_M >> h1_mb",passDTF);
    t->Draw("df_Bp_M >> h2_mb", df_status);

    gStyle->SetOptStat(0);

    // B+ mass
    TCanvas c;
    c.cd();
    h2_mb->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2_mb->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_mb->SetTitle("");
    h1_mb->SetLineColor(kBlack);
    h1_mb->SetFillColorAlpha(kBlack, 0.25);
    h2_mb->SetLineColor(kBlue);
    h2_mb->SetFillColorAlpha(kBlue, 0.25);
    h2_mb->DrawNormalized("same");
    h1_mb->DrawNormalized("same");

    TLegend* leg = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(h1_mb,"Nominal DTF","f");
    leg->AddEntry(h2_mb,"New fit","f");
    leg->Draw("same");

    int bin1_1 = h1_mb->FindFirstBinAbove(h1_mb->GetMaximum()/2);
    int bin2_1 = h1_mb->FindLastBinAbove(h1_mb->GetMaximum()/2);
    double fwhm_1 = h1_mb->GetBinCenter(bin2_1) - h1_mb->GetBinCenter(bin1_1);

    int bin1_2 = h2_mb->FindFirstBinAbove(h2_mb->GetMaximum()/2);
    int bin2_2 = h2_mb->FindLastBinAbove(h2_mb->GetMaximum()/2);
    double fwhm_2 = h2_mb->GetBinCenter(bin2_2) - h2_mb->GetBinCenter(bin1_2);

    double sig_1 = fwhm_1/2.355;
    double sig_2 = fwhm_2/2.355;

    int bin1_max = h1_mb->GetMaximumBin();
    int bin2_max = h2_mb->GetMaximumBin();
    double m1_peak = h1_mb->GetBinCenter(bin1_max);
    double m2_peak = h1_mb->GetBinCenter(bin2_max);

    cout << "Nominal DTF : B+ mass has resolution of " << sig_1 << " MeV" << endl;
    cout << "Standalone fit : B+ mass has resolution of " << sig_2 << " MeV" << endl;

    cout << "Nominal DTF: B+ mass peak is " << m1_peak << " MeV" << endl;
    cout << "Standalone fit: B+ mass peak is " << m2_peak << " MeV" << endl;

    c.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass.gif");
    c.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass.pdf");

    // B+ mass error
    TH1D* h1_mberr = new TH1D("h1_mberr", "h1_mberr", 100, 0, 5000);
    TH1D* h2_mberr = new TH1D("h2_mberr", "h2_mberr", 100, 0, 5000);
    t->Draw("Bp_dtf_12_MERR >> h1_mberr",passDTF);
    t->Draw("df_Bp_MERR >> h2_mberr", df_status);

    TCanvas c1;
    c1.cd();
    h2_mberr->GetXaxis()->SetTitle("#Delta m_{B} (MeV)");
    h2_mberr->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_mberr->SetTitle("");
    h1_mberr->SetLineColor(kBlack);
    h1_mberr->SetFillColorAlpha(kBlack, 0.25);
    h2_mberr->SetLineColor(kBlue);
    h2_mberr->SetFillColorAlpha(kBlue, 0.25);
    h2_mberr->DrawNormalized("same");
    h1_mberr->DrawNormalized("same");

    TLegend* leg1 = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h1_mberr,"Nominal DTF","f");
    leg1->AddEntry(h2_mberr,"New fit","f");
    leg1->Draw("same");
    c1.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_err.gif");
    c1.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_err.pdf");

    // B+ mass pull
    gStyle->SetOptStat(1);
    TH1D* h1_pull = new TH1D("h1_pull", "h1_pull", 100, -10, 10);
    TH1D* h2_pull = new TH1D("h2_pull", "h2_pull", 100, -10, 10);
    t->Draw("(Bp_dtf_12_M-5279)/Bp_dtf_12_MERR >> h1_pull",passDTF);
    t->Draw("(df_Bp_M-5279)/df_Bp_MERR >> h2_pull", df_status);

    TCanvas c2;
    c2.cd();
    h1_pull->GetXaxis()->SetTitle("Pull: (m_{B}-m_{B}^{PDG})/#Delta m_{B}");
    h1_pull->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_pull->SetTitle("");
    h1_pull->SetLineColor(kBlack);
    h1_pull->SetFillColorAlpha(kBlack, 0.25);
    h2_pull->SetLineColor(kBlue);
    h2_pull->SetFillColorAlpha(kBlue, 0.25);

    h1_pull->DrawNormalized();
    h2_pull->DrawNormalized("same");

    TLegend* leg2 = new TLegend(0.2, 0.6, 0.5, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(h1_pull,"Nominal DTF","f");
    leg2->AddEntry(h2_pull,"New fit","f");
    leg2->Draw("same");
    c2.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_pull.gif");
    c2.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_pull.pdf");

    // Bias in antineutrino PZ
    TH1D* h1_bias1 = new TH1D("h1_bias1", "h1_bias1", 100, -100000, 100000);
    TH1D* h2_bias1 = new TH1D("h2_bias1", "h2_bias1", 100, -100000, 100000);
    t->Draw("Bp_dtf_12_tauminus_nu_tau_PZ - (taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z) >> h1_bias1",passDTF);
    t->Draw("df_antinutau_PZ - (taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z) >> h2_bias1", df_status);

    TCanvas c3;
    c3.cd();
    h1_bias1->GetXaxis()->SetTitle("Fit - TRUE antineutrino PZ (MeV)");
    h1_bias1->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_bias1->SetTitle("");
    h1_bias1->SetLineColor(kBlack);
    h1_bias1->SetFillColorAlpha(kBlack, 0.25);
    h2_bias1->SetLineColor(kBlue);
    h2_bias1->SetFillColorAlpha(kBlue, 0.25);
    h1_bias1->DrawNormalized("same");
    h2_bias1->DrawNormalized("same");

    TLegend* leg3 = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.03);
    leg3->AddEntry(h1_bias1,"Nominal DTF","f");
    leg3->AddEntry(h2_bias1,"New fit","f");
    leg3->Draw("same");
    c3.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_antinutau_PZ_bias.gif");
    c3.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_antinutau_PZ_bias.pdf");

    // Bias in neutrino PZ
    TH1D* h1_bias2 = new TH1D("h1_bias2", "h1_bias2", 100, -100000, 100000);
    TH1D* h2_bias2 = new TH1D("h2_bias2", "h2_bias2", 100, -100000, 100000);
    t->Draw("Bp_dtf_12_tauminus_0_nu_tau_PZ - (taum_TRUEP_Z - taum_pi1_TRUEP_Z - taum_pi2_TRUEP_Z - taum_pi3_TRUEP_Z) >> h1_bias2",passDTF);
    t->Draw("df_nutau_PZ - (taum_TRUEP_Z - taum_pi1_TRUEP_Z - taum_pi2_TRUEP_Z - taum_pi3_TRUEP_Z) >> h2_bias2", df_status);

    TCanvas c4;
    c4.cd();
    h1_bias2->GetXaxis()->SetTitle("Fit - TRUE neutrino PZ (MeV)");
    h1_bias2->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_bias2->SetTitle("");
    h1_bias2->SetLineColor(kBlack);
    h1_bias2->SetFillColorAlpha(kBlack, 0.25);
    h2_bias2->SetLineColor(kBlue);
    h2_bias2->SetFillColorAlpha(kBlue, 0.25);
    h1_bias2->DrawNormalized("same");
    h2_bias2->DrawNormalized("same");

    TLegend* leg4 = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.03);
    leg4->AddEntry(h1_bias2,"Nominal DTF","f");
    leg4->AddEntry(h2_bias2,"New fit","f");
    leg4->Draw("same");
    c4.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_nutau_PZ_bias.gif");
    c4.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/Bmass_nutau_PZ_bias.pdf");

    // tau+ mass
    TH1D* h1_taup = new TH1D("h1_mtaup", "h1_mtaup", 100, 1600, 1900);
    TH1D* h2_taup = new TH1D("h2_mtaup", "h2_mtaup", 100, 1600, 1900);

    t->Draw("Bp_dtf_12_tauminus_M >> h1_mtaup",passDTF);
    t->Draw("df_taup_M  >> h2_mtaup", df_status);

    TCanvas cs5;
    cs5.cd();
    h2_taup->GetXaxis()->SetTitle("m_{#tau^{+}} (MeV)");
    h2_taup->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_taup->SetTitle("");
    h1_taup->SetLineColor(kBlack);
    h1_taup->SetFillColorAlpha(kBlack, 0.25);
    h2_taup->SetLineColor(kBlue);
    h2_taup->SetFillColorAlpha(kBlue, 0.25);
    h2_taup->DrawNormalized("same");
    h1_taup->DrawNormalized("same");

    TLegend* legs5 = new TLegend(0.7, 0.6, 0.85, 0.85);
    legs5->SetBorderSize(0);
    legs5->SetTextSize(0.03);
    legs5->AddEntry(h1_taup,"Nominal DTF","f");
    legs5->AddEntry(h2_taup,"New fit","f");
    legs5->Draw("same");
    cs5.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/taup.gif");
    cs5.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/taup.pdf");

    // Energy conservation in DV1
    TH1D* h1_E1 = new TH1D("h1_E1", "h1_E1", 100, -500, 500);
    TH1D* h2_E1 = new TH1D("h2_E1", "h2_E1", 100, -500, 500);

    t->Draw("Bp_dtf_12_tauminus_PE - Bp_dtf_12_tauminus_piplus_PE - Bp_dtf_12_tauminus_piplus_0_PE - Bp_dtf_12_tauminus_piplus_1_PE - Bp_dtf_12_tauminus_nu_tau_PE >> h1_E1",passDTF);
    t->Draw(" sqrt( pow(1776.86,2) + pow(df_taup_PX,2) + pow(df_taup_PY,2) + pow(df_taup_PZ,2) ) - taup_pi1_PE - taup_pi2_PE - taup_pi3_PE - sqrt( pow(df_antinutau_PX,2) + pow(df_antinutau_PY,2) + pow(df_antinutau_PZ,2) ) >> h2_E1", df_status);

    TCanvas c5;
    c5.cd();
    h2_E1->GetXaxis()->SetTitle("#Delta E in DV1 (MeV)");
    h2_E1->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_E1->SetTitle("");
    h1_E1->SetLineColor(kBlack);
    h1_E1->SetFillColorAlpha(kBlack, 0.25);
    h2_E1->SetLineColor(kBlue);
    h2_E1->SetFillColorAlpha(kBlue, 0.25);
    h2_E1->DrawNormalized("same");
    h1_E1->DrawNormalized("same");

    TLegend* leg5 = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg5->SetBorderSize(0);
    leg5->SetTextSize(0.03);
    leg5->AddEntry(h1_E1,"Nominal DTF","f");
    leg5->AddEntry(h2_E1,"New fit","f");
    leg5->Draw("same");
    c5.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/E_conservation_DV1.gif");
    c5.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/E_conservation_DV1.pdf");

    // Energy conservation in DV2
    TH1D* h1_E2 = new TH1D("h1_E2", "h1_E2", 100, -500, 500);
    TH1D* h2_E2 = new TH1D("h2_E2", "h2_E2", 100, -500, 500);

    t->Draw("Bp_dtf_12_tauminus_0_PE - Bp_dtf_12_tauminus_0_piplus_PE - Bp_dtf_12_tauminus_0_piplus_0_PE - Bp_dtf_12_tauminus_0_piplus_1_PE - Bp_dtf_12_tauminus_0_nu_tau_PE >> h1_E2",passDTF);
    t->Draw(" sqrt( pow(1776.86,2) + pow(df_taum_PX,2) + pow(df_taum_PY,2) + pow(df_taum_PZ,2) ) - taum_pi1_PE - taum_pi2_PE - taum_pi3_PE - sqrt( pow(df_nutau_PX,2) + pow(df_nutau_PY,2) + pow(df_nutau_PZ,2) ) >> h2_E2", df_status);

    TCanvas c6;
    c6.cd();
    h2_E2->GetXaxis()->SetTitle("#Delta E in DV2 (MeV)");
    h2_E2->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_E2->SetTitle("");
    h1_E2->SetLineColor(kBlack);
    h1_E2->SetFillColorAlpha(kBlack, 0.25);
    h2_E2->SetLineColor(kBlue);
    h2_E2->SetFillColorAlpha(kBlue, 0.25);
    h2_E2->DrawNormalized("same");
    h1_E2->DrawNormalized("same");

    TLegend* leg6 = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg6->SetBorderSize(0);
    leg6->SetTextSize(0.03);
    leg6->AddEntry(h1_E2,"Nominal DTF","f");
    leg6->AddEntry(h2_E2,"New fit","f");
    leg6->Draw("same");
    c6.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/E_conservation_DV2.gif");
    c6.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/E_conservation_DV2.pdf");

    // Energy conservation in BV
    TH1D* h1_E = new TH1D("h1_E", "h1_E", 100, -500, 500);
    TH1D* h2_E = new TH1D("h2_E", "h2_E", 100, -500, 500);

    t->Draw("Bp_dtf_12_PE - Bp_dtf_12_tauminus_PE - Bp_dtf_12_tauminus_0_PE - Kp_PE >> h1_E",passDTF);
    t->Draw(" sqrt( pow(df_Bp_M,2) + pow(df_Bp_PX,2) + pow(df_Bp_PY,2) + pow(df_Bp_PZ,2) ) - sqrt( pow(1776.86,2) + pow(df_taum_PX,2) + pow(df_taum_PY,2) + pow(df_taum_PZ,2) ) - sqrt( pow(1776.86,2) + pow(df_taup_PX,2) + pow(df_taup_PY,2) + pow(df_taup_PZ,2) ) - Kp_PE  >> h2_E", df_status);

    TCanvas c7;
    c7.cd();
    h1_E->GetXaxis()->SetTitle("#Delta E in BV (MeV)");
    h1_E->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_E->SetTitle("");
    h1_E->SetLineColor(kBlack);
    h1_E->SetFillColorAlpha(kBlack, 0.25);
    h2_E->SetLineColor(kBlue);
    h2_E->SetFillColorAlpha(kBlue, 0.25);
    h1_E->DrawNormalized("same");
    h2_E->DrawNormalized("same");

    TLegend* leg7 = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg7->SetBorderSize(0);
    leg7->SetTextSize(0.03);
    leg7->AddEntry(h1_E,"Nominal DTF","f");
    leg7->AddEntry(h2_E,"New fit","f");
    leg7->Draw("same");
    c7.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/E_conservation_BV.gif");
    c7.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/E_conservation_BV.pdf");

    // 3-momentum conservation in DV1
    TH1D* h1_P1x = new TH1D("h1_P1x", "h1_P1x", 100, -10, 10);
    TH1D* h2_P1x = new TH1D("h2_P1x", "h2_P1x", 100, -10, 10);
    TH1D* h1_P1y = new TH1D("h1_P1y", "h1_P1y", 100, -10, 10);
    TH1D* h2_P1y = new TH1D("h2_P1y", "h2_P1y", 100, -10, 10);
    TH1D* h1_P1z = new TH1D("h1_P1z", "h1_P1z", 100, -10, 10);
    TH1D* h2_P1z = new TH1D("h2_P1z", "h2_P1z", 100, -10, 10);

    t->Draw("Bp_dtf_12_tauminus_PX - Bp_dtf_12_tauminus_piplus_PX - Bp_dtf_12_tauminus_piplus_0_PX - Bp_dtf_12_tauminus_piplus_1_PX - Bp_dtf_12_tauminus_nu_tau_PX >> h1_P1x",passDTF);
    t->Draw("Bp_dtf_12_tauminus_PY - Bp_dtf_12_tauminus_piplus_PY - Bp_dtf_12_tauminus_piplus_0_PY - Bp_dtf_12_tauminus_piplus_1_PY - Bp_dtf_12_tauminus_nu_tau_PY >> h1_P1y",passDTF);
    t->Draw("Bp_dtf_12_tauminus_PZ - Bp_dtf_12_tauminus_piplus_PZ - Bp_dtf_12_tauminus_piplus_0_PZ - Bp_dtf_12_tauminus_piplus_1_PZ - Bp_dtf_12_tauminus_nu_tau_PZ >> h1_P1z",passDTF);

    t->Draw("df_taup_PX - df_m_7 - df_m_10 - df_m_13 - df_antinutau_PX  >> h2_P1x", df_status);
    t->Draw("df_taup_PY - df_m_8 - df_m_11 - df_m_14 - df_antinutau_PY  >> h2_P1y", df_status);
    t->Draw("df_taup_PZ - df_m_9 - df_m_12 - df_m_15 - df_antinutau_PZ  >> h2_P1z", df_status);

    TCanvas c8x;
    c8x.cd();
    h1_P1x->GetXaxis()->SetTitle("#Delta PX in DV1 (MeV)");
    h1_P1x->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_P1x->SetTitle("");
    h1_P1x->SetLineColor(kBlack);
    h1_P1x->SetFillColorAlpha(kBlack, 0.25);
    h2_P1x->SetLineColor(kBlue);
    h2_P1x->SetFillColorAlpha(kBlue, 0.25);
    h1_P1x->DrawNormalized("same");
    h2_P1x->DrawNormalized("same");

    TLegend* leg8x = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg8x->SetBorderSize(0);
    leg8x->SetTextSize(0.03);
    leg8x->AddEntry(h1_P1x,"Nominal DTF","f");
    leg8x->AddEntry(h2_P1x,"New fit","f");
    leg8x->Draw("same");
    c8x.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PX_conservation_DV1.gif");
    c8x.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PX_conservation_DV1.pdf");

    TCanvas c8y;
    c8y.cd();
    h1_P1y->GetXaxis()->SetTitle("#Delta PY in DV1 (MeV)");
    h1_P1y->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_P1y->SetTitle("");
    h1_P1y->SetLineColor(kBlack);
    h1_P1y->SetFillColorAlpha(kBlack, 0.25);
    h2_P1y->SetLineColor(kBlue);
    h2_P1y->SetFillColorAlpha(kBlue, 0.25);
    h1_P1y->DrawNormalized("same");
    h2_P1y->DrawNormalized("same");

    TLegend* leg8y = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg8y->SetBorderSize(0);
    leg8y->SetTextSize(0.03);
    leg8y->AddEntry(h1_P1y,"Nominal DTF","f");
    leg8y->AddEntry(h2_P1y,"New fit","f");
    leg8y->Draw("same");
    c8y.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PY_conservation_DV1.gif");
    c8y.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PY_conservation_DV1.pdf");

    TCanvas c8z;
    c8z.cd();
    h1_P1z->GetXaxis()->SetTitle("#Delta PZ in DV1 (MeV)");
    h1_P1z->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_P1z->SetTitle("");
    h1_P1z->SetLineColor(kBlack);
    h1_P1z->SetFillColorAlpha(kBlack, 0.25);
    h2_P1z->SetLineColor(kBlue);
    h2_P1z->SetFillColorAlpha(kBlue, 0.25);
    h1_P1z->DrawNormalized("same");
    h2_P1z->DrawNormalized("same");

    TLegend* leg8z = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg8z->SetBorderSize(0);
    leg8z->SetTextSize(0.03);
    leg8z->AddEntry(h1_P1z,"Nominal DTF","f");
    leg8z->AddEntry(h2_P1z,"New fit","f");
    leg8z->Draw("same");
    c8z.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PZ_conservation_DV1.gif");
    c8z.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PZ_conservation_DV1.pdf");

    // 3-momentum conservation in DV2
    TH1D* h1_P2x = new TH1D("h1_P2x", "h1_P2x", 100, -10, 10);
    TH1D* h2_P2x = new TH1D("h2_P2x", "h2_P2x", 100, -10, 10);
    TH1D* h1_P2y = new TH1D("h1_P2y", "h1_P2y", 100, -10, 10);
    TH1D* h2_P2y = new TH1D("h2_P2y", "h2_P2y", 100, -10, 10);
    TH1D* h1_P2z = new TH1D("h1_P2z", "h1_P2z", 100, -10, 10);
    TH1D* h2_P2z = new TH1D("h2_P2z", "h2_P2z", 100, -10, 10);

    t->Draw("Bp_dtf_12_tauminus_0_PX - Bp_dtf_12_tauminus_0_piplus_PX - Bp_dtf_12_tauminus_0_piplus_0_PX - Bp_dtf_12_tauminus_0_piplus_1_PX - Bp_dtf_12_tauminus_0_nu_tau_PX >> h1_P2x",passDTF);
    t->Draw("Bp_dtf_12_tauminus_0_PY - Bp_dtf_12_tauminus_0_piplus_PY - Bp_dtf_12_tauminus_0_piplus_0_PY - Bp_dtf_12_tauminus_0_piplus_1_PY - Bp_dtf_12_tauminus_0_nu_tau_PY >> h1_P2y",passDTF);
    t->Draw("Bp_dtf_12_tauminus_0_PZ - Bp_dtf_12_tauminus_0_piplus_PZ - Bp_dtf_12_tauminus_0_piplus_0_PZ - Bp_dtf_12_tauminus_0_piplus_1_PZ - Bp_dtf_12_tauminus_0_nu_tau_PZ >> h1_P2z",passDTF);

    t->Draw("df_taum_PX - df_m_19 - df_m_22 - df_m_25 - df_nutau_PX  >> h2_P2x", df_status);
    t->Draw("df_taum_PY - df_m_20 - df_m_23 - df_m_26 - df_nutau_PY  >> h2_P2y", df_status);
    t->Draw("df_taum_PZ - df_m_21 - df_m_24 - df_m_27 - df_nutau_PZ  >> h2_P2z", df_status);

    TCanvas c9x;
    c9x.cd();
    h1_P2x->GetXaxis()->SetTitle("#Delta PX in DV2 (MeV)");
    h1_P2x->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_P2x->SetTitle("");
    h1_P2x->SetLineColor(kBlack);
    h1_P2x->SetFillColorAlpha(kBlack, 0.25);
    h2_P2x->SetLineColor(kBlue);
    h2_P2x->SetFillColorAlpha(kBlue, 0.25);
    h1_P2x->DrawNormalized("same");
    h2_P2x->DrawNormalized("same");

    TLegend* leg9x = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg9x->SetBorderSize(0);
    leg9x->SetTextSize(0.03);
    leg9x->AddEntry(h1_P2x,"Nominal DTF","f");
    leg9x->AddEntry(h2_P2x,"New fit","f");
    leg9x->Draw("same");
    c9x.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PX_conservation_DV2.gif");
    c9x.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PX_conservation_DV2.pdf");

    TCanvas c9y;
    c9y.cd();
    h1_P2y->GetXaxis()->SetTitle("#Delta PY in DV2 (MeV)");
    h1_P2y->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_P2y->SetTitle("");
    h1_P2y->SetLineColor(kBlack);
    h1_P2y->SetFillColorAlpha(kBlack, 0.25);
    h2_P2y->SetLineColor(kBlue);
    h2_P2y->SetFillColorAlpha(kBlue, 0.25);
    h1_P2y->DrawNormalized("same");
    h2_P2y->DrawNormalized("same");

    TLegend* leg9y = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg9y->SetBorderSize(0);
    leg9y->SetTextSize(0.03);
    leg9y->AddEntry(h1_P2y,"Nominal DTF","f");
    leg9y->AddEntry(h2_P2y,"New fit","f");
    leg9y->Draw("same");
    c9y.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PY_conservation_DV2.gif");
    c9y.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PY_conservation_DV2.pdf");

    TCanvas c9z;
    c9z.cd();
    h1_P2z->GetXaxis()->SetTitle("#Delta PZ in DV2 (MeV)");
    h1_P2z->GetYaxis()->SetTitle("Events / (100) MeV");
    h1_P2z->SetTitle("");
    h1_P2z->SetLineColor(kBlack);
    h1_P2z->SetFillColorAlpha(kBlack, 0.25);
    h2_P2z->SetLineColor(kBlue);
    h2_P2z->SetFillColorAlpha(kBlue, 0.25);
    h1_P2z->DrawNormalized("same");
    h2_P2z->DrawNormalized("same");

    TLegend* leg9z = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg9z->SetBorderSize(0);
    leg9z->SetTextSize(0.03);
    leg9z->AddEntry(h1_P2z,"Nominal DTF","f");
    leg9z->AddEntry(h2_P2z,"New fit","f");
    leg9z->Draw("same");
    c9z.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PZ_conservation_DV2.gif");
    c9z.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PZ_conservation_DV2.pdf");

    // 3-momentum conservation in BV
    TH1D* h1_Px = new TH1D("h1_Px", "h1_Px", 100, -10, 10);
    TH1D* h2_Px = new TH1D("h2_Px", "h2_Px", 100, -10, 10);
    TH1D* h1_Py = new TH1D("h1_Py", "h1_Py", 100, -10, 10);
    TH1D* h2_Py = new TH1D("h2_Py", "h2_Py", 100, -10, 10);
    TH1D* h1_Pz = new TH1D("h1_Pz", "h1_Pz", 100, -10, 10);
    TH1D* h2_Pz = new TH1D("h2_Pz", "h2_Pz", 100, -10, 10);

    t->Draw("Bp_dtf_12_PX - Bp_dtf_12_tauminus_PX - Bp_dtf_12_tauminus_0_PX - Kp_PX >> h1_Px",passDTF);
    t->Draw("Bp_dtf_12_PY - Bp_dtf_12_tauminus_PY - Bp_dtf_12_tauminus_0_PY - Kp_PY >> h1_Py",passDTF);
    t->Draw("Bp_dtf_12_PZ - Bp_dtf_12_tauminus_PZ - Bp_dtf_12_tauminus_0_PZ - Kp_PZ >> h1_Pz",passDTF);

    t->Draw("df_Bp_PX - Kp_PX - df_taup_PX - df_taum_PX  >> h2_Px", df_status);
    t->Draw("df_Bp_PY - Kp_PY - df_taup_PY - df_taum_PY  >> h2_Py", df_status);
    t->Draw("df_Bp_PZ - Kp_PZ - df_taup_PZ - df_taum_PZ  >> h2_Pz", df_status);

    TCanvas c10x;
    c10x.cd();
    h2_Px->GetXaxis()->SetTitle("#Delta PX in BV (MeV)");
    h2_Px->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_Px->SetTitle("");
    h1_Px->SetLineColor(kBlack);
    h1_Px->SetFillColorAlpha(kBlack, 0.25);
    h2_Px->SetLineColor(kBlue);
    h2_Px->SetFillColorAlpha(kBlue, 0.25);
    h2_Px->DrawNormalized("same");
    h1_Px->DrawNormalized("same");

    TLegend* leg10x = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg10x->SetBorderSize(0);
    leg10x->SetTextSize(0.03);
    leg10x->AddEntry(h1_Px,"Nominal DTF","f");
    leg10x->AddEntry(h2_Px,"New fit","f");
    leg10x->Draw("same");
    c10x.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PX_conservation_BV.gif");
    c10x.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PX_conservation_BV.pdf");

    TCanvas c10y;
    c10y.cd();
    h2_Py->GetXaxis()->SetTitle("#Delta PY in BV (MeV)");
    h2_Py->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_Py->SetTitle("");
    h1_Py->SetLineColor(kBlack);
    h1_Py->SetFillColorAlpha(kBlack, 0.25);
    h2_Py->SetLineColor(kBlue);
    h2_Py->SetFillColorAlpha(kBlue, 0.25);
    h2_Py->DrawNormalized("same");
    h1_Py->DrawNormalized("same");

    TLegend* leg10y = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg10y->SetBorderSize(0);
    leg10y->SetTextSize(0.03);
    leg10y->AddEntry(h1_Py,"Nominal DTF","f");
    leg10y->AddEntry(h2_Py,"New fit","f");
    leg10y->Draw("same");
    c10y.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PY_conservation_BV.gif");
    c10y.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PY_conservation_BV.pdf");

    TCanvas c10z;
    c10z.cd();
    h2_Pz->GetXaxis()->SetTitle("#Delta PZ in BV (MeV)");
    h2_Pz->GetYaxis()->SetTitle("Events / (100) MeV");
    h2_Pz->SetTitle("");
    h1_Pz->SetLineColor(kBlack);
    h1_Pz->SetFillColorAlpha(kBlack, 0.25);
    h2_Pz->SetLineColor(kBlue);
    h2_Pz->SetFillColorAlpha(kBlue, 0.25);
    h2_Pz->DrawNormalized("same");
    h1_Pz->DrawNormalized("same");

    TLegend* leg10z = new TLegend(0.7, 0.6, 0.85, 0.85);
    leg10z->SetBorderSize(0);
    leg10z->SetTextSize(0.03);
    leg10z->AddEntry(h1_Pz,"Nominal DTF","f");
    leg10z->AddEntry(h2_Pz,"New fit","f");
    leg10z->Draw("same");
    c10z.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PZ_conservation_BV.gif");
    c10z.SaveAs("/panfs/felician/B2Ktautau/ROOT_Sim/2018/Plots/PZ_conservation_BV.pdf");

}