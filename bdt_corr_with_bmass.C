

void bdt_corr_with_bmass()
{
    gStyle->SetOptStat(0);

    Double_t bdt_cuts[] = {0., 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95};
    Int_t colors[] = {kBlack+0, kBlue+0, kMagenta+0, kRed+0, kGreen+1, kCyan+1, kViolet+1, kOrange-3, kRed+3};
    Int_t N = sizeof(bdt_cuts)/sizeof(bdt_cuts[0]);

    ///////////////////////////////////////////////////////////////////////////////// 3pi3pi MC ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    std::vector<TH1D*> histos_isolation_3pi3pi;
    std::vector<TH1D*> histos_kinematic_3pi3pi;
    TLegend *leg_isolation_3pi3pi = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *leg_kinematic_3pi3pi = new TLegend(0.7, 0.7, 0.9, 0.9);

    for(int i = 0; i < N; i++)
    {
        histos_isolation_3pi3pi.push_back( new TH1D( Form("h_iso_3pi3pi_%i",i), Form("h_iso_3pi3pi_%i",i), 100, 4000, 8000) );
        histos_kinematic_3pi3pi.push_back( new TH1D( Form("h_kin_3pi3pi_%i",i), Form("h_kin_3pi3pi_%i",i), 100, 4000, 8000) );

        t_3pi3pi_2016->Draw(Form("df_Bp_M >> h_iso_3pi3pi_%i",i), Form("(df_status==0) && (BDT1 > %f)",bdt_cuts[i]));
        t_3pi3pi_2016->Draw(Form("df_Bp_M >> h_kin_3pi3pi_%i",i), Form("(df_status==0) && (BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_3pi3pi[i]->SetLineColor(colors[i]);
        histos_kinematic_3pi3pi[i]->SetLineColor(colors[i]);

        if(i == N-1)
        {
            histos_isolation_3pi3pi[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_3pi3pi[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_3pi3pi[i]->SetTitle("Isolation MVA (3#pi3#pi MC) 2016-2018");

            histos_kinematic_3pi3pi[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_3pi3pi[i]->GetYaxis()->SetTitle("BDT2");
            histos_kinematic_3pi3pi[i]->SetTitle("Topological MVA (3#pi3#pi MC) 2016-2018");
        }
    }

    TCanvas c_3pi3pi_isolation;
    c_3pi3pi_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_3pi3pi[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_isolation_3pi3pi[N-1-i]->DrawNormalized("same");
        }
        leg_isolation_3pi3pi->AddEntry(Form("h_iso_3pi3pi_%i",N-1-i), Form("BDT1 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_isolation_3pi3pi->Draw("same");
    c_3pi3pi_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_3pi3pi_MC.pdf");

    TCanvas c_3pi3pi_kinematic;
    c_3pi3pi_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_3pi3pi[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_3pi3pi[N-1-i]->DrawNormalized("same");
        }
        leg_kinematic_3pi3pi->AddEntry(Form("h_kin_3pi3pi_%i",N-1-i), Form("BDT2 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_kinematic_3pi3pi->Draw("same");
    c_3pi3pi_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_3pi3pi_MC.pdf");

    TH2D* h2D_isolation_3pi3pi = new TH2D("h2d_iso_3pi3pi", "h2d_iso_3pi3pi", 100, 4000, 8000, 100, 0, 1);
    TH2D* h2D_kinematic_3pi3pi = new TH2D("h2d_kin_3pi3pi", "h2d_kin_3pi3pi", 100, 4000, 8000, 100, 0, 1);

    t_3pi3pi_2016->Draw("BDT1 : df_Bp_M >> h2d_iso_3pi3pi", "df_status==0");
    t_3pi3pi_2016->Draw("BDT2 : df_Bp_M >> h2d_kin_3pi3pi", "df_status==0");

    TCanvas c1_3pi3pi_isolation;
    c1_3pi3pi_isolation.cd();
    c1_3pi3pi_isolation.SetLogz();
    h2D_isolation_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_3pi3pi->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_3pi3pi->SetTitle(Form("Isolation MVA (3#pi3#pi MC) 2016-2018 (corr = %.3f)",h2D_isolation_3pi3pi->GetCorrelationFactor()));
    h2D_isolation_3pi3pi->Draw("COLZ");
    c1_3pi3pi_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_3pi3pi_MC.pdf");

    TCanvas c1_3pi3pi_kinematic;
    c1_3pi3pi_kinematic.cd();
    c1_3pi3pi_kinematic.SetLogz();
    h2D_kinematic_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_3pi3pi->GetYaxis()->SetTitle("BDT2");
    h2D_kinematic_3pi3pi->SetTitle(Form("Topological MVA (3#pi3#pi MC) 2016-2018 (corr = %.3lf)",h2D_kinematic_3pi3pi->GetCorrelationFactor()));
    h2D_kinematic_3pi3pi->Draw("COLZ");
    c1_3pi3pi_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT2_vs_Bmass_3pi3pi_MC.pdf");

    TCanvas c2_3pi3pi_isolation;
    c2_3pi3pi_isolation.cd();
    TH1D* h_profile_isolation_3pi3pi = (TH1D*)h2D_isolation_3pi3pi->ProfileX();
    h_profile_isolation_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_isolation_3pi3pi->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_isolation_3pi3pi->SetTitle(Form("Isolation MVA (3#pi3#pi MC) 2016-2018 (corr = %.3f)",h2D_isolation_3pi3pi->GetCorrelationFactor()));
    h_profile_isolation_3pi3pi->Draw();
    c2_3pi3pi_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_profile_3pi3pi_MC.pdf");

    TCanvas c2_3pi3pi_kinematic;
    c2_3pi3pi_kinematic.cd();
    TH1D* h_profile_kinematic_3pi3pi = (TH1D*)h2D_kinematic_3pi3pi->ProfileX();
    h_profile_kinematic_3pi3pi->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_kinematic_3pi3pi->GetYaxis()->SetTitle("Average BDT2 per m_{B} bin");
    h_profile_kinematic_3pi3pi->SetTitle(Form("Topological MVA (3#pi3#pi MC) 2016-2018 (corr = %.3f)",h2D_kinematic_3pi3pi->GetCorrelationFactor()));
    h_profile_kinematic_3pi3pi->Draw();
    c2_3pi3pi_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_profile_3pi3pi_MC.pdf");

    Double_t isolation_peaks_3pi3pi[N], kinematic_peaks_3pi3pi[N];
    Double_t isolation_peaks_error_3pi3pi[N], kinematic_peaks_error_3pi3pi[N];

    for(int i = 0; i < N; i++)
    {
        isolation_peaks_3pi3pi[i] = histos_isolation_3pi3pi[i]->GetBinCenter(histos_isolation_3pi3pi[i]->GetMaximumBin());
        kinematic_peaks_3pi3pi[i] = histos_kinematic_3pi3pi[i]->GetBinCenter(histos_kinematic_3pi3pi[i]->GetMaximumBin());

        isolation_peaks_error_3pi3pi[i] = 20; // half of the bin width
        kinematic_peaks_error_3pi3pi[i] = 20;

    }

    TGraph* g_isolation_3pi3pi = new TGraphErrors(N, bdt_cuts, isolation_peaks_3pi3pi, 0, isolation_peaks_error_3pi3pi);
    TGraph* g_kinematic_3pi3pi = new TGraphErrors(N, bdt_cuts, kinematic_peaks_3pi3pi, 0, kinematic_peaks_error_3pi3pi);

    TCanvas c3_3pi3pi_isolation;
    c3_3pi3pi_isolation.cd();
    c3_3pi3pi_isolation.SetGrid();
    c3_3pi3pi_isolation.GetPad(0)->SetLeftMargin(0.15);
    g_isolation_3pi3pi->SetTitle("Isolation MVA (3#pi3#pi MC) 2016-2018 ; BDT1 > x; X-position of B^{+} mass peak (MeV)");
    g_isolation_3pi3pi->SetMarkerStyle(8);
    g_isolation_3pi3pi->GetXaxis()->SetLimits(-0.1,1.1);
    g_isolation_3pi3pi->Draw("AP");
    c3_3pi3pi_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt1_3pi3pi_MC.pdf");

    TCanvas c3_3pi3pi_kinematic;
    c3_3pi3pi_kinematic.cd();
    c3_3pi3pi_kinematic.SetGrid();
    c3_3pi3pi_kinematic.GetPad(0)->SetLeftMargin(0.15);
    g_kinematic_3pi3pi->SetTitle("Topological MVA (3#pi3#pi MC) 2016-2018 ; BDT2 > x; X-position of B^{+} mass peak (MeV)");
    g_kinematic_3pi3pi->SetMarkerStyle(8);
    g_kinematic_3pi3pi->GetXaxis()->SetLimits(-0.1,1.1);
    g_kinematic_3pi3pi->Draw("AP");
    c3_3pi3pi_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt2_3pi3pi_MC.pdf");

    ///////////// 3pi3pi pi0 MC
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

    ////// 3pi3pi 2pi0 MC
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

    TH1D* h_3pi3pi_iso = new TH1D("h_3pi3pi_iso", "h_3pi3pi_iso", 100, 0, 1);
    TH1D* h_3pi3pipi0_iso = new TH1D("h_3pi3pipi0_iso", "h_3pi3pipi0_iso", 100, 0, 1);
    TH1D* h_3pi3pi2pi0_iso = new TH1D("h_3pi3pi2pi0_iso", "h_3pi3pi2pi0_iso", 100, 0, 1);

    t_3pi3pi_2016->Draw("BDT1 >> h_3pi3pi_iso", "df_status==0");
    t_3pi3pipi0_2016->Draw("BDT1 >> h_3pi3pipi0_iso", "df_status==0");
    t_3pi3pi2pi0_2016->Draw("BDT1 >> h_3pi3pi2pi0_iso", "df_status==0");

    TCanvas c4;
    c4.cd();
    h_3pi3pi_iso->GetXaxis()->SetTitle("BDT1");
    h_3pi3pi_iso->GetYaxis()->SetTitle("Normalised entries / (0.025)");
    h_3pi3pi_iso->SetTitle("Isolation MVA 2016-2018");
    h_3pi3pi_iso->SetLineColor(kBlue);
    h_3pi3pipi0_iso->SetLineColor(kOrange-3);
    h_3pi3pi2pi0_iso->SetLineColor(kGreen+1);
    h_3pi3pi_iso->DrawNormalized();
    h_3pi3pipi0_iso->DrawNormalized("same");
    h_3pi3pi2pi0_iso->DrawNormalized("same");
    TLegend *leg_comp_iso = new TLegend(0.1, 0.7, 0.4, 0.9);
    leg_comp_iso->AddEntry(h_3pi3pi_iso, "3#pi3#pi MC", "lp");
    leg_comp_iso->AddEntry(h_3pi3pipi0_iso, "3#pi3#pi #pi^{0} MC", "lp");
    leg_comp_iso->AddEntry(h_3pi3pi2pi0_iso, "3#pi3#pi 2#pi^{0} MC", "lp");
    leg_comp_iso->Draw("same");
    c4.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_mc_components.pdf");

    TH1D* h_3pi3pi_kin = new TH1D("h_3pi3pi_kin", "h_3pi3pi_kin", 100, 0, 1);
    TH1D* h_3pi3pipi0_kin = new TH1D("h_3pi3pipi0_kin", "h_3pi3pipi0_kin", 100, 0, 1);
    TH1D* h_3pi3pi2pi0_kin = new TH1D("h_3pi3pi2pi0_kin", "h_3pi3pi2pi0_kin", 100, 0, 1);

    t_3pi3pi_2016->Draw("BDT2 >> h_3pi3pi_kin", "df_status==0");
    t_3pi3pipi0_2016->Draw("BDT2 >> h_3pi3pipi0_kin", "df_status==0");
    t_3pi3pi2pi0_2016->Draw("BDT2 >> h_3pi3pi2pi0_kin", "df_status==0");

    TCanvas c5;
    c5.cd();
    h_3pi3pi_kin->SetLineColor(kBlue);
    h_3pi3pipi0_kin->SetLineColor(kOrange-3);
    h_3pi3pi2pi0_kin->SetLineColor(kGreen+1);
    h_3pi3pi_kin->GetXaxis()->SetTitle("BDT2");
    h_3pi3pi_kin->GetYaxis()->SetTitle("Normalised entries / (0.025)");
    h_3pi3pi_kin->SetTitle("Topological MVA 2016-2018");
    h_3pi3pi_kin->DrawNormalized();
    h_3pi3pipi0_kin->DrawNormalized("same");
    h_3pi3pi2pi0_kin->DrawNormalized("same");
    TLegend *leg_comp_kin = new TLegend(0.1, 0.7, 0.4, 0.9);
    leg_comp_kin->AddEntry(h_3pi3pi_kin, "3#pi3#pi MC", "lp");
    leg_comp_kin->AddEntry(h_3pi3pipi0_kin, "3#pi3#pi #pi^{0} MC", "lp");
    leg_comp_kin->AddEntry(h_3pi3pi2pi0_kin, "3#pi3#pi 2#pi^{0} MC", "lp");
    leg_comp_kin->Draw("same");
    c5.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_mc_components.pdf");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////////////////// RS data /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFileCollection *fc_rs_data_2016 = new TFileCollection("fc_rs_data_2016", "fc_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt");
    TFileCollection *fc_rs_data_2017 = new TFileCollection("fc_rs_data_2017", "fc_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt");
    TFileCollection *fc_rs_data_2018 = new TFileCollection("fc_rs_data_2018", "fc_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt");

    TChain *t_rs_data_2016 = new TChain("DecayTree");
    TChain *t_rs_data_2017 = new TChain("DecayTree");
    TChain *t_rs_data_2018 = new TChain("DecayTree");

    t_rs_data_2016->AddFileInfoList((TCollection*)fc_rs_data_2016->GetList());
    t_rs_data_2017->AddFileInfoList((TCollection*)fc_rs_data_2017->GetList());
    t_rs_data_2018->AddFileInfoList((TCollection*)fc_rs_data_2018->GetList());

    t_rs_data_2016->GetEntries();
    t_rs_data_2017->GetEntries();
    t_rs_data_2018->GetEntries();

    t_rs_data_2016->Add(t_rs_data_2017);
    t_rs_data_2016->Add(t_rs_data_2018);

    TFileCollection *fc1_rs_data_2016 = new TFileCollection("fc1_rs_data_2016", "fc1_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_2/bdt_output.txt");
    TFileCollection *fc1_rs_data_2017 = new TFileCollection("fc1_rs_data_2017", "fc1_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_2/bdt_output.txt");
    TFileCollection *fc1_rs_data_2018 = new TFileCollection("fc1_rs_data_2018", "fc1_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_2/bdt_output.txt");

    TChain *t1_rs_data_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_rs_data_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_rs_data_2018 = new TChain("XGBoost/DecayTree");

    t1_rs_data_2016->AddFileInfoList((TCollection*)fc1_rs_data_2016->GetList());
    t1_rs_data_2017->AddFileInfoList((TCollection*)fc1_rs_data_2017->GetList());
    t1_rs_data_2018->AddFileInfoList((TCollection*)fc1_rs_data_2018->GetList());

    t1_rs_data_2016->GetEntries();
    t1_rs_data_2017->GetEntries();
    t1_rs_data_2018->GetEntries();

    t1_rs_data_2016->Add(t1_rs_data_2017);
    t1_rs_data_2016->Add(t1_rs_data_2018);

    t_rs_data_2016->AddFriend(t1_rs_data_2016);

    std::vector<TH1D*> histos_isolation_rs_data;
    std::vector<TH1D*> histos_kinematic_rs_data;
    TLegend *leg_isolation_rs_data = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *leg_kinematic_rs_data = new TLegend(0.7, 0.7, 0.9, 0.9);

    for(int i = 0; i < N; i++)
    {
        histos_isolation_rs_data.push_back( new TH1D( Form("h_iso_rs_data_%i",i), Form("h_iso_rs_data_%i",i), 100, 4000, 8000) );
        histos_kinematic_rs_data.push_back( new TH1D( Form("h_kin_rs_data_%i",i), Form("h_kin_rs_data_%i",i), 100, 4000, 8000) );

        t_rs_data_2016->Draw(Form("df_Bp_M >> h_iso_rs_data_%i",i), Form("(df_status==0) && ((df_Bp_M < 4700) || (df_Bp_M > 5800)) && (BDT1 > %f)",bdt_cuts[i]));
        t_rs_data_2016->Draw(Form("df_Bp_M >> h_kin_rs_data_%i",i), Form("(df_status==0) && ((df_Bp_M < 4700) || (df_Bp_M > 5800)) && (BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_rs_data[i]->SetLineColor(colors[i]);
        histos_kinematic_rs_data[i]->SetLineColor(colors[i]);

        if(i == 0)
        {
            histos_isolation_rs_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_rs_data[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_rs_data[i]->SetTitle("Isolation MVA (RS data) 2016-2018");

            histos_kinematic_rs_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_rs_data[i]->GetYaxis()->SetTitle("BDT2");
            histos_kinematic_rs_data[i]->SetTitle("Topological MVA (RS data) 2016-2018");
        }
    }

    TCanvas c_rs_data_isolation;
    c_rs_data_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_rs_data[i]->DrawNormalized();
        }
        else
        {
            histos_isolation_rs_data[i]->DrawNormalized("same");
        }
        leg_isolation_rs_data->AddEntry(Form("h_iso_rs_data_%i",i), Form("BDT1 > %.2f",bdt_cuts[i]), "lp");
    }
    leg_isolation_rs_data->Draw("same");
    c_rs_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_rs_data_MC.pdf");

    TCanvas c_rs_data_kinematic;
    c_rs_data_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_rs_data[i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_rs_data[i]->DrawNormalized("same");
        }
        leg_kinematic_rs_data->AddEntry(Form("h_kin_rs_data_%i",i), Form("BDT2 > %.2f",bdt_cuts[i]), "lp");
    }
    leg_kinematic_rs_data->Draw("same");
    c_rs_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_rs_data_MC.pdf");

    TH2D* h2D_isolation_rs_data = new TH2D("h2d_iso_rs_data", "h2d_iso_rs_data", 100, 4000, 8000, 100, 0, 1);
    TH2D* h2D_kinematic_rs_data = new TH2D("h2d_kin_rs_data", "h2d_kin_rs_data", 100, 4000, 8000, 100, 0, 1);

    t_rs_data_2016->Draw("BDT1 : df_Bp_M >> h2d_iso_rs_data", "df_status==0");
    t_rs_data_2016->Draw("BDT2 : df_Bp_M >> h2d_kin_rs_data", "df_status==0");

    TCanvas c1_rs_data_isolation;
    c1_rs_data_isolation.cd();
    c1_rs_data_isolation.SetLogz();
    h2D_isolation_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_rs_data->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_rs_data->SetTitle(Form("Isolation MVA (RS data) 2016-2018 (corr = %.3f)",h2D_isolation_rs_data->GetCorrelationFactor()));
    h2D_isolation_rs_data->Draw("COLZ");
    c1_rs_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_rs_data_MC.pdf");

    TCanvas c1_rs_data_kinematic;
    c1_rs_data_kinematic.cd();
    c1_rs_data_kinematic.SetLogz();
    h2D_kinematic_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_rs_data->GetYaxis()->SetTitle("BDT2");
    h2D_kinematic_rs_data->SetTitle(Form("Topological MVA (RS data) 2016-2018 (corr = %.3lf)",h2D_kinematic_rs_data->GetCorrelationFactor()));
    h2D_kinematic_rs_data->Draw("COLZ");
    c1_rs_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT2_vs_Bmass_rs_data_MC.pdf");

    TCanvas c2_rs_data_isolation;
    c2_rs_data_isolation.cd();
    TH1D* h_profile_isolation_rs_data = (TH1D*)h2D_isolation_rs_data->ProfileX();
    h_profile_isolation_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_isolation_rs_data->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_isolation_rs_data->SetTitle(Form("Isolation MVA (RS data) 2016-2018 (corr = %.3f)",h2D_isolation_rs_data->GetCorrelationFactor()));
    h_profile_isolation_rs_data->Draw();
    c2_rs_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_profile_rs_data_MC.pdf");

    TCanvas c2_rs_data_kinematic;
    c2_rs_data_kinematic.cd();
    TH1D* h_profile_kinematic_rs_data = (TH1D*)h2D_kinematic_rs_data->ProfileX();
    h_profile_kinematic_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_kinematic_rs_data->GetYaxis()->SetTitle("Average BDT2 per m_{B} bin");
    h_profile_kinematic_rs_data->SetTitle(Form("Topological MVA (RS data) 2016-2018 (corr = %.3f)",h2D_kinematic_rs_data->GetCorrelationFactor()));
    h_profile_kinematic_rs_data->Draw();
    c2_rs_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_profile_rs_data_MC.pdf");

    Double_t isolation_peaks_rs_data[N], kinematic_peaks_rs_data[N];
    Double_t isolation_peaks_error_rs_data[N], kinematic_peaks_error_rs_data[N];

    for(int i = 0; i < N; i++)
    {
        isolation_peaks_rs_data[i] = histos_isolation_rs_data[i]->GetBinCenter(histos_isolation_rs_data[i]->GetMaximumBin());
        kinematic_peaks_rs_data[i] = histos_kinematic_rs_data[i]->GetBinCenter(histos_kinematic_rs_data[i]->GetMaximumBin());

        isolation_peaks_error_rs_data[i] = 20; // half of the bin width
        kinematic_peaks_error_rs_data[i] = 20;

    }

    TGraph* g_isolation_rs_data = new TGraphErrors(N, bdt_cuts, isolation_peaks_rs_data, 0, isolation_peaks_error_rs_data);
    TGraph* g_kinematic_rs_data = new TGraphErrors(N, bdt_cuts, kinematic_peaks_rs_data, 0, kinematic_peaks_error_rs_data);

    TCanvas c3_rs_data_isolation;
    c3_rs_data_isolation.cd();
    c3_rs_data_isolation.SetGrid();
    c3_rs_data_isolation.GetPad(0)->SetLeftMargin(0.15);
    g_isolation_rs_data->SetTitle("Isolation MVA (RS data) 2016-2018 ; BDT1 > x; X-position of B^{+} mass peak (MeV)");
    g_isolation_rs_data->SetMarkerStyle(8);
    g_isolation_rs_data->GetXaxis()->SetLimits(-0.1,1.1);
    g_isolation_rs_data->Draw("AP");
    c3_rs_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt1_rs_data_MC.pdf");

    TCanvas c3_rs_data_kinematic;
    c3_rs_data_kinematic.cd();
    c3_rs_data_kinematic.SetGrid();
    c3_rs_data_kinematic.GetPad(0)->SetLeftMargin(0.15);
    g_kinematic_rs_data->SetTitle("Topological MVA (RS data) 2016-2018 ; BDT2 > x; X-position of B^{+} mass peak (MeV)");
    g_kinematic_rs_data->SetMarkerStyle(8);
    g_kinematic_rs_data->GetXaxis()->SetLimits(-0.1,1.1);
    g_kinematic_rs_data->Draw("AP");
    c3_rs_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt2_rs_data_MC.pdf");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////// WS data /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFileCollection *fc_ws_data_2016 = new TFileCollection("fc_ws_data_2016", "fc_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt");
    TFileCollection *fc_ws_data_2017 = new TFileCollection("fc_ws_data_2017", "fc_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt");
    TFileCollection *fc_ws_data_2018 = new TFileCollection("fc_ws_data_2018", "fc_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt");

    TChain *t_ws_data_2016 = new TChain("DecayTree");
    TChain *t_ws_data_2017 = new TChain("DecayTree");
    TChain *t_ws_data_2018 = new TChain("DecayTree");

    t_ws_data_2016->AddFileInfoList((TCollection*)fc_ws_data_2016->GetList());
    t_ws_data_2017->AddFileInfoList((TCollection*)fc_ws_data_2017->GetList());
    t_ws_data_2018->AddFileInfoList((TCollection*)fc_ws_data_2018->GetList());

    t_ws_data_2016->GetEntries();
    t_ws_data_2017->GetEntries();
    t_ws_data_2018->GetEntries();

    t_ws_data_2016->Add(t_ws_data_2017);
    t_ws_data_2016->Add(t_ws_data_2018);

    TFileCollection *fc1_ws_data_2016 = new TFileCollection("fc1_ws_data_2016", "fc1_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_3/bdt_output.txt");
    TFileCollection *fc1_ws_data_2017 = new TFileCollection("fc1_ws_data_2017", "fc1_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_3/bdt_output.txt");
    TFileCollection *fc1_ws_data_2018 = new TFileCollection("fc1_ws_data_2018", "fc1_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_3/bdt_output.txt");

    TChain *t1_ws_data_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_ws_data_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_ws_data_2018 = new TChain("XGBoost/DecayTree");

    t1_ws_data_2016->AddFileInfoList((TCollection*)fc1_ws_data_2016->GetList());
    t1_ws_data_2017->AddFileInfoList((TCollection*)fc1_ws_data_2017->GetList());
    t1_ws_data_2018->AddFileInfoList((TCollection*)fc1_ws_data_2018->GetList());

    t1_ws_data_2016->GetEntries();
    t1_ws_data_2017->GetEntries();
    t1_ws_data_2018->GetEntries();

    t1_ws_data_2016->Add(t1_ws_data_2017);
    t1_ws_data_2016->Add(t1_ws_data_2018);

    t_ws_data_2016->AddFriend(t1_ws_data_2016);

    std::vector<TH1D*> histos_isolation_ws_data;
    std::vector<TH1D*> histos_kinematic_ws_data;
    TLegend *leg_isolation_ws_data = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *leg_kinematic_ws_data = new TLegend(0.7, 0.7, 0.9, 0.9);

    for(int i = 0; i < N; i++)
    {
        histos_isolation_ws_data.push_back( new TH1D( Form("h_iso_ws_data_%i",i), Form("h_iso_ws_data_%i",i), 100, 4000, 8000) );
        histos_kinematic_ws_data.push_back( new TH1D( Form("h_kin_ws_data_%i",i), Form("h_kin_ws_data_%i",i), 100, 4000, 8000) );

        t_ws_data_2016->Draw(Form("df_Bp_M >> h_iso_ws_data_%i",i), Form("(df_status==0) && (BDT1 > %f)",bdt_cuts[i]));
        t_ws_data_2016->Draw(Form("df_Bp_M >> h_kin_ws_data_%i",i), Form("(df_status==0) && (BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_ws_data[i]->SetLineColor(colors[i]);
        histos_kinematic_ws_data[i]->SetLineColor(colors[i]);

        if(i == 0)
        {
            histos_isolation_ws_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_ws_data[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_ws_data[i]->SetTitle("Isolation MVA (WS data) 2016-2018");

            histos_kinematic_ws_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_ws_data[i]->GetYaxis()->SetTitle("BDT2");
            histos_kinematic_ws_data[i]->SetTitle("Topological MVA (WS data) 2016-2018");
        }
    }

    TCanvas c_ws_data_isolation;
    c_ws_data_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_ws_data[i]->DrawNormalized();
        }
        else
        {
            histos_isolation_ws_data[i]->DrawNormalized("same");
        }
        leg_isolation_ws_data->AddEntry(Form("h_iso_ws_data_%i",i), Form("BDT1 > %.2f",bdt_cuts[i]), "lp");
    }
    leg_isolation_ws_data->Draw("same");
    c_ws_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_ws_data_MC.pdf");

    TCanvas c_ws_data_kinematic;
    c_ws_data_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_ws_data[i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_ws_data[i]->DrawNormalized("same");
        }
        leg_kinematic_ws_data->AddEntry(Form("h_kin_ws_data_%i",i), Form("BDT2 > %.2f",bdt_cuts[i]), "lp");
    }
    leg_kinematic_ws_data->Draw("same");
    c_ws_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_ws_data_MC.pdf");

    TH2D* h2D_isolation_ws_data = new TH2D("h2d_iso_ws_data", "h2d_iso_ws_data", 100, 4000, 8000, 100, 0, 1);
    TH2D* h2D_kinematic_ws_data = new TH2D("h2d_kin_ws_data", "h2d_kin_ws_data", 100, 4000, 8000, 100, 0, 1);

    t_ws_data_2016->Draw("BDT1 : df_Bp_M >> h2d_iso_ws_data", "df_status==0");
    t_ws_data_2016->Draw("BDT2 : df_Bp_M >> h2d_kin_ws_data", "df_status==0");

    TCanvas c1_ws_data_isolation;
    c1_ws_data_isolation.cd();
    c1_ws_data_isolation.SetLogz();
    h2D_isolation_ws_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_ws_data->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_ws_data->SetTitle(Form("Isolation MVA (WS data) 2016-2018 (corr = %.3f)",h2D_isolation_ws_data->GetCorrelationFactor()));
    h2D_isolation_ws_data->Draw("COLZ");
    c1_ws_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_ws_data_MC.pdf");

    TCanvas c1_ws_data_kinematic;
    c1_ws_data_kinematic.cd();
    c1_ws_data_kinematic.SetLogz();
    h2D_kinematic_ws_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_ws_data->GetYaxis()->SetTitle("BDT2");
    h2D_kinematic_ws_data->SetTitle(Form("Topological MVA (WS data) 2016-2018 (corr = %.3lf)",h2D_kinematic_ws_data->GetCorrelationFactor()));
    h2D_kinematic_ws_data->Draw("COLZ");
    c1_ws_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT2_vs_Bmass_ws_data_MC.pdf");

    TCanvas c2_ws_data_isolation;
    c2_ws_data_isolation.cd();
    TH1D* h_profile_isolation_ws_data = (TH1D*)h2D_isolation_ws_data->ProfileX();
    h_profile_isolation_ws_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_isolation_ws_data->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_isolation_ws_data->SetTitle(Form("Isolation MVA (WS data) 2016-2018 (corr = %.3f)",h2D_isolation_ws_data->GetCorrelationFactor()));
    h_profile_isolation_ws_data->Draw();
    c2_ws_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_profile_ws_data_MC.pdf");

    TCanvas c2_ws_data_kinematic;
    c2_ws_data_kinematic.cd();
    TH1D* h_profile_kinematic_ws_data = (TH1D*)h2D_kinematic_ws_data->ProfileX();
    h_profile_kinematic_ws_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_kinematic_ws_data->GetYaxis()->SetTitle("Average BDT2 per m_{B} bin");
    h_profile_kinematic_ws_data->SetTitle(Form("Topological MVA (WS data) 2016-2018 (corr = %.3f)",h2D_kinematic_ws_data->GetCorrelationFactor()));
    h_profile_kinematic_ws_data->Draw();
    c2_ws_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_profile_ws_data_MC.pdf");

    Double_t isolation_peaks_ws_data[N], kinematic_peaks_ws_data[N];
    Double_t isolation_peaks_error_ws_data[N], kinematic_peaks_error_ws_data[N];

    for(int i = 0; i < N; i++)
    {
        isolation_peaks_ws_data[i] = histos_isolation_ws_data[i]->GetBinCenter(histos_isolation_ws_data[i]->GetMaximumBin());
        kinematic_peaks_ws_data[i] = histos_kinematic_ws_data[i]->GetBinCenter(histos_kinematic_ws_data[i]->GetMaximumBin());

        isolation_peaks_error_ws_data[i] = 20; // half of the bin width
        kinematic_peaks_error_ws_data[i] = 20;

    }

    TGraph* g_isolation_ws_data = new TGraphErrors(N, bdt_cuts, isolation_peaks_ws_data, 0, isolation_peaks_error_ws_data);
    TGraph* g_kinematic_ws_data = new TGraphErrors(N, bdt_cuts, kinematic_peaks_ws_data, 0, kinematic_peaks_error_ws_data);

    TCanvas c3_ws_data_isolation;
    c3_ws_data_isolation.cd();
    c3_ws_data_isolation.SetGrid();
    c3_ws_data_isolation.GetPad(0)->SetLeftMargin(0.15);
    g_isolation_ws_data->SetTitle("Isolation MVA (WS data) 2016-2018 ; BDT1 > x; X-position of B^{+} mass peak (MeV)");
    g_isolation_ws_data->SetMarkerStyle(8);
    g_isolation_ws_data->GetXaxis()->SetLimits(-0.1,1.1);
    g_isolation_ws_data->Draw("AP");
    c3_ws_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt1_ws_data_MC.pdf");

    TCanvas c3_ws_data_kinematic;
    c3_ws_data_kinematic.cd();
    c3_ws_data_kinematic.SetGrid();
    c3_ws_data_kinematic.GetPad(0)->SetLeftMargin(0.15);
    g_kinematic_ws_data->SetTitle("Topological MVA (WS data) 2016-2018 ; BDT2 > x; X-position of B^{+} mass peak (MeV)");
    g_kinematic_ws_data->SetMarkerStyle(8);
    g_kinematic_ws_data->GetXaxis()->SetLimits(-0.1,1.1);
    g_kinematic_ws_data->Draw("AP");
    c3_ws_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt2_ws_data_MC.pdf");

    // h_3pi3pi_iso
    TH1D* h_rs_data_iso = new TH1D("h_rs_data_iso", "h_rs_data_iso", 100, 0, 1);
    TH1D* h_ws_data_iso = new TH1D("h_ws_data_iso", "h_ws_data_iso", 100, 0, 1);

    t_rs_data_2016->Draw("BDT1 >> h_rs_data_iso", "df_status==0");
    t_ws_data_2016->Draw("BDT1 >> h_ws_data_iso", "df_status==0");

    TCanvas c6;
    c6.cd();
    h_ws_data_iso->GetXaxis()->SetTitle("BDT1");
    h_ws_data_iso->GetYaxis()->SetTitle("Normalised entries / (0.025)");
    h_ws_data_iso->SetTitle("Isolation MVA 2016-2018");
    h_3pi3pi_iso->SetLineColor(kBlue);
    h_rs_data_iso->SetLineColor(kBlack);
    h_ws_data_iso->SetLineColor(kRed);
    h_ws_data_iso->DrawNormalized();
    h_rs_data_iso->DrawNormalized("same");
    h_3pi3pi_iso->DrawNormalized("same");
    TLegend *leg1_comp_iso = new TLegend(0.3, 0.7, 0.6, 0.9);
    leg1_comp_iso->AddEntry(h_3pi3pi_iso, "3#pi3#pi MC", "lp");
    leg1_comp_iso->AddEntry(h_rs_data_iso, "RS data", "lp");
    leg1_comp_iso->AddEntry(h_ws_data_iso, "WS data", "lp");
    leg1_comp_iso->Draw("same");
    c6.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_mc_data.pdf");

    TH1D* h_rs_data_kin = new TH1D("h_rs_data_kin", "h_rs_data_kin", 100, 0, 1);
    TH1D* h_ws_data_kin = new TH1D("h_ws_data_kin", "h_ws_data_kin", 100, 0, 1);

    t_rs_data_2016->Draw("BDT2 >> h_rs_data_kin", "df_status==0");
    t_ws_data_2016->Draw("BDT2 >> h_ws_data_kin", "df_status==0");

    TCanvas c7;
    c7.cd();
    h_ws_data_kin->GetXaxis()->SetTitle("BDT2");
    h_ws_data_kin->GetYaxis()->SetTitle("Normalised entries / (0.025)");
    h_ws_data_kin->SetTitle("Topological MVA 2016-2018");
    h_3pi3pi_kin->SetLineColor(kBlue);
    h_rs_data_kin->SetLineColor(kBlack);
    h_ws_data_kin->SetLineColor(kRed);
    h_ws_data_kin->DrawNormalized();
    h_rs_data_kin->DrawNormalized("same");
    h_3pi3pi_kin->DrawNormalized("same");
    TLegend *leg1_comp_kin = new TLegend(0.3, 0.7, 0.6, 0.9);
    leg1_comp_kin->AddEntry(h_3pi3pi_kin, "3#pi3#pi MC", "lp");
    leg1_comp_kin->AddEntry(h_rs_data_kin, "RS data", "lp");
    leg1_comp_kin->AddEntry(h_ws_data_kin, "WS data", "lp");
    leg1_comp_kin->Draw("same");
    c7.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_mc_data.pdf");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////////////////// Normalisation channel: MC (loose rectangular cuts) /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFileCollection *fc_norm_mc_2016 = new TFileCollection("fc_norm_mc_2016", "fc_norm_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_71/pre_sel_tree.txt");
    TFileCollection *fc_norm_mc_2017 = new TFileCollection("fc_norm_mc_2017", "fc_norm_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_71/pre_sel_tree.txt");
    TFileCollection *fc_norm_mc_2018 = new TFileCollection("fc_norm_mc_2018", "fc_norm_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_71/pre_sel_tree.txt");

    TChain *t_norm_mc_2016 = new TChain("DecayTree");
    TChain *t_norm_mc_2017 = new TChain("DecayTree");
    TChain *t_norm_mc_2018 = new TChain("DecayTree");

    t_norm_mc_2016->AddFileInfoList((TCollection*)fc_norm_mc_2016->GetList());
    t_norm_mc_2017->AddFileInfoList((TCollection*)fc_norm_mc_2017->GetList());
    t_norm_mc_2018->AddFileInfoList((TCollection*)fc_norm_mc_2018->GetList());

    t_norm_mc_2016->GetEntries();
    t_norm_mc_2017->GetEntries();
    t_norm_mc_2018->GetEntries();

    t_norm_mc_2016->Add(t_norm_mc_2017);
    t_norm_mc_2016->Add(t_norm_mc_2018);

    TFileCollection *fc1_norm_mc_2016 = new TFileCollection("fc1_norm_mc_2016", "fc1_norm_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_71/bdt_output.txt");
    TFileCollection *fc1_norm_mc_2017 = new TFileCollection("fc1_norm_mc_2017", "fc1_norm_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_71/bdt_output.txt");
    TFileCollection *fc1_norm_mc_2018 = new TFileCollection("fc1_norm_mc_2018", "fc1_norm_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_71/bdt_output.txt");

    TChain *t1_norm_mc_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_norm_mc_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_norm_mc_2018 = new TChain("XGBoost/DecayTree");

    t1_norm_mc_2016->AddFileInfoList((TCollection*)fc1_norm_mc_2016->GetList());
    t1_norm_mc_2017->AddFileInfoList((TCollection*)fc1_norm_mc_2017->GetList());
    t1_norm_mc_2018->AddFileInfoList((TCollection*)fc1_norm_mc_2018->GetList());

    t1_norm_mc_2016->GetEntries();
    t1_norm_mc_2017->GetEntries();
    t1_norm_mc_2018->GetEntries();

    t1_norm_mc_2016->Add(t1_norm_mc_2017);
    t1_norm_mc_2016->Add(t1_norm_mc_2018);

    t_norm_mc_2016->AddFriend(t1_norm_mc_2016);

    std::vector<TH1D*> histos_isolation_norm_mc;
    std::vector<TH1D*> histos_kinematic_norm_mc;
    TLegend *leg_isolation_norm_mc = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *leg_kinematic_norm_mc = new TLegend(0.7, 0.7, 0.9, 0.9);

    for(int i = 0; i < N; i++)
    {
        histos_isolation_norm_mc.push_back( new TH1D( Form("h_iso_norm_mc_%i",i), Form("h_iso_norm_mc_%i",i), 100, 5235, 5355) );
        histos_kinematic_norm_mc.push_back( new TH1D( Form("h_kin_norm_mc_%i",i), Form("h_kin_norm_mc_%i",i), 100, 5235, 5355) );

        t_norm_mc_2016->Draw(Form("Bp_dtf_M[0] >> h_iso_norm_mc_%i",i), Form("(BDT1 > %f)",bdt_cuts[i]));
        t_norm_mc_2016->Draw(Form("Bp_dtf_M[0] >> h_kin_norm_mc_%i",i), Form("(BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_norm_mc[i]->SetLineColor(colors[i]);
        histos_kinematic_norm_mc[i]->SetLineColor(colors[i]);

        if(i == N-1)
        {
            histos_isolation_norm_mc[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_norm_mc[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_norm_mc[i]->SetTitle("Isolation MVA (Norm. MC) 2016-2018");

            histos_kinematic_norm_mc[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_norm_mc[i]->GetYaxis()->SetTitle("BDT2");
            histos_kinematic_norm_mc[i]->SetTitle("Topological MVA (Norm. MC) 2016-2018");
        }
    }

    TCanvas c_norm_mc_isolation;
    c_norm_mc_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_norm_mc[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_isolation_norm_mc[N-1-i]->DrawNormalized("same");
        }
        leg_isolation_norm_mc->AddEntry(Form("h_iso_norm_mc_%i",N-1-i), Form("BDT1 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_isolation_norm_mc->Draw("same");
    c_norm_mc_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_norm_MC.pdf");

    TCanvas c_norm_mc_kinematic;
    c_norm_mc_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_norm_mc[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_norm_mc[N-1-i]->DrawNormalized("same");
        }
        leg_kinematic_norm_mc->AddEntry(Form("h_kin_norm_mc_%i",N-1-i), Form("BDT2 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_kinematic_norm_mc->Draw("same");
    c_norm_mc_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_norm_MC.pdf");

    TH2D* h2D_isolation_norm_mc = new TH2D("h2d_iso_norm_mc", "h2d_iso_norm_mc", 100, 5235, 5355, 100, 0, 1);
    TH2D* h2D_kinematic_norm_mc = new TH2D("h2d_kin_norm_mc", "h2d_kin_norm_mc", 100, 5235, 5355, 100, 0, 1);

    t_norm_mc_2016->Draw("BDT1 : Bp_dtf_M[0] >> h2d_iso_norm_mc");
    t_norm_mc_2016->Draw("BDT2 : Bp_dtf_M[0] >> h2d_kin_norm_mc");

    TCanvas c1_norm_mc_isolation;
    c1_norm_mc_isolation.cd();
    c1_norm_mc_isolation.SetLogz();
    h2D_isolation_norm_mc->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_norm_mc->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_norm_mc->SetTitle(Form("Isolation MVA (Norm. MC) 2016-2018 (corr = %.3f)",h2D_isolation_norm_mc->GetCorrelationFactor()));
    h2D_isolation_norm_mc->Draw("COLZ");
    c1_norm_mc_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_norm_MC.pdf");

    TCanvas c1_norm_mc_kinematic;
    c1_norm_mc_kinematic.cd();
    c1_norm_mc_kinematic.SetLogz();
    h2D_kinematic_norm_mc->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_norm_mc->GetYaxis()->SetTitle("BDT2");
    h2D_kinematic_norm_mc->SetTitle(Form("Topological MVA (Norm. MC) 2016-2018 (corr = %.3lf)",h2D_kinematic_norm_mc->GetCorrelationFactor()));
    h2D_kinematic_norm_mc->Draw("COLZ");
    c1_norm_mc_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT2_vs_Bmass_norm_MC.pdf");

    TCanvas c2_norm_mc_isolation;
    c2_norm_mc_isolation.cd();
    TH1D* h_profile_isolation_norm_mc = (TH1D*)h2D_isolation_norm_mc->ProfileX();
    h_profile_isolation_norm_mc->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_isolation_norm_mc->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_isolation_norm_mc->SetTitle(Form("Isolation MVA (Norm. MC) 2016-2018 (corr = %.3f)",h2D_isolation_norm_mc->GetCorrelationFactor()));
    h_profile_isolation_norm_mc->Draw();
    c2_norm_mc_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_profile_norm_MC.pdf");

    TCanvas c2_norm_mc_kinematic;
    c2_norm_mc_kinematic.cd();
    TH1D* h_profile_kinematic_norm_mc = (TH1D*)h2D_kinematic_norm_mc->ProfileX();
    h_profile_kinematic_norm_mc->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_kinematic_norm_mc->GetYaxis()->SetTitle("Average BDT2 per m_{B} bin");
    h_profile_kinematic_norm_mc->SetTitle(Form("Topological MVA (Norm. MC) 2016-2018 (corr = %.3f)",h2D_kinematic_norm_mc->GetCorrelationFactor()));
    h_profile_kinematic_norm_mc->Draw();
    c2_norm_mc_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_profile_norm_MC.pdf");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////// Normalisation channel: data (loose rectangular cuts) /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFileCollection *fc_norm_data_2016 = new TFileCollection("fc_norm_data_2016", "fc_norm_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_81/pre_sel_tree.txt");
    TFileCollection *fc_norm_data_2017 = new TFileCollection("fc_norm_data_2017", "fc_norm_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_81/pre_sel_tree.txt");
    TFileCollection *fc_norm_data_2018 = new TFileCollection("fc_norm_data_2018", "fc_norm_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_81/pre_sel_tree.txt");

    TChain *t_norm_data_2016 = new TChain("DecayTree");
    TChain *t_norm_data_2017 = new TChain("DecayTree");
    TChain *t_norm_data_2018 = new TChain("DecayTree");

    t_norm_data_2016->AddFileInfoList((TCollection*)fc_norm_data_2016->GetList());
    t_norm_data_2017->AddFileInfoList((TCollection*)fc_norm_data_2017->GetList());
    t_norm_data_2018->AddFileInfoList((TCollection*)fc_norm_data_2018->GetList());

    t_norm_data_2016->GetEntries();
    t_norm_data_2017->GetEntries();
    t_norm_data_2018->GetEntries();

    t_norm_data_2016->Add(t_norm_data_2017);
    t_norm_data_2016->Add(t_norm_data_2018);

    TFileCollection *fc1_norm_data_2016 = new TFileCollection("fc1_norm_data_2016", "fc1_norm_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_81/bdt_output.txt");
    TFileCollection *fc1_norm_data_2017 = new TFileCollection("fc1_norm_data_2017", "fc1_norm_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_81/bdt_output.txt");
    TFileCollection *fc1_norm_data_2018 = new TFileCollection("fc1_norm_data_2018", "fc1_norm_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_81/bdt_output.txt");

    TChain *t1_norm_data_2016 = new TChain("XGBoost/DecayTree");
    TChain *t1_norm_data_2017 = new TChain("XGBoost/DecayTree");
    TChain *t1_norm_data_2018 = new TChain("XGBoost/DecayTree");

    t1_norm_data_2016->AddFileInfoList((TCollection*)fc1_norm_data_2016->GetList());
    t1_norm_data_2017->AddFileInfoList((TCollection*)fc1_norm_data_2017->GetList());
    t1_norm_data_2018->AddFileInfoList((TCollection*)fc1_norm_data_2018->GetList());

    t1_norm_data_2016->GetEntries();
    t1_norm_data_2017->GetEntries();
    t1_norm_data_2018->GetEntries();

    t1_norm_data_2016->Add(t1_norm_data_2017);
    t1_norm_data_2016->Add(t1_norm_data_2018);

    t_norm_data_2016->AddFriend(t1_norm_data_2016);

    std::vector<TH1D*> histos_isolation_norm_data;
    std::vector<TH1D*> histos_kinematic_norm_data;
    TLegend *leg_isolation_norm_data = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *leg_kinematic_norm_data = new TLegend(0.7, 0.7, 0.9, 0.9);

    for(int i = 0; i < N; i++)
    {
        histos_isolation_norm_data.push_back( new TH1D( Form("h_iso_norm_data_%i",i), Form("h_iso_norm_data_%i",i), 100, 5235, 5355) );
        histos_kinematic_norm_data.push_back( new TH1D( Form("h_kin_norm_data_%i",i), Form("h_kin_norm_data_%i",i), 100, 5235, 5355) );

        t_norm_data_2016->Draw(Form("Bp_dtf_M[0] >> h_iso_norm_data_%i",i), Form("(BDT1 > %f)",bdt_cuts[i]));
        t_norm_data_2016->Draw(Form("Bp_dtf_M[0] >> h_kin_norm_data_%i",i), Form("(BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_norm_data[i]->SetLineColor(colors[i]);
        histos_kinematic_norm_data[i]->SetLineColor(colors[i]);

        if(i == N-1)
        {
            histos_isolation_norm_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_norm_data[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_norm_data[i]->SetTitle("Isolation MVA (Norm. data) 2016-2018");

            histos_kinematic_norm_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_norm_data[i]->GetYaxis()->SetTitle("BDT2");
            histos_kinematic_norm_data[i]->SetTitle("Topological MVA (Norm. data) 2016-2018");
        }
    }

    TCanvas c_norm_data_isolation;
    c_norm_data_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_norm_data[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_isolation_norm_data[N-1-i]->DrawNormalized("same");
        }
        leg_isolation_norm_data->AddEntry(Form("h_iso_norm_data_%i",N-1-i), Form("BDT1 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_isolation_norm_data->Draw("same");
    c_norm_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_norm_data.pdf");

    TCanvas c_norm_data_kinematic;
    c_norm_data_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_norm_data[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_norm_data[N-1-i]->DrawNormalized("same");
        }
        leg_kinematic_norm_data->AddEntry(Form("h_kin_norm_data_%i",N-1-i), Form("BDT2 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_kinematic_norm_data->Draw("same");
    c_norm_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_norm_data.pdf");

    TH2D* h2D_isolation_norm_data = new TH2D("h2d_iso_norm_data", "h2d_iso_norm_data", 100, 5235, 5355, 100, 0, 1);
    TH2D* h2D_kinematic_norm_data = new TH2D("h2d_kin_norm_data", "h2d_kin_norm_data", 100, 5235, 5355, 100, 0, 1);

    t_norm_data_2016->Draw("BDT1 : Bp_dtf_M[0] >> h2d_iso_norm_data");
    t_norm_data_2016->Draw("BDT2 : Bp_dtf_M[0] >> h2d_kin_norm_data");

    TCanvas c1_norm_data_isolation;
    c1_norm_data_isolation.cd();
    c1_norm_data_isolation.SetLogz();
    h2D_isolation_norm_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_norm_data->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_norm_data->SetTitle(Form("Isolation MVA (Norm. data) 2016-2018 (corr = %.3f)",h2D_isolation_norm_data->GetCorrelationFactor()));
    h2D_isolation_norm_data->Draw("COLZ");
    c1_norm_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_norm_data.pdf");

    TCanvas c1_norm_data_kinematic;
    c1_norm_data_kinematic.cd();
    c1_norm_data_kinematic.SetLogz();
    h2D_kinematic_norm_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_norm_data->GetYaxis()->SetTitle("BDT2");
    h2D_kinematic_norm_data->SetTitle(Form("Topological MVA (Norm. data) 2016-2018 (corr = %.3lf)",h2D_kinematic_norm_data->GetCorrelationFactor()));
    h2D_kinematic_norm_data->Draw("COLZ");
    c1_norm_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT2_vs_Bmass_norm_data.pdf");

    TCanvas c2_norm_data_isolation;
    c2_norm_data_isolation.cd();
    TH1D* h_profile_isolation_norm_data = (TH1D*)h2D_isolation_norm_data->ProfileX();
    h_profile_isolation_norm_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_isolation_norm_data->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_isolation_norm_data->SetTitle(Form("Isolation MVA (Norm. data) 2016-2018 (corr = %.3f)",h2D_isolation_norm_data->GetCorrelationFactor()));
    h_profile_isolation_norm_data->Draw();
    c2_norm_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_profile_norm_data.pdf");

    TCanvas c2_norm_data_kinematic;
    c2_norm_data_kinematic.cd();
    TH1D* h_profile_kinematic_norm_data = (TH1D*)h2D_kinematic_norm_data->ProfileX();
    h_profile_kinematic_norm_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_kinematic_norm_data->GetYaxis()->SetTitle("Average BDT2 per m_{B} bin");
    h_profile_kinematic_norm_data->SetTitle(Form("Topological MVA (Norm. data) 2016-2018 (corr = %.3f)",h2D_kinematic_norm_data->GetCorrelationFactor()));
    h_profile_kinematic_norm_data->Draw();
    c2_norm_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_profile_norm_data.pdf");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////// Normalisation channel: background M > 5320 MeV (loose rectangular cuts) /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<TH1D*> histos_isolation_norm_bkg;
    std::vector<TH1D*> histos_kinematic_norm_bkg;
    TLegend *leg_isolation_norm_bkg = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *leg_kinematic_norm_bkg = new TLegend(0.7, 0.7, 0.9, 0.9);

    for(int i = 0; i < N; i++)
    {
        histos_isolation_norm_bkg.push_back( new TH1D( Form("h_iso_norm_bkg_%i",i), Form("h_iso_norm_bkg_%i",i), 100, 5320, 5355) );
        histos_kinematic_norm_bkg.push_back( new TH1D( Form("h_kin_norm_bkg_%i",i), Form("h_kin_norm_bkg_%i",i), 100, 5320, 5355) );

        t_norm_data_2016->Draw(Form("Bp_dtf_M[0] >> h_iso_norm_bkg_%i",i), Form("(BDT1 > %f) && (Bp_dtf_M[0] > 5320)",bdt_cuts[i]));
        t_norm_data_2016->Draw(Form("Bp_dtf_M[0] >> h_kin_norm_bkg_%i",i), Form("(BDT2 > %f) && (Bp_dtf_M[0] > 5320)",bdt_cuts[i]));

        histos_isolation_norm_bkg[i]->SetLineColor(colors[i]);
        histos_kinematic_norm_bkg[i]->SetLineColor(colors[i]);

        if(i == N-1)
        {
            histos_isolation_norm_bkg[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_norm_bkg[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_norm_bkg[i]->SetTitle("Isolation MVA (Norm. bkg) 2016-2018");

            histos_kinematic_norm_bkg[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_norm_bkg[i]->GetYaxis()->SetTitle("BDT2");
            histos_kinematic_norm_bkg[i]->SetTitle("Topological MVA (Norm. bkg) 2016-2018");
        }
    }

    TCanvas c_norm_bkg_isolation;
    c_norm_bkg_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_norm_bkg[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_isolation_norm_bkg[N-1-i]->DrawNormalized("same");
        }
        leg_isolation_norm_bkg->AddEntry(Form("h_iso_norm_bkg_%i",N-1-i), Form("BDT1 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_isolation_norm_bkg->Draw("same");
    c_norm_bkg_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_norm_bkg.pdf");

    TCanvas c_norm_bkg_kinematic;
    c_norm_bkg_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_norm_bkg[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_norm_bkg[N-1-i]->DrawNormalized("same");
        }
        leg_kinematic_norm_bkg->AddEntry(Form("h_kin_norm_bkg_%i",N-1-i), Form("BDT2 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_kinematic_norm_bkg->Draw("same");
    c_norm_bkg_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_norm_bkg.pdf");

    TH2D* h2D_isolation_norm_bkg = new TH2D("h2d_iso_norm_bkg", "h2d_iso_norm_bkg", 100, 5320, 5355, 100, 0, 1);
    TH2D* h2D_kinematic_norm_bkg = new TH2D("h2d_kin_norm_bkg", "h2d_kin_norm_bkg", 100, 5320, 5355, 100, 0, 1);

    t_norm_data_2016->Draw("BDT1 : Bp_dtf_M[0] >> h2d_iso_norm_bkg", "Bp_dtf_M[0] > 5320");
    t_norm_data_2016->Draw("BDT2 : Bp_dtf_M[0] >> h2d_kin_norm_bkg", "Bp_dtf_M[0] > 5320");

    TCanvas c1_norm_bkg_isolation;
    c1_norm_bkg_isolation.cd();
    c1_norm_bkg_isolation.SetLogz();
    h2D_isolation_norm_bkg->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_norm_bkg->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_norm_bkg->SetTitle(Form("Isolation MVA (Norm. bkg) 2016-2018 (corr = %.3f)",h2D_isolation_norm_bkg->GetCorrelationFactor()));
    h2D_isolation_norm_bkg->Draw("COLZ");
    c1_norm_bkg_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_norm_bkg.pdf");

    TCanvas c1_norm_bkg_kinematic;
    c1_norm_bkg_kinematic.cd();
    c1_norm_bkg_kinematic.SetLogz();
    h2D_kinematic_norm_bkg->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_norm_bkg->GetYaxis()->SetTitle("BDT2");
    h2D_kinematic_norm_bkg->SetTitle(Form("Topological MVA (Norm. bkg) 2016-2018 (corr = %.3lf)",h2D_kinematic_norm_bkg->GetCorrelationFactor()));
    h2D_kinematic_norm_bkg->Draw("COLZ");
    c1_norm_bkg_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT2_vs_Bmass_norm_bkg.pdf");

    TCanvas c2_norm_bkg_isolation;
    c2_norm_bkg_isolation.cd();
    TH1D* h_profile_isolation_norm_bkg = (TH1D*)h2D_isolation_norm_bkg->ProfileX();
    h_profile_isolation_norm_bkg->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_isolation_norm_bkg->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_isolation_norm_bkg->SetTitle(Form("Isolation MVA (Norm. bkg) 2016-2018 (corr = %.3f)",h2D_isolation_norm_bkg->GetCorrelationFactor()));
    h_profile_isolation_norm_bkg->Draw();
    c2_norm_bkg_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_profile_norm_bkg.pdf");

    TCanvas c2_norm_bkg_kinematic;
    c2_norm_bkg_kinematic.cd();
    TH1D* h_profile_kinematic_norm_bkg = (TH1D*)h2D_kinematic_norm_bkg->ProfileX();
    h_profile_kinematic_norm_bkg->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_profile_kinematic_norm_bkg->GetYaxis()->SetTitle("Average BDT2 per m_{B} bin");
    h_profile_kinematic_norm_bkg->SetTitle(Form("Topological MVA (Norm. bkg) 2016-2018 (corr = %.3f)",h2D_kinematic_norm_bkg->GetCorrelationFactor()));
    h_profile_kinematic_norm_bkg->Draw();
    c2_norm_bkg_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_profile_norm_bkg.pdf");

    TH1D* h_norm_mc_iso = new TH1D("h_norm_mc_iso", "h_norm_mc_iso", 100, 0, 1);
    TH1D* h_norm_data_iso = new TH1D("h_norm_data_iso", "h_norm_data_iso", 100, 0, 1);
    TH1D* h_norm_bkg_iso = new TH1D("h_norm_bkg_iso", "h_norm_bkg_iso", 100, 0, 1);

    t_norm_mc_2016->Draw("BDT1 >> h_norm_mc_iso");
    t_norm_data_2016->Draw("BDT1 >> h_norm_data_iso");
    t_norm_data_2016->Draw("BDT1 >> h_norm_bkg_iso", "Bp_dtf_M[0] > 5320");

    TCanvas c8;
    c8.cd();
    h_norm_bkg_iso->GetXaxis()->SetTitle("BDT1");
    h_norm_bkg_iso->GetYaxis()->SetTitle("Normalised entries / (0.025)");
    h_norm_bkg_iso->SetTitle("Isolation MVA 2016-2018");
    h_norm_mc_iso->SetLineColor(kBlue);
    h_norm_data_iso->SetLineColor(kBlack);
    h_norm_bkg_iso->SetLineColor(kRed);
    h_norm_bkg_iso->DrawNormalized();
    h_norm_data_iso->DrawNormalized("same");
    h_norm_mc_iso->DrawNormalized("same");
    TLegend *leg2_comp_iso = new TLegend(0.3, 0.7, 0.6, 0.9);
    leg2_comp_iso->AddEntry(h_3pi3pi_iso, "Norm. MC", "lp");
    leg2_comp_iso->AddEntry(h_rs_data_iso, "Norm. data", "lp");
    leg2_comp_iso->AddEntry(h_ws_data_iso, "Norm. bkg.", "lp");
    leg2_comp_iso->Draw("same");
    c8.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt1_norm_mc_data.pdf");

    TH1D* h_norm_mc_kin = new TH1D("h_norm_mc_kin", "h_norm_mc_kin", 100, 0, 1);
    TH1D* h_norm_data_kin = new TH1D("h_norm_data_kin", "h_norm_data_kin", 100, 0, 1);
    TH1D* h_norm_bkg_kin = new TH1D("h_norm_bkg_kin", "h_norm_bkg_kin", 100, 0, 1);

    t_norm_mc_2016->Draw("BDT2 >> h_norm_mc_kin");
    t_norm_data_2016->Draw("BDT2 >> h_norm_data_kin");
    t_norm_data_2016->Draw("BDT2 >> h_norm_bkg_kin", "Bp_dtf_M[0] > 5320");

    TCanvas c9;
    c9.cd();
    h_norm_bkg_kin->GetXaxis()->SetTitle("BDT2");
    h_norm_bkg_kin->GetYaxis()->SetTitle("Normalised entries / (0.025)");
    h_norm_bkg_kin->SetTitle("Topological MVA 2016-2018");
    h_norm_mc_kin->SetLineColor(kBlue);
    h_norm_data_kin->SetLineColor(kBlack);
    h_norm_bkg_kin->SetLineColor(kRed);
    h_norm_bkg_kin->DrawNormalized();
    h_norm_data_kin->DrawNormalized("same");
    h_norm_mc_kin->DrawNormalized("same");
    TLegend *leg2_comp_kin = new TLegend(0.3, 0.7, 0.6, 0.9);
    leg2_comp_kin->AddEntry(h_3pi3pi_kin, "Norm. MC", "lp");
    leg2_comp_kin->AddEntry(h_rs_data_kin, "Norm. data", "lp");
    leg2_comp_kin->AddEntry(h_ws_data_kin, "Norm. bkg.", "lp");
    leg2_comp_kin->Draw("same");
    c9.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/bdt2_norm_mc_data.pdf");

}
