

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
        histos_isolation_3pi3pi.push_back( new TH1D( Form("h_iso_%i",i), Form("h_iso_%i",i), 100, 4000, 8000) );
        histos_kinematic_3pi3pi.push_back( new TH1D( Form("h_kin_%i",i), Form("h_kin_%i",i), 100, 4000, 8000) );

        t_3pi3pi_2016->Draw(Form("df_Bp_M >> h_iso_%i",i), Form("(df_status==0) && (BDT1 > %f)",bdt_cuts[i]));
        t_3pi3pi_2016->Draw(Form("df_Bp_M >> h_kin_%i",i), Form("(df_status==0) && (BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_3pi3pi[i]->SetLineColor(colors[i]);
        histos_kinematic_3pi3pi[i]->SetLineColor(colors[i]);

        if(i == N-1)
        {
            histos_isolation_3pi3pi[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_3pi3pi[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_3pi3pi[i]->SetTitle("Isolation MVA (3#pi3#pi MC) 2016-2018");

            histos_kinematic_3pi3pi[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_3pi3pi[i]->GetYaxis()->SetTitle("BDT1");
            histos_kinematic_3pi3pi[i]->SetTitle("Kinematic MVA (3#pi3#pi MC) 2016-2018");
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
        leg_isolation_3pi3pi->AddEntry(Form("h_iso_%i",N-1-i), Form("BDT1 > %.2f",bdt_cuts[N-1-i]), "lp");
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
        leg_kinematic_3pi3pi->AddEntry(Form("h_kin_%i",N-1-i), Form("BDT2 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_kinematic_3pi3pi->Draw("same");
    c_3pi3pi_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_3pi3pi_MC.pdf");

    TH2D* h2D_isolation_3pi3pi = new TH2D("h2d_iso", "h2d_iso", 100, 4000, 8000, 100, 0, 1);
    TH2D* h2D_kinematic_3pi3pi = new TH2D("h2d_kin", "h2d_kin", 100, 4000, 8000, 100, 0, 1);

    t_3pi3pi_2016->Draw("BDT1 : df_Bp_M >> h2d_iso", "df_status==0");
    t_3pi3pi_2016->Draw("BDT2 : df_Bp_M >> h2d_kin", "df_status==0");

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
    h2D_kinematic_3pi3pi->GetYaxis()->SetTitle("BDT1");
    h2D_kinematic_3pi3pi->SetTitle(Form("Kinematic MVA (3#pi3#pi MC) 2016-2018 (corr = %.3lf)",h2D_kinematic_3pi3pi->GetCorrelationFactor()));
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
    h_profile_kinematic_3pi3pi->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_kinematic_3pi3pi->SetTitle(Form("Kinematic MVA (3#pi3#pi MC) 2016-2018 (corr = %.3f)",h2D_kinematic_3pi3pi->GetCorrelationFactor()));
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
    g_kinematic_3pi3pi->SetTitle("Kinematic MVA (3#pi3#pi MC) 2016-2018 ; BDT2 > x; X-position of B^{+} mass peak (MeV)");
    g_kinematic_3pi3pi->SetMarkerStyle(8);
    g_kinematic_3pi3pi->GetXaxis()->SetLimits(-0.1,1.1);
    g_kinematic_3pi3pi->Draw("AP");
    c3_3pi3pi_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt2_3pi3pi_MC.pdf");
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
        histos_isolation_rs_data.push_back( new TH1D( Form("h_iso_%i",i), Form("h_iso_%i",i), 100, 4000, 8000) );
        histos_kinematic_rs_data.push_back( new TH1D( Form("h_kin_%i",i), Form("h_kin_%i",i), 100, 4000, 8000) );

        t_rs_data_2016->Draw(Form("df_Bp_M >> h_iso_%i",i), Form("(df_status==0) && ((df_Bp_M < 4700) || (df_Bp_M > 5800)) && (BDT1 > %f)",bdt_cuts[i]));
        t_rs_data_2016->Draw(Form("df_Bp_M >> h_kin_%i",i), Form("(df_status==0) && ((df_Bp_M < 4700) || (df_Bp_M > 5800)) && (BDT2 > %f)",bdt_cuts[i]));

        histos_isolation_rs_data[i]->SetLineColor(colors[i]);
        histos_kinematic_rs_data[i]->SetLineColor(colors[i]);

        if(i == N-1)
        {
            histos_isolation_rs_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_isolation_rs_data[i]->GetYaxis()->SetTitle("BDT1");
            histos_isolation_rs_data[i]->SetTitle("Isolation MVA (RS data) 2016-2018");

            histos_kinematic_rs_data[i]->GetXaxis()->SetTitle("m_{B} (MeV)");
            histos_kinematic_rs_data[i]->GetYaxis()->SetTitle("BDT1");
            histos_kinematic_rs_data[i]->SetTitle("Kinematic MVA (RS data) 2016-2018");
        }
    }

    TCanvas c_rs_data_isolation;
    c_rs_data_isolation.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_isolation_rs_data[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_isolation_rs_data[N-1-i]->DrawNormalized("same");
        }
        leg_isolation_rs_data->AddEntry(Form("h_iso_%i",N-1-i), Form("BDT1 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_isolation_rs_data->Draw("same");
    c_rs_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_isolation_bdt_cut_rs_data_MC.pdf");

    TCanvas c_rs_data_kinematic;
    c_rs_data_kinematic.cd();
    for(int i = 0; i < N; i++)
    {
        if(i == 0)
        {
            histos_kinematic_rs_data[N-1-i]->DrawNormalized();
        }
        else
        {
            histos_kinematic_rs_data[N-1-i]->DrawNormalized("same");
        }
        leg_kinematic_rs_data->AddEntry(Form("h_kin_%i",N-1-i), Form("BDT2 > %.2f",bdt_cuts[N-1-i]), "lp");
    }
    leg_kinematic_rs_data->Draw("same");
    c_rs_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/Bmass_kinematic_bdt_cut_rs_data_MC.pdf");

    TH2D* h2D_isolation_rs_data = new TH2D("h2d_iso", "h2d_iso", 100, 4000, 8000, 100, 0, 1);
    TH2D* h2D_kinematic_rs_data = new TH2D("h2d_kin", "h2d_kin", 100, 4000, 8000, 100, 0, 1);

    t_rs_data_2016->Draw("BDT1 : df_Bp_M >> h2d_iso", "df_status==0");
    t_rs_data_2016->Draw("BDT2 : df_Bp_M >> h2d_kin", "df_status==0");

    TCanvas c1_rs_data_isolation;
    c1_rs_data_isolation.cd();
    c1_rs_data_isolation.SetLogz();
    h2D_isolation_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_isolation_rs_data->GetYaxis()->SetTitle("BDT1");
    h2D_isolation_rs_data->SetTitle(Form("Isolation MVA (3#pi3#pi MC) 2016-2018 (corr = %.3f)",h2D_isolation_rs_data->GetCorrelationFactor()));
    h2D_isolation_rs_data->Draw("COLZ");
    c1_rs_data_isolation.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/BDT1_vs_Bmass_rs_data_MC.pdf");

    TCanvas c1_rs_data_kinematic;
    c1_rs_data_kinematic.cd();
    c1_rs_data_kinematic.SetLogz();
    h2D_kinematic_rs_data->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2D_kinematic_rs_data->GetYaxis()->SetTitle("BDT1");
    h2D_kinematic_rs_data->SetTitle(Form("Kinematic MVA (RS data) 2016-2018 (corr = %.3lf)",h2D_kinematic_rs_data->GetCorrelationFactor()));
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
    h_profile_kinematic_rs_data->GetYaxis()->SetTitle("Average BDT1 per m_{B} bin");
    h_profile_kinematic_rs_data->SetTitle(Form("Kinematic MVA (RS data) 2016-2018 (corr = %.3f)",h2D_kinematic_rs_data->GetCorrelationFactor()));
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
    g_kinematic_rs_data->SetTitle("Kinematic MVA (RS data) 2016-2018 ; BDT2 > x; X-position of B^{+} mass peak (MeV)");
    g_kinematic_rs_data->SetMarkerStyle(8);
    g_kinematic_rs_data->GetXaxis()->SetLimits(-0.1,1.1);
    g_kinematic_rs_data->Draw("AP");
    c3_rs_data_kinematic.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_correlation/mass_peak_pos_vs_bdt2_rs_data_MC.pdf");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}
