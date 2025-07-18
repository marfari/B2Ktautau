void draw_mass_plot(Int_t entry, TTree* t1, TTree* t2, TTree* t3, TString name, TCut mass_vetoe_cuts, Int_t Npar, Int_t isWrong);
Int_t n_choose_k(Int_t n, Int_t k);

void mass_vetoes_plots()
{
    // Input files
    // 3pi3pi MC
    TFile* f1 = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_1/post_sel_tree_bdt1_0.6_bdt2_0.8.root");
    TTree* t1 = (TTree*)f1->Get("DecayTree");  

    cout << "3pi3pi MC" << endl;
    cout << t1->GetEntries() << endl;

    // RS data
    TFile* f2 = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt1_0.6_bdt2_0.8.root");
    TTree* t2 = (TTree*)f2->Get("DecayTree");  

    cout << "RS data" << endl;
    cout << t2->GetEntries() << endl;

    // WS data
    TFile* f3 = new TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0.6_bdt2_0.8.root");
    TTree* t3 = (TTree*)f3->Get("DecayTree");      

    cout << "WS data" << endl;
    cout << t3->GetEntries() << endl;

    Int_t N2 = n_choose_k(7,2);
    Int_t N3 = n_choose_k(7,3);
    Int_t N4 = n_choose_k(7,4);
    Int_t N5 = n_choose_k(7,5);
    Int_t N6 = n_choose_k(7,6);

    // Invariant mass plots (correct mass hypothesis) 
    TCanvas* c_2[N2], c_3[N3], c_4[N4], c_5[N5], c_6[N6];
    TH1D* h_2_mc[N2], h_3_mc[N3], h_4_mc[N4], h_5_mc[N5], h_6_mc[N6];
    TH1D* h_2_rs_data[N2], h_3_rs_data[N3], h_4_rs_data[N4], h_5_rs_data[N5], h_6_rs_data[N6];
    TH1D* h_2_ws_data[N2], h_3_ws_data[N3], h_4_ws_data[N4], h_5_ws_data[N5], h_6_ws_data[N6];

    Int_t a2 = 0;
    Int_t a3 = 0;
    Int_t a4 = 0;
    Int_t a5 = 0;
    Int_t a6 = 0;

    Int_t b2 = 0;
    Int_t b3 = 0;
    Int_t b4 = 0;
    Int_t b5 = 0;
    Int_t b6 = 0;

    TString costheta_pi2 = "(Kp_PX*taup_pi2_PX + Kp_PY*taup_pi2_PY + Kp_PZ*taup_pi2_PZ)/(TMath::Sqrt(pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2))*TMath::Sqrt(pow(taup_pi2_PX,2) + pow(taup_pi2_PY,2) + pow(taup_pi2_PZ,2)))";
    // TString costheta_pi4 = "TMath::Acos(  (Kp_PX*taum_pi1_PX + Kp_PY*taum_pi1_PY + Kp_PZ*taum_pi1_PZ)/(TMath::Sqrt(pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2))*TMath::Sqrt(pow(taum_pi1_PX,2) + pow(taum_pi1_PY,2) + pow(taum_pi1_PZ,2)))  )";
    // TString costheta_pi6 = "TMath::Acos(  (Kp_PX*taum_pi3_PX + Kp_PY*taum_pi3_PY + Kp_PZ*taum_pi3_PZ)/(TMath::Sqrt(pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2))*TMath::Sqrt(pow(taum_pi3_PX,2) + pow(taum_pi3_PY,2) + pow(taum_pi3_PZ,2)))  )";

    gStyle->SetOptStat(0);

    TString Cx_taup = "(df_DV1y-df_BVy)*df_Kp_PZ - (df_DV1z-df_BVz)*df_Kp_PY";
    TString Cy_taup = "(df_DV1z-df_BVz)*df_Kp_PX - (df_DV1x-df_BVx)*df_Kp_PZ";
    TString Cz_taup = "(df_DV1x-df_BVx)*df_Kp_PY - (df_DV1y-df_BVy)*df_Kp_PX";
    TString C_taup = "TMath::Sqrt( pow("+Cx_taup+",2) + pow("+Cy_taup+",2) + pow("+Cz_taup+",2) )";
    TString IP_taup = "2*"+C_taup+"/(TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) + pow(df_Kp_PZ,2) ))";
    
    TString Cx_taum = "(df_DV2y-df_BVy)*df_Kp_PZ - (df_DV2z-df_BVz)*df_Kp_PY";
    TString Cy_taum = "(df_DV2z-df_BVz)*df_Kp_PX - (df_DV2x-df_BVx)*df_Kp_PZ";
    TString Cz_taum = "(df_DV2x-df_BVx)*df_Kp_PY - (df_DV2y-df_BVy)*df_Kp_PX";
    TString C_taum = "TMath::Sqrt( pow("+Cx_taum+",2) + pow("+Cy_taum+",2) + pow("+Cz_taum+",2) )";
    TString IP_taum = "2*"+C_taum+"/(TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) + pow(df_Kp_PZ,2) ))";

    TString IP_taup_cut = "("+IP_taup+" > 0.25) || (TMath::Abs(Bp_M02-892) > 50)";
    TString IP_taum_cut = "("+IP_taum+" > 0.25) || ((TMath::Abs(Bp_M04-892) > 50) && (TMath::Abs(Bp_M06-892) > 50))";
    
    TString IP_2_particles_cut = "("+IP_taup_cut+") &&  ("+IP_taum_cut+")";
    TString IP_4_particles_cut = "( (("+IP_taup+" > 0.15) && ("+IP_taum+" > 0.15)) || ((TMath::Abs(Bp_M0124-1864.84) > 30) && (TMath::Abs(Bp_M0126-1864.84) > 30) && (TMath::Abs(Bp_M0146-1864.84) > 30) && (TMath::Abs(Bp_M0234-1864.84) > 30) && (TMath::Abs(Bp_M0236-1864.84) > 30) && (TMath::Abs(Bp_M0245-1864.84) > 30) && (TMath::Abs(Bp_M0256-1864.84) > 30)) )";

    TCut mass_vetoe_cuts = "";
    mass_vetoe_cuts += "(TMath::Abs(Bp_M02-1864.84) > 30) && (TMath::Abs(Bp_M04-1864.84) > 30) && (TMath::Abs(Bp_M06-1864.84) > 30)"; // 2 particles: D0
    mass_vetoe_cuts += "(TMath::Abs(Bp_M046-1869.66) > 20)"; // 3 particles: D+
    mass_vetoe_cuts += "(TMath::Abs(Bp_M0456-1864.84) > 30)"; // 4 particles D0
    mass_vetoe_cuts += "(TMath::Abs(Bp_M01246-2010.26) > 20) && (TMath::Abs(Bp_M02346-2010.26) > 20) && (TMath::Abs(Bp_M02456-2010.26) > 20)"; // 5 particles (D*-)
    mass_vetoe_cuts += IP_2_particles_cut;
    mass_vetoe_cuts += IP_4_particles_cut;

    // TCanvas c;
    // c.cd();
    // TH1D* h_sig = new TH1D("h_sig", "h_sig", 50, 0.99, 1);
    // TH1D* h_rs = new TH1D("h_rs", "h_rs", 50, 0.99, 1);
    // t1->Draw(costheta_pi2+" >> h_sig");
    // t2->Draw(costheta_pi2+" >> h_rs");
    // h_sig->GetXaxis()->SetTitle("cos(#theta_{K,#pi^{2}})");
    // h_sig->GetYaxis()->SetTitle("Normalized entries / (0.04)");
    // h_sig->SetTitle("Events in K*0 mass region (within 50 MeV)");
    // h_rs->SetLineColor(kBlack);
    // h_sig->SetLineColor(kBlue);
    // h_sig->DrawNormalized();
    // h_rs->DrawNormalized("same");
    // c.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/costheta2.pdf");

    // TH1D* pkpi_02_mc = new TH1D("pkpi_02_mc", "pkpi_02_mc", 100, -10, 10);
    // TH1D* pkpi_04_mc = new TH1D("pkpi_04_mc", "pkpi_04_mc", 100, -10, 10);
    // TH1D* pkpi_06_mc = new TH1D("pkpi_06_mc", "pkpi_06_mc", 100, -10, 10);
    // TH1D* pkpi_02_rs = new TH1D("pkpi_02_rs", "pkpi_02_rs", 100, -10, 10);
    // TH1D* pkpi_04_rs = new TH1D("pkpi_04_rs", "pkpi_04_rs", 100, -10, 10);
    // TH1D* pkpi_06_rs = new TH1D("pkpi_06_rs", "pkpi_06_rs", 100, -10, 10);

    // t1->Draw("(Kp_PX + taup_pi2_PX)/(Kp_PY + taup_pi2_PY) - (df_DV1x - df_BVx)/(df_DV1y - df_BVy) >> pkpi_02_mc");
    // t1->Draw("(Kp_PX + taum_pi1_PX)/(Kp_PY + taum_pi1_PY) - (df_DV2x - df_BVx)/(df_DV2y - df_BVy) >> pkpi_04_mc");
    // t1->Draw("(Kp_PX + taum_pi3_PX)/(Kp_PY + taum_pi3_PY) - (df_DV2x - df_BVx)/(df_DV2y - df_BVy) >> pkpi_06_mc");

    // t2->Draw("(Kp_PX + taup_pi2_PX)/(Kp_PY + taup_pi2_PY) - (df_DV1x - df_BVx)/(df_DV1y - df_BVy) >> pkpi_02_rs");
    // t2->Draw("(Kp_PX + taum_pi1_PX)/(Kp_PY + taum_pi1_PY) - (df_DV2x - df_BVx)/(df_DV2y - df_BVy) >> pkpi_04_rs");
    // t2->Draw("(Kp_PX + taum_pi3_PX)/(Kp_PY + taum_pi3_PY) - (df_DV2x - df_BVx)/(df_DV2y - df_BVy) >> pkpi_06_rs");

    // TCanvas B2;
    // B2.cd();
    // pkpi_02_mc->GetXaxis()->SetTitle("#vec{p}_{Kpi} parallel to DV-BV constraint");
    // pkpi_02_mc->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // pkpi_02_mc->SetTitle("02: x/y");
    // pkpi_02_rs->SetLineColor(kBlack);
    // pkpi_02_mc->SetLineColor(kBlue);
    // pkpi_02_mc->DrawNormalized();
    // pkpi_02_rs->DrawNormalized("same");
    // B2.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pkpi_02.pdf");

    // TCanvas B4;
    // B4.cd();
    // pkpi_04_rs->GetXaxis()->SetTitle("#vec{p}_{Kpi} parallel to DV-BV constraint");
    // pkpi_04_rs->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // pkpi_04_rs->SetTitle("04: x/y");
    // pkpi_04_rs->SetLineColor(kBlack);
    // pkpi_04_mc->SetLineColor(kBlue);
    // pkpi_04_rs->DrawNormalized();
    // pkpi_04_mc->DrawNormalized("same");
    // B4.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pkpi_04.pdf");

    // TCanvas B6;
    // B6.cd();
    // pkpi_06_rs->GetXaxis()->SetTitle("#vec{p}_{Kpi} parallel to DV-BV constraint");
    // pkpi_06_rs->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // pkpi_06_rs->SetTitle("06: x/y");
    // pkpi_06_rs->SetLineColor(kBlack);
    // pkpi_06_mc->SetLineColor(kBlue);
    // pkpi_06_rs->DrawNormalized();
    // pkpi_06_mc->DrawNormalized("same");
    // B6.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pkpi_06.pdf");

    // TH1D* pkpi_02_mc_1 = new TH1D("pkpi_02_mc_1", "pkpi_02_mc_1", 100, -50, 50);
    // TH1D* pkpi_04_mc_1 = new TH1D("pkpi_04_mc_1", "pkpi_04_mc_1", 100, -50, 50);
    // TH1D* pkpi_06_mc_1 = new TH1D("pkpi_06_mc_1", "pkpi_06_mc_1", 100, -50, 50);
    // TH1D* pkpi_02_rs_1 = new TH1D("pkpi_02_rs_1", "pkpi_02_rs_1", 100, -50, 50);
    // TH1D* pkpi_04_rs_1 = new TH1D("pkpi_04_rs_1", "pkpi_04_rs_1", 100, -50, 50);
    // TH1D* pkpi_06_rs_1 = new TH1D("pkpi_06_rs_1", "pkpi_06_rs_1", 100, -50, 50);

    // t1->Draw("(Kp_PZ + taup_pi2_PZ)/(Kp_PY + taup_pi2_PY) - (df_DV1z - df_BVz)/(df_DV1y - df_BVy) >> pkpi_02_mc_1");
    // t1->Draw("(Kp_PZ + taum_pi1_PZ)/(Kp_PY + taum_pi1_PY) - (df_DV2z - df_BVz)/(df_DV2y - df_BVy) >> pkpi_04_mc_1");
    // t1->Draw("(Kp_PZ + taum_pi3_PZ)/(Kp_PY + taum_pi3_PY) - (df_DV2z - df_BVz)/(df_DV2y - df_BVy) >> pkpi_06_mc_1");

    // t2->Draw("(Kp_PZ + taup_pi2_PZ)/(Kp_PY + taup_pi2_PY) - (df_DV1z - df_BVz)/(df_DV1y - df_BVy) >> pkpi_02_rs_1");
    // t2->Draw("(Kp_PZ + taum_pi1_PZ)/(Kp_PY + taum_pi1_PY) - (df_DV2z - df_BVz)/(df_DV2y - df_BVy) >> pkpi_04_rs_1");
    // t2->Draw("(Kp_PZ + taum_pi3_PZ)/(Kp_PY + taum_pi3_PY) - (df_DV2z - df_BVz)/(df_DV2y - df_BVy) >> pkpi_06_rs_1");

    // TCanvas B1;
    // B1.cd();
    // pkpi_02_mc_1->GetXaxis()->SetTitle("#vec{p}_{Kpi} parallel to DV-BV constraint");
    // pkpi_02_mc_1->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // pkpi_02_mc_1->SetTitle("02: z/y");
    // pkpi_02_rs_1->SetLineColor(kBlack);
    // pkpi_02_mc_1->SetLineColor(kBlue);
    // pkpi_02_mc_1->DrawNormalized();
    // pkpi_02_rs_1->DrawNormalized("same");
    // B1.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pkpi_02_1.pdf");

    // TCanvas B3;
    // B3.cd();
    // pkpi_04_rs_1->GetXaxis()->SetTitle("#vec{p}_{Kpi} parallel to DV-BV constraint");
    // pkpi_04_rs_1->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // pkpi_04_rs_1->SetTitle("04: z/y");
    // pkpi_04_rs_1->SetLineColor(kBlack);
    // pkpi_04_mc_1->SetLineColor(kBlue);
    // pkpi_04_rs_1->DrawNormalized();
    // pkpi_04_mc_1->DrawNormalized("same");
    // B3.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pkpi_04_1.pdf");

    // TCanvas B5;
    // B5.cd();
    // pkpi_06_rs_1->GetXaxis()->SetTitle("#vec{p}_{Kpi} parallel to DV-BV constraint");
    // pkpi_06_rs_1->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // pkpi_06_rs_1->SetTitle("06: z/y");
    // pkpi_06_rs_1->SetLineColor(kBlack);
    // pkpi_06_mc_1->SetLineColor(kBlue);
    // pkpi_06_rs_1->DrawNormalized();
    // pkpi_06_mc_1->DrawNormalized("same");
    // B5.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pkpi_06_1.pdf");

    // TH1D* p_Ktaup_mc = new TH1D("p_Ktaup_mc", "p_Ktaup_mc", 100, -50, 50);
    // TH1D* p_Ktaup_rs = new TH1D("p_Ktaup_rs", "p_Ktaup_rs", 100, -50, 50);
    // TH1D* p_Ktaum_mc = new TH1D("p_Ktaum_mc", "p_Ktaum_mc", 100, -50, 50);
    // TH1D* p_Ktaum_rs = new TH1D("p_Ktaum_rs", "p_Ktaum_rs", 100, -50, 50);

    // t1->Draw("(Kp_PZ + df_taup_PZ)/(Kp_PY + df_taup_PY) - (df_DV1z - df_BVz)/(df_DV1y - df_BVy) >> p_Ktaup_mc");
    // t1->Draw("(Kp_PZ + df_taum_PZ)/(Kp_PY + df_taum_PY) - (df_DV2z - df_BVz)/(df_DV2y - df_BVy) >> p_Ktaum_mc");
    // t2->Draw("(Kp_PZ + df_taup_PZ)/(Kp_PY + df_taup_PY) - (df_DV1z - df_BVz)/(df_DV1y - df_BVy) >> p_Ktaup_rs");
    // t2->Draw("(Kp_PZ + df_taum_PZ)/(Kp_PY + df_taum_PY) - (df_DV2z - df_BVz)/(df_DV2y - df_BVy) >> p_Ktaum_rs");

    // TCanvas C1;
    // C1.cd();
    // p_Ktaup_mc->GetXaxis()->SetTitle("#vec{p}_{Ktau+} parallel to DV-BV constraint");
    // p_Ktaup_mc->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // p_Ktaup_mc->SetTitle("z/y");
    // p_Ktaup_rs->SetLineColor(kBlack);
    // p_Ktaup_mc->SetLineColor(kBlue);
    // p_Ktaup_mc->DrawNormalized();
    // p_Ktaup_rs->DrawNormalized("same");
    // C1.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pktaup.pdf");

    // TCanvas C2;
    // C2.cd();
    // p_Ktaum_mc->GetXaxis()->SetTitle("#vec{p}_{Ktau-} parallel to DV-BV constraint");
    // p_Ktaum_mc->GetYaxis()->SetTitle("Normalized entries / (0.2)");
    // p_Ktaum_mc->SetTitle("z/y");
    // p_Ktaum_rs->SetLineColor(kBlack);
    // p_Ktaum_mc->SetLineColor(kBlue);
    // p_Ktaum_mc->DrawNormalized();
    // p_Ktaum_rs->DrawNormalized("same");
    // C2.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/pktaum.pdf");

    // TH1D* BV_DV1_distance_mc = new TH1D("BV_DV1_distance_mc", "BV_DV1_distance_mc", 100, 0, 20);
    // TH1D* BV_DV2_distance_mc = new TH1D("BV_DV2_distance_mc", "BV_DV2_distance_mc", 100, 0, 20);
    // TH1D* BV_DV1_distance_rs = new TH1D("BV_DV1_distance_rs", "BV_DV1_distance_rs", 100, 0, 20);
    // TH1D* BV_DV2_distance_rs = new TH1D("BV_DV2_distance_rs", "BV_DV2_distance_rs", 100, 0, 20);

    // t1->Draw("TMath::Sqrt(pow(df_BVx - df_DV1x,2) + pow(df_BVy-df_DV1y,2) + pow(df_BVz-df_DV1z,2)) >> BV_DV1_distance_mc", "(TMath::Abs(Bp_M02-892) < 40) && (TMath::Abs(Bp_M04-892) < 40) && (TMath::Abs(Bp_M06-892) < 40)");
    // t1->Draw("TMath::Sqrt(pow(df_BVx - df_DV2x,2) + pow(df_BVy-df_DV2y,2) + pow(df_BVz-df_DV2z,2)) >> BV_DV2_distance_mc", "(TMath::Abs(Bp_M02-892) < 40) && (TMath::Abs(Bp_M04-892) < 40) && (TMath::Abs(Bp_M06-892) < 40)");
    // t2->Draw("TMath::Sqrt(pow(df_BVx - df_DV1x,2) + pow(df_BVy-df_DV1y,2) + pow(df_BVz-df_DV1z,2)) >> BV_DV1_distance_rs", "(TMath::Abs(Bp_M02-892) < 40) && (TMath::Abs(Bp_M04-892) < 40) && (TMath::Abs(Bp_M06-892) < 40)");
    // t2->Draw("TMath::Sqrt(pow(df_BVx - df_DV2x,2) + pow(df_BVy-df_DV2y,2) + pow(df_BVz-df_DV2z,2)) >> BV_DV2_distance_rs", "(TMath::Abs(Bp_M02-892) < 40) && (TMath::Abs(Bp_M04-892) < 40) && (TMath::Abs(Bp_M06-892) < 40)");

    // TCanvas D1;
    // D1.cd();
    // BV_DV1_distance_mc->GetXaxis()->SetTitle("BV-DV1 distance (mm)");
    // BV_DV1_distance_mc->GetYaxis()->SetTitle("Normalized entries / (0.02 mm)");
    // BV_DV1_distance_mc->SetTitle("Events in K^{*0} mass region (< 40 MeV)");
    // BV_DV1_distance_rs->SetLineColor(kBlack);
    // BV_DV1_distance_mc->SetLineColor(kBlue);
    // BV_DV1_distance_mc->DrawNormalized();
    // BV_DV1_distance_rs->DrawNormalized("same");
    // D1.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/BV_DV1_distance.pdf");

    // TCanvas D2;
    // D2.cd();
    // BV_DV2_distance_mc->GetXaxis()->SetTitle("BV-DV2 distance (mm)");
    // BV_DV2_distance_mc->GetYaxis()->SetTitle("Normalized entries / (0.02 mm)");
    // BV_DV2_distance_mc->SetTitle("Events in K^{*0} mass region (< 40 MeV)");
    // BV_DV2_distance_rs->SetLineColor(kBlack);
    // BV_DV2_distance_mc->SetLineColor(kBlue);
    // BV_DV2_distance_mc->DrawNormalized();
    // BV_DV2_distance_rs->DrawNormalized("same");
    // D2.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/BV_DV2_distance.pdf");


    // Double_t xmin = 800;
    // Double_t xmax = 1100;
    // Double_t nbins = 60;

    // TCanvas c_phys;
    // c_phys.cd();
    // TH1D* h_phys_1 = new TH1D("h_phys_1", "h_phys_1", nbins, xmin, xmax);
    // TH1D* h_phys_2 = new TH1D("h_phys_2", "h_phys_2", nbins, xmin, xmax);
    // TH1D* h_phys_3 = new TH1D("h_phys_3", "h_phys_3", nbins, xmin, xmax);
    // TH1D* h_phys_4 = new TH1D("h_phys_4", "h_phys_4", nbins, xmin, xmax);
    // t2->Draw("Bp_M02 >> h_phys_1", mass_vetoe_cuts);
    // t2->Draw("Bp_M02 >> h_phys_2", "BDT1 > 0.8"+mass_vetoe_cuts);
    // t2->Draw("Bp_M02 >> h_phys_3", "BDT1 > 0.9"+mass_vetoe_cuts);
    // t2->Draw("Bp_M02 >> h_phys_4", "BDT1 > 0.95"+mass_vetoe_cuts);
    // h_phys_1->SetLineColor(kBlack);
    // h_phys_2->SetLineColor(kBlue);
    // h_phys_3->SetLineColor(kGreen+1);
    // h_phys_4->SetLineColor(kRed);
    // h_phys_4->GetXaxis()->SetTitle("Bp_M02");
    // h_phys_4->GetYaxis()->SetTitle(Form("Normalized entries / (%.2f MeV)", (xmax-xmin)/nbins));
    // h_phys_4->SetTitle("Effect of physics BDT");
    // h_phys_4->DrawNormalized();
    // h_phys_3->DrawNormalized("same");
    // h_phys_2->DrawNormalized("same");
    // h_phys_1->DrawNormalized("same");
    // TLegend* leg_phys = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg_phys->AddEntry(h_phys_1, "BDT1 > 0.6", "lp");
    // leg_phys->AddEntry(h_phys_2, "BDT1 > 0.8", "lp");
    // leg_phys->AddEntry(h_phys_3, "BDT1 > 0.9", "lp");
    // leg_phys->AddEntry(h_phys_4, "BDT1 > 0.95", "lp");
    // leg_phys->Draw("same");
    // c_phys.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/Bp_M02_phys_BDT.pdf");

    // TCanvas c1_phys;
    // c1_phys.cd();
    // TH1D* h1_phys_1 = new TH1D("h1_phys_1", "h1_phys_1", nbins, xmin, xmax);
    // TH1D* h1_phys_2 = new TH1D("h1_phys_2", "h1_phys_2", nbins, xmin, xmax);
    // TH1D* h1_phys_3 = new TH1D("h1_phys_3", "h1_phys_3", nbins, xmin, xmax);
    // TH1D* h1_phys_4 = new TH1D("h1_phys_4", "h1_phys_4", nbins, xmin, xmax);
    // t2->Draw("Bp_M04 >> h1_phys_1", mass_vetoe_cuts);
    // t2->Draw("Bp_M04 >> h1_phys_2", "BDT1 > 0.8"+mass_vetoe_cuts);
    // t2->Draw("Bp_M04 >> h1_phys_3", "BDT1 > 0.9"+mass_vetoe_cuts);
    // t2->Draw("Bp_M04 >> h1_phys_4", "BDT1 > 0.95"+mass_vetoe_cuts);
    // h1_phys_1->SetLineColor(kBlack);
    // h1_phys_2->SetLineColor(kBlue);
    // h1_phys_3->SetLineColor(kGreen+1);
    // h1_phys_4->SetLineColor(kRed);
    // h1_phys_1->GetXaxis()->SetTitle("Bp_M04");
    // h1_phys_1->GetYaxis()->SetTitle(Form("Normalized entries / (%.2f MeV)", (xmax-xmin)/nbins));
    // h1_phys_1->SetTitle("Effect of physics BDT");
    // h1_phys_1->DrawNormalized();
    // h1_phys_3->DrawNormalized("same");
    // h1_phys_2->DrawNormalized("same");
    // h1_phys_4->DrawNormalized("same");
    // TLegend* leg1_phys = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg1_phys->AddEntry(h1_phys_1, "BDT1 > 0.6", "lp");
    // leg1_phys->AddEntry(h1_phys_2, "BDT1 > 0.8", "lp");
    // leg1_phys->AddEntry(h1_phys_3, "BDT1 > 0.9", "lp");
    // leg1_phys->AddEntry(h1_phys_4, "BDT1 > 0.95", "lp");
    // leg1_phys->Draw("same");
    // c1_phys.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/Bp_M04_phys_BDT.pdf");

    // TCanvas c2_phys;
    // c2_phys.cd();
    // TH1D* h2_phys_1 = new TH1D("h2_phys_1", "h1_phys_1", nbins, xmin, xmax);
    // TH1D* h2_phys_2 = new TH1D("h2_phys_2", "h1_phys_2", nbins, xmin, xmax);
    // TH1D* h2_phys_3 = new TH1D("h2_phys_3", "h1_phys_3", nbins, xmin, xmax);
    // TH1D* h2_phys_4 = new TH1D("h2_phys_4", "h1_phys_4", nbins, xmin, xmax);
    // t2->Draw("Bp_M06 >> h2_phys_1", mass_vetoe_cuts);
    // t2->Draw("Bp_M06 >> h2_phys_2", "BDT1 > 0.8"+mass_vetoe_cuts);
    // t2->Draw("Bp_M06 >> h2_phys_3", "BDT1 > 0.9"+mass_vetoe_cuts);
    // t2->Draw("Bp_M06 >> h2_phys_4", "BDT1 > 0.95"+mass_vetoe_cuts);
    // h2_phys_1->SetLineColor(kBlack);
    // h2_phys_2->SetLineColor(kBlue);
    // h2_phys_3->SetLineColor(kGreen+1);
    // h2_phys_4->SetLineColor(kRed);
    // h2_phys_4->GetXaxis()->SetTitle("Bp_M06");
    // h2_phys_4->GetYaxis()->SetTitle(Form("Normalized entries / (%.2f MeV)", (xmax-xmin)/nbins));
    // h2_phys_4->SetTitle("Effect of physics BDT");
    // h2_phys_4->DrawNormalized();
    // h2_phys_3->DrawNormalized("same");
    // h2_phys_2->DrawNormalized("same");
    // h2_phys_1->DrawNormalized("same");
    // TLegend* leg2_phys = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg2_phys->AddEntry(h2_phys_1, "BDT1 > 0.6", "lp");
    // leg2_phys->AddEntry(h2_phys_2, "BDT1 > 0.8", "lp");
    // leg2_phys->AddEntry(h2_phys_3, "BDT1 > 0.9", "lp");
    // leg2_phys->AddEntry(h2_phys_4, "BDT1 > 0.95", "lp");
    // leg2_phys->Draw("same");
    // c2_phys.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/Bp_M06_phys_BDT.pdf");

    // TCanvas c_comb;
    // c_comb.cd();
    // TH1D* h_comb_1 = new TH1D("h_comb_1", "h_comb_1", nbins, xmin, xmax);
    // TH1D* h_comb_2 = new TH1D("h_comb_2", "h_comb_2", nbins, xmin, xmax);
    // TH1D* h_comb_3 = new TH1D("h_comb_3", "h_comb_3", nbins, xmin, xmax);
    // TH1D* h_comb_4 = new TH1D("h_comb_4", "h_comb_4", nbins, xmin, xmax);
    // t2->Draw("Bp_M02 >> h_comb_1", mass_vetoe_cuts);
    // t2->Draw("Bp_M02 >> h_comb_2", "BDT2 > 0.9"+mass_vetoe_cuts);
    // t2->Draw("Bp_M02 >> h_comb_3", "BDT2 > 0.95"+mass_vetoe_cuts);
    // t2->Draw("Bp_M02 >> h_comb_4", "BDT2 > 0.97"+mass_vetoe_cuts);
    // h_comb_1->SetLineColor(kBlack);
    // h_comb_2->SetLineColor(kBlue);
    // h_comb_3->SetLineColor(kGreen+1);
    // h_comb_4->SetLineColor(kRed);
    // h_comb_3->GetXaxis()->SetTitle("Bp_M02");
    // h_comb_3->GetYaxis()->SetTitle(Form("Normalized entries / (%.2f MeV)", (xmax-xmin)/nbins));
    // h_comb_3->SetTitle("Effect of combinatorial BDT");
    // h_comb_3->DrawNormalized();
    // h_comb_1->DrawNormalized("same");
    // h_comb_2->DrawNormalized("same");
    // h_comb_4->DrawNormalized("same");
    // TLegend* leg_comb = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg_comb->AddEntry(h_comb_1, "BDT2 > 0.8", "lp");
    // leg_comb->AddEntry(h_comb_2, "BDT2 > 0.9", "lp");
    // leg_comb->AddEntry(h_comb_3, "BDT2 > 0.95", "lp");
    // leg_comb->AddEntry(h_comb_4, "BDT2 > 0.97", "lp");
    // leg_comb->Draw("same");
    // c_comb.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/Bp_M02_comb_BDT.pdf");

    // TCanvas c1_comb;
    // c1_comb.cd();
    // TH1D* h1_comb_1 = new TH1D("h1_comb_1", "h1_comb_1", nbins, xmin, xmax);
    // TH1D* h1_comb_2 = new TH1D("h1_comb_2", "h1_comb_2", nbins, xmin, xmax);
    // TH1D* h1_comb_3 = new TH1D("h1_comb_3", "h1_comb_3", nbins, xmin, xmax);
    // TH1D* h1_comb_4 = new TH1D("h1_comb_4", "h1_comb_4", nbins, xmin, xmax);
    // t2->Draw("Bp_M04 >> h1_comb_1", mass_vetoe_cuts);
    // t2->Draw("Bp_M04 >> h1_comb_2", "BDT2 > 0.9"+mass_vetoe_cuts);
    // t2->Draw("Bp_M04 >> h1_comb_3", "BDT2 > 0.95"+mass_vetoe_cuts);
    // t2->Draw("Bp_M04 >> h1_comb_4", "BDT2 > 0.97"+mass_vetoe_cuts);
    // h1_comb_1->SetLineColor(kBlack);
    // h1_comb_2->SetLineColor(kBlue);
    // h1_comb_3->SetLineColor(kGreen+1);
    // h1_comb_4->SetLineColor(kRed);
    // h1_comb_1->GetXaxis()->SetTitle("Bp_M04");
    // h1_comb_1->GetYaxis()->SetTitle(Form("Normalized entries / (%.2f MeV)", (xmax-xmin)/nbins));
    // h1_comb_1->SetTitle("Effect of combinatorial BDT");
    // h1_comb_1->DrawNormalized();
    // h1_comb_3->DrawNormalized("same");
    // h1_comb_2->DrawNormalized("same");
    // h1_comb_4->DrawNormalized("same");
    // TLegend* leg1_comb = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg1_comb->AddEntry(h1_comb_1, "BDT2 > 0.8", "lp");
    // leg1_comb->AddEntry(h1_comb_2, "BDT2 > 0.9", "lp");
    // leg1_comb->AddEntry(h1_comb_3, "BDT2 > 0.95", "lp");
    // leg1_comb->AddEntry(h1_comb_4, "BDT2 > 0.97", "lp");
    // leg1_comb->Draw("same");
    // c1_comb.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/Bp_M04_comb_BDT.pdf");

    // TCanvas c2_comb;
    // c2_comb.cd();
    // TH1D* h2_comb_1 = new TH1D("h2_comb_1", "h1_comb_1", nbins, xmin, xmax);
    // TH1D* h2_comb_2 = new TH1D("h2_comb_2", "h1_comb_2", nbins, xmin, xmax);
    // TH1D* h2_comb_3 = new TH1D("h2_comb_3", "h1_comb_3", nbins, xmin, xmax);
    // TH1D* h2_comb_4 = new TH1D("h2_comb_4", "h1_comb_4", nbins, xmin, xmax);
    // t2->Draw("Bp_M06 >> h2_comb_1", mass_vetoe_cuts);
    // t2->Draw("Bp_M06 >> h2_comb_2", "BDT2 > 0.9"+mass_vetoe_cuts);
    // t2->Draw("Bp_M06 >> h2_comb_3", "BDT2 > 0.95"+mass_vetoe_cuts);
    // t2->Draw("Bp_M06 >> h2_comb_4", "BDT2 > 0.97"+mass_vetoe_cuts);
    // h2_comb_1->SetLineColor(kBlack);
    // h2_comb_2->SetLineColor(kBlue);
    // h2_comb_3->SetLineColor(kGreen+1);
    // h2_comb_4->SetLineColor(kRed);
    // h2_comb_4->GetXaxis()->SetTitle("Bp_M06");
    // h2_comb_4->GetYaxis()->SetTitle(Form("Normalized entries / (%.2f MeV)", (xmax-xmin)/nbins));
    // h2_comb_4->SetTitle("Effect of combinatorial BDT");
    // h2_comb_4->DrawNormalized();
    // h2_comb_3->DrawNormalized("same");
    // h2_comb_2->DrawNormalized("same");
    // h2_comb_1->DrawNormalized("same");
    // TLegend* leg2_comb = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg2_comb->AddEntry(h2_comb_1, "BDT2 > 0.8", "lp");
    // leg2_comb->AddEntry(h2_comb_2, "BDT2 > 0.9", "lp");
    // leg2_comb->AddEntry(h2_comb_3, "BDT2 > 0.95", "lp");
    // leg2_comb->AddEntry(h2_comb_4, "BDT2 > 0.97", "lp");
    // leg2_comb->Draw("same");
    // c2_comb.SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/Bp_M06_comb_BDT.pdf");

    Double_t N_sig_num = t1->GetEntries(mass_vetoe_cuts);
    Double_t N_rs_num = t2->GetEntries(mass_vetoe_cuts);
    Double_t N_ws_num = t3->GetEntries(mass_vetoe_cuts);

    cout << "Signal eff = " << N_sig_num/(t1->GetEntries())*100 << endl;
    cout << "RS data eff = " << N_rs_num/(t2->GetEntries())*100 << endl;
    cout << "WS data eff = " << N_ws_num/(t3->GetEntries())*100 << endl;

    for(int i1 = 0; i1 < 7; i1++)
    {
        for(int i2 = 1; i2 < 7; i2++)
        {   
            if(i2 > i1)
            {
                TString name2 = Form("Bp_M%i%i",i1,i2);
                draw_mass_plot(a2, t1, t2, t3, name2, mass_vetoe_cuts, 2, 0);
                TString name2_massH1 = Form("Bp_M%i%i_massH1",i1,i2);
                draw_mass_plot(a2, t1, t2, t3, name2_massH1, mass_vetoe_cuts, 2, 1);
                if(i1 == 0)
                {
                    TString name2_massH2 = Form("Bp_M%i%i_massH2",i1,i2);
                    draw_mass_plot(b2, t1, t2, t3, name2_massH2, mass_vetoe_cuts, 2, 2);
                    b2 += 1;
                }
                for(int i3 = 1; i3 < 7; i3++)
                {
                    if(i3 > i2)
                    {
                        TString name3 = Form("Bp_M%i%i%i",i1,i2,i3);
                        draw_mass_plot(a3, t1, t2, t3, name3, mass_vetoe_cuts, 3, 0);
                        TString name3_massH1 = Form("Bp_M%i%i%i_massH1",i1,i2,i3);
                        draw_mass_plot(a3, t1, t2, t3, name3_massH1, mass_vetoe_cuts, 3, 1);
                        if(i1 == 0)
                        {
                            TString name3_massH2 = Form("Bp_M%i%i%i_massH2",i1,i2,i3);
                            draw_mass_plot(b3, t1, t2, t3, name3_massH2, mass_vetoe_cuts, 3, 2);
                            b3 += 1;
                        }
                        for(int i4 = 1; i4 < 7; i4++)
                        {
                            if(i4 > i3)
                            {
                                TString name4 = Form("Bp_M%i%i%i%i",i1,i2,i3,i4);
                                draw_mass_plot(a4, t1, t2, t3, name4, mass_vetoe_cuts, 4, 0);
                                TString name4_massH1 = Form("Bp_M%i%i%i%i_massH1",i1,i2,i3,i4);
                                draw_mass_plot(a4, t1, t2, t3, name4_massH1, mass_vetoe_cuts, 4, 1);
                                if(i1 == 0)
                                {
                                    TString name4_massH2 = Form("Bp_M%i%i%i%i_massH2",i1,i2,i3,i4);
                                    draw_mass_plot(b4, t1, t2, t3, name4_massH2, mass_vetoe_cuts, 4, 2);
                                    b4 += 1;
                                }
                                for(int i5 = 1; i5 < 7; i5++)
                                {
                                    if(i5 > i4)
                                    {
                                        TString name5 = Form("Bp_M%i%i%i%i%i",i1,i2,i3,i4,i5);
                                        draw_mass_plot(a5, t1, t2, t3, name5, mass_vetoe_cuts, 5, 0);
                                        TString name5_massH1 = Form("Bp_M%i%i%i%i%i_massH1",i1,i2,i3,i4,i5);
                                        draw_mass_plot(a5, t1, t2, t3, name5_massH1, mass_vetoe_cuts, 5, 1);
                                        if(i1 ==  0)
                                        {
                                            TString name5_massH2 = Form("Bp_M%i%i%i%i%i_massH2",i1,i2,i3,i4,i5);
                                            draw_mass_plot(b5, t1, t2, t3, name5_massH2, mass_vetoe_cuts, 5, 2);
                                            b5 += 1;

                                        }
                                        for(int i6 = 1; i6 < 7; i6++)
                                        {   
                                            if(i6 > i5)
                                            {
                                                TString name6 = Form("Bp_M%i%i%i%i%i%i",i1,i2,i3,i4,i5,i6);
                                                draw_mass_plot(a6, t1, t2, t3, name6, mass_vetoe_cuts, 6, 0);
                                                TString name6_massH1 = Form("Bp_M%i%i%i%i%i%i_massH1",i1,i2,i3,i4,i5,i6);
                                                draw_mass_plot(a6, t1, t2, t3, name6_massH1, mass_vetoe_cuts, 6, 1);
                                                if(i1 == 0)
                                                {
                                                    TString name6_massH2 = Form("Bp_M%i%i%i%i%i%i_massH2",i1,i2,i3,i4,i5,i6);
                                                    draw_mass_plot(b6, t1, t2, t3, name6_massH2, mass_vetoe_cuts, 6, 2);
                                                    b6 += 1;
                                                }
                                                a6 += 1;
                                            }
                                        }
                                        a5 += 1;
                                    }
                                }
                                a4 +=1;
                            }
                        }   
                        a3 += 1;
                    }
                }        
                a2 += 1;
            }
        }
    }
    
}

void draw_mass_plot(Int_t entry, TTree* t1, TTree* t2, TTree* t3, TString name, TCut mass_vetoe_cuts, Int_t Npar, Int_t isWrong)
{
    gStyle->SetOptStat(0);
    // MC vs RS data vs WS data
    TCanvas c;
    c.cd();

    Int_t nbins = 80;

    std::vector<Double_t> x_min_vec;
    x_min_vec.push_back(t1->GetMinimum(name));
    x_min_vec.push_back(t2->GetMinimum(name));
    x_min_vec.push_back(t3->GetMinimum(name));

    std::vector<Double_t> x_max_vec;
    x_max_vec.push_back(t1->GetMaximum(name));
    x_max_vec.push_back(t2->GetMaximum(name));
    x_max_vec.push_back(t3->GetMaximum(name));

    Double_t x_min = x_min_vec[ std::distance( x_min_vec.begin(), std::min_element(x_min_vec.begin(), x_min_vec.end()) ) ];
    Double_t x_max = x_max_vec[ std::distance( x_max_vec.begin(), std::max_element(x_max_vec.begin(), x_max_vec.end()) ) ];

    if(Npar == 2)
    {
        if((name == "Bp_M02") || (name == "Bp_M04") || (name == "Bp_M06"))
        {
            x_min = 500;
            x_max = 3500;
            nbins = 80;
        }
    }
    else if(Npar == 3)
    {
        if(name == "Bp_M046")
        {
            x_min = 500;
            x_max = 3500;
            nbins = 100;
        }
    }
    else if(Npar == 4)
    {
        if(name == "Bp_M0456")
        {
            x_min = 1000;
            x_max = 4000;
            nbins = 50;
        }
    }

    TH1D* h_mc = new TH1D(Form("h_%i_mc_%i_%i",Npar,entry,isWrong), Form("h_%i_mc_%i_%i",Npar,entry,isWrong), nbins, x_min, x_max);
    TH1D* h_rs_data = new TH1D(Form("h_%i_rs_data_%i_%i",Npar,entry,isWrong), Form("h_%i_rs_data_%i_%i",Npar,entry,isWrong), nbins, x_min, x_max);
    TH1D* h_ws_data = new TH1D(Form("h_%i_ws_data_%i_%i",Npar,entry,isWrong), Form("h_%i_ws_data_%i_%i",Npar,entry,isWrong), nbins, x_min, x_max);
    TH1D* h1_mc = new TH1D(Form("h1_%i_mc_%i_%i",Npar,entry,isWrong), Form("h1_%i_mc_%i_%i",Npar,entry,isWrong), nbins, x_min, x_max);
    TH1D* h1_rs_data = new TH1D(Form("h1_%i_rs_data_%i_%i",Npar,entry,isWrong), Form("h1_%i_rs_data_%i_%i",Npar,entry,isWrong), nbins, x_min, x_max);

    t1->Draw(name+Form(" >> h_%i_mc_%i_%i",Npar,entry,isWrong)); 
    t2->Draw(name+Form(" >> h_%i_rs_data_%i_%i",Npar,entry,isWrong));
    t1->Draw(name+Form(" >> h1_%i_mc_%i_%i",Npar,entry,isWrong), mass_vetoe_cuts); 
    t2->Draw(name+Form(" >> h1_%i_rs_data_%i_%i",Npar,entry,isWrong), mass_vetoe_cuts);
    // t3->Draw(name+Form(" >> h_%i_ws_data_%i",Npar,entry));

    h_mc->SetLineColor(kBlue);
    h_rs_data->SetLineColor(kBlack);
    h_ws_data->SetLineColor(kRed);
    h1_mc->SetLineColor(kAzure+1);
    h1_rs_data->SetLineColor(kGreen+1);

    h_mc->GetXaxis()->SetTitle(name+" (MeV)");
    h_mc->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/nbins  ));

    h_rs_data->GetXaxis()->SetTitle(name+" (MeV)");
    h_rs_data->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/nbins  ));

    h1_mc->GetXaxis()->SetTitle(name+" (MeV)");
    h1_mc->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/nbins  ));

    h1_rs_data->GetXaxis()->SetTitle(name+" (MeV)");
    h1_rs_data->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/nbins  ));

    // h_ws_data->GetXaxis()->SetTitle(name+" (MeV)");
    // h_ws_data->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/nbins  ));

    if(isWrong)
    {
        h_mc->SetTitle("Wrong mass hypothesis");
        h_rs_data->SetTitle("Wrong mass hyopthesis");
        h1_mc->SetTitle("Wrong mass hypothesis");
        h1_rs_data->SetTitle("Wrong mass hyopthesis");
    }
    else
    {
        h_mc->SetTitle("Correct mass hypothesis");
        h_rs_data->SetTitle("Correct mass hyopthesis");
        h1_mc->SetTitle("Correct mass hypothesis");
        h1_rs_data->SetTitle("Correct mass hyopthesis");
    }

    h_mc->Scale(1./h_mc->Integral());
    h_rs_data->Scale(1./h_rs_data->Integral());
    h1_mc->Scale(1./h1_mc->Integral());
    h1_rs_data->Scale(1./h1_rs_data->Integral());
    // h_ws_data->Scale(1./h_ws_data->Integral());

    if( (h_rs_data->GetMaximum() >= h_mc->GetMaximum()) && (h_rs_data->GetMaximum() >= h1_mc->GetMaximum()) &&  (h_rs_data->GetMaximum() >= h1_rs_data->GetMaximum()) )
    {
        h_mc->GetYaxis()->SetRangeUser(0, 1.1*(h_rs_data->GetMaximum()));
    }
    else if( (h1_rs_data->GetMaximum() >= h_mc->GetMaximum()) && (h1_rs_data->GetMaximum() >= h1_mc->GetMaximum()) &&  (h1_rs_data->GetMaximum() >= h_rs_data->GetMaximum()) )
    {
        h_mc->GetYaxis()->SetRangeUser(0, 1.1*(h1_rs_data->GetMaximum()));
    }
    else if( (h1_mc->GetMaximum() >= h_mc->GetMaximum()) && (h1_mc->GetMaximum() >= h_rs_data->GetMaximum()) &&  (h1_mc->GetMaximum() >= h1_rs_data->GetMaximum()) )
    {
        h_mc->GetYaxis()->SetRangeUser(0, 1.1*(h1_mc->GetMaximum()));
    }
    else
    {
        h_mc->GetYaxis()->SetRangeUser(0, 1.1*(h_mc->GetMaximum()));
    }

    h_mc->Draw("HIST");
    h1_mc->Draw("HIST same");
    h_rs_data->Draw("HIST same");
    h1_rs_data->Draw("HIST same");
    // h_ws_data->Draw("HIST same");

    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_mc, "MC", "lp");
    leg->AddEntry(h_rs_data, "RS data", "lp");
    leg->AddEntry(h1_mc, "MC (after cut)", "lp");
    leg->AddEntry(h1_rs_data, "RS data (after cut)", "lp");
    // leg->AddEntry(h_ws_data, "WS data", "lp");
    leg->Draw("same");

    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/%i_particles/",Npar)+name+".pdf");
}


Int_t n_choose_k(Int_t n, Int_t k)
{
    return TMath::Factorial(n)/( TMath::Factorial(k)*TMath::Factorial(n-k) );
}
