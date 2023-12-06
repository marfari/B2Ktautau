
void exact_constraints( int year, int species, TString RECO_files, TString FIT_files )
{
    TCut df_status = "(df_status==0) && (df_F_tolerance < pow(10,-6))";
    TCut passDTF = "Bp_dtf_12_status==0";

    Bool_t isMC = false;
    if( (species == 10) || (species == 11) || (species == 12) || (species == 1) )
    {
        isMC = true;
    }

    Double_t tolerance = pow(10,-6);

    // TTree from DaVinci
    TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    // TTree with results from the standalone fitter
    TFileCollection* fc1 = new TFileCollection("fc1", "fc1", FIT_files);
    TChain* t1 = new TChain("DecayTree");
    t1->AddFileInfoList((TCollection*)fc1->GetList());
    t->AddFriend(t1);

    int dimM = 22; // number of measured parameters
    int dimX = 19; // number of unknown parameters
    int dimC = 20; // number of exact constraints

    std::vector<TH1D> histos;
    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        histos.push_back(TH1D(Form("h_%i",i), Form("h_%i",i), 100, -tolerance, tolerance) );
    }

    for(int i = 0; i < dimM+dimX+dimC; i++)
    {
        t->Draw(Form("df_F_%i >> h_%i",i,i),df_status);
    }

    TString name[] = {
        "dL/dPVx (MeV/mm)",
        "dL/dPVy (MeV/mm)",
        "dL/dPVz (MeV/mm)",
        "dL/dDV1x (MeV/mm)",
        "dL/dDV1y (MeV/mm)",
        "dL/DV1z (MeV/mm)",
        "dL/dp3pi1x ()",
        "dL/dp3pi1y ()",
        "dL/dp3pi1z ()",
        "dL/dE3pi1 ()",
        "dL/dDV2x (MeV/mm)",
        "dL/dDV2y (MeV/mm)",
        "dL/dDV2z (MeV/mm)",
        "dL/dp3pi2x ()",
        "dL/dp3pi2y ()",
        "dL/dp3pi2z ()",
        "dL/dE3pi2 ()",
        "dL/dRPx (MeV/mm)",
        "dL/dRPy (MeV/mm)",
        "dL/dpKx ()",
        "dL/dpKy ()",
        "dL/dpKz ()",
        "dL/dBVx (MeV/mm)",
        "dL/dBVy (MeV/mm)",
        "dL/dBVz (MeV/mm)",
        "dL/dpBx ()",
        "dL/dpBy ()",
        "dL/dpBz ()",
        "dL/dMB_squared (1/MeV)",
        "dL/dptau1x ()",
        "dL/dptau1y ()",
        "dL/dptau1z ()",
        "dL/dpnu1x ()",
        "dL/dpnu1y ()",
        "dL/dpnu1z ()",
        "dL/dptau2x ()",
        "dL/dptau2y ()",
        "dL/dptau2z ()",
        "dL/dpnu2x ()",
        "dL/dpnu2y ()",
        "dL/dpnu2z ()",
        "pB must point back to PV (x,y) (MeV/mm)",
        "pB must point back to PV (x,z) (MeV/mm)",
        "ptau1 must point back to BV (x,y) (MeV/mm)",
        "ptau1 must point back to BV (x,z) (MeV/mm)",
        "Momentum conservation in DV1 (x) (MeV)",
        "Momentum conservation in DV1 (y) (MeV)",
        "Momentum conservation in DV1 (z) (MeV)",
        "Energy conservation in DV1 (MeV)",
        "ptau2 must point back to BV (x,y) (MeV/mm)",
        "ptau2 must point back to BV (x,z) (MeV/mm)",
        "Momentum conservation in DV2 (x) (MeV)",
        "Momentum conservation in DV2 (y) (MeV)",
        "Momentum conservation in DV2 (z) (MeV)",
        "Energy conservation in DV2 (MeV)",
        "BV must lie in K trajectory (x,y) (MeV/mm)",
        "BV must lie in K trajectory (x,z) (MeV/mm)",
        "Momentum conservation in BV (x) (MeV)",
        "Momentum conservation in BV (y) (MeV)",
        "Momentum conservation in BV (z) (MeV)",
        "Energy conservation in BV (MeV)"
    };

    for(int i = 0; i < dimM+dimX+dimC;i++)
    {
        TCanvas c = new TCanvas();
        histos[i].GetXaxis()->SetTitle(name[i]);
        histos[i].GetYaxis()->SetTitle( TString::Format("Events /(%g)",(histos[i].GetXaxis()->GetXmax() - histos[i].GetXaxis()->GetXmin())/100) );
        histos[i].SetTitle(name[i]);
        histos[i].Draw();
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/F%i.gif",year,species,i));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/F%i.pdf",year,species,i));
    }

    // B+ mass
    TH1D* h1_mb = new TH1D("h1_mb", "h1_mb", 100, 4000, 8000);
    TH1D* h2_mb = new TH1D("h2_mb", "h2_mb", 100, 4000, 8000);

    t->Draw("Bp_dtf_12_M >> h1_mb",passDTF);
    t->Draw("df_Bp_M >> h2_mb", df_status);

    gStyle->SetOptStat(0);

    TCanvas c1;
    c1.cd();
    h2_mb->GetXaxis()->SetTitle("m_{B} (MeV)");
    h2_mb->GetYaxis()->SetTitle( TString::Format("Events /(%g)",(h2_mb->GetXaxis()->GetXmax() - h2_mb->GetXaxis()->GetXmin())/100) );
    h2_mb->SetTitle("");
    h1_mb->SetLineColor(kBlack);
    h1_mb->SetFillColorAlpha(kBlack, 0.25);
    h2_mb->SetLineColor(kBlue);
    h2_mb->SetFillColorAlpha(kBlue, 0.25);
    h2_mb->DrawNormalized("same");
    h1_mb->DrawNormalized("same");

    TLegend* leg = new TLegend(0.6, 0.6, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(h1_mb,Form("DTF: m %0.2f w %0.2f", h1_mb->GetMean(), h1_mb->GetStdDev()),"f");
    leg->AddEntry(h2_mb,Form("Fitter: m %0.2f w %0.2f", h2_mb->GetMean(), h2_mb->GetStdDev()),"f");
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

    double DTF_num = t->GetEntries(passDTF);
    double fitter_num = t->GetEntries(df_status);
    double num_entries = t->GetEntries();

    double DTF_pass = DTF_num/num_entries;
    double fitter_pass = fitter_num/num_entries;

    cout << "Nominal DTF : passing rate of " << DTF_pass*100 << " \% " << endl;
    cout << "Standalone fit : passing rate of " << fitter_pass*100 << " \% " << endl;

    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass.gif",year,species));
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass.pdf",year,species));

    // B+ mass 2D histogram
    TH2D* h2 = new TH2D("h2", "h2", 100, 2000, 10000, 100, 2000, 10000);

    t->Draw("df_Bp_M : Bp_dtf_12_M >> h2",passDTF+df_status);

    TCanvas c1_2;
    c1_2.cd();
    h2->GetXaxis()->SetTitle("DTF mass (MeV)");
    h2->GetYaxis()->SetTitle("Fitter mass (MeV)");
    h2->SetTitle("2D plot for events that pass DTF and Fitter");
    h2->Draw("COLZ");
    c1_2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass_2D.gif",year,species));
    c1_2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass_2D.pdf",year,species));

    // B+ mass error
    TH1D* h1_mberr = new TH1D("h1_mberr", "h1_mberr", 100, 0, 1000);
    TH1D* h2_mberr = new TH1D("h2_mberr", "h2_mberr", 100, 0, 1000);
    t->Draw("Bp_dtf_12_MERR >> h1_mberr",passDTF);
    t->Draw("df_Bp_MERR >> h2_mberr", df_status);

    TCanvas c2;
    c2.cd();
    h2_mberr->GetXaxis()->SetTitle("#Delta m_{B} (MeV)");
    h2_mberr->GetYaxis()->SetTitle( TString::Format("Events /(%g)",(h2_mberr->GetXaxis()->GetXmax() - h2_mberr->GetXaxis()->GetXmin())/100) );
    h2_mberr->SetTitle("");
    h1_mberr->SetLineColor(kBlack);
    h1_mberr->SetFillColorAlpha(kBlack, 0.25);
    h2_mberr->SetLineColor(kBlue);
    h2_mberr->SetFillColorAlpha(kBlue, 0.25);
    h2_mberr->DrawNormalized("same");
    h1_mberr->DrawNormalized("same");

    TLegend* leg1 = new TLegend(0.6, 0.6, 0.85, 0.85);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h1_mberr,Form("DTF: m %0.2f w %0.2f",h1_mberr->GetMean(),h1_mberr->GetStdDev()),"f");
    leg1->AddEntry(h2_mberr,Form("Fitter: m %0.2f w %0.2f",h2_mberr->GetMean(),h2_mberr->GetStdDev()),"f");
    leg1->Draw("same");
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass_err.gif",year,species));
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass_err.pdf",year,species));

    // B+ mass pull
    if(isMC)
    {
        TH1D* h1_pull = new TH1D("h1_pull", "h1_pull", 100, -10, 10);
        TH1D* h2_pull = new TH1D("h2_pull", "h2_pull", 100, -10, 10);

        t->Draw("(Bp_dtf_12_M-5279)/Bp_dtf_12_MERR >> h1_pull",passDTF);
        t->Draw("(df_Bp_M - 5279)/df_Bp_MERR >> h2_pull", df_status);

        TCanvas c3;
        c3.cd();
        h1_pull->GetXaxis()->SetTitle("(m_{B}-m_{B}^{PDG})/#Delta m_{B}");
        h1_pull->GetYaxis()->SetTitle( TString::Format("Events /(%g)",(h1_pull->GetXaxis()->GetXmax() - h1_pull->GetXaxis()->GetXmin())/100) );
        h1_pull->SetTitle("");
        h1_pull->SetLineColor(kBlack);
        h1_pull->SetFillColorAlpha(kBlack, 0.25);
        h2_pull->SetLineColor(kBlue);
        h2_pull->SetFillColorAlpha(kBlue, 0.25);
        h1_pull->DrawNormalized();
        h2_pull->DrawNormalized("same");

        TLegend* leg2 = new TLegend(0.6, 0.6, 0.85, 0.85);
        leg2->SetBorderSize(0);
        leg2->SetTextSize(0.03);
        leg2->AddEntry(h1_pull,Form("DTF: m %0.2f w %0.2f",h1_pull->GetMean(),h1_pull->GetStdDev()),"f");
        leg2->AddEntry(h2_pull,Form("Fitter: m %0.2f w %0.2f",h2_pull->GetMean(),h2_pull->GetStdDev()),"f");
        leg2->Draw("same");
        c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass_pull.gif",year,species));
        c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/Bmass_pull.pdf",year,species));
    }

    Double_t range[] = {
        0.2, // PVx
        0.2, // PVy
        1,  // PVz
        2, // DV1x
        2, // DV1y
        20, // DV1z
        100, // p3pi1x
        100, // p3pi1y
        1000, // p3pi1z
        1000, // E3pi1
        2, // DV2x
        2, // DV2y
        20, // DV2z
        100, // p3pi2x
        100, // p3pi2y
        1000, // p3pi2z
        1000, // E3pi2
        1, // RPx
        1, // RPy
        200, // pKx
        200, // pKy
        1000, // pKz
        0.4, // BVx
        0.4, // BVy
        10, // BVz
        10000, // pBx
        10000, // pBy
        250000, // pBz
        5*pow(10,7), // mB^2
        10000, // ptau1x
        10000, // ptau1y
        200000, // ptau1z
        20000, // pnu1x
        20000, // pnu1y
        200000, // pnu1z
        10000, // ptau2x
        10000, // ptau2y
        200000, // ptau2z
        20000, // pnu2x
        20000, // pnu2y
        200000 // pnu2z
    };

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
        -10000, // pnu1x
        -10000, // pnu1y
        -100, // pnu1z
        -20000, // ptau2x
        -20000, // ptau2y
        0, // ptau2z
        -10000, // pnu2x
        -10000, // pnu2y
        -100 // pnu2z
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
        10000, // pnu1x
        10000, // pnu1y
        100000, // pnu1z
        20000, // ptau2x
        20000, // ptau2y
        500000, // ptau2z
        10000, // pnu2x
        10000, // pnu2y
        100000 // pnu2z
    };

    // Bias in fit parameters
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
        "df_Bp_M2",
        "df_taup_PX",
        "df_taup_PY",
        "df_taup_PZ",
        "df_antinutau_PX",
        "df_antinutau_PY",
        "df_antinutau_PZ",
        "df_taum_PX",
        "df_taum_PY",
        "df_taum_PZ",
        "df_nutau_PX",
        "df_nutau_PY",
        "df_nutau_PZ"
    };

    TString DTF_name[] = {
        "Bp_dtf_12_PV_X",
        "Bp_dtf_12_PV_Y",
        "Bp_dtf_12_PV_Z",
        "Bp_dtf_12_taup_ENDVERTEX_X",
        "Bp_dtf_12_taup_ENDVERTEX_Y",
        "Bp_dtf_12_taup_ENDVERTEX_Z",
        "(Bp_dtf_12_tauminus_piplus_PX + Bp_dtf_12_tauminus_piplus_0_PX + Bp_dtf_12_tauminus_piplus_1_PX)",
        "(Bp_dtf_12_tauminus_piplus_PY + Bp_dtf_12_tauminus_piplus_0_PY + Bp_dtf_12_tauminus_piplus_1_PY)",
        "(Bp_dtf_12_tauminus_piplus_PZ + Bp_dtf_12_tauminus_piplus_0_PZ + Bp_dtf_12_tauminus_piplus_1_PZ)",
        "(Bp_dtf_12_tauminus_piplus_PE + Bp_dtf_12_tauminus_piplus_0_PE + Bp_dtf_12_tauminus_piplus_1_PE)",
        "Bp_dtf_12_taum_ENDVERTEX_X",
        "Bp_dtf_12_taum_ENDVERTEX_Y",
        "Bp_dtf_12_taum_ENDVERTEX_Z",
        "(Bp_dtf_12_tauminus_0_piplus_PX + Bp_dtf_12_tauminus_0_piplus_0_PX + Bp_dtf_12_tauminus_0_piplus_1_PX)",
        "(Bp_dtf_12_tauminus_0_piplus_PY + Bp_dtf_12_tauminus_0_piplus_0_PY + Bp_dtf_12_tauminus_0_piplus_1_PY)",
        "(Bp_dtf_12_tauminus_0_piplus_PZ + Bp_dtf_12_tauminus_0_piplus_0_PZ + Bp_dtf_12_tauminus_0_piplus_1_PZ)",
        "(Bp_dtf_12_tauminus_0_piplus_PE + Bp_dtf_12_tauminus_0_piplus_0_PE + Bp_dtf_12_tauminus_0_piplus_1_PE)",
        "Bp_dtf_12_Bp_ENDVERTEX_X",
        "Bp_dtf_12_Bp_ENDVERTEX_Y",
        "Bp_dtf_12_Kplus_PX",
        "Bp_dtf_12_Kplus_PY",
        "Bp_dtf_12_Kplus_PZ",
        "Bp_dtf_12_Bp_ENDVERTEX_X",
        "Bp_dtf_12_Bp_ENDVERTEX_Y",
        "Bp_dtf_12_Bp_ENDVERTEX_Z",
        "Bp_dtf_12_PX",
        "Bp_dtf_12_PY",
        "Bp_dtf_12_PZ",
        "pow(Bp_dtf_12_M,2)",
        "Bp_dtf_12_tauminus_PX",
        "Bp_dtf_12_tauminus_PY",
        "Bp_dtf_12_tauminus_PZ",
        "(Bp_dtf_12_tauminus_PX - Bp_dtf_12_tauminus_piplus_PX - Bp_dtf_12_tauminus_piplus_0_PX - Bp_dtf_12_tauminus_piplus_1_PX)",
        "(Bp_dtf_12_tauminus_PY - Bp_dtf_12_tauminus_piplus_PY - Bp_dtf_12_tauminus_piplus_0_PY - Bp_dtf_12_tauminus_piplus_1_PY)",
        "(Bp_dtf_12_tauminus_PZ - Bp_dtf_12_tauminus_piplus_PZ - Bp_dtf_12_tauminus_piplus_0_PZ - Bp_dtf_12_tauminus_piplus_1_PZ)",
        "Bp_dtf_12_tauminus_0_PX",
        "Bp_dtf_12_tauminus_0_PY",
        "Bp_dtf_12_tauminus_0_PZ",
        "(Bp_dtf_12_tauminus_0_PX - Bp_dtf_12_tauminus_0_piplus_PX - Bp_dtf_12_tauminus_0_piplus_0_PX - Bp_dtf_12_tauminus_0_piplus_1_PX)",
        "(Bp_dtf_12_tauminus_0_PY - Bp_dtf_12_tauminus_0_piplus_PY - Bp_dtf_12_tauminus_0_piplus_0_PY - Bp_dtf_12_tauminus_0_piplus_1_PY)",
        "(Bp_dtf_12_tauminus_0_PZ - Bp_dtf_12_tauminus_0_piplus_PZ - Bp_dtf_12_tauminus_0_piplus_0_PZ - Bp_dtf_12_tauminus_0_piplus_1_PZ)"
    };

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
        "pow(5279.34,2)",
        "taup_TRUEP_X",
        "taup_TRUEP_Y",
        "taup_TRUEP_Z",
        "(taup_TRUEP_X - taup_pi1_TRUEP_X - taup_pi2_TRUEP_X - taup_pi3_TRUEP_X)",
        "(taup_TRUEP_Y - taup_pi1_TRUEP_Y - taup_pi2_TRUEP_Y - taup_pi3_TRUEP_Y)",
        "(taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z)",
        "taum_TRUEP_X",
        "taum_TRUEP_Y",
        "taum_TRUEP_Z",
        "(taum_TRUEP_X - taum_pi1_TRUEP_X - taum_pi2_TRUEP_X - taum_pi3_TRUEP_X)",
        "(taum_TRUEP_Y - taum_pi1_TRUEP_Y - taum_pi2_TRUEP_Y - taum_pi3_TRUEP_Y)",
        "(taum_TRUEP_Z - taum_pi1_TRUEP_Z - taum_pi2_TRUEP_Z - taum_pi3_TRUEP_Z)"
    };

    TString name1[] = {
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
        "MB^{2} (MeV^{2})",
        "ptau1x (MeV)",
        "ptau1y (MeV)",
        "ptau1z (MeV)",
        "pnu1x (MeV)",
        "pnu1y (MeV)",
        "pnu1z (MeV)",
        "ptau2x (MeV)",
        "ptau2y (MeV)",
        "ptau2z (MeV)",
        "pnu2x (MeV)",
        "pnu2y (MeV)",
        "pnu2z (MeV)"
    };

    std::vector<TH1D> df_histos;
    std::vector<TH1D> DTF_histos;
    if(isMC)
    {
        for(int i = 0; i < dimM+dimX; i++)
        {
            df_histos.push_back( TH1D(Form("df_h_%i",i), Form("df_h_%i",i), 100, -range[i], range[i]) );
            DTF_histos.push_back( TH1D(Form("DTF_h_%i",i), Form("DTF_h_%i",i), 100, -range[i], range[i]) );
        }

        for(int i = 0; i < dimM+dimX; i++)
        {
            t->Draw(df_name[i]+" - "+true_name[i]+Form(" >> df_h_%i",i),df_status);
            t->Draw(DTF_name[i]+" - "+true_name[i]+Form(" >> DTF_h_%i",i),df_status);
        }

        for(int i = 0; i < dimM+dimX;i++)
        {
            TCanvas c3 = new TCanvas();
            c3.cd();
            df_histos[i].GetXaxis()->SetTitle(name1[i]);
            df_histos[i].GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(df_histos[i].GetXaxis()->GetXmax() - df_histos[i].GetXaxis()->GetXmin())/100) );
            df_histos[i].SetTitle("Fit - true");
            DTF_histos[i].GetXaxis()->SetTitle(name1[i]);
            DTF_histos[i].GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(DTF_histos[i].GetXaxis()->GetXmax() - DTF_histos[i].GetXaxis()->GetXmin())/100) );
            DTF_histos[i].SetTitle("Fit - true");

            DTF_histos[i].SetLineColor(kBlack);
            DTF_histos[i].SetFillColorAlpha(kBlack, 0.25);
            df_histos[i].SetLineColor(kBlue);
            df_histos[i].SetFillColorAlpha(kBlue, 0.25);

            if( (df_histos[i].GetMaximum())/(df_histos[i].Integral()) > (DTF_histos[i].GetMaximum())/(DTF_histos[i].Integral()) )
            {
                df_histos[i].DrawNormalized();
                DTF_histos[i].DrawNormalized("same");
            }
            else
            {
                DTF_histos[i].DrawNormalized();
                df_histos[i].DrawNormalized("same");
            }

            TLegend* leg3 = new TLegend(0.6, 0.6, 0.85, 0.85);
            leg3->SetBorderSize(0);
            leg3->SetTextSize(0.03);
            leg3->AddEntry(&DTF_histos[i],Form("DTF: m %0.2f w %0.2f",DTF_histos[i].GetMean(),DTF_histos[i].GetStdDev()),"f");
            leg3->AddEntry(&df_histos[i],Form("Fitter: m %0.2f w %0.2f",df_histos[i].GetMean(),df_histos[i].GetStdDev()),"f");
            leg3->Draw("same");

            c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/var_%i.gif",year,species,i));
            c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/var_%i.pdf",year,species,i));
        }
    }
    else
    {
        for(int i = 0; i < dimM+dimX; i++)
        {
            df_histos.push_back( TH1D(Form("df_h_%i",i), Form("df_h_%i",i), 100, var_min[i], var_max[i]) );
            DTF_histos.push_back( TH1D(Form("DTF_h_%i",i), Form("DTF_h_%i",i), 100, var_min[i], var_max[i]) );
        }

        for(int i = 0; i < dimM+dimX; i++)
        {
            t->Draw(df_name[i]+Form(" >> df_h_%i",i),df_status);
            t->Draw(DTF_name[i]+Form(" >> DTF_h_%i",i),df_status);
        }

        for(int i = 0; i < dimM+dimX;i++)
        {
            TCanvas c3 = new TCanvas();
            c3.cd();
            df_histos[i].GetXaxis()->SetTitle(name1[i]);
            df_histos[i].GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(df_histos[i].GetXaxis()->GetXmax() - df_histos[i].GetXaxis()->GetXmin())/100) );
            df_histos[i].SetTitle("");
            DTF_histos[i].GetXaxis()->SetTitle(name1[i]);
            DTF_histos[i].GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(DTF_histos[i].GetXaxis()->GetXmax() - DTF_histos[i].GetXaxis()->GetXmin())/100) );
            DTF_histos[i].SetTitle("");

            DTF_histos[i].SetLineColor(kBlack);
            DTF_histos[i].SetFillColorAlpha(kBlack, 0.25);
            df_histos[i].SetLineColor(kBlue);
            df_histos[i].SetFillColorAlpha(kBlue, 0.25);

            if( (df_histos[i].GetMaximum())/(df_histos[i].Integral()) > (DTF_histos[i].GetMaximum())/(DTF_histos[i].Integral()) )
            {
                df_histos[i].DrawNormalized();
                DTF_histos[i].DrawNormalized("same");
            }
            else
            {
                DTF_histos[i].DrawNormalized();
                df_histos[i].DrawNormalized("same");
            }

            TLegend* leg3 = new TLegend(0.6, 0.6, 0.85, 0.85);
            leg3->SetBorderSize(0);
            leg3->SetTextSize(0.03);
            leg3->AddEntry(&DTF_histos[i],Form("DTF: m %0.2f w %0.2f",DTF_histos[i].GetMean(),DTF_histos[i].GetStdDev()),"f");
            leg3->AddEntry(&df_histos[i],Form("Fitter: m %0.2f w %0.2f",df_histos[i].GetMean(),df_histos[i].GetStdDev()),"f");
            leg3->Draw("same");

            c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/var_%i.gif",year,species,i));
            c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/var_%i.pdf",year,species,i));
        } 
    }

    // Fitter pull distributions
    TString df_err_name[] = {
        "df_PVx_err",
        "df_PVy_err",
        "df_PVz_err",
        "df_DV1x_err",
        "df_DV1y_err",
        "df_DV1z_err",
        "df_3p1_PX_err",
        "df_3p1_PY_err",
        "df_3p1_PZ_err",
        "df_3p1_PE_err",
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
        "df_Bp_M2_err",
        "df_taup_PX_err",
        "df_taup_PY_err",
        "df_taup_PZ_err",
        "df_antinutau_PX_err",
        "df_antinutau_PY_err",
        "df_antinutau_PZ_err",
        "df_taum_PX_err",
        "df_taum_PY_err",
        "df_taum_PZ_err",
        "df_nutau_PX_err",
        "df_nutau_PY_err",
        "df_nutau_PZ_err"
    };

    std::vector<TH1D> df_histos1;
    for(int i = 0; i < dimM+dimX; i++)
    {
        df_histos1.push_back( TH1D(Form("df_h1_%i",i), Form("df_h1_%i",i), 100, -10, 10) );
    }

    if(isMC)
    {
        for(int i = 0; i < dimM+dimX; i++)
        {
            t->Draw("("+df_name[i]+" - "+true_name[i]+") / "+df_err_name[i]+Form(" >> df_h1_%i",i),df_status);
        }

        for(int i = 0; i < dimM+dimX;i++)
        {
            TCanvas c4 = new TCanvas();
            c4.cd();
            df_histos1[i].GetXaxis()->SetTitle("Pull of "+name1[i]);
            df_histos1[i].GetYaxis()->SetTitle( TString::Format("Events /(%g)",(df_histos1[i].GetXaxis()->GetXmax() - df_histos1[i].GetXaxis()->GetXmin())/100) );
            df_histos1[i].SetTitle("");

            df_histos1[i].SetLineColor(kBlue);
            df_histos1[i].SetFillColorAlpha(kBlue, 0.25);

            df_histos1[i].Draw();

            c4.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/pull_%i.gif",year,species,i));
            c4.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/exact_constraints/201%i/Species_%i/pull_%i.pdf",year,species,i));
        }
    }
}
