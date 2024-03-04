
void compare( int year, TString FILES1, TString FILES2 )
{
    // TCut df_status = "df_status==0";

    int dimM = 22; // number of measured parameters
    int dimX = 23; // number of unknown parameters

    // TTree from DaVinci
    TFileCollection* fc1 = new TFileCollection("f1", "f1", FILES1);
    TFileCollection* fc2 = new TFileCollection("f2", "f2", FILES2, 1);

    TChain* t1 = new TChain("ntuple/DecayTree");
    TChain* t2 = new TChain("ntuple_SS/DecayTree");

    t1->AddFileInfoList((TCollection*)fc1->GetList());
    t2->AddFileInfoList((TCollection*)fc2->GetList(),10000);

    TCut truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
    TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
    TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
    TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
    TCut trigger = L0_trigger+HLT1_trigger+HLT2_trigger;

    TString var_names[] = {"Bp_VTXISONUMVTX_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B",
                           "Bp_VTXISONUMVTX_taup", "Bp_VTXISODCHI2ONETRACK_taup", "Bp_VTXISODCHI2MASSONETRACK_taup", "Bp_VTXISODCHI2TWOTRACK_taup", "Bp_VTXISODCHI2MASSTWOTRACK_taup",
                           "Bp_VTXISONUMVTX_taum", "Bp_VTXISODCHI2ONETRACK_taum", "Bp_VTXISODCHI2MASSONETRACK_taum", "Bp_VTXISODCHI2TWOTRACK_taum", "Bp_VTXISODCHI2MASSTWOTRACK_taum",
                           "Bp_CC_05_MULT_B", "Bp_CC_05_SPT_B", "Bp_CC_05_VPT_B", "Bp_CC_05_PX_B", "Bp_CC_05_PY_B", "Bp_CC_05_PZ_B", "Bp_CC_05_PASYM_B", "Bp_CC_05_PTASYM_B", "Bp_CC_05_PXASYM_B", "Bp_CC_05_PYASYM_B", "Bp_CC_05_PZASYM_B", "Bp_CC_05_DELTAETA_B", "Bp_CC_05_DELTAPHI_B", "Bp_CC_05_PX_B", "Bp_CC_05_IT_B", "Bp_CC_05_MAXPT_Q_B", "Bp_CC_05_MAXPT_PT_B", "Bp_CC_05_MAXPT_PX_B", "Bp_CC_05_MAXPT_PY_B", "Bp_CC_05_MAXPT_PZ_B", "Bp_CC_05_MAXPT_PE_B", "Bp_NC_05_MULT_B", "Bp_NC_05_SPT_B", "Bp_NC_05_VPT_B", "Bp_NC_05_PX_B", "Bp_NC_05_PY_B", "Bp_NC_05_PZ_B", "Bp_NC_05_PASYM_B", "Bp_NC_05_PTASYM_B", "Bp_NC_05_PXASYM_B", "Bp_NC_05_PYASYM_B", "Bp_NC_05_PZASYM_B", "Bp_NC_05_DELTAETA_B", "Bp_NC_05_DELTAPHI_B", "Bp_NC_05_IT_B", "Bp_NC_05_MAXPT_PT_B", "Bp_NC_05_MAXPT_PX_B", "Bp_NC_05_MAXPT_PY_B", "Bp_NC_05_MAXPT_PZ_B", "Bp_CCNC_05_IT_B",
                           "Bp_CC_10_MULT_B", "Bp_CC_10_SPT_B", "Bp_CC_10_VPT_B", "Bp_CC_10_PX_B", "Bp_CC_10_PY_B", "Bp_CC_10_PZ_B", "Bp_CC_10_PASYM_B", "Bp_CC_10_PTASYM_B", "Bp_CC_10_PXASYM_B", "Bp_CC_10_PYASYM_B", "Bp_CC_10_PZASYM_B", "Bp_CC_10_DELTAETA_B", "Bp_CC_10_DELTAPHI_B", "Bp_CC_10_PX_B", "Bp_CC_10_IT_B", "Bp_CC_10_MAXPT_Q_B", "Bp_CC_10_MAXPT_PT_B", "Bp_CC_10_MAXPT_PX_B", "Bp_CC_10_MAXPT_PY_B", "Bp_CC_10_MAXPT_PZ_B", "Bp_CC_10_MAXPT_PE_B", "Bp_NC_10_MULT_B", "Bp_NC_10_SPT_B", "Bp_NC_10_VPT_B", "Bp_NC_10_PX_B", "Bp_NC_10_PY_B", "Bp_NC_10_PZ_B", "Bp_NC_10_PASYM_B", "Bp_NC_10_PTASYM_B", "Bp_NC_10_PXASYM_B", "Bp_NC_10_PYASYM_B", "Bp_NC_10_PZASYM_B", "Bp_NC_10_DELTAETA_B", "Bp_NC_10_DELTAPHI_B", "Bp_NC_10_IT_B", "Bp_NC_10_MAXPT_PT_B", "Bp_NC_10_MAXPT_PX_B", "Bp_NC_10_MAXPT_PY_B", "Bp_NC_10_MAXPT_PZ_B", "Bp_CCNC_10_IT_B",
                           "Bp_CC_15_MULT_B", "Bp_CC_15_SPT_B", "Bp_CC_15_VPT_B", "Bp_CC_15_PX_B", "Bp_CC_15_PY_B", "Bp_CC_15_PZ_B", "Bp_CC_15_PASYM_B", "Bp_CC_15_PTASYM_B", "Bp_CC_15_PXASYM_B", "Bp_CC_15_PYASYM_B", "Bp_CC_15_PZASYM_B", "Bp_CC_15_DELTAETA_B", "Bp_CC_15_DELTAPHI_B", "Bp_CC_15_PX_B", "Bp_CC_15_IT_B", "Bp_CC_15_MAXPT_Q_B", "Bp_CC_15_MAXPT_PT_B", "Bp_CC_15_MAXPT_PX_B", "Bp_CC_15_MAXPT_PY_B", "Bp_CC_15_MAXPT_PZ_B", "Bp_CC_15_MAXPT_PE_B", "Bp_NC_15_MULT_B", "Bp_NC_15_SPT_B", "Bp_NC_15_VPT_B", "Bp_NC_15_PX_B", "Bp_NC_15_PY_B", "Bp_NC_15_PZ_B", "Bp_NC_15_PASYM_B", "Bp_NC_15_PTASYM_B", "Bp_NC_15_PXASYM_B", "Bp_NC_15_PYASYM_B", "Bp_NC_15_PZASYM_B", "Bp_NC_15_DELTAETA_B", "Bp_NC_15_DELTAPHI_B", "Bp_NC_15_IT_B", "Bp_NC_15_MAXPT_PT_B", "Bp_NC_15_MAXPT_PX_B", "Bp_NC_15_MAXPT_PY_B", "Bp_NC_15_MAXPT_PZ_B", "Bp_CCNC_15_IT_B",
                           "Bp_CC_20_MULT_B", "Bp_CC_20_SPT_B", "Bp_CC_20_VPT_B", "Bp_CC_20_PX_B", "Bp_CC_20_PY_B", "Bp_CC_20_PZ_B", "Bp_CC_20_PASYM_B", "Bp_CC_20_PTASYM_B", "Bp_CC_20_PXASYM_B", "Bp_CC_20_PYASYM_B", "Bp_CC_20_PZASYM_B", "Bp_CC_20_DELTAETA_B", "Bp_CC_20_DELTAPHI_B", "Bp_CC_20_PX_B", "Bp_CC_20_IT_B", "Bp_CC_20_MAXPT_Q_B", "Bp_CC_20_MAXPT_PT_B", "Bp_CC_20_MAXPT_PX_B", "Bp_CC_20_MAXPT_PY_B", "Bp_CC_20_MAXPT_PZ_B", "Bp_CC_20_MAXPT_PE_B", "Bp_NC_20_MULT_B", "Bp_NC_20_SPT_B", "Bp_NC_20_VPT_B", "Bp_NC_20_PX_B", "Bp_NC_20_PY_B", "Bp_NC_20_PZ_B", "Bp_NC_20_PASYM_B", "Bp_NC_20_PTASYM_B", "Bp_NC_20_PXASYM_B", "Bp_NC_20_PYASYM_B", "Bp_NC_20_PZASYM_B", "Bp_NC_20_DELTAETA_B", "Bp_NC_20_DELTAPHI_B", "Bp_NC_20_IT_B", "Bp_NC_20_MAXPT_PT_B", "Bp_NC_20_MAXPT_PX_B", "Bp_NC_20_MAXPT_PY_B", "Bp_NC_20_MAXPT_PZ_B", "Bp_CCNC_20_IT_B",
                           "Bp_VTXISOBDTHARDFIRSTVALUE_B", "Bp_VTXISOBDTHARDSECONDVALUE_B", "Bp_VTXISOBDTHARDTHIRDVALUE_B",
                           "Bp_TRKISOBDTFIRSTVALUE_K", "Bp_TRKISOBDTSECONDVALUE_K", "Bp_TRKISOBDTTHIRDVALUE_K",
                           "Bp_TRKISOBDTFIRSTVALUE_taup_pi1", "Bp_TRKISOBDTSECONDVALUE_taup_pi1", "Bp_TRKISOBDTTHIRDVALUE_taup_pi1",
                           "Bp_TRKISOBDTFIRSTVALUE_taup_pi2", "Bp_TRKISOBDTSECONDVALUE_taup_pi2", "Bp_TRKISOBDTTHIRDVALUE_taup_pi2",
                           "Bp_TRKISOBDTFIRSTVALUE_taup_pi3", "Bp_TRKISOBDTSECONDVALUE_taup_pi3", "Bp_TRKISOBDTTHIRDVALUE_taup_pi3",
                           "Bp_TRKISOBDTFIRSTVALUE_taum_pi1", "Bp_TRKISOBDTSECONDVALUE_taum_pi1", "Bp_TRKISOBDTTHIRDVALUE_taum_pi1",
                           "Bp_TRKISOBDTFIRSTVALUE_taum_pi2", "Bp_TRKISOBDTSECONDVALUE_taum_pi2", "Bp_TRKISOBDTTHIRDVALUE_taum_pi2",
                           "Bp_TRKISOBDTFIRSTVALUE_taum_pi3", "Bp_TRKISOBDTSECONDVALUE_taum_pi3", "Bp_TRKISOBDTTHIRDVALUE_taum_pi3",
                           "Kp_ProbNNk", "taup_pi1_ProbNNk", "taup_pi2_ProbNNk", "taup_pi3_ProbNNk", "taum_pi1_ProbNNk", "taum_pi2_ProbNNk", "taum_pi3_ProbNNk",
                           "Kp_TRGHOSTPROB"
                           };

    Double_t var_min[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -12000, 0, -1.05, 0, -6000, -6000, 0, 0, 0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -1.05, 0, -10000, -10000, 0, -1.05,
                          0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -12000, 0, -1.05, 0, -6000, -6000, 0, 0, 0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -1.05, 0, -10000, -10000, 0, -1.05,
                          0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -12000, 0, -1.05, 0, -6000, -6000, 0, 0, 0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -1.05, 0, -10000, -10000, 0, -1.05,
                          0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -12000, 0, -1.05, 0, -6000, -6000, 0, 0, 0, 0, 0, -20000, -20000, 0, -1, -1, -1, -1, -1, -1, 0, -1.05, 0, -10000, -10000, 0, -1.05,
                          -1.05, -1.05, -1.05, 
                          -0.5, -0.5, -0.5,
                          -0.5, -0.5, -0.5,
                          -0.5, -0.5, -0.5,
                          -0.5, -0.5, -0.5,
                          -0.5, -0.5, -0.5,
                          -0.5, -0.5, -0.5,
                          -0.5, -0.5, -0.5,
                          0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0
                         };
    
    Double_t var_max[] = {10, 100, 15000, 20, 15000, 10, 30, 5000, 100, 5000, 10, 100, 5000, 100, 5000,
                          20, 50000, 50000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 12000, 1.05, 1.05, 10000, 6000, 6000, 100000, 100000, 20, 20000, 20000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 1.05, 10000, 10000, 10000, 100000, 1.05,
                          20, 50000, 50000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 12000, 1.05, 1.05, 10000, 6000, 6000, 100000, 100000, 20, 20000, 20000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 1.05, 10000, 10000, 10000, 100000, 1.05,
                          20, 50000, 50000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 12000, 1.05, 1.05, 10000, 6000, 6000, 100000, 100000, 20, 20000, 20000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 1.05, 10000, 10000, 10000, 100000, 1.05,
                          20, 50000, 50000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 12000, 1.05, 1.05, 10000, 6000, 6000, 100000, 100000, 20, 20000, 20000, 20000, 20000, 100000, 1.05, 1.05, 1.05, 1.05, 1.05, 6, 3.2, 1.05, 10000, 10000, 10000, 100000, 1.05,                          
                          1, 1, 1, 
                          0.5, 0.5, 0.5,
                          0.5, 0.5, 0.5,
                          0.5, 0.5, 0.5,
                          0.5, 0.5, 0.5,
                          0.5, 0.5, 0.5,
                          0.5, 0.5, 0.5,
                          0.5, 0.5, 0.5,
                          1.05, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                          0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3
                         };

    Int_t N = sizeof(var_names)/sizeof(var_names[0]);
    gStyle->SetOptStat(0);
    
    std::vector<TH1D> h1, h2;
    for(int i = 0; i < N; i++)
    {
        h1.push_back( TH1D( Form("h1_%i",i), Form("h1_%i",i), 100, var_min[i], var_max[i] ) );
        h2.push_back( TH1D( Form("h2_%i",i), Form("h2_%i",i), 100, var_min[i], var_max[i] ) );

        t1->Draw(var_names[i]+Form(" >> h1_%i",i),truthMatch+trigger,"",50000);
        t2->Draw(var_names[i]+Form(" >> h2_%i",i),trigger,"",50000);

        TCanvas c = new TCanvas();
        c.cd();
        h1[i].GetXaxis()->SetTitle(var_names[i]);
        h1[i].GetYaxis()->SetTitle("Normalised entries (100 bins)");
        h1[i].SetTitle("MC (truth-match + trigger)");
        h2[i].GetXaxis()->SetTitle(var_names[i]);
        h2[i].GetYaxis()->SetTitle("Normalised entries (100 bins)");
        h2[i].SetTitle("WS data (trigger)");

        h2[i].SetLineColor(kRed);
        h2[i].SetFillColorAlpha(kRed, 0.25);
        h1[i].SetLineColor(kBlue);
        h1[i].SetFillColorAlpha(kBlue, 0.25);


        if( (h1[i].GetMaximum())/(h1[i].Integral()) > (h2[i].GetMaximum())/(h2[i].Integral()) )
        {
            h1[i].DrawNormalized();
            h2[i].DrawNormalized("same");
        }
        else
        {
            h2[i].DrawNormalized();
            h1[i].DrawNormalized("same");
        }

        TLegend* leg = new TLegend(0.6, 0.6, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.03);
        leg->AddEntry(&h1[i],Form("MC: m %0.2f w %0.2f",h1[i].GetMean(),h1[i].GetStdDev()),"f");
        leg->AddEntry(&h2[i],Form("WS data: m %0.2f w %0.2f",h2[i].GetMean(),h2[i].GetStdDev()),"f");
        leg->Draw("same");

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/comparisons/201%i/TruthMatch_trigger_MC_vs_WS_data/var_%i.gif",year,i));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/comparisons/201%i/TruthMatch_trigger_MC_vs_WS_data/var_%i.pdf",year,i));
    }
    
}