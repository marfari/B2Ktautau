
using namespace std;

#define compare 2
// 0 = compares signal data and signal MC
// 1 = compares tail and core of signal MC
// 2 = compares signal MC with BuDD background MC

void comparison(){

      TString year = "2016";

      TString folder_name;
      if(compare == 0){folder_name = "Data_MC";}
      else if(compare == 1){folder_name = "Core_tail";}
      else if(compare == 2){folder_name = "Signal_vs_BuDD";}

      TString name1, name2;
      if(compare == 0){
            name1 = "MC";
            name2 = "Data";
      }
      else if(compare == 1){
            name1 = "Core";
            name2 = "Tail";
      }
      else if(compare == 2){
            name1 = "Signal MC";
            name2 = "BuDD bkg. MC";
      }

      TFileCollection *fc1; // signal data files on the grid
      if(compare == 0){
            fc1 = new TFileCollection("data_signal", "data_signal", "data_2016_MagUp.txt", 1);
            fc1->AddFromFile("data_2016_MagDown.txt", 1);
      }

      TChain *t1 = new TChain("ntuple/DecayTree");
      if((compare == 0) || (compare == 1) || (compare == 2)){t1->Add("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");}
      TChain *t2 = new TChain("ntuple/DecayTree");
      if(compare == 2){t2->Add("/panfs/felician/BDD_bkg/ROOT_Sim/2016/BuDD_mc_2016.root");}
      else if (compare == 0){t2->AddFileInfoList((TCollection*)fc1->GetList());}

      TString cut1;
      if( (compare == 0) || (compare == 2)){cut1 = "(abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (Bp_ConsBp_0_status == 0)";}
      else if(compare == 1){cut1 = "(abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (Bp_ConsBp_0_status == 0) && (abs(Bp_ConsBp_0_M - 5279) < 300)";}
      TString cut2;
      if(compare == 0){cut2 = "(Bp_ConsBp_0_status == 0)";}
      else if(compare == 1){cut2 =  "(abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (Bp_ConsBp_0_status == 0) && (Bp_ConsBp_0_M > 6500)";}
      else if(compare == 2){cut2 = "(abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 20213) && (abs(taum_TRUEID) == 20213) && (Bp_ConsBp_0_status == 0)";}

      TString variables[]  = {"Bp_BPVDIRA", "Bp_CONEDELTAETA_B", "Bp_CONEDELTAETA_K", "Bp_CONEDELTAETA_tau1", "Bp_CONEDELTAETA_tau2", "Bp_CONEMULT_B", "Bp_CONEMULT_K", "Bp_CONEMULT_tau1", "Bp_CONEMULT_tau2", "Bp_CONEPASYM_B", "Bp_CONEPASYM_K", "Bp_CONEPASYM_tau1", "Bp_CONEPASYM_tau2", "Bp_CONEPTASYM_B", "Bp_CONEPTASYM_K", "Bp_CONEPTASYM_tau1", "Bp_CONEPTASYM_tau2",
                              "Bp_ConsBp_0_Kplus_PE", "Bp_ConsBp_0_Kplus_PX", "Bp_ConsBp_0_Kplus_PY", "Bp_ConsBp_0_Kplus_PZ", "Bp_ConsBp_0_M", "Bp_ConsBp_0_MERR", "Bp_ConsBp_0_PE", "Bp_ConsBp_0_PX", "Bp_ConsBp_0_PY", "Bp_ConsBp_0_PZ", "Bp_ConsBp_0_PV_X", "Bp_ConsBp_0_PV_Y", "Bp_ConsBp_0_PV_Z", "Bp_ConsBp_0_ctau", "Bp_ConsBp_0_ctauErr", "Bp_ConsBp_0_decayLength", "Bp_ConsBp_0_decayLengthErr", 
                              "Bp_ConsBp_0_tauminus_0_M", "Bp_ConsBp_0_tauminus_0_MERR", "Bp_ConsBp_0_tauminus_0_PE", "Bp_ConsBp_0_tauminus_0_PX", "Bp_ConsBp_0_tauminus_0_PY", "Bp_ConsBp_0_tauminus_0_PZ", "Bp_ConsBp_0_tauminus_0_ctau", "Bp_ConsBp_0_tauminus_0_ctauErr", "Bp_ConsBp_0_tauminus_0_decayLength", "Bp_ConsBp_0_tauminus_0_decayLengthErr",
                              "Bp_ConsBp_0_tauminus_0_nu_tau_PE", "Bp_ConsBp_0_tauminus_0_nu_tau_PX", "Bp_ConsBp_0_tauminus_0_nu_tau_PY", "Bp_ConsBp_0_tauminus_0_nu_tau_PZ", "Bp_ConsBp_0_tauminus_0_piplus_0_PE", "Bp_ConsBp_0_tauminus_0_piplus_0_PX", "Bp_ConsBp_0_tauminus_0_piplus_0_PY", "Bp_ConsBp_0_tauminus_0_piplus_0_PZ", "Bp_ConsBp_0_tauminus_0_piplus_1_PE", "Bp_ConsBp_0_tauminus_0_piplus_1_PX", "Bp_ConsBp_0_tauminus_0_piplus_1_PY", "Bp_ConsBp_0_tauminus_0_piplus_1_PZ", "Bp_ConsBp_0_tauminus_0_piplus_PE", "Bp_ConsBp_0_tauminus_0_piplus_PX", "Bp_ConsBp_0_tauminus_0_piplus_PY", "Bp_ConsBp_0_tauminus_0_piplus_PZ", 
                              "Bp_ConsBp_0_tauminus_M", "Bp_ConsBp_0_tauminus_MERR", "Bp_ConsBp_0_tauminus_PE", "Bp_ConsBp_0_tauminus_PX", "Bp_ConsBp_0_tauminus_PY", "Bp_ConsBp_0_tauminus_PZ", "Bp_ConsBp_0_tauminus_ctau", "Bp_ConsBp_0_tauminus_ctauErr", "Bp_ConsBp_0_tauminus_decayLength", "Bp_ConsBp_0_tauminus_decayLengthErr",
                              "Bp_ConsBp_0_tauminus_nu_tau_PE", "Bp_ConsBp_0_tauminus_nu_tau_PX", "Bp_ConsBp_0_tauminus_nu_tau_PY", "Bp_ConsBp_0_tauminus_nu_tau_PZ", "Bp_ConsBp_0_tauminus_piplus_0_PE", "Bp_ConsBp_0_tauminus_piplus_0_PX", "Bp_ConsBp_0_tauminus_piplus_0_PY", "Bp_ConsBp_0_tauminus_piplus_0_PZ", "Bp_ConsBp_0_tauminus_piplus_1_PE", "Bp_ConsBp_0_tauminus_piplus_1_PX", "Bp_ConsBp_0_tauminus_piplus_1_PY", "Bp_ConsBp_0_tauminus_piplus_1_PZ", "Bp_ConsBp_0_tauminus_piplus_PE", "Bp_ConsBp_0_tauminus_piplus_PX", "Bp_ConsBp_0_tauminus_piplus_PY", "Bp_ConsBp_0_tauminus_piplus_PZ",
                              "Bp_DIRA_OWNPV", "Bp_ENDVERTEX_CHI2", "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_XERR", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_YERR", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_ZERR", "Bp_ETA", "Bp_FDCHI2_OWNPV", "Bp_FD_OWNPV", "Bp_M", "Bp_M01", "Bp_M02", "Bp_M03", "Bp_M04", "Bp_M05", "Bp_M06", "Bp_M012", "Bp_M013", "Bp_M014", "Bp_M015", "Bp_M016", "Bp_M023", "Bp_M024", "Bp_M025", "Bp_M026", "Bp_M034", "Bp_M035", "Bp_M036", "Bp_M045", "Bp_M046", "Bp_M056", "Bp_MCORR", "Bp_OWNPV_CHI2", "Bp_OWNPV_X", "Bp_OWNPV_XERR", "Bp_OWNPV_Y", "Bp_OWNPV_YERR", "Bp_OWNPV_Z", "Bp_OWNPV_ZERR",
                              "Bp_PE", "Bp_PHI", "Bp_PX", "Bp_PY", "Bp_PZ", "Bp_VCHI2PDOF", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISONUMVTX_B", 
                              "Kp_ETA", "Kp_M", "Kp_PE", "Kp_PHI", "Kp_PIDK", "Kp_PX", "Kp_PY", "Kp_PZ",  
                              "taum_DIRA_ORIVX", "taum_DIRA_OWNPV", "taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_X", "taum_ENDVERTEX_Y", "taum_ENDVERTEX_Z", "taum_ENDVERTEX_XERR", "taum_ENDVERTEX_YERR", "taum_ENDVERTEX_ZERR", "taum_FDCHI2_ORIVX", "taum_FDCHI2_OWNPV", "taum_FD_ORIVX", "taum_FD_OWNPV", 
                              "taum_M", "taum_M12", "taum_M13", "taum_M23", "taum_PE", "taum_PX", "taum_PY", "taum_PZ", 
                              "taup_DIRA_ORIVX", "taup_DIRA_OWNPV", "taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_X", "taup_ENDVERTEX_Y", "taup_ENDVERTEX_Z", "taup_ENDVERTEX_XERR", "taup_ENDVERTEX_YERR", "taup_ENDVERTEX_ZERR", "taup_FDCHI2_ORIVX", "taup_FDCHI2_OWNPV", "taup_FD_ORIVX", "taup_FD_OWNPV", 
                              "taup_M", "taup_M12", "taup_M13", "taup_M23", "taup_PE", "taup_PX", "taup_PY", "taup_PZ"
                            };

      // TString variables[]  = {"B^{+} DIRA to PV (^{#circ})", "B^{+} cone #Delta #eta", "K^{+} cone #Delta #eta", "#tau^{+} cone #Delta #eta", "#tau^{-} cone #Delta #eta", "B^{+} cone multiplicity", "K^{+} cone multiplicity", "#tau^{+} cone multiplicity", "#tau^{-} cone multiplicity", "B^{+} cone p^{asym}", "K^{+} cone p^{asym}", "#tau^{+} cone p^{asym}", "#tau^{-} cone p^{asym}", "B^{+} cone p_{T}^{asym}", "K^{+} cone p_{T}^{asym}", "#tau^{+} cone p_{T}^{asym}", "#tau^{-} cone p_{T}^{asym}",
      //                         "DTF K^{+} E}", "DTF K^{+} p^{x}", "DTF K^{+} p^{y}", "DTF K^{+} p^{z}", "DTF B^{+} M", "DTF B^{+} #Delta M", "DTF B^{+} E", "DTF B^{+} p^{x}", "DTF B^{+} p^{y}", "DTF B^{+} p^{z}", "DTF PV^{x}", "DTF PV^{y}", "DTF PV^{z}", "DTF B^{+} c#tau", "DTF B^{+} #Delta(c#tau)", "DTF B^{+} decay length", "DTF B^{+} error on decay length", 
      //                         "DTF #tau^{-} M", "DTF #tau^{-} #Delta M", "DTF #tau^{-} E", "DTF #tau^{-} p^{x}", "DTF #tau^{-} p^{y}", "DTF #tau^{-} p^{z}", "DTF #tau^{-} c#tau", "DTF #tau^{-} #Delta(c#tau)", "DTF #tau^{-} decay length", "DTF #tau^{-} error on decay length",
      //                         "DTF #nu_{#tau} E", "DTF #nu_{#tau} p^{x}", "DTF #nu_{#tau} p^{y}", "DTF #nu_{#tau} p^{z}", "DTF taum_pi1 E", "DTF taum_pi1 p^{x}", "DTF taum_pi1 p^{y}", "DTF taum_pi1 p^{z}", "DTF taum_pi2 E", "DTF taum_pi2 p^{x}", "DTF taum_pi2 p^{y}", "DTF taum_pi2 p^{z}", "DTF taum_pi3 E", "DTF taum_pi3 p^{x}", "DTF taum_pi3 p^{y}", "DTF taum_pi3 p^{z}", 
      //                         "DTF #tau^{+} M", "DTF #tau^{+} #DELTA M", "DTF #tau^{+} E", "DTF #tau^{+} p^{x}", "DTF #tau^{+} p^{y}", "DTF #tau^{+} p^{z}", "DTF #tau^{+} c#tau", "DTF #tau^{+} #Delta(c#tau)", "DTF #tau^{+} decay length", "DTF #tau^{+} error on decay length",
      //                         "DTF #bar{#nu}_{#tau} E", "DTF #bar{#nu}_{#tau} p^{x}", "DTF #bar{#nu}_{#tau} p^{y}", "DTF #bar{#nu}_{#tau} p^{z}", "DTF taup_pi1 E", "DTF taup_pi1 p^{x}", "DTF taup_pi1 p^{y}", "DTF taup_pi1 p^{z}", "DTF taup_pi2 E", "DTF taup_pi2 p^{x}", "DTF taup_pi2 p^{y}", "DTF taup_pi2 p^{z}", "DTF taup_pi3 E", "DTF taup_pi3 p^{x}", "DTF taup_pi3 p^{y}", "DTF taup_pi3 p^{z}",
      //                         "Reco. B^{+} DIRA to PV", "Reco. BV #chi^{2}", "Reco. BV^{x}", "Reco. #Delta(BV^{x})", "Reco. BV^{y}", "Reco. #Delta(BV^{y})", "Reco. BV^{z}", "Reco. #Delta(BV^{z})", "Reco. B^{+} #eta", "Reco. B^{+} FD to PV #chi^{2}", "Reco. B^{+} FD to PV", "Reco. B^{+} M", "Reco. M(taup_pi1 taup_pi2)", "Reco. M(taup_pi1 taup_pi3)", "Reco. M(taup_pi1 taum_pi1)", "Reco. M(taup_pi1 taum_pi2)", "Reco. M(taup_pi1 taum_pi3)", "Reco. M(taup_pi1 K)", "Reco. M(taup_pi2 taup_pi3)", "Reco. M(taup_pi2 taum_pi1)", "Reco. M(taup_pi2 taum_pi2)", "Reco. M(taup_pi2 taum_pi3)", "Reco. M(taup_pi2 K)", "Reco. M(taup_pi3 taum_pi1)", "Reco. M(taup_pi3 taum_pi2)", "Reco. M(taup_pi3 taum_pi3)", "Reco. M(taup_pi3 K)", "Reco. M(taum_pi1 taum_pi2)", "Reco. M(taum_pi1 taum_pi3)", "Reco. M(taum_pi1 K)", "Reco. M(taum_pi2 taum_pi3)", "M(taum_pi2 K)", "Reco. M(taum_pi3 K)", "Reco. B^{+} M^{corr}", "Reco. PV #chi^{2}", "Reco. PV^{x}", "Reco. #Delta(PV^{x})", "Reco. PV^{y}", "Reco. #Delta(PV^{y})", "Reco. PV^{z}", "Reco. #Delta(PV^{z})",
      //                         "Reco. B^{+} E", "Reco. B^{+} #phi", "Reco. B^{+} p^{x}", "Reco. B^{+} p^{y}", "Reco. B^{+} p^{z}", "Bp_VCHI2PDOF", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISONUMVTX_B", 
      //                         "Reco. K^{+} #eta", "Reco. K^{+} M", "Reco. K^{+} E", "Reco. K^{+} #phi", "Reco. K^{+} PID", "Reco. K^{+} p^{x}", "Reco. K^{+} p^{y}", "Reco. K^{+} p^{z}", 
      //                         "taum_DIRA_ORIVX", "taum_DIRA_OWNPV", "taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_X", "taum_ENDVERTEX_Y", "taum_ENDVERTEX_Z", "taum_ENDVERTEX_XERR", "taum_ENDVERTEX_YERR", "taum_ENDVERTEX_ZERR", "taum_FDCHI2_ORIVX", "taum_FDCHI2_OWNPV", "taum_FD_ORIVX", "taum_FD_OWNPV", 
      //                         "taum_M", "taum_M12", "taum_M13", "taum_M23", "taum_PE", "taum_PX", "taum_PY", "taum_PZ", 
      //                         "taup_DIRA_ORIVX", "taup_DIRA_OWNPV", "taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_X", "taup_ENDVERTEX_Y", "taup_ENDVERTEX_Z", "taup_ENDVERTEX_XERR", "taup_ENDVERTEX_YERR", "taup_ENDVERTEX_ZERR", "taup_FDCHI2_ORIVX", "taup_FDCHI2_OWNPV", "taup_FD_ORIVX", "taup_FD_OWNPV", 
      //                         "taup_M", "taup_M12", "taup_M13", "taup_M23", "taup_PE", "taup_PX", "taup_PY", "taup_PZ"
      //                       };

      std::vector<TCanvas *> c;
      gStyle->SetOptStat(0);

      t1->SetLineColor(kBlue);
      t1->SetFillColorAlpha(kBlue, 0.25);
      t2->SetLineColor(kRed);
      t2->SetFillColorAlpha(kRed, 0.25);

      int n_var =  sizeof(variables) / sizeof(variables[0]);

      for(int i = 0; i < n_var; i++){
            c.push_back(new TCanvas());
            c[i]->cd();

            t1->Draw(variables[i],cut1,"norm, HIST");
            t2->Draw(variables[i],cut2,"norm, HIST, SAME");
            TH1D *h = (TH1D*)gROOT->FindObject("htemp");
            h->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin())/h->GetNbinsX()) );
            h->SetTitle(variables[i]);
            gPad->Update();

            TLegend* leg = new TLegend(0.6, 0.6, 0.75, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->AddEntry(t1,name1);
            leg->AddEntry(t2,name2);
            leg->Draw("same");

            c[i]->SaveAs(folder_name+"/"+variables[i]+".png");
            c[i]->SaveAs(folder_name+"/"+variables[i]+".pdf");
      }

      return 0;
}