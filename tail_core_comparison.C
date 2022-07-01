

using namespace std;

void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram(RooRealVar var, TTree* t, int n, bool core);

using namespace std;

void tail_core_comparison(){

      TString year = "2016";

      TFile* fin_mc = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");
      TTree* t_mc = (TTree*)fin_mc->Get("ntuple/DecayTree");
 
 
      TString variables[] = {"Bp_ENDVERTEX_CHI2", "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_XERR", "Bp_ENDVERTEX_YERR", "Bp_ENDVERTEX_ZERR", "Bp_DIRA_OWNPV", "Bp_FDCHI2_OWNPV", "Bp_FD_OWNPV", "Bp_IPCHI2_OWNPV", "Bp_IP_OWNPV", "Bp_OWNPV_CHI2", "Bp_OWNPV_X", "Bp_OWNPV_Y", "Bp_OWNPV_Z", "Bp_OWNPV_XERR", "Bp_OWNPV_YERR", "Bp_OWNPV_ZERR", "Bp_ConsBp_M", "Bp_M", "Bp_PE", "Bp_PT", "Bp_PX", "Bp_PY", "Bp_PZ", 
                             "Kp_PE", "Kp_PT", "Kp_PX", "Kp_PY", "Kp_PZ", 
                             "taup_pim0_PE", "taup_pim0_PT", "taup_pim0_PX", "taup_pim0_PY", "taup_pim0_PZ",
                             "taum_pim0_PE", "taum_pim0_PT", "taum_pim0_PX", "taum_pim0_PY", "taum_pim0_PZ",
                             "taum_pim1_PE", "taum_pim1_PT", "taum_pim1_PX", "taum_pim1_PY", "taum_pim1_PZ",
                             "taup_pip0_PE", "taup_pip0_PT", "taup_pip0_PX", "taup_pip0_PY", "taup_pip0_PZ",
                             "taup_pip1_PE", "taup_pip1_PT", "taup_pip1_PX", "taup_pip1_PY", "taup_pip1_PZ",
                             "taum_pip0_PE", "taum_pip0_PT", "taum_pip0_PX", "taum_pip0_PY", "taum_pip0_PZ",
                             "taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_X", "taum_ENDVERTEX_Y", "taum_ENDVERTEX_Z", "taum_ENDVERTEX_XERR", "taum_ENDVERTEX_YERR", "taum_ENDVERTEX_ZERR", "taum_DIRA_ORIVX", "taum_DIRA_OWNPV", "taum_FDCHI2_ORIVX", "taum_FDCHI2_OWNPV", "taum_FD_ORIVX", "taum_FD_OWNPV", "taum_IPCHI2_OWNPV", "taum_IP_OWNPV", "taum_M", "taum_PE", "taum_PT", "taum_PX", "taum_PY", "taum_PZ",
                             "taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_X", "taup_ENDVERTEX_Y", "taup_ENDVERTEX_Z", "taup_ENDVERTEX_XERR", "taup_ENDVERTEX_YERR", "taup_ENDVERTEX_ZERR", "taup_DIRA_ORIVX", "taup_DIRA_OWNPV", "taup_FDCHI2_ORIVX", "taup_FDCHI2_OWNPV", "taup_FD_ORIVX", "taup_FD_OWNPV", "taup_IPCHI2_OWNPV", "taup_IP_OWNPV", "taup_M", "taup_PE", "taup_PT", "taup_PX", "taup_PY", "taup_PZ"};

      TString x_labels[] = {"BV #chi^{2}", "BVx (mm)", "BVy (mm)", "BVz (mm)", "BVxerr (mm)", "BVyerr (mm)", "BVzerr (mm)", "B^{+} DIRA to PV (rad)", "B^{+} FD to PV #chi^{2}", "B^{+} FD to PV (mm)", "B^{+} IP to PV #chi^{2}", "B^{+} IP to PV (mm)", "PV #chi^{2}", "PVx (mm)", "PVy (mm)", "PVz (mm)", "PVxerr (mm)", "PVyerr (mm)", "PVzerr (mm)", "m_{B^{+}}^{DTF} (MeV)", "B^{+} m (MeV)", "B^{+} E (MeV)", "B^{+} p_{T} (MeV)", "B^{+} p_{x} (MeV)", "B^{+} p_{y} (MeV)", "B^{+} p_{z} (MeV)",
                            "K^{+} E (MeV)", "K^{+} p_{T} (MeV)", "K^{+} p_{x} (MeV)", "K^{+} p_{y} (MeV)", "K^{+} p_{z} (MeV)",
                            "#tau^{+} #pi^{-} E (MeV)", "#tau^{+} #pi^{-} p_{T} (MeV)", "#tau^{+} #pi^{-} p_{x} (MeV)", "#tau^{+} #pi^{-} p_{y} (MeV)", "#tau^{+} #pi^{-} p_{z} (MeV)",
                            "#tau^{-} #pi^{-}_{1} E (MeV)", "#tau^{-} #pi^{-}_{1} p_{T} (MeV)", "#tau^{-} #pi^{-}_{1} p_{x} (MeV)", "#tau^{-} #pi^{-}_{1} p_{y} (MeV)", "#tau^{-} #pi^{-}_{1} p_{z} (MeV)",
                            "#tau^{-} #pi^{-}_{2} E (MeV)", "#tau^{-} #pi^{-}_{2} p_{T} (MeV)", "#tau^{-} #pi^{-}_{2} p_{x} (MeV)", "#tau^{-} #pi^{-}_{2} p_{y} (MeV)", "#tau^{-} #pi^{-}_{2} p_{z} (MeV)",
                            "#tau^{+} #pi^{+}_{1} E (MeV)", "#tau^{+} #pi^{+}_{1} p_{T} (MeV)", "#tau^{+} #pi^{+}_{1} p_{x} (MeV)", "#tau^{+} #pi^{+}_{1} p_{y} (MeV)", "#tau^{+} #pi^{+}_{1} p_{z} (MeV)",
                            "#tau^{+} #pi^{+}_{2} E (MeV)", "#tau^{+} #pi^{+}_{2} p_{T} (MeV)", "#tau^{+} #pi^{+}_{2} p_{x} (MeV)", "#tau^{+} #pi^{+}_{2} p_{y} (MeV)", "#tau^{+} #pi^{+}_{2} p_{z} (MeV)", 
                            "#tau^{-} #pi^{+} E (MeV)", "#tau^{-} #pi^{+} p_{T} (MeV)", "#tau^{-} #pi^{+} p_{x} (MeV)", "#tau^{-} #pi^{+} p_{y} (MeV)", "#tau^{-} #pi^{+} p_{z} (MeV)", 
                            "#tau^{-} DV #chi^{2}", "#tau^{-} DVx (mm)", "#tau^{-} DVy (mm)", "#tau^{-} DVz (mm)", "#tau^{-} DVxerr (mm)", "#tau^{-} DVyerr (mm)", "#tau^{-} DVzerr (mm)", "#tau^{-} DIRA to BV (rad)", "#tau^{-} DIRA to PV (rad)", "#tau^{-} FD to BV #chi^{2}", "#tau^{-} FD to PV #chi^{2}", "#tau^{-} FD to BV (mm)", "#tau^{-} FD to PV (mm)", "#tau^{-} IP to PV #chi^{2}", "#tau^{-} IP to PV (mm)", "#tau^{-} m (MeV)", "#tau^{-} E (MeV)", "#tau^{-} p_{T} (MeV)", "#tau^{-} p_{x} (MeV)", "#tau^{-} p_{y} (MeV)", "#tau^{-} p_{z} (MeV)",
                            "#tau^{+} DV #chi^{2}", "#tau^{+} DVx (mm)", "#tau^{+} DVy (mm)", "#tau^{+} DVz (mm)", "#tau^{+} DVxerr (mm)", "#tau^{+} DVyerr (mm)", "#tau^{+} DVzerr (mm)", "#tau^{+} DIRA to BV (rad)", "#tau^{+} DIRA to PV (rad)", "#tau^{+} FD to BV #chi^{2}", "#tau^{+} FD to PV #chi^{2}", "#tau^{+} FD to BV (mm)", "#tau^{+} FD to PV (mm)", "#tau^{+} IP to PV #chi^{2}", "#tau^{+} IP to PV (mm)", "#tau^{+} m (MeV)", "#tau^{+} E (MeV)", "#tau^{+} p_{T} (MeV)", "#tau^{+} p_{x} (MeV)", "#tau^{+} p_{y} (MeV)", "#tau^{+} p_{z} (MeV)"};

      std::vector<TH1D*> histos_core;
      std::vector<TH1D*> histos_tail;
      int n_bins = 80;
      int n_var = sizeof(variables)/sizeof(variables[0]);

      RooWorkspace* ws = new RooWorkspace("ws");
      set_up_workspace_variables(*ws);

      std::vector<TString> names;

      for(int i=0; i<n_var; i++){
            histos_core.push_back(create_histogram((*ws->var(variables[i])), t_mc, n_bins, true));
            histos_tail.push_back(create_histogram((*ws->var(variables[i])), t_mc, n_bins, false));
            names.push_back(TString(variables[i]));
      }

      for(int i=0; i<n_var; i++){

            TCanvas a = new TCanvas();
            a.SetName(TString::Format("a_%i",i));
            a.SetTitle(" Core vs Tail - "+year);
            a.cd();
            histos_core[i]->SetYTitle( TString::Format("Candidates / (%g)",((ws->var(variables[i]))->getMax()-(ws->var(variables[i]))->getMin())/n_bins) );
            histos_core[i]->SetXTitle(x_labels[i]);
            histos_core[i]->SetStats(0);
            histos_core[i]->SetStats(0);
            histos_core[i]->SetTitle(" ");

            gStyle->SetOptStat(0);

            //normalization
            histos_tail[i]->Scale(1/histos_tail[i]->Integral());
            histos_core[i]->Scale(1/histos_core[i]->Integral());
            histos_core[i]->GetYaxis()->SetRangeUser(0.1*histos_core[i]->GetMinimum(),1.1*histos_core[i]->GetMaximum());
            histos_tail[i]->Draw("HIST");
            histos_core[i]->Draw("HIST same");    

            //y axis: maximum and minimum
            if((histos_core[i]->GetMaximum() > histos_tail[i]->GetMaximum())){
                  histos_tail[i]->GetYaxis()->SetRangeUser(0.1*histos_core[i]->GetMinimum(), 1.1*histos_core[i]->GetMaximum());}
            else if((histos_tail[i]->GetMaximum() > histos_core[i]->GetMaximum())){
                  histos_tail[i]->GetYaxis()->SetRangeUser(0.1*histos_tail[i]->GetMinimum(), 1.1*histos_tail[i]->GetMaximum());}

            //--TRATIO--//
            auto rp = new TRatioPlot(histos_tail[i], histos_core[i], "divsym");
            a.SetTicks(0, 1);
            rp->SetH1DrawOpt("HIST");
            rp->SetH2DrawOpt("HIST");
            rp->Draw();
            rp->GetLowerRefYaxis()->SetTitle("Core / Tail");
            rp->GetUpperRefYaxis()->SetTitle( TString::Format("Candidates / (%g)",((*ws->var(variables[i])).getMax()-(*ws->var(variables[i])).getMin())/n_bins) );
            rp->GetLowerRefXaxis()->SetTitle(x_labels[i]);
            rp->GetUpperRefXaxis()->SetTitle(x_labels[i]);
            rp->GetLowerRefGraph()->SetLineColor(kBlack); 
            rp->GetLowerRefGraph()->SetMarkerColor(kBlack);
            a.Update();
        
            TLegend* leg;
            leg = new TLegend(0.7, 0.8, 0.85, 0.9);
            leg->SetBorderSize(0);
            leg->AddEntry(histos_core[i], "Core", "lp");
            leg->AddEntry(histos_tail[i], "Tail", "lp");
            leg->SetTextSize(0.03);
            leg->Draw("same");

            a.SaveAs("Core_tail/core_vs_tail_"+names[i]+"_"+year+".pdf");
            a.SaveAs("Core_tail/core_vs_tail_"+names[i]+"_"+year+".gif");
      }

      return 0;

}

TH1D* create_histogram(RooRealVar var, TTree* t, int n, bool core){

      TH1D* h = new TH1D(var.GetName(), var.GetName(), n, var.getMin(), var.getMax());
      
      TString name_string = TString(var.GetName()) + ">>htemp(" + Form("%d",n) +"," + Form("%lf", var.getMin()) + "," + Form("%lf", var.getMax()) + ")";

      TString cut_core = "(Bp_ConsBp_status == 0) && (abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521) && (abs(Bp_ConsBp_M - 5279) < 300)";
      TString cut_tail = "(Bp_ConsBp_status == 0) &&(abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521) && (Bp_ConsBp_M > 6500)";

      if(core){t->Draw(name_string,cut_core);}
      else{t->Draw(name_string,cut_tail);}
      
      h = (TH1D*)gDirectory->Get("htemp")->Clone();
      h->SetTitle("");
      //h->SetMarkerStyle(29);,
      if(core){
            //h->SetMarkerColor(kBlue);
            h->SetLineColor(kBlue);
      }
      else{
            //h->SetMarkerColor(kRed);
            h->SetLineColor(kRed);
      }
      h->SetMarkerSize(1);
      return h;

}

void set_up_workspace_variables(RooWorkspace& w){

      float Bp_ENDVERTEX_CHI2_min, Bp_ENDVERTEX_CHI2_max;
      float Bp_ENDVERTEX_X_min, Bp_ENDVERTEX_X_max;
      float Bp_ENDVERTEX_Y_min, Bp_ENDVERTEX_Y_max;
      float Bp_ENDVERTEX_Z_min, Bp_ENDVERTEX_Z_max;
      float Bp_ENDVERTEX_XERR_min, Bp_ENDVERTEX_XERR_max;
      float Bp_ENDVERTEX_YERR_min, Bp_ENDVERTEX_YERR_max;
      float Bp_ENDVERTEX_ZERR_min, Bp_ENDVERTEX_ZERR_max;

      float Bp_DIRA_OWNPV_min, Bp_DIRA_OWNPV_max;
      float Bp_FDCHI2_OWNPV_min, Bp_FDCHI2_OWNPV_max;
      float Bp_FD_OWNPV_min, Bp_FD_OWNPV_max;
      float Bp_IPCHI2_OWNPV_min, Bp_IPCHI2_OWNPV_max;
      float Bp_IP_OWNPV_min, Bp_IP_OWNPV_max;
      float Bp_OWNPV_CHI2_min, Bp_OWNPV_CHI2_max;
      float Bp_OWNPV_X_min, Bp_OWNPV_X_max;
      float Bp_OWNPV_Y_min, Bp_OWNPV_Y_max;
      float Bp_OWNPV_Z_min, Bp_OWNPV_Z_max;
      float Bp_OWNPV_XERR_min, Bp_OWNPV_XERR_max;
      float Bp_OWNPV_YERR_min, Bp_OWNPV_YERR_max;
      float Bp_OWNPV_ZERR_min, Bp_OWNPV_ZERR_max;

      float Bp_ConsBp_M_min, Bp_ConsBp_M_max;
      float Bp_M_min, Bp_M_max;
      float Bp_PE_min, Bp_PE_max;
      float Bp_PT_min, Bp_PT_max;
      float Bp_PX_min, Bp_PX_max;
      float Bp_PY_min, Bp_PY_max;
      float Bp_PZ_min, Bp_PZ_max;

      float Kp_PE_min, Kp_PE_max;
      float Kp_PT_min, Kp_PT_max;
      float Kp_PX_min, Kp_PX_max;
      float Kp_PY_min, Kp_PY_max;
      float Kp_PZ_min, Kp_PZ_max;

      float taup_pim0_PE_min, taup_pim0_PE_max;
      float taup_pim0_PT_min, taup_pim0_PT_max;
      float taup_pim0_PX_min, taup_pim0_PX_max;
      float taup_pim0_PY_min, taup_pim0_PY_max;
      float taup_pim0_PZ_min, taup_pim0_PZ_max;

      float taum_pim0_PE_min, taum_pim0_PE_max;
      float taum_pim0_PT_min, taum_pim0_PT_max;
      float taum_pim0_PX_min, taum_pim0_PX_max;
      float taum_pim0_PY_min, taum_pim0_PY_max;
      float taum_pim0_PZ_min, taum_pim0_PZ_max;

      float taum_pim1_PE_min, taum_pim1_PE_max;
      float taum_pim1_PT_min, taum_pim1_PT_max;
      float taum_pim1_PX_min, taum_pim1_PX_max;
      float taum_pim1_PY_min, taum_pim1_PY_max;
      float taum_pim1_PZ_min, taum_pim1_PZ_max;

      float taup_pip0_PE_min, taup_pip0_PE_max;
      float taup_pip0_PT_min, taup_pip0_PT_max;
      float taup_pip0_PX_min, taup_pip0_PX_max;
      float taup_pip0_PY_min, taup_pip0_PY_max;
      float taup_pip0_PZ_min, taup_pip0_PZ_max;

      float taup_pip1_PE_min, taup_pip1_PE_max;
      float taup_pip1_PT_min, taup_pip1_PT_max;
      float taup_pip1_PX_min, taup_pip1_PX_max;
      float taup_pip1_PY_min, taup_pip1_PY_max;
      float taup_pip1_PZ_min, taup_pip1_PZ_max;

      float taum_pip0_PE_min, taum_pip0_PE_max;
      float taum_pip0_PT_min, taum_pip0_PT_max;
      float taum_pip0_PX_min, taum_pip0_PX_max;
      float taum_pip0_PY_min, taum_pip0_PY_max;
      float taum_pip0_PZ_min, taum_pip0_PZ_max;

      float taum_ENDVERTEX_CHI2_min, taum_ENDVERTEX_CHI2_max;
      float taum_ENDVERTEX_X_min, taum_ENDVERTEX_X_max;
      float taum_ENDVERTEX_Y_min, taum_ENDVERTEX_Y_max;
      float taum_ENDVERTEX_Z_min, taum_ENDVERTEX_Z_max;
      float taum_ENDVERTEX_XERR_min, taum_ENDVERTEX_XERR_max;
      float taum_ENDVERTEX_YERR_min, taum_ENDVERTEX_YERR_max;
      float taum_ENDVERTEX_ZERR_min, taum_ENDVERTEX_ZERR_max;

      float taum_DIRA_ORIVX_min, taum_DIRA_ORIVX_max;
      float taum_DIRA_OWNPV_min, taum_DIRA_OWNPV_max;
      float taum_FDCHI2_ORIVX_min, taum_FDCHI2_ORIVX_max;
      float taum_FDCHI2_OWNPV_min, taum_FDCHI2_OWNPV_max;
      float taum_FD_ORIVX_min, taum_FD_ORIVX_max;
      float taum_FD_OWNPV_min, taum_FD_OWNPV_max;
      float taum_IPCHI2_OWNPV_min, taum_IPCHI2_OWNPV_max;
      float taum_IP_OWNPV_min, taum_IP_OWNPV_max;

      float taum_M_min, taum_M_max;
      float taum_PE_min, taum_PE_max;
      float taum_PT_min, taum_PT_max;
      float taum_PX_min, taum_PX_max;
      float taum_PY_min, taum_PY_max;
      float taum_PZ_min, taum_PZ_max;

      float taup_ENDVERTEX_CHI2_min, taup_ENDVERTEX_CHI2_max;
      float taup_ENDVERTEX_X_min, taup_ENDVERTEX_X_max;
      float taup_ENDVERTEX_Y_min, taup_ENDVERTEX_Y_max;
      float taup_ENDVERTEX_Z_min, taup_ENDVERTEX_Z_max;
      float taup_ENDVERTEX_XERR_min, taup_ENDVERTEX_XERR_max;
      float taup_ENDVERTEX_YERR_min, taup_ENDVERTEX_YERR_max;
      float taup_ENDVERTEX_ZERR_min, taup_ENDVERTEX_ZERR_max;

      float taup_DIRA_ORIVX_min, taup_DIRA_ORIVX_max;
      float taup_DIRA_OWNPV_min, taup_DIRA_OWNPV_max;
      float taup_FDCHI2_ORIVX_min, taup_FDCHI2_ORIVX_max;
      float taup_FDCHI2_OWNPV_min, taup_FDCHI2_OWNPV_max;
      float taup_FD_ORIVX_min, taup_FD_ORIVX_max;
      float taup_FD_OWNPV_min, taup_FD_OWNPV_max;
      float taup_IPCHI2_OWNPV_min, taup_IPCHI2_OWNPV_max;
      float taup_IP_OWNPV_min, taup_IP_OWNPV_max;

      float taup_M_min, taup_M_max;
      float taup_PE_min, taup_PE_max;
      float taup_PT_min, taup_PT_max;
      float taup_PX_min, taup_PX_max;
      float taup_PY_min, taup_PY_max;
      float taup_PZ_min, taup_PZ_max;

      // B+ decay vertex
      Bp_ENDVERTEX_CHI2_min = 0.;
      Bp_ENDVERTEX_CHI2_max = 2000.;

      Bp_ENDVERTEX_X_min = -8.;
      Bp_ENDVERTEX_X_max = 10.;

      Bp_ENDVERTEX_Y_min = -8.;
      Bp_ENDVERTEX_Y_max = 10.;

      Bp_ENDVERTEX_Z_min = -150.;
      Bp_ENDVERTEX_Z_max = 250.;

      Bp_ENDVERTEX_XERR_min = 0.;
      Bp_ENDVERTEX_XERR_max = 0.15;

      Bp_ENDVERTEX_YERR_min = 0.;
      Bp_ENDVERTEX_YERR_max = 0.15;

      Bp_ENDVERTEX_ZERR_min = 0.;
      Bp_ENDVERTEX_ZERR_max = 4.;

      // Primary vertex (B+)
      Bp_DIRA_OWNPV_min = -1.5;
      Bp_DIRA_OWNPV_max = 1.5;

      Bp_FDCHI2_OWNPV_min = 0.;
      Bp_FDCHI2_OWNPV_max = 60000;

      Bp_FD_OWNPV_min = 0.;
      Bp_FD_OWNPV_max = 110.;

      Bp_IPCHI2_OWNPV_min = 0.;
      Bp_IPCHI2_OWNPV_max = 1500.;

      Bp_IP_OWNPV_min = 0.;
      Bp_IP_OWNPV_max = 1.5;

      Bp_OWNPV_CHI2_min = 0.;
      Bp_OWNPV_CHI2_max = 150.;

      Bp_OWNPV_X_min = 0.65;
      Bp_OWNPV_X_max = 1.05;

      Bp_OWNPV_Y_min = -0.4;
      Bp_OWNPV_Y_max = 0.;

      Bp_OWNPV_Z_min = -150.;
      Bp_OWNPV_Z_max = 150.;

      Bp_OWNPV_XERR_min = 0.;
      Bp_OWNPV_XERR_max = 0.06;

      Bp_OWNPV_YERR_min = 0.;
      Bp_OWNPV_YERR_max = 0.06;

      Bp_OWNPV_ZERR_min = 0.;
      Bp_OWNPV_ZERR_max = 0.5;

      // B+ kinematics
      Bp_ConsBp_M_min = 4000;
      Bp_ConsBp_M_max = 8000;

      Bp_M_min = 3000;
      Bp_M_max = 6000;
      
      Bp_PE_min = 0.;
      Bp_PE_max = 500000;

      Bp_PT_min = 0.;
      Bp_PT_max = 50000;

      Bp_PX_min = -40000.;
      Bp_PX_max = 40000;

      Bp_PY_min = -50000.;
      Bp_PY_max = 40000;

      Bp_PZ_min = 0.;
      Bp_PZ_max = 500000;

      // K+ kinematics
      Kp_PE_min = 0.;
      Kp_PE_max = 250000;

      Kp_PT_min = 0.;
      Kp_PT_max = 16000;

      Kp_PX_min = -10000.;
      Kp_PX_max = 15000;

      Kp_PY_min = -10000.;
      Kp_PY_max = 10000;

      Kp_PZ_min = 0.;
      Kp_PZ_max = 250000;

      // taup_pim0 
      taup_pim0_PE_min = 0.;
      taup_pim0_PE_max = 120000;

      taup_pim0_PT_min = 0.;
      taup_pim0_PT_max = 10000;

      taup_pim0_PX_min = -7000.;
      taup_pim0_PX_max = 7000;

      taup_pim0_PY_min = -7000.;
      taup_pim0_PY_max = 7000;

      taup_pim0_PZ_min = 0.;
      taup_pim0_PZ_max = 120000;

      // taum_pim0
      taum_pim0_PE_min = 0.;
      taum_pim0_PE_max = 120000;

      taum_pim0_PT_min = 0.;
      taum_pim0_PT_max = 10000;

      taum_pim0_PX_min = -7000.;
      taum_pim0_PX_max = 7000;

      taum_pim0_PY_min = -6000.;
      taum_pim0_PY_max = 6000;

      taum_pim0_PZ_min = 0.;
      taum_pim0_PZ_max = 120000;

      // taum_pim1
      taum_pim1_PE_min = 0.;
      taum_pim1_PE_max = 120000;

      taum_pim1_PT_min = 0.;
      taum_pim1_PT_max = 10000;

      taum_pim1_PX_min = -7000.;
      taum_pim1_PX_max = 7000;

      taum_pim1_PY_min = -7000.;
      taum_pim1_PY_max = 7000;

      taum_pim1_PZ_min = 0.;
      taum_pim1_PZ_max = 120000;

      // taup_pip0
      taup_pip0_PE_min = 0.;
      taup_pip0_PE_max = 120000;

      taup_pip0_PT_min = 0.;
      taup_pip0_PT_max = 10000;

      taup_pip0_PX_min = -7000.;
      taup_pip0_PX_max = 7000;

      taup_pip0_PY_min = -7000.;
      taup_pip0_PY_max = 7000;

      taup_pip0_PZ_min = 0.;
      taup_pip0_PZ_max = 120000;

      // taup_pip1
      taup_pip1_PE_min = 0.;
      taup_pip1_PE_max = 120000;

      taup_pip1_PT_min = 0.;
      taup_pip1_PT_max = 10000;

      taup_pip1_PX_min = -7000.;
      taup_pip1_PX_max = 7000;

      taup_pip1_PY_min = -7000.;
      taup_pip1_PY_max = 7000;

      taup_pip1_PZ_min = 0.;
      taup_pip1_PZ_max = 120000;

      // taum_pip0
      taum_pip0_PE_min = 0.;
      taum_pip0_PE_max = 120000;

      taum_pip0_PT_min = 0.;
      taum_pip0_PT_max = 10000;

      taum_pip0_PX_min = -7000.;
      taum_pip0_PX_max = 7000;

      taum_pip0_PY_min = -7000.;
      taum_pip0_PY_max = 7000;

      taum_pip0_PZ_min = 0.;
      taum_pip0_PZ_max = 120000;      

      // tau- decay vertex
      taum_ENDVERTEX_CHI2_min = 0.;
      taum_ENDVERTEX_CHI2_max = 20.;

      taum_ENDVERTEX_X_min = -7.;
      taum_ENDVERTEX_X_max = 7.;

      taum_ENDVERTEX_Y_min = -7.;
      taum_ENDVERTEX_Y_max = 7.;

      taum_ENDVERTEX_Z_min = -150.;
      taum_ENDVERTEX_Z_max = 300.;

      taum_ENDVERTEX_XERR_min = 0.;
      taum_ENDVERTEX_XERR_max = 0.20;

      taum_ENDVERTEX_YERR_min = 0.;
      taum_ENDVERTEX_YERR_max = 0.25;

      taum_ENDVERTEX_ZERR_min = 0.;
      taum_ENDVERTEX_ZERR_max = 6.;

      // tau- PV
      taum_DIRA_ORIVX_min = -1.5;
      taum_DIRA_ORIVX_max = 1.5;

      taum_DIRA_OWNPV_min = 0.996;
      taum_DIRA_OWNPV_max = 1.;

      taum_FDCHI2_ORIVX_min = 0.;
      taum_FDCHI2_ORIVX_max = 4000;

      taum_FDCHI2_OWNPV_min = 0.;
      taum_FDCHI2_OWNPV_max = 60000;

      taum_FD_ORIVX_min = 0.;
      taum_FD_ORIVX_max = 50.;

      taum_FD_OWNPV_min = 0.;
      taum_FD_OWNPV_max = 200.;

      taum_IPCHI2_OWNPV_min = 0.;
      taum_IPCHI2_OWNPV_max = 8000;

      taum_IP_OWNPV_min = 0.;
      taum_IP_OWNPV_max = 5.;

      // tau- kineamtics
      taum_M_min = 400.;
      taum_M_max = 1800.;

      taum_PE_min = 0.;
      taum_PE_max = 250000;

      taum_PT_min = 0.;
      taum_PT_max = 15000;

      taum_PX_min = -15000.;
      taum_PX_max = 15000;

      taum_PY_min = -15000.;
      taum_PY_max = 15000;

      taum_PZ_min = 0.;
      taum_PZ_max = 250000;  

      // tau+ decay vertex   
      taup_ENDVERTEX_CHI2_min = 0.;
      taup_ENDVERTEX_CHI2_max = 20.;

      taup_ENDVERTEX_X_min = -6.;
      taup_ENDVERTEX_X_max = 8.;

      taup_ENDVERTEX_Y_min = -6.;
      taup_ENDVERTEX_Y_max = 6.;

      taup_ENDVERTEX_Z_min = -150;
      taup_ENDVERTEX_Z_max = 250;

      taup_ENDVERTEX_XERR_min = 0.;
      taup_ENDVERTEX_XERR_max = 0.2;

      taup_ENDVERTEX_YERR_min = 0.;
      taup_ENDVERTEX_YERR_max = 0.25;

      taup_ENDVERTEX_ZERR_min = 0.;
      taup_ENDVERTEX_ZERR_max = 5.;

      // tau+ PV
      taup_DIRA_ORIVX_min = -1.5;
      taup_DIRA_ORIVX_max = 1.5;

      taup_DIRA_OWNPV_min = 0.996;
      taup_DIRA_OWNPV_max = 1.;

      taup_FDCHI2_ORIVX_min = 0.;
      taup_FDCHI2_ORIVX_max = 4000;

      taup_FDCHI2_OWNPV_min = 0.;
      taup_FDCHI2_OWNPV_max = 60000;

      taup_FD_ORIVX_min = 0.;
      taup_FD_ORIVX_max = 60;

      taup_FD_OWNPV_min = 0.;
      taup_FD_OWNPV_max = 200;

      taup_IPCHI2_OWNPV_min = 0.;
      taup_IPCHI2_OWNPV_max = 6000;

      taup_IP_OWNPV_min = 0.;
      taup_IP_OWNPV_max = 3.;

      // tau+ kinematics
      taup_M_min = 400;
      taup_M_max = 1800;

      taup_PE_min = 0.;
      taup_PE_max = 200000;

      taup_PT_min = 0.;
      taup_PT_max = 15000;

      taup_PX_min = -15000.;
      taup_PX_max = 15000;

      taup_PY_min = -15000.;
      taup_PY_max = 15000;

      taup_PZ_min = 0.;
      taup_PZ_max = 200000;  

      // B+ decay vertex
      RooRealVar Bp_ENDVERTEX_CHI2("Bp_ENDVERTEX_CHI2", "Bp_ENDVERTEX_CHI2", Bp_ENDVERTEX_CHI2_min, Bp_ENDVERTEX_CHI2_max);
      RooRealVar Bp_ENDVERTEX_X("Bp_ENDVERTEX_X", "Bp_ENDVERTEX_X", Bp_ENDVERTEX_X_min, Bp_ENDVERTEX_X_max, "mm");
      RooRealVar Bp_ENDVERTEX_Y("Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Y", Bp_ENDVERTEX_Y_min, Bp_ENDVERTEX_Y_max, "mm");
      RooRealVar Bp_ENDVERTEX_Z("Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_Z", Bp_ENDVERTEX_Z_min, Bp_ENDVERTEX_Z_max, "mm");
      RooRealVar Bp_ENDVERTEX_XERR("Bp_ENDVERTEX_XERR", "Bp_ENDVERTEX_XERR", Bp_ENDVERTEX_XERR_min, Bp_ENDVERTEX_XERR_max, "mm");
      RooRealVar Bp_ENDVERTEX_YERR("Bp_ENDVERTEX_YERR", "Bp_ENDVERTEX_YERR", Bp_ENDVERTEX_YERR_min, Bp_ENDVERTEX_YERR_max, "mm");
      RooRealVar Bp_ENDVERTEX_ZERR("Bp_ENDVERTEX_ZERR", "Bp_ENDVERTEX_ZERR", Bp_ENDVERTEX_ZERR_min, Bp_ENDVERTEX_ZERR_max, "mm");

      // Primary vertex (B+)
      RooRealVar Bp_DIRA_OWNPV("Bp_DIRA_OWNPV", "Bp_DIRA_OWNPV", Bp_DIRA_OWNPV_min, Bp_DIRA_OWNPV_max, "rad");
      RooRealVar Bp_FDCHI2_OWNPV("Bp_FDCHI2_OWNPV", "Bp_FDCHI2_OWNPV", Bp_FDCHI2_OWNPV_min, Bp_FDCHI2_OWNPV_max);
      RooRealVar Bp_FD_OWNPV("Bp_FD_OWNPV", "Bp_FD_OWNPV", Bp_FD_OWNPV_min, Bp_FD_OWNPV_max, "mm");
      RooRealVar Bp_IPCHI2_OWNPV("Bp_IPCHI2_OWNPV", "Bp_IPCHI2_OWNPV", Bp_IPCHI2_OWNPV_min, Bp_IPCHI2_OWNPV_max);
      RooRealVar Bp_IP_OWNPV("Bp_IP_OWNPV", "Bp_IP_OWNPV", Bp_IP_OWNPV_min, Bp_IP_OWNPV_max, "mm");
      RooRealVar Bp_OWNPV_CHI2("Bp_OWNPV_CHI2", "Bp_OWNPV_CHI2", Bp_OWNPV_CHI2_min, Bp_OWNPV_CHI2_max);
      RooRealVar Bp_OWNPV_X("Bp_OWNPV_X", "Bp_OWNPV_X", Bp_OWNPV_X_min, Bp_OWNPV_X_max, "mm");
      RooRealVar Bp_OWNPV_Y("Bp_OWNPV_Y", "Bp_OWNPV_Y", Bp_OWNPV_Y_min, Bp_OWNPV_Y_max, "mm");
      RooRealVar Bp_OWNPV_Z("Bp_OWNPV_Z", "Bp_OWNPV_Z", Bp_OWNPV_Z_min, Bp_OWNPV_Z_max, "mm");
      RooRealVar Bp_OWNPV_XERR("Bp_OWNPV_XERR", "Bp_OWNPV_XERR", Bp_OWNPV_XERR_min, Bp_OWNPV_XERR_max, "mm");
      RooRealVar Bp_OWNPV_YERR("Bp_OWNPV_YERR", "Bp_OWNPV_YERR", Bp_OWNPV_YERR_min, Bp_OWNPV_YERR_max, "mm");
      RooRealVar Bp_OWNPV_ZERR("Bp_OWNPV_ZERR", "Bp_OWNPV_ZERR", Bp_OWNPV_ZERR_min, Bp_OWNPV_ZERR_max, "mm");

      // B+ kineamtics
      RooRealVar Bp_ConsBp_M("Bp_ConsBp_M", "Bp_ConsBp_M", Bp_ConsBp_M_min, Bp_ConsBp_M_max, "MeV");
      RooRealVar Bp_M("Bp_M", "Bp_M", Bp_M_min, Bp_M_max, "MeV");
      RooRealVar Bp_PE("Bp_PE", "Bp_PE", Bp_PE_min, Bp_PE_max, "MeV");
      RooRealVar Bp_PT("Bp_PT", "Bp_PT", Bp_PT_min, Bp_PT_max, "MeV");
      RooRealVar Bp_PX("Bp_PX", "Bp_PX", Bp_PX_min, Bp_PX_max, "MeV");
      RooRealVar Bp_PY("Bp_PY", "Bp_PY", Bp_PY_min, Bp_PY_max, "MeV");
      RooRealVar Bp_PZ("Bp_PZ", "Bp_PZ", Bp_PZ_min, Bp_PZ_max, "MeV");

      // K+ kinematics
      RooRealVar Kp_PE("Kp_PE", "Kp_PE", Kp_PE_min, Kp_PE_max, "MeV");
      RooRealVar Kp_PT("Kp_PT", "Kp_PT", Kp_PT_min, Kp_PT_max, "MeV");
      RooRealVar Kp_PX("Kp_PX", "Kp_PX", Kp_PX_min, Kp_PX_max, "MeV");
      RooRealVar Kp_PY("Kp_PY", "Kp_PY", Kp_PY_min, Kp_PY_max, "MeV");
      RooRealVar Kp_PZ("Kp_PZ", "Kp_PZ", Kp_PZ_min, Kp_PZ_max, "MeV");

      // taup_pim0
      RooRealVar taup_pim0_PE("taup_pim0_PE", "taup_pim0_PE", taup_pim0_PE_min, taup_pim0_PE_max, "MeV");
      RooRealVar taup_pim0_PT("taup_pim0_PT", "taup_pim0_PT", taup_pim0_PT_min, taup_pim0_PT_max, "MeV");
      RooRealVar taup_pim0_PX("taup_pim0_PX", "taup_pim0_PX", taup_pim0_PX_min, taup_pim0_PX_max, "MeV");
      RooRealVar taup_pim0_PY("taup_pim0_PY", "taup_pim0_PY", taup_pim0_PY_min, taup_pim0_PY_max, "MeV");
      RooRealVar taup_pim0_PZ("taup_pim0_PZ", "taup_pim0_PZ", taup_pim0_PZ_min, taup_pim0_PZ_max, "MeV");

      // taum_pim0
      RooRealVar taum_pim0_PE("taum_pim0_PE", "taum_pim0_PE", taum_pim0_PE_min, taum_pim0_PE_max, "MeV");
      RooRealVar taum_pim0_PT("taum_pim0_PT", "taum_pim0_PT", taum_pim0_PT_min, taum_pim0_PT_max, "MeV");
      RooRealVar taum_pim0_PX("taum_pim0_PX", "taum_pim0_PX", taum_pim0_PX_min, taum_pim0_PX_max, "MeV");
      RooRealVar taum_pim0_PY("taum_pim0_PY", "taum_pim0_PY", taum_pim0_PY_min, taum_pim0_PY_max, "MeV");
      RooRealVar taum_pim0_PZ("taum_pim0_PZ", "taum_pim0_PZ", taum_pim0_PZ_min, taum_pim0_PZ_max, "MeV");

      // taum_pim1
      RooRealVar taum_pim1_PE("taum_pim1_PE", "taum_pim1_PE", taum_pim1_PE_min, taum_pim1_PE_max, "MeV");
      RooRealVar taum_pim1_PT("taum_pim1_PT", "taum_pim1_PT", taum_pim1_PT_min, taum_pim1_PT_max, "MeV");
      RooRealVar taum_pim1_PX("taum_pim1_PX", "taum_pim1_PX", taum_pim1_PX_min, taum_pim1_PX_max, "MeV");
      RooRealVar taum_pim1_PY("taum_pim1_PY", "taum_pim1_PY", taum_pim1_PY_min, taum_pim1_PY_max, "MeV");
      RooRealVar taum_pim1_PZ("taum_pim1_PZ", "taum_pim1_PZ", taum_pim1_PZ_min, taum_pim1_PZ_max, "MeV");

      // taup_pip0
      RooRealVar taup_pip0_PE("taup_pip0_PE", "taup_pip0_PE", taup_pip0_PE_min, taup_pip0_PE_max, "MeV");
      RooRealVar taup_pip0_PT("taup_pip0_PT", "taup_pip0_PT", taup_pip0_PT_min, taup_pip0_PT_max, "MeV");
      RooRealVar taup_pip0_PX("taup_pip0_PX", "taup_pip0_PX", taup_pip0_PX_min, taup_pip0_PX_max, "MeV");
      RooRealVar taup_pip0_PY("taup_pip0_PY", "taup_pip0_PY", taup_pip0_PY_min, taup_pip0_PY_max, "MeV");
      RooRealVar taup_pip0_PZ("taup_pip0_PZ", "taup_pip0_PZ", taup_pip0_PZ_min, taup_pip0_PZ_max, "MewV");

      // taup_pip1
      RooRealVar taup_pip1_PE("taup_pip1_PE", "taup_pip1_PE", taup_pip1_PE_min, taup_pip1_PE_max, "MeV");
      RooRealVar taup_pip1_PT("taup_pip1_PT", "taup_pip1_PT", taup_pip1_PT_min, taup_pip1_PT_max, "MeV");
      RooRealVar taup_pip1_PX("taup_pip1_PX", "taup_pip1_PX", taup_pip1_PX_min, taup_pip1_PX_max, "MeV");
      RooRealVar taup_pip1_PY("taup_pip1_PY", "taup_pip1_PY", taup_pip1_PY_min, taup_pip1_PY_max, "MeV");
      RooRealVar taup_pip1_PZ("taup_pip1_PZ", "taup_pip1_PZ", taup_pip1_PZ_min, taup_pip1_PZ_max, "MeV");

      // taum_pip0
      RooRealVar taum_pip0_PE("taum_pip0_PE", "taum_pip0_PE", taum_pip0_PE_min, taum_pip0_PE_max, "MeV");
      RooRealVar taum_pip0_PT("taum_pip0_PT", "taum_pip0_PT", taum_pip0_PT_min, taum_pip0_PT_max, "MeV");
      RooRealVar taum_pip0_PX("taum_pip0_PX", "taum_pip0_PX", taum_pip0_PX_min, taum_pip0_PX_max, "MeV");
      RooRealVar taum_pip0_PY("taum_pip0_PY", "taum_pip0_PY", taum_pip0_PY_min, taum_pip0_PY_max, "MeV");
      RooRealVar taum_pip0_PZ("taum_pip0_PZ", "taum_pip0_PZ", taum_pip0_PZ_min, taum_pip0_PZ_max, "MeV");

      // tau+ decay vertex
      RooRealVar taup_ENDVERTEX_CHI2("taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_CHI2", taup_ENDVERTEX_CHI2_min, taup_ENDVERTEX_CHI2_max);
      RooRealVar taup_ENDVERTEX_X("taup_ENDVERTEX_X", "taup_ENDVERTEX_X", taup_ENDVERTEX_X_min, taup_ENDVERTEX_X_max, "mm");
      RooRealVar taup_ENDVERTEX_Y("taup_ENDVERTEX_Y", "taup_ENDVERTEX_Y", taup_ENDVERTEX_Y_min, taup_ENDVERTEX_Y_max, "mm");
      RooRealVar taup_ENDVERTEX_Z("taup_ENDVERTEX_Z", "taup_ENDVERTEX_Z", taup_ENDVERTEX_Z_min, taup_ENDVERTEX_Z_max, "mm");
      RooRealVar taup_ENDVERTEX_XERR("taup_ENDVERTEX_XERR", "taup_ENDVERTEX_XERR", taup_ENDVERTEX_XERR_min, taup_ENDVERTEX_XERR_max, "mm");
      RooRealVar taup_ENDVERTEX_YERR("taup_ENDVERTEX_YERR", "taup_ENDVERTEX_YERR", taup_ENDVERTEX_YERR_min, taup_ENDVERTEX_YERR_max, "mm");
      RooRealVar taup_ENDVERTEX_ZERR("taup_ENDVERTEX_ZERR", "taup_ENDVERTEX_ZERR", taup_ENDVERTEX_ZERR_min, taup_ENDVERTEX_ZERR_max, "mm");

      // tau+ PV
      RooRealVar taup_DIRA_ORIVX("taup_DIRA_ORIVX", "taup_DIRA_ORIVX", taup_DIRA_ORIVX_min, taup_DIRA_ORIVX_max, "rad");
      RooRealVar taup_DIRA_OWNPV("taup_DIRA_OWNPV", "taup_DIRA_OWNPV", taup_DIRA_OWNPV_min, taup_DIRA_OWNPV_max, "rad");
      RooRealVar taup_FDCHI2_ORIVX("taup_FDCHI2_ORIVX", "taup_FDCHI2_ORIVX", taup_FDCHI2_ORIVX_min, taup_FDCHI2_ORIVX_max);
      RooRealVar taup_FDCHI2_OWNPV("taup_FDCHI2_OWNPV", "taup_FDCHI2_OWNPV", taup_FDCHI2_OWNPV_min, taup_FDCHI2_OWNPV_max);
      RooRealVar taup_FD_ORIVX("taup_FD_ORIVX", "taup_FD_ORIVX", taup_FD_ORIVX_min, taup_FD_ORIVX_max, "mm");
      RooRealVar taup_FD_OWNPV("taup_FD_OWNPV", "taup_FD_OWNPV", taup_FD_OWNPV_min, taup_FD_OWNPV_max, "mm");
      RooRealVar taup_IPCHI2_OWNPV("taup_IPCHI2_OWNPV", "taup_IPCHI2_OWNPV", taup_IPCHI2_OWNPV_min, taup_IPCHI2_OWNPV_max);
      RooRealVar taup_IP_OWNPV("taup_IP_OWNPV", "taup_IP_OWNPV", taup_IP_OWNPV_min, taup_IP_OWNPV_max, "mm");

      // tau+ kinematics
      RooRealVar taup_M("taup_M", "taup_M", taup_M_min, taup_M_max, "MeV");
      RooRealVar taup_PE("taup_PE", "taup_PE", taup_PE_min, taup_PE_max, "MeV");
      RooRealVar taup_PT("taup_PT", "taup_PT", taup_PT_min, taup_PT_max, "MeV");
      RooRealVar taup_PX("taup_PX", "taup_PX", taup_PX_min, taup_PX_max, "MeV");
      RooRealVar taup_PY("taup_PY", "taup_PY", taup_PY_min, taup_PY_max, "MeV");
      RooRealVar taup_PZ("taup_PZ", "taup_PZ", taup_PZ_min, taup_PZ_max, "MeV");

      // tau- decay vertex
      RooRealVar taum_ENDVERTEX_CHI2("taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_CHI2", taum_ENDVERTEX_CHI2_min, taum_ENDVERTEX_CHI2_max);
      RooRealVar taum_ENDVERTEX_X("taum_ENDVERTEX_X", "taum_ENDVERTEX_X", taum_ENDVERTEX_X_min, taum_ENDVERTEX_X_max, "mm");
      RooRealVar taum_ENDVERTEX_Y("taum_ENDVERTEX_Y", "taum_ENDVERTEX_Y", taum_ENDVERTEX_Y_min, taum_ENDVERTEX_Y_max, "mm");
      RooRealVar taum_ENDVERTEX_Z("taum_ENDVERTEX_Z", "taum_ENDVERTEX_Z", taum_ENDVERTEX_Z_min, taum_ENDVERTEX_Z_max, "mm");
      RooRealVar taum_ENDVERTEX_XERR("taum_ENDVERTEX_XERR", "taum_ENDVERTEX_XERR", taum_ENDVERTEX_XERR_min, taum_ENDVERTEX_XERR_max, "mm");
      RooRealVar taum_ENDVERTEX_YERR("taum_ENDVERTEX_YERR", "taum_ENDVERTEX_YERR", taum_ENDVERTEX_YERR_min, taum_ENDVERTEX_YERR_max, "mm");
      RooRealVar taum_ENDVERTEX_ZERR("taum_ENDVERTEX_ZERR", "taum_ENDVERTEX_ZERR", taum_ENDVERTEX_ZERR_min, taum_ENDVERTEX_ZERR_max, "mm");

      // tau- PV
      RooRealVar taum_DIRA_ORIVX("taum_DIRA_ORIVX", "taum_DIRA_ORIVX", taum_DIRA_ORIVX_min, taum_DIRA_ORIVX_max, "rad");
      RooRealVar taum_DIRA_OWNPV("taum_DIRA_OWNPV", "taum_DIRA_OWNPV", taum_DIRA_OWNPV_min, taum_DIRA_OWNPV_max, "rad");
      RooRealVar taum_FDCHI2_ORIVX("taum_FDCHI2_ORIVX", "taum_FDCHI2_ORIVX", taum_FDCHI2_ORIVX_min, taum_FDCHI2_ORIVX_max);
      RooRealVar taum_FDCHI2_OWNPV("taum_FDCHI2_OWNPV", "taum_FDCHI2_OWNPV", taum_FDCHI2_OWNPV_min, taum_FDCHI2_OWNPV_max);
      RooRealVar taum_FD_ORIVX("taum_FD_ORIVX", "taum_FD_ORIVX", taum_FD_ORIVX_min, taum_FD_ORIVX_max, "mm");
      RooRealVar taum_FD_OWNPV("taum_FD_OWNPV", "taum_FD_OWNPV", taum_FD_OWNPV_min, taum_FD_OWNPV_max, "mm");
      RooRealVar taum_IPCHI2_OWNPV("taum_IPCHI2_OWNPV", "taum_IPCHI2_OWNPV", taum_IPCHI2_OWNPV_min, taum_IPCHI2_OWNPV_max);
      RooRealVar taum_IP_OWNPV("taum_IP_OWNPV", "taum_IP_OWNPV", taum_IP_OWNPV_min, taum_IP_OWNPV_max, "mm");

      // tau- kinematics
      RooRealVar taum_M("taum_M", "taum_M", taum_M_min, taum_M_max, "MeV");
      RooRealVar taum_PE("taum_PE", "taum_PE", taum_PE_min, taum_PE_max, "MeV");
      RooRealVar taum_PT("taum_PT", "taum_PT", taum_PT_min, taum_PT_max, "MeV");
      RooRealVar taum_PX("taum_PX", "taum_PX", taum_PX_min, taum_PX_max, "MeV");
      RooRealVar taum_PY("taum_PY", "taum_PY", taum_PY_min, taum_PY_max, "MeV");
      RooRealVar taum_PZ("taum_PZ", "taum_PZ", taum_PZ_min, taum_PZ_max, "MeV");

      // B+ decay vertex
      w.import(Bp_ENDVERTEX_CHI2);
      w.import(Bp_ENDVERTEX_X);
      w.import(Bp_ENDVERTEX_Y);
      w.import(Bp_ENDVERTEX_Z);
      w.import(Bp_ENDVERTEX_XERR);
      w.import(Bp_ENDVERTEX_YERR);
      w.import(Bp_ENDVERTEX_ZERR);

      // Primary vertex (B+)
      w.import(Bp_DIRA_OWNPV);
      w.import(Bp_FDCHI2_OWNPV);
      w.import(Bp_FD_OWNPV);
      w.import(Bp_IPCHI2_OWNPV);
      w.import(Bp_IP_OWNPV);
      w.import(Bp_OWNPV_CHI2);
      w.import(Bp_OWNPV_X);
      w.import(Bp_OWNPV_Y);
      w.import(Bp_OWNPV_Z);
      w.import(Bp_OWNPV_XERR);
      w.import(Bp_OWNPV_YERR);
      w.import(Bp_OWNPV_ZERR);

      // B+ kinematics
      w.import(Bp_ConsBp_M);
      w.import(Bp_M);
      w.import(Bp_PE);
      w.import(Bp_PT);
      w.import(Bp_PX);
      w.import(Bp_PY);
      w.import(Bp_PZ);

      // K+ kinematics
      w.import(Kp_PE);
      w.import(Kp_PT);
      w.import(Kp_PX);
      w.import(Kp_PY);
      w.import(Kp_PZ);

      // taup_pim0
      w.import(taup_pim0_PE);
      w.import(taup_pim0_PT);
      w.import(taup_pim0_PX);
      w.import(taup_pim0_PY);
      w.import(taup_pim0_PZ);

      // taum_pim0
      w.import(taum_pim0_PE);
      w.import(taum_pim0_PT);
      w.import(taum_pim0_PX);
      w.import(taum_pim0_PY);
      w.import(taum_pim0_PZ);

      // taum_pim1
      w.import(taum_pim1_PE);
      w.import(taum_pim1_PT);
      w.import(taum_pim1_PX);
      w.import(taum_pim1_PY);
      w.import(taum_pim1_PZ);

      // taup_pip0
      w.import(taup_pip0_PE);
      w.import(taup_pip0_PT);
      w.import(taup_pip0_PX);
      w.import(taup_pip0_PY);
      w.import(taup_pip0_PZ);

      // taup_pip1
      w.import(taup_pip1_PE);
      w.import(taup_pip1_PT);
      w.import(taup_pip1_PX);
      w.import(taup_pip1_PY);
      w.import(taup_pip1_PZ);

      // taum_pip0
      w.import(taum_pip0_PE);
      w.import(taum_pip0_PT);
      w.import(taum_pip0_PX);
      w.import(taum_pip0_PY);
      w.import(taum_pip0_PZ);

      // tau- decay vertex
      w.import(taum_ENDVERTEX_CHI2);
      w.import(taum_ENDVERTEX_X);
      w.import(taum_ENDVERTEX_Y);
      w.import(taum_ENDVERTEX_Z);
      w.import(taum_ENDVERTEX_XERR);
      w.import(taum_ENDVERTEX_YERR);
      w.import(taum_ENDVERTEX_ZERR);

      // tau- PV
      w.import(taum_DIRA_ORIVX);
      w.import(taum_DIRA_OWNPV);
      w.import(taum_FDCHI2_ORIVX);
      w.import(taum_FDCHI2_OWNPV);
      w.import(taum_FD_ORIVX);
      w.import(taum_FD_OWNPV);
      w.import(taum_IPCHI2_OWNPV);
      w.import(taum_IP_OWNPV);

      // tau- kinematics
      w.import(taum_M);
      w.import(taum_PE);
      w.import(taum_PT);
      w.import(taum_PX);
      w.import(taum_PY);
      w.import(taum_PZ);

      // tau+ decay vertex
      w.import(taup_ENDVERTEX_CHI2);
      w.import(taup_ENDVERTEX_X);
      w.import(taup_ENDVERTEX_Y);
      w.import(taup_ENDVERTEX_Z);
      w.import(taup_ENDVERTEX_XERR);
      w.import(taup_ENDVERTEX_YERR);
      w.import(taup_ENDVERTEX_ZERR);

      // tau+ PV
      w.import(taup_DIRA_ORIVX);
      w.import(taup_DIRA_OWNPV);
      w.import(taup_FDCHI2_ORIVX);
      w.import(taup_FDCHI2_OWNPV);
      w.import(taup_FD_ORIVX);
      w.import(taup_FD_OWNPV);
      w.import(taup_IPCHI2_OWNPV);
      w.import(taup_IP_OWNPV);

      // tau+ kinematics
      w.import(taup_M);
      w.import(taup_PE);
      w.import(taup_PT);
      w.import(taup_PX);
      w.import(taup_PY);
      w.import(taup_PZ);

      return 0;
}
