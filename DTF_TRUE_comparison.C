using namespace std;

void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram(RooRealVar var, TTree* t, int n, bool core);

using namespace std;

void DTF_TRUE_comparison(){

      TString year = "2016";

      TFile* fin_mc = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");
      TTree* t_mc = (TTree*)fin_mc->Get("ntuple/DecayTree");
 
 
      TString variables_DTF[] = {"Bp_ConsBp_Bp_ENDVERTEX_X[Bp_ConsBp_nIters-1]", "Bp_ConsBp_Bp_ENDVERTEX_Y[Bp_ConsBp_nIters-1]", "Bp_ConsBp_Bp_ENDVERTEX_Z[Bp_ConsBp_nIters-1]", "Bp_ConsBp_PV_X", "Bp_ConsBp_PV_Y", "Bp_ConsBp_PV_Z", "Bp_ConsBp_PE", "Bp_ConsBp_PX", "Bp_ConsBp_PY", "Bp_ConsBp_PZ",
                                 "Bp_ConsBp_Kplus_PE", "Bp_ConsBp_Kplus_PX", "Bp_ConsBp_Kplus_PY", "Bp_ConsBp_Kplus_PZ",
                                 "Bp_ConsBp_Taum_ENDVERTEX_X[Bp_ConsBp_nIters-1]", "Bp_ConsBp_Taum_ENDVERTEX_Y[Bp_ConsBp_nIters-1]", "Bp_ConsBp_Taum_ENDVERTEX_Z[Bp_ConsBp_nIters-1]", "Bp_ConsBp_tauminus_0_PE", "Bp_ConsBp_tauminus_0_PX", "Bp_ConsBp_tauminus_0_PY", "Bp_ConsBp_tauminus_0_PZ", 
                                 "Bp_ConsBp_tauminus_0_piplus_0_PE", "Bp_ConsBp_tauminus_0_piplus_0_PX", "Bp_ConsBp_tauminus_0_piplus_0_PY", "Bp_ConsBp_tauminus_0_piplus_0_PZ",
                                 "Bp_ConsBp_tauminus_0_piplus_1_PE", "Bp_ConsBp_tauminus_0_piplus_1_PX", "Bp_ConsBp_tauminus_0_piplus_1_PY", "Bp_ConsBp_tauminus_0_piplus_1_PZ",
                                 "Bp_ConsBp_tauminus_0_piplus_PE", "Bp_ConsBp_tauminus_0_piplus_PX", "Bp_ConsBp_tauminus_0_piplus_PY", "Bp_ConsBp_tauminus_0_piplus_PZ",
                                 "Bp_ConsBp_Taup_ENDVERTEX_X[Bp_ConsBp_nIters-1]", "Bp_ConsBp_Taup_ENDVERTEX_Y[Bp_ConsBp_nIters-1]", "Bp_ConsBp_Taup_ENDVERTEX_Z[Bp_ConsBp_nIters-1]", "Bp_ConsBp_tauminus_PE", "Bp_ConsBp_tauminus_PX", "Bp_ConsBp_tauminus_PY", "Bp_ConsBp_tauminus_PZ",
                                 "Bp_ConsBp_tauminus_piplus_0_PE", "Bp_ConsBp_tauminus_piplus_0_PX", "Bp_ConsBp_tauminus_piplus_0_PY", "Bp_ConsBp_tauminus_piplus_0_PZ",
                                 "Bp_ConsBp_tauminus_piplus_1_PE", "Bp_ConsBp_tauminus_piplus_1_PX", "Bp_ConsBp_tauminus_piplus_1_PY", "Bp_ConsBp_tauminus_piplus_1_PZ",
                                 "Bp_ConsBp_tauminus_piplus_PE", "Bp_ConsBp_tauminus_piplus_PX", "Bp_ConsBp_tauminus_piplus_PY", "Bp_ConsBp_tauminus_piplus_PZ"};

      TString variables_TRUE[] = {"Bp_TRUEENDVERTEX_X", "Bp_TRUEENDVERTEX_Y", "Bp_TRUEENDVERTEX_Z", "Bp_TRUEORIGINVERTEX_X", "Bp_TRUEORIGINVERTEX_Y", "Bp_TRUEORIGINVERTEX_Z", "Bp_TRUEP_E", "Bp_TRUEP_X", "Bp_TRUEP_Y", "Bp_TRUEP_Z",
                                  "Kp_TRUEP_E", "Kp_TRUEP_X", "Kp_TRUEP_Y", "Kp_TRUEP_Z", 
                                  "taum_TRUEENDVERTEX_X", "taum_TRUEENDVERTEX_Y", "taum_TRUEENDVERTEX_Z", "taum_TRUEP_E", "taum_TRUEP_X", "taum_TRUEP_Y", "taum_TRUEP_Z",
                                  "taum_pim0_TRUEP_E", "taum_pim0_TRUEP_X", "taum_pim0_TRUEP_Y", "taum_pim0_TRUEP_Z",
                                  "taum_pim1_TRUEP_E", "taum_pim1_TRUEP_X", "taum_pim1_TRUEP_Y", "taum_pim1_TRUEP_Z",
                                  "taum_pip0_TRUEP_E", "taum_pip0_TRUEP_X", "taum_pip0_TRUEP_Y", "taum_pip0_TRUEP_Z",
                                  "taup_TRUEENDVERTEX_X", "taup_TRUEENDVERTEX_Y", "taup_TRUEENDVERTEX_Z", "taup_TRUEP_E", "taup_TRUEP_X", "taup_TRUEP_Y", "taup_TRUEP_Z",
                                  "taup_pim0_TRUEP_E", "taup_pim0_TRUEP_X", "taup_pim0_TRUEP_Y", "taup_pim0_TRUEP_Z",
                                  "taup_pip0_TRUEP_E", "taup_pip0_TRUEP_X", "taup_pip0_TRUEP_Y", "taup_pip0_TRUEP_Z",
                                  "taup_pip1_TRUEP_E", "taup_pip1_TRUEP_X", "taup_pip1_TRUEP_Y", "taup_pip1_TRUEP_Z"};

      double h_min[] = {-40., -50., -40., -40., -40., -80., -120000., -25000., -25000., -120000., -1000., -400., -400., -1000., -40., -40., -40., -120000., -25000., -25000., -120000., -120000., -25000., -25000., -120000., -1000., -1000., -1000., -1000., -120000., -25000., -25000., -120000., -40., -40., -40., -120000., -25000., -25000., -120000., -1000., -200., -200., -1000., -1000., -1000., -1000., -1000., -1000., -1000., -1000., -1000.};
      double h_max[] = {40., 50., 40., 40., 40., 80., 120000., 25000., 25000., 120000., 1000., 400., 400., 1000., 40., 40., 40., 120000., 25000., 25000., 120000., 120000., 25000., 25000., 120000., 1000., 1000., 1000., 1000., 120000., 25000., 25000., 120000., 40., 40., 40., 120000., 25000., 25000., 120000., 1000., 200., 200., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.};

                        
      std::vector<TH1D*> histos;
      int n_bins = 80;
      int n_var = sizeof(variables_DTF)/sizeof(variables_DTF[0]);

      std::vector<TString> names;

      TFile* fout = new TFile("Mass_resolution_visualisation/DTF_vs_TRUE.root", "recreate");
      fout->cd();

      TString cut = "(Bp_ConsBp_status == 0) && (abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
      for(int i = 0; i < n_var; ++i){
            histos.push_back(new TH1D(Form("h_%.i",i), Form("h_%.i", i), 80, h_min[i], h_max[i]));  
            t_mc->Draw(variables_DTF[i]+"-"+variables_TRUE[i]+Form(" >> h_%.i",i), cut);

            histos[i]->GetXaxis()->SetTitle(variables_DTF[i]+"-"+variables_TRUE[i]);
            histos[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(histos[i]->GetXaxis()->GetXmax() - histos[i]->GetXaxis()->GetXmin())/histos[i]->GetNbinsX()) );
            histos[i]->Write();
      }
      fout->Close();

      return 0;

}

