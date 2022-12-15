
using namespace std;

#define compare 0
// 0 = compares signal data and signal MC
// 1 = compares tail and core of signal MC
// 2 = compares signal MC with BuDD background MC
// 3 = compares signal MC and BuDDK cocktail background MC
// 4 = compares signal MC and BdDDK cocktail background MC

#define oneDplots 1
// 0 = does not make 1D plots
// 1 = makes 1D plots

void comparison(){

      TString year = "2018";

      TString folder_name;
      if(compare == 0){folder_name = "Data_MC_"+year;}
      else if(compare == 1){folder_name = "Core_tail_"+year;}
      else if(compare == 2){folder_name = "Signal_vs_BuDD_"+year;}
      else if(compare == 3){folder_name = "Signal_vs_BuDDK_"+year;}
      else if(compare == 4){folder_name = "Signal_vs_BdDDK_"+year;}

      TString name1, name2;
      if(compare == 0){
            name1 = "Simulated signal";
            name2 = "Data";
      }
      else if(compare == 1){
            name1 = "Core";
            name2 = "Tail";
      }
      else if(compare == 2){
            name1 = "Signal MC";
            name2 = "BuDD MC";
      }
      else if(compare == 3){
            name1 = "Signal MC";
            name2 = "BuDDK cocktail MC";
      }
      else if(compare == 4){
            name1 = "Signal MC";
            name2 = "BdDDK cocktail MC";
      }

      TFileCollection *fc1; // signal data files on the grid
      if(compare == 0){
            fc1 = new TFileCollection("data_signal", "data_signal", "data_"+year+"_MagUp.txt", 5);
            fc1->AddFromFile("data_"+year+"_MagDown.txt", 5);
      }

      TChain *t1 = new TChain("ntuple/DecayTree");
      if((compare == 0) || (compare == 1) || (compare == 2) || (compare == 3) || (compare == 4)){t1->Add("/panfs/felician/B2Ktautau/ROOT_Sim/"+year+"/mc_"+year+".root");}
      TChain *t2 = new TChain("ntuple/DecayTree");
      if(compare == 2){t2->Add("/panfs/felician/BDD_bkg/ROOT_Sim/"+year+"/BuDD_mc_"+year+".root");}
      else if(compare == 0){t2->AddFileInfoList((TCollection*)fc1->GetList());}
      else if(compare == 3){t2->Add("/panfs/felician/BuDDK_cocktail/ROOT_Sim/"+year+"/mc_"+year+"_BuDDK_cocktail.root");}
      else if(compare == 4){t2->Add("/panfs/felician/BdDDK_cocktail/ROOT_Sim/"+year+"/mc_"+year+"_BdDDK_cocktail.root");}

      // Additional variables
      TString tau_separation = "sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";

      TString theta = " acos( ( Bp_ConsBp_0_PX*(Bp_ENDVERTEX_X - Bp_OWNPV_X) + Bp_ConsBp_0_PY*(Bp_ENDVERTEX_Y - Bp_OWNPV_Y) + Bp_ConsBp_0_PZ*(Bp_ENDVERTEX_Z - Bp_OWNPV_Z) )/( sqrt(pow(Bp_ConsBp_0_PX,2) + pow(Bp_ConsBp_0_PY,2) + pow(Bp_ConsBp_0_PZ,2))*sqrt(pow(Bp_ENDVERTEX_X-Bp_OWNPV_X,2) + pow(Bp_ENDVERTEX_Y-Bp_OWNPV_Y,2) + pow(Bp_ENDVERTEX_Z-Bp_OWNPV_Z,2) ) ) )  ";

      TString deltaBpPx = "Bp_ConsBp_0_PX-Bp_PX";
      TString deltaBpPy = "Bp_ConsBp_0_PY-Bp_PY";
      TString deltaBpPz = "Bp_ConsBp_0_PZ-Bp_PZ";

      TString cut1;
      if( (compare == 0) || (compare == 2) || (compare == 3) || (compare == 4)){cut1 = "(abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (Bp_ConsBp_0_status == 0)";}
      else if(compare == 1){cut1 = "(abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (Bp_ConsBp_0_status == 0) && (abs(Bp_ConsBp_0_M - 5279) < 300)";}
      cut1+="&& (Bp_ConsBp_0_tauminus_0_decayLength > 0.5) && (Bp_ConsBp_0_tauminus_decayLength > 0.5) && (Bp_ConsBp_0_decayLength > 5) && ("+tau_separation+" > 0.5) && (Bp_ConsBp_0_chi2 < 20)";
      //cut1+="&& (Bp_ConsBp_0_tauminus_0_decayLength > 0.5) && (Bp_ConsBp_0_tauminus_decayLength > 0.5) && (Bp_ConsBp_0_decayLength > 5) && ( sqrt(pow(Bp_ConsBp_0_tauminus_0_PE + Bp_ConsBp_0_tauminus_PE,2) - pow(Bp_ConsBp_0_tauminus_0_PX + Bp_ConsBp_0_tauminus_PX,2) - pow(Bp_ConsBp_0_tauminus_0_PY + Bp_ConsBp_0_tauminus_PY,2) - pow(Bp_ConsBp_0_tauminus_0_PZ + Bp_ConsBp_0_tauminus_PZ,2)) < 4000) && ( sqrt(pow(Kp_PE + Bp_ConsBp_0_tauminus_0_PE,2) - pow(Kp_PX + Bp_ConsBp_0_tauminus_0_PX,2) - pow(Kp_PY + Bp_ConsBp_0_tauminus_0_PY,2) - pow(Kp_PZ + Bp_ConsBp_0_tauminus_0_PZ,2)) < 5500)";
      TString cut2;
      if(compare == 0){cut2 = "(Bp_ConsBp_0_status == 0)";}
      else if(compare == 1){cut2 =  "(abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (Bp_ConsBp_0_status == 0) && (Bp_ConsBp_0_M > 6500)";}
      else if(compare == 2){cut2 = "(abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_TRUEID) == 20213) && (abs(taum_TRUEID) == 20213) && (Bp_ConsBp_0_status == 0)";}
      else if(compare == 3){cut2 = "(Bp_ConsBp_0_status[0] == 0) && (abs(Bp_TRUEID) == 521) && (abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID == 211) && abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) != 521) && (abs(taum_TRUEID) != 521)";}
      else if(compare == 4){cut2 = "(Bp_ConsBp_0_status[0] == 0) && (abs(Bp_TRUEID) == 511) && (abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) != 511) && (abs(taum_TRUEID) != 511)";}
      cut2+="&& (Bp_ConsBp_0_tauminus_0_decayLength > 0.5) && (Bp_ConsBp_0_tauminus_decayLength > 0.5) && (Bp_ConsBp_0_decayLength > 5) && ("+tau_separation+" > 0.5) && (Bp_ConsBp_0_chi2 < 20)";
      //cut2+="&& (Bp_ConsBp_0_tauminus_0_decayLength > 0.5) && (Bp_ConsBp_0_tauminus_decayLength > 0.5) && (Bp_ConsBp_0_decayLength > 5) && ( sqrt(pow(Bp_ConsBp_0_tauminus_0_PE + Bp_ConsBp_0_tauminus_PE,2) - pow(Bp_ConsBp_0_tauminus_0_PX + Bp_ConsBp_0_tauminus_PX,2) - pow(Bp_ConsBp_0_tauminus_0_PY + Bp_ConsBp_0_tauminus_PY,2) - pow(Bp_ConsBp_0_tauminus_0_PZ + Bp_ConsBp_0_tauminus_PZ,2)) < 4000) && ( sqrt(pow(Kp_PE + Bp_ConsBp_0_tauminus_0_PE,2) - pow(Kp_PX + Bp_ConsBp_0_tauminus_0_PX,2) - pow(Kp_PY + Bp_ConsBp_0_tauminus_0_PY,2) - pow(Kp_PZ + Bp_ConsBp_0_tauminus_0_PZ,2)) < 5500)";

      TString variables[] = {deltaBpPx, deltaBpPy, deltaBpPz, tau_separation, theta, "Bp_ConsBp_0_chi2", "Bp_ConsBp_0_nIter", "Bp_CC_ANGLE_B", "Bp_CC_MULT_B", "Bp_CC_SPT_B", "Bp_CC_VPT_B", "Bp_CC_PX_B", "Bp_CC_PY_B", "Bp_CC_PZ_B", "Bp_CC_PASYM_B", "Bp_CC_PTASYM_B", "Bp_CC_PXASYM_B", "Bp_CC_PYASYM_B", "Bp_CC_PZASYM_B", "Bp_CC_DELTAETA_B", "Bp_CC_DELTAPHI_B", "Bp_CC_PX_B", "Bp_CC_IT_B", "Bp_CC_MAXPT_Q_B", "Bp_CC_MAXPT_PT_B", "Bp_CC_MAXPT_PX_B", "Bp_CC_MAXPT_PY_B", "Bp_CC_MAXPT_PZ_B", "Bp_CC_MAXPT_PE_B", "Bp_NC_ANGLE_B", "Bp_NC_MULT_B", "Bp_NC_SPT_B", "Bp_NC_VPT_B", "Bp_NC_PX_B", "Bp_NC_PY_B", "Bp_NC_PZ_B", "Bp_NC_PASYM_B", "Bp_NC_PTASYM_B", "Bp_NC_PXASYM_B", "Bp_NC_PYASYM_B", "Bp_NC_PZASYM_B", "Bp_NC_DELTAETA_B", "Bp_NC_DELTAPHI_B", "Bp_NC_IT_B", "Bp_NC_MAXPT_PT_B", "Bp_NC_MAXPT_PX_B", "Bp_NC_MAXPT_PY_B", "Bp_NC_MAXPT_PZ_B", "Bp_CCNC_IT_B",
                  "Bp_CC_ANGLE_K", "Bp_CC_MULT_K", "Bp_CC_SPT_K", "Bp_CC_VPT_K", "Bp_CC_PX_K", "Bp_CC_PY_K", "Bp_CC_PZ_K", "Bp_CC_PASYM_K", "Bp_CC_PTASYM_K", "Bp_CC_PXASYM_K", "Bp_CC_PYASYM_K", "Bp_CC_PZASYM_K", "Bp_CC_DELTAETA_K", "Bp_CC_DELTAPHI_K", "Bp_CC_PX_K", "Bp_CC_IT_K", "Bp_CC_MAXPT_Q_K", "Bp_CC_MAXPT_PT_K", "Bp_CC_MAXPT_PX_K", "Bp_CC_MAXPT_PY_K", "Bp_CC_MAXPT_PZ_K", "Bp_CC_MAXPT_PE_K", "Bp_NC_ANGLE_K", "Bp_NC_MULT_K", "Bp_NC_SPT_K", "Bp_NC_VPT_K", "Bp_NC_PX_K", "Bp_NC_PY_K", "Bp_NC_PZ_K", "Bp_NC_PASYM_K", "Bp_NC_PTASYM_K", "Bp_NC_PXASYM_K", "Bp_NC_PYASYM_K", "Bp_NC_PZASYM_K", "Bp_NC_DELTAETA_K", "Bp_NC_DELTAPHI_K", "Bp_NC_IT_K", "Bp_NC_MAXPT_PT_K", "Bp_NC_MAXPT_PX_K", "Bp_NC_MAXPT_PY_K", "Bp_NC_MAXPT_PZ_K", "Bp_CCNC_IT_K",
                  "Bp_CC_ANGLE_tau1", "Bp_CC_MULT_tau1", "Bp_CC_SPT_tau1", "Bp_CC_VPT_tau1", "Bp_CC_PX_tau1", "Bp_CC_PY_tau1", "Bp_CC_PZ_tau1", "Bp_CC_PASYM_tau1", "Bp_CC_PTASYM_tau1", "Bp_CC_PXASYM_tau1", "Bp_CC_PYASYM_tau1", "Bp_CC_PZASYM_tau1", "Bp_CC_DELTAETA_tau1", "Bp_CC_DELTAPHI_tau1", "Bp_CC_PX_tau1", "Bp_CC_IT_tau1", "Bp_CC_MAXPT_Q_tau1", "Bp_CC_MAXPT_PT_tau1", "Bp_CC_MAXPT_PX_tau1", "Bp_CC_MAXPT_PY_tau1", "Bp_CC_MAXPT_PZ_tau1", "Bp_CC_MAXPT_PE_tau1", "Bp_NC_ANGLE_tau1", "Bp_NC_MULT_tau1", "Bp_NC_SPT_tau1", "Bp_NC_VPT_tau1", "Bp_NC_PX_tau1", "Bp_NC_PY_tau1", "Bp_NC_PZ_tau1", "Bp_NC_PASYM_tau1", "Bp_NC_PTASYM_tau1", "Bp_NC_PXASYM_tau1", "Bp_NC_PYASYM_tau1", "Bp_NC_PZASYM_tau1", "Bp_NC_DELTAETA_tau1", "Bp_NC_DELTAPHI_tau1", "Bp_NC_IT_tau1", "Bp_NC_MAXPT_PT_tau1", "Bp_NC_MAXPT_PX_tau1", "Bp_NC_MAXPT_PY_tau1", "Bp_NC_MAXPT_PZ_tau1", "Bp_CCNC_IT_tau1",
                  "Bp_CC_ANGLE_tau2", "Bp_CC_MULT_tau2", "Bp_CC_SPT_tau2", "Bp_CC_VPT_tau2", "Bp_CC_PX_tau2", "Bp_CC_PY_tau2", "Bp_CC_PZ_tau2", "Bp_CC_PASYM_tau2", "Bp_CC_PTASYM_tau2", "Bp_CC_PXASYM_tau2", "Bp_CC_PYASYM_tau2", "Bp_CC_PZASYM_tau2", "Bp_CC_DELTAETA_tau2", "Bp_CC_DELTAPHI_tau2", "Bp_CC_PX_tau2", "Bp_CC_IT_tau2", "Bp_CC_MAXPT_Q_tau2", "Bp_CC_MAXPT_PT_tau2", "Bp_CC_MAXPT_PX_tau2", "Bp_CC_MAXPT_PY_tau2", "Bp_CC_MAXPT_PZ_tau2", "Bp_CC_MAXPT_PE_tau2", "Bp_NC_ANGLE_tau2", "Bp_NC_MULT_tau2", "Bp_NC_SPT_tau2", "Bp_NC_VPT_tau2", "Bp_NC_PX_tau2", "Bp_NC_PY_tau2", "Bp_NC_PZ_tau2", "Bp_NC_PASYM_tau2", "Bp_NC_PTASYM_tau2", "Bp_NC_PXASYM_tau2", "Bp_NC_PYASYM_tau2", "Bp_NC_PZASYM_tau2", "Bp_NC_DELTAETA_tau2", "Bp_NC_DELTAPHI_tau2", "Bp_NC_IT_tau2", "Bp_NC_MAXPT_PT_tau2", "Bp_NC_MAXPT_PX_tau2", "Bp_NC_MAXPT_PY_tau2", "Bp_NC_MAXPT_PZ_tau2", "Bp_CCNC_IT_tau2",
                  "Bp_ConsBp_0_Kplus_PE", "Bp_ConsBp_0_Kplus_PX", "Bp_ConsBp_0_Kplus_PY", "Bp_ConsBp_0_Kplus_PZ", "Bp_ConsBp_0_M", "Bp_ConsBp_0_MERR", "Bp_ConsBp_0_PE", "Bp_ConsBp_0_PX", "Bp_ConsBp_0_PY", "Bp_ConsBp_0_PZ", "Bp_ConsBp_0_PV_X", "Bp_ConsBp_0_PV_Y", "Bp_ConsBp_0_PV_Z", "Bp_ConsBp_0_ctau", "Bp_ConsBp_0_ctauErr", "Bp_ConsBp_0_decayLength", "Bp_ConsBp_0_decayLengthErr", 
                  "Bp_ConsBp_0_tauminus_0_M", "Bp_ConsBp_0_tauminus_0_MERR", "Bp_ConsBp_0_tauminus_0_PE", "Bp_ConsBp_0_tauminus_0_PX", "Bp_ConsBp_0_tauminus_0_PY", "Bp_ConsBp_0_tauminus_0_PZ", "Bp_ConsBp_0_tauminus_0_ctau", "Bp_ConsBp_0_tauminus_0_ctauErr", "Bp_ConsBp_0_tauminus_0_decayLength", "Bp_ConsBp_0_tauminus_0_decayLengthErr",
                  "Bp_ConsBp_0_tauminus_0_nu_tau_PE", "Bp_ConsBp_0_tauminus_0_nu_tau_PX", "Bp_ConsBp_0_tauminus_0_nu_tau_PY", "Bp_ConsBp_0_tauminus_0_nu_tau_PZ", "Bp_ConsBp_0_tauminus_0_piplus_0_PE", "Bp_ConsBp_0_tauminus_0_piplus_0_PX", "Bp_ConsBp_0_tauminus_0_piplus_0_PY", "Bp_ConsBp_0_tauminus_0_piplus_0_PZ", "Bp_ConsBp_0_tauminus_0_piplus_1_PE", "Bp_ConsBp_0_tauminus_0_piplus_1_PX", "Bp_ConsBp_0_tauminus_0_piplus_1_PY", "Bp_ConsBp_0_tauminus_0_piplus_1_PZ", "Bp_ConsBp_0_tauminus_0_piplus_PE", "Bp_ConsBp_0_tauminus_0_piplus_PX", "Bp_ConsBp_0_tauminus_0_piplus_PY", "Bp_ConsBp_0_tauminus_0_piplus_PZ",  
                  "Bp_ConsBp_0_tauminus_M", "Bp_ConsBp_0_tauminus_MERR", "Bp_ConsBp_0_tauminus_PE", "Bp_ConsBp_0_tauminus_PX", "Bp_ConsBp_0_tauminus_PY", "Bp_ConsBp_0_tauminus_PZ", "Bp_ConsBp_0_tauminus_ctau", "Bp_ConsBp_0_tauminus_ctauErr", "Bp_ConsBp_0_tauminus_decayLength", "Bp_ConsBp_0_tauminus_decayLengthErr",
                  "Bp_ConsBp_0_tauminus_nu_tau_PE", "Bp_ConsBp_0_tauminus_nu_tau_PX", "Bp_ConsBp_0_tauminus_nu_tau_PY", "Bp_ConsBp_0_tauminus_nu_tau_PZ", "Bp_ConsBp_0_tauminus_piplus_0_PE", "Bp_ConsBp_0_tauminus_piplus_0_PX", "Bp_ConsBp_0_tauminus_piplus_0_PY", "Bp_ConsBp_0_tauminus_piplus_0_PZ", "Bp_ConsBp_0_tauminus_piplus_1_PE", "Bp_ConsBp_0_tauminus_piplus_1_PX", "Bp_ConsBp_0_tauminus_piplus_1_PY", "Bp_ConsBp_0_tauminus_piplus_1_PZ", "Bp_ConsBp_0_tauminus_piplus_PE", "Bp_ConsBp_0_tauminus_piplus_PX", "Bp_ConsBp_0_tauminus_piplus_PY", "Bp_ConsBp_0_tauminus_piplus_PZ", 
                  "Bp_DIRA_OWNPV", "Bp_ENDVERTEX_CHI2", "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_XERR", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_YERR", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_ZERR", "Bp_ETA", "Bp_FDCHI2_OWNPV", "Bp_FD_OWNPV", "Bp_M", "Bp_M01", "Bp_M02", "Bp_M03", "Bp_M04", "Bp_M05", "Bp_M06", "Bp_M012", "Bp_M013", "Bp_M014", "Bp_M015", "Bp_M016", "Bp_M023", "Bp_M024", "Bp_M025", "Bp_M026", "Bp_M034", "Bp_M035", "Bp_M036", "Bp_M045", "Bp_M046", "Bp_M056", "Bp_MCORR", "Bp_OWNPV_CHI2", "Bp_OWNPV_X", "Bp_OWNPV_XERR", "Bp_OWNPV_Y", "Bp_OWNPV_YERR", "Bp_OWNPV_Z", "Bp_OWNPV_ZERR",
                  "Bp_PE", "Bp_PHI", "Bp_PX", "Bp_PY", "Bp_PZ", "Bp_PT", "Bp_VCHI2PDOF", 
                  "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISONUMVTX_B", 
                  "Bp_VTXISODCHI2MASSONETRACK_K", "Bp_VTXISODCHI2MASSTWOTRACK_K", "Bp_VTXISODCHI2ONETRACK_K", "Bp_VTXISODCHI2TWOTRACK_K", "Bp_VTXISONUMVTX_K", 
                  "Bp_VTXISODCHI2MASSONETRACK_tau1", "Bp_VTXISODCHI2MASSTWOTRACK_tau1", "Bp_VTXISODCHI2ONETRACK_tau1", "Bp_VTXISODCHI2TWOTRACK_tau1", "Bp_VTXISONUMVTX_tau1", 
                  "Bp_VTXISODCHI2MASSONETRACK_tau2", "Bp_VTXISODCHI2MASSTWOTRACK_tau2", "Bp_VTXISODCHI2ONETRACK_tau2", "Bp_VTXISODCHI2TWOTRACK_tau2", "Bp_VTXISONUMVTX_tau2", 
                  "Kp_ETA", "Kp_M", "Kp_PE", "Kp_PHI", "Kp_PIDK", "Kp_PX", "Kp_PY", "Kp_PZ", "Kp_PT",
                  "taum_DIRA_ORIVX", "taum_DIRA_OWNPV", "taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_X", "taum_ENDVERTEX_Y", "taum_ENDVERTEX_Z", "taum_ENDVERTEX_XERR", "taum_ENDVERTEX_YERR", "taum_ENDVERTEX_ZERR", "taum_FDCHI2_ORIVX", "taum_FDCHI2_OWNPV", "taum_FD_ORIVX", "taum_FD_OWNPV", 
                  "taum_M", "taum_M12", "taum_M13", "taum_M23", "taum_PE", "taum_PX", "taum_PY", "taum_PZ", "taum_PT",
                  "taup_DIRA_ORIVX", "taup_DIRA_OWNPV", "taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_X", "taup_ENDVERTEX_Y", "taup_ENDVERTEX_Z", "taup_ENDVERTEX_XERR", "taup_ENDVERTEX_YERR", "taup_ENDVERTEX_ZERR", "taup_FDCHI2_ORIVX", "taup_FDCHI2_OWNPV", "taup_FD_ORIVX", "taup_FD_OWNPV", 
                  "taup_M", "taup_M12", "taup_M13", "taup_M23", "taup_PE", "taup_PX", "taup_PY", "taup_PZ", "taup_PT",
                  };

      // other years
      // TString variables[] =  {"Bp_BPVDIRA", "Bp_CONEDELTAETA_B", "Bp_CONEDELTAETA_K", "Bp_CONEDELTAETA_tau1", "Bp_CONEDELTAETA_tau2", "Bp_CONEMULT_B", "Bp_CONEMULT_K", "Bp_CONEMULT_tau1", "Bp_CONEMULT_tau2", "Bp_CONEPASYM_B", "Bp_CONEPASYM_K", "Bp_CONEPASYM_tau1", "Bp_CONEPASYM_tau2", "Bp_CONEPTASYM_B", "Bp_CONEPTASYM_K", "Bp_CONEPTASYM_tau1", "Bp_CONEPTASYM_tau2",
      //             "Bp_ConsBp_0_Kplus_PE", "Bp_ConsBp_0_Kplus_PX", "Bp_ConsBp_0_Kplus_PY", "Bp_ConsBp_0_Kplus_PZ", "Bp_ConsBp_0_M", "Bp_ConsBp_0_MERR", "Bp_ConsBp_0_PE", "Bp_ConsBp_0_PX", "Bp_ConsBp_0_PY", "Bp_ConsBp_0_PZ", "Bp_ConsBp_0_PV_X", "Bp_ConsBp_0_PV_Y", "Bp_ConsBp_0_PV_Z", "Bp_ConsBp_0_ctau", "Bp_ConsBp_0_ctauErr", "Bp_ConsBp_0_decayLength", "Bp_ConsBp_0_decayLengthErr", 
      //             "Bp_ConsBp_0_tauminus_0_M", "Bp_ConsBp_0_tauminus_0_MERR", "Bp_ConsBp_0_tauminus_0_PE", "Bp_ConsBp_0_tauminus_0_PX", "Bp_ConsBp_0_tauminus_0_PY", "Bp_ConsBp_0_tauminus_0_PZ", "Bp_ConsBp_0_tauminus_0_ctau", "Bp_ConsBp_0_tauminus_0_ctauErr", "Bp_ConsBp_0_tauminus_0_decayLength", "Bp_ConsBp_0_tauminus_0_decayLengthErr",
      //             "Bp_ConsBp_0_tauminus_0_nu_tau_PE", "Bp_ConsBp_0_tauminus_0_nu_tau_PX", "Bp_ConsBp_0_tauminus_0_nu_tau_PY", "Bp_ConsBp_0_tauminus_0_nu_tau_PZ", "Bp_ConsBp_0_tauminus_0_piplus_0_PE", "Bp_ConsBp_0_tauminus_0_piplus_0_PX", "Bp_ConsBp_0_tauminus_0_piplus_0_PY", "Bp_ConsBp_0_tauminus_0_piplus_0_PZ", "Bp_ConsBp_0_tauminus_0_piplus_1_PE", "Bp_ConsBp_0_tauminus_0_piplus_1_PX", "Bp_ConsBp_0_tauminus_0_piplus_1_PY", "Bp_ConsBp_0_tauminus_0_piplus_1_PZ", "Bp_ConsBp_0_tauminus_0_piplus_PE", "Bp_ConsBp_0_tauminus_0_piplus_PX", "Bp_ConsBp_0_tauminus_0_piplus_PY", "Bp_ConsBp_0_tauminus_0_piplus_PZ",  
      //             "Bp_ConsBp_0_tauminus_M", "Bp_ConsBp_0_tauminus_MERR", "Bp_ConsBp_0_tauminus_PE", "Bp_ConsBp_0_tauminus_PX", "Bp_ConsBp_0_tauminus_PY", "Bp_ConsBp_0_tauminus_PZ", "Bp_ConsBp_0_tauminus_ctau", "Bp_ConsBp_0_tauminus_ctauErr", "Bp_ConsBp_0_tauminus_decayLength", "Bp_ConsBp_0_tauminus_decayLengthErr",
      //             "Bp_ConsBp_0_tauminus_nu_tau_PE", "Bp_ConsBp_0_tauminus_nu_tau_PX", "Bp_ConsBp_0_tauminus_nu_tau_PY", "Bp_ConsBp_0_tauminus_nu_tau_PZ", "Bp_ConsBp_0_tauminus_piplus_0_PE", "Bp_ConsBp_0_tauminus_piplus_0_PX", "Bp_ConsBp_0_tauminus_piplus_0_PY", "Bp_ConsBp_0_tauminus_piplus_0_PZ", "Bp_ConsBp_0_tauminus_piplus_1_PE", "Bp_ConsBp_0_tauminus_piplus_1_PX", "Bp_ConsBp_0_tauminus_piplus_1_PY", "Bp_ConsBp_0_tauminus_piplus_1_PZ", "Bp_ConsBp_0_tauminus_piplus_PE", "Bp_ConsBp_0_tauminus_piplus_PX", "Bp_ConsBp_0_tauminus_piplus_PY", "Bp_ConsBp_0_tauminus_piplus_PZ", 
      //             "Bp_DIRA_OWNPV", "Bp_ENDVERTEX_CHI2", "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_XERR", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_YERR", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_ZERR", "Bp_ETA", "Bp_FDCHI2_OWNPV", "Bp_FD_OWNPV", "Bp_M", "Bp_M01", "Bp_M02", "Bp_M03", "Bp_M04", "Bp_M05", "Bp_M06", "Bp_M012", "Bp_M013", "Bp_M014", "Bp_M015", "Bp_M016", "Bp_M023", "Bp_M024", "Bp_M025", "Bp_M026", "Bp_M034", "Bp_M035", "Bp_M036", "Bp_M045", "Bp_M046", "Bp_M056", "Bp_MCORR", "Bp_OWNPV_CHI2", "Bp_OWNPV_X", "Bp_OWNPV_XERR", "Bp_OWNPV_Y", "Bp_OWNPV_YERR", "Bp_OWNPV_Z", "Bp_OWNPV_ZERR",
      //             "Bp_PE", "Bp_PHI", "Bp_PX", "Bp_PY", "Bp_PZ", "Bp_PT", "Bp_VCHI2PDOF", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISONUMVTX_B", 
      //             "Kp_ETA", "Kp_M", "Kp_PE", "Kp_PHI", "Kp_PIDK", "Kp_PX", "Kp_PY", "Kp_PZ", "Kp_PT",
      //             "taum_DIRA_ORIVX", "taum_DIRA_OWNPV", "taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_X", "taum_ENDVERTEX_Y", "taum_ENDVERTEX_Z", "taum_ENDVERTEX_XERR", "taum_ENDVERTEX_YERR", "taum_ENDVERTEX_ZERR", "taum_FDCHI2_ORIVX", "taum_FDCHI2_OWNPV", "taum_FD_ORIVX", "taum_FD_OWNPV", 
      //             "taum_M", "taum_M12", "taum_M13", "taum_M23", "taum_PE", "taum_PX", "taum_PY", "taum_PZ", "taum_PT",
      //             "taup_DIRA_ORIVX", "taup_DIRA_OWNPV", "taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_X", "taup_ENDVERTEX_Y", "taup_ENDVERTEX_Z", "taup_ENDVERTEX_XERR", "taup_ENDVERTEX_YERR", "taup_ENDVERTEX_ZERR", "taup_FDCHI2_ORIVX", "taup_FDCHI2_OWNPV", "taup_FD_ORIVX", "taup_FD_OWNPV", 
      //             "taup_M", "taup_M12", "taup_M13", "taup_M23", "taup_PE", "taup_PX", "taup_PY", "taup_PZ", "taup_PT",
      //             };

      TString titles[]  = {"Area of triangle PV-DV1-DV2 (mm^{2})", "Impact parameter PV to K^{+} traj. (mm)", "Impact parameter DV1 to K^{+} traj. (mm)", "Impact parameter DV2 to K^{+} traj. (mm)", "Distance btw BV and DV1 (mm)", "Distance btw BV and DV2 (mm)", "Distance btw DV1 and DV2 (mm)", "Angle btw K^{+} momentum and surface formed by PV-DV1-DV2 (^{#circ})",
                              "Bp_BPVDIRA", "B^{+} cone #Delta #eta", "K^{+} cone #Delta #eta", "#tau^{+} cone #Delta #eta", "#tau^{-} cone #Delta #eta", "B^{+} cone multiplicity", "K^{+} cone multiplicity", "#tau^{+} cone multiplicity", "#tau^{-} cone multiplicity", "B^{+} cone p asymmetry", "K^{+} cone p asymmetry", "#tau^{+} cone p asymmetry", "#tau^{-} cone p asymmetry", "B^{+} cone pT asymmetry", "K^{+} cone pT asymmetry", "#tau^{+} cone pT asymmetry", "#tau^{-} cone pT asymmetry",
                              "DTF K^{+} E (MeV)", "DTF K^{+} Px (MeV)", "DTF K^{+} Py (MeV)", "DTF K^{+} Pz (MeV)", "DTF B^{+} M (MeV)", "DTF B^{+} #Delta M (MeV)", "DTF B^{+} E (MeV)", "DTF B^{+} Px (MeV)", "DTF B^{+} Py (MeV)", "DTF B^{+} Pz (MeV)", "DTF PVx (mm)", "DTF PVy (mm)", "DTF PVz (mm)", "DTF B^{+} c#tau (mm)", "DTF B^{+} #Delta(c#tau) (mm)", "DTF B^{+} decay length (mm)", "DTF B^{+} decay length error (mm)", 
                              "DTF #tau^{-} M (MeV)", "DTF #tau^{-} #Delta M (MeV)", "DTF #tau^{-} E (MeV)", "DTF #tau^{-} Px (MeV)", "DTF #tau^{-} Py (MeV)", "DTF #tau^{-} Pz (MeV)", "DTF #tau^{-} c#tau (mm)", "DTF #tau^{-} #Delta(c#tau) (mm)", "DTF #tau^{-} decay length (mm)", "DTF #tau^{-} decay length error (mm)",
                              "DTF #nu_{#tau} E (MeV)", "DTF #nu_{#tau} Px (MeV)", "DTF #nu_{#tau} Py (MeV)", "DTF #nu_{#tau} Pz (MeV)", "DTF #tau^{-} #pi_1 E (MeV)", "DTF #tau^{-} #pi_1 Px (MeV)", "DTF #tau^{-} #pi_1 Py (MeV)", "DTF #tau^{-} #pi_1 Pz (MeV)", "DTF #tau^{-} #pi_2 E (MeV)", "DTF #tau^{-} #pi_2 Px (MeV)", "DTF #tau^{-} #pi_2 Py (MeV)", "DTF #tau^{-} #pi_2 Pz (MeV)", "DTF #tau^{-} #pi_3 E (MeV)", "DTF #tau^{-} #pi_3 Px (MeV)", "DTF #tau^{-} #pi_3 Py (MeV)", "DTF #tau^{-} #pi_3 Pz (MeV)",  
                              "DTF #tau^{+} M (MeV)", "DTF #tau^{+} #Delta M (MeV)", "DTF #tau^{+} E (MeV)", "DTF #tau^{+} Px (MeV)", "DTF #tau^{+} Py (MeV)", "DTF #tau^{+} Pz (MeV)", "DTF #tau^{+} c#tau (mm)", "DTF #tau^{+} #Delta(c#tau) (mm)", "DTF #tau^{+} decay length (mm)", "DTF #tau^{+} decay length error (mm)",
                              "DTF #bar{#nu}_{#tau} E (MeV)", "DTF #bar{#nu}_{#tau} Px (MeV)", "DTF #bar{#nu}_{#tau} Py (MeV)", "DTF #bar{#nu}_{#tau} Pz (MeV)", "DTF #tau^{+} #pi_1 E (MeV)", "DTF #tau^{+} #pi_1 Px (MeV)", "DTF #tau^{+} #pi_1 Py (MeV)", "DTF #tau^{+} #pi_1 Pz (MeV)", "DTF #tau^{+} #pi_2 E (MeV)", "DTF #tau^{+} #pi_2 Px (MeV)", "DTF #tau^{+} #pi_2 Py (MeV)", "DTF #tau^{+} #pi_2 Pz (MeV)", "DTF #tau^{+} #pi_3 E (MeV)", "DTF #tau^{+} #pi_3 Px (MeV)", "DTF #tau^{+} #pi_3 Py (MeV)", "DTF #tau^{+} #pi_3 Pz (MeV)", 
                              "B^{+} DIRA to PV", "BV #chi^{2}", "BVx (mm)", "BVx error (mm)", "BVy (mm)", "BVy error (mm)", "BVz (mm)", "BVz error (mm)", "B^{+} #eta", "B^{+} FD to PV #chi^{2}", "B^{+} FD to BV (mm)", "B^{+} visible mass (MeV)", "M(K^{+} #pi^{+}) (MeV) (#tau^{-})", "M(K^{+} #pi^{+}) (MeV) (#tau^{-})", "M(K^{+} #pi^{-}) (MeV) (#tau^{-})", "M(K^{+} #pi^{+}) (MeV) (#tau^{+})", "M(K^{+} #pi^{+}) (MeV) (#tau^{+})", "M(K^{+} #pi^{-}) (MeV) (#tau^{+})", "M(K^{+} #pi^{+} #pi^{-}) (MeV) (same #tau^{-})", "M(K^{+} #pi^{+} #pi^{-}) (MeV) (same #tau^{-})", "M(K^{+} #pi^{+} #pi^{+}) (MeV) (different #tau)", "M(K^{+} #pi^{+} #pi^{+}) (MeV) (different #tau)", "M(K^{+} #pi^{+} #pi^{-}) (MeV) (different #tau)", "M(K^{+} #pi^{-} #pi^{-}) (MeV) (same #tau^{-})", "M(K^{+} #pi^{-} #pi^{+}) (MeV) (different #tau)", "M(K^{+} #pi^{-} #pi^{+}) (MeV) (different #tau)", "M(K^{+} #pi^{-} #pi^{-}) (MeV) (different #tau)", "M(K^{+} #pi^{-} #pi^{+}) (MeV) (different #tau)", "M(K^{+} #pi^{-} #pi^{+}) (MeV) (different #tau)", "M(K^{+} #pi^{-} #pi^{-}) (MeV) (different #tau)", "M(K^{+} #pi^{+} #pi^{+}) (MeV) (same #tau^{+})", "M(K^{+} #pi^{+} #pi^{-}) (MeV) (same #tau^{+})", "M(K^{+} #pi^{+} #pi^{-}) (MeV) (same #tau^{+})", "B^{+} corrected mass (MeV)", "PV #chi^{2}", "PVx (mm)", "PVx error (mm)", "PVy (mm)", "PVy error (mm)", "PVz (mm)", "PVz error (mm)",
                              "B^{+} E (MeV)", "B^{+} #phi (rad)", "B^{+} Px (MeV)", "B^{+} Py (MeV)", "B^{+} Pz (MeV)", "B^{+} pT (MeV)", "Bp_VCHI2PDOF", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISONUMVTX_B", 
                              "K^{+} #eta", "K^{+} M (MeV)", "K^{+} E (MeV)", "K^{+} #phi", "Kp_PIDK", "K^{+} Px (MeV)", "K^{+} Py (MeV)", "K^{+} Pz (MeV)", "K^{+} pT (MeV)",
                              "#tau^{-} DIRA to BV", "#tau^{-} DIRA to PV", "#tau^{-} DV #chi^{2}", "#tau^{-} DVx (mm)", "#tau^{-} DVy (mm)", "#tau^{-} DVz (mm)", "#tau^{-} DVx error (mm)", "#tau^{-} DVy error (mm)", "#tau^{-} DVz error (mm)", "#tau^{-} FD to BV #chi^{2}", "#tau^{-} FD to PV #chi^{2}", "#tau^{-} FD to BV (mm)", "#tau^{-} FD to PV (mm)", 
                              "#tau^{-} M (MeV)", "#tau^{-} M12 (MeV)", "#tau^{-} M13 (MeV)", "#tau^{-} M23 (MeV)", "#tau^{-} E (MeV)", "#tau^{-} Px (MeV)", "#tau^{-} Py (MeV)", "#tau^{-} Pz (MeV)", "#tau^{-} pT (MeV)",
                              "#tau^{+} DIRA to BV", "#tau^{+} DIRA to PV", "#tau^{+} DV #chi^{2}", "#tau^{+} DVx (mm)", "#tau^{+} DVy (mm)", "#tau^{+} DVz (mm)", "#tau^{+} DVx error (mm)", "#tau^{+} DVy error (mm)", "#tau^{+} DVz error (mm)", "#tau^{+} FD to BV #chi^{2}", "#tau^{+} FD to PV #chi^{2}", "#tau^{+} FD to BV (mm)", "#tau^{+} FD to PV (mm)", 
                              "#tau^{+} M (MeV)", "#tau^{+} M12 (MeV)", "#tau^{+} M13 (MeV)", "#tau^{+} M23 (MeV)", "#tau^{+} E (MeV)", "#tau^{+} Px (MeV)", "#tau^{+} Py (MeV)", "#tau^{+} Pz (MeV)", "#tau^{+} pT (MeV)",
                            };


      std::vector<TCanvas*> c;
      gStyle->SetOptStat(0);

      t1->SetLineColor(kBlue);
      t1->SetFillColorAlpha(kBlue, 0.25);
      t2->SetLineColor(kRed);
      t2->SetFillColorAlpha(kRed, 0.25);

      int n_var =  sizeof(variables) / sizeof(variables[0]);

      std::vector<TH1D*> h1;
      std::vector<TH1D*> h2;

      TString var_name;

      for(int i = 0; i < n_var; i++){
            double xmin;  
            double xmax;

            if(variables[i] == "Bp_ConsBp_0_ctau"){
                  xmin = 0;
                  xmax = 50;
            }
            else if(variables[i] == "Bp_ConsBp_0_ctauErr"){
                  xmin = 0;
                  xmax = 5; 
            }
            else if(variables[i] == "Bp_ConsBp_0_decayLength"){
                  xmin = 0;
                  xmax = 250;        
            }
            else if(variables[i] == "Bp_ConsBp_0_decayLengthErr"){
                  xmin = 0;
                  xmax = 10;              
            }
            else if(variables[i] == "Bp_ConsBp_0_M"){
                  xmin = 4000;
                  xmax = 8000;  
            }
            else if(variables[i] == "Bp_ConsBp_0_MERR"){
                  xmin = 0;
                  xmax = 1000;                   
            }
            else if(variables[i] == "Bp_ConsBp_0_PE"){
                  xmin = 0;
                  xmax = 1000000;                     
            }
            else if(variables[i] == "Bp_ConsBp_0_PX"){
                  xmin = -50000;
                  xmax = 50000;
            }
            else if(variables[i] == "Bp_ConsBp_0_PY"){
                  xmin = -50000;
                  xmax = 50000;
            }
            else if(variables[i] == "Bp_ConsBp_0_PZ"){
                  xmin = 0;
                  xmax = 1000000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_ctau"){
                  xmin = 0;
                  xmax = 3;                  
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_ctauErr"){
                  xmin = 0;
                  xmax = 0.3;                    
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_decayLength"){
                  xmin = 0;
                  xmax = 80;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_decayLengthErr"){
                  xmin = 0;
                  xmax = 6;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_M"){
                  xmin = 1400;
                  xmax = 2000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_MERR"){
                  xmin = 0;
                  xmax = 100;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_nu_tau_PE"){
                  xmin = 0;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_nu_tau_PX"){
                  xmin = -5000;
                  xmax = 5000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_nu_tau_PY"){
                  xmin = -5000;
                  xmax = 5000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_nu_tau_PZ"){
                  xmin = 0;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_PE"){
                  xmin = 0;
                  xmax = 100000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_PX"){
                  xmin = -15000;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_PY"){
                  xmin = -15000;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_0_PZ"){
                  xmin = 0;
                  xmax = 100000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_ctau"){
                  xmin = -2;
                  xmax = 4;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_ctauErr"){
                  xmin = 0;
                  xmax = 0.3;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_decayLength"){
                  xmin = 0;
                  xmax = 100;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_decayLengthErr"){
                  xmin = 0;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_M"){
                  xmin = 1400;
                  xmax = 2000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_MERR"){
                  xmin = 0;
                  xmax = 100;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_nu_tau_PE"){
                  xmin = 0;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_nu_tau_PX"){
                  xmin = -5000;
                  xmax = 5000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_nu_tau_PY"){
                  xmin = -5000;
                  xmax = 5000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_nu_tau_PZ"){
                  xmin = 0;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_PE"){
                  xmin = 0;
                  xmax = 100000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_PX"){
                  xmin = -15000;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_PY"){
                  xmin = -15000;
                  xmax = 15000;
            }
            else if(variables[i] == "Bp_ConsBp_0_tauminus_PZ"){
                  xmin = 0;
                  xmax = 1000000;
            }
            else if(variables[i] == "Bp_ENDVERTEX_CHI2"){
                  xmin = 0;
                  xmax = 50;
            }
            else if(variables[i] == "Bp_FDCHI2_OWNPV"){
                  xmin = 0;
                  xmax = 150000;
            }
            else if(variables[i] == "Bp_MCORR"){
                  xmin = 2000;
                  xmax = 8000;
            }
            else if(variables[i] == "Bp_VCHI2PDOF"){
                  xmin = 0;
                  xmax = 300;
            }
            else if(variables[i] == "Bp_VTXISODCHI2MASSONETRACK_B"){
                  xmin = 0;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_VTXISODCHI2MASSTWOTRACK_B"){
                  xmin = 0;
                  xmax = 30000;
            }
            else if(variables[i] == "Bp_VTXISODCHI2ONETRACK_B"){
                  xmin = 0;
                  xmax = 5000;
            }
            else if(variables[i] == "Bp_VTXISODCHI2TWOTRACK_B"){
                  xmin = 0;
                  xmax = 10000;
            }
            else if(variables[i] == "Bp_VTXISONUMVTX_B"){
                  xmin = 0;
                  xmax = 1.05;
            }
            else if(variables[i] == "Kp_M"){
                  xmin = 490;
                  xmax = 500;
            }
            else if(variables[i] == "taum_DIRA_OWNPV"){
                  xmin = 0.98;
                  xmax = 1.02;
            }
            else if(variables[i] == "taum_FD_ORIVX"){
                  xmin = 0;
                  xmax = 50;
            }
            else if(variables[i] == "taum_FD_OWNPV"){
                  xmin = 0;
                  xmax = 150;
            }
            else if(variables[i] == "taum_FDCHI2_ORIVX"){
                  xmin = 0;
                  xmax = 100;
            }
            else if(variables[i] == "taum_FDCHI2_OWNPV"){
                  xmin = 0;
                  xmax = 60000;
            }
            else if(variables[i] == "taup_DIRA_OWNPV"){
                  xmin = 0.98;
                  xmax = 1.02;
            }
            else if(variables[i] == "taup_ENDVERTEX_XERR"){
                  xmin = 0;
                  xmax = 0.15;
            }
            else if(variables[i] == "taup_FD_ORIVX"){
                  xmin = 0;
                  xmax = 100;
            }
            else if(variables[i] == "taup_FD_OWNPV"){
                  xmin = 0;
                  xmax = 150;
            }
            else if(variables[i] == "taup_FDCHI2_ORIVX"){
                  xmin = 0;
                  xmax = 100;
            }
            else if(variables[i] == "taup_FDCHI2_OWNPV"){
                  xmin = 0;
                  xmax = 150000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PE_B"){
                  xmin = 0;
                  xmax = 1000000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PE_K"){
                  xmin = 0;
                  xmax = 1000000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PE_tau1"){
                  xmin = 0;
                  xmax = 1000000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PE_tau2"){
                  xmin = 0;
                  xmax = 1000000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PT_B"){
                  xmin = 0;
                  xmax = 50000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PT_K"){
                  xmin = 0;
                  xmax = 50000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PT_tau1"){
                  xmin = 0;
                  xmax = 50000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PT_tau2"){
                  xmin = 0;
                  xmax = 50000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PX_B"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PX_K"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PX_tau1"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PX_tau2"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PY_B"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PY_K"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PY_tau1"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PY_tau2"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PZ_B"){
                  xmin = 0;
                  xmax = 500000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PZ_K"){
                  xmin = 0;
                  xmax = 500000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PZ_tau1"){
                  xmin = 0;
                  xmax = 500000;
            }
            else if(variables[i] == "Bp_CC_MAXPT_PZ_tau2"){
                  xmin = 0;
                  xmax = 500000;
            }
            else if(variables[i] == "Bp_CC_PX_B"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PX_K"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PX_tau1"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PX_tau2"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PXASYM_B"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PXASYM_K"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PXASYM_tau1"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PXASYM_tau2"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PY_B"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PY_K"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PY_tau1"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PY_tau2"){
                  xmin = -20000;
                  xmax = 20000;
            }
            else if(variables[i] == "Bp_CC_PYASYM_B"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PYASYM_K"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PYASYM_tau1"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PYASYM_tau2"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_CC_PZ_B"){
                  xmin = 0;
                  xmax = 500000;
            }
           else if(variables[i] == "Bp_CC_PZ_K"){
                  xmin = 0;
                  xmax = 500000;
            }
           else if(variables[i] == "Bp_CC_PZ_tau1"){
                  xmin = 0;
                  xmax = 500000;
            }
           else if(variables[i] == "Bp_CC_PZ_tau2"){
                  xmin = 0;
                  xmax = 500000;
            }
           else if(variables[i] == "Bp_CC_SPT_B"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_SPT_K"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_SPT_tau1"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_SPT_tau2"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_VPT_B"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_VPT_K"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_VPT_tau1"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_CC_VPT_tau2"){
                  xmin = 0;
                  xmax = 50000;
            }
           else if(variables[i] == "Bp_NC_PXASYM_B"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PXASYM_K"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PXASYM_tau1"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PXASYM_tau2"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PYASYM_B"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PYASYM_K"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PYASYM_tau1"){
                  xmin = -10;
                  xmax = 10;
            }
           else if(variables[i] == "Bp_NC_PYASYM_tau2"){
                  xmin = -10;
                  xmax = 10;
            }
            else if(variables[i] == "Bp_ConsBp_0_chi2"){
                  xmin = 0;
                  xmax = 50;
            }
            else if(variables[i] == "Bp_ConsBp_0_nIter"){
                  xmin = 0;
                  xmax = 25;
            }
            else if(variables[i] == theta){
                  xmin = 0;
                  xmax = 0.05;
            }
            else if(variables[i] == deltaBpPx){
                  xmin = -25000;
                  xmax = 25000;
            }
            else if(variables[i] == deltaBpPy){
                  xmin = -50000;
                  xmax = 50000;
            }
            else if(variables[i] == deltaBpPz){
                  xmin = 0;
                  xmax = 2000000;
            }
            else if(variables[i] == "taum_DIRA_ORIVX"){
                  xmin = 0.9;
                  xmax = 1.0;
            }
            else if(variables[i] == "taup_DIRA_ORIVX"){
                  xmin = 0.9;
                  xmax = 1.0;
            }
            else if(variables[i] == "Bp_DIRA_OOWNPV"){
                  xmin = 0.9;
                  xmax = 1.0;
            }
            else{
                  xmin = t2->GetMinimum(variables[i]);
                  xmax = t2->GetMaximum(variables[i]);
            }

            if(variables[i] == tau_separation){var_name = "tau_separation";}
            else if(variables[i] == theta){var_name = "theta";}
            else if(variables[i] == deltaBpPx){var_name = "DTF_minus_visible_Bp_Px";}
            else if(variables[i] == deltaBpPy){var_name = "DTF_minus_visible_Bp_Py";}
            else if(variables[i] == deltaBpPz){var_name = "DTF_minus_visible_Bp_Pz";}
            else{var_name = variables[i];}

            h1.push_back(new TH1D(Form("h1_%.i",i), Form("h1_%.i",i), 100, xmin - 0.1*abs(xmin), xmax + 0.1*abs(xmax)));
            h2.push_back(new TH1D(Form("h2_%.i",i), Form("h2_%.i",i), 100, xmin - 0.1*abs(xmin), xmax + 0.1*abs(xmax)));

            t1->Draw(variables[i]+Form(" >> h1_%.i",i),cut1);
            t2->Draw(variables[i]+Form(" >> h2_%.i",i),cut2);

            c.push_back(new TCanvas());
            c[i]->cd();

            h1[i]->SetFillColorAlpha(kBlue,0.25);
            h2[i]->SetFillColorAlpha(kRed,0.25);

            h1[i]->SetTitle(var_name);
            h2[i]->SetTitle(var_name);
            h1[i]->GetXaxis()->SetTitle(var_name);
            h1[i]->GetYaxis()->SetTitle( "Normalised entries" );
            h2[i]->GetXaxis()->SetTitle(var_name);
            h2[i]->GetYaxis()->SetTitle( "Normalised entries" );

            if( (h1[i]->GetMaximum())/(h1[i]->Integral()) > (h2[i]->GetMaximum())/(h2[i]->Integral()) ){
                  h1[i]->DrawNormalized("HIST");
                  h2[i]->DrawNormalized("HIST same");
            }
            else{
                  h2[i]->DrawNormalized("HIST");
                  h1[i]->DrawNormalized("HIST same");
            }

            TLegend* leg = new TLegend(0.7, 0.6, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->AddEntry(h1[i],name1,"f");
            leg->AddEntry(h2[i],name2,"f");
            leg->Draw("same");

            c[i]->SaveAs(folder_name+"/"+var_name+".png");
            c[i]->SaveAs(folder_name+"/"+var_name+".pdf");
      }

      return 0;
}