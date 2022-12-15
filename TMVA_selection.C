#define type 2
// 0 -> taus (kinematic)
// 1 -> taus (isolation)
// 2 -> B (kinematic)
// 3 -> B (isolation)
// 4 -> tau+ (isolation)
// 5 -> tau- (isolation)

#define background_proxy 1
// 0 -> use RS data, with M > 6.5 GeV
// 1 -> use WS data

#define year 8

void TMVA_selection(){
    // Load TMVA libraries
    TMVA::Tools::Instance(); 

    // Create output file for TMVA
    TString name_type;
    if(type == 0){name_type = "taus_kinematic";}
    else if(type == 1){name_type = "taus_isolation";}
    else if(type == 2){name_type = "B_kinematic";}
    else if(type == 3){name_type = "B_isolation";}
    else if(type == 4){name_type = "tau1_isolation";}
    else if(type == 5){name_type = "tau2_isolation";}
    TString name_bkg_proxy;
    if(background_proxy == 0){name_bkg_proxy = "RS_data_right_sideband";}
    else if(background_proxy == 1){name_bkg_proxy = "WS_data";}
    TString outfileName = Form("TMVA/tmva_output_201%i_", year)+name_type+"_"+name_bkg_proxy+".root";
    TFile* outputFile = TFile::Open( outfileName, "RECREATE"); 
    outputFile->cd();

    // Create TMVA factory
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",outputFile,"V:!Silent:Color:Transformations=I:DrawProgressBar:AnalysisType=Classification"); 
    TMVA::DataLoader *dataloader = new TMVA::DataLoader(Form("tmva_dataset_201%i_",year)+name_type+"_"+name_bkg_proxy);

    // Open file and retrieve tries
    TFileCollection* fc = new TFileCollection("RS_data", "RS_data", Form("data_201%i_MagUp.txt",year));
    fc->AddFromFile(Form("data_201%i_MagDown.txt", year));

    TChain* bkg_tree;
    if(background_proxy == 0){bkg_tree = new TChain("ntuple/DecayTree");}
    else if(background_proxy == 1){bkg_tree = new TChain("ntuple_SS/DecayTree");}
    bkg_tree->AddFileInfoList((TCollection*)fc->GetList());

    TChain* sig_tree = new TChain("DecayTree");
    sig_tree->Add(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_truth_matched.root",year,year));

    // Pre-selections
    //TCut truth_matching = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
    TCut pass_DTF = "(Bp_ConsBp_0_status[0]==0)";
    TCut sig_mass = "(Bp_ConsBp_0_M[0]>4500) && (Bp_ConsBp_0_M[0]<6000)";
    TCut bkg_mass;
    if(background_proxy==0){bkg_mass = "(Bp_ConsBp_0_M>6500)";}
    else if(background_proxy==1){bkg_mass = "";}
    TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
    TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
    TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
    TCut trigger_lines = L0_trigger+HLT1_trigger+HLT2_trigger;
    TCut common_selections = "(Bp_ConsBp_0_decayLength[0]>5) && (Bp_ConsBp_0_tauminus_0_decayLength[0]>0.5) && (Bp_ConsBp_0_tauminus_decayLength[0]>0.5) && (taup_DIRA_ORIVX > 0.99) && (taum_DIRA_ORIVX > 0.99) && (Bp_DIRA_OWNPV > 0.99) && (sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) > 0.5) && (Bp_ConsBp_0_chi2 < 20)";

    TCut mycuts = pass_DTF+sig_mass+trigger_lines+common_selections;
    TCut mycutb = pass_DTF+bkg_mass+trigger_lines+common_selections;

    // Add variables to the training
    if(type == 0){ // taus kinematic
        dataloader->AddVariable("taum_M", "taum_M", "MeV", 'F', 800, 1600);
        dataloader->AddVariable("taum_M12", "taum_M12", "MeV", 'F', 200, 1400);
        dataloader->AddVariable("taum_M13", "taum_M13", "MeV", 'F', 200, 1400);
        dataloader->AddVariable("taum_M23", "taum_M23", "MeV", 'F', 200, 1400);
        dataloader->AddVariable("taup_M", "taup_M", "MeV", 'F', 800, 1600);
        dataloader->AddVariable("taup_M12", "taup_M12", "MeV", 'F', 200, 1400);
        dataloader->AddVariable("taup_M13", "taup_M13", "MeV", 'F', 200, 1400);
        dataloader->AddVariable("taup_M23", "taup_M23", "MeV", 'F', 200, 1400);  
        dataloader->AddVariable("taum_ENDVERTEX_CHI2", "taum_ENDVERTEX_CHI2", "", 'F', 0, 20);
        dataloader->AddVariable("taup_ENDVERTEX_CHI2", "taup_ENDVERTEX_CHI2", "", 'F', 0, 20);
        dataloader->AddVariable("Bp_ConsBp_0_tauminus_0_decayLength[0]", "DTF_taum_decayLength", "mm", 'F', 0, 40);
        dataloader->AddVariable("Bp_ConsBp_0_tauminus_0_decayLengthErr[0]", "DTF_taum_decayLengthError", "mm", 'F', 0, 6);
        dataloader->AddVariable("Bp_ConsBp_0_tauminus_decayLength[0]", "DTF_taup_decayLength", "mm", 'F', 0, 40);
        dataloader->AddVariable("Bp_ConsBp_0_tauminus_decayLengthErr[0]", "DTF_taup_decayLengthError", "mm", 'F', 0, 6);
        dataloader->AddVariable("log(1-taum_DIRA_ORIVX)", "log(1-taum_DIRA_ORIVX)", "", 'F', -20, -2);
        dataloader->AddVariable("log(1-taum_DIRA_OWNPV)", "log(1-taum_DIRA_OWNPV)", "", 'F', -16, -5);
        dataloader->AddVariable("log(1-taup_DIRA_ORIVX)", "log(1-taup_DIRA_ORIVX)", "", 'F', -20, -2);
        dataloader->AddVariable("log(1-taup_DIRA_OWNPV)", "log(1-taup_DIRA_OWNPV)", "", 'F', -16, -5);
        dataloader->AddVariable("min(taum_pi1_IPCHI2_OWNPV, min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))", "min(taum_pi1_IPCHI2_OWNPV, min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))", "", 'F', 0, 3000);
        dataloader->AddVariable("min(taup_pi1_IPCHI2_OWNPV, min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))", "min(taup_pi1_IPCHI2_OWNPV, min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))", "", 'F', 0, 3000);
        TString E = "Bp_ConsBp_0_tauminus_0_PE[0] + Bp_ConsBp_0_tauminus_PE[0]";
        TString Px = "Bp_ConsBp_0_tauminus_0_PX[0] + Bp_ConsBp_0_tauminus_PX[0]";
        TString Py = "Bp_ConsBp_0_tauminus_0_PY[0] + Bp_ConsBp_0_tauminus_PY[0]";
        TString Pz = "Bp_ConsBp_0_tauminus_0_PZ[0] + Bp_ConsBp_0_tauminus_PZ[0]";
        TString q2 = "abs(pow("+E+",2) - pow("+Px+",2) - pow("+Py+",2) - pow("+Pz+",2))";
        dataloader->AddVariable(q2, "q^{2}", "MeV^{2}", 'F', 0, 100000000);
    }
    else if((type == 1) && (year == 6)){ // taus isolation, year 2016
        dataloader->AddVariable("Bp_CONEDELTAETA_tau1", "CONEDELTAETA_tau1", "", 'F');
        dataloader->AddVariable("Bp_CONEMULT_tau1", "CONEMULT_tau1", "", 'I');
        dataloader->AddVariable("Bp_CONEPASYM_tau1", "CONEPASYM_tau1", "", 'F');
        dataloader->AddVariable("Bp_CONEPTASYM_tau1", "CONEPTASYM_tau1", "", 'F');
        dataloader->AddVariable("Bp_CONEDELTAETA_tau2", "CONEDELTAETA_tau2", "", 'F');
        dataloader->AddVariable("Bp_CONEMULT_tau2", "CONEMULT_tau2", "", 'I');
        dataloader->AddVariable("Bp_CONEPASYM_tau2", "CONEPASYM_tau2", "", 'F');
        dataloader->AddVariable("Bp_CONEPTASYM_tau2", "CONEPTASYM_tau2", "", 'F');
    }
    else if((type == 1) && (year == 8)){ // taus isolation, year 2018
        dataloader->AddVariable("Bp_CC_IT_tau1", "CC_IT_tau1", "", 'F', 0.1, 1);
        dataloader->AddVariable("Bp_CC_MAXPT_PE_tau1", "CC_MAXPT_PE_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_PT_tau1", "CC_MAXPT_PT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_Q_tau1", "CC_MAXPT_Q_tau1", "", 'I');
        dataloader->AddVariable("Bp_CC_MULT_tau1", "CC_MULT_tau1", "", 'I');
        dataloader->AddVariable("Bp_CC_PZ_tau1", "CC_PZ_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_PZASYM_tau1", "CC_PZASYM_tau1", "", 'F');
        dataloader->AddVariable("Bp_CC_SPT_tau1", "CC_SPT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_VPT_tau1", "CC_VPT_tau1", "MeV", 'F');

        dataloader->AddVariable("Bp_NC_DELTAETA_tau1", "NC_DELTAETA_tau1", "", 'F');
        dataloader->AddVariable("Bp_NC_DELTAPHI_tau1", "NC_DELTAPHI_tau1", "", 'F');
        dataloader->AddVariable("Bp_NC_IT_tau1", "NC_IT_tau1", "", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PT_tau1", "NC_MAXPT_PT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PZ_tau1", "NC_MAXPT_PZ_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MULT_tau1", "NC_MULT_tau1", "", 'I');
        dataloader->AddVariable("Bp_NC_PTASYM_tau1", "NC_PTASYM_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_PZ_tau1", "NC_PZ_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_SPT_tau1", "NC_SPT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_VPT_tau1", "NC_VPT_tau1", "MeV", 'F');

        dataloader->AddVariable("Bp_VTXISODCHI2ONETRACK_tau1", "VTXISODCHI2ONETRACK_tau1", "", 'F');
        dataloader->AddVariable("Bp_VTXISODCHI2TWOTRACK_tau1", "VTXISODCHI2TWOTRACK_tau1", "", 'F');
        dataloader->AddVariable("Bp_VTXISONUMVTX_tau1", "VTXISONUMVTX_tau1", "", 'I');

        dataloader->AddVariable("Bp_CC_IT_tau2", "CC_IT_tau2", "", 'F', 0.1, 1);
        dataloader->AddVariable("Bp_CC_MAXPT_PE_tau2", "CC_MAXPT_PE_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_PT_tau2", "CC_MAXPT_PT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_Q_tau2", "CC_MAXPT_Q_tau2", "", 'I');
        dataloader->AddVariable("Bp_CC_MULT_tau2", "CC_MULT_tau2", "", 'I');
        dataloader->AddVariable("Bp_CC_PZ_tau2", "CC_PZ_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_PZASYM_tau2", "CC_PZASYM_tau2", "", 'F');
        dataloader->AddVariable("Bp_CC_SPT_tau2", "CC_SPT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_VPT_tau2", "CC_VPT_tau2", "MeV", 'F');

        dataloader->AddVariable("Bp_NC_DELTAETA_tau2", "NC_DELTAETA_tau2", "", 'F');
        dataloader->AddVariable("Bp_NC_DELTAPHI_tau2", "NC_DELTAPHI_tau2", "", 'F');
        dataloader->AddVariable("Bp_NC_IT_tau2", "NC_IT_tau2", "", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PT_tau2", "NC_MAXPT_PT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PZ_tau2", "NC_MAXPT_PZ_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MULT_tau2", "NC_MULT_tau2", "", 'I');
        dataloader->AddVariable("Bp_NC_PTASYM_tau2", "NC_PTASYM_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_PZ_tau2", "NC_PZ_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_SPT_tau2", "NC_SPT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_VPT_tau2", "NC_VPT_tau2", "MeV", 'F');

        dataloader->AddVariable("Bp_VTXISODCHI2ONETRACK_tau2", "VTXISODCHI2ONETRACK_tau2", "", 'F');
        dataloader->AddVariable("Bp_VTXISODCHI2TWOTRACK_tau2", "VTXISODCHI2TWOTRACK_tau2", "", 'F');
        dataloader->AddVariable("Bp_VTXISONUMVTX_tau2", "VTXISONUMVTX_tau2", "", 'I');
    }
    else if(type==2){ // B kinematic
        dataloader->AddVariable("Bp_ConsBp_0_PX[0]", "DTF_Bp_PX", "MeV", 'F');
        dataloader->AddVariable("Bp_ConsBp_0_PY[0]", "DTF_Bp_PY", "MeV", 'F');
        dataloader->AddVariable("Bp_ConsBp_0_PZ[0]", "DTF_Bp_PZ", "MeV", 'F');
        dataloader->AddVariable("Bp_ConsBp_0_decayLengthErr[0]", "DTF_Bp_decayLength", "mm", 'F');
        dataloader->AddVariable("Bp_ENDVERTEX_CHI2", "Bp_ENDVERTEX_CHI2", "mm", 'F', 0, 100);
        dataloader->AddVariable("Bp_ENDVERTEX_X", "BVx", "mm", 'F', -8, 8);
        dataloader->AddVariable("Bp_ENDVERTEX_XERR", "BVx_error", "mm", 'F', 0, 0.15);
        dataloader->AddVariable("Bp_ENDVERTEX_Y", "BVy", "mm", 'F', -8, 8);
        dataloader->AddVariable("Bp_ENDVERTEX_YERR", "BVy_error", "mm", 'F', 0, 0.15);
        dataloader->AddVariable("Bp_ENDVERTEX_ZERR", "BVz_errpr", "mm", 'F', 0, 5);
        dataloader->AddVariable("Bp_ETA", "Bp_eta", "", 'F', 0, 5);
        dataloader->AddVariable("Bp_FD_OWNPV", "Bp_FD_OWNPV", "mm", 'F', 0, 80);
        dataloader->AddVariable("Bp_M01", "Bp_M01", "MeV", 'F', 0, 3000);
        dataloader->AddVariable("Bp_M02", "Bp_M02", "MeV", 'F', 0, 3000);
        dataloader->AddVariable("Bp_M03", "Bp_M03", "MeV", 'F', 0, 3000);
        dataloader->AddVariable("Bp_M04", "Bp_M04", "MeV", 'F', 0, 3000);
        dataloader->AddVariable("Bp_M05", "Bp_M05", "MeV", 'F', 0, 3000);
        dataloader->AddVariable("Bp_M06", "Bp_M06", "MeV", 'F', 0, 3000);
        dataloader->AddVariable("Bp_M012", "Bp_M012", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M013", "Bp_M013", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M014", "Bp_M014", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M015", "Bp_M015", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_M016", "Bp_M016", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_M023", "Bp_M023", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M024", "Bp_M024", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M025", "Bp_M025", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M026", "Bp_M026", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M034", "Bp_M034", "MeV", 'F', 0, 3500);
        dataloader->AddVariable("Bp_M035", "Bp_M035", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_M036", "Bp_M036", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_M045", "Bp_M045", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_M046", "Bp_M046", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_M056", "Bp_M056", "MeV", 'F', 0, 4000);
        dataloader->AddVariable("Bp_PT", "Bp_PT", "MeV", 'F', 0, 40000);
    }
    else if((type==4) && (year==8)){ // tau+ isolation
        dataloader->AddVariable("Bp_CC_IT_tau1", "CC_IT_tau1", "", 'F', 0.1, 1);
        dataloader->AddVariable("Bp_CC_MAXPT_PE_tau1", "CC_MAXPT_PE_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_PT_tau1", "CC_MAXPT_PT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_Q_tau1", "CC_MAXPT_Q_tau1", "", 'I');
        dataloader->AddVariable("Bp_CC_MULT_tau1", "CC_MULT_tau1", "", 'I');
        dataloader->AddVariable("Bp_CC_PZ_tau1", "CC_PZ_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_PZASYM_tau1", "CC_PZASYM_tau1", "", 'F');
        dataloader->AddVariable("Bp_CC_SPT_tau1", "CC_SPT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_VPT_tau1", "CC_VPT_tau1", "MeV", 'F');

        dataloader->AddVariable("Bp_NC_DELTAETA_tau1", "NC_DELTAETA_tau1", "", 'F');
        dataloader->AddVariable("Bp_NC_DELTAPHI_tau1", "NC_DELTAPHI_tau1", "", 'F');
        dataloader->AddVariable("Bp_NC_IT_tau1", "NC_IT_tau1", "", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PT_tau1", "NC_MAXPT_PT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PZ_tau1", "NC_MAXPT_PZ_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MULT_tau1", "NC_MULT_tau1", "", 'I');
        dataloader->AddVariable("Bp_NC_PTASYM_tau1", "NC_PTASYM_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_PZ_tau1", "NC_PZ_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_SPT_tau1", "NC_SPT_tau1", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_VPT_tau1", "NC_VPT_tau1", "MeV", 'F');

        dataloader->AddVariable("Bp_VTXISODCHI2ONETRACK_tau1", "VTXISODCHI2ONETRACK_tau1", "", 'F');
        dataloader->AddVariable("Bp_VTXISODCHI2TWOTRACK_tau1", "VTXISODCHI2TWOTRACK_tau1", "", 'F');
        dataloader->AddVariable("Bp_VTXISONUMVTX_tau1", "VTXISONUMVTX_tau1", "", 'I');
    }
    else if((type==5) && (year==8)){ // tau- isolation
        dataloader->AddVariable("Bp_CC_IT_tau2", "CC_IT_tau2", "", 'F', 0.1, 1);
        dataloader->AddVariable("Bp_CC_MAXPT_PE_tau2", "CC_MAXPT_PE_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_PT_tau2", "CC_MAXPT_PT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_MAXPT_Q_tau2", "CC_MAXPT_Q_tau2", "", 'I');
        dataloader->AddVariable("Bp_CC_MULT_tau2", "CC_MULT_tau2", "", 'I');
        dataloader->AddVariable("Bp_CC_PZ_tau2", "CC_PZ_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_PZASYM_tau2", "CC_PZASYM_tau2", "", 'F');
        dataloader->AddVariable("Bp_CC_SPT_tau2", "CC_SPT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_CC_VPT_tau2", "CC_VPT_tau2", "MeV", 'F');

        dataloader->AddVariable("Bp_NC_DELTAETA_tau2", "NC_DELTAETA_tau2", "", 'F');
        dataloader->AddVariable("Bp_NC_DELTAPHI_tau2", "NC_DELTAPHI_tau2", "", 'F');
        dataloader->AddVariable("Bp_NC_IT_tau2", "NC_IT_tau2", "", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PT_tau2", "NC_MAXPT_PT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MAXPT_PZ_tau2", "NC_MAXPT_PZ_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_MULT_tau2", "NC_MULT_tau2", "", 'I');
        dataloader->AddVariable("Bp_NC_PTASYM_tau2", "NC_PTASYM_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_PZ_tau2", "NC_PZ_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_SPT_tau2", "NC_SPT_tau2", "MeV", 'F');
        dataloader->AddVariable("Bp_NC_VPT_tau2", "NC_VPT_tau2", "MeV", 'F');

        dataloader->AddVariable("Bp_VTXISODCHI2ONETRACK_tau2", "VTXISODCHI2ONETRACK_tau2", "", 'F');
        dataloader->AddVariable("Bp_VTXISODCHI2TWOTRACK_tau2", "VTXISODCHI2TWOTRACK_tau2", "", 'F');
        dataloader->AddVariable("Bp_VTXISONUMVTX_tau2", "VTXISONUMVTX_tau2", "", 'I');
    }

    // Set the event weights per tree
    double sigWeight = 1.0;
    double bkgWeight = 1.0;

    // Add signal and background trees to the factory
    dataloader->AddSignalTree(sig_tree,sigWeight);
    dataloader->AddBackgroundTree(bkg_tree,bkgWeight);

    dataloader->AddCut(TCut("(abs(Bp_ConsBp_0_PX[0]) < 50000) && (abs(Bp_ConsBp_0_PY[0]) < 50000) && (Bp_ConsBp_0_PZ[0] > 0) && (Bp_ConsBp_0_PZ[0] < 800000) && (Bp_ConsBp_0_decayLengthErr[0] < 8) && (Bp_ConsBp_0_decayLengthErr[0] > 0)"), "Signal");
    dataloader->AddCut(TCut("(abs(Bp_ConsBp_0_PX[0]) < 50000) && (abs(Bp_ConsBp_0_PY[0]) < 50000) && (Bp_ConsBp_0_PZ[0] > 0) && (Bp_ConsBp_0_PZ[0] < 800000) && (Bp_ConsBp_0_decayLengthErr[0] < 8) && (Bp_ConsBp_0_decayLengthErr[0] > 0)"), "Background");

    // Prepare factory for training with cuts definned for signal and background; Samples will be split: half for training and half for testing; Randomly split
    dataloader->PrepareTrainingAndTestTree(mycuts,mycutb,"random"); 

    // Define TMVA method we want to use (BDT, neural network...); Specify options for the method (#trees to train, maximum depth of each tree in BDT)
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "NTrees=400:MaxDepth=2");
    //factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN", "");

    // Train and test the TMVA method
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    // Once this is done, a .root file is created with the results of the training. We can see:
    // - performance of each variable (signal vs backrgound)
    // - BDT responde for signal and background after training 
    // If you look at the overtraining test, then you can see the BDT response of the test and training samples superimposed.
    // If they do not match, then the BDT has been overtrained, which usually means you don't have enough events in you training samples, and the BDT will not perform as well as you expect. 
    outputFile->Close();

    cout << "==> Wrote root file: " << outputFile->GetName() << endl;
    cout << "==> TMVAClassification is done!" << endl;

    delete factory;
    delete dataloader;

    if (!gROOT->IsBatch()){TMVA::TMVAGui( outfileName );} 
    // Use this script in order to run the various individual macros
    // that plot the output of TMVA (e.g. running TMVAClassification.C),
    // stored in the file "TMVA.root"
}