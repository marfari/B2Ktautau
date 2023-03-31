// (tmva_type)

void TMVA_training(int year, int tmva_type, int background_proxy, TString DATA_files, TString MC_truthMatched_files){
    // Load TMVA libraries
    TMVA::Tools::Instance(); 

    // Create output file for TMVA
    TString outfileName = Form("/panfs/felician/B2Ktautau/workflow/TMVA_training/201%i/tmva_output_type%i_bkgProxy%i.root", year, tmva_type, background_proxy);
    TFile* outputFile = TFile::Open( outfileName, "RECREATE"); 
    outputFile->cd();

    // Create TMVA factory
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",outputFile,"V:!Silent:Color:Transformations=I:DrawProgressBar:AnalysisType=Classification"); 
    TMVA::DataLoader *dataloader = new TMVA::DataLoader(Form("tmva_dataset_type%i_bkgProxy%i",tmva_type,background_proxy));

    // Open file and retrieve data trees (only uses the 1 data grid files for training to save time/space)
    TFileCollection* fc_data = new TFileCollection("data", "data", DATA_files, 1); // takes only the 1st grid file
    TChain* bkg_tree;
    if(background_proxy == 0){bkg_tree = new TChain("ntuple/DecayTree");}
    else if(background_proxy == 1){bkg_tree = new TChain("ntuple_SS/DecayTree");}
    bkg_tree->AddFileInfoList((TCollection*)fc_data->GetList());

    // Open file and retrieve MC trees (truth-matched, 3pi 3pi)
    TFileCollection* fc_mc = new TFileCollection("MC", "MC", MC_truthMatched_files);
    TChain* sig_tree = new TChain("DecayTree");
    sig_tree->AddFileInfoList((TCollection*)fc_mc->GetList());

    // Pre-selections
    //TCut truth_matching = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
    TCut pass_DTF = "(Bp_ConsBp_seq_12_status[0]==0)";
    TCut sig_mass = "(Bp_ConsBp_seq_12_M[0]>4700) && (Bp_ConsBp_seq_12_M[0]<5800)";
    TCut bkg_mass;
    if(background_proxy==0){bkg_mass = "(Bp_ConsBp_seq_12_M[0]>6500)";}
    else if(background_proxy==1){bkg_mass = "(Bp_ConsBp_seq_12_M[0]>4700) && (Bp_ConsBp_seq_12_M[0]<5800)";}
    TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
    TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
    TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
    TCut trigger_lines = L0_trigger+HLT1_trigger+HLT2_trigger;
    TCut common_selections = "(Bp_ConsBp_seq_12_decayLength[0]>5) && (Bp_DIRA_OWNPV>0.99) && (max(taum_M, taup_M) < 1600) && (min(taum_M, taup_M) > 800) && (max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1) < 6)";

    TCut mycuts = pass_DTF+sig_mass+trigger_lines+common_selections;
    TCut mycutb = pass_DTF+bkg_mass+trigger_lines+common_selections;

    cout << "# signal events after pre-selections =  " << sig_tree->GetEntries(mycuts) << endl;
    cout << "# background events after pre-selections = " << bkg_tree->GetEntries(mycutb) << endl;

    // Add variables to the training
    if(tmva_type == 0){ // taus kinematic
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
        dataloader->AddVariable("Bp_ConsBp_seq_12_tauminus_0_decayLength[0]", "DTF_taum_decayLength", "mm", 'F', 0, 40);
        dataloader->AddVariable("Bp_ConsBp_seq_12_tauminus_0_decayLengthErr[0]", "DTF_taum_decayLengthError", "mm", 'F', 0, 6);
        dataloader->AddVariable("Bp_ConsBp_seq_12_tauminus_decayLength[0]", "DTF_taup_decayLength", "mm", 'F', 0, 40);
        dataloader->AddVariable("Bp_ConsBp_seq_12_tauminus_decayLengthErr[0]", "DTF_taup_decayLengthError", "mm", 'F', 0, 6);
        dataloader->AddVariable("log(1-taum_DIRA_ORIVX)", "log(1-taum_DIRA_ORIVX)", "", 'F', -20, -2);
        dataloader->AddVariable("log(1-taum_DIRA_OWNPV)", "log(1-taum_DIRA_OWNPV)", "", 'F', -16, -5);
        dataloader->AddVariable("log(1-taup_DIRA_ORIVX)", "log(1-taup_DIRA_ORIVX)", "", 'F', -20, -2);
        dataloader->AddVariable("log(1-taup_DIRA_OWNPV)", "log(1-taup_DIRA_OWNPV)", "", 'F', -16, -5);
        dataloader->AddVariable("min(taum_pi1_IPCHI2_OWNPV, min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))", "min(taum_pi1_IPCHI2_OWNPV, min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))", "", 'F', 0, 3000);
        dataloader->AddVariable("min(taup_pi1_IPCHI2_OWNPV, min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))", "min(taup_pi1_IPCHI2_OWNPV, min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))", "", 'F', 0, 3000);
        TString E = "Bp_ConsBp_seq_12_tauminus_0_PE[0] + Bp_ConsBp_seq_12_tauminus_PE[0]";
        TString Px = "Bp_ConsBp_seq_12_tauminus_0_PX[0] + Bp_ConsBp_seq_12_tauminus_PX[0]";
        TString Py = "Bp_ConsBp_seq_12_tauminus_0_PY[0] + Bp_ConsBp_seq_12_tauminus_PY[0]";
        TString Pz = "Bp_ConsBp_seq_12_tauminus_0_PZ[0] + Bp_ConsBp_seq_12_tauminus_PZ[0]";
        TString q2 = "abs(pow("+E+",2) - pow("+Px+",2) - pow("+Py+",2) - pow("+Pz+",2))";
        dataloader->AddVariable(q2, "q^{2}", "MeV^{2}", 'F', 0, 100000000);
    }
    else if((tmva_type == 1) && (year == 6)){ // taus isolation, year 2016
        dataloader->AddVariable("Bp_CONEDELTAETA_tau1", "CONEDELTAETA_tau1", "", 'F');
        dataloader->AddVariable("Bp_CONEMULT_tau1", "CONEMULT_tau1", "", 'I');
        dataloader->AddVariable("Bp_CONEPASYM_tau1", "CONEPASYM_tau1", "", 'F');
        dataloader->AddVariable("Bp_CONEPTASYM_tau1", "CONEPTASYM_tau1", "", 'F');
        dataloader->AddVariable("Bp_CONEDELTAETA_tau2", "CONEDELTAETA_tau2", "", 'F');
        dataloader->AddVariable("Bp_CONEMULT_tau2", "CONEMULT_tau2", "", 'I');
        dataloader->AddVariable("Bp_CONEPASYM_tau2", "CONEPASYM_tau2", "", 'F');
        dataloader->AddVariable("Bp_CONEPTASYM_tau2", "CONEPTASYM_tau2", "", 'F');
    }
    else if((tmva_type == 1) && (year == 8)){ // taus isolation, year 2018
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
    else if(tmva_type==2){ // B kinematic
        // DTF
        dataloader->AddVariable("Bp_ConsBp_seq_12_MERR[0]", "DTF_Bp_MERR", "MeV", 'F');
        dataloader->AddVariable("Bp_ConsBp_seq_12_PX[0]", "DTF_Bp_PX", "MeV", 'F');
        dataloader->AddVariable("Bp_ConsBp_seq_12_PZ[0]", "DTF_Bp_PZ", "MeV", 'F');
        dataloader->AddVariable("Bp_ConsBp_seq_12_decayLength[0]/Bp_ConsBp_seq_12_decayLengthErr[0]", "DTF_Bp_decayLength_sig", "", 'F');
        // Visible
        dataloader->AddVariable("Bp_ENDVERTEX_X", "BVx", "mm", 'F', -8, 8);
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
    else if((tmva_type==4) && (year==8)){ // tau+ isolation
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
    else if((tmva_type==5) && (year==8)){ // tau- isolation
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
    else if(tmva_type == 6) // mix
    {
        // tau (kinematics)
        dataloader->AddVariable("max(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", "max(DVm_chi2, DVp_chi2)", "", 'F', 0, 25);
        dataloader->AddVariable("min(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", "min(DVm_chi2, DVp_chi2)", "", 'F');
        dataloader->AddVariable("max(taum_M, taup_M)", "max(taum_M, taup_M)", "MeV", 'F');
        dataloader->AddVariable("max(taum_M12, taup_M12)", "max(taum_M12, taup_M12)", "MeV", 'F');
        dataloader->AddVariable("max(taum_M23, taup_M23)", "max(taum_M23, taup_M23)", "MeV", 'F');
        dataloader->AddVariable("min(taum_M, taup_M)", "min(taum_M, taup_M)", "MeV", 'F');
        dataloader->AddVariable("min(taum_M12, taup_M12)", "min(taum_M12, taup_M12)", "MeV", 'F');
        dataloader->AddVariable("min(taum_M23, taup_M23)", "min(taum_M23, taup_M23)", "MeV", 'F');
        dataloader->AddVariable("min( log(min(taum_pi1_IPCHI2_OWNPV, min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV)) ), log(min(taup_pi1_IPCHI2_OWNPV, min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) )", "min(log_taum_daughters_min_IP_CHI2_PV,log_taup_daughters_min_IP_CHI2_PV)", "", 'F');
        dataloader->AddVariable("log( sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) )", "log(taus_separation_3D)", "", 'F');
        dataloader->AddVariable("max( log(taup_FDCHI2_ORIVX), log(taum_FDCHI2_ORIVX) )", "max( log(taup_FDCHI2_ORIVX),log(taum_FDCHI2_ORIVX) )", "", 'F', 0, 1000);
        dataloader->AddVariable("min( log(taup_FDCHI2_ORIVX), log(taum_FDCHI2_ORIVX) )", "min( log(taup_FDCHI2_ORIVX),log(taum_FDCHI2_ORIVX) )", "", 'F', 0, 200);
        dataloader->AddVariable("max( log(1-abs(taup_DIRA_ORIVX))*sign(taup_DIRA_ORIVX), log(1-abs(taum_DIRA_ORIVX))*sign(taum_DIRA_ORIVX) )", "max(log(1-|taup_DIRA_ORIVX|)*sign,log(1-|taum_DIRA_ORIVX|)*sign)", "", 'F');
        dataloader->AddVariable("min( log(1-abs(taup_DIRA_ORIVX))*sign(taup_DIRA_ORIVX), log(1-abs(taum_DIRA_ORIVX))*sign(taum_DIRA_ORIVX) )", "min(log(1-|taup_DIRA_ORIVX|)*sign,log(1-|taum_DIRA_ORIVX|)*sign)", "", 'F');

        // Impact parameter of tau DVs wrt to K+ trajectory
        TString Cx_taup = "(taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY";
        TString Cy_taup = "(taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ";
        TString Cz_taup = "(taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX";
        TString C_taup = "sqrt( pow("+Cx_taup+",2) + pow("+Cy_taup+",2) + pow("+Cz_taup+",2) )";
        TString IP_taup = "2*"+C_taup+"/(sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";

        TString Cx_taum = "(taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY";
        TString Cy_taum = "(taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ";
        TString Cz_taum = "(taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX";
        TString C_taum = "sqrt( pow("+Cx_taum+",2) + pow("+Cy_taum+",2) + pow("+Cz_taum+",2) )";
        TString IP_taum = "2*"+C_taum+"/(sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";

        dataloader->AddVariable("max("+IP_taup+", "+IP_taum+")", "max(IP_taup,IP_taum)", "mm", 'F');
        dataloader->AddVariable("min("+IP_taup+", "+IP_taum+")", "min(IP_taup,IP_taum)", "mm", 'F');

        // Area of the triangle formed by the PV and the tau+ and tau- decay vertices
        TString a = "sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) )";
        TString b = "sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) )";
        TString c = "sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";
        TString s = "("+a+"+"+b+"+"+c+")*0.5";
        TString A = "sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";
        dataloader->AddVariable(A, "aread_triangle_PV_DVp_DVm", "mm^{2}", 'F');

        // Angle between K+ momentum and surface formed by PV, tau+ and tau- decay vertices
        TString v1x = "Bp_OWNPV_X - taup_ENDVERTEX_X";
        TString v1y = "Bp_OWNPV_Y - taup_ENDVERTEX_Y";
        TString v1z = "Bp_OWNPV_Z - taup_ENDVERTEX_Z";
        TString v2x = "taum_ENDVERTEX_X - taup_ENDVERTEX_X";
        TString v2y = "taum_ENDVERTEX_Y - taup_ENDVERTEX_Y";
        TString v2z = "taum_ENDVERTEX_Z - taup_ENDVERTEX_Z";
        TString nx = v2y+"*"+v1z+" - "+v2z+"*"+v1y;
        TString ny = v2z+"*"+v1x+" - "+v2x+"*"+v1z;
        TString nz = v2x+"*"+v1y+" - "+v2y+"*"+v1x;
        TString theta = "acos( ("+nx+"*Kp_PX + "+ny+"*Kp_PY + "+nz+"*Kp_PZ)/(sqrt( pow("+nx+",2) + pow("+ny+",2) + pow("+nz+",2)  )*sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) )) )";
        dataloader->AddVariable(theta, "angle between Kp momentum and PV-DVm-DVp surface", "", 'F');

        // tau (isolation)
        dataloader->AddVariable("max(Bp_CC_IT_tau1, Bp_CC_IT_tau2)", "max(CC_IT_tau1,CC_IT_tau2)", "", 'F');
        dataloader->AddVariable("max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2)", "max(CC_MAXPT_PT_tau1,CC_MAXPT_PT_tau2)", "MeV", 'F');
        dataloader->AddVariable("max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2)", "max(CC_MULT_tau1,CC_MULT_tau2)", "", 'F');
        dataloader->AddVariable("min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1)", "min(CC_PZ_tau1,CC_PZ_tau2)", "MeV", 'F');
        dataloader->AddVariable("min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1)", "min(CC_PZASYM_tau1,CC_PZASYM_tau2)", "", 'F');
        dataloader->AddVariable("min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1)", "min(CC_SPT_tau1,CC_SPT_tau2)", "MeV", 'F');
        dataloader->AddVariable("min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1)", "min(CC_VPT_tau1,CC_VPT_tau2)", "MeV", 'F');
        dataloader->AddVariable("max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2)", "max(NC_DELTAETA_tau1,NC_DELTAETA_tau2)", "", 'F', -1, 5);
        dataloader->AddVariable("max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2)", "max(NC_DELTAPHI_tau1,NC_DELTAPHI_tau2)", "", 'F', -1, 5);
        dataloader->AddVariable("max(Bp_NC_IT_tau1,Bp_NC_IT_tau2)", "max(NC_IT_tau1,NC_IT_tau2)", "", 'F', 0., 1.5);
        dataloader->AddVariable("max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2)", "max(NC_MAXPT_PT_tau1,NC_MAXPT_PT_tau2)", "MeV", 'F');
        dataloader->AddVariable("max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2)", "max(NC_MAXPT_PZ_tau1,NC_MAXPT_PZ_tau2)", "MeV", 'F');
        dataloader->AddVariable("min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1)", "min(NC_MULT_tau1,NC_MULT_tau2)", "", 'F');
        dataloader->AddVariable("min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1)", "min(NC_PTASYM_tau1,NC_PTASYM_tau2)", "MeV", 'F', -2, 3);
        dataloader->AddVariable("min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1)", "min(NC_PZ_tau1,NC_PZ_tau)2", "MeV", 'F');
        dataloader->AddVariable("min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1)", "min(NC_SPT_tau1,NC_SPT_tau2)", "MeV", 'F');
        dataloader->AddVariable("min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1)", "min(NC_VPT_tau1,NC_VPT_tau2)", "MeV", 'F');

        dataloader->AddVariable("max(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", "max(VTXISODCHI2ONETRACK_tau1, VTXISODCHI2ONETRACK_tau2)", "", 'F', 0, 150);
        dataloader->AddVariable("min(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", "min(VTXISODCHI2ONETRACK_tau1, VTXISODCHI2ONETRACK_tau2)", "", 'F', 0, 150);
        dataloader->AddVariable("max(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", "max(VTXISODCHI2TWOTRACK_tau1, VTXISODCHI2TWOTRACK_tau2)", "", 'F', 0, 100);
        dataloader->AddVariable("min(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", "min(VTXISODCHI2TWOTRACK_tau1, VTXISODCHI2TWOTRACK_tau2)", "", 'F', 0, 100);
        dataloader->AddVariable("max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", "max(VTXISONUMVTX_tau1,VTXISONUMVTX_tau2)", "", 'F', 0, 15);
        dataloader->AddVariable("min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", "min(VTXISONUMVTX_tau1,VTXISONUMVTX_tau2)", "", 'F', 0, 15);

        // kaon
        dataloader->AddVariable("log(1-Kp_ProbNNk)", "log(1-Kp_ProbNNk)", "", 'F');

        // B meson
        dataloader->AddVariable("Bp_PT", "Bp_PT", "MeV", 'F');
        dataloader->AddVariable("Bp_FDCHI2_OWNPV", "Bp_FDCHI2_OWNPV", "", 'F');

        // others
        dataloader->AddVariable("Bp_ConsBp_seq_12_chi2[0]", "DTF #chi^{2}", "", 'F');
        // IP of PV wrt to K+ trajectory
        TString Cx_PV = "(Bp_OWNPV_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (Bp_OWNPV_Z-Bp_ENDVERTEX_Z)*Kp_PY";
        TString Cy_PV = "(Bp_OWNPV_Z-Bp_ENDVERTEX_Z)*Kp_PX - (Bp_OWNPV_X-Bp_ENDVERTEX_X)*Kp_PZ";
        TString Cz_PV = "(Bp_OWNPV_X-Bp_ENDVERTEX_X)*Kp_PY - (Bp_OWNPV_Y-Bp_ENDVERTEX_Y)*Kp_PX";
        TString C_PV = "sqrt( pow("+Cx_PV+",2) + pow("+Cy_PV+",2) + pow("+Cz_PV+",2) )";
        TString IP_PV = "2*"+C_PV+"/(sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";
        //dataloader->AddVariable(IP_PV, "IP_PV", "mm", 'F');
        
    }

    // Set the event weights per tree
    double sigWeight = 1.0;
    double bkgWeight = 1.0;

    // Add signal and background trees to the factory
    dataloader->AddSignalTree(sig_tree,sigWeight);
    dataloader->AddBackgroundTree(bkg_tree,bkgWeight);

    if(tmva_type == 2)
    {
        dataloader->AddCut(TCut("(abs(Bp_ConsBp_seq_12_PX[0]) < 50000) && (Bp_ConsBp_seq_12_PZ[0] > 0) && (Bp_ConsBp_seq_12_PZ[0] < 800000) && (Bp_ConsBp_seq_12_MERR > 0) && (Bp_ConsBp_seq_12_MERR < 1500)"), "Signal");
        dataloader->AddCut(TCut("(abs(Bp_ConsBp_seq_12_PX[0]) < 50000) && (Bp_ConsBp_seq_12_PZ[0] > 0) && (Bp_ConsBp_seq_12_PZ[0] < 800000) && (Bp_ConsBp_seq_12_MERR > 0) && (Bp_ConsBp_seq_12_MERR < 1500)"), "Background");
    }
    else if(tmva_type == 6)
    {
        dataloader->AddCut(TCut("(max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2) > 0) && (max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2) > -1) && (max(Bp_NC_IT_tau1,Bp_NC_IT_tau2) > 0) && (min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1) > 0) && (min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1) > -1)"), "Signal");
        dataloader->AddCut(TCut("(max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2) > 0) && (max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2) > -1) && (max(Bp_NC_IT_tau1,Bp_NC_IT_tau2) > 0) && (min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1) > 0) && (min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1) > -1)"), "Background");
    }


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