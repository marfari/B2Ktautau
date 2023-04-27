// (tmva_type)

void TMVA_training(int year, int tmva_type, int background_proxy, TString RS_DATA_files, TString WS_DATA_files, TString SIGNAL_MC_files){
    // Load TMVA libraries
    TMVA::Tools::Instance(); 

    // Create output file for TMVA
    TString outfileName = Form("/panfs/felician/B2Ktautau/workflow/TMVA_training/201%i/tmva_output_type%i_bkgProxy%i.root", year, tmva_type, background_proxy);
    TFile* outputFile = TFile::Open( outfileName, "RECREATE"); 
    outputFile->cd();

    // Create TMVA factory
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",outputFile,"V:!Silent:Color:Transformations=I:DrawProgressBar:AnalysisType=Classification"); 
    TMVA::DataLoader *dataloader = new TMVA::DataLoader(Form("tmva_dataset_type%i_bkgProxy%i",tmva_type,background_proxy));

    TFileCollection* fc_data;
    if(background_proxy == 0) // RS data (after pre-selections)
    {
        fc_data = new TFileCollection("data", "data", RS_DATA_files, 1); // takes only the 1st root file
    } 
    else if(background_proxy == 1) // WS data (after pre-selections)
    {
        fc_data = new TFileCollection("data", "data", WS_DATA_files, 1); // takes only the 1st root file
    }
    TChain* bkg_tree = new TChain("DecayTree");
    bkg_tree->AddFileInfoList((TCollection*)fc_data->GetList());

    TFileCollection* fc_mc = new TFileCollection("MC", "MC", SIGNAL_MC_files); // Signal MC (component=0) files (after pre-selections)
    TChain* sig_tree = new TChain("DecayTree");
    sig_tree->AddFileInfoList((TCollection*)fc_mc->GetList());

    // Pre-selections
    TCut sig_mass = "(Bp_ConsBp_seq_12_M[0]>4700) && (Bp_ConsBp_seq_12_M[0]<5800)";
    TCut bkg_mass;
    if(background_proxy==0){bkg_mass = "(Bp_ConsBp_seq_12_M[0]>6500)";}
    else if(background_proxy==1){bkg_mass = "(Bp_ConsBp_seq_12_M[0]>4700) && (Bp_ConsBp_seq_12_M[0]<5800)";}

    TCut mycuts = sig_mass;
    TCut mycutb = bkg_mass;

    cout << "# signal events after pre-selections in signal region =  " << sig_tree->GetEntries(mycuts) << endl;
    cout << "# background events after pre-selections in signal region = " << bkg_tree->GetEntries(mycutb) << endl;

    // Add variables to the training
    if(tmva_type == 6) // mix
    {
        if(year == 8)
        {
            // tau (kinematics)
            dataloader->AddVariable("TMath::Max(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", "Max(DVm_chi2,DVp_chi2)", "", 'F');
            dataloader->AddVariable("TMath::Min(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", "Min(DVm_chi2,DVp_chi2)", "", 'F');
            dataloader->AddVariable("TMath::Max(taum_M, taup_M)", "Max(taum_M,taup_M)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(taum_M12, taup_M12)", "Max(taum_M12,taup_M12)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(taum_M23, taup_M23)", "Max(taum_M23,taup_M23)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(taum_M, taup_M)", "Min(taum_M,taup_M)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(taum_M12, taup_M12)", "Min(taum_M12,taup_M12)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(taum_M23, taup_M23)", "Min(taum_M23,taup_M23)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min( TMath::Log( TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV)) ), TMath::Log( TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV)) ) )", "Min(log_taum_daughters_min_IP_CHI2_PV,log_taup_daughters_min_IP_CHI2_PV)", "", 'F');
            dataloader->AddVariable("TMath::Log( TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) )", "Log(taus_separation_3D)", "", 'F');
            dataloader->AddVariable("TMath::Max( TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX) )", "Max( Log(taup_FDCHI2_ORIVX),Log(taum_FDCHI2_ORIVX) )", "", 'F');
            dataloader->AddVariable("TMath::Min( TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX) )", "Min( Log(taup_FDCHI2_ORIVX),Log(taum_FDCHI2_ORIVX) )", "", 'F');
            dataloader->AddVariable("TMath::Max( TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) )", "Max(Log(1-|taup_DIRA_ORIVX|)*sign,Log(1-|taum_DIRA_ORIVX|)*sign)", "", 'F');
            dataloader->AddVariable("TMath::Min( TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) )", "Min(Log(1-|taup_DIRA_ORIVX|)*sign,Log(1-|taum_DIRA_ORIVX|)*sign)", "", 'F');

            // Impact parameter of tau DVs wrt to K+ trajectory
            TString Cx_taup = "(taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY";
            TString Cy_taup = "(taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ";
            TString Cz_taup = "(taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX";
            TString C_taup = "TMath::Sqrt( pow("+Cx_taup+",2) + pow("+Cy_taup+",2) + pow("+Cz_taup+",2) )";
            TString IP_taup = "2*"+C_taup+"/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";

            TString Cx_taum = "(taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY";
            TString Cy_taum = "(taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ";
            TString Cz_taum = "(taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX";
            TString C_taum = "TMath::Sqrt( pow("+Cx_taum+",2) + pow("+Cy_taum+",2) + pow("+Cz_taum+",2) )";
            TString IP_taum = "2*"+C_taum+"/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";

            dataloader->AddVariable("TMath::Max("+IP_taup+", "+IP_taum+")", "Max(IP_taup,IP_taum)", "mm", 'F');
            dataloader->AddVariable("TMath::Min("+IP_taup+", "+IP_taum+")", "Min(IP_taup,IP_taum)", "mm", 'F');

            // Area of the triangle formed by the PV and the tau+ and tau- decay vertices
            TString a = "TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) )";
            TString b = "TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) )";
            TString c = "TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";
            TString s = "("+a+"+"+b+"+"+c+")*0.5";
            TString A = "TMath::Sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";
            // dataloader->AddVariable(A, "area_triangle_PV_DVp_DVm", "mm^{2}", 'F');

            // tau (isolation)
            dataloader->AddVariable("TMath::Max(Bp_CC_IT_tau1, Bp_CC_IT_tau2)", "Max(CC_IT_tau1,CC_IT_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2)", "Max(CC_MAXPT_PT_tau1,CC_MAXPT_PT_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2)", "Max(CC_MULT_tau1,CC_MULT_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1)", "Min(CC_PZ_tau1,CC_PZ_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1)", "Min(CC_PZASYM_tau1,CC_PZASYM_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1)", "Min(CC_SPT_tau1,CC_SPT_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1)", "Min(CC_VPT_tau1,CC_VPT_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2)", "Max(NC_DELTAETA_tau1,NC_DELTAETA_tau2)", "", 'F', -1, 5);
            dataloader->AddVariable("TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2)", "Max(NC_DELTAPHI_tau1,NC_DELTAPHI_tau2)", "", 'F', -1, 5);
            dataloader->AddVariable("TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2)", "Max(NC_IT_tau1,NC_IT_tau2)", "", 'F', 0., 1.5);
            dataloader->AddVariable("TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2)", "Max(NC_MAXPT_PT_tau1,NC_MAXPT_PT_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2)", "Max(NC_MAXPT_PZ_tau1,NC_MAXPT_PZ_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1)", "Min(NC_MULT_tau1,NC_MULT_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1)", "Min(NC_PTASYM_tau1,NC_PTASYM_tau2)", "MeV", 'F', -2, 3);
            dataloader->AddVariable("TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1)", "Min(NC_PZ_tau1,NC_PZ_tau)2", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1)", "Min(NC_SPT_tau1,NC_SPT_tau2)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1)", "Min(NC_VPT_tau1,NC_VPT_tau2)", "MeV", 'F');

            dataloader->AddVariable("TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", "Max(VTXISODCHI2ONETRACK_tau1, VTXISODCHI2ONETRACK_tau2)", "", 'F', 0, 150);
            dataloader->AddVariable("TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", "Min(VTXISODCHI2ONETRACK_tau1, VTXISODCHI2ONETRACK_tau2)", "", 'F', 0, 150);
            dataloader->AddVariable("TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", "Max(VTXISODCHI2TWOTRACK_tau1, VTXISODCHI2TWOTRACK_tau2)", "", 'F', 0, 100);
            dataloader->AddVariable("TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", "Min(VTXISODCHI2TWOTRACK_tau1, VTXISODCHI2TWOTRACK_tau2)", "", 'F', 0, 100);
            dataloader->AddVariable("TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", "Max(VTXISONUMVTX_tau1,VTXISONUMVTX_tau2)", "", 'F', 0, 15);
            dataloader->AddVariable("TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", "Min(VTXISONUMVTX_tau1,VTXISONUMVTX_tau2)", "", 'F', 0, 15);

            // kaon
            dataloader->AddVariable("TMath::Log(1-Kp_ProbNNk)", "Log(1-Kp_ProbNNk)", "", 'F');

            // B meson
            dataloader->AddVariable("Bp_PT", "Bp_PT", "MeV", 'F');

            // others
            dataloader->AddVariable("Bp_ConsBp_seq_12_chi2[0]", "DTF #chi^{2}", "", 'F');
            // IP of PV wrt to K+ trajectory
            TString Cx_PV = "(Bp_OWNPV_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (Bp_OWNPV_Z-Bp_ENDVERTEX_Z)*Kp_PY";
            TString Cy_PV = "(Bp_OWNPV_Z-Bp_ENDVERTEX_Z)*Kp_PX - (Bp_OWNPV_X-Bp_ENDVERTEX_X)*Kp_PZ";
            TString Cz_PV = "(Bp_OWNPV_X-Bp_ENDVERTEX_X)*Kp_PY - (Bp_OWNPV_Y-Bp_ENDVERTEX_Y)*Kp_PX";
            TString C_PV = "TMath::Sqrt( pow("+Cx_PV+",2) + pow("+Cy_PV+",2) + pow("+Cz_PV+",2) )";
            TString IP_PV = "2*"+C_PV+"/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";
            //dataloader->AddVariable(IP_PV, "IP_PV", "mm", 'F');
        }
        else
        {
            // tau (kinematics)
            dataloader->AddVariable("TMath::Max(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", "Max(DVm_chi2,DVp_chi2)", "", 'F');
            dataloader->AddVariable("TMath::Min(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", "Min(DVm_chi2,DVp_chi2)", "", 'F');
            dataloader->AddVariable("TMath::Max(taum_M, taup_M)", "Max(taum_M,taup_M)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(taum_M12, taup_M12)", "Max(taum_M12,taup_M12)", "MeV", 'F');
            dataloader->AddVariable("TMath::Max(taum_M23, taup_M23)", "Max(taum_M23,taup_M23)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(taum_M, taup_M)", "Min(taum_M,taup_M)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(taum_M12, taup_M12)", "Min(taum_M12,taup_M12)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min(taum_M23, taup_M23)", "Min(taum_M23,taup_M23)", "MeV", 'F');
            dataloader->AddVariable("TMath::Min( TMath::Log( TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV)) ), TMath::Log( TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV)) ) )", "Min(log_taum_daughters_min_IP_CHI2_PV,log_taup_daughters_min_IP_CHI2_PV)", "", 'F');
            dataloader->AddVariable("TMath::Log( TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) )", "Log(taus_separation_3D)", "", 'F');
            dataloader->AddVariable("TMath::Max( TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX) )", "Max( Log(taup_FDCHI2_ORIVX),Log(taum_FDCHI2_ORIVX) )", "", 'F');
            dataloader->AddVariable("TMath::Min( TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX) )", "Min( Log(taup_FDCHI2_ORIVX),Log(taum_FDCHI2_ORIVX) )", "", 'F');
            dataloader->AddVariable("TMath::Max( TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) )", "Max(Log(1-|taup_DIRA_ORIVX|)*sign,Log(1-|taum_DIRA_ORIVX|)*sign)", "", 'F');
            dataloader->AddVariable("TMath::Min( TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) )", "Min(Log(1-|taup_DIRA_ORIVX|)*sign,Log(1-|taum_DIRA_ORIVX|)*sign)", "", 'F');

            // Impact parameter of tau DVs wrt to K+ trajectory
            TString Cx_taup = "(taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY";
            TString Cy_taup = "(taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ";
            TString Cz_taup = "(taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX";
            TString C_taup = "TMath::Sqrt( pow("+Cx_taup+",2) + pow("+Cy_taup+",2) + pow("+Cz_taup+",2) )";
            TString IP_taup = "2*"+C_taup+"/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";

            TString Cx_taum = "(taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY";
            TString Cy_taum = "(taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ";
            TString Cz_taum = "(taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX";
            TString C_taum = "TMath::Sqrt( pow("+Cx_taum+",2) + pow("+Cy_taum+",2) + pow("+Cz_taum+",2) )";
            TString IP_taum = "2*"+C_taum+"/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";

            dataloader->AddVariable("TMath::Max("+IP_taup+", "+IP_taum+")", "Max(IP_taup,IP_taum)", "mm", 'F');
            dataloader->AddVariable("TMath::Min("+IP_taup+", "+IP_taum+")", "Min(IP_taup,IP_taum)", "mm", 'F');

            // Area of the triangle formed by the PV and the tau+ and tau- decay vertices
            TString a = "TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) )";
            TString b = "TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) )";
            TString c = "TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";
            TString s = "("+a+"+"+b+"+"+c+")*0.5";
            TString A = "TMath::Sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";
            // dataloader->AddVariable(A, "area_triangle_PV_DVp_DVm", "mm^{2}", 'F');

            // tau (isolation)
            dataloader->AddVariable("TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", "Max(CONEMULT_tau1,CONEMULT_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", "Max(CONEPASYM_tau1,CONEPASYM_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", "Max(CONEPTASYM_tau1,CONEPTASYM_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", "Max(CONEDELTAETA_tau1,CONEDELTAETA_tau2)", "", 'F');
            
            dataloader->AddVariable("TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", "Max(CONEMULT_tau1,CONEMULT_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", "Max(CONEPASYM_tau1,CONEPASYM_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", "Max(CONEPTASYM_tau1,CONEPTASYM_tau2)", "", 'F');
            dataloader->AddVariable("TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", "Max(CONEDELTAETA_tau1,CONEDELTAETA_tau2)", "", 'F');

            // kaon
            dataloader->AddVariable("TMath::Log(1-Kp_ProbNNk)", "Log(1-Kp_ProbNNk)", "", 'F');

            // B meson
            dataloader->AddVariable("Bp_PT", "Bp_PT", "MeV", 'F');

            // others
            dataloader->AddVariable("Bp_ConsBp_seq_12_chi2[0]", "DTF #chi^{2}", "", 'F');
            // IP of PV wrt to K+ trajectory
            TString Cx_PV = "(Bp_OWNPV_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (Bp_OWNPV_Z-Bp_ENDVERTEX_Z)*Kp_PY";
            TString Cy_PV = "(Bp_OWNPV_Z-Bp_ENDVERTEX_Z)*Kp_PX - (Bp_OWNPV_X-Bp_ENDVERTEX_X)*Kp_PZ";
            TString Cz_PV = "(Bp_OWNPV_X-Bp_ENDVERTEX_X)*Kp_PY - (Bp_OWNPV_Y-Bp_ENDVERTEX_Y)*Kp_PX";
            TString C_PV = "TMath::Sqrt( pow("+Cx_PV+",2) + pow("+Cy_PV+",2) + pow("+Cz_PV+",2) )";
            TString IP_PV = "2*"+C_PV+"/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ))";
            //dataloader->AddVariable(IP_PV, "IP_PV", "mm", 'F');
        }
    }

    // Set the event weights per tree
    double sigWeight = 1.0;
    double bkgWeight = 1.0;

    // Add signal and background trees to the factory
    dataloader->AddSignalTree(sig_tree,sigWeight);
    dataloader->AddBackgroundTree(bkg_tree,bkgWeight);

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