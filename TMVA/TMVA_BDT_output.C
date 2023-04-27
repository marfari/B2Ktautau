
void TMVA_BDT_output(int year, int species, int tmva_type, int background_proxy, TString MC_files, TString MC_0_files, TString MC_1_files, TString MC_2_files, TString RS_DATA_files, TString WS_DATA_files, TString weights,  int line)
{
    // creates a root file containing the BDT response 

    TFileCollection* fc; 
    TFile* target = new TFile(Form("/panfs/felician/B2Ktautau/workflow/TMVA_response/201%i/Species_%i/BDT_output_tmva_type%i_bkgProxy%i_%i.root",year,species,tmva_type,background_proxy,line), "RECREATE");

    if(species == 0) // MC (all components)
    {   
        fc = new TFileCollection("fc", "fc", MC_files, 1, line);
    }
    else if(species == 10) // MC (3pi3pi)
    {
        fc = new TFileCollection("fc", "fc", MC_0_files, 1, line);
    }
    else if(species == 11) // MC (3pi3pi pi0)
    {
        fc = new TFileCollection("fc", "fc", MC_1_files, 1, line);
    }
    else if(species == 12) // MC (3pi3pi 2pi0)
    {
        fc = new TFileCollection("fc", "fc", MC_2_files, 1, line);
    }
    else if(species == 1) // RS data
    {
        fc = new TFileCollection("fc", "fc", RS_DATA_files, 1, line);
    }
    else if(species == 2) // WS data
    {
        fc = new TFileCollection("fc", "fc", WS_DATA_files, 1, line);
    }
    else
    {
        cout << "Species is wrong. Try either 0 for 'MC', 1 for 'RS_data' or 2 for 'WS_data'" << endl;
        return;
    }
    TChain* dataTree = new TChain("DecayTree");
    dataTree->AddFileInfoList((TCollection*)fc->GetList());

    TTree* tree = new TTree("DecayTree", "DecayTree");

    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "V:Color:!Silent" ); 

    Float_t var8[42];
    Float_t var[27];

    if(year == 8)
    {
        // tau (kinematics)
        reader->AddVariable("TMath::Max(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2)", &var8[0]);
        reader->AddVariable("TMath::Min(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2)", &var8[1]);
        reader->AddVariable("TMath::Max(taum_M,taup_M)", &var8[2]);
        reader->AddVariable("TMath::Max(taum_M12,taup_M12)", &var8[3]);
        reader->AddVariable("TMath::Max(taum_M23,taup_M23)", &var8[4]);
        reader->AddVariable("TMath::Min(taum_M,taup_M)", &var8[5]);
        reader->AddVariable("TMath::Min(taum_M12,taup_M12)", &var8[6]);
        reader->AddVariable("TMath::Min(taum_M23,taup_M23)", &var8[7]);
        reader->AddVariable("TMath::Min( TMath::Log(TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV)) ), "
                                        "TMath::Log(TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) )", &var8[8]);
        reader->AddVariable("TMath::Log( TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + "
                                                    "pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + "
                                                    "pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) )", &var8[9]);
        reader->AddVariable("TMath::Max(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX))", &var8[10]);
        reader->AddVariable("TMath::Min(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX))", &var8[11]);
        reader->AddVariable("TMath::Max(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX),"
                                    " TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX))", &var8[12]);
        reader->AddVariable("TMath::Min(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX),"
                                    " TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX))", &var8[13]);
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
        reader->AddVariable("TMath::Max("+IP_taup+", "+IP_taum+")", &var8[14]);
        reader->AddVariable("TMath::Min("+IP_taup+", "+IP_taum+")", &var8[15]);
        // Area of the triangle formed by the PV and the tau+ and tau- decay vertices
        TString a = "TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) )";
        TString b = "TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) )";
        TString c = "TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";
        TString s = "("+a+"+"+b+"+"+c+")*0.5";
        TString A = "TMath::Sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";
        // reader->AddVariable(A, &var8[16]);
        // tau (isolation)
        reader->AddVariable("TMath::Max(Bp_CC_IT_tau1, Bp_CC_IT_tau2)", &var8[16]);
        reader->AddVariable("TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2)", &var8[17]);
        reader->AddVariable("TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2)", &var8[18]);
        reader->AddVariable("TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1)", &var8[19]);
        reader->AddVariable("TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1)", &var8[20]);
        reader->AddVariable("TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1)", &var8[21]);
        reader->AddVariable("TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1)", &var8[22]);
        reader->AddVariable("TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2)", &var8[23]);
        reader->AddVariable("TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2)", &var8[24]);
        reader->AddVariable("TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2)", &var8[25]);
        reader->AddVariable("TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2)", &var8[26]);
        reader->AddVariable("TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2)", &var8[27]);
        reader->AddVariable("TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1)", &var8[28]);
        reader->AddVariable("TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1)", &var8[29]);
        reader->AddVariable("TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1)", &var8[30]);
        reader->AddVariable("TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1)", &var8[31]);
        reader->AddVariable("TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1)", &var8[32]);
        reader->AddVariable("TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", &var8[33]);
        reader->AddVariable("TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", &var8[34]);
        reader->AddVariable("TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", &var8[35]);
        reader->AddVariable("TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", &var8[36]);
        reader->AddVariable("TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", &var8[37]);
        reader->AddVariable("TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", &var8[38]);
        // kaon
        reader->AddVariable("TMath::Log(1-Kp_ProbNNk)", &var8[39]);
        // B meson
        reader->AddVariable("Bp_PT", &var8[40]);
        // others
        reader->AddVariable("Bp_ConsBp_seq_12_chi2[0]", &var8[41]);
    }
    else
    {
        // tau (kinematics)
        reader->AddVariable("TMath::Max(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2)", &var[0]);
        reader->AddVariable("TMath::Min(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2)", &var[1]);
        reader->AddVariable("TMath::Max(taum_M,taup_M)", &var[2]);
        reader->AddVariable("TMath::Max(taum_M12,taup_M12)", &var[3]);
        reader->AddVariable("TMath::Max(taum_M23,taup_M23)", &var[4]);
        reader->AddVariable("TMath::Min(taum_M,taup_M)", &var[5]);
        reader->AddVariable("TMath::Min(taum_M12,taup_M12)", &var[6]);
        reader->AddVariable("TMath::Min(taum_M23,taup_M23)", &var[7]);
        reader->AddVariable("TMath::Min( TMath::Log(TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV)) ), "
                                        "TMath::Log(TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) )", &var[8]);
        reader->AddVariable("TMath::Log( TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + "
                                                    "pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + "
                                                    "pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) )", &var[9]);
        reader->AddVariable("TMath::Max(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX))", &var[10]);
        reader->AddVariable("TMath::Min(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX))", &var[11]);
        reader->AddVariable("TMath::Max(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX),"
                                    " TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX))", &var[12]);
        reader->AddVariable("TMath::Min(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX),"
                                    " TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX))", &var[13]);
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
        reader->AddVariable("TMath::Max("+IP_taup+", "+IP_taum+")", &var[14]);
        reader->AddVariable("TMath::Min("+IP_taup+", "+IP_taum+")", &var[15]);
        // Area of the triangle formed by the PV and the tau+ and tau- decay vertices
        TString a = "TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) )";
        TString b = "TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) )";
        TString c = "TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";
        TString s = "("+a+"+"+b+"+"+c+")*0.5";
        TString A = "TMath::Sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";
        // reader->AddVariable(A, &var[16]);
        // tau (isolation)
        reader->AddVariable("TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", &var[16]);
        reader->AddVariable("TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", &var[17]);
        reader->AddVariable("TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", &var[18]);
        reader->AddVariable("TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", &var[19]);

        reader->AddVariable("TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", &var[20]);
        reader->AddVariable("TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", &var[21]);
        reader->AddVariable("TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", &var[22]);
        reader->AddVariable("TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", &var[23]);
        // kaon
        reader->AddVariable("TMath::Log(1-Kp_ProbNNk)", &var[24]);
        // B meson
        reader->AddVariable("Bp_PT", &var[25]);
        // others
        reader->AddVariable("Bp_ConsBp_seq_12_chi2[0]", &var[26]);
    }

    reader->BookMVA("BDT method", Form("/panfs/felician/B2Ktautau/workflow/TMVA_training/201%i/tmva_dataset_type%i_bkgProxy%i/weights/TMVAClassification_BDT.weights.xml",year,tmva_type,background_proxy));  

    // training features
    Double_t taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2, taum_M, taum_M12, taum_M13, taum_M23, taup_M, taup_M12, taup_M13, taup_M23;
    Double_t taum_pi1_IPCHI2_OWNPV, taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV, taup_pi1_IPCHI2_OWNPV, taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV;
    Double_t taup_ENDVERTEX_X, taup_ENDVERTEX_Y, taup_ENDVERTEX_Z, taum_ENDVERTEX_X, taum_ENDVERTEX_Y, taum_ENDVERTEX_Z;
    Double_t taup_FDCHI2_ORIVX, taum_FDCHI2_ORIVX, taup_DIRA_ORIVX, taum_DIRA_ORIVX;
    Double_t Bp_ENDVERTEX_X, Bp_ENDVERTEX_Y, Bp_ENDVERTEX_Z, Bp_ENDVERTEX_CHI2, Bp_PT;
    Double_t Kp_PX, Kp_PY, Kp_PZ, Kp_P, Kp_PT, Kp_ProbNNk;

    Double_t Bp_CC_IT_tau1, Bp_CC_MAXPT_PT_tau1, Bp_CC_PZ_tau2, Bp_CC_PZASYM_tau2, Bp_CC_SPT_tau2, Bp_CC_VPT_tau2;
    Double_t Bp_CC_IT_tau2, Bp_CC_MAXPT_PT_tau2, Bp_CC_PZ_tau1, Bp_CC_PZASYM_tau1, Bp_CC_SPT_tau1, Bp_CC_VPT_tau1;
    Double_t Bp_CC_MULT_tau1, Bp_NC_MULT_tau2, Bp_VTXISONUMVTX_tau2; 
    Double_t Bp_CC_MULT_tau2, Bp_NC_MULT_tau1, Bp_VTXISONUMVTX_tau1;
    Double_t Bp_NC_DELTAETA_tau1, Bp_NC_DELTAPHI_tau1, Bp_NC_IT_tau1, Bp_NC_MAXPT_PT_tau1, Bp_NC_MAXPT_PZ_tau1, Bp_NC_PTASYM_tau2, Bp_NC_PZ_tau2, Bp_NC_SPT_tau2, Bp_NC_VPT_tau2;
    Double_t Bp_NC_DELTAETA_tau2, Bp_NC_DELTAPHI_tau2, Bp_NC_IT_tau2, Bp_NC_MAXPT_PT_tau2, Bp_NC_MAXPT_PZ_tau2, Bp_NC_PTASYM_tau1, Bp_NC_PZ_tau1, Bp_NC_SPT_tau1, Bp_NC_VPT_tau1;
    Double_t Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2ONETRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1;    

    Double_t Bp_CONEMULT_tau1, Bp_CONEMULT_tau2, Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2, Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2, Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2;

    Float_t Bp_ConsBp_seq_12_chi2[100];
    Double_t Bp_OWNPV_X, Bp_OWNPV_Y, Bp_OWNPV_Z, Bp_FDCHI2_OWNPV;
    
    dataTree->SetBranchAddress("taum_ENDVERTEX_CHI2", &taum_ENDVERTEX_CHI2);
    dataTree->SetBranchAddress("taup_ENDVERTEX_CHI2", &taup_ENDVERTEX_CHI2);
    dataTree->SetBranchAddress("taum_M", &taum_M);
    dataTree->SetBranchAddress("taum_M12", &taum_M12);
    dataTree->SetBranchAddress("taum_M13", &taum_M13);
    dataTree->SetBranchAddress("taum_M23", &taum_M23);
    dataTree->SetBranchAddress("taup_M", &taup_M);
    dataTree->SetBranchAddress("taup_M12", &taup_M12);
    dataTree->SetBranchAddress("taup_M13", &taup_M13);
    dataTree->SetBranchAddress("taup_M23", &taup_M23);
    dataTree->SetBranchAddress("taum_pi1_IPCHI2_OWNPV", &taum_pi1_IPCHI2_OWNPV);
    dataTree->SetBranchAddress("taum_pi2_IPCHI2_OWNPV", &taum_pi2_IPCHI2_OWNPV);
    dataTree->SetBranchAddress("taum_pi3_IPCHI2_OWNPV", &taum_pi3_IPCHI2_OWNPV);
    dataTree->SetBranchAddress("taup_pi1_IPCHI2_OWNPV", &taup_pi1_IPCHI2_OWNPV);
    dataTree->SetBranchAddress("taup_pi2_IPCHI2_OWNPV", &taup_pi2_IPCHI2_OWNPV);
    dataTree->SetBranchAddress("taup_pi3_IPCHI2_OWNPV", &taup_pi3_IPCHI2_OWNPV);
    dataTree->SetBranchAddress("taup_ENDVERTEX_X", &taup_ENDVERTEX_X);
    dataTree->SetBranchAddress("taup_ENDVERTEX_Y", &taup_ENDVERTEX_Y);
    dataTree->SetBranchAddress("taup_ENDVERTEX_Z", &taup_ENDVERTEX_Z);
    dataTree->SetBranchAddress("taum_ENDVERTEX_X", &taum_ENDVERTEX_X);
    dataTree->SetBranchAddress("taum_ENDVERTEX_Y", &taum_ENDVERTEX_Y);
    dataTree->SetBranchAddress("taum_ENDVERTEX_Z", &taum_ENDVERTEX_Z);
    dataTree->SetBranchAddress("taup_FDCHI2_ORIVX", &taup_FDCHI2_ORIVX);
    dataTree->SetBranchAddress("taum_FDCHI2_ORIVX", &taum_FDCHI2_ORIVX);
    dataTree->SetBranchAddress("taup_DIRA_ORIVX", &taup_DIRA_ORIVX);
    dataTree->SetBranchAddress("taum_DIRA_ORIVX", &taum_DIRA_ORIVX);
    dataTree->SetBranchAddress("Bp_ENDVERTEX_X", &Bp_ENDVERTEX_X);
    dataTree->SetBranchAddress("Bp_ENDVERTEX_Y", &Bp_ENDVERTEX_Y);
    dataTree->SetBranchAddress("Bp_ENDVERTEX_Z", &Bp_ENDVERTEX_Z);
    dataTree->SetBranchAddress("Kp_PX", &Kp_PX);
    dataTree->SetBranchAddress("Kp_PY", &Kp_PY);
    dataTree->SetBranchAddress("Kp_PZ", &Kp_PZ);
    dataTree->SetBranchAddress("Kp_P", &Kp_P);
    if(year == 8)
    {
        dataTree->SetBranchAddress("Bp_CC_IT_tau1", &Bp_CC_IT_tau1);
        dataTree->SetBranchAddress("Bp_CC_IT_tau2", &Bp_CC_IT_tau2);
        dataTree->SetBranchAddress("Bp_CC_MAXPT_PT_tau1", &Bp_CC_MAXPT_PT_tau1);
        dataTree->SetBranchAddress("Bp_CC_MAXPT_PT_tau2", &Bp_CC_MAXPT_PT_tau2);
        dataTree->SetBranchAddress("Bp_CC_MULT_tau1", &Bp_CC_MULT_tau1);
        dataTree->SetBranchAddress("Bp_CC_MULT_tau2", &Bp_CC_MULT_tau2);
        dataTree->SetBranchAddress("Bp_CC_PZ_tau1", &Bp_CC_PZ_tau1);
        dataTree->SetBranchAddress("Bp_CC_PZ_tau2", &Bp_CC_PZ_tau2);
        dataTree->SetBranchAddress("Bp_CC_PZASYM_tau1", &Bp_CC_PZASYM_tau1);
        dataTree->SetBranchAddress("Bp_CC_PZASYM_tau2", &Bp_CC_PZASYM_tau2);
        dataTree->SetBranchAddress("Bp_CC_SPT_tau1", &Bp_CC_SPT_tau1);
        dataTree->SetBranchAddress("Bp_CC_SPT_tau2", &Bp_CC_SPT_tau2);
        dataTree->SetBranchAddress("Bp_CC_VPT_tau1", &Bp_CC_VPT_tau1);
        dataTree->SetBranchAddress("Bp_CC_VPT_tau2", &Bp_CC_VPT_tau2);
        dataTree->SetBranchAddress("Bp_NC_DELTAETA_tau1", &Bp_NC_DELTAETA_tau1);
        dataTree->SetBranchAddress("Bp_NC_DELTAETA_tau2", &Bp_NC_DELTAETA_tau2);
        dataTree->SetBranchAddress("Bp_NC_DELTAPHI_tau1", &Bp_NC_DELTAPHI_tau1);
        dataTree->SetBranchAddress("Bp_NC_DELTAPHI_tau2", &Bp_NC_DELTAPHI_tau2);
        dataTree->SetBranchAddress("Bp_NC_IT_tau1", &Bp_NC_IT_tau1);
        dataTree->SetBranchAddress("Bp_NC_IT_tau2", &Bp_NC_IT_tau2);
        dataTree->SetBranchAddress("Bp_NC_MAXPT_PT_tau1", &Bp_NC_MAXPT_PT_tau1);
        dataTree->SetBranchAddress("Bp_NC_MAXPT_PT_tau2", &Bp_NC_MAXPT_PT_tau2);
        dataTree->SetBranchAddress("Bp_NC_MAXPT_PZ_tau1", &Bp_NC_MAXPT_PZ_tau1);
        dataTree->SetBranchAddress("Bp_NC_MAXPT_PZ_tau2", &Bp_NC_MAXPT_PZ_tau2);
        dataTree->SetBranchAddress("Bp_NC_MULT_tau1", &Bp_NC_MULT_tau1);
        dataTree->SetBranchAddress("Bp_NC_MULT_tau2", &Bp_NC_MULT_tau2);
        dataTree->SetBranchAddress("Bp_NC_PTASYM_tau1", &Bp_NC_PTASYM_tau1);
        dataTree->SetBranchAddress("Bp_NC_PTASYM_tau2", &Bp_NC_PTASYM_tau2);
        dataTree->SetBranchAddress("Bp_NC_PZ_tau1", &Bp_NC_PZ_tau1);
        dataTree->SetBranchAddress("Bp_NC_PZ_tau2", &Bp_NC_PZ_tau2);
        dataTree->SetBranchAddress("Bp_NC_SPT_tau1", &Bp_NC_SPT_tau1);
        dataTree->SetBranchAddress("Bp_NC_SPT_tau2", &Bp_NC_SPT_tau2);
        dataTree->SetBranchAddress("Bp_NC_VPT_tau1", &Bp_NC_VPT_tau1);
        dataTree->SetBranchAddress("Bp_NC_VPT_tau2", &Bp_NC_VPT_tau2);
        dataTree->SetBranchAddress("Bp_VTXISODCHI2ONETRACK_tau1", &Bp_VTXISODCHI2ONETRACK_tau1);
        dataTree->SetBranchAddress("Bp_VTXISODCHI2ONETRACK_tau2", &Bp_VTXISODCHI2ONETRACK_tau2);
        dataTree->SetBranchAddress("Bp_VTXISODCHI2TWOTRACK_tau1", &Bp_VTXISODCHI2TWOTRACK_tau1);
        dataTree->SetBranchAddress("Bp_VTXISODCHI2TWOTRACK_tau2", &Bp_VTXISODCHI2TWOTRACK_tau2);
        dataTree->SetBranchAddress("Bp_VTXISONUMVTX_tau1", &Bp_VTXISONUMVTX_tau1);
        dataTree->SetBranchAddress("Bp_VTXISONUMVTX_tau2", &Bp_VTXISONUMVTX_tau2);
    }
    else
    {
        dataTree->SetBranchAddress("Bp_CONEMULT_tau1", &Bp_CONEMULT_tau1);
        dataTree->SetBranchAddress("Bp_CONEMULT_tau2", &Bp_CONEMULT_tau2);
        dataTree->SetBranchAddress("Bp_CONEPASYM_tau1", &Bp_CONEPASYM_tau1);
        dataTree->SetBranchAddress("Bp_CONEPASYM_tau2", &Bp_CONEPASYM_tau2);
        dataTree->SetBranchAddress("Bp_CONEPTASYM_tau1", &Bp_CONEPTASYM_tau1);
        dataTree->SetBranchAddress("Bp_CONEPTASYM_tau2", &Bp_CONEPTASYM_tau2);
        dataTree->SetBranchAddress("Bp_CONEDELTAETA_tau1", &Bp_CONEDELTAETA_tau1);
        dataTree->SetBranchAddress("Bp_CONEDELTAETA_tau2", &Bp_CONEDELTAETA_tau2);
    }
    dataTree->SetBranchAddress("Kp_PT", &Kp_PT);
    dataTree->SetBranchAddress("Kp_ProbNNk", &Kp_ProbNNk);
    dataTree->SetBranchAddress("Bp_ENDVERTEX_CHI2", &Bp_ENDVERTEX_CHI2);
    dataTree->SetBranchAddress("Bp_PT", &Bp_PT);
    dataTree->SetBranchAddress("Bp_ConsBp_seq_12_chi2", Bp_ConsBp_seq_12_chi2);
    dataTree->SetBranchAddress("Bp_OWNPV_X", &Bp_OWNPV_X);
    dataTree->SetBranchAddress("Bp_OWNPV_Y", &Bp_OWNPV_Y);
    dataTree->SetBranchAddress("Bp_OWNPV_Z", &Bp_OWNPV_Z);
    dataTree->SetBranchAddress("Bp_FDCHI2_OWNPV", &Bp_FDCHI2_OWNPV);

    Float_t BDT_response;
    tree->Branch("BDT_response", &BDT_response);

    int counter = 0;

    for(int i = 0; i < dataTree->GetEntries(); i++)
    {
        dataTree->GetEntry(i);

        if(year == 8)
        {
            var8[0] = TMath::Max(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2);
            var8[1] = TMath::Min(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2);
            var8[2] = TMath::Max(taum_M,taup_M);
            var8[3] = TMath::Max(taum_M12,taup_M12);
            var8[4] = TMath::Max(taum_M23,taup_M23);
            var8[5] = TMath::Min(taum_M,taup_M);
            var8[6] = TMath::Min(taum_M12,taup_M12);
            var8[7] = TMath::Min(taum_M23,taup_M23);
            var8[8] = TMath::Min( TMath::Log(TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))), 
                                TMath::Log(TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) );
            var8[9] = TMath::Log( TMath::Sqrt( pow(taup_ENDVERTEX_X-taum_ENDVERTEX_X,2) + 
                                            pow(taup_ENDVERTEX_Y-taum_ENDVERTEX_Y,2) + 
                                            pow(taup_ENDVERTEX_Z-taum_ENDVERTEX_Z,2) )  );
            var8[10] = TMath::Max(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX));
            var8[11] = TMath::Min(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX));
            var8[12] = TMath::Max(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), 
                                TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX));
            var8[13] = TMath::Min(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), 
                                TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX));
            Double_t Cx_p = (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY;
            Double_t Cy_p = (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ;
            Double_t Cz_p = (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX;
            Double_t C_p = TMath::Sqrt( pow(Cx_p,2) + pow(Cy_p,2) + pow(Cz_p,2) );
            Double_t IP_p = 2*C_p/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ));
            Double_t Cx_m = (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY;
            Double_t Cy_m = (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ;
            Double_t Cz_m = (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX;
            Double_t C_m = TMath::Sqrt( pow(Cx_m,2) + pow(Cy_m,2) + pow(Cz_m,2) );
            Double_t IP_m = 2*C_m/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ));
            var8[14] = TMath::Max(IP_p, IP_m);
            var8[15] = TMath::Min(IP_p, IP_m);
            Double_t a = TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) );
            Double_t b = TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) );
            Double_t c = TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) );
            Double_t s = (a+b+c)*0.5;
            Double_t A = TMath::Sqrt(s*(s-a)*(s-b)*(s-c));
            // var8[16] = A;
            var8[16] = TMath::Max(Bp_CC_IT_tau1,Bp_CC_IT_tau2);
            var8[17] = TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2);
            var8[18] = TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2);
            var8[19] = TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1);
            var8[20] = TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1);
            var8[21] = TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1);
            var8[22] = TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1);
            var8[23] = TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2);
            var8[24] = TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2);
            var8[25] = TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2);
            var8[26] = TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2);
            var8[27] = TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2);
            var8[28] = TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1);
            var8[29] = TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1);
            var8[30] = TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1);
            var8[31] = TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1);
            var8[32] = TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1);
            var8[33] = TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1,Bp_VTXISODCHI2ONETRACK_tau2);
            var8[34] = TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1,Bp_VTXISODCHI2ONETRACK_tau2);
            var8[35] = TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2,Bp_VTXISODCHI2TWOTRACK_tau1);
            var8[36] = TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2,Bp_VTXISODCHI2TWOTRACK_tau1);
            var8[37] = TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1);
            var8[38] = TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1);
            var8[39] = TMath::Log(1-Kp_ProbNNk);
            var8[40] = Bp_PT;
            var8[41] = Bp_ConsBp_seq_12_chi2[0];
        }
        else
        {
            var[0] = TMath::Max(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2);
            var[1] = TMath::Min(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2);
            var[2] = TMath::Max(taum_M,taup_M);
            var[3] = TMath::Max(taum_M12,taup_M12);
            var[4] = TMath::Max(taum_M23,taup_M23);
            var[5] = TMath::Min(taum_M,taup_M);
            var[6] = TMath::Min(taum_M12,taup_M12);
            var[7] = TMath::Min(taum_M23,taup_M23);
            var[8] = TMath::Min( TMath::Log(TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))), 
                                TMath::Log(TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) );
            var[9] = TMath::Log( TMath::Sqrt( pow(taup_ENDVERTEX_X-taum_ENDVERTEX_X,2) + 
                                            pow(taup_ENDVERTEX_Y-taum_ENDVERTEX_Y,2) + 
                                            pow(taup_ENDVERTEX_Z-taum_ENDVERTEX_Z,2) )  );
            var[10] = TMath::Max(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX));
            var[11] = TMath::Min(TMath::Log(taup_FDCHI2_ORIVX), TMath::Log(taum_FDCHI2_ORIVX));
            var[12] = TMath::Max(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), 
                                TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX));
            var[13] = TMath::Min(TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), 
                                TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX));
            Double_t Cx_p = (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY;
            Double_t Cy_p = (taup_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ;
            Double_t Cz_p = (taup_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taup_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX;
            Double_t C_p = TMath::Sqrt( pow(Cx_p,2) + pow(Cy_p,2) + pow(Cz_p,2) );
            Double_t IP_p = 2*C_p/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ));
            Double_t Cx_m = (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PZ - (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PY;
            Double_t Cy_m = (taum_ENDVERTEX_Z-Bp_ENDVERTEX_Z)*Kp_PX - (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PZ;
            Double_t Cz_m = (taum_ENDVERTEX_X-Bp_ENDVERTEX_X)*Kp_PY - (taum_ENDVERTEX_Y-Bp_ENDVERTEX_Y)*Kp_PX;
            Double_t C_m = TMath::Sqrt( pow(Cx_m,2) + pow(Cy_m,2) + pow(Cz_m,2) );
            Double_t IP_m = 2*C_m/(TMath::Sqrt( pow(Kp_PX,2) + pow(Kp_PY,2) + pow(Kp_PZ,2) ));
            var[14] = TMath::Max(IP_p, IP_m);
            var[15] = TMath::Min(IP_p, IP_m);
            Double_t a = TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) );
            Double_t b = TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) );
            Double_t c = TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) );
            Double_t s = (a+b+c)*0.5;
            Double_t A = TMath::Sqrt(s*(s-a)*(s-b)*(s-c));
            // var[16] = A;
            var[16] = TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2);
            var[17] = TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2);
            var[18] = TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2);
            var[19] = TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2);
            var[20] = TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2);
            var[21] = TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2);
            var[22] = TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2);
            var[23] = TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2);
            var[24] = TMath::Log(1-Kp_ProbNNk);
            var[25] = Bp_PT;
            var[26] = Bp_ConsBp_seq_12_chi2[0];
        }

        BDT_response=reader->EvaluateMVA("BDT method");

        tree->Fill();

        if ( i > 1.0*counter*(dataTree->GetEntries())/100 ) {
            cout<<counter<<"%"<<endl;
            counter += 10;
        }
    }

    tree->Write();
    target->Close();

    return;
}