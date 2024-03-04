
void TMVA_response(int year, int species, int tmva_method, int tmva_type, 
                   int background_proxy, TString MC_files, 
                   TString MC_0_files, TString MC_1_files, TString MC_2_files, 
                   TString RS_DATA_files, TString WS_DATA_files, TString weights, int line)
{
    // creates a root file containing the BDT response 

    TFileCollection* fc;
    TFile* target = new TFile(Form("/panfs/felician/B2Ktautau/workflow/TMVA_response/201%i/Species_%i/tmva_method%i_type%i_bkgProxy%i/%i.root",year,species,tmva_method,tmva_type,background_proxy,line), "RECREATE");

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

    int num_vars;
    if(tmva_type == 0)
    {
        num_vars = 17;
    }
    else if(tmva_type == 1)
    {
        if(year == 8)
        {
            num_vars = 31;
        }
        else
        {
            num_vars = 8;
        }
    }
    else if(tmva_type == 6)
    {
        num_vars = 34;
    }

    Float_t var[num_vars];

    // Add variables to the reader
    if(tmva_type == 0)
    {
        reader->AddVariable("TMath::Max( TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) )", &var[0]);
        reader->AddVariable("TMath::Min( TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) )", &var[1]);
            
        reader->AddVariable("TMath::Max(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", &var[2]);
        reader->AddVariable("TMath::Min(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2)", &var[3]);

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

        reader->AddVariable("TMath::Max("+IP_taup+", "+IP_taum+")", &var[4]);
        reader->AddVariable("TMath::Min("+IP_taup+", "+IP_taum+")", &var[5]);
        
        // Area of the triangle formed by the PV and the tau+ and tau- decay vertices
        TString a = "TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) )";
        TString b = "TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) )";
        TString c = "TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) )";
        TString s = "("+a+"+"+b+"+"+c+")*0.5";
        TString A = "TMath::Sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";
        reader->AddVariable("A := "+A, &var[6]);

        reader->AddVariable("TMath::Max(taum_M, taup_M)", &var[7]);
        reader->AddVariable("TMath::Max(taum_M12, taup_M12)", &var[8]);
        reader->AddVariable("TMath::Max(taum_M23, taup_M23)", &var[9]);
        reader->AddVariable("TMath::Min(taum_M, taup_M)", &var[10]);
        reader->AddVariable("TMath::Min(taum_M12, taup_M12)", &var[11]);
        reader->AddVariable("TMath::Min(taum_M23, taup_M23)", &var[12]);

        reader->AddVariable("TMath::Max(taup_PT, taum_PT)", &var[13]);
        reader->AddVariable("TMath::Min(taup_PT, taum_PT)", &var[14]);
        reader->AddVariable("TMath::Max(taup_PZ, taum_PZ)", &var[15]);
        reader->AddVariable("TMath::Min(taup_PZ, taum_PZ)", &var[16]);
    }
    else if(tmva_type == 1)
    {
        if(year == 8)
        {
            reader->AddVariable("TMath::Max(Bp_CC_IT_tau1, Bp_CC_IT_tau2)", &var[0]);
            reader->AddVariable("TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2)", &var[1]);
            reader->AddVariable("TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2)", &var[2]);
            reader->AddVariable("TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1)", &var[3]);
            reader->AddVariable("TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1)", &var[4]);
            reader->AddVariable("TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1)", &var[5]);
            reader->AddVariable("TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1)", &var[6]);
            reader->AddVariable("TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2)", &var[7]);
            reader->AddVariable("TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2)", &var[8]);
            reader->AddVariable("TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2)", &var[9]);
            reader->AddVariable("TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2)", &var[10]);
            reader->AddVariable("TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2)", &var[11]);
            reader->AddVariable("TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1)", &var[12]);
            reader->AddVariable("TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1)", &var[13]);
            reader->AddVariable("TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1)", &var[14]);
            reader->AddVariable("TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1)", &var[15]);
            reader->AddVariable("TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1)", &var[16]);

            reader->AddVariable("TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", &var[17]);
            reader->AddVariable("TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", &var[18]);
            reader->AddVariable("TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", &var[19]);
            reader->AddVariable("TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", &var[20]);
            reader->AddVariable("TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", &var[21]);
            reader->AddVariable("TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", &var[22]);
        
            reader->AddVariable("TMath::Max(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2)", &var[23]);
            reader->AddVariable("TMath::Max(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2)", &var[24]);
            reader->AddVariable("TMath::Max(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2)", &var[25]);
            reader->AddVariable("TMath::Max(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2)", &var[26]);

            reader->AddVariable("TMath::Min(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2)", &var[27]);
            reader->AddVariable("TMath::Min(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2)", &var[28]);
            reader->AddVariable("TMath::Min(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2)", &var[29]);
            reader->AddVariable("TMath::Min(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2)", &var[30]);
        }
        else
        {
            reader->AddVariable("TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", &var[0]);
            reader->AddVariable("TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", &var[1]);
            reader->AddVariable("TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", &var[2]);
            reader->AddVariable("TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", &var[3]);
            
            reader->AddVariable("TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", &var[4]);
            reader->AddVariable("TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", &var[5]);
            reader->AddVariable("TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", &var[6]);
            reader->AddVariable("TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", &var[7]);
        }
    }
    else if(tmva_type == 6)
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
        reader->AddVariable("TMath::Min( TMath::Log10(TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV)) ), "
                                        "TMath::Log10(TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) )", &var[8]);
        reader->AddVariable("TMath::Log10( TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + "
                                                    "pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + "
                                                    "pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) ) )", &var[9]);
        reader->AddVariable("TMath::Max(TMath::Log10(taup_FDCHI2_ORIVX), TMath::Log10(taum_FDCHI2_ORIVX))", &var[10]);
        reader->AddVariable("TMath::Min(TMath::Log10(taup_FDCHI2_ORIVX), TMath::Log10(taum_FDCHI2_ORIVX))", &var[11]);
        reader->AddVariable("TMath::Max(TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX),"
                                    " TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX))", &var[12]);
        reader->AddVariable("TMath::Min(TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX),"
                                    " TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX))", &var[13]);
        
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
        reader->AddVariable(A, &var[16]);
        
        // kaon
        reader->AddVariable("TMath::Log10(1-Kp_ProbNNk)", &var[17]);
        
        // B meson
        reader->AddVariable("Bp_PT", &var[18]);
        
        // others
        reader->AddVariable("Bp_dtf_12_chi2[0]", &var[19]);
        
        if(year == 8) // tau isolation
        {
            // reader->AddVariable("TMath::Max(Bp_CC_IT_tau1, Bp_CC_IT_tau2)", &var[20]);
            // reader->AddVariable("TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2)", &var[21]);
            // reader->AddVariable("TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2)", &var[22]);
            // reader->AddVariable("TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1)", &var[23]);
            // reader->AddVariable("TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1)", &var[24]);
            // reader->AddVariable("TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1)", &var[25]);
            // reader->AddVariable("TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1)", &var[26]);
            // reader->AddVariable("TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2)", &var[27]);
            // reader->AddVariable("TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2)", &var[28]);
            // reader->AddVariable("TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2)", &var[29]);
            // reader->AddVariable("TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2)", &var[30]);
            // reader->AddVariable("TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2)", &var[31]);
            // reader->AddVariable("TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1)", &var[32]);
            // reader->AddVariable("TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1)", &var[33]);
            // reader->AddVariable("TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1)", &var[34]);
            // reader->AddVariable("TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1)", &var[35]);
            // reader->AddVariable("TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1)", &var[36]);
            reader->AddVariable("TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", &var[20]);
            reader->AddVariable("TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1, Bp_VTXISODCHI2ONETRACK_tau2)", &var[21]);
            reader->AddVariable("TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", &var[22]);
            reader->AddVariable("TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2, Bp_VTXISODCHI2TWOTRACK_tau1)", &var[23]);
            reader->AddVariable("TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", &var[24]);
            reader->AddVariable("TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1)", &var[25]);
            reader->AddVariable("TMath::Max(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2)", &var[26]);
            reader->AddVariable("TMath::Max(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2)", &var[27]);
            reader->AddVariable("TMath::Max(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2)", &var[28]);
            reader->AddVariable("TMath::Max(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2)", &var[29]);
            reader->AddVariable("TMath::Min(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2)", &var[30]);
            reader->AddVariable("TMath::Min(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2)", &var[31]);
            reader->AddVariable("TMath::Min(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2)", &var[32]);
            reader->AddVariable("TMath::Min(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2)", &var[33]);
        }
        else
        {
            reader->AddVariable("TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", &var[20]);
            reader->AddVariable("TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", &var[21]);
            reader->AddVariable("TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", &var[22]);
            reader->AddVariable("TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", &var[23]);
            reader->AddVariable("TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2)", &var[24]);
            reader->AddVariable("TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2)", &var[25]);
            reader->AddVariable("TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2)", &var[26]);
            reader->AddVariable("TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2)", &var[27]);
        }
    }
    reader->BookMVA("TMVA method", weights);  

    // Set branch address
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

    Double_t Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2, Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2, Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2;

    Double_t Bp_CONEMULT_tau1, Bp_CONEMULT_tau2, Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2, Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2, Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2;

    Float_t Bp_dtf_12_chi2[100];
    Double_t Bp_OWNPV_X, Bp_OWNPV_Y, Bp_OWNPV_Z, Bp_FDCHI2_OWNPV;
    Double_t taup_PT, taup_PZ, taum_PT, taum_PZ;
    
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
        dataTree->SetBranchAddress("Bp_CC_PASYM_tau1", &Bp_CC_PASYM_tau1);
        dataTree->SetBranchAddress("Bp_CC_PASYM_tau2", &Bp_CC_PASYM_tau2);
        dataTree->SetBranchAddress("Bp_CC_PTASYM_tau1", &Bp_CC_PTASYM_tau1);
        dataTree->SetBranchAddress("Bp_CC_PTASYM_tau2", &Bp_CC_PTASYM_tau2);
        dataTree->SetBranchAddress("Bp_CC_DELTAETA_tau1", &Bp_CC_DELTAETA_tau1);
        dataTree->SetBranchAddress("Bp_CC_DELTAETA_tau2", &Bp_CC_DELTAETA_tau2);
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
    dataTree->SetBranchAddress("Bp_dtf_12_chi2", Bp_dtf_12_chi2);
    dataTree->SetBranchAddress("Bp_OWNPV_X", &Bp_OWNPV_X);
    dataTree->SetBranchAddress("Bp_OWNPV_Y", &Bp_OWNPV_Y);
    dataTree->SetBranchAddress("Bp_OWNPV_Z", &Bp_OWNPV_Z);
    dataTree->SetBranchAddress("Bp_FDCHI2_OWNPV", &Bp_FDCHI2_OWNPV);
    dataTree->SetBranchAddress("taup_PT", &taup_PT);
    dataTree->SetBranchAddress("taup_PZ", &taup_PZ);
    dataTree->SetBranchAddress("taum_PT", &taum_PT);
    dataTree->SetBranchAddress("taum_PZ", &taum_PZ);

    Float_t TMVA_response = 0;
    tree->Branch("BDT_response", &TMVA_response);

    int counter = 0;

    for(int i = 0; i < dataTree->GetEntries(); i++)
    {
        dataTree->GetEntry(i);

        if(tmva_type == 0)
        {
            var[0] = TMath::Max( TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) );
            var[1] = TMath::Min( TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) );
            var[2] = TMath::Max(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2);
            var[3] = TMath::Min(taum_ENDVERTEX_CHI2, taup_ENDVERTEX_CHI2);
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
            var[4] = TMath::Max(IP_p, IP_m);
            var[5] = TMath::Min(IP_p, IP_m);
            Double_t a = TMath::Sqrt( pow(Bp_OWNPV_X - taup_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taup_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taup_ENDVERTEX_Z,2) );
            Double_t b = TMath::Sqrt( pow(Bp_OWNPV_X - taum_ENDVERTEX_X,2) + pow(Bp_OWNPV_Y - taum_ENDVERTEX_Y,2) + pow(Bp_OWNPV_Z - taum_ENDVERTEX_Z,2) );
            Double_t c = TMath::Sqrt( pow(taup_ENDVERTEX_X - taum_ENDVERTEX_X,2) + pow(taup_ENDVERTEX_Y - taum_ENDVERTEX_Y,2) + pow(taup_ENDVERTEX_Z - taum_ENDVERTEX_Z,2) );
            Double_t s = (a+b+c)*0.5;
            Double_t A = TMath::Sqrt(s*(s-a)*(s-b)*(s-c));
            var[6] = A;
            var[7] = TMath::Max(taum_M, taup_M);
            var[8] = TMath::Max(taum_M12, taup_M12);
            var[9] = TMath::Max(taum_M23, taup_M23);
            var[10] = TMath::Min(taum_M, taup_M);
            var[11] = TMath::Min(taum_M12, taup_M12);
            var[12] = TMath::Min(taum_M23, taup_M23);
            var[13] = TMath::Max(taup_PT, taum_PT);
            var[14] = TMath::Min(taup_PT, taum_PT);
            var[15] = TMath::Max(taup_PZ, taum_PZ);
            var[16] = TMath::Min(taup_PZ, taum_PZ);
        }
        else if(tmva_type == 1)
        {
            if(year == 8)
            {
                var[0] = TMath::Max(Bp_CC_IT_tau1,Bp_CC_IT_tau2);
                var[1] = TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2);
                var[2] = TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2);
                var[3] = TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1);
                var[4] = TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1);
                var[5] = TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1);
                var[6] = TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1);
                var[7] = TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2);
                var[8] = TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2);
                var[9] = TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2);
                var[10] = TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2);
                var[11] = TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2);
                var[12] = TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1);
                var[13] = TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1);
                var[14] = TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1);
                var[15] = TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1);
                var[16] = TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1);
                var[17] = TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1,Bp_VTXISODCHI2ONETRACK_tau2);
                var[18] = TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1,Bp_VTXISODCHI2ONETRACK_tau2);
                var[19] = TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2,Bp_VTXISODCHI2TWOTRACK_tau1);
                var[20] = TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2,Bp_VTXISODCHI2TWOTRACK_tau1);
                var[21] = TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1);
                var[22] = TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1);
                var[23] = TMath::Max(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2);
                var[24] = TMath::Max(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2);
                var[25] = TMath::Max(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2);
                var[26] = TMath::Max(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2);
                var[27] = TMath::Min(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2);
                var[28] = TMath::Min(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2);
                var[29] = TMath::Min(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2);
                var[30] = TMath::Min(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2);

            }
            else
            {
                var[0] = TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2);
                var[1] = TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2);
                var[2] = TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2);
                var[3] = TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2);
                var[4] = TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2);
                var[5] = TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2);
                var[6] = TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2);
                var[7] = TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2);
            }
        }
        else if(tmva_type == 6)
        {
            var[0] = TMath::Max(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2);
            var[1] = TMath::Min(taum_ENDVERTEX_CHI2,taup_ENDVERTEX_CHI2);
            var[2] = TMath::Max(taum_M,taup_M);
            var[3] = TMath::Max(taum_M12,taup_M12);
            var[4] = TMath::Max(taum_M23,taup_M23);
            var[5] = TMath::Min(taum_M,taup_M);
            var[6] = TMath::Min(taum_M12,taup_M12);
            var[7] = TMath::Min(taum_M23,taup_M23);
            var[8] = TMath::Min( TMath::Log10(TMath::Min(taum_pi1_IPCHI2_OWNPV, TMath::Min(taum_pi2_IPCHI2_OWNPV, taum_pi3_IPCHI2_OWNPV))), 
                                TMath::Log10(TMath::Min(taup_pi1_IPCHI2_OWNPV, TMath::Min(taup_pi2_IPCHI2_OWNPV, taup_pi3_IPCHI2_OWNPV))) );
            var[9] = TMath::Log10( TMath::Sqrt( pow(taup_ENDVERTEX_X-taum_ENDVERTEX_X,2) + 
                                            pow(taup_ENDVERTEX_Y-taum_ENDVERTEX_Y,2) + 
                                            pow(taup_ENDVERTEX_Z-taum_ENDVERTEX_Z,2) )  );
            var[10] = TMath::Max(TMath::Log10(taup_FDCHI2_ORIVX), TMath::Log10(taum_FDCHI2_ORIVX));
            var[11] = TMath::Min(TMath::Log10(taup_FDCHI2_ORIVX), TMath::Log10(taum_FDCHI2_ORIVX));
            var[12] = TMath::Max(TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), 
                                TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX));
            var[13] = TMath::Min(TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), 
                                TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX));
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
            var[16] = A;
            var[17] = TMath::Log10(1-Kp_ProbNNk);
            var[18] = Bp_PT;
            var[19] = Bp_dtf_12_chi2[0];
            if(year == 8)
            {
                // var[20] = TMath::Max(Bp_CC_IT_tau1,Bp_CC_IT_tau2);
                // var[21] = TMath::Max(Bp_CC_MAXPT_PT_tau1,Bp_CC_MAXPT_PT_tau2);
                // var[22] = TMath::Max(Bp_CC_MULT_tau1,Bp_CC_MULT_tau2);
                // var[23] = TMath::Min(Bp_CC_PZ_tau2,Bp_CC_PZ_tau1);
                // var[24] = TMath::Min(Bp_CC_PZASYM_tau2,Bp_CC_PZASYM_tau1);
                // var[25] = TMath::Min(Bp_CC_SPT_tau2,Bp_CC_SPT_tau1);
                // var[26] = TMath::Min(Bp_CC_VPT_tau2,Bp_CC_VPT_tau1);
                // var[27] = TMath::Max(Bp_NC_DELTAETA_tau1,Bp_NC_DELTAETA_tau2);
                // var[28] = TMath::Max(Bp_NC_DELTAPHI_tau1,Bp_NC_DELTAPHI_tau2);
                // var[29] = TMath::Max(Bp_NC_IT_tau1,Bp_NC_IT_tau2);
                // var[30] = TMath::Max(Bp_NC_MAXPT_PT_tau1,Bp_NC_MAXPT_PT_tau2);
                // var[31] = TMath::Max(Bp_NC_MAXPT_PZ_tau1,Bp_NC_MAXPT_PZ_tau2);
                // var[32] = TMath::Min(Bp_NC_MULT_tau2,Bp_NC_MULT_tau1);
                // var[33] = TMath::Min(Bp_NC_PTASYM_tau2,Bp_NC_PTASYM_tau1);
                // var[34] = TMath::Min(Bp_NC_PZ_tau2,Bp_NC_PZ_tau1);
                // var[35] = TMath::Min(Bp_NC_SPT_tau2,Bp_NC_SPT_tau1);
                // var[36] = TMath::Min(Bp_NC_VPT_tau2,Bp_NC_VPT_tau1);
                var[20] = TMath::Max(Bp_VTXISODCHI2ONETRACK_tau1,Bp_VTXISODCHI2ONETRACK_tau2);
                var[21] = TMath::Min(Bp_VTXISODCHI2ONETRACK_tau1,Bp_VTXISODCHI2ONETRACK_tau2);
                var[22] = TMath::Max(Bp_VTXISODCHI2TWOTRACK_tau2,Bp_VTXISODCHI2TWOTRACK_tau1);
                var[23] = TMath::Min(Bp_VTXISODCHI2TWOTRACK_tau2,Bp_VTXISODCHI2TWOTRACK_tau1);
                var[24] = TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1);
                var[25] = TMath::Min(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1);
                var[26] = TMath::Max(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2);
                var[27] = TMath::Max(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2);
                var[28] = TMath::Max(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2);
                var[29] = TMath::Max(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2);
                var[30] = TMath::Min(Bp_CC_MULT_tau1, Bp_CC_MULT_tau2);
                var[31] = TMath::Min(Bp_CC_PASYM_tau1, Bp_CC_PASYM_tau2);
                var[32] = TMath::Min(Bp_CC_PTASYM_tau1, Bp_CC_PTASYM_tau2);
                var[33] = TMath::Min(Bp_CC_DELTAETA_tau1, Bp_CC_DELTAETA_tau2);

            }
            else
            {
                var[20] = TMath::Max(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2);
                var[21] = TMath::Max(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2);
                var[22] = TMath::Max(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2);
                var[23] = TMath::Max(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2);
                var[24] = TMath::Min(Bp_CONEMULT_tau1, Bp_CONEMULT_tau2);
                var[25] = TMath::Min(Bp_CONEPASYM_tau1, Bp_CONEPASYM_tau2);
                var[26] = TMath::Min(Bp_CONEPTASYM_tau1, Bp_CONEPTASYM_tau2);
                var[27] = TMath::Min(Bp_CONEDELTAETA_tau1, Bp_CONEDELTAETA_tau2);
            }
        }

        TMVA_response=reader->EvaluateMVA("TMVA method");
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