
Double_t eps_error(Double_t Num, Double_t Den);
TCut DD_cut(TString name, Int_t mother_ID);
TCut truthMatch_cocktailMC(TString name);

void create_pre_sel_tree(int year, int species, int line, bool createTable)
{   
    // TString path = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/Num_entries/%i.txt",year,species,line);
    // std::ofstream file(path);
    Bool_t isKtautauMC = false;
    Bool_t isBuDDKp_cocktail = false;
    Bool_t isBdDDKp_cocktail = false;
    Bool_t isBsDDKp_cocktail = false;
    Bool_t isBuDDK0_cocktail = false;
    Bool_t isBdDDK0_cocktail = false;
    Bool_t isBuDD_cocktail = false;
    Bool_t isBdDD_cocktail = false;
    Bool_t isBsDD_cocktail = false;

    if( (species == 10) || (species == 11) || (species == 12) || (species == 1) )
    {
        isKtautauMC = true;
    }
    if( (species == 100) || (species == 101) || (species == 102) )
    {
        isBuDDKp_cocktail = true;
    }
    if(species == 110)
    {
        isBdDDKp_cocktail = true;
    }
    if(species == 120)
    {
        isBsDDKp_cocktail = true;
    }
    if(species == 130)
    {
        isBuDDK0_cocktail = true;
    }
    if( (species == 140) || (species == 141) )
    {
        isBdDDK0_cocktail = true;
    }
    if( (species == 150) || (species == 151) )
    {
        isBuDD_cocktail = true;
    }
    if( (species == 160) || (species == 161) || (species == 162) || (species == 163) )
    {
        isBdDD_cocktail = true;
    }
    if( (species == 170) || (species == 171) || (species == 172) || (species == 173) )
    {
        isBsDD_cocktail = true;
    }

    Bool_t is_cocktailMC = false;
    if(isBuDDKp_cocktail || isBdDDKp_cocktail || isBsDDKp_cocktail || isBuDDK0_cocktail || isBdDDK0_cocktail || isBuDD_cocktail || isBdDD_cocktail || isBsDD_cocktail)
    {
        is_cocktailMC = true;
    }

    TString FILES;
    if(isKtautauMC) // Ktautau MC
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_1/pid_corr.txt",year);
    }
    else if(species == 4) // DDK MC
    {
        FILES = Form("Files_on_grid/MC_DpDmK_201%i.txt",year);
    }
    else if((species == 2) || (species == 3)) // Ktautau data
    {
        FILES = Form("Files_on_grid/data_201%i.txt",year);
    }
    else if((species == 5) || (species == 6)) // D+D-K data
    {
        FILES = Form("Files_on_grid/data_DDK_201%i.txt",year);
    }
    else if( (species == 7) || (species == 71) || (species == 74) )
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_7/pid_corr.txt",year);
    }
    else if((species == 8) || (species == 81) || (species == 82) || (species == 83) || (species == 84) || (species == 85))
    {
        FILES = Form("Files_on_grid/data_D0Dps_201%i.txt",year);
    }
    else if(species == 9) // D0D0K MC
    {
        FILES = Form("Files_on_grid/MC_D0D0K_201%i.txt",year);
    }
    else if((species == 0) || (species == -1)) // D0D0K data
    {
        FILES = Form("Files_on_grid/data_D0D0K_201%i.txt",year);
    }
    else if(isBuDDKp_cocktail)
    {
        FILES = Form("Files_on_grid/MC_201%i_BuDDKp_cocktail.txt",year);
    }
    else if(isBdDDKp_cocktail)
    {
        FILES = Form("Files_on_grid/MC_201%i_BdDDKp_cocktail.txt",year);   
    }
    else if(isBsDDKp_cocktail)
    {
        FILES = Form("Files_on_grid/MC_201%i_BsDDKp_cocktail.txt",year);  
    }
    else if(isBuDDK0_cocktail)
    {
        FILES = Form("Files_on_grid/MC_201%i_BuDDK0_cocktail.txt",year);  
    }
    else if(isBdDDK0_cocktail)
    {
        FILES = Form("/panfs/felician/BdDDK0_cocktail/root_files_201%i.txt",year); 
    }
    else if(isBuDD_cocktail)
    {
        FILES = Form("Files_on_grid/MC_201%i_BuDD_cocktail.txt",year);
    }
    else if(isBdDD_cocktail)
    {
        FILES = Form("/panfs/felician/BdDD_cocktail/root_files_201%i.txt",year);
    }
    else if(isBsDD_cocktail)
    {
        FILES = Form("/panfs/felician/BsDD_cocktail/root_files_201%i.txt",year);
    }

    TString MC_component_file;
    if(isKtautauMC)
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/%i.root",year,line);
    }

    /////////////////////////////////////////////////// TRUTH MATCH CUTS ///////////////////////////////////////////////////////////////////////
    // Other MC:
    TCut truthMatch;
    if(isKtautauMC){
        truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
    } 
    else if(species == 4) // D+D-K+ MC
    {
        truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(Dp_K_TRUEID) == 321) && (abs(Dp_pi1_TRUEID) == 211) && (abs(Dp_pi2_TRUEID) == 211) && (abs(Dm_K_TRUEID) == 321) && (abs(Dm_pi1_TRUEID) == 211) && (abs(Dm_pi2_TRUEID) == 211) && (abs(Dp_TRUEID) == 411) && (abs(Dm_TRUEID) == 411) && (abs(Bp_TRUEID) == 521)";
    }
    else if((species == 7) || (species == 74) ) // D0D+s MC
    {
        truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(Dsp_K1_TRUEID) == 321) && (abs(Dsp_K2_TRUEID) == 321) && (abs(D0bar_pi_TRUEID) == 211) && (abs(Dsp_pi_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(Dsp_TRUEID) == 431) && (abs(Bp_TRUEID) == 521)";
    }
    else if(species == 9) // D0D0K MC
    {
        truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(D0bar_pi1_TRUEID) == 211) && (abs(D0bar_pi2_TRUEID) == 211) && (abs(D0bar_pi3_TRUEID) == 211) && (abs(D0_K_TRUEID) == 321) && (abs(D0_pi1_TRUEID) == 211) && (abs(D0_pi2_TRUEID) == 211) && (abs(D0_pi3_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(D0_TRUEID) == 421) && (abs(Kp_TRUEID) == 321) && (abs(Bp_TRUEID) == 521)";
    }

    ///////////////////////////////////////////////////////// TRIGGER ////////////////////////////////////////////////////////////////////////////////////////
    TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
    TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
    TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
    TCut trigger = L0_trigger+HLT1_trigger+HLT2_trigger;

    ///////////////////////////////////////////////////////// Rectangular cuts /////////////////////////////////////////////////////////////////////////////////
    TCut others;
    std::vector<TCut> other_cuts;
    if(isKtautauMC || (species == 2) || (species == 3) || is_cocktailMC)
    {
        other_cuts.push_back("(taup_M > 750) && (taup_M < 1650)");
        other_cuts.push_back("(taum_M > 750) && (taum_M < 1650)");
        other_cuts.push_back("(Bp_VTXISODCHI2MASSONETRACK_B > 3600)");
        other_cuts.push_back("(Bp_VTXISOBDTHARDFIRSTVALUE_B < 0)");
        other_cuts.push_back("(Bp_BPVVD > 4)");
        other_cuts.push_back("(TMath::Max(Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup,Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum) > -0.1)");
    }
    if((species == 4) || (species == 5) || (species == 6)) // D+D-K+
    {   
        other_cuts.push_back("(Dp_K_ProbNNk > 0.1) && (Dm_K_ProbNNk > 0.1)");
        other_cuts.push_back("(TMath::Min( TMath::Log10(1-TMath::Abs(Dp_DIRA_ORIVX))*TMath::Sign(1,Dp_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(Dm_DIRA_ORIVX))*TMath::Sign(1,Dm_DIRA_ORIVX) ) < -2.)"); 
        other_cuts.push_back("(TMath::Log10(TMath::Abs(1-Bp_DIRA_OWNPV))*TMath::Sign(1,Bp_DIRA_OWNPV) < -4)");
        other_cuts.push_back("(Bp_ENDVERTEX_CHI2 < 40)");
        other_cuts.push_back("(abs(Dp_M-1869.66) < 50) && (abs(Dm_M < 1869.66) < 50)");
        other_cuts.push_back("(Bp_FD_OWNPV > 2)");
    }
    if((species == 7) || (species == 8) ) // D0barD+s
    { 
        other_cuts.push_back("Bp_FD_OWNPV < 80");
        other_cuts.push_back("(D0bar_M > 1820) && (D0bar_M < 1910)");
        other_cuts.push_back("(Dsp_M > 1930) && (Dsp_M < 2000)");
        if(species == 7)
        {
            other_cuts.push_back("D0bar_K_ProbNNk_pidgen_default > 0.55");
            other_cuts.push_back("Dsp_K1_ProbNNk_pidgen_default > 0.55");
            other_cuts.push_back("Dsp_K2_ProbNNk_pidgen_default > 0.55");
        }
        else
        {
            other_cuts.push_back("D0bar_K_ProbNNk > 0.55");
            other_cuts.push_back("Dsp_K1_ProbNNk > 0.55");
            other_cuts.push_back("Dsp_K2_ProbNNk > 0.55");
        }

    }
    if(species == 81) // D0bar (K+K-pi+)
    {
        other_cuts.push_back("Bp_FD_OWNPV < 80");
        other_cuts.push_back("D0bar_K_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K1_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K2_ProbNNk > 0.55");
        other_cuts.push_back("(D0bar_M > 1820) && (D0bar_M < 1910)");
        other_cuts.push_back("(Dsp_M < 1930) || (Dsp_M > 2000)"); // D+s sidebands
    }
    if(species == 82) // (K+pi-) D+s
    {
        other_cuts.push_back("Bp_FD_OWNPV < 80");
        other_cuts.push_back("D0bar_K_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K1_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K2_ProbNNk > 0.55");
        other_cuts.push_back("(D0bar_M < 1820) || (D0bar_M > 1910)");
        other_cuts.push_back("(Dsp_M > 1930) && (Dsp_M < 2000)");
    }
    if(species == 83)
    {
        other_cuts.push_back("Bp_FD_OWNPV < 80");
        other_cuts.push_back("D0bar_K_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K1_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K2_ProbNNk > 0.55");
        other_cuts.push_back("(D0bar_M < 1820) || (D0bar_M > 1910)");
        other_cuts.push_back("(Dsp_M < 1930) || (Dsp_M > 2000)");
    }
    if((species == 9) || (species == 0) || (species == -1)) // D0D0K MC, RS data, WS data
    {
        other_cuts.push_back("Kp_PT > 1000");
        other_cuts.push_back("Kp_ProbNNk > 0.55");
        other_cuts.push_back("D0bar_K_ProbNNk > 0.55");
        other_cuts.push_back("D0_K_ProbNNk > 0.55");
        other_cuts.push_back("(D0bar_M > 1840) && (D0bar_M < 1890)");
        other_cuts.push_back("(D0_M > 1840) && (D0_M < 1890)");
    }

    Int_t N = other_cuts.size();
    for(int i = 0; i < N; i++)
    {
        others += other_cuts[i];
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFileCollection* fc = new TFileCollection("fc", "fc", FILES, 1, line);
    TChain* t;
    if(isKtautauMC || (species == 7) || (species == 74) ) // Ktautau MC + normalisation channel (PID corrected)
    {
        t = new TChain("DecayTree");
    }
    else if((species == 3) || (species == 6) || (species == -1)) // WS data (Ktautau,D+D-K+,D0D0K)
    {
        t = new TChain("ntuple_SS/DecayTree");
    }
    else // MC and RS data
    {
        t = new TChain("ntuple/DecayTree");
    }
    t->AddFileInfoList((TCollection*)fc->GetList());

    TCut MC_component = "";
    if(isKtautauMC)
    {
        t->AddFriend("DecayTree", MC_component_file);
    
        if(species == 10){MC_component = "component == 0";} // 3pi 3pi
        else if(species == 11){MC_component = "component == 1";} // 3pi 3pi pi0
        else if(species == 12){MC_component = "component == 2";} // 3pi 3pi 2pi0
        else if(species == 1){MC_component = "";} // all 3 components
    }

    TCut pre_selections;
    if(isKtautauMC) // Ktautau MC
    {
        pre_selections = MC_component+truthMatch+trigger+others;
    }
    else if((species == 2) || (species == 3)) // Ktautau data
    {
        pre_selections = trigger+others;
    }
    else if((isBuDDKp_cocktail) || (isBdDDKp_cocktail) || (isBsDDKp_cocktail) || (isBuDDK0_cocktail) || (isBdDDK0_cocktail) || (isBuDD_cocktail) || (isBdDD_cocktail) || (isBsDD_cocktail))
    {
        pre_selections = truthMatch+trigger+others;
    }
    else if((species == 4) || (species == 7) || (species == 9) ) // D+D-K+, D0barD+s, D0D0K MC
    {
        pre_selections = truthMatch+trigger+others+"(Bp_dtf_status[0]==0)";
    }
    else if(species == 74)
    {
        pre_selections = truthMatch+trigger+"(Bp_dtf_status[0]==0)";
    }
    else if((species == 5) || (species == 6) || (species == 8) || (species == 81) || (species == 82) || (species == 83) || (species == 0) || (species == -1)) // D+D-K+ and D0barD+s data
    {
        pre_selections = trigger+others+"(Bp_dtf_status[0]==0)";
    }
    else if(species == 84)
    {
        pre_selections = trigger+"(Bp_dtf_status[0]==0)";
    }
    else if(species == 85)
    {
        pre_selections = "";
    }
    else if(is_cocktailMC)
    {
        pre_selections = trigger+others;
    }

    TString fout_name = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/%i.root",year,species,line);
    TFile* fout = new TFile(fout_name, "RECREATE");   
    fout->cd();
    if(is_cocktailMC)
    {
        TString names1[] = {"BuD0D0Kp", "BuD0starD0Kp", "BuD0D0starKp", "BuD0starD0starKp"}; // B+ -> D0 D0 K+
        TString names2[] = {"BuDpDmKp", "BuDpstarDmKp", "BuDpDmstarKp", "BuDpstarDmstarKp"}; // B+ -> D+ D- K+
        TString names3[] = {"BuDsDsKp", "BuDsstarDsKp", "BuDsDsstarKp", "BuDsstarDsstarKp"}; // B+ -> Ds+ Ds- K+
        TString names4[] = {"BdDmD0Kp", "BdDmstarD0Kp", "BdDmD0starKp", "BdDmstarD0starKp"}; // B0 -> D- D0 K+
        TString names5[] = {"BsDsD0Kp", "BsDsstarD0Kp", "BsDsD0starKp", "BsDsstarD0starKp"}; // Bs -> Ds- D0 K+
        TString names6[] = {"BuD0DpK0", "BuD0starDpK0", "BuD0DpstarK0", "BuD0starDpstarK0"}; // B+ -> D0 D+ K0
        TString names7[] = {"BdDpDmK0", "BdDpstarDmK0", "BdDpDmstarK0", "BdDpstarDmstarK0"}; // B0 -> D+ D- K0
        TString names8[] = {"BdD0D0K0", "BdD0starD0K0", "BdD0D0starK0", "BdD0starD0starK0"}; // B0 -> D0 D0 K0
        TString names9[] = {"BuD0Ds", "BuD0starDs", "BuD0Dsstar", "BuD0starDsstar"}; // B+ -> D0 Ds+
        TString names10[] = {"BuD0Dp", "BuD0starDp", "BuD0Dpstar", "BuD0starDpstar"}; // B+ -> D0 D+
        TString names11[] = {"BdD0D0", "BdD0starD0", "BdD0D0star", "BdD0starD0star"}; // B0 -> D0 D0
        TString names12[] = {"BdDpDm", "BdDpstarDm", "BdDpDmstar", "BdDpstarDmstar"}; // B0 -> D+ D-
        TString names13[] = {"BdDpDs", "BdDpstarDs", "BdDpDsstar", "BdDpstarDsstar"}; // B0 -> D- Ds+
        TString names14[] = {"BdDsDs", "BdDsstarDs", "BdDsDsstar", "BdDsstarDsstar"}; // B0 -> Ds+ Ds-
        TString names15[] = {"BsDsDs", "BsDsstarDs", "BsDsDsstar", "BsDsstarDsstar"}; // Bs -> Ds+ Ds-
        TString names16[] = {"BsDpDs", "BsDpstarDs", "BsDpDsstar", "BsDpstarDsstar"}; // Bs -> D- Ds+
        TString names17[] = {"BsDpDm", "BsDpstarDm", "BsDpDmstar", "BsDpstarDmstar"}; // Bs -> D+ D-
        TString names18[] = {"BsD0D0", "BsD0starD0", "BsD0D0star", "BsD0starD0star"}; // Bs -> D0 D0

        for(int i = 0; i < 4; i++)
        {
            TString name;
            if(species == 100){name = names1[i];} // B+ -> D0 D0 K+
            else if(species == 101){name = names2[i];} // B+ -> D+ D- K+
            else if(species == 102){name = names3[i];} // B+ -> Ds+ Ds- K+
            else if(species == 110){name = names4[i];} // B0 -> D- D0 K+
            else if(species == 120){name = names5[i];} // Bs -> Ds- D0 K+
            else if(species == 130){name = names6[i];} // B+ -> D0 D+ K0
            else if(species == 140){name = names7[i];} // B0 -> D+ D- K0
            else if(species == 141){name = names8[i];} // B0 -> D0 D0 K0
            else if(species == 150){name = names9[i];} // B+ -> D0 Ds+
            else if(species == 151){name = names10[i];} // B+ -> D0 D+
            else if(species == 160){name = names11[i];} // B0 -> D0 D0
            else if(species == 161){name = names12[i];} // B0 -> D+ D-
            else if(species == 162){name = names13[i];} // B0 -> D- Ds+
            else if(species == 163){name = names14[i];} // B0 -> Ds+ Ds-
            else if(species == 170){name = names15[i];} // Bs -> Ds+ Ds-
            else if(species == 171){name = names16[i];} // Bs -> D- Ds+
            else if(species == 172){name = names17[i];} // Bs -> D+ D-
            else if(species == 173){name = names18[i];} // Bs -> D0 D0

            fout->mkdir(name);
            TTree* t_pre_sel = (TTree*)t->CopyTree(pre_selections+truthMatch_cocktailMC(name));
            fout->cd(name);
            t_pre_sel->Write();
        }
    }
    else
    {
        TTree* t_pre_sel = (TTree*)t->CopyTree(pre_selections);
        t_pre_sel->Write();
    }
    fout->Close();

    cout << "Finished successfully" << endl;

    // if(createTable) // creates the pre-selection efficiency tables
    // {
    //     cout << "Creating pre-selection efficiency tables" << endl;
    //     Bool_t isMC = false;
    //     if((species == 1) || (species == 10) || (species == 11) || (species == 12) || (species == 4) || (species == 7) || (species == 9))
    //     {
    //         isMC = true;
    //     }

    //     TCut mass = ""; // background definition
    //     // if(species == 8) 
    //     // {
    //     //     mass = "Bp_dtf_M[0] > 5350";

    //     // }

    //     TFileCollection* fc_all;
    //     if(isMC)
    //     {
    //         fc_all = new TFileCollection("fc_all", "fc_all", FILES);
    //     }
    //     else
    //     {
    //         fc_all = new TFileCollection("fc_all", "fc_all", FILES);
    //     }
         
    //     TChain* t_reco;
    //     if((species == 3) || (species == 6) || (species == -1)) // WS data (Ktautau,D+D-K+,D0D0K)
    //     {
    //         t_reco = new TChain("ntuple_SS/DecayTree");
    //     }
    //     else // MC and RS data
    //     {
    //         t_reco = new TChain("ntuple/DecayTree");
    //     }
    //     t_reco->AddFileInfoList((TCollection*)fc_all->GetList());

    //     // Reconstruction efficiency
    //     cout << "Reconstruction efficiency" << endl;
    //     TChain* t_gen;
    //     TChain* t_gen1;
    //     TChain* t_gen2;
    //     TChain* t_gen3;
    //     Double_t eps_reco = 1.;
    //     Double_t N_reco = t_reco->GetEntries(mass);
    //     Double_t N_gen;
    //     if(isMC)
    //     {
    //         if(species == 10)
    //         {
    //             t_gen = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
    //         }
    //         else if(species == 11)
    //         {
    //             t_gen = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
    //             t_gen1 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
    //         }
    //         else if(species == 12)
    //         {
    //             t_gen = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
    //         }
    //         else if(species == 1)
    //         {
    //             t_gen1 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
    //             t_gen2 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
    //             t_gen3 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
    //             t_gen = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
    //         }
    //         else
    //         {
    //             t_gen = new TChain("mc_ntuple/MCDecayTree");
    //         }
    //         t_gen->AddFileInfoList((TCollection*)fc_all->GetList());

    //         if(species == 11)
    //         {
    //             t_gen1->AddFileInfoList((TCollection*)fc_all->GetList());
    //         }
    //         if(species == 1)
    //         {
    //             t_gen1->AddFileInfoList((TCollection*)fc_all->GetList());
    //             t_gen2->AddFileInfoList((TCollection*)fc_all->GetList());
    //             t_gen3->AddFileInfoList((TCollection*)fc_all->GetList());
    //         }

    //         if(isKtautauMC)
    //         {
    //             TFileCollection* fc_comp = new TFileCollection("fc_comp", "fc_comp", Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/mc_components.txt",year));
    //             TChain* t_comp = new TChain("DecayTree");
    //             t_comp->AddFileInfoList((TCollection*)fc_comp->GetList());

    //             t_reco->AddFriend(t_comp);
    //             N_reco = t_reco->GetEntries(MC_component+truthMatch);
    //         }
    //         else
    //         {
    //             N_reco = t_reco->GetEntries(truthMatch);
    //         }

    //         if(species == 11)
    //         {
    //             N_gen = t_gen->GetEntries();
    //             N_gen += t_gen1->GetEntries();
    //         }
    //         else if(species == 1)
    //         {
    //             N_gen = t_gen->GetEntries();
    //             N_gen += t_gen1->GetEntries();
    //             N_gen += t_gen2->GetEntries();
    //             N_gen += t_gen3->GetEntries();
    //         }
    //         else
    //         {
    //             N_gen = t_gen->GetEntries();
    //         }

    //         eps_reco = N_reco / N_gen;
    //     }

    //     // Trigger
    //     cout << "Trigger efficiency" << endl;
    //     Double_t N_L0 = t_reco->GetEntries(MC_component+mass+truthMatch+L0_trigger);
    //     Double_t N_L0_HLT1 = t_reco->GetEntries(MC_component+mass+truthMatch+L0_trigger+HLT1_trigger);
    //     Double_t N_trigger = t_reco->GetEntries(MC_component+mass+truthMatch+trigger);
    //     Double_t eps_L0 = N_L0 / N_reco;
    //     Double_t eps_L0_HLT1 = N_L0_HLT1 / N_reco; 
    //     Double_t eps_trigger = N_trigger / N_reco;

    //     TTree* t_intermediate;
    //     TChain* t_fit;
    //     if(addFit)
    //     {
    //         t_intermediate = (TTree*)t_reco->CopyTree(MC_component+truthMatch+trigger);

    //         TFileCollection* fc_fit = new TFileCollection("fc_fit", "fc_fit", FILES_fit);
    //         t_fit = new TChain("DecayTree");
    //         t_fit->AddFileInfoList((TCollection*)fc_fit->GetList());

    //         t_intermediate->AddFriend(t_fit);
    //     }
    
    //     // DTF
    //     cout << "DTF efficiency" << endl;
    //     Double_t N_pass;
    //     if(addFit)
    //     {
    //         N_pass = t_intermediate->GetEntries(passDTF);
    //     }
    //     else
    //     {
    //         N_pass = t_reco->GetEntries(mass+truthMatch+trigger+passDTF);
    //     }
    //     Double_t eps_dtf = N_pass / N_trigger;

    //     // Pre-selections
    //     cout << "Rectangular cuts efficiency" << endl;
    //     std::vector<Double_t> eps_others;
    //     std::vector<Double_t> N_others;
    //     for(int i = 0; i < N; i++)
    //     {
    //         if(addFit)
    //         {
    //             N_others.push_back( t_intermediate->GetEntries(passDTF+other_cuts[i]) );
    //         }
    //         else
    //         {
    //             N_others.push_back( t_reco->GetEntries(mass+truthMatch+trigger+passDTF+other_cuts[i]) );
    //         }
    //         eps_others.push_back( N_others[i] / N_pass );
    //     }

    //     Double_t N_others_all;
    //     if(addFit)
    //     {
    //         N_others_all = t_intermediate->GetEntries(passDTF+others);
    //     }
    //     else
    //     {
    //         N_others_all = t_reco->GetEntries(mass+truthMatch+trigger+passDTF+others);
    //     }
    //     Double_t eps_others_all = N_others_all / N_pass;

    //     file << " \\begin{table}[!htbp]" << std::endl;
    //     file << " \\centering " << std::endl;
    //     file << " \\begin{tabular}{|c|c|}" << std::endl;
    //     file << " \\hline" << std::endl;
    //     file << " & " << "Efficiency & " << " \\\\ " << std::endl;
    //     file << "\\hline" << std::endl;
    //     if(isMC)
    //     {
    //         file << "Reconstruction & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_reco*100,  eps_error(N_reco,N_gen) ) << "\\% \\\\" << std::endl;
    //     }
    //     file << "L0 & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_L0*100, eps_error(N_L0, N_reco) ) << "\\% \\\\" << std::endl;
    //     file << "L0 + HLT1 & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_L0_HLT1*100, eps_error(N_L0_HLT1, N_reco) ) << "\\% \\\\" << std::endl;
    //     file << "Trigger & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_trigger*100, eps_error(N_trigger, N_reco) ) << "\\% \\\\" << std::endl;
    //     file << "Pass DTF & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_dtf*100, eps_error(N_pass, N_trigger) ) << "\\% \\\\" << std::endl;
    //     for(int i = 0; i < N; i++)
    //     {
    //         file << other_cuts[i]+" & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_others[i]*100, eps_error(N_others[i], N_pass) ) << "\\% \\\\" << std::endl;
    //     }
    //     file << "All & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_others_all*100, eps_error(N_others_all, N_pass) ) << "\\% \\\\" << std::endl;
    //     file << "\\hline" << std::endl;
    //     file << "\\end{tabular}" << std::endl;
    //     file << "\\end{table}" << std::endl;
    //     file.close();
    // }
}

Double_t eps_error(Double_t Num, Double_t Den)
{
    return (Num/Den)*sqrt( 1./Num + 1./Den )*100;
}

TCut truthMatch_cocktailMC(TString name)
{
    map<TString, TCut> theMap;

    TCut final_state = "(abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(Kp_TRUEID ) == 321)";

    // BuDDKp
    // B+ -> D0D0K+
    TCut BuD0D0Kp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("D0_D0",521); 
    TCut BuD0starD0Kp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("D0star_D0",521);
    TCut BuD0D0starKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("D0_D0star",521);
    TCut BuD0starD0starKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("D0star_D0star",521);
    theMap["BuD0D0Kp"] = BuD0D0Kp;
    theMap["BuD0starD0Kp"] = BuD0starD0Kp;
    theMap["BuD0D0starKp"] = BuD0D0starKp;
    theMap["BuD0starD0starKp"] = BuD0starD0starKp;

    // B+ -> D+D-K+
    TCut BuDpDmKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Dp_Dm",521);
    TCut BuDpstarDmKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Dpstar_Dm",521);
    TCut BuDpDmstarKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Dp_Dmstar",521);
    TCut BuDpstarDmstarKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Dpstar_Dmstar",521);
    theMap["BuDpDmKp"] = BuDpDmKp;
    theMap["BuDpstarDmKp"] = BuDpstarDmKp;
    theMap["BuDpDmstarKp"] = BuDpDmstarKp;
    theMap["BuDpstarDmstarKp"] = BuDpstarDmstarKp;

    // B+ -> Ds+Ds-K+
    TCut BuDsDsKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Ds_Ds",521);
    TCut BuDsstarDsKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Dsstar_Ds",521);
    TCut BuDsDsstarKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Ds_Dsstar",521);
    TCut BuDsstarDsstarKp = "(abs(Bp_TRUEID) == 521)"+final_state+DD_cut("Dsstar_Dsstar",521);
    theMap["BuDsDsKp"] = BuDsDsKp;
    theMap["BuDsstarDsKp"] = BuDsstarDsKp;
    theMap["BuDsDsstarKp"] = BuDsDsstarKp;
    theMap["BuDsstarDsstarKp"] = BuDsstarDsstarKp;

    // BdDDKp
    // B0 -> D-D0K+
    TCut BdDmD0Kp = "(abs(Bp_TRUEID) == 511)"+final_state+DD_cut("Dp_D0",511);
    TCut BdDmstarD0Kp = "(abs(Bp_TRUEID) == 511)"+final_state+DD_cut("Dpstar_D0",511);
    TCut BdDmD0starKp = "(abs(Bp_TRUEID) == 511)"+final_state+DD_cut("Dp_D0star",511);
    TCut BdDmstarD0starKp = "(abs(Bp_TRUEID) == 511)"+final_state+DD_cut("Dpstar_D0star",511);
    theMap["BdDmD0Kp"] = BdDmD0Kp;
    theMap["BdDmstarD0Kp"] = BdDmstarD0Kp;
    theMap["BdDmD0starKp"] = BdDmD0starKp;
    theMap["BdDmstarD0starKp"] = BdDmstarD0starKp;

    // BsDDKp
    // Bs -> Ds-D0K+
    TCut BsDsD0Kp = "(abs(Bp_TRUEID) == 531)"+final_state+DD_cut("D0_Ds",531);
    TCut BsDsstarD0Kp = "(abs(Bp_TRUEID) == 531)"+final_state+DD_cut("D0_Dsstar",531);
    TCut BsDsD0starKp = "(abs(Bp_TRUEID) == 531)"+final_state+DD_cut("D0star_Ds",531);
    TCut BsDsstarD0starKp = "(abs(Bp_TRUEID) == 531)"+final_state+DD_cut("D0star_Dsstar",531);
    theMap["BsDsD0Kp"] = BsDsD0Kp;
    theMap["BsDsstarD0Kp"] = BsDsstarD0Kp;
    theMap["BsDsD0starKp"] = BsDsD0starKp;
    theMap["BsDsstarD0starKp"] = BsDsstarD0starKp;

    // BuDDK0
    // B+ -> D0D+K0
    TCut BuD0DpK0 = final_state+DD_cut("Dp_D0",521);
    TCut BuD0starDpK0 = final_state+DD_cut("Dp_D0star",521);
    TCut BuD0DpstarK0 = final_state+DD_cut("Dpstar_D0",521);
    TCut BuD0starDpstarK0 = final_state+DD_cut("Dpstar_D0star",521);
    theMap["BuD0DpK0"] = BuD0DpK0;
    theMap["BuD0starDpK0"] = BuD0starDpK0;
    theMap["BuD0DpstarK0"] = BuD0DpstarK0;
    theMap["BuD0starDpstarK0"] = BuD0starDpstarK0;

    // BdDDK0
    // B0 -> D+D-K0
    TCut BdDpDmK0 = final_state+DD_cut("Dp_Dm",511);
    TCut BdDpstarDmK0 = final_state+DD_cut("Dpstar_Dm",511);
    TCut BdDpDmstarK0 = final_state+DD_cut("Dp_Dmstar",511);
    TCut BdDpstarDmstarK0 = final_state+DD_cut("Dpstar_Dmstar",511);
    theMap["BdDpDmK0"] = BdDpDmK0;
    theMap["BdDpstarDmK0"] = BdDpstarDmK0;
    theMap["BdDpDmstarK0"] = BdDpDmstarK0;
    theMap["BdDpstarDmstarK0"] = BdDpstarDmstarK0;

    // B0 -> D0D0K0
    TCut BdD0D0K0 = final_state+DD_cut("D0_D0",511);
    TCut BdD0starD0K0 = final_state+DD_cut("D0star_D0",511);
    TCut BdD0D0starK0 = final_state+DD_cut("D0_D0star",511);
    TCut BdD0starD0starK0 = final_state+DD_cut("D0star_D0star",511);
    theMap["BdD0D0K0"] = BdD0D0K0;
    theMap["BdD0starD0K0"] = BdD0starD0K0;
    theMap["BdD0D0starK0"] = BdD0D0starK0;
    theMap["BdD0starD0starK0"] = BdD0starD0starK0;

    // BuDD
    // B+ -> D0Ds+
    TCut BuD0Ds = final_state+DD_cut("D0_Ds",521);
    TCut BuD0starDs = final_state+DD_cut("D0star_Ds",521);
    TCut BuD0Dsstar = final_state+DD_cut("D0_Dsstar",521);
    TCut BuD0starDsstar = final_state+DD_cut("D0star_Dsstar",521);
    theMap["BuD0Ds"] = BuD0Ds;
    theMap["BuD0starDs"] = BuD0starDs;
    theMap["BuD0Dsstar"] = BuD0Dsstar;
    theMap["BuD0starDsstar"] = BuD0starDsstar;

    // B+ -> D0D+
    TCut BuD0Dp = final_state+DD_cut("Dp_D0",521);
    TCut BuD0starDp = final_state+DD_cut("Dp_D0star",521);
    TCut BuD0Dpstar = final_state+DD_cut("Dpstar_D0",521);
    TCut BuD0starDpstar = final_state+DD_cut("Dpstar_D0star",521);
    theMap["BuD0Dp"] = BuD0Dp;
    theMap["BuD0starDp"] = BuD0starDp;
    theMap["BuD0Dpstar"] = BuD0Dpstar;
    theMap["BuD0starDpstar"] = BuD0starDpstar;

    // BdDD
    // B0 -> D0D0
    TCut BdD0D0 = final_state+DD_cut("D0_D0",511);
    TCut BdD0starD0 = final_state+DD_cut("D0star_D0",511);
    TCut BdD0D0star = final_state+DD_cut("D0_D0star",511);
    TCut BdD0starD0star = final_state+DD_cut("D0star_D0star",511);
    theMap["BdD0D0"] = BdD0D0;
    theMap["BdD0starD0"] = BdD0starD0;
    theMap["BdD0D0star"] = BdD0D0star;
    theMap["BdD0starD0star"] = BdD0starD0star;

    // B0 -> D+D-
    TCut BdDpDm = final_state+DD_cut("Dp_Dm",511);
    TCut BdDpstarDm = final_state+DD_cut("Dpstar_Dm",511);
    TCut BdDpDmstar = final_state+DD_cut("Dp_Dmstar",511);
    TCut BdDpstarDmstar = final_state+DD_cut("Dpstar_Dmstar",511);
    theMap["BdDpDm"] = BdDpDm;
    theMap["BdDpstarDm"] = BdDpstarDm;
    theMap["BdDpDmstar"] = BdDpDmstar;
    theMap["BdDpstarDmstar"] = BdDpstarDmstar;

    // B0 -> D-Ds+
    TCut BdDpDs = final_state+DD_cut("Dp_Ds",511);
    TCut BdDpstarDs = final_state+DD_cut("Dpstar_Ds",511);
    TCut BdDpDsstar = final_state+DD_cut("Dp_Dsstar",511);
    TCut BdDpstarDsstar = final_state+DD_cut("Dpstar_Dsstar",511);
    theMap["BdDpDs"] = BdDpDs;
    theMap["BdDpstarDs"] = BdDpstarDs;
    theMap["BdDpDsstar"] = BdDpDsstar;
    theMap["BdDpstarDsstar"] = BdDpstarDsstar;

    // B0 -> Ds+Ds-
    TCut BdDsDs = final_state+DD_cut("Ds_Ds",511);
    TCut BdDsstarDs = final_state+DD_cut("Dsstar_Ds",511);
    TCut BdDsDsstar = final_state+DD_cut("Ds_Dsstar",511);
    TCut BdDsstarDsstar = final_state+DD_cut("Dsstar_Dsstar",511);
    theMap["BdDsDs"] = BdDsDs;
    theMap["BdDsstarDs"] = BdDsstarDs;
    theMap["BdDsDsstar"] = BdDsDsstar;
    theMap["BdDsstarDsstar"] = BdDsstarDsstar;

    // Bs -> DD
    // Bs -> Ds+Ds-
    TCut BsDsDs = final_state+DD_cut("Ds_Ds",531);
    TCut BsDsstarDs = final_state+DD_cut("Dsstar_Ds",531);
    TCut BsDsDsstar = final_state+DD_cut("Ds_Dsstar",531);
    TCut BsDsstarDsstar = final_state+DD_cut("Dsstar_Dsstar",531);
    theMap["BsDsDs"] = BsDsDs;
    theMap["BsDsstarDs"] = BsDsstarDs;
    theMap["BsDsDsstar"] = BsDsDsstar;
    theMap["BsDsstarDsstar"] = BsDsstarDsstar;

    // Bs -> D-Ds+
    TCut BsDpDs = final_state+DD_cut("Dp_Ds",531);
    TCut BsDpstarDs = final_state+DD_cut("Dpstar_Ds",531);
    TCut BsDpDsstar = final_state+DD_cut("Dp_Dsstar",531);
    TCut BsDpstarDsstar = final_state+DD_cut("Dsstar_Dsstar",531);
    theMap["BsDpDs"] = BsDpDs;
    theMap["BsDpstarDs"] = BsDpstarDs;
    theMap["BsDpDsstar"] = BsDpDsstar;
    theMap["BsDpstarDsstar"] = BsDpstarDsstar;

    // Bs -> D+D-
    TCut BsDpDm = final_state+DD_cut("Dp_Dm",531);
    TCut BsDpstarDm = final_state+DD_cut("Dpstar_Dm",531);
    TCut BsDpDmstar = final_state+DD_cut("Dp_Dmstar",531);
    TCut BsDpstarDmstar = final_state+DD_cut("Dpstar_Dmstar",531);
    theMap["BsDpDm"] = BsDpDm;
    theMap["BsDpstarDm"] = BsDpstarDm;
    theMap["BsDpDmstar"] = BsDpDmstar;
    theMap["BsDpstarDmstar"] = BsDpstarDmstar;

    // Bs -> D0D0
    TCut BsD0D0 = final_state+DD_cut("D0_D0",531);
    TCut BsD0starD0 = final_state+DD_cut("D0star_D0",531);
    TCut BsD0D0star = final_state+DD_cut("D0_D0star",531);
    TCut BsD0starD0star = final_state+DD_cut("D0star_D0star",531);
    theMap["BsD0D0"] = BsD0D0;
    theMap["BsD0starD0"] = BsD0starD0;
    theMap["BsD0D0star"] = BsD0D0star;
    theMap["BsD0starD0star"] = BsD0starD0star;

    return theMap[name];
}

TCut DD_cut(TString name, Int_t mother_ID)
{
    // Ds+ Ds-
    TCut Ds_Ds = Form("(abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i)",mother_ID,mother_ID);
    TCut Dsstar_Ds = Form("((taup_TRUEID == 431) && (taup_MC_MOTHER_ID == 433) &&  (taup_MC_GD_MOTHER_ID == %i) && (taum_TRUEID == -431) && (taum_MC_MOTHER_ID == %i)) || ((taum_TRUEID == 431) && (taum_MC_MOTHER_ID == 433) &&  (taum_MC_GD_MOTHER_ID == %i) && (taup_TRUEID == -431) && (taup_MC_MOTHER_ID == %i)) || ((taup_TRUEID == -431) && (taup_MC_MOTHER_ID == -433) &&  (taup_MC_GD_MOTHER_ID == -%i) && (taum_TRUEID == 431) && (taum_MC_MOTHER_ID == -%i)) || ((taum_TRUEID == -431) && (taum_MC_MOTHER_ID == -433) &&  (taum_MC_GD_MOTHER_ID == -%i) && (taup_TRUEID == 431) && (taup_MC_MOTHER_ID == -%i))", mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID);
    TCut Ds_Dsstar = Form("((taup_TRUEID == 431) && (taup_MC_MOTHER_ID == %i) && (taum_TRUEID == -431) && (taum_MC_MOTHER_ID == -433) && (taum_MC_GD_MOTHER_ID == %i)) || ((taum_TRUEID == 431) && (taum_MC_MOTHER_ID == %i) && (taup_TRUEID == -431) && (taup_MC_MOTHER_ID == -433) && (taup_MC_GD_MOTHER_ID == %i)) || ((taup_TRUEID == -431) && (taup_MC_MOTHER_ID == -%i) && (taum_TRUEID == 431) && (taum_MC_MOTHER_ID == 433) && (taum_MC_GD_MOTHER_ID == -%i)) || ((taum_TRUEID == -431) && (taum_MC_MOTHER_ID == -%i) && (taup_TRUEID == 431) && (taup_MC_MOTHER_ID == 433) && (taup_MC_GD_MOTHER_ID == -%i))", mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID);
    TCut Dsstar_Dsstar = Form("(abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)==%i)",mother_ID,mother_ID);

    // D0bar D0
    TCut D0_D0 = Form("(abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i)",mother_ID,mother_ID);
    TCut D0star_D0 = Form("((taup_TRUEID == -421) && (taup_MC_MOTHER_ID == -423) && (taup_MC_GD_MOTHER_ID == %i) && (taum_TRUEID == 421) && (taum_MC_MOTHER_ID == %i)) || ((taum_TRUEID == -421) && (taum_MC_MOTHER_ID == -423) && (taum_MC_GD_MOTHER_ID == %i) && (taup_TRUEID == 421) && (taup_MC_MOTHER_ID == %i)) || ((taup_TRUEID == 421) && (taup_MC_MOTHER_ID == 423) && (taup_MC_GD_MOTHER_ID == -%i) && (taum_TRUEID == -421) && (taum_MC_MOTHER_ID == -%i)) || ((taum_TRUEID == 421) && (taum_MC_MOTHER_ID == 423) && (taum_MC_GD_MOTHER_ID == -%i) && (taup_TRUEID == -421) && (taup_MC_MOTHER_ID == -%i))", mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID);
    TCut D0_D0star = Form("((taup_TRUEID == -421) && (taup_MC_MOTHER_ID == %i) && (taum_TRUEID == 421) && (taum_MC_MOTHER_ID == 423) && (taum_MC_GD_MOTHER_ID == %i)) || ((taum_TRUEID == -421) && (taum_MC_MOTHER_ID == %i) && (taup_TRUEID == 421) && (taup_MC_MOTHER_ID == 423) && (taup_MC_GD_MOTHER_ID == %i)) || ((taup_TRUEID == 421) && (taup_MC_MOTHER_ID == -%i) && (taum_TRUEID == -421) && (taum_MC_MOTHER_ID == -423) && (taum_MC_GD_MOTHER_ID == -%i)) || ((taum_TRUEID == 421) && (taum_MC_MOTHER_ID == -%i) && (taup_TRUEID == -421) && (taup_MC_MOTHER_ID == -423) && (taup_MC_GD_MOTHER_ID == -%i))", mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID, mother_ID);
    TCut D0star_D0star = Form("(abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)==%i)",mother_ID,mother_ID);

    // D+ D-
    TCut Dp_Dm = Form("(abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i)",mother_ID,mother_ID);
    TCut Dpstar_Dm = Form("(((taup_TRUEID == 411) && (taum_TRUEID == -411) && (taup_MC_MOTHER_ID == 413) && (taup_MC_GD_MOTHER_ID == %i) && (taum_MC_MOTHER_ID == %i)) || ((taup_TRUEID == 413) && (taum_TRUEID == -411) && (taup_MC_MOTHER_ID == %i) && (taum_MC_MOTHER_ID == %i))) || (((taum_TRUEID == 411) && (taup_TRUEID == -411) && (taum_MC_MOTHER_ID == 413) && (taum_MC_GD_MOTHER_ID == %i) && (taup_MC_MOTHER_ID == %i)) || ((taum_TRUEID == 413) && (taup_TRUEID == -411) && (taum_MC_MOTHER_ID == %i) && (taup_MC_MOTHER_ID == %i))) || (((taup_TRUEID == -411) && (taum_TRUEID == 411) && (taup_MC_MOTHER_ID == -413) && (taup_MC_GD_MOTHER_ID == -%i) && (taum_MC_MOTHER_ID == -%i)) || ((taup_TRUEID == -413) && (taum_TRUEID == 411) && (taup_MC_MOTHER_ID == -%i) && (taum_MC_MOTHER_ID == -%i))) || (((taum_TRUEID == -411) && (taup_TRUEID == 411) && (taum_MC_MOTHER_ID == -413) && (taum_MC_GD_MOTHER_ID == -%i) && (taup_MC_MOTHER_ID == -%i)) || ((taum_TRUEID == -413) && (taup_TRUEID == 411) && (taum_MC_MOTHER_ID == -%i) && (taup_MC_MOTHER_ID == -%i)))",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dp_Dmstar = Form("(((taup_TRUEID == 411) && (taum_TRUEID == -411) && (taup_MC_MOTHER_ID == %i) && (taum_MC_MOTHER_ID == -413) && (taum_MC_GD_MOTHER_ID == %i)) || ((taup_TRUEID == 411) && (taum_TRUEID == -413) && (taup_MC_MOTHER_ID == %i) && (taum_MC_MOTHER_ID == %i))) || (((taum_TRUEID == 411) && (taup_TRUEID == -411) && (taum_MC_MOTHER_ID == %i) && (taup_MC_MOTHER_ID == -413) && (taup_MC_GD_MOTHER_ID == %i)) || ((taum_TRUEID == 411) && (taup_TRUEID == -413) && (taum_MC_MOTHER_ID == %i) && (taup_MC_MOTHER_ID == %i))) || (((taup_TRUEID == -411) && (taum_TRUEID == 411) && (taup_MC_MOTHER_ID == -%i) && (taum_MC_MOTHER_ID == 413) && (taum_MC_GD_MOTHER_ID == -%i)) || ((taup_TRUEID == -411) && (taum_TRUEID == 413) && (taup_MC_MOTHER_ID == -%i) && (taum_MC_MOTHER_ID == -%i))) || (((taum_TRUEID == -411) && (taup_TRUEID == 411) && (taum_MC_MOTHER_ID == -%i) && (taup_MC_MOTHER_ID == 413) && (taup_MC_GD_MOTHER_ID == -%i)) || ((taum_TRUEID == -411) && (taup_TRUEID == 413) && (taum_MC_MOTHER_ID == -%i) && (taup_MC_MOTHER_ID == -%i)))",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dpstar_Dmstar = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);

    // D0 Ds
    TCut D0_Ds = Form("( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut D0star_Ds = Form("( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut D0_Dsstar = Form("( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut D0star_Dsstar = Form("( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);

    // D+ D0
    TCut Dp_D0 = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dpstar_D0 = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dp_D0star = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dpstar_D0star = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);

    // D+ Ds
    TCut Dp_Ds = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dpstar_Ds = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==413) &&  (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dp_Dsstar = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID);
    TCut Dpstar_Dsstar = Form("( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==%i) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)==%i) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)==%i) &&  (abs(taum_MC_MOTHER_ID)==%i) )",mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID,mother_ID);

    if(name == "Ds_Ds") // Ds+Ds-
    {
        return Ds_Ds;
    }
    else if(name == "Dsstar_Ds")
    {
        return Dsstar_Ds;
    }
    else if(name == "Ds_Dsstar")
    {
        return Ds_Dsstar;
    }
    else if(name == "Dsstar_Dsstar")
    {
        return Dsstar_Dsstar;
    }
    else if(name == "D0_D0") // D0barD0
    {
        return D0_D0;
    }
    else if(name == "D0star_D0")
    {
        return D0star_D0;
    }
    else if(name == "D0_D0star")
    {
        return D0_D0star;
    }
    else if(name == "D0star_D0star")
    {
        return D0star_D0star;
    }
    else if(name == "Dp_Dm") // D+D-
    {
        return Dp_Dm;
    }
    else if(name == "Dpstar_Dm")
    {
        return Dpstar_Dm;
    }
    else if(name == "Dp_Dmstar")
    {
        return Dp_Dmstar;
    }
    else if(name == "Dp_Dmstar")
    {
        return Dp_Dmstar;
    }
    else if(name == "Dpstar_Dmstar")
    {
        return Dpstar_Dmstar;
    }
    else if(name == "D0_Ds") // D0barDs+
    {
        return D0_Ds;
    }
    else if(name == "D0star_Ds")
    {
        return D0star_Ds;
    }
    else if(name == "D0_Dsstar")
    {
        return D0_Dsstar;
    }
    else if(name == "D0star_Dsstar")
    {
        return D0star_Dsstar;
    }
    else if(name == "Dp_D0") // D+D+
    {
        return Dp_D0;
    }
    else if(name == "Dpstar_D0")
    {
        return Dpstar_D0;
    }
    else if(name == "Dp_D0star")
    {
        return Dp_D0star;
    }
    else if(name == "Dpstar_D0star")
    {
        return Dpstar_D0star;
    }
    else if(name == "Dp_Ds") // D+Ds-
    {
        return Dp_Ds;
    }
    else if(name == "Dpstar_Ds")
    {
        return Dpstar_Ds;
    }
    else if(name == "Dp_Dsstar")
    {
        return Dp_Dsstar;
    }
    else if(name == "Dpstar_Dsstar")
    {
        return Dpstar_Dsstar;
    }
}