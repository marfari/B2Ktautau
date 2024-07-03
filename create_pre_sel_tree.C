
Double_t eps_error(Double_t Num, Double_t Den);

Bool_t addFit = false;

void create_pre_sel_tree(int year, int species, int line, bool createTable)
{   
    TString path;
    if(createTable)
    {
        path = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/pre_sel_table_species_%i.tex",year,species);
    }
    else
    {
        path = "";
    }
    std::ofstream file(path);

    Bool_t isKtautauMC = false;
    if( (species == 10) || (species == 11) || (species == 12) || (species == 1) )
    {
        isKtautauMC = true;
    }

    if(isKtautauMC || (species == 2) || (species == 3))
    {
        addFit = true;
    }

    TString FILES;
    if(isKtautauMC) // Ktautau MC
    {
        FILES = Form("Files_on_grid/MC_201%i.txt",year);
    }
    else if(species == 4) // DDK MC
    {
        FILES = Form("Files_on_grid/MC_DDK_201%i.txt",year);
    }
    else if((species == 2) || (species == 3)) // Ktautau data
    {
        FILES = Form("Files_on_grid/data_201%i.txt",year);
    }
    else if((species == 5) || (species == 6)) // DDK data
    {
        FILES = Form("Files_on_grid/data_DDK_201%i.txt",year);
    }
    else if( (species == 7) || (species == 71))
    {
        FILES = Form("Files_on_grid/MC_D0Dps_201%i.txt",year);
    }
    else if((species == 8) || (species == 81) || (species == 82) || (species == 83))
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

    TString MC_component_file;
    if(isKtautauMC)
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/%i.root",year,line);
    }

    TString FILES_fit;
    if(addFit)
    {
        FILES_fit = Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/fit_results.txt",year,species);
    }

    // Pre-selection cuts
    // MC truth-match
    TCut truthMatch;
    if(isKtautauMC){
        truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
    } 
    else if(species == 4) // D+D-K+ MC
    {
        truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(Dp_K_TRUEID) == 321) && (abs(Dp_pi1_TRUEID) == 211) && (abs(Dp_pi2_TRUEID) == 211) && (abs(Dm_K_TRUEID) == 321) && (abs(Dm_pi1_TRUEID) == 211) && (abs(Dm_pi2_TRUEID) == 211) && (abs(Dp_TRUEID) == 411) && (abs(Dm_TRUEID) == 411) && (abs(Bp_TRUEID) == 521)";
    }
    else if(species == 7) // D0D+s MC
    {
        truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(Dsp_K1_TRUEID) == 321) && (abs(Dsp_K2_TRUEID) == 321) && (abs(D0bar_pi_TRUEID) == 211) && (abs(Dsp_pi_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(Dsp_TRUEID) == 431) && (abs(Bp_TRUEID) == 521)";
    }
    else if(species == 9) // D0D0K MC
    {
        truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(D0bar_pi1_TRUEID) == 211) && (abs(D0bar_pi2_TRUEID) == 211) && (abs(D0bar_pi3_TRUEID) == 211) && (abs(D0_K_TRUEID) == 321) && (abs(D0_pi1_TRUEID) == 211) && (abs(D0_pi2_TRUEID) == 211) && (abs(D0_pi3_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(D0_TRUEID) == 421) && (abs(Kp_TRUEID) == 321) && (abs(Bp_TRUEID) == 521)";
    }

    // Trigger cuts:
    TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
    TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
    TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
    TCut trigger = L0_trigger+HLT1_trigger+HLT2_trigger;

    // pass fitter
    TCut passDTF;
    if(addFit)
    {
        passDTF = "df_status==0";
    }
    else
    {
        passDTF = "Bp_dtf_status[0]==0";
    }

    TCut others;
    std::vector<TCut> other_cuts;
    if(isKtautauMC || (species == 2) || (species == 3))
    {
        other_cuts.push_back("sqrt( pow(df_BVx - df_PVx,2) + pow(df_BVy - df_PVy,2) + pow(df_BVz - df_PVz,2) ) > 2");
        other_cuts.push_back("(taup_M > 800) && (taup_M < 1650)");
        other_cuts.push_back("(taum_M > 800) && (taum_M < 1650)");
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
    if((species == 7) || (species == 8)) // D0barD+s
    { 
        other_cuts.push_back("Bp_FD_OWNPV < 80");
        other_cuts.push_back("D0bar_K_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K1_ProbNNk > 0.55");
        other_cuts.push_back("Dsp_K2_ProbNNk > 0.55");
        other_cuts.push_back("(D0bar_M > 1820) && (D0bar_M < 1910)");
        other_cuts.push_back("(Dsp_M > 1930) && (Dsp_M < 2000)");
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

    TFileCollection* fc = new TFileCollection("fc", "fc", FILES, 1, line);
    TChain* t;

    if((species == 3) || (species == 6) || (species == -1)) // WS data (Ktautau,D+D-K+,D0D0K)
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
        pre_selections = MC_component+truthMatch+trigger+passDTF+others;
    }
    else if((species == 4) || (species == 7) || (species == 9)) // D+D-K+, D0barD+s, D0D0K MC
    {
        pre_selections = truthMatch+trigger+passDTF+others;
    }
    else if(species == 71) // D0barD+s un-TM MC
    {
        pre_selections = trigger+passDTF+others;
    }
    else if((species == 2) || (species == 3)) // Ktautau data
    {
        pre_selections = trigger+passDTF+others;
    }
    else if((species == 5) || (species == 6) || (species == 8) || (species == 81) || (species == 82) || (species == 83) || (species == 0) || (species == -1)) // D+D-K+ and D0barD+s data
    {
        pre_selections = trigger+passDTF+others;
    }


    TString fout_name = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/%i.root",year,species,line);
    if(addFit)
    {
        TTree* t_intermediate;
        if(isKtautauMC)
        {
            t_intermediate = (TTree*)t->CopyTree(MC_component+truthMatch+trigger);
        }
        else
        {
            t_intermediate = (TTree*)t->CopyTree(trigger);   
        }

        // TFileCollection* fc_fit = new TFileCollection("fc_fit", "fc_fit", FILES_fit, 1, line);
        // t_fit = new TChain("DecayTree");
        // t_fit->AddFileInfoList((TCollection*)fc_fit->GetList());
        TFile* f_fit = new TFile(Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/%i.root",year,species,line));
        TTree* t_fit = (TTree*)f_fit->Get("DecayTree");

        // cout << t_intermediate->GetEntries() << endl;
        // cout << t_fit->GetEntries() << endl;

        t_intermediate->AddFriend(t_fit, "fit");
        // t_fit->AddFriend(t_intermediate);

        // TTree* t_meas_pre_sel = (TTree*)t_intermediate->CopyTree(pre_selections);
        // TTree* t_fit_pre_sel = (TTree*)t_fit->CopyTree(pre_selections);
        // t_meas_pre_sel->AddFriend(t_fit_pre_sel);

        ROOT::RDataFrame df(*t_intermediate);
        df.Filter(pre_selections.GetTitle()).Snapshot("DecayTree", fout_name);
    }
    else
    {
        TFile* fout = new TFile(fout_name, "RECREATE");   
        fout->cd();
        TTree* t_pre_sel = (TTree*)t->CopyTree(pre_selections);
        t_pre_sel->Write();
        fout->Close();
    }

    // I want a TTree with t_intermediate and t_fit branches

    // TChain* t_fit_pre_sel;
    // TChain* t_pre_sel_measured;
    // if(addFit)
    // {   
    //     t_fit_pre_sel = (TChain*)t_fit->CopyTree(pre_selections);
    //     t_fit_pre_sel->SetName("DecayTree_fit");
    //     t_fit_pre_sel->SetTitle("DecayTree_fit");

    //     t_pre_sel_measured = (TChain*)t_intermediate->CopyTree(pre_selections);
    //     t_pre_sel_measured->SetName("DecayTree_measured");
    //     t_pre_sel_measured->SetTitle("DecayTree_measured");
    // }

    // TTree* t_pre_sel;
    // if(addFit)
    // {
    //     // t_fit_pre_sel->AddFriend(t_pre_sel_measured);
    //     ROOT::RDataFrame df(*t_fit_pre_sel);
    //     df.Snapshot("DecayTree", fout_name);

    //     // t_pre_sel_measured->Add(t_fit_pre_sel);
    //     // t_pre_sel_measured->Write();

    //     // t_pre_sel_measured->AddFriend(t_fit_pre_sel);
    //     // t_pre_sel_measured->Write();
    // }


    cout << "Finished successfully" << endl;

    if(createTable) // creates the pre-selection efficiency tables
    {
        cout << "Creating pre-selection efficiency tables" << endl;
        Bool_t isMC = false;
        if((species == 1) || (species == 10) || (species == 11) || (species == 12) || (species == 4) || (species == 7) || (species == 9))
        {
            isMC = true;
        }

        TCut mass = ""; // background definition
        // if(species == 8) 
        // {
        //     mass = "Bp_dtf_M[0] > 5350";

        // }

        TFileCollection* fc_all;
        if(isMC)
        {
            fc_all = new TFileCollection("fc_all", "fc_all", FILES);
        }
        else
        {
            fc_all = new TFileCollection("fc_all", "fc_all", FILES);
        }
         
        TChain* t_reco;
        if((species == 3) || (species == 6) || (species == -1)) // WS data (Ktautau,D+D-K+,D0D0K)
        {
            t_reco = new TChain("ntuple_SS/DecayTree");
        }
        else // MC and RS data
        {
            t_reco = new TChain("ntuple/DecayTree");
        }
        t_reco->AddFileInfoList((TCollection*)fc_all->GetList());

        // Reconstruction efficiency
        cout << "Reconstruction efficiency" << endl;
        TChain* t_gen;
        TChain* t_gen1;
        TChain* t_gen2;
        TChain* t_gen3;
        Double_t eps_reco = 1.;
        Double_t N_reco = t_reco->GetEntries(mass);
        Double_t N_gen;
        if(isMC)
        {
            if(species == 10)
            {
                t_gen = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
            }
            else if(species == 11)
            {
                t_gen = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
                t_gen1 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
            }
            else if(species == 12)
            {
                t_gen = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
            }
            else if(species == 1)
            {
                t_gen1 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
                t_gen2 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
                t_gen3 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
                t_gen = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
            }
            else
            {
                t_gen = new TChain("mc_ntuple/MCDecayTree");
            }
            t_gen->AddFileInfoList((TCollection*)fc_all->GetList());

            if(species == 11)
            {
                t_gen1->AddFileInfoList((TCollection*)fc_all->GetList());
            }
            if(species == 1)
            {
                t_gen1->AddFileInfoList((TCollection*)fc_all->GetList());
                t_gen2->AddFileInfoList((TCollection*)fc_all->GetList());
                t_gen3->AddFileInfoList((TCollection*)fc_all->GetList());
            }

            if(isKtautauMC)
            {
                TFileCollection* fc_comp = new TFileCollection("fc_comp", "fc_comp", Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/mc_components.txt",year));
                TChain* t_comp = new TChain("DecayTree");
                t_comp->AddFileInfoList((TCollection*)fc_comp->GetList());

                t_reco->AddFriend(t_comp);
                N_reco = t_reco->GetEntries(MC_component+truthMatch);
            }
            else
            {
                N_reco = t_reco->GetEntries(truthMatch);
            }

            if(species == 11)
            {
                N_gen = t_gen->GetEntries();
                N_gen += t_gen1->GetEntries();
            }
            else if(species == 1)
            {
                N_gen = t_gen->GetEntries();
                N_gen += t_gen1->GetEntries();
                N_gen += t_gen2->GetEntries();
                N_gen += t_gen3->GetEntries();
            }
            else
            {
                N_gen = t_gen->GetEntries();
            }

            eps_reco = N_reco / N_gen;
        }

        // Trigger
        cout << "Trigger efficiency" << endl;
        Double_t N_L0 = t_reco->GetEntries(MC_component+mass+truthMatch+L0_trigger);
        Double_t N_L0_HLT1 = t_reco->GetEntries(MC_component+mass+truthMatch+L0_trigger+HLT1_trigger);
        Double_t N_trigger = t_reco->GetEntries(MC_component+mass+truthMatch+trigger);
        Double_t eps_L0 = N_L0 / N_reco;
        Double_t eps_L0_HLT1 = N_L0_HLT1 / N_reco; 
        Double_t eps_trigger = N_trigger / N_reco;

        TTree* t_intermediate;
        TChain* t_fit;
        if(addFit)
        {
            t_intermediate = (TTree*)t_reco->CopyTree(MC_component+truthMatch+trigger);

            TFileCollection* fc_fit = new TFileCollection("fc_fit", "fc_fit", FILES_fit);
            t_fit = new TChain("DecayTree");
            t_fit->AddFileInfoList((TCollection*)fc_fit->GetList());

            t_intermediate->AddFriend(t_fit);
        }
    
        // DTF
        cout << "DTF efficiency" << endl;
        Double_t N_pass;
        if(addFit)
        {
            N_pass = t_intermediate->GetEntries(passDTF);
        }
        else
        {
            N_pass = t_reco->GetEntries(mass+truthMatch+trigger+passDTF);
        }
        Double_t eps_dtf = N_pass / N_trigger;

        // Pre-selections
        cout << "Rectangular cuts efficiency" << endl;
        std::vector<Double_t> eps_others;
        std::vector<Double_t> N_others;
        for(int i = 0; i < N; i++)
        {
            if(addFit)
            {
                N_others.push_back( t_intermediate->GetEntries(passDTF+other_cuts[i]) );
            }
            else
            {
                N_others.push_back( t_reco->GetEntries(mass+truthMatch+trigger+passDTF+other_cuts[i]) );
            }
            eps_others.push_back( N_others[i] / N_pass );
        }

        Double_t N_others_all;
        if(addFit)
        {
            N_others_all = t_intermediate->GetEntries(passDTF+others);
        }
        else
        {
            N_others_all = t_reco->GetEntries(mass+truthMatch+trigger+passDTF+others);
        }
        Double_t eps_others_all = N_others_all / N_pass;

        file << " \\begin{table}[!htbp]" << std::endl;
        file << " \\centering " << std::endl;
        file << " \\begin{tabular}{|c|c|}" << std::endl;
        file << " \\hline" << std::endl;
        file << " & " << "Efficiency & " << " \\\\ " << std::endl;
        file << "\\hline" << std::endl;
        if(isMC)
        {
            file << "Reconstruction & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_reco*100,  eps_error(N_reco,N_gen) ) << "\\% \\\\" << std::endl;
        }
        file << "L0 & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_L0*100, eps_error(N_L0, N_reco) ) << "\\% \\\\" << std::endl;
        file << "L0 + HLT1 & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_L0_HLT1*100, eps_error(N_L0_HLT1, N_reco) ) << "\\% \\\\" << std::endl;
        file << "Trigger & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_trigger*100, eps_error(N_trigger, N_reco) ) << "\\% \\\\" << std::endl;
        file << "Pass DTF & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_dtf*100, eps_error(N_pass, N_trigger) ) << "\\% \\\\" << std::endl;
        for(int i = 0; i < N; i++)
        {
            file << other_cuts[i]+" & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_others[i]*100, eps_error(N_others[i], N_pass) ) << "\\% \\\\" << std::endl;
        }
        file << "All & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_others_all*100, eps_error(N_others_all, N_pass) ) << "\\% \\\\" << std::endl;
        file << "\\hline" << std::endl;
        file << "\\end{tabular}" << std::endl;
        file << "\\end{table}" << std::endl;
        file.close();
    }
}

Double_t eps_error(Double_t Num, Double_t Den)
{
    return (Num/Den)*sqrt( 1./Num + 1./Den )*100;
}