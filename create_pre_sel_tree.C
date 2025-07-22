
Double_t eps_error(Double_t Num, Double_t Den);
TCut DD_cut(TString name, Int_t mother_ID);
TCut truthMatch_cocktailMC(TString name);

void create_pre_sel_tree(int year, int species, int line)
{   
    // TString path = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/Num_entries/%i.txt",year,species,line);
    // std::ofstream file(path);
    Bool_t is_cocktailMC = false;
    if((species == 100) || (species == 110) || (species == 120) || (species == 130) || (species == 140) || (species == 150) || (species == 160) || (species == 170)){is_cocktailMC = true;}

    Bool_t is_KtautauMC = false;
    if((species == 1) || (species == 10))
    {
        is_KtautauMC = true;
    }

    Bool_t is_KtautauData = false;
    if((species == 2) || (species == 3))
    {
        is_KtautauData = true;
    }

    TString FILES;
    if(is_KtautauMC) // Ktautau MC
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_1/pid_corr.txt",year);
    }
    else if(is_KtautauData) // Ktautau data
    {
        FILES = Form("Files_on_grid/data_201%i.txt",year);
    }
    else if((species == 100) || (species == 1000))
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_100/pid_corr.txt",year);
    }
    else if((species == 110) || (species == 1100))
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_110/pid_corr.txt",year);
    }
    else if((species == 120) || (species == 1200))
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_120/pid_corr.txt",year);
    }
    else if((species == 130) || (species == 1300))
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_130/pid_corr.txt",year);
    }
    else if((species == 150) || (species == 1500))
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_150/pid_corr.txt",year);
    }
    else if( (species == 7) || (species == 71) || (species == 72) ) // DDs MC
    {
        FILES = Form("/panfs/felician/B2Ktautau/workflow/PIDCalib/201%i/Species_7/pid_corr.txt",year);
    }
    else if((species == 8) || (species == 81)) // DDs data
    {
        FILES = Form("Files_on_grid/data_D0Dps_201%i.txt",year);
    }
    else if(species == 4) // DDK MC
    {
        FILES = Form("Files_on_grid/MC_DpDmK_201%i.txt",year);
    }
    else if((species == 5) || (species == 6)) // D+D-K data
    {
        FILES = Form("Files_on_grid/data_DDK_201%i.txt",year);
    }
    else if(species == 9) // D0D0K MC
    {
        FILES = Form("Files_on_grid/MC_D0D0K_201%i.txt",year);
    }
    else if((species == 0) || (species == -1)) // D0D0K data
    {
        FILES = Form("Files_on_grid/data_D0D0K_201%i.txt",year);
    }

    TString MC_component_file = "";
    if(is_KtautauMC)
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/mc_components.txt",year);
    }
    else if((species == 100) || (species == 1000))
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201%i/Species_100/cocktail_mc_components.txt",year);
    }
    else if((species == 110) || (species == 1100))
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201%i/Species_110/cocktail_mc_components.txt",year);
    }
    else if((species == 120) || (species == 1200))
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201%i/Species_120/cocktail_mc_components.txt",year);
    }
    else if((species == 130) || (species == 1300))
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201%i/Species_130/cocktail_mc_components.txt",year);
    }
    else if((species == 150) || (species == 1500))
    {
        MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201%i/Species_150/cocktail_mc_components.txt",year);
    }

    // Pre-selection TTree:
    // Kttautau + cocktail MCs: (truth-match) + trigger + rectangular cuts
    // Normalisation mode: (truth-match) + trigger + rectangular cuts + pass DTF + fit region

    TFile* t_selections = new TFile(Form("/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_%i/selections.root",species));
    TObjString* truthMatch_exp = (TObjString*)t_selections->Get("truthMatch");
    TObjString* trigger_exp = (TObjString*)t_selections->Get("trigger");
    TObjString* rectangular_cuts_exp = (TObjString*)t_selections->Get("rectangular_cuts");
    TObjString* pass_mass_fit_exp = (TObjString*)t_selections->Get("pass_mass_fit");
    TObjString* fit_range_exp = (TObjString*)t_selections->Get("fit_range");

    TCut truthMatch(truthMatch_exp->GetString().Data());
    TCut trigger(trigger_exp->GetString().Data());
    TCut rectangular_cuts(rectangular_cuts_exp->GetString().Data());
    TCut pass_mass_fit(pass_mass_fit_exp->GetString().Data());
    TCut fit_range(fit_range_exp->GetString().Data());

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCut pre_selections;
    if(species == 1) // Ktautau MC
    {
        pre_selections = truthMatch+trigger+rectangular_cuts;
    }
    else if((species == 2) || (species == 3) || (species == 10) || is_cocktailMC) // Ktautau data,  Ktautau MC w/o TM, cocktail MCs w/o TM
    {
        pre_selections = trigger+rectangular_cuts;
    }
    else if((species == 1000) || (species == 1100) || (species == 1200) || (species == 1300) || (species == 1500)) // Cocktail MCs w/ TM cuts
    {
        pre_selections = truthMatch;
    }
    else if(species == 7) // DDs MC
    {
        pre_selections = truthMatch+trigger+rectangular_cuts+pass_mass_fit+fit_range;
    }
    else if((species == 8) || (species == 72)) // DDs data and DDs MC w/o TM cuts
    {
        pre_selections = trigger+rectangular_cuts+pass_mass_fit+fit_range;
    }
    else if(species == 71) // D0bar Ds+ MC without rectangular cuts for the MVA validation
    {
        TCut stripping_cuts = "(D0bar_pi_ProbNNpi_pidgen_default > 0.55) && (Dsp_pi_ProbNNpi_pidgen_default > 0.55)";
        pre_selections = stripping_cuts+truthMatch+trigger+pass_mass_fit;
    }
    else if(species == 81) // D0bar Ds+ data without the rectangular cuts for MVA validation
    {
        pre_selections = trigger+"(Bp_dtf_status[0]==0)";
    }
    else if((species == 4) || (species == 9) ) // D+D-K+, D0barD+s, D0D0K MC
    {
        pre_selections = truthMatch+trigger+rectangular_cuts+pass_mass_fit;
    }
    else if((species == 5) || (species == 6) || (species == 0) || (species == -1)) // D+D-K+ and D0barD+s data
    {
        pre_selections = trigger+rectangular_cuts+pass_mass_fit;
    }

    TFileCollection* fc = new TFileCollection("fc", "fc", FILES, 1, line);
    TChain* t;
    if(is_KtautauMC || (species == 7) || (species == 71) || (species == 72) || is_cocktailMC || (species == 1000) || (species == 1100) || (species == 1200) || (species == 1300) || (species == 1500)) // Ktautau MC + normalisation channel (PID corrected)
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


    if(is_KtautauMC || is_cocktailMC)
    {
        TFileCollection* fc1 = new TFileCollection("fc1", "fc1", MC_component_file, 1, line);
        TChain* t1 = new TChain("DecayTree");
        t1->AddFileInfoList((TCollection*)fc1->GetList());

        Int_t Nraw_entries = t->GetEntries();
        Int_t Nbranch_entries = t1->GetEntries();

        cout << "From grid: " << Nraw_entries << endl;
        cout << "MC compo: " << Nbranch_entries << endl;
        if(Nraw_entries != Nbranch_entries)
        {
            cout << "Wrong number of entries" << endl;
            return;
        }
        t->AddFriend(t1, "MC_component");
    }

    TString fout_name = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/%i.root",year,species,line);
    ROOT::RDataFrame df(*t);
    df.Filter(pre_selections.GetTitle()).Snapshot("DecayTree", fout_name);
    
    cout << "Finished successfully" << endl;
}

