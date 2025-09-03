
void create_post_selection_tree(Int_t species, Double_t BDT)
{
    Bool_t isKtautau = false;
    if((species == 1) || (species == 10) || (species == 2) || (species == 3))
    {
        isKtautau = true;
    }

    Bool_t is_cocktailMC = false;
    if((species == 100) || (species == 110) || (species == 120) || (species == 130) || (species == 150))
    {
        is_cocktailMC = true;
    }

    // Selections    
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

    TCut all_selections = trigger+rectangular_cuts+pass_mass_fit+fit_range;

    if(isKtautau || is_cocktailMC)
    {
        TObjString* mass_vetoes_exp = (TObjString*)t_selections->Get("mass_vetoes");
        TCut mass_vetoes(mass_vetoes_exp->GetString().Data());

        all_selections += mass_vetoes;
    }
    if((species == 1) || (species == 7))
    {
        all_selections += truthMatch;
    }
    if((species == 72) || (species == 8))
    {
        TCut best_cand("is_best_cand == 1");
        all_selections += best_cand;
    }

    TCut bdt_cuts = Form("BDT > %f",BDT);
    if(isKtautau || is_cocktailMC)
    {
        all_selections += bdt_cuts;
    }
    
    // Pre-selection tree
    TFileCollection *fc1_2016 = new TFileCollection("fc1_2016", "fc1_2016", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_%i/pre_sel_tree.txt",species));
    TFileCollection *fc1_2017 = new TFileCollection("fc1_2017", "fc1_2017", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_%i/pre_sel_tree.txt",species));
    TFileCollection *fc1_2018 = new TFileCollection("fc1_2018", "fc1_2018", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_%i/pre_sel_tree.txt",species));

    TChain* t1_2016 = new TChain("DecayTree");
    TChain* t1_2017 = new TChain("DecayTree");
    TChain* t1_2018 = new TChain("DecayTree");

    t1_2016->AddFileInfoList((TCollection*)fc1_2016->GetList());
    t1_2017->AddFileInfoList((TCollection*)fc1_2017->GetList());
    t1_2018->AddFileInfoList((TCollection*)fc1_2018->GetList());

    cout << "Pre-selection files" << endl;
    cout << t1_2016->GetEntries() << endl;
    cout << t1_2017->GetEntries() << endl;
    cout << t1_2018->GetEntries() << endl;

    t1_2016->Add(t1_2017);
    t1_2016->Add(t1_2018);

    // Fit results
    if(isKtautau || is_cocktailMC)
    {
        TFileCollection *fc2_2016 = new TFileCollection("fc2_2016", "fc2_2016", Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_%i/fit_results.txt",species));
        TFileCollection *fc2_2017 = new TFileCollection("fc2_2017", "fc2_2017", Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_%i/fit_results.txt",species));
        TFileCollection *fc2_2018 = new TFileCollection("fc2_2018", "fc2_2018", Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_%i/fit_results.txt",species));

        TChain* t2_2016 = new TChain("DecayTree");
        TChain* t2_2017 = new TChain("DecayTree");
        TChain* t2_2018 = new TChain("DecayTree");

        t2_2016->AddFileInfoList((TCollection*)fc2_2016->GetList());
        t2_2017->AddFileInfoList((TCollection*)fc2_2017->GetList());
        t2_2018->AddFileInfoList((TCollection*)fc2_2018->GetList());

        cout << "Fit results" << endl;
        cout << t2_2016->GetEntries() << endl;
        cout << t2_2017->GetEntries() << endl;
        cout << t2_2018->GetEntries() << endl;

        t2_2016->Add(t2_2017);
        t2_2016->Add(t2_2018);

        t1_2016->AddFriend(t2_2016, "gsl");

        // Mass combinations
        TFileCollection *fc5_2016 = new TFileCollection("fc5_2016", "fc4_2016", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_%i/invariant_mass_tree.txt",species));
        TFileCollection *fc5_2017 = new TFileCollection("fc5_2017", "fc4_2017", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_%i/invariant_mass_tree.txt",species));
        TFileCollection *fc5_2018 = new TFileCollection("fc5_2018", "fc4_2018", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_%i/invariant_mass_tree.txt",species));

        TChain* t5_2016 = new TChain("DecayTree");
        TChain* t5_2017 = new TChain("DecayTree");
        TChain* t5_2018 = new TChain("DecayTree");

        t5_2016->AddFileInfoList((TCollection*)fc5_2016->GetList());
        t5_2017->AddFileInfoList((TCollection*)fc5_2017->GetList());
        t5_2018->AddFileInfoList((TCollection*)fc5_2018->GetList());

        cout << "Mass combinations" << endl;
        cout << t5_2016->GetEntries() << endl;
        cout << t5_2017->GetEntries() << endl;
        cout << t5_2018->GetEntries() << endl;

        t5_2016->Add(t5_2017);
        t5_2016->Add(t5_2018);

        t1_2016->AddFriend(t5_2016, "mass");
    }

    // BDT response
    TFileCollection *fc3_2016 = new TFileCollection("fc3_2016", "fc3_2016", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_%i/bdt_output.txt",species));
    TFileCollection *fc3_2017 = new TFileCollection("fc3_2017", "fc3_2017", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_%i/bdt_output.txt",species));
    TFileCollection *fc3_2018 = new TFileCollection("fc3_2018", "fc3_2018", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_%i/bdt_output.txt",species));

    TChain* t3_2016 = new TChain("DecayTree");
    TChain* t3_2017 = new TChain("DecayTree");
    TChain* t3_2018 = new TChain("DecayTree");

    t3_2016->AddFileInfoList((TCollection*)fc3_2016->GetList());
    t3_2017->AddFileInfoList((TCollection*)fc3_2017->GetList());
    t3_2018->AddFileInfoList((TCollection*)fc3_2018->GetList());

    cout << "BDT response" << endl;
    cout << t3_2016->GetEntries() << endl;
    cout << t3_2017->GetEntries() << endl;
    cout << t3_2018->GetEntries() << endl;

    t3_2016->Add(t3_2017);
    t3_2016->Add(t3_2018);

    t1_2016->AddFriend(t3_2016, "bdt");

    if((species == 72) || (species == 8)) 
    {
        // Best candidate 
        TFileCollection *fc6_2016 = new TFileCollection("fc6_2016", "fc6_2016", Form("/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_%i/multiple_events.txt",species));
        TFileCollection *fc6_2017 = new TFileCollection("fc6_2017", "fc6_2017", Form("/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_%i/multiple_events.txt",species));
        TFileCollection *fc6_2018 = new TFileCollection("fc6_2018", "fc6_2018", Form("/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_%i/multiple_events.txt",species));

        TChain* t6_2016 = new TChain("DecayTree");
        TChain* t6_2017 = new TChain("DecayTree");
        TChain* t6_2018 = new TChain("DecayTree");

        t6_2016->AddFileInfoList((TCollection*)fc6_2016->GetList());
        t6_2017->AddFileInfoList((TCollection*)fc6_2017->GetList());
        t6_2018->AddFileInfoList((TCollection*)fc6_2018->GetList());

        cout << "Best candidate" << endl;
        cout << t6_2016->GetEntries() << endl;
        cout << t6_2017->GetEntries() << endl;
        cout << t6_2018->GetEntries() << endl;

        t6_2016->Add(t6_2017);
        t6_2016->Add(t6_2018);

        t1_2016->AddFriend(t6_2016, "best_cand");
    }

    ROOT::RDataFrame df(*t1_2016);
    df.Filter(all_selections.GetTitle()).Snapshot("DecayTree", Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_%i/post_sel_tree_bdt_%.4f.root",species,BDT));
    
}