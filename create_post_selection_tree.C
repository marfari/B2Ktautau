

void create_post_selection_tree(Int_t species)
{
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

    // BDT response
    TFileCollection *fc3_2016 = new TFileCollection("fc3_2016", "fc3_2016", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_%i/bdt_output.txt",species));
    TFileCollection *fc3_2017 = new TFileCollection("fc3_2017", "fc3_2017", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_%i/bdt_output.txt",species));
    TFileCollection *fc3_2018 = new TFileCollection("fc3_2018", "fc3_2018", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_%i/bdt_output.txt",species));

    TChain* t3_2016 = new TChain("XGBoost/DecayTree");
    TChain* t3_2017 = new TChain("XGBoost/DecayTree");
    TChain* t3_2018 = new TChain("XGBoost/DecayTree");

    t3_2016->AddFileInfoList((TCollection*)fc3_2016->GetList());
    t3_2017->AddFileInfoList((TCollection*)fc3_2017->GetList());
    t3_2018->AddFileInfoList((TCollection*)fc3_2018->GetList());

    cout << "BDT response" << endl;
    cout << t3_2016->GetEntries() << endl;
    cout << t3_2017->GetEntries() << endl;
    cout << t3_2018->GetEntries() << endl;

    t3_2016->Add(t3_2017);
    t3_2016->Add(t3_2018);

    // Mass combinations
    TFileCollection *fc4_2016 = new TFileCollection("fc4_2016", "fc4_2016", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_%i/invariant_mass_tree.txt",species));
    TFileCollection *fc4_2017 = new TFileCollection("fc4_2017", "fc4_2017", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_%i/invariant_mass_tree.txt",species));
    TFileCollection *fc4_2018 = new TFileCollection("fc4_2018", "fc4_2018", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_%i/invariant_mass_tree.txt",species));

    TChain* t4_2016 = new TChain("DecayTree");
    TChain* t4_2017 = new TChain("DecayTree");
    TChain* t4_2018 = new TChain("DecayTree");

    t4_2016->AddFileInfoList((TCollection*)fc4_2016->GetList());
    t4_2017->AddFileInfoList((TCollection*)fc4_2017->GetList());
    t4_2018->AddFileInfoList((TCollection*)fc4_2018->GetList());

    cout << "Mass combinations" << endl;
    cout << t4_2016->GetEntries() << endl;
    cout << t4_2017->GetEntries() << endl;
    cout << t4_2018->GetEntries() << endl;

    t4_2016->Add(t4_2017);
    t4_2016->Add(t4_2018);

    t1_2016->AddFriend(t2_2016, "gsl");
    t1_2016->AddFriend(t3_2016, "bdt");
    t1_2016->AddFriend(t4_2016, "mass");

    TCut pass_fitter = "(df_status==0)";
    TCut bdt_cuts = "(BDT1 > 0.968) && (BDT2 > 0.953)";
    TCut cuts = pass_fitter+bdt_cuts;

    ROOT::RDataFrame df(*t1_2016);
    df.Filter(cuts.GetTitle()).Snapshot("DecayTree", Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_%i/post_sel_tree.root",species));

}