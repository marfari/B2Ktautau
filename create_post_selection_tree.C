

void create_post_selection_tree(Int_t year, Int_t species, Int_t line)
{
    // Pre-selection tree
    TFileCollection *fc1 = new TFileCollection("fc1", "fc1", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/pre_sel_tree.txt",year,species), 1, line);
    TChain* t1 = new TChain("DecayTree");
    t1->AddFileInfoList((TCollection*)fc1->GetList());

    // Fit results
    TFileCollection *fc2 = new TFileCollection("fc2", "fc2", Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/fit_results.txt",year,species), 1, line);
    TChain* t2 = new TChain("DecayTree");
    t2->AddFileInfoList((TCollection*)fc2->GetList());

    // BDT response
    TFileCollection *fc3 = new TFileCollection("fc3", "fc3", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/201%i/Species_%i/bdt_output.txt",year,species), 1, line);
    TChain* t3 = new TChain("XGBoost/DecayTree");
    t3->AddFileInfoList((TCollection*)fc3->GetList());

    // Mass combinations
    TFileCollection *fc4 = new TFileCollection("fc4", "fc4", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201%i/Species_%i/invariant_mass_tree.txt",year,species), 1, line);
    TChain* t4 = new TChain("DecayTree");
    t4->AddFileInfoList((TCollection*)fc4->GetList());

    t1->AddFriend(t2, "gsl");
    t1->AddFriend(t3, "bdt");
    t1->AddFriend(t4, "mass");

    TCut pass_fitter = "(df_status==0)";
    TCut bdt_cuts = "(BDT1 > 0.968) && (BDT2 > 0.953)";
    TCut cuts = pass_fitter+bdt_cuts;

    ROOT::RDataFrame df(*t1);
    df.Filter(cuts.GetTitle()).Snapshot("DecayTree", Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/201%i/Species_%i/%i.root",year,species,line));

}