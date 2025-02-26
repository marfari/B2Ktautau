void create_table(Int_t species, TCut pass_fitter, TCut fit_region, TCut bdt_cuts, TCut mass_vetoes);
Double_t eps_error(Double_t Num, Double_t Den);

void create_post_selection_tree(Int_t species, Double_t BDT1, Double_t BDT2, bool createTable)
{
    // Creates post selection tree fot Ktautau
    TCut pass_fitter = "(df_status==0)";
    TCut fit_region = "(df_Bp_M > 4000) && (df_Bp_M < 8000)";
    TCut mass_vetoes = "";
    if((BDT1 != 0.8) && (BDT2 != 0.8)) // no BDT cut
    {   
        cout << "Applying mass vetoe cuts" << endl;
        mass_vetoes += "(TMath::Abs(Bp_M02-1864.84) > 20) && (TMath::Abs(Bp_M04-1864.84) > 20) && (TMath::Abs(Bp_M06-1864.84) > 20)"; // 2 particles: D0
        mass_vetoes += "(TMath::Abs(Bp_M046-1869.66) > 20)"; // 3 particles: D+
        mass_vetoes += "(TMath::Abs(Bp_M0456-1864.84) > 30)"; // 4 particles: D0
        mass_vetoes += "(TMath::Abs(Bp_M02456-2010.26) > 30)"; // 5 particle: D*-
    }
    TCut bdt_cuts = Form("(BDT1 > %f) && (BDT2 > %f)",BDT1,BDT2);
    TCut cuts = pass_fitter+fit_region+mass_vetoes+bdt_cuts;
    if(species != 1)
    {
        cuts += "is_best_cand == 1";
    }

    if(createTable)
    {
        create_table(species, pass_fitter, fit_region, bdt_cuts, mass_vetoes);
        return;
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
    
    if(species != 1) // If it is not the truth-match MC
    {
        // Best candidate 
        TFileCollection *fc5_2016 = new TFileCollection("fc5_2016", "fc5_2016", Form("/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_%i/multiple_events.txt",species));
        TFileCollection *fc5_2017 = new TFileCollection("fc5_2017", "fc5_2017", Form("/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_%i/multiple_events.txt",species));
        TFileCollection *fc5_2018 = new TFileCollection("fc5_2018", "fc5_2018", Form("/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_%i/multiple_events.txt",species));

        TChain* t5_2016 = new TChain("DecayTree");
        TChain* t5_2017 = new TChain("DecayTree");
        TChain* t5_2018 = new TChain("DecayTree");

        t5_2016->AddFileInfoList((TCollection*)fc5_2016->GetList());
        t5_2017->AddFileInfoList((TCollection*)fc5_2017->GetList());
        t5_2018->AddFileInfoList((TCollection*)fc5_2018->GetList());

        cout << "Best candidate" << endl;
        cout << t5_2016->GetEntries() << endl;
        cout << t5_2017->GetEntries() << endl;
        cout << t5_2018->GetEntries() << endl;

        t5_2016->Add(t5_2017);
        t5_2016->Add(t5_2018);

        t1_2016->AddFriend(t5_2016, "best_cand");
    }

    ROOT::RDataFrame df(*t1_2016);
    if((BDT1 == 0) && (BDT2 == 0))
    {
        df.Filter(cuts.GetTitle()).Snapshot("DecayTree", Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_%i/post_sel_tree_bdt1_%.0f_bdt2_%.0f.root",species,BDT1,BDT2));
    }
    else
    {
        df.Filter(cuts.GetTitle()).Snapshot("DecayTree", Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_%i/post_sel_tree_bdt1_%.1f_bdt2_%.1f.root",species,BDT1,BDT2));
    }

}

void create_table(Int_t species, TCut pass_fitter, TCut fit_region, TCut bdt_cuts, TCut mass_vetoes)
{
    Bool_t isBuDDKp_cocktail = false;
    Bool_t isBdDDKp_cocktail = false;
    Bool_t isBsDDKp_cocktail = false;
    Bool_t isBuDDK0_cocktail = false;
    Bool_t isBdDDK0_cocktail = false;
    Bool_t isBuDD_cocktail = false;
    Bool_t isBdDD_cocktail = false;
    Bool_t isBsDD_cocktail = false;

    if(species == 100)
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
    if(species == 140)
    {
        isBdDDK0_cocktail = true;
    }
    if(species == 150)
    {
        isBuDD_cocktail = true;
    }
    if(species == 160)
    {
        isBdDD_cocktail = true;
    }
    if(species == 170)
    {
        isBsDD_cocktail = true;
    }

    Bool_t is_cocktailMC = false;
    if(isBuDDKp_cocktail || isBdDDKp_cocktail || isBsDDKp_cocktail || isBuDDK0_cocktail || isBdDDK0_cocktail || isBuDD_cocktail || isBdDD_cocktail || isBsDD_cocktail)
    {
        is_cocktailMC = true;
    }

    for(int year = 6; year < 9; year++)
    {
        cout << Form("Creating pre-selection efficiency tables for year %i", year) << endl;
        std::ofstream file(Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_%i/post_sel_table_201%i.tex",species,year));

        // Pre-sel tree
        TFileCollection *fc1 = new TFileCollection("fc1", "fc1", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/pre_sel_tree.txt",year,species));
        TChain* t1 = new TChain("DecayTree");
        t1->AddFileInfoList((TCollection*)fc1->GetList());

        // GSL tree
        TFileCollection *fc2 = new TFileCollection("fc2", "fc2", Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/fit_results.txt",year,species));
        TChain* t2 = new TChain("DecayTree");
        t2->AddFileInfoList((TCollection*)fc2->GetList());


        // BDT tree
        TFileCollection *fc3 = new TFileCollection("fc3", "fc3", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/201%i/Species_%i/bdt_output.txt",year,species));
        TChain* t3 = new TChain("XGBoost/DecayTree");
        t3->AddFileInfoList((TCollection*)fc3->GetList());


        // Invariant mass tree
        TFileCollection *fc4 = new TFileCollection("fc4", "fc4", Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201%i/Species_%i/invariant_mass_tree.txt",year,species));
        TChain* t4 = new TChain("DecayTree");
        t4->AddFileInfoList((TCollection*)fc4->GetList());

        t1->AddFriend(t2);
        t1->AddFriend(t3);
        t1->AddFriend(t4);

        // Gen tree
        TChain *t5;
        TChain *t5_3pi3pi, *t5_3pi3pi_pi0, *t5_3pi3pi2pi0;
        TChain *t5_DD, *t5_DstarD, *t5_DDstar, *t5_DstarDstar;
        if((species == 1) || is_cocktailMC)
        {
            TFileCollection *fc5;
            if(species == 1)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i.txt",year));

                // 3pi3pi
                t5_3pi3pi = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
                t5_3pi3pi->AddFileInfoList((TCollection*)fc5->GetList());
                
                // 3pi3pi pi0
                t5_3pi3pi_pi0 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
                TChain* t5_3pi_pi0_3pi = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
                t5_3pi3pi_pi0->AddFileInfoList((TCollection*)fc5->GetList());
                t5_3pi_pi0_3pi->AddFileInfoList((TCollection*)fc5->GetList());
                t5_3pi3pi_pi0->GetEntries();
                t5_3pi_pi0_3pi->GetEntries();
                t5_3pi3pi_pi0->Add(t5_3pi_pi0_3pi);
                
                // 3pi3pi 2pi0
                t5_3pi3pi2pi0 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
                t5_3pi3pi2pi0->AddFileInfoList((TCollection*)fc5->GetList());

                // All
                t5 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
                t5->AddFileInfoList((TCollection*)fc5->GetList());

                t5->GetEntries();
                t5->Add(t5_3pi3pi_pi0);
                t5->Add(t5_3pi3pi2pi0);
            }
            else if(isBuDDKp_cocktail)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i_BuDDKp_cocktail.txt",year));
                if(species == 100)
                {
                    t5_DD = new TChain("mc_ntuple_BuD0D0Kp/MCDecayTree");
                    t5_DstarD = new TChain("mc_ntuple_BuD0starD0Kp/MCDecayTree");
                    t5_DDstar = new TChain("mc_ntuple_BuD0D0starKp/MCDecayTree");
                    t5_DstarDstar = new TChain("mc_ntuple_BuD0starD0starKp/MCDecayTree");
                    t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
                }
                else if(species == 101)
                {
                    t5_DD = new TChain("mc_ntuple_BuDpDmKp/MCDecayTree");
                    t5_DstarD = new TChain("mc_ntuple_BuDpstarDmKp/MCDecayTree");
                    t5_DDstar = new TChain("mc_ntuple_BuDpDmstarKp/MCDecayTree");
                    t5_DstarDstar = new TChain("mc_ntuple_BuDpstarDmstarKp/MCDecayTree");
                    t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
                }
                else if(species == 102)
                {
                    t5_DD = new TChain("mc_ntuple_BuDsDsKp/MCDecayTree");
                    t5_DstarD = new TChain("mc_ntuple_BuDsstarDsKp/MCDecayTree");
                    t5_DDstar = new TChain("mc_ntuple_BuDsDsstarKp/MCDecayTree");
                    t5_DstarDstar = new TChain("mc_ntuple_BuDsstarDsstarKp/MCDecayTree");
                    t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
                }

            }
            else if(isBdDDKp_cocktail)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i_BdDDKp_cocktail.txt",year));
                t5_DD = new TChain("mc_ntuple_BdDmD0Kp/MCDecayTree");
                t5_DstarD = new TChain("mc_ntuple_BdDmstarD0Kp/MCDecayTree");
                t5_DDstar = new TChain("mc_ntuple_BdDmD0starKp/MCDecayTree");
                t5_DstarDstar = new TChain("mc_ntuple_BdDmstarD0starKp/MCDecayTree");
                t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
            }
            else if(isBsDDKp_cocktail)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i_BsDDKp_cocktail.txt",year));
                t5_DD = new TChain("mc_ntuple_BsDsD0Kp/MCDecayTree");
                t5_DstarD = new TChain("mc_ntuple_BsDsstarD0Kp/MCDecayTree");
                t5_DDstar = new TChain("mc_ntuple_BsDsD0starKp/MCDecayTree");
                t5_DstarDstar = new TChain("mc_ntuple_BsDsstarD0starKp/MCDecayTree");
                t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
            }
            else if(isBuDDK0_cocktail)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i_BuDDK0_cocktail.txt",year));
                t5_DD = new TChain("mc_ntuple_BuD0DpK0/MCDecayTree");
                t5_DstarD = new TChain("mc_ntuple_BuD0starDpK0/MCDecayTree");
                t5_DDstar = new TChain("mc_ntuple_BuD0DpstarK0/MCDecayTree");
                t5_DstarDstar = new TChain("mc_ntuple_BuD0starDpstarK0/MCDecayTree");
                t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
            }
            else if(isBuDD_cocktail)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i_BuDD_cocktail.txt",year));
                if(species == 150)
                {
                    t5_DD = new TChain("mc_ntuple_BuD0Ds/MCDecayTree");
                    t5_DstarD = new TChain("mc_ntuple_BuD0starDs/MCDecayTree");
                    t5_DDstar = new TChain("mc_ntuple_BuD0Dsstar/MCDecayTree");
                    t5_DstarDstar = new TChain("mc_ntuple_BuD0starDsstar/MCDecayTree");
                    t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
                }
                else if(species == 151)
                {
                    t5_DD = new TChain("mc_ntuple_BuD0Dp/MCDecayTree");
                    t5_DstarD = new TChain("mc_ntuple_BuD0starDp/MCDecayTree");
                    t5_DDstar = new TChain("mc_ntuple_BuD0Dpstar/MCDecayTree");
                    t5_DstarDstar = new TChain("mc_ntuple_BuD0starDpstar/MCDecayTree");
                    t5_DD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarD->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DDstar->AddFileInfoList((TCollection*)fc5->GetList());
                    t5_DstarDstar->AddFileInfoList((TCollection*)fc5->GetList());
                }
            }
        }

        Double_t N_initial;
        Double_t N_initial_3pi3pi, N_initial_3pi3pipi0, N_initial_3pi3pi2pi0;
        if(species == 1)
        {
            N_initial_3pi3pi = t5_3pi3pi->GetEntries();
            N_initial_3pi3pipi0 = t5_3pi3pi_pi0->GetEntries();
            N_initial_3pi3pi2pi0 = t5_3pi3pi2pi0->GetEntries();
            N_initial = t5->GetEntries();
        }

        Double_t N_presel;
        Double_t N_presel_3pi3pi, N_presel_3pi3pipi0, N_presel_3pi3pi2pi0;
        if(species == 1)
        {
            N_presel_3pi3pi = t1->GetEntries("component==0");
            N_presel_3pi3pipi0 = t1->GetEntries("component==1");
            N_presel_3pi3pi2pi0 = t1->GetEntries("component==2");
        }
        N_presel = t1->GetEntries();

        Double_t N_gsl;
        Double_t N_gsl_3pi3pi, N_gsl_3pi3pipi0, N_gsl_3pi3pi2pi0;
        Double_t eps_gsl;
        Double_t eps_gsl_3pi3pi, eps_gsl_3pi3pipi0, eps_gsl_3pi3pi2pi0;
        if(species == 1)
        {
            N_gsl_3pi3pi = t1->GetEntries(pass_fitter+"component==0");
            N_gsl_3pi3pipi0 = t1->GetEntries(pass_fitter+"component==1");
            N_gsl_3pi3pi2pi0 = t1->GetEntries(pass_fitter+"component==2");

            eps_gsl_3pi3pi = N_gsl_3pi3pi/N_presel_3pi3pi;
            eps_gsl_3pi3pipi0 = N_gsl_3pi3pipi0/N_presel_3pi3pipi0;
            eps_gsl_3pi3pi2pi0 = N_gsl_3pi3pi2pi0/N_presel_3pi3pi2pi0;
        }
        N_gsl = t1->GetEntries(pass_fitter);
        eps_gsl = N_gsl/N_presel;

        Double_t N_fit_region;
        Double_t N_fit_region_3pi3pi, N_fit_region_3pi3pipi0, N_fit_region_3pi3pi2pi0;
        Double_t eps_fit_region;
        Double_t eps_fit_region_3pi3pi, eps_fit_region_3pi3pipi0, eps_fit_region_3pi3pi2pi0;
        if(species == 1)
        {
            N_fit_region_3pi3pi = t1->GetEntries(pass_fitter+fit_region+"component==0");
            N_fit_region_3pi3pipi0 = t1->GetEntries(pass_fitter+fit_region+"component==1");
            N_fit_region_3pi3pi2pi0 = t1->GetEntries(pass_fitter+fit_region+"component==2");

            eps_fit_region_3pi3pi = N_fit_region_3pi3pi/N_gsl_3pi3pi;
            eps_fit_region_3pi3pipi0 = N_fit_region_3pi3pipi0/N_gsl_3pi3pipi0;
            eps_fit_region_3pi3pi2pi0 = N_fit_region_3pi3pi2pi0/N_gsl_3pi3pi2pi0;
        }
        N_fit_region = t1->GetEntries(pass_fitter+fit_region);
        eps_fit_region = N_fit_region/N_gsl;

        Double_t N_mass_vetoes;
        Double_t N_mass_vetoes_3pi3pi, N_mass_vetoes_3pi3pipi0, N_mass_vetoes_3pi3pi2pi0;
        Double_t eps_mass;
        Double_t eps_mass_3pi3pi, eps_mass_3pi3pipi0, eps_mass_3pi3pi2pi0;
        if(species == 1)
        {
            N_mass_vetoes_3pi3pi = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+"component==0");
            N_mass_vetoes_3pi3pipi0 = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+"component==1");
            N_mass_vetoes_3pi3pi2pi0 = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+"component==2");

            eps_mass_3pi3pi = N_mass_vetoes_3pi3pi/N_fit_region_3pi3pi;
            eps_mass_3pi3pipi0 = N_mass_vetoes_3pi3pipi0/N_fit_region_3pi3pipi0;
            eps_mass_3pi3pi2pi0 = N_mass_vetoes_3pi3pi2pi0/N_fit_region_3pi3pi2pi0;
        }
        N_mass_vetoes = t1->GetEntries(pass_fitter+fit_region+mass_vetoes);
        eps_mass = N_mass_vetoes/N_fit_region;

        Double_t N_bdt;
        Double_t N_bdt_3pi3pi, N_bdt_3pi3pipi0, N_bdt_3pi3pi2pi0;
        Double_t eps_bdt;
        Double_t eps_bdt_3pi3pi, eps_bdt_3pi3pipi0, eps_bdt_3pi3pi2pi0;
        if(species == 1)
        {
            N_bdt_3pi3pi = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+bdt_cuts+"component==0");
            N_bdt_3pi3pipi0 = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+bdt_cuts+"component==1");
            N_bdt_3pi3pi2pi0 = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+bdt_cuts+"component==2");

            eps_bdt_3pi3pi = N_bdt_3pi3pi/N_mass_vetoes_3pi3pi;
            eps_bdt_3pi3pipi0 = N_bdt_3pi3pipi0/N_mass_vetoes_3pi3pipi0;
            eps_bdt_3pi3pi2pi0 = N_bdt_3pi3pi2pi0/N_mass_vetoes_3pi3pi2pi0; 
        }
        N_bdt = t1->GetEntries(pass_fitter+fit_region+mass_vetoes+bdt_cuts);
        eps_bdt = N_bdt/N_mass_vetoes;

        Double_t eps_post_strip;
        Double_t eps_post_strip_3pi3pi, eps_post_strip_3pi3pipi0, eps_post_strip_3pi3pi2pi0;
        if(species == 1)
        {
            eps_post_strip_3pi3pi = N_bdt_3pi3pi/N_initial_3pi3pi;
            eps_post_strip_3pi3pipi0 = N_bdt_3pi3pipi0/N_initial_3pi3pipi0;
            eps_post_strip_3pi3pi2pi0 = N_bdt_3pi3pi2pi0/N_initial_3pi3pi2pi0;
        }
        eps_post_strip = N_bdt/N_initial;

        Double_t acc_error, strip_error;
        Double_t eps_acc, eps_strip;
        if((species == 1) && (year == 6))
        {
            eps_acc = (5.318/100);
            eps_strip = (1.028/100);
            acc_error = (0.010/100);
            strip_error = (0.003/100);
        }
        else if((species == 1) && (year == 7))
        {
            eps_acc = (5.321/100);
            eps_strip = (1.056/100);
            acc_error = (0.011/100);
            strip_error = (0.003/100);
        }
        else if((species == 1) && (year == 8))
        {
            eps_acc = (5.330/100);
            eps_strip = (1.050/100);
            acc_error = (0.010/100);
            strip_error = (0.002/100);
        }

        Double_t eps_total;
        Double_t eps_total_3pi3pi, eps_total_3pi3pipi0, eps_total_3pi3pi2pi0;
        if(species == 1)
        {
            eps_total_3pi3pi = eps_post_strip_3pi3pi*eps_acc*eps_strip;
            eps_total_3pi3pipi0 = eps_post_strip_3pi3pipi0*eps_acc*eps_strip;
            eps_total_3pi3pi2pi0 = eps_post_strip_3pi3pi2pi0*eps_acc*eps_strip;
        }
        eps_total = eps_post_strip*eps_acc*eps_strip;

        file << " \\begin{table}[!htbp]" << std::endl;
        file << " \\centering " << std::endl;
        file << " \\begin{tabular}{|c|c|}" << std::endl;
        file << " \\hline" << std::endl;
        file << " Cut & " << "Efficiency " << " \\\\ " << std::endl;
        file << "\\hline" << std::endl;

        if(species == 1)
        {
            file << Form("Pass GSL (3pi3pi) & %.2lf $\\pm$ %.2lf \\%% ", eps_gsl_3pi3pi*100,  eps_error(N_gsl_3pi3pi,N_presel_3pi3pi) ) << " \\\\ " << std::endl;
            file << Form("Pass GSL (3pi3pi\\_pi0) & %.2lf $\\pm$ %.2lf \\%% ", eps_gsl_3pi3pipi0*100,  eps_error(N_gsl_3pi3pipi0,N_presel_3pi3pipi0) ) << " \\\\ " << std::endl;
            file << Form("Pass GSL (3pi3pi\\_2pi0) & %.2lf $\\pm$ %.2lf \\%% ", eps_gsl_3pi3pi2pi0*100,  eps_error(N_gsl_3pi3pi2pi0,N_presel_3pi3pi2pi0) ) << " \\\\ " << std::endl;
            file << Form("Pass GSL (All) & %.2lf $\\pm$ %.2lf \\%% ", eps_gsl*100,  eps_error(N_gsl,N_presel) ) << " \\\\ \\hline" << std::endl;

            file << Form("Fit region (3pi3pi) & %.2lf $\\pm$ %.2lf \\%% ", eps_fit_region_3pi3pi*100,  eps_error(N_fit_region_3pi3pi,N_gsl_3pi3pi) ) << " \\\\" << std::endl;
            file << Form("Fit region (3pi3pi\\_pi0) & %.2lf $\\pm$ %.2lf \\%% ", eps_fit_region_3pi3pipi0*100,  eps_error(N_fit_region_3pi3pipi0,N_gsl_3pi3pipi0) ) << " \\\\" << std::endl;
            file << Form("Fit region (3pi3pi\\_2pi0) & %.2lf $\\pm$ %.2lf \\%% ", eps_fit_region_3pi3pi2pi0*100,  eps_error(N_fit_region_3pi3pi2pi0,N_gsl_3pi3pi2pi0) ) << " \\\\" << std::endl;
            file << Form("Fit region (All) & %.2lf $\\pm$ %.2lf \\%% ", eps_fit_region*100,  eps_error(N_fit_region,N_gsl) ) << " \\\\ \\hline" << std::endl;

            file << Form("Mass vetoes (3pi3pi) & %.2lf $\\pm$ %.2lf \\%% ", eps_mass_3pi3pi*100, eps_error(N_mass_vetoes_3pi3pi, N_fit_region_3pi3pi) ) << " \\\\" << std::endl;
            file << Form("Mass vetoes (3pi3pi\\_pi0) & %.2lf $\\pm$ %.2lf \\%% ", eps_mass_3pi3pipi0*100, eps_error(N_mass_vetoes_3pi3pipi0, N_fit_region_3pi3pipi0) ) << " \\\\" << std::endl;
            file << Form("Mass vetoes (3pi3pi\\_2pi0) & %.2lf $\\pm$ %.2lf \\%% ", eps_mass_3pi3pi2pi0*100, eps_error(N_mass_vetoes_3pi3pi2pi0, N_fit_region_3pi3pi2pi0) ) << " \\\\" << std::endl;
            file << Form("Mass vetoes (All) & %.2lf $\\pm$ %.2lf \\%% ", eps_mass*100, eps_error(N_mass_vetoes, N_fit_region) ) << " \\\\ \\hline" << std::endl;

            // file << "BDTs (3pi3pi) & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_bdt_3pi3pi*100, eps_error(N_bdt_3pi3pi,N_mass_vetoes_3pi3pi) ) << " \\\\ " << std::endl;
            // file << "BDTs (3pi3pi\\_pi0) & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_bdt_3pi3pipi0*100, eps_error(N_bdt_3pi3pipi0,N_mass_vetoes_3pi3pipi0) ) << " \\\\ " << std::endl;
            // file << "BDTs (3pi3pi\\_2pi0) & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_bdt_3pi3pi2pi0*100, eps_error(N_bdt_3pi3pi2pi0,N_mass_vetoes_3pi3pi2pi0) ) << " \\\\ " << std::endl;
            // file << "BDTs (All) & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_bdt*100, eps_error(N_bdt,N_mass_vetoes) ) << " \\\\ \\hline" << std::endl;
            
            Double_t post_strip_error_3pi3pi = eps_error(N_mass_vetoes_3pi3pi, N_initial_3pi3pi)/100.;
            Double_t post_strip_error_3pi3pipi0 = eps_error(N_mass_vetoes_3pi3pipi0, N_initial_3pi3pipi0)/100.;
            Double_t post_strip_error_3pi3pi2pi0 = eps_error(N_mass_vetoes_3pi3pi2pi0, N_initial_3pi3pi2pi0)/100.;
            Double_t post_strip_error = eps_error(N_mass_vetoes, N_initial)/100.;

            file << Form("Total (3pi3pi) & %.6lf $\\pm$ %.6lf \\%% ", eps_total_3pi3pi*100, TMath::Sqrt( pow(eps_acc*eps_strip,2)*pow(post_strip_error_3pi3pi,2) + pow(eps_post_strip_3pi3pi*eps_strip,2)*pow(acc_error,2) + pow(eps_post_strip_3pi3pi*eps_acc,2)*pow(strip_error,2) )*100 ) << " \\\\ " << std::endl;
            file << Form("Total (3pi3pi\\_pi0) & %.6lf $\\pm$ %.6lf \\%% ", eps_total_3pi3pipi0*100, TMath::Sqrt( pow(eps_acc*eps_strip,2)*pow(post_strip_error_3pi3pipi0,2) + pow(eps_post_strip_3pi3pipi0*eps_strip,2)*pow(acc_error,2) + pow(eps_post_strip_3pi3pipi0*eps_acc,2)*pow(strip_error,2) )*100 ) << " \\\\ " << std::endl;
            file << Form("Total (3pi3pi\\_2pi0) & %.6lf $\\pm$ %.6lf \\%% ", eps_total_3pi3pi2pi0*100, TMath::Sqrt( pow(eps_acc*eps_strip,2)*pow(post_strip_error_3pi3pi2pi0,2) + pow(eps_post_strip_3pi3pi2pi0*eps_strip,2)*pow(acc_error,2) + pow(eps_post_strip_3pi3pi2pi0*eps_acc,2)*pow(strip_error,2) )*100 ) << " \\\\ " << std::endl;
            file << Form("Total (All) & %.6lf $\\pm$ %.6lf \\%% ", eps_total*100, TMath::Sqrt( pow(eps_acc*eps_strip,2)*pow(post_strip_error,2) + pow(eps_post_strip*eps_strip,2)*pow(acc_error,2) + pow(eps_post_strip*eps_acc,2)*pow(strip_error,2) )*100 ) << " \\\\ " << std::endl;

        }
        else if((species == 2) || (species == 3))
        {
            file << Form("Pass GSL & %.2lf $\\pm$ %.2lf \\%% ", eps_gsl*100,  eps_error(N_gsl,N_presel) ) << " \\\\ \\hline " << std::endl;
            file << Form("Fit region & %.2lf $\\pm$ %.2lf \\%% ", eps_fit_region*100,  eps_error(N_fit_region,N_gsl) ) << " \\\\ \\hline " << std::endl;
            file << Form("Mass vetoes & %.2lf $\\pm$ %.2lf \\%% ", eps_mass*100, eps_error(N_mass_vetoes, N_fit_region) ) << " \\\\" << std::endl;
            // file << "BDTs & " << Form("%.6lf $\\pm$ %.6lf \\% & ", eps_bdt*100, eps_error(N_bdt,N_mass_vetoes) ) << " \\\\ " << std::endl;
        }

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