void create_table(Int_t species, TCut pass_fitter, TCut bdt_cuts, TCut mass_vetoes);
Double_t eps_error(Double_t Num, Double_t Den);

void create_post_selection_tree(Int_t species, bool createTable)
{
    // Creates post selection tree fot Ktautau
    TCut pass_fitter = "(df_status==0)";
    TCut bdt_cuts = "(BDT1 > 0.968) && (BDT2 > 0.953)"; // these cuts were optmised with the rough sensitivity study
    TCut mass_vetoes = "(TMath::Abs(Bp_M0456 - 1864.84) > 25)";
    TCut cuts = pass_fitter+bdt_cuts+mass_vetoes;

    if(createTable)
    {
        create_table(species, pass_fitter, bdt_cuts, mass_vetoes);
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

    ROOT::RDataFrame df(*t1_2016);
    df.Filter(cuts.GetTitle()).Snapshot("DecayTree", Form("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_%i/post_sel_tree.root",species));

}

void create_table(Int_t species, TCut pass_fitter, TCut bdt_cuts, TCut mass_vetoes)
{
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
        TChain *t5_DD, *t5_DstarD, *t5_DDstar, *t5_DstarDstar;
        if(isKtautauMC || is_cocktailMC)
        {
            TFileCollection *fc5;
            if(isKtautauMC)
            {
                fc5 = new TFileCollection("fc5", "fc5", Form("Files_on_grid/MC_201%i.txt",year));

                if(species == 10)
                {
                    t5 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
                    t5->AddFileInfoList((TCollection*)fc5->GetList());
                }
                else if(species == 11)
                {
                    t5 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
                    TChain* t6 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
                    t5->AddFileInfoList((TCollection*)fc5->GetList());
                    t6->AddFileInfoList((TCollection*)fc5->GetList());
                    t5->GetEntries();
                    t6->GetEntries();
                    t5->Add(t6);
                }
                else if(species == 12)
                {
                    t5 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
                    t5->AddFileInfoList((TCollection*)fc5->GetList());
                }
                else if(species == 1)
                {
                    t5 = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
                    TChain* t6 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
                    TChain* t7 = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
                    TChain* t8 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
                    t5->AddFileInfoList((TCollection*)fc5->GetList());
                    t6->AddFileInfoList((TCollection*)fc5->GetList());
                    t7->AddFileInfoList((TCollection*)fc5->GetList());
                    t8->AddFileInfoList((TCollection*)fc5->GetList());
                    t5->GetEntries();
                    t6->GetEntries();
                    t7->GetEntries();
                    t8->GetEntries();
                    t5->Add(t6);
                    t5->Add(t7);
                    t5->Add(t8);
                }
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
        if(isKtautauMC)
        {
            N_initial = t5->GetEntries();
        }
        Double_t N_presel = t1->GetEntries();
        Double_t N_gsl = t1->GetEntries(pass_fitter);
        Double_t N_bdt = t1->GetEntries(pass_fitter+bdt_cuts);
        Double_t N_mass_vetoes = t1->GetEntries(pass_fitter+bdt_cuts+mass_vetoes);

        Double_t eps_gsl = N_gsl/N_presel;
        Double_t eps_bdt = N_bdt/N_gsl;
        Double_t eps_mass = N_mass_vetoes/N_bdt;
        Double_t eps_post_strip = N_mass_vetoes/N_initial;

        Double_t acc_error, strip_error;
        Double_t eps_acc, eps_strip;
        if(isKtautauMC && (year == 6))
        {
            eps_acc = (5.318/100);
            eps_strip = (1.028/100);
            acc_error = (0.010/100);
            strip_error = (0.003/100);
        }
        else if(isKtautauMC && (year == 7))
        {
            eps_acc = (5.321/100);
            eps_strip = (1.056/100);
            acc_error = (0.011/100);
            strip_error = (0.003/100);
        }
        else if(isKtautauMC && (year == 8))
        {
            eps_acc = (5.330/100);
            eps_strip = (1.050/100);
            acc_error = (0.010/100);
            strip_error = (0.002/100);
        }

        Double_t eps_total;
        if(isKtautauMC)
        {
            eps_total = eps_post_strip*eps_acc*eps_strip;
        }

        file << " \\begin{table}[!htbp]" << std::endl;
        file << " \\centering " << std::endl;
        file << " \\begin{tabular}{|c|c|}" << std::endl;
        file << " \\hline" << std::endl;
        file << " Cut & " << "Efficiency " << " \\\\ " << std::endl;
        file << "\\hline" << std::endl;

        if(isKtautauMC || (species == 2) || (species == 3))
        {
            file << "Pass GSL & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_gsl*100,  eps_error(N_gsl,N_presel) ) << " \\\\ \\hline " << std::endl;
            
            if(isKtautauMC)
            {
                file << "BDTs & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_bdt*100, eps_error(N_bdt,N_gsl) ) << " \\\\ \\hline " << std::endl;
                file << "Mass vetoes & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_mass*100, eps_error(N_mass_vetoes, N_bdt) ) << " \\\\ \\hline" << std::endl;
                Double_t post_strip_error = eps_error(N_mass_vetoes, N_initial)/100.;
                file << "Total & " << Form("%.6lf $\\pm$ %.6lf \\% & ", eps_total*100, TMath::Sqrt( pow(eps_acc*eps_strip,2)*pow(post_strip_error,2) + pow(eps_post_strip*eps_strip,2)*pow(acc_error,2) + pow(eps_post_strip*eps_acc,2)*pow(strip_error,2) )*100 ) << " \\\\ " << std::endl;
            }
            else
            {
                file << "BDTs & " << Form("%.4lf $\\pm$ %.4lf \\% & ", eps_bdt*100, eps_error(N_bdt,N_gsl) ) << " \\\\ \\hline " << std::endl;
                file << "Mass vetoes & " << Form("%.2lf $\\pm$ %.2lf \\% & ", eps_mass*100, eps_error(N_mass_vetoes, N_bdt) ) << " \\\\ " << std::endl;
            }
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