
void create_cocktail_MC_component_branch(int year, int species, int line)
{
    struct event
    {
        UInt_t runNumber;
        ULong64_t eventNumber;
    };

    TString MC_files;
    if(species == 100)
    {
        MC_files = Form("Files_on_grid/MC_201%i_BuDDKp_cocktail.txt",year);
    }
    else if(species == 110)
    {
        MC_files = Form("Files_on_grid/MC_201%i_BdDDKp_cocktail.txt",year);
    }
    else if(species == 120)
    {
        MC_files = Form("Files_on_grid/MC_201%i_BsDDKp_cocktail.txt",year);   
    }
    else if(species == 130)
    {
        MC_files = Form("Files_on_grid/MC_201%i_BuDDK0_cocktail.txt",year);   
    }
    else if(species == 150)
    {
        MC_files = Form("Files_on_grid/MC_201%i_BuDD_cocktail.txt",year);   
    }   
    else
    {
        cout << "Wrong species name for cocktail MC" << endl;
        return;
    }

    Int_t N;
    if(species == 100){N = 4*3;}
    else if(species == 150){N = 4*2;}
    else{N = 4;}

    TString names[N];
    if(species == 100){ 
        names[0] = "mc_ntuple_BuD0D0Kp"; names[1] = "mc_ntuple_BuD0starD0Kp"; names[2] = "mc_ntuple_BuD0D0starKp"; names[3] = "mc_ntuple_BuD0starD0starKp"; 
        names[4] = "mc_ntuple_BuDpDmKp"; names[5] = "mc_ntuple_BuDpstarDmKp"; names[6] = "mc_ntuple_BuDpDmstarKp"; names[7] = "mc_ntuple_BuDpstarDmstarKp";
        names[8] = "mc_ntuple_BuDsDsKp"; names[9] = "mc_ntuple_BuDsstarDsKp"; names[10] = "mc_ntuple_BuDsDsstarKp"; names[11] = "mc_ntuple_BuDsstarDsstarKp";
    }   
    else if(species == 110){ names[0] = "mc_ntuple_BdDmD0Kp"; names[1] = "mc_ntuple_BdDmstarD0Kp"; names[2] = "mc_ntuple_BdDmD0starKp"; names[3] = "mc_ntuple_BdDmstarD0starKp"; }
    else if(species == 120){ names[0] = "mc_ntuple_BsDsD0Kp"; names[1] = "mc_ntuple_BsDsstarD0Kp"; names[2] = "mc_ntuple_BsDsD0starKp"; names[3] = "mc_ntuple_BsDsstarD0starKp"; }
    else if(species == 130){ names[0] = "mc_ntuple_BuD0DpK0"; names[1] = "mc_ntuple_BuD0starDpK0"; names[2] = "mc_ntuple_BuD0DpstarK0"; names[3] = "mc_ntuple_BuD0starDpstarK0"; }
    else if(species == 150){ 
        names[0] = "mc_ntuple_BuD0Ds"; names[1] = "mc_ntuple_BuD0starDs"; names[2] = "mc_ntuple_BuD0Dsstar"; names[3] = "mc_ntuple_BuD0starDsstar"; 
        names[4] = "mc_ntuple_BuD0Dp"; names[5] = "mc_ntuple_BuD0starDp"; names[6] = "mc_ntuple_BuD0Dpstar"; names[7] = "mc_ntuple_BuD0starDpstar"; 
    }

    TFileCollection* fc = new TFileCollection("MC", "MC", MC_files, 1, line);
    TChain* t_reco = new TChain("ntuple/DecayTree");
    t_reco->AddFileInfoList((TCollection*)fc->GetList());

    std::vector<TChain*> t_gen;
    for(int i = 0; i < N; i++)
    {
        t_gen.push_back( new TChain(names[i]+"/MCDecayTree") );
        t_gen[i]->AddFileInfoList((TCollection*)fc->GetList());
    }

    UInt_t runNumber_gen[N];
    ULong64_t eventNumber_gen[N];
    UInt_t gen_entries[N];
    std::vector< std::vector<event> > all_gen_events;

    for(int i = 0; i < N; i++)
    {
        t_gen[i]->SetBranchAddress("runNumber", &runNumber_gen[i]);
        t_gen[i]->SetBranchAddress("eventNumber", &eventNumber_gen[i]);
        gen_entries[i] = t_gen[i]->GetEntries();

        cout << Form("Looping over %i-th GEN tree",i) << endl;
        cout << gen_entries[i] << endl;

        std::vector<event> events_gen;
        for(int evt = 0; evt < gen_entries[i]; evt++)
        {
            t_gen[i]->GetEntry(evt);

            event e;
            e.runNumber = runNumber_gen[i];
            e.eventNumber = eventNumber_gen[i];

            events_gen.push_back(e);
        }
        all_gen_events.push_back(events_gen);

    }

    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201%i/Species_%i/%i.root",year,species,line), "RECREATE");
    TTree* newtree = new TTree("DecayTree", "DecayTree");

    Int_t component = -999;
    Int_t decay = -999;
    newtree->Branch("component", &component);
    newtree->Branch("species", &decay);

    cout << "Looping over RECO" << endl;
    UInt_t runNumber;
    ULong64_t eventNumber;
    t_reco->SetBranchAddress("runNumber", &runNumber);
    t_reco->SetBranchAddress("eventNumber", &eventNumber);

    UInt_t reco_entries = t_reco->GetEntries();
    cout << reco_entries << endl;
    for(int i = 0; i < reco_entries; i++)
    {
        component = -999;
        decay = -999;

        t_reco->GetEntry(i);

        for(int i = 0; i < N; i++)
        {
            for(event e: all_gen_events[i])
            {
                if(e.runNumber == runNumber && e.eventNumber == eventNumber)
                {
                    component = i%4;
                    decay = species + N/4 - int(N/4 - i/4);
                    break;
                }
            }
            if(component == i%4) continue;
        }

        newtree->Fill();

    }
    newtree->Write();
    fout->Close();

    cout << "Finished successfully" << endl;
}