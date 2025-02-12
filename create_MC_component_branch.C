
void create_MC_component_branch(int year, int line)
{
    struct event
    {
        UInt_t runNumber;
        ULong64_t eventNumber;
    };

    TString MC_files = Form("Files_on_grid/MC_201%i.txt",year);
    TFileCollection* fc = new TFileCollection("MC", "MC", MC_files, 1, line);

    TChain* t_reco = new TChain("ntuple/DecayTree");
    TChain* t_gen_3pi_3pi = new TChain("mc_ntuple_3pi_3pi/MCDecayTree");
    TChain* t_gen_3pi_3pipi0 = new TChain("mc_ntuple_3pi_3pipi0/MCDecayTree");
    TChain* t_gen_3pipi0_3pi = new TChain("mc_ntuple_3pipi0_3pi/MCDecayTree");
    TChain* t_gen_3pi_3pi_2pi0 = new TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree");
    t_reco->AddFileInfoList((TCollection*)fc->GetList());
    t_gen_3pi_3pi->AddFileInfoList((TCollection*)fc->GetList());
    t_gen_3pi_3pipi0->AddFileInfoList((TCollection*)fc->GetList());
    t_gen_3pipi0_3pi->AddFileInfoList((TCollection*)fc->GetList());
    t_gen_3pi_3pi_2pi0->AddFileInfoList((TCollection*)fc->GetList());

    t_gen_3pi_3pipi0->GetEntries();
    t_gen_3pipi0_3pi->GetEntries();
    t_gen_3pi_3pipi0->Add(t_gen_3pipi0_3pi);

    // TChain* t_gen_3pi_3pi_pi0 =  new TChain("t_gen_3pi_3pi_pi0", "t_gen_3pi_3pi_pi0");
    // t_gen_3pi_3pi_pi0->Add(t_gen_3pi_3pipi0);
    // t_gen_3pi_3pi_pi0->Add(t_gen_3pipi0_3pi);

    // 3pi 3pi (GEN)
    cout << "Looping over GEN 3pi 3pi" << endl;
    std::vector<event> events_gen_3pi3pi;
    UInt_t runNumber_3pi_3pi;
    ULong64_t eventNumber_3pi_3pi;
    t_gen_3pi_3pi->SetBranchAddress("runNumber", &runNumber_3pi_3pi);
    t_gen_3pi_3pi->SetBranchAddress("eventNumber", &eventNumber_3pi_3pi);

    UInt_t gen_3pi_3pi_entries = t_gen_3pi_3pi->GetEntries();
    for(int i = 0; i < gen_3pi_3pi_entries; i++)
    {
        t_gen_3pi_3pi->GetEntry(i);

        event e;
        e.runNumber = runNumber_3pi_3pi;
        e.eventNumber = eventNumber_3pi_3pi;

        events_gen_3pi3pi.push_back(e);
    }

    // 3pi 3pi pi0 (GEN)
    cout << "Looping over GEN 3pi 3pi pi0" << endl;
    std::vector<event> events_gen_3pi3pipi0;
    UInt_t runNumber_3pi_3pi_pi0;
    ULong64_t eventNumber_3pi_3pi_pi0;
    t_gen_3pi_3pipi0->SetBranchAddress("runNumber", &runNumber_3pi_3pi_pi0);
    t_gen_3pi_3pipi0->SetBranchAddress("eventNumber", &eventNumber_3pi_3pi_pi0);

    UInt_t gen_3pi_3pi_pi0_entries = t_gen_3pi_3pipi0->GetEntries();
    for(int i = 0; i < gen_3pi_3pi_pi0_entries; i++)
    {
        t_gen_3pi_3pipi0->GetEntry(i);

        event e;
        e.runNumber = runNumber_3pi_3pi_pi0;
        e.eventNumber = eventNumber_3pi_3pi_pi0;

        events_gen_3pi3pipi0.push_back(e);
    }

    // 3pi 3pi pi0 (GEN)
    cout << "Looping over GEN 3pi 3pi 2pi0" << endl;
    std::vector<event> events_gen_3pi3pi2pi0;
    UInt_t runNumber_3pi_3pi_2pi0;
    ULong64_t eventNumber_3pi_3pi_2pi0;
    t_gen_3pi_3pi_2pi0->SetBranchAddress("runNumber", &runNumber_3pi_3pi_2pi0);
    t_gen_3pi_3pi_2pi0->SetBranchAddress("eventNumber", &eventNumber_3pi_3pi_2pi0);

    UInt_t gen_3pi_3pi_2pi0_entries = t_gen_3pi_3pi_2pi0->GetEntries();
    for(int i = 0; i < gen_3pi_3pi_2pi0_entries; i++)
    {
        t_gen_3pi_3pi_2pi0->GetEntry(i);

        event e;
        e.runNumber = runNumber_3pi_3pi_2pi0;
        e.eventNumber = eventNumber_3pi_3pi_2pi0;

        events_gen_3pi3pi2pi0.push_back(e);
    }
    
    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/%i.root",year,line), "RECREATE");
    fout->cd();
    TTree * newtree = new TTree("DecayTree", "DecayTree");
    Int_t component = -999;
    newtree->Branch("component", &component);

    cout << "Looping over RECO" << endl;
    UInt_t runNumber;
    ULong64_t eventNumber;
    t_reco->SetBranchAddress("runNumber", &runNumber);
    t_reco->SetBranchAddress("eventNumber", &eventNumber);

    UInt_t reco_entries = t_reco->GetEntries();
    for(int i = 0; i < reco_entries; i++)
    {
        component = -999;

        t_reco->GetEntry(i);
        for(event e: events_gen_3pi3pi)
        {
            if(e.runNumber == runNumber && e.eventNumber == eventNumber)
            {
                component = 0;
                newtree->Fill();
                break;
            }
        }
        if(component == 0) continue;

        for(event e: events_gen_3pi3pipi0)
        {
            if(e.runNumber == runNumber && e.eventNumber == eventNumber)
            {
                component = 1;
                newtree->Fill();
                break;
            }
        }
        if(component == 1) continue;

        for(event e: events_gen_3pi3pi2pi0)
        {
            if(e.runNumber == runNumber && e.eventNumber == eventNumber)
            {
                component = 2;
                newtree->Fill();
                break;
            }
        }
        if(component == 2) continue;

        newtree->Fill();

    }
    newtree->Write();
    fout->Close();

    cout << "Finished successfully" << endl;
}