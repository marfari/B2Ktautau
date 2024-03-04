using namespace RooFit;
using namespace RooStats;

void createDataset(int year, int species, TString FILES)
{
    TFileCollection *fc = new TFileCollection("fc", "fc", FILES);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    RooRealVar mass("mass", "mass", 0, 100000); // put the range large enough to avoid cutting events in the mass distribution
    RooDataSet* data = new RooDataSet("data", "data", mass);

    //float DTF_mass;
    Float_t DTF_mass[100];
    Float_t DTF_status[100];
    if((species == 4) || (species == 5) || (species == 6))
    {
        t->SetBranchAddress("Bp_dtf_DDK_M", DTF_mass);
        t->SetBranchAddress("Bp_dtf_DDK_status", DTF_status);
    }

    int counter = 0;
    cout << "Looping over data" << endl;
    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);

        if((species == 4) || (species == 5) || (species == 6))
        {
            if( (DTF_status[0] == 0) && (DTF_mass[0] > 5180) && (DTF_mass[0] < 5600) )
            {
                mass.setVal(DTF_mass[0]);
                data->add(mass);
            }
        }
        
        if ( i > 1.0*counter*(t->GetEntries())/100 ) {
            cout<<counter<<"%"<<endl;
            counter += 10;
        }
    }
    cout << "Finished looping over data" << endl;

    RooWorkspace* w = new RooWorkspace("w", "w");
    w->import(*data);

    TFile* fout = new TFile( Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_%i/mass_dataset.root",year,species), "RECREATE"); 
    fout->cd();
    w->Write();
    fout->Close();

    return;
}

