using namespace std;

void compute_B_yield(int year, TString RS_DATA_files, TString WS_DATA_files)
{
    // files after pre-selections
    TFileCollection* fc_rs_data = new TFileCollection("RS_DATA", "RS_DATA", RS_DATA_files);
    TFileCollection* fc_ws_data = new TFileCollection("WS_DATA", "RS_DATA", WS_DATA_files);

    TChain* t_rs_data = new TChain("DecayTree");
    TChain* t_ws_data = new TChain("DecayTree");

    t_rs_data->AddFileInfoList((TCollection*)fc_rs_data->GetList());
    t_ws_data->AddFileInfoList((TCollection*)fc_ws_data->GetList());

    TCut signal = "(Bp_ConsBp_seq_12_M[0]>4700) && (Bp_ConsBp_seq_12_M[0]<5800)";  
    TCut right = "(Bp_ConsBp_seq_12_M[0]>5800) && (Bp_ConsBp_seq_12_M[0]<8000)";  

    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201%i/bkg_yield.root",year), "RECREATE");
    fout->cd();                                  
    TTree* t_bkg_yield = new TTree("DecayTree", "DecayTree");
    Int_t RS_signal;
    t_bkg_yield->Branch("RS_signal", &RS_signal);

    Int_t RS_right = t_rs_data->GetEntries(right);
    Int_t WS_right = t_ws_data->GetEntries(right);
    Int_t WS_signal = t_ws_data->GetEntries(signal);
    RS_signal = (RS_right/WS_right)*WS_signal;
    t_bkg_yield->Fill();

    t_bkg_yield->Write();
    fout->Close();

    return;
}