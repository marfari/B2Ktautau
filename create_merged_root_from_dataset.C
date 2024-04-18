
void create_merged_root_from_dataset(Int_t species, Bool_t getSignal)
{
    RooDataSet* data;
    if((species == 5) && getSignal) // only signal component of species 5
    {
        TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_%i/splot_result.root",species));
        RooWorkspace* w = (RooWorkspace*)f->Get("w_out");
        data = (RooDataSet*)w->data("dataWithSWeights");
        // data = new RooDataSet{data_sw->GetName(), data_sw->GetTitle(), data_sw, *data_sw->get(), nullptr, "n_signal_sw"};
        data->Print();
    }
    else
    {
        TFile* f1 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/2016/Species_%i/mass_dataset.root",species));
        TFile* f2 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/2017/Species_%i/mass_dataset.root",species));
        TFile* f3 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/2018/Species_%i/mass_dataset.root",species));

        RooWorkspace* w1 = (RooWorkspace*)f1->Get("w");
        RooWorkspace* w2 = (RooWorkspace*)f2->Get("w");
        RooWorkspace* w3 = (RooWorkspace*)f3->Get("w");

        data = (RooDataSet*)w1->data("data");    
        RooDataSet* data2 = (RooDataSet*)w2->data("data");    
        RooDataSet* data3 = (RooDataSet*)w3->data("data");    
        data->append(*data2);
        data->append(*data3);

    }
    TFile* fout;
    if(getSignal)
    {
        fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/species_%i_signal_tree.root",species), "RECREATE");
    } 
    else
    {
        fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/species_%i_tree.root",species), "RECREATE");
    }
    TTree* t = (TTree*)data->GetClonedTree();
    t->SetName("DatasetTree");
    
    fout->cd();
    t->Write();
    fout->Close();
}