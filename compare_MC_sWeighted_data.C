using namespace RooStats;
using namespace RooFit;
using namespace std;

TH1D* create_histogram_mc(RooRealVar var, TTree* t, TString weight, Int_t i, Double_t var_min, Double_t var_max);
void addWeight( TTree* t, TString WEIGHT_file, Int_t year, Int_t species);
void compute_2Dweights(RooWorkspace* w, RooDataSet* data_mc, RooDataSet data_signal_sw, TString var1_name, TString var2_name, Int_t year, Int_t species);
Double_t reduced_chi2(TH1D* h_mc, TH1D* h_data);

Int_t N = 20;

void compare_MC_sWeighted_data(Int_t year, Int_t species, Bool_t applyWeight)
{
    TString MC_files, SPLOT_files;
    if(species == 7)
    {
        MC_files = Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_7/mass_dataset.root",year);
        SPLOT_files = Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201%i/Species_8/splot_result.root",year);
    }

    TFile* f = new TFile(SPLOT_files);
    RooWorkspace* w = (RooWorkspace*)f->Get("w_out");
    RooDataSet* data_sw = (RooDataSet*)w->data("dataWithSWeights");
    RooDataSet data_signal_sw{data_sw->GetName(), data_sw->GetTitle(), data_sw, *data_sw->get(), nullptr, "n_signal_sw"};

    TFile* f1 = new TFile(MC_files);
    RooWorkspace* w1 = (RooWorkspace*)f1->Get("w");
    RooDataSet* data_mc = (RooDataSet*)w1->data("data");

    RooArgSet vars = w->allVars();

    RooRealVar* mass = w->var("mass");

    // DDs 
    RooRealVar* mu = w->var("mu");
    RooRealVar* sigma = w->var("sig");
    RooRealVar* alpha = w->var("alpha");
    RooRealVar* n = w->var("n");
    RooRealVar* frac = w->var("frac");
    RooRealVar* sigma1 = w->var("sig1");
    RooRealVar* alpha1 = w->var("alpha1");
    RooRealVar* lambda = w->var("lambda");
    RooRealVar* n_signal = w->var("n_signal");
    RooRealVar* n_combinatorial = w->var("n_combinatorial");
    RooRealVar* L_n_combinatorial = w->var("L_n_combinatorial");
    RooRealVar* L_n_signal = w->var("L_n_signal");
    RooRealVar* n_combinatorial_sw = w->var("n_combinatorial_sw");
    RooRealVar* n_signal_sw = w->var("n_signal_sw");
    vars.remove(*mu);
    vars.remove(*sigma);
    vars.remove(*alpha);
    vars.remove(*n);
    vars.remove(*frac);
    vars.remove(*sigma1);
    vars.remove(*alpha1);
    vars.remove(*lambda);
    vars.remove(*n_signal);
    vars.remove(*n_combinatorial);
    vars.remove(*L_n_signal);
    vars.remove(*L_n_combinatorial);
    vars.remove(*n_combinatorial_sw);
    vars.remove(*n_signal_sw);

    Int_t n_vars = vars.getSize();

    gStyle->SetOptStat(0);

    TString weight = "";
    TChain* t;
    if(applyWeight)
    {
        // Get TTree
        TFileCollection* fc = new TFileCollection("fc", "fc", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/pre_sel_tree.txt",year,species));
        t = new TChain("DecayTree");
        t->AddFileInfoList((TCollection*)fc->GetList());

        // Create TTree with weight branch
        TFile* f2 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201%i/DDs_correction/Species_%i/tree_with_weights.root",year,species));
        TTree* tw = (TTree*)f2->Get("DecayTree");

        t->AddFriend(tw);
        weight = "weights";

        cout << tw->GetEntries() << endl;
        cout << t->GetEntries() << endl;
    }

    std::vector<double> red_chi2_values;

    for(int i = 0; i < n_vars; i++)
    {
        RooRealVar* var = (RooRealVar*)vars[i];
        TString name = var->GetName();

        TString iso_vars[] = {"Bp_BSTAUTAUTAUISOBDTFIRSTVALUETAUP", "Bp_BSTAUTAUTAUISOBDTSECONDVALUETAUP", "Bp_BSTAUTAUTAUISOBDTTHIRDVALUETAUP", "Bp_BSTAUTAUTAUISOBDTFIRSTVALUETAUM", "Bp_BSTAUTAUTAUISOBDTSECONDVALUETAUM", "Bp_BSTAUTAUTAUISOBDTTHIRDVALUETAUM",
                    "Bp_BSTAUTAUTAUISOFIRSTVALUETAUP", "Bp_BSTAUTAUTAUISOSECONDVALUETAUP", "Bp_BSTAUTAUTAUISOFIRSTVALUETAUM", "Bp_BSTAUTAUTAUISOSECONDVALUETAUM",
                    "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIM", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUPPIM", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIM", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIP", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUMPIP", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIP", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM2", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM2", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM2",
                    "Bp_BSTAUTAUTRKISOFIRSTVALUETAUPPIM", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUMPIP", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUMPIM2"};
     

        Double_t var_min, var_max;
        if(std::find(std::begin(iso_vars), std::end(iso_vars), name) != std::end(iso_vars))
        {
            var_min = -2;
            var_max = 2;
        }   
        else
        {
            var_min = var->getMin();
            var_max = var->getMax();
        }

        TH1D* h_sig_MC;
        TH1D* h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*var,Binning(N,var_min,var_max));
        if(applyWeight)
        {
            h_sig_MC = create_histogram_mc(*var, t, weight, i, var_min, var_max);
        }
        else
        {
            h_sig_MC = (TH1D*)data_mc->createHistogram(Form("h_mc_%i",i),*var,Binning(N,var_min,var_max));
        }

        Double_t chi2, red_chi2;
        Int_t ndf, igood;
        if(applyWeight)
        {
            h_sig_MC->Chi2TestX(h_sig_data, chi2, ndf, igood, "UW");
        }
        else
        {
            h_sig_MC->Chi2TestX(h_sig_data, chi2, ndf, igood, "WW");
        }

        if(ndf != 0)
        {
            red_chi2 = chi2/ndf;
        }
        else
        {
            red_chi2 = 0;
        }
        red_chi2_values.push_back(red_chi2);

        TCanvas c;
        c.cd();

        h_sig_MC->GetXaxis()->SetTitle(var->GetName());
        h_sig_data->GetXaxis()->SetTitle(var->GetName());

        h_sig_MC->GetYaxis()->SetTitle( TString::Format("Normalised entries / (%g)",(var->getMax() - var->getMin())/N) );
        h_sig_data->GetYaxis()->SetTitle( TString::Format("Normalised entries / (%g)",(var->getMax() - var->getMin())/N) );

        h_sig_MC->SetMarkerColor(kOrange+7);
        h_sig_MC->SetLineColor(kOrange+7);

        h_sig_data->SetMarkerColor(kBlue);
        h_sig_data->SetLineColor(kBlue);

        h_sig_data->Scale(1/h_sig_data->Integral());
        h_sig_MC->Scale(1/h_sig_MC->Integral());

        if((h_sig_data->GetMaximum() > h_sig_MC->GetMaximum()))
        {
            h_sig_MC->GetYaxis()->SetRangeUser(0, 1.3*h_sig_data->GetMaximum());
            h_sig_data->GetYaxis()->SetRangeUser(0, 1.3*h_sig_data->GetMaximum());
        }
        else if((h_sig_MC->GetMaximum() > h_sig_data->GetMaximum()))
        {
            h_sig_MC->GetYaxis()->SetRangeUser(0, 1.3*h_sig_MC->GetMaximum());
            h_sig_data->GetYaxis()->SetRangeUser(0, 1.3*h_sig_MC->GetMaximum());
        }
        
        h_sig_MC->SetTitle(Form("#chi^{2} = %.3lf",red_chi2));
        h_sig_data->SetTitle("");
        h_sig_MC->Draw();
        h_sig_data->Draw("same");

        //--TRATIO--//
        auto rp = new TRatioPlot(h_sig_data, h_sig_MC, "diffsig");
        c.SetTicks(0, 1);
        rp->SetH1DrawOpt("E");
        rp->Draw();
        rp->GetLowerRefYaxis()->SetTitle("(Data - MC) / #sigma");
        rp->GetUpperRefYaxis()->SetTitle( TString::Format("Normalised entries / (%g)",(var->getMax() - var->getMin())/N) );
        c.Update();

        TLegend* leg;
        leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg->AddEntry(h_sig_data, "Signal in data", "lfp");
        leg->AddEntry(h_sig_MC, "Signal in MC", "lfp");
        leg->SetTextSize(0.03);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw("same");

        if(applyWeight)
        {
            c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_reweighted/var_%i.pdf",year,species,i));            
        }
        else
        {
            c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/var_%i.pdf",year,species,i));  
        }
    }

    std::vector<int> idx(red_chi2_values.size());
    iota(idx.begin(), idx.end(), 0);
    stable_sort(idx.begin(), idx.end(), [&red_chi2_values](int i1, int i2) {return red_chi2_values[i1] < red_chi2_values[i2];});

    // std::sort(red_chi2_values.begin(),red_chi2_values.end()); //Sorting the vector using greater<int>() function

    for(int i = 0; i < n_vars; i++)
    {
        RooRealVar* var = (RooRealVar*)vars[idx[i]];
        TString name = var->GetName();

        cout << name <<  " chi2/ndf = " << red_chi2_values[idx[i]] << endl;
    }

    // 2D Weight map histogram
    // if(!applyWeight)
    // {
    //     if(species == 5)
    //     {
    //         compute_2Dweights(w, data_mc, data_signal_sw, "nTracks", "nSPDHits", year, species);
    //     }
    //     if(species == 8)
    //     {
    //         compute_2Dweights(w, data_mc, data_signal_sw, "Bp_PT", "nSPDHits", year, species);
    //     }
    // }
}

TH1D* create_histogram_mc(RooRealVar var, TTree* t, TString weight, Int_t i, Double_t var_min, Double_t var_max)
{    
    TString name = var.GetName();

    if(name == "mass")
    {
        name = "Bp_dtf_M";
    }
    
    TH1D* h = new TH1D(Form("h_MC_%i",i), Form("h_MC_%i",i), N, var_min, var_max);
    t->Draw( name+Form(" >> h_MC_%i",i) , weight);

    return h;
}

void addWeight( TTree* t, TString WEIGHT_file, Int_t year, Int_t species)
{   
    TFile* f = new TFile(WEIGHT_file, "READ");
    TH2D* h_wei = (TH2D*)f->Get("weight");

    Int_t num_entries  = t->GetEntries();

    Int_t n_bins_x = h_wei->GetNbinsX();
    Int_t n_bins_y = h_wei->GetNbinsY();

    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/tree_with_weights.root",year,species), "RECREATE");
    TTree* tw = new TTree("tw", "tw");
    Double_t weight_value = 0.; 
    tw->Branch("weight", &weight_value);

    Double_t varX, varY;
    t->SetBranchAddress("Bp_ETA", &varX);
    t->SetBranchAddress("nSPDHits", &varY);

    for(int i = 0; i < num_entries; i++)
    {
        t->GetEntry(i);

        weight_value = h_wei->GetBinContent( h_wei->FindBin(varX,varY) );
        tw->Fill();
    }

    fout->cd();
    tw->Write();
    fout->Close();
}

void compute_2Dweights(RooWorkspace* w, RooDataSet* data_mc, RooDataSet data_signal_sw, TString var1_name, TString var2_name, Int_t year, Int_t species)
{
    RooRealVar* var1 =  (RooRealVar*)w->var(var1_name);
    RooRealVar* var2 =  (RooRealVar*)w->var(var2_name);

    Double_t var1_min = var1->getMin();
    Double_t var1_max = var1->getMax();
    Double_t var2_min = var2->getMin();
    Double_t var2_max = var2->getMax();

    TH2D* h_2D_sig_data = new TH2D("h_2D_sig_data", "h_2D_sig_data", N, var1_min, var1_max, N, var2_min, var2_max);
    TH2D* h_2D_sig_MC = new TH2D("h_2D_sig_MC", "h_2D_sig_MC", N, var1_min, var1_max, N, var2_min, var2_max);
    TH2D* h_2D_wei;

    h_2D_sig_MC = (TH2D*)data_mc->createHistogram("h_2D_sig_MC", *var1, Binning(N,var1_min,var1_max), YVar(*var2, Binning(N,var2_min,var2_max)));
    h_2D_sig_data = (TH2D*)data_signal_sw.createHistogram("h_2D_sig_data", *var1, Binning(N,var1_min,var1_max), YVar(*var2, Binning(N,var2_min,var2_max)) );

    h_2D_sig_data->Scale(1./h_2D_sig_data->Integral());
    h_2D_sig_MC->Scale(1./h_2D_sig_MC->Integral());

    TCanvas c1;
    c1.cd();
    h_2D_sig_MC->GetXaxis()->SetTitle(var1_name);
    h_2D_sig_MC->GetYaxis()->SetTitle(var2_name);
    h_2D_sig_MC->SetTitle("MC");
    h_2D_sig_MC->Draw("COLZ");
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/MC_unweighted/2D_MC.gif",year,species));
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/MC_unweighted/2D_MC.pdf",year,species));

    TCanvas c2;
    c2.cd();
    h_2D_sig_data->GetXaxis()->SetTitle(var1_name);
    h_2D_sig_data->GetYaxis()->SetTitle(var2_name);
    h_2D_sig_data->SetTitle("Data (sp)");
    h_2D_sig_data->Draw("COLZ");
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/MC_unweighted/2D_data.gif",year,species));
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/MC_unweighted/2D_data.pdf",year,species));

    h_2D_wei->SetDefaultSumw2(kTRUE);
    h_2D_wei = (TH2D*)h_2D_sig_data->Clone("weight");
    h_2D_wei->Divide(h_2D_sig_MC);

    TCanvas c3;
    c3.cd();
    h_2D_wei->GetXaxis()->SetTitle(var1_name);
    h_2D_wei->GetYaxis()->SetTitle(var2_name);
    h_2D_wei->SetTitle("Data (sp) / MC");
    h_2D_wei->Draw("COLZ");
    c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/MC_unweighted/2D_weights.gif",year,species));
    c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/MC_unweighted/2D_weights.pdf",year,species));

    TFile* f_wei = new TFile(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/Species_%i/weights.root",year,species), "RECREATE");
    f_wei->cd();
    h_2D_wei->Write();
    f_wei->Close();
}