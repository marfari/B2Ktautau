using namespace RooStats;
using namespace RooFit;
using namespace std;

TH1D* create_histogram_mc(RooRealVar var, TTree* t, TString weight, Int_t i, Double_t var_min, Double_t var_max);
// void addWeight( TTree* t, TString WEIGHT_file, Int_t year, Int_t species);
void compute_1Dweights(RooWorkspace* w, RooDataSet* data_mc, RooDataSet data_signal_sw, TString var_name, Int_t year, Int_t species);
void compute_2Dweights(RooWorkspace* w, RooDataSet* data_mc, RooDataSet data_signal_sw, TString var1_name, TString var2_name, Int_t year, Int_t species);
Double_t reduced_chi2(TH1D* h_mc, TH1D* h_data);

Int_t N = 50;

void compare_MC_sWeighted_data(Int_t year, Int_t species, Bool_t applyWeight)
{
    TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201%i/Species_8/splot_result.root",year));
    RooWorkspace* w = (RooWorkspace*)f->Get("w_out");
    RooDataSet* data_sw = (RooDataSet*)w->data("dataWithSWeights");
    RooDataSet data_signal_sw{data_sw->GetName(), data_sw->GetTitle(), data_sw, *data_sw->get(), nullptr, "n_signal_sw"};

    TFile* f1 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_7/mass_dataset.root",year));
    RooWorkspace* w1 = (RooWorkspace*)f1->Get("w");
    RooDataSet* data_mc = (RooDataSet*)w1->data("data");

    RooArgSet vars = w1->allVars();
    RooRealVar* mass = w->var("mass");

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

        TFileCollection* fc1 = new TFileCollection("fc1", "fc1", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/201%i/Species_%i/bdt_output.txt",year,species));
        TChain* t1 = new TChain("XGBoost/DecayTree");
        t1->AddFileInfoList((TCollection*)fc1->GetList());

        // Get TTree with weight branch
        TFile* f2 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201%i/DDs_correction/Species_%i/tree_with_weights.root",year,species));
        // TFile* f2 = new TFile(Form("/panfs/felician/B2Ktautau/workflow/addWeight/DDs_correction/201%i/Species_%i/tree_with_weight_1D.root",year,species));
        TTree* tw = (TTree*)f2->Get("DecayTree");


        t->AddFriend(tw);
        t->AddFriend(t1);
        weight = "weights";

        cout << tw->GetEntries() << endl;
        cout << t->GetEntries() << endl;
        cout << t1->GetEntries() << endl;
    }
    
    std::vector<double> red_chi2_values;

    for(int i = 0; i < n_vars; i++)
    {
        RooRealVar* var = (RooRealVar*)vars[i];
        TString name = var->GetName();

        TString ipchi2_vars[] = {"D0bar_K_IPCHI2_OWNPV", "D0bar_pi_IPCHI2_OWNPV", "Dsp_K1_IPCHI2_OWNPV", "Dsp_K2_IPCHI2_OWNPV", "Dsp_pi_IPCHI2_OWNPV"};
        TString pid_vars_1[] = {"D0bar_K_ProbNNk", "Dsp_K1_ProbNNk", "Dsp_K2_ProbNNk", "D0bar_K_ProbNNk_pidgen_default", "Dsp_K1_ProbNNk_pidgen_default", "Dsp_K2_ProbNNk_pidgen_default"};
        TString pid_vars_2[] = {"D0bar_pi_ProbNNpi", "Dsp_pi_ProbNNpi", "D0bar_pi_ProbNNpi_pidgen_default", "Dsp_pi_ProbNNpi_pidgen_default"};
        TString pid_vars_3[] = {"D0bar_K_ProbNNpi", "Dsp_K1_ProbNNpi", "Dsp_K2_ProbNNpi", "D0bar_K_ProbNNpi_pidgen_default", "Dsp_K1_ProbNNpi_pidgen_default", "Dsp_K2_ProbNNpi_pidgen_default"};
        TString pid_vars_4[] = {"D0bar_pi_ProbNNk", "Dsp_pi_ProbNNk", "D0bar_pi_ProbNNk_pidgen_default", "Dsp_pi_ProbNNk_pidgen_default"};

        Double_t var_min, var_max;
        if(std::find(std::begin(ipchi2_vars), std::end(ipchi2_vars), name) != std::end(ipchi2_vars))
        {
            var_min = 0;
            var_max = 15000;
        }
        else if(std::find(std::begin(pid_vars_1), std::end(pid_vars_1), name) != std::end(pid_vars_1))
        {
            var_min = 0.98;
            var_max = 1;
        }
        else if(std::find(std::begin(pid_vars_2), std::end(pid_vars_2), name) != std::end(pid_vars_2))
        {
            var_min = 0.94;
            var_max = 1.;
        }
        else if(std::find(std::begin(pid_vars_3), std::end(pid_vars_3), name) != std::end(pid_vars_3))
        {
            var_min = 0.;
            var_max = 0.04;
        }
        else if(std::find(std::begin(pid_vars_4), std::end(pid_vars_4), name) != std::end(pid_vars_4))
        {
            var_min = 0.;
            var_max = 0.002;
        }
        else if(name == "D0bar_K_IPCHI2_OWNPV")
        {
            var_max = 20000;
        }
        else if(name == "Dsp_pi_IPCHI2_OWNPV")
        {
            var_max = 10000;
        }
        else if(name == "Bp_VTXISODCHI2ONETRACK_B")
        {
            var_max = 3000;
        }        
        else if(name == "Bp_VTXISODCHI2MASSONETRACK_B")
        {
            var_max = 20000;
        }
        else if(name == "Bp_VTXISODCHI2TWOTRACK_B")
        {
            var_max = 3000;
        }
        else if(name == "Bp_VTXISODCHI2MASSTWOTRACK_B")
        {
            var_max = 25000;
        }
        else if(name == "Bp_VTXISODCHI2ONETRACK_D0bar")
        {
            var_max = 1000;
        }
        else if(name == "Bp_VTXISODCHI2MASSONETRACK_D0bar")
        {
            var_max = 6000;
        }
        else if(name == "Bp_VTXISODCHI2TWOTRACK_D0bar")
        {
            var_min = 0;
            var_max = 1000;
        }
        else if(name == "Bp_VTXISODCHI2MASSTWOTRACK_D0bar")
        {
            var_max = 6000;
        }
        else if(name == "Bp_VTXISODCHI2ONETRACK_Dsp")
        {
            var_max = 1000;
        }
        else if(name == "Bp_VTXISODCHI2MASSONETRACK_Dsp")
        {
            var_max = 6500;
        }
        else if(name == "Bp_VTXISODCHI2TWOTRACK_Dsp")
        {
            var_min = 0;
            var_max = 2000;
        }
        else if(name == "Bp_VTXISODCHI2MASSTWOTRACK_Dsp")
        {
            var_max = 6500;
        }
        else if((name == "Bp_CC_05_SPT_B") || (name == "Bp_CC_05_VPT_B") || (name == "Bp_NC_05_SPT_B") || (name == "Bp_NC_05_VPT_B") || (name == "Bp_CC_05_SPT_D0bar") || (name == "Bp_CC_05_VPT_D0bar") || (name == "Bp_NC_05_SPT_D0bar") || (name == "Bp_NC_05_VPT_D0bar") || (name == "Bp_CC_05_SPT_Dsp") || (name == "Bp_CC_05_VPT_Dsp") || (name == "Bp_NC_05_SPT_Dsp") || (name == "Bp_NC_05_VPT_Dsp"))
        {
            var_max = 10000;
        }
        else if((name == "Bp_CC_05_PX_B") || (name == "Bp_CC_05_PY_B") || (name == "Bp_CC_05_PX_D0bar") || (name == "Bp_CC_05_PY_D0bar") || (name == "Bp_CC_05_PX_Dsp") || (name == "Bp_CC_05_PY_Dsp"))
        {
            var_min = -10000;
            var_max = 10000;
        }
        else if((name == "Bp_CC_05_PZ_B") || (name == "Bp_CC_05_PZ_D0bar") || (name == "Bp_CC_05_PZ_Dsp"))
        {
            var_max = 100000;
        }
        else if((name == "Bp_CC_05_PXASYM_B") || (name == "Bp_CC_05_PYASYM_B") || (name == "Bp_NC_05_PXASYM_B") || (name == "Bp_NC_05_PYASYM_B") || (name == "Bp_CC_05_PXASYM_D0bar") || (name == "Bp_CC_05_PYASYM_D0bar") || (name == "Bp_NC_05_PXASYM_D0bar") || (name == "Bp_NC_05_PYASYM_D0bar") || (name == "Bp_CC_05_PXASYM_Dsp") || (name == "Bp_CC_05_PYASYM_Dsp") || (name == "Bp_NC_05_PXASYM_Dsp") || (name == "Bp_NC_05_PYASYM_Dsp"))
        {
            var_min = -1;
        }
        else if((name == "Bp_CC_05_MAXPT_PX_B") || (name == "Bp_CC_05_MAXPT_PY_B") || (name == "Bp_NC_05_PX_B") || (name == "Bp_NC_05_PY_B") || (name == "Bp_NC_05_MAXPT_PX_B") || (name == "Bp_NC_05_MAXPT_PY_B") || (name == "Bp_CC_05_MAXPT_PX_D0bar") || (name == "Bp_CC_05_MAXPT_PY_D0bar") || (name == "Bp_NC_05_PX_D0bar") || (name == "Bp_NC_05_PY_D0bar") || (name == "Bp_NC_05_MAXPT_PX_D0bar") || (name == "Bp_NC_05_MAXPT_PY_D0bar") || (name == "Bp_CC_05_MAXPT_PX_Dsp") || (name == "Bp_CC_05_MAXPT_PY_Dsp") || (name == "Bp_NC_05_PX_Dsp") || (name == "Bp_NC_05_PY_Dsp") || (name == "Bp_NC_05_MAXPT_PX_Dsp") || (name == "Bp_NC_05_MAXPT_PY_Dsp"))
        {
            var_min = -5000;
            var_max = 5000;
        }
        else if((name == "Bp_CC_05_MAXPT_PZ_B") || (name == "Bp_CC_05_MAXPT_PE_B") || (name == "Bp_NC_05_PZ_B") || (name == "Bp_NC_05_MAXPT_PZ_B") || (name == "Bp_NC_05_MAXPT_PE_B") || (name == "Bp_CC_05_MAXPT_PZ_D0bar") || (name == "Bp_CC_05_MAXPT_PE_D0bar") || (name == "Bp_NC_05_PZ_D0bar") || (name == "Bp_NC_05_MAXPT_PZ_D0bar") || (name == "Bp_NC_05_MAXPT_PE_D0bar") || (name == "Bp_CC_05_MAXPT_PZ_Dsp") || (name == "Bp_CC_05_MAXPT_PE_Dsp") || (name == "Bp_NC_05_PZ_Dsp") || (name == "Bp_NC_05_MAXPT_PZ_Dsp") || (name == "Bp_NC_05_MAXPT_PE_Dsp"))
        {
            var_max = 50000;
        }
        else if((name == "Bp_NC_05_MULT_B") || (name == "Bp_NC_05_MULT_D0bar") || (name == "Bp_NC_05_MULT_Dsp"))
        {
            var_min = 0;
        }
        else if((name == "Bp_NC_05_PASYM_B") || (name == "Bp_NC_05_PTASYM_B") || (name == "Bp_NC_05_PZASYM_B") || (name == "Bp_NC_05_PASYM_D0bar") || (name == "Bp_NC_05_PTASYM_D0bar") || (name == "Bp_NC_05_PZASYM_D0bar") || (name == "Bp_NC_05_PASYM_Dsp") || (name == "Bp_NC_05_PTASYM_Dsp") || (name == "Bp_NC_05_PZASYM_Dsp"))
        {
            var_min = -1;
            var_max = 1;
        }
        else if((name == "Bp_NC_05_DELTAETA_B") || (name == "Bp_NC_05_DELTAETA_D0bar") || (name == "Bp_NC_05_DELTAETA_Dsp"))
        {
            var_min = -1;
            var_max = 5;
        }
        else if((name == "Bp_NC_05_DELTAPHI_B") || (name == "Bp_NC_05_DELTAPHI_D0bar") || (name == "Bp_NC_05_DELTAPHI_Dsp"))
        {
            var_min = 0;
            var_max = 3.14;
        }
        else if((name == "Bp_NC_05_IT_B") || (name == "Bp_CCNC_05_IT_B") || (name == "Bp_NC_05_IT_D0bar") || (name == "Bp_CCNC_05_IT_D0bar") || (name == "Bp_NC_05_IT_Dsp") || (name == "Bp_CCNC_05_IT_Dsp"))
        {
            var_min = 0.1;
            var_max = 1;
        }
        else if((name == "Bp_NC_05_PTASYM_B") || (name == "Bp_NC_05_PTASYM_D0nar") || (name == "Bp_NC_05_PTASYM_Dsp"))
        {
            var_min = -0.4;
        }
        else if((name == "Bp_Bstautau_ISOBDTSECONDVALUE_taup") || (name == "Bp_Bstautau_ISOBDTSECONDVALUE_taum"))
        {
            var_min = -0.5;
            var_max = 0;
        }
        else if((name == "Bp_Bstautau_ISOBDTTHIRDVALUE_taup") || (name == "Bp_Bstautau_ISOBDTTHIRDVALUE_taum"))
        {
            var_min = -0.5;
            var_max = 0.5;
        }
        else if(name == "BDT2")
        {
            var_min = 0.8;
        }
        else
        {
            var_min = var->getMin();
            var_max = var->getMax();
        }

        TH1D* h_sig_MC;
        if(applyWeight)
        {
            h_sig_MC = create_histogram_mc(*var, t, weight, i, var_min, var_max);
        }
        else
        {
            h_sig_MC = (TH1D*)data_mc->createHistogram(Form("h_mc_%i",i),*var,Binning(N,var_min,var_max));
        }
         
        TH1D* h_sig_data; 

        if( strcmp(var->GetName(), "D0bar_K_ProbNNk_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("D0bar_K_ProbNNk")),Binning(N,var_min,var_max));
        }
        else if( strcmp(var->GetName(), "D0bar_pi_ProbNNk_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("D0bar_pi_ProbNNk")),Binning(N,var_min,var_max));
        }
        else if( strcmp( var->GetName(), "Dsp_K1_ProbNNk_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("Dsp_K1_ProbNNk")),Binning(N,var_min,var_max));
        }
        else if( strcmp(var->GetName(), "Dsp_K2_ProbNNk_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("Dsp_K2_ProbNNk")),Binning(N,var_min,var_max));
        }
        else if( strcmp( var->GetName(), "Dsp_pi_ProbNNk_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("Dsp_pi_ProbNNk")),Binning(N,var_min,var_max));
        } 
        else if( strcmp(var->GetName(), "D0bar_K_ProbNNpi_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("D0bar_K_ProbNNpi")),Binning(N,var_min,var_max));
        }
        else if( strcmp(var->GetName(), "D0bar_pi_ProbNNpi_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("D0bar_pi_ProbNNpi")),Binning(N,var_min,var_max));
        }
        else if( strcmp(var->GetName(), "Dsp_K1_ProbNNpi_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("Dsp_K1_ProbNNpi")),Binning(N,var_min,var_max));
        }
        else if( strcmp(var->GetName(), "Dsp_K2_ProbNNpi_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("Dsp_K2_ProbNNpi")),Binning(N,var_min,var_max));
        }
        else if( strcmp(var->GetName(), "Dsp_pi_ProbNNpi_pidgen_default") == 0 )
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*(w->var("Dsp_pi_ProbNNpi")),Binning(N,var_min,var_max));
        }
        else
        {
            h_sig_data = (TH1D*)data_signal_sw.createHistogram(Form("h_data_%i",i),*var,Binning(N,var_min,var_max));
        }

        Double_t chi2, red_chi2;
        Int_t ndf, igood;
        if(applyWeight)
        {
            h_sig_MC->Chi2TestX(h_sig_data, chi2, ndf, igood, "WW");
        }
        else
        {
            h_sig_MC->Chi2TestX(h_sig_data, chi2, ndf, igood, "UW");
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
        auto rp = new TRatioPlot(h_sig_data, h_sig_MC, "divsym");
        c.SetTicks(0, 1);
        rp->SetH1DrawOpt("E");
        rp->Draw();
        rp->GetLowerRefYaxis()->SetTitle("Data (sp) / MC");
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

    // 1D weight histogram
    if(!applyWeight)
    {
        compute_1Dweights(w, data_mc, data_signal_sw, "nSPDHits", year, species);
    }

    // 2D Weight map histogram
    // if(!applyWeight)
    // {
    //     compute_2Dweights(w, data_mc, data_signal_sw, "nTracks", "nSPDHits", year, species);
    // }
}

TH1D* create_histogram_mc(RooRealVar var, TTree* t, TString weight, Int_t i, Double_t var_min, Double_t var_max)
{    
    TString name = var.GetName();

    if(name == "mass")
    {
        name = "Bp_dtf_M[0]";
    }
    
    TH1D* h = new TH1D(Form("h_MC_%i",i), Form("h_MC_%i",i), N, var_min, var_max);
    t->Draw( name+Form(" >> h_MC_%i",i) , weight+"*((Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355))");


    return h;
}

// void addWeight( TTree* t, TString WEIGHT_file, Int_t year, Int_t species)
// {   
//     TFile* f = new TFile(WEIGHT_file, "READ");
//     TH2D* h_wei = (TH2D*)f->Get("weight");

//     Int_t num_entries  = t->GetEntries();

//     Int_t n_bins_x = h_wei->GetNbinsX();
//     Int_t n_bins_y = h_wei->GetNbinsY();

//     TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/tree_with_weights.root",year,species), "RECREATE");
//     TTree* tw = new TTree("tw", "tw");
//     Double_t weight_value = 0.; 
//     tw->Branch("weight", &weight_value);

//     Double_t varX, varY;
//     t->SetBranchAddress("Bp_ETA", &varX);
//     t->SetBranchAddress("nSPDHits", &varY);

//     for(int i = 0; i < num_entries; i++)
//     {
//         t->GetEntry(i);

//         weight_value = h_wei->GetBinContent( h_wei->FindBin(varX,varY) );
//         tw->Fill();
//     }

//     fout->cd();
//     tw->Write();
//     fout->Close();
// }

void compute_1Dweights(RooWorkspace* w, RooDataSet* data_mc, RooDataSet data_signal_sw, TString var_name, Int_t year, Int_t species)
{
    RooRealVar* var =  (RooRealVar*)w->var(var_name);

    Double_t var_min = var->getMin();
    Double_t var_max = var->getMax();

    // Int_t n_bins = 9;
    // // Double_t bins[5] = {0, 100, 150, 200, 500};
    // Double_t bins[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450};

    // TH1D* h_sig_data = new TH1D("h_sig_data", "h_sig_data", n_bins, bins);
    // TH1D* h_sig_MC = new TH1D("h_sig_MC", "h_sig_MC", n_bins, bins);
    TH1D* h_wei;

    // RooBinning x_bins(n_bins, bins, "x_bins");

    // h_sig_MC = (TH1D*)data_mc->createHistogram("h_sig_MC", *var, Binning(x_bins));
    // h_sig_data = (TH1D*)data_signal_sw.createHistogram("h_sig_data", *var, Binning(x_bins));
    TH1D* h_sig_MC = (TH1D*)data_mc->createHistogram("h_sig_MC", *var, Binning(N,var_min,var_max));
    TH1D* h_sig_data = (TH1D*)data_signal_sw.createHistogram("h_sig_data", *var, Binning(N,var_min,var_max));

    h_sig_data->Scale(1./h_sig_data->Integral());
    h_sig_MC->Scale(1./h_sig_MC->Integral());

    TCanvas c1;
    c1.cd();
    h_sig_MC->GetXaxis()->SetTitle(var_name);
    h_sig_MC->GetYaxis()->SetTitle(Form("Entries / (%i bins)",N));
    h_sig_MC->SetTitle("MC");
    h_sig_MC->Draw("COLZ");
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/1D_MC.pdf",year,species));

    TCanvas c2;
    c2.cd();
    h_sig_data->GetXaxis()->SetTitle(var_name);
    h_sig_data->GetYaxis()->SetTitle(Form("Entries / (%i bins)",N));
    h_sig_data->SetTitle("Data (sp)");
    h_sig_data->Draw("COLZ");
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/1D_data.pdf",year,species));

    h_wei->SetDefaultSumw2(kTRUE);
    h_wei = (TH1D*)h_sig_data->Clone("weight");
    h_wei->Divide(h_sig_MC);

    TCanvas c3;
    c3.cd();
    h_wei->GetXaxis()->SetTitle(var_name);
    h_wei->GetYaxis()->SetTitle(Form("Entries / (%i bins)",N));
    h_wei->SetTitle("Data (sp) / MC");
    h_wei->Draw("COLZ");
    c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/1D_weights.pdf",year,species));

    TFile* f_wei = new TFile(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/weights.root",year,species), "RECREATE");
    f_wei->cd();
    h_wei->Write();
    f_wei->Close();
}

void compute_2Dweights(RooWorkspace* w, RooDataSet* data_mc, RooDataSet data_signal_sw, TString var1_name, TString var2_name, Int_t year, Int_t species)
{
    RooRealVar* var1 =  (RooRealVar*)w->var(var1_name);
    RooRealVar* var2 =  (RooRealVar*)w->var(var2_name);

    Double_t var1_min = var1->getMin();
    Double_t var1_max = var1->getMax();
    Double_t var2_min = var2->getMin();
    Double_t var2_max = var2->getMax();

    // Int_t N = 10;

    Double_t x_bins[] = {0,50,100,150,200,400};
    Double_t y_bins[] = {0,100,200,300,500};

    TH2D* h_2D_sig_data = new TH2D("h_2D_sig_data", "h_2D_sig_data", N, x_bins, N, y_bins);
    TH2D* h_2D_sig_MC = new TH2D("h_2D_sig_MC", "h_2D_sig_MC", N, x_bins, N, y_bins);
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
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/2D_MC.pdf",year,species));

    TCanvas c2;
    c2.cd();
    h_2D_sig_data->GetXaxis()->SetTitle(var1_name);
    h_2D_sig_data->GetYaxis()->SetTitle(var2_name);
    h_2D_sig_data->SetTitle("Data (sp)");
    h_2D_sig_data->Draw("COLZ");
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/2D_data.pdf",year,species));

    h_2D_wei->SetDefaultSumw2(kTRUE);
    h_2D_wei = (TH2D*)h_2D_sig_data->Clone("weight");
    h_2D_wei->Divide(h_2D_sig_MC);

    TCanvas c3;
    c3.cd();
    h_2D_wei->GetXaxis()->SetTitle(var1_name);
    h_2D_wei->GetYaxis()->SetTitle(var2_name);
    h_2D_wei->SetTitle("Data (sp) / MC");
    h_2D_wei->Draw("COLZ");
    c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/MC_unweighted/2D_weights.pdf",year,species));

    TFile* f_wei = new TFile(Form("/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201%i/DDs_correction/Species_%i/weights.root",year,species), "RECREATE");
    f_wei->cd();
    h_2D_wei->Write();
    f_wei->Close();
}