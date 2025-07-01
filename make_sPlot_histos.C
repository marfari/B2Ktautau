using namespace RooStats;
using namespace RooFit;
using namespace std;

Int_t nbins = 50;

void make_sPlot_histos(Int_t species)
{
    TFile* f = new TFile("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_8/mass_fit_result.root");
    RooWorkspace* w = (RooWorkspace*)f->Get("w_out");

    // 1) Get RooDataSet + RooAbsPdf + RooRealVars
    RooDataSet* data = (RooDataSet*)w->data("data");
    RooAbsPdf* model = w->pdf("model");
    RooAbsPdf* bkg_pdf = w->pdf("exp");

    RooRealVar* mass = w->var("mass");
    RooArgSet vars = w->allVars();

    RooRealVar *mu, *sigma, *sigma1, *alpha, *alpha1, *n, *frac, *lambda, *n_signal, *n_combinatorial;
    mu = w->var("mu");
    sigma = w->var("sig");
    alpha = w->var("alpha");
    n = w->var("n");
    frac = w->var("frac");
    sigma1 = w->var("sig1");
    alpha1 = w->var("alpha1");
    lambda = w->var("lambda");
    n_signal = w->var("n_signal");
    n_combinatorial = w->var("n_combinatorial");
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
    
    Int_t n_vars = vars.getSize();

    // Sideband subtraction
    // 1) Divide RooDataSet in 2 regions: (signal peak region) + (sideband regions)
    Double_t mass_min = mass->getMin();
    Double_t mass_max = mass->getMax();
    Double_t peak_min = 5235.;
    Double_t peak_max = 5320.;

    RooDataSet* peak_data = (RooDataSet*)data->reduce(Form("(mass > %f) && (mass < %f)",peak_min,peak_max));
    RooDataSet* sideband_data = (RooDataSet*)data->reduce(Form("((mass > %f) && (mass < %f)) || ((mass > %f) && (mass < %f))",mass_min, peak_min, peak_max, mass_max));

    // 2) Compute r by integrating the bkg PDF
    mass->setRange("L", mass_min, peak_min);
    mass->setRange("P", peak_min, peak_max);
    mass->setRange("R", peak_max, mass_max);

    RooAbsReal* L = bkg_pdf->createIntegral(*mass, *mass, "L");
    RooAbsReal* P = bkg_pdf->createIntegral(*mass, *mass, "P");
    RooAbsReal* R = bkg_pdf->createIntegral(*mass, *mass, "R");

    Double_t r = P->getVal()/( L->getVal() + R->getVal() );

    // 3) Create sidbeand subtraction histograms: signal, background and total in the same plot
    std::vector<TH1D*> sideband_histos;
    for(int i = 0; i < n_vars; i++)
    {
        RooRealVar* var = (RooRealVar*)vars[i];
        TString var_name = var->GetName();
        Double_t x_min = var->getMin();
        Double_t x_max = var->getMax(); 

        if(var_name == "D0bar_K_IPCHI2_OWNPV")
        {
            x_max = 20000;
        }
        else if(var_name == "Dsp_pi_IPCHI2_OWNPV")
        {
            x_max = 10000;
        }
        else if(var_name == "Bp_VTXISODCHI2ONETRACK_B")
        {
            x_max = 3000;
        }        
        else if(var_name == "Bp_VTXISODCHI2MASSONETRACK_B")
        {
            x_max = 20000;
        }
        else if(var_name == "Bp_VTXISODCHI2TWOTRACK_B")
        {
            x_max = 3000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSTWOTRACK_B")
        {
            x_max = 25000;
        }
        else if(var_name == "Bp_VTXISODCHI2ONETRACK_D0bar")
        {
            x_max = 1000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSONETRACK_D0bar")
        {
            x_max = 6000;
        }
        else if(var_name == "Bp_VTXISODCHI2TWOTRACK_D0bar")
        {
            x_min = 0;
            x_max = 1000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSTWOTRACK_D0bar")
        {
            x_max = 6000;
        }
        else if(var_name == "Bp_VTXISODCHI2ONETRACK_Dsp")
        {
            x_max = 1000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSONETRACK_Dsp")
        {
            x_max = 6500;
        }
        else if(var_name == "Bp_VTXISODCHI2TWOTRACK_Dsp")
        {
            x_min = 0;
            x_max = 2000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSTWOTRACK_Dsp")
        {
            x_max = 6500;
        }
        else if((var_name == "Bp_CC_05_SPT_B") || (var_name == "Bp_CC_05_VPT_B") || (var_name == "Bp_NC_05_SPT_B") || (var_name == "Bp_NC_05_VPT_B") || (var_name == "Bp_CC_05_SPT_D0bar") || (var_name == "Bp_CC_05_VPT_D0bar") || (var_name == "Bp_NC_05_SPT_D0bar") || (var_name == "Bp_NC_05_VPT_D0bar") || (var_name == "Bp_CC_05_SPT_Dsp") || (var_name == "Bp_CC_05_VPT_Dsp") || (var_name == "Bp_NC_05_SPT_Dsp") || (var_name == "Bp_NC_05_VPT_Dsp"))
        {
            x_max = 10000;
        }
        else if((var_name == "Bp_CC_05_PX_B") || (var_name == "Bp_CC_05_PY_B") || (var_name == "Bp_CC_05_PX_D0bar") || (var_name == "Bp_CC_05_PY_D0bar") || (var_name == "Bp_CC_05_PX_Dsp") || (var_name == "Bp_CC_05_PY_Dsp"))
        {
            x_min = -10000;
            x_max = 10000;
        }
        else if((var_name == "Bp_CC_05_PZ_B") || (var_name == "Bp_CC_05_PZ_D0bar") || (var_name == "Bp_CC_05_PZ_Dsp"))
        {
            x_max = 100000;
        }
        else if((var_name == "Bp_CC_05_PXASYM_B") || (var_name == "Bp_CC_05_PYASYM_B") || (var_name == "Bp_NC_05_PXASYM_B") || (var_name == "Bp_NC_05_PYASYM_B") || (var_name == "Bp_CC_05_PXASYM_D0bar") || (var_name == "Bp_CC_05_PYASYM_D0bar") || (var_name == "Bp_NC_05_PXASYM_D0bar") || (var_name == "Bp_NC_05_PYASYM_D0bar") || (var_name == "Bp_CC_05_PXASYM_Dsp") || (var_name == "Bp_CC_05_PYASYM_Dsp") || (var_name == "Bp_NC_05_PXASYM_Dsp") || (var_name == "Bp_NC_05_PYASYM_Dsp"))
        {
            x_min = -1;
        }
        else if((var_name == "Bp_CC_05_MAXPT_PX_B") || (var_name == "Bp_CC_05_MAXPT_PY_B") || (var_name == "Bp_NC_05_PX_B") || (var_name == "Bp_NC_05_PY_B") || (var_name == "Bp_NC_05_MAXPT_PX_B") || (var_name == "Bp_NC_05_MAXPT_PY_B") || (var_name == "Bp_CC_05_MAXPT_PX_D0bar") || (var_name == "Bp_CC_05_MAXPT_PY_D0bar") || (var_name == "Bp_NC_05_PX_D0bar") || (var_name == "Bp_NC_05_PY_D0bar") || (var_name == "Bp_NC_05_MAXPT_PX_D0bar") || (var_name == "Bp_NC_05_MAXPT_PY_D0bar") || (var_name == "Bp_CC_05_MAXPT_PX_Dsp") || (var_name == "Bp_CC_05_MAXPT_PY_Dsp") || (var_name == "Bp_NC_05_PX_Dsp") || (var_name == "Bp_NC_05_PY_Dsp") || (var_name == "Bp_NC_05_MAXPT_PX_Dsp") || (var_name == "Bp_NC_05_MAXPT_PY_Dsp"))
        {
            x_min = -5000;
            x_max = 5000;
        }
        else if((var_name == "Bp_CC_05_MAXPT_PZ_B") || (var_name == "Bp_CC_05_MAXPT_PE_B") || (var_name == "Bp_NC_05_PZ_B") || (var_name == "Bp_NC_05_MAXPT_PZ_B") || (var_name == "Bp_NC_05_MAXPT_PE_B") || (var_name == "Bp_CC_05_MAXPT_PZ_D0bar") || (var_name == "Bp_CC_05_MAXPT_PE_D0bar") || (var_name == "Bp_NC_05_PZ_D0bar") || (var_name == "Bp_NC_05_MAXPT_PZ_D0bar") || (var_name == "Bp_NC_05_MAXPT_PE_D0bar") || (var_name == "Bp_CC_05_MAXPT_PZ_Dsp") || (var_name == "Bp_CC_05_MAXPT_PE_Dsp") || (var_name == "Bp_NC_05_PZ_Dsp") || (var_name == "Bp_NC_05_MAXPT_PZ_Dsp") || (var_name == "Bp_NC_05_MAXPT_PE_Dsp"))
        {
            x_max = 50000;
        }
        else if((var_name == "Bp_NC_05_MULT_B") || (var_name == "Bp_NC_05_MULT_D0bar") || (var_name == "Bp_NC_05_MULT_Dsp"))
        {
            x_min = 0;
        }
        else if((var_name == "Bp_NC_05_PASYM_B") || (var_name == "Bp_NC_05_PTASYM_B") || (var_name == "Bp_NC_05_PZASYM_B") || (var_name == "Bp_NC_05_PASYM_D0bar") || (var_name == "Bp_NC_05_PTASYM_D0bar") || (var_name == "Bp_NC_05_PZASYM_D0bar") || (var_name == "Bp_NC_05_PASYM_Dsp") || (var_name == "Bp_NC_05_PTASYM_Dsp") || (var_name == "Bp_NC_05_PZASYM_Dsp"))
        {
            x_min = -1;
            x_max = 1;
        }
        else if((var_name == "Bp_NC_05_DELTAETA_B") || (var_name == "Bp_NC_05_DELTAETA_D0bar") || (var_name == "Bp_NC_05_DELTAETA_Dsp"))
        {
            x_min = -1;
            x_max = 5;
        }
        else if((var_name == "Bp_NC_05_DELTAPHI_B") || (var_name == "Bp_NC_05_DELTAPHI_D0bar") || (var_name == "Bp_NC_05_DELTAPHI_Dsp"))
        {
            x_min = 0;
            x_max = 3.14;
        }
        else if((var_name == "Bp_NC_05_IT_B") || (var_name == "Bp_CCNC_05_IT_B") || (var_name == "Bp_NC_05_IT_D0bar") || (var_name == "Bp_CCNC_05_IT_D0bar") || (var_name == "Bp_NC_05_IT_Dsp") || (var_name == "Bp_CCNC_05_IT_Dsp"))
        {
            x_min = 0.1;
            x_max = 1;
        }
        else if((var_name == "Bp_NC_05_PTASYM_B") || (var_name == "Bp_NC_05_PTASYM_D0nar") || (var_name == "Bp_NC_05_PTASYM_Dsp"))
        {
            x_min = -0.4;
        }
        else if((var_name == "Bp_Bstautau_ISOBDTSECONDVALUE_taup") || (var_name == "Bp_Bstautau_ISOBDTSECONDVALUE_taum"))
        {
            x_min = -0.5;
            x_max = 0;
        }
        else if((var_name == "Bp_Bstautau_ISOBDTTHIRDVALUE_taup") || (var_name == "Bp_Bstautau_ISOBDTTHIRDVALUE_taum"))
        {
            x_min = -0.5;
            x_max = 0.5;
        }
        else if(var_name == "BDT2")
        {
            x_min = 0.8;
        }
        else
        {
            x_min = var->getMin();
            x_max = var->getMax();
        }

        TH1D* h_peak = (TH1D*)peak_data->createHistogram(Form("h_peak_%i",i), *var, Binning(nbins, x_min, x_max));
        TH1D* h_sig = new TH1D(*h_peak); // h_sig = h_peak
        TH1D* h_sideband = (TH1D*)sideband_data->createHistogram(Form("h_sideband_%i",i), *var, Binning(nbins, x_min, x_max));
        h_sig->Add(h_sideband, -r); // h_sig = h_peak - r*h_sideband
        h_sideband->Scale(r);

        gStyle->SetOptStat(0);
        h_peak->SetMarkerColor(kBlack);
        h_peak->SetLineColor(kBlack);
        h_sig->SetMarkerColor(kBlue);
        h_sig->SetLineColor(kBlue);
        h_sideband->SetMarkerColor(kRed);
        h_sideband->SetLineColor(kRed);
        h_peak->SetTitle("");
        h_sideband->SetTitle("");
        h_sig->SetTitle("");

        h_peak->GetYaxis()->SetRangeUser(0, 1.3*h_peak->GetMaximum());

        TCanvas c;
        c.cd();

        h_peak->Draw();
        h_sig->Draw("same");
        h_sideband->Draw("same");

        TLegend *leg = new TLegend(0.7, 0.9, 0.9, 1.0);
        leg->AddEntry(h_peak, "Total", "lp");
        leg->AddEntry(h_sig, "Signal", "lp");
        leg->AddEntry(h_sideband, "Background", "lp");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw("same");
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_%i/Sideband_plots/var_%i.pdf",species,i));

        sideband_histos.push_back(h_sig);
    }

    // sPlot
    // The sPlot technique requires the model parameters (other than the yields) to be fixed
    mu->setConstant();
    sigma->setConstant();
    alpha->setConstant();
    n->setConstant();
    frac->setConstant();
    sigma1->setConstant();
    alpha1->setConstant();
    lambda->setConstant();
    
    RooMsgService::instance().setSilentMode(true);

    // Fit the model to the data
    RooStats::SPlot sData{"sData", "An SPlot", *data, model, RooArgList(*n_signal, *n_combinatorial)};

    std::cout << "\n\nThe dataset after creating sWeights:\n";
    data->Print();

    // Check that our weights have the desired properties
    std::cout << "\n\n------------------------------------------\n\nCheck SWeights:" << std::endl;
    std::cout << std::endl
                << "Yield of signal is\t" << n_signal->getVal() << ".  From sWeights it is "
                << sData.GetYieldFromSWeight("n_signal") << std::endl;

    std::cout << "Yield of background is\t" << n_combinatorial->getVal() << ".  From sWeights it is "
                << sData.GetYieldFromSWeight("n_combinatorial") << std::endl
                << std::endl;

    for (Int_t i = 0; i < 10; i++) {
        std::cout << "signal Weight for event " << i << std::right << std::setw(12) << sData.GetSWeight(i, "n_signal") << "  background Weight"
                << std::setw(12) << sData.GetSWeight(i, "n_combinatorial") << "  Total Weight" << std::setw(12) << sData.GetSumOfEventSWeight(i)
                << std::endl;
    }
    std::cout << std::endl;

    // Import new dataset with sWeights
    std::cout << "import new dataset with sWeights" << std::endl;
    w->import(*data, Rename("dataWithSWeights"));

    RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);

    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_%i/splot_result.root",species), "RECREATE");
    fout->cd();
    w->Write();
    fout->Close();

    // Make splot plots
    std::cout << "make plots" << std::endl;

    // Create weighted datasets from sWeighted dataset
    RooDataSet* data_sw = (RooDataSet*)w->data("dataWithSWeights");
    RooDataSet data_signal_sw{data_sw->GetName(), data_sw->GetTitle(), data_sw, *data_sw->get(), nullptr, "n_signal_sw"};
    RooDataSet data_background_sw{data_sw->GetName(), data_sw->GetTitle(), data_sw, *data_sw->get(), nullptr, "n_combinatorial_sw"};

    // Plot total vs signal vs background distributions in the same plot
    std::vector<TH1D*> splot_histos; 
    for(int i = 0; i < n_vars; i++)
    {
        RooRealVar* var = (RooRealVar*)vars[i];
        TString var_name = var->GetName();

        TCanvas c;
        c.cd();
        Double_t x_min = var->getMin();
        Double_t x_max = var->getMax();

        if(var_name == "D0bar_K_IPCHI2_OWNPV")
        {
            x_max = 20000;
        }
        else if(var_name == "Dsp_pi_IPCHI2_OWNPV")
        {
            x_max = 10000;
        }
        else if(var_name == "Bp_VTXISODCHI2ONETRACK_B")
        {
            x_max = 3000;
        }        
        else if(var_name == "Bp_VTXISODCHI2MASSONETRACK_B")
        {
            x_max = 20000;
        }
        else if(var_name == "Bp_VTXISODCHI2TWOTRACK_B")
        {
            x_max = 3000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSTWOTRACK_B")
        {
            x_max = 25000;
        }
        else if(var_name == "Bp_VTXISODCHI2ONETRACK_D0bar")
        {
            x_max = 1000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSONETRACK_D0bar")
        {
            x_max = 6000;
        }
        else if(var_name == "Bp_VTXISODCHI2TWOTRACK_D0bar")
        {
            x_min = 0;
            x_max = 1000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSTWOTRACK_D0bar")
        {
            x_max = 6000;
        }
        else if(var_name == "Bp_VTXISODCHI2ONETRACK_Dsp")
        {
            x_max = 1000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSONETRACK_Dsp")
        {
            x_max = 6500;
        }
        else if(var_name == "Bp_VTXISODCHI2TWOTRACK_Dsp")
        {
            x_min = 0;
            x_max = 2000;
        }
        else if(var_name == "Bp_VTXISODCHI2MASSTWOTRACK_Dsp")
        {
            x_max = 6500;
        }
        else if((var_name == "Bp_CC_05_SPT_B") || (var_name == "Bp_CC_05_VPT_B") || (var_name == "Bp_NC_05_SPT_B") || (var_name == "Bp_NC_05_VPT_B") || (var_name == "Bp_CC_05_SPT_D0bar") || (var_name == "Bp_CC_05_VPT_D0bar") || (var_name == "Bp_NC_05_SPT_D0bar") || (var_name == "Bp_NC_05_VPT_D0bar") || (var_name == "Bp_CC_05_SPT_Dsp") || (var_name == "Bp_CC_05_VPT_Dsp") || (var_name == "Bp_NC_05_SPT_Dsp") || (var_name == "Bp_NC_05_VPT_Dsp"))
        {
            x_max = 10000;
        }
        else if((var_name == "Bp_CC_05_PX_B") || (var_name == "Bp_CC_05_PY_B") || (var_name == "Bp_CC_05_PX_D0bar") || (var_name == "Bp_CC_05_PY_D0bar") || (var_name == "Bp_CC_05_PX_Dsp") || (var_name == "Bp_CC_05_PY_Dsp"))
        {
            x_min = -10000;
            x_max = 10000;
        }
        else if((var_name == "Bp_CC_05_PZ_B") || (var_name == "Bp_CC_05_PZ_D0bar") || (var_name == "Bp_CC_05_PZ_Dsp"))
        {
            x_max = 100000;
        }
        else if((var_name == "Bp_CC_05_PXASYM_B") || (var_name == "Bp_CC_05_PYASYM_B") || (var_name == "Bp_NC_05_PXASYM_B") || (var_name == "Bp_NC_05_PYASYM_B") || (var_name == "Bp_CC_05_PXASYM_D0bar") || (var_name == "Bp_CC_05_PYASYM_D0bar") || (var_name == "Bp_NC_05_PXASYM_D0bar") || (var_name == "Bp_NC_05_PYASYM_D0bar") || (var_name == "Bp_CC_05_PXASYM_Dsp") || (var_name == "Bp_CC_05_PYASYM_Dsp") || (var_name == "Bp_NC_05_PXASYM_Dsp") || (var_name == "Bp_NC_05_PYASYM_Dsp"))
        {
            x_min = -1;
        }
        else if((var_name == "Bp_CC_05_MAXPT_PX_B") || (var_name == "Bp_CC_05_MAXPT_PY_B") || (var_name == "Bp_NC_05_PX_B") || (var_name == "Bp_NC_05_PY_B") || (var_name == "Bp_NC_05_MAXPT_PX_B") || (var_name == "Bp_NC_05_MAXPT_PY_B") || (var_name == "Bp_CC_05_MAXPT_PX_D0bar") || (var_name == "Bp_CC_05_MAXPT_PY_D0bar") || (var_name == "Bp_NC_05_PX_D0bar") || (var_name == "Bp_NC_05_PY_D0bar") || (var_name == "Bp_NC_05_MAXPT_PX_D0bar") || (var_name == "Bp_NC_05_MAXPT_PY_D0bar") || (var_name == "Bp_CC_05_MAXPT_PX_Dsp") || (var_name == "Bp_CC_05_MAXPT_PY_Dsp") || (var_name == "Bp_NC_05_PX_Dsp") || (var_name == "Bp_NC_05_PY_Dsp") || (var_name == "Bp_NC_05_MAXPT_PX_Dsp") || (var_name == "Bp_NC_05_MAXPT_PY_Dsp"))
        {
            x_min = -5000;
            x_max = 5000;
        }
        else if((var_name == "Bp_CC_05_MAXPT_PZ_B") || (var_name == "Bp_CC_05_MAXPT_PE_B") || (var_name == "Bp_NC_05_PZ_B") || (var_name == "Bp_NC_05_MAXPT_PZ_B") || (var_name == "Bp_NC_05_MAXPT_PE_B") || (var_name == "Bp_CC_05_MAXPT_PZ_D0bar") || (var_name == "Bp_CC_05_MAXPT_PE_D0bar") || (var_name == "Bp_NC_05_PZ_D0bar") || (var_name == "Bp_NC_05_MAXPT_PZ_D0bar") || (var_name == "Bp_NC_05_MAXPT_PE_D0bar") || (var_name == "Bp_CC_05_MAXPT_PZ_Dsp") || (var_name == "Bp_CC_05_MAXPT_PE_Dsp") || (var_name == "Bp_NC_05_PZ_Dsp") || (var_name == "Bp_NC_05_MAXPT_PZ_Dsp") || (var_name == "Bp_NC_05_MAXPT_PE_Dsp"))
        {
            x_max = 50000;
        }
        else if((var_name == "Bp_NC_05_MULT_B") || (var_name == "Bp_NC_05_MULT_D0bar") || (var_name == "Bp_NC_05_MULT_Dsp"))
        {
            x_min = 0;
        }
        else if((var_name == "Bp_NC_05_PASYM_B") || (var_name == "Bp_NC_05_PTASYM_B") || (var_name == "Bp_NC_05_PZASYM_B") || (var_name == "Bp_NC_05_PASYM_D0bar") || (var_name == "Bp_NC_05_PTASYM_D0bar") || (var_name == "Bp_NC_05_PZASYM_D0bar") || (var_name == "Bp_NC_05_PASYM_Dsp") || (var_name == "Bp_NC_05_PTASYM_Dsp") || (var_name == "Bp_NC_05_PZASYM_Dsp"))
        {
            x_min = -1;
            x_max = 1;
        }
        else if((var_name == "Bp_NC_05_DELTAETA_B") || (var_name == "Bp_NC_05_DELTAETA_D0bar") || (var_name == "Bp_NC_05_DELTAETA_Dsp"))
        {
            x_min = -1;
            x_max = 5;
        }
        else if((var_name == "Bp_NC_05_DELTAPHI_B") || (var_name == "Bp_NC_05_DELTAPHI_D0bar") || (var_name == "Bp_NC_05_DELTAPHI_Dsp"))
        {
            x_min = 0;
            x_max = 3.14;
        }
        else if((var_name == "Bp_NC_05_IT_B") || (var_name == "Bp_CCNC_05_IT_B") || (var_name == "Bp_NC_05_IT_D0bar") || (var_name == "Bp_CCNC_05_IT_D0bar") || (var_name == "Bp_NC_05_IT_Dsp") || (var_name == "Bp_CCNC_05_IT_Dsp"))
        {
            x_min = 0.1;
            x_max = 1;
        }
        else if((var_name == "Bp_NC_05_PTASYM_B") || (var_name == "Bp_NC_05_PTASYM_D0nar") || (var_name == "Bp_NC_05_PTASYM_Dsp"))
        {
            x_min = -0.4;
        }
        else if((var_name == "Bp_Bstautau_ISOBDTSECONDVALUE_taup") || (var_name == "Bp_Bstautau_ISOBDTSECONDVALUE_taum"))
        {
            x_min = -0.5;
            x_max = 0;
        }
        else if((var_name == "Bp_Bstautau_ISOBDTTHIRDVALUE_taup") || (var_name == "Bp_Bstautau_ISOBDTTHIRDVALUE_taum"))
        {
            x_min = -0.5;
            x_max = 0.5;
        }
        else if(var_name == "BDT2")
        {
            x_min = 0.8;
        }
        else
        {
            x_min = var->getMin();
            x_max = var->getMax();
        }

        TH1D* h_tot = (TH1D*)data_sw->createHistogram("tot",*var,Binning(nbins,x_min,x_max));
        TH1D* h_sig = (TH1D*)data_signal_sw.createHistogram("sig",*var,Binning(nbins,x_min,x_max));
        TH1D* h_bkg = (TH1D*)data_background_sw.createHistogram("bkg",*var,Binning(nbins,x_min,x_max));

        gStyle->SetOptStat(0);
        h_tot->SetLineColor(kBlack);
        h_sig->SetLineColor(kBlue);
        h_bkg->SetLineColor(kRed);

        h_tot->GetXaxis()->SetTitle(var_name);
        h_tot->GetYaxis()->SetTitle( TString::Format("Events /(%g)",(var->getMax() - var->getMin())/float(nbins)) );
        h_tot->SetTitle("");

        h_tot->Draw("E");
        h_sig->Draw("E same");
        h_bkg->Draw("E same");

        TLegend *leg1 = new TLegend(0.7, 0.9, 0.9, 1.0);
        leg1->AddEntry(h_tot, "Total", "lp");
        leg1->AddEntry(h_sig, "Signal", "lp");
        leg1->AddEntry(h_bkg, "Background", "lp");
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->Draw("same");
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_%i/Splot_plots/var_%i.pdf",species,i));

        splot_histos.push_back(h_sig);
    }

    // Sideband vs splot
    for(int i = 0; i < n_vars; i++)
    {
        RooRealVar* var = (RooRealVar*)vars[i];
        TString var_name = var->GetName();

        sideband_histos[i]->SetMarkerColor(kPink+7);
        sideband_histos[i]->SetLineColor(kPink+7);
        splot_histos[i]->SetMarkerColor(kGreen+2);
        splot_histos[i]->SetLineColor(kGreen+2);

        TCanvas c;
        c.cd();

        sideband_histos[i]->SetTitle("");
        splot_histos[i]->SetTitle("");

        sideband_histos[i]->Scale(1/sideband_histos[i]->Integral());
        splot_histos[i]->Scale(1/splot_histos[i]->Integral());

        if((sideband_histos[i]->GetMaximum() > splot_histos[i]->GetMaximum()))
        {
            sideband_histos[i]->GetYaxis()->SetRangeUser(0, 1.3*sideband_histos[i]->GetMaximum());
            sideband_histos[i]->Draw();
            splot_histos[i]->Draw("same");
        }
        else if((splot_histos[i]->GetMaximum() > sideband_histos[i]->GetMaximum()))
        {
            splot_histos[i]->GetYaxis()->SetRangeUser(0, 1.3*splot_histos[i]->GetMaximum());
            splot_histos[i]->Draw();
            sideband_histos[i]->Draw("same");
        }

        auto rp = new TRatioPlot(sideband_histos[i], splot_histos[i], "diffsig");
        c.SetTicks(0, 1);
        rp->SetH1DrawOpt("E");
        rp->Draw();
        rp->GetLowerRefYaxis()->SetTitle("Pull");
        rp->GetUpperRefYaxis()->SetTitle( TString::Format("Normalised entries / (%g)",(var->getMax() - var->getMin())/nbins) );
        c.Update();

        TLegend *leg2 = new TLegend(0.6, 0.8, 0.9, 0.9);
        leg2->AddEntry(sideband_histos[i], "Sideband subtraction", "lp");
        leg2->AddEntry(splot_histos[i], "sPlot", "lp");
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetTextSize(0.03);
        leg2->Draw("same");
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_%i/Splot_vs_sideband_plots/var_%i.pdf",species,i));
    }

}

