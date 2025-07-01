using namespace RooStats;
using namespace RooFit;
using namespace std;

void validate_fit(RooAbsPdf* model, RooRealVar* mass, RooRealVar* var, Int_t species);

void fit_mass(Int_t species)
{

    TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/Species_%i/mass_dataset.root",species));
    RooWorkspace* w = (RooWorkspace*)f->Get("w");
    RooDataSet *data = (RooDataSet*)w->data("data");

    RooRealVar* mass = new RooRealVar("mass", "mass", 5235, 5355, "MeV");
   
    RooWorkspace* w_out;
    RooFitResult* fit;
    RooFitResult* fit1;    
    if(species == 72)
    {
        Double_t nbins = 30;
        mass->setBins(nbins);

        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar mu1("mu1", "#mu1", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 7., 0.5, 100., "MeV");
        RooRealVar sig1("sig1", "#sigma_{1}", 3, 0.5, 100., "MeV");
        RooRealVar sig2("sig2", "#sigma_{2}", 10, 0.5, 500., "MeV");

        RooRealVar alpha("alpha", "#alpha", 2, 0., 20.);
        RooRealVar alpha1("alpha1", "#alpha1", -2, -20., 5.);
        RooRealVar n("n", "n", 10, 0, 300);
        RooRealVar frac("frac", "f", 0.5, 0, 1);
        // n.setConstant();

        RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *mass, mu, sig, alpha, n);
        RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *mass, mu, sig1, alpha1, n);

        RooGaussian* g = new RooGaussian("g", "g", *mass, mu, sig);
        RooGaussian* g1 = new RooGaussian("g1", "g1", *mass, mu, sig1);
        RooGaussian* g2 = new RooGaussian("g2", "g2", *mass, mu, sig2);

        RooRealVar xi("xi", "xi", 0., -1. ,1.);
        RooRealVar ro1("ro1", "ro1", -0.05, -1., 1.);
        RooRealVar ro2("ro2", "ro2", 0.05, -1., 1.);
        xi.setConstant();

        RooBukinPdf* bukin = new RooBukinPdf("bukin", "bukin", *mass, mu, sig, xi, ro1, ro2); 

        RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*cb2), frac); 
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*g1,*g2), frac);
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*bukin,*g1), frac);
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*g1), frac); 

        Double_t n_signal_initial = data->sumEntries();
        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0.,2*(data->sumEntries()));

        RooExtendPdf* model = new RooExtendPdf("model", "model", *sig_pdf, *n_signal);

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Extended(kTRUE), RooFit::Save());
        fit->Print();

        TCanvas c("c", "c", 2000,1500);
        c.SetTitle("");

        TPad *p1 = new TPad("p1","p1",0.,0.27,1.,1.);
        p1->SetTitle("");
        p1->SetBorderMode(1);
        p1->SetFrameBorderMode(0);
        p1->SetBorderSize(2);
        p1->SetBottomMargin(0.10);
        p1->Draw();

        TPad *p2 = new TPad("p2","p2",0.,0.075,1.,0.25);
        p2->SetTitle("");
        p2->SetTopMargin(0.);
        p2->SetBottomMargin(0.2);
        p2->SetBorderMode(1);
        p2->Draw();

        p1->cd();
        RooPlot* massframe = mass->frame();
        data->plotOn(massframe, RooFit::Name("Data"));
        model->plotOn(massframe, RooFit::Name("Fit"), RooFit::LineColor(kBlue), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("CB1"), RooFit::Components("cb1"), RooFit::LineColor(kCyan), RooFit::LineStyle(3), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("CB2"), RooFit::Components("cb2"), RooFit::LineColor(kGreen), RooFit::LineStyle(3), RooFit::LineWidth(2));
        model->paramOn(massframe, RooFit::Layout(0.6,0.85, 0.9));
        massframe->SetXTitle("m_{B^{+}} (MeV)");
        massframe->Draw();

        double n_float_params = fit->floatParsFinal().getSize();
        double chis = massframe->chiSquare("Fit", "Data", n_float_params);
        TLatex* tex = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %f",chis));//%.3lf
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.035);
        tex->Draw("same");

        Double_t ndf = nbins - n_float_params;
        Double_t pvalue =  TMath::Prob(chis*ndf, ndf);
        TLatex* tex1 = new TLatex(0.13, 0.75, Form("p = %f",pvalue));//%.3lf
        tex1->SetNDC(kTRUE);
        tex1->SetTextFont(42);
        tex1->SetTextSize(0.035);
        tex1->Draw("same");

        TLegend* leg = new TLegend(0.5, 0.6, 0.6, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Fit","Fit");
        leg->AddEntry("CB1","CB1");
        leg->AddEntry("CB2","CB2");
        leg->Draw("same");

        RooHist* pull_hist = massframe->pullHist("Data","Fit");
        RooPlot* pull_plot = mass->frame(Title(""));

        pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
        pull_plot->SetTitle("");

        pull_plot->GetXaxis()->SetTitle("");
        pull_plot->GetXaxis()->SetTitleSize(0.15);
        pull_plot->GetXaxis()->SetTitleOffset(0.9);
        pull_plot->GetXaxis()->SetLabelSize(0.15);
        pull_plot->GetXaxis()->SetLabelOffset(0.01);
        pull_plot->GetXaxis()->SetTickLength(0.13);

        pull_plot->GetYaxis()->SetTitle("Pull");
        pull_plot->GetYaxis()->SetTitleSize(0.13);
        pull_plot->GetYaxis()->SetTitleOffset(0.18);
        pull_plot->GetYaxis()->SetLabelSize(0.13);
        pull_plot->GetYaxis()->SetLabelOffset(0.005);
        pull_plot->GetYaxis()->SetNdivisions(305);

        p2->cd();
        pull_plot->Draw();
        TLine* l = new TLine(5235,0, 5355,0);
        l->Draw("same");
        l->SetLineColor(kBlue);
        l->SetLineStyle(2);
        l->Draw("same");

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit.pdf",species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        if(fit_status == 0)
        {
            validate_fit(model, mass, n_signal, species);
        }
    }
    else if(species == 8)
    {
        Double_t nbins = 30;
        mass->setBins(nbins);

        // Signal PDF    
        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 7., 0.5, 100., "MeV");
        RooRealVar sig1("sig1", "#sigma_{1}", 3, 0.5, 100., "MeV");
        RooRealVar n("n", "n", 10, 0, 300);
        // n.setConstant();

        TFile* fr = new TFile("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_72/mass_fit_result.root");
        
        RooFitResult* fitres = (RooFitResult*)fr->Get("fitresult_model_data");
        RooRealVar* n_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("n"));
        RooRealVar* frac_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("frac"));
        RooRealVar* sig1_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("sig1"));
        RooRealVar* alpha_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("alpha"));
        RooRealVar* alpha1_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("alpha1"));

        n_mc->setConstant();
        frac_mc->setConstant();
        sig1_mc->setConstant();
        alpha_mc->setConstant();
        alpha1_mc->setConstant();

        RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *mass, mu, sig, *alpha_mc, *n_mc);
        RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *mass, mu, sig1, *alpha1_mc, *n_mc);

        RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*cb2), *frac_mc);

        Double_t mass_peak = 5279.34;
        Double_t n_signal_initial = data->sumEntries(TString::Format("(abs(mass-%g)<20)", mass_peak));
        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0.,2*(data->sumEntries()));

        // Background PDF
        RooRealVar a1("a1", "a1", 0.5, 0, 5);
        RooPolynomial *poly = new RooPolynomial("poly", "poly", *mass, RooArgList(a1));

        RooRealVar lambda("lambda", "lambda", -0.1, -10, 10);
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);

        Double_t n_combinatorial_initial = data->sumEntries() - n_signal_initial;
        RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial", "n_combinatorial", n_combinatorial_initial, 0.,2*(data->sumEntries()));

        // Charmless B (D+s sideband)
        // TFile* fr = new TFile("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_81/mass_fit_result.root");
        // RooFitResult* fitres = (RooFitResult*)fr->Get("fitresult_model_data");
        // RooRealVar* a1_charmless = (RooRealVar*)fitres->floatParsFinal().find(Form("a1"));
        // RooRealVar* frac_charmless = (RooRealVar*)fitres->floatParsFinal().find(Form("frac"));
        // RooRealVar* lambda_charmless = (RooRealVar*)fitres->floatParsFinal().find(Form("lambda"));
        // a1_charmless->setConstant();
        // frac_charmless->setConstant();
        // lambda_charmless->setConstant();

        // RooRealVar* N_charmless = new RooRealVar("N_charmless", "N_charmless", 0, data->sumEntries());

        // RooPolynomial* poly_charmless = new RooPolynomial("poly_charmless", "poly_charmless", *mass, RooArgList(*a1_charmless));
        // RooExponential* exp_charmless = new RooExponential("exp_charmless", "exp_charmless", *mass, *lambda_charmless);

        // RooAddPdf* pdf_charmless = new RooAddPdf("pdf_charmless", "pdf_charmless", RooArgList(*exp_charmless,*poly_charmless), *frac_charmless);

        // Total PDF
        RooAddPdf* model = new RooAddPdf("model", "model", RooArgList(*sig_pdf,*exp), RooArgList(*n_signal,*n_combinatorial));

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Extended(kTRUE), RooFit::Save());
        fit->Print();

        TCanvas c("c", "c", 2000,1500);
        c.SetTitle("");

        TPad *p1 = new TPad("p1","p1",0.,0.27,1.,1.);
        p1->SetTitle("");
        p1->SetBorderMode(1);
        p1->SetFrameBorderMode(0);
        p1->SetBorderSize(2);
        p1->SetBottomMargin(0.10);
        p1->Draw();

        TPad *p2 = new TPad("p2","p2",0.,0.075,1.,0.25);
        p2->SetTitle("");
        p2->SetTopMargin(0.);
        p2->SetBottomMargin(0.2);
        p2->SetBorderMode(1);
        p2->Draw();

        p1->cd();
        RooPlot* massframe = mass->frame();
        data->plotOn(massframe, RooFit::Name("Data"));
        model->plotOn(massframe, RooFit::Name("Fit"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("Signal"),RooFit::Components("sig_pdf"), RooFit::LineColor(kBlue), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("CB1"),RooFit::Components("cb1"), RooFit::LineColor(kBlue), RooFit::LineStyle(3), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("CB2"),RooFit::Components("cb2"), RooFit::LineColor(kCyan), RooFit::LineStyle(3), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("Comb. bkg."),RooFit::Components("exp"), RooFit::LineColor(kGreen+1), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->paramOn(massframe, RooFit::Layout(0.6,0.85, 0.9));
        massframe->SetXTitle("m_{B^{+}} (MeV)");
        massframe->Draw();

        double n_float_params = fit->floatParsFinal().getSize();
        double chis = massframe->chiSquare("Fit", "Data", n_float_params);
        TLatex* tex = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %f",chis));//%.3lf
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.035);
        tex->Draw("same");

        Double_t ndf = nbins - n_float_params;
        Double_t pvalue =  TMath::Prob(chis*ndf, ndf);
        TLatex* tex1 = new TLatex(0.13, 0.75, Form("p = %f",pvalue));//%.3lf
        tex1->SetNDC(kTRUE);
        tex1->SetTextFont(42);
        tex1->SetTextSize(0.035);
        tex1->Draw("same");

        TLegend* leg = new TLegend(0.13, 0.3, 0.3, 0.7);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Fit","Fit");
        leg->AddEntry("Signal","Signal");
        leg->AddEntry("Comb. bkg.","Comb. bkg.");
        leg->AddEntry("CB1","CB1");
        leg->AddEntry("CB2","CB2");
        leg->Draw("same");

        RooHist* pull_hist = massframe->pullHist("Data","Fit");
        RooPlot* pull_plot = mass->frame(Title(""));

        pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
        pull_plot->SetTitle("");

        pull_plot->GetXaxis()->SetTitle("");
        pull_plot->GetXaxis()->SetTitleSize(0.15);
        pull_plot->GetXaxis()->SetTitleOffset(0.9);
        pull_plot->GetXaxis()->SetLabelSize(0.15);
        pull_plot->GetXaxis()->SetLabelOffset(0.01);
        pull_plot->GetXaxis()->SetTickLength(0.13);

        pull_plot->GetYaxis()->SetTitle("Pull");
        pull_plot->GetYaxis()->SetTitleSize(0.13);
        pull_plot->GetYaxis()->SetTitleOffset(0.18);
        pull_plot->GetYaxis()->SetLabelSize(0.13);
        pull_plot->GetYaxis()->SetLabelOffset(0.005);
        pull_plot->GetYaxis()->SetNdivisions(305);

        p2->cd();
        pull_plot->Draw();
        TLine* l = new TLine(5235,0, 5355,0);
        l->Draw("same");
        l->SetLineColor(kBlue);
        l->SetLineStyle(2);
        l->Draw("same");

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit.pdf",species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        if(fit_status == 0)
        {
            validate_fit(model, mass, n_signal, species);
        }

    }

    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit_result.root",species),"RECREATE");
    fout->cd();
    w_out->Write();
    fit->Write();
    fout->Close();

}

void validate_fit(RooAbsPdf* model, RooRealVar* mass, RooRealVar* var, Int_t species)
{
    // Validate fit
    RooMCStudy* mcstudy = new RooMCStudy(*model, *mass, Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));
    mcstudy->generateAndFit(5000);

    RooPlot* framesPull, *framesParam;
    framesPull = mcstudy->plotPull(*var,FrameBins(200));//FrameRange(-5,5)
    framesPull->SetTitle("");
    framesParam = mcstudy->plotParam(*var,FrameBins(50));
    framesParam->SetTitle("");

    TGraph* h = static_cast<TGraph*>(framesPull->getObject(0));

    gStyle->SetOptFit(0111);

    TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

    gPad->SetLeftMargin(0.15);

    c_pull->cd();
    h->SetTitle("");
    h->Draw();
    c_pull->Update();
    h->Fit("gaus","","",-20,20);
    h->GetFunction("gaus")->SetLineColor(4);
    h->GetFunction("gaus")->SetLineWidth(5);
    h->GetXaxis()->SetTitle("Pull");
    h->GetYaxis()->SetTitle("Toy MCs");
    h->GetXaxis()->SetRangeUser(-10,10);
    h->Draw("same");

    TCanvas* c_params = new TCanvas("params", "params", 900, 800);
    c_params->cd();
    framesParam->GetYaxis()->SetTitleOffset(1.4);
    framesParam->Draw();

    c_pull->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/pulls_poisson.pdf",species));
    c_params->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/params_poisson.pdf",species));
    
}