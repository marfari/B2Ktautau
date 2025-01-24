using namespace RooStats;
using namespace RooFit;
using namespace std;

void validate_fit(Bool_t fit_status, RooAbsPdf* model, RooRealVar* mass, RooRealVar* var, Int_t species, Int_t year);

void fit_mass(Int_t year, Int_t species)
{
    Bool_t isKtautauMC = false;
    if((species == 1) || (species == 10) || (species == 11) || (species == 12))
    {
        isKtautauMC = true;
    }

    RooDataSet* data;
    if((species == 81) || (species == 5)) // merge 3 years because of low statistics
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
    else // fit per year
    {
        TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_%i/mass_dataset.root",year,species));
        RooWorkspace* w = (RooWorkspace*)f->Get("w");
        data = (RooDataSet*)w->data("data");
    }

    RooRealVar* mass;
    if((species == 4) || (species == 5)) // D+D-K+
    {
        mass = new RooRealVar("mass", "mass", 5180, 5600, "MeV");
    }
    else if((species == 7) || (species == 71) || (species == 8) || (species == 81)) // D0bar D+s
    {
        mass = new RooRealVar("mass", "mass", 5235, 5355, "MeV");
    }
    else if((species == 9) || (species == 0))
    {
        mass = new RooRealVar("mass", "mass", 5240, 5340, "MeV");
    }
    else if(isKtautauMC || (species == 2) || (species == 3))
    {
        mass = new RooRealVar("df_Bp_M", "df_Bp_M", 4000, 8000, "MeV");
    }
   
    RooWorkspace* w_out;
    RooFitResult* fit;
    RooFitResult* fit1;
    if(species == 5)
    {
        mass->setBins(84);
        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 17, 0.5, 30., "MeV");
        RooRealVar frac("frac", "frac", 0., 1.);

        RooGaussian* gauss = new RooGaussian("gauss", "gauss", *mass, mu, sig);

        RooRealVar lambda("lambda", "#lambda", -0.01,-5.,1.0);
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);

        Double_t mass_peak = 5279.34;
        Double_t n_signal_initial = data->sumEntries(TString::Format("(abs(mass-%g)<25)", mass_peak));
        Double_t n_combinatorial_initial = data->sumEntries() - n_signal_initial;
        cout << "n_signal_initial = " << n_signal_initial << endl;
        cout << "n_combinatorial_initial = " << n_combinatorial_initial << endl;

        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0.,data->sumEntries());
        RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial", "n_combinatorial", n_combinatorial_initial, 0., (data->sumEntries())*2);

        RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*gauss,*exp),RooArgList(*n_signal,*n_combinatorial));

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Save(), RooFit::Extended(kTRUE));
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
        model->plotOn(massframe, RooFit::Name("Fit"), RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("Signal"), RooFit::Components("gauss"), RooFit::LineColor(kBlue), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("Background"), RooFit::Components("exp"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->paramOn(massframe, RooFit::Layout(0.55,0.85, 0.9));
        massframe->SetXTitle("m_{B^{+}} (MeV)");
        massframe->Draw();

        double n_float_params = fit->floatParsFinal().getSize();
        double chis = massframe->chiSquare("Fit", "Data", n_float_params);
        TLatex* tex = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %f",chis));//%.3lf
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.035);
        tex->Draw("same");

        TLegend* leg = new TLegend(0.4, 0.6, 0.5, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Signal","Signal");
        leg->AddEntry("Background","Background");
        leg->AddEntry("Fit","Fit");
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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit.gif",species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit.pdf",species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if(species == 0)
    {
        Int_t nbins = 80;
        mass->setBins(nbins);
        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 17, 0.5, 30., "MeV");
        RooRealVar sig2("sig2", "#sigma2", 7, 0.5, 30., "MeV");
        RooRealVar frac("frac", "frac", 0., 1.);

        RooGaussian* gauss = new RooGaussian("gauss", "gauss", *mass, mu, sig);
        RooGaussian* gauss2 = new RooGaussian("gauss2", "gauss2", *mass, mu, sig2);
        RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*gauss,*gauss2), frac);
    
        RooRealVar alpha("alpha", "#alpha", 2, 0., 20.);
        RooRealVar n("n", "n", 2, 0, 300);
        n.setConstant();
        RooCBShape* cb = new RooCBShape("cb", "cb", *mass, mu, sig, alpha, n);

        RooRealVar lambda("lambda", "#lambda", -0.01,-5.,1.0);
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);

        RooRealVar a1("a1", "a1", 0.05, -10, 10);
        RooPolynomial* poly = new RooPolynomial("poly", "poly", *mass, RooArgList(a1));

        Double_t mass_peak = 5279.34;
        Double_t n_signal_initial = data->sumEntries(TString::Format("(abs(mass-%g)<25)", mass_peak));
        Double_t n_combinatorial_initial = data->sumEntries() - n_signal_initial;
        cout << "n_signal_initial = " << n_signal_initial << endl;
        cout << "n_combinatorial_initial = " << n_combinatorial_initial << endl;

        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0.,data->sumEntries());
        RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial", "n_combinatorial", n_combinatorial_initial, 0., (data->sumEntries())*2);

        RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*gauss,*exp),RooArgList(*n_signal,*n_combinatorial));

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Save(), RooFit::Extended(kTRUE));
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
        model->plotOn(massframe, RooFit::Name("Fit"), RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("Signal"), RooFit::Components("gauss"), RooFit::LineColor(kBlue), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->plotOn(massframe, RooFit::Name("Background"), RooFit::Components("exp"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->paramOn(massframe, RooFit::Layout(0.55,0.85, 0.9));
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

        TLegend* leg = new TLegend(0.4, 0.6, 0.5, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Signal","Signal");
        leg->AddEntry("Background","Background");
        leg->AddEntry("Fit","Fit");
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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if(species == 4)
    {
        mass->setRange("fit_range", 5180, 5600);

        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar mu1("mu1", "#mu1", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 7., 0.5, 100., "MeV");
        RooRealVar sig1("sig1", "#sigma_{1}", 3, 0.5, 100., "MeV");
        RooRealVar sig2("sig2", "#sigma_{2}", 10, 0.5, 500., "MeV");

        RooRealVar alpha("alpha", "#alpha", 0., 20.);
        RooRealVar n("n", "n", 0, 300);
        RooRealVar frac("frac", "f", 0.3, 0, 1);
        RooRealVar frac1("frac1", "f1", 0.3, 0, 1);

        RooCBShape* cb = new RooCBShape("cb", "cb", *mass, mu, sig, alpha, n);

        RooGaussian* g = new RooGaussian("g", "g", *mass, mu, sig);
        RooGaussian* g1 = new RooGaussian("g1", "g1", *mass, mu, sig1);
        RooGaussian* g2 = new RooGaussian("g2", "g2", *mass, mu, sig2);

        RooRealVar xi("xi", "xi", -100. ,100.);
        RooRealVar ro1("ro1", "ro1", -5., 5.);
        RooRealVar ro2("ro2", "ro2", -5., 5.);

        RooBukinPdf* pdf = new RooBukinPdf("pdf", "pdf", *mass, mu, sig, xi, ro1, ro2); 
        // RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(*bu,*g1), frac);

        Double_t mass_peak = 5279.34;
        Double_t n_signal_initial = data->sumEntries();
        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0, (data->sumEntries())*2);

        RooExtendPdf* model = new RooExtendPdf("model", "model", *pdf, *n_signal);

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Range("fit_range"));
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
        model->paramOn(massframe, RooFit::Layout(0.55,0.85, 0.9));
        massframe->SetXTitle("m_{B^{+}} (MeV)");
        massframe->Draw();

        double n_float_params = fit->floatParsFinal().getSize();
        double chis = massframe->chiSquare("Fit", "Data", n_float_params);
        TLatex* tex = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %f",chis));//%.3lf
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.035);
        tex->Draw("same");

        TLegend* leg = new TLegend(0.4, 0.6, 0.5, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Fit","Fit");
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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if(species == 7)
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
        RooRealVar n("n", "n", 50, 0, 300);
        RooRealVar frac("frac", "f", 0.5, 0, 1);
        n.setConstant();

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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if(species == 81)
    {
        Double_t nbins = 46;
        mass->setBins(nbins);
        RooRealVar a1("a1", "a1", 0.05, -10, 10);
        RooRealVar a2("a2", "a2", -100, 100);
        RooRealVar a3("a3", "a3", 0, 100);

        RooRealVar lambda("lambda", "lambda", -0.05, -10, 5);
        RooRealVar lambda1("lambda1", "lambda1", -10, 5);

        RooRealVar frac("frac", "frac", 0.5, 0, 1);
        RooRealVar frac1("frac1", "frac1", 0.3, 0, 1);

        RooPolynomial* poly = new RooPolynomial("poly", "poly", *mass, RooArgList(a1));
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);
        RooExponential* exp1 = new RooExponential("exp1", "exp1", *mass, lambda1);

        RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(*exp,*poly), frac);

        Double_t n_initial = data->sumEntries();
        RooRealVar* N = new RooRealVar("N", "N", n_initial, 0.,2*(data->sumEntries()));

        RooExtendPdf* model = new RooExtendPdf("model", "model", *pdf, *N);

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
        TLine* l = new TLine(5235,0, 5350,0);
        l->Draw("same");
        l->SetLineColor(kBlue);
        l->SetLineStyle(2);
        l->Draw("same");

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit.gif",species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit.pdf",species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);
    
        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, N, species, year);
    }
    else if(species == 8)
    {
        Double_t nbins = 30;
        mass->setBins(nbins);

        // Signal PDF
        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar mu1("mu1", "#mu1", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 7., 0.5, 100., "MeV");
        RooRealVar sig1("sig1", "#sigma_{1}", 3, 0.5, 100., "MeV");
        RooRealVar sig2("sig2", "#sigma_{2}", 10, 0.5, 500., "MeV");

        RooRealVar alpha("alpha", "#alpha", 2, 0., 20.);
        RooRealVar alpha1("alpha1", "#alpha1", -2, -20., 5.);
        RooRealVar n("n", "n", 50, 0, 300);
        n.setConstant();

        RooGaussian* g = new RooGaussian("g", "g", *mass, mu, sig);
        RooGaussian* g2 = new RooGaussian("g2", "g2", *mass, mu, sig2);

        TFile* fr = new TFile(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_7/mass_fit_result.root",year));
        RooFitResult* fitres = (RooFitResult*)fr->Get("fitresult_model_data");
        RooRealVar* frac_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("frac"));
        RooRealVar* alpha_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("alpha"));
        RooRealVar* alpha1_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("alpha1"));
        RooRealVar* sig1_mc = (RooRealVar*)fitres->floatParsFinal().find(Form("sig1"));
        frac_mc->setConstant();
        alpha_mc->setConstant();
        alpha1_mc->setConstant();
        sig1_mc->setConstant();

        RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *mass, mu, sig, *alpha_mc, n);
        RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *mass, mu, *sig1_mc, *alpha1_mc, n);

        RooRealVar xi( "xi", "xi", 0., -1, 1 );
        xi.setConstant();

        RooRealVar ro1("ro1", "ro1", -0.05, -1., 1.);
        RooRealVar ro2("ro2", "ro2", 0.05, -1., 1.);
        RooRealVar frac("frac", "frac", 0.5, 0, 1);

        RooBukinPdf* bukin = new RooBukinPdf("bukin", "bukin", *mass, mu, sig, xi, ro1, ro2); 
        RooGaussian* g1 = new RooGaussian("g1", "g1", *mass, mu, *sig1_mc);

        RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*cb2), *frac_mc); 
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*g1,*g2), frac);
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*bukin,*g1), frac);

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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if((species == 85) || (species == 86) || (species == 87))
    {
        Double_t nbins = 30;
        mass->setBins(nbins);

        // Signal PDF
        RooRealVar mu("mu", "#mu", 5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 7., 0.5, 100., "MeV");

        RooGaussian* sig_pdf = new RooGaussian("sig_pdf", "sig_pdf", *mass, mu, sig);

        Double_t mass_peak = 5279.34;
        Double_t n_signal_initial = data->sumEntries(TString::Format("(abs(mass-%g)<20)", mass_peak));
        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0.,2*(data->sumEntries()));

        // Background PDF
        RooRealVar lambda("lambda", "lambda", -0.1, -10, 10);
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);

        Double_t n_combinatorial_initial = data->sumEntries() - n_signal_initial;
        RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial", "n_combinatorial", n_combinatorial_initial, 0.,2*(data->sumEntries()));

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
        // leg->Draw("same");

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


        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if(species == 71)
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
        RooRealVar n("n", "n", 50, 0, 300);
        RooRealVar frac("frac", "f", 0.5, 0, 1);
        n.setConstant();

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
        RooRealVar frac1("frac1", "frac1", 0.5, 0, 1);

        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*cb2), frac); 
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*g1,*g2), frac);
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*bukin,*g1), frac);
        // RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*g1), frac); 
        RooAddPdf* sig_pdf = new RooAddPdf("sig_pdf", "sig_pdf", RooArgList(*cb1,*cb2,*g2), RooArgList(frac,frac1)); 

        RooRealVar a1("a1", "a1", 0.5, 0, 5);
        RooPolynomial *poly = new RooPolynomial("poly", "poly", *mass, RooArgList(a1));

        RooRealVar lambda("lambda", "lambda", -0.1, -10, 10);
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);
        RooAddPdf* tot_pdf = new RooAddPdf("tot_pdf", "tot_pdf", RooArgList(*sig_pdf,*poly), frac1);

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
        model->plotOn(massframe, RooFit::Name("Exp"), RooFit::Components("poly"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
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
        leg->AddEntry("Poly","Poly");
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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
        // validate_fit(fit_status, model, mass, n_signal, species, year);
    }
    else if(isKtautauMC)
    {
        mass->setRange("fit_range", 4000, 8000);
        mass->setBins(40);

        RooRealVar mu("mu", "#mu", 5279.41,5240., 5310., "MeV");
        RooRealVar sig("sig", "#sigma", 200., 0.5, 900., "MeV");
        RooRealVar sig1("sig1", "#sigma1", 50., 0.5, 900., "MeV");
        RooRealVar sig2("sig2", "#sigma1", 300., 0.5, 900., "MeV");

        RooRealVar alpha("alpha", "#alpha", 0., 20.);
        RooRealVar alpha1("alpha1", "#alpha1", -20., 0.);
        RooRealVar n("n", "n", 0, 30);
        RooRealVar frac("frac", "frac", 0.5, 0, 1);
        n.setConstant();

        RooGaussian* g = new RooGaussian("g", "g", *mass, mu, sig);
        RooGaussian* g1 = new RooGaussian("g1", "g1", *mass, mu, sig1);
        RooCBShape* cb = new RooCBShape("cb", "cb", *mass, mu, sig, alpha, n);
        RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *mass, mu, sig1, alpha1, n);
        
        RooAddPdf* pdf2 = new RooAddPdf("pdf2", "pdf2", RooArgList(*g,*cb), frac);
        RooAddPdf* pdf4 = new RooAddPdf("pdf4", "pdf4", RooArgList(*cb,*cb1), frac);
        RooAddPdf* pdf5 = new RooAddPdf("pdf5", "pdf5", RooArgList(*g,*cb1), frac);
        RooAddPdf* pdf6 = new RooAddPdf("pdf6", "pdf6", RooArgList(*g,*g1), frac);

        RooRealVar xi("xi", "xi", -100. ,100.);
        RooRealVar ro1("ro1", "ro1", -5., 5.);
        RooRealVar ro2("ro2", "ro2", -5., 5.);
        RooRealVar ro3("ro3", "ro3", -5., 5.);
        RooBukinPdf* buk = new RooBukinPdf("buk", "buk", *mass, mu, sig, xi, ro1, ro2); 
        RooBukinPdf* buk1 = new RooBukinPdf("buk1", "buk1", *mass, mu, sig1, xi, ro1, ro3); 
        RooBukinPdf* buk2 = new RooBukinPdf("buk2", "buk2", *mass, mu, sig2, xi, ro1, ro3); 

        RooBifurGauss* bg = new RooBifurGauss("bg", "bg", *mass, mu, sig, sig1);

        RooAddPdf* pdf3 = new RooAddPdf("pdf3", "pdf3", RooArgList(*g1,*buk), frac);
        RooAddPdf* pdf7 = new RooAddPdf("pdf7", "pdf7", RooArgList(*cb1,*buk), frac);
        RooAddPdf* pdf8 = new RooAddPdf("pdf8", "pdf8", RooArgList(*buk,*buk1), frac);
        RooAddPdf* pdf9 = new RooAddPdf("pdf9", "pdf9", RooArgList(*bg,*buk2), frac);

        Double_t mass_peak = 5279.41;
        Double_t n_signal_initial = data->sumEntries();
        RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0, (data->sumEntries())*2);

        RooExtendPdf* model = new RooExtendPdf("model", "model", *pdf3, *n_signal);

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Range("fit_range"));
        fit->Print();

        TF1 * f = model->asTF( RooArgList(*mass) );
        Double_t max = f->GetMaximum();
        Double_t x_max = f->GetMaximumX();
        cout << "max = " << max << endl;
        cout << "x_max = " << x_max << endl;

        double fwhm_left = f->GetX(max/2,4000, x_max);
        double fwhm_right = f->GetX(max/2, x_max, 8000);
        cout << "fwhm_left = " << fwhm_left << endl;
        cout << "fwhm_right = " << fwhm_right << endl;
        double fwhm = fwhm_right - fwhm_left;
        cout << "FWHM = " << fwhm << endl;
        cout << "sig = " << fwhm/2.355 << endl;

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
        // model->plotOn(massframe, RooFit::Name("Gauss"), RooFit::Components("g"), RooFit::LineColor(kGreen), RooFit::LineStyle(1), RooFit::LineWidth(2));
        // model->plotOn(massframe, RooFit::Name("Bukin"), RooFit::Components("buk1"), RooFit::LineColor(kCyan), RooFit::LineStyle(1), RooFit::LineWidth(2));
        model->paramOn(massframe, RooFit::Layout(0.55,0.85, 0.9));
        massframe->SetXTitle("m_{B^{+}} (MeV)");
        massframe->Draw();

        double n_float_params = fit->floatParsFinal().getSize();
        double chis = massframe->chiSquare("Fit", "Data", n_float_params);
        TLatex* tex = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %f",chis));//%.3lf
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.035);
        tex->Draw("same");

        TLegend* leg = new TLegend(0.4, 0.6, 0.5, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Fit","Fit");
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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
    }
    else if(species == 2)
    {
        mass->setRange("left_sideband", 4000, 4700);
        mass->setRange("right_sideband", 5810, 8000);
        mass->setBins(40);

        RooRealVar lambda("lambda", "lambda", -0.1, -10, 10);
        // RooRealVar lambda1("lambda1", "lambda1", -0.1, -10, 10);
        RooRealVar frac("frac", "frac", 0.5, 0., 1.);
        RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);
        // RooExponential* exp1 = new RooExponential("exp1", "exp1", *mass, lambda1);
        // RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(*exp,*exp1), frac);

        RooRealVar mu("mu", "#mu", 4700., 4800., "MeV");
        RooRealVar sig("sig", "#sigma", 7., 0.5, 100., "MeV");
        RooRealVar xi("xi", "xi", -100. ,100.);
        RooRealVar ro1("ro1", "ro1", -5., 5.);
        RooRealVar ro2("ro2", "ro2", -5., 5.);

        RooBukinPdf* bukin = new RooBukinPdf("bukin", "bukin", *mass, mu, sig, xi, ro1, ro2); 
        RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(*exp,*bukin), frac);

        Double_t n_bkg_initial = data->sumEntries();
        RooRealVar* n_bkg = new RooRealVar("n_bkg", "n_bkg", n_bkg_initial, 0.,2*(data->sumEntries()));

        RooExtendPdf* model = new RooExtendPdf("model", "model", *pdf, *n_bkg);

        fit = model->fitTo(*data, RooFit::Minos(true), RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Range("left_sideband,right_sideband"));
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
        model->paramOn(massframe, RooFit::Layout(0.55,0.85, 0.9));
        massframe->SetXTitle("m_{B^{+}} (MeV)");
        massframe->Draw();

        double n_float_params = fit->floatParsFinal().getSize();
        double chis = massframe->chiSquare("Fit", "Data", n_float_params);
        TLatex* tex = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %f",chis));//%.3lf
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.035);
        tex->Draw("same");

        TLegend* leg = new TLegend(0.4, 0.6, 0.5, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->AddEntry("Data","Data");
        leg->AddEntry("Fit","Fit");
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

        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.gif",year,species));
        c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit.pdf",year,species));

        w_out = new RooWorkspace("w_out");
        w_out->import(*data);
        w_out->import(*model);

        Int_t fit_status = fit->status();
    }

    if(species == 5)
    {
        TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/mass_fit_result.root",species),"RECREATE");
        fout->cd();
        w_out->Write();
        fit->Write();
        fout->Close();
    }
    else
    {
        TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/mass_fit_result.root",year,species),"RECREATE");
        fout->cd();
        w_out->Write();
        fit->Write();
        fout->Close();
    } 

}

void validate_fit(Bool_t fit_status, RooAbsPdf* model, RooRealVar* mass, RooRealVar* var, Int_t species, Int_t year)
{
    // Validate fit
    if(fit_status == 0)
    {
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

        if((species == 81) || (species == 5))
        {
            c_pull->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/pulls_poisson.gif",species));
            c_pull->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/pulls_poisson.pdf",species));
            c_params->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/params_poisson.gif",species));
            c_params->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_%i/params_poisson.pdf",species));
        }
        else
        {
            c_pull->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/pulls_poisson.gif",year,species));
            c_pull->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/pulls_poisson.pdf",year,species));
            c_params->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/params_poisson.gif",year,species));
            c_params->SaveAs(Form("/panfs/felician/B2Ktautau/workflow/fit_mass/201%i/Species_%i/params_poisson.pdf",year,species));
        }

    }  
}