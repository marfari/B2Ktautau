using namespace RooStats;
using namespace RooFit;
using namespace std;

void fit_mass(Int_t year, Int_t species)
{
    TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_%i/mass_dataset.root",year,species));
    RooWorkspace* w = (RooWorkspace*)f->Get("w");
    RooDataSet* data = (RooDataSet*)w->data("data");
    RooRealVar* mass;
    if((species == 4) || (species == 5))
    {
        mass = new RooRealVar("mass", "mass", 5180, 5600, "MeV");
    }

    RooRealVar mu("mu", "mu", 5240., 5310.);
    RooRealVar sig("sig", "sig", 17, 0.5, 30.);
    RooRealVar frac("frac", "frac", 0., 1.);

    RooGaussian* gauss = new RooGaussian("gauss", "gauss", *mass, mu, sig);

    RooRealVar lambda("lambda", "lambda", -1,-5.,1.0);
    RooExponential* exp = new RooExponential("exp", "exp", *mass, lambda);

    Double_t mass_peak = 5279;
    Double_t n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<25", mass_peak));
    Double_t n_combinatorial_initial = data->sumEntries() - n_signal_initial;

    RooRealVar* n_signal = new RooRealVar("n_signal", "n_signal", n_signal_initial, 0.,(data->sumEntries())*2);
    RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial", "n_combinatorial", n_combinatorial_initial, 0., data->sumEntries());

    RooAddPdf model("model","model",RooArgList(*gauss,*exp),RooArgList(*n_signal,*n_combinatorial));

    RooFitResult* fit = model.fitTo(*data, RooFit::Minos(true), RooFit::Save(), RooFit::Extended(kTRUE));
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
    model.plotOn(massframe, RooFit::Name("Fit"), RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::LineWidth(2));
    model.plotOn(massframe, RooFit::Name("Signal"), RooFit::Components("gauss"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
    model.plotOn(massframe, RooFit::Name("Background"), RooFit::Components("exp"), RooFit::LineColor(kBlue), RooFit::LineStyle(1), RooFit::LineWidth(2));
    model.paramOn(massframe,Layout(0.60,0.99,0.8));
    massframe->SetXTitle("m_{B^{+}} (MeV)");
    massframe->Draw();

    double n_float_params = fit->floatParsFinal().getSize();
    double chis = massframe->chiSquare("Fit", "Data", n_float_params);
    TLatex* tex = new TLatex(0.15, 0.75, Form("#chi^{2}/ndf = %f",chis));//%.3lf
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetTextSize(0.035);
    tex->Draw("same");

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

}