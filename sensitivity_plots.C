std::vector<double> sensitivity_plots_it(int year, int nprob, int bdt_it);

void sensitivity_plots(int year, int nprob)
{
    double BDT_cut[nprob], S[nprob], B[nprob], sig_bdt_eff[nprob], bkg_bdt_eff[nprob], sensitivity[nprob];
    double epsilon = 0.01;

    for(int i = 0; i < nprob; i++)
    {
        std::vector<double> result = sensitivity_plots_it(year, nprob, i);

        BDT_cut[i] = result[0];
        S[i] = result[1];
        B[i] = result[2];
        sig_bdt_eff[i] = result[3];
        bkg_bdt_eff[i] = result[4];

        if(result[3] < epsilon)
        {
            continue;
        }
        else
        {
            sensitivity[i] = result[5];
        }

    }

    TGraph* epsilon_bdt_S = new TGraph(nprob, BDT_cut, sig_bdt_eff);
    TGraph* epsilon_bdt_B = new TGraph(nprob, BDT_cut, bkg_bdt_eff);

    TCanvas c;
    c.cd();
    gPad->SetGrid(1,1);
    epsilon_bdt_S->SetMarkerStyle(8);
    epsilon_bdt_S->SetMarkerColor(kBlue);
    epsilon_bdt_B->SetMarkerStyle(8);
    epsilon_bdt_B->SetMarkerColor(kRed);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(epsilon_bdt_S);
    mg->Add(epsilon_bdt_B);
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("BDT cut: bdt > x");
    mg->GetYaxis()->SetTitle("BDT efficiency");
    mg->SetTitle("");
    TLegend *leg = new TLegend(0.7, 0.7, 0.9,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(epsilon_bdt_S, "Signal BDT eff", "lp");
    leg->AddEntry(epsilon_bdt_B, "Background BDT eff", "lp");
    leg->Draw("same");
    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/bdt_efficiency.gif",year,nprob));
    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/bdt_efficiency.pdf",year,nprob));

    TGraph* g_limit = new TGraph(nprob, BDT_cut, sensitivity);

    TCanvas c1;
    c1.cd();
    gPad->SetGrid(1,1);
    g_limit->SetMarkerStyle(8);
    g_limit->SetMarkerColor(kBlack);
    g_limit->Draw("AP");
    g_limit->GetXaxis()->SetTitle("BDT cut: bdt > x");
    g_limit->GetYaxis()->SetTitle("Sensitivity");
    g_limit->SetTitle("");
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/sensitivity.gif",year,nprob));
    c1.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/sensitivity.pdf",year,nprob));

    TGraph* sig_yield = new TGraph(nprob, BDT_cut, S);
    TGraph* bkg_yield = new TGraph(nprob, BDT_cut, B);

    TCanvas c2;
    c2.cd();
    gPad->SetGrid(1,1);
    sig_yield->SetMarkerStyle(8);
    sig_yield->SetMarkerColor(kBlue);
    sig_yield->Draw("AP");
    sig_yield->GetXaxis()->SetTitle("BDT cut: bdt > x");
    sig_yield->GetYaxis()->SetTitle("Signal Yield");
    sig_yield->SetTitle("");
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/signal_yield.gif",year,nprob));
    c2.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/signal_yield.pdf",year,nprob));

    TCanvas c3;
    c3.cd();
    gPad->SetGrid(1,1);
    bkg_yield->SetMarkerStyle(8);
    bkg_yield->SetMarkerColor(kRed);
    bkg_yield->Draw("AP");
    bkg_yield->GetXaxis()->SetTitle("BDT cut: bdt > x");
    bkg_yield->GetYaxis()->SetTitle("Background Yield");
    bkg_yield->SetTitle("");
    c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/background_yield.gif",year,nprob));
    c3.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201%i/N_prob_%i/background_yield.pdf",year,nprob));

    for(int i = 0; i < nprob; i++)
    {
        if(sig_bdt_eff[i] < epsilon)
        {
            sensitivity[i] = 100;
        }
    }

    double bdt_cut_limit, upper_limit;
    for(int i = 0; i < nprob; i++)
    {
        if(i == 0)
        {
            bdt_cut_limit = BDT_cut[0];
            upper_limit = sensitivity[0];
        }
        else
        {
            if( (sensitivity[i] < upper_limit))
            {
                upper_limit = sensitivity[i];
                bdt_cut_limit = BDT_cut[i];
            }
        }
    }

    cout << Form("bdt > %f minimises 95 percent upper limit on B(Ktautau) < %f", bdt_cut_limit, upper_limit) << endl; 

    return;
}

std::vector<double> sensitivity_plots_it(int year, int nprob, int bdt_it)
{
    TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/sensitivity/201%i/N_prob_%i/BDT_cut_%i.root",year,nprob,bdt_it));
    TTree* t = (TTree*)f->Get("DecayTree");

    double BDT_cut, S, B, sig_bdt_eff, bkg_bdt_eff, sensitivity;
    t->SetBranchAddress("BDT_cut", &BDT_cut);
    t->SetBranchAddress("S", &S);
    t->SetBranchAddress("B", &B);
    t->SetBranchAddress("sig_bdt_eff", &sig_bdt_eff);
    t->SetBranchAddress("bkg_bdt_eff", &bkg_bdt_eff);
    t->SetBranchAddress("sensitivity", &sensitivity);
    t->GetEntry(0);

    std::vector<double> result;
    result.push_back(BDT_cut);
    result.push_back(S);
    result.push_back(B);
    result.push_back(sig_bdt_eff);
    result.push_back(bkg_bdt_eff);
    result.push_back(sensitivity);

    return result;
}