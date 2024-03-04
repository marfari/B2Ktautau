
#define year1 6
#define year2 8

void compare_years()
{

    TFile* f1 = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_truth_matched.root",year1,year1));
    TFile* f2 = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_truth_matched.root",year2,year2));

    TTree* t1 = (TTree*)f1->Get("DecayTree");
    TTree* t2 = (TTree*)f2->Get("DecayTree");

    float mb1, mb2, status1, status2;

    t1->SetBranchAddress("Bp_ConsBp_seq_M", &mb1);
    t1->SetBranchAddress("Bp_ConsBp_seq_status", &status1);
    t2->SetBranchAddress("Bp_ConsBp_seq_M", &mb2);
    t2->SetBranchAddress("Bp_ConsBp_seq_status", &status2);

    TH1D* h1 = new TH1D("h1", "h1", 80, 4000, 8000);
    TH1D* h2 = new TH1D("h2", "h2", 80, 4000, 8000);

    for(int i = 0; i < t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if(status1 == 0)
        {
            h1->Fill(mb1);
        }
    }

    for(int i = 0; i < t2->GetEntries(); i++)
    {
        t2->GetEntry(i);

        if(status2 == 0)
        {
            h2->Fill(mb2);
        }
    }

    TCanvas c;
    c.cd();
    h1->SetLineColor(kBlue);
    h1->SetFillColorAlpha(kBlue, 0.25);
    h2->SetLineColor(kRed);
    h2->SetFillColorAlpha(kRed, 0.25);

    h1->SetTitle(" ");
    h1->GetXaxis()->SetTitle("m_{B^{+}}^{DTF} (MeV)");
    h1->GetYaxis()->SetTitle(Form("Normalised entries (%i bins)",80));
    h2->SetTitle(" ");
    h2->GetXaxis()->SetTitle("m_{B^{+}}^{DTF} (MeV)");
    h2->GetYaxis()->SetTitle(Form("Normalised entries (%i bins)",80));
    

    if( (h1->GetMaximum())/(h1->Integral()) > (h2->GetMaximum())/(h2->Integral()) )
    {
        h1->DrawNormalized();
        h2->DrawNormalized("same");
    }
    else
    {
        h2->DrawNormalized();
        h1->DrawNormalized("same");
    }

    TLegend* leg = new TLegend(0.6, 0.6, 0.75, 0.85);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(h1,"2016","f");
    leg->AddEntry(h2,"2018","f");
    leg->Draw("same");

   int bin1_1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
   int bin2_1 = h1->FindLastBinAbove(h1->GetMaximum()/2);

   int bin1_2 = h2->FindFirstBinAbove(h2->GetMaximum()/2);
   int bin2_2 = h2->FindLastBinAbove(h2->GetMaximum()/2);

    cout << "FWHM 2016 = " << h1->GetBinCenter(bin2_1) - h1->GetBinCenter(bin1_1) << endl;
    cout << "FWHM 2018 = " << h2->GetBinCenter(bin2_2) - h2->GetBinCenter(bin1_2) << endl;

    c.SaveAs(Form("/home/felician/B2Ktautau/Compare_years/201%i_201%i_Bp_ConsBp_seq_M.gif",year1,year2));
    c.SaveAs(Form("/home/felician/B2Ktautau/Compare_years/201%i_201%i_Bp_ConsBp_seq_M.pdf",year1,year2));

}