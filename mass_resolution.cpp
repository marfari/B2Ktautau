using namespace std;

void mass_resolution(){

    TFile* fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");
    TTree* t = (TTree*)fin->Get("ntuple/DecayTree");

    Double_t taup_FD, taum_FD, taup_FD_chi2, taum_FD_chi2, taup_IP_chi2, taum_IP_chi2;
    Float_t status, Bmass, Pnu1x, Pnu1y, Pnu1z, Pnu2x, Pnu2y, Pnu2z;
    Int_t Kp_TRUEID, taup_pip0_TRUEID, taup_pim0_TRUEID, taup_pip1_TRUEID, taum_pim0_TRUEID, taum_pip0_TRUEID, taum_pim1_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;
    Double_t DV1x, DV1y, DV1z, DV2x, DV2y, DV2z, BVx, BVy, BVz;
    Double_t taup_pip0_TRUEPZ, taup_pim0_TRUEPZ, taup_pip1_TRUEPZ, taum_pim0_TRUEPZ, taum_pip0_TRUEPZ, taum_pim1_TRUEPZ, taup_TRUEPZ, taum_TRUEPZ;

    t->SetBranchAddress("taup_TRUEP_Z", &taup_TRUEPZ);
    t->SetBranchAddress("taum_TRUEP_Z", &taum_TRUEPZ);
    t->SetBranchAddress("taup_pip0_TRUEP_Z", &taup_pip0_TRUEPZ);
    t->SetBranchAddress("taup_pim0_TRUEP_Z", &taup_pim0_TRUEPZ);
    t->SetBranchAddress("taup_pip1_TRUEP_Z", &taup_pip1_TRUEPZ);
    t->SetBranchAddress("taum_pim0_TRUEP_Z", &taum_pim0_TRUEPZ);
    t->SetBranchAddress("taum_pip0_TRUEP_Z", &taum_pip0_TRUEPZ);
    t->SetBranchAddress("taum_pim1_TRUEP_Z", &taum_pim1_TRUEPZ);    

    t->SetBranchAddress("Bp_ConsBp_M",&Bmass);
    t->SetBranchAddress("taup_FD_ORIVX",&taup_FD);
    t->SetBranchAddress("taum_FD_ORIVX",&taum_FD); 
    t->SetBranchAddress("taup_FDCHI2_ORIVX",&taup_FD_chi2);
    t->SetBranchAddress("taum_FDCHI2_ORIVX",&taum_FD_chi2);
    t->SetBranchAddress("taup_IPCHI2_OWNPV",&taup_IP_chi2);
    t->SetBranchAddress("taum_IPCHI2_OWNPV",&taum_IP_chi2);

    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PX",&Pnu1x);
    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PY",&Pnu1y);
    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PZ",&Pnu1z);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PX",&Pnu2x);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PY",&Pnu2y);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PZ",&Pnu2z);

    t->SetBranchAddress("Kp_TRUEID", &Kp_TRUEID);
    t->SetBranchAddress("Bp_TRUEID", &Bp_TRUEID);
    t->SetBranchAddress("taup_TRUEID", &taup_TRUEID);
    t->SetBranchAddress("taum_TRUEID", &taum_TRUEID);
    t->SetBranchAddress("taup_pip0_TRUEID", &taup_pip0_TRUEID);
    t->SetBranchAddress("taup_pim0_TRUEID", &taup_pim0_TRUEID);
    t->SetBranchAddress("taup_pip1_TRUEID", &taup_pip1_TRUEID);
    t->SetBranchAddress("taum_pim0_TRUEID", &taum_pim0_TRUEID);
    t->SetBranchAddress("taum_pip0_TRUEID", &taum_pip0_TRUEID);
    t->SetBranchAddress("taum_pim1_TRUEID", &taum_pim1_TRUEID);
    t->SetBranchAddress("Bp_ConsBp_status",&status);

    t->SetBranchAddress("taup_ENDVERTEX_X",&DV1x);
    t->SetBranchAddress("taup_ENDVERTEX_Y",&DV1y);
    t->SetBranchAddress("taup_ENDVERTEX_Z",&DV1z);
    t->SetBranchAddress("taum_ENDVERTEX_X",&DV2x);
    t->SetBranchAddress("taum_ENDVERTEX_Y",&DV2y);
    t->SetBranchAddress("taum_ENDVERTEX_Z",&DV2z);

    t->SetBranchAddress("Bp_ENDVERTEX_X",&BVx);
    t->SetBranchAddress("Bp_ENDVERTEX_Y",&BVy);
    t->SetBranchAddress("Bp_ENDVERTEX_Z",&BVz);

    TH1D* histo = new TH1D("histo", "histo", 80, 4000., 8000.);
    TH1D* histo1 = new TH1D("histo1", "histo1", 80, 4000., 8000.);
    TH1D* histo2 = new TH1D("histo2", "histo2", 80, 4000., 8000.);
    TH1D* histo3 = new TH1D("histo3", "histo3", 80, 4000., 8000.);

    TH1D* h_nutau = new TH1D("h_nutau", "h_nutau", 80, -120000., 120000);
    TH1D* h_antinutau = new TH1D("h_antinutau", "h_antinutau", 80, -120000., 120000);

    double P1;
    double P2;

    double nutau_TRUEPZ;
    double antinutau_TRUEPZ;

    for(int i = 0; i<t->GetEntries(); i++){
        t->GetEntry(i);

        P1 = sqrt(pow(Pnu1x,2) + pow(Pnu1y,2) + pow(Pnu1z,2));
        P2 = sqrt(pow(Pnu2x,2) + pow(Pnu2y,2) + pow(Pnu2z,2));

        nutau_TRUEPZ = taum_TRUEPZ - (taum_pim0_TRUEPZ + taum_pip0_TRUEPZ + taum_pim1_TRUEPZ);
        antinutau_TRUEPZ = taup_TRUEPZ - (taup_pip0_TRUEPZ + taup_pim0_TRUEPZ + taup_pip1_TRUEPZ);

        if( (status == 0) && ((abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521))){
            histo->Fill(Bmass);

            h_nutau->Fill(Pnu2z - nutau_TRUEPZ);
            h_antinutau->Fill(Pnu1z - antinutau_TRUEPZ);

            if(((P1 < 5000) && (P2 < 5000))){
                histo1->Fill(Bmass);
            }
            if((taup_FD_chi2 > 10) && (taum_FD_chi2 > 10)){
                histo2->Fill(Bmass);
            }
            if((taup_FD > 2.) && (taum_FD > 2.)){
                histo3->Fill(Bmass);
            }

        }
    }

    TLegend* leg;
    leg = new TLegend(0.5, 0.7, 0.85, 0.85);
    leg->AddEntry(histo, "No cut", "lp");
    //leg->AddEntry(histo1, "Pnu1,2 < 5000", "lp");
    leg->AddEntry(histo2, "taup,m_FD_chi2 > 10");
    leg->AddEntry(histo3, "taup,m_FD > 2 mm");
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);

    // normalization
    histo->Scale(1/histo->Integral());
    histo1->Scale(1/histo1->Integral());
    histo2->Scale(1/histo2->Integral());
    histo3->Scale(1/histo3->Integral());

    //y axis: maximum and minimum
    double y[] = {histo->GetMaximum(), histo1->GetMaximum(), histo2->GetMaximum(), histo3->GetMaximum()};
    histo->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);
    histo1->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);
    histo2->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);
    histo3->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);

    TCanvas c;
    c.cd();
    gStyle->SetOptStat(0);
    histo->GetXaxis()->SetTitle("m_{B} (MeV)");
    histo->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(histo->GetXaxis()->GetXmax() - histo->GetXaxis()->GetXmin())/histo->GetNbinsX()));
    histo->SetTitle(" ");
    histo->SetLineColor(kBlack);
    histo1->SetLineColor(kBlue);
    histo2->SetLineColor(kRed);
    histo3->SetLineColor(kBlue);
    histo->Draw("HIST");
    //histo1->Draw("HIST same");
    histo2->Draw("HIST same");
    histo3->Draw("HIST same");
    leg->Draw("same");
    c.SaveAs("Mass_resolution_visualisation/mass.gif");
    c.SaveAs("Mass_resolution_visualisation/mass.pdf");


    TCanvas c1;
    c1.cd();
    h_nutau->GetXaxis()->SetTitle("nutau_DTF_PZ - nutau_TRUE_PZ (MeV)");
    h_nutau->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau->GetXaxis()->GetXmax() - h_nutau->GetXaxis()->GetXmin())/h_nutau->GetNbinsX()) );
    h_nutau->SetTitle(" ");
    h_nutau->Draw();
    c1.SaveAs("Mass_resolution_visualisation/nutau_DTF_TRUE.gif");
    c1.SaveAs("Mass_resolution_visualisation/nutau_DTF_TRUE.pdf");

    TCanvas c2;
    c2.cd();
    h_antinutau->GetXaxis()->SetTitle("antinutau_DTF_PZ - antinutau_TRUE_PZ (MeV)");
    h_antinutau->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau->GetXaxis()->GetXmax() - h_antinutau->GetXaxis()->GetXmin())/h_antinutau->GetNbinsX()) );
    h_antinutau->SetTitle(" ");
    h_antinutau->Draw();
    c2.SaveAs("Mass_resolution_visualisation/antinutau_DTF_TRUE.gif");
    c2.SaveAs("Mass_resolution_visualisation/antinutau_DTF_TRUE.pdf");

}