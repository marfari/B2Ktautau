#define year 8
#define pnu0_delta 5
// 0 = 10 GeV (0.01)
// 1 = 1 TeV (0.01)
// 2 = 1 GeV (0.01)
// 3 = 100 TeV (0.01)
// 4 = 10 GeV (0.001)
// 5 = official MC

void bias_DTF()
{
    // TFile* f;
    // if(pnu0_delta == 0){f = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_10GeV.root",year,year));}
    // else if(pnu0_delta == 1){f = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_1TeV.root",year,year));}
    // else if(pnu0_delta == 2){f = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_1GeV.root",year,year));}
    // else if(pnu0_delta == 3){f = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_100TeV.root",year,year));}
    // else if(pnu0_delta == 4){f = new TFile(Form("/panfs/felician/B2Ktautau/ROOT_Sim/201%i/mc_201%i_10GeV_0.001.root",year,year));}

    TFileCollection *fc = new TFileCollection("MC", "MC", Form("Files_on_grid/MC_201%i_all_MagUp.txt",year));
    fc->AddFromFile(Form("Files_on_grid/MC_201%i_all_MagDown.txt",year));

    //TTree* t = (TTree*)f->Get("ntuple/DecayTree");
    TChain* t = new TChain("ntuple/DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    int Bp_TRUEID, taup_TRUEID, taum_TRUEID, Kp_TRUEID, taup_pi1_TRUEID, taup_pi2_TRUEID, taup_pi3_TRUEID, taum_pi1_TRUEID, taum_pi2_TRUEID, taum_pi3_TRUEID;
    double taup_TRUEP_Z, taup_pi1_TRUEP_Z, taup_pi2_TRUEP_Z, taup_pi3_TRUEP_Z;
    double taum_TRUEP_Z, taum_pi1_TRUEP_Z, taum_pi2_TRUEP_Z, taum_pi3_TRUEP_Z;
    float status_0, status_1, status_2, status_012, status_021, status_102, status_120, status_201, status_210;
    float M_0, M_1, M_2, M_012, M_021, M_102, M_120, M_201, M_210;
    float antinutau_PZ_0, antinutau_PZ_1, antinutau_PZ_2, antinutau_PZ_012, antinutau_PZ_021, antinutau_PZ_102, antinutau_PZ_120, antinutau_PZ_201, antinutau_PZ_210;
    float nutau_PZ_0, nutau_PZ_1, nutau_PZ_2, nutau_PZ_012, nutau_PZ_021, nutau_PZ_102, nutau_PZ_120, nutau_PZ_201, nutau_PZ_210;
    float MERR_0, MERR_1, MERR_2, MERR_012, MERR_021, MERR_102, MERR_120, MERR_201, MERR_210;
    float MTau_0, MTau_1, MTau_2, MTau_012, MTau_021, MTau_102, MTau_120, MTau_201, MTau_210;
    float MTauERR_0, MTauERR_1, MTauERR_2, MTauERR_012, MTauERR_021, MTauERR_102, MTauERR_120, MTauERR_201, MTauERR_210;
    int whichFitter_012, whichFitter_021, whichFitter_102, whichFitter_120, whichFitter_201, whichFitter_210;

    float status[] = {status_0, status_1, status_2, status_012, status_021, status_102, status_120, status_201, status_210};
    float M[] = {M_0, M_1, M_2, M_012, M_021, M_102, M_120, M_201, M_210};
    float antinutau_PZ[] = {antinutau_PZ_0, antinutau_PZ_1, antinutau_PZ_2, antinutau_PZ_012, antinutau_PZ_021, antinutau_PZ_102, antinutau_PZ_120, antinutau_PZ_201, antinutau_PZ_210};
    float nutau_PZ[] = {nutau_PZ_0, nutau_PZ_1, nutau_PZ_2, nutau_PZ_012, nutau_PZ_021, nutau_PZ_102, nutau_PZ_120, nutau_PZ_201, nutau_PZ_210};
    float MERR[] = {MERR_0, MERR_1, MERR_2, MERR_012, MERR_021, MERR_102, MERR_120, MERR_201, MERR_210};
    float MTau[] = {MTau_0, MTau_1, MTau_2, MTau_012, MTau_021, MTau_102, MTau_120, MTau_201, MTau_210};
    float MTauERR[] = {MTauERR_0, MTauERR_1, MTauERR_2, MTauERR_012, MTauERR_021, MTauERR_102, MTauERR_120, MTauERR_201, MTauERR_210};
    int whichFitter[] = {whichFitter_012, whichFitter_021, whichFitter_102, whichFitter_120, whichFitter_201, whichFitter_210};

    double M_PDG = 5279.25; // MeV (PDG average)
    double MTau_PDG = 1776.86; 
    int sequence[] = {0, 1, 2, 12, 21, 102, 120, 201, 210};
    
    t->SetBranchAddress("Bp_TRUEID", &Bp_TRUEID);
    t->SetBranchAddress("taup_TRUEID", &taup_TRUEID);
    t->SetBranchAddress("taum_TRUEID", &taum_TRUEID);
    t->SetBranchAddress("Kp_TRUEID", &Kp_TRUEID);
    t->SetBranchAddress("taup_pi1_TRUEID", &taup_pi1_TRUEID);
    t->SetBranchAddress("taup_pi2_TRUEID", &taup_pi2_TRUEID);
    t->SetBranchAddress("taup_pi3_TRUEID", &taup_pi3_TRUEID);
    t->SetBranchAddress("taum_pi1_TRUEID", &taum_pi1_TRUEID);
    t->SetBranchAddress("taum_pi2_TRUEID", &taum_pi2_TRUEID);
    t->SetBranchAddress("taum_pi3_TRUEID", &taum_pi3_TRUEID);

    t->SetBranchAddress("taup_TRUEP_Z", &taup_TRUEP_Z);
    t->SetBranchAddress("taup_pi1_TRUEP_Z", &taup_pi1_TRUEP_Z);
    t->SetBranchAddress("taup_pi2_TRUEP_Z", &taup_pi2_TRUEP_Z);
    t->SetBranchAddress("taup_pi3_TRUEP_Z", &taup_pi3_TRUEP_Z);
    t->SetBranchAddress("taum_TRUEP_Z", &taum_TRUEP_Z);
    t->SetBranchAddress("taum_pi1_TRUEP_Z", &taum_pi1_TRUEP_Z);
    t->SetBranchAddress("taum_pi2_TRUEP_Z", &taum_pi2_TRUEP_Z);
    t->SetBranchAddress("taum_pi3_TRUEP_Z", &taum_pi3_TRUEP_Z);

    std::vector<TH1D*> h_pz, h_pz_sig, h_mass_pull, h_mass_pull_sig, h_tau_mass_pull, h_tau_mass_pull_sig, h_MERR_0, h_MERR_1, h_MERR_2;
    std::vector<TH2D*> h_corr_bias_PZ;

    for(int i = 0; i < size(sequence); i++)
    {
        int seq = sequence[i];
        if((seq == 0) || (seq == 1) || (seq == 2))
        {
            t->SetBranchAddress(Form("Bp_ConsBp_%i_status", seq), &status[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_%i_M", seq), &M[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_%i_MERR", seq), &MERR[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_%i_tauminus_nu_tau_PZ", seq), &antinutau_PZ[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_%i_tauminus_0_nu_tau_PZ", seq), &nutau_PZ[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_%i_tauminus_M", seq), &MTau[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_%i_tauminus_MERR", seq), &MTauERR[i]);
        }
        else
        {
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_status", seq), &status[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_M", seq), &M[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_MERR", seq), &MERR[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_tauminus_nu_tau_PZ", seq), &antinutau_PZ[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_tauminus_0_nu_tau_PZ", seq), &nutau_PZ[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_tauminus_M", seq), &MTau[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_tauminus_MERR", seq), &MTauERR[i]);
            t->SetBranchAddress(Form("Bp_ConsBp_seq_%i_whichFitter", seq), &whichFitter[i-3]);
        }

        h_pz.push_back(new TH1D(Form("h_pz_%i",seq), Form("h_pz_%i",seq), 100, -100000, 100000));
        h_pz_sig.push_back(new TH1D(Form("h_pz_sig_%i",seq), Form("h_pz_sig_%i",seq), 100, -100000, 100000));
        h_mass_pull.push_back(new TH1D(Form("h_mass_pull_%i",seq), Form("h_mass_pull_%i",seq), 100, -10, 10));
        h_mass_pull_sig.push_back(new TH1D(Form("h_mass_pull_sig_%i",seq), Form("h_mass_pull_sig_%i",seq), 100, -10, 10));
        h_tau_mass_pull.push_back(new TH1D(Form("h_tau_mass_pull_%i",seq), Form("h_tau_mass_pull_%i",seq), 100, -5, 5));
        h_tau_mass_pull_sig.push_back(new TH1D(Form("h_tau_mass_pull_sig_%i",seq), Form("h_tau_mass_pull_sig_%i",seq), 100, -5, 5));

        if(!((seq == 0) || (seq == 1) || (seq == 2)))
        {
            h_MERR_0.push_back(new TH1D(Form("h_MERR_0_%i",seq), Form("h_MERR_0_%i",seq), 100, 0, 2000));
            h_MERR_1.push_back(new TH1D(Form("h_MERR_1_%i",seq), Form("h_MERR_1_%i",seq), 100, 0, 2000));
            h_MERR_2.push_back(new TH1D(Form("h_MERR_2_%i",seq), Form("h_MERR_2_%i",seq), 100, 0, 2000));

            h_corr_bias_PZ.push_back(new TH2D(Form("h_corr_bias_PZ_%i",seq), Form("h_corr_bias_PZ_%i",seq), 100, -100000, 100000, 100, -100000, 100000));
        }
    }

    int counter = 0;

    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);

        for(int j = 0; j < size(sequence); j++)
        {
            int seq = sequence[j];

            if((abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521))
            {// truth matched

                if(status[j] == 0)
                {// pass DTF
                    h_pz[j]->Fill(antinutau_PZ[j] - (taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z));
                    h_mass_pull[j]->Fill((M[j] - M_PDG)/(MERR[j]));
                    h_tau_mass_pull[j]->Fill((MTau[j] - MTau_PDG)/(MTauERR[j]));

                    if(!((seq == 0) || (seq == 1) || (seq == 2)))
                    {
                        if(whichFitter[j-3] == 0)
                        {
                            h_MERR_0[j-3]->Fill(MERR[j-3]);
                        }
                        else if(whichFitter[j-3] == 1)
                        {
                            h_MERR_1[j-3]->Fill(MERR[j-3]);
                        }
                        else if(whichFitter[j-3] == 2)
                        {
                            h_MERR_2[j-3]->Fill(MERR[j-3]);
                        }

                        h_corr_bias_PZ[j-3]->Fill(antinutau_PZ[j-3] - (taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z), nutau_PZ[j-3] - (taum_TRUEP_Z - taum_pi1_TRUEP_Z - taum_pi2_TRUEP_Z - taum_pi3_TRUEP_Z)  );
                    }
                    
                    if((M[j] > 4500) && (M[j] < 6000))
                    {// inside signal region
                        h_pz_sig[j]->Fill(antinutau_PZ[j] - (taup_TRUEP_Z - taup_pi1_TRUEP_Z - taup_pi2_TRUEP_Z - taup_pi3_TRUEP_Z));
                        h_mass_pull_sig[j]->Fill((M[j] - M_PDG)/(MERR[j]));
                        h_tau_mass_pull_sig[j]->Fill((MTau[j] - MTau_PDG)/(MTauERR[j]));
                    }
                }
            }
        }

        // status display
        if ( i> 1.0*counter*(t->GetEntries())/100 ) {
            cout<<counter<<"%"<<endl;
            counter += 10;
        }
    }

    for(int i = 0; i < size(sequence); i++)
    {
        int seq = sequence[i];

        TCanvas c1 = new TCanvas();
        c1.cd();
        h_pz[i]->SetFillColorAlpha(kBlue-9, 0.5);
        h_pz[i]->GetXaxis()->SetTitle("antinutau DTF - TRUE PZ (MeV)");
        h_pz[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_pz[i]->GetXaxis()->GetXmax() - h_pz[i]->GetXaxis()->GetXmin())/h_pz[i]->GetNbinsX()) );
        h_pz[i]->SetTitle(Form("Method %i", seq));
        h_pz[i]->Draw();
        if(pnu0_delta == 0)
        {// 10 GeV
            c1.SaveAs(Form("Bias_DTF/201%i/10GeV/antinutau_PZ_%i.gif",year,seq));
            c1.SaveAs(Form("Bias_DTF/201%i/10GeV/antinutau_PZ_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 1)
        {// 1 TeV
            c1.SaveAs(Form("Bias_DTF/201%i/1TeV/antinutau_PZ_%i.gif",year,seq));
            c1.SaveAs(Form("Bias_DTF/201%i/1TeV/antinutau_PZ_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 2)
        {// 1 GeV
            c1.SaveAs(Form("Bias_DTF/201%i/1GeV/antinutau_PZ_%i.gif",year,seq));
            c1.SaveAs(Form("Bias_DTF/201%i/1GeV/antinutau_PZ_%i.pdf",year,seq));       
        }
        else if(pnu0_delta == 3)
        {// 100 TeV
            c1.SaveAs(Form("Bias_DTF/201%i/100TeV/antinutau_PZ_%i.gif",year,seq));
            c1.SaveAs(Form("Bias_DTF/201%i/100TeV/antinutau_PZ_%i.pdf",year,seq));       
        }
        else if(pnu0_delta == 4)
        {
            c1.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/antinutau_PZ_%i.gif",year,seq));
            c1.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/antinutau_PZ_%i.pdf",year,seq));     
        }
        else if(pnu0_delta == 5)
        {
            c1.SaveAs(Form("Bias_DTF/201%i/OfficialMC/antinutau_PZ_%i.gif",year,seq));
            c1.SaveAs(Form("Bias_DTF/201%i/OfficialMC/antinutau_PZ_%i.pdf",year,seq));        
        }
        c1.Update();

        TCanvas c2 = new TCanvas();
        c2.cd();
        h_pz_sig[i]->SetFillColorAlpha(kBlue-9, 0.5);
        h_pz_sig[i]->GetXaxis()->SetTitle("antinutau DTF - TRUE PZ (MeV)");
        h_pz_sig[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_pz_sig[i]->GetXaxis()->GetXmax() - h_pz_sig[i]->GetXaxis()->GetXmin())/h_pz_sig[i]->GetNbinsX()) );
        h_pz_sig[i]->SetTitle(Form("Method %i", seq));
        h_pz_sig[i]->Draw();
        if(pnu0_delta == 0)
        {// 10 GeV
            c2.SaveAs(Form("Bias_DTF/201%i/10GeV/antinutau_PZ_sig_%i.gif",year,seq));
            c2.SaveAs(Form("Bias_DTF/201%i/10GeV/antinutau_PZ_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 1)
        {// 1 TeV
            c2.SaveAs(Form("Bias_DTF/201%i/1TeV/antinutau_PZ_sig_%i.gif",year,seq));
            c2.SaveAs(Form("Bias_DTF/201%i/1TeV/antinutau_PZ_sig_%i.pdf",year,seq));  
        }
        else if(pnu0_delta == 2)
        {// 1 GeV
            c2.SaveAs(Form("Bias_DTF/201%i/1GeV/antinutau_PZ_sig_%i.gif",year,seq));
            c2.SaveAs(Form("Bias_DTF/201%i/1GeV/antinutau_PZ_sig_%i.pdf",year,seq));     
        }
        else if(pnu0_delta == 3)
        {// 100 TeV
            c2.SaveAs(Form("Bias_DTF/201%i/100TeV/antinutau_PZ_sig_%i.gif",year,seq));
            c2.SaveAs(Form("Bias_DTF/201%i/100TeV/antinutau_PZ_sig_%i.pdf",year,seq));      
        }
        else if(pnu0_delta == 4)
        {
            c2.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/antinutau_PZ_sig_%i.gif",year,seq));
            c2.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/antinutau_PZ_sig_%i.pdf",year,seq));     
        }
        else if(pnu0_delta == 5)
        {
            c2.SaveAs(Form("Bias_DTF/201%i/OfficialMC/antinutau_PZ_sig_%i.gif",year,seq));
            c2.SaveAs(Form("Bias_DTF/201%i/OfficialMC/antinutau_PZ_sig_%i.pdf",year,seq));
        }
        c2.Update();

        TCanvas c3 = new TCanvas();
        c3.cd();
        h_mass_pull[i]->SetFillColorAlpha(kBlue-9, 0.5);
        h_mass_pull[i]->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/ #Delta m^{DTF}");
        h_mass_pull[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass_pull[i]->GetXaxis()->GetXmax() - h_mass_pull[i]->GetXaxis()->GetXmin())/h_mass_pull[i]->GetNbinsX()) );
        h_mass_pull[i]->SetTitle(Form("Method %i", seq));
        h_mass_pull[i]->Draw();
        if(pnu0_delta == 0)
        {// 10 GeV
            c3.SaveAs(Form("Bias_DTF/201%i/10GeV/Bp_mass_pull_%i.gif",year,seq));
            c3.SaveAs(Form("Bias_DTF/201%i/10GeV/Bp_mass_pull_%i.pdf",year,seq));     
        }
        else if(pnu0_delta == 1)
        {// 1 TeV
            c3.SaveAs(Form("Bias_DTF/201%i/1TeV/Bp_mass_pull_%i.gif",year,seq));
            c3.SaveAs(Form("Bias_DTF/201%i/1TeV/Bp_mass_pull_%i.pdf",year,seq));     
        }
        else if(pnu0_delta == 2)
        {// 1 GeV
            c3.SaveAs(Form("Bias_DTF/201%i/1GeV/Bp_mass_pull_%i.gif",year,seq));
            c3.SaveAs(Form("Bias_DTF/201%i/1GeV/Bp_mass_pull_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 3)
        {// 100 TeV
            c3.SaveAs(Form("Bias_DTF/201%i/100TeV/Bp_mass_pull_%i.gif",year,seq));
            c3.SaveAs(Form("Bias_DTF/201%i/100TeV/Bp_mass_pull_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 4)
        {
            c3.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Bp_mass_pull_%i.gif",year,seq));
            c3.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Bp_mass_pull_%i.pdf",year,seq)); 
        }
        else if(pnu0_delta == 5)
        {
            c3.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Bp_mass_pull_%i.gif",year,seq));
            c3.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Bp_mass_pull_%i.pdf",year,seq)); 
        }
        c3.Update();

        TCanvas c4 = new TCanvas();
        c4.cd();
        h_mass_pull_sig[i]->SetFillColorAlpha(kBlue-9, 0.5);
        h_mass_pull_sig[i]->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/ #Delta m^{DTF}");
        h_mass_pull_sig[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass_pull_sig[i]->GetXaxis()->GetXmax() - h_mass_pull_sig[i]->GetXaxis()->GetXmin())/h_mass_pull_sig[i]->GetNbinsX()) );
        h_mass_pull_sig[i]->SetTitle(Form("Method %i", seq));
        h_mass_pull_sig[i]->Draw();
        if(pnu0_delta == 0)
        {// 10 GeV
            c4.SaveAs(Form("Bias_DTF/201%i/10GeV/Bp_mass_pull_sig_%i.gif",year,seq));
            c4.SaveAs(Form("Bias_DTF/201%i/10GeV/Bp_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 1)
        {// 1 TeV
            c4.SaveAs(Form("Bias_DTF/201%i/1TeV/Bp_mass_pull_sig_%i.gif",year,seq));
            c4.SaveAs(Form("Bias_DTF/201%i/1TeV/Bp_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 2)
        {// 1 GeV
            c4.SaveAs(Form("Bias_DTF/201%i/1GeV/Bp_mass_pull_sig_%i.gif",year,seq));
            c4.SaveAs(Form("Bias_DTF/201%i/1GeV/Bp_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 3)
        {// 100 TeV
            c4.SaveAs(Form("Bias_DTF/201%i/100TeV/Bp_mass_pull_sig_%i.gif",year,seq));
            c4.SaveAs(Form("Bias_DTF/201%i/100TeV/Bp_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 4)
        {
            c4.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Bp_mass_pull_sig_%i.gif",year,seq));
            c4.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Bp_mass_pull_sig_%i.pdf",year,seq)); 
        }
        else if(pnu0_delta == 5)
        {
            c4.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Bp_mass_pull_sig_%i.gif",year,seq));
            c4.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Bp_mass_pull_sig_%i.pdf",year,seq));  
        }
        c4.Update();

        TCanvas c5 = new TCanvas();
        c5.cd();
        h_tau_mass_pull[i]->SetFillColorAlpha(kBlue-9, 0.5);
        h_tau_mass_pull[i]->GetXaxis()->SetTitle("(m_{#tau}^{DTF} - m_{#tau}^{PDG})/ #Delta m_{#tau}^{DTF}");
        h_tau_mass_pull[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_tau_mass_pull[i]->GetXaxis()->GetXmax() - h_tau_mass_pull[i]->GetXaxis()->GetXmin())/h_tau_mass_pull[i]->GetNbinsX()) );
        h_tau_mass_pull[i]->SetTitle(Form("Method %i", seq));
        h_tau_mass_pull[i]->Draw();
        if(pnu0_delta == 0)
        {// 10 GeV
            c5.SaveAs(Form("Bias_DTF/201%i/10GeV/Taup_mass_pull_%i.gif",year,seq));
            c5.SaveAs(Form("Bias_DTF/201%i/10GeV/Taup_mass_pull_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 1)
        {// 1 TeV
            c5.SaveAs(Form("Bias_DTF/201%i/1TeV/Taup_mass_pull_%i.gif",year,seq));
            c5.SaveAs(Form("Bias_DTF/201%i/1TeV/Taup_mass_pull_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 2)
        {// 1 GeV
            c5.SaveAs(Form("Bias_DTF/201%i/1GeV/Taup_mass_pull_%i.gif",year,seq));
            c5.SaveAs(Form("Bias_DTF/201%i/1GeV/Taup_mass_pull_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 3)
        {// 100 TeV
            c5.SaveAs(Form("Bias_DTF/201%i/100TeV/Taup_mass_pull_%i.gif",year,seq));
            c5.SaveAs(Form("Bias_DTF/201%i/100TeV/Taup_mass_pull_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 4)
        {
            c5.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Taup_mass_pull_%i.gif",year,seq));
            c5.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Taup_mass_pull_%i.pdf",year,seq));  
        }
        else if(pnu0_delta == 5)
        {
            c5.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Taup_mass_pull_%i.gif",year,seq));
            c5.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Taup_mass_pull_%i.pdf",year,seq));  
        }
        c5.Update();

        TCanvas c6 = new TCanvas();
        c6.cd();
        h_tau_mass_pull_sig[i]->SetFillColorAlpha(kBlue-9, 0.5);
        h_tau_mass_pull_sig[i]->GetXaxis()->SetTitle("(m_{#tau}^{DTF} - m_{#tau}^{PDG})/ #Delta m_{#tau}^{DTF}");
        h_tau_mass_pull_sig[i]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_tau_mass_pull_sig[i]->GetXaxis()->GetXmax() - h_tau_mass_pull_sig[i]->GetXaxis()->GetXmin())/h_tau_mass_pull_sig[i]->GetNbinsX()) );
        h_tau_mass_pull_sig[i]->SetTitle(Form("Method %i", seq));
        h_tau_mass_pull_sig[i]->Draw();
        if(pnu0_delta == 0)
        {// 10 GeV
            c6.SaveAs(Form("Bias_DTF/201%i/10GeV/Taup_mass_pull_sig_%i.gif",year,seq));
            c6.SaveAs(Form("Bias_DTF/201%i/10GeV/Taup_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 1)
        {// 1 TeV
            c6.SaveAs(Form("Bias_DTF/201%i/1TeV/Taup_mass_pull_sig_%i.gif",year,seq));
            c6.SaveAs(Form("Bias_DTF/201%i/1TeV/Taup_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 2)
        {// 1 GeV
            c6.SaveAs(Form("Bias_DTF/201%i/1GeV/Taup_mass_pull_sig_%i.gif",year,seq));
            c6.SaveAs(Form("Bias_DTF/201%i/1GeV/Taup_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 3)
        {// 100 TeV
            c6.SaveAs(Form("Bias_DTF/201%i/100TeV/Taup_mass_pull_sig_%i.gif",year,seq));
            c6.SaveAs(Form("Bias_DTF/201%i/100TeV/Taup_mass_pull_sig_%i.pdf",year,seq));
        }
        else if(pnu0_delta == 4)
        {
            c6.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Taup_mass_pull_sig_%i.gif",year,seq));
            c6.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Taup_mass_pull_sig_%i.pdf",year,seq));   
        }
        else if(pnu0_delta == 5)
        {
            c6.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Taup_mass_pull_sig_%i.gif",year,seq));
            c6.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Taup_mass_pull_sig_%i.pdf",year,seq));     
        }
        c6.Update();

        if(!((seq == 0) || (seq == 1) || (seq == 2)))
        {
            TCanvas c7 = new TCanvas();
            c7.cd();
            h_MERR_0[i-3]->SetLineColor(kBlack);
            h_MERR_1[i-3]->SetLineColor(kBlue);
            h_MERR_2[i-3]->SetLineColor(kRed);

            h_MERR_0[i-3]->SetFillColorAlpha(kBlack,0.25);
            h_MERR_1[i-3]->SetFillColorAlpha(kBlue, 0.25);
            h_MERR_2[i-3]->SetFillColorAlpha(kRed, 0.25);

            if(seq == 12)
            {
                h_MERR_2[i-3]->GetXaxis()->SetTitle("DTF B^{+} MERR (MeV)");
                h_MERR_2[i-3]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_MERR_0[i-3]->GetXaxis()->GetXmax() - h_MERR_0[i-3]->GetXaxis()->GetXmin())/h_MERR_0[i-3]->GetNbinsX()) );
                h_MERR_2[i-3]->SetTitle("Method 012"); 

                h_MERR_2[i-3]->DrawNormalized();
                h_MERR_1[i-3]->DrawNormalized("same");
                h_MERR_0[i-3]->DrawNormalized("same");   
            }
            else if(seq == 21)
            {
                h_MERR_1[i-3]->GetXaxis()->SetTitle("DTF B^{+} MERR (MeV)");
                h_MERR_1[i-3]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_MERR_0[i-3]->GetXaxis()->GetXmax() - h_MERR_0[i-3]->GetXaxis()->GetXmin())/h_MERR_0[i-3]->GetNbinsX()) );
                h_MERR_1[i-3]->SetTitle("Method 021"); 

                h_MERR_1[i-3]->DrawNormalized();
                h_MERR_2[i-3]->DrawNormalized("same");
                h_MERR_0[i-3]->DrawNormalized("same"); 
            }
            else if(seq == 102)
            {
                h_MERR_2[i-3]->GetXaxis()->SetTitle("DTF B^{+} MERR (MeV)");
                h_MERR_2[i-3]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_MERR_0[i-3]->GetXaxis()->GetXmax() - h_MERR_0[i-3]->GetXaxis()->GetXmin())/h_MERR_0[i-3]->GetNbinsX()) );
                h_MERR_2[i-3]->SetTitle("Method 102"); 

                h_MERR_2[i-3]->DrawNormalized();
                h_MERR_0[i-3]->DrawNormalized("same");
                h_MERR_1[i-3]->DrawNormalized("same"); 
            }
            else if(seq == 120)
            {
                h_MERR_0[i-3]->GetXaxis()->SetTitle("DTF B^{+} MERR (MeV)");
                h_MERR_0[i-3]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_MERR_0[i-3]->GetXaxis()->GetXmax() - h_MERR_0[i-3]->GetXaxis()->GetXmin())/h_MERR_0[i-3]->GetNbinsX()) );
                h_MERR_0[i-3]->SetTitle("Method 120"); 

                h_MERR_0[i-3]->DrawNormalized();
                h_MERR_2[i-3]->DrawNormalized("same");
                h_MERR_1[i-3]->DrawNormalized("same"); 
            }
            else if(seq == 201)
            {
                h_MERR_1[i-3]->GetXaxis()->SetTitle("DTF B^{+} MERR (MeV)");
                h_MERR_1[i-3]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_MERR_0[i-3]->GetXaxis()->GetXmax() - h_MERR_0[i-3]->GetXaxis()->GetXmin())/h_MERR_0[i-3]->GetNbinsX()) );
                h_MERR_1[i-3]->SetTitle("Method 201"); 

                h_MERR_1[i-3]->DrawNormalized();
                h_MERR_0[i-3]->DrawNormalized("same");
                h_MERR_2[i-3]->DrawNormalized("same"); 
            }
            else if(seq == 210)
            {
                h_MERR_0[i-3]->GetXaxis()->SetTitle("DTF B^{+} MERR (MeV)");
                h_MERR_0[i-3]->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_MERR_0[i-3]->GetXaxis()->GetXmax() - h_MERR_0[i-3]->GetXaxis()->GetXmin())/h_MERR_0[i-3]->GetNbinsX()) );
                h_MERR_0[i-3]->SetTitle("Method 210"); 

                h_MERR_0[i-3]->DrawNormalized();
                h_MERR_1[i-3]->DrawNormalized("same");
                h_MERR_2[i-3]->DrawNormalized("same"); 
            }

            TLegend* leg = new TLegend(0.6, 0.6, 0.75, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->AddEntry(h_MERR_0[i-3],"Strategy 0","f");
            leg->AddEntry(h_MERR_1[i-3],"Strategy 1","f");
            leg->AddEntry(h_MERR_2[i-3],"Strategy 2","f");
            leg->Draw("same");

            if(pnu0_delta == 0)
            {// 10 GeV
                c7.SaveAs(Form("Bias_DTF/201%i/10GeV/Bp_MERR_whichFitter_%i.gif",year,seq));
                c7.SaveAs(Form("Bias_DTF/201%i/10GeV/Bp_MERR_whichFitter_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 1)
            {// 1 TeV
                c7.SaveAs(Form("Bias_DTF/201%i/1TeV/Bp_MERR_whichFitter_%i.gif",year,seq));
                c7.SaveAs(Form("Bias_DTF/201%i/1TeV/Bp_MERR_whichFitter_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 2)
            {// 1 GeV
                c7.SaveAs(Form("Bias_DTF/201%i/1GeV/Bp_MERR_whichFitter_%i.gif",year,seq));
                c7.SaveAs(Form("Bias_DTF/201%i/1GeV/Bp_MERR_whichFitter_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 3)
            {// 100 TeV
                c7.SaveAs(Form("Bias_DTF/201%i/100TeV/Bp_MERR_whichFitter_%i.gif",year,seq));
                c7.SaveAs(Form("Bias_DTF/201%i/100TeV/Bp_MERR_whichFitter_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 4)
            {
                c7.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Bp_MERR_whichFitter_%i.gif",year,seq));
                c7.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/Bp_MERR_whichFitter_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 5)
            {
                c7.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Bp_MERR_whichFitter_%i.gif",year,seq));
                c7.SaveAs(Form("Bias_DTF/201%i/OfficialMC/Bp_MERR_whichFitter_%i.pdf",year,seq));  
            }
            c7.Update();

            TCanvas c8;
            c8.cd();
            h_corr_bias_PZ[i-3]->GetXaxis()->SetTitle("antinutau DTF - TRUE PZ (MeV)");
            h_corr_bias_PZ[i-3]->GetYaxis()->SetTitle("nutau DTF - TRUE PZ (MeV)");
            h_corr_bias_PZ[i-3]->SetTitle(Form("Method %i", seq));
            h_corr_bias_PZ[i-3]->Draw("COL");
            if(pnu0_delta == 0)
            {// 10 GeV
                c8.SaveAs(Form("Bias_DTF/201%i/10GeV/antinutau_nutau_bias_PZ_%i.gif",year,seq));
                c8.SaveAs(Form("Bias_DTF/201%i/10GeV/antinutau_nutau_bias_PZ_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 1)
            {// 1 TeV
                c8.SaveAs(Form("Bias_DTF/201%i/1TeV/antinutau_nutau_bias_PZ_%i.gif",year,seq));
                c8.SaveAs(Form("Bias_DTF/201%i/1TeV/antinutau_nutau_bias_PZ_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 2)
            {// 1 GeV
                c8.SaveAs(Form("Bias_DTF/201%i/1GeV/antinutau_nutau_bias_PZ_%i.gif",year,seq));
                c8.SaveAs(Form("Bias_DTF/201%i/1GeV/antinutau_nutau_bias_PZ_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 3)
            {// 100 TeV
                c8.SaveAs(Form("Bias_DTF/201%i/100TeV/antinutau_nutau_bias_PZ_%i.gif",year,seq));
                c8.SaveAs(Form("Bias_DTF/201%i/100TeV/antinutau_nutau_bias_PZ_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 4)
            {
                c8.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/antinutau_nutau_bias_PZ_%i.gif",year,seq));
                c8.SaveAs(Form("Bias_DTF/201%i/10GeV_0.001/antinutau_nutau_bias_PZ_%i.pdf",year,seq));
            }
            else if(pnu0_delta == 5)
            {
                c8.SaveAs(Form("Bias_DTF/201%i/OfficialMC/antinutau_nutau_bias_PZ_%i.gif",year,seq));
                c8.SaveAs(Form("Bias_DTF/201%i/OfficialMC/antinutau_nutau_bias_PZ_%i.pdf",year,seq));   
            }
            c8.Update();
        }
    }
}