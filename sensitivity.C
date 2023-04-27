using namespace std;

vector<double> range(double min, double max, size_t N);

void sensitivity(int year, int nprob, int bdt_it, TString MC_files, TString RS_DATA_files, TString MC_BDT_files, TString RS_DATA_BDT_files, TString BKG_YIELD_file, TString LUMI_file)
{
    vector<double> bdts = range(-1, 1, nprob);
    double bdt_cut = bdts[bdt_it];

    // all these files have already passed the pre-selection cuts
    TFileCollection* fc_mc = new TFileCollection("MC", "MC", MC_files);
    TFileCollection* fc_rs_data = new TFileCollection("RS_DATA", "RS_DATA", RS_DATA_files);
    TFileCollection* fc_mc_bdt = new TFileCollection("MC_BDT", "MC_BDT", MC_BDT_files);
    TFileCollection* fc_rs_data_bdt = new TFileCollection("RS_DATA_BDT", "RS_DATA_BDT", RS_DATA_BDT_files);

    TChain* t_mc = new TChain("DecayTree");
    TChain* t_rs_data = new TChain("DecayTree");
    TChain* t_mc_bdt = new TChain("DecayTree");
    TChain* t_rs_data_bdt = new TChain("DecayTree");

    t_mc->AddFileInfoList((TCollection*)fc_mc->GetList());
    t_rs_data->AddFileInfoList((TCollection*)fc_rs_data->GetList());
    t_mc_bdt->AddFileInfoList((TCollection*)fc_mc_bdt->GetList());
    t_rs_data_bdt->AddFileInfoList((TCollection*)fc_rs_data_bdt->GetList());

    t_mc->AddFriend(t_mc_bdt);
    t_rs_data->AddFriend(t_rs_data_bdt);

    // External parameters
    //double sigma = 86.6*pow(10,-6); // cross-section (fixed)
    double sigma = 2*0.412*560*pow(10,-6);
    double B_3pi_nu = 9.31*pow(10,-2); // from PDG (fixed)
    double B_3pi_pi0_nu = 4.62*pow(10,-2); // from PDF (fixed)

    TFile* fc_lumi = new TFile(LUMI_file);
    TTree* t_lumi = (TTree*)fc_lumi->Get("DecayTree");

    double L;
    t_lumi->SetBranchAddress("lumi_in_subset", &L);
    t_lumi->GetEntry(0);
    L *= pow(10,15);
    cout << "lumi = " << L << endl;

    TFile* fc_yield = new TFile(BKG_YIELD_file);
    TTree* t_yield = (TTree*)fc_yield->Get("DecayTree");

    Int_t B_pre_sel;
    t_yield->SetBranchAddress("RS_signal", &B_pre_sel);
    t_yield->GetEntry(0);
    cout << "RS signal yield = " << B_pre_sel << endl;

    // Efficiencies
    // Acceptance efficiency (MC)
    double sig_acc_eff;
    if(year == 6){sig_acc_eff = 5.335*pow(10,-2);}
    else if(year == 7){sig_acc_eff = 5.315*pow(10,-2);}
    else if(year == 8){sig_acc_eff = 5.330*pow(10,-2);}

    // Stripping efficiency (MC)
    double sig_strip_eff;
    if(year == 6){sig_strip_eff = 1.256*pow(10,-2);}
    else if(year == 7){sig_strip_eff = 1.283*pow(10,-2);}
    else if(year == 8){sig_strip_eff = 1.028*pow(10,-2);}

    // Pre-selection efficiency (MC)
    double sig_presel_eff;
    if(year == 7){sig_presel_eff = 17.2*pow(10,-2);}
    else if(year == 8){sig_presel_eff = 17*pow(10,-2);}

    // BDT efficiency (MC)
    double sig_bdt_eff;
    // BDT efficiency (RS data)
    double bkg_bdt_eff;

    // vector<double> bdt_cut = range(-1, 1, nprob);
    // for(int i = 0; i < nprob; i++)
    // {
    //     cout << bdt_cut[i] << endl;
    // }

    double sig_num, sig_den, bkg_num, bkg_den;
    double B, S, epsilon_S, sensitivity;

    double bdt_values;

    sig_den = t_mc->GetEntries();
    bkg_den = t_rs_data->GetEntries("(Bp_ConsBp_seq_12_M[0] > 5800) && (Bp_ConsBp_seq_12_M[0] < 8000)"); // RS data, right sideband

    // BDT efficiencies
    TCut BDT_cut = Form("(BDT_response > %f)", bdt_cut);

    // Signal BDT efficiency
    sig_num = t_mc->GetEntries(BDT_cut);
    sig_bdt_eff = sig_num/sig_den;

    // Background BDT efficiency 
    bkg_num = t_rs_data->GetEntries("(Bp_ConsBp_seq_12_M[0] > 5800) && (Bp_ConsBp_seq_12_M[0] < 8000)"+BDT_cut);
    bkg_bdt_eff = bkg_num/bkg_den;

    B = B_pre_sel*bkg_bdt_eff;
    S = 0.5*( pow(1.64,2) + sqrt(pow(1.64,4) + 4*pow(1.64,2)*B) ); // 95% CL upper limit -> 1.64 sigma

    //epsilon_S = sig_strip_eff*sig_presel_eff*sig_bdt_eff;
    epsilon_S = sig_acc_eff*sig_strip_eff*sig_presel_eff*sig_bdt_eff;

    if(epsilon_S != 0)
    {
        sensitivity = (1/(L*sigma))*(S/(epsilon_S*pow(B_3pi_nu + B_3pi_pi0_nu,2))); 
    }
    else
    {
        sensitivity = 100000000; // 1/0 indetermination
    }

    cout << "BDT > " << bdt_cut << endl;
    cout << "Signal yield = " << S << endl;
    cout << "Background yield = " << B << endl;
    cout << "Signal efficiency = " << sig_bdt_eff << endl;
    cout << "Background efficiency = " << bkg_bdt_eff << endl;
    cout << "Sensitivity = " << sensitivity << endl;

    TFile* fout = new TFile( Form("/panfs/felician/B2Ktautau/workflow/sensitivity/201%i/N_prob_%i/BDT_cut_%i.root",year,nprob,bdt_it), "RECREATE");
    fout->cd();           
    TTree* tout = new TTree("DecayTree", "DecayTree");
    tout->Branch("BDT_cut", &bdt_cut);
    tout->Branch("S", &S);
    tout->Branch("B", &B);
    tout->Branch("sig_bdt_eff", &sig_bdt_eff);
    tout->Branch("bkg_bdt_eff", &bkg_bdt_eff);
    tout->Branch("sensitivity", &sensitivity);
    tout->Fill();
    fout->Write();
    fout->Close(); 

    return;
}

vector<double> range(double min, double max, size_t N) {
    vector<double> range;
    double delta = (max-min)/double(N-1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}