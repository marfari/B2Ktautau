
void plot(TFile* fout, TH1D* h_mass, std::vector<double> m, TString title, TString xlabel, TString year, int i);
void plot1(TFile* fout, TH1D* h_mass, TString title, TString xlabel, TString year, int i);
void plot_comparison(TFile* fout, TH1D* h_mass, TH1D* h_mass2, TH1D* h_mass3, std::vector<double> m, TString title, TString xlabel, TString year, int i);

TString year = "2016";

bool isData = true;

void backgrounds(){

    TFile* fin;
    if(isData){fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Data/"+year+"/data_"+year+".root");}
    else{fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/"+year+"/mc_"+year+".root");}
    TTree* tin = (TTree*)fin->Get("ntuple/DecayTree");

    TFile* fout;
    if(isData){fout = new TFile("Backgrounds/bkg_"+year+".root","recreate");}
    else{fout = new TFile("Backgrounds/bkg_"+year+"_MC.root","recreate");}

    double Ek, Pkx, Pky, Pkz, Epi11, Ppi11x, Ppi11y, Ppi11z, Epi12, Ppi12x, Ppi12y, Ppi12z, Epi13, Ppi13x, Ppi13y, Ppi13z, Epi21, Ppi21x, Ppi21y, Ppi21z, Epi22, Ppi22x, Ppi22y, Ppi22z, Epi23, Ppi23x, Ppi23y, Ppi23z;
    double taup_FD_chi2, taum_FD_chi2, taup_ENDVERTEX_CHI2, taum_ENDVERTEX_CHI2;
    float status;
    int Kp_TRUEID, taup_pip0_TRUEID, taup_pim0_TRUEID, taup_pip1_TRUEID, taum_pim0_TRUEID, taum_pip0_TRUEID, taum_pim1_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;

    tin->SetBranchAddress("Bp_ConsBp_status",&status);
    tin->SetBranchAddress("Kp_PE",&Ek);
    tin->SetBranchAddress("Kp_PX",&Pkx);
    tin->SetBranchAddress("Kp_PY",&Pky);
    tin->SetBranchAddress("Kp_PZ",&Pkz);
    tin->SetBranchAddress("taup_pim0_PE",&Epi11);
    tin->SetBranchAddress("taup_pim0_PX",&Ppi11x);
    tin->SetBranchAddress("taup_pim0_PY",&Ppi11y);
    tin->SetBranchAddress("taup_pim0_PZ",&Ppi11z);
    tin->SetBranchAddress("taup_pip0_PE",&Epi12);
    tin->SetBranchAddress("taup_pip0_PX",&Ppi12x);
    tin->SetBranchAddress("taup_pip0_PY",&Ppi12y);
    tin->SetBranchAddress("taup_pip0_PZ",&Ppi12z);
    tin->SetBranchAddress("taup_pip1_PE",&Epi13);
    tin->SetBranchAddress("taup_pip1_PX",&Ppi13x);
    tin->SetBranchAddress("taup_pip1_PY",&Ppi13y);
    tin->SetBranchAddress("taup_pip1_PZ",&Ppi13z);
    tin->SetBranchAddress("taum_pim0_PE",&Epi21);
    tin->SetBranchAddress("taum_pim0_PX",&Ppi21x);
    tin->SetBranchAddress("taum_pim0_PY",&Ppi21y);
    tin->SetBranchAddress("taum_pim0_PZ",&Ppi21z);
    tin->SetBranchAddress("taum_pim1_PE",&Epi22);
    tin->SetBranchAddress("taum_pim1_PX",&Ppi22x);
    tin->SetBranchAddress("taum_pim1_PY",&Ppi22y);
    tin->SetBranchAddress("taum_pim1_PZ",&Ppi22z);
    tin->SetBranchAddress("taum_pip0_PE",&Epi23);
    tin->SetBranchAddress("taum_pip0_PX",&Ppi23x);
    tin->SetBranchAddress("taum_pip0_PY",&Ppi23y);
    tin->SetBranchAddress("taum_pip0_PZ",&Ppi23z);
    tin->SetBranchAddress("taup_FDCHI2_ORIVX",&taup_FD_chi2);
    tin->SetBranchAddress("taum_FDCHI2_ORIVX",&taum_FD_chi2);
    tin->SetBranchAddress("taup_ENDVERTEX_CHI2",&taup_ENDVERTEX_CHI2);
    tin->SetBranchAddress("taum_ENDVERTEX_CHI2",&taum_ENDVERTEX_CHI2);

    if(!(isData)){
        tin->SetBranchAddress("Kp_TRUEID", &Kp_TRUEID);
        tin->SetBranchAddress("Bp_TRUEID", &Bp_TRUEID);
        tin->SetBranchAddress("taup_TRUEID", &taup_TRUEID);
        tin->SetBranchAddress("taum_TRUEID", &taum_TRUEID);
        tin->SetBranchAddress("taup_pip0_TRUEID", &taup_pip0_TRUEID);
        tin->SetBranchAddress("taup_pim0_TRUEID", &taup_pim0_TRUEID);
        tin->SetBranchAddress("taup_pip1_TRUEID", &taup_pip1_TRUEID);
        tin->SetBranchAddress("taum_pim0_TRUEID", &taum_pim0_TRUEID);
        tin->SetBranchAddress("taum_pip0_TRUEID", &taum_pip0_TRUEID);
        tin->SetBranchAddress("taum_pim1_TRUEID", &taum_pim1_TRUEID);
    }

    // Masses in MeV
    double m_D0 = 1864.84; 
    double m_D0_err = 0.05; 

    double m_Dpm = 1869.66;
    double m_Dpm_err = 0.05;

    double m_phi = 1019.461;
    double m_phi_err = 0.016;

    double m_Kstar = 895.5;
    double m_Kstar_err = 0.8;

    double m_K = 493.677;
    double m_K_err = 0.016;

    double m_Ds_pm = 1968.35;
    double m_Ds_pm_err = 0.07;

    double m_ro = 775.11; // (tau decays and e+e-)
    double m_ro_err = 0.34; 

    double m_w = 782.66;
    double m_w_err = 0.13;

    double m_Ks = 497.611;
    double m_Ks_err = 0.013;

    double m_Dstar = 2006.85;
    double m_Dstar_err = 0.05;

    double m_Dstar1 = 2010.26;
    double m_Dstar1_err = 0.05;


    // (3pi)+
    TH1D* h_3pip = new TH1D("h_3pip", "h_3pip", 80, 400, 2000);
    TH1D* h_3pim = new TH1D("h_3pim", "h_3pim", 80, 400, 2000);

    // (pi+ pi-)
    TH1D* h_pip_pim1 = new TH1D("h_pip_pim1", "h_pip_pim1", 100, 200, 2000);
    TH1D* h_pip_pim2 = new TH1D("h_pip_pim2", "h_pip_pim2", 100, 200, 2000);
    TH1D* h_pip_pim3 = new TH1D("h_pip_pim3", "h_pip_pim3", 100, 200, 2000);
    TH1D* h_pip_pim4 = new TH1D("h_pip_pim4", "h_pip_pim4", 100, 200, 2000);

    // Kpi
    TH1D* h_Kpi_1 = new TH1D("h_K_pi_1", "h_K_pi_1", 100, 400, 2500);
    TH1D* h_Kpi_2 = new TH1D("h_K_pi_2", "h_K_pi_2", 100, 400, 2500);
    TH1D* h_Kpi_3 = new TH1D("h_K_pi_3", "h_K_pi_3", 100, 400, 2500);

    // K(2pi)+
    TH1D* h_K2pip1 = new TH1D("h_K2pip1", "h_K2pip1", 100, 400, 3000);
    TH1D* h_K2pip2 = new TH1D("h_K2pip2", "h_K2pip2", 100, 400, 3000);
    TH1D* h_K2pip3 = new TH1D("h_K2pip3", "h_K2pip3", 100, 400, 3000);
    
    // 2(pi+pi-)
    TH1D* h_2pip_2pim1 = new TH1D("h_2pip_2pim1", "h_2pip_2pim1", 100, 400, 3000);
    TH1D* h_2pip_2pim2 = new TH1D("h_2pip_2pim2", "h_2pip_2pim2", 100, 400, 3000);
    TH1D* h_2pip_2pim3 = new TH1D("h_2pip_2pim3", "h_2pip_2pim3", 100, 400, 3000);
    TH1D* h_2pip_2pim4 = new TH1D("h_2pip_2pim4", "h_2pip_2pim4", 100, 400, 3000);

    // K(2pi)-pi+
    TH1D* h_K_2pim_pip1 = new TH1D("h_K_2pim_pip1", "h_K_2pim_pip1", 100, 400, 3500);
    TH1D* h_K_2pim_pip2 = new TH1D("h_K_2pim_pip2", "h_K_2pim_pip2", 100, 400, 3500);
    TH1D* h_K_2pim_pip3 = new TH1D("h_K_2pim_pip3", "h_K_2pim_pip3", 100, 400, 3500);
    TH1D* h_K_2pim_pip4 = new TH1D("h_K_2pim_pip4", "h_K_2pim_pip4", 100, 400, 3500);
    TH1D* h_K_2pim_pip5 = new TH1D("h_K_2pim_pip5", "h_K_2pim_pip5", 100, 400, 3500);
    TH1D* h_K_2pim_pip6 = new TH1D("h_K_2pim_pip6", "h_K_2pim_pip6", 100, 400, 3500);
    TH1D* h_K_2pim_pip7 = new TH1D("h_K_2pim_pip7", "h_K_2pim_pip7", 100, 400, 3500);
    TH1D* h_K_2pim_pip8 = new TH1D("h_K_2pim_pip8", "h_K_2pim_pip8", 100, 400, 3500);
    TH1D* h_K_2pim_pip9 = new TH1D("h_K_2pim_pip9", "h_K_2pim_pip9", 100, 400, 3500);

    //Kpi+pi-
    TH1D* h_K_pip_pim1 = new TH1D("h_K_pip_pim1", "h_K_pip_pim1", 100, 400, 3500);
    TH1D* h_K_pip_pim2 = new TH1D("h_K_pip_pim2", "h_K_pip_pim2", 100, 400, 3500);
    TH1D* h_K_pip_pim3 = new TH1D("h_K_pip_pim3", "h_K_pip_pim3", 100, 400, 3500);
    TH1D* h_K_pip_pim4 = new TH1D("h_K_pip_pim4", "h_K_pip_pim4", 100, 400, 3500);
    TH1D* h_K_pip_pim5 = new TH1D("h_K_pip_pim5", "h_K_pip_pim5", 100, 400, 3500);
    TH1D* h_K_pip_pim6 = new TH1D("h_K_pip_pim6", "h_K_pip_pim6", 100, 400, 3500);
    TH1D* h_K_pip_pim7 = new TH1D("h_K_pip_pim7", "h_K_pip_pim7", 100, 400, 3500);
    TH1D* h_K_pip_pim8 = new TH1D("h_K_pip_pim8", "h_K_pip_pim8", 100, 400, 3500);
    TH1D* h_K_pip_pim9 = new TH1D("h_K_pip_pim9", "h_K_pip_pim9", 100, 400, 3500);

    // K(2pi-)(pi+pi-)
    TH1D* h_K_2pim_pip_pim1 = new TH1D("h_K_2pim_pip_pim1", "h_K_2pim_pip_pim1", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim2 = new TH1D("h_K_2pim_pip_pim2", "h_K_2pim_pip_pim2", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim3 = new TH1D("h_K_2pim_pip_pim3", "h_K_2pim_pip_pim3", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim4 = new TH1D("h_K_2pim_pip_pim4", "h_K_2pim_pip_pim4", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim5 = new TH1D("h_K_2pim_pip_pim5", "h_K_2pim_pip_pim5", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim6 = new TH1D("h_K_2pim_pip_pim6", "h_K_2pim_pip_pim6", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim7 = new TH1D("h_K_2pim_pip_pim7", "h_K_2pim_pip_pim7", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim8 = new TH1D("h_K_2pim_pip_pim8", "h_K_2pim_pip_pim8", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim9 = new TH1D("h_K_2pim_pip_pim9", "h_K_2pim_pip_pim9", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim10 = new TH1D("h_K_2pim_pip_pim10", "h_K_2pim_pip_pim10", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim11 = new TH1D("h_K_2pim_pip_pim11", "h_K_2pim_pip_pim11", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim12 = new TH1D("h_K_2pim_pip_pim12", "h_K_2pim_pip_pim12", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim13 = new TH1D("h_K_2pim_pip_pim13", "h_K_2pim_pip_pim13", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim14 = new TH1D("h_K_2pim_pip_pim14", "h_K_2pim_pip_pim14", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim15 = new TH1D("h_K_2pim_pip_pim15", "h_K_2pim_pip_pim15", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim16 = new TH1D("h_K_2pim_pip_pim16", "h_K_2pim_pip_pim16", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim17 = new TH1D("h_K_2pim_pip_pim17", "h_K_2pim_pip_pim17", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim18 = new TH1D("h_K_2pim_pip_pim18", "h_K_2pim_pip_pim18", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim19 = new TH1D("h_K_2pim_pip_pim19", "h_K_2pim_pip_pim19", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim20 = new TH1D("h_K_2pim_pip_pim20", "h_K_2pim_pip_pim20", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim21 = new TH1D("h_K_2pim_pip_pim21", "h_K_2pim_pip_pim21", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim22 = new TH1D("h_K_2pim_pip_pim22", "h_K_2pim_pip_pim22", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim23 = new TH1D("h_K_2pim_pip_pim23", "h_K_2pim_pip_pim23", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim24 = new TH1D("h_K_2pim_pip_pim24", "h_K_2pim_pip_pim24", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim25 = new TH1D("h_K_2pim_pip_pim25", "h_K_2pim_pip_pim25", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim26 = new TH1D("h_K_2pim_pip_pim26", "h_K_2pim_pip_pim26", 100, 1000, 4500);
    TH1D* h_K_2pim_pip_pim27 = new TH1D("h_K_2pim_pip_pim27", "h_K_2pim_pip_pim27", 100, 1000, 4500);

    // 2(pi+pi-)pi+
    TH1D* h_2pippim_pip1 = new TH1D("h_2pippim_pip1", "h_2pippim_pip1", 100, 1000, 4500);
    TH1D* h_2pippim_pip2 = new TH1D("h_2pippim_pip2", "h_2pippim_pip2", 100, 1000, 4500);
    TH1D* h_2pippim_pip3 = new TH1D("h_2pippim_pip3", "h_2pippim_pip3", 100, 1000, 4500);

    // (2pi)+pi-
    TH1D* h_2pip_pim1 = new TH1D("h_2pip_pim1", "h_2pip_pim1", 100, 400, 3500);
    TH1D* h_2pip_pim2 = new TH1D("h_2pip_pim2", "h_2pip_pim2", 100, 400, 3500);
    TH1D* h_2pip_pim3 = new TH1D("h_2pip_pim3", "h_2pip_pim3", 100, 400, 3500);
    TH1D* h_2pip_pim4 = new TH1D("h_2pip_pim4", "h_2pip_pim4", 100, 400, 3500);
    TH1D* h_2pip_pim5 = new TH1D("h_2pip_pim5", "h_2pip_pim5", 100, 400, 3500);
    TH1D* h_2pip_pim6 = new TH1D("h_2pip_pim6", "h_2pip_pim6", 100, 400, 3500);
    TH1D* h_2pip_pim7 = new TH1D("h_2pip_pim7", "h_2pip_pim7", 100, 400, 3500);
    TH1D* h_2pip_pim8 = new TH1D("h_2pip_pim8", "h_2pip_pim8", 100, 400, 3500);
    TH1D* h_2pip_pim9 = new TH1D("h_2pip_pim9", "h_2pip_pim9", 100, 400, 3500);

    // pi+pi- (misID)
    TH1D* h_pip_pim_ID1 = new TH1D("h_pip_pim_ID1", "h_pip_pim_ID1", 100, 400, 3000);
    TH1D* h_pip_pim_ID2 = new TH1D("h_pip_pim_ID2", "h_pip_pim_ID2", 100, 400, 3000);
    TH1D* h_pip_pim_ID3 = new TH1D("h_pip_pim_ID3", "h_pip_pim_ID3", 100, 400, 3000);
    TH1D* h_pip_pim_ID4 = new TH1D("h_pip_pim_ID4", "h_pip_pim_ID4", 100, 400, 3000);
    TH1D* h_pip_pim_ID5 = new TH1D("h_pip_pim_ID5", "h_pip_pim_ID5", 100, 400, 3000);
    TH1D* h_pip_pim_ID6 = new TH1D("h_pip_pim_ID6", "h_pip_pim_ID6", 100, 400, 3000);
    TH1D* h_pip_pim_ID7 = new TH1D("h_pip_pim_ID7", "h_pip_pim_ID7", 100, 400, 3000);
    TH1D* h_pip_pim_ID8 = new TH1D("h_pip_pim_ID8", "h_pip_pim_ID8", 100, 400, 3000);
    TH1D* h_pip_pim_ID9 = new TH1D("h_pip_pim_ID9", "h_pip_pim_ID9", 100, 400, 3000);

    // K+(pi+pi-) (misID)
    TH1D* h_K_2pi_ID1 = new TH1D("h_K_2pi_ID1", "h_K_2pi_ID1", 100, 400, 3500);
    TH1D* h_K_2pi_ID2 = new TH1D("h_K_2pi_ID2", "h_K_2pi_ID2", 100, 400, 3500);
    TH1D* h_K_2pi_ID3 = new TH1D("h_K_2pi_ID3", "h_K_2pi_ID3", 100, 400, 3500);
    TH1D* h_K_2pi_ID4 = new TH1D("h_K_2pi_ID4", "h_K_2pi_ID4", 100, 400, 3500);
    TH1D* h_K_2pi_ID5 = new TH1D("h_K_2pi_ID5", "h_K_2pi_ID5", 100, 400, 3500);
    TH1D* h_K_2pi_ID6 = new TH1D("h_K_2pi_ID6", "h_K_2pi_ID6", 100, 400, 3500);
    TH1D* h_K_2pi_ID7 = new TH1D("h_K_2pi_ID7", "h_K_2pi_ID7", 100, 400, 3500);
    TH1D* h_K_2pi_ID8 = new TH1D("h_K_2pi_ID8", "h_K_2pi_ID8", 100, 400, 3500);
    TH1D* h_K_2pi_ID9 = new TH1D("h_K_2pi_ID9", "h_K_2pi_ID9", 100, 400, 3500);

    // pi+2pi- (misID)
    TH1D* h_pip_2pim_ID1 = new TH1D("h_pip_2pim_ID1", "h_pip_2pim_ID1", 100, 400, 3500);
    TH1D* h_pip_2pim_ID2 = new TH1D("h_pip_2pim_ID2", "h_pip_2pim_ID2", 100, 400, 3500);
    TH1D* h_pip_2pim_ID3 = new TH1D("h_pip_2pim_ID3", "h_pip_2pim_ID3", 100, 400, 3500);
    TH1D* h_pip_2pim_ID4 = new TH1D("h_pip_2pim_ID4", "h_pip_2pim_ID4", 100, 400, 3500);
    TH1D* h_pip_2pim_ID5 = new TH1D("h_pip_2pim_ID5", "h_pip_2pim_ID5", 100, 400, 3500);
    TH1D* h_pip_2pim_ID6 = new TH1D("h_pip_2pim_ID6", "h_pip_2pim_ID6", 100, 400, 3500);
    TH1D* h_pip_2pim_ID7 = new TH1D("h_pip_2pim_ID7", "h_pip_2pim_ID7", 100, 400, 3500);
    TH1D* h_pip_2pim_ID8 = new TH1D("h_pip_2pim_ID8", "h_pip_2pim_ID8", 100, 400, 3500);
    TH1D* h_pip_2pim_ID9 = new TH1D("h_pip_2pim_ID9", "h_pip_2pim_ID9", 100, 400, 3500);

    // pi+(pi+pi-) (misiD)
    TH1D* h_pip_2pi_ID1 = new TH1D("h_pip_2pi_ID1", "h_pip_2pi_ID1", 100, 400, 3500);
    TH1D* h_pip_2pi_ID2 = new TH1D("h_pip_2pi_ID2", "h_pip_2pi_ID2", 100, 400, 3500);
    TH1D* h_pip_2pi_ID3 = new TH1D("h_pip_2pi_ID3", "h_pip_2pi_ID3", 100, 400, 3500);
    TH1D* h_pip_2pi_ID4 = new TH1D("h_pip_2pi_ID4", "h_pip_2pi_ID4", 100, 400, 3500);
    TH1D* h_pip_2pi_ID5 = new TH1D("h_pip_2pi_ID5", "h_pip_2pi_ID5", 100, 400, 3500);
    TH1D* h_pip_2pi_ID6 = new TH1D("h_pip_2pi_ID6", "h_pip_2pi_ID6", 100, 400, 3500);
    TH1D* h_pip_2pi_ID7 = new TH1D("h_pip_2pi_ID7", "h_pip_2pi_ID7", 100, 400, 3500);
    TH1D* h_pip_2pi_ID8 = new TH1D("h_pip_2pi_ID8", "h_pip_2pi_ID8", 100, 400, 3500);
    TH1D* h_pip_2pi_ID9 = new TH1D("h_pip_2pi_ID9", "h_pip_2pi_ID9", 100, 400, 3500);

    // 2pi+2pi- (misID)
    TH1D* h_2pip_2pim_ID1 = new TH1D("h_2pip_2pim_ID1", "h_2pip_2pim_ID1", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID2 = new TH1D("h_2pip_2pim_ID2", "h_2pip_2pim_ID2", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID3 = new TH1D("h_2pip_2pim_ID3", "h_2pip_2pim_ID3", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID4 = new TH1D("h_2pip_2pim_ID4", "h_2pip_2pim_ID4", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID5 = new TH1D("h_2pip_2pim_ID5", "h_2pip_2pim_ID5", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID6 = new TH1D("h_2pip_2pim_ID6", "h_2pip_2pim_ID6", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID7 = new TH1D("h_2pip_2pim_ID7", "h_2pip_2pim_ID7", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID8 = new TH1D("h_2pip_2pim_ID8", "h_2pip_2pim_ID8", 100, 400, 3500);
    TH1D* h_2pip_2pim_ID9 = new TH1D("h_2pip_2pim_ID9", "h_2pip_2pim_ID9", 100, 400, 3500);

    // pi+
    TH1D* h_pip1 = new TH1D("h_pip1", "h_pip1", 80, 0, 200);
    TH1D* h_pip2 = new TH1D("h_pip2", "h_pip2", 80, 0, 200);
    TH1D* h_pip3 = new TH1D("h_pip3", "h_pip3", 80, 0, 200);

    // 2pi (same tau)
    TH1D* h_taup_2pi_1 = new TH1D("h_taup_2pi_1", "h_taup_2pi_1", 100, 200, 2000);
    TH1D* h_taup_2pi_2 = new TH1D("h_taup_2pi_2", "h_taup_2pi_2", 100, 200, 2000);
    TH1D* h_taup_2pi_3 = new TH1D("h_taup_2pi_3", "h_taup_2pi_3", 100, 200, 2000);
    TH1D* h_taum_2pi_1 = new TH1D("h_taum_2pi_1", "h_taum_2pi_1", 100, 200, 2000);
    TH1D* h_taum_2pi_2 = new TH1D("h_taum_2pi_2", "h_taum_2pi_2", 100, 200, 2000);
    TH1D* h_taum_2pi_3 = new TH1D("h_taum_2pi_3", "h_taum_2pi_3", 100, 200, 2000);

    // 2pi (different tau)
    TH1D* h_2pi_1 = new TH1D("h_2pi_1", "h_2pi_1", 100, 200, 2000);
    TH1D* h_2pi_2 = new TH1D("h_2pi_2", "h_2pi_2", 100, 200, 2000);
    TH1D* h_2pi_3 = new TH1D("h_2pi_3", "h_2pi_3", 100, 200, 2000);
    TH1D* h_2pi_4 = new TH1D("h_2pi_4", "h_2pi_4", 100, 200, 2000);
    TH1D* h_2pi_5 = new TH1D("h_2pi_5", "h_2pi_5", 100, 200, 2000);
    TH1D* h_2pi_6 = new TH1D("h_2pi_6", "h_2pi_6", 100, 200, 2000);
    TH1D* h_2pi_7 = new TH1D("h_2pi_7", "h_2pi_7", 100, 200, 2000);
    TH1D* h_2pi_8 = new TH1D("h_2pi_8", "h_2pi_8", 100, 200, 2000);
    TH1D* h_2pi_9 = new TH1D("h_2pi_9", "h_2pi_9", 100, 200, 2000);


    for(int i = 0; i < tin->GetEntries(); ++i){
        tin->GetEntry(i);

        if( (isData && (status == 0)) || ( (!(isData)) && (status == 0) && (abs(Kp_TRUEID) == 321) && (abs(taup_pip0_TRUEID) == 211) && (abs(taup_pim0_TRUEID) == 211) && (abs(taup_pip1_TRUEID) == 211) && (abs(taum_pim0_TRUEID) == 211) && (abs(taum_pip0_TRUEID) == 211) && (abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521) ) ){
            ROOT::Math::PxPyPzEVector Pk(Pkx, Pky, Pkz, Ek);
            ROOT::Math::PxPyPzEVector Ppi11(Ppi11x, Ppi11y, Ppi11z, Epi11);
            ROOT::Math::PxPyPzEVector Ppi12(Ppi12x, Ppi12y, Ppi12z, Epi12);
            ROOT::Math::PxPyPzEVector Ppi13(Ppi13x, Ppi13y, Ppi13z, Epi13);
            ROOT::Math::PxPyPzEVector Ppi21(Ppi21x, Ppi21y, Ppi21z, Epi21);
            ROOT::Math::PxPyPzEVector Ppi22(Ppi22x, Ppi22y, Ppi22z, Epi22);
            ROOT::Math::PxPyPzEVector Ppi23(Ppi23x, Ppi23y, Ppi23z, Epi23);

            double m_kaon = 493.677; // pion misidentified as a kaon
            double Epi11_K = sqrt(pow(m_kaon,2) + pow(Ppi11x,2) + pow(Ppi11y,2) + pow(Ppi11z,2));
            double Epi12_K = sqrt(pow(m_kaon,2) + pow(Ppi12x,2) + pow(Ppi12y,2) + pow(Ppi12z,2));
            double Epi13_K = sqrt(pow(m_kaon,2) + pow(Ppi13x,2) + pow(Ppi13y,2) + pow(Ppi13z,2));
            double Epi21_K = sqrt(pow(m_kaon,2) + pow(Ppi21x,2) + pow(Ppi21y,2) + pow(Ppi21z,2));
            double Epi22_K = sqrt(pow(m_kaon,2) + pow(Ppi22x,2) + pow(Ppi22y,2) + pow(Ppi22z,2));
            double Epi23_K = sqrt(pow(m_kaon,2) + pow(Ppi23x,2) + pow(Ppi23y,2) + pow(Ppi23z,2));

            ROOT::Math::PxPyPzEVector Ppi11_K(Ppi11x, Ppi11y, Ppi11z, Epi11_K);
            ROOT::Math::PxPyPzEVector Ppi12_K(Ppi12x, Ppi12y, Ppi12z, Epi12_K);
            ROOT::Math::PxPyPzEVector Ppi13_K(Ppi13x, Ppi13y, Ppi13z, Epi13_K);
            ROOT::Math::PxPyPzEVector Ppi21_K(Ppi21x, Ppi21y, Ppi21z, Epi21_K);
            ROOT::Math::PxPyPzEVector Ppi22_K(Ppi22x, Ppi22y, Ppi22z, Epi22_K);
            ROOT::Math::PxPyPzEVector Ppi23_K(Ppi23x, Ppi23y, Ppi23z, Epi23_K);

            // (3pi)^+
            ROOT::Math::PxPyPzEVector P_3pip = Ppi11 + Ppi12 + Ppi13;
            ROOT::Math::PxPyPzEVector P_3pim = Ppi21 + Ppi22 + Ppi23;

            // (pi+pi-)
            ROOT::Math::PxPyPzEVector P_pip_pim1 = Ppi12 + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_pim2 = Ppi12 + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_pim3 = Ppi13 + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_pim4 = Ppi13 + Ppi22;

            // Kpi
            ROOT::Math::PxPyPzEVector P_Kpi_1 = Pk + Ppi11;
            ROOT::Math::PxPyPzEVector P_Kpi_2 = Pk + Ppi21;
            ROOT::Math::PxPyPzEVector P_Kpi_3 = Pk + Ppi22;

            // K(2pi)+
            ROOT::Math::PxPyPzEVector P_K2pip1 = Pk + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K2pip2 = Pk + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K2pip3 = Pk + Ppi21 + Ppi22;

            // 2(pi+pi-)
            ROOT::Math::PxPyPzEVector P_2pip_2pim1 = Ppi11 + Ppi21 + Ppi12 + Ppi13;
            ROOT::Math::PxPyPzEVector P_2pip_2pim2 = Ppi11 + Ppi21 + Ppi12 + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pip_2pim3 = Ppi11 + Ppi22 + Ppi12 + Ppi13;
            ROOT::Math::PxPyPzEVector P_2pip_2pim4 = Ppi11 + Ppi22 + Ppi12 + Ppi23;

            // K(2pi)-pi+
            ROOT::Math::PxPyPzEVector P_K_2pim_pip1 = Pk + Ppi11 + Ppi21 + Ppi12;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip2 = Pk + Ppi11 + Ppi21 + Ppi13;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip3 = Pk + Ppi11 + Ppi21 + Ppi23;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip4 = Pk + Ppi11 + Ppi22 + Ppi12;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip5 = Pk + Ppi11 + Ppi22 + Ppi13;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip6 = Pk + Ppi11 + Ppi22 + Ppi23;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip7 = Pk + Ppi21 + Ppi22 + Ppi12;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip8 = Pk + Ppi21 + Ppi22 + Ppi13;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip9 = Pk + Ppi21 + Ppi22 + Ppi23;

            // Kpi+pi-
            ROOT::Math::PxPyPzEVector P_K_pip_pim1 = Pk + Ppi12 + Ppi11;
            ROOT::Math::PxPyPzEVector P_K_pip_pim2 = Pk + Ppi12 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_pip_pim3 = Pk + Ppi12 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_pip_pim4 = Pk + Ppi13 + Ppi11;
            ROOT::Math::PxPyPzEVector P_K_pip_pim5 = Pk + Ppi13 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_pip_pim6 = Pk + Ppi13 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_pip_pim7 = Pk + Ppi23 + Ppi11;
            ROOT::Math::PxPyPzEVector P_K_pip_pim8 = Pk + Ppi23 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_pip_pim9 = Pk + Ppi23 + Ppi22;

            // K(2pi)-pi+pi-
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim1 = Pk + Ppi12 + Ppi11 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim2 = Pk + Ppi12 + Ppi21 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim3 = Pk + Ppi12 + Ppi22 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim4 = Pk + Ppi13 + Ppi11 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim5 = Pk + Ppi13 + Ppi21 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim6 = Pk + Ppi13 + Ppi22 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim7 = Pk + Ppi23 + Ppi11 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim8 = Pk + Ppi23 + Ppi21 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim9 = Pk + Ppi23 + Ppi22 + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim10 = Pk + Ppi12 + Ppi11 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim11 = Pk + Ppi12 + Ppi21 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim12 = Pk + Ppi12 + Ppi22 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim13 = Pk + Ppi13 + Ppi11 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim14 = Pk + Ppi13 + Ppi21 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim15 = Pk + Ppi13 + Ppi22 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim16 = Pk + Ppi23 + Ppi11 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim17 = Pk + Ppi23 + Ppi21 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim18 = Pk + Ppi23 + Ppi22 + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim19 = Pk + Ppi12 + Ppi11 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim20 = Pk + Ppi12 + Ppi21 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim21 = Pk + Ppi12 + Ppi22 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim22 = Pk + Ppi13 + Ppi11 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim23 = Pk + Ppi13 + Ppi21 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim24 = Pk + Ppi13 + Ppi22 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim25 = Pk + Ppi23 + Ppi11 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim26 = Pk + Ppi23 + Ppi21 + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pim_pip_pim27 = Pk + Ppi23 + Ppi22 + Ppi21 + Ppi22;

            // 2(pi+pi-)pi+
            ROOT::Math::PxPyPzEVector P_2pippim_pip1 = Ppi11 + Ppi12 + Ppi21 + Ppi23 + Ppi12;
            ROOT::Math::PxPyPzEVector P_2pippim_pip2 = Ppi11 + Ppi12 + Ppi21 + Ppi23 + Ppi13;
            ROOT::Math::PxPyPzEVector P_2pippim_pip3 = Ppi11 + Ppi12 + Ppi21 + Ppi23 + Ppi23;

            // (2pi)+pi-
            ROOT::Math::PxPyPzEVector P_2pip_pim1 = Ppi12 + Ppi13 + Ppi11;
            ROOT::Math::PxPyPzEVector P_2pip_pim2 = Ppi12 + Ppi13 + Ppi21;
            ROOT::Math::PxPyPzEVector P_2pip_pim3 = Ppi12 + Ppi13 + Ppi22;
            ROOT::Math::PxPyPzEVector P_2pip_pim4 = Ppi12 + Ppi23 + Ppi11;
            ROOT::Math::PxPyPzEVector P_2pip_pim5 = Ppi12 + Ppi23 + Ppi21;
            ROOT::Math::PxPyPzEVector P_2pip_pim6 = Ppi12 + Ppi23 + Ppi22;
            ROOT::Math::PxPyPzEVector P_2pip_pim7 = Ppi23 + Ppi13 + Ppi11;
            ROOT::Math::PxPyPzEVector P_2pip_pim8 = Ppi23 + Ppi13 + Ppi21;
            ROOT::Math::PxPyPzEVector P_2pip_pim9 = Ppi23 + Ppi13 + Ppi22;

            // pi+pi- (misID)
            ROOT::Math::PxPyPzEVector P_pip_pim_ID1 = Ppi12_K + Ppi11;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID2 = Ppi12_K + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID3 = Ppi12_K + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID4 = Ppi13_K + Ppi11;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID5 = Ppi13_K + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID6 = Ppi13_K + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID7 = Ppi23_K + Ppi11;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID8 = Ppi23_K + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_pim_ID9 = Ppi23_K + Ppi22;

            // K+(pi+pi-) (misID)
            ROOT::Math::PxPyPzEVector P_K_2pi_ID1 = Pk + Ppi12_K + Ppi11;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID2 = Pk + Ppi12_K + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID3 = Pk + Ppi12_K + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID4 = Pk + Ppi13_K + Ppi11;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID5 = Pk + Ppi13_K + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID6 = Pk + Ppi13_K + Ppi22;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID7 = Pk + Ppi23_K + Ppi11;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID8 = Pk + Ppi23_K + Ppi21;
            ROOT::Math::PxPyPzEVector P_K_2pi_ID9 = Pk + Ppi23_K + Ppi22;    

            // pi+2pi- (misID)
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID1 = Ppi12_K + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID2 = Ppi13_K + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID3 = Ppi23_K + Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID4 = Ppi12_K + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID5 = Ppi13_K + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID6 = Ppi23_K + Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID7 = Ppi12_K + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID8 = Ppi13_K + Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_pip_2pim_ID9 = Ppi23_K + Ppi21 + Ppi22;

            // pi+(pi+pi-) (misID)
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID1 = Ppi11 + Ppi12 + Ppi12_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID2 = Ppi11 + Ppi12 + Ppi13_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID3 = Ppi11 + Ppi12 + Ppi23_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID4 = Ppi11 + Ppi13 + Ppi12_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID5 = Ppi11 + Ppi13 + Ppi13_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID6 = Ppi11 + Ppi13 + Ppi23_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID7 = Ppi11 + Ppi23 + Ppi12_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID8 = Ppi11 + Ppi23 + Ppi13_K;
            ROOT::Math::PxPyPzEVector P_pip_2pi_ID9 = Ppi11 + Ppi23 + Ppi23_K;

            // 2pi+2pi- (misID)
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID1 = Ppi11 + Ppi21 + Ppi12_K + Ppi13;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID2 = Ppi11 + Ppi21 + Ppi12_K + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID3 = Ppi11 + Ppi21 + Ppi13_K + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID4 = Ppi11 + Ppi22 + Ppi12_K + Ppi13;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID5 = Ppi11 + Ppi22 + Ppi12_K + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID6 = Ppi11 + Ppi22 + Ppi13_K + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID7 = Ppi21 + Ppi22 + Ppi12_K + Ppi13;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID8 = Ppi21 + Ppi22 + Ppi12_K + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pip_2pim_ID9 = Ppi21 + Ppi22 + Ppi13_K + Ppi23;

            // pi+
            ROOT::Math::PxPyPzEVector P_pip1 = Ppi11;
            ROOT::Math::PxPyPzEVector P_pip2 = Ppi12;
            ROOT::Math::PxPyPzEVector P_pip3 = Ppi13;

            // 2pi (same tau)
            ROOT::Math::PxPyPzEVector P_taup_2pi_1 = Ppi11 + Ppi12;
            ROOT::Math::PxPyPzEVector P_taup_2pi_2 = Ppi12 + Ppi13;
            ROOT::Math::PxPyPzEVector P_taup_2pi_3 = Ppi11 + Ppi13;
            ROOT::Math::PxPyPzEVector P_taum_2pi_1 = Ppi21 + Ppi22;
            ROOT::Math::PxPyPzEVector P_taum_2pi_2 = Ppi21 + Ppi23;
            ROOT::Math::PxPyPzEVector P_taum_2pi_3 = Ppi22 + Ppi23;

            // 2pi (different tau)
            ROOT::Math::PxPyPzEVector P_2pi_1 = Ppi11 + Ppi21;
            ROOT::Math::PxPyPzEVector P_2pi_2 = Ppi11 + Ppi22;
            ROOT::Math::PxPyPzEVector P_2pi_3 = Ppi11 + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pi_4 = Ppi12 + Ppi21;
            ROOT::Math::PxPyPzEVector P_2pi_5 = Ppi12 + Ppi22;
            ROOT::Math::PxPyPzEVector P_2pi_6 = Ppi12 + Ppi23;
            ROOT::Math::PxPyPzEVector P_2pi_7 = Ppi13 + Ppi21;
            ROOT::Math::PxPyPzEVector P_2pi_8 = Ppi13 + Ppi22;
            ROOT::Math::PxPyPzEVector P_2pi_9 = Ppi13 + Ppi23;

            // (3pi)^+
            double m_3pip = P_3pip.M();
            double m_3pim = P_3pim.M();

            // (pi+pi-)
            double m_pip_pim1 = P_pip_pim1.M();
            double m_pip_pim2 = P_pip_pim2.M();
            double m_pip_pim3 = P_pip_pim3.M();
            double m_pip_pim4 = P_pip_pim4.M();

            // Kpi
            double m_Kpi_1 = P_Kpi_1.M();
            double m_Kpi_2 = P_Kpi_2.M();
            double m_Kpi_3 = P_Kpi_3.M();

            // K(2pi)+
            double m_K2pip1 = P_K2pip1.M();
            double m_K2pip2 = P_K2pip2.M();
            double m_K2pip3 = P_K2pip3.M();

            // 2(pi+pi-)
            double m_2pip_2pim1 = P_2pip_2pim1.M();
            double m_2pip_2pim2 = P_2pip_2pim2.M();
            double m_2pip_2pim3 = P_2pip_2pim3.M();
            double m_2pip_2pim4 = P_2pip_2pim4.M();

            // K(2pi)-pi+
            double m_K_2pim_pip1 = P_K_2pim_pip1.M();
            double m_K_2pim_pip2 = P_K_2pim_pip2.M();
            double m_K_2pim_pip3 = P_K_2pim_pip3.M();
            double m_K_2pim_pip4 = P_K_2pim_pip4.M();
            double m_K_2pim_pip5 = P_K_2pim_pip5.M();
            double m_K_2pim_pip6 = P_K_2pim_pip6.M();
            double m_K_2pim_pip7 = P_K_2pim_pip7.M();
            double m_K_2pim_pip8 = P_K_2pim_pip8.M();
            double m_K_2pim_pip9 = P_K_2pim_pip9.M();

            // Kpi+pi-
            double m_K_pip_pim1 = P_K_pip_pim1.M();
            double m_K_pip_pim2 = P_K_pip_pim2.M();
            double m_K_pip_pim3 = P_K_pip_pim3.M();
            double m_K_pip_pim4 = P_K_pip_pim4.M();
            double m_K_pip_pim5 = P_K_pip_pim5.M();
            double m_K_pip_pim6 = P_K_pip_pim6.M();
            double m_K_pip_pim7 = P_K_pip_pim7.M();
            double m_K_pip_pim8 = P_K_pip_pim8.M();
            double m_K_pip_pim9 = P_K_pip_pim9.M();

            // K(2pi)-pi+pi-
            double m_K_2pim_pip_pim1 = P_K_2pim_pip_pim1.M();
            double m_K_2pim_pip_pim2 = P_K_2pim_pip_pim2.M();
            double m_K_2pim_pip_pim3 = P_K_2pim_pip_pim3.M();
            double m_K_2pim_pip_pim4 = P_K_2pim_pip_pim4.M();
            double m_K_2pim_pip_pim5 = P_K_2pim_pip_pim5.M();
            double m_K_2pim_pip_pim6 = P_K_2pim_pip_pim6.M();
            double m_K_2pim_pip_pim7 = P_K_2pim_pip_pim7.M();
            double m_K_2pim_pip_pim8 = P_K_2pim_pip_pim8.M();
            double m_K_2pim_pip_pim9 = P_K_2pim_pip_pim9.M();
            double m_K_2pim_pip_pim10 = P_K_2pim_pip_pim10.M();
            double m_K_2pim_pip_pim11 = P_K_2pim_pip_pim11.M();
            double m_K_2pim_pip_pim12 = P_K_2pim_pip_pim12.M();
            double m_K_2pim_pip_pim13 = P_K_2pim_pip_pim13.M();
            double m_K_2pim_pip_pim14 = P_K_2pim_pip_pim14.M();
            double m_K_2pim_pip_pim15 = P_K_2pim_pip_pim15.M();
            double m_K_2pim_pip_pim16 = P_K_2pim_pip_pim16.M();
            double m_K_2pim_pip_pim17 = P_K_2pim_pip_pim17.M();
            double m_K_2pim_pip_pim18 = P_K_2pim_pip_pim18.M();
            double m_K_2pim_pip_pim19 = P_K_2pim_pip_pim19.M();
            double m_K_2pim_pip_pim20 = P_K_2pim_pip_pim20.M();
            double m_K_2pim_pip_pim21 = P_K_2pim_pip_pim21.M();
            double m_K_2pim_pip_pim22 = P_K_2pim_pip_pim22.M();
            double m_K_2pim_pip_pim23 = P_K_2pim_pip_pim23.M();
            double m_K_2pim_pip_pim24 = P_K_2pim_pip_pim24.M();
            double m_K_2pim_pip_pim25 = P_K_2pim_pip_pim25.M();
            double m_K_2pim_pip_pim26 = P_K_2pim_pip_pim26.M();
            double m_K_2pim_pip_pim27 = P_K_2pim_pip_pim27.M();

            // 2(pi+pi-)pi+
            double m_2pippim_pip1 = P_2pippim_pip1.M();
            double m_2pippim_pip2 = P_2pippim_pip2.M();
            double m_2pippim_pip3 = P_2pippim_pip3.M();

            // (2pi)+pi-
            double m_2pip_pim1 = P_2pip_pim1.M();
            double m_2pip_pim2 = P_2pip_pim2.M();
            double m_2pip_pim3 = P_2pip_pim3.M();
            double m_2pip_pim4 = P_2pip_pim4.M();
            double m_2pip_pim5 = P_2pip_pim5.M();
            double m_2pip_pim6 = P_2pip_pim6.M();
            double m_2pip_pim7 = P_2pip_pim7.M();
            double m_2pip_pim8 = P_2pip_pim8.M();
            double m_2pip_pim9 = P_2pip_pim9.M();

            // pi+pi- (misID)
            double m_pip_pim_ID1 = P_pip_pim_ID1.M();
            double m_pip_pim_ID2 = P_pip_pim_ID2.M();
            double m_pip_pim_ID3 = P_pip_pim_ID3.M();
            double m_pip_pim_ID4 = P_pip_pim_ID4.M();
            double m_pip_pim_ID5 = P_pip_pim_ID5.M();
            double m_pip_pim_ID6 = P_pip_pim_ID6.M();
            double m_pip_pim_ID7 = P_pip_pim_ID7.M();
            double m_pip_pim_ID8 = P_pip_pim_ID8.M();
            double m_pip_pim_ID9 = P_pip_pim_ID9.M();

            // K+(pi+pi-) (misID)
            double m_K_2pi_ID1 = P_K_2pi_ID1.M();
            double m_K_2pi_ID2 = P_K_2pi_ID2.M();
            double m_K_2pi_ID3 = P_K_2pi_ID3.M();
            double m_K_2pi_ID4 = P_K_2pi_ID4.M();
            double m_K_2pi_ID5 = P_K_2pi_ID5.M();
            double m_K_2pi_ID6 = P_K_2pi_ID6.M();
            double m_K_2pi_ID7 = P_K_2pi_ID7.M();
            double m_K_2pi_ID8 = P_K_2pi_ID8.M();
            double m_K_2pi_ID9 = P_K_2pi_ID9.M();

            // pi+2pi- (misID)
            double m_pip_2pim_ID1 = P_pip_2pim_ID1.M();
            double m_pip_2pim_ID2 = P_pip_2pim_ID2.M();
            double m_pip_2pim_ID3 = P_pip_2pim_ID3.M();
            double m_pip_2pim_ID4 = P_pip_2pim_ID4.M();
            double m_pip_2pim_ID5 = P_pip_2pim_ID5.M();
            double m_pip_2pim_ID6 = P_pip_2pim_ID6.M();
            double m_pip_2pim_ID7 = P_pip_2pim_ID7.M();
            double m_pip_2pim_ID8 = P_pip_2pim_ID8.M();
            double m_pip_2pim_ID9 = P_pip_2pim_ID9.M();

            // pi+(pi+pi-)
            double m_pip_2pi_ID1 = P_pip_2pi_ID1.M();
            double m_pip_2pi_ID2 = P_pip_2pi_ID2.M();
            double m_pip_2pi_ID3 = P_pip_2pi_ID3.M();
            double m_pip_2pi_ID4 = P_pip_2pi_ID4.M();
            double m_pip_2pi_ID5 = P_pip_2pi_ID5.M();
            double m_pip_2pi_ID6 = P_pip_2pi_ID6.M();
            double m_pip_2pi_ID7 = P_pip_2pi_ID7.M();
            double m_pip_2pi_ID8 = P_pip_2pi_ID8.M();
            double m_pip_2pi_ID9 = P_pip_2pi_ID9.M();

            // 2pi+2pi- (misID)
            double m_2pip_2pim_ID1 = P_2pip_2pim_ID1.M();
            double m_2pip_2pim_ID2 = P_2pip_2pim_ID2.M();
            double m_2pip_2pim_ID3 = P_2pip_2pim_ID3.M();
            double m_2pip_2pim_ID4 = P_2pip_2pim_ID4.M();
            double m_2pip_2pim_ID5 = P_2pip_2pim_ID5.M();
            double m_2pip_2pim_ID6 = P_2pip_2pim_ID6.M();
            double m_2pip_2pim_ID7 = P_2pip_2pim_ID7.M();
            double m_2pip_2pim_ID8 = P_2pip_2pim_ID8.M();
            double m_2pip_2pim_ID9 = P_2pip_2pim_ID9.M();

            // pi+
            double m_pip1 = P_pip1.M();
            double m_pip2 = P_pip2.M();
            double m_pip3 = P_pip3.M();

            // 2pi (same tau)
            double m_taup_2pi_1 = P_taup_2pi_1.M();
            double m_taup_2pi_2 = P_taup_2pi_2.M();
            double m_taup_2pi_3 = P_taup_2pi_3.M();
            double m_taum_2pi_1 = P_taum_2pi_1.M();
            double m_taum_2pi_2 = P_taum_2pi_2.M();
            double m_taum_2pi_3 = P_taum_2pi_3.M();

            // 2pi (different tau)
            double m_2pi_1 = P_2pi_1.M();
            double m_2pi_2 = P_2pi_2.M();
            double m_2pi_3 = P_2pi_3.M();
            double m_2pi_4 = P_2pi_4.M();
            double m_2pi_5 = P_2pi_5.M();
            double m_2pi_6 = P_2pi_6.M();
            double m_2pi_7 = P_2pi_7.M();
            double m_2pi_8 = P_2pi_8.M();
            double m_2pi_9 = P_2pi_9.M();

            // (3pi)^+
            h_3pip->Fill(m_3pip);
            h_3pim->Fill(m_3pim);

            // (pi+pi-)
            h_pip_pim1->Fill(m_pip_pim1);
            h_pip_pim2->Fill(m_pip_pim2);
            h_pip_pim3->Fill(m_pip_pim3);
            h_pip_pim4->Fill(m_pip_pim4);

            // Kpi
            h_Kpi_1->Fill(m_Kpi_1);
            h_Kpi_2->Fill(m_Kpi_2);
            h_Kpi_3->Fill(m_Kpi_3);

            // K(2pi)+
            h_K2pip1->Fill(m_K2pip1);
            h_K2pip2->Fill(m_K2pip2);
            h_K2pip3->Fill(m_K2pip3);

            // 2(pi+pi-)
            h_2pip_2pim1->Fill(m_2pip_2pim1);
            h_2pip_2pim2->Fill(m_2pip_2pim2);
            h_2pip_2pim3->Fill(m_2pip_2pim3);
            h_2pip_2pim4->Fill(m_2pip_2pim4);

            // K(2pi)+pi-
            h_K_2pim_pip1->Fill(m_K_2pim_pip1);
            h_K_2pim_pip2->Fill(m_K_2pim_pip2);
            h_K_2pim_pip3->Fill(m_K_2pim_pip3);
            h_K_2pim_pip4->Fill(m_K_2pim_pip4);
            h_K_2pim_pip5->Fill(m_K_2pim_pip5);
            h_K_2pim_pip6->Fill(m_K_2pim_pip6);
            h_K_2pim_pip7->Fill(m_K_2pim_pip7);
            h_K_2pim_pip8->Fill(m_K_2pim_pip8);
            h_K_2pim_pip9->Fill(m_K_2pim_pip9);

            // Kpi+pi-
            h_K_pip_pim1->Fill(m_K_pip_pim1);
            h_K_pip_pim2->Fill(m_K_pip_pim2);
            h_K_pip_pim3->Fill(m_K_pip_pim3);
            h_K_pip_pim4->Fill(m_K_pip_pim4);
            h_K_pip_pim5->Fill(m_K_pip_pim5);
            h_K_pip_pim6->Fill(m_K_pip_pim6);
            h_K_pip_pim7->Fill(m_K_pip_pim7);
            h_K_pip_pim8->Fill(m_K_pip_pim8);
            h_K_pip_pim9->Fill(m_K_pip_pim9);

            // K(2pi)-pi+pi-
            h_K_2pim_pip_pim1->Fill(m_K_2pim_pip_pim1);
            h_K_2pim_pip_pim2->Fill(m_K_2pim_pip_pim2);
            h_K_2pim_pip_pim3->Fill(m_K_2pim_pip_pim3);
            h_K_2pim_pip_pim4->Fill(m_K_2pim_pip_pim4);
            h_K_2pim_pip_pim5->Fill(m_K_2pim_pip_pim5);
            h_K_2pim_pip_pim6->Fill(m_K_2pim_pip_pim6);
            h_K_2pim_pip_pim7->Fill(m_K_2pim_pip_pim7);
            h_K_2pim_pip_pim8->Fill(m_K_2pim_pip_pim8);
            h_K_2pim_pip_pim9->Fill(m_K_2pim_pip_pim9);
            h_K_2pim_pip_pim10->Fill(m_K_2pim_pip_pim10);
            h_K_2pim_pip_pim11->Fill(m_K_2pim_pip_pim11);
            h_K_2pim_pip_pim12->Fill(m_K_2pim_pip_pim12);
            h_K_2pim_pip_pim13->Fill(m_K_2pim_pip_pim13);
            h_K_2pim_pip_pim14->Fill(m_K_2pim_pip_pim14);
            h_K_2pim_pip_pim15->Fill(m_K_2pim_pip_pim15);
            h_K_2pim_pip_pim16->Fill(m_K_2pim_pip_pim16);
            h_K_2pim_pip_pim17->Fill(m_K_2pim_pip_pim17);
            h_K_2pim_pip_pim18->Fill(m_K_2pim_pip_pim18);
            h_K_2pim_pip_pim19->Fill(m_K_2pim_pip_pim19);
            h_K_2pim_pip_pim20->Fill(m_K_2pim_pip_pim20);
            h_K_2pim_pip_pim21->Fill(m_K_2pim_pip_pim21);
            h_K_2pim_pip_pim22->Fill(m_K_2pim_pip_pim22);
            h_K_2pim_pip_pim23->Fill(m_K_2pim_pip_pim23);
            h_K_2pim_pip_pim24->Fill(m_K_2pim_pip_pim24);
            h_K_2pim_pip_pim25->Fill(m_K_2pim_pip_pim25);
            h_K_2pim_pip_pim26->Fill(m_K_2pim_pip_pim26);
            h_K_2pim_pip_pim27->Fill(m_K_2pim_pip_pim27);

            // 2(pi+pi-)pi+
            h_2pippim_pip1->Fill(m_2pippim_pip1);
            h_2pippim_pip2->Fill(m_2pippim_pip2);
            h_2pippim_pip3->Fill(m_2pippim_pip3);

            // (2pi)+pi-
            h_2pip_pim1->Fill(m_2pip_pim1);
            h_2pip_pim2->Fill(m_2pip_pim2);
            h_2pip_pim3->Fill(m_2pip_pim3);
            h_2pip_pim4->Fill(m_2pip_pim4);
            h_2pip_pim5->Fill(m_2pip_pim5);
            h_2pip_pim6->Fill(m_2pip_pim6);
            h_2pip_pim7->Fill(m_2pip_pim7);
            h_2pip_pim8->Fill(m_2pip_pim8);
            h_2pip_pim9->Fill(m_2pip_pim9);

            // pi+pi- (misID)
            h_pip_pim_ID1->Fill(m_pip_pim_ID1);
            h_pip_pim_ID2->Fill(m_pip_pim_ID2);
            h_pip_pim_ID3->Fill(m_pip_pim_ID3);
            h_pip_pim_ID4->Fill(m_pip_pim_ID4);
            h_pip_pim_ID5->Fill(m_pip_pim_ID5);
            h_pip_pim_ID6->Fill(m_pip_pim_ID6);
            h_pip_pim_ID7->Fill(m_pip_pim_ID7);
            h_pip_pim_ID8->Fill(m_pip_pim_ID8);
            h_pip_pim_ID9->Fill(m_pip_pim_ID9);

            // K+(pi+pi-) (misID)
            h_K_2pi_ID1->Fill(m_K_2pi_ID1);
            h_K_2pi_ID2->Fill(m_K_2pi_ID2);
            h_K_2pi_ID3->Fill(m_K_2pi_ID3);
            h_K_2pi_ID4->Fill(m_K_2pi_ID4);
            h_K_2pi_ID5->Fill(m_K_2pi_ID5);
            h_K_2pi_ID6->Fill(m_K_2pi_ID6);
            h_K_2pi_ID7->Fill(m_K_2pi_ID7);
            h_K_2pi_ID8->Fill(m_K_2pi_ID8);
            h_K_2pi_ID9->Fill(m_K_2pi_ID9);

            // pi+2pi- (misiD)
            h_pip_2pim_ID1->Fill(m_pip_2pim_ID1);
            h_pip_2pim_ID2->Fill(m_pip_2pim_ID2);
            h_pip_2pim_ID3->Fill(m_pip_2pim_ID3);
            h_pip_2pim_ID4->Fill(m_pip_2pim_ID4);
            h_pip_2pim_ID5->Fill(m_pip_2pim_ID5);
            h_pip_2pim_ID6->Fill(m_pip_2pim_ID6);
            h_pip_2pim_ID7->Fill(m_pip_2pim_ID7);
            h_pip_2pim_ID8->Fill(m_pip_2pim_ID8);
            h_pip_2pim_ID9->Fill(m_pip_2pim_ID9);

            // pi+(pi+pi-) (misiD)
            h_pip_2pi_ID1->Fill(m_pip_2pi_ID1);
            h_pip_2pi_ID2->Fill(m_pip_2pi_ID2);
            h_pip_2pi_ID3->Fill(m_pip_2pi_ID3);
            h_pip_2pi_ID4->Fill(m_pip_2pi_ID4);
            h_pip_2pi_ID5->Fill(m_pip_2pi_ID5);
            h_pip_2pi_ID6->Fill(m_pip_2pi_ID6);
            h_pip_2pi_ID7->Fill(m_pip_2pi_ID7);
            h_pip_2pi_ID8->Fill(m_pip_2pi_ID8);
            h_pip_2pi_ID9->Fill(m_pip_2pi_ID9);

            // 2pi+2pi- (misiD)
            h_2pip_2pim_ID1->Fill(m_2pip_2pim_ID1);
            h_2pip_2pim_ID2->Fill(m_2pip_2pim_ID2);
            h_2pip_2pim_ID3->Fill(m_2pip_2pim_ID3);
            h_2pip_2pim_ID4->Fill(m_2pip_2pim_ID4);
            h_2pip_2pim_ID5->Fill(m_2pip_2pim_ID5);
            h_2pip_2pim_ID6->Fill(m_2pip_2pim_ID6);
            h_2pip_2pim_ID7->Fill(m_2pip_2pim_ID7);
            h_2pip_2pim_ID8->Fill(m_2pip_2pim_ID8);
            h_2pip_2pim_ID9->Fill(m_2pip_2pim_ID9);

            // pi+
            h_pip1->Fill(m_pip1);
            h_pip2->Fill(m_pip2);
            h_pip3->Fill(m_pip3);

            // 2pi (same tau)
            h_taup_2pi_1->Fill(m_taup_2pi_1);
            h_taup_2pi_2->Fill(m_taup_2pi_2);
            h_taup_2pi_3->Fill(m_taup_2pi_3);
            h_taum_2pi_1->Fill(m_taum_2pi_1);
            h_taum_2pi_2->Fill(m_taum_2pi_2);
            h_taum_2pi_3->Fill(m_taum_2pi_3);

            // 2pi (different tau)
            h_2pi_1->Fill(m_2pi_1);
            h_2pi_2->Fill(m_2pi_2);
            h_2pi_3->Fill(m_2pi_3);
            h_2pi_4->Fill(m_2pi_4);
            h_2pi_5->Fill(m_2pi_5);
            h_2pi_6->Fill(m_2pi_6);
            h_2pi_7->Fill(m_2pi_7);
            h_2pi_8->Fill(m_2pi_8);
            h_2pi_9->Fill(m_2pi_9);

        }
    }

    fout->cd();
    // (3pi)^+
    std::vector mass_3pi = {m_K, m_Dpm};
    plot(fout, h_3pip, mass_3pi, "(3#pi)^{-} invariant mass", "m(3#pi) (MeV)", year, 1);
    plot(fout, h_3pip, mass_3pi, "(3#pi)^{-} invariant mass", "m(3#pi) (MeV)", year, 2);

    // (pi+pi-)
    std::vector mass_pip_pim = {m_ro, m_w, m_Ks};
    plot(fout, h_pip_pim1, mass_pip_pim, "#pi^{+} #pi^{-} invariant mass nº1", "m(#pi^{+} #pi^{-}) (MeV)", year, 3);
    plot(fout, h_pip_pim2, mass_pip_pim, "#pi^{+} #pi^{-} invariant mass nº2", "m(#pi^{+} #pi^{-}) (MeV)", year, 4);
    plot(fout, h_pip_pim3, mass_pip_pim, "#pi^{+} #pi^{-} invariant mass nº3", "m(#pi^{+} #pi^{-}) (MeV)", year, 5);
    plot(fout, h_pip_pim4, mass_pip_pim, "#pi^{+} #pi^{-} invariant mass nº4", "m(#pi^{+} #pi^{-}) (MeV)", year, 6);

    // Kpi
    std::vector mass_Kpi = {m_Kstar, m_D0};
    plot(fout, h_Kpi_1, mass_Kpi, "K^{+} #pi^{-} invariant mass nº1", "m(K#pi) (MeV)", year, 7);   
    plot(fout, h_Kpi_2, mass_Kpi, "K^{+} #pi^{-} invariant mass nº2", "m(K#pi) (MeV)", year, 8);
    plot(fout, h_Kpi_3, mass_Kpi, "K^{+} #pi^{-} invariant mass nº3", "m(K#pi) (MeV)", year, 9); 

    // K(2pi)+
    std::vector mass_K2pip = {m_Dpm};
    plot(fout, h_K2pip1, mass_K2pip, "K^{+} (2#pi)^{-} invariant mass nº1", "m(K2#pi) (MeV)", year, 10);
    plot(fout, h_K2pip2, mass_K2pip, "K^{+} (2#pi)^{-} invariant mass nº2", "m(K2#pi) (MeV)", year, 11);
    plot(fout, h_K2pip3, mass_K2pip, "K^{+} (2#pi)^{-} invariant mass nº3", "m(K2#pi) (MeV)", year, 12);

    // 2(pi+pi-)
    std::vector mass_2pip_2pim = {m_D0};
    plot(fout, h_2pip_2pim1, mass_2pip_2pim, "2(#pi^{+} #pi^{-}) invariant mass nº1", "m(2(#pi^{+}#pi^{-})) (MeV)", year, 13);
    plot(fout, h_2pip_2pim2, mass_2pip_2pim, "2(#pi^{+} #pi^{-}) invariant mass nº2", "m(2(#pi^{+}#pi^{-})) (MeV)", year, 14);
    plot(fout, h_2pip_2pim3, mass_2pip_2pim, "2(#pi^{+} #pi^{-}) invariant mass nº3", "m(2(#pi^{+}#pi^{-})) (MeV)", year, 15);
    plot(fout, h_2pip_2pim4, mass_2pip_2pim, "2(#pi^{+} #pi^{-}) invariant mass nº4", "m(2(#pi^{+}#pi^{-})) (MeV)", year, 16);

    // K(2pi)+pi-
    std::vector mass_K_2pim_pip = {m_D0};
    plot(fout, h_K_2pim_pip1, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº1", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 17);
    plot(fout, h_K_2pim_pip2, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº2", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 18);
    plot(fout, h_K_2pim_pip3, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº3", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 19);
    plot(fout, h_K_2pim_pip4, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº4", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 20);
    plot(fout, h_K_2pim_pip5, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº5", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 21);
    plot(fout, h_K_2pim_pip6, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº6", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 22);
    plot(fout, h_K_2pim_pip7, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº7", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 23);
    plot(fout, h_K_2pim_pip8, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº8", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 24);
    plot(fout, h_K_2pim_pip9, mass_K_2pim_pip, "K^{+} 2#pi^{-} #pi^{+} invariant mass nº9", "m(K(2#pi))^{-} #pi^{+})) (MeV)", year, 25);

    // Kpi+pi-
    std::vector mass_K_pip_pim = {m_Ds_pm, m_Dstar, m_Dstar1};
    plot(fout, h_K_pip_pim1, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº1", "m(K pi^{+} #pi^{-}) (MeV)", year, 26);
    plot(fout, h_K_pip_pim2, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº2", "m(K pi^{+} #pi^{-}) (MeV)", year, 27);
    plot(fout, h_K_pip_pim3, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº3", "m(K pi^{+} #pi^{-}) (MeV)", year, 28);
    plot(fout, h_K_pip_pim4, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº4", "m(K pi^{+} #pi^{-}) (MeV)", year, 29);
    plot(fout, h_K_pip_pim5, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº5", "m(K pi^{+} #pi^{-}) (MeV)", year, 30);
    plot(fout, h_K_pip_pim6, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº6", "m(K pi^{+} #pi^{-}) (MeV)", year, 31);
    plot(fout, h_K_pip_pim7, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº7", "m(K pi^{+} #pi^{-}) (MeV)", year, 32);
    plot(fout, h_K_pip_pim8, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº8", "m(K pi^{+} #pi^{-}) (MeV)", year, 33);
    plot(fout, h_K_pip_pim9, mass_K_pip_pim, "K^{+} #pi^{+} #pi^{-} invariant mass nº9", "m(K pi^{+} #pi^{-}) (MeV)", year, 34);

    // K(2pi-)(pi+pi-)
    std::vector mass_K_2pim_pip_pim = {m_Ds_pm, m_Dstar, m_Dstar1};
    plot(fout, h_K_2pim_pip_pim1, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº1", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 35);
    plot(fout, h_K_2pim_pip_pim2, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº2", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 36);
    plot(fout, h_K_2pim_pip_pim3, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº3", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 37);
    plot(fout, h_K_2pim_pip_pim4, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº4", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 38);
    plot(fout, h_K_2pim_pip_pim5, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº5", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 39);
    plot(fout, h_K_2pim_pip_pim6, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº6", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 40);
    plot(fout, h_K_2pim_pip_pim7, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº7", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 41);
    plot(fout, h_K_2pim_pip_pim8, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº8", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 42);
    plot(fout, h_K_2pim_pip_pim9, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº9", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 43);
    plot(fout, h_K_2pim_pip_pim10, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº10", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 44);
    plot(fout, h_K_2pim_pip_pim11, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº11", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 45);
    plot(fout, h_K_2pim_pip_pim12, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº12", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 46);
    plot(fout, h_K_2pim_pip_pim13, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº13", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 47);
    plot(fout, h_K_2pim_pip_pim14, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº14", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 48);
    plot(fout, h_K_2pim_pip_pim15, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº15", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 49);
    plot(fout, h_K_2pim_pip_pim16, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº16", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 50);
    plot(fout, h_K_2pim_pip_pim17, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº17", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 51);
    plot(fout, h_K_2pim_pip_pim18, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº18", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 52);
    plot(fout, h_K_2pim_pip_pim19, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº19", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 53);
    plot(fout, h_K_2pim_pip_pim20, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº20", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 54);
    plot(fout, h_K_2pim_pip_pim21, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº21", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 55);
    plot(fout, h_K_2pim_pip_pim22, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº22", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 56);
    plot(fout, h_K_2pim_pip_pim23, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº23", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 57);
    plot(fout, h_K_2pim_pip_pim24, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº24", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 58);
    plot(fout, h_K_2pim_pip_pim25, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº25", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 59);
    plot(fout, h_K_2pim_pip_pim26, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº26", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 60);
    plot(fout, h_K_2pim_pip_pim27, mass_K_2pim_pip_pim, "K^{+} (2#pi)^{-} #pi^{+} #pi^{-} invariant mass nº27", "m(K (2#pi)^{-} pi^{+} #pi^{-}) (MeV)", year, 61);

    // 2(pi+pi-)pi+
    std::vector mass_2pippim_pip = {m_Ds_pm, m_Dstar, m_Dstar1};
    plot(fout, h_2pippim_pip1, mass_2pippim_pip, "2(#pi^{+} #pi^{-}) #pi^{+} invariant mass nº1", "m(2(#pi^{+} #pi^{-}) pi^{+}) (MeV)", year, 62);
    plot(fout, h_2pippim_pip2, mass_2pippim_pip, "2(#pi^{+} #pi^{-}) #pi^{+} invariant mass nº2", "m(2(#pi^{+} #pi^{-}) pi^{+}) (MeV)", year, 63);
    plot(fout, h_2pippim_pip3, mass_2pippim_pip, "2(#pi^{+} #pi^{-}) #pi^{+} invariant mass nº3", "m(2(#pi^{+} #pi^{-}) pi^{+}) (MeV)", year, 64);

    // (2pi)+pi-
    std::vector mass_2pip_pim = {m_Ds_pm};
    plot(fout, h_2pip_pim1, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº1", "m((2#pi)^{+} pi^{-}) (MeV)", year, 65);
    plot(fout, h_2pip_pim2, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº2", "m((2#pi)^{+} pi^{-}) (MeV)", year, 66);
    plot(fout, h_2pip_pim3, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº3", "m((2#pi)^{+} pi^{-}) (MeV)", year, 67);
    plot(fout, h_2pip_pim4, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº4", "m((2#pi)^{+} pi^{-}) (MeV)", year, 68);
    plot(fout, h_2pip_pim5, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº5", "m((2#pi)^{+} pi^{-}) (MeV)", year, 69);
    plot(fout, h_2pip_pim6, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº6", "m((2#pi)^{+} pi^{-}) (MeV)", year, 70);
    plot(fout, h_2pip_pim7, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº7", "m((2#pi)^{+} pi^{-}) (MeV)", year, 71);
    plot(fout, h_2pip_pim8, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº8", "m((2#pi)^{+} pi^{-}) (MeV)", year, 72);
    plot(fout, h_2pip_pim9, mass_2pip_pim, "(2#pi)^{+} #pi^{-} invariant mass nº9", "m((2#pi)^{+} pi^{-}) (MeV)", year, 73);

    // pi+pi- (misID)
    std::vector mass_pip_pim_ID = {m_phi, m_Kstar, m_D0};
    plot(fout, h_pip_pim_ID1, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº1", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 74);    
    plot(fout, h_pip_pim_ID2, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº2", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 75);    
    plot(fout, h_pip_pim_ID3, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº3", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 76);    
    plot(fout, h_pip_pim_ID4, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº4", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 77);    
    plot(fout, h_pip_pim_ID5, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº5", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 78);    
    plot(fout, h_pip_pim_ID6, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº6", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 79);    
    plot(fout, h_pip_pim_ID7, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº7", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 80);    
    plot(fout, h_pip_pim_ID8, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº8", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 81);    
    plot(fout, h_pip_pim_ID9, mass_pip_pim_ID, "#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº9", "m(#pi^{+}#pi^{-})^{misID} (MeV)", year, 82);    

    // K+(pi+pi-) (misID)
    std::vector mass_K_2pi_ID = {m_Ds_pm};
    plot(fout, h_K_2pi_ID1, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº1", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 83);    
    plot(fout, h_K_2pi_ID2, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº2", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 84);    
    plot(fout, h_K_2pi_ID3, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº3", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 85);    
    plot(fout, h_K_2pi_ID4, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº4", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 86);    
    plot(fout, h_K_2pi_ID5, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº5", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 87);    
    plot(fout, h_K_2pi_ID6, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº6", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 88);    
    plot(fout, h_K_2pi_ID7, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº7", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 89);    
    plot(fout, h_K_2pi_ID8, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº8", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 90);    
    plot(fout, h_K_2pi_ID9, mass_K_2pi_ID, "K^{+}#pi^{+}#pi^{-} (K<->#pi mis ID) invariant mass nº9", "m(K^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 91);    

    // pi+2pi- (misID)
    std::vector mass_pip_2pim_ID = {m_Dpm};
    plot(fout, h_pip_2pim_ID1, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº1", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 92);    
    plot(fout, h_pip_2pim_ID2, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº2", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 93);    
    plot(fout, h_pip_2pim_ID3, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº3", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 94);    
    plot(fout, h_pip_2pim_ID4, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº4", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 95);    
    plot(fout, h_pip_2pim_ID5, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº5", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 96);    
    plot(fout, h_pip_2pim_ID6, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº6", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 97);    
    plot(fout, h_pip_2pim_ID7, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº7", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 98);    
    plot(fout, h_pip_2pim_ID8, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº8", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 99);    
    plot(fout, h_pip_2pim_ID9, mass_pip_2pim_ID, "#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº9", "m(#pi^{+}2#pi^{-})^{misID} (MeV)", year, 100);    

    // pi+(pi+pi-) (misID)
    std::vector mass_pip_2pi_ID = {m_Dstar, m_Dstar1};
    plot(fout, h_pip_2pi_ID1, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº1", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 101);    
    plot(fout, h_pip_2pi_ID2, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº2", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 102);    
    plot(fout, h_pip_2pi_ID3, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº3", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 103);    
    plot(fout, h_pip_2pi_ID4, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº4", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 104);    
    plot(fout, h_pip_2pi_ID5, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº5", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 105);    
    plot(fout, h_pip_2pi_ID6, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº6", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 106);    
    plot(fout, h_pip_2pi_ID7, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº7", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 107);    
    plot(fout, h_pip_2pi_ID8, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº8", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 108);    
    plot(fout, h_pip_2pi_ID9, mass_pip_2pi_ID, "#pi^{+}(#pi^{+}#pi^{-}) (K<->#pi mis ID) invariant mass nº9", "m(#pi^{+}#pi^{+}#pi^{-})^{misID} (MeV)", year, 109);    

    // 2pi+2pi- (misID)
    std::vector mass_2pip_2pim_ID = {m_D0};    
    plot(fout, h_2pip_2pim_ID1, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº1", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 110);
    plot(fout, h_2pip_2pim_ID2, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº2", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 111);
    plot(fout, h_2pip_2pim_ID3, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº3", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 112);
    plot(fout, h_2pip_2pim_ID4, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº4", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 113);
    plot(fout, h_2pip_2pim_ID5, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº5", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 114);
    plot(fout, h_2pip_2pim_ID6, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº6", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 115);
    plot(fout, h_2pip_2pim_ID7, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº7", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 116);
    plot(fout, h_2pip_2pim_ID8, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº8", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 117);
    plot(fout, h_2pip_2pim_ID9, mass_2pip_2pim_ID, "2#pi^{+}2#pi^{-} (K<->#pi mis ID) invariant mass nº9", "m(2#pi^{+}2#pi^{-})^{misID} (MeV)", year, 118);

    // pi+
    std::vector mass_pip = {m_K};   
    plot(fout, h_pip1, mass_pip, "#pi^{+} invariant mass nº1", "m(#pi^{+}) (MeV)", year, 119);
    plot(fout, h_pip2, mass_pip, "#pi^{+} invariant mass nº2", "m(#pi^{+}) (MeV)", year, 120);
    plot(fout, h_pip3, mass_pip, "#pi^{+} invariant mass nº3", "m(#pi^{+}) (MeV)", year, 121);

    // 2pi (same tau)
    plot1(fout, h_taup_2pi_1, "2#pi (from #tau^{+}) invariant mass nº1", "m(#pi^{+} #pi^{-}) (MeV)", year, 122);
    plot1(fout, h_taup_2pi_2, "2#pi (from #tau^{+}) invariant mass nº2", "m(#pi^{+} #pi^{+}) (MeV)", year, 123);
    plot1(fout, h_taup_2pi_3, "2#pi (from #tau^{+}) invariant mass nº3", "m(#pi^{+} #pi^{-}) (MeV)", year, 124);
    plot1(fout, h_taum_2pi_1, "2#pi (from #tau^{-}) invariant mass nº1", "m(#pi^{-} #pi^{-}) (MeV)", year, 125);
    plot1(fout, h_taum_2pi_2, "2#pi (from #tau^{-}) invariant mass nº2", "m(#pi^{+} #pi^{-}) (MeV)", year, 126);
    plot1(fout, h_taum_2pi_3, "2#pi (from #tau^{-}) invariant mass nº3", "m(#pi^{+} #pi^{-}) (MeV)", year, 127);

    // 2pi (different tau)
    plot1(fout, h_2pi_1, "2#pi (from different #tau) invariant mass nº1", "m(#pi^{-} #pi^{-}) (MeV)", year, 128);
    plot1(fout, h_2pi_2, "2#pi (from different #tau) invariant mass nº2", "m(#pi^{-} #pi^{-}) (MeV)", year, 129);
    plot1(fout, h_2pi_3, "2#pi (from different #tau) invariant mass nº3", "m(#pi^{-} #pi^{+}) (MeV)", year, 130);
    plot1(fout, h_2pi_4, "2#pi (from different #tau) invariant mass nº4", "m(#pi^{+} #pi^{-}) (MeV)", year, 131);
    plot1(fout, h_2pi_5, "2#pi (from different #tau) invariant mass nº5", "m(#pi^{+} #pi^{-}) (MeV)", year, 132);
    plot1(fout, h_2pi_6, "2#pi (from different #tau) invariant mass nº6", "m(#pi^{+} #pi^{+}) (MeV)", year, 133);
    plot1(fout, h_2pi_7, "2#pi (from different #tau) invariant mass nº7", "m(#pi^{+} #pi^{-}) (MeV)", year, 134);
    plot1(fout, h_2pi_8, "2#pi (from different #tau) invariant mass nº8", "m(#pi^{+} #pi^{-}) (MeV)", year, 135);
    plot1(fout, h_2pi_9, "2#pi (from different #tau) invariant mass nº9", "m(#pi^{+} #pi^{+}) (MeV)", year, 136);

    fout->Close();
}

void plot(TFile* fout, TH1D* h_mass, std::vector<double> m, TString title, TString xlabel, TString year, int i){

    h_mass->SetLineColor(kBlack);
    h_mass->SetLineWidth(2);
    h_mass->GetXaxis()->SetTitle(xlabel);
    h_mass->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass->SetTitle(title+" ("+year+")");

    TCanvas c;
    c.cd();
    c.SetTitle(Form("c_%i",i));
    c.SetName(Form("c_%i",i));
    c.Update();
    h_mass->Draw();

    int size_m = sizeof(m)/sizeof(m[0]);
    double ymin = h_mass->GetMinimum();
    double ymax = h_mass->GetMaximum();
    for(int i = 0; i < size_m; ++i){
        TLine* line = new TLine(m[i], ymin, m[i], ymax);
        line->SetLineColor(kBlue);
        line->SetLineWidth(2);
        line->Draw("same");
    }
    c.Write();

}

void plot1(TFile* fout, TH1D* h_mass, TString title, TString xlabel, TString year, int i){

    h_mass->SetLineColor(kBlack);
    h_mass->SetLineWidth(2);
    h_mass->GetXaxis()->SetTitle(xlabel);
    h_mass->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass->SetTitle(title+" ("+year+")");

    TCanvas c;
    c.cd();
    c.SetTitle(Form("c_%i",i));
    c.SetName(Form("c_%i",i));
    c.Update();
    h_mass->Draw();

    c.Write();

}

void plot_comparison(TFile* fout, TH1D* h_mass, TH1D* h_mass2, TH1D* h_mass3, std::vector<double> m, TString title, TString xlabel, TString year, int i){

    h_mass->SetLineColor(kBlack);
    h_mass->SetLineWidth(2);
    h_mass->GetXaxis()->SetTitle(xlabel);
    h_mass->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass->SetTitle(title+" ("+year+")");

    h_mass2->SetLineColor(kRed);
    h_mass2->SetLineWidth(2);

    h_mass3->SetLineColor(kGreen);
    h_mass3->SetLineWidth(2);

    TLegend* leg;
    leg = new TLegend(0.6, 0.7, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->AddEntry(h_mass, "No cut", "lp");
    leg->AddEntry(h_mass2, "#tau vtx #chi^{2} < 4", "lp");
    leg->AddEntry(h_mass3, "#tau vtx #chi^{2} < 9", "lp");

    TCanvas c;
    c.cd();
    c.SetTitle(Form("c_%i",i));
    c.SetName(Form("c_%i",i));
    gStyle->SetOptStat(0);
    h_mass->DrawNormalized();
    h_mass2->DrawNormalized("same");
    h_mass3->DrawNormalized("same");

    int size_m = sizeof(m)/sizeof(m[0]);
    double ymin = h_mass->GetMinimum();
    double ymax = h_mass->GetMaximum();
    for(int i = 0; i < size_m; ++i){
        TLine* line = new TLine(m[i], ymin, m[i], ymax);
        line->SetLineColor(kBlue);
        line->SetLineWidth(2);
        line->Draw("same");
    }
    leg->Draw("same");

    c.Write();

}