RooRealVar* fit(RooDataSet* data, RooRealVar mass, TString cut, TString name, TString folder_naeme);

using namespace RooStats;
using namespace RooFit;
using namespace std;

std::ofstream file;

bool make_fit = false;
bool myInit = false;

void mass_resolution(){

    TFile* fin;
    if(myInit){fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016_myInit.root");}
    else{fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");}
    TTree* t = (TTree*)fin->Get("ntuple/DecayTree");
    TFile* f_kstar = new TFile("/panfs/avenkate/KstTauTau/2016_MD/Stripping/ntuple_TM.root");
    TTree* t_kstar = (TTree*)f_kstar->Get("DecayTreeTuple");

    TString folder_name;
    if(myInit){folder_name = "Mass_resolution_visualisation_myInit/";}
    else{folder_name = "Mass_resolution_visualisation/";}

    Double_t Bmass_visible, taup_FD, taum_FD, taup_FD_chi2, taum_FD_chi2, taup_IP_chi2, taum_IP_chi2;
    Float_t status, Bmass, Bmass_err, Enu1, Pnu1x, Pnu1y, Pnu1z, Enu2, Pnu2x, Pnu2y, Pnu2z, DTF_EB, DTF_PBx, DTF_PBy, DTF_PBz;
    Int_t Kp_TRUEID, taup_pip0_TRUEID, taup_pim0_TRUEID, taup_pip1_TRUEID, taum_pim0_TRUEID, taum_pip0_TRUEID, taum_pim1_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;
    Double_t DV1x, DV1y, DV1z, DV2x, DV2y, DV2z, BVx, BVy, BVz, PVx, PVy, PVz, KVx, KVy, KVz;
    Double_t taup_pip0_TRUEPE, taup_pim0_TRUEPE, taup_pip1_TRUEPE, taum_pim0_TRUEPE, taum_pip0_TRUEPE, taum_pim1_TRUEPE, taup_TRUEPE, taum_TRUEPE;
    Double_t taup_pip0_TRUEPX, taup_pim0_TRUEPX, taup_pip1_TRUEPX, taum_pim0_TRUEPX, taum_pip0_TRUEPX, taum_pim1_TRUEPX, taup_TRUEPX, taum_TRUEPX;
    Double_t taup_pip0_TRUEPY, taup_pim0_TRUEPY, taup_pip1_TRUEPY, taum_pim0_TRUEPY, taum_pip0_TRUEPY, taum_pim1_TRUEPY, taup_TRUEPY, taum_TRUEPY;
    Double_t taup_pip0_TRUEPZ, taup_pim0_TRUEPZ, taup_pip1_TRUEPZ, taum_pim0_TRUEPZ, taum_pip0_TRUEPZ, taum_pim1_TRUEPZ, taup_TRUEPZ, taum_TRUEPZ;
    Float_t taum_DTF_M, taup_DTF_M, taum_DTF_Merr, taup_DTF_Merr, DTF_chi2, DTF_nIter;
    Double_t Kp_PT, taup_PT, taum_PT, Pkx, Pky, Pkz;
    Float_t M_kstar, status_kstar, Enu1_s, Pnu1x_s, Pnu1y_s, Pnu1z_s, Enu2_s, Pnu2x_s, Pnu2y_s, Pnu2z_s;
    Double_t taup_pip0_TRUEPE_s, taup_pim0_TRUEPE_s, taup_pip1_TRUEPE_s, taum_pim0_TRUEPE_s, taum_pip0_TRUEPE_s, taum_pim1_TRUEPE_s, taup_TRUEPE_s, taum_TRUEPE_s;
    Double_t taup_pip0_TRUEPX_s, taup_pim0_TRUEPX_s, taup_pip1_TRUEPX_s, taum_pim0_TRUEPX_s, taum_pip0_TRUEPX_s, taum_pim1_TRUEPX_s, taup_TRUEPX_s, taum_TRUEPX_s;
    Double_t taup_pip0_TRUEPY_s, taup_pim0_TRUEPY_s, taup_pip1_TRUEPY_s, taum_pim0_TRUEPY_s, taum_pip0_TRUEPY_s, taum_pim1_TRUEPY_s, taup_TRUEPY_s, taum_TRUEPY_s;
    Double_t taup_pip0_TRUEPZ_s, taup_pim0_TRUEPZ_s, taup_pip1_TRUEPZ_s, taum_pim0_TRUEPZ_s, taum_pip0_TRUEPZ_s, taum_pim1_TRUEPZ_s, taup_TRUEPZ_s, taum_TRUEPZ_s;

    t->SetBranchAddress("taup_TRUEP_E", &taup_TRUEPE);
    t->SetBranchAddress("taum_TRUEP_E", &taum_TRUEPE);
    t->SetBranchAddress("taup_pip0_TRUEP_E", &taup_pip0_TRUEPE);
    t->SetBranchAddress("taup_pim0_TRUEP_E", &taup_pim0_TRUEPE);
    t->SetBranchAddress("taup_pip1_TRUEP_E", &taup_pip1_TRUEPE);
    t->SetBranchAddress("taum_pim0_TRUEP_E", &taum_pim0_TRUEPE);
    t->SetBranchAddress("taum_pip0_TRUEP_E", &taum_pip0_TRUEPE);
    t->SetBranchAddress("taum_pim1_TRUEP_E", &taum_pim1_TRUEPE); 

    t->SetBranchAddress("taup_TRUEP_X", &taup_TRUEPX);
    t->SetBranchAddress("taum_TRUEP_X", &taum_TRUEPX);
    t->SetBranchAddress("taup_pip0_TRUEP_X", &taup_pip0_TRUEPX);
    t->SetBranchAddress("taup_pim0_TRUEP_X", &taup_pim0_TRUEPX);
    t->SetBranchAddress("taup_pip1_TRUEP_X", &taup_pip1_TRUEPX);
    t->SetBranchAddress("taum_pim0_TRUEP_X", &taum_pim0_TRUEPX);
    t->SetBranchAddress("taum_pip0_TRUEP_X", &taum_pip0_TRUEPX);
    t->SetBranchAddress("taum_pim1_TRUEP_X", &taum_pim1_TRUEPX); 

    t->SetBranchAddress("taup_TRUEP_Y", &taup_TRUEPY);
    t->SetBranchAddress("taum_TRUEP_Y", &taum_TRUEPY);
    t->SetBranchAddress("taup_pip0_TRUEP_Y", &taup_pip0_TRUEPY);
    t->SetBranchAddress("taup_pim0_TRUEP_Y", &taup_pim0_TRUEPY);
    t->SetBranchAddress("taup_pip1_TRUEP_Y", &taup_pip1_TRUEPY);
    t->SetBranchAddress("taum_pim0_TRUEP_Y", &taum_pim0_TRUEPY);
    t->SetBranchAddress("taum_pip0_TRUEP_Y", &taum_pip0_TRUEPY);
    t->SetBranchAddress("taum_pim1_TRUEP_Y", &taum_pim1_TRUEPY); 

    t->SetBranchAddress("taup_TRUEP_Z", &taup_TRUEPZ);
    t->SetBranchAddress("taum_TRUEP_Z", &taum_TRUEPZ);
    t->SetBranchAddress("taup_pip0_TRUEP_Z", &taup_pip0_TRUEPZ);
    t->SetBranchAddress("taup_pim0_TRUEP_Z", &taup_pim0_TRUEPZ);
    t->SetBranchAddress("taup_pip1_TRUEP_Z", &taup_pip1_TRUEPZ);
    t->SetBranchAddress("taum_pim0_TRUEP_Z", &taum_pim0_TRUEPZ);
    t->SetBranchAddress("taum_pip0_TRUEP_Z", &taum_pip0_TRUEPZ);
    t->SetBranchAddress("taum_pim1_TRUEP_Z", &taum_pim1_TRUEPZ);    

    t->SetBranchAddress("Bp_ConsBp_M",&Bmass);
    t->SetBranchAddress("Bp_M",&Bmass_visible);
    t->SetBranchAddress("Bp_ConsBp_MERR",&Bmass_err);
    t->SetBranchAddress("taup_FD_ORIVX",&taup_FD);
    t->SetBranchAddress("taum_FD_ORIVX",&taum_FD); 
    t->SetBranchAddress("taup_FDCHI2_ORIVX",&taup_FD_chi2);
    t->SetBranchAddress("taum_FDCHI2_ORIVX",&taum_FD_chi2);
    t->SetBranchAddress("taup_IPCHI2_OWNPV",&taup_IP_chi2);
    t->SetBranchAddress("taum_IPCHI2_OWNPV",&taum_IP_chi2);

    // B+ kinematics
    t->SetBranchAddress("Bp_ConsBp_PE",&DTF_EB);
    t->SetBranchAddress("Bp_ConsBp_PX",&DTF_PBx);
    t->SetBranchAddress("Bp_ConsBp_PY",&DTF_PBy);
    t->SetBranchAddress("Bp_ConsBp_PZ",&DTF_PBz);

    // DTF neutrino kinematics
    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PE",&Enu1);
    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PX",&Pnu1x);
    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PY",&Pnu1y);
    t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PZ",&Pnu1z);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PE",&Enu2);
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
    t->SetBranchAddress("Bp_OWNPV_X",&PVx);
    t->SetBranchAddress("Bp_OWNPV_Y",&PVy);
    t->SetBranchAddress("Bp_OWNPV_Z",&PVz);
    t->SetBranchAddress("refPoint_X",&KVx);
    t->SetBranchAddress("refPoint_Y",&KVy);
    t->SetBranchAddress("refPoint_Z",&KVz);

    t->SetBranchAddress("Kp_PX",&Pkx);
    t->SetBranchAddress("Kp_PY",&Pky);
    t->SetBranchAddress("Kp_PZ",&Pkz);

    t->SetBranchAddress("Bp_ConsBp_tauminus_M",&taup_DTF_M);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_M",&taum_DTF_M);
    t->SetBranchAddress("Bp_ConsBp_tauminus_MERR",&taup_DTF_Merr);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_MERR",&taum_DTF_Merr);
    t->SetBranchAddress("Bp_ConsBp_chi2",&DTF_chi2);
    t->SetBranchAddress("Bp_ConsBp_nIter",&DTF_nIter);

    t->SetBranchAddress("Kp_PT",&Kp_PT);
    t->SetBranchAddress("taup_PT",&taup_PT);
    t->SetBranchAddress("taum_PT",&taum_PT);

    t_kstar->SetBranchAddress("B0_kttdtf0_M",&M_kstar);
    t_kstar->SetBranchAddress("B0_kttdtf0_status",&status_kstar);

    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_nu_tau_PE",&Enu1_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_nu_tau_PX",&Pnu1x_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_nu_tau_PY",&Pnu1y_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_nu_tau_PZ",&Pnu1z_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_0_nu_tau_PE",&Enu2_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_0_nu_tau_PX",&Pnu2x_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_0_nu_tau_PY",&Pnu2y_s);
    t_kstar->SetBranchAddress("B0_kttdtf0_tauminus_0_nu_tau_PZ",&Pnu2z_s);

    t_kstar->SetBranchAddress("Taup_TRUEP_E",&taup_TRUEPE_s);
    t_kstar->SetBranchAddress("Taup_TRUEP_X",&taup_TRUEPX_s);
    t_kstar->SetBranchAddress("Taup_TRUEP_Y",&taup_TRUEPY_s);
    t_kstar->SetBranchAddress("Taup_TRUEP_Z",&taup_TRUEPZ_s);
    t_kstar->SetBranchAddress("Taum_TRUEP_E",&taum_TRUEPE_s);
    t_kstar->SetBranchAddress("Taum_TRUEP_X",&taum_TRUEPX_s);
    t_kstar->SetBranchAddress("Taum_TRUEP_Y",&taum_TRUEPY_s);
    t_kstar->SetBranchAddress("Taum_TRUEP_Z",&taum_TRUEPZ_s);

    t_kstar->SetBranchAddress("Taup_pi1_TRUEP_E",&taup_pip0_TRUEPE_s);
    t_kstar->SetBranchAddress("Taup_pi1_TRUEP_X",&taup_pip0_TRUEPX_s);
    t_kstar->SetBranchAddress("Taup_pi1_TRUEP_Y",&taup_pip0_TRUEPY_s);
    t_kstar->SetBranchAddress("Taup_pi1_TRUEP_Z",&taup_pip0_TRUEPZ_s);
    t_kstar->SetBranchAddress("Taup_pi2_TRUEP_E",&taup_pim0_TRUEPE_s);
    t_kstar->SetBranchAddress("Taup_pi2_TRUEP_X",&taup_pim0_TRUEPX_s);
    t_kstar->SetBranchAddress("Taup_pi2_TRUEP_Y",&taup_pim0_TRUEPY_s);
    t_kstar->SetBranchAddress("Taup_pi2_TRUEP_Z",&taup_pim0_TRUEPZ_s);
    t_kstar->SetBranchAddress("Taup_pi3_TRUEP_E",&taup_pip1_TRUEPE_s);
    t_kstar->SetBranchAddress("Taup_pi3_TRUEP_X",&taup_pip1_TRUEPX_s);
    t_kstar->SetBranchAddress("Taup_pi3_TRUEP_Y",&taup_pip1_TRUEPY_s);
    t_kstar->SetBranchAddress("Taup_pi3_TRUEP_Z",&taup_pip1_TRUEPZ_s);

    t_kstar->SetBranchAddress("Taum_pi1_TRUEP_E",&taum_pim0_TRUEPE_s);
    t_kstar->SetBranchAddress("Taum_pi1_TRUEP_X",&taum_pim0_TRUEPX_s);
    t_kstar->SetBranchAddress("Taum_pi1_TRUEP_Y",&taum_pim0_TRUEPY_s);
    t_kstar->SetBranchAddress("Taum_pi1_TRUEP_Z",&taum_pim0_TRUEPZ_s);
    t_kstar->SetBranchAddress("Taum_pi2_TRUEP_E",&taum_pip0_TRUEPE_s);
    t_kstar->SetBranchAddress("Taum_pi2_TRUEP_X",&taum_pip0_TRUEPX_s);
    t_kstar->SetBranchAddress("Taum_pi2_TRUEP_Y",&taum_pip0_TRUEPY_s);
    t_kstar->SetBranchAddress("Taum_pi2_TRUEP_Z",&taum_pip0_TRUEPZ_s);
    t_kstar->SetBranchAddress("Taum_pi3_TRUEP_E",&taum_pim1_TRUEPE_s);
    t_kstar->SetBranchAddress("Taum_pi3_TRUEP_X",&taum_pim1_TRUEPX_s);
    t_kstar->SetBranchAddress("Taum_pi3_TRUEP_Y",&taum_pim1_TRUEPY_s);
    t_kstar->SetBranchAddress("Taum_pi3_TRUEP_Z",&taum_pim1_TRUEPZ_s);

    TH1D* h_nutau_E = new TH1D("h_nutau_E", "h_nutau_E", 80, -120000., 120000);
    TH1D* h_antinutau_E = new TH1D("h_antinutau_E", "h_antinutau_E", 80, -120000., 120000);
    TH1D* h_nutau_x = new TH1D("h_nutau_x", "h_nutau_x", 80, -25000., 25000);
    TH1D* h_antinutau_x = new TH1D("h_antinutau_x", "h_antinutau_x", 80, -25000., 25000);
    TH1D* h_nutau_y = new TH1D("h_nutau_y", "h_nutau_y", 80, -25000., 25000);
    TH1D* h_antinutau_y = new TH1D("h_antinutau_y", "h_antinutau_y", 80, -25000., 25000);
    TH1D* h_nutau_z = new TH1D("h_nutau_z", "h_nutau_z", 80, -120000., 120000);
    TH1D* h_antinutau_z = new TH1D("h_antinutau_z", "h_antinutau_z", 80, -120000., 120000);
    TH1D* h_nutau_z_core = new TH1D("h_nutau_z_core", "h_nutau_z_core", 80, -120000., 120000);
    TH1D* h_antinutau_z_core = new TH1D("h_antinutau_z_core", "h_antinutau_z_core", 80, -120000., 120000);
    TH1D* h_nutau_z_tail = new TH1D("h_nutau_z_tail", "h_nutau_z_tail", 80, -120000., 120000);
    TH1D* h_antinutau_z_tail = new TH1D("h_antinutau_z_tail", "h_antinutau_z_tail", 80, -120000., 120000);

    TH1D* h_nutau_E_s = new TH1D("h_nutau_E_s", "h_nutau_E_s", 80, -120000., 120000);
    TH1D* h_antinutau_E_s = new TH1D("h_antinutau_E_s", "h_antinutau_E_s", 80, -120000., 120000);
    TH1D* h_nutau_x_s = new TH1D("h_nutau_x_s", "h_nutau_x_s", 80, -25000., 25000);
    TH1D* h_antinutau_x_s = new TH1D("h_antinutau_x_s", "h_antinutau_x_s", 80, -25000., 25000);
    TH1D* h_nutau_y_s = new TH1D("h_nutau_y_s", "h_nutau_y_s", 80, -25000., 25000);
    TH1D* h_antinutau_y_s = new TH1D("h_antinutau_y_s", "h_antinutau_y_s", 80, -25000., 25000);
    TH1D* h_nutau_z_s = new TH1D("h_nutau_z_s", "h_nutau_z_s", 80, -25000., 25000);
    TH1D* h_antinutau_z_s = new TH1D("h_antinutau_z_s", "h_antinutau_z_s", 80, -120000., 120000);

    TH1D* h_mass = new TH1D("h_mass", "h_mass", 80, 4000, 8000);
    TH1D* h_pull = new TH1D("h_pull", "h_pull", 80, -10., 10);
    TH1D* h_pull_taup = new TH1D("h_pull_taup", "h_pull_taup", 80, -5, 5);
    TH1D* h_pull_taum = new TH1D("h_pull_taum", "h_pull_taum", 80, -5, 5);

    TH1D* h_mass_best_cut = new TH1D("h_mass_best_cut", "h_mass_best_cut", 80, 4000, 8000);
    TH1D* h_mass_2ndbest_cut = new TH1D("h_mass_2ndbest_cut", "h_mass_2ndbest_cut", 80, 4000, 8000);
    TH1D* h_mass_s = new TH1D("h_mass_s", "h_mass_s", 80, 4000, 8000);
    TH1D* h_mass_new = new TH1D("h_mass_new", "h_mass_new", 80, 4000, 8000);

    TH1D* h_angle_core = new TH1D("h_angle_core", "h_angle_core", 40, 0., 10);
    TH1D* h_angle_tail = new TH1D("h_angle_tail", "h_angle_tail", 40, 0, 10);
    TH1D* h_angle_fail = new TH1D("h_angle_fail", "h_angle_fail", 40, 0, 10);

    TH2D* h_corr1_PX_PZ = new TH2D("h_corr1_PX_PZ", "h_corr1_PX_PZ", 200, -4000, 4000, 200, 0, 40000);
    TH2D* h_corr2_PX_PZ = new TH2D("h_corr2_PX_PZ", "h_corr2_PX_PZ", 200, -4000, 4000, 200, 0, 40000);
    TH2D* h_corr1_PY_PZ = new TH2D("h_corr1_PY_PZ", "h_corr1_PY_PZ", 200, -4000, 4000, 200, 0, 40000);
    TH2D* h_corr2_PY_PZ = new TH2D("h_corr2_PY_PZ", "h_corr2_PY_PZ", 200, -4000, 4000, 200, 0, 40000);
    TH2D* h_corr1_M_PZ = new TH2D("h_corr1_M_PZ", "h_corr1_M_PZ", 200, 4000, 8000, 200, 0, 40000);
    TH2D* h_corr2_M_PZ = new TH2D("h_corr2_M_PZ", "h_corr2_M_PZ", 200, 4000, 8000, 200, 0, 40000);

    TH2D* h_taup_FD_nuPz_core = new TH2D("h_taup_FD_nuPz_core", "h_taup_FD_nuPz_core", 200, 0, 10, 200, -20000, 50000);
    TH2D* h_taum_FD_nuPz_core = new TH2D("h_taum_FD_nuPz_core", "h_taum_FD_nuPz_core", 200, 0, 10, 200, -20000, 50000);
    TH2D* h_taup_FD_nuPz_tail = new TH2D("h_taup_FD_nuPz_tail", "h_taup_FD_nuPz_tail", 200, 0, 10, 200, -20000, 50000);
    TH2D* h_taum_FD_nuPz_tail = new TH2D("h_taum_FD_nuPz_tail", "h_taum_FD_nuPz_tail", 200, 0, 10, 200, -20000, 50000);

    TH2D* h_m_merr_pass = new TH2D("h_m_merr_pass", "h_m_merr_pass", 80, 0, 1000, 40, 4000, 8000);
    TH2D* h_m_merr_fail = new TH2D("h_m_merr_fail", "h_m_merr_fail", 80, 0, 1000, 40, 4000, 8000);

    TH2D* h_pull_pnu1z = new TH2D("h_pull_pnu1z", "h_pull_pnu1z", 80, -10, 10, 80, -120000, 120000);
    TH2D* h_pull_pnu2z = new TH2D("h_pull_pnu2z", "h_pull_pnu2z", 80, -10, 10, 80, -120000, 120000);

    TH2D* h_DTF_taup_taum_M = new TH2D("h_DTF_taup_taum_M", "h_DTF_taup_taum_M", 80, -5000, 5000, 80, -5000, 5000);
    TH2D* h_DTF_taup_Bmass_M = new TH2D("h_DTF_taup_Bmass_M", "h_DTF_taup_Bmass_M", 80, 4000, 8000, 80, -1000, 3000);
    TH2D* h_DTF_taum_Bmass_M = new TH2D("h_DTF_taum_Bmass_M", "h_DTF_taum_Bmass_M", 80, 4000, 8000, 80, -1000, 3000);

    TH1D* h_DTF_pass_mass = new TH1D("h_DTF_pass_mass", "h_DTF_pass_mass", 80, 2000, 8000);
    TH1D* h_DTF_fail_mass = new TH1D("h_DTF_fail_mass", "h_DTF_fail_mass", 80, 2000, 8000);
    TH1D* h_visible_mass = new TH1D("h_visible_mass", "h_visible_mass", 80, 2000, 8000);


    double P1;
    double P2;

    double nutau_TRUEPE;
    double antinutau_TRUEPE;
    double nutau_TRUEPX;
    double antinutau_TRUEPX;
    double nutau_TRUEPY;
    double antinutau_TRUEPY;
    double nutau_TRUEPZ;
    double antinutau_TRUEPZ;

    RooRealVar mass("mass", "mass", 4000, 8000, "MeV");
    RooDataSet* data1 = new RooDataSet("data1", "data1", mass);
    RooDataSet* data2 = new RooDataSet("data2", "data2", mass);
    RooDataSet* data3 = new RooDataSet("data3", "data3", mass);
    RooDataSet* data4 = new RooDataSet("data4", "data4", mass);
    RooDataSet* data5 = new RooDataSet("data5", "data5", mass);
    RooDataSet* data6 = new RooDataSet("data6", "data6", mass);
    RooDataSet* data7 = new RooDataSet("data7", "data7", mass);
    RooDataSet* data7_1 = new RooDataSet("data7_1", "data7_1", mass);
    RooDataSet* data8 = new RooDataSet("data8", "data8", mass);
    RooDataSet* data9 = new RooDataSet("data9", "data9", mass);
    RooDataSet* data10 = new RooDataSet("data10", "data10", mass);
    RooDataSet* data10_1 = new RooDataSet("data10_1", "data10_1", mass);
    RooDataSet* data11 = new RooDataSet("data11", "data11", mass);
    RooDataSet* data12 = new RooDataSet("data12", "data12", mass);
    RooDataSet* data13 = new RooDataSet("data13", "data13", mass);
    RooDataSet* data14 = new RooDataSet("data14", "data14", mass);

    for(int i = 0; i<t->GetEntries(); i++){
        t->GetEntry(i);
        mass.setVal(Bmass);

        // Triangle area 
        ROOT::Math::XYZPoint PV(PVx,PVy,PVz);
        ROOT::Math::XYZPoint DV1(DV1x,DV1y,DV1z);
        ROOT::Math::XYZPoint DV2(DV2x,DV2y,DV2z);

        double a = sqrt((PV-DV1).Mag2());
        double b = sqrt((PV-DV2).Mag2());
        double c = sqrt((DV2-DV1).Mag2());
        double s = (a+b+c)/2.;
        double A = sqrt(s*(s-a)*(s-b)*(s-c));

        // IPs
        ROOT::Math::XYZPoint KV(KVx, KVy, KVz);
        ROOT::Math::XYZVector Pk(Pkx, Pky, Pkz);

        ROOT::Math::XYZVector KV1(KV.x() + Pk.x(), KV.y() + Pk.y(), KV.z() + Pk.z()); 
        ROOT::Math::XYZVector A1(DV1.x() - KV.x(), DV1.y() - KV.y(), DV1.z() - KV.z());
        ROOT::Math::XYZVector A2(DV2.x() - KV.x(), DV2.y() - KV.y(), DV2.z() - KV.z());
        ROOT::Math::XYZVector A3(PV.x() - KV.x(), PV.y()-KV.y(), PV.z()-KV.z());

        double IP1 = (2*sqrt( (A1.Cross(Pk)).Mag2() ))/sqrt(Pk.Mag2());
        double IP2 = (2*sqrt( (A2.Cross(Pk)).Mag2() ))/sqrt(Pk.Mag2());
        double IP3 = (2*sqrt( (A3.Cross(Pk)).Mag2() ))/sqrt(Pk.Mag2());

        // Angle between K+ trajectory and PV-DV1-DV2 triangle
        ROOT::Math::XYZVector v1 = PV-DV1;
        ROOT::Math::XYZVector v2 = DV2-DV1;
        ROOT::Math::XYZVector n = v2.Cross(v1);

        double theta = acos( (n.Dot(Pk))/(sqrt(n.Mag2())*sqrt(Pk.Mag2())) );

        // Magnitude of DTF neutrino momentum
        P1 = sqrt(pow(Pnu1x,2) + pow(Pnu1y,2) + pow(Pnu1z,2));
        P2 = sqrt(pow(Pnu2x,2) + pow(Pnu2y,2) + pow(Pnu2z,2));

        // DTF vs TRUE nu momentum
        nutau_TRUEPE = taum_TRUEPE - (taum_pim0_TRUEPE + taum_pip0_TRUEPE + taum_pim1_TRUEPE);
        antinutau_TRUEPE = taup_TRUEPE - (taup_pip0_TRUEPE + taup_pim0_TRUEPE + taup_pip1_TRUEPE);
        nutau_TRUEPX = taum_TRUEPX - (taum_pim0_TRUEPX + taum_pip0_TRUEPX + taum_pim1_TRUEPX);
        antinutau_TRUEPX = taup_TRUEPX - (taup_pip0_TRUEPX + taup_pim0_TRUEPX + taup_pip1_TRUEPX);
        nutau_TRUEPY = taum_TRUEPY - (taum_pim0_TRUEPY + taum_pip0_TRUEPY + taum_pim1_TRUEPY);
        antinutau_TRUEPY = taup_TRUEPY - (taup_pip0_TRUEPY + taup_pim0_TRUEPY + taup_pip1_TRUEPY);
        nutau_TRUEPZ = taum_TRUEPZ - (taum_pim0_TRUEPZ + taum_pip0_TRUEPZ + taum_pim1_TRUEPZ);
        antinutau_TRUEPZ = taup_TRUEPZ - (taup_pip0_TRUEPZ + taup_pim0_TRUEPZ + taup_pip1_TRUEPZ);

        // New Bmass (with all DTF kinematics except for neutrino Pz, which is the TRUE value)
        DTF_PBx -= (Pnu1x + Pnu2x);
        DTF_PBx += (nutau_TRUEPX + antinutau_TRUEPX);
        DTF_PBy -= (Pnu1y + Pnu2y);
        DTF_PBy += (nutau_TRUEPY + antinutau_TRUEPY);
        DTF_PBz -= (Pnu1z + Pnu2z);
        DTF_PBz += (nutau_TRUEPZ + antinutau_TRUEPZ);

        DTF_EB -= (Enu1 + Enu2);
        DTF_EB += ( sqrt(pow(antinutau_TRUEPX,2) + pow(antinutau_TRUEPY,2) + pow(antinutau_TRUEPZ,2)) + sqrt(pow(nutau_TRUEPX,2) + pow(nutau_TRUEPY,2) + pow(nutau_TRUEPZ,2)) );

        // Tau separation in XY plane
        double Dxy = sqrt(pow(DV1x - DV2x,2) + pow(DV1y - DV2y,2));

        double Bmass_new = sqrt(pow(DTF_EB,2) - pow(DTF_PBx,2) - pow(DTF_PBy,2) - pow(DTF_PBz,2));

        if( (abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)){

            h_visible_mass->Fill(Bmass_visible);

            if(status == 0){
                h_mass->Fill(Bmass);
                h_mass_new->Fill(Bmass_new);
                h_DTF_pass_mass->Fill(Bmass);

                h_nutau_E->Fill(Enu2 - nutau_TRUEPE);
                h_antinutau_E->Fill(Enu1 - antinutau_TRUEPE);
                h_nutau_x->Fill(Pnu2x - nutau_TRUEPX);
                h_antinutau_x->Fill(Pnu1x - antinutau_TRUEPX);
                h_nutau_y->Fill(Pnu2y - nutau_TRUEPY);
                h_antinutau_y->Fill(Pnu1y - antinutau_TRUEPY);
                h_nutau_z->Fill(Pnu2z - nutau_TRUEPZ);
                h_antinutau_z->Fill(Pnu1z - antinutau_TRUEPZ);

                h_corr1_PX_PZ->Fill(Pnu1x, Pnu1z);
                h_corr2_PX_PZ->Fill(Pnu2x, Pnu2z);
                h_corr1_PY_PZ->Fill(Pnu1y, Pnu1z);
                h_corr2_PY_PZ->Fill(Pnu2y, Pnu2z);
                h_corr1_M_PZ->Fill(Bmass, Pnu1z);
                h_corr2_M_PZ->Fill(Bmass, Pnu2z);

                h_pull->Fill((Bmass - 5279)/Bmass_err);
                h_pull_taup->Fill((taup_DTF_M-1776.86)/taup_DTF_Merr);
                h_pull_taum->Fill((taum_DTF_M-1776.86)/taum_DTF_Merr);

                h_pull_pnu1z->Fill((Bmass - 5279)/Bmass_err,Pnu1z-antinutau_TRUEPZ);
                h_pull_pnu2z->Fill((Bmass - 5279)/Bmass_err,Pnu2z-nutau_TRUEPZ);

                h_m_merr_pass->Fill(Bmass_err,Bmass);

                h_DTF_taup_taum_M->Fill(taup_DTF_M, taum_DTF_M);
                h_DTF_taup_Bmass_M->Fill(Bmass, taup_DTF_M);
                h_DTF_taum_Bmass_M->Fill(Bmass, taum_DTF_M);

                if(abs(Bmass - 5279) < 300){
                    h_nutau_z_core->Fill(Pnu2z - nutau_TRUEPZ);
                    h_antinutau_z_core->Fill(Pnu1z - antinutau_TRUEPZ);

                    h_angle_core->Fill(abs((0.5*M_PI-theta)*(180/M_PI)));

                    h_taup_FD_nuPz_core->Fill(taup_FD, Pnu1z);
                    h_taum_FD_nuPz_core->Fill(taum_FD, Pnu2z);

                }
                if(Bmass > 6500){
                    h_nutau_z_tail->Fill(Pnu2z - nutau_TRUEPZ);
                    h_antinutau_z_tail->Fill(Pnu1z - antinutau_TRUEPZ);

                    h_angle_tail->Fill(abs((0.5*M_PI-theta)*(180/M_PI)));

                    h_taup_FD_nuPz_tail->Fill(taup_FD, Pnu1z);
                    h_taum_FD_nuPz_tail->Fill(taum_FD, Pnu2z);

                }

                if((Bmass > 4000) && (Bmass < 8000)){
                    data1->add(mass);

                    if(A > 1.){data2->add(mass);}
                    if((IP1 > 0.4) && (IP2 > 0.4)){data3->add(mass);}
                    if( ((taup_DTF_M > 1750) && (taup_DTF_M < 1800)) && ((taum_DTF_M > 1750) && (taum_DTF_M < 1800)) ){
                        data4->add(mass);
                        h_mass_2ndbest_cut->Fill(Bmass);
                    }
                    if(((P1 < 5000) && (P2 < 5000))){
                        data5->add(mass);
                        h_mass_best_cut->Fill(Bmass);
                    }
                    if((taup_FD > 1.) && (taum_FD > 1.)){data6->add(mass);}
                    if((taup_FD_chi2 > 4) && (taum_FD_chi2 > 4)){data7->add(mass);}
                    if((taup_FD_chi2 > 9) && (taum_FD_chi2 > 9)){data7_1->add(mass);}
                    if(DTF_chi2 < 15){data8->add(mass);}
                    if(DTF_nIter < 5){data9->add(mass);}
                    if(abs(DV1z - DV2z) > 1.5){data10->add(mass);}
                    if(Dxy > 0.15){data10_1->add(mass);}
                    if(Kp_PT > 2000){data11->add(mass);}
                    if((taum_PT > 3000) && (taup_PT > 3000)){data12->add(mass);}
                    if((DV1z > (BVz - 0.2) && (DV2z > (BVz - 0.2)))){data13->add(mass);}
                    if(abs(abs((0.5*M_PI-theta)*(180/M_PI))) > 0.3 ){data14->add(mass);}

                }

            }
            else{
                h_DTF_fail_mass->Fill(Bmass);
                h_angle_fail->Fill(abs((0.5*M_PI-theta)*(180/M_PI)));

                h_m_merr_fail->Fill(Bmass_err,Bmass);
            }
        }
    }

    RooRealVar mass_kstar("mass_kstar", "mass_kstar", 4000, 8000, "MeV");
    RooDataSet* data_kstar = new RooDataSet("data_kstar", "data_kstar", mass_kstar);

    double nutau_TRUEPE_s;
    double antinutau_TRUEPE_s;
    double nutau_TRUEPX_s;
    double antinutau_TRUEPX_s;
    double nutau_TRUEPY_s;
    double antinutau_TRUEPY_s;
    double nutau_TRUEPZ_s;
    double antinutau_TRUEPZ_s;

    for(int i = 0; i < t_kstar->GetEntries(); i++){
        t_kstar->GetEntry(i);

        // DTF vs TRUE nu momentum
        nutau_TRUEPE_s = taum_TRUEPE_s - (taum_pim0_TRUEPE_s + taum_pip0_TRUEPE_s + taum_pim1_TRUEPE_s);
        antinutau_TRUEPE_s = taup_TRUEPE_s - (taup_pip0_TRUEPE_s + taup_pim0_TRUEPE_s + taup_pip1_TRUEPE_s);
        nutau_TRUEPX_s = taum_TRUEPX_s - (taum_pim0_TRUEPX_s + taum_pip0_TRUEPX_s + taum_pim1_TRUEPX_s);
        antinutau_TRUEPX_s = taup_TRUEPX_s - (taup_pip0_TRUEPX_s + taup_pim0_TRUEPX_s + taup_pip1_TRUEPX_s);
        nutau_TRUEPY_s = taum_TRUEPY_s - (taum_pim0_TRUEPY_s + taum_pip0_TRUEPY_s + taum_pim1_TRUEPY_s);
        antinutau_TRUEPY_s = taup_TRUEPY_s - (taup_pip0_TRUEPY_s + taup_pim0_TRUEPY_s + taup_pip1_TRUEPY_s);
        nutau_TRUEPZ_s = taum_TRUEPZ_s - (taum_pim0_TRUEPZ_s + taum_pip0_TRUEPZ_s + taum_pim1_TRUEPZ_s);
        antinutau_TRUEPZ_s = taup_TRUEPZ_s - (taup_pip0_TRUEPZ_s + taup_pim0_TRUEPZ_s + taup_pip1_TRUEPZ_s);

        if(status_kstar == 0){

            h_mass_s->Fill(M_kstar);

            if((M_kstar > 4000) && (M_kstar < 8000)){
                mass_kstar.setVal(M_kstar);
                data_kstar->add(mass_kstar);

            h_nutau_E_s->Fill(Enu2_s - nutau_TRUEPE_s);
            h_antinutau_E_s->Fill(Enu1_s - antinutau_TRUEPE_s);
            h_nutau_x_s->Fill(Pnu2x_s - nutau_TRUEPX_s);
            h_antinutau_x_s->Fill(Pnu1x_s - antinutau_TRUEPX_s);
            h_nutau_y_s->Fill(Pnu2y_s - nutau_TRUEPY_s);
            h_antinutau_y_s->Fill(Pnu1y_s - antinutau_TRUEPY_s);
            h_nutau_z_s->Fill(Pnu2z_s - nutau_TRUEPZ_s);
            h_antinutau_z_s->Fill(Pnu1z_s - antinutau_TRUEPZ_s);

            }
        }
    }

    if(make_fit){

        RooRealVar* sigma1 = fit(data1, mass, "No cut", "Bmass_cut1", folder_name);
        RooRealVar* sigma2 = fit(data2, mass, "Area of PV-DV1-DV2 triangle > 1 mm^{2}", "Bmass_cut2", folder_name);
        RooRealVar* sigma3 = fit(data3, mass, "IP_{#tau} > 0.2 mm", "Bmass_cut3", folder_name);
        RooRealVar* sigma4 = fit(data4, mass, "DTF #tau mass #in [1750,1800] MeV", "Bmass_cut4", folder_name);
        RooRealVar* sigma5 = fit(data5, mass, "DTF P_{#nu} < 5000 MeV", "Bmass_cut5", folder_name);
        RooRealVar* sigma6 = fit(data6, mass, "#tau FD to BV > 1 mm", "Bmass_cut6", folder_name);
        RooRealVar* sigma7 = fit(data7, mass, "#tau FD to BV #chi^{2} > 4", "Bmass_cut7", folder_name);
        RooRealVar* sigma7_1 = fit(data7_1, mass, "#tau FD to BV #chi^{2} > 9", "Bmass_cut7_1", folder_name);
        RooRealVar* sigma8 = fit(data8, mass, "DTF #chi^{2} < 15", "Bmass_cut8", folder_name);
        RooRealVar* sigma9 = fit(data9, mass, "Number of DTF iterations < 5", "Bmass_cut9", folder_name);
        RooRealVar* sigma10 = fit(data10, mass, "#tau Z separation > 1.5 mm", "Bmass_cut10", folder_name);
        RooRealVar* sigma10_1 = fit(data10_1, mass, "#tau XY separation > 0.15 mm", "Bmass_cut10_1", folder_name);
        RooRealVar* sigma11 = fit(data11, mass, "K^{+} p_{T} > 2000 MeV", "Bmass_cut11", folder_name);
        RooRealVar* sigma12 = fit(data12, mass, "#tau p_{T} > 3000 MeV", "Bmass_cut12", folder_name);
        RooRealVar* sigma13 = fit(data13, mass, "#tau DVz > (BVz - 0.2 mm)", "Bmass_cut13", folder_name);
        RooRealVar* sigma14 = fit(data14, mass, "abs(#theta) > 0.45 #circ", "Bmass_cut14", folder_name);
        RooRealVar* sigma_kstar = fit(data_kstar, mass_kstar, "K^{*} mass", "Kstar_mass", folder_name);

        // FWHM + efficiency
        TString out = folder_name+"fwhm_eff_resolution.tex";
        file.open(out);  

        if(!file.is_open())
            cout << "Output file not opened!";

        file << "\\begin{table}" << std::endl;
        file << Form("\\caption{Resolution and efficiency of $B^{+}$ mass distribution for different cuts. Initial mass resolution is %.2lf $\\pm$ %.2lf MeV.}", sigma1->getVal(), sigma1->getError()) << std::endl;
        file << "\\centering" << std::endl;
        file << "\\begin{tabular}{|c|c|c|}" << std::endl;
        file << "\\hline" << std::endl; 

        double n1 = data1->sumEntries();
        double n2 = data2->sumEntries();
        double n3 = data3->sumEntries();
        double n4 = data4->sumEntries();
        double n5 = data5->sumEntries();
        double n6 = data6->sumEntries();
        double n7 = data7->sumEntries();
        double n7_1 = data7_1->sumEntries();
        double n8 = data8->sumEntries();
        double n9 = data9->sumEntries();
        double n10 = data10->sumEntries();
        double n10_1 = data10_1->sumEntries();
        double n11 = data11->sumEntries();
        double n12 = data12->sumEntries();
        double n13 = data13->sumEntries();
        double n14 = data14->sumEntries();

        file << "Cut & " <<  "Resolution (MeV) & " << "Efficiency (\\%)" << "\\\\ \\hline" << std::endl;   
        file << "Area of PV-DV1-DV2 triangle $> 1\\,$mm$^2$ & " <<  Form("%.2lf $\\pm$ %.2lf & ", sigma2->getVal(), sigma2->getError()) << Form("%.2lf ", (n2/n1)*100.) << "\\\\" << std::endl;    
        file << "$IP_{\\tau^{\\pm}} > 0.2$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma3->getVal(), sigma3->getError()) << Form("%.2lf ", (n3/n1)*100.) << "\\\\" << std::endl;    
        file << "DTF $\\tau^{\\pm}$ mass $\\in [1750,1800]$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma4->getVal(), sigma4->getError()) << Form("%.2lf ", (n4/n1)*100.) << "\\\\" << std::endl;    
        file << "DTF $P_{\\nu} < 5000$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma5->getVal(), sigma5->getError()) << Form("%.2lf ", (n5/n1)*100.) << "\\\\" << std::endl;    
        file << "$\\tau^{\\pm}$ FD to BV$ > 1$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma6->getVal(), sigma6->getError()) << Form("%.2lf ", (n6/n1)*100.) << "\\\\" << std::endl;    
        file << "$\\tau^{\\pm}$ FD to BV $\\chi^2 > 4$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma7->getVal(), sigma7->getError()) << Form("%.2lf ", (n7/n1)*100.) << "\\\\" << std::endl;    
        file << "$\\tau^{\\pm}$ FD to BV $\\chi^2 > 9$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma7_1->getVal(), sigma7_1->getError()) << Form("%.2lf ", (n7_1/n1)*100.) << "\\\\" << std::endl;    
        file << "DTF $\\chi^2 < 15$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma8->getVal(), sigma8->getError()) << Form("%.2lf ", (n8/n1)*100.) << "\\\\" << std::endl;    
        file << "Number of DTF iterations $< 5$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma9->getVal(), sigma9->getError()) << Form("%.2lf ", (n9/n1)*100.) << "\\\\" << std::endl;    
        file << "$\\tau$ z separation $> 1.5$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma10->getVal(), sigma10->getError()) << Form("%.2lf ", (n10/n1)*100.) << "\\\\" << std::endl;    
        file << "$\\tau$ xy separation $> 0.15$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma10_1->getVal(), sigma10_1->getError()) << Form("%.2lf ", (n10_1/n1)*100.) << "\\\\" << std::endl;    
        file << "$K^{+} p_{T} > 2000$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma11->getVal(), sigma11->getError()) << Form("%.2lf ", (n11/n1)*100.) << "\\\\" << std::endl;    
        file << "$\\tau^{\\pm} p_{T} > 3000$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma12->getVal(), sigma12->getError()) << Form("%.2lf ", (n12/n1)*100.) << "\\\\" << std::endl;  
        file << "$\\tau^{\\pm}$ DVz $>$ (BVz-0.2\\,mm) & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma13->getVal(), sigma13->getError()) << Form("%.2lf ", (n13/n1)*100.) << "\\\\" << std::endl;  
        file << "$abs(\\theta) > 0.45 ^{\\circ}$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma14->getVal(), sigma14->getError()) << Form("%.2lf ", (n14/n1)*100.) << "\\\\ \\hline" << std::endl; 

        file << "\\end{tabular}" << std::endl;
        file << "\\end{table}" << std::endl;
        file.close();
    }   

    TCanvas c;
    c.cd();
    h_mass_new->GetXaxis()->SetTitle("mass_new (MeV)");
    h_mass_new->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass_new->GetXaxis()->GetXmax() - h_mass_new->GetXaxis()->GetXmin())/h_mass_new->GetNbinsX()) );
    h_mass_new->SetTitle(" ");
    h_mass_new->Draw();
    c.SaveAs(folder_name+"mass_new.gif");
    c.SaveAs(folder_name+"mass_new.pdf");

    TCanvas c1E;
    c1E.cd();
    h_nutau_E->GetXaxis()->SetTitle("nutau_DTF_PE - nutau_TRUE_PE (MeV)");
    h_nutau_E->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_E->GetXaxis()->GetXmax() - h_nutau_E->GetXaxis()->GetXmin())/h_nutau_E->GetNbinsX()) );
    h_nutau_E->SetTitle(" ");
    h_nutau_E->Draw();
    c1E.SaveAs(folder_name+"nutau_E_DTF_TRUE.gif");
    c1E.SaveAs(folder_name+"nutau_E_DTF_TRUE.pdf");

    TCanvas c2E;
    c2E.cd();
    h_antinutau_E->GetXaxis()->SetTitle("antinutau_DTF_PE - antinutau_TRUE_PE (MeV)");
    h_antinutau_E->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_E->GetXaxis()->GetXmax() - h_antinutau_E->GetXaxis()->GetXmin())/h_antinutau_E->GetNbinsX()) );
    h_antinutau_E->SetTitle(" ");
    h_antinutau_E->Draw();
    c2E.SaveAs(folder_name+"antinutau_E_DTF_TRUE.gif");
    c2E.SaveAs(folder_name+"antinutau_E_DTF_TRUE.pdf");

    TCanvas c3E;
    c3E.cd();
    h_nutau_E_s->GetXaxis()->SetTitle("nutau_DTF_PE - nutau_TRUE_PE (MeV)");
    h_nutau_E_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_E_s->GetXaxis()->GetXmax() - h_nutau_E_s->GetXaxis()->GetXmin())/h_nutau_E_s->GetNbinsX()) );
    h_nutau_E_s->SetTitle(" ");
    h_nutau_E_s->Draw();
    c3E.SaveAs(folder_name+"nutau_E_DTF_TRUE_s.gif");
    c3E.SaveAs(folder_name+"nutau_E_DTF_TRUE_s.pdf");

    TCanvas c4E;
    c4E.cd();
    h_antinutau_E_s->GetXaxis()->SetTitle("antinutau_DTF_PE - antinutau_TRUE_PE (MeV)");
    h_antinutau_E_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_E_s->GetXaxis()->GetXmax() - h_antinutau_E_s->GetXaxis()->GetXmin())/h_antinutau_E_s->GetNbinsX()) );
    h_antinutau_E_s->SetTitle(" ");
    h_antinutau_E_s->Draw();
    c4E.SaveAs(folder_name+"antinutau_E_DTF_TRUE_s.gif");
    c4E.SaveAs(folder_name+"antinutau_E_DTF_TRUE_s.pdf");

    TCanvas c1x;
    c1x.cd();
    h_nutau_x->GetXaxis()->SetTitle("nutau_DTF_PX - nutau_TRUE_PX (MeV)");
    h_nutau_x->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_x->GetXaxis()->GetXmax() - h_nutau_x->GetXaxis()->GetXmin())/h_nutau_x->GetNbinsX()) );
    h_nutau_x->SetTitle(" ");
    h_nutau_x->Draw();
    c1x.SaveAs(folder_name+"nutau_x_DTF_TRUE.gif");
    c1x.SaveAs(folder_name+"nutau_x_DTF_TRUE.pdf");

    TCanvas c2x;
    c2x.cd();
    h_antinutau_x->GetXaxis()->SetTitle("antinutau_DTF_PX - antinutau_TRUE_PX (MeV)");
    h_antinutau_x->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_x->GetXaxis()->GetXmax() - h_antinutau_x->GetXaxis()->GetXmin())/h_antinutau_x->GetNbinsX()) );
    h_antinutau_x->SetTitle(" ");
    h_antinutau_x->Draw();
    c2x.SaveAs(folder_name+"antinutau_x_DTF_TRUE.gif");
    c2x.SaveAs(folder_name+"antinutau_x_DTF_TRUE.pdf");

    TCanvas c3x;
    c3x.cd();
    h_nutau_x_s->GetXaxis()->SetTitle("nutau_DTF_PX - nutau_TRUE_PX (MeV)");
    h_nutau_x_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_x_s->GetXaxis()->GetXmax() - h_nutau_x_s->GetXaxis()->GetXmin())/h_nutau_x_s->GetNbinsX()) );
    h_nutau_x_s->SetTitle(" ");
    h_nutau_x_s->Draw();
    c3x.SaveAs(folder_name+"nutau_x_DTF_TRUE_s.gif");
    c3x.SaveAs(folder_name+"nutau_x_DTF_TRUE_s.pdf");

    TCanvas c4x;
    c4x.cd();
    h_antinutau_x_s->GetXaxis()->SetTitle("antinutau_DTF_PX - antinutau_TRUE_PX (MeV)");
    h_antinutau_x_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_x_s->GetXaxis()->GetXmax() - h_antinutau_x_s->GetXaxis()->GetXmin())/h_antinutau_x_s->GetNbinsX()) );
    h_antinutau_x_s->SetTitle(" ");
    h_antinutau_x_s->Draw();
    c4x.SaveAs(folder_name+"antinutau_x_DTF_TRUE_s.gif");
    c4x.SaveAs(folder_name+"antinutau_x_DTF_TRUE_s.pdf");

    TCanvas c1y;
    c1y.cd();
    h_nutau_y->GetXaxis()->SetTitle("nutau_DTF_PY - nutau_TRUE_PY (MeV)");
    h_nutau_y->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_y->GetXaxis()->GetXmax() - h_nutau_y->GetXaxis()->GetXmin())/h_nutau_y->GetNbinsX()) );
    h_nutau_y->SetTitle(" ");
    h_nutau_y->Draw();
    c1y.SaveAs(folder_name+"nutau_y_DTF_TRUE.gif");
    c1y.SaveAs(folder_name+"nutau_y_DTF_TRUE.pdf");

    TCanvas c2y;
    c2y.cd();
    h_antinutau_y->GetXaxis()->SetTitle("antinutau_DTF_PY - antinutau_TRUE_PY (MeV)");
    h_antinutau_y->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_y->GetXaxis()->GetXmax() - h_antinutau_y->GetXaxis()->GetXmin())/h_antinutau_y->GetNbinsX()) );
    h_antinutau_y->SetTitle(" ");
    h_antinutau_y->Draw();
    c2y.SaveAs(folder_name+"antinutau_y_DTF_TRUE.gif");
    c2y.SaveAs(folder_name+"antinutau_y_DTF_TRUE.pdf");
    
    TCanvas c3y;
    c3y.cd();
    h_nutau_y_s->GetXaxis()->SetTitle("nutau_DTF_PY - nutau_TRUE_PY (MeV)");
    h_nutau_y_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_y_s->GetXaxis()->GetXmax() - h_nutau_y_s->GetXaxis()->GetXmin())/h_nutau_y_s->GetNbinsX()) );
    h_nutau_y_s->SetTitle(" ");
    h_nutau_y_s->Draw();
    c3y.SaveAs(folder_name+"nutau_y_DTF_TRUE_s.gif");
    c3y.SaveAs(folder_name+"nutau_y_DTF_TRUE_s.pdf");

    TCanvas c4y;
    c4y.cd();
    h_antinutau_y_s->GetXaxis()->SetTitle("antinutau_DTF_PY - antinutau_TRUE_PY (MeV)");
    h_antinutau_y_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_y_s->GetXaxis()->GetXmax() - h_antinutau_y_s->GetXaxis()->GetXmin())/h_antinutau_y_s->GetNbinsX()) );
    h_antinutau_y_s->SetTitle(" ");
    h_antinutau_y_s->Draw();
    c4y.SaveAs(folder_name+"antinutau_y_DTF_TRUE_s.gif");
    c4y.SaveAs(folder_name+"antinutau_y_DTF_TRUE_s.pdf");

    TCanvas c1z;
    c1z.cd();
    h_nutau_z->GetXaxis()->SetTitle("nutau_DTF_PZ - nutau_TRUE_PZ (MeV)");
    h_nutau_z->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_z->GetXaxis()->GetXmax() - h_nutau_z->GetXaxis()->GetXmin())/h_nutau_z->GetNbinsX()) );
    h_nutau_z->SetTitle(" ");
    h_nutau_z->Draw();
    c1z.SaveAs(folder_name+"nutau_z_DTF_TRUE.gif");
    c1z.SaveAs(folder_name+"nutau_z_DTF_TRUE.pdf");

    TCanvas c2z;
    c2z.cd();
    h_antinutau_z->GetXaxis()->SetTitle("antinutau_DTF_PZ - antinutau_TRUE_PZ (MeV)");
    h_antinutau_z->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_z->GetXaxis()->GetXmax() - h_antinutau_z->GetXaxis()->GetXmin())/h_antinutau_z->GetNbinsX()) );
    h_antinutau_z->SetTitle(" ");
    h_antinutau_z->Draw();
    c2z.SaveAs(folder_name+"antinutau_z_DTF_TRUE.gif");
    c2z.SaveAs(folder_name+"antinutau_z_DTF_TRUE.pdf");

    TCanvas c3z;
    c3z.cd();
    h_nutau_z_s->GetXaxis()->SetTitle("nutau_DTF_PZ - nutau_TRUE_PZ (MeV)");
    h_nutau_z_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_z_s->GetXaxis()->GetXmax() - h_nutau_z_s->GetXaxis()->GetXmin())/h_nutau_z_s->GetNbinsX()) );
    h_nutau_z_s->SetTitle(" ");
    h_nutau_z_s->Draw();
    c3z.SaveAs(folder_name+"nutau_z_DTF_TRUE_s.gif");
    c3z.SaveAs(folder_name+"nutau_z_DTF_TRUE_s.pdf");

    TCanvas c4z;
    c4z.cd();
    h_antinutau_z_s->GetXaxis()->SetTitle("antinutau_DTF_PZ - antinutau_TRUE_PZ (MeV)");
    h_antinutau_z_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_z_s->GetXaxis()->GetXmax() - h_antinutau_z_s->GetXaxis()->GetXmin())/h_antinutau_z_s->GetNbinsX()) );
    h_antinutau_z_s->SetTitle(" ");
    h_antinutau_z_s->Draw();
    c4z.SaveAs(folder_name+"antinutau_z_DTF_TRUE_s.gif");
    c4z.SaveAs(folder_name+"antinutau_z_DTF_TRUE_s.pdf");

    gStyle->SetOptStat(0);
    TCanvas c5z;
    c5z.cd();
    h_nutau_z_core->GetXaxis()->SetTitle("nutau_DTF_PZ - nutau_TRUE_PZ (MeV)");
    h_nutau_z_core->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_z_core->GetXaxis()->GetXmax() - h_nutau_z_core->GetXaxis()->GetXmin())/h_nutau_z_core->GetNbinsX()) );
    h_nutau_z_core->SetTitle(" ");
    h_nutau_z_core->SetLineColor(kBlue);
    h_nutau_z_tail->SetLineColor(kRed);
    //normalization
    h_nutau_z_core->Scale(1/h_nutau_z_core->Integral());
    h_nutau_z_tail->Scale(1/h_nutau_z_tail->Integral());    
    h_nutau_z_core->Draw("HIST");
    h_nutau_z_tail->Draw("HIST same");
    TLegend* leg1;
    leg1 = new TLegend(0.7, 0.8, 0.85, 0.85);
    leg1->SetBorderSize(0);
    leg1->AddEntry(h_nutau_z_core, "Core", "lp");
    leg1->AddEntry(h_nutau_z_tail, "Tail", "lp");
    leg1->SetTextSize(0.03);
    leg1->Draw("same");
    c5z.SaveAs(folder_name+"nutau_DTF_z_TRUE.gif");
    c5z.SaveAs(folder_name+"nutau_DTF_z_TRUE.pdf");

    TCanvas c6z;
    c6z.cd();
    h_antinutau_z_core->GetXaxis()->SetTitle("antinutau_DTF_PZ - antinutau_TRUE_PZ (MeV)");
    h_antinutau_z_core->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_z_core->GetXaxis()->GetXmax() - h_antinutau_z_core->GetXaxis()->GetXmin())/h_antinutau_z_core->GetNbinsX()) );
    h_antinutau_z_core->SetTitle(" ");
    h_antinutau_z_core->SetLineColor(kBlue);
    h_antinutau_z_tail->SetLineColor(kRed);
    //normalization
    h_antinutau_z_core->Scale(1/h_antinutau_z_core->Integral());
    h_antinutau_z_tail->Scale(1/h_antinutau_z_tail->Integral());
    h_antinutau_z_core->Draw("HIST");
    h_antinutau_z_tail->Draw("HIST same");
    TLegend* leg2;
    leg2 = new TLegend(0.7, 0.8, 0.85, 0.85);
    leg2->SetBorderSize(0);
    leg2->AddEntry(h_antinutau_z_core, "Core", "lp");
    leg2->AddEntry(h_antinutau_z_tail, "Tail", "lp");
    leg2->SetTextSize(0.03);
    leg2->Draw("same");
    c6z.SaveAs(folder_name+"antinutau_z_DTF_TRUE.gif");
    c6z.SaveAs(folder_name+"antinutau_z_DTF_TRUE.pdf");

    TCanvas c1;
    c1.cd();
    h_mass_s->SetTitle(" ");
    h_mass_s->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass_s->GetXaxis()->SetTitle("m_{B} (MeV)");
    h_mass->SetLineColor(kBlue);
    h_mass_s->SetLineColor(kRed);
    h_mass_best_cut->SetLineColor(kBlue);
    h_mass_s->DrawNormalized();
    //h_mass->DrawNormalized("same");
    h_mass_best_cut->DrawNormalized("same");
    TLegend* legs;
    legs = new TLegend(0.7, 0.8, 0.85, 0.85);
    legs->SetBorderSize(0);
    legs->AddEntry(h_mass, "KpTauTau Pnu < 5000 MeV", "lp");
    legs->AddEntry(h_mass_s, "KstarTauTau", "lp");
    legs->SetTextSize(0.03);
    legs->Draw("same");
    c1.SaveAs(folder_name+"Kplus_Kstar_mass.gif");
    c1.SaveAs(folder_name+"Kplus_Kstar_mass.pdf");

    TCanvas c2;
    c2.cd();
    h_angle_tail->GetXaxis()->SetTitle("|#theta| (#circ)");
    h_angle_tail->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_angle_core->GetXaxis()->GetXmax() - h_angle_core->GetXaxis()->GetXmin())/h_angle_core->GetNbinsX()) );
    h_angle_tail->SetTitle(" ");
    h_angle_core->SetLineColor(kBlack);
    h_angle_tail->SetLineColor(kBlue);
    h_angle_fail->SetLineColor(kRed);
    h_angle_tail->DrawNormalized();
    h_angle_fail->DrawNormalized("same");
    h_angle_core->DrawNormalized("same");
    TLegend* leg3;
    leg3 = new TLegend(0.6, 0.7, 0.85, 0.85);
    leg3->SetBorderSize(0);
    leg3->AddEntry(h_angle_core, Form("Core mean = %.3lf #circ",h_angle_core->GetMean()), "lp");
    leg3->AddEntry(h_angle_tail, Form("Tail mean = %.3lf #circ",h_angle_tail->GetMean()), "lp");
    leg3->AddEntry(h_angle_fail, Form("Fail mean = %.3lf #circ",h_angle_fail->GetMean()), "lp");
    leg3->SetTextSize(0.03);
    leg3->Draw("same");
    c2.SaveAs(folder_name+"angle.gif");
    c2.SaveAs(folder_name+"angle.pdf");

    TCanvas c3;
    c3.cd();
    h_corr1_PX_PZ->GetXaxis()->SetTitle("Antineutrino p_{x} (MeV)");
    h_corr1_PX_PZ->GetYaxis()->SetTitle("Antineutrino p_{z} (MeV)");
    h_corr1_PX_PZ->SetTitle("");
    h_corr1_PX_PZ->Draw("COLZ");
    c3.SaveAs(folder_name+"corr1_PX_PZ.gif");
    c3.SaveAs(folder_name+"corr1_PX_PZ.pdf");

    TCanvas c4;
    c4.cd();
    h_corr2_PX_PZ->GetXaxis()->SetTitle("Neutrino p_{x} (MeV)");
    h_corr2_PX_PZ->GetYaxis()->SetTitle("Neutrino p_{z} (MeV)");
    h_corr2_PX_PZ->SetTitle("");
    h_corr2_PX_PZ->Draw("COLZ");
    c4.SaveAs(folder_name+"corr2_PX_PZ.gif");
    c4.SaveAs(folder_name+"corr2_PX_PZ.pdf");

    TCanvas c5;
    c5.cd();
    h_corr1_PY_PZ->GetXaxis()->SetTitle("Antineutrino p_{y} (MeV)");
    h_corr1_PY_PZ->GetYaxis()->SetTitle("Antineutrino p_{z} (MeV)");
    h_corr1_PY_PZ->SetTitle("");
    h_corr1_PY_PZ->Draw("COLZ");
    c5.SaveAs(folder_name+"corr1_PY_PZ.gif");
    c5.SaveAs(folder_name+"corr1_PY_PZ.pdf");

    TCanvas c6;
    c6.cd();
    h_corr2_PY_PZ->GetXaxis()->SetTitle("Neutrino p_{y} (MeV)");
    h_corr2_PY_PZ->GetYaxis()->SetTitle("Neutrino p_{z} (MeV)");
    h_corr2_PY_PZ->SetTitle("");
    h_corr2_PY_PZ->Draw("COLZ");
    c6.SaveAs(folder_name+"corr2_PY_PZ.gif");
    c6.SaveAs(folder_name+"corr2_PY_PZ.pdf");

    TCanvas c7;
    c7.cd();
    h_corr1_M_PZ->GetXaxis()->SetTitle("Antineutrino m_{B} (MeV)");
    h_corr1_M_PZ->GetYaxis()->SetTitle("Antineutrino p_{z} (MeV)");
    h_corr1_M_PZ->SetTitle("");
    h_corr1_M_PZ->Draw("COLZ");
    c7.SaveAs(folder_name+"corr1_M_PZ.gif");
    c7.SaveAs(folder_name+"corr1_M_PZ.pdf");

    TCanvas c8;
    c8.cd();
    h_corr2_M_PZ->GetXaxis()->SetTitle("Neutrino m_{B} (MeV)");
    h_corr2_M_PZ->GetYaxis()->SetTitle("Neutrino p_{z} (MeV)");
    h_corr2_M_PZ->SetTitle("");
    h_corr2_M_PZ->Draw("COLZ");
    c8.SaveAs(folder_name+"corr2_M_PZ.gif");
    c8.SaveAs(folder_name+"corr2_M_PZ.pdf");

    TCanvas b;
    b.cd();
    h_visible_mass->SetLineColor(kRed);
    h_DTF_pass_mass->SetLineColor(kBlue);
    h_DTF_pass_mass->GetXaxis()->SetTitle("m_{B^{+}} (MeV)");
    h_DTF_pass_mass->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_DTF_pass_mass->GetXaxis()->GetXmax() - h_DTF_pass_mass->GetXaxis()->GetXmin())/h_DTF_pass_mass->GetNbinsX()) );
    h_DTF_pass_mass->DrawNormalized();
    h_visible_mass->DrawNormalized("same");
    TLegend* l;
    l = new TLegend(0.7, 0.6, 0.85, 0.85);
    l->SetBorderSize(0);
    l->AddEntry(h_visible_mass, "Visible", "lp");
    l->AddEntry(h_DTF_pass_mass, "DTF (pass)", "lp");
    l->Draw("same");
    b.SaveAs(folder_name+"Bmass_visible_DTF.gif");
    b.SaveAs(folder_name+"Bmass_visible_DTF.pdf");

    TCanvas b1;
    b1.cd();
    h_DTF_fail_mass->SetLineColor(kRed);
    h_DTF_pass_mass->SetLineColor(kBlue);
    h_DTF_pass_mass->GetXaxis()->SetTitle("m_{B^{+}} (MeV)");
    h_DTF_pass_mass->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_DTF_pass_mass->GetXaxis()->GetXmax() - h_DTF_pass_mass->GetXaxis()->GetXmin())/h_DTF_pass_mass->GetNbinsX()) );
    h_DTF_pass_mass->GetXaxis()->SetRangeUser(4000, 8000);
    h_DTF_pass_mass->DrawNormalized();
    h_DTF_fail_mass->DrawNormalized("same");
    TLegend* l1;
    l1 = new TLegend(0.6, 0.6, 0.85, 0.85);
    l1->SetBorderSize(0);
    l1->AddEntry(h_DTF_pass_mass, Form("DTF (pass) : %.0lf events", h_DTF_pass_mass->GetEntries()), "lp");
    l1->AddEntry(h_DTF_fail_mass, Form("DTF (fail) : %.0lf events", h_DTF_fail_mass->GetEntries()), "lp");
    l1->Draw("same");
    b1.SaveAs(folder_name+"Bmass_DTF_pass_fail.gif");
    b1.SaveAs(folder_name+"Bmass_DTF_pass_fail.pdf");

    TCanvas b2;
    b2.cd();
    h_mass->SetLineColor(kRed);
    h_mass_s->SetLineColor(kBlue);
    h_mass->SetFillColorAlpha(kRed,0.25);
    h_mass_s->SetFillColorAlpha(kBlue,0.25);
    h_mass->GetXaxis()->SetTitle("m_{B^{+}} (MeV)");
    h_mass->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass_s->DrawNormalized();
    h_mass->DrawNormalized("same");
    TLegend* l2;
    l2 = new TLegend(0.6, 0.6, 0.85, 0.85);
    l2->SetBorderSize(0);
    l2->AddEntry(h_mass, "B^{+} -> K^{+}  #tau #tau", "lp");
    l2->AddEntry(h_mass_s, "B^{0} -> K^{*0} #tau #tau", "lp");
    l2->Draw("same");
    b2.SaveAs(folder_name+"Bmass_KpKstar.gif");
    b2.SaveAs(folder_name+"Bmass_KpKstar.pdf");

    TCanvas b3;
    b3.cd();
    h_mass_best_cut->SetLineColor(kRed);
    h_mass_s->SetLineColor(kBlue);
    h_mass_best_cut->SetFillColorAlpha(kRed,0.25);
    h_mass_s->SetFillColorAlpha(kBlue,0.25);
    h_mass_best_cut->GetXaxis()->SetTitle("m_{B^{+}} (MeV)");
    h_mass_best_cut->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass_s->DrawNormalized();
    h_mass_best_cut->DrawNormalized("same");
    TLegend* l3;
    l3 = new TLegend(0.6, 0.5, 0.85, 0.85);
    l3->SetBorderSize(0);
    l3->AddEntry(h_mass, "B^{+} -> K^{+}  #tau #tau, P_{#nu}^{DTF} < 5000 MeV", "lp");
    l3->AddEntry(h_mass_s, "B^{0} -> K^{*0} #tau #tau", "lp");
    l3->Draw("same");
    b3.SaveAs(folder_name+"Bmass_KpKstar1.gif");
    b3.SaveAs(folder_name+"Bmass_KpKstar1.pdf");

    TCanvas b4;
    b4.cd();
    h_mass_2ndbest_cut->SetLineColor(kRed);
    h_mass_s->SetLineColor(kBlue);
    h_mass_2ndbest_cut->SetFillColorAlpha(kRed,0.25);
    h_mass_s->SetFillColorAlpha(kBlue,0.25);
    h_mass_2ndbest_cut->GetXaxis()->SetTitle("m_{B^{+}} (MeV)");
    h_mass_2ndbest_cut->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass->GetXaxis()->GetXmax() - h_mass->GetXaxis()->GetXmin())/h_mass->GetNbinsX()) );
    h_mass_s->DrawNormalized();
    h_mass_2ndbest_cut->DrawNormalized("same");
    TLegend* l4;
    l4 = new TLegend(0.6, 0.5, 0.85, 0.85);
    l4->SetBorderSize(0);
    l4->AddEntry(h_mass, "B^{+} -> K^{+}  #tau #tau, m_{#tau}^{DTF} #in [1750,1800] MeV", "lp");
    l4->AddEntry(h_mass_s, "B^{0} -> K^{*0} #tau #tau", "lp");
    l4->Draw("same");
    b4.SaveAs(folder_name+"Bmass_KpKstar2.gif");
    b4.SaveAs(folder_name+"Bmass_KpKstar2.pdf");

    gStyle->SetOptFit(1011);
    TCanvas c9;
    c9.cd();
    h_pull->Fit("gaus","","",-1.5,4);
    h_pull->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/#Delta m^{DTF}");
    h_pull->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_pull->GetXaxis()->GetXmax() - h_pull->GetXaxis()->GetXmin())/h_pull->GetNbinsX()) );
    h_pull->SetTitle("Pull distribution");
    h_pull->Draw();
    c9.SaveAs(folder_name+"pull.gif");
    c9.SaveAs(folder_name+"pull.pdf");

    gStyle->SetOptFit(0);
    TCanvas c10;
    c10.cd();
    h_taup_FD_nuPz_core->GetXaxis()->SetTitle("#tau^{+} FD (mm)");
    h_taup_FD_nuPz_core->GetYaxis()->SetTitle("Anti-neutrino p_{z} (MeV)");
    h_taup_FD_nuPz_core->SetTitle("Core");
    h_taup_FD_nuPz_core->Draw("COL");
    c10.SaveAs(folder_name+"taup_FD_nuPz_core.gif");
    c10.SaveAs(folder_name+"taup_FD_nuPz_core.pdf");

    TCanvas c11;
    c11.cd();
    h_taup_FD_nuPz_tail->GetXaxis()->SetTitle("#tau^{+} FD (mm)");
    h_taup_FD_nuPz_tail->GetYaxis()->SetTitle("Anti-neutrino p_{z} (MeV)");
    h_taup_FD_nuPz_tail->SetTitle("Tail");
    h_taup_FD_nuPz_tail->Draw("COL");
    c11.SaveAs(folder_name+"taup_FD_nuPz_tail.gif");
    c11.SaveAs(folder_name+"taup_FD_nuPz_tail.pdf");

    TCanvas c12;
    c12.cd();
    h_taum_FD_nuPz_core->GetXaxis()->SetTitle("#tau^{-} FD (mm)");
    h_taum_FD_nuPz_core->GetYaxis()->SetTitle("Neutrino p_{z} (MeV)");
    h_taum_FD_nuPz_core->SetTitle("Core");
    h_taum_FD_nuPz_core->Draw("COL");
    c12.SaveAs(folder_name+"taum_FD_nuPz_core.gif");
    c12.SaveAs(folder_name+"taum_FD_nuPz_core.pdf");

    TCanvas c13;
    c13.cd();
    h_taum_FD_nuPz_tail->GetXaxis()->SetTitle("#tau^{-} FD (mm)");
    h_taum_FD_nuPz_tail->GetYaxis()->SetTitle("Neutrino p_{z} (MeV)");
    h_taum_FD_nuPz_tail->SetTitle("Tail");
    h_taum_FD_nuPz_tail->Draw("COL");
    c12.SaveAs(folder_name+"taum_FD_nuPz_tail.gif");
    c12.SaveAs(folder_name+"taum_FD_nuPz_tail.pdf");

    TCanvas c14;
    c14.cd();
    h_m_merr_pass->GetXaxis()->SetTitle("#Delta m_{B} (MeV)");
    h_m_merr_pass->GetYaxis()->SetTitle("m_{B} (MeV)");
    h_m_merr_pass->SetTitle("Pass");
    h_m_merr_pass->Draw("COL");
    c14.SaveAs(folder_name+"m_merr_pass.gif");
    c14.SaveAs(folder_name+"m_merr_pass.pdf");

    TCanvas c16;
    c16.cd();
    h_m_merr_fail->GetXaxis()->SetTitle("#Delta m_{B} (MeV)");
    h_m_merr_fail->GetYaxis()->SetTitle("m_{B} (MeV)");
    h_m_merr_fail->SetTitle("Fail");
    h_m_merr_fail->Draw("COL");
    c16.SaveAs(folder_name+"m_merr_fail.gif");
    c16.SaveAs(folder_name+"m_merr_fail.pdf");

    TCanvas c17;
    c17.cd();
    h_pull_pnu1z->GetYaxis()->SetTitle("Anti-neutrino p_{z}^{DTF} - p_{z}^{TRUE} (MeV)");
    h_pull_pnu1z->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/#Delta m^{DTF}");
    h_pull_pnu1z->Draw("COL");
    c17.SaveAs(folder_name+"pull_pnu1z.gif");
    c17.SaveAs(folder_name+"pull_pnu1z.pdf");

    TCanvas c18;
    c18.cd();
    h_pull_pnu2z->GetYaxis()->SetTitle("Neutrino p_{z}^{DTF} - p_{z}^{TRUE} (MeV)");
    h_pull_pnu2z->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/#Delta m^{DTF}");
    h_pull_pnu2z->Draw("COL");
    c18.SaveAs(folder_name+"pull_pnu2z.gif");
    c18.SaveAs(folder_name+"pull_pnu2z.pdf");

    TCanvas c19;
    c19.cd();
    h_pull_taup->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/#Delta m^{DTF}");
    h_pull_taup->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_pull_taup->GetXaxis()->GetXmax() - h_pull_taup->GetXaxis()->GetXmin())/h_pull_taup->GetNbinsX()) );
    h_pull_taup->SetTitle("Pull distribution (#tau^{+})");
    h_pull_taup->Draw();
    c19.SaveAs(folder_name+"pull_taup.gif");
    c19.SaveAs(folder_name+"pull_taup.pdf");

    TCanvas c20;
    c20.cd();
    h_pull_taum->GetXaxis()->SetTitle("(m^{DTF} - m^{PDG})/#Delta m^{DTF}");
    h_pull_taum->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_pull_taum->GetXaxis()->GetXmax() - h_pull_taum->GetXaxis()->GetXmin())/h_pull_taum->GetNbinsX()) );
    h_pull_taum->SetTitle("Pull distribution (#tau^{-})");
    h_pull_taum->Draw();
    c20.SaveAs(folder_name+"pull_taum.gif");
    c20.SaveAs(folder_name+"pull_taum.pdf");

    TCanvas c21;
    c21.cd();
    h_DTF_taup_taum_M->GetYaxis()->SetTitle("DTF #tau^{-} mass (MeV)");
    h_DTF_taup_taum_M->GetXaxis()->SetTitle("DTF #tau^{+} mass (MeV)");
    h_DTF_taup_taum_M->Draw("COL");
    c21.SaveAs(folder_name+"DTF_taup_taum_M.gif");
    c21.SaveAs(folder_name+"DTF_taup_taum_M.pdf");

    TCanvas c22;
    c22.cd();
    h_DTF_taup_Bmass_M->GetXaxis()->SetTitle("DTF B^{+} mass (MeV)");
    h_DTF_taup_Bmass_M->GetYaxis()->SetTitle("DTF #tau^{+} mass (MeV)");
    h_DTF_taup_Bmass_M->Draw("COL");
    c22.SaveAs(folder_name+"DTF_taup_Bmass_M.gif");
    c22.SaveAs(folder_name+"DTF_taup_Bmass_M.pdf");

    TCanvas c23;
    c23.cd();
    h_DTF_taum_Bmass_M->GetXaxis()->SetTitle("DTF B^{+} mass (MeV)");
    h_DTF_taum_Bmass_M->GetYaxis()->SetTitle("DTF #tau^{-} mass (MeV)");
    h_DTF_taum_Bmass_M->Draw("COL");
    c23.SaveAs(folder_name+"DTF_taum_Bmass_M.gif");
    c23.SaveAs(folder_name+"DTF_taum_Bmass_M.pdf");

}

RooRealVar* fit(RooDataSet* data, RooRealVar mass, TString cut, TString name, TString folder_name){

    RooRealVar xp("xp", "xp", 5000., 6000.);
    RooRealVar sigp("sigp", "sigp", 0., 500.);
    RooRealVar xi("xi", "xi", -100. ,100.);
    RooRealVar ro1("ro1", "ro1", -5., 5.);
    RooRealVar ro2("ro2", "ro2", -5., 5.);
    RooBukinPdf* pdf = new RooBukinPdf("pdf", "pdf", mass, xp, sigp, xi, ro1, ro2); 

    RooFitResult* fit = pdf->fitTo(*data, Minos(true), Save());
    fit->Print();
    RooRealVar* sigma = (RooRealVar*)fit->floatParsFinal().find("sigp");

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
    RooPlot* massframe = mass.frame(Title(cut));
    data->plotOn(massframe, RooFit::Name("Data"));
    pdf->plotOn(massframe, RooFit::Name("Fit"), LineColor(kRed), LineStyle(1), LineWidth(2));
    pdf->paramOn(massframe,Layout(0.60,0.99,0.8));
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
    RooPlot* pull_plot = mass.frame(Title(""));

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

    c.SaveAs(folder_name+name+".gif");
    c.SaveAs(folder_name+name+".pdf");

    return sigma;
}