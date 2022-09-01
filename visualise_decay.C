#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;

ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint,  ROOT::Math::XYZPoint thePoint, bool invFlag);
void plot(ROOT::Math::XYZPoint PV_t, ROOT::Math::XYZPoint KV_t, ROOT::Math::XYZPoint BV_t, ROOT::Math::XYZPoint BV_true_t, ROOT::Math::XYZPoint DV1_t, ROOT::Math::XYZPoint DV2_t, ROOT::Math::XYZPoint DV1err, ROOT::Math::XYZPoint DV2err, ROOT::Math::XYZVector Xerr1_t, ROOT::Math::XYZVector Xerr2_t, ROOT::Math::XYZVector Yerr1_t, ROOT::Math::XYZVector Yerr2_t, ROOT::Math::XYZVector Pb_t, ROOT::Math::XYZVector Pnu1_t, ROOT::Math::XYZVector P3pi11_t, ROOT::Math::XYZVector P3pi12_t, ROOT::Math::XYZVector P3pi13_t, ROOT::Math::XYZVector Pnu2_t, ROOT::Math::XYZVector P3pi21_t, ROOT::Math::XYZVector P3pi22_t, ROOT::Math::XYZVector P3pi23_t, TString name, TFile* fout, int i, float DTF_chi2);
void plot_comparison(TH1D* histo_core, TH1D* histo_tail, TString name, TString unit, int n_bins);

bool DTF_neutrino = false;
bool all_true = false;

// DTF_neutrino = false, all_true = false -> nuTRUE
// DTF_neutrino = true, all_true = false -> nuDTF
// all_true = true, DTF_neutrino = false -> allTRUE

void visualise_decay(){
    if(all_true){DTF_neutrino = false;}
    if(DTF_neutrino){all_true = false;}

    TFile* fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");
    TTree* t = (TTree*)fin->Get("ntuple/DecayTree");

    TFile* fout;
    if(DTF_neutrino){fout = new TFile("Mass_resolution_visualisation/MC_2016_kinematics_nuDTF.root","recreate");}
    else if(all_true){fout = new TFile("Mass_resolution_visualisation/MC_2016_kinematics_allTRUE.root","recreate");}
    else{fout = new TFile("Mass_resolution_visualisation/MC_2016_kinematics_nuTRUE.root","recreate");}
    
    double n_bins = 100;
    TH1D* h_pass_core = new TH1D("h_pass_core", "h_pass_core", n_bins, 0., 15);
    TH1D* h_pass_tail = new TH1D("h_pass_tail", "h_pass_tail", n_bins, 0., 15);
    TH1D* h_fail = new TH1D("h_fail", "h_fail", n_bins, 0., 15);

    TH1D* h_PV_core = new TH1D("h_PV_IP_core", "h_PV_IP_core", n_bins, 0, 6);
    TH1D* h_DV1_core = new TH1D("h_DV1_IP_core", "h_DV1_IP_core", n_bins, 0, 2);
    TH1D* h_DV2_core = new TH1D("h_DV2_IP_core", "h_DV2_IP_core", n_bins, 0, 2);
    TH1D* h_PV_tail = new TH1D("h_PV_IP_tail", "h_PV_IP_tail", n_bins, 0, 6);
    TH1D* h_DV1_tail = new TH1D("h_DV1_IP_tail", "h_DV1_IP_tail", n_bins, 0, 2);
    TH1D* h_DV2_tail = new TH1D("h_DV2_IP_tail", "h_DV2_IP_tail", n_bins, 0, 2);
    TH1D* h_PV_fail = new TH1D("h_PV_IP_fail", "h_PV_IP_fail", n_bins, 0, 6);
    TH1D* h_DV1_fail = new TH1D("h_DV1_IP_fail", "h_DV1_IP_fail", n_bins, 0, 2);
    TH1D* h_DV2_fail = new TH1D("h_DV2_IP_fail", "h_DV2_IP_fail", n_bins, 0, 2);

    // Histos for comparins core with tail
    TH1D* h_taup_FDCHI2_BV_core = new TH1D("h_taup_FDCHI2_BV_core", "h_taup_FDCHI2_core", n_bins, 0., 50.);
    TH1D* h_taup_FDCHI2_BV_tail = new TH1D("h_taup_FDCHI2_BV_tail", "h_taup_FDCHI2_tail", n_bins, 0., 50.);
    TH1D* h_taum_FDCHI2_BV_core = new TH1D("h_taum_FDCHI2_BV_core", "h_taum_FDCHI2_core", n_bins, 0., 50.);
    TH1D* h_taum_FDCHI2_BV_tail = new TH1D("h_taum_FDCHI2_BV_tail", "h_taum_FDCHI2_tail", n_bins, 0., 50.);  
    TH1D* h_taup_FD_BV_core = new TH1D("h_taup_FD_BV_core", "h_taup_FD_core", n_bins, 0., 50.);
    TH1D* h_taup_FD_BV_tail = new TH1D("h_taup_FD_BV_tail", "h_taup_FD_tail", n_bins, 0., 50.);
    TH1D* h_taum_FD_BV_core = new TH1D("h_taum_FD_BV_core", "h_taum_FD_core", n_bins, 0., 50.);
    TH1D* h_taum_FD_BV_tail = new TH1D("h_taum_FD_BV_tail", "h_taum_FD_tail", n_bins, 0., 50.); 
    TH1D* h_taup_sep_core = new TH1D("h_taup_sep_core", "h_taup_sep_core", n_bins, 0., 30.);
    TH1D* h_taup_sep_tail = new TH1D("h_taup_sep_tail", "h_taup_sep_tail", n_bins, 0., 30.);
    TH1D* h_taup_DIRA_BV_core = new TH1D("h_taup_DIRA_BV_core", "h_taup_DIRA_BV_core", n_bins, -1.5, 1.5);
    TH1D* h_taup_DIRA_BV_tail = new TH1D("h_taup_DIRA_BV_tail", "h_taup_DIRA_BV_tail", n_bins, -1.5, 1.5);
    TH1D* h_nutau_P_core = new TH1D("h_nutau_P_core", "h_nutau_P_core", n_bins, 0., 150000);
    TH1D* h_nutau_P_tail = new TH1D("h_nutau_P_tail", "h_nutau_P_tail", n_bins, 0., 150000);
    TH1D* h_nutau0_P_core = new TH1D("h_nutau0_P_core", "h_nutau0_P_core", n_bins, 0., 150000);
    TH1D* h_nutau0_P_tail = new TH1D("h_nutau0_P_tail", "h_nutau0_P_tail", n_bins, 0., 150000);

    Float_t Bmass;
    Double_t P3pi11x, P3pi12x, P3pi13x, P3pi21x, P3pi22x, P3pi23x, P3pi11y, P3pi12y, P3pi13y, P3pi21y, P3pi22y, P3pi23y, P3pi11z, P3pi12z, P3pi13z, P3pi21z, P3pi22z, P3pi23z, P3pi2x, P3pi2y, P3pi2z, Pkx, Pky, Pkz;
    Double_t KVx, KVy, KVz, PVx, PVy, PVz, BVx, BVy, BVz;
    Double_t DV1x, DV1y, DV1z, DV2x, DV2y, DV2z, DV2xerr, DV1xerr, DV1yerr, DV1zerr, DV2yerr, DV2zerr;
    Double_t taup_FD_BV_chi2, taup_FD_BV, taum_FD_BV_chi2, taum_FD_BV, taup_DIRA_BV, taup_IP_PV_chi2, taup_IP_PV;
    Float_t status, Pnu1x, Pnu1y, Pnu1z, Pnu2x, Pnu2y, Pnu2z;
    Int_t Kp_TRUEID, taup_pip0_TRUEID, taup_pim0_TRUEID, taup_pip1_TRUEID, taum_pim0_TRUEID, taum_pip0_TRUEID, taum_pim1_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;
    Double_t BVx_true, BVy_true, BVz_true;
    Double_t taup_TRUEPX, taum_TRUEPX, taup_TRUEPY, taum_TRUEPY, taup_TRUEPZ, taum_TRUEPZ;
    Float_t DTF_chi2, DTF_ndf, DTF_taup_M, DTF_taum_M;

    t->SetBranchAddress("taup_FDCHI2_ORIVX",&taup_FD_BV_chi2);
    t->SetBranchAddress("taup_FD_ORIVX",&taup_FD_BV);
    t->SetBranchAddress("taum_FDCHI2_ORIVX",&taum_FD_BV_chi2);
    t->SetBranchAddress("taum_FD_ORIVX",&taum_FD_BV);
    t->SetBranchAddress("taup_DIRA_ORIVX",&taup_DIRA_BV);
    t->SetBranchAddress("taup_IPCHI2_OWNPV",&taup_IP_PV_chi2);
    t->SetBranchAddress("taup_IP_OWNPV",&taup_IP_PV_chi2);

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

    t->SetBranchAddress("Bp_ENDVERTEX_X",&BVx);
    t->SetBranchAddress("Bp_ENDVERTEX_Y",&BVy);
    t->SetBranchAddress("Bp_ENDVERTEX_Z",&BVz);
    t->SetBranchAddress("Bp_TRUEENDVERTEX_X",&BVx_true);
    t->SetBranchAddress("Bp_TRUEENDVERTEX_Y",&BVy_true);
    t->SetBranchAddress("Bp_TRUEENDVERTEX_Z",&BVz_true);

    t->SetBranchAddress("refPoint_X",&KVx);
    t->SetBranchAddress("refPoint_Y",&KVy);
    t->SetBranchAddress("refPoint_Z",&KVz);

    t->SetBranchAddress("Bp_ConsBp_tauminus_M",&DTF_taup_M);
    t->SetBranchAddress("Bp_ConsBp_tauminus_0_M",&DTF_taum_M);

    if(DTF_neutrino){
      t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PX",&Pnu1x);
      t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PY",&Pnu1y);
      t->SetBranchAddress("Bp_ConsBp_tauminus_nu_tau_PZ",&Pnu1z);
      t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PX",&Pnu2x);
      t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PY",&Pnu2y);
      t->SetBranchAddress("Bp_ConsBp_tauminus_0_nu_tau_PZ",&Pnu2z);
    }
    else{
      t->SetBranchAddress("taup_TRUEP_X", &taup_TRUEPX);
      t->SetBranchAddress("taum_TRUEP_X", &taum_TRUEPX);
      t->SetBranchAddress("taup_pip0_TRUEP_X",&P3pi11x);
      t->SetBranchAddress("taup_pim0_TRUEP_X",&P3pi12x);
      t->SetBranchAddress("taup_pip1_TRUEP_X",&P3pi13x);
      t->SetBranchAddress("taum_pim0_TRUEP_X",&P3pi21x);
      t->SetBranchAddress("taum_pip0_TRUEP_X",&P3pi22x);
      t->SetBranchAddress("taum_pim1_TRUEP_X",&P3pi23x);


      t->SetBranchAddress("taup_TRUEP_Y", &taup_TRUEPY);
      t->SetBranchAddress("taum_TRUEP_Y", &taum_TRUEPY);
      t->SetBranchAddress("taup_pip0_TRUEP_Y",&P3pi11y);
      t->SetBranchAddress("taup_pim0_TRUEP_Y",&P3pi12y);
      t->SetBranchAddress("taup_pip1_TRUEP_Y",&P3pi13y);
      t->SetBranchAddress("taum_pim0_TRUEP_Y",&P3pi21y);
      t->SetBranchAddress("taum_pip0_TRUEP_Y",&P3pi22y);
      t->SetBranchAddress("taum_pim1_TRUEP_Y",&P3pi23y);

      t->SetBranchAddress("taup_TRUEP_Z", &taup_TRUEPZ);
      t->SetBranchAddress("taum_TRUEP_Z", &taum_TRUEPZ);
      t->SetBranchAddress("taup_pip0_TRUEP_Z",&P3pi11z);
      t->SetBranchAddress("taup_pim0_TRUEP_Z",&P3pi12z);
      t->SetBranchAddress("taup_pip1_TRUEP_Z",&P3pi13z);
      t->SetBranchAddress("taum_pim0_TRUEP_Z",&P3pi21z);
      t->SetBranchAddress("taum_pip0_TRUEP_Z",&P3pi22z);
      t->SetBranchAddress("taum_pim1_TRUEP_Z",&P3pi23z);  
    }

    t->SetBranchAddress("Bp_ConsBp_M",&Bmass);
    t->SetBranchAddress("Bp_ConsBp_chi2",&DTF_chi2);
    t->SetBranchAddress("Bp_ConsBp_nDOF",&DTF_ndf);

    if(all_true){
      t->SetBranchAddress("Bp_TRUEORIGINVERTEX_X",&PVx);
      t->SetBranchAddress("Bp_TRUEORIGINVERTEX_Y",&PVy);
      t->SetBranchAddress("Bp_TRUEORIGINVERTEX_Z",&PVz);

      t->SetBranchAddress("taup_TRUEENDVERTEX_X",&DV1x);
      t->SetBranchAddress("taup_TRUEENDVERTEX_Y",&DV1y);
      t->SetBranchAddress("taup_TRUEENDVERTEX_Z",&DV1z);

      t->SetBranchAddress("taum_TRUEENDVERTEX_X",&DV2x);
      t->SetBranchAddress("taum_TRUEENDVERTEX_Y",&DV2y);
      t->SetBranchAddress("taum_TRUEENDVERTEX_Z",&DV2z);

      DV1xerr = 0.;
      DV1yerr = 0.;
      DV1zerr = 0.;

      DV2xerr = 0.;
      DV2yerr = 0.;
      DV2zerr = 0.;

      t->SetBranchAddress("Kp_TRUEP_X",&Pkx);
      t->SetBranchAddress("Kp_TRUEP_Y",&Pky);
      t->SetBranchAddress("Kp_TRUEP_Z",&Pkz);

    }
    else{
      t->SetBranchAddress("Bp_OWNPV_X",&PVx);
      t->SetBranchAddress("Bp_OWNPV_Y",&PVy);
      t->SetBranchAddress("Bp_OWNPV_Z",&PVz);

      t->SetBranchAddress("taup_pip0_PX",&P3pi11x);
      t->SetBranchAddress("taup_pim0_PX",&P3pi12x);
      t->SetBranchAddress("taup_pip1_PX",&P3pi13x);

      t->SetBranchAddress("taup_pip0_PY",&P3pi11y);
      t->SetBranchAddress("taup_pim0_PY",&P3pi12y);
      t->SetBranchAddress("taup_pip1_PY",&P3pi13y);

      t->SetBranchAddress("taup_pip0_PZ",&P3pi11z);
      t->SetBranchAddress("taup_pim0_PZ",&P3pi12z);
      t->SetBranchAddress("taup_pip1_PZ",&P3pi13z);

      t->SetBranchAddress("taum_pim0_PX",&P3pi21x);
      t->SetBranchAddress("taum_pip0_PX",&P3pi22x);
      t->SetBranchAddress("taum_pim1_PX",&P3pi23x);

      t->SetBranchAddress("taum_pim0_PY",&P3pi21y);
      t->SetBranchAddress("taum_pip0_PY",&P3pi22y);
      t->SetBranchAddress("taum_pim1_PY",&P3pi23y);

      t->SetBranchAddress("taum_pim0_PZ",&P3pi21z);
      t->SetBranchAddress("taum_pip0_PZ",&P3pi22z);
      t->SetBranchAddress("taum_pim1_PZ",&P3pi23z);

      t->SetBranchAddress("taup_ENDVERTEX_X",&DV1x);
      t->SetBranchAddress("taup_ENDVERTEX_Y",&DV1y);
      t->SetBranchAddress("taup_ENDVERTEX_Z",&DV1z);

      t->SetBranchAddress("taum_ENDVERTEX_X",&DV2x);
      t->SetBranchAddress("taum_ENDVERTEX_Y",&DV2y);
      t->SetBranchAddress("taum_ENDVERTEX_Z",&DV2z);

      t->SetBranchAddress("taup_ENDVERTEX_XERR",&DV1xerr);
      t->SetBranchAddress("taup_ENDVERTEX_YERR",&DV1yerr);
      t->SetBranchAddress("taup_ENDVERTEX_ZERR",&DV1zerr);

      t->SetBranchAddress("taum_ENDVERTEX_XERR",&DV2xerr);
      t->SetBranchAddress("taum_ENDVERTEX_YERR",&DV2yerr);
      t->SetBranchAddress("taum_ENDVERTEX_ZERR",&DV2zerr);

      t->SetBranchAddress("Kp_PX",&Pkx);
      t->SetBranchAddress("Kp_PY",&Pky);
      t->SetBranchAddress("Kp_PZ",&Pkz);

    }

    double mpeak = 5279; // from TTree 

    for(int i = 0; i < t->GetEntries(); i++){
      t->GetEntry(i);

      if(!(DTF_neutrino)){
        Pnu1x = taup_TRUEPX - (P3pi11x + P3pi12x + P3pi13x);
        Pnu2x = taum_TRUEPX - (P3pi21x + P3pi22x + P3pi23x);
        Pnu1y = taup_TRUEPY - (P3pi11y + P3pi12y + P3pi13y);
        Pnu2y = taum_TRUEPY - (P3pi21y + P3pi22y + P3pi23y);    
        Pnu1z = taup_TRUEPZ - (P3pi11z + P3pi12z + P3pi13z);
        Pnu2z = taum_TRUEPZ - (P3pi21z + P3pi22z + P3pi23z);
      }

      // Initialise points and vectors
      ROOT::Math::XYZPoint KV(KVx,KVy,KVz);
      ROOT::Math::XYZPoint BV(BVx, BVy, BVz);
      ROOT::Math::XYZPoint BV_true(BVx_true, BVy_true, BVz_true);
      ROOT::Math::XYZPoint PV(PVx,PVy,PVz);
      ROOT::Math::XYZPoint DV1(DV1x,DV1y,DV1z);
      ROOT::Math::XYZPoint DV2(DV2x,DV2y,DV2z);
      ROOT::Math::XYZPoint DV1err(DV1xerr,DV1yerr,DV1zerr);
      ROOT::Math::XYZPoint DV2err(DV2xerr,DV2yerr,DV2zerr);
      ROOT::Math::XYZVector Pnu1(Pnu1x,Pnu1y,Pnu1z);
      ROOT::Math::XYZVector Pnu2(Pnu2x,Pnu2y,Pnu2z);
      ROOT::Math::XYZVector Pk(Pkx,Pky,Pkz);
      ROOT::Math::XYZVector P3pi11(P3pi11x,P3pi11y,P3pi11z);
      ROOT::Math::XYZVector P3pi12(P3pi12x,P3pi12y,P3pi12z);
      ROOT::Math::XYZVector P3pi13(P3pi13x,P3pi13y,P3pi13z);
      ROOT::Math::XYZVector P3pi21(P3pi21x,P3pi21y,P3pi21z);
      ROOT::Math::XYZVector P3pi22(P3pi22x,P3pi22y,P3pi22z);
      ROOT::Math::XYZVector P3pi23(P3pi23x,P3pi23y,P3pi23z);
      ROOT::Math::XYZVector Xerr1(2*DV1err.x(), 0., 0.);
      ROOT::Math::XYZVector Xerr2(2*DV2err.x(), 0., 0.);
      ROOT::Math::XYZVector Yerr1(0., 2*DV1err.y(), 0.);
      ROOT::Math::XYZVector Yerr2(0., 2*DV2err.y(), 0.);

      // Transform x and y coordinates to plane perpendicular to K+ trajectory (z-axis points along K+ trajectory)
      ROOT::Math::XYZPoint KV_t = makeTransformation_point(Pk, KV, KV, false);
      ROOT::Math::XYZPoint BV_t = makeTransformation_point(Pk, KV, BV, false);
      ROOT::Math::XYZPoint BV_true_t = makeTransformation_point(Pk, KV, BV_true, false);
      ROOT::Math::XYZPoint PV_t = makeTransformation_point(Pk, KV, PV, false);
      ROOT::Math::XYZPoint DV1_t = makeTransformation_point(Pk, KV, DV1, false);
      ROOT::Math::XYZPoint DV2_t = makeTransformation_point(Pk, KV, DV2, false);
      ROOT::Math::XYZVector Pnu1_t = makeTransformation_vec(Pk, KV, Pnu1, false);
      ROOT::Math::XYZVector Pnu2_t = makeTransformation_vec(Pk, KV, Pnu2, false);
      ROOT::Math::XYZVector Pk_t = makeTransformation_vec(Pk, KV, Pk, false);
      ROOT::Math::XYZVector P3pi11_t = makeTransformation_vec(Pk, KV, P3pi11, false);
      ROOT::Math::XYZVector P3pi12_t = makeTransformation_vec(Pk, KV, P3pi12, false);
      ROOT::Math::XYZVector P3pi13_t = makeTransformation_vec(Pk, KV, P3pi13, false);
      ROOT::Math::XYZVector P3pi21_t = makeTransformation_vec(Pk, KV, P3pi21, false);
      ROOT::Math::XYZVector P3pi22_t = makeTransformation_vec(Pk, KV, P3pi22, false);
      ROOT::Math::XYZVector P3pi23_t = makeTransformation_vec(Pk, KV, P3pi23, false);
      ROOT::Math::XYZVector Xerr1_t = makeTransformation_vec(Pk, KV, Xerr1, false);
      ROOT::Math::XYZVector Xerr2_t = makeTransformation_vec(Pk, KV, Xerr2, false);
      ROOT::Math::XYZVector Yerr1_t = makeTransformation_vec(Pk, KV, Yerr1, false);
      ROOT::Math::XYZVector Yerr2_t = makeTransformation_vec(Pk, KV, Yerr2, false);
      ROOT::Math::XYZVector Pb_t = Pk_t + P3pi11_t + P3pi12_t + P3pi13_t + P3pi21_t + P3pi22_t + P3pi23_t + Pnu1_t + Pnu2_t; 

      double a = sqrt((PV-DV1).Mag2());
      double b = sqrt((PV-DV2).Mag2());
      double c = sqrt((DV2-DV1).Mag2());
      double s = (a+b+c)/2.;
      double A = sqrt(s*(s-a)*(s-b)*(s-c));

      ROOT::Math::XYZVector KV1(KV.x() + Pk.x(), KV.y() + Pk.y(), KV.z() + Pk.z()); 
      ROOT::Math::XYZVector A1(DV1.x() - KV.x(), DV1.y() - KV.y(), DV1.z() - KV.z());
      ROOT::Math::XYZVector A2(DV2.x() - KV.x(), DV2.y() - KV.y(), DV2.z() - KV.z());
      ROOT::Math::XYZVector A3(PV.x() - KV.x(), PV.y()-KV.y(), PV.z()-KV.z());

      double IP1 = (2*sqrt( (A1.Cross(Pk)).Mag2() ))/sqrt(Pk.Mag2());
      double IP2 = (2*sqrt( (A2.Cross(Pk)).Mag2() ))/sqrt(Pk.Mag2());
      double IP3 = (2*sqrt( (A3.Cross(Pk)).Mag2() ))/sqrt(Pk.Mag2());

      double tau_distance;
      float DTF_norm_chi2 = DTF_chi2/DTF_ndf;

      if( (abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)){ // truth match

        tau_distance = sqrt(pow(DV1.x() - DV2.x(),2) + pow(DV1.y() - DV2.y(),2) + pow(DV1.z() - DV2.z(),2));

        if(status == 0){ // pass DTF

          if((DTF_taup_M != 1776.86)){plot(PV_t, KV_t, BV_t, BV_true_t, DV1_t, DV2_t, DV1err, DV2err, Xerr1_t, Xerr2_t, Yerr1_t, Yerr2_t, Pb_t, Pnu1_t, P3pi11_t, P3pi12_t, P3pi13_t, Pnu2_t, P3pi21_t, P3pi22_t, P3pi23_t, "pass_taup_no_PDG", fout, i, DTF_norm_chi2);}
          if((DTF_taum_M != 1776.86)){plot(PV_t, KV_t, BV_t, BV_true_t, DV1_t, DV2_t, DV1err, DV2err, Xerr1_t, Xerr2_t, Yerr1_t, Yerr2_t, Pb_t, Pnu1_t, P3pi11_t, P3pi12_t, P3pi13_t, Pnu2_t, P3pi21_t, P3pi22_t, P3pi23_t, "pass_taum_no_PDG", fout, i, DTF_norm_chi2);}

          if( abs(Bmass - mpeak) < 300 ){ // core
            plot(PV_t, KV_t, BV_t, BV_true_t, DV1_t, DV2_t, DV1err, DV2err, Xerr1_t, Xerr2_t, Yerr1_t, Yerr2_t, Pb_t, Pnu1_t, P3pi11_t, P3pi12_t, P3pi13_t, Pnu2_t, P3pi21_t, P3pi22_t, P3pi23_t, "pass_core", fout, i, DTF_norm_chi2);
            
            h_pass_core->Fill(A);
            h_DV1_core->Fill(IP1);
            h_DV2_core->Fill(IP2);
            h_PV_core->Fill(IP3);

            h_taup_FDCHI2_BV_core->Fill(taup_FD_BV_chi2);
            h_taup_FD_BV_core->Fill(taup_FD_BV);
            h_taum_FDCHI2_BV_core->Fill(taum_FD_BV_chi2);
            h_taum_FD_BV_core->Fill(taum_FD_BV);
            h_taup_sep_core->Fill(tau_distance);
            h_taup_DIRA_BV_core->Fill(taup_DIRA_BV);
            h_nutau_P_core->Fill(sqrt(pow(Pnu1x,2) + pow(Pnu1y,2) + pow(Pnu1z,2)));
            h_nutau0_P_core->Fill(sqrt(pow(Pnu2x,2) + pow(Pnu2y,2) + pow(Pnu2z,2)));

          }
          else if(Bmass > 6500.){ // tails
            plot(PV_t, KV_t, BV_t, BV_true_t, DV1_t, DV2_t, DV1err, DV2err, Xerr1_t, Xerr2_t, Yerr1_t, Yerr2_t, Pb_t, Pnu1_t, P3pi11_t, P3pi12_t, P3pi13_t, Pnu2_t, P3pi21_t, P3pi22_t, P3pi23_t, "pass_tail", fout, i, DTF_norm_chi2);
            
            h_pass_tail->Fill(A);
            h_DV1_tail->Fill(IP1);
            h_DV2_tail->Fill(IP2);
            h_PV_tail->Fill(IP3);

            h_taup_FDCHI2_BV_tail->Fill(taup_FD_BV_chi2);
            h_taup_FD_BV_tail->Fill(taup_FD_BV);
            h_taum_FDCHI2_BV_tail->Fill(taum_FD_BV_chi2);
            h_taum_FD_BV_tail->Fill(taum_FD_BV);
            h_taup_sep_tail->Fill(tau_distance);
            h_taup_DIRA_BV_tail->Fill(taup_DIRA_BV);
            h_nutau_P_tail->Fill(sqrt(pow(Pnu1x,2) + pow(Pnu1y,2) + pow(Pnu1z,2)));
            h_nutau0_P_tail->Fill(sqrt(pow(Pnu2x,2) + pow(Pnu2y,2) + pow(Pnu2z,2)));

          }
        }
        else{ // fail DTF
          plot(PV_t, KV_t, BV_t, BV_true_t, DV1_t, DV2_t, DV1err, DV2err, Xerr1_t, Xerr2_t, Yerr1_t, Yerr2_t, Pb_t, Pnu1_t, P3pi11_t, P3pi12_t, P3pi13_t, Pnu2_t, P3pi21_t, P3pi22_t, P3pi23_t, "fail", fout, i, DTF_norm_chi2);
          
          h_fail->Fill(A);
          h_DV1_fail->Fill(IP1);
          h_DV2_fail->Fill(IP2);
          h_PV_fail->Fill(IP3);
          
        }
    }
  }

  plot_comparison(h_taup_FDCHI2_BV_core, h_taup_FDCHI2_BV_tail, "taup FD to BV chi2", " ", n_bins);
  plot_comparison(h_taup_FD_BV_core, h_taup_FDCHI2_BV_tail, "taup FD to BV", "(mm)", n_bins);
  plot_comparison(h_taum_FDCHI2_BV_core, h_taum_FDCHI2_BV_tail, "taum FD to BV chi2", " ", n_bins);
  plot_comparison(h_taum_FD_BV_core, h_taum_FDCHI2_BV_tail, "taum FD to BV", "(mm)", n_bins);
  plot_comparison(h_taup_sep_core, h_taup_sep_tail, "distance between tau vertices", "(mm)", n_bins);
  plot_comparison(h_taup_DIRA_BV_core, h_taup_DIRA_BV_tail, "taup DIRA to BV", "(rad)", n_bins);
  plot_comparison(h_nutau_P_core, h_nutau_P_tail, "antinutau momentum", "(MeV)", n_bins);
  plot_comparison(h_nutau0_P_core, h_nutau0_P_tail, "nutau momentum", "(MeV)", n_bins);

  TCanvas c2;
  c2.cd();
  gStyle->SetOptStat(0);
  h_pass_core->SetLineColor(kBlack);
  h_pass_tail->SetLineColor(kBlue);
  h_fail->SetLineColor(kRed);
  h_pass_core->SetTitle("Triangle area (PV, #tau^{+} and #tau^{-} vertices");
  h_pass_core->GetXaxis()->SetTitle("DV1-DV2-PV traingle area (mm^{2})");
  h_pass_core->GetYaxis()->SetTitle( TString::Format("Candidates / (%g)",(h_pass_core->GetXaxis()->GetXmax() - h_pass_core->GetXaxis()->GetXmin())/n_bins) );

  // normalization
  h_pass_core->Scale(1/h_PV_core->Integral());
  h_pass_tail->Scale(1/h_PV_tail->Integral());
  h_fail->Scale(1/h_PV_fail->Integral());

  //y axis: maximum and minimum
  double y[] = {h_pass_core->GetMaximum(), h_pass_tail->GetMaximum(), h_fail->GetMaximum()};
  h_pass_core->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);
  h_pass_tail->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);
  h_fail->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);

  TLegend* leg = new TLegend(0.3, 0.7, 0.9, 0.9);
  leg->AddEntry(h_pass_core, TString::Format("Pass core ; mean = %f",h_pass_core->GetMean()));
  leg->AddEntry(h_pass_tail, TString::Format("Pass tail ; mean = %f",h_pass_tail->GetMean()));
  leg->AddEntry(h_fail, TString::Format("Fail ; mean = %f",h_fail->GetMean()));

  h_pass_core->Draw("HIST");
  h_pass_tail->Draw("HIST same");
  h_fail->Draw("HIST same");
  leg->Draw("same");
  c2.SaveAs("Mass_resolution_visualisation/triangle_area.gif");
  c2.SaveAs("Mass_resolution_visualisation/triangle_area.pdf");

  TCanvas c4;
  c4.cd();
  gStyle->SetOptStat(0);
  h_DV1_core->SetLineColor(kBlack);
  h_DV1_tail->SetLineColor(kBlue);
  h_DV1_fail->SetLineColor(kRed);
  h_DV1_core->SetTitle("Distance from DV1 to K^{+} traj.");
  h_DV1_core->GetXaxis()->SetTitle("IP (mm)");
  h_DV1_core->GetYaxis()->SetTitle( TString::Format("Candidates / (%g)",(h_DV1_core->GetXaxis()->GetXmax() - h_DV1_core->GetXaxis()->GetXmin())/n_bins) );
 
  //y axis: maximum and minimum
  double y1[] = {h_DV1_core->GetMaximum(), h_DV1_tail->GetMaximum(), h_DV1_fail->GetMaximum()};
  h_DV1_core->GetYaxis()->SetRangeUser(0., y1[TMath::LocMax(sizeof(y1)/sizeof(y1[0]),y1)]);
  h_DV1_tail->GetYaxis()->SetRangeUser(0., y1[TMath::LocMax(sizeof(y1)/sizeof(y1[0]),y1)]);
  h_DV1_fail->GetYaxis()->SetRangeUser(0., y1[TMath::LocMax(sizeof(y1)/sizeof(y1[0]),y1)]);

  TLegend* leg1 = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg1->AddEntry(h_DV1_core, TString::Format("Pass core ; mean = %f",h_DV1_core->GetMean()));
  leg1->AddEntry(h_DV1_tail, TString::Format("Pass tail ; mean = %f",h_DV1_tail->GetMean()));
  leg1->AddEntry(h_DV1_fail, TString::Format("Fail ; mean = %f",h_DV1_fail->GetMean()));

  // normalization
  h_DV1_core->Scale(1/h_DV1_core->Integral());
  h_DV1_tail->Scale(1/h_DV1_tail->Integral());
  h_DV1_fail->Scale(1/h_DV1_fail->Integral());

  h_DV1_core->Draw("HIST");
  h_DV1_tail->Draw("HIST same");
  h_DV1_fail->Draw("HIST same");
  leg1->Draw("same");
  c4.SaveAs("Mass_resolution_visualisation/DV1_IP.gif");
  c4.SaveAs("Mass_resolution_visualisation/DV1_IP.pdf");

  TCanvas c5;
  c5.cd();
  gStyle->SetOptStat(0);
  h_DV2_core->SetLineColor(kBlack);
  h_DV2_tail->SetLineColor(kBlue);
  h_DV2_fail->SetLineColor(kRed);
  h_DV2_core->SetTitle("Distance from DV2 to K^{+} traj.");
  h_DV2_core->GetXaxis()->SetTitle("IP (mm)");
  h_DV2_core->GetYaxis()->SetTitle( TString::Format("Candidates / (%g)",(h_DV2_core->GetXaxis()->GetXmax() - h_DV2_core->GetXaxis()->GetXmin())/n_bins) );
 
  //y axis: maximum and minimum
  double y2[] = {h_DV2_core->GetMaximum(), h_DV2_tail->GetMaximum(), h_DV2_fail->GetMaximum()};
  h_DV2_core->GetYaxis()->SetRangeUser(0., y2[TMath::LocMax(sizeof(y2)/sizeof(y2[0]),y2)]);
  h_DV2_tail->GetYaxis()->SetRangeUser(0., y2[TMath::LocMax(sizeof(y2)/sizeof(y2[0]),y2)]);
  h_DV2_fail->GetYaxis()->SetRangeUser(0., y2[TMath::LocMax(sizeof(y2)/sizeof(y2[0]),y2)]);

  TLegend* leg2 = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg2->AddEntry(h_DV2_core, TString::Format("Pass core ; mean = %f",h_DV2_core->GetMean()));
  leg2->AddEntry(h_DV2_tail, TString::Format("Pass tail ; mean = %f",h_DV2_tail->GetMean()));
  leg2->AddEntry(h_DV2_fail, TString::Format("Fail ; mean = %f",h_DV2_fail->GetMean()));

  // normalisation
  h_DV2_core->Scale(1/h_DV2_core->Integral());
  h_DV2_tail->Scale(1/h_DV2_tail->Integral());
  h_DV2_fail->Scale(1/h_DV2_fail->Integral());

  h_DV2_core->Draw("HIST");
  h_DV2_tail->Draw("HIST same");
  h_DV2_fail->Draw("HIST same");
  leg2->Draw("same");
  c5.SaveAs("Mass_resolution_visualisation/DV2_IP.gif");
  c5.SaveAs("Mass_resolution_visualisation/DV2_IP.pdf");

  TCanvas c6;
  c6.cd();
  gStyle->SetOptStat(0);
  h_PV_core->SetLineColor(kBlack);
  h_PV_tail->SetLineColor(kBlue);
  h_PV_fail->SetLineColor(kRed);
  h_PV_core->SetTitle("Distance from PV to K^{+} traj.");
  h_PV_core->GetXaxis()->SetTitle("IP (mm)");
  h_PV_core->GetYaxis()->SetTitle( TString::Format("Candidates / (%g)",(h_PV_core->GetXaxis()->GetXmax() - h_PV_core->GetXaxis()->GetXmin())/n_bins) );
 
  //y axis: maximum and minimum
  double y3[] = {h_PV_core->GetMaximum(), h_PV_tail->GetMaximum(), h_PV_fail->GetMaximum()};
  h_PV_core->GetYaxis()->SetRangeUser(0., y3[TMath::LocMax(sizeof(y3)/sizeof(y3[0]),y3)]);
  h_PV_tail->GetYaxis()->SetRangeUser(0., y3[TMath::LocMax(sizeof(y3)/sizeof(y3[0]),y3)]);
  h_PV_fail->GetYaxis()->SetRangeUser(0., y3[TMath::LocMax(sizeof(y3)/sizeof(y3[0]),y3)]);

  TLegend* leg3 = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg3->AddEntry(h_PV_core, TString::Format("Pass core ; mean = %f",h_PV_core->GetMean()));
  leg3->AddEntry(h_PV_tail, TString::Format("Pass tail ; mean = %f",h_PV_tail->GetMean()));
  leg3->AddEntry(h_PV_fail, TString::Format("Fail ; mean = %f",h_PV_fail->GetMean()));

  // normalization
  h_PV_core->Scale(1/h_PV_core->Integral());
  h_PV_tail->Scale(1/h_PV_tail->Integral());
  h_PV_fail->Scale(1/h_PV_fail->Integral());

  h_PV_core->Draw("HIST");
  h_PV_tail->Draw("HIST same");
  h_PV_fail->Draw("HIST same");
  leg3->Draw("same");
  c6.SaveAs("Mass_resolution_visualisation/PV_IP.gif");
  c6.SaveAs("Mass_resolution_visualisation/PV_IP.pdf");

  return;
}

void plot_comparison(TH1D* histo_core, TH1D* histo_tail, TString name, TString unit, int n_bins){

  TCanvas c;
  c.cd();
  gStyle->SetOptStat(0);
  histo_core->SetLineColor(kBlue);
  histo_tail->SetLineColor(kRed);
  histo_core->SetTitle(name);
  histo_core->GetXaxis()->SetTitle(name + " " + unit);
  histo_core->GetYaxis()->SetTitle( TString::Format("Candidates / (%g)",(histo_core->GetXaxis()->GetXmax() - histo_core->GetXaxis()->GetXmin())/n_bins) );

  //y axis: maximum and minimum
  double y[] = {histo_core->GetMaximum(), histo_tail->GetMaximum()};
  histo_core->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);
  histo_tail->GetYaxis()->SetRangeUser(0., y[TMath::LocMax(sizeof(y)/sizeof(y[0]),y)]);

  TLegend* leg = new TLegend(0.6, 0.8, 0.9, 0.9);
  leg->AddEntry(histo_core, TString::Format("Core ; mean = %f",histo_core->GetMean()));
  leg->AddEntry(histo_tail, TString::Format("Tail ; mean = %f",histo_tail->GetMean()));

  // normalization
  histo_core->Scale(1/histo_core->Integral());
  histo_tail->Scale(1/histo_tail->Integral());

  histo_core->Draw("HIST");
  histo_tail->Draw("HIST same");
  leg->Draw("same");
  c.SaveAs("Mass_resolution_visualisation/"+name+"_comparison.gif");
  c.SaveAs("Mass_resolution_visualisation/"+name+"_comparison.pdf");

}

void plot(ROOT::Math::XYZPoint PV_t, ROOT::Math::XYZPoint KV_t, ROOT::Math::XYZPoint BV_t, ROOT::Math::XYZPoint BV_true_t, ROOT::Math::XYZPoint DV1_t, ROOT::Math::XYZPoint DV2_t, ROOT::Math::XYZPoint DV1err, ROOT::Math::XYZPoint DV2err,  ROOT::Math::XYZVector Xerr1_t, ROOT::Math::XYZVector Xerr2_t, ROOT::Math::XYZVector Yerr1_t, ROOT::Math::XYZVector Yerr2_t, ROOT::Math::XYZVector Pb_t, ROOT::Math::XYZVector Pnu1_t, ROOT::Math::XYZVector P3pi11_t, ROOT::Math::XYZVector P3pi12_t, ROOT::Math::XYZVector P3pi13_t, ROOT::Math::XYZVector Pnu2_t, ROOT::Math::XYZVector P3pi21_t, ROOT::Math::XYZVector P3pi22_t, ROOT::Math::XYZVector P3pi23_t, TString name, TFile* fout, int i, float DTF_chi2){

        double r1x = 0.5*sqrt(pow(Xerr1_t.x(),2) + pow(Xerr1_t.y(),2) + pow(Xerr1_t.z(),2));
        double r2x = 0.5*sqrt(pow(Xerr2_t.x(),2) + pow(Xerr2_t.y(),2) + pow(Xerr2_t.z(),2));
        double r1y = 0.5*sqrt(pow(Yerr1_t.x(),2) + pow(Yerr1_t.y(),2) + pow(Yerr1_t.z(),2));
        double r2y = 0.5*sqrt(pow(Yerr2_t.x(),2) + pow(Yerr2_t.y(),2) + pow(Yerr2_t.z(),2));
        double theta1 = acos(Xerr1_t.x()/sqrt(pow(Xerr1_t.x(),2) + pow(Xerr1_t.y(),2) + pow(Xerr1_t.z(),2)));
        double theta2 = acos(Xerr2_t.x()/sqrt(pow(Xerr2_t.x(),2) + pow(Xerr2_t.y(),2) + pow(Xerr2_t.z(),2)));
        theta1 *= 180/M_PI;
        theta2 *= 180/M_PI;

        TEllipse* e1 = new TEllipse(DV1_t.x(), DV1_t.y(), r1x, r1y, 0., 360., theta1);
        TEllipse* e2 = new TEllipse(DV2_t.x(), DV2_t.y(), r2x, r2y, 0., 360., theta2);
        e1->SetLineColor(kBlue);
        e2->SetLineColor(kRed);
        e1->SetLineStyle(7);
        e2->SetLineStyle(7);
        e1->SetLineWidth(2);
        e2->SetLineWidth(2);
        e1->SetFillColorAlpha(kBlue,0.5);
        e2->SetFillColorAlpha(kRed,0.5);

        double r = 0.0001;
        double Pb_lineX[2] = {0., Pb_t.x()*r};
        double Pb_lineY[2] = {0., Pb_t.y()*r};
        double Pnu1_lineX[2] = {DV1_t.x(), DV1_t.x() + Pnu1_t.x()*r};
        double Pnu1_lineY[2] = {DV1_t.y(), DV1_t.y() + Pnu1_t.y()*r};
        double Pnu2_lineX[2] = {DV2_t.x(), DV2_t.x() + Pnu2_t.x()*r};
        double Pnu2_lineY[2] = {DV2_t.y(), DV2_t.x() + Pnu2_t.y()*r};
        double P3pi11_lineX[2] = {DV1_t.x(), DV1_t.x() + P3pi11_t.x()*r};
        double P3pi11_lineY[2] = {DV1_t.y(), DV1_t.y() + P3pi11_t.y()*r};
        double P3pi12_lineX[2] = {DV1_t.x(), DV1_t.x() + P3pi12_t.x()*r};
        double P3pi12_lineY[2] = {DV1_t.y(), DV1_t.y() + P3pi12_t.y()*r};
        double P3pi13_lineX[2] = {DV1_t.x(), DV1_t.x() + P3pi13_t.x()*r};
        double P3pi13_lineY[2] = {DV1_t.y(), DV1_t.y() + P3pi13_t.y()*r};
        double P3pi21_lineX[2] = {DV2_t.x(), DV2_t.x() + P3pi21_t.x()*r};
        double P3pi21_lineY[2] = {DV2_t.y(), DV2_t.y() + P3pi21_t.y()*r};
        double P3pi22_lineX[2] = {DV2_t.x(), DV2_t.x() + P3pi22_t.x()*r};
        double P3pi22_lineY[2] = {DV2_t.y(), DV2_t.y() + P3pi22_t.y()*r};
        double P3pi23_lineX[2] = {DV2_t.x(), DV2_t.x() + P3pi23_t.x()*r};
        double P3pi23_lineY[2] = {DV2_t.y(), DV2_t.y() + P3pi23_t.y()*r};
        double PV_pointX[1] = {PV_t.x()};
        double PV_pointY[1] = {PV_t.y()};

        TGraph* Pb_line = new TGraph(2, Pb_lineX, Pb_lineY);
        TGraph* Pnu1_line = new TGraph(2, Pnu1_lineX, Pnu1_lineY);
        TGraph* Pnu2_line = new TGraph(2, Pnu2_lineX, Pnu2_lineY);
        TGraph* P3pi11_line = new TGraph(2, P3pi11_lineX, P3pi11_lineY);
        TGraph* P3pi12_line = new TGraph(2, P3pi12_lineX, P3pi12_lineY);
        TGraph* P3pi13_line = new TGraph(2, P3pi13_lineX, P3pi13_lineY);
        TGraph* P3pi21_line = new TGraph(2, P3pi21_lineX, P3pi21_lineY);
        TGraph* P3pi22_line = new TGraph(2, P3pi22_lineX, P3pi22_lineY);
        TGraph* P3pi23_line = new TGraph(2, P3pi23_lineX, P3pi23_lineY);

        TGraph* DV1_point = new TGraph();
        TGraph* DV2_point = new TGraph();
        TGraph* PV_point = new TGraph();
        TGraph* origin = new TGraph();
        DV1_point->SetPoint(0, DV1_t.x(), DV1_t.y());
        DV2_point->SetPoint(0, DV2_t.x(), DV2_t.y());
        PV_point->SetPoint(0, PV_t.x(), PV_t.y());
        origin->SetPoint(0, 0., 0.);
        DV1_point->SetMarkerColor(kBlue);
        DV2_point->SetMarkerColor(kRed);
        PV_point->SetMarkerColor(kBlack);
        origin->SetMarkerColor(kGray);
        DV1_point->SetMarkerStyle(8);
        DV2_point->SetMarkerStyle(8);
        PV_point->SetMarkerStyle(8);
        origin->SetMarkerStyle(8);

        TCanvas c = new TCanvas();
        c.SetName(TString::Format(name+"_c%i",i));
        c.SetWindowSize(1250,500);
        c.Divide(2,1);

        // XY PLANE
        c.cd(1);

        double X[] = {0., PV_t.x(), DV1_t.x(), DV2_t.x(), DV1_t.x() + P3pi11_t.x()*r, DV1_t.x() + P3pi12_t.x()*r, DV1_t.x() + P3pi13_t.x()*r, DV1_t.x() + P3pi21_t.x()*r, DV1_t.x() + P3pi22_t.x()*r, DV1_t.x() + P3pi23_t.x()*r, DV1_t.x() + Pnu1_t.x()*r, DV2_t.x() + Pnu2_t.x()*r, Pb_t.x()*r};
        double Y[] = {0., PV_t.y(), DV1_t.y(), DV2_t.y(), DV1_t.y() + P3pi11_t.y()*r, DV1_t.y() + P3pi12_t.y()*r, DV1_t.y() + P3pi13_t.y()*r, DV1_t.y() + P3pi21_t.y()*r, DV1_t.y() + P3pi22_t.y()*r, DV1_t.y() + P3pi23_t.y()*r, DV1_t.y() + Pnu1_t.y()*r, DV2_t.y() + Pnu2_t.y()*r, Pb_t.y()*r};
        double Xmin = X[TMath::LocMin(sizeof(X)/sizeof(X[0]),X)] - 0.1;
        double Xmax = X[TMath::LocMax(sizeof(X)/sizeof(X[0]),X)] + 0.1;
        double Ymin = Y[TMath::LocMin(sizeof(Y)/sizeof(Y[0]),Y)] - 0.1;
        double Ymax = Y[TMath::LocMax(sizeof(Y)/sizeof(Y[0]),Y)] + 0.1;

        gPad->DrawFrame(Xmin,Ymin,Xmax,Ymax);
        gPad->SetGrid(1,1);
        gPad->Draw();

        Pb_line->SetLineColor(kGray);
        Pnu1_line->SetLineColor(kBlue);
        P3pi11_line->SetLineColor(kBlue);
        P3pi12_line->SetLineColor(kBlue);
        P3pi13_line->SetLineColor(kBlue);
        Pnu2_line->SetLineColor(kRed);
        P3pi21_line->SetLineColor(kRed);
        P3pi22_line->SetLineColor(kRed);
        P3pi23_line->SetLineColor(kRed);
        Pb_line->SetLineWidth(2);
        Pnu1_line->SetLineWidth(2);
        P3pi11_line->SetLineWidth(2);
        P3pi12_line->SetLineWidth(2);
        P3pi13_line->SetLineWidth(2);
        Pnu2_line->SetLineWidth(2);
        P3pi21_line->SetLineWidth(2);
        P3pi22_line->SetLineWidth(2);
        P3pi23_line->SetLineWidth(2);
        Pnu1_line->SetLineStyle(2);
        Pnu2_line->SetLineStyle(2);

        TMultiGraph* mg = new TMultiGraph();
        mg->Add(Pb_line,"L");
        mg->Add(Pnu1_line,"L");
        mg->Add(P3pi11_line,"L");
        mg->Add(P3pi12_line,"L");
        mg->Add(P3pi13_line,"L");
        mg->Add(Pnu2_line,"L");
        mg->Add(P3pi21_line,"L");
        mg->Add(P3pi22_line,"L");
        mg->Add(P3pi23_line,"L");
        mg->Add(DV1_point,"AP");
        mg->Add(DV2_point,"AP");
        mg->Add(PV_point,"AP");
        mg->Add(origin,"AP");

        mg->Draw("L");
        mg->SetTitle("Projection on plane transverse to K^{+} trajectory");
        //gPad->SetTitle("Projection on plane transverse to K^{+} trajectory");
        mg->GetXaxis()->SetTitle("x (mm)");
        mg->GetYaxis()->SetTitle("y (mm)");

        e1->Draw("same");
        e2->Draw("same");

        // Z AXIS
        c.cd(2);
        gPad->SetGrid(1,1);
        gPad->Draw();
        //gPad->DrawFrame(Zmin, -0.1 , Zmax, 0.1);

        double DV1err_lineZ[2] = {DV1_t.z() - DV1err.z(), DV1_t.z() + DV1err.z()};
        double DV2err_lineZ[2] = {DV2_t.z() - DV2err.z() , DV2_t.z() + DV2err.z()};
        double P3pi1_lineZ[2] = {DV1_t.z(), DV1_t.z() + P3pi11_t.z() + P3pi12_t.z() + P3pi13_t.z()};
        double P3pi2_lineZ[2] = {DV2_t.z(), DV2_t.z() + P3pi21_t.z() + P3pi22_t.z() + P3pi23_t.z()};
        double zero2[2] = {0.,0.};

        TGraph* DV1z_point = new TGraph();
        TGraph* DV2z_point = new TGraph();
        TGraph* PVz_point = new TGraph();
        TGraph* BVz_point = new TGraph();
        TGraph* BVz_true_point = new TGraph();
        DV1z_point->SetPoint(0, DV1_t.z(), 0.);
        DV2z_point->SetPoint(0, DV2_t.z(), 0.);
        PVz_point->SetPoint(0, PV_t.z(), 0.);
        BVz_point->SetPoint(0, BV_t.z(), 0.);
        BVz_true_point->SetPoint(0, BV_true_t.z(), 0.);

        DV1z_point->SetMarkerColor(kBlue);
        DV2z_point->SetMarkerColor(kRed);
        PVz_point->SetMarkerColor(kBlack);
        BVz_point->SetMarkerColor(kGray);
        BVz_true_point->SetMarkerColor(kGray+2);
        DV1z_point->SetMarkerStyle(8);
        DV2z_point->SetMarkerStyle(8);
        PVz_point->SetMarkerStyle(8);
        BVz_point->SetMarkerStyle(8);
        BVz_true_point->SetMarkerStyle(8);

        TGraph* DV1err_line = new TGraph(2, DV1err_lineZ, zero2);
        TGraph* DV2err_line = new TGraph(2, DV2err_lineZ, zero2);

        DV1err_line->SetLineColor(kBlue);
        DV2err_line->SetLineColor(kRed);
        DV1err_line->SetLineWidth(2);
        DV2err_line->SetLineWidth(2);
        DV1err_line->SetLineStyle(7);
        DV2err_line->SetLineStyle(7);

        TMultiGraph* mg1 = new TMultiGraph();
        mg1->Add(DV1z_point,"AP");
        mg1->Add(DV2z_point,"AP");
        mg1->Add(PVz_point,"AP");
        mg1->Add(BVz_point,"AP");
        mg1->Add(BVz_true_point,"AP");
        mg1->Add(DV1err_line,"L");
        mg1->Add(DV2err_line,"L");

        TLatex* tex = new TLatex(0.3, 0.8, Form("DTF #chi^{2}/ndf = %.3lf",DTF_chi2));
        tex->SetNDC(kTRUE);
        tex->SetTextFont(42);
        tex->SetTextSize(0.045);

        TLatex* tex1 = new TLatex(0.3, 0.7, Form("Anti-neutrino P_{z} = %.3lf MeV",Pnu1_t.z()));
        tex1->SetNDC(kTRUE);
        tex1->SetTextFont(42);
        tex1->SetTextSize(0.045);

        TLatex* tex2 = new TLatex(0.3, 0.6, Form("Neutrino P_{z} = %.3lf MeV",Pnu2_t.z()));
        tex2->SetNDC(kTRUE);
        tex2->SetTextFont(42);
        tex2->SetTextSize(0.045);

        TLatex* tex3 = new TLatex(0.3, 0.95, Form("Event number : %i",i));
        tex3->SetNDC(kTRUE);
        tex3->SetTextFont(42);
        tex3->SetTextSize(0.045);

        mg1->Draw();
        tex->Draw("same");
        tex1->Draw("same");
        tex2->Draw("same");
        tex3->Draw("same");
        //mg1->SetTitle("Along K^{+} trajectory");
        mg1->GetXaxis()->SetTitle("z (mm)");

        fout->cd();
        c.Write();

        return 0;
}

ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag){
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Rotation about original Z axis to bring p_K into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(p_K.x()/p_K.y());
  //Clockwise rotation about new X axis to align Z axis with p_K
  double theta = -1 * TMath::ATan(p_K.Rho() * (p_K.y()/TMath::Abs(p_K.y())) /p_K.z());

  ROOT::Math::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  ROOT::Math::Transform3D myTransformation     = ROOT::Math::Transform3D(myRotation) * ROOT::Math::Transform3D(myShift);
  ROOT::Math::Transform3D myTransformation_inv = myTransformation.Inverse();

  ROOT::Math::XYZVector theVector_t;
  if(invFlag)
    theVector_t = myTransformation_inv * theVector;
  else 
    theVector_t = myTransformation * theVector;
  
  return theVector_t;
}

ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint,  ROOT::Math::XYZPoint thePoint, bool invFlag){
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Rotation about original Z axis to bring p_K into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(p_K.x()/p_K.y());
  //Clockwise rotation about new X axis to align Z axis with p_K
  double theta = -1 * TMath::ATan(p_K.Rho() * (p_K.y()/TMath::Abs(p_K.y())) /p_K.z());

  ROOT::Math::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object
  ROOT::Math::Transform3D myTransformation     = ROOT::Math::Transform3D(myRotation) * ROOT::Math::Transform3D(myShift);
  ROOT::Math::Transform3D myTransformation_inv = myTransformation.Inverse();

  ROOT::Math::XYZPoint thePoint_t;
  if(invFlag)
    thePoint_t = myTransformation_inv * thePoint;
  else 
    thePoint_t = myTransformation * thePoint;
  
  return thePoint_t;
}