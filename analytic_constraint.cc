#include <TH1.h>
#include <TH1D.h>
#include <TMath.h>
#include <RooAbsPdf.h>
#include <iomanip>

using namespace std;

std::vector<double> rotate(double x, double y, double z, double theta, double phi);
std::vector<double> translate(double x, double y, double z, double Px, double Py, double Pz);
std::vector<double> rotate_inverse(double xp, double yp, double zp, double theta, double phi);
std::vector<double> translate_inverse(double x, double y, double z, double Px, double Py, double Pz);

void analytic_constraint(){

  TFile* f = new TFile("./RapidSim/B2Ktautau_tree.root");
  TTree* t = (TTree*)f->Get("DecayTree");

  double PVx, PVy, PVz, DV1x, DV1y, DV1z, DV2x, DV2y, DV2z, Ek, Pkx, Pky, Pkz, Epi_11, Epi_12, Epi_13, Epi_21, Epi_22, Epi_23, Ppix_11, Ppix_12, Ppix_13, Ppix_21, Ppix_22, Ppix_23, Ppiy_11, Ppiy_12, Ppiy_13, Ppiy_21, Ppiy_22, Ppiy_23, Ppiz_11, Ppiz_12, Ppiz_13, Ppiz_21, Ppiz_22, Ppiz_23, B0_M, MCorr, m3pi1_tree, m3pi2_tree, BVgenX, BVgenY, BVgenZ;

  t->SetBranchAddress("Bp_0_origX_TRUE",&PVx);
  t->SetBranchAddress("Bp_0_origY_TRUE",&PVy);
  t->SetBranchAddress("Bp_0_origZ_TRUE",&PVz);
  t->SetBranchAddress("Bp_0_vtxX_TRUE",&BVgenX);
  t->SetBranchAddress("Bp_0_vtxY_TRUE",&BVgenY);
  t->SetBranchAddress("Bp_0_vtxZ_TRUE",&BVgenZ);
  t->SetBranchAddress("taup_0_vtxX_TRUE",&DV1x);
  t->SetBranchAddress("taup_0_vtxY_TRUE",&DV1y);
  t->SetBranchAddress("taup_0_vtxZ_TRUE",&DV1z);
  t->SetBranchAddress("taum_0_vtxX_TRUE",&DV2x);
  t->SetBranchAddress("taum_0_vtxY_TRUE",&DV2y);
  t->SetBranchAddress("taum_0_vtxZ_TRUE",&DV2z);
  t->SetBranchAddress("Kp_0_E_TRUE",&Ek);
  t->SetBranchAddress("Kp_0_PX_TRUE",&Pkx);
  t->SetBranchAddress("Kp_0_PY_TRUE",&Pky);
  t->SetBranchAddress("Kp_0_PZ_TRUE",&Pkz);
  t->SetBranchAddress("pip_0_E_TRUE",&Epi_11);
  t->SetBranchAddress("pip_0_PX_TRUE",&Ppix_11);
  t->SetBranchAddress("pip_0_PY_TRUE",&Ppiy_11);
  t->SetBranchAddress("pip_0_PZ_TRUE",&Ppiz_11);
  t->SetBranchAddress("pim_0_E_TRUE",&Epi_12);
  t->SetBranchAddress("pim_0_PX_TRUE",&Ppix_12);
  t->SetBranchAddress("pim_0_PY_TRUE",&Ppiy_12);
  t->SetBranchAddress("pim_0_PZ_TRUE",&Ppiz_12);
  t->SetBranchAddress("pip_1_E_TRUE",&Epi_13);
  t->SetBranchAddress("pip_1_PX_TRUE",&Ppix_13);
  t->SetBranchAddress("pip_1_PY_TRUE",&Ppiy_13);
  t->SetBranchAddress("pip_1_PZ_TRUE",&Ppiz_13);
  t->SetBranchAddress("pip_2_E_TRUE",&Epi_21);
  t->SetBranchAddress("pip_2_PX_TRUE",&Ppix_21);
  t->SetBranchAddress("pip_2_PY_TRUE",&Ppiy_21);
  t->SetBranchAddress("pip_2_PZ_TRUE",&Ppiz_21);
  t->SetBranchAddress("pim_1_E_TRUE",&Epi_22);
  t->SetBranchAddress("pim_1_PX_TRUE",&Ppix_22);
  t->SetBranchAddress("pim_1_PY_TRUE",&Ppiy_22);
  t->SetBranchAddress("pim_1_PZ_TRUE",&Ppiz_22);
  t->SetBranchAddress("pim_2_E_TRUE",&Epi_23);
  t->SetBranchAddress("pim_2_PX_TRUE",&Ppix_23);
  t->SetBranchAddress("pim_2_PY_TRUE",&Ppiy_23);
  t->SetBranchAddress("pim_2_PZ_TRUE",&Ppiz_23);
  t->SetBranchAddress("Bp_0_M",&B0_M);
  t->SetBranchAddress("MCorr",&MCorr);
  t->SetBranchAddress("pip_0_pim_0_pip_1_M_TRUE",&m3pi1_tree);
  t->SetBranchAddress("pip_2_pim_1_pim_2_M_TRUE",&m3pi2_tree);

  // Unknown parameters  
  TH1D* h_BVx = new TH1D("BVx", "BVx", 100, -5, 5);
  TH1D* h_BVy = new TH1D("BVy", "BVy", 100, -5, 5);
  TH1D* h_BVz = new TH1D("BVz", "BVz", 100, -5, 5);
  TH1D* h_BVz1 = new TH1D("BVz1", "BVz1", 100, -5, 5);
  TH1D* h_BVz2 = new TH1D("BVz2", "BVz2", 100, -5, 5);
  TH1D* h_Etau1 = new TH1D("Etau1", "Etau1", 100, 0, 50);
  TH1D* h_Ptaux1 = new TH1D("Ptaux1", "Ptaux1", 100, -5, 5);
  TH1D* h_Ptauy1 = new TH1D("Ptauy1", "Ptauy1", 100, -5, 5);
  TH1D* h_Ptauz1 = new TH1D("Ptauz1", "Ptauz1", 100, 0, 50);
  TH1D* h_Etau2 = new TH1D("Etau2", "Etau2", 100, 0, 50);
  TH1D* h_Ptaux2 = new TH1D("Ptaux2", "Ptaux2", 100, -5, 5);
  TH1D* h_Ptauy2 = new TH1D("Ptauy2", "Ptauy2", 100, -5, 5);
  TH1D* h_Ptauz2 = new TH1D("Ptauz2", "Ptauz2", 100, 0, 50);
  TH1D* h_Enu1 = new TH1D("Enu1", "Enu1", 100, 0, 20);
  TH1D* h_Pnu1 = new TH1D("Pnu1", "Pnu1", 100, 0, 20);
  TH1D* h_Pnux1 = new TH1D("Pnux1", "Pnux1", 100, -5, 5);
  TH1D* h_Pnuy1 = new TH1D("Pnuy1", "Pnuy1", 100, -5, 5);
  TH1D* h_Pnuz1 = new TH1D("Pnuz1", "Pnuz1", 100, 0, 20);
  TH1D* h_Enu2 = new TH1D("Enu2", "Enu2", 100, 0, 20);
  TH1D* h_Pnu2 = new TH1D("Pnu2", "Pnu2", 100, 0, 20);
  TH1D* h_Pnux2 = new TH1D("Pnux2", "Pnux2", 100, -5, 5);
  TH1D* h_Pnuy2 = new TH1D("Pnuy2", "Pnuy2", 100, -5, 5);
  TH1D* h_Pnuz2 = new TH1D("Pnuz2", "Pnuz2", 100, 0, 20);
  TH1D* h_Eb = new TH1D("Eb", "Eb", 100, 0, 100);
  TH1D* h_Pbx = new TH1D("Pbx", "Pbx", 100, -6, 6);
  TH1D* h_Pby = new TH1D("Pby", "Pby", 100, -6, 6);
  TH1D* h_Pbz = new TH1D("Pbz", "Pbz", 100, 0, 100);

  TH1D* h_BP_M = new TH1D("B0_M", "B0_M", 100, 2, 8);
  TH1D* h_MCorr = new TH1D("MCorr", "MCorr", 100, 2, 8);
  TH1D* h_Mb = new TH1D("Mb", "Mb", 100, 2, 8);
  TH1D* h_Mb1 = new TH1D("Mb1", "Mb1", 100, 2, 8);
  TH1D* h_Pb = new TH1D("Pb", "Pb", 100, 0, 100);
  TH1D* h_Tb = new TH1D("Tb", "Tb", 100, 0, 6); 

  TH1D* h_m3pi1 = new TH1D("m3pi1", "m3pi1", 100, 0, 2);
  TH1D* h_m3pi2 = new TH1D("m3pi2", "m3pi2", 100, 0, 2);
  TH1D* h_m3pi1_tree = new TH1D("m3pi1_tree", "m3pi1_tree", 100, 0, 2);
  TH1D* h_m3pi2_tree = new TH1D("m3pi2_tree", "m3pi2_tree", 100, 0, 2);
  TH1D* h_E3pi1 = new TH1D("E3pi1", "E3pi1", 100, 0, 50);
  TH1D* h_E3pi2 = new TH1D("E3pi2", "E3pi2", 100, 0, 50);
  TH1D* h_P3pi1x = new TH1D("P3pi1x", "P3pi1x", 100, -5, 5);
  TH1D* h_P3pi1y = new TH1D("P3pi1y", "P3pi1y", 100, -4, 7);
  TH1D* h_P3pi1z = new TH1D("P3pi1z", "P3pi1z", 100, 0, 70);
  TH1D* h_P3pi2x = new TH1D("P3pi2x", "P3pi2x", 100, -5, 5);
  TH1D* h_P3pi2y = new TH1D("P3pi2y", "P3pi2y", 100, -4, 7);
  TH1D* h_P3pi2z = new TH1D("P3pi2z", "P3pi2z", 100, 0, 70);

  TH2D* h_ng_hm = new TH2D("ng_hm", "ng_hm", 100, -10, 10, 100, -10, 10);
  TH2D* h_og_im = new TH2D("og_im", "og_im", 100, 0, 100, 100, 0, 100);
  TH1D* h_sqrt = new TH1D("sqrt", "sqrt", 100, -30, 30);
  TH1D* h_sqrt2 = new TH1D("sqrt2", "sqrt2", 100, -30, 30);
  TH2D* h_2d = new TH2D("2d", "2d", 100, 0, 10, 100, 0, 10);
  TH2D* h_hg_nm = new TH2D("hg_nm", "hg_nm", 100, -10, 10, 100, -10, 10);
  TH2D* h_g_m = new TH2D("g_m", "g_m", 100, -50, 50, 100, -50, 50);
  TH2D* h_g_n = new TH2D("g_n", "g_n", 100, -50, 50, 100, -50, 50);
  TH2D* h_g_o = new TH2D("g_o", "g_o", 100, -50, 50, 100, -50, 50);
  TH2D* h_h_m = new TH2D("h_m", "h_m", 100, -50, 50, 100, -50, 50);
  TH2D* h_h_n = new TH2D("h_n", "h_n", 100, -50, 50, 100, -50, 50);
  TH2D* h_h_o = new TH2D("h_o", "h_o", 100, -50, 50, 100, -50, 50);
  TH2D* h_i_m = new TH2D("i_m", "i_m", 100, -50, 50, 100, -50, 50);
  TH2D* h_i_n = new TH2D("i_n", "i_n", 100, -50, 50, 100, -50, 50);
  TH2D* h_i_o = new TH2D("i_o", "i_o", 100, -50, 50, 100, -50, 50);

  TH1D* h_Etau_E3pi1 = new TH1D("Etau_E3pi1", "Etau_E3pi1", 100, -5, 5);
  TH1D* h_Etau_E3pi2 = new TH1D("Etau_E3pi2", "Etau_E3pi2", 100, -5, 5);

  double mtau = 1.77686; // in GeV

  double passed = 0;
  double total = 0;

  for(int evt = 0; evt < t->GetEntries(); evt++){
    t->GetEntry(evt);

    h_m3pi1_tree->Fill(m3pi1_tree);
    h_m3pi2_tree->Fill(m3pi2_tree);
    h_MCorr->Fill(MCorr);
    h_BP_M->Fill(B0_M);

    // Change reference frame (z-axis in K+ trajectory)

    double Pk = sqrt(pow(Pkx,2) + pow(Pky,2) + pow(Pkz,2));
    double Pk2d = sqrt(pow(Pkx,2) + pow(Pky,2));
    double theta = atan(Pky/Pkx);
    double phi = atan(Pk2d/Pkz);

    std::vector<double> PV_trans = translate(PVx, PVy, PVz, BVgenX, BVgenY, BVgenZ);
    std::vector<double> PV = rotate(PV_trans[0], PV_trans[1], PV_trans[2], theta, phi);
    std::vector<double> DV1_trans = translate(DV1x, DV1y, DV1z, BVgenX, BVgenY, BVgenZ);
    std::vector<double> DV1 = rotate(DV1_trans[0], DV1_trans[1], DV1_trans[2], theta, phi);
    std::vector<double> DV2_trans = translate(DV2x, DV2y, DV2z, BVgenX, BVgenY, BVgenZ);
    std::vector<double> DV2 = rotate(DV2_trans[0], DV2_trans[1], DV2_trans[2], theta, phi);

    std::vector<double> Ppi_11 = rotate(Ppix_11, Ppiy_11, Ppiz_11, theta, phi);
    std::vector<double> Ppi_12 = rotate(Ppix_12, Ppiy_12, Ppiz_12, theta, phi);
    std::vector<double> Ppi_13 = rotate(Ppix_13, Ppiy_13, Ppiz_13, theta, phi);
    std::vector<double> Ppi_21 = rotate(Ppix_21, Ppiy_21, Ppiz_21, theta, phi);
    std::vector<double> Ppi_22 = rotate(Ppix_22, Ppiy_22, Ppiz_22, theta, phi);
    std::vector<double> Ppi_23 = rotate(Ppix_23, Ppiy_23, Ppiz_23, theta, phi);

    double PVxp = PV[0];
    double PVyp = PV[1];
    double PVzp = PV[2];
    double DV1xp = DV1[0];
    double DV1yp = DV1[1];
    double DV1zp = DV1[2];
    double DV2xp = DV2[0];
    double DV2yp = DV2[1];
    double DV2zp = DV2[2];
    double Ppixp_11 = Ppi_11[0];
    double Ppiyp_11 = Ppi_11[1];
    double Ppizp_11 = Ppi_11[2];
    double Ppixp_12 = Ppi_12[0];
    double Ppiyp_12 = Ppi_12[1];
    double Ppizp_12 = Ppi_12[2];
    double Ppixp_13 = Ppi_13[0];
    double Ppiyp_13 = Ppi_13[1];
    double Ppizp_13 = Ppi_13[2];
    double Ppixp_21 = Ppi_21[0];
    double Ppiyp_21 = Ppi_21[1];
    double Ppizp_21 = Ppi_21[2];
    double Ppixp_22 = Ppi_22[0];
    double Ppiyp_22 = Ppi_22[1];
    double Ppizp_22 = Ppi_22[2];
    double Ppixp_23 = Ppi_23[0];
    double Ppiyp_23 = Ppi_23[1];
    double Ppizp_23 = Ppi_23[2];

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Convert distances to GeV-1 (not really necessary because only ratios of distances appear)
    // Distances in RapidSim are in cm
   
    PVxp *= (1/0.1975)*pow(10,13);
    PVyp *= (1/0.1975)*pow(10,13);
    PVzp *= (1/0.1975)*pow(10,13);

    DV1xp *= (1/0.1975)*pow(10,13);
    DV1yp *= (1/0.1975)*pow(10,13);
    DV1zp *= (1/0.1975)*pow(10,13);

    DV2xp *= (1/0.1975)*pow(10,13);
    DV2yp *= (1/0.1975)*pow(10,13);
    DV2zp *= (1/0.1975)*pow(10,13);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double p3pi1x = Ppixp_11 + Ppixp_12 + Ppixp_13;
    double p3pi1y = Ppiyp_11 + Ppiyp_12 + Ppiyp_13;
    double p3pi1z = Ppizp_11 + Ppizp_12 + Ppizp_13;
    double p3pi2x = Ppixp_21 + Ppixp_22 + Ppixp_23;
    double p3pi2y = Ppiyp_21 + Ppiyp_22 + Ppiyp_23;
    double p3pi2z = Ppizp_21 + Ppizp_22 + Ppizp_23;

    double E3pi1 = Epi_11 + Epi_12 + Epi_13;
    double E3pi2 = Epi_21 + Epi_22 + Epi_23;
    double p3pi1 = sqrt(pow(p3pi1x,2) + pow(p3pi1y,2) + pow(p3pi1z,2));
    double p3pi2 = sqrt(pow(p3pi2x,2) + pow(p3pi2y,2) + pow(p3pi2z,2));

    double m3pi1 = sqrt(pow(E3pi1,2) - pow(p3pi1,2));
    double m3pi2 = sqrt(pow(E3pi2,2) - pow(p3pi2,2));

    h_m3pi1->Fill(m3pi1); // same as generated (ok)
    h_m3pi2->Fill(m3pi2); // same as generates (ok)
 
    h_E3pi1->Fill(E3pi1);
    h_E3pi2->Fill(E3pi2);
    h_P3pi1x->Fill(p3pi1x);
    h_P3pi1y->Fill(p3pi1y);
    h_P3pi1z->Fill(p3pi1z);
    h_P3pi2x->Fill(p3pi2x);
    h_P3pi2y->Fill(p3pi2y);
    h_P3pi2z->Fill(p3pi2z);

    double a1 = DV1yp/DV1xp;
    double a2 = DV2yp/DV2xp;

    double b = (a1*PVxp - PVyp)/(PVyp - a2*PVxp);
    double c1 = b*((DV2zp - DV1zp)/DV2xp);
    double c2 = b*(DV1xp/DV2xp);
    double d1 = ( (DV1zp - PVzp)*(1 + b) + PVxp*c1 )/( DV1xp*(1 + b) - PVxp*(1 + c2) );
    double d2 = (PVxp*Pk)/( DV1xp*(1 + b) - PVxp*(1 + c2) );
    double e1 = c1 + c2*d1;
    double e2 = c2*d2;

    // expressions are equivalent to the ones found in https://lphe.epfl.ch/publications/theses/these.an.pdf, with a different notation
    double g = 1 + pow(a1,2) + pow(d1,2) - pow((p3pi1x + a1*p3pi1y + d1*p3pi1z)/E3pi1, 2);
    double h = 2*d1*d2 - ((p3pi1x + a1*p3pi1y + d1*p3pi1z)*(pow(mtau,2) + pow(m3pi1,2) + 2*d2*p3pi1z))/(pow(E3pi1,2));
    double i = pow(mtau,2) + pow(d2,2) - pow( (pow(mtau,2) + pow(m3pi1,2) + 2*d2*p3pi1z)/(2*E3pi1), 2);

    double m = pow(b,2) + pow(b,2)*pow(a2,2) + pow(e1,2) - pow( (b*p3pi2x + a2*b*p3pi2y + e1*p3pi2z)/E3pi2, 2);
    double n = 2*e1*e2 - ((b*p3pi2x + a2*b*p3pi2y + e1*p3pi2z)*(pow(mtau,2) + pow(m3pi2,2) + 2*e2*p3pi2z))/(pow(E3pi2,2));
    double o = pow(mtau,2) + pow(e2,2) - pow( (pow(mtau,2) + pow(m3pi2,2) + 2*e2*p3pi2z )/(2*E3pi2), 2);

    h_ng_hm->Fill(n*g,h*m);
    h_og_im->Fill(o*g,i*m);
    h_sqrt->Fill(pow(h,2) - 4*i*g);
    h_sqrt2->Fill(pow(n,2) - 4*o*m);
    h_hg_nm->Fill(h/g, n/m);
    h_g_m->Fill(g,m);
    h_g_n->Fill(g,n);
    h_g_o->Fill(g,o);
    h_h_m->Fill(h,m);
    h_h_n->Fill(h,n);
    h_h_o->Fill(h,o);
    h_i_m->Fill(i,m);
    h_i_n->Fill(i,n);
    h_i_o->Fill(i,o);  

    double y1 = pow((g*o - i*m)/(h*m - g*n),2);
    double y2 = (h*o - i*n)/(h*m - g*n);
    h_2d->Fill(y1,y2);

    total += 1;
    double Ptaux1 = 0;
    double r = 0.001;

    double x1p = (-h + sqrt(pow(h,2) - 4*g*i))/(2*g);
    double x1m = (-h - sqrt(pow(h,2) - 4*g*i))/(2*g);
    double x2p = (-n + sqrt(pow(n,2) - 4*m*o))/(2*m); 
    double x2m = (-n - sqrt(pow(n,2) - 4*m*o))/(2*m);

    if( abs(m*pow(x1p,2) + n*x1p + o) < r ){Ptaux1 = x1p;}
    else if( abs(m*pow(x1m,2) + n*x1m + o) < r ){Ptaux1 = x1m;}
    else if( abs(g*pow(x2p,2) + h*x2p + i) < r ){Ptaux1 = x2p;}
    else if( abs(g*pow(x2m,2) + h*x2m + i) < r ){Ptaux1 = x2m;}
    else{continue;}

    passed += 1;

    
    //cout << g*pow(Ptaux1,2) + h*Ptaux1 + i << endl;

    //double Ptaux1 = (g*o - i*m)/(h*m - g*n);
    double Ptauy1 = a1*Ptaux1;
    double Ptauz1 = d1*Ptaux1 + d2;

    double Ptaux2 = b*Ptaux1;
    double Ptauy2 = a2*b*Ptaux1;
    double Ptauz2 = e1*Ptaux1 + e2;
    //cout << "Ptauz2 1 = " << e1*Ptaux1 + e2 << endl; // (same ok)
    //cout << "Ptauz2 2 = " << c1*Ptaux1 + c2*Ptauz1 << endl; // (same ok)

    double Pnux1 = Ptaux1 - p3pi1x;
    double Pnuy1 = Ptauy1 - p3pi1y;
    double Pnuz1 = Ptauz1 - p3pi1z;

    double Pnux2 = Ptaux2 - p3pi2x;
    double Pnuy2 = Ptauy2 - p3pi2y;
    double Pnuz2 = Ptauz2 - p3pi2z;

    double Pbx = Ptaux1 + Ptaux2;
    double Pby = Ptauy1 + Ptauy2;
    double Pbz = Ptauz1 + Ptauz2 + Pk;

    double BVz = DV1zp - DV1xp*(Ptauz1/Ptaux1);
    double BVz1 = DV2zp - DV2xp*(Ptauz2/Ptaux2);
    double BVz2 = PVzp - PVxp*((Ptauz1 + Ptauz2 + Pk)/(Ptaux1 + Ptaux2));

    // Make reference frame transformation back to LHCb 

    std::vector<double> PTAU1 = rotate_inverse(Ptaux1, Ptauy1, Ptauz1, theta, phi);
    std::vector<double> PTAU2 = rotate_inverse(Ptaux2, Ptauy2, Ptauz2, theta, phi);
    std::vector<double> PNU1 = rotate_inverse(Pnux1, Pnuy1, Pnuz1, theta, phi);
    std::vector<double> PNU2 = rotate_inverse(Pnux2, Pnuy2, Pnuz2, theta, phi);
    std::vector<double> PB = rotate_inverse(Pbx, Pby, Pbz, theta, phi);

    std::vector<double> BV_rot = rotate_inverse(0, 0, BVz, theta, phi);
    std::vector<double> BV = translate_inverse(BV_rot[0], BV_rot[1], BV_rot[2], BVgenX, BVgenY, BVgenZ);
    std::vector<double> BV_rot1 = rotate_inverse(0, 0, BVz1, theta, phi);
    std::vector<double> BV1 = translate_inverse(BV_rot1[0], BV_rot1[1], BV_rot1[2], BVgenX, BVgenY, BVgenZ);
    std::vector<double> BV_rot2 = rotate_inverse(0, 0, BVz2, theta, phi);
    std::vector<double> BV2 = translate_inverse(BV_rot2[0], BV_rot2[1], BV_rot2[2], BVgenX, BVgenY, BVgenZ);

    double PtauX1 = PTAU1[0];
    double PtauY1 = PTAU1[1];
    double PtauZ1 = PTAU1[2];
    double PtauX2 = PTAU2[0];
    double PtauY2 = PTAU2[1];
    double PtauZ2 = PTAU2[2];
    double PnuX1 = PNU1[0];
    double PnuY1 = PNU1[1];
    double PnuZ1 = PNU1[2];
    double PnuX2 = PNU2[0];
    double PnuY2 = PNU2[1];
    double PnuZ2 = PNU2[2];
    double PbX = PB[0];
    double PbY = PB[1];
    double PbZ = PB[2];
    double BVX = BV[0];
    double BVY = BV[1];
    double BVZ = BV[2];
    double BVX1 = BV1[0];
    double BVY1 = BV1[1];
    double BVZ1 = BV1[2];
    double BVX2 = BV2[0];
    double BVY2 = BV2[1];
    double BVZ2 = BV2[2];

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // convert to cm
    double BVX_cm = BVX*0.1975*pow(10,-13);
    double BVY_cm = BVY*0.1975*pow(10,-13);
    double BVZ_cm = BVZ*0.1975*pow(10,-13);
    double BVZ1_cm = BVZ1*0.1975*pow(10,-13);
    double BVZ2_cm = BVZ2*0.1975*pow(10,-13);

    // Scalars
    double Ptau1 = sqrt(pow(PtauX1,2) + pow(PtauY1,2) + pow(PtauZ1,2));
    double Ptau2 = sqrt(pow(PtauX2,2) + pow(PtauY2,2) + pow(PtauZ2,2));
    double Pnu1 = sqrt(pow(PnuX1,2) + pow(PnuY1,2) + pow(PnuZ1,2));
    double Pnu2 = sqrt(pow(PnuX2,2) + pow(PnuY2,2) + pow(PnuZ2,2));
    double Pb = sqrt(pow(PbX,2) + pow(PbY,2) + pow(PbZ,2));

    double Etau1 = sqrt(pow(mtau,2) + pow(Ptau1,2));
    double Etau2 = sqrt(pow(mtau,2) + pow(Ptau2,2));
    double Enu1 = Etau1 - E3pi1;
    double Enu2 = Etau2 - E3pi2;
    double Eb = Etau1 + Etau2 + Ek;
    double Mb = sqrt(pow(Eb,2) - pow(Pb,2));

    // B+ lifetime
    double Dx = BVX - PVx;
    double Dy = BVY - PVy;
    double Dz = BVz - PVz;
    double D = sqrt( pow(Dx,2) + pow(Dy,2) + pow(Dz,2) ); // in GeV^-1
    double Tb = (Mb/Pb)*D*6.59*pow(10,-13); // in ps

    h_Ptaux1->Fill(PtauX1);
    h_Ptauy1->Fill(PtauY1);
    h_Ptauz1->Fill(PtauZ1);
    h_Etau1->Fill(Etau1);
    h_Ptaux2->Fill(PtauX2);
    h_Ptauy2->Fill(PtauY2);
    h_Ptauz2->Fill(PtauZ2);
    h_Etau2->Fill(Etau2);
    h_Pnux1->Fill(PnuX1);
    h_Pnuy1->Fill(PnuY1);
    h_Pnuz1->Fill(PnuZ1);
    h_Enu1->Fill(Enu1);
    h_Pnu1->Fill(Pnu1);
    h_Pnux2->Fill(PnuX2);
    h_Pnuy2->Fill(PnuY2);
    h_Pnuz2->Fill(PnuZ2);
    h_Enu2->Fill(Enu2);
    h_Pnu2->Fill(Pnu2);
    h_Pbx->Fill(PbX);
    h_Pby->Fill(PbY);
    h_Pbz->Fill(PbZ);
    h_Pb->Fill(Pb);
    h_Eb->Fill(Eb);
    h_Mb->Fill(Mb);
    h_Mb1->Fill(Mb);
    h_BVx->Fill(10*BVX_cm);
    h_BVy->Fill(10*BVY_cm);
    h_BVz->Fill(10*BVZ_cm);
    h_BVz1->Fill(10*BVZ1_cm);
    h_BVz2->Fill(10*BVZ2_cm);
    h_Tb->Fill(Tb);

  }

  cout << (passed/total)*100 << " % of events passed" << endl;

  h_Tb->Scale(1/h_Tb->Integral());
  h_Tb->Fit("expo");

  TCanvas a1;
  a1.cd();
  h_ng_hm->GetXaxis()->SetTitle("ng");
  h_ng_hm->GetYaxis()->SetTitle("hm");
  h_ng_hm->Draw();
  a1.SaveAs("./Plots/ng_hm.gif");
  a1.SaveAs("./Plots/ng_hm.pdf");

  TCanvas a2;
  a2.cd();
  h_og_im->GetXaxis()->SetTitle("og");
  h_og_im->GetYaxis()->SetTitle("im");
  h_og_im->Draw();
  a2.SaveAs("./Plots/og_im.gif");
  a2.SaveAs("./Plots/og_im.pdf");

  TCanvas a3;
  a3.cd();
  h_sqrt->GetXaxis()->SetTitle("h^{2} - 4ig");
  h_sqrt->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_sqrt->GetXaxis()->GetXmax() - h_sqrt->GetXaxis()->GetXmin())/h_sqrt->GetNbinsX()) );
  h_sqrt->Draw();
  a3.SaveAs("./Plots/sqrt.gif");
  a3.SaveAs("./Plots/sqrt.pdf");

  TCanvas a4;
  a4.cd();
  h_sqrt2->GetXaxis()->SetTitle("n^{2} - 4om");
  h_sqrt2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_sqrt2->GetXaxis()->GetXmax() - h_sqrt2->GetXaxis()->GetXmin())/h_sqrt2->GetNbinsX()) );
  h_sqrt2->Draw();
  a4.SaveAs("./Plots/sqrt2.gif");
  a4.SaveAs("./Plots/sqrt2.pdf");

  TCanvas a5;
  a5.cd();
  h_Etau_E3pi1->GetXaxis()->SetTitle("E_{#tau 1} - E_{3#pi 1} [GeV]");
  h_Etau_E3pi1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Etau_E3pi1->GetXaxis()->GetXmax() - h_Etau_E3pi1->GetXaxis()->GetXmin())/h_Etau_E3pi1->GetNbinsX()) );
  h_Etau_E3pi1->Draw();  
  a5.SaveAs("./Plots/Etau_E3pi1.gif");
  a5.SaveAs("./Plots/Etau_E3pi1.pdf");

  TCanvas a6;
  a6.cd();
  h_Etau_E3pi2->GetXaxis()->SetTitle("E_{#tau 2} - E_{3#pi 2} [GeV]");
  h_Etau_E3pi2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Etau_E3pi2->GetXaxis()->GetXmax() - h_Etau_E3pi2->GetXaxis()->GetXmin())/h_Etau_E3pi2->GetNbinsX()) );
  h_Etau_E3pi2->Draw();
  a6.SaveAs("./Plots/Etau_E3pi2.gif");
  a6.SaveAs("./Plots/Etau_E3pi2.pdf");

  TCanvas a7;
  a7.cd();
  h_2d->GetXaxis()->SetTitle("y1");
  h_2d->GetYaxis()->SetTitle("y2");
  h_2d->Draw();
  a7.SaveAs("./Plots/x_x2.gif");
  a7.SaveAs("./Plots/x_x2.pdf");  

  TCanvas a8;
  a8.cd();
  h_hg_nm->GetXaxis()->SetTitle("h/g");
  h_hg_nm->GetYaxis()->SetTitle("n/m");
  h_hg_nm->Draw();
  a8.SaveAs("./Plots/hg_nm.gif");
  a8.SaveAs("./Plots/hg_nm.pdf");

  TCanvas a9;
  a9.cd();
  h_g_m->SetTitle("g vs m ");
  h_g_m->GetXaxis()->SetTitle("g [GeV^{-2}]");
  h_g_m->GetYaxis()->SetTitle("m [GeV^{-2}]");
  h_g_m->Draw();
  a9.SaveAs("./Plots/g_m.gif");
  a9.SaveAs("./Plots/g_m.pdf");

  TCanvas a9_1;
  a9_1.cd();
  h_g_n->SetTitle("g vs n");
  h_g_n->GetXaxis()->SetTitle("g [GeV^{-2}]");
  h_g_n->GetYaxis()->SetTitle("n [GeV^{-1}]");
  h_g_n->Draw();
  a9_1.SaveAs("./Plots/g_n.gif");
  a9_1.SaveAs("./Plots/g_n.pdf");

  TCanvas a9_2;
  a9_2.cd();
  h_g_o->SetTitle("g vs o");
  h_g_o->GetXaxis()->SetTitle("g [GeV^{-2}]");
  h_g_o->GetYaxis()->SetTitle("o");
  h_g_o->Draw();
  a9_2.SaveAs("./Plots/g_o.gif");
  a9_2.SaveAs("./Plots/g_o.pdf");

  TCanvas a10;
  a10.cd();
  h_h_m->SetTitle("h vs m");
  h_h_m->GetXaxis()->SetTitle("h [GeV^{-1}]");
  h_h_m->GetYaxis()->SetTitle("m [GeV^{-2}]");
  h_h_m->Draw();
  a10.SaveAs("./Plots/h_m.gif");
  a10.SaveAs("./Plots/h_m.pdf");

  TCanvas a10_1;
  a10_1.cd();
  h_h_n->SetTitle("h vs n");
  h_h_n->GetXaxis()->SetTitle("h [GeV^{-1}]");
  h_h_n->GetYaxis()->SetTitle("n [GeV^{-1}]");
  h_h_n->Draw();
  a10_1.SaveAs("./Plots/h_n.gif");
  a10_1.SaveAs("./Plots/h_n.pdf");

  TCanvas a10_2;
  a10_2.cd();
  h_h_o->SetTitle("h vs o");
  h_h_o->GetXaxis()->SetTitle("h [GeV^{-1}]");
  h_h_o->GetYaxis()->SetTitle("o");
  h_h_o->Draw();
  a10_2.SaveAs("./Plots/h_o.gif");
  a10_2.SaveAs("./Plots/h_o.pdf");

  TCanvas a11;
  a11.cd();
  h_i_m->SetTitle("i vs m");
  h_i_m->GetXaxis()->SetTitle("i");
  h_i_m->GetYaxis()->SetTitle("m [GeV^{-2}]");
  h_i_m->Draw();
  a11.SaveAs("./Plots/i_m.gif");
  a11.SaveAs("./Plots/i_m.pdf");

  TCanvas a11_1;
  a11_1.cd();
  h_i_n->SetTitle("i vs n");
  h_i_n->GetXaxis()->SetTitle("i");
  h_i_n->GetYaxis()->SetTitle("n [GeV^{-1}]");
  h_i_n->Draw();
  a11_1.SaveAs("./Plots/i_n.gif");
  a11_1.SaveAs("./Plots/i_n.pdf");

  TCanvas a11_2;
  a11_2.cd();
  h_i_o->SetTitle("i vs o");
  h_i_o->GetXaxis()->SetTitle("i");
  h_i_o->GetYaxis()->SetTitle("o");
  h_i_o->Draw();
  a11_2.SaveAs("./Plots/i_o.gif");
  a11_2.SaveAs("./Plots/i_o.pdf");

  TCanvas b;
  b.cd();
  h_m3pi1->SetFillColorAlpha(kBlue, 0.35);
  h_m3pi1_tree->SetFillColorAlpha(kRed, 0.35);
  h_m3pi1->GetXaxis()->SetTitle("M_{3#pi1} [GeV]");
  h_m3pi1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_m3pi1->GetXaxis()->GetXmax() - h_m3pi1->GetXaxis()->GetXmin())/h_m3pi1->GetNbinsX()) );
  h_m3pi1->Draw();
  h_m3pi1_tree->Draw("same");
  b.SaveAs("./Plots/m3pi1.gif");
  b.SaveAs("./Plots/m3pi1.pdf");

  TCanvas b1;
  b1.cd();
  h_m3pi2->SetFillColorAlpha(kBlue, 0.35);
  h_m3pi2_tree->SetFillColorAlpha(kRed, 0.35);
  h_m3pi2->GetXaxis()->SetTitle("M_{3#pi2} [GeV]");
  h_m3pi2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_m3pi2->GetXaxis()->GetXmax() - h_m3pi2->GetXaxis()->GetXmin())/h_m3pi2->GetNbinsX()) );
  h_m3pi2->Draw();
  h_m3pi2_tree->Draw("same");
  b1.SaveAs("./Plots/m3pi2.gif");
  b1.SaveAs("./Plots/m3pi2.pdf");

  TCanvas b2;
  //b2.SetLogy();
  b2.cd();
  h_Pb->GetXaxis()->SetTitle("P_{B} [GeV]");
  h_Pb->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pb->GetXaxis()->GetXmax() - h_Pb->GetXaxis()->GetXmin())/h_Pb->GetNbinsX()) );
  h_Pb->Draw();
  b2.SaveAs("./Plots/Pb.gif");
  b2.SaveAs("./Plots/Pb.pdf");

  TCanvas b3;
  b3.cd();
  h_E3pi1->GetXaxis()->SetTitle("E_{3#pi1} [GeV]");
  h_E3pi1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_E3pi1->GetXaxis()->GetXmax() - h_E3pi1->GetXaxis()->GetXmin())/h_E3pi1->GetNbinsX()) );
  h_E3pi1->Draw();
  b3.SaveAs("./Plots/E3pi1.gif");
  b3.SaveAs("./Plots/E3pi1.pdf");

  TCanvas b4;
  b4.cd();
  h_P3pi1x->GetXaxis()->SetTitle("P_{3#pi1x} [GeV]");
  h_P3pi1x->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_P3pi1x->GetXaxis()->GetXmax() - h_P3pi1x->GetXaxis()->GetXmin())/h_P3pi1x->GetNbinsX()) );
  h_P3pi1x->Draw();
  b4.SaveAs("./Plots/P3pi1x.gif");
  b4.SaveAs("./Plots/P3pi1x.pdf");

  TCanvas b5;
  b5.cd();
  h_P3pi1y->GetXaxis()->SetTitle("P_{3#pi1y} [GeV]");
  h_P3pi1y->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_P3pi1y->GetXaxis()->GetXmax() - h_P3pi1y->GetXaxis()->GetXmin())/h_P3pi1y->GetNbinsX()) );
  h_P3pi1y->Draw();
  b5.SaveAs("./Plots/P3pi1y.gif");
  b5.SaveAs("./Plots/P3pi1y.pdf");

  TCanvas b6;
  b6.cd();
  h_P3pi1z->GetXaxis()->SetTitle("P_{3#pi1z} [GeV]");
  h_P3pi1z->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_P3pi1z->GetXaxis()->GetXmax() - h_P3pi1z->GetXaxis()->GetXmin())/h_P3pi1z->GetNbinsX()) );
  h_P3pi1z->Draw();
  b6.SaveAs("./Plots/P3pi1z.gif");
  b6.SaveAs("./Plots/P3pi1z.pdf");

  TCanvas b7;
  b7.cd();
  h_E3pi2->GetXaxis()->SetTitle("E_{3#pi2} [GeV]");
  h_E3pi2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_E3pi2->GetXaxis()->GetXmax() - h_E3pi2->GetXaxis()->GetXmin())/h_E3pi2->GetNbinsX()) );
  h_E3pi2->Draw();
  b7.SaveAs("./Plots/E3pi2.gif");
  b7.SaveAs("./Plots/E3pi2.pdf");

  TCanvas b8;
  b8.cd();
  h_P3pi2x->GetXaxis()->SetTitle("P_{3#pi2x} [GeV]");
  h_P3pi2x->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_P3pi2x->GetXaxis()->GetXmax() - h_P3pi2x->GetXaxis()->GetXmin())/h_P3pi2x->GetNbinsX()) );
  h_P3pi2x->Draw();
  b8.SaveAs("./Plots/P3pi2x.gif");
  b8.SaveAs("./Plots/P3pi2x.pdf");

  TCanvas b9;
  b9.cd();
  h_P3pi2y->GetXaxis()->SetTitle("P_{3#pi2y} [GeV]");
  h_P3pi2y->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_P3pi2y->GetXaxis()->GetXmax() - h_P3pi2y->GetXaxis()->GetXmin())/h_P3pi2y->GetNbinsX()) );
  h_P3pi2y->Draw();
  b9.SaveAs("./Plots/P3pi2y.gif");
  b9.SaveAs("./Plots/P3pi2y.pdf");

  TCanvas b10;
  b10.cd();
  h_P3pi2z->GetXaxis()->SetTitle("P_{3#pi2z} [GeV]");
  h_P3pi2z->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_P3pi2z->GetXaxis()->GetXmax() - h_P3pi2z->GetXaxis()->GetXmin())/h_P3pi2z->GetNbinsX()) );
  h_P3pi2z->Draw();
  b10.SaveAs("./Plots/P3pi2z.gif");
  b10.SaveAs("./Plots/P3pi2z.pdf");

  TLegend *leg = new TLegend (0.7,0.5,0.88,0.68);
  leg->AddEntry(h_BP_M, "Visible" ,"fp");
  leg->AddEntry(h_MCorr, "Corrected" ,"fp");
  leg->AddEntry(h_Mb, "Analytically reco." ,"fp");

  TCanvas c;
  c.cd();
  gStyle->SetOptStat(0);
  h_BP_M->SetTitle("");
  h_MCorr->SetTitle("");
  h_Mb->SetTitle("");
  h_Mb->Scale(1/h_Mb->Integral());
  h_BP_M->Scale(1/h_BP_M->Integral());
  h_MCorr->Scale(1/h_MCorr->Integral());
  h_BP_M->SetFillColorAlpha(kRed, 0.35);
  h_MCorr->SetFillColorAlpha(kGreen, 0.35);
  h_Mb->SetFillColorAlpha(kBlue, 0.35);
  h_MCorr->GetXaxis()->SetTitle("M_{B} [GeV]");
  h_MCorr->GetYaxis()->SetTitle( "Normalized entries" );
  h_MCorr->GetYaxis()->SetRangeUser(0,0.1);
  h_MCorr->Draw("HIST");
  h_BP_M->Draw("HIST same");
  h_Mb->Draw("HISTsame");
  leg->Draw("same");
  c.SaveAs("./Plots/B_mass.gif");
  c.SaveAs("./Plots/B_mass.pdf");
  gStyle->SetOptStat(1);

  TCanvas c0;
  c0.cd();
  h_Mb1->GetXaxis()->SetTitle("M_{B} [GeV]");
  h_Mb1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Mb1->GetXaxis()->GetXmax() - h_Mb1->GetXaxis()->GetXmin())/h_Mb1->GetNbinsX()) );
  h_Mb1->Draw();
  c0.SaveAs("./Plots/B_my_mass.gif");
  c0.SaveAs("./Plots/B_my_mass.pdf");

  TCanvas c1;
  //c1.SetLogy();
  c1.cd();
  h_BVz->GetXaxis()->SetTitle("BVz [mm]");
  h_BVz->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_BVz->GetXaxis()->GetXmax() - h_BVz->GetXaxis()->GetXmin())/h_BVz->GetNbinsX()) );
  h_BVz->SetFillColorAlpha(kBlue, 0.35);
  h_BVz1->SetFillColorAlpha(kRed, 0.35);
  h_BVz2->SetFillColorAlpha(kGreen, 0.35);
  h_BVz->Draw();
  h_BVz1->Draw("same");
  h_BVz2->Draw("same");
  c1.SaveAs("./Plots/BVz.gif");
  c1.SaveAs("./Plots/BVz.pdf");
 
  TCanvas c1_2;
  c1_2.cd();
  h_BVx->GetXaxis()->SetTitle("BVx [mm]");
  h_BVx->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_BVx->GetXaxis()->GetXmax() - h_BVx->GetXaxis()->GetXmin())/h_BVx->GetNbinsX()) );
  h_BVx->Draw();
  c1_2.SaveAs("./Plots/BVx.gif");
  c1_2.SaveAs("./Plots/BVx.pdf");

  TCanvas c1_3;
  c1_3.cd();
  h_BVy->GetXaxis()->SetTitle("BVy [mm]");
  h_BVy->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_BVy->GetXaxis()->GetXmax() - h_BVy->GetXaxis()->GetXmin())/h_BVy->GetNbinsX()) );
  h_BVy->Draw();
  c1_3.SaveAs("./Plots/BVy.gif");
  c1_3.SaveAs("./Plots/BVy.pdf");

  TCanvas c2;
  c2.cd();
  h_Ptaux1->GetXaxis()->SetTitle("Px_{#tau1} [GeV]");
  h_Ptaux1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptaux1->GetXaxis()->GetXmax() - h_Ptaux1->GetXaxis()->GetXmin())/h_Ptaux1->GetNbinsX()) );
  h_Ptaux1->Draw();
  c2.SaveAs("./Plots/Ptau1x.gif");
  c2.SaveAs("./Plots/Ptau1x.pdf");

  TCanvas c3;
  c3.cd();
  h_Ptauy1->GetXaxis()->SetTitle("Py_{#tau1} [GeV]");
  h_Ptauy1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptauy1->GetXaxis()->GetXmax() - h_Ptauy1->GetXaxis()->GetXmin())/h_Ptauy1->GetNbinsX()) );
  h_Ptauy1->Draw();
  c3.SaveAs("./Plots/Ptau1y.gif");
  c3.SaveAs("./Plots/Ptau1y.pdf");
  
  TCanvas c4;
  //c4.SetLogy();
  c4.cd();
  h_Ptauz1->GetXaxis()->SetTitle("Pz_{#tau1} [GeV]");
  h_Ptauz1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptauz1->GetXaxis()->GetXmax() - h_Ptauz1->GetXaxis()->GetXmin())/h_Ptauz1->GetNbinsX()) );
  h_Ptauz1->Draw();
  c4.SaveAs("./Plots/Ptau1z.gif");
  c4.SaveAs("./Plots/Ptau1z.pdf");

  TCanvas c5;
  //c5.SetLogy();
  c5.cd();
  h_Etau1->GetXaxis()->SetTitle("E_{#tau1} [GeV]");
  h_Etau1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Etau1->GetXaxis()->GetXmax() - h_Etau1->GetXaxis()->GetXmin())/h_Etau1->GetNbinsX()) );
  h_Etau1->Draw();
  c5.SaveAs("./Plots/Etau1.gif");
  c5.SaveAs("./Plots/Etau1.pdf");

  TCanvas c6;
  c6.cd();
  h_Ptaux2->GetXaxis()->SetTitle("Px_{#tau2} [GeV]");
  h_Ptaux2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptaux2->GetXaxis()->GetXmax() - h_Ptaux2->GetXaxis()->GetXmin())/h_Ptaux2->GetNbinsX()) );
  h_Ptaux2->Draw();
  c6.SaveAs("./Plots/Ptau2x.gif");
  c6.SaveAs("./Plots/Ptau2x.pdf");

  TCanvas c7;
  c7.cd();
  h_Ptauy2->GetXaxis()->SetTitle("Py_{#tau2} [GeV]");
  h_Ptauy2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptauy2->GetXaxis()->GetXmax() - h_Ptauy2->GetXaxis()->GetXmin())/h_Ptauy2->GetNbinsX()) );
  h_Ptauy2->Draw();
  c7.SaveAs("./Plots/Ptau2y.gif");
  c7.SaveAs("./Plots/Ptau2y.pdf");

  TCanvas c8;
  //c8.SetLogy();
  c8.cd();
  h_Ptauz2->GetXaxis()->SetTitle("Pz_{#tau2} [GeV]");
  h_Ptauz2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptauz2->GetXaxis()->GetXmax() - h_Ptauz2->GetXaxis()->GetXmin())/h_Ptauz2->GetNbinsX()) );
  h_Ptauz2->Draw();
  c8.SaveAs("./Plots/Ptau2z.gif");
  c8.SaveAs("./Plots/Ptau2z.pdf");

  TCanvas c9;
  //c9.SetLogy();
  c9.cd();
  h_Etau2->GetXaxis()->SetTitle("E_{#tau2} [GeV]");
  h_Etau2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Etau2->GetXaxis()->GetXmax() - h_Etau2->GetXaxis()->GetXmin())/h_Etau2->GetNbinsX()) );
  h_Etau2->Draw();
  c9.SaveAs("./Plots/Etau2.gif");
  c9.SaveAs("./Plots/Etau2.pdf");

  TCanvas c10;
  //c10.SetLogy();
  c10.cd();
  h_Pnux1->GetXaxis()->SetTitle("Px_{#nu1} [GeV]");
  h_Pnux1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pnux1->GetXaxis()->GetXmax() - h_Pnux1->GetXaxis()->GetXmin())/h_Pnux1->GetNbinsX()) );
  h_Pnux1->Draw();
  c10.SaveAs("./Plots/Pnu1x.gif");
  c10.SaveAs("./Plots/Pnu1x.pdf");

  TCanvas c11;
  //c11.SetLogy();
  c11.cd();
  h_Pnuy1->GetXaxis()->SetTitle("Py_{#nu1} [GeV]");
  h_Pnuy1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pnuy1->GetXaxis()->GetXmax() - h_Pnuy1->GetXaxis()->GetXmin())/h_Pnuy1->GetNbinsX()) );
  h_Pnuy1->Draw();
  c11.SaveAs("./Plots/Pnu1y.gif");
  c11.SaveAs("./Plots/Pnu1y.pdf");

  TCanvas c12;
  //c12.SetLogy();
  c12.cd();
  h_Pnuz1->GetXaxis()->SetTitle("Pz_{#nu1} [GeV]");
  h_Pnuz1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pnuz1->GetXaxis()->GetXmax() - h_Pnuz1->GetXaxis()->GetXmin())/h_Pnuz1->GetNbinsX()) );
  h_Pnuz1->Draw();
  c12.SaveAs("./Plots/Pnu1z.gif");
  c12.SaveAs("./Plots/Pnu1z.pdf");

  TCanvas c13;
  c13.cd();
  h_Enu1->SetFillColorAlpha(kRed, 0.35);
  h_Pnu1->SetFillColorAlpha(kBlue, 0.35);
  h_Enu1->GetXaxis()->SetTitle("E_{#nu1} [GeV]");
  h_Enu1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Enu1->GetXaxis()->GetXmax() - h_Enu1->GetXaxis()->GetXmin())/h_Enu1->GetNbinsX()) );
  h_Enu1->Draw();
  h_Pnu1->Draw("same");
  c13.SaveAs("./Plots/Enu1.gif");
  c13.SaveAs("./Plots/Enu1.pdf");

  TCanvas c14;
  //c14.SetLogy();
  c14.cd();
  h_Pnux2->GetXaxis()->SetTitle("Px_{#nu2} [GeV]");
  h_Pnux2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pnux2->GetXaxis()->GetXmax() - h_Pnux2->GetXaxis()->GetXmin())/h_Pnux2->GetNbinsX()) );
  h_Pnux2->Draw();
  c14.SaveAs("./Plots/Pnu2x.gif");
  c14.SaveAs("./Plots/Pnu2x.pdf");

  TCanvas c15;
  //c15.SetLogy();
  c15.cd();
  h_Pnuy2->GetXaxis()->SetTitle("Py_{#nu2} [GeV]");
  h_Pnuy2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pnuy2->GetXaxis()->GetXmax() - h_Pnuy2->GetXaxis()->GetXmin())/h_Pnuy2->GetNbinsX()) );
  h_Pnuy2->Draw();
  c15.SaveAs("./Plots/Pnu2y.gif");
  c15.SaveAs("./Plots/Pnu2y.pdf");

  TCanvas c16;
  //c16.SetLogy();
  c16.cd();
  h_Pnuz2->GetXaxis()->SetTitle("Pz_{#nu2} [GeV]");
  h_Pnuz2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pnuz2->GetXaxis()->GetXmax() - h_Pnuz2->GetXaxis()->GetXmin())/h_Pnuz2->GetNbinsX()) );
  h_Pnuz2->Draw();
  c16.SaveAs("./Plots/Pnu2z.gif");
  c16.SaveAs("./Plots/Pnu2z.pdf");

  TCanvas c17;
  //c17.SetLogy();
  c17.cd();
  h_Enu2->SetFillColorAlpha(kRed, 0.35);
  h_Pnu2->SetFillColorAlpha(kBlue, 0.35);
  h_Enu2->GetXaxis()->SetTitle("E_{#nu2} [GeV]");
  h_Enu2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Enu2->GetXaxis()->GetXmax() - h_Enu2->GetXaxis()->GetXmin())/h_Enu2->GetNbinsX()) );
  h_Enu2->Draw();
  h_Pnu2->Draw("same");
  c17.SaveAs("./Plots/Enu2.gif");
  c17.SaveAs("./Plots/Enu2.pdf");

  TCanvas c18;
  c18.cd();
  h_Pbx->GetXaxis()->SetTitle("Px_{B} [GeV]");
  h_Pbx->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pbx->GetXaxis()->GetXmax() - h_Pbx->GetXaxis()->GetXmin())/h_Pbx->GetNbinsX()) );
  h_Pbx->Draw();
  c18.SaveAs("./Plots/Pbx.gif");
  c18.SaveAs("./Plots/Pbx.pdf");

  TCanvas c19;
  c19.cd();
  h_Pby->GetXaxis()->SetTitle("Py_{B} [GeV]");
  h_Pby->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pby->GetXaxis()->GetXmax() - h_Pby->GetXaxis()->GetXmin())/h_Pby->GetNbinsX()) );
  h_Pby->Draw();
  c19.SaveAs("./Plots/Pby.gif");
  c19.SaveAs("./Plots/Pby.pdf");

  TCanvas c20;
  //c20.SetLogy();
  c20.cd();
  h_Pbz->GetXaxis()->SetTitle("Pz_{B} [GeV]");
  h_Pbz->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pbz->GetXaxis()->GetXmax() - h_Pbz->GetXaxis()->GetXmin())/h_Pbz->GetNbinsX()) );
  h_Pbz->Draw();
  c20.SaveAs("./Plots/Pbz.gif");
  c20.SaveAs("./Plots/Pbz.pdf");

  TCanvas c21;
  //c21.SetLogy();
  c21.cd();
  h_Eb->GetXaxis()->SetTitle("E_{B} [GeV]");
  h_Eb->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Eb->GetXaxis()->GetXmax() - h_Eb->GetXaxis()->GetXmin())/h_Eb->GetNbinsX()) );
  h_Eb->Draw();
  c21.SaveAs("./Plots/Eb.gif");
  c21.SaveAs("./Plots/Eb.pdf");

  TCanvas a;
  a.cd();
  gStyle->SetOptFit(1111);
  h_Tb->SetTitle("B^{+} lifetime");
  h_Tb->GetXaxis()->SetTitle("#tau_{B} [ps]");
  h_Tb->GetYaxis()->SetTitle("p(#tau_{B})");
  h_Tb->Draw();
  a.SaveAs("./Plots/Tb.gif");
  a.SaveAs("./Plots/Tb.pdf");

  return;
 
}


std::vector<double> rotate(double x, double y, double z, double theta, double phi){

  std::vector<double> rotated;

  double xp = cos(phi)*cos(theta)*x + cos(phi)*sin(theta)*y + sin(phi)*z;   
  double yp = -sin(theta)*x + cos(theta)*y;
  double zp = sin(phi)*cos(theta)*x + sin(phi)*sin(theta)*y - cos(phi)*z;

  rotated.push_back(xp);
  rotated.push_back(yp);
  rotated.push_back(zp);

  return rotated;

}

std::vector<double> translate(double x, double y, double z, double Px, double Py, double Pz){

  std::vector<double> translated;

  double xp = x - Px;    
  double yp = y - Py;
  double zp = z - Pz;

  translated.push_back(xp);
  translated.push_back(yp);
  translated.push_back(zp);

  return translated;

}

std::vector<double> translate_inverse(double xp, double yp, double zp, double Px, double Py, double Pz){

  std::vector<double> translated;

  double x = xp + Px;
  double y = yp + Py;
  double z = zp + Pz;

  translated.push_back(x);
  translated.push_back(y);
  translated.push_back(z);

  return translated;

}

std::vector<double> rotate_inverse(double xp, double yp, double zp, double theta, double phi){

  std::vector<double> rotated;

  double x = cos(phi)*cos(theta)*xp - sin(theta)*yp + sin(phi)*cos(theta)*zp;
  double y = cos(phi)*sin(theta)*xp + cos(theta)*yp + sin(phi)*sin(theta)*zp;
  double z = sin(phi)*xp - cos(phi)*zp;

  rotated.push_back(x);
  rotated.push_back(y);
  rotated.push_back(z);

  return rotated;

}


