#include <TH1.h>
#include <TH1D.h>

void analytic_constraint(){

  TFile* f = new TFile("./RapidSim/B2Ktautau_tree.root");
  TTree* t = (TTree*)f->Get("DecayTree");

  // Known parameters
  double PVx, PVy, PVz, DV1x, DV1y, DV1z, DV2x, DV2y, DV2z, Ek, Pkx, Pky, Pkz, Epi_11, Epi_12, Epi_13, Epi_21, Epi_22, Epi_23, Ppix_11, Ppix_12, Ppix_13, Ppix_21, Ppix_22, Ppix_23, Ppiy_11, Ppiy_12, Ppiy_13, Ppiy_21, Ppiy_22, Ppiy_23, Ppiz_11, Ppiz_12, Ppiz_13, Ppiz_21, Ppiz_22, Ppiz_23, B0_M, MCorr, m3pi1_tree, m3pi2_tree;

  double BVX, BVY, BVZ, PBX, PBY, PBZ, mtau1, mtau2;

  t->SetBranchAddress("Bp_0_PX",&PBX);
  t->SetBranchAddress("Bp_0_PY",&PBY);
  t->SetBranchAddress("Bp_0_PZ",&PBZ);
  t->SetBranchAddress("Bp_0_vtxX",&BVX);
  t->SetBranchAddress("Bp_0_vtxY",&BVY);
  t->SetBranchAddress("Bp_0_vtxZ",&BVZ);
  t->SetBranchAddress("Bp_0_origX",&PVx);
  t->SetBranchAddress("Bp_0_origY",&PVy);
  t->SetBranchAddress("Bp_0_origZ",&PVz);
  t->SetBranchAddress("taup_0_vtxX",&DV1x);
  t->SetBranchAddress("taup_0_vtxY",&DV1y);
  t->SetBranchAddress("taup_0_vtxZ",&DV1z);
  t->SetBranchAddress("taum_0_vtxX",&DV2x);
  t->SetBranchAddress("taum_0_vtxY",&DV2y);
  t->SetBranchAddress("taum_0_vtxZ",&DV2z);
  t->SetBranchAddress("Kp_0_E",&Ek);
  t->SetBranchAddress("Kp_0_PX",&Pkx);
  t->SetBranchAddress("Kp_0_PY",&Pky);
  t->SetBranchAddress("Kp_0_PZ",&Pkz);
  t->SetBranchAddress("pip_0_E",&Epi_11);
  t->SetBranchAddress("pip_0_PX",&Ppix_11);
  t->SetBranchAddress("pip_0_PY",&Ppiy_11);
  t->SetBranchAddress("pip_0_PZ",&Ppiz_11);
  t->SetBranchAddress("pim_0_E",&Epi_12);
  t->SetBranchAddress("pim_0_PX",&Ppix_12);
  t->SetBranchAddress("pim_0_PY",&Ppiy_12);
  t->SetBranchAddress("pim_0_PZ",&Ppiz_12);
  t->SetBranchAddress("pip_1_E",&Epi_13);
  t->SetBranchAddress("pip_1_PX",&Ppix_13);
  t->SetBranchAddress("pip_1_PY",&Ppiy_13);
  t->SetBranchAddress("pip_1_PZ",&Ppiz_13);
  t->SetBranchAddress("pip_2_E",&Epi_21);
  t->SetBranchAddress("pip_2_PX",&Ppix_21);
  t->SetBranchAddress("pip_2_PY",&Ppiy_21);
  t->SetBranchAddress("pip_2_PZ",&Ppiz_21);
  t->SetBranchAddress("pim_1_E",&Epi_22);
  t->SetBranchAddress("pim_1_PX",&Ppix_22);
  t->SetBranchAddress("pim_1_PY",&Ppiy_22);
  t->SetBranchAddress("pim_1_PZ",&Ppiz_22);
  t->SetBranchAddress("pim_2_E",&Epi_23);
  t->SetBranchAddress("pim_2_PX",&Ppix_23);
  t->SetBranchAddress("pim_2_PY",&Ppiy_23);
  t->SetBranchAddress("pim_2_PZ",&Ppiz_23);
  t->SetBranchAddress("Bp_0_M",&B0_M);
  t->SetBranchAddress("MCorr",&MCorr);
  t->SetBranchAddress("pip_0_pim_0_pip_1_M",&m3pi1_tree);
  t->SetBranchAddress("pip_2_pim_1_pim_2_M",&m3pi2_tree);
  t->SetBranchAddress("taup_0_M_TRUE",&mtau1);
  t->SetBranchAddress("taum_0_M_TRUE",&mtau2);

  // Unknown parameters  
  TH1D* h_BVz = new TH1D("BVz", "BVz", 100, -1, 1);
  TH1D* h_BVz1 = new TH1D("BVz1", "BVz1", 100, -1, 1);
  TH1D* h_BVz2 = new TH1D("BVz2", "BVz2", 100, -1, 1);
  TH1D* h_Etau1 = new TH1D("Etau1", "Etau1", 100, 0, 150);
  TH1D* h_Ptaux1 = new TH1D("Ptaux1", "Ptaux1", 100, -2, 2);
  TH1D* h_Ptauy1 = new TH1D("Ptauy1", "Ptauy1", 100, -2, 2);
  TH1D* h_Ptauz1 = new TH1D("Ptauz1", "Ptauz1", 100, 0, 150);
  TH1D* h_Etau2 = new TH1D("Etau2", "Etau2", 100, 0, 150);
  TH1D* h_Ptaux2 = new TH1D("Ptaux2", "Ptaux2", 100, -2, 2);
  TH1D* h_Ptauy2 = new TH1D("Ptauy2", "Ptauy2", 100, -2, 2);
  TH1D* h_Ptauz2 = new TH1D("Ptauz2", "Ptauz2", 100, 0, 150);
  TH1D* h_Enu1 = new TH1D("Enu1", "Enu1", 100, 0, 200);
  TH1D* h_Pnux1 = new TH1D("Pnux1", "Pnux1", 100, -10, 10);
  TH1D* h_Pnuy1 = new TH1D("Pnuy1", "Pnuy1", 100, -10, 10);
  TH1D* h_Pnuz1 = new TH1D("Pnuz1", "Pnuz1", 100, -200, 5);
  TH1D* h_Enu2 = new TH1D("Enu2", "Enu2", 100, 0, 500);
  TH1D* h_Pnux2 = new TH1D("Pnux2", "Pnux2", 100, -10, 10);
  TH1D* h_Pnuy2 = new TH1D("Pnuy2", "Pnuy2", 100, -10, 10);
  TH1D* h_Pnuz2 = new TH1D("Pnuz2", "Pnuz2", 100, -200, 5);
  TH1D* h_Eb = new TH1D("Eb", "Eb", 100, 0, 300);
  TH1D* h_Pbx = new TH1D("Pbx", "Pbx", 100, -2, 2);
  TH1D* h_Pby = new TH1D("Pby", "Pby", 100, -2, 2);
  TH1D* h_Pbz = new TH1D("Pbz", "Pbz", 100, 0, 1500);

  TH1D* h_BP_M = new TH1D("B0_M", "B0_M", 100, 2, 8);
  TH1D* h_MCorr = new TH1D("MCorr", "MCorr", 100, 2, 8);
  TH1D* h_Mb = new TH1D("Mb", "Mb", 100, 2, 8);
  TH1D* h_Mb1 = new TH1D("Mb1", "Mb1", 100, 2, 20);
  TH1D* h_Pb = new TH1D("Pb", "Pb", 100, 0, 300);
 
  TH1D* h_m3pi1 = new TH1D("m3pi1", "m3pi1", 100, 0, 2);
  TH1D* h_m3pi2 = new TH1D("m3pi2", "m3pi2", 100, 0, 2);
  TH1D* h_m3pi1_tree = new TH1D("m3pi1_tree", "m3pi1_tree", 100, 0, 2);
  TH1D* h_m3pi2_tree = new TH1D("m3pi2_tree", "m3pi2_tree", 100, 0, 2);

  TH1D* h_Mtau1 = new TH1D("Mtau1", "Mtau1", 100, 1.5, 2.);
  TH1D* h_Mtau2 = new TH1D("Mtau2", "Mtau2", 100, 1.5, 2.);
  TH1D* h_Mk = new TH1D("Mk", "Mk", 100, 0, 1.);

  TH1D* h_Btau = new TH1D("Btau", "Btau", 100, 0, 6.);
  TH1D* h_Btau1 = new TH1D("Btau1", "Btau1", 100, 0, 6.);
  TH1D* h_BD = new TH1D("BD", "BD", 100, 0, 10);

  //double mtau = 1.77686; // in GeV

  for(int evt = 0; evt < t->GetEntries(); evt++){
    t->GetEntry(evt);

    h_m3pi1_tree->Fill(m3pi1_tree);
    h_m3pi2_tree->Fill(m3pi2_tree);

    // B+ lifetime
    double Dx = BVX - PVx;
    double Dy = BVY - PVy;
    double Dz = BVZ - PVz;
    double D = sqrt(pow(Dx,2) + pow(Dy,2) + pow(Dz,2));
    double P = sqrt(pow(PBX,2) + pow(PBY,2) + pow(PBZ,2));
    double T = (5.279/P)*D;

    h_Btau->Fill(T);
    h_BD->Fill(D);
    // [t] = ps -> [d] = 0.1975/6.59 cm = (1/6.59)*pow(10,13) GeV^-1

    // Change reference frame (z-axis along K momentum) + move PV to (1,0,0) 

    double Pk = sqrt(pow(Pkx,2) + pow(Pky,2) + pow(Pkz,2));
    double Pk2d = sqrt(pow(Pkx,2) + pow(Pky,2));

    double Mk = sqrt(pow(Ek,2) - pow(Pk,2));
    h_Mk->Fill(Mk);

    h_MCorr->Fill(MCorr);
    h_BP_M->Fill(B0_M);

    double PVxp = -(Pky/Pk2d)*PVx + (Pkx/Pk2d)*PVy;
    double PVyp = -((Pkz*Pkx)/(Pk*Pk2d))*PVx - ((Pky*Pkz)/(Pk*Pk2d))*PVy + (Pk2d/Pk)*PVz;
    double PVzp = (Pkx*PVx + Pky*PVy + Pkz*PVz)/Pk;
    double DV1xp = -(Pky/Pk2d)*DV1x + (Pkx/Pk2d)*DV1y;
    double DV1yp = -((Pkz*Pkx)/(Pk*Pk2d))*DV1x - ((Pky*Pkz)/(Pk*Pk2d))*DV1y + (Pk2d/Pk)*DV1z;
    double DV1zp = (Pkx*DV1x + Pky*DV1y + Pkz*DV1z)/Pk; 
    double DV2xp = -(Pky/Pk2d)*DV2x + (Pkx/Pk2d)*DV2y;
    double DV2yp = -((Pkz*Pkx)/(Pk*Pk2d))*DV2x - ((Pky*Pkz)/(Pk*Pk2d))*DV2y + (Pk2d/Pk)*DV2z;
    double DV2zp = (Pkx*DV2x + Pky*DV2y + Pkz*DV2z)/Pk;
    double Ppixp_11 = -(Pky/Pk2d)*Ppix_11 + (Pkx/Pk2d)*Ppiy_11;
    double Ppiyp_11 = -((Pkz*Pkx)/(Pk*Pk2d))*Ppix_11 - ((Pky*Pkz)/(Pk*Pk2d))*Ppiy_11 + (Pk2d/Pk)*Ppiz_11;
    double Ppizp_11 = (Pkx*Ppix_11 + Pky*Ppiy_11 + Pkz*Ppiz_11)/Pk;
    double Ppixp_12 = -(Pky/Pk2d)*Ppix_12 + (Pkx/Pk2d)*Ppiy_12;
    double Ppiyp_12 = -((Pkz*Pkx)/(Pk*Pk2d))*Ppix_12 - ((Pky*Pkz)/(Pk*Pk2d))*Ppiy_12 + (Pk2d/Pk)*Ppiz_12;
    double Ppizp_12 = (Pkx*Ppix_12 + Pky*Ppiy_12 + Pkz*Ppiz_12)/Pk;
    double Ppixp_13 = -(Pky/Pk2d)*Ppix_13 + (Pkx/Pk2d)*Ppiy_13;
    double Ppiyp_13 = -((Pkz*Pkx)/(Pk*Pk2d))*Ppix_13 - ((Pky*Pkz)/(Pk*Pk2d))*Ppiy_13 + (Pk2d/Pk)*Ppiz_13;
    double Ppizp_13 = (Pkx*Ppix_13 + Pky*Ppiy_13 + Pkz*Ppiz_13)/Pk;
    double Ppixp_21 = -(Pky/Pk2d)*Ppix_21 + (Pkx/Pk2d)*Ppiy_21;
    double Ppiyp_21 = -((Pkz*Pkx)/(Pk*Pk2d))*Ppix_21 - ((Pky*Pkz)/(Pk*Pk2d))*Ppiy_21 + (Pk2d/Pk)*Ppiz_21;
    double Ppizp_21 = (Pkx*Ppix_21 + Pky*Ppiy_21 + Pkz*Ppiz_21)/Pk;
    double Ppixp_22 = -(Pky/Pk2d)*Ppix_22 + (Pkx/Pk2d)*Ppiy_22;
    double Ppiyp_22 = -((Pkz*Pkx)/(Pk*Pk2d))*Ppix_22 - ((Pky*Pkz)/(Pk*Pk2d))*Ppiy_22 + (Pk2d/Pk)*Ppiz_22;
    double Ppizp_22 = (Pkx*Ppix_22 + Pky*Ppiy_22 + Pkz*Ppiz_22)/Pk;
    double Ppixp_23 = -(Pky/Pk2d)*Ppix_23 + (Pkx/Pk2d)*Ppiy_23;
    double Ppiyp_23 = -((Pkz*Pkx)/(Pk*Pk2d))*Ppix_23 - ((Pky*Pkz)/(Pk*Pk2d))*Ppiy_23 + (Pk2d/Pk)*Ppiz_23;
    double Ppizp_23 = (Pkx*Ppix_23 + Pky*Ppiy_23 + Pkz*Ppiz_23)/Pk;
    double Pkxp = -(Pky/Pk2d)*Pkx + (Pkx/Pk2d)*Pky;
    double Pkyp = -((Pkz*Pkx)/(Pk*Pk2d))*Pkx - ((Pky*Pkz)/(Pk*Pk2d))*Pky + (Pk2d/Pk)*Pkz;
    double Pkzp = (Pkx*Pkx + Pky*Pky + Pkz*Pkz)/Pk;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Convert distances to GeV-1

    PVxp *= (1/6.59)*pow(10,13);
    PVyp *= (1/6.59)*pow(10,13);
    PVzp *= (1/6.59)*pow(10,13);

    DV1xp *= (1/6.59)*pow(10,13);
    DV1yp *= (1/6.59)*pow(10,13);
    DV1zp *= (1/6.59)*pow(10,13);

    DV2xp *= (1/6.59)*pow(10,13);
    DV2yp *= (1/6.59)*pow(10,13);
    DV2zp *= (1/6.59)*pow(10,13);

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

    double a1 = DV1yp/DV1xp;
    double a2 = DV2yp/DV2xp;

    double b = (a1*PVxp - PVyp)/(PVyp - a2*PVxp);
    double c1 = b*(DV2zp - DV1zp)/DV2xp;
    double c2 = b*DV1xp/DV2xp;
    double d1 = ( (DV1zp - PVzp)*(1+b) + PVxp*c1 )/( DV1xp*(1+b) - PVxp*(1+c2) );
    double d2 = (PVxp*Pk)/( DV1xp*(1+b) - PVxp*(1+c2) );
    double e1 = c1 + c2*d1;
    double e2 = c2*d2;

    double f1 = p3pi1x + a1*p3pi1y + d1*p3pi1z;
    double f2 = d2*p3pi1z;
    double j1 = b*p3pi2x + a2*b*p3pi2y + e1*p3pi2z;
    double j2 = e2*p3pi2z;

    double g = 1 + pow(a1,2) + pow(d1,2) - pow(f1,2)/pow(E3pi1,2);
    double h = 2*d1*d2 - ( pow(mtau1,2) + pow(m3pi1,2) )*(f1/pow(E3pi1,2)) - (2*f1*f2)/pow(E3pi1,2);
    double i = pow(mtau1,2) + pow(d2,2) - pow( pow(mtau1,2) + pow(m3pi1,2) ,2)/(4*pow(E3pi1,2)) - ( pow(mtau1,2) + pow(m3pi1,2) )*(f2/pow(E3pi1,2)) - pow(f2,2)/pow(E3pi1,2);

    double m = pow(b,2)*(1 + pow(a2,2)) + pow(e1,2) - pow(j1,2)/pow(E3pi2,2);
    double n = 2*e1*e2 - ( pow(mtau2,2) + pow(m3pi2,2) )*(j1/pow(E3pi2,2)) - (2*j1*j2)/pow(E3pi2,2);
    double o = pow(mtau2,2) + pow(e2,2) - pow( pow(mtau2,2) + pow(m3pi2,2) ,2)/(4*pow(E3pi2,2)) - ( pow(mtau2,2) + pow(m3pi2,2) )*(j2/pow(E3pi2,2)) - pow(j2,2)/pow(E3pi2,2);

    double Ptaux1 = (i*m - o*g)/(g*n - h*m);
    double Ptauy1 = a1*Ptaux1;
    double Ptauz1 = d1*Ptaux1 + d2;
    double Ptau1 = sqrt(pow(Ptaux1,2) + pow(Ptauy1,2) + pow(Ptauz1,2));
    double Etau1 = sqrt(pow(mtau1,2) + pow(Ptau1,2));

    h_Ptaux1->Fill(Ptaux1);
    h_Ptauy1->Fill(Ptauy1);
    h_Ptauz1->Fill(Ptauz1);
    h_Etau1->Fill(Etau1);

    double Ptaux2 = b*Ptaux1;
    double Ptauy2 = a2*b*Ptaux1;
    double Ptauz2 = e1*Ptaux1 + e2;
    //cout << "Ptauz2 1 = " << e1*Ptaux1 + e2 << endl; // (same ok)
    //cout << "Ptauz2 2 = " << c1*Ptaux1 + c2*Ptauz1 << endl; // (same ok)
    double Ptau2 = sqrt(pow(Ptaux2,2) + pow(Ptauy2,2) + pow(Ptauz2,2)); 
    double Etau2 = sqrt(pow(mtau2,2) + pow(Ptau2,2));

    h_Ptaux2->Fill(Ptaux2);
    h_Ptauy2->Fill(Ptauy2);
    h_Ptauz2->Fill(Ptauz2);
    h_Etau2->Fill(Etau2);

    double Pnux1 = Ptaux1 - p3pi1x;
    double Pnuy1 = Ptauy1 - p3pi1y;
    double Pnuz1 = Ptauz1 - p3pi1z;
    double Enu1 = sqrt(pow(Pnux1,2) + pow(Pnuy1,2) + pow(Pnuz1,2));

    h_Pnux1->Fill(Pnux1);
    h_Pnuy1->Fill(Pnuy1);
    h_Pnuz1->Fill(Pnuz1);
    h_Enu1->Fill(Enu1);

    double Pnux2 = Ptaux2 - p3pi2x;
    double Pnuy2 = Ptauy2 - p3pi2y;
    double Pnuz2 = Ptauz2 - p3pi2z;
    double Enu2 = sqrt(pow(Pnux2,2) + pow(Pnuy2,2) + pow(Pnuz2,2));

    h_Pnux2->Fill(Pnux2);
    h_Pnuy2->Fill(Pnuy2);
    h_Pnuz2->Fill(Pnuz2);
    h_Enu2->Fill(Enu2);

    double Pbx = Ptaux1 + Ptaux2;
    double Pby = Ptauy1 + Ptauy2;
    double Pbz = Ptauz1 + Ptauz2 + Pk;
    double Pb = sqrt(pow(Pbx,2) + pow(Pby,2) + pow(Pbz,2));
    double Eb = Etau1 + Etau2 + Ek;
    double Mb = sqrt(pow(Eb,2) - pow(Pb,2));
    cout << Mb << endl;

    h_Pbx->Fill(Pbx);
    h_Pby->Fill(Pby);
    h_Pbz->Fill(Pbz);
    h_Pb->Fill(Pb);
    h_Eb->Fill(Eb);
    h_Mb->Fill(Mb);
    h_Mb1->Fill(Mb);

    double BVz = DV1zp - DV1xp*(Ptauz1/Ptaux1);
    double BVz1 = DV2zp - DV2xp*(Ptauz2/Ptaux2);
    double BVz2 = PVzp - PVxp*((Ptauz1 + Ptauz2 + Pk)/(Ptaux1 + Ptaux2));

    // convert to mm
    double BVz_mm = BVz*0.1975*pow(10,-12);
    double BVz1_mm = BVz1*0.1975*pow(10,-12);
    double BVz2_mm = BVz2*0.1975*pow(10,-12);

    //cout << "BVz 1 = " << BVz_mm << endl; // (same ok)
    //cout << "BVz 2 = " << BVz1_mm << endl; // (same ok)
    //cout << "BVz 3 = " << BVz2_mm << endl; // (same ok)

    h_BVz->Fill(BVz_mm);
    h_BVz1->Fill(BVz1_mm);
    h_BVz2->Fill(BVz2_mm);

    // Analitically reconstructed lifetime
    double Dxp = -PVxp*0.1975*pow(10,-12);
    double Dyp = -PVyp*0.1975*pow(10,-12);
    double Dzp = (BVz2 - PVzp)*0.1975*pow(10,-12);
    double Dp = sqrt(pow(Dxp,2) + pow(Dyp,2) + pow(Dzp,2));
    double Tp = (Mb/Pb)*Dp;
    h_Btau1->Fill(Tp);

  }

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
  b2.SetLogy();
  b2.cd();
  h_Pb->GetXaxis()->SetTitle("P_{B} [GeV]");
  h_Pb->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Pb->GetXaxis()->GetXmax() - h_Pb->GetXaxis()->GetXmin())/h_Pb->GetNbinsX()) );
  h_Pb->Draw();
  b2.SaveAs("./Plots/Pb.gif");
  b2.SaveAs("./Plots/Pb.pdf");

  TLegend *leg = new TLegend (0.7,0.5,0.88,0.68);
  leg->AddEntry(h_BP_M, "Visible" ,"fp");
  leg->AddEntry(h_MCorr, "Corrected" ,"fp");
  leg->AddEntry(h_Mb, "Analytically reco." ,"fp");

  TCanvas c;
  c.cd();
  h_BP_M->SetTitle("");
  h_MCorr->SetTitle("");
  h_Mb->SetTitle("");
  h_BP_M->SetFillColorAlpha(kRed, 0.35);
  h_MCorr->SetFillColorAlpha(kGreen, 0.35);
  h_Mb->SetFillColorAlpha(kBlue, 0.35);
  h_MCorr->GetXaxis()->SetTitle("M_{B} [GeV]");
  h_MCorr->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_MCorr->GetXaxis()->GetXmax() - h_MCorr->GetXaxis()->GetXmin())/h_MCorr->GetNbinsX()) );
  h_MCorr->Draw();
  h_BP_M->Draw("same");
  h_Mb->Draw("same");
  leg->Draw("same");
  c.SaveAs("./Plots/B_mass.gif");
  c.SaveAs("./Plots/B_mass.pdf");

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
  c4.SetLogy();
  c4.cd();
  h_Ptauz1->GetXaxis()->SetTitle("Pz_{#tau1} [GeV]");
  h_Ptauz1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptauz1->GetXaxis()->GetXmax() - h_Ptauz1->GetXaxis()->GetXmin())/h_Ptauz1->GetNbinsX()) );
  h_Ptauz1->Draw();
  c4.SaveAs("./Plots/Ptau1z.gif");
  c4.SaveAs("./Plots/Ptau1z.pdf");

  TCanvas c5;
  c5.SetLogy();
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
  c8.SetLogy();
  c8.cd();
  h_Ptauz2->GetXaxis()->SetTitle("Pz_{#tau2} [GeV]");
  h_Ptauz2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Ptauz2->GetXaxis()->GetXmax() - h_Ptauz2->GetXaxis()->GetXmin())/h_Ptauz2->GetNbinsX()) );
  h_Ptauz2->Draw();
  c8.SaveAs("./Plots/Ptau2z.gif");
  c8.SaveAs("./Plots/Ptau2z.pdf");

  TCanvas c9;
  c9.SetLogy();
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
  h_Enu1->GetXaxis()->SetTitle("E_{#nu1} [GeV]");
  h_Enu1->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Enu1->GetXaxis()->GetXmax() - h_Enu1->GetXaxis()->GetXmin())/h_Enu1->GetNbinsX()) );
  h_Enu1->Draw();
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
  h_Enu2->GetXaxis()->SetTitle("E_{#nu2} [GeV]");
  h_Enu2->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Enu2->GetXaxis()->GetXmax() - h_Enu2->GetXaxis()->GetXmin())/h_Enu2->GetNbinsX()) );
  h_Enu2->Draw();
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
  c21.SetLogy();
  c21.cd();
  h_Eb->GetXaxis()->SetTitle("E_{B} [GeV]");
  h_Eb->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Eb->GetXaxis()->GetXmax() - h_Eb->GetXaxis()->GetXmin())/h_Eb->GetNbinsX()) );
  h_Eb->Draw();
  c21.SaveAs("./Plots/Eb.gif");
  c21.SaveAs("./Plots/Eb.pdf");

  TCanvas d;
  d.cd();
  h_Btau->SetFillColorAlpha(kRed, 0.35);
  h_Btau1->SetFillColorAlpha(kGreen, 0.35);
  h_Btau->GetXaxis()->SetTitle("#tau_{B} [ps]");
  h_Btau->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_Btau->GetXaxis()->GetXmax() - h_Btau->GetXaxis()->GetXmin())/h_Btau->GetNbinsX()) );
  h_Btau->Draw();
  h_Btau1->Draw("same");
  d.SaveAs("./Plots/Tb.gif");
  d.SaveAs("./Plots/Tb.pdf");

  TCanvas d1;
  d1.cd();
  h_BD->GetXaxis()->SetTitle("D_{B} [mm]");
  h_BD->GetYaxis()->SetTitle( TString::Format("Events / (%g)",(h_BD->GetXaxis()->GetXmax() - h_BD->GetXaxis()->GetXmin())/h_BD->GetNbinsX()) );
  h_BD->Draw();
  d1.SaveAs("./Plots/Db.gif");
  d1.SaveAs("./Plots/Db.pdf");

}
