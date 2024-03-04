ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);
ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint,  ROOT::Math::XYZPoint thePoint, bool invFlag);

void analytic_reconstruction(int year, TString RECO_files, int species){

  Bool_t isMC = false;
  if( (species == 10) || (species == 11) || (species == 12) || (species == 1) )
  {
    isMC = true;
  }

  // is isMC uses the true quantities for m 

  TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files); 
  TChain* t = new TChain("DecayTree");
  t->AddFileInfoList((TCollection*)fc->GetList());

  TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/analytical_reco/201%i/ana_reco_vars_%i.root",year,species), "RECREATE");
  fout->cd();
  TTree* tout = new TTree("DecayTree", "DecayTree");

  // Known parameters
  Double_t Ek, Pkx, Pky, Pkz; // kaon 4-momentum
  Double_t Epi1_1, Ppi1_x_1, Ppi1_y_1, Ppi1_z_1; // 6 pions 4-momentum
  Double_t Epi1_2, Ppi1_x_2, Ppi1_y_2, Ppi1_z_2;
  Double_t Epi1_3, Ppi1_x_3, Ppi1_y_3, Ppi1_z_3;
  Double_t Epi2_1, Ppi2_x_1, Ppi2_y_1, Ppi2_z_1;
  Double_t Epi2_2, Ppi2_x_2, Ppi2_y_2, Ppi2_z_2;
  Double_t Epi2_3, Ppi2_x_3, Ppi2_y_3, Ppi2_z_3;
  Double_t DV1x, DV1y, DV1z; // tau+ decay vertex
  Double_t DV2x, DV2y, DV2z; // tau- decay vertex
  Double_t PVx, PVy, PVz; // primary vertex
  Double_t Rx, Ry, Rz; // point in K+ trajectory

  // Unknown parameters
  Double_t Enu1, Pnu1x, Pnu1y, Pnu1z, nu1_M; // Antineutrino 4-momentum
  Double_t Enu2, Pnu2x, Pnu2y, Pnu2z, nu2_M; // Neutrino 4-momentum
  Double_t Etau1, Ptau1x, Ptau1y, Ptau1z, tau1_M; // tau+ 4-momentum
  Double_t Etau2, Ptau2x, Ptau2y, Ptau2z, tau2_M; // tau- 4-momentum
  Double_t Eb, Pbx, Pby, Pbz, Mb; // B+ 4-momentum
  Double_t BVx, BVy, BVz; // B+ decay vertex
  Double_t Enu1_TRUE, Pnu1x_TRUE, Pnu1y_TRUE, Pnu1z_TRUE; // Antineutrino 4-momentum
  Double_t Enu2_TRUE, Pnu2x_TRUE, Pnu2y_TRUE, Pnu2z_TRUE; // Neutrino 4-momentum
  Double_t Etau1_TRUE, Ptau1x_TRUE, Ptau1y_TRUE, Ptau1z_TRUE; // tau+ 4-momentum
  Double_t Etau2_TRUE, Ptau2x_TRUE, Ptau2y_TRUE, Ptau2z_TRUE; // tau- 4-momentum
  Double_t Eb_TRUE, Pbx_TRUE, Pby_TRUE, Pbz_TRUE; // B+ 4-momentum
  Double_t BVx_TRUE, BVy_TRUE, BVz_TRUE; // B+ decay vertex
  Double_t L1x, L1y, L1z;
  Double_t L2x, L2y, L2z;
  Double_t Lx, Ly, Lz;
  Double_t sol1, sol2;
  Double_t a1, a2, b, c, d, e, f, g, h, i, j, x1, x2, p1, p2, p3, q1, q2, q3;

  // truth information
  Double_t Epi1_1_true, Epi1_2_true, Epi1_3_true;
  Double_t Epi2_1_true, Epi2_2_true, Epi2_3_true;
  Double_t Ppi1x_1_true, Ppi1x_2_true, Ppi1x_3_true;
  Double_t Ppi2x_1_true, Ppi2x_2_true, Ppi2x_3_true;
  Double_t Ppi1y_1_true, Ppi1y_2_true, Ppi1y_3_true;
  Double_t Ppi2y_1_true, Ppi2y_2_true, Ppi2y_3_true;
  Double_t Ppi1z_1_true, Ppi1z_2_true, Ppi1z_3_true;
  Double_t Ppi2z_1_true, Ppi2z_2_true, Ppi2z_3_true;

  // Truth-match
  Int_t Kp_TRUEID, taup_pi1_TRUEID, taup_pi2_TRUEID, taup_pi3_TRUEID, taum_pi1_TRUEID, taum_pi2_TRUEID, taum_pi3_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;

  // Input: measured quantities (PV,DV1,P3pi1,DV2,P3pi2,RP,pK)
  if(isMC)
  {
    // Primary vertex
    t->SetBranchAddress("Bp_TRUEORIGINVERTEX_X", &PVx);
    t->SetBranchAddress("Bp_TRUEORIGINVERTEX_Y", &PVy);
    t->SetBranchAddress("Bp_TRUEORIGINVERTEX_Z", &PVz);
    // tau+ decay vertex:
    t->SetBranchAddress("taup_TRUEENDVERTEX_X", &DV1x);
    t->SetBranchAddress("taup_TRUEENDVERTEX_Y", &DV1y);
    t->SetBranchAddress("taup_TRUEENDVERTEX_Z", &DV1z);
    // tau+ pions kinematics:
    t->SetBranchAddress("taup_pi1_TRUEP_E", &Epi1_1);
    t->SetBranchAddress("taup_pi1_TRUEP_X", &Ppi1_x_1);
    t->SetBranchAddress("taup_pi1_TRUEP_Y", &Ppi1_y_1);
    t->SetBranchAddress("taup_pi1_TRUEP_Z", &Ppi1_z_1);
    t->SetBranchAddress("taup_pi2_TRUEP_E", &Epi1_2);
    t->SetBranchAddress("taup_pi2_TRUEP_X", &Ppi1_x_2);
    t->SetBranchAddress("taup_pi2_TRUEP_Y", &Ppi1_y_2);
    t->SetBranchAddress("taup_pi2_TRUEP_Z", &Ppi1_z_2);
    t->SetBranchAddress("taup_pi3_TRUEP_E", &Epi1_3);
    t->SetBranchAddress("taup_pi3_TRUEP_X", &Ppi1_x_3);
    t->SetBranchAddress("taup_pi3_TRUEP_Y", &Ppi1_y_3);
    t->SetBranchAddress("taup_pi3_TRUEP_Z", &Ppi1_z_3);
    // tau- decay vertex:
    t->SetBranchAddress("taum_TRUEENDVERTEX_X", &DV2x);
    t->SetBranchAddress("taum_TRUEENDVERTEX_Y", &DV2y);
    t->SetBranchAddress("taum_TRUEENDVERTEX_Z", &DV2z);
    // tau- pions kinematics:
    t->SetBranchAddress("taum_pi1_TRUEP_E", &Epi2_1);
    t->SetBranchAddress("taum_pi1_TRUEP_X", &Ppi2_x_1);
    t->SetBranchAddress("taum_pi1_TRUEP_Y", &Ppi2_y_1);
    t->SetBranchAddress("taum_pi1_TRUEP_Z", &Ppi2_z_1);
    t->SetBranchAddress("taum_pi2_TRUEP_E", &Epi2_2);
    t->SetBranchAddress("taum_pi2_TRUEP_X", &Ppi2_x_2);
    t->SetBranchAddress("taum_pi2_TRUEP_Y", &Ppi2_y_2);
    t->SetBranchAddress("taum_pi2_TRUEP_Z", &Ppi2_z_2);
    t->SetBranchAddress("taum_pi3_TRUEP_E", &Epi2_3);
    t->SetBranchAddress("taum_pi3_TRUEP_X", &Ppi2_x_3);
    t->SetBranchAddress("taum_pi3_TRUEP_Y", &Ppi2_y_3);
    t->SetBranchAddress("taum_pi3_TRUEP_Z", &Ppi2_z_3);
    // B+ (true) decay vertex lies in K+ trajectory:
    t->SetBranchAddress("Bp_TRUEENDVERTEX_X", &Rx); 
    t->SetBranchAddress("Bp_TRUEENDVERTEX_Y", &Ry);
    t->SetBranchAddress("Bp_TRUEENDVERTEX_Z", &Rz);
    // K+ kinematics
    t->SetBranchAddress("Kp_TRUEP_E", &Ek);
    t->SetBranchAddress("Kp_TRUEP_X", &Pkx);
    t->SetBranchAddress("Kp_TRUEP_Y", &Pky);
    t->SetBranchAddress("Kp_TRUEP_Z", &Pkz);
  }
  else
  {
    // K+ kinematics
    t->SetBranchAddress("Kp_PE", &Ek);
    t->SetBranchAddress("Kp_PX", &Pkx);
    t->SetBranchAddress("Kp_PY", &Pky);
    t->SetBranchAddress("Kp_PZ", &Pkz);
    // tau+ pions kinematics
    t->SetBranchAddress("taup_pi1_PE", &Epi1_1);
    t->SetBranchAddress("taup_pi1_PX", &Ppi1_x_1);
    t->SetBranchAddress("taup_pi1_PY", &Ppi1_y_1);
    t->SetBranchAddress("taup_pi1_PZ", &Ppi1_z_1);
    t->SetBranchAddress("taup_pi2_PE", &Epi1_2);
    t->SetBranchAddress("taup_pi2_PX", &Ppi1_x_2);
    t->SetBranchAddress("taup_pi2_PY", &Ppi1_y_2);
    t->SetBranchAddress("taup_pi2_PZ", &Ppi1_z_2);
    t->SetBranchAddress("taup_pi3_PE", &Epi1_3);
    t->SetBranchAddress("taup_pi3_PX", &Ppi1_x_3);
    t->SetBranchAddress("taup_pi3_PY", &Ppi1_y_3);
    t->SetBranchAddress("taup_pi3_PZ", &Ppi1_z_3);
    // tau- pions kinematics
    t->SetBranchAddress("taum_pi1_PE", &Epi2_1);
    t->SetBranchAddress("taum_pi1_PX", &Ppi2_x_1);
    t->SetBranchAddress("taum_pi1_PY", &Ppi2_y_1);
    t->SetBranchAddress("taum_pi1_PZ", &Ppi2_z_1);
    t->SetBranchAddress("taum_pi2_PE", &Epi2_2);
    t->SetBranchAddress("taum_pi2_PX", &Ppi2_x_2);
    t->SetBranchAddress("taum_pi2_PY", &Ppi2_y_2);
    t->SetBranchAddress("taum_pi2_PZ", &Ppi2_z_2);
    t->SetBranchAddress("taum_pi3_PE", &Epi2_3);
    t->SetBranchAddress("taum_pi3_PX", &Ppi2_x_3);
    t->SetBranchAddress("taum_pi3_PY", &Ppi2_y_3);
    t->SetBranchAddress("taum_pi3_PZ", &Ppi2_z_3);
    // tau+ decay vertex
    t->SetBranchAddress("taup_ENDVERTEX_X", &DV1x);
    t->SetBranchAddress("taup_ENDVERTEX_Y", &DV1y);
    t->SetBranchAddress("taup_ENDVERTEX_Z", &DV1z);
    // tau- decay vertex
    t->SetBranchAddress("taum_ENDVERTEX_X", &DV2x);
    t->SetBranchAddress("taum_ENDVERTEX_Y", &DV2y);
    t->SetBranchAddress("taum_ENDVERTEX_Z", &DV2z);
    // reference point in K+ trajectory
    t->SetBranchAddress("refPoint_X", &Rx); 
    t->SetBranchAddress("refPoint_Y", &Ry);
    t->SetBranchAddress("refPoint_Z", &Rz);
    // primary vertex
    t->SetBranchAddress("Bp_OWNPV_X", &PVx);
    t->SetBranchAddress("Bp_OWNPV_Y", &PVy);
    t->SetBranchAddress("Bp_OWNPV_Z", &PVz);
  }

  // Saves the unmeasured quantities (BV,pB,mB,Ptau1,Pnu1,Ptau2,Pnu2)
  // B+ decay vertex 
  tout->Branch("BVx", &BVx);
  tout->Branch("BVy", &BVy);
  tout->Branch("BVz", &BVz);
  // B+ 4-momentum + mass
  tout->Branch("Bp_PE", &Eb);
  tout->Branch("Bp_PX", &Pbx);
  tout->Branch("Bp_PY", &Pby);
  tout->Branch("Bp_PZ", &Pbz);
  tout->Branch("Bp_M", &Mb);
  // tau+ 4-momentum
  tout->Branch("taup_PE", &Etau1);
  tout->Branch("taup_PX", &Ptau1x);
  tout->Branch("taup_PY", &Ptau1y);
  tout->Branch("taup_PZ", &Ptau1z);
  tout->Branch("taup_M", &tau1_M);
  // antineutrino 4-momentum
  tout->Branch("antinutau_PE", &Enu1);
  tout->Branch("antinutau_PX", &Pnu1x);
  tout->Branch("antinutau_PY", &Pnu1y);
  tout->Branch("antinutau_PZ", &Pnu1z);
  tout->Branch("antinutau_M", &nu1_M);
  // tau- 4-momentum
  tout->Branch("taum_PE", &Etau2);
  tout->Branch("taum_PX", &Ptau2x);
  tout->Branch("taum_PY", &Ptau2y);
  tout->Branch("taum_PZ", &Ptau2z);
  tout->Branch("taum_M", &tau2_M);
  // nutau 4-momentum
  tout->Branch("nutau_PE", &Enu2);
  tout->Branch("nutau_PX", &Pnu2x);
  tout->Branch("nutau_PY", &Pnu2y);
  tout->Branch("nutau_PZ", &Pnu2z);
  tout->Branch("nutau_M", &nu2_M);

  // scale factors
  tout->Branch("L1x", &L1x);
  tout->Branch("L1y", &L1y);
  tout->Branch("L1z", &L1z);
  tout->Branch("L2x", &L2x);
  tout->Branch("L2y", &L2y);
  tout->Branch("L2z", &L2z);
  tout->Branch("Lx", &Lx);
  tout->Branch("Ly", &Ly);
  tout->Branch("Lz", &Lz);
  tout->Branch("sol1", &sol1);
  tout->Branch("sol2", &sol2);

  tout->Branch("a1", &a1);
  tout->Branch("a2", &a2);
  tout->Branch("b", &b);
  tout->Branch("c", &c);
  tout->Branch("d", &d);
  tout->Branch("e", &e);
  tout->Branch("f", &f);
  tout->Branch("g", &g);
  tout->Branch("h", &h);
  tout->Branch("i", &i);
  tout->Branch("j", &j);
  tout->Branch("x1", &x1);
  tout->Branch("x2", &x2);
  tout->Branch("p1", &p1);
  tout->Branch("p2", &p2);
  tout->Branch("p3", &p3);
  tout->Branch("q1", &q1);
  tout->Branch("q2", &q2);
  tout->Branch("q3", &q3);

  UInt_t num_entries = t->GetEntries();
  for(int evt = 0; evt < num_entries; evt++)
  {
    t->GetEntry(evt);

    ROOT::Math::XYZPoint PV(PVx, PVy, PVz);
    ROOT::Math::XYZPoint DV1(DV1x, DV1y, DV1z);
    ROOT::Math::XYZVector P3pi1(Ppi1_x_1 + Ppi1_x_2 + Ppi1_x_3, Ppi1_y_1 + Ppi1_y_2 + Ppi1_y_3, Ppi1_z_1 + Ppi1_z_2 + Ppi1_z_3);
    ROOT::Math::XYZPoint DV2(DV2x, DV2y, DV2z);
    ROOT::Math::XYZVector P3pi2(Ppi2_x_1 + Ppi2_x_2 + Ppi2_x_3, Ppi2_y_1 + Ppi2_y_2 + Ppi2_y_3, Ppi2_z_1 + Ppi2_z_2 + Ppi2_z_3);
    ROOT::Math::XYZPoint R(Rx, Ry, Rz);
    ROOT::Math::XYZVector Pk(Pkx, Pky, Pkz);

    // Make transformation from LHCb to ANA frame (z-axis along K+ trajectory)
    ROOT::Math::XYZPoint PV_t = makeTransformation_point(Pk, R, PV, false);
    ROOT::Math::XYZPoint DV1_t = makeTransformation_point(Pk, R, DV1, false);
    ROOT::Math::XYZVector P3pi1_t = makeTransformation_vec(Pk, R, P3pi1, false);
    ROOT::Math::XYZPoint DV2_t = makeTransformation_point(Pk, R, DV2, false);
    ROOT::Math::XYZVector P3pi2_t = makeTransformation_vec(Pk, R, P3pi2, false);
    ROOT::Math::XYZPoint R_t = makeTransformation_point(Pk, R, R, false);
    ROOT::Math::XYZVector Pk_t = makeTransformation_vec(Pk, R, Pk, false);

    a1 = (DV1_t.y())/(DV1_t.x());
    a2 = (DV2_t.y())/(DV2_t.x());
    b = (PV_t.y() - a1*PV_t.x())/(a2*PV_t.x() - PV_t.y());
    c = b*(DV1_t.x())/(DV2_t.x());
    d = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()));
    e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
    f = ( PV_t.x()*sqrt(Pk_t.Mag2()) )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
    g = c*e + d;
    h = f*c;
    i = DV1_t.z() - e*DV1_t.x();
    j = f*DV1_t.x();

    Double_t mtau = 1776.86; // MeV (PDG, 2023)
    Double_t E3pi1 = Epi1_1 + Epi1_2 + Epi1_3;
    Double_t E3pi2 = Epi2_1 + Epi2_2 + Epi2_3;
    Double_t m3pi1 = sqrt( pow(E3pi1,2) - P3pi1_t.Mag2() );  
    Double_t m3pi2 = sqrt( pow(E3pi2,2) - P3pi2_t.Mag2() );

    x1 = P3pi1_t.x() + a1*P3pi1_t.y() + e*P3pi1_t.z();
    x2 = b*P3pi2_t.x() + a2*b*P3pi2_t.y() + g*P3pi2_t.z();

    p1 = 1 + pow(a1,2) + pow(e,2) - pow(x1/E3pi1,2);
    p2 = 2*e*f - ( pow(mtau,2) + pow(m3pi1,2) + 2*f*P3pi1_t.z() )*(x1/pow(E3pi1,2) );
    p3 = pow(mtau,2) + pow(f,2) - pow( ( pow(mtau,2) + pow(m3pi1,2) + 2*f*P3pi1_t.z() )/(2*E3pi1), 2);
    q1 = pow(b,2) + pow(a2*b,2) + pow(g,2) - pow(x2/E3pi2,2);
    q2 = 2*g*h - ( pow(mtau,2) + pow(m3pi2,2) + 2*h*P3pi2_t.z() )*(x2/pow(E3pi2,2) );
    q3 = pow(mtau,2) + pow(h,2) - pow( ( pow(mtau,2) + pow(m3pi2,2) + 2*h*P3pi2_t.z() )/(2*E3pi2),2 );

    Double_t Ptau1x_t = (p1*q3 - p3*q1)/(p2*q1 - p1*q2);
    Double_t Ptau1y_t = a1*Ptau1x_t;
    Double_t Ptau1z_t = e*Ptau1x_t + f;
    Double_t Ptau2x_t = b*Ptau1x_t;
    Double_t Ptau2y_t = a2*b*Ptau1x_t;
    Double_t Ptau2z_t = g*Ptau1x_t + h;
    Double_t BVz_t = i - j*(1/Ptau1x_t);

    sol1 = p1*pow(Ptau1x_t,2) + p2*Ptau1x_t + p3;
    sol2 = q1*pow(Ptau1x_t,2) + q2*Ptau1x_t + q3;

    ROOT::Math::XYZPoint BV_t(0, 0, BVz_t);
    ROOT::Math::XYZVector Ptau1_t(Ptau1x_t, Ptau1y_t, Ptau1z_t);
    ROOT::Math::XYZVector Ptau2_t(Ptau2x_t, Ptau2y_t, Ptau2z_t); 
    ROOT::Math::XYZVector Pb_t = Ptau1_t + Ptau2_t + Pk_t;
    ROOT::Math::XYZVector Pnu1_t = Ptau1_t - P3pi1_t;
    ROOT::Math::XYZVector Pnu2_t = Ptau2_t - P3pi2_t;

    // make transformation back to LHCb frame
    ROOT::Math::XYZPoint BV = makeTransformation_point(Pk, R, BV_t, true);
    ROOT::Math::XYZVector Pb = makeTransformation_vec(Pk, R, Pb_t, true);
    ROOT::Math::XYZVector Ptau1 = makeTransformation_vec(Pk, R, Ptau1_t, true);
    ROOT::Math::XYZVector Pnu1 = makeTransformation_vec(Pk, R, Pnu1_t, true);
    ROOT::Math::XYZVector Ptau2 = makeTransformation_vec(Pk, R, Ptau2_t, true);
    ROOT::Math::XYZVector Pnu2 = makeTransformation_vec(Pk, R, Pnu2_t, true);

    L1x = Ptau1_t.x()/DV1_t.x();
    L1y = Ptau1_t.y()/DV1_t.y();
    L1z = Ptau1_t.z()/(DV1_t.z() - BV_t.z());
    L2x = Ptau2_t.x()/DV2_t.x();
    L2y = Ptau2_t.y()/DV2_t.y();
    L2z = Ptau2_t.z()/(DV2_t.z() - BV_t.z());
    Lx = -Pb_t.x()/PV_t.x();
    Ly = -Pb_t.y()/PV_t.y();
    Lz = Pb_t.z()/(BV_t.z() - PV_t.z());

    BVx = BV.x();
    BVy = BV.y();
    BVz = BV.z();
    Pbx = Pb.x();
    Pby = Pb.y();
    Pbz = Pb.z();
    Eb = Etau1 + Etau2 + Ek;
    Mb = sqrt( pow(Eb,2) - Pb.Mag2() );
    Ptau1x = Ptau1.x();
    Ptau1y = Ptau1.y();
    Ptau1z = Ptau1.z();
    Etau1 = sqrt( pow(mtau,2) + Ptau1_t.Mag2() );
    Pnu1x = Pnu1.x();
    Pnu1y = Pnu1.y();
    Pnu1z = Pnu1.z();
    Enu1 = sqrt(Pnu1_t.Mag2());
    Ptau2x = Ptau2.x();
    Ptau2y = Ptau2.y();
    Ptau2z = Ptau2.z();
    Etau2 = sqrt( pow(mtau,2) + Ptau2_t.Mag2() );
    Pnu2x = Pnu2.x();
    Pnu2y = Pnu2.y();
    Pnu2z = Pnu2.z();
    Enu2 = sqrt(Pnu2_t.Mag2());
    tau1_M = sqrt( pow(Etau1,2) - Ptau1_t.Mag2() );
    tau2_M = sqrt( pow(Etau2,2) - Ptau2_t.Mag2() );
    nu1_M = sqrt( pow(Enu1,2) - Pnu1_t.Mag2() );
    nu2_M = sqrt( pow(Enu2,2) - Pnu2_t.Mag2() );
    
    tout->Fill();
  }

  tout->Write();
  fout->Close();

}

ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector p_K, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag){
  // invFlag makes inverse transformation
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
  // invFlag makes inverse transformation
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


