RooRealVar* fit(RooDataSet* data, RooRealVar mass, TString cut, TString name);

using namespace RooStats;
using namespace RooFit;
using namespace std;

std::ofstream file;

void mass_resolution(){

    TFile* fin = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2016/mc_2016.root");
    TTree* t = (TTree*)fin->Get("ntuple/DecayTree");
    TFile* f_kstar = new TFile("/panfs/avenkate/KstTauTau/2016_MD/Stripping/ntuple_TM.root");
    TTree* t_kstar = (TTree*)f_kstar->Get("DecayTreeTuple");

    Double_t taup_FD, taum_FD, taup_FD_chi2, taum_FD_chi2, taup_IP_chi2, taum_IP_chi2;
    Float_t status, Bmass, Enu1, Pnu1x, Pnu1y, Pnu1z, Enu2, Pnu2x, Pnu2y, Pnu2z, DTF_EB, DTF_PBx, DTF_PBy, DTF_PBz;
    Int_t Kp_TRUEID, taup_pip0_TRUEID, taup_pim0_TRUEID, taup_pip1_TRUEID, taum_pim0_TRUEID, taum_pip0_TRUEID, taum_pim1_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;
    Double_t DV1x, DV1y, DV1z, DV2x, DV2y, DV2z, BVx, BVy, BVz, PVx, PVy, PVz, KVx, KVy, KVz;
    Double_t taup_pip0_TRUEPE, taup_pim0_TRUEPE, taup_pip1_TRUEPE, taum_pim0_TRUEPE, taum_pip0_TRUEPE, taum_pim1_TRUEPE, taup_TRUEPE, taum_TRUEPE;
    Double_t taup_pip0_TRUEPX, taup_pim0_TRUEPX, taup_pip1_TRUEPX, taum_pim0_TRUEPX, taum_pip0_TRUEPX, taum_pim1_TRUEPX, taup_TRUEPX, taum_TRUEPX;
    Double_t taup_pip0_TRUEPY, taup_pim0_TRUEPY, taup_pip1_TRUEPY, taum_pim0_TRUEPY, taum_pip0_TRUEPY, taum_pim1_TRUEPY, taup_TRUEPY, taum_TRUEPY;
    Double_t taup_pip0_TRUEPZ, taup_pim0_TRUEPZ, taup_pip1_TRUEPZ, taum_pim0_TRUEPZ, taum_pip0_TRUEPZ, taum_pim1_TRUEPZ, taup_TRUEPZ, taum_TRUEPZ;
    Float_t taum_DTF_M, taup_DTF_M, DTF_chi2, DTF_nIter;
    Double_t Kp_PT, taup_PT, taum_PT, Pkx, Pky, Pkz;
    Float_t M_kstar, status_kstar;

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
    t->SetBranchAddress("Bp_ConsBp_chi2",&DTF_chi2);
    t->SetBranchAddress("Bp_ConsBp_nIter",&DTF_nIter);

    t->SetBranchAddress("Kp_PT",&Kp_PT);
    t->SetBranchAddress("taup_PT",&taup_PT);
    t->SetBranchAddress("taum_PT",&taum_PT);

    t_kstar->SetBranchAddress("B0_kttdtf0_M",&M_kstar);
    t_kstar->SetBranchAddress("B0_kttdtf0_status",&status_kstar);

    TH1D* h_nutau_E = new TH1D("h_nutau_E", "h_nutau_E", 80, -120000., 120000);
    TH1D* h_antinutau_E = new TH1D("h_antinutau_E", "h_antinutau_E", 80, -120000., 120000);
    TH1D* h_nutau_x = new TH1D("h_nutau_x", "h_nutau_x", 80, -120000., 120000);
    TH1D* h_antinutau_x = new TH1D("h_antinutau_x", "h_antinutau_x", 80, -120000., 120000);
    TH1D* h_nutau_y = new TH1D("h_nutau_y", "h_nutau_y", 80, -120000., 120000);
    TH1D* h_antinutau_y = new TH1D("h_antinutau_y", "h_antinutau_y", 80, -120000., 120000);
    TH1D* h_nutau_z_core = new TH1D("h_nutau_z_core", "h_nutau_z_core", 80, -120000., 120000);
    TH1D* h_antinutau_z_core = new TH1D("h_antinutau_z_core", "h_antinutau_z_core", 80, -120000., 120000);
    TH1D* h_nutau_z_tail = new TH1D("h_nutau_z_tail", "h_nutau_z_tail", 80, -120000., 120000);
    TH1D* h_antinutau_z_tail = new TH1D("h_antinutau_z_tail", "h_antinutau_z_tail", 80, -120000., 120000);

    TH1D* h_mass_new = new TH1D("h_mass_new", "h_mass_new", 80, 4000, 8000);

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
    RooDataSet* data8 = new RooDataSet("data8", "data8", mass);
    RooDataSet* data9 = new RooDataSet("data9", "data9", mass);
    RooDataSet* data10 = new RooDataSet("data10", "data10", mass);
    RooDataSet* data11 = new RooDataSet("data11", "data11", mass);
    RooDataSet* data12 = new RooDataSet("data12", "data12", mass);
    RooDataSet* data13 = new RooDataSet("data13", "data13", mass);

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
        DTF_PBz -= (Pnu1z + Pnu2z);
        DTF_PBz += (nutau_TRUEPZ + antinutau_TRUEPZ);

        DTF_EB -= (Enu1 + Enu2);
        DTF_EB += ( sqrt(pow(Pnu1x,2) + pow(Pnu1y,2) + pow(antinutau_TRUEPZ,2)) + sqrt(pow(Pnu2x,2) + pow(Pnu2y,2) + pow(nutau_TRUEPZ,2)) );

        double Bmass_new = sqrt(pow(DTF_EB,2) - pow(DTF_PBx,2) - pow(DTF_PBy,2) - pow(DTF_PBz,2));

        if( (status == 0) && ((abs(Kp_TRUEID) == 321 && abs(taup_pip0_TRUEID) == 211 && abs(taup_pim0_TRUEID) == 211 && abs(taup_pip1_TRUEID) == 211 && abs(taum_pim0_TRUEID) == 211 && abs(taum_pip0_TRUEID) == 211 && abs(taum_pim1_TRUEID) == 211) && (abs(taup_TRUEID) == 15 && abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521))){

            h_mass_new->Fill(Bmass_new);

            h_nutau_E->Fill(Enu2 - nutau_TRUEPE);
            h_antinutau_E->Fill(Enu1 - antinutau_TRUEPE);
            h_nutau_x->Fill(Pnu2x - nutau_TRUEPX);
            h_antinutau_x->Fill(Pnu1x - antinutau_TRUEPX);
            h_nutau_y->Fill(Pnu2y - nutau_TRUEPY);
            h_antinutau_y->Fill(Pnu1y - antinutau_TRUEPY);

            if(abs(Bmass - 5279) < 300){
                h_nutau_z_core->Fill(Pnu2z - nutau_TRUEPZ);
                h_antinutau_z_core->Fill(Pnu1z - antinutau_TRUEPZ);
            }
            if(Bmass > 6500){
                h_nutau_z_tail->Fill(Pnu2z - nutau_TRUEPZ);
                h_antinutau_z_tail->Fill(Pnu1z - antinutau_TRUEPZ);
            }

            if((Bmass > 4000) && (Bmass < 8000)){
                data1->add(mass);

                if(A > 1.){data2->add(mass);}
                if((IP1 > 0.4) && (IP2 > 0.4)){data3->add(mass);}
                if( ((taup_DTF_M > 1750) && (taup_DTF_M < 1800)) && ((taum_DTF_M > 1750) && (taum_DTF_M < 1800)) ){data4->add(mass);}
                if(((P1 < 5000) && (P2 < 5000))){data5->add(mass);}
                if((taup_FD > 1.) && (taum_FD > 1.)){data6->add(mass);}
                if((taup_FD_chi2 > 5) && (taum_FD_chi2 > 5)){data7->add(mass);}
                if(DTF_chi2 < 15){data8->add(mass);}
                if(DTF_nIter < 5){data9->add(mass);}
                if(abs(DV1z - DV2z) > 1){data10->add(mass);}
                if(Kp_PT > 2000){data11->add(mass);}
                if((taum_PT > 3000) && (taup_PT > 3000)){data12->add(mass);}
                if((Pnu2z - nutau_TRUEPZ < 0) && (Pnu1z - antinutau_TRUEPZ < 0)){data13->add(mass);}

            }

        }
    }

    RooRealVar mass_kstar("mass_kstar", "mass_kstar", 4000, 8000, "MeV");
    RooDataSet* data_kstar = new RooDataSet("data_kstar", "data_kstar", mass_kstar);

    for(int i = 0; i < t_kstar->GetEntries(); i++){
        t_kstar->GetEntry(i);

        if(status_kstar == 0){
            if((M_kstar > 4000) && (M_kstar < 8000)){
                mass_kstar.setVal(M_kstar);
                data_kstar->add(mass_kstar);
            }
        }
    }

    RooRealVar* sigma1 = fit(data1, mass, "No cut", "Bmass_cut1");
    RooRealVar* sigma2 = fit(data2, mass, "Area of PV-DV1-DV2 triangle > 1 mm^{2}", "Bmass_cut2");
    RooRealVar* sigma3 = fit(data3, mass, "IP_{#tau} > 0.2 mm", "Bmass_cut3");
    RooRealVar* sigma4 = fit(data4, mass, "DTF #tau mass #in [1750,1800] MeV", "Bmass_cut4");
    RooRealVar* sigma5 = fit(data5, mass, "DTF P_{#nu} < 5000 MeV", "Bmass_cut5");
    RooRealVar* sigma6 = fit(data6, mass, "#tau FD to BV > 1 mm", "Bmass_cut6");
    RooRealVar* sigma7 = fit(data7, mass, "#tau FD to BV #chi^{2} > 5", "Bmass_cut7");
    RooRealVar* sigma8 = fit(data8, mass, "DTF #chi^{2} < 15", "Bmass_cut8");
    RooRealVar* sigma9 = fit(data9, mass, "Number of DTF iterations < 5", "Bmass_cut9");
    RooRealVar* sigma10 = fit(data10, mass, "#tau Z separation > 1 mm", "Bmass_cut10");
    RooRealVar* sigma11 = fit(data11, mass, "K^{+} p_{T} > 2000 MeV", "Bmass_cut11");
    RooRealVar* sigma12 = fit(data12, mass, "#tau p_{T} > 3000 MeV", "Bmass_cut12");
    RooRealVar* sigma13 = fit(data13, mass, "Difference between DTF and TRUE nu Pz < 0 MeV", "Bmass_cut13");
    RooRealVar* sigma_kstar = fit(data_kstar, mass_kstar, "K^{*} mass", "Kstar_mass");

    // FWHM + efficiency
    std::string out = "Mass_resolution_visualisation/fwhm_eff_resolution.tex";
    file.open(out);  

    if(!file.is_open())
        cout << "Output file not opened!";

    file << "\\begin{table}" << std::endl;
    file << Form("\\caption{FWHM and efficiency of $B^{+}$ mass distribution for different cuts. Initial mass resolution is %.2lf $\\pm$ %.2lf MeV.}", sigma1->getVal(), sigma1->getError()) << std::endl;
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
    double n8 = data8->sumEntries();
    double n9 = data9->sumEntries();
    double n10 = data10->sumEntries();
    double n11 = data11->sumEntries();
    double n12 = data12->sumEntries();
    double n13 = data13->sumEntries();

    file << "Cut & " <<  "Resolution (MeV) & " << "Efficiency (\\%)" << "\\\\ \\hline" << std::endl;   
    file << "Area of PV-DV1-DV2 triangle $> 1\\,$mm$^2$ & " <<  Form("%.2lf $\\pm$ %.2lf & ", sigma2->getVal(), sigma2->getError()) << Form("%.2lf ", (n2/n1)*100.) << "\\\\" << std::endl;    
    file << "$IP_{\\tau^{\\pm}} > 0.2$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma3->getVal(), sigma3->getError()) << Form("%.2lf ", (n3/n1)*100.) << "\\\\" << std::endl;    
    file << "DTF $\\tau^{\\pm}$ mass $\\in [1750,1800]$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma4->getVal(), sigma4->getError()) << Form("%.2lf ", (n4/n1)*100.) << "\\\\" << std::endl;    
    file << "DTF $P_{\\nu} < 5000$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma5->getVal(), sigma5->getError()) << Form("%.2lf ", (n5/n1)*100.) << "\\\\" << std::endl;    
    file << "$\\tau^{\\pm}$ FD to BV$ > 1$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma6->getVal(), sigma6->getError()) << Form("%.2lf ", (n6/n1)*100.) << "\\\\" << std::endl;    
    file << "$\\tau^{\\pm}$ FD to BV $\\chi^2 > 5$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma7->getVal(), sigma7->getError()) << Form("%.2lf ", (n7/n1)*100.) << "\\\\" << std::endl;    
    file << "DTF $\\chi^2 < 15$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma8->getVal(), sigma8->getError()) << Form("%.2lf ", (n8/n1)*100.) << "\\\\" << std::endl;    
    file << "Number of DTF iterations $< 5$ & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma9->getVal(), sigma9->getError()) << Form("%.2lf ", (n9/n1)*100.) << "\\\\" << std::endl;    
    file << "$\\tau$ z separation $> 1$\\,mm & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma10->getVal(), sigma10->getError()) << Form("%.2lf ", (n10/n1)*100.) << "\\\\" << std::endl;    
    file << "$K^{+} p_{T} > 2000$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma11->getVal(), sigma11->getError()) << Form("%.2lf ", (n11/n1)*100.) << "\\\\" << std::endl;    
    file << "$\\tau^{\\pm} p_{T} > 3000$\\,MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma12->getVal(), sigma12->getError()) << Form("%.2lf ", (n12/n1)*100.) << "\\\\" << std::endl;  
    file << "Difference btw DTF and TRUE $\\nu$ Pz $< 0$ MeV & " <<  Form("%.2lf $\\pm$ %.2lf & ",  sigma13->getVal(), sigma13->getError()) << Form("%.2lf ", (n13/n1)*100.) << "\\\\ \\hline" << std::endl;    

    file << "\\end{tabular}" << std::endl;
    file << "\\end{table}" << std::endl;
    file.close();

    TCanvas c;
    c.cd();
    h_mass_new->GetXaxis()->SetTitle("mass_new (MeV)");
    h_mass_new->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_mass_new->GetXaxis()->GetXmax() - h_mass_new->GetXaxis()->GetXmin())/h_mass_new->GetNbinsX()) );
    h_mass_new->SetTitle(" ");
    h_mass_new->Draw();
    c.SaveAs("Mass_resolution_visualisation/mass_new.gif");
    c.SaveAs("Mass_resolution_visualisation/mass_new.pdf");

    TCanvas c1E;
    c1E.cd();
    h_nutau_E->GetXaxis()->SetTitle("nutau_DTF_PE - nutau_TRUE_PE (MeV)");
    h_nutau_E->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_E->GetXaxis()->GetXmax() - h_nutau_E->GetXaxis()->GetXmin())/h_nutau_E->GetNbinsX()) );
    h_nutau_E->SetTitle(" ");
    h_nutau_E->Draw();
    c1E.SaveAs("Mass_resolution_visualisation/nutau_E_DTF_TRUE.gif");
    c1E.SaveAs("Mass_resolution_visualisation/nutau_E_DTF_TRUE.pdf");

    TCanvas c2E;
    c2E.cd();
    h_antinutau_E->GetXaxis()->SetTitle("antinutau_DTF_PE - antinutau_TRUE_PE (MeV)");
    h_antinutau_E->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_E->GetXaxis()->GetXmax() - h_antinutau_E->GetXaxis()->GetXmin())/h_antinutau_E->GetNbinsX()) );
    h_antinutau_E->SetTitle(" ");
    h_antinutau_E->Draw();
    c2E.SaveAs("Mass_resolution_visualisation/antinutau_E_DTF_TRUE.gif");
    c2E.SaveAs("Mass_resolution_visualisation/antinutau_E_DTF_TRUE.pdf");

    TCanvas c1x;
    c1x.cd();
    h_nutau_x->GetXaxis()->SetTitle("nutau_DTF_PX - nutau_TRUE_PX (MeV)");
    h_nutau_x->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_x->GetXaxis()->GetXmax() - h_nutau_x->GetXaxis()->GetXmin())/h_nutau_x->GetNbinsX()) );
    h_nutau_x->SetTitle(" ");
    h_nutau_x->Draw();
    c1x.SaveAs("Mass_resolution_visualisation/nutau_x_DTF_TRUE.gif");
    c1x.SaveAs("Mass_resolution_visualisation/nutau_x_DTF_TRUE.pdf");

    TCanvas c2x;
    c2x.cd();
    h_antinutau_x->GetXaxis()->SetTitle("antinutau_DTF_PX - antinutau_TRUE_PX (MeV)");
    h_antinutau_x->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_x->GetXaxis()->GetXmax() - h_antinutau_x->GetXaxis()->GetXmin())/h_antinutau_x->GetNbinsX()) );
    h_antinutau_x->SetTitle(" ");
    h_antinutau_x->Draw();
    c2x.SaveAs("Mass_resolution_visualisation/antinutau_x_DTF_TRUE.gif");
    c2x.SaveAs("Mass_resolution_visualisation/antinutau_x_DTF_TRUE.pdf");

    TCanvas c1y;
    c1y.cd();
    h_nutau_y->GetXaxis()->SetTitle("nutau_DTF_PY - nutau_TRUE_PY (MeV)");
    h_nutau_y->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_nutau_y->GetXaxis()->GetXmax() - h_nutau_y->GetXaxis()->GetXmin())/h_nutau_y->GetNbinsX()) );
    h_nutau_y->SetTitle(" ");
    h_nutau_y->Draw();
    c1y.SaveAs("Mass_resolution_visualisation/nutau_y_DTF_TRUE.gif");
    c1y.SaveAs("Mass_resolution_visualisation/nutau_y_DTF_TRUE.pdf");

    TCanvas c2y;
    c2y.cd();
    h_antinutau_y->GetXaxis()->SetTitle("antinutau_DTF_PY - antinutau_TRUE_PY (MeV)");
    h_antinutau_y->GetYaxis()->SetTitle( TString::Format("Events / (%g MeV)",(h_antinutau_y->GetXaxis()->GetXmax() - h_antinutau_y->GetXaxis()->GetXmin())/h_antinutau_y->GetNbinsX()) );
    h_antinutau_y->SetTitle(" ");
    h_antinutau_y->Draw();
    c2y.SaveAs("Mass_resolution_visualisation/antinutau_y_DTF_TRUE.gif");
    c2y.SaveAs("Mass_resolution_visualisation/antinutau_y_DTF_TRUE.pdf");

    gStyle->SetOptStat(0);
    TCanvas c1z;
    c1z.cd();
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
    c1z.SaveAs("Mass_resolution_visualisation/nutau_DTF_z_TRUE.gif");
    c1z.SaveAs("Mass_resolution_visualisation/nutau_DTF_z_TRUE.pdf");

    TCanvas c2z;
    c2z.cd();
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
    c2z.SaveAs("Mass_resolution_visualisation/antinutau_z_DTF_TRUE.gif");
    c2z.SaveAs("Mass_resolution_visualisation/antinutau_z_DTF_TRUE.pdf");

}

RooRealVar* fit(RooDataSet* data, RooRealVar mass, TString cut, TString name){

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

    c.SaveAs("Mass_resolution_visualisation/"+name+".gif");
    c.SaveAs("Mass_resolution_visualisation/"+name+".pdf");

    return sigma;
}
