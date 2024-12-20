#include<fstream>
#include <iostream>   // std::cout
#include <string>
#include<algorithm>
#include <typeinfo>
#include <map>
#include <typeinfo>

bool is_digits(const std::string &str);
bool isInside(const std::string & str, char c);
std::map<TString,map<TString,std::vector<double>>>table_reader(TString year, TString Magnet);
Double_t fraction_zero( Double_t E , Int_t region);
Double_t energy_overlap( Double_t E , Double_t X , Double_t Y , Int_t region1 , Int_t region2 );
void ReadFile( FILE * fp , Double_t * array , Int_t * nbins );
void ReadTables(int year, Bool_t isMagUp);

void L0Hadron_TOS_trigger_corr(int year, int species, int line){

  TString files;
  if(species == 1)
  {
    files = Form("Files_on_grid/MC_201%i.txt",year);
  }
  else if(species == 7)
  {
    files = Form("Files_on_grid/MC_D0Dps_201%i.txt",year);
  }
  else
  {
    cout << "Not a valid species name" << endl;
    return;
  }
  TFileCollection* fc = new TFileCollection("fc", "fc", files, 1, line);
  TChain* t = new TChain("ntuple/DecayTree");
  t->AddFileInfoList((TCollection*)fc->GetList());

  TString MC_component_file;
  if(species == 1)
  {
      MC_component_file = Form("/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201%i/%i.root",year,line);
      t->AddFriend("DecayTree", MC_component_file);
  }

  string filename;
  ifstream file(files);
  if(file.is_open())
  {
    for(int i = 0; i < line; i++)
    {
      getline(file,filename);
    }
    file.close();
  }
  else
  {
    cout << "Unable to open file" << endl;
  }

  Bool_t isMagUp = filename.find("MagUp") != string::npos;

  gStyle->SetOptStat(0);
  
  // Set branch address
  // Ktautau
  Double_t Kp_L0Calo_HCAL_realET, taup_pi1_L0Calo_HCAL_realET, taup_pi2_L0Calo_HCAL_realET, taup_pi3_L0Calo_HCAL_realET, taum_pi1_L0Calo_HCAL_realET, taum_pi2_L0Calo_HCAL_realET, taum_pi3_L0Calo_HCAL_realET;
  Int_t Kp_L0Calo_HCAL_region, taup_pi1_L0Calo_HCAL_region, taup_pi2_L0Calo_HCAL_region, taup_pi3_L0Calo_HCAL_region, taum_pi1_L0Calo_HCAL_region, taum_pi2_L0Calo_HCAL_region, taum_pi3_L0Calo_HCAL_region;
  Double_t Kp_L0Calo_HCAL_xProjection, taup_pi1_L0Calo_HCAL_xProjection, taup_pi2_L0Calo_HCAL_xProjection, taup_pi3_L0Calo_HCAL_xProjection, taum_pi1_L0Calo_HCAL_xProjection, taum_pi2_L0Calo_HCAL_xProjection, taum_pi3_L0Calo_HCAL_xProjection;
  Double_t Kp_L0Calo_HCAL_yProjection, taup_pi1_L0Calo_HCAL_yProjection, taup_pi2_L0Calo_HCAL_yProjection, taup_pi3_L0Calo_HCAL_yProjection, taum_pi1_L0Calo_HCAL_yProjection, taum_pi2_L0Calo_HCAL_yProjection, taum_pi3_L0Calo_HCAL_yProjection;
  // DDs
  Double_t D0bar_K_L0Calo_HCAL_realET, D0bar_pi_L0Calo_HCAL_realET, Dsp_K1_L0Calo_HCAL_realET, Dsp_K2_L0Calo_HCAL_realET, Dsp_pi_L0Calo_HCAL_realET;
  Int_t D0bar_K_L0Calo_HCAL_region, D0bar_pi_L0Calo_HCAL_region, Dsp_K1_L0Calo_HCAL_region, Dsp_K2_L0Calo_HCAL_region, Dsp_pi_L0Calo_HCAL_region;
  Double_t D0bar_K_L0Calo_HCAL_xProjection, D0bar_pi_L0Calo_HCAL_xProjection, Dsp_K1_L0Calo_HCAL_xProjection, Dsp_K2_L0Calo_HCAL_xProjection, Dsp_pi_L0Calo_HCAL_xProjection;
  Double_t D0bar_K_L0Calo_HCAL_yProjection, D0bar_pi_L0Calo_HCAL_yProjection, Dsp_K1_L0Calo_HCAL_yProjection, Dsp_K2_L0Calo_HCAL_yProjection, Dsp_pi_L0Calo_HCAL_yProjection;

  Int_t Kp_TRUEID, taup_pi1_TRUEID, taup_pi2_TRUEID, taup_pi3_TRUEID, taum_pi1_TRUEID, taum_pi2_TRUEID, taum_pi3_TRUEID, taup_TRUEID, taum_TRUEID, Bp_TRUEID;
  Int_t D0bar_K_TRUEID, D0bar_pi_TRUEID, Dsp_K1_TRUEID, Dsp_K2_TRUEID, Dsp_pi_TRUEID, D0bar_TRUEID, Dsp_TRUEID;
  Int_t mc_component;

  Bool_t Bp_L0HadronDecision_TOS, Bp_L0HadronDecision_TIS, Bp_L0MuonDecision_TIS, Bp_L0ElectronDecision_TIS, Bp_L0PhotonDecision_TIS, Bp_Hlt1TrackMVADecision_TOS, Bp_Hlt1TwoTrackMVADecision_TOS, Bp_Hlt2Topo2BodyDecision_TOS, Bp_Hlt2Topo3BodyDecision_TOS, Bp_Hlt2Topo4BodyDecision_TOS;

  if(species == 1)
  {
    t->SetBranchAddress("Kp_TRUEID", &Kp_TRUEID);
    t->SetBranchAddress("taup_pi1_TRUEID", &taup_pi1_TRUEID);
    t->SetBranchAddress("taup_pi2_TRUEID", &taup_pi2_TRUEID);
    t->SetBranchAddress("taup_pi3_TRUEID", &taup_pi3_TRUEID);
    t->SetBranchAddress("taum_pi1_TRUEID", &taum_pi1_TRUEID);
    t->SetBranchAddress("taum_pi2_TRUEID", &taum_pi2_TRUEID);
    t->SetBranchAddress("taum_pi3_TRUEID", &taum_pi3_TRUEID);
    t->SetBranchAddress("Kp_L0Calo_HCAL_realET", &Kp_L0Calo_HCAL_realET);
    t->SetBranchAddress("taup_pi1_L0Calo_HCAL_realET", &taup_pi1_L0Calo_HCAL_realET);
    t->SetBranchAddress("taup_pi2_L0Calo_HCAL_realET", &taup_pi2_L0Calo_HCAL_realET);
    t->SetBranchAddress("taup_pi3_L0Calo_HCAL_realET", &taup_pi3_L0Calo_HCAL_realET);
    t->SetBranchAddress("taum_pi1_L0Calo_HCAL_realET", &taum_pi1_L0Calo_HCAL_realET);
    t->SetBranchAddress("taum_pi2_L0Calo_HCAL_realET", &taum_pi2_L0Calo_HCAL_realET);
    t->SetBranchAddress("taum_pi3_L0Calo_HCAL_realET", &taum_pi3_L0Calo_HCAL_realET);
    t->SetBranchAddress("Kp_L0Calo_HCAL_region", &Kp_L0Calo_HCAL_region);
    t->SetBranchAddress("taup_pi1_L0Calo_HCAL_region", &taup_pi1_L0Calo_HCAL_region);
    t->SetBranchAddress("taup_pi2_L0Calo_HCAL_region", &taup_pi2_L0Calo_HCAL_region);
    t->SetBranchAddress("taup_pi3_L0Calo_HCAL_region", &taup_pi3_L0Calo_HCAL_region);
    t->SetBranchAddress("taum_pi1_L0Calo_HCAL_region", &taum_pi1_L0Calo_HCAL_region);
    t->SetBranchAddress("taum_pi2_L0Calo_HCAL_region", &taum_pi2_L0Calo_HCAL_region);
    t->SetBranchAddress("taum_pi3_L0Calo_HCAL_region", &taum_pi3_L0Calo_HCAL_region);
    t->SetBranchAddress("Kp_L0Calo_HCAL_xProjection", &Kp_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("taup_pi1_L0Calo_HCAL_xProjection", &taup_pi1_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("taup_pi2_L0Calo_HCAL_xProjection", &taup_pi2_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("taup_pi3_L0Calo_HCAL_xProjection", &taup_pi3_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("taum_pi1_L0Calo_HCAL_xProjection", &taum_pi1_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("taum_pi2_L0Calo_HCAL_xProjection", &taum_pi2_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("taum_pi3_L0Calo_HCAL_xProjection", &taum_pi3_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("Kp_L0Calo_HCAL_yProjection", &Kp_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("taup_pi1_L0Calo_HCAL_yProjection", &taup_pi1_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("taup_pi2_L0Calo_HCAL_yProjection", &taup_pi2_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("taup_pi3_L0Calo_HCAL_yProjection", &taup_pi3_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("taum_pi1_L0Calo_HCAL_yProjection", &taum_pi1_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("taum_pi2_L0Calo_HCAL_yProjection", &taum_pi2_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("taum_pi3_L0Calo_HCAL_yProjection", &taum_pi3_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("Kp_TRUEID", &Kp_TRUEID);
    t->SetBranchAddress("taup_pi1_TRUEID", &taup_pi1_TRUEID);
    t->SetBranchAddress("taup_pi2_TRUEID", &taup_pi2_TRUEID);
    t->SetBranchAddress("taup_pi3_TRUEID", &taup_pi3_TRUEID);
    t->SetBranchAddress("taum_pi1_TRUEID", &taum_pi1_TRUEID);
    t->SetBranchAddress("taum_pi2_TRUEID", &taum_pi2_TRUEID);
    t->SetBranchAddress("taum_pi3_TRUEID", &taum_pi3_TRUEID);
    t->SetBranchAddress("taup_TRUEID", &taup_TRUEID);
    t->SetBranchAddress("taum_TRUEID", &taum_TRUEID);
    t->SetBranchAddress("component", &mc_component);
  }
  else
  {
    t->SetBranchAddress("D0bar_K_TRUEID", &D0bar_K_TRUEID);
    t->SetBranchAddress("D0bar_pi_TRUEID", &D0bar_pi_TRUEID);
    t->SetBranchAddress("Dsp_K1_TRUEID", &Dsp_K1_TRUEID);
    t->SetBranchAddress("Dsp_K2_TRUEID", &Dsp_K2_TRUEID);
    t->SetBranchAddress("Dsp_pi_TRUEID", &Dsp_pi_TRUEID);
    t->SetBranchAddress("D0bar_K_L0Calo_HCAL_realET", &D0bar_K_L0Calo_HCAL_realET);
    t->SetBranchAddress("D0bar_pi_L0Calo_HCAL_realET", &D0bar_pi_L0Calo_HCAL_realET);
    t->SetBranchAddress("Dsp_K1_L0Calo_HCAL_realET", &Dsp_K1_L0Calo_HCAL_realET);
    t->SetBranchAddress("Dsp_K2_L0Calo_HCAL_realET", &Dsp_K2_L0Calo_HCAL_realET);
    t->SetBranchAddress("Dsp_pi_L0Calo_HCAL_realET", &Dsp_pi_L0Calo_HCAL_realET);
    t->SetBranchAddress("D0bar_K_L0Calo_HCAL_region", &D0bar_K_L0Calo_HCAL_region);
    t->SetBranchAddress("D0bar_pi_L0Calo_HCAL_region", &D0bar_pi_L0Calo_HCAL_region);
    t->SetBranchAddress("Dsp_K1_L0Calo_HCAL_region", &Dsp_K1_L0Calo_HCAL_region);
    t->SetBranchAddress("Dsp_K2_L0Calo_HCAL_region", &Dsp_K2_L0Calo_HCAL_region);
    t->SetBranchAddress("Dsp_pi_L0Calo_HCAL_region", &Dsp_pi_L0Calo_HCAL_region);
    t->SetBranchAddress("D0bar_K_L0Calo_HCAL_xProjection", &D0bar_K_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("D0bar_pi_L0Calo_HCAL_xProjection", &D0bar_pi_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("Dsp_K1_L0Calo_HCAL_xProjection", &Dsp_K1_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("Dsp_K2_L0Calo_HCAL_xProjection", &Dsp_K2_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("Dsp_pi_L0Calo_HCAL_xProjection", &Dsp_pi_L0Calo_HCAL_xProjection);
    t->SetBranchAddress("D0bar_K_L0Calo_HCAL_yProjection", &D0bar_K_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("D0bar_pi_L0Calo_HCAL_yProjection", &D0bar_pi_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("Dsp_K1_L0Calo_HCAL_yProjection", &Dsp_K1_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("Dsp_K2_L0Calo_HCAL_yProjection", &Dsp_K2_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("Dsp_pi_L0Calo_HCAL_yProjection", &Dsp_pi_L0Calo_HCAL_yProjection);
    t->SetBranchAddress("D0bar_K_TRUEID", &D0bar_K_TRUEID);
    t->SetBranchAddress("D0bar_pi_TRUEID", &D0bar_pi_TRUEID);
    t->SetBranchAddress("Dsp_K1_TRUEID", &Dsp_K1_TRUEID);
    t->SetBranchAddress("Dsp_K2_TRUEID", &Dsp_K2_TRUEID);
    t->SetBranchAddress("Dsp_pi_TRUEID", &Dsp_pi_TRUEID);
    t->SetBranchAddress("D0bar_TRUEID", &D0bar_TRUEID);
    t->SetBranchAddress("Dsp_TRUEID", &Dsp_TRUEID);
  }
  t->SetBranchAddress("Bp_TRUEID", &Bp_TRUEID);
  t->SetBranchAddress("Bp_L0HadronDecision_TOS", &Bp_L0HadronDecision_TOS);
  t->SetBranchAddress("Bp_L0HadronDecision_TIS", &Bp_L0HadronDecision_TIS);
  t->SetBranchAddress("Bp_L0MuonDecision_TIS", &Bp_L0MuonDecision_TIS);
  t->SetBranchAddress("Bp_L0ElectronDecision_TIS", &Bp_L0ElectronDecision_TIS);
  t->SetBranchAddress("Bp_L0PhotonDecision_TIS", &Bp_L0PhotonDecision_TIS);
  t->SetBranchAddress("Bp_Hlt1TrackMVADecision_TOS", &Bp_Hlt1TrackMVADecision_TOS);
  t->SetBranchAddress("Bp_Hlt1TwoTrackMVADecision_TOS", &Bp_Hlt1TwoTrackMVADecision_TOS);
  t->SetBranchAddress("Bp_Hlt2Topo2BodyDecision_TOS", &Bp_Hlt2Topo2BodyDecision_TOS);
  t->SetBranchAddress("Bp_Hlt2Topo3BodyDecision_TOS", &Bp_Hlt2Topo3BodyDecision_TOS);
  t->SetBranchAddress("Bp_Hlt2Topo4BodyDecision_TOS", &Bp_Hlt2Topo4BodyDecision_TOS);

  std::map<TString, map < TString, std::vector<double> > > mytable;
  cout<<"reading table"<<endl;
  ReadTables(year, isMagUp);// now we can use the run2 method for run1
  cout<<"Done"<<endl;

  Int_t nDau;
  if(species == 1)
  {
    nDau = 7;
  }
  else
  {
    nDau = 5;
  }

  std::vector<Int_t> pid, region;
  std::vector<Double_t> x, y, ET;
  std::vector<Double_t> eff_tot, effUp_tot, effDn_tot;
  std::vector<Int_t> comp_tot;
  std::vector<Bool_t> HadronTOS, HadronTIS, ElectronTIS, PhotonTIS, MuonTIS, TrackMVATOS, TwoTrackMVATOS, Topo2TOS, Topo3TOS, Topo4TOS;

  Int_t nentries = t->GetEntries();
  for(int evt = 0; evt < nentries; evt++)
  {
    if(evt%10000==0)std::cout<<"entry "<<evt<<"/"<<nentries<<endl;
    t->GetEntry(evt);

    Bool_t truthMatch;
    if(species == 1)
    {
      truthMatch = (abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521);
    }
    else
    {
      truthMatch = (abs(D0bar_K_TRUEID) == 321) && (abs(Dsp_K1_TRUEID) == 321) && (abs(Dsp_K2_TRUEID) == 321) && (abs(D0bar_pi_TRUEID) == 211) && (abs(Dsp_pi_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(Dsp_TRUEID) == 431) && (abs(Bp_TRUEID) == 521);
    }

    if(truthMatch)
    {
      Int_t Ktautau_pid[7] = { Kp_TRUEID, taup_pi1_TRUEID, taup_pi2_TRUEID, taup_pi3_TRUEID, taum_pi1_TRUEID, taum_pi2_TRUEID, taum_pi3_TRUEID };
      Int_t Ktautau_region[7] = { Kp_L0Calo_HCAL_region, taup_pi1_L0Calo_HCAL_region, taup_pi2_L0Calo_HCAL_region, taup_pi3_L0Calo_HCAL_region, taum_pi1_L0Calo_HCAL_region, taum_pi2_L0Calo_HCAL_region, taum_pi3_L0Calo_HCAL_region };
      Double_t Ktautau_x[7] = { Kp_L0Calo_HCAL_xProjection, taup_pi1_L0Calo_HCAL_xProjection, taup_pi2_L0Calo_HCAL_xProjection, taup_pi3_L0Calo_HCAL_xProjection, taum_pi1_L0Calo_HCAL_xProjection, taum_pi2_L0Calo_HCAL_xProjection, taum_pi3_L0Calo_HCAL_xProjection };
      Double_t Ktautau_y[7] = { Kp_L0Calo_HCAL_yProjection, taup_pi1_L0Calo_HCAL_yProjection, taup_pi2_L0Calo_HCAL_yProjection, taup_pi3_L0Calo_HCAL_yProjection, taum_pi1_L0Calo_HCAL_yProjection, taum_pi2_L0Calo_HCAL_yProjection, taum_pi3_L0Calo_HCAL_yProjection };
      Double_t Ktautau_ET[7] = { Kp_L0Calo_HCAL_realET, taup_pi1_L0Calo_HCAL_realET, taup_pi2_L0Calo_HCAL_realET, taup_pi3_L0Calo_HCAL_realET, taum_pi1_L0Calo_HCAL_realET, taum_pi2_L0Calo_HCAL_realET, taum_pi3_L0Calo_HCAL_realET };

      Int_t DDs_pid[5] = { D0bar_K_TRUEID, D0bar_pi_TRUEID, Dsp_K1_TRUEID, Dsp_K2_TRUEID, Dsp_pi_TRUEID };
      Int_t DDs_region[5] = { D0bar_K_L0Calo_HCAL_region, D0bar_pi_L0Calo_HCAL_region, Dsp_K1_L0Calo_HCAL_region, Dsp_K2_L0Calo_HCAL_region, Dsp_pi_L0Calo_HCAL_region };
      Double_t DDs_x[5] = { D0bar_K_L0Calo_HCAL_xProjection, D0bar_pi_L0Calo_HCAL_xProjection, Dsp_K1_L0Calo_HCAL_xProjection, Dsp_K2_L0Calo_HCAL_xProjection, Dsp_pi_L0Calo_HCAL_xProjection };
      Double_t DDs_y[5] = { D0bar_K_L0Calo_HCAL_yProjection, D0bar_pi_L0Calo_HCAL_yProjection, Dsp_K1_L0Calo_HCAL_yProjection, Dsp_K2_L0Calo_HCAL_yProjection, Dsp_pi_L0Calo_HCAL_yProjection };
      Double_t DDs_ET[5] = { D0bar_K_L0Calo_HCAL_realET, D0bar_pi_L0Calo_HCAL_realET, Dsp_K1_L0Calo_HCAL_realET, Dsp_K2_L0Calo_HCAL_realET, Dsp_pi_L0Calo_HCAL_realET };

      for(int i = 0; i < nDau; i++)
      {
        if(species == 1)
        {
          pid.push_back(Ktautau_pid[i]);
          region.push_back(Ktautau_region[i]);
          x.push_back(Ktautau_x[i]);
          y.push_back(Ktautau_y[i]);
          ET.push_back(Ktautau_ET[i]);
        }
        else
        {
          pid.push_back(DDs_pid[i]);
          region.push_back(DDs_region[i]);
          x.push_back(DDs_x[i]);
          y.push_back(DDs_y[i]);
          ET.push_back(DDs_ET[i]);
        }
      }

      // Find particle ID and area
      char particle_ID[20][10];
      TString _particle_ID[20];
      char area[20][10];
      for(int i = 0; i < nDau; i++)
      {
        if(pid[i] == 211)
        {
          snprintf( particle_ID[i] , sizeof particle_ID[i] , "%s" , "piplus" );
          _particle_ID[ i ] = "piplus";
        }
        else if(pid[i] == -211)
        {
          snprintf( particle_ID[i] , sizeof particle_ID[i] , "%s" , "piminus" ) ;
          _particle_ID[ i ] = "piminus";
        }
        else if(pid[i] == 321)
        {
          snprintf( particle_ID[i] , sizeof particle_ID[i] , "%s" , "Kplus" ) ;
          _particle_ID[ i ] = "Kplus";
        }
        else if(pid[i] == -321)
        {
          snprintf( particle_ID[i] , sizeof particle_ID[i] , "%s" , "Kminus" ) ;
          _particle_ID[ i ] = "Kminus";
        }
        else
        {
          snprintf( particle_ID[i] , sizeof particle_ID[i] , "%s" , "" ) ;
          _particle_ID[i] = "";
        }
        if(region[i] == 0) snprintf(area[i], sizeof area[i] , "%s" , "outer");
        else if(region[i] == 1) snprintf(area[i] , sizeof area[i] , "%s" , "inner");
        else snprintf(area[i], sizeof area[i], "%s" , "");
      }
      Double_t fraction_zeros[100] ;
      Double_t energy_overlaps[100][100] ;
      Double_t E0[100] ;
      Double_t betarel[100] ;
      Double_t eff[nDau+1] ;
      Double_t eff_Up[nDau+1];
      Double_t eff_Dn[nDau+1];

      // Fill fraction of tracks that do not reach HCAL, overlap tables, E0 tables and betarel tables
      for(int i = 0; i < nDau; i++)
      {
        eff[i] = 0.;
        eff_Up[i] = 0.;
        eff_Dn[i] = 0;
        fraction_zeros[i] = fraction_zero(ET[i], region[i]);
        for(int j = 0; j < nDau; j++)
        {
          if(i == j)
          {
            energy_overlaps[i][j] = 0;
          }
          else
          {
            energy_overlaps[i][j] = energy_overlap(ET[j], x[j]-x[i], y[j] - y[i], region[i], region[j]);
          }
        }
        if(strcmp(area[ i ], "") == 0) continue ; // test if area[i] is empty, if yes (e.g. the particle doesn't go in the HCAL) start next iteration in the for loop (on i here)
      
        char histoname[100];
        Int_t bin_number;
        snprintf(histoname, sizeof histoname , "E0_%s" , area[i]); // this comes from the histograms read with the function 'ReadtaTables'
        bin_number = ((TH2D*)gDirectory->Get(histoname))->FindBin(x[i], y[i]);
        E0[i] = ((TH1D*)gDirectory->Get(histoname))->GetBinContent(bin_number);
        snprintf(histoname, sizeof histoname, "betarel_%s", area[i]) ;
        bin_number = ((TH2D*)gDirectory->Get(histoname))->FindBin(x[i], y[i]);
        betarel[i] = ((TH1D*)gDirectory->Get(histoname))->GetBinContent(bin_number);
      }

      // for loop system with bitset test in order to consider each overlap combination, we have 4 tracks, for one given track (i) we have to question the overlap with:
      // only the track itself (=0 cf energy_overlaps def) (1 config) , each of the 3 other tracks individualy (3 configs), each combination of 2 of the 3 other tracks (3 config), the 3 other tracks together (1 config)
      // the fraction_zeros is also computed because used in the final computation
      // all of this in order to obtain the efficiency attached to each track (containing overlap, zero frac and betarel)

      char histoname[100];
      Int_t bin_number;
      for(int i = 0; i < nDau; i++)
      {
        eff[i] = 0. ;
        if ( strcmp( area[i] , "" ) == 0 ) continue ;
        if ( strcmp( particle_ID[i] , "" ) == 0 ) continue ;

        for ( int j = 0 ; j < ( 1 << nDau ) ; ++j )// 1<<ndau=2^x -> loop over bitset to condider the 16 possible combinations of tracks deposite/nodeposite into the HCAL
        {
          std::bitset< 32 > bj(j);
          if ( bj.test(i) ) // continue only if the reference track is turning on ((depositing into HCAL)) 
          {
            Double_t fraction_zeros_temp[ 100 ] ;
            Double_t energies_temp[ 100 ] ;
            for ( int k = 0 ; k < nDau ; ++k )//loop over the track considered to apply the overlap (including the reference track itself)
            {
              if ( bj.test( k ) ) 
              {//define operations if the k track is turning on or not (depositing into HCAL)
                fraction_zeros_temp[k] = 1. - fraction_zeros[k] ;
                energies_temp[k] = energy_overlaps[i][k] ;
              }
              else 
              {
                fraction_zeros_temp[k] = fraction_zeros[k] ;
                energies_temp[k] = 0. ;
              }
            }

            fraction_zeros_temp[i] = 1. ; // fix frac 0 of the i track itself to 1 (neutral element of the product)
            energies_temp[i] = ET[ i ]  ; // fix energy overlap of the i track itself to the transverse energy of the i track itself
            snprintf( histoname , sizeof histoname , "eff_%s_%s" , area[i] , particle_ID[i] ) ;//take the efficiency table in function of ET as histoname
            bin_number = ((TH1D*)gDirectory->Get( histoname ))->FindBin( betarel[ i ]*std::accumulate( energies_temp , energies_temp + nDau , 0 ) + 4*E0[i] );// find efficiency bin (ET) corresponding to ((sum_overlaps * betarel) +E0) for i track, E0 has to be considered because efficiency maps are build corrected from E0 maps to be general but E0 is there in the HCAL
            if (bin_number > ((TH1D*)gDirectory->Get( histoname ))->GetNbinsX()) bin_number = bin_number - 1 ;
            //std::accumulate( fraction_zeros_temp , fraction_zeros_temp + nDau , 1. , std::multiplies< double >() ) is equal to 1*product_{k,0}^{nDau(=4)}(fraction_zeros_temp[k])
            eff[i] += std::accumulate( fraction_zeros_temp , fraction_zeros_temp + nDau , 1. , std::multiplies< double >() )*((TH1D*)gDirectory->Get( histoname ))->GetBinContent( bin_number ) ;//final efficiency on the track i is incremented eff+=(product_empty * eff[efficiency bin previously determined])
            // determine eff envelop (only last loop bin_number is keeped -> the knowledge of all loops is in eff[i], bin number must also contain memories or maybe the ET deviation is considered small ?)
            eff_Up[i] = eff[i] + ((TH1D*)gDirectory->Get( histoname ))->GetBinError( bin_number );
            eff_Dn[i] = eff[i] - ((TH1D*)gDirectory->Get( histoname ))->GetBinError( bin_number );

          }
        }
      }

      // now combine all efficiencies to get the mother one (total efficiency TOS and noTOS of the event from combination from eff of each track)
      eff[nDau] = 0.; // eff at nDau = eff from all tracks together
      eff_Up[nDau] = 0.;
      eff_Dn[nDau] = 0.;

      if(species == 1)
      {
        for ( int j = 1 ; j < ( 1 << nDau ) ; ++j )//similar to previous loop system in order to consider tot efficiency = product of all combination
        {
          std::bitset< 7 > bj(j);

          Double_t eff_temp[100];
          Double_t eff_Up_temp[100];
          Double_t eff_Dn_temp[100];
          for ( int k = 0 ; k < nDau ; ++k )
          {
            if ( bj.test(k) )
            {
              eff_temp[k] = eff[k];
              eff_Up_temp[k] = eff_Up[k];
              eff_Dn_temp[k] = eff_Dn[k];
            }
            else
            {
              eff_temp[k] = 1. - eff[k];
              eff_Up_temp[k] = 1. - eff_Up[k];
              eff_Dn_temp[k] = 1. - eff_Dn[k];
            }
          }
          //cout<< bj <<endl;
          //cout<< std::accumulate( eff_temp , eff_temp + nDau , 1. , std::multiplies< double >() ) <<endl;
          eff[nDau] += std::accumulate( eff_temp , eff_temp + nDau , 1. , std::multiplies< double >() ) ;
          eff_Up[nDau] += std::accumulate( eff_Up_temp , eff_Up_temp + nDau , 1. , std::multiplies<double>() );
          eff_Dn[nDau] += std::accumulate( eff_Dn_temp , eff_Dn_temp + nDau , 1. , std::multiplies<double>() ); 
        }
      }
      else
      {
        for ( int j = 1 ; j < ( 1 << nDau ) ; ++j )//similar to previous loop system in order to consider tot efficiency = product of all combination
        {
          std::bitset< 5 > bj(j);

          Double_t eff_temp[100];
          Double_t eff_Up_temp[100];
          Double_t eff_Dn_temp[100];
          for ( int k = 0 ; k < nDau ; ++k )
          {
            if ( bj.test(k) )
            {
              eff_temp[k] = eff[k];
              eff_Up_temp[k] = eff_Up[k];
              eff_Dn_temp[k] = eff_Dn[k];
            }
            else
            {
              eff_temp[k] = 1. - eff[k];
              eff_Up_temp[k] = 1. - eff_Up[k];
              eff_Dn_temp[k] = 1. - eff_Dn[k];
            }
          }
          //cout<< bj <<endl;
          //cout<< std::accumulate( eff_temp , eff_temp + nDau , 1. , std::multiplies< double >() ) <<endl;
          eff[nDau] += std::accumulate( eff_temp , eff_temp + nDau , 1. , std::multiplies< double >() ) ;
          eff_Up[nDau] += std::accumulate( eff_Up_temp , eff_Up_temp + nDau , 1. , std::multiplies<double>() );
          eff_Dn[nDau] += std::accumulate( eff_Dn_temp , eff_Dn_temp + nDau , 1. , std::multiplies<double>() ); 
        }
      }
      
      // eff[nDau] is the TOS L0Hadron corrected efficiency 
      // cout << "eff TOS event " << evt << " eff=" << eff[nDau] << endl;

      eff_tot.push_back(eff[nDau]);
      effUp_tot.push_back(eff_Up[nDau]);
      effDn_tot.push_back(eff_Dn[nDau]);
      if(species == 1)
      {
        comp_tot.push_back(mc_component);
      }

      HadronTOS.push_back(Bp_L0HadronDecision_TOS);
      HadronTIS.push_back(Bp_L0HadronDecision_TIS);
      ElectronTIS.push_back(Bp_L0ElectronDecision_TIS);
      MuonTIS.push_back(Bp_L0MuonDecision_TIS);
      PhotonTIS.push_back(Bp_L0PhotonDecision_TIS);
      TrackMVATOS.push_back(Bp_Hlt1TrackMVADecision_TOS);
      TwoTrackMVATOS.push_back(Bp_Hlt1TwoTrackMVADecision_TOS);
      Topo2TOS.push_back(Bp_Hlt2Topo2BodyDecision_TOS);
      Topo3TOS.push_back(Bp_Hlt2Topo3BodyDecision_TOS);
      Topo4TOS.push_back(Bp_Hlt2Topo4BodyDecision_TOS);

    }

  }

  Double_t L0Hadron_TOS_eff, L0Hadron_TOS_eff_errUp, L0Hadron_TOS_eff_errDown;
  Int_t mc_comp;
  Bool_t H_TOS, H_TIS, E_TIS, P_TIS, M_TIS, HLT1_one_TOS, HLT1_two_TOS, HLT2_2_TOS, HLT2_3_TOS, HLT2_4_TOS;
  TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/TriggerCorrection/201%i/Species_%i/%i.root",year,species,line), "RECREATE");
  TTree* tout = new TTree("DecayTree", "DecayTree");
  tout->Branch("L0Hadron_TOS_eff_corr", &L0Hadron_TOS_eff);
  tout->Branch("L0Hadron_TOS_eff_corr_errUp", &L0Hadron_TOS_eff_errUp);
  tout->Branch("L0Hadron_TOS_eff_corr_errDown", &L0Hadron_TOS_eff_errDown);
  if(species == 1)
  {
    tout->Branch("component", &mc_comp);
  }
  tout->Branch("Bp_L0HadronDecision_TOS", &H_TOS);
  tout->Branch("Bp_L0HadronDecision_TIS", &H_TIS);
  tout->Branch("Bp_L0MuonDecision_TIS", &M_TIS);
  tout->Branch("Bp_L0ElectronDecision_TIS", &E_TIS);
  tout->Branch("Bp_L0PhotonDecision_TIS", &P_TIS);
  tout->Branch("Bp_Hlt1TrackMVADecision_TOS", &HLT1_one_TOS);
  tout->Branch("Bp_Hlt1TwoTrackMVADecision_TOS", &HLT1_two_TOS);
  tout->Branch("Bp_Hlt2Topo2BodyDecision_TOS", &HLT2_2_TOS);
  tout->Branch("Bp_Hlt2Topo3BodyDecision_TOS", &HLT2_3_TOS);
  tout->Branch("Bp_Hlt2Topo4BodyDecision_TOS", &HLT2_4_TOS);

  cout << "Filling TTree" << endl;
  Int_t N = eff_tot.size();
  for(int i = 0; i < N; i++)
  {
    L0Hadron_TOS_eff = eff_tot[i];
    L0Hadron_TOS_eff_errUp = effUp_tot[i];
    L0Hadron_TOS_eff_errDown = effDn_tot[i];
    if(species == 1)
    {
      mc_comp = comp_tot[i];
    }
    H_TOS = HadronTOS[i];
    H_TIS = HadronTIS[i];
    E_TIS = ElectronTIS[i];
    P_TIS = PhotonTIS[i];
    M_TIS = MuonTIS[i];
    HLT1_one_TOS = TrackMVATOS[i];
    HLT1_two_TOS = TwoTrackMVATOS[i];
    HLT2_2_TOS = Topo2TOS[i];
    HLT2_3_TOS = Topo3TOS[i];
    HLT2_4_TOS = Topo4TOS[i];

    tout->Fill();
  }

  fout->cd();
  tout->Write();
  fout->Close();
}

Double_t fraction_zero( Double_t E , Int_t region)
{
  Double_t A, B ;

  if ( region == 1 ) {
    A = 0.075924 ;
    B = 0.000651308 ;}
  else if ( region == 0 )
  {
    A = 0.233658;
    B = 0.0007411413;
  }
  else
    return 1. ;

  return ( A * TMath::Exp( -B*E ) ) ;
}

Double_t energy_overlap( Double_t E , Double_t X , Double_t Y , Int_t region1 , Int_t region2 )
{
  // region = 1 pour inner, 0 pour outer
  Double_t sigma1, sigma2, c , a1 ;

  if ( region2 == 1 )
  {
    sigma1 = 64. ;
    sigma2 = 202. ;
    a1 = 333. / ( 333. + 103. ) ;
  } else if ( region2 == 0 )
  {
    sigma1 = 100. ;
    sigma2 = 240. ;
    a1 = 238. / (238.+62.) ;
  } else return 0;

  if ( region1 == 1 ) c = 131. ;
  else if ( region1 == 0 ) c = 262. ;
  else return 0. ;

  Double_t A1, B1, C1, D1, A2, B2, C2, D2 ;

  A1 = TMath::Erf( ( X + c ) / ( TMath::Sqrt2() * sigma1 ) ) ;
  B1 = TMath::Erf( ( c - X ) / ( TMath::Sqrt2() * sigma1 ) ) ;
  C1 = TMath::Erf( ( Y + c ) / ( TMath::Sqrt2() * sigma1 ) ) ;
  D1 = TMath::Erf( ( c - Y ) / ( TMath::Sqrt2() * sigma1 ) ) ;
  A2 = TMath::Erf( ( X + c ) / ( TMath::Sqrt2() * sigma2 ) ) ;
  B2 = TMath::Erf( ( c - X ) / ( TMath::Sqrt2() * sigma2 ) ) ;
  C2 = TMath::Erf( ( Y + c ) / ( TMath::Sqrt2() * sigma2 ) ) ;
  D2 = TMath::Erf( ( c - Y ) / ( TMath::Sqrt2() * sigma2 ) ) ;

  return E * ( a1 * ( A1 * C1 + B1 * C1 + A1 * D1 + B1 * D1 ) + ( 1- a1 ) * ( A2 * C2 + B2 * C2 + A2 * D2 + B2 * D2 ) ) / 4. ;
}

void ReadFile( FILE * fp , Double_t * array , Int_t * nbins )
{
  int i = 0 ;
  while (!feof(fp))// feof say if we are at the end of the file
  {
    if (fscanf(fp,"%lf",&array[i]) != 1)//fscanf allow to extract data from a file and put them into a variable (here acces the input file, extract one double number and put it into an array)
    {
      continue;
    }
    i++;
  }
  if ( nbins ) *nbins = i ;
}

void ReadTables(int year, Bool_t isMagUp)
{
  FILE *fptr ;

  Double_t eff_bins[1000];
  Int_t nbins_eff ;
  fptr = fopen( "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/eff_inner_HCAL_bins.txt" , "r" ) ;
  ReadFile( fptr , eff_bins , &nbins_eff ) ;
  fclose( fptr ) ;

  cout<<"Read file inner HCAL bins Done"<<endl;

  TString used_year = Form("201%i",year);
  cout << used_year << endl;
  TString Magnet;
  if(isMagUp)
  {
    Magnet = "MU";
  }
  else
  {
    Magnet = "MD";
  }

  // Tables of efficiency
  char area[2][10] = { "inner" , "outer" } ;
  char species[4][10] = { "Kminus" , "Kplus" , "piplus" , "piminus" } ;
  //  char magnet[3] = Magnet ;
  //  char year[5] = Year ;
  char filename[1000];

  Double_t eff[ 1000 ] ;
  Double_t eff_err[ 1000 ] ;
  int i, j, k ;
  for ( i = 0 ; i < 2 ; ++i ) {
    for ( j = 0 ; j < 4 ; ++j ) {
      snprintf( filename , sizeof filename , "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/eff_%s_%s_%s_%s.txt" , area[i] , species[j] , used_year.Data() , Magnet.Data() ) ;//function to build path to file at each iteration
      fptr = fopen( filename , "r" ) ;
      ReadFile( fptr , eff , 0 ) ;
      fclose( fptr ) ;
      //cout<<"Read file eff "<<area[i]<<" "<< species[j]<<" "<< used_year.Data()<<" "<< Magnet.Data() <<" Done"<<endl;
      
      snprintf( filename , sizeof filename , "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/eff_%s_%s_%s_%s_err.txt" , area[i] , species[j] , used_year.Data() , Magnet.Data() ) ;
      fptr = fopen( filename , "r" ) ;
      ReadFile( fptr , eff_err , 0 ) ;
      fclose( fptr ) ;
      snprintf( filename , sizeof filename , "eff_%s_%s" , area[i] , species[j] ) ;
      TH1D* h = new TH1D( filename , filename , nbins_eff - 1 , eff_bins ) ;
      h->SetCanExtend(false) ;
      for ( k = 0 ; k < nbins_eff-1 ; ++k )
      {
        h->SetBinContent( k+1 , eff[ k ] ) ;
        h->SetBinError( k+1 , eff_err[ k ] ) ;
      }
    }
  }

  // Tables of E0
  Double_t E0[ 10000 ] ;
  Int_t nbins_xy ;
  snprintf( filename , sizeof filename , "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/E0_bb_inner_%s_%s.txt" , used_year.Data(), Magnet.Data() ) ;
  fptr = fopen( filename , "r" ) ;
  ReadFile( fptr , E0 , &nbins_xy ) ;
  fclose( fptr ) ;
  //cout<<"Read file E0 bb inner"<< used_year.Data()<<" "<< Magnet.Data() <<" Done"<<endl;

  Double_t bins_x[ 1000 ] , bins_y[ 1000 ] ;
  Int_t nbins_hcal_x , nbins_hcal_y ;
  fptr = fopen( "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/2d_hcal_inner_bins_x.txt" , "r" ) ;
  ReadFile( fptr , bins_x , &nbins_hcal_x ) ;
  fclose( fptr ) ;
  fptr = fopen( "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/2d_hcal_inner_bins_y.txt" , "r" ) ;
  ReadFile( fptr , bins_y , &nbins_hcal_y ) ;
  fclose( fptr ) ;

  snprintf( filename , sizeof filename , "E0_%s" , "inner" ) ;
  Double_t E0_2D[nbins_hcal_x-1][nbins_hcal_y-1] ;
  memcpy( E0_2D , E0 , nbins_xy * sizeof( Double_t ) ) ;//function to pass from flat E0 to 2D 
  TH2D* hh = new TH2D( filename , filename , nbins_hcal_x - 1 , bins_x , nbins_hcal_y - 1 , bins_y ) ;
  for ( i = 0 ; i < nbins_hcal_x - 1 ; ++i )
    for ( j = 0 ; j < nbins_hcal_y - 1 ; ++j )
      hh->SetBinContent( i+1 , j+1 , E0_2D[ i ][ j ] ) ;
  hh->SetCanExtend(false) ;

  Double_t E0_outer[ 10000 ] ;
  Int_t nbins_xy_outer ;
  snprintf( filename , sizeof filename , "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/E0_bb_outer_%s_%s.txt" , used_year.Data(), Magnet.Data() ) ;
  fptr = fopen( filename , "r" ) ;
  ReadFile( fptr , E0_outer , &nbins_xy_outer ) ;
  fclose( fptr ) ;
  //cout<<"Read file E0 bb outer"<< used_year.Data()<<" "<< Magnet.Data() <<" Done"<<endl;

  Double_t bins_x_outer[ 1000 ] , bins_y_outer[ 1000 ] ;
  Int_t nbins_hcal_x_outer , nbins_hcal_y_outer ;
  fptr = fopen( "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/2d_hcal_outer_bins_x.txt" , "r" ) ;
  ReadFile( fptr , bins_x_outer , &nbins_hcal_x_outer ) ;
  fclose( fptr ) ;
  fptr = fopen( "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/2d_hcal_outer_bins_y.txt" , "r" ) ;
  ReadFile( fptr , bins_y_outer , &nbins_hcal_y_outer ) ;
  fclose( fptr ) ;

  snprintf( filename , sizeof filename , "E0_%s" , "outer" ) ;
  Double_t E0_2D_outer[nbins_hcal_x_outer-1][nbins_hcal_y_outer-1] ;
  memcpy( E0_2D_outer , E0_outer , nbins_xy_outer * sizeof( Double_t ) ) ;
  hh = new TH2D( filename , filename , nbins_hcal_x_outer - 1 , bins_x_outer ,
                 nbins_hcal_y_outer - 1 , bins_y_outer ) ;
  for ( i = 0 ; i < nbins_hcal_x_outer - 1 ; ++i )
    for ( j = 0 ; j < nbins_hcal_y_outer - 1 ; ++j )
      hh->SetBinContent( i+1 , j+1 , E0_2D_outer[ i ][ j ] ) ;
  hh->SetCanExtend(false) ;

  // tables of betarel
  Double_t betarel[ 10000 ] ;
  snprintf( filename , sizeof filename , "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/betarel_inner_%s_%s.txt" , used_year.Data(), Magnet.Data() ) ;
  fptr = fopen( filename , "r" ) ;
  ReadFile( fptr , betarel , &nbins_xy ) ;
  fclose( fptr ) ;

  snprintf( filename , sizeof filename , "betarel_%s" , "inner" ) ;
  Double_t betarel_2D[nbins_hcal_x-1][nbins_hcal_y-1] ;
  memcpy( betarel_2D , betarel , nbins_xy * sizeof( Double_t ) ) ;
  hh = new TH2D( filename , filename , nbins_hcal_x - 1 , bins_x , nbins_hcal_y - 1 , bins_y ) ;
  for ( i = 0 ; i < nbins_hcal_x - 1 ; ++i )
    for ( j = 0 ; j < nbins_hcal_y - 1 ; ++j )
      hh->SetBinContent( i+1 , j+1 , betarel_2D[ i ][ j ] ) ;
  hh->SetCanExtend(false) ;

  Double_t betarel_outer[ 10000 ] ;
  snprintf( filename , sizeof filename , "/home/felician/B2Ktautau/TriggerCorr/L0Hadron_corr_inputs/betarel_outer_%s_%s.txt" , used_year.Data(), Magnet.Data() ) ;
  fptr = fopen( filename , "r" ) ;
  ReadFile( fptr , betarel_outer , &nbins_xy_outer ) ;
  fclose( fptr ) ;

  snprintf( filename , sizeof filename , "betarel_%s" , "outer" ) ;
  Double_t betarel_2D_outer[nbins_hcal_x_outer-1][nbins_hcal_y_outer-1] ;
  memcpy( betarel_2D_outer , betarel_outer , nbins_xy_outer * sizeof( Double_t ) ) ;
  hh = new TH2D( filename , filename , nbins_hcal_x_outer - 1 , bins_x_outer ,
                 nbins_hcal_y_outer - 1 , bins_y_outer ) ;
  for ( i = 0 ; i < nbins_hcal_x_outer - 1  ; ++i )
    for ( j = 0 ; j < nbins_hcal_y_outer - 1 ; ++j )
      hh->SetBinContent( i+1 , j+1 , betarel_2D_outer[ i ][ j ] ) ;
  hh->SetCanExtend(false) ;

  return ;
}

bool is_digits(const std::string &str)
{
    return std::all_of(str.begin(), str.end(), ::isdigit); // C++11
}

bool isInside(const std::string & str, char c)
{
    return str.find(c) != std::string::npos;
}