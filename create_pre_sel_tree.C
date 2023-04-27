
void create_pre_sel_tree(int year, int species, int component, TString MC_files, TString MC_component_file, TString DATA_files, int line){

    TFileCollection* fc;
    TChain* t_reco;
    if(species == 0) // MC
    {
        fc = new TFileCollection("fc", "fc", MC_files, 1, line);
        t_reco = new TChain("ntuple/DecayTree");
    }
    else if((species == 1) || (species == 2))
    {
        fc = new TFileCollection("fc", "fc", DATA_files, 1, line);
        if(species == 1) // RS data
        {
            t_reco = new TChain("ntuple/DecayTree");
        }
        else // WS data
        {
            t_reco = new TChain("ntuple_SS/DecayTree");
        }
    }
    else
    {
        cout << "Not a valid species value. Try 0 (MC), 1 (RS data) or 2 (WS data)" << endl;
        return;
    }

    t_reco->AddFileInfoList((TCollection*)fc->GetList());

    TCut MC_component;
    if(species == 0)
    {
        t_reco->AddFriend("DecayTree", MC_component_file);
    
        if(component == 0){MC_component = "component == 0";} // 3pi 3pi
        else if(component == 1){MC_component = "component == 1";} // 3pi 3pi pi0
        else if(component == 2){MC_component = "component == 2";} // 3pi 3pi 2pi0
        else if(component == -1){MC_component = "";} // all 3 components
        else
        {
            cout << "Not a valid component value. Try 0 (3pi3pi), 1 (3pi3pipi0), 2 (3pi3pi2pi0) or -1 (all)" << endl;
            return;
        }
    }
    
    TCut truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)";
    TCut pass_DTF = "(Bp_ConsBp_seq_12_status[0]==0)";
    TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
    TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
    TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
    TCut trigger = L0_trigger+HLT1_trigger+HLT2_trigger;
    TCut others;
    if(year == 8)
    {
        others = "(Bp_ConsBp_seq_12_decayLength[0]>5) && "
                    "(Bp_DIRA_OWNPV>0.99) && "
                    "(TMath::Max(taum_M, taup_M) < 1600) && "
                    "(TMath::Min(taum_M, taup_M) > 800) && "
                    "(TMath::Max(Bp_VTXISONUMVTX_tau2,Bp_VTXISONUMVTX_tau1) < 6) && "
                    "((TMath::Min( TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), "
                                  "TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) ) < 0)) ";
    }
    else
    {
        others = "(Bp_ConsBp_seq_12_decayLength[0]>5) && "
                    "(Bp_DIRA_OWNPV>0.99) && "
                    "(TMath::Max(taum_M, taup_M) < 1600) && "
                    "(TMath::Min(taum_M, taup_M) > 800) && "
                    "((TMath::Min( TMath::Log(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), "
                                  "TMath::Log(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX) ) < 0)) ";
    }
    
    TCut pre_selections;
    if(species == 0) // MC
    {
        pre_selections = truthMatch+pass_DTF+trigger+others+MC_component;
    }
    else // data
    {
        pre_selections = pass_DTF+trigger+others;
    }

    TFile* fout;
    if(species == 0)
    {
        fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201%i/Species_%i/Component_%i/%i.root",year,species,component,line), "RECREATE");
    }
    else
    {
        fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201%i/Species_%i/%i.root",year,species,line), "RECREATE");   
    }

    fout->cd();
    TTree* t_pre_sel = (TTree*)t_reco->CopyTree(pre_selections);
    t_pre_sel->Write();
    fout->Close();

    cout << "Finished successfully" << endl;
}