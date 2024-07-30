using namespace RooFit;
using namespace RooStats;

void store_vars_in_workspace(RooWorkspace* w, TTree* t, Int_t n_vars, TString variables[n_vars], Int_t species, Bool_t isKtautauMC);

// D+D-K+ variables
TString variables_DDK[] = {"mass",
                    "Bp_FD_OWNPV", "Dp_FD_ORIVX", "Dm_FD_ORIVX", 
                    "Bp_PT", "Bp_ETA", "Dp_PT", "Dm_PT", "Kp_PT", "Kp_ETA",
                    "Dp_PX", "Dp_PY", "Dp_PZ", "Dm_PX", "Dm_PY", "Dm_PZ", "Kp_PX", "Kp_PY", "Kp_PZ", "Bp_PX", "Bp_PY", "Bp_PZ",
                    "Kp_IPCHI2_OWNPV", "Dp_K_IPCHI2_OWNPV", "Dp_pi1_IPCHI2_OWNPV", "Dp_pi2_IPCHI2_OWNPV", "Dm_K_IPCHI2_OWNPV", "Dm_pi1_IPCHI2_OWNPV", "Dm_pi2_IPCHI2_OWNPV",
                    "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
                    "Dp_ENDVERTEX_X", "Dp_ENDVERTEX_Y", "Dp_ENDVERTEX_Z", "Dp_ENDVERTEX_CHI2",
                    "Dm_ENDVERTEX_X", "Dm_ENDVERTEX_Y", "Dm_ENDVERTEX_Z", "Dm_ENDVERTEX_CHI2",
                    "Dm_M", "Dp_M",
                    "Dp_K_ProbNNk", "Dm_K_ProbNNk", "Kp_ProbNNk", "Dp_pi1_ProbNNk", "Dp_pi2_ProbNNk", "Dm_pi1_ProbNNk", "Dm_pi2_ProbNNk",
                    "Dp_K_ProbNNpi", "Dm_K_ProbNNpi", "Kp_ProbNNpi", "Dp_pi1_ProbNNpi", "Dp_pi2_ProbNNpi", "Dm_pi1_ProbNNpi", "Dm_pi2_ProbNNpi",
                    "Bp_VTXISONUMVTX_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B",
                    "Bp_VTXISONUMVTX_Dp", "Bp_VTXISODCHI2ONETRACK_Dp", "Bp_VTXISODCHI2MASSONETRACK_Dp", "Bp_VTXISODCHI2TWOTRACK_Dp", "Bp_VTXISODCHI2MASSTWOTRACK_Dp",
                    "Bp_VTXISONUMVTX_Dm", "Bp_VTXISODCHI2ONETRACK_Dm", "Bp_VTXISODCHI2MASSONETRACK_Dm", "Bp_VTXISODCHI2TWOTRACK_Dm", "Bp_VTXISODCHI2MASSTWOTRACK_Dm",
                    "Bp_CC_05_MULT_B", "Bp_CC_05_SPT_B", "Bp_CC_05_VPT_B", "Bp_CC_05_PY_B", "Bp_CC_05_PZ_B", "Bp_CC_05_PASYM_B", "Bp_CC_05_PTASYM_B", "Bp_CC_05_PXASYM_B", "Bp_CC_05_PYASYM_B", "Bp_CC_05_PZASYM_B", "Bp_CC_05_DELTAETA_B", "Bp_CC_05_DELTAPHI_B", "Bp_CC_05_PX_B", "Bp_CC_05_IT_B", "Bp_CC_05_MAXPT_Q_B", "Bp_CC_05_MAXPT_PT_B", "Bp_CC_05_MAXPT_PX_B", "Bp_CC_05_MAXPT_PY_B", "Bp_CC_05_MAXPT_PZ_B", "Bp_CC_05_MAXPT_PE_B", "Bp_NC_05_MULT_B", "Bp_NC_05_SPT_B", "Bp_NC_05_VPT_B", "Bp_NC_05_PX_B", "Bp_NC_05_PY_B", "Bp_NC_05_PZ_B", "Bp_NC_05_PASYM_B", "Bp_NC_05_PTASYM_B", "Bp_NC_05_PXASYM_B", "Bp_NC_05_PYASYM_B", "Bp_NC_05_PZASYM_B", "Bp_NC_05_DELTAETA_B", "Bp_NC_05_DELTAPHI_B", "Bp_NC_05_IT_B", "Bp_NC_05_MAXPT_PT_B", "Bp_NC_05_MAXPT_PX_B", "Bp_NC_05_MAXPT_PY_B", "Bp_NC_05_MAXPT_PZ_B", "Bp_CCNC_05_IT_B",
                    "Bp_VTXISOBDTHARDFIRSTVALUE_B", "Bp_VTXISOBDTHARDSECONDVALUE_B", "Bp_VTXISOBDTHARDTHIRDVALUE_B",
                    "Bp_TRKISOBDTFIRSTVALUE_K", "Bp_TRKISOBDTSECONDVALUE_K", "Bp_TRKISOBDTTHIRDVALUE_K",
                    "Bp_TRKISOBDTFIRSTVALUE_Dp_K", "Bp_TRKISOBDTSECONDVALUE_Dp_K", "Bp_TRKISOBDTTHIRDVALUE_Dp_K",
                    "Bp_TRKISOBDTFIRSTVALUE_Dp_pi1", "Bp_TRKISOBDTSECONDVALUE_Dp_pi1", "Bp_TRKISOBDTTHIRDVALUE_Dp_pi1",
                    "Bp_TRKISOBDTFIRSTVALUE_Dp_pi2", "Bp_TRKISOBDTSECONDVALUE_Dp_pi2", "Bp_TRKISOBDTTHIRDVALUE_Dp_pi2",
                    "Bp_TRKISOBDTFIRSTVALUE_Dm_K", "Bp_TRKISOBDTSECONDVALUE_Dm_K", "Bp_TRKISOBDTTHIRDVALUE_Dm_K",
                    "Bp_TRKISOBDTFIRSTVALUE_Dm_pi1", "Bp_TRKISOBDTSECONDVALUE_Dm_pi1", "Bp_TRKISOBDTTHIRDVALUE_Dm_pi1",
                    "Bp_TRKISOBDTFIRSTVALUE_Dm_pi2", "Bp_TRKISOBDTSECONDVALUE_Dm_pi2", "Bp_TRKISOBDTTHIRDVALUE_Dm_pi2",
                    "Kp_TRGHOSTPROB", "nPVs", "nTracks", "nLongTracks", "nDownstreamTracks", "nUpstreamTracks", "nVeloTracks", "nTTracks", 
                    "nVeloClusters", "nITClusters", "nSPDHits"
                    };

// D0bar D+s variables
TString variables_DDs[] = {"mass",
                    "Bp_FD_OWNPV", "D0bar_FD_ORIVX", "Dsp_FD_ORIVX", 
                    "Bp_PT", "Bp_ETA", "D0bar_PT", "Dsp_PT",
                    "D0bar_PX", "D0bar_PY", "D0bar_PZ", "Dsp_PX", "Dsp_PY", "Dsp_PZ", "Bp_PX", "Bp_PY", "Bp_PZ",
                    "D0bar_K_IPCHI2_OWNPV", "D0bar_pi_IPCHI2_OWNPV", "Dsp_K1_IPCHI2_OWNPV", "Dsp_K2_IPCHI2_OWNPV", "Dsp_pi_IPCHI2_OWNPV",
                    "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
                    "D0bar_ENDVERTEX_X", "D0bar_ENDVERTEX_Y", "D0bar_ENDVERTEX_Z", "D0bar_ENDVERTEX_CHI2",
                    "Dsp_ENDVERTEX_X", "Dsp_ENDVERTEX_Y", "Dsp_ENDVERTEX_Z", "Dsp_ENDVERTEX_CHI2",
                    "D0bar_M", "Dsp_M",
                    "D0bar_K_ProbNNk", "D0bar_pi_ProbNNk", "Dsp_K1_ProbNNk", "Dsp_K2_ProbNNk", "Dsp_pi_ProbNNk",
                    "D0bar_K_ProbNNpi", "D0bar_pi_ProbNNpi", "Dsp_K1_ProbNNpi", "Dsp_K2_ProbNNpi", "Dsp_pi_ProbNNpi",
                    "D0bar_K_TRGHOSTPROB", "D0bar_pi_TRGHOSTPROB", "Dsp_K1_TRGHOSTPROB", "Dsp_K2_TRGHOSTPROB", "Dsp_pi_TRGHOSTPROB",
                    "Bp_BSTAUTAUTAUISOBDTFIRSTVALUETAUP", "Bp_BSTAUTAUTAUISOBDTSECONDVALUETAUP", "Bp_BSTAUTAUTAUISOBDTTHIRDVALUETAUP", "Bp_BSTAUTAUTAUISOBDTFIRSTVALUETAUM", "Bp_BSTAUTAUTAUISOBDTSECONDVALUETAUM", "Bp_BSTAUTAUTAUISOBDTTHIRDVALUETAUM",
                    "Bp_BSTAUTAUTAUISOFIRSTVALUETAUP", "Bp_BSTAUTAUTAUISOSECONDVALUETAUP", "Bp_BSTAUTAUTAUISOFIRSTVALUETAUM", "Bp_BSTAUTAUTAUISOSECONDVALUETAUM",
                    "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIM", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUPPIM", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIM", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIP", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUMPIP", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIP", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM2", "Bp_BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM2", "Bp_BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM2",
                    "Bp_BSTAUTAUTRKISOFIRSTVALUETAUPPIM", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUPPIP1", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUPPIP2", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUMPIP", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUMPIM1", "Bp_BSTAUTAUTRKISOFIRSTVALUETAUMPIM2",
                    "nPVs", "nTracks", "nLongTracks", "nDownstreamTracks", "nUpstreamTracks", "nVeloTracks", "nTTracks", "nVeloClusters", "nITClusters", "nSPDHits"
                    };

// D0D0K
TString variables_D0D0K[] = {"mass",
                    "Bp_FD_OWNPV", "D0bar_FD_ORIVX", "D0_FD_ORIVX", 
                    "Bp_PT", "Bp_ETA", "D0bar_PT", "D0_PT",
                    "D0bar_PX", "D0bar_PY", "D0bar_PZ", "D0_PX", "D0_PY", "D0_PZ", "Bp_PX", "Bp_PY", "Bp_PZ",
                    "D0bar_K_IPCHI2_OWNPV", "D0bar_pi1_IPCHI2_OWNPV", "D0bar_pi2_IPCHI2_OWNPV", "D0bar_pi3_IPCHI2_OWNPV", "D0_K_IPCHI2_OWNPV", "D0_pi1_IPCHI2_OWNPV", "D0_pi2_IPCHI2_OWNPV", "D0_pi3_IPCHI2_OWNPV",
                    "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
                    "D0bar_ENDVERTEX_X", "D0bar_ENDVERTEX_Y", "D0bar_ENDVERTEX_Z", "D0bar_ENDVERTEX_CHI2",
                    "D0_ENDVERTEX_X", "D0_ENDVERTEX_Y", "D0_ENDVERTEX_Z", "D0_ENDVERTEX_CHI2",
                    "D0bar_M", "D0_M",
                    "D0bar_K_ProbNNk", "D0bar_pi1_ProbNNk", "D0bar_pi2_ProbNNk", "D0bar_pi3_ProbNNk", "D0_K_ProbNNk", "D0_pi1_ProbNNk", "D0_pi2_ProbNNk", "D0_pi3_ProbNNk",
                    "D0bar_K_ProbNNpi", "D0bar_pi1_ProbNNpi", "D0bar_pi2_ProbNNpi", "D0bar_pi3_ProbNNpi", "D0_K_ProbNNpi", "D0_pi1_ProbNNpi", "D0_pi2_ProbNNpi", "D0_pi3_ProbNNpi",
                    "D0bar_K_TRGHOSTPROB", "D0bar_pi1_TRGHOSTPROB", "D0bar_pi2_TRGHOSTPROB", "D0bar_pi3_TRGHOSTPROB", "D0_K_TRGHOSTPROB", "D0_pi1_TRGHOSTPROB", "D0_pi2_TRGHOSTPROB", "D0_pi3_TRGHOSTPROB",
                    "Bp_CONEANGLE_15", "Bp_CONEMULT_15", "Bp_CONEPTASYM_15",
                    "Bp_CONEANGLE_17", "Bp_CONEMULT_17", "Bp_CONEPTASYM_17",
                    "Bp_CONEANGLE_10", "Bp_CONEMULT_10", "Bp_CONEPTASYM_10",
                    "nPVs", "nTracks", "nLongTracks", "nDownstreamTracks", "nUpstreamTracks", "nVeloTracks", "nTTracks", "nVeloClusters", "nITClusters", "nSPDHits"
                    };
    
// Ktautau variables
TString variables_Ktautau[] = {"mass"};

void createDataset(int year, int species, TString FILES)
{
    Bool_t isKtautauMC = false;
    if((species == 1) || (species == 10) || (species == 11) || (species == 12))
    {
        isKtautauMC = true;
    }

    TFileCollection *fc;
    if(isKtautauMC)
    {
        fc = new TFileCollection("fc", "fc", Form("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201%i/Species_%i/fit_results.txt",year,species));
    }
    else
    {
        fc = new TFileCollection("fc", "fc", FILES);
    }
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    // Require events to pass DTF
    TTree* t_cut = (TTree*)t->CopyTree("");

    // Flatten the DTF B+ mass (vector branch)
    Double_t theMass;
    TBranch* bmass = t_cut->Branch("mass", &theMass);

    Float_t DTF_mass[100];
    Double_t fit_mass;
    if(isKtautauMC)
    {
        t_cut->SetBranchAddress("df_Bp_M", &fit_mass);
    }
    else
    {
        t_cut->SetBranchAddress("Bp_dtf_M", DTF_mass);
    }

    for(int i = 0; i < t_cut->GetEntries(); i++)
    {
        t_cut->GetEntry(i);
        if(isKtautauMC)
        {
            theMass = fit_mass;
        }
        else
        {
            theMass = DTF_mass[0];
        }
        bmass->Fill();
    }

    // Create workspace
    RooWorkspace* w = new RooWorkspace("w");

    // Setup workspace variables
    Int_t n_vars;
    RooArgList arg_list("arg_list");
    if(isKtautauMC || (species == 2) || (species == 3))
    {
        n_vars = sizeof(variables_Ktautau)/sizeof(variables_Ktautau[0]);
        store_vars_in_workspace(w, t_cut, n_vars, variables_Ktautau, species, isKtautauMC);
        for(int i = 0; i < n_vars; i++)
        {
            arg_list.add(*w->var(variables_Ktautau[i]));
        }
    }
    else if((species == 4) || (species == 5) || (species == 6))
    {
        n_vars = sizeof(variables_DDK)/sizeof(variables_DDK[0]);
        store_vars_in_workspace(w, t_cut, n_vars, variables_DDK, species, isKtautauMC);
        for(int i = 0; i < n_vars; i++)
        {
            arg_list.add(*w->var(variables_DDK[i]));
        }   
    }
    else if((species == 7) || (species == 8))
    {
        n_vars = sizeof(variables_DDs)/sizeof(variables_DDs[0]);
        store_vars_in_workspace(w, t_cut, n_vars, variables_DDs, species, isKtautauMC);
        for(int i = 0; i < n_vars; i++)
        {
            arg_list.add(*w->var(variables_DDs[i]));
        }
    }
    else if((species == 9) || (species == 0) || (species == -1))
    {
        n_vars = sizeof(variables_D0D0K)/sizeof(variables_D0D0K[0]);
        store_vars_in_workspace(w, t_cut, n_vars, variables_D0D0K, species, isKtautauMC);
        for(int i = 0; i < n_vars; i++)
        {
            arg_list.add(*w->var(variables_D0D0K[i]));
        }
    }

    // Add variables to a RooDataSet
    RooDataSet* data;
    data = new RooDataSet("data","data",t_cut,arg_list);
    
    w->import(*data);
    w->Print();

    TFile* fout = new TFile( Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_%i/mass_dataset.root",year,species), "RECREATE"); 
    fout->cd();
    w->Write();
    fout->Close();
}

void store_vars_in_workspace(RooWorkspace* w, TTree* t, Int_t n_vars, TString variables[n_vars], Int_t species, Bool_t isKtautauMC)
{    
    std::vector<RooRealVar*> vars;
    for(int i = 0; i < n_vars; i++)
    {
        if(variables[i] == "mass")
        {
            if((species == 4) || (species == 5) || (species == 6)) // D+D-K+
            {
                vars.push_back( new RooRealVar(variables[i], variables[i], 5180, 5600 ) ); 
            }
            else if((species == 7) || (species == 71) || (species == 8) || (species == 81)) // D0bar D+s
            {
                vars.push_back( new RooRealVar(variables[i], variables[i], 5235, 5355 ) ); 
            }
            else if( (species == 9) || (species == 0) || (species == -1) ) // D0D0K
            {
                vars.push_back( new RooRealVar(variables[i], variables[i], 5200, 5600 ) );
            }
            else if(isKtautauMC)
            {
                vars.push_back( new RooRealVar(variables[i], variables[i], 4000, 8000 ) );
            }
        }
        else
        {
            vars.push_back( new RooRealVar(variables[i], variables[i], t->GetMinimum(variables[i]), t->GetMaximum(variables[i]) ) );
        }
        w->import(*vars[i]);
    }
}