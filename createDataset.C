using namespace RooFit;
using namespace RooStats;

void store_vars_in_workspace(RooWorkspace* w, TTree* t, Int_t n_vars, TString variables[n_vars]);

// D0bar D+s variables
TString variables[] = {"mass", "BDT1", "BDT2",
                    "Bp_FD_OWNPV", "D0bar_FD_ORIVX", "Dsp_FD_ORIVX", 
                    "Bp_P", "Bp_PT", "Bp_ETA", "D0bar_PT", "Dsp_PT",
                    "D0bar_PX", "D0bar_PY", "D0bar_PZ", "Dsp_PX", "Dsp_PY", "Dsp_PZ", "Bp_PX", "Bp_PY", "Bp_PZ",
                    "D0bar_K_IPCHI2_OWNPV", "D0bar_pi_IPCHI2_OWNPV", "Dsp_K1_IPCHI2_OWNPV", "Dsp_K2_IPCHI2_OWNPV", "Dsp_pi_IPCHI2_OWNPV",
                    "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
                    "D0bar_ENDVERTEX_X", "D0bar_ENDVERTEX_Y", "D0bar_ENDVERTEX_Z", "D0bar_ENDVERTEX_CHI2",
                    "Dsp_ENDVERTEX_X", "Dsp_ENDVERTEX_Y", "Dsp_ENDVERTEX_Z", "Dsp_ENDVERTEX_CHI2",
                    "D0bar_M", "Dsp_M",
                    "D0bar_K_ProbNNk", "D0bar_pi_ProbNNk", "Dsp_K1_ProbNNk", "Dsp_K2_ProbNNk", "Dsp_pi_ProbNNk",
                    "D0bar_K_ProbNNpi", "D0bar_pi_ProbNNpi", "Dsp_K1_ProbNNpi", "Dsp_K2_ProbNNpi", "Dsp_pi_ProbNNpi",
                    "Bp_VTXISONUMVTX_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B",
                    "Bp_VTXISONUMVTX_D0bar", "Bp_VTXISODCHI2ONETRACK_D0bar", "Bp_VTXISODCHI2MASSONETRACK_D0bar", "Bp_VTXISODCHI2TWOTRACK_D0bar", "Bp_VTXISODCHI2MASSTWOTRACK_D0bar",
                    "Bp_VTXISONUMVTX_Dsp", "Bp_VTXISODCHI2ONETRACK_Dsp", "Bp_VTXISODCHI2MASSONETRACK_Dsp", "Bp_VTXISODCHI2TWOTRACK_Dsp", "Bp_VTXISODCHI2MASSTWOTRACK_Dsp",
                    "Bp_CC_05_MULT_B", "Bp_CC_05_SPT_B", "Bp_CC_05_VPT_B", "Bp_CC_05_PX_B", "Bp_CC_05_PY_B", "Bp_CC_05_PZ_B", "Bp_CC_05_PASYM_B", "Bp_CC_05_PTASYM_B", "Bp_CC_05_PXASYM_B", "Bp_CC_05_PYASYM_B", "Bp_CC_05_PZASYM_B", "Bp_CC_05_DELTAETA_B", "Bp_CC_05_DELTAPHI_B", "Bp_CC_05_IT_B", "Bp_CC_05_MAXPT_Q_B", "Bp_CC_05_MAXPT_PT_B", "Bp_CC_05_MAXPT_PX_B", "Bp_CC_05_MAXPT_PY_B", "Bp_CC_05_MAXPT_PZ_B", "Bp_CC_05_MAXPT_PE_B", 
                    "Bp_NC_05_MULT_B", "Bp_NC_05_SPT_B", "Bp_NC_05_VPT_B", "Bp_NC_05_PX_B", "Bp_NC_05_PY_B", "Bp_NC_05_PZ_B", "Bp_NC_05_PASYM_B", "Bp_NC_05_PTASYM_B", "Bp_NC_05_PXASYM_B", "Bp_NC_05_PYASYM_B", "Bp_NC_05_PZASYM_B", "Bp_NC_05_DELTAETA_B", "Bp_NC_05_DELTAPHI_B", "Bp_NC_05_IT_B", "Bp_NC_05_MAXPT_PT_B", "Bp_NC_05_MAXPT_PX_B", "Bp_NC_05_MAXPT_PY_B", "Bp_NC_05_MAXPT_PZ_B", "Bp_CCNC_05_IT_B",
                    "Bp_CC_05_MULT_D0bar", "Bp_CC_05_SPT_D0bar", "Bp_CC_05_VPT_D0bar", "Bp_CC_05_PX_D0bar", "Bp_CC_05_PY_D0bar", "Bp_CC_05_PZ_D0bar", "Bp_CC_05_PASYM_D0bar", "Bp_CC_05_PTASYM_D0bar", "Bp_CC_05_PXASYM_D0bar", "Bp_CC_05_PYASYM_D0bar", "Bp_CC_05_PZASYM_D0bar", "Bp_CC_05_DELTAETA_D0bar", "Bp_CC_05_DELTAPHI_D0bar", "Bp_CC_05_IT_D0bar", "Bp_CC_05_MAXPT_Q_D0bar", "Bp_CC_05_MAXPT_PT_D0bar", "Bp_CC_05_MAXPT_PX_D0bar", "Bp_CC_05_MAXPT_PY_D0bar", "Bp_CC_05_MAXPT_PZ_D0bar", "Bp_CC_05_MAXPT_PE_D0bar", 
                    "Bp_NC_05_MULT_D0bar", "Bp_NC_05_SPT_D0bar", "Bp_NC_05_VPT_D0bar", "Bp_NC_05_PX_D0bar", "Bp_NC_05_PY_D0bar", "Bp_NC_05_PZ_D0bar", "Bp_NC_05_PASYM_D0bar", "Bp_NC_05_PTASYM_D0bar", "Bp_NC_05_PXASYM_D0bar", "Bp_NC_05_PYASYM_D0bar", "Bp_NC_05_PZASYM_D0bar", "Bp_NC_05_DELTAETA_D0bar", "Bp_NC_05_DELTAPHI_D0bar", "Bp_NC_05_IT_D0bar", "Bp_NC_05_MAXPT_PT_D0bar", "Bp_NC_05_MAXPT_PX_D0bar", "Bp_NC_05_MAXPT_PY_D0bar", "Bp_NC_05_MAXPT_PZ_D0bar", "Bp_CCNC_05_IT_D0bar",
                    "Bp_CC_05_MULT_Dsp", "Bp_CC_05_SPT_Dsp", "Bp_CC_05_VPT_Dsp", "Bp_CC_05_PX_Dsp", "Bp_CC_05_PY_Dsp", "Bp_CC_05_PZ_Dsp", "Bp_CC_05_PASYM_Dsp", "Bp_CC_05_PTASYM_Dsp", "Bp_CC_05_PXASYM_Dsp", "Bp_CC_05_PYASYM_Dsp", "Bp_CC_05_PZASYM_Dsp", "Bp_CC_05_DELTAETA_Dsp", "Bp_CC_05_DELTAPHI_Dsp", "Bp_CC_05_IT_Dsp", "Bp_CC_05_MAXPT_Q_Dsp", "Bp_CC_05_MAXPT_PT_Dsp", "Bp_CC_05_MAXPT_PX_Dsp", "Bp_CC_05_MAXPT_PY_Dsp", "Bp_CC_05_MAXPT_PZ_Dsp", "Bp_CC_05_MAXPT_PE_Dsp", 
                    "Bp_NC_05_MULT_Dsp", "Bp_NC_05_SPT_Dsp", "Bp_NC_05_VPT_Dsp", "Bp_NC_05_PX_Dsp", "Bp_NC_05_PY_Dsp", "Bp_NC_05_PZ_Dsp", "Bp_NC_05_PASYM_Dsp", "Bp_NC_05_PTASYM_Dsp", "Bp_NC_05_PXASYM_Dsp", "Bp_NC_05_PYASYM_Dsp", "Bp_NC_05_PZASYM_Dsp", "Bp_NC_05_DELTAETA_Dsp", "Bp_NC_05_DELTAPHI_Dsp", "Bp_NC_05_IT_Dsp", "Bp_NC_05_MAXPT_PT_Dsp", "Bp_NC_05_MAXPT_PX_Dsp", "Bp_NC_05_MAXPT_PY_Dsp", "Bp_NC_05_MAXPT_PZ_Dsp", "Bp_CCNC_05_IT_Dsp",
                    "Bp_VTXISOBDTHARDFIRSTVALUE_B", "Bp_VTXISOBDTHARDSECONDVALUE_B", "Bp_VTXISOBDTHARDTHIRDVALUE_B",
                    "Bp_VTXISOBDTHARDFIRSTVALUE_D0bar", "Bp_VTXISOBDTHARDSECONDVALUE_D0bar", "Bp_VTXISOBDTHARDTHIRDVALUE_D0bar",
                    "Bp_VTXISOBDTHARDFIRSTVALUE_Dsp", "Bp_VTXISOBDTHARDSECONDVALUE_Dsp", "Bp_VTXISOBDTHARDTHIRDVALUE_Dsp",
                    "Bp_TRKISOBDTFIRSTVALUE_D0bar_K", "Bp_TRKISOBDTSECONDVALUE_D0bar_K", "Bp_TRKISOBDTTHIRDVALUE_D0bar_K",
                    "Bp_TRKISOBDTFIRSTVALUE_D0bar_pi", "Bp_TRKISOBDTSECONDVALUE_D0bar_pi", "Bp_TRKISOBDTTHIRDVALUE_D0bar_pi",
                    "Bp_TRKISOBDTFIRSTVALUE_Dsp_K1", "Bp_TRKISOBDTSECONDVALUE_Dsp_K1", "Bp_TRKISOBDTTHIRDVALUE_Dsp_K1",
                    "Bp_TRKISOBDTFIRSTVALUE_Dsp_K2", "Bp_TRKISOBDTSECONDVALUE_Dsp_K2", "Bp_TRKISOBDTTHIRDVALUE_Dsp_K2",
                    "Bp_TRKISOBDTFIRSTVALUE_Dsp_pi", "Bp_TRKISOBDTSECONDVALUE_Dsp_pi", "Bp_TRKISOBDTTHIRDVALUE_Dsp_pi",
                    "Bp_Bstautau_ISOBDTFIRSTVALUE_taup", "Bp_Bstautau_ISOBDTSECONDVALUE_taup", "Bp_Bstautau_ISOBDTTHIRDVALUE_taup", "Bp_Bstautau_ISOBDTFIRSTVALUE_taum", "Bp_Bstautau_ISOBDTSECONDVALUE_taum", "Bp_Bstautau_ISOBDTTHIRDVALUE_taum",
                    "nPVs", "nTracks", "nLongTracks", "nDownstreamTracks", "nUpstreamTracks", "nVeloTracks", "nTTracks", "nVeloClusters", "nITClusters", "nSPDHits",
                    }; 

void createDataset(int year, int species)
{
    cout << "Opening .root files" << endl;

    TChain* t, *t1;

    if(year != -1) // w/o best candidate selection
    {
        TFileCollection *fc = new TFileCollection("fc", "fc", Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/pre_sel_tree.txt",year,species));
        t = new TChain("DecayTree");
        t->AddFileInfoList((TCollection*)fc->GetList());
        
        TFileCollection *fc1 = new TFileCollection("fc1", "fc1", Form("/panfs/felician/B2Ktautau/workflow/sklearn_response/201%i/Species_%i/bdt_output.txt",year,species));
        t1 = new TChain("XGBoost/DecayTree");
        t1->AddFileInfoList((TCollection*)fc1->GetList());

        Int_t Nentries = t->GetEntries();
        Int_t N1entries = t1->GetEntries();
        if(Nentries != N1entries)
        {
            cout << "Wrong number of entries" << endl;
            return;
        }
        else
        {
            t->AddFriend(t1);
        }
    } 
    else // w/ best candidate selection
    {   
        TFile* f = new TFile(Form("/panfs/felician/B2Ktautau/workflow/multiple_events/Species_%i/single_candidate_post_sel.root",species));
        t = (TChain*)f->Get("DecayTree");
    }

    TFile* fout = new TFile( Form("/panfs/felician/B2Ktautau/workflow/create_dataset/201%i/Species_%i/mass_dataset.root",year,species), "RECREATE"); 
    fout->cd();

    cout << "Creating TTree" << endl;
    TTree* t_cut = (TTree*)t->CopyTree("");
    
    // Flatten the DTF B+ mass (vector branch)
    Double_t theMass;
    TBranch* bmass = t_cut->Branch("mass", &theMass);

    Float_t DTF_mass[100];
    t_cut->SetBranchAddress("Bp_dtf_M", DTF_mass);
    cout << "Filling branches" << endl;
    for(int i = 0; i < t_cut->GetEntries(); i++)
    {
        t_cut->GetEntry(i);
        theMass = DTF_mass[0];
        bmass->Fill();
    }

    // BDT variables
    Double_t the_bdt1, the_bdt2;
    TBranch* bdt1 = t_cut->Branch("BDT1", &the_bdt1);
    TBranch* bdt2 = t_cut->Branch("BDT2", &the_bdt2);

    Double_t BDT1, BDT2;
    if(year != -1)
    {
        t1->SetBranchAddress("BDT1", &BDT1);
        t1->SetBranchAddress("BDT2", &BDT2);
    
        for(int i = 0; i < t1->GetEntries(); i++)
        {
            t1->GetEntry(i);
            the_bdt1 = BDT1;
            the_bdt2 = BDT2;
            bdt1->Fill();
            bdt2->Fill();
        }
    }

    // Create workspace
    RooWorkspace* w = new RooWorkspace("w");

    // Setup workspace variables
    cout << "Adding variables to workspace" << endl;
    Int_t n_vars;
    RooArgList arg_list("arg_list");

    n_vars = sizeof(variables)/sizeof(variables[0]);
    store_vars_in_workspace(w, t_cut, n_vars, variables);
    for(int i = 0; i < n_vars; i++)
    {
        arg_list.add(*w->var(variables[i]));
    }

    if( (species == 7) || (species == 71) || (species == 72) )
    {
        TString variables1[] =  {"D0bar_K_ProbNNk_pidgen_default", "D0bar_pi_ProbNNk_pidgen_default", "Dsp_K1_ProbNNk_pidgen_default", "Dsp_K2_ProbNNk_pidgen_default", "Dsp_pi_ProbNNk_pidgen_default",
            "D0bar_K_ProbNNpi_pidgen_default", "D0bar_pi_ProbNNpi_pidgen_default", "Dsp_K1_ProbNNpi_pidgen_default", "Dsp_K2_ProbNNpi_pidgen_default", "Dsp_pi_ProbNNpi_pidgen_default"};

        Int_t n_vars1;
        n_vars1 = sizeof(variables1)/sizeof(variables1[0]);
        store_vars_in_workspace(w, t_cut, n_vars1, variables1);
        for(int i = 0; i < n_vars1; i++)
        {
            arg_list.add(*w->var(variables1[i]));
        }
    }
    
    // Add variables to a RooDataSet
    cout << "Creating RooDataSet from TTree" << endl;
    RooDataSet* data = new RooDataSet("data","data",t_cut,arg_list);
    w->import(*data);
    w->Print();

    w->Write();
    fout->Close();
}

void store_vars_in_workspace(RooWorkspace* w, TTree* t, Int_t n_vars, TString variables[n_vars])
{    
    std::vector<RooRealVar*> vars;
    for(int i = 0; i < n_vars; i++)
    {
        if(variables[i] == "mass")
        {
            vars.push_back( new RooRealVar(variables[i], variables[i], 5235, 5355 ) ); 
        }
        else
        {
            vars.push_back( new RooRealVar(variables[i], variables[i], t->GetMinimum(variables[i]), t->GetMaximum(variables[i]) ) );
        }
        w->import(*vars[i]);
    }
}