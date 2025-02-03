Int_t n_choose_k(Int_t n, Int_t k);

void mass_vetoes_plots()
{
    // Invariant mass files
    // 3pi3pi MC
    TFileCollection *fc1_2016 = new TFileCollection("fc1_2016", "fc1_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_10/invariant_mass_tree.txt");
    TFileCollection *fc1_2017 = new TFileCollection("fc1_2017", "fc1_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_10/invariant_mass_tree.txt");
    TFileCollection *fc1_2018 = new TFileCollection("fc1_2018", "fc1_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_10/invariant_mass_tree.txt");

    TChain* t1_2016 = new TChain("DecayTree");
    TChain* t1_2017 = new TChain("DecayTree");
    TChain* t1_2018 = new TChain("DecayTree");

    t1_2016->AddFileInfoList((TCollection*)fc1_2016->GetList());
    t1_2017->AddFileInfoList((TCollection*)fc1_2017->GetList());
    t1_2018->AddFileInfoList((TCollection*)fc1_2018->GetList());

    // RS data
    TFileCollection* fc2_2016 = new TFileCollection("fc2_2016", "fc2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_2/invariant_mass_tree.txt", 100);
    TFileCollection* fc2_2017 = new TFileCollection("fc2_2017", "fc2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_2/invariant_mass_tree.txt", 100);
    TFileCollection* fc2_2018 = new TFileCollection("fc2_2018", "fc2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_2/invariant_mass_tree.txt", 100);

    TChain* t2_2016 = new TChain("DecayTree");
    TChain* t2_2017 = new TChain("DecayTree");
    TChain* t2_2018 = new TChain("DecayTree");

    t2_2016->AddFileInfoList((TCollection*)fc2_2016->GetList());
    t2_2017->AddFileInfoList((TCollection*)fc2_2017->GetList());
    t2_2018->AddFileInfoList((TCollection*)fc2_2018->GetList());

    // WS data
    TFileCollection* fc3_2016 = new TFileCollection("fc3_2016", "fc3_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_3/invariant_mass_tree.txt", 100);
    TFileCollection* fc3_2017 = new TFileCollection("fc3_2017", "fc3_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_3/invariant_mass_tree.txt", 100);
    TFileCollection* fc3_2018 = new TFileCollection("fc3_2018", "fc3_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_3/invariant_mass_tree.txt", 100);

    TChain* t3_2016 = new TChain("DecayTree");
    TChain* t3_2017 = new TChain("DecayTree");
    TChain* t3_2018 = new TChain("DecayTree");

    t3_2016->AddFileInfoList((TCollection*)fc3_2016->GetList());
    t3_2017->AddFileInfoList((TCollection*)fc3_2017->GetList());
    t3_2018->AddFileInfoList((TCollection*)fc3_2018->GetList());

    // BDT response files
    // 3pi3pi MC
    TFileCollection* fc4_2016 = new TFileCollection("fc4_2016", "fc4_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_10/bdt_output.txt");
    TFileCollection* fc4_2017 = new TFileCollection("fc4_2017", "fc4_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_10/bdt_output.txt");
    TFileCollection* fc4_2018 = new TFileCollection("fc4_2018", "fc4_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_10/bdt_output.txt");

    TChain* t4_2016 = new TChain("XGBoost/DecayTree");
    TChain* t4_2017 = new TChain("XGBoost/DecayTree");
    TChain* t4_2018 = new TChain("XGBoost/DecayTree");

    t4_2016->AddFileInfoList((TCollection*)fc4_2016->GetList());
    t4_2017->AddFileInfoList((TCollection*)fc4_2017->GetList());
    t4_2018->AddFileInfoList((TCollection*)fc4_2018->GetList());

    cout << "3pi3pi MC" << endl;
    cout << t1_2016->GetEntries() << "  " << t1_2017->GetEntries() << "  " << t1_2018->GetEntries() << endl;
    cout << t4_2016->GetEntries() << "  " << t4_2017->GetEntries() << "  " << t4_2018->GetEntries() << endl;

    t1_2016->Add(t1_2017);
    t1_2016->Add(t1_2018);

    t4_2016->Add(t4_2017);
    t4_2016->Add(t4_2018);

    // cout << t1_2016->GetEntries() << endl;
    // cout << t4_2016->GetEntries() << endl;
    t1_2016->AddFriend(t4_2016);

    // RS data
    TFileCollection* fc5_2016 = new TFileCollection("fc5_2016", "fc5_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_2/bdt_output.txt", 100);
    TFileCollection* fc5_2017 = new TFileCollection("fc5_2017", "fc5_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_2/bdt_output.txt", 100);
    TFileCollection* fc5_2018 = new TFileCollection("fc5_2018", "fc5_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_2/bdt_output.txt", 100);

    TChain* t5_2016 = new TChain("XGBoost/DecayTree");
    TChain* t5_2017 = new TChain("XGBoost/DecayTree");
    TChain* t5_2018 = new TChain("XGBoost/DecayTree");

    t5_2016->AddFileInfoList((TCollection*)fc5_2016->GetList());
    t5_2017->AddFileInfoList((TCollection*)fc5_2017->GetList());
    t5_2018->AddFileInfoList((TCollection*)fc5_2018->GetList());

    cout << "RS data" << endl;
    cout << t2_2016->GetEntries() << "  " << t2_2017->GetEntries() << "  " << t2_2018->GetEntries() << endl;
    cout << t5_2016->GetEntries() << "  " << t5_2017->GetEntries() << "  " << t5_2018->GetEntries() << endl;

    t2_2016->Add(t2_2017);
    t2_2016->Add(t2_2018);

    t5_2016->Add(t5_2017);
    t5_2016->Add(t5_2018);

    // cout << t2_2016->GetEntries() << endl;
    // cout << t5_2016->GetEntries() << endl;
    t2_2016->AddFriend(t5_2016);

    // WS data
    TFileCollection* fc6_2016 = new TFileCollection("fc6_2016", "fc6_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_3/bdt_output.txt", 100);
    TFileCollection* fc6_2017 = new TFileCollection("fc6_2017", "fc6_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_3/bdt_output.txt", 100);
    TFileCollection* fc6_2018 = new TFileCollection("fc6_2018", "fc6_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_3/bdt_output.txt", 100);

    TChain* t6_2016 = new TChain("XGBoost/DecayTree");
    TChain* t6_2017 = new TChain("XGBoost/DecayTree");
    TChain* t6_2018 = new TChain("XGBoost/DecayTree");

    t6_2016->AddFileInfoList((TCollection*)fc6_2016->GetList());
    t6_2017->AddFileInfoList((TCollection*)fc6_2017->GetList());
    t6_2018->AddFileInfoList((TCollection*)fc6_2018->GetList());

    cout << "WS data" << endl;
    cout << t3_2016->GetEntries() << "  " << t3_2017->GetEntries() << "  " << t3_2018->GetEntries() << endl;
    cout << t6_2016->GetEntries() << "  " << t6_2017->GetEntries() << "  " << t6_2018->GetEntries() << endl;

    t3_2016->Add(t3_2017);
    t3_2016->Add(t3_2018);

    t6_2016->Add(t6_2017);
    t6_2016->Add(t6_2018);

    // cout << t3_2016->GetEntries() << endl;
    // cout << t6_2016->GetEntries() << endl;
    t3_2016->AddFriend(t6_2016);

    Int_t N2 = n_choose_k(7,2);
    Int_t N3 = n_choose_k(7,3);
    Int_t N4 = n_choose_k(7,4);
    Int_t N5 = n_choose_k(7,5);
    Int_t N6 = n_choose_k(7,6);

    // Invariant mass plots (correct mass hypothesis) 
    TCanvas* c_2[N2], c_3[N3], c_4[N4], c_5[N5], c_6[N6];
    TH1D* h_2_mc[N2], h_3_mc[N3], h_4_mc[N4], h_5_mc[N5], h_6_mc[N6];
    TH1D* h_2_rs_data[N2], h_3_rs_data[N3], h_4_rs_data[N4], h_5_rs_data[N5], h_6_rs_data[N6];
    TH1D* h_2_ws_data[N2], h_3_ws_data[N3], h_4_ws_data[N4], h_5_ws_data[N5], h_6_ws_data[N6];

    Int_t a2 = 0;
    Int_t a3 = 0;
    Int_t a4 = 0;
    Int_t a5 = 0;
    Int_t a6 = 0;

    TCut bdt_cut = "(BDT1 > 0.968) && (BDT2 > 0.953)";

    for(int i1 = 0; i1 < 7; i1++)
    {
        for(int i2 = 0; i2 < i1; i2++)
        {   
            TString name2 = Form("Bp_M%i%i",i1,i2);

            c_2[a2] = new TCanvas(Form("c_2_%i",a2), Form("c_2_%i",a2));
            c_2[a2]->cd();

            Double_t x_min = t1_2016->GetMinimum(name2);
            Double_t x_max = t1_2016->GetMaximum(name2);

            h_2_mc[a2] = new TH1D(Form("h_2_mc_%i",a2), Form("h_2_mc_%i",a2), 100, x_min, x_max );
            h_2_rs_data[a2] = new TH1D(Form("h_2_rs_data_%i",a2), Form("h_2_rs_data_%i",a2), 100, x_min, x_max );
            h_2_ws_data[a2] =  new TH1D(Form("h_2_ws_data_%i",a2), Form("h_2_ws_data_%i",a2), 100, x_min, x_max );

            t1_2016->Draw(name2+Form(" >> h_2_mc_%i",a2), bdt_cut); 
            t2_2016->Draw(name2+Form(" >> h_2_rs_data_%i",a2));
            t3_2016->Draw(name2+Form(" >> h_2_ws_data_%i",a2));

            h_2_mc[a2]->SetLineColor(kBlue);
            h_2_rs_data[a2]->SetLineColor(kBlack);
            h_2_ws_data[a2]->SetLineColor(kRed);

            h_2_mc[a2]->GetXaxis()->SetTitle(name2+" (MeV)");
            h_2_mc[a2]->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/100.  ));
            h_2_mc[a2]->SetTitle("");

            h_2_mc[a2]->DrawNormalized();
            h_2_rs_data[a2]->DrawNormalized("same");
            h_2_ws_data[a2]->DrawNormalized("same");

            c_2[a2]->SaveAs("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/2_particles/"+name2+".pdf");

            for(int i3 = 0; i3 < i2; i3++)
            {
                // tout->Branch(Form("Bp_M%i%i%i",i1,i2,i3), &mass_comb_3[a3]);
                // tout->Branch(Form("Bp_M%i%i%i_massH",i1,i2,i3), &mass_comb_3_massH[a3]);
                for(int i4 = 0; i4 < i3; i4++)
                {
                    // tout->Branch(Form("Bp_M%i%i%i%i",i1,i2,i3,i4), &mass_comb_4[a4]);
                    // tout->Branch(Form("Bp_M%i%i%i%i_massH",i1,i2,i3,i4), &mass_comb_4_massH[a4]);
                    for(int i5 = 0; i5 < i4; i5++)
                    {
                        // tout->Branch(Form("Bp_M%i%i%i%i%i",i1,i2,i3,i4,i5), &mass_comb_5[a5]);
                        // tout->Branch(Form("Bp_M%i%i%i%i%i_massH",i1,i2,i3,i4,i5), &mass_comb_5_massH[a5]);
                        for(int i6 = 0; i6 < i5; i6++)
                        {   
                            // tout->Branch(Form("Bp_M%i%i%i%i%i%i",i1,i2,i3,i4,i5,i6), &mass_comb_6[a6]);
                            // tout->Branch(Form("Bp_M%i%i%i%i%i%i_massH",i1,i2,i3,i4,i5,i6), &mass_comb_6_massH[a6]);
                        }
                        a2 += 1;
                    }
                    a3 +=1;
                }   
                a4 += 1;
            }        
            a5 += 1;
        }
        a6 += 1;
    }


}

Int_t n_choose_k(Int_t n, Int_t k)
{
    return TMath::Factorial(n)/( TMath::Factorial(k)*TMath::Factorial(n-k) );
}
