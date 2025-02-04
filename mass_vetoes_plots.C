void draw_mass_plot(Int_t entry, TChain* t1_2016, TChain* t2_2016, TChain* t3_2016, TString name, TCut mass_vetoe_cuts, Int_t Npar);
Int_t n_choose_k(Int_t n, Int_t k);

void mass_vetoes_plots()
{
    // Input files
    // 3pi3pi MC
    cout << "3pi3pi MC" << endl;
    TFileCollection* fc1_2016 = new TFileCollection("fc1_2016", "fc1_2016", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2016/Species_10/post_sel_tree.txt");
    TFileCollection* fc1_2017 = new TFileCollection("fc1_2017", "fc1_2017", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2017/Species_10/post_sel_tree.txt");
    TFileCollection* fc1_2018 = new TFileCollection("fc1_2018", "fc1_2018", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2018/Species_10/post_sel_tree.txt");

    TChain* t1_2016 = new TChain("DecayTree");
    TChain* t1_2017 = new TChain("DecayTree");
    TChain* t1_2018 = new TChain("DecayTree");

    t1_2016->AddFileInfoList((TCollection*)fc1_2016->GetList());
    t1_2017->AddFileInfoList((TCollection*)fc1_2017->GetList());
    t1_2018->AddFileInfoList((TCollection*)fc1_2018->GetList());

    cout << "2016: " << t1_2016->GetEntries() << endl;
    cout << "2017: " << t1_2017->GetEntries() << endl;
    cout << "2018: " << t1_2018->GetEntries() << endl;

    t1_2016->Add(t1_2017);
    t1_2016->Add(t1_2018);

    // RS data
    cout << "RS data" << endl;
    TFileCollection* fc2_2016 = new TFileCollection("fc2_2016", "fc2_2016", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2016/Species_2/post_sel_tree.txt", 10);
    TFileCollection* fc2_2017 = new TFileCollection("fc2_2017", "fc2_2017", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2017/Species_2/post_sel_tree.txt", 10);
    TFileCollection* fc2_2018 = new TFileCollection("fc2_2018", "fc2_2018", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2018/Species_2/post_sel_tree.txt", 10);

    TChain* t2_2016 = new TChain("DecayTree");
    TChain* t2_2017 = new TChain("DecayTree");
    TChain* t2_2018 = new TChain("DecayTree");

    t2_2016->AddFileInfoList((TCollection*)fc2_2016->GetList());
    t2_2017->AddFileInfoList((TCollection*)fc2_2017->GetList());
    t2_2018->AddFileInfoList((TCollection*)fc2_2018->GetList());

    cout << "2016: " << t2_2016->GetEntries() << endl;
    cout << "2017: " << t2_2017->GetEntries() << endl;
    cout << "2018: " << t2_2018->GetEntries() << endl;

    t2_2016->Add(t2_2017);
    t2_2016->Add(t2_2018);

    // WS data
    cout << "WS data" << endl;
    TFileCollection* fc3_2016 = new TFileCollection("fc3_2016", "fc3_2016", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2016/Species_3/post_sel_tree.txt", 10);
    TFileCollection* fc3_2017 = new TFileCollection("fc3_2017", "fc3_2017", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2017/Species_3/post_sel_tree.txt", 10);
    TFileCollection* fc3_2018 = new TFileCollection("fc3_2018", "fc3_2018", "/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/2018/Species_3/post_sel_tree.txt", 10);

    TChain* t3_2016 = new TChain("DecayTree");
    TChain* t3_2017 = new TChain("DecayTree");
    TChain* t3_2018 = new TChain("DecayTree");

    t3_2016->AddFileInfoList((TCollection*)fc3_2016->GetList());
    t3_2017->AddFileInfoList((TCollection*)fc3_2017->GetList());
    t3_2018->AddFileInfoList((TCollection*)fc3_2018->GetList());

    cout << "2016: " << t3_2016->GetEntries() << endl;
    cout << "2017: " << t3_2017->GetEntries() << endl;
    cout << "2018: " << t3_2018->GetEntries() << endl;

    t3_2016->Add(t3_2017);
    t3_2016->Add(t3_2018);

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

    TCut mass_vetoe_cuts = "";

    for(int i1 = 0; i1 < 7; i1++)
    {
        for(int i2 = 1; i2 < 7; i2++)
        {   
            if(i2 > i1)
            {
                TString name2 = Form("Bp_M%i%i",i1,i2);
                draw_mass_plot(a2, t1_2016, t2_2016, t3_2016, name2, mass_vetoe_cuts, 2);
                for(int i3 = 1; i3 < 7; i3++)
                {
                    if(i3 > i2)
                    {
                        TString name3 = Form("Bp_M%i%i%i",i1,i2,i3);
                        draw_mass_plot(a3, t1_2016, t2_2016, t3_2016, name3, mass_vetoe_cuts, 3);
                        for(int i4 = 1; i4 < 7; i4++)
                        {
                            if(i4 > i3)
                            {
                                TString name4 = Form("Bp_M%i%i%i%i",i1,i2,i3,i4);
                                draw_mass_plot(a4, t1_2016, t2_2016, t3_2016, name4, mass_vetoe_cuts, 4);
                                for(int i5 = 1; i5 < 7; i5++)
                                {
                                    if(i5 > i4)
                                    {
                                        TString name5 = Form("Bp_M%i%i%i%i%i",i1,i2,i3,i4,i5);
                                        draw_mass_plot(a5, t1_2016, t2_2016, t3_2016, name5, mass_vetoe_cuts, 5);
                                        for(int i6 = 1; i6 < 7; i6++)
                                        {   
                                            if(i6 > i5)
                                            {
                                                TString name6 = Form("Bp_M%i%i%i%i%i%i",i1,i2,i3,i4,i5,i6);
                                                draw_mass_plot(a6, t1_2016, t2_2016, t3_2016, name6, mass_vetoe_cuts, 6);
                                                a6 += 1;
                                            }
                                        }
                                        a5 += 1;
                                    }
                                }
                                a4 +=1;
                            }
                        }   
                        a3 += 1;
                    }
                }        
                a2 += 1;
            }
        }
    }
}

void draw_mass_plot(Int_t entry, TChain* t1_2016, TChain* t2_2016, TChain* t3_2016, TString name, TCut mass_vetoe_cuts, Int_t Npar)
{
    gStyle->SetOptStat(0);
    // MC vs RS data vs WS data
    TCanvas c;
    c.cd();

    // Double_t x_min_vec[3] = { t1_2016->GetMinimum(name), t2_2016->GetMinimum(name), t3_2016->GetMinimum(name) };
    // Double_t x_max_vec[3] = { t1_2016->GetMaximum(name), t2_2016->GetMaximum(name), t3_2016->GetMaximum(name) };

    std::vector<Double_t> x_min_vec;
    x_min_vec.push_back(t1_2016->GetMinimum(name));
    x_min_vec.push_back(t2_2016->GetMinimum(name));
    x_min_vec.push_back(t3_2016->GetMinimum(name));

    std::vector<Double_t> x_max_vec;
    x_max_vec.push_back(t1_2016->GetMaximum(name));
    x_max_vec.push_back(t2_2016->GetMaximum(name));
    x_max_vec.push_back(t3_2016->GetMaximum(name));

    Double_t x_min = x_min_vec[ std::distance( x_min_vec.begin(), std::min_element(x_min_vec.begin(), x_min_vec.end()) ) ];
    Double_t x_max = x_max_vec[ std::distance( x_max_vec.begin(), std::max_element(x_max_vec.begin(), x_max_vec.end()) ) ];

    TH1D* h_mc = new TH1D(Form("h_%i_mc_%i",Npar,entry), Form("h_%i_mc_%i",Npar,entry), 100, x_min, x_max);
    TH1D* h_rs_data = new TH1D(Form("h_%i_rs_data_%i",Npar,entry), Form("h_%i_rs_data_%i",Npar,entry), 100, x_min, x_max);
    TH1D* h_ws_data = new TH1D(Form("h_%i_ws_data_%i",Npar,entry), Form("h_%i_ws_data_%i",Npar,entry), 100, x_min, x_max);

    t1_2016->Draw(name+Form(" >> h_%i_mc_%i",Npar,entry)); 
    t2_2016->Draw(name+Form(" >> h_%i_rs_data_%i",Npar,entry));
    t3_2016->Draw(name+Form(" >> h_%i_ws_data_%i",Npar,entry));

    h_mc->SetLineColor(kBlue);
    h_rs_data->SetLineColor(kBlack);
    h_ws_data->SetLineColor(kRed);

    h_mc->GetXaxis()->SetTitle(name+" (MeV)");
    h_mc->GetYaxis()->SetTitle(Form( "Normalized entries / (%.2f MeV)", (x_max - x_min)/100.  ));
    h_mc->SetTitle("");

    h_mc->DrawNormalized();
    h_rs_data->DrawNormalized("same");
    h_ws_data->DrawNormalized("same");

    c.SaveAs(Form("/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/%i_particles/",Npar)+name+".pdf");
}

Int_t n_choose_k(Int_t n, Int_t k)
{
    return TMath::Factorial(n)/( TMath::Factorial(k)*TMath::Factorial(n-k) );
}
