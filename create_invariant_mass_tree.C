Int_t n_choose_k(Int_t n, Int_t k);

void create_invariant_mass_tree(Int_t year, Int_t species, Int_t line)
{
    TString FILES = Form("/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201%i/Species_%i/pre_sel_tree.txt", year, species);
    TFileCollection *fc = new TFileCollection("fc", "fc", FILES, 1, line);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    Double_t branch_values[7][4];

    TString branch_names[7][4] = 
    {
        {"taup_pi1_PX", "taup_pi1_PY", "taup_pi1_PZ", "taup_pi1_PE"},
        {"taup_pi2_PX", "taup_pi2_PY", "taup_pi2_PZ", "taup_pi2_PE"},
        {"taup_pi3_PX", "taup_pi3_PY", "taup_pi3_PZ", "taup_pi3_PE"},
        {"taum_pi1_PX", "taum_pi1_PY", "taum_pi1_PZ", "taum_pi1_PE"},
        {"taum_pi2_PX", "taum_pi2_PY", "taum_pi2_PZ", "taum_pi2_PE"},
        {"taum_pi3_PX", "taum_pi3_PY", "taum_pi3_PZ", "taum_pi3_PE"},
        {"Kp_PX", "Kp_PY", "Kp_PZ", "Kp_PE"}
    };

    for(int i = 0; i < 7; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            t->SetBranchAddress(branch_names[i][j], &branch_values[i][j]);
        }
    }

    TFile* fout = new TFile(Form("/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201%i/Species_%i/%i.root",year,species,line), "RECREATE");
    TTree* tout = new TTree("DecayTree", "DecayTree");

    Int_t N2 = n_choose_k(7,2);
    Int_t N3 = n_choose_k(7,3);
    Int_t N4 = n_choose_k(7,4);
    Int_t N5 = n_choose_k(7,5);
    Int_t N6 = n_choose_k(7,6);

    Double_t mass_comb_2[N2];
    Double_t mass_comb_3[N3];
    Double_t mass_comb_4[N4];
    Double_t mass_comb_5[N5];
    Double_t mass_comb_6[N6];

    Double_t mass_comb_2_massH[N2];
    Double_t mass_comb_3_massH[N3];
    Double_t mass_comb_4_massH[N4];
    Double_t mass_comb_5_massH[N5];
    Double_t mass_comb_6_massH[N6];

    Int_t a2 = 0;
    Int_t a3 = 0;
    Int_t a4 = 0;
    Int_t a5 = 0;
    Int_t a6 = 0;
   
    for(int i1 = 0; i1 < 7; i1++)
    {
        // 2 particle combination
        for(int i2 = 0; i2 < i1; i2++)
        {   
            tout->Branch(Form("Bp_M%i%i",i1,i2), &mass_comb_2[a2]);   
            tout->Branch(Form("Bp_M%i%i_massH",i1,i2), &mass_comb_2_massH[a2]);                     
            for(int i3 = 0; i3 < i2; i3++)
            {
                tout->Branch(Form("Bp_M%i%i%i",i1,i2,i3), &mass_comb_3[a3]);
                tout->Branch(Form("Bp_M%i%i%i_massH",i1,i2,i3), &mass_comb_3_massH[a3]);
                for(int i4 = 0; i4 < i3; i4++)
                {
                    tout->Branch(Form("Bp_M%i%i%i%i",i1,i2,i3,i4), &mass_comb_4[a4]);
                    tout->Branch(Form("Bp_M%i%i%i%i_massH",i1,i2,i3,i4), &mass_comb_4_massH[a4]);
                    for(int i5 = 0; i5 < i4; i5++)
                    {
                        tout->Branch(Form("Bp_M%i%i%i%i%i",i1,i2,i3,i4,i5), &mass_comb_5[a5]);
                        tout->Branch(Form("Bp_M%i%i%i%i%i_massH",i1,i2,i3,i4,i5), &mass_comb_5_massH[a5]);
                        for(int i6 = 0; i6 < i5; i6++)
                        {   
                            tout->Branch(Form("Bp_M%i%i%i%i%i%i",i1,i2,i3,i4,i5,i6), &mass_comb_6[a6]);
                            tout->Branch(Form("Bp_M%i%i%i%i%i%i_massH",i1,i2,i3,i4,i5,i6), &mass_comb_6_massH[a6]);
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

    Double_t mkaon = 493.677;
    Double_t mpion = 139.57039;

    Int_t num_entries = t->GetEntries();
    for(int evt = 0 ; evt < num_entries; evt++)
    {
        t->GetEntry(evt);

        Int_t b2 = 0;
        Int_t b3 = 0;
        Int_t b4 = 0;
        Int_t b5 = 0;
        Int_t b6 = 0;

        // Loop over mass combinations
        for(int i1 = 0; i1 < 7; i1++)
        {
            Bool_t isKaon = false;
            if(i1 == 6)
            {
                isKaon = true;
            }

            for(int i2 = 0; i2 < i1; i2++)
            {
                ROOT::Math::PxPyPzEVector P1( branch_values[i1][0], branch_values[i1][1], branch_values[i1][2], branch_values[i1][3] );
                ROOT::Math::PxPyPzEVector P2( branch_values[i2][0], branch_values[i2][1], branch_values[i2][2], branch_values[i2][3] );

                Double_t M2_squared = (P1+P2).M2();
                if(M2_squared > 0)
                {
                    mass_comb_2[b2] = TMath::Sqrt( M2_squared );
                }
                else
                {
                    mass_comb_2[b2] = -TMath::Sqrt( TMath::Abs(M2_squared) );
                }

                Double_t M2_squared_massH = isKaon ? M2_squared - pow(mkaon,2) + pow(mpion,2) : M2_squared - pow(mpion,2) + pow(mkaon,2);
                if(M2_squared_massH > 0)
                {
                    mass_comb_2_massH[b2] = TMath::Sqrt(M2_squared_massH);
                }
                else
                {
                    mass_comb_2_massH[b2] = -TMath::Sqrt( TMath::Abs(M2_squared_massH) );
                }

                for(int i3 = 0; i3 < i2; i3++)
                {
                    ROOT::Math::PxPyPzEVector P3( branch_values[i3][0], branch_values[i3][1], branch_values[i3][2], branch_values[i3][3] );

                    Double_t M3_squared = (P1+P2+P3).M2();
                    if(M3_squared > 0)
                    {
                        mass_comb_3[b3] = TMath::Sqrt(M3_squared);
                    }
                    else
                    {
                        mass_comb_3[b3] = -TMath::Sqrt( TMath::Abs(M3_squared) );
                    }

                    Double_t M3_squared_massH = isKaon ? M3_squared - pow(mkaon,2) + pow(mpion,2) : M3_squared - pow(mpion,2) + pow(mkaon,2);
                    if(M3_squared_massH > 0)
                    {
                        mass_comb_3_massH[b3] = TMath::Sqrt(M3_squared_massH);
                    }
                    else
                    {
                        mass_comb_3_massH[b3] = -TMath::Sqrt( TMath::Abs(M3_squared_massH) );
                    }

                    for(int i4 = 0; i4 < i3; i4++)
                    {
                        ROOT::Math::PxPyPzEVector P4( branch_values[i4][0], branch_values[i4][1], branch_values[i4][2], branch_values[i4][3] );

                        Double_t M4_squared = (P1+P2+P3+P4).M2();
                        if(M4_squared > 0)
                        {
                            mass_comb_4[b4] = TMath::Sqrt(M4_squared);
                        }
                        else
                        {
                            mass_comb_4[b4] = -TMath::Sqrt( TMath::Abs(M4_squared) );
                        }

                        Double_t M4_squared_massH = isKaon ? M4_squared - pow(mkaon,2) + pow(mpion,2) : M4_squared - pow(mpion,2) + pow(mkaon,2);
                        if(M4_squared_massH > 0)
                        {
                            mass_comb_4_massH[b4] = TMath::Sqrt(M4_squared_massH);
                        }
                        else
                        {
                            mass_comb_4_massH[b4] = -TMath::Sqrt( TMath::Abs(M4_squared_massH) );
                        }

                        for(int i5 = 0; i5 < i4; i5++)
                        {
                            ROOT::Math::PxPyPzEVector P5( branch_values[i5][0], branch_values[i5][1], branch_values[i5][2], branch_values[i5][3] );

                            Double_t M5_squared = (P1+P2+P3+P4+P5).M2();
                            if(M5_squared > 0)
                            {
                                mass_comb_5[b5] = TMath::Sqrt(M5_squared);
                            }
                            else
                            {
                                mass_comb_5[b5] = -TMath::Sqrt( TMath::Abs(M5_squared) );
                            }

                            Double_t M5_squared_massH = isKaon ? M5_squared - pow(mkaon,2) + pow(mpion,2) : M5_squared - pow(mpion,2) + pow(mkaon,2);
                            if(M5_squared_massH > 0)
                            {
                                mass_comb_5_massH[b5] = TMath::Sqrt(M5_squared_massH);
                            }
                            else
                            {
                                mass_comb_5_massH[b5] = -TMath::Sqrt( TMath::Abs(M5_squared_massH) );
                            }

                            for(int i6 = 0; i6 < i5; i6++)
                            {
                                ROOT::Math::PxPyPzEVector P6( branch_values[i6][0], branch_values[i6][1], branch_values[i6][2], branch_values[i6][3] );
                                
                                Double_t M6_squared = (P1+P2+P3+P4+P5+P6).M2();

                                if(M6_squared > 0)
                                {
                                    mass_comb_6[b6] = TMath::Sqrt(M6_squared);
                                }
                                else
                                {
                                    mass_comb_6[b6] = -TMath::Sqrt( TMath::Abs(M6_squared) );
                                }

                                Double_t M6_squared_massH = isKaon ? M6_squared - pow(mkaon,2) + pow(mpion,2) : M6_squared - pow(mpion,2) + pow(mkaon,2);
                                if(M6_squared_massH > 0)
                                {
                                    mass_comb_6_massH[b6] = TMath::Sqrt(M6_squared_massH);
                                }
                                else
                                {
                                    mass_comb_6_massH[b6] = -TMath::Sqrt( TMath::Abs(M6_squared_massH) );
                                }
                            }   
                            b2 += 1;
                        }
                        b3 += 1;
                    }
                    b4 += 1;
                }
                b5 += 1;
            }
            b6 += 1;
        }
    
        tout->Fill();
    }

    fout->cd();
    tout->Write();
    fout->Close();

}

Int_t n_choose_k(Int_t n, Int_t k)
{
    return TMath::Factorial(n)/( TMath::Factorial(k)*TMath::Factorial(n-k) );
}
