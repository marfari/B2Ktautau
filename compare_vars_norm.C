
using namespace std;

void compare_vars_norm()
{
      /////////////////////////////////////////////////////////////////////////////// Input files /////////////////////////////////////////////////////////////////////////////////////////////////////
      TString FILES_MC = "/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_7/pid_corr.txt"; 
      TString FILES_data = "Files_on_grid/data_D0Dps_2016.txt"; 

      TFileCollection *fc = new TFileCollection("fc", "fc", FILES_MC, 10);
      TChain *t = new TChain("DecayTree");
      t->AddFileInfoList((TCollection*)fc->GetList());

      TFileCollection *fc1 = new TFileCollection("fc1", "fc1", FILES_data, 100);
      TChain *t1 = new TChain("ntuple/DecayTree");
      t1->AddFileInfoList((TCollection*)fc1->GetList());
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      cout << t->GetEntries() << endl;
      cout << t1->GetEntries() << endl;

      TCut truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(Dsp_K1_TRUEID) == 321) && (abs(Dsp_K2_TRUEID) == 321) && (abs(D0bar_pi_TRUEID) == 211) && (abs(Dsp_pi_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(Dsp_TRUEID) == 431) && (abs(Bp_TRUEID) == 521)";

      TCut L0_trigger = "(Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1))";
      TCut HLT1_trigger = "(Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1)";
      TCut HLT2_trigger = "(Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1)";
      TCut trigger = L0_trigger+HLT1_trigger+HLT2_trigger;

      // Some fancy variables
      // IP of tau DV to K+ trajectory
      TString Cx_taup = "(df_DV1y-df_BVy)*df_Kp_PZ - (df_DV1z-df_BVz)*df_Kp_PY";
      TString Cy_taup = "(df_DV1z-df_BVz)*df_Kp_PX - (df_DV1x-df_BVx)*df_Kp_PZ";
      TString Cz_taup = "(df_DV1x-df_BVx)*df_Kp_PY - (df_DV1y-df_BVy)*df_Kp_PX";
      TString C_taup = "TMath::Sqrt( pow("+Cx_taup+",2) + pow("+Cy_taup+",2) + pow("+Cz_taup+",2) )";
      TString IP_taup = "2*"+C_taup+"/(TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) + pow(df_Kp_PZ,2) ))";

      TString Cx_taum = "(df_DV2y-df_BVy)*df_Kp_PZ - (df_DV2z-df_BVz)*df_Kp_PY";
      TString Cy_taum = "(df_DV2z-df_BVz)*df_Kp_PX - (df_DV2x-df_BVx)*df_Kp_PZ";
      TString Cz_taum = "(df_DV2x-df_BVx)*df_Kp_PY - (df_DV2y-df_BVy)*df_Kp_PX";
      TString C_taum = "TMath::Sqrt( pow("+Cx_taum+",2) + pow("+Cy_taum+",2) + pow("+Cz_taum+",2) )";
      TString IP_taum = "2*"+C_taum+"/(TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) + pow(df_Kp_PZ,2) ))";

      // Area of the DV1-DV2-PV triangle
      TString a = "TMath::Sqrt( pow(df_PVx - df_DV1x,2) + pow(df_PVy - df_DV1y,2) + pow(df_PVz - df_DV1z,2) )";
      TString b = "TMath::Sqrt( pow(df_PVx - df_DV2x,2) + pow(df_PVy - df_DV2y,2) + pow(df_PVz - df_DV2z,2) )";
      TString c = "TMath::Sqrt( pow(df_DV1x - df_DV2x,2) + pow(df_DV1y - df_DV2y,2) + pow(df_DV1z - df_DV2z,2) )";
      TString s = "("+a+"+"+b+"+"+c+")*0.5";
      TString A = "TMath::Sqrt("+s+"*("+s+"-"+a+")*("+s+"-"+b+")*("+s+"-"+c+"))";

      // B+ FD to PV
      TString Bp_FD_PV = "TMath::Sqrt( pow(df_BVx - df_PVx,2) + pow(df_BVy - df_PVy,2) + pow(df_BVz - df_PVz,2) )";

      // taus FD to Bv
      TString taup_FD_BV = "TMath::Sqrt( pow(df_DV1x - df_BVx,2) + pow(df_DV1y - df_BVy,2) + pow(df_DV1z - df_BVz,2) )";
      TString taum_FD_BV = "TMath::Sqrt( pow(df_DV2x - df_BVx,2) + pow(df_DV2y - df_BVy,2) + pow(df_DV2z - df_BVz,2) )";

      // distance between the two tau DVs
      TString DV1_DV2_distance = "TMath::Sqrt( pow(df_DV1x - df_DV2x,2) + pow(df_DV1y - df_DV2y,2) + pow(df_DV1z - df_DV2z,2) )";

      // For Ktautau:
      TString variables[] = {"Bp_BPVVD", "D0bar_M", "Dsp_M", "D0bar_K_ProbNNk", "Dsp_K1_ProbNNk", "Dsp_K2_ProbNNk"

                        // "TMath::Max(Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taup,Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taum)", "TMath::Max(Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup,Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum)", "TMath::Max(Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup,Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum)",
                        // "TMath::Max( TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX)  )",
                        // "TMath::Min( TMath::Log10(1-TMath::Abs(taup_DIRA_ORIVX))*TMath::Sign(1,taup_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(taum_DIRA_ORIVX))*TMath::Sign(1,taum_DIRA_ORIVX)  )",
                        // "Kp_P", "Kp_PT", "Bp_P", "Bp_PT", "Bp_BPVVD",
                        // "taup_M", "taum_M", "taum_M12", "taup_M12", "taum_M23", "taup_M23",
                        // "df_PVx", "df_PVy", "df_PVz", "df_BVx", "df_BVy", "df_BVz", "df_DV1x", "df_DV1y", "df_DV1z", "df_DV2x", "df_DV2y", "df_DV2z",
                        // "df_taup_PX", "df_taup_PY", "df_taup_PZ", "TMath::Sqrt( pow(df_taup_PX,2) + pow(df_taup_PY,2) )",
                        // "df_taum_PX", "df_taum_PY", "df_taum_PZ", "TMath::Sqrt( pow(df_taum_PX,2) + pow(df_taum_PY,2) )",
                        // "df_Bp_PX", "df_Bp_PY", "df_Bp_PZ", "TMath::Sqrt( pow(df_Bp_PX,2) + pow(df_Bp_PY,2) )", 
                        // "df_Kp_PX", "df_Kp_PY", "df_Kp_PZ", "TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) )",
                        // "df_antinutau_PX", "df_antinutau_PY", "df_antinutau_PZ", "TMath::Sqrt( pow(df_antinutau_PX,2) + pow(df_antinutau_PY,2) )",
                        // "df_nutau_PX", "df_nutau_PY", "df_nutau_PZ", "TMath::Sqrt( pow(df_nutau_PX,2) + pow(df_nutau_PY,2) )",
                        // "df_chi2", "taup_ENDVERTEX_CHI2", "taum_ENDVERTEX_CHI2",
                        // Bp_FD_PV, taup_FD_BV, taum_FD_BV,
                        // "Kp_ProbNNk_pidgen_default", "taup_pi1_ProbNNk_pidgen_default", "taup_pi2_ProbNNk_pidgen_default", "taup_pi3_ProbNNk_pidgen_default", "taum_pi1_ProbNNk_pidgen_default", "taum_pi2_ProbNNk_pidgen_default", "taum_pi3_ProbNNk_pidgen_default",
                        // "Kp_ProbNNpi_pidgen_default", "taup_pi1_ProbNNpi_pidgen_default", "taup_pi2_ProbNNpi_pidgen_default", "taup_pi3_ProbNNpi_pidgen_default", "taum_pi1_ProbNNpi_pidgen_default", "taum_pi2_ProbNNpi_pidgen_default", "taum_pi3_ProbNNpi_pidgen_default",
                        // "Bp_ETA", "Bp_PHI", "Kp_PHI", "taup_pi1_PHI", "taup_pi2_PHI", "taup_pi3_PHI", "taum_pi1_PHI", "taum_pi2_PHI", "taum_pi3_PHI", "Kp_ETA", "taup_pi1_ETA", "taup_pi2_ETA", "taup_pi3_ETA", "taum_pi1_ETA", "taum_pi2_ETA", "taum_pi3_ETA",
                        // "Kp_IPCHI2_OWNPV", "taup_pi1_IPCHI2_OWNPV", "taup_pi2_IPCHI2_OWNPV", "taup_pi3_IPCHI2_OWNPV", "taum_pi1_IPCHI2_OWNPV", "taum_pi2_IPCHI2_OWNPV", "taum_pi3_IPCHI2_OWNPV",
                        // "Bp_VTXISONUMVTX_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B",
                        // "Bp_VTXISONUMVTX_taup", "Bp_VTXISODCHI2ONETRACK_taup", "Bp_VTXISODCHI2MASSONETRACK_taup", "Bp_VTXISODCHI2TWOTRACK_taup", "Bp_VTXISODCHI2MASSTWOTRACK_taup",
                        // "Bp_VTXISONUMVTX_taum", "Bp_VTXISODCHI2ONETRACK_taum", "Bp_VTXISODCHI2MASSONETRACK_taum", "Bp_VTXISODCHI2TWOTRACK_taum", "Bp_VTXISODCHI2MASSTWOTRACK_taum",
                        // "Bp_CC_05_MULT_B", "Bp_CC_05_SPT_B", "Bp_CC_05_VPT_B", "Bp_CC_05_PX_B", "Bp_CC_05_PY_B", "Bp_CC_05_PZ_B", "Bp_CC_05_PASYM_B", "Bp_CC_05_PTASYM_B", "Bp_CC_05_PXASYM_B", "Bp_CC_05_PYASYM_B", "Bp_CC_05_PZASYM_B", "Bp_CC_05_DELTAETA_B", "Bp_CC_05_DELTAPHI_B", "Bp_CC_05_IT_B", "Bp_CC_05_MAXPT_Q_B", "Bp_CC_05_MAXPT_PT_B", "Bp_CC_05_MAXPT_PX_B", "Bp_CC_05_MAXPT_PY_B", "Bp_CC_05_MAXPT_PZ_B", "Bp_CC_05_MAXPT_PE_B", "Bp_NC_05_MULT_B", "Bp_NC_05_SPT_B", "Bp_NC_05_VPT_B", "Bp_NC_05_PX_B", "Bp_NC_05_PY_B", "Bp_NC_05_PZ_B", "Bp_NC_05_PASYM_B", "Bp_NC_05_PTASYM_B", "Bp_NC_05_PXASYM_B", "Bp_NC_05_PYASYM_B", "Bp_NC_05_PZASYM_B", "Bp_NC_05_DELTAETA_B", "Bp_NC_05_DELTAPHI_B", "Bp_NC_05_IT_B", "Bp_NC_05_MAXPT_PT_B", "Bp_NC_05_MAXPT_PX_B", "Bp_NC_05_MAXPT_PY_B", "Bp_NC_05_MAXPT_PZ_B", "Bp_CCNC_05_IT_B",
                        // "Bp_CC_05_MULT_K", "Bp_CC_05_SPT_K", "Bp_CC_05_VPT_K", "Bp_CC_05_PX_K", "Bp_CC_05_PY_K", "Bp_CC_05_PZ_K", "Bp_CC_05_PASYM_K", "Bp_CC_05_PTASYM_K", "Bp_CC_05_PXASYM_K", "Bp_CC_05_PYASYM_K", "Bp_CC_05_PZASYM_K", "Bp_CC_05_DELTAETA_K", "Bp_CC_05_DELTAPHI_K", "Bp_CC_05_IT_K", "Bp_CC_05_MAXPT_Q_K", "Bp_CC_05_MAXPT_PT_K", "Bp_CC_05_MAXPT_PX_K", "Bp_CC_05_MAXPT_PY_K", "Bp_CC_05_MAXPT_PZ_K", "Bp_CC_05_MAXPT_PE_K", "Bp_NC_05_MULT_K", "Bp_NC_05_SPT_K", "Bp_NC_05_VPT_K", "Bp_NC_05_PX_K", "Bp_NC_05_PY_K", "Bp_NC_05_PZ_K", "Bp_NC_05_PASYM_K", "Bp_NC_05_PTASYM_K", "Bp_NC_05_PXASYM_K", "Bp_NC_05_PYASYM_K", "Bp_NC_05_PZASYM_K", "Bp_NC_05_DELTAETA_K", "Bp_NC_05_DELTAPHI_K", "Bp_NC_05_IT_K", "Bp_NC_05_MAXPT_PT_K", "Bp_NC_05_MAXPT_PX_K", "Bp_NC_05_MAXPT_PY_K", "Bp_NC_05_MAXPT_PZ_K", "Bp_CCNC_05_IT_K",
                        // "Bp_CC_05_MULT_taup", "Bp_CC_05_SPT_taup", "Bp_CC_05_VPT_taup", "Bp_CC_05_PX_taup", "Bp_CC_05_PY_taup", "Bp_CC_05_PZ_taup", "Bp_CC_05_PASYM_taup", "Bp_CC_05_PTASYM_taup", "Bp_CC_05_PXASYM_taup", "Bp_CC_05_PYASYM_taup", "Bp_CC_05_PZASYM_taup", "Bp_CC_05_DELTAETA_taup", "Bp_CC_05_DELTAPHI_taup", "Bp_CC_05_IT_taup", "Bp_CC_05_MAXPT_Q_taup", "Bp_CC_05_MAXPT_PT_taup", "Bp_CC_05_MAXPT_PX_taup", "Bp_CC_05_MAXPT_PY_taup", "Bp_CC_05_MAXPT_PZ_taup", "Bp_CC_05_MAXPT_PE_taup", "Bp_NC_05_MULT_taup", "Bp_NC_05_SPT_taup", "Bp_NC_05_VPT_taup", "Bp_NC_05_PX_taup", "Bp_NC_05_PY_taup", "Bp_NC_05_PZ_taup", "Bp_NC_05_PASYM_taup", "Bp_NC_05_PTASYM_taup", "Bp_NC_05_PXASYM_taup", "Bp_NC_05_PYASYM_taup", "Bp_NC_05_PZASYM_taup", "Bp_NC_05_DELTAETA_taup", "Bp_NC_05_DELTAPHI_taup", "Bp_NC_05_IT_taup", "Bp_NC_05_MAXPT_PT_taup", "Bp_NC_05_MAXPT_PX_taup", "Bp_NC_05_MAXPT_PY_taup", "Bp_NC_05_MAXPT_PZ_taup", "Bp_CCNC_05_IT_taup",
                        // "Bp_CC_05_MULT_taum", "Bp_CC_05_SPT_taum", "Bp_CC_05_VPT_taum", "Bp_CC_05_PX_taum", "Bp_CC_05_PY_taum", "Bp_CC_05_PZ_taum", "Bp_CC_05_PASYM_taum", "Bp_CC_05_PTASYM_taum", "Bp_CC_05_PXASYM_taum", "Bp_CC_05_PYASYM_taum", "Bp_CC_05_PZASYM_taum", "Bp_CC_05_DELTAETA_taum", "Bp_CC_05_DELTAPHI_taum", "Bp_CC_05_IT_taum", "Bp_CC_05_MAXPT_Q_taum", "Bp_CC_05_MAXPT_PT_taum", "Bp_CC_05_MAXPT_PX_taum", "Bp_CC_05_MAXPT_PY_taum", "Bp_CC_05_MAXPT_PZ_taum", "Bp_CC_05_MAXPT_PE_taum", "Bp_NC_05_MULT_taum", "Bp_NC_05_SPT_taum", "Bp_NC_05_VPT_taum", "Bp_NC_05_PX_taum", "Bp_NC_05_PY_taum", "Bp_NC_05_PZ_taum", "Bp_NC_05_PASYM_taum", "Bp_NC_05_PTASYM_taum", "Bp_NC_05_PXASYM_taum", "Bp_NC_05_PYASYM_taum", "Bp_NC_05_PZASYM_taum", "Bp_NC_05_DELTAETA_taum", "Bp_NC_05_DELTAPHI_taum", "Bp_NC_05_IT_taum", "Bp_NC_05_MAXPT_PT_taum", "Bp_NC_05_MAXPT_PX_taum", "Bp_NC_05_MAXPT_PY_taum", "Bp_NC_05_MAXPT_PZ_taum", "Bp_CCNC_05_IT_taum",
                        // "Bp_VTXISOBDTHARDFIRSTVALUE_B", "Bp_VTXISOBDTHARDSECONDVALUE_B", "Bp_VTXISOBDTHARDTHIRDVALUE_B",
                        // "Bp_VTXISOBDTHARDFIRSTVALUE_taup", "Bp_VTXISOBDTHARDSECONDVALUE_taup", "Bp_VTXISOBDTHARDTHIRDVALUE_taup",
                        // "Bp_VTXISOBDTHARDFIRSTVALUE_taum", "Bp_VTXISOBDTHARDSECONDVALUE_taum", "Bp_VTXISOBDTHARDTHIRDVALUE_taum",
                        // "Bp_TRKISOBDTFIRSTVALUE_K", "Bp_TRKISOBDTSECONDVALUE_K", "Bp_TRKISOBDTTHIRDVALUE_K",
                        // "Bp_TRKISOBDTFIRSTVALUE_taup_pi1", "Bp_TRKISOBDTSECONDVALUE_taup_pi1", "Bp_TRKISOBDTTHIRDVALUE_taup_pi1",
                        // "Bp_TRKISOBDTFIRSTVALUE_taup_pi2", "Bp_TRKISOBDTSECONDVALUE_taup_pi2", "Bp_TRKISOBDTTHIRDVALUE_taup_pi2",
                        // "Bp_TRKISOBDTFIRSTVALUE_taup_pi3", "Bp_TRKISOBDTSECONDVALUE_taup_pi3", "Bp_TRKISOBDTTHIRDVALUE_taup_pi3",
                        // "Bp_TRKISOBDTFIRSTVALUE_taum_pi1", "Bp_TRKISOBDTSECONDVALUE_taum_pi1", "Bp_TRKISOBDTTHIRDVALUE_taum_pi1",
                        // "Bp_TRKISOBDTFIRSTVALUE_taum_pi2", "Bp_TRKISOBDTSECONDVALUE_taum_pi2", "Bp_TRKISOBDTTHIRDVALUE_taum_pi2",
                        // "Bp_TRKISOBDTFIRSTVALUE_taum_pi3", "Bp_TRKISOBDTSECONDVALUE_taum_pi3", "Bp_TRKISOBDTTHIRDVALUE_taum_pi3",
                        // "Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taup", "Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup", "Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup", "Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taum", "Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum", "Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum",
                        // "taup_AMAXDOCA", "taup_AMINDOCA", "taup_DOCACHI2MAX",
                        // "taum_AMAXDOCA", "taum_AMINDOCA", "taum_DOCACHI2MAX",
                        // "taup_CosTheta", "taum_CosTheta", "Kp_CosTheta", "taup_pi1_CosTheta", "taup_pi2_CosTheta", "taup_pi3_CosTheta", "taum_pi1_CosTheta", "taum_pi2_CosTheta", "taum_pi3_CosTheta"
                        };
      // For D+D-K:
      // TString variables[] = {"Dp_M", "Dm_M",
      //                      "TMath::Max( TMath::Log10(1-TMath::Abs(Dp_DIRA_ORIVX))*TMath::Sign(1,Dp_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(Dm_DIRA_ORIVX))*TMath::Sign(1,Dm_DIRA_ORIVX) )",
      //                      "TMath::Min( TMath::Log10(1-TMath::Abs(Dp_DIRA_ORIVX))*TMath::Sign(1,Dp_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(Dm_DIRA_ORIVX))*TMath::Sign(1,Dm_DIRA_ORIVX) )",
      //                      "TMath::Max(Dm_M, Dp_M)", "TMath::Max(Dm_M12, Dp_M12)", "TMath::Max(Dm_M23, Dp_M23)", "TMath::Min(Dm_M, Dp_M)", "TMath::Min(Dm_M12, Dp_M12)", "TMath::Min(Dm_M23, Dp_M23)",
      //                      "Bp_FD_OWNPV", "Dp_FD_ORIVX", "Dm_FD_ORIVX",
      //                      "TMath::Log10(TMath::Abs(1-Bp_DIRA_OWNPV))*TMath::Sign(1,Bp_DIRA_OWNPV)",
      //                      "Kp_IPCHI2_OWNPV", "Dp_K_IPCHI2_OWNPV", "Dp_pi1_IPCHI2_OWNPV", "Dp_pi2_IPCHI2_OWNPV", "Dm_K_IPCHI2_OWNPV", "Dm_pi1_IPCHI2_OWNPV", "Dm_pi2_IPCHI2_OWNPV",
      //                      "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
      //                      "Dp_ENDVERTEX_X", "Dp_ENDVERTEX_Y", "Dp_ENDVERTEX_Z", "Dp_ENDVERTEX_CHI2",
      //                      "Dm_ENDVERTEX_X", "Dm_ENDVERTEX_Y", "Dm_ENDVERTEX_Z", "Dm_ENDVERTEX_CHI2",
      //                      "Dp_K_ProbNNk", "Dm_K_ProbNNk", "Kp_ProbNNk", "Dp_pi1_ProbNNk", "Dp_pi2_ProbNNk", "Dm_pi1_ProbNNk", "Dm_pi2_ProbNNk",
      //                      "Dp_K_ProbNNpi", "Dm_K_ProbNNpi", "Kp_ProbNNpi", "Dp_pi1_ProbNNpi", "Dp_pi2_ProbNNpi", "Dm_pi1_ProbNNpi", "Dm_pi2_ProbNNpi",
      //                      "Dp_PX", "Dp_PY", "Dp_PZ", "Dm_PX", "Dm_PY", "Dm_PZ", "Kp_PX", "Kp_PY", "Kp_PZ",
      //                      "Bp_VTXISONUMVTX_B", "Bp_VTXISODCHI2ONETRACK_B", "Bp_VTXISODCHI2MASSONETRACK_B", "Bp_VTXISODCHI2TWOTRACK_B", "Bp_VTXISODCHI2MASSTWOTRACK_B",
      //                      "Bp_VTXISONUMVTX_Dp", "Bp_VTXISODCHI2ONETRACK_Dp", "Bp_VTXISODCHI2MASSONETRACK_Dp", "Bp_VTXISODCHI2TWOTRACK_Dp", "Bp_VTXISODCHI2MASSTWOTRACK_Dp",
      //                      "Bp_VTXISONUMVTX_Dm", "Bp_VTXISODCHI2ONETRACK_Dm", "Bp_VTXISODCHI2MASSONETRACK_Dm", "Bp_VTXISODCHI2TWOTRACK_Dm", "Bp_VTXISODCHI2MASSTWOTRACK_Dm",
      //                      "Bp_CC_05_MULT_B", "Bp_CC_05_SPT_B", "Bp_CC_05_VPT_B", "Bp_CC_05_PX_B", "Bp_CC_05_PY_B", "Bp_CC_05_PZ_B", "Bp_CC_05_PASYM_B", "Bp_CC_05_PTASYM_B", "Bp_CC_05_PXASYM_B", "Bp_CC_05_PYASYM_B", "Bp_CC_05_PZASYM_B", "Bp_CC_05_DELTAETA_B", "Bp_CC_05_DELTAPHI_B", "Bp_CC_05_IT_B", "Bp_CC_05_MAXPT_Q_B", "Bp_CC_05_MAXPT_PT_B", "Bp_CC_05_MAXPT_PX_B", "Bp_CC_05_MAXPT_PY_B", "Bp_CC_05_MAXPT_PZ_B", "Bp_CC_05_MAXPT_PE_B", "Bp_NC_05_MULT_B", "Bp_NC_05_SPT_B", "Bp_NC_05_VPT_B", "Bp_NC_05_PX_B", "Bp_NC_05_PY_B", "Bp_NC_05_PZ_B", "Bp_NC_05_PASYM_B", "Bp_NC_05_PTASYM_B", "Bp_NC_05_PXASYM_B", "Bp_NC_05_PYASYM_B", "Bp_NC_05_PZASYM_B", "Bp_NC_05_DELTAETA_B", "Bp_NC_05_DELTAPHI_B", "Bp_NC_05_IT_B", "Bp_NC_05_MAXPT_PT_B", "Bp_NC_05_MAXPT_PX_B", "Bp_NC_05_MAXPT_PY_B", "Bp_NC_05_MAXPT_PZ_B", "Bp_CCNC_05_IT_B",
      //                      "Bp_CC_10_MULT_B", "Bp_CC_10_SPT_B", "Bp_CC_10_VPT_B", "Bp_CC_10_PX_B", "Bp_CC_10_PY_B", "Bp_CC_10_PZ_B", "Bp_CC_10_PASYM_B", "Bp_CC_10_PTASYM_B", "Bp_CC_10_PXASYM_B", "Bp_CC_10_PYASYM_B", "Bp_CC_10_PZASYM_B", "Bp_CC_10_DELTAETA_B", "Bp_CC_10_DELTAPHI_B", "Bp_CC_10_IT_B", "Bp_CC_10_MAXPT_Q_B", "Bp_CC_10_MAXPT_PT_B", "Bp_CC_10_MAXPT_PX_B", "Bp_CC_10_MAXPT_PY_B", "Bp_CC_10_MAXPT_PZ_B", "Bp_CC_10_MAXPT_PE_B", "Bp_NC_10_MULT_B", "Bp_NC_10_SPT_B", "Bp_NC_10_VPT_B", "Bp_NC_10_PX_B", "Bp_NC_10_PY_B", "Bp_NC_10_PZ_B", "Bp_NC_10_PASYM_B", "Bp_NC_10_PTASYM_B", "Bp_NC_10_PXASYM_B", "Bp_NC_10_PYASYM_B", "Bp_NC_10_PZASYM_B", "Bp_NC_10_DELTAETA_B", "Bp_NC_10_DELTAPHI_B", "Bp_NC_10_IT_B", "Bp_NC_10_MAXPT_PT_B", "Bp_NC_10_MAXPT_PX_B", "Bp_NC_10_MAXPT_PY_B", "Bp_NC_10_MAXPT_PZ_B", "Bp_CCNC_10_IT_B",
      //                      "Bp_CC_15_MULT_B", "Bp_CC_15_SPT_B", "Bp_CC_15_VPT_B", "Bp_CC_15_PX_B", "Bp_CC_15_PY_B", "Bp_CC_15_PZ_B", "Bp_CC_15_PASYM_B", "Bp_CC_15_PTASYM_B", "Bp_CC_15_PXASYM_B", "Bp_CC_15_PYASYM_B", "Bp_CC_15_PZASYM_B", "Bp_CC_15_DELTAETA_B", "Bp_CC_15_DELTAPHI_B", "Bp_CC_15_IT_B", "Bp_CC_15_MAXPT_Q_B", "Bp_CC_15_MAXPT_PT_B", "Bp_CC_15_MAXPT_PX_B", "Bp_CC_15_MAXPT_PY_B", "Bp_CC_15_MAXPT_PZ_B", "Bp_CC_15_MAXPT_PE_B", "Bp_NC_15_MULT_B", "Bp_NC_15_SPT_B", "Bp_NC_15_VPT_B", "Bp_NC_15_PX_B", "Bp_NC_15_PY_B", "Bp_NC_15_PZ_B", "Bp_NC_15_PASYM_B", "Bp_NC_15_PTASYM_B", "Bp_NC_15_PXASYM_B", "Bp_NC_15_PYASYM_B", "Bp_NC_15_PZASYM_B", "Bp_NC_15_DELTAETA_B", "Bp_NC_15_DELTAPHI_B", "Bp_NC_15_IT_B", "Bp_NC_15_MAXPT_PT_B", "Bp_NC_15_MAXPT_PX_B", "Bp_NC_15_MAXPT_PY_B", "Bp_NC_15_MAXPT_PZ_B", "Bp_CCNC_15_IT_B",
      //                      "Bp_CC_20_MULT_B", "Bp_CC_20_SPT_B", "Bp_CC_20_VPT_B", "Bp_CC_20_PX_B", "Bp_CC_20_PY_B", "Bp_CC_20_PZ_B", "Bp_CC_20_PASYM_B", "Bp_CC_20_PTASYM_B", "Bp_CC_20_PXASYM_B", "Bp_CC_20_PYASYM_B", "Bp_CC_20_PZASYM_B", "Bp_CC_20_DELTAETA_B", "Bp_CC_20_DELTAPHI_B", "Bp_CC_20_IT_B", "Bp_CC_20_MAXPT_Q_B", "Bp_CC_20_MAXPT_PT_B", "Bp_CC_20_MAXPT_PX_B", "Bp_CC_20_MAXPT_PY_B", "Bp_CC_20_MAXPT_PZ_B", "Bp_CC_20_MAXPT_PE_B", "Bp_NC_20_MULT_B", "Bp_NC_20_SPT_B", "Bp_NC_20_VPT_B", "Bp_NC_20_PX_B", "Bp_NC_20_PY_B", "Bp_NC_20_PZ_B", "Bp_NC_20_PASYM_B", "Bp_NC_20_PTASYM_B", "Bp_NC_20_PXASYM_B", "Bp_NC_20_PYASYM_B", "Bp_NC_20_PZASYM_B", "Bp_NC_20_DELTAETA_B", "Bp_NC_20_DELTAPHI_B", "Bp_NC_20_IT_B", "Bp_NC_20_MAXPT_PT_B", "Bp_NC_20_MAXPT_PX_B", "Bp_NC_20_MAXPT_PY_B", "Bp_NC_20_MAXPT_PZ_B", "Bp_CCNC_20_IT_B",
      //                      "Bp_VTXISOBDTHARDFIRSTVALUE_B", "Bp_VTXISOBDTHARDSECONDVALUE_B", "Bp_VTXISOBDTHARDTHIRDVALUE_B",
      //                      "Bp_TRKISOBDTFIRSTVALUE_K", "Bp_TRKISOBDTSECONDVALUE_K", "Bp_TRKISOBDTTHIRDVALUE_K",
      //                      "Bp_TRKISOBDTFIRSTVALUE_Dp_K", "Bp_TRKISOBDTSECONDVALUE_Dp_K", "Bp_TRKISOBDTTHIRDVALUE_Dp_K",
      //                      "Bp_TRKISOBDTFIRSTVALUE_Dp_pi1", "Bp_TRKISOBDTSECONDVALUE_Dp_pi1", "Bp_TRKISOBDTTHIRDVALUE_Dp_pi1",
      //                      "Bp_TRKISOBDTFIRSTVALUE_Dp_pi2", "Bp_TRKISOBDTSECONDVALUE_Dp_pi2", "Bp_TRKISOBDTTHIRDVALUE_Dp_pi2",
      //                      "Bp_TRKISOBDTFIRSTVALUE_Dm_K", "Bp_TRKISOBDTSECONDVALUE_Dm_K", "Bp_TRKISOBDTTHIRDVALUE_Dm_K",
      //                      "Bp_TRKISOBDTFIRSTVALUE_Dm_pi1", "Bp_TRKISOBDTSECONDVALUE_Dm_pi1", "Bp_TRKISOBDTTHIRDVALUE_Dm_pi1",
      //                      "Bp_TRKISOBDTFIRSTVALUE_Dm_pi2", "Bp_TRKISOBDTSECONDVALUE_Dm_pi2", "Bp_TRKISOBDTTHIRDVALUE_Dm_pi2",
      //                      "Kp_TRGHOSTPROB"
      //                      };

      // For D0D+s:
      // TString variables[] = {"D0bar_M", "Dsp_M",
      //                      "TMath::Log10(1-TMath::Abs(D0bar_DIRA_ORIVX))*TMath::Sign(1,D0bar_DIRA_ORIVX)",
      //                      "TMath::Log10(1-TMath::Abs(Dsp_DIRA_ORIVX))*TMath::Sign(1,Dsp_DIRA_ORIVX)",
      //                      "D0bar_M12", "Dsp_M12", "D0bar_M23", "Dsp_M23",
      //                      "Bp_FD_OWNPV", "D0bar_FD_ORIVX", "Dsp_FD_ORIVX",
      //                      "TMath::Log10(TMath::Abs(1-Bp_DIRA_OWNPV))*TMath::Sign(1,Bp_DIRA_OWNPV)",
      //                      "D0bar_K_IPCHI2_OWNPV", "D0bar_pi_IPCHI2_OWNPV", "Dsp_K1_IPCHI2_OWNPV", "Dsp_K2_IPCHI2_OWNPV", "Dsp_pi_IPCHI2_OWNPV",
      //                      "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
      //                      "D0bar_ENDVERTEX_X", "D0bar_ENDVERTEX_Y", "D0bar_ENDVERTEX_Z", "D0bar_ENDVERTEX_CHI2",
      //                      "Dsp_ENDVERTEX_X", "Dsp_ENDVERTEX_Y", "Dsp_ENDVERTEX_Z", "Dsp_ENDVERTEX_CHI2",
      //                      "D0bar_K_ProbNNk", "Dsp_K1_ProbNNk", "Dsp_K2_ProbNNk", "D0bar_pi_ProbNNk", "Dsp_pi_ProbNNk",
      //                      "D0bar_K_ProbNNpi", "Dsp_K1_ProbNNpi", "Dsp_K2_ProbNNpi", "D0bar_pi_ProbNNpi", "Dsp_pi_ProbNNpi",
      //                      "D0bar_PX", "D0bar_PY", "D0bar_PZ", "Dsp_PX", "Dsp_PY", "Dsp_PZ", "Bp_PX", "Bp_PY", "Bp_PZ",
      //                      "D0bar_K_TRGHOSTPROB", "Dsp_K1_TRGHOSTPROB", "Dsp_K2_TRGHOSTPROB", 
      //                      "nTracks", "nSPDHits"
      //                      };

      // For D0D0K
      // TString variables[] = {"D0bar_M", "D0_M",
      //                      "TMath::Log10(1-TMath::Abs(D0bar_DIRA_ORIVX))*TMath::Sign(1,D0bar_DIRA_ORIVX)",
      //                      "TMath::Log10(1-TMath::Abs(D0_DIRA_ORIVX))*TMath::Sign(1,D0_DIRA_ORIVX)",
      //                      "D0bar_M12", "D0_M12", "D0_M23", "D0_M23",
      //                      "Bp_FD_OWNPV", "D0bar_FD_ORIVX", "D0_FD_ORIVX",
      //                      "TMath::Log10(TMath::Abs(1-Bp_DIRA_OWNPV))*TMath::Sign(1,Bp_DIRA_OWNPV)",
      //                      "D0bar_K_IPCHI2_OWNPV", "D0bar_pi1_IPCHI2_OWNPV", "D0bar_pi2_IPCHI2_OWNPV", "D0bar_pi3_IPCHI2_OWNPV", 
      //                      "D0_K_IPCHI2_OWNPV", "D0_pi1_IPCHI2_OWNPV", "D0_pi2_IPCHI2_OWNPV", "D0_pi3_IPCHI2_OWNPV",
      //                      "Bp_ENDVERTEX_X", "Bp_ENDVERTEX_Y", "Bp_ENDVERTEX_Z", "Bp_ENDVERTEX_CHI2",
      //                      "D0bar_ENDVERTEX_X", "D0bar_ENDVERTEX_Y", "D0bar_ENDVERTEX_Z", "D0bar_ENDVERTEX_CHI2",
      //                      "D0_ENDVERTEX_X", "D0_ENDVERTEX_Y", "D0_ENDVERTEX_Z", "D0_ENDVERTEX_CHI2",
      //                      "Kp_ProbNNk", "D0bar_K_ProbNNk", "D0_K_ProbNNk", "D0bar_pi1_ProbNNk", "D0bar_pi2_ProbNNk", "D0bar_pi3_ProbNNk", "D0_pi1_ProbNNk", "D0_pi2_ProbNNk", "D0_pi3_ProbNNk",
      //                      "Kp_ProbNNpi", "D0bar_K_ProbNNpi", "D0_K_ProbNNpi", "D0bar_pi1_ProbNNpi", "D0bar_pi2_ProbNNpi", "D0bar_pi3_ProbNNpi", "D0_pi1_ProbNNpi", "D0_pi2_ProbNNpi", "D0_pi3_ProbNNpi",
      //                      "D0bar_PX", "D0bar_PY", "D0bar_PZ", "D0_PX", "D0_PY", "D0_PZ", "Bp_PX", "Bp_PY", "Bp_PZ",
      //                      "Kp_TRGHOSTPROB", "D0bar_K_TRGHOSTPROB", "D0_K_TRGHOSTPROB", 
      //                      "nTracks", "nSPDHits",
      //                      "Bp_CONEMULT_15", "Bp_CONEPTASYM_15", 
      //                      "Bp_CONEMULT_17", "Bp_CONEPTASYM_17", 
      //                      "Bp_CONEMULT_10", "Bp_CONEPTASYM_10"
      //                      };

      int n_var =  sizeof(variables) / sizeof(variables[0]);

      std::vector<TCanvas*> C;
      gStyle->SetOptStat(0);

      std::vector<TH1D*> h1;
      std::vector<TH1D*> h2;

      TString folder_name = "/panfs/felician/B2Ktautau/workflow/compare_vars/2016/DDs_sig_vs_bkg_rect_cuts";

      TString name1 = "MC";
      TString name2 = "Data (right)";

      for(int i = 0; i < n_var; i++){
            double x1_min = t->GetMinimum(variables[i]);
            double x1_max = t->GetMaximum(variables[i]);
            double x2_min = t1->GetMinimum(variables[i]);
            double x2_max = t1->GetMaximum(variables[i]);

            double X_min, X_max;
            if(x1_min < x2_min){X_min = x1_min;}
            else{X_min = x2_min;}
            if(x1_max > x2_max){X_max = x1_max;}
            else{X_max = x2_max;}

            if((variables[i] == "df_Bp_PX") || (variables[i] == "df_Bp_PY") || (variables[i] == "Bp_CC_20_PX_B") || (variables[i] == "Bp_CC_20_PY_B") )
            {
                  X_min = -50000;
                  X_max = 50000;
            }
            else if((variables[i] == "df_Bp_PZ") || (variables[i] == "df_Bp_PE") || (variables[i] == "Bp_CC_20_PZ_B") )
            {
                  X_min = 0;
                  X_max = 800000;
            }
            else if((variables[i] == "df_BVx") || (variables[i] == "df_BVy"))
            {
                  X_min = -10;
                  X_max = 10;
            }
            else if(variables[i] == "df_BVz")
            {
                  X_min = -150;
                  X_max = 150;
            }
            else if((variables[i] == "df_taup_PX") || (variables[i] == "df_taup_PY") || (variables[i] == "df_taum_PX") || (variables[i] == "df_taum_PY") || (variables[i] == "df_antinutau_PX") || (variables[i] == "df_antinutau_PY") || (variables[i] == "df_nutau_PX") || (variables[i] == "df_nutau_PY"))
            {
                  X_min = -20000;
                  X_max = 20000;
            }
            else if((variables[i] == "df_taup_PZ") || (variables[i] == "df_taup_PE") || (variables[i] == "df_taum_PZ") || (variables[i] == "df_taum_PE") || (variables[i] == "df_antinutau_PZ") || (variables[i] == "df_antinutau_PE") || (variables[i] == "df_nutau_PZ") || (variables[i] == "df_nutau_PE"))
            {
                  X_min = 0;
                  X_max = 300000;
            }
            else if( (variables[i] == "Bp_MIPCHI2DV") || (variables[i] == "Kp_MIPCHI2DV") || (variables[i] == "taup_pi1_MIPCHI2DV") || (variables[i] == "taup_pi2_MIPCHI2DV") || (variables[i] == "taup_pi3_MIPCHI2DV") || (variables[i] == "taum_pi1_MIPCHI2DV") || (variables[i] == "taum_pi2_MIPCHI2DV") || (variables[i] == "taum_pi3_MIPCHI2DV") )
            {
                  X_max = 1000;
            }
            else if(variables[i] == "Bp_BPVVDCHI2")
            {
                  X_max = 60;
            }
            else if( (variables[i] == "Dp_M") || (variables[i] == "Dm_M"))
            {
                  X_min = 1800;
                  X_max = 1950;
            }
            else if( (variables[i] == "Bp_dtf_12_M") || (variables[i] == "df_Bp_M") )
            {
                  X_min = 4000;
                  X_max = 8000;
            }
            else if( (variables[i] == "Bp_MCORR") )
            {
                  X_min = 2000;
                  X_max = 8000;
            }
            else if( (variables[i] == "Kp_IPCHI2_OWNPV") || (variables[i] == "taup_pi1_IPCHI2_OWNPV") || (variables[i] == "taup_pi2_IPCHI2_OWNPV") || (variables[i] == "taup_pi2_IPCHI2_OWNPV") || (variables[i] == "taum_pi1_IPCHI2_OWNPV") || (variables[i] == "taum_pi2_IPCHI2_OWNPV") || (variables[i] == "taum_pi3_IPCHI2_OWNPV") )
            {
                  X_max = 10000;
            }
            else if( (variables[i] == "Bp_VTXISODCHI2ONETRACK_B") || (variables[i] == "Bp_VTXISODCHI2TWOTRACK_B") )
            {
                  X_max = 200;
            }
            else if( (variables[i] == "Bp_VTXISODCHI2ONETRACK_taup") || (variables[i] == "Bp_VTXISODCHI2ONETRACK_taum") || (variables[i] == "Bp_VTXISODCHI2TWOTRACK_taup") || (variables[i] == "Bp_VTXISODCHI2TWOTRACK_taum") )
            {
                  X_max = 25;
            }
            else if( (variables[i] == "Bp_CC_20_SPT_B") || (variables[i] == "Bp_CC_20_VPT_B") )
            {
                  X_max = 25000;
            }
            else if( (variables[i] == "Bp_CC_20_PXASYM_B") || (variables[i] == "Bp_CC_20_PYASYM_B") || (variables[i] == "Bp_NC_20_PXASYM_B") || (variables[i] == "Bp_NC_20_PYASYM_B") )
            {
                  X_min = -5;
                  X_max = 5;
            }
            else if( variables[i] == "Bp_CC_20_MAXPT_PT_B" )
            {
                  X_max = 20000;
            }
            else if( (variables[i] == "Bp_CC_20_MAXPT_PX_B") || (variables[i] == "Bp_CC_20_MAXPT_PY_B") )
            {
                  X_min = -20000;
                  X_max = 20000;
            }
            else if( (variables[i] == "Bp_CC_20_MAXPT_PZ_B") || (variables[i] == "Bp_CC_20_MAXPT_PE_B") )
            {
                  X_max = 200000;
            }
            else if(variables[i] == "Bp_TAUCHI2")
            {
                  X_min = 0;
                  X_max = 1000;
            }
            else if( (variables[i] == "taup_TAUCHI2") || (variables[i] == "taum_TAUCHI2") )
            {
                  X_min = 0;
                  X_max = 25;
            }
            else if( (variables[i] == "Bp_VTXISODCHI2MASSONETRACK_B") )
            {
                  X_max = 10000;
            }
            else if( (variables[i] == "Bp_VTXISODCHI2MASSTWOTRACK_B") )
            {
                  X_max = 20000;
            }
            else if( (variables[i] == "Bp_VTXISODCHI2MASSONETRACK_taup") || (variables[i] == "Bp_VTXISODCHI2MASSONETRACK_taum") )
            {
                  X_max = 5000;
            }
            else if(variables[i] == "df_chi2")
            {
                  X_min = 0;
                  X_max = 100;
            }
            else if((variables[i] == "taup_M") || (variables[i] == "taum_M"))
            {
                  X_min = 400;
            }
            
            h1.push_back(new TH1D(Form("h1_%.i",i), Form("h1_%.i",i), 100, X_min, X_max));
            h2.push_back(new TH1D(Form("h2_%.i",i), Form("h2_%.i",i), 100, X_min, X_max));

            if((variables[i] == "D0bar_K_ProbNNk") || (variables[i] == "Dsp_K1_ProbNNk") || (variables[i] == "Dsp_K2_ProbNNk"))
            {
                  t->Draw(variables[i]+"_pidgen_default"+Form(" >> h1_%.i",i), truthMatch+trigger);
            }
            else
            {
                  t->Draw(variables[i]+Form(" >> h1_%.i",i), truthMatch+trigger);
            }
            t1->Draw(variables[i]+Form(" >> h2_%.i",i), trigger+"Bp_dtf_M[0] > 5350");

            C.push_back(new TCanvas());
            C[i]->cd();

            // 3pi3pi MC
            h1[i]->SetMarkerColor(kBlue);
            h1[i]->SetLineColor(kBlue);
            h1[i]->SetMarkerStyle(2);

            // RS data
            h2[i]->SetMarkerColor(kRed);
            h2[i]->SetLineColor(kRed);
            h2[i]->SetMarkerStyle(2);

            TString names[] = {"B^{+} FD from PV (mm)", "m_{#bar{D}^{0}} (MeV)", "m_{D^{+}_{s}} (MeV)", "D0bar_K_ProbNNk", "Dsp_K1_ProbNNk", "Dsp_K2_ProbNNk"};

            h1[i]->SetTitle("");
            h2[i]->SetTitle("");
            h1[i]->GetXaxis()->SetTitle(names[i]);
            h1[i]->GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(h1[i]->GetXaxis()->GetXmax() - h1[i]->GetXaxis()->GetXmin())/100) );
            h2[i]->GetXaxis()->SetTitle(names[i]);
            h2[i]->GetYaxis()->SetTitle( TString::Format("Normalised entries /(%g)",(h2[i]->GetXaxis()->GetXmax() - h2[i]->GetXaxis()->GetXmin())/100) );

            h1[i]->Scale(1/h1[i]->Integral());
            h2[i]->Scale(1/h2[i]->Integral());

            double ymax = TMath::Max(1.3*h1[i]->GetMaximum(), 1.3*h2[i]->GetMaximum());
            h1[i]->GetYaxis()->SetRangeUser(0, ymax);

            h1[i]->Draw();
            h2[i]->Draw("same");

            double x_cut1, x_cut2;
            if(variables[i] == "Bp_BPVVD")
            {
                  x_cut1 = 4;
                  x_cut2 = 80;
            }
            if(variables[i] == "D0bar_M")
            {
                  x_cut1 = 1820;
                  x_cut2 = 1910;
            }
            if(variables[i] == "Dsp_M")
            {
                  x_cut1 = 1930;
                  x_cut2 = 2000;
            }
            if((variables[i] == "D0bar_K_ProbNNk") || (variables[i] == "Dsp_K1_ProbNNk") || (variables[i] == "Dsp_K2_ProbNNk"))
            {
                  x_cut1 = 0.55;
                  x_cut2 = x_cut1;
            }

            TLine* line1 = new TLine(x_cut1, 0, x_cut1, ymax);
            TLine* line2 = new TLine(x_cut2, 0, x_cut2, ymax);

            line1->SetLineColor(kBlack);
            line1->SetLineStyle(9);
            line1->SetLineWidth(2);
      
            line2->SetLineColor(kBlack);
            line2->SetLineStyle(9);
            line2->SetLineWidth(2);

            line1->Draw("same");
            line2->Draw("same");

            TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.05);
            leg->AddEntry(h1[i],name1,"p");
            leg->AddEntry(h2[i],name2,"p");
            leg->Draw("same");

            C[i]->SaveAs(folder_name+"/"+Form("var_%i",i)+".pdf");
      }
      
}