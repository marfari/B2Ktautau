import sys
import ROOT


def main(argv):
    ROOT.gStyle.SetOptStat(0)

    ###################################################### Tight truth-matched Ktautau MC #####################################################################
    cut_Ktautau = "(df_status == 0) && (df_Bp_M > 4000) && (df_Bp_M < 8000)"
    
    # pre-sel tree
    fc_BuKtautau_mc_2016 = ROOT.TFileCollection("fc_BuKtautau_mc_2016", "fc_BuKtautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1/pre_sel_tree.txt")
    fc_BuKtautau_mc_2017 = ROOT.TFileCollection("fc_BuKtautau_mc_2017", "fc_BuKtautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1/pre_sel_tree.txt")
    fc_BuKtautau_mc_2018 = ROOT.TFileCollection("fc_BuKtautau_mc_2018", "fc_BuKtautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1/pre_sel_tree.txt")

    t_BuKtautau_mc_2016 = ROOT.TChain("DecayTree")
    t_BuKtautau_mc_2017 = ROOT.TChain("DecayTree")
    t_BuKtautau_mc_2018 = ROOT.TChain("DecayTree")

    t_BuKtautau_mc_2016.AddFileInfoList(fc_BuKtautau_mc_2016.GetList())
    t_BuKtautau_mc_2017.AddFileInfoList(fc_BuKtautau_mc_2017.GetList())
    t_BuKtautau_mc_2018.AddFileInfoList(fc_BuKtautau_mc_2018.GetList())

    t_BuKtautau_mc_2016.Add(t_BuKtautau_mc_2017)
    t_BuKtautau_mc_2016.Add(t_BuKtautau_mc_2018)

    # gsl 
    fc1_BuKtautau_mc_2016 = ROOT.TFileCollection("fc1_BuKtautau_mc_2016", "fc1_BuKtautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1/fit_results.txt")
    fc1_BuKtautau_mc_2017 = ROOT.TFileCollection("fc1_BuKtautau_mc_2017", "fc1_BuKtautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1/fit_results.txt")
    fc1_BuKtautau_mc_2018 = ROOT.TFileCollection("fc1_BuKtautau_mc_2018", "fc1_BuKtautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1/fit_results.txt")

    t1_BuKtautau_mc_2016 = ROOT.TChain("DecayTree")
    t1_BuKtautau_mc_2017 = ROOT.TChain("DecayTree")
    t1_BuKtautau_mc_2018 = ROOT.TChain("DecayTree")

    t1_BuKtautau_mc_2016.AddFileInfoList(fc1_BuKtautau_mc_2016.GetList())
    t1_BuKtautau_mc_2017.AddFileInfoList(fc1_BuKtautau_mc_2017.GetList())
    t1_BuKtautau_mc_2018.AddFileInfoList(fc1_BuKtautau_mc_2018.GetList())

    t1_BuKtautau_mc_2016.Add(t1_BuKtautau_mc_2017)
    t1_BuKtautau_mc_2016.Add(t1_BuKtautau_mc_2018)

    # sklearn
    fc2_BuKtautau_mc_2016 = ROOT.TFileCollection("fc2_BuKtautau_mc_2016", "fc2_BuKtautau_mc_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_1/bdt_output.txt")
    fc2_BuKtautau_mc_2017 = ROOT.TFileCollection("fc2_BuKtautau_mc_2017", "fc2_BuKtautau_mc_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_1/bdt_output.txt")
    fc2_BuKtautau_mc_2018 = ROOT.TFileCollection("fc2_BuKtautau_mc_2018", "fc2_BuKtautau_mc_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_1/bdt_output.txt")

    t2_BuKtautau_mc_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuKtautau_mc_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuKtautau_mc_2018 = ROOT.TChain("XGBoost/DecayTree")

    t2_BuKtautau_mc_2016.AddFileInfoList(fc2_BuKtautau_mc_2016.GetList())
    t2_BuKtautau_mc_2017.AddFileInfoList(fc2_BuKtautau_mc_2017.GetList())
    t2_BuKtautau_mc_2018.AddFileInfoList(fc2_BuKtautau_mc_2018.GetList())

    t2_BuKtautau_mc_2016.Add(t2_BuKtautau_mc_2017)
    t2_BuKtautau_mc_2016.Add(t2_BuKtautau_mc_2018)

    t_BuKtautau_mc_2016.AddFriend(t1_BuKtautau_mc_2016)
    t_BuKtautau_mc_2016.AddFriend(t2_BuKtautau_mc_2016)

    ######################################################################################################################################################
    
    ###################################################### Loose truth-matched cocktail MCs #####################################################################
    cut_cocktail_mc = "(df_status == 0) && (df_Bp_M > 4000) && (df_Bp_M < 8000) && (taup_BKGCAT <= 60) && (taum_BKGCAT <= 60) && (Bp_BKGCAT <= 60)"

    # pre-sel tree
    fc_BuDDKp_mc_2016 = ROOT.TFileCollection("fc_BuDDKp_mc_2016", "fc_BuDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_100/pre_sel_tree.txt")
    fc_BdDDKp_mc_2016 = ROOT.TFileCollection("fc_BdDDKp_mc_2016", "fc_BdDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_110/pre_sel_tree.txt")
    fc_BsDDKp_mc_2016 = ROOT.TFileCollection("fc_BsDDKp_mc_2016", "fc_BsDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_120/pre_sel_tree.txt")
    fc_BuDDK0_mc_2016 = ROOT.TFileCollection("fc_BuDDK0_mc_2016", "fc_BuDDK0_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_130/pre_sel_tree.txt")
    fc_BuDD_mc_2016 = ROOT.TFileCollection("fc_BuDD_mc_2016", "fc_BuDD_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_150/pre_sel_tree.txt")

    fc_BuDDKp_mc_2017 = ROOT.TFileCollection("fc_BuDDKp_mc_2017", "fc_BuDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_100/pre_sel_tree.txt")
    fc_BdDDKp_mc_2017 = ROOT.TFileCollection("fc_BdDDKp_mc_2017", "fc_BdDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_110/pre_sel_tree.txt")
    fc_BsDDKp_mc_2017 = ROOT.TFileCollection("fc_BsDDKp_mc_2017", "fc_BsDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_120/pre_sel_tree.txt")
    fc_BuDDK0_mc_2017 = ROOT.TFileCollection("fc_BuDDK0_mc_2017", "fc_BuDDK0_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_130/pre_sel_tree.txt")
    fc_BuDD_mc_2017 = ROOT.TFileCollection("fc_BuDD_mc_2017", "fc_BuDD_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_150/pre_sel_tree.txt")

    fc_BuDDKp_mc_2018 = ROOT.TFileCollection("fc_BuDDKp_mc_2018", "fc_BuDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_100/pre_sel_tree.txt")
    fc_BdDDKp_mc_2018 = ROOT.TFileCollection("fc_BdDDKp_mc_2018", "fc_BdDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_110/pre_sel_tree.txt")
    fc_BsDDKp_mc_2018 = ROOT.TFileCollection("fc_BsDDKp_mc_2018", "fc_BsDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_120/pre_sel_tree.txt")
    fc_BuDDK0_mc_2018 = ROOT.TFileCollection("fc_BuDDK0_mc_2018", "fc_BuDDK0_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_130/pre_sel_tree.txt")
    fc_BuDD_mc_2018 = ROOT.TFileCollection("fc_BuDD_mc_2018", "fc_BuDD_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_150/pre_sel_tree.txt")

    t_BuDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t_BdDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t_BsDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t_BuDDK0_mc_2016 = ROOT.TChain("DecayTree")
    t_BuDD_mc_2016 = ROOT.TChain("DecayTree")

    t_BuDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t_BdDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t_BsDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t_BuDDK0_mc_2017 = ROOT.TChain("DecayTree")
    t_BuDD_mc_2017 = ROOT.TChain("DecayTree")

    t_BuDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t_BdDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t_BsDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t_BuDDK0_mc_2018 = ROOT.TChain("DecayTree")
    t_BuDD_mc_2018 = ROOT.TChain("DecayTree")

    t_BuDDKp_mc_2016.AddFileInfoList(fc_BuDDKp_mc_2016.GetList())
    t_BdDDKp_mc_2016.AddFileInfoList(fc_BdDDKp_mc_2016.GetList())
    t_BsDDKp_mc_2016.AddFileInfoList(fc_BsDDKp_mc_2016.GetList())
    t_BuDDK0_mc_2016.AddFileInfoList(fc_BuDDK0_mc_2016.GetList())
    t_BuDD_mc_2016.AddFileInfoList(fc_BuDD_mc_2016.GetList())

    t_BuDDKp_mc_2017.AddFileInfoList(fc_BuDDKp_mc_2017.GetList())
    t_BdDDKp_mc_2017.AddFileInfoList(fc_BdDDKp_mc_2017.GetList())
    t_BsDDKp_mc_2017.AddFileInfoList(fc_BsDDKp_mc_2017.GetList())
    t_BuDDK0_mc_2017.AddFileInfoList(fc_BuDDK0_mc_2017.GetList())
    t_BuDD_mc_2017.AddFileInfoList(fc_BuDD_mc_2017.GetList())

    t_BuDDKp_mc_2018.AddFileInfoList(fc_BuDDKp_mc_2018.GetList())
    t_BdDDKp_mc_2018.AddFileInfoList(fc_BdDDKp_mc_2018.GetList())
    t_BsDDKp_mc_2018.AddFileInfoList(fc_BsDDKp_mc_2018.GetList())
    t_BuDDK0_mc_2018.AddFileInfoList(fc_BuDDK0_mc_2018.GetList())
    t_BuDD_mc_2018.AddFileInfoList(fc_BuDD_mc_2018.GetList())

    t_BuDDKp_mc_2016.Add(t_BuDDKp_mc_2017)
    t_BuDDKp_mc_2016.Add(t_BuDDKp_mc_2018)

    t_BdDDKp_mc_2016.Add(t_BdDDKp_mc_2017)
    t_BdDDKp_mc_2016.Add(t_BdDDKp_mc_2018)

    t_BsDDKp_mc_2016.Add(t_BsDDKp_mc_2017)
    t_BsDDKp_mc_2016.Add(t_BsDDKp_mc_2018)

    t_BuDDK0_mc_2016.Add(t_BuDDK0_mc_2017)
    t_BuDDK0_mc_2016.Add(t_BuDDK0_mc_2018)

    t_BuDD_mc_2016.Add(t_BuDD_mc_2017)
    t_BuDD_mc_2016.Add(t_BuDD_mc_2018)

    # gsl 
    fc1_BuDDKp_mc_2016 = ROOT.TFileCollection("fc1_BuDDKp_mc_2016", "fc1_BuDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_100/fit_results.txt")
    fc1_BdDDKp_mc_2016 = ROOT.TFileCollection("fc1_BdDDKp_mc_2016", "fc1_BdDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_110/fit_results.txt")
    fc1_BsDDKp_mc_2016 = ROOT.TFileCollection("fc1_BsDDKp_mc_2016", "fc1_BsDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_120/fit_results.txt")
    fc1_BuDDK0_mc_2016 = ROOT.TFileCollection("fc1_BuDDK0_mc_2016", "fc1_BuDDK0_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_130/fit_results.txt")
    fc1_BuDD_mc_2016 = ROOT.TFileCollection("fc1_BuDD_mc_2016", "fc1_BuDD_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_150/fit_results.txt")

    fc1_BuDDKp_mc_2017 = ROOT.TFileCollection("fc1_BuDDKp_mc_2017", "fc1_BuDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_100/fit_results.txt")
    fc1_BdDDKp_mc_2017 = ROOT.TFileCollection("fc1_BdDDKp_mc_2017", "fc1_BdDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_110/fit_results.txt")
    fc1_BsDDKp_mc_2017 = ROOT.TFileCollection("fc1_BsDDKp_mc_2017", "fc1_BsDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_120/fit_results.txt")
    fc1_BuDDK0_mc_2017 = ROOT.TFileCollection("fc1_BuDDK0_mc_2017", "fc1_BuDDK0_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_130/fit_results.txt")
    fc1_BuDD_mc_2017 = ROOT.TFileCollection("fc1_BuDD_mc_2017", "fc1_BuDD_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_150/fit_results.txt")

    fc1_BuDDKp_mc_2018 = ROOT.TFileCollection("fc1_BuDDKp_mc_2018", "fc1_BuDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_100/fit_results.txt")
    fc1_BdDDKp_mc_2018 = ROOT.TFileCollection("fc1_BdDDKp_mc_2018", "fc1_BdDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_110/fit_results.txt")
    fc1_BsDDKp_mc_2018 = ROOT.TFileCollection("fc1_BsDDKp_mc_2018", "fc1_BsDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_120/fit_results.txt")
    fc1_BuDDK0_mc_2018 = ROOT.TFileCollection("fc1_BuDDK0_mc_2018", "fc1_BuDDK0_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_130/fit_results.txt")
    fc1_BuDD_mc_2018 = ROOT.TFileCollection("fc1_BuDD_mc_2018", "fc1_BuDD_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_150/fit_results.txt")

    t1_BuDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t1_BdDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t1_BsDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t1_BuDDK0_mc_2016 = ROOT.TChain("DecayTree")
    t1_BuDD_mc_2016 = ROOT.TChain("DecayTree")

    t1_BuDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t1_BdDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t1_BsDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t1_BuDDK0_mc_2017 = ROOT.TChain("DecayTree")
    t1_BuDD_mc_2017 = ROOT.TChain("DecayTree")

    t1_BuDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t1_BdDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t1_BsDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t1_BuDDK0_mc_2018 = ROOT.TChain("DecayTree")
    t1_BuDD_mc_2018 = ROOT.TChain("DecayTree")

    t1_BuDDKp_mc_2016.AddFileInfoList(fc1_BuDDKp_mc_2016.GetList())
    t1_BdDDKp_mc_2016.AddFileInfoList(fc1_BdDDKp_mc_2016.GetList())
    t1_BsDDKp_mc_2016.AddFileInfoList(fc1_BsDDKp_mc_2016.GetList())
    t1_BuDDK0_mc_2016.AddFileInfoList(fc1_BuDDK0_mc_2016.GetList())
    t1_BuDD_mc_2016.AddFileInfoList(fc1_BuDD_mc_2016.GetList())

    t1_BuDDKp_mc_2017.AddFileInfoList(fc1_BuDDKp_mc_2017.GetList())
    t1_BdDDKp_mc_2017.AddFileInfoList(fc1_BdDDKp_mc_2017.GetList())
    t1_BsDDKp_mc_2017.AddFileInfoList(fc1_BsDDKp_mc_2017.GetList())
    t1_BuDDK0_mc_2017.AddFileInfoList(fc1_BuDDK0_mc_2017.GetList())
    t1_BuDD_mc_2017.AddFileInfoList(fc1_BuDD_mc_2017.GetList())

    t1_BuDDKp_mc_2018.AddFileInfoList(fc1_BuDDKp_mc_2018.GetList())
    t1_BdDDKp_mc_2018.AddFileInfoList(fc1_BdDDKp_mc_2018.GetList())
    t1_BsDDKp_mc_2018.AddFileInfoList(fc1_BsDDKp_mc_2018.GetList())
    t1_BuDDK0_mc_2018.AddFileInfoList(fc1_BuDDK0_mc_2018.GetList())
    t1_BuDD_mc_2018.AddFileInfoList(fc1_BuDD_mc_2018.GetList())

    t1_BuDDKp_mc_2016.Add(t1_BuDDKp_mc_2017)
    t1_BuDDKp_mc_2016.Add(t1_BuDDKp_mc_2018)

    t1_BdDDKp_mc_2016.Add(t1_BdDDKp_mc_2017)
    t1_BdDDKp_mc_2016.Add(t1_BdDDKp_mc_2018)

    t1_BsDDKp_mc_2016.Add(t1_BsDDKp_mc_2017)
    t1_BsDDKp_mc_2016.Add(t1_BsDDKp_mc_2018)

    t1_BuDDK0_mc_2016.Add(t1_BuDDK0_mc_2017)
    t1_BuDDK0_mc_2016.Add(t1_BuDDK0_mc_2018)

    t1_BuDD_mc_2016.Add(t1_BuDD_mc_2017)
    t1_BuDD_mc_2016.Add(t1_BuDD_mc_2018)

    # sklearn
    fc2_BuDDKp_mc_2016 = ROOT.TFileCollection("fc2_BuDDKp_mc_2016", "fc2_BuDDKp_mc_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_100/bdt_output.txt")
    fc2_BdDDKp_mc_2016 = ROOT.TFileCollection("fc2_BdDDKp_mc_2016", "fc2_BdDDKp_mc_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_110/bdt_output.txt")
    fc2_BsDDKp_mc_2016 = ROOT.TFileCollection("fc2_BsDDKp_mc_2016", "fc2_BsDDKp_mc_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_120/bdt_output.txt")
    fc2_BuDDK0_mc_2016 = ROOT.TFileCollection("fc2_BuDDK0_mc_2016", "fc2_BuDDK0_mc_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_130/bdt_output.txt")
    fc2_BuDD_mc_2016 = ROOT.TFileCollection("fc2_BuDD_mc_2016", "fc2_BuDD_mc_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_150/bdt_output.txt")

    fc2_BuDDKp_mc_2017 = ROOT.TFileCollection("fc2_BuDDKp_mc_2017", "fc2_BuDDKp_mc_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_100/bdt_output.txt")
    fc2_BdDDKp_mc_2017 = ROOT.TFileCollection("fc2_BdDDKp_mc_2017", "fc2_BdDDKp_mc_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_110/bdt_output.txt")
    fc2_BsDDKp_mc_2017 = ROOT.TFileCollection("fc2_BsDDKp_mc_2017", "fc2_BsDDKp_mc_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_120/bdt_output.txt")
    fc2_BuDDK0_mc_2017 = ROOT.TFileCollection("fc2_BuDDK0_mc_2017", "fc2_BuDDK0_mc_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_130/bdt_output.txt")
    fc2_BuDD_mc_2017 = ROOT.TFileCollection("fc2_BuDD_mc_2017", "fc2_BuDD_mc_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_150/bdt_output.txt")

    fc2_BuDDKp_mc_2018 = ROOT.TFileCollection("fc2_BuDDKp_mc_2018", "fc2_BuDDKp_mc_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_100/bdt_output.txt")
    fc2_BdDDKp_mc_2018 = ROOT.TFileCollection("fc2_BdDDKp_mc_2018", "fc2_BdDDKp_mc_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_110/bdt_output.txt")
    fc2_BsDDKp_mc_2018 = ROOT.TFileCollection("fc2_BsDDKp_mc_2018", "fc2_BsDDKp_mc_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_120/bdt_output.txt")
    fc2_BuDDK0_mc_2018 = ROOT.TFileCollection("fc2_BuDDK0_mc_2018", "fc2_BuDDK0_mc_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_130/bdt_output.txt")
    fc2_BuDD_mc_2018 = ROOT.TFileCollection("fc2_BuDD_mc_2018", "fc2_BuDD_mc_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_150/bdt_output.txt")

    t2_BuDDKp_mc_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_BdDDKp_mc_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_BsDDKp_mc_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuDDK0_mc_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuDD_mc_2016 = ROOT.TChain("XGBoost/DecayTree")

    t2_BuDDKp_mc_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_BdDDKp_mc_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_BsDDKp_mc_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuDDK0_mc_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuDD_mc_2017 = ROOT.TChain("XGBoost/DecayTree")

    t2_BuDDKp_mc_2018 = ROOT.TChain("XGBoost/DecayTree")
    t2_BdDDKp_mc_2018 = ROOT.TChain("XGBoost/DecayTree")
    t2_BsDDKp_mc_2018 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuDDK0_mc_2018 = ROOT.TChain("XGBoost/DecayTree")
    t2_BuDD_mc_2018 = ROOT.TChain("XGBoost/DecayTree")

    t2_BuDDKp_mc_2016.AddFileInfoList(fc2_BuDDKp_mc_2016.GetList())
    t2_BdDDKp_mc_2016.AddFileInfoList(fc2_BdDDKp_mc_2016.GetList())
    t2_BsDDKp_mc_2016.AddFileInfoList(fc2_BsDDKp_mc_2016.GetList())
    t2_BuDDK0_mc_2016.AddFileInfoList(fc2_BuDDK0_mc_2016.GetList())
    t2_BuDD_mc_2016.AddFileInfoList(fc2_BuDD_mc_2016.GetList())

    t2_BuDDKp_mc_2017.AddFileInfoList(fc2_BuDDKp_mc_2017.GetList())
    t2_BdDDKp_mc_2017.AddFileInfoList(fc2_BdDDKp_mc_2017.GetList())
    t2_BsDDKp_mc_2017.AddFileInfoList(fc2_BsDDKp_mc_2017.GetList())
    t2_BuDDK0_mc_2017.AddFileInfoList(fc2_BuDDK0_mc_2017.GetList())
    t2_BuDD_mc_2017.AddFileInfoList(fc2_BuDD_mc_2017.GetList())

    t2_BuDDKp_mc_2018.AddFileInfoList(fc2_BuDDKp_mc_2018.GetList())
    t2_BdDDKp_mc_2018.AddFileInfoList(fc2_BdDDKp_mc_2018.GetList())
    t2_BsDDKp_mc_2018.AddFileInfoList(fc2_BsDDKp_mc_2018.GetList())
    t2_BuDDK0_mc_2018.AddFileInfoList(fc2_BuDDK0_mc_2018.GetList())
    t2_BuDD_mc_2018.AddFileInfoList(fc2_BuDD_mc_2018.GetList())

    t2_BuDDKp_mc_2016.Add(t2_BuDDKp_mc_2017)
    t2_BuDDKp_mc_2016.Add(t2_BuDDKp_mc_2018)

    t2_BdDDKp_mc_2016.Add(t2_BdDDKp_mc_2017)
    t2_BdDDKp_mc_2016.Add(t2_BdDDKp_mc_2018)

    t2_BsDDKp_mc_2016.Add(t2_BsDDKp_mc_2017)
    t2_BsDDKp_mc_2016.Add(t2_BsDDKp_mc_2018)

    t2_BuDDK0_mc_2016.Add(t2_BuDDK0_mc_2017)
    t2_BuDDK0_mc_2016.Add(t2_BuDDK0_mc_2018)

    t2_BuDD_mc_2016.Add(t2_BuDD_mc_2017)
    t2_BuDD_mc_2016.Add(t2_BuDD_mc_2018)


    t_BuDDKp_mc_2016.AddFriend(t1_BuDDKp_mc_2016)
    t_BuDDKp_mc_2016.AddFriend(t2_BuDDKp_mc_2016)

    t_BdDDKp_mc_2016.AddFriend(t1_BdDDKp_mc_2016)
    t_BdDDKp_mc_2016.AddFriend(t2_BdDDKp_mc_2016)

    t_BsDDKp_mc_2016.AddFriend(t1_BsDDKp_mc_2016)
    t_BsDDKp_mc_2016.AddFriend(t2_BsDDKp_mc_2016)

    t_BuDDK0_mc_2016.AddFriend(t1_BuDDK0_mc_2016)
    t_BuDDK0_mc_2016.AddFriend(t2_BuDDK0_mc_2016)

    t_BuDD_mc_2016.AddFriend(t1_BuDD_mc_2016)
    t_BuDD_mc_2016.AddFriend(t2_BuDD_mc_2016)
    ######################################################################################################################################################

    ###################################################### RS and WS data #####################################################################    
    # pre-sel tree
    fc_rs_data_2016 = ROOT.TFileCollection("fc_rs_data_2016", "fc_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_2/pre_sel_tree.txt", 10)
    fc_rs_data_2017 = ROOT.TFileCollection("fc_rs_data_2017", "fc_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_2/pre_sel_tree.txt", 10)
    fc_rs_data_2018 = ROOT.TFileCollection("fc_rs_data_2018", "fc_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_2/pre_sel_tree.txt", 10)

    fc_ws_data_2016 = ROOT.TFileCollection("fc_ws_data_2016", "fc_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 10)
    fc_ws_data_2017 = ROOT.TFileCollection("fc_ws_data_2017", "fc_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 10)
    fc_ws_data_2018 = ROOT.TFileCollection("fc_ws_data_2018", "fc_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 10)

    t_rs_data_2016 = ROOT.TChain("DecayTree")
    t_rs_data_2017 = ROOT.TChain("DecayTree")
    t_rs_data_2018 = ROOT.TChain("DecayTree")

    t_ws_data_2016 = ROOT.TChain("DecayTree")
    t_ws_data_2017 = ROOT.TChain("DecayTree")
    t_ws_data_2018 = ROOT.TChain("DecayTree")

    t_rs_data_2016.AddFileInfoList(fc_rs_data_2016.GetList())
    t_rs_data_2017.AddFileInfoList(fc_rs_data_2017.GetList())
    t_rs_data_2018.AddFileInfoList(fc_rs_data_2018.GetList())

    t_ws_data_2016.AddFileInfoList(fc_ws_data_2016.GetList())
    t_ws_data_2017.AddFileInfoList(fc_ws_data_2017.GetList())
    t_ws_data_2018.AddFileInfoList(fc_ws_data_2018.GetList())

    t_rs_data_2016.Add(t_rs_data_2017)
    t_rs_data_2016.Add(t_rs_data_2018)

    t_ws_data_2016.Add(t_ws_data_2017)
    t_ws_data_2016.Add(t_ws_data_2018)

    # gsl 
    fc1_rs_data_2016 = ROOT.TFileCollection("fc1_rs_data_2016", "fc1_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt", 10)
    fc1_rs_data_2017 = ROOT.TFileCollection("fc1_rs_data_2017", "fc1_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt", 10)
    fc1_rs_data_2018 = ROOT.TFileCollection("fc1_rs_data_2018", "fc1_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt", 10)

    fc1_ws_data_2016 = ROOT.TFileCollection("fc1_ws_data_2016", "fc1_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 10)
    fc1_ws_data_2017 = ROOT.TFileCollection("fc1_ws_data_2017", "fc1_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 10)
    fc1_ws_data_2018 = ROOT.TFileCollection("fc1_ws_data_2018", "fc1_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 10)

    t1_rs_data_2016 = ROOT.TChain("DecayTree")
    t1_rs_data_2017 = ROOT.TChain("DecayTree")
    t1_rs_data_2018 = ROOT.TChain("DecayTree")

    t1_ws_data_2016 = ROOT.TChain("DecayTree")
    t1_ws_data_2017 = ROOT.TChain("DecayTree")
    t1_ws_data_2018 = ROOT.TChain("DecayTree")

    t1_rs_data_2016.AddFileInfoList(fc1_rs_data_2016.GetList())
    t1_rs_data_2017.AddFileInfoList(fc1_rs_data_2017.GetList())
    t1_rs_data_2018.AddFileInfoList(fc1_rs_data_2018.GetList())

    t1_ws_data_2016.AddFileInfoList(fc1_ws_data_2016.GetList())
    t1_ws_data_2017.AddFileInfoList(fc1_ws_data_2017.GetList())
    t1_ws_data_2018.AddFileInfoList(fc1_ws_data_2018.GetList())

    t1_rs_data_2016.Add(t1_rs_data_2017)
    t1_rs_data_2016.Add(t1_rs_data_2018)

    t1_ws_data_2016.Add(t1_ws_data_2017)
    t1_ws_data_2016.Add(t1_ws_data_2018)

    # sklearn
    fc2_rs_data_2016 = ROOT.TFileCollection("fc2_rs_data_2016", "fc2_rs_data_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_2/bdt_output.txt", 10)
    fc2_rs_data_2017 = ROOT.TFileCollection("fc2_rs_data_2017", "fc2_rs_data_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_2/bdt_output.txt", 10)
    fc2_rs_data_2018 = ROOT.TFileCollection("fc2_rs_data_2018", "fc2_rs_data_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_2/bdt_output.txt", 10)

    fc2_ws_data_2016 = ROOT.TFileCollection("fc2_ws_data_2016", "fc2_ws_data_2016", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_3/bdt_output.txt", 10)
    fc2_ws_data_2017 = ROOT.TFileCollection("fc2_ws_data_2017", "fc2_ws_data_2017", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_3/bdt_output.txt", 10)
    fc2_ws_data_2018 = ROOT.TFileCollection("fc2_ws_data_2018", "fc2_ws_data_2018", f"/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_3/bdt_output.txt", 10)

    t2_rs_data_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_rs_data_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_rs_data_2018 = ROOT.TChain("XGBoost/DecayTree")

    t2_ws_data_2016 = ROOT.TChain("XGBoost/DecayTree")
    t2_ws_data_2017 = ROOT.TChain("XGBoost/DecayTree")
    t2_ws_data_2018 = ROOT.TChain("XGBoost/DecayTree")

    t2_rs_data_2016.AddFileInfoList(fc2_rs_data_2016.GetList())
    t2_rs_data_2017.AddFileInfoList(fc2_rs_data_2017.GetList())
    t2_rs_data_2018.AddFileInfoList(fc2_rs_data_2018.GetList())

    t2_ws_data_2016.AddFileInfoList(fc2_ws_data_2016.GetList())
    t2_ws_data_2017.AddFileInfoList(fc2_ws_data_2017.GetList())
    t2_ws_data_2018.AddFileInfoList(fc2_ws_data_2018.GetList())

    t2_rs_data_2016.Add(t2_rs_data_2017)
    t2_rs_data_2016.Add(t2_rs_data_2018)

    t2_ws_data_2016.Add(t2_ws_data_2017)
    t2_ws_data_2016.Add(t2_ws_data_2018)


    t_rs_data_2016.AddFriend(t1_rs_data_2016)
    t_rs_data_2016.AddFriend(t2_rs_data_2016)

    t_ws_data_2016.AddFriend(t1_ws_data_2016)
    t_ws_data_2016.AddFriend(t2_ws_data_2016)
    ######################################################################################################################################################

    ### PLOTS: 3 MC components
    # Physics BDT
    h_3pi3pi_mc_phys = ROOT.TH1D("h_3pi3pi_mc_phys", "h_3pi3pi_mc_phys", 100, 0, 1)
    h_3pi3pipi0_mc_phys = ROOT.TH1D("h_3pi3pipi0_mc_phys", "h_3pi3pipi0_mc_phys", 100, 0, 1)
    h_3pi3pi2pi0_mc_phys = ROOT.TH1D("h_3pi3pi2pi0_mc_phys", "h_3pi3pi2pi0_mc_phys", 100, 0, 1)

    t_BuKtautau_mc_2016.Draw("BDT1 >> h_3pi3pi_mc_phys", cut_Ktautau+" && (component==0)")
    t_BuKtautau_mc_2016.Draw("BDT1 >> h_3pi3pipi0_mc_phys", cut_Ktautau+" && (component==1)")
    t_BuKtautau_mc_2016.Draw("BDT1 >> h_3pi3pi2pi0_mc_phys", cut_Ktautau+" && (component==2)")

    c = ROOT.TCanvas()
    c.cd()
    h_3pi3pi_mc_phys.GetXaxis().SetTitle("BDT1")
    h_3pi3pi_mc_phys.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_3pi3pi_mc_phys.SetTitle("Physics BDT")
    h_3pi3pi_mc_phys.SetLineColor(4)
    h_3pi3pipi0_mc_phys.SetLineColor(8)
    h_3pi3pi2pi0_mc_phys.SetLineColor(2)
    h_3pi3pi_mc_phys.DrawNormalized()
    h_3pi3pipi0_mc_phys.DrawNormalized("same")
    h_3pi3pi2pi0_mc_phys.DrawNormalized("same")
    phys_leg = ROOT.TLegend(0.1, 0.7, 0.4, 0.89)
    phys_leg.AddEntry(h_3pi3pi_mc_phys, "3#pi3#pi MC", "lp")
    phys_leg.AddEntry(h_3pi3pipi0_mc_phys, "3#pi3#pi #pi^{0} MC", "lp")
    phys_leg.AddEntry(h_3pi3pi2pi0_mc_phys, "3#pi3#pi 2#pi^{0} MC", "lp")
    phys_leg.SetBorderSize(0)
    phys_leg.Draw("same")
    c.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/phys_ktautau_mc_components.pdf")

    # Combinatorial BDT
    h_3pi3pi_mc_comb = ROOT.TH1D("h_3pi3pi_mc_comb", "h_3pi3pi_mc_comb", 100, 0, 1)
    h_3pi3pipi0_mc_comb = ROOT.TH1D("h_3pi3pipi0_mc_comb", "h_3pi3pipi0_mc_comb", 100, 0, 1)
    h_3pi3pi2pi0_mc_comb = ROOT.TH1D("h_3pi3pi2pi0_mc_comb", "h_3pi3pi2pi0_mc_comb", 100, 0, 1)

    t_BuKtautau_mc_2016.Draw("BDT2 >> h_3pi3pi_mc_comb", cut_Ktautau+" && (component==0)")
    t_BuKtautau_mc_2016.Draw("BDT2 >> h_3pi3pipi0_mc_comb", cut_Ktautau+" && (component==1)")
    t_BuKtautau_mc_2016.Draw("BDT2 >> h_3pi3pi2pi0_mc_comb", cut_Ktautau+" && (component==2)")

    c1 = ROOT.TCanvas()
    c1.cd()
    h_3pi3pi_mc_comb.GetXaxis().SetTitle("BDT2")
    h_3pi3pi_mc_comb.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_3pi3pi_mc_comb.SetTitle("Combinatorial BDT")
    h_3pi3pi_mc_comb.SetLineColor(4)
    h_3pi3pipi0_mc_comb.SetLineColor(8)
    h_3pi3pi2pi0_mc_comb.SetLineColor(2)
    h_3pi3pi_mc_comb.DrawNormalized()
    h_3pi3pipi0_mc_comb.DrawNormalized("same")
    h_3pi3pi2pi0_mc_comb.DrawNormalized("same")
    comb_leg = ROOT.TLegend(0.1, 0.7, 0.4, 0.89)
    comb_leg.AddEntry(h_3pi3pi_mc_comb, "3#pi3#pi MC", "lp")
    comb_leg.AddEntry(h_3pi3pipi0_mc_comb, "3#pi3#pi #pi^{0} MC", "lp")
    comb_leg.AddEntry(h_3pi3pi2pi0_mc_comb, "3#pi3#pi 2#pi^{0} MC", "lp")
    comb_leg.SetBorderSize(0)
    comb_leg.Draw("same")
    c1.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/comb_ktautau_mc_components.pdf")

    ### PLOTS: cocktail MCs
    # Physics BDT
    h_BuDDKp_mc_phys = ROOT.TH1D("h_BuDDKp_mc_phys", "h_BuDDKp_mc_phys", 100, 0, 1)
    h_BdDDKp_mc_phys = ROOT.TH1D("h_BdDDKp_mc_phys", "h_BdDDKp_mc_phys", 100, 0, 1)
    h_BsDDKp_mc_phys = ROOT.TH1D("h_BsDDKp_mc_phys", "h_BsDDKp_mc_phys", 100, 0, 1)
    h_BuDDK0_mc_phys = ROOT.TH1D("h_BuDDK0_mc_phys", "h_BuDDK0_mc_phys", 100, 0, 1)
    h_BuDD_mc_phys = ROOT.TH1D("h_BuDD_mc_phys", "h_BuDD_mc_phys", 100, 0, 1)

    t_BuDDKp_mc_2016.Draw("BDT1 >> h_BuDDKp_mc_phys", cut_cocktail_mc)
    t_BdDDKp_mc_2016.Draw("BDT1 >> h_BdDDKp_mc_phys", cut_cocktail_mc)
    t_BsDDKp_mc_2016.Draw("BDT1 >> h_BsDDKp_mc_phys", cut_cocktail_mc)
    t_BuDDK0_mc_2016.Draw("BDT1 >> h_BuDDK0_mc_phys", cut_cocktail_mc)
    t_BuDD_mc_2016.Draw("BDT1 >> h_BuDD_mc_phys", cut_cocktail_mc)

    c2 = ROOT.TCanvas()
    c2.cd()
    h_BsDDKp_mc_phys.GetXaxis().SetTitle("BDT1")
    h_BsDDKp_mc_phys.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_BsDDKp_mc_phys.SetTitle("Physics BDT")
    h_BuDDKp_mc_phys.SetLineColor(4)
    h_BdDDKp_mc_phys.SetLineColor(8)
    h_BsDDKp_mc_phys.SetLineColor(2)
    h_BuDDK0_mc_phys.SetLineColor(51)
    h_BuDD_mc_phys.SetLineColor(94)
    h_BsDDKp_mc_phys.DrawNormalized()
    h_BdDDKp_mc_phys.DrawNormalized("same")
    h_BuDDKp_mc_phys.DrawNormalized("same")
    h_BuDDK0_mc_phys.DrawNormalized("same")
    h_BuDD_mc_phys.DrawNormalized("same")
    phys_leg_1 = ROOT.TLegend(0.4, 0.6, 0.7, 0.89)
    phys_leg_1.AddEntry(h_BuDDKp_mc_phys, "B^{+} #rightarrow D D K^{+}", "lp")
    phys_leg_1.AddEntry(h_BdDDKp_mc_phys, "B^{0} #rightarrow D D K^{+}", "lp")
    phys_leg_1.AddEntry(h_BsDDKp_mc_phys, "B^{0}_{s} #rightarrow D D K^{+}", "lp")
    phys_leg_1.AddEntry(h_BuDDK0_mc_phys, "B^{+} #rightarrow D D K^{0}", "lp")
    phys_leg_1.AddEntry(h_BuDD_mc_phys, "B^{+} #rightarrow D D ", "lp")
    phys_leg_1.SetBorderSize(0)
    phys_leg_1.Draw("same")
    c2.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/phys_cocktail_mcs.pdf")

    # Combinatorial BDT
    h_BuDDKp_mc_comb = ROOT.TH1D("h_BuDDKp_mc_comb", "h_BuDDKp_mc_comb", 100, 0, 1)
    h_BdDDKp_mc_comb = ROOT.TH1D("h_BdDDKp_mc_comb", "h_BdDDKp_mc_comb", 100, 0, 1)
    h_BsDDKp_mc_comb = ROOT.TH1D("h_BsDDKp_mc_comb", "h_BsDDKp_mc_comb", 100, 0, 1)
    h_BuDDK0_mc_comb = ROOT.TH1D("h_BuDDK0_mc_comb", "h_BuDDK0_mc_comb", 100, 0, 1)
    h_BuDD_mc_comb = ROOT.TH1D("h_BuDD_mc_comb", "h_BuDD_mc_comb", 100, 0, 1)

    t_BuDDKp_mc_2016.Draw("BDT2 >> h_BuDDKp_mc_comb", cut_cocktail_mc)
    t_BdDDKp_mc_2016.Draw("BDT2 >> h_BdDDKp_mc_comb", cut_cocktail_mc)
    t_BsDDKp_mc_2016.Draw("BDT2 >> h_BsDDKp_mc_comb", cut_cocktail_mc)
    t_BuDDK0_mc_2016.Draw("BDT2 >> h_BuDDK0_mc_comb", cut_cocktail_mc)
    t_BuDD_mc_2016.Draw("BDT2 >> h_BuDD_mc_comb", cut_cocktail_mc)

    c3 = ROOT.TCanvas()
    c3.cd()
    h_BuDD_mc_comb.GetXaxis().SetTitle("BDT2")
    h_BuDD_mc_comb.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_BuDD_mc_comb.SetTitle("Combinatorial BDT")
    h_BuDDKp_mc_comb.SetLineColor(4)
    h_BdDDKp_mc_comb.SetLineColor(8)
    h_BsDDKp_mc_comb.SetLineColor(2)
    h_BuDDK0_mc_comb.SetLineColor(51)
    h_BuDD_mc_comb.SetLineColor(94)
    h_BuDD_mc_comb.DrawNormalized()
    h_BdDDKp_mc_comb.DrawNormalized("same")
    h_BsDDKp_mc_comb.DrawNormalized("same")
    h_BuDDK0_mc_comb.DrawNormalized("same")
    h_BuDDKp_mc_comb.DrawNormalized("same")
    comb_leg_1 = ROOT.TLegend(0.4, 0.6, 0.7, 0.89)
    comb_leg_1.AddEntry(h_BuDDKp_mc_comb, "B^{+} #rightarrow D D K^{+}", "lp")
    comb_leg_1.AddEntry(h_BdDDKp_mc_comb, "B^{0} #rightarrow D D K^{+}", "lp")
    comb_leg_1.AddEntry(h_BsDDKp_mc_comb, "B^{0}_{s} #rightarrow D D K^{+}", "lp")
    comb_leg_1.AddEntry(h_BuDDK0_mc_comb, "B^{+} #rightarrow D D K^{0}", "lp")
    comb_leg_1.AddEntry(h_BuDD_mc_comb, "B^{+} #rightarrow D D ", "lp")
    comb_leg_1.SetBorderSize(0)
    comb_leg_1.Draw("same")
    c3.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/comb_cocktail_mcs.pdf")

    # PLOTS: RS and WS data
    # Physics BDT
    h_rs_data_phys = ROOT.TH1D("h_rs_data_phys", "h_rs_data_phys", 100, 0, 1)
    h_ws_data_phys = ROOT.TH1D("h_ws_data_phys", "h_ws_data_phys", 100, 0, 1)

    t_rs_data_2016.Draw("BDT1 >> h_rs_data_phys", cut_Ktautau)
    t_ws_data_2016.Draw("BDT1 >> h_ws_data_phys", cut_Ktautau)

    c4 = ROOT.TCanvas()
    c4.cd()
    h_ws_data_phys.GetXaxis().SetTitle("BDT1")
    h_ws_data_phys.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_ws_data_phys.SetTitle("Physics BDT")
    h_rs_data_phys.SetLineColor(4)
    h_ws_data_phys.SetLineColor(8)
    h_ws_data_phys.DrawNormalized()
    h_rs_data_phys.DrawNormalized("same")
    phys_leg_2 = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
    phys_leg_2.AddEntry(h_rs_data_phys, "RS data", "lp")
    phys_leg_2.AddEntry(h_ws_data_phys, "WS data", "lp")
    phys_leg_2.SetBorderSize(0)
    phys_leg_2.Draw("same")
    c4.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/phys_ktautau_data.pdf")

    # Combinatorial BDT
    h_rs_data_comb = ROOT.TH1D("h_rs_data_comb", "h_rs_data_comb", 100, 0, 1)
    h_ws_data_comb = ROOT.TH1D("h_ws_data_comb", "h_ws_data_comb", 100, 0, 1)

    t_rs_data_2016.Draw("BDT2 >> h_rs_data_comb", cut_Ktautau)
    t_ws_data_2016.Draw("BDT2 >> h_ws_data_comb", cut_Ktautau)

    c5 = ROOT.TCanvas()
    c5.cd()
    h_ws_data_comb.GetXaxis().SetTitle("BDT2")
    h_ws_data_comb.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_ws_data_comb.SetTitle("Combinatorial BDT")
    h_rs_data_comb.SetLineColor(4)
    h_ws_data_comb.SetLineColor(8)
    h_ws_data_comb.DrawNormalized()
    h_rs_data_comb.DrawNormalized("same")
    comb_leg_2 = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
    comb_leg_2.AddEntry(h_rs_data_comb, "RS data", "lp")
    comb_leg_2.AddEntry(h_ws_data_comb, "WS data", "lp")
    comb_leg_2.SetBorderSize(0)
    comb_leg_2.Draw("same")
    c5.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/comb_ktautau_data.pdf")

    # PLOTS: Signal vs background
    # Physics BDT
    h_ktautau_mc_phys = ROOT.TH1D("h_ktautau_mc_phys", "h_ktautau_mc_phys", 100, 0, 1)
    t_BuKtautau_mc_2016.Draw("BDT1 >> h_ktautau_mc_phys", cut_Ktautau)

    h_cocktail_mc_phys = ROOT.TH1D("h_cocktail_mc_phys", "h_cocktail_mc_phys", 100, 0, 1)
    h_cocktail_mc_phys.Sumw2()
    h_cocktail_mc_phys.Add(h_BuDDKp_mc_phys)
    h_cocktail_mc_phys.Add(h_BdDDKp_mc_phys)
    h_cocktail_mc_phys.Add(h_BsDDKp_mc_phys)
    h_cocktail_mc_phys.Add(h_BuDDK0_mc_phys)
    h_cocktail_mc_phys.Add(h_BuDD_mc_phys)

    c5 = ROOT.TCanvas()
    c5.cd()
    h_cocktail_mc_phys.GetXaxis().SetTitle("BDT1")
    h_cocktail_mc_phys.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_cocktail_mc_phys.SetTitle("Physics BDT")
    h_ktautau_mc_phys.SetLineColor(4)
    h_cocktail_mc_phys.SetLineColor(2)
    h_cocktail_mc_phys.DrawNormalized("hist")
    h_ktautau_mc_phys.DrawNormalized("same")
    phys_leg_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    phys_leg_3.AddEntry(h_ktautau_mc_phys, "Signal", "lp")
    phys_leg_3.AddEntry(h_cocktail_mc_phys, "Background", "lp")
    phys_leg_3.SetBorderSize(0)
    phys_leg_3.SetBorderSize(0)
    phys_leg_3.Draw("same")
    c5.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/phys_sig_vs_bkg.pdf")

    # Combinatorial BDT
    h_ktautau_mc_comb = ROOT.TH1D("h_ktautau_mc_comb", "h_ktautau_mc_comb", 100, 0, 1)
    t_BuKtautau_mc_2016.Draw("BDT2 >> h_ktautau_mc_comb", cut_Ktautau)

    c6 = ROOT.TCanvas()
    c6.cd()
    h_ws_data_comb.GetXaxis().SetTitle("BDT2")
    h_ws_data_comb.GetYaxis().SetTitle("Normalized entries / (0.01)")
    h_ws_data_comb.SetTitle("Combinatorial BDT")
    h_ktautau_mc_comb.SetLineColor(4)
    h_ws_data_comb.SetLineColor(2)
    h_ws_data_comb.DrawNormalized()
    h_ktautau_mc_comb.DrawNormalized("same")
    comb_leg_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    comb_leg_3.AddEntry(h_ktautau_mc_comb, "Signal", "lp")
    comb_leg_3.AddEntry(h_ws_data_comb, "Background", "lp")
    comb_leg_3.SetBorderSize(0)
    comb_leg_3.Draw("same")
    c6.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/comb_sig_vs_bkg.pdf")

    ### B+ mass evolution w/ BDT cut


if __name__ == "__main__":
    main(sys.argv)