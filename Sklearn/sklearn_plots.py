import sys
import ROOT


def main(argv):
    name = "response_merged"

    ROOT.gStyle.SetOptStat(0)

    ###################################################### Tight truth-matched Ktautau MC #####################################################################    
    t_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_1/selections.root")
    ktautau_mc_truthMatch_exp = t_sig.Get("truthMatch")
    trigger_exp = t_sig.Get("trigger")
    ktautau_mc_rectangular_cuts_exp = t_sig.Get("rectangular_cuts")
    pass_mass_fit_exp = t_sig.Get("pass_mass_fit")
    fit_range_exp = t_sig.Get("fit_range")
    mass_vetoes_exp = t_sig.Get("mass_vetoes")

    ktautau_mc_truthMatch = ROOT.TCut(ktautau_mc_truthMatch_exp.GetString().Data())
    trigger = ROOT.TCut(trigger_exp.GetString().Data())
    ktautau_mc_rectangular_cuts = ROOT.TCut(ktautau_mc_rectangular_cuts_exp.GetString().Data())
    pass_mass_fit = ROOT.TCut(pass_mass_fit_exp.GetString().Data())
    fit_range = ROOT.TCut(fit_range_exp.GetString().Data())
    mass_vetoes = ROOT.TCut(mass_vetoes_exp.GetString().Data())

    t_bkg = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_3/selections.root")
    ktautau_data_rectangular_cuts_exp = t_bkg.Get("rectangular_cuts")
    ktautau_data_rectangular_cuts = ROOT.TCut(ktautau_data_rectangular_cuts_exp.GetString().Data())
    
    cut_Ktautau_mc = (ktautau_mc_truthMatch+trigger+ktautau_mc_rectangular_cuts+pass_mass_fit+fit_range+mass_vetoes).GetTitle()
    cut_Ktautau_data = (trigger+ktautau_data_rectangular_cuts+pass_mass_fit+fit_range+mass_vetoes).GetTitle()

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
    fc2_BuKtautau_mc_2016 = ROOT.TFileCollection("fc2_BuKtautau_mc_2016", "fc2_BuKtautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_1/bdt_output.txt")
    fc2_BuKtautau_mc_2017 = ROOT.TFileCollection("fc2_BuKtautau_mc_2017", "fc2_BuKtautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_1/bdt_output.txt")
    fc2_BuKtautau_mc_2018 = ROOT.TFileCollection("fc2_BuKtautau_mc_2018", "fc2_BuKtautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_1/bdt_output.txt")

    t2_BuKtautau_mc_2016 = ROOT.TChain("DecayTree")
    t2_BuKtautau_mc_2017 = ROOT.TChain("DecayTree")
    t2_BuKtautau_mc_2018 = ROOT.TChain("DecayTree")

    t2_BuKtautau_mc_2016.AddFileInfoList(fc2_BuKtautau_mc_2016.GetList())
    t2_BuKtautau_mc_2017.AddFileInfoList(fc2_BuKtautau_mc_2017.GetList())
    t2_BuKtautau_mc_2018.AddFileInfoList(fc2_BuKtautau_mc_2018.GetList())

    t2_BuKtautau_mc_2016.Add(t2_BuKtautau_mc_2017)
    t2_BuKtautau_mc_2016.Add(t2_BuKtautau_mc_2018)

    # mass tree
    fc3_BuKtautau_mc_2016 = ROOT.TFileCollection("fc3_BuKtautau_mc_2016", "fc3_BuKtautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1/invariant_mass_tree.txt")
    fc3_BuKtautau_mc_2017 = ROOT.TFileCollection("fc3_BuKtautau_mc_2017", "fc3_BuKtautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1/invariant_mass_tree.txt")
    fc3_BuKtautau_mc_2018 = ROOT.TFileCollection("fc3_BuKtautau_mc_2018", "fc3_BuKtautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1/invariant_mass_tree.txt")

    t3_BuKtautau_mc_2016 = ROOT.TChain("DecayTree")
    t3_BuKtautau_mc_2017 = ROOT.TChain("DecayTree")
    t3_BuKtautau_mc_2018 = ROOT.TChain("DecayTree")

    t3_BuKtautau_mc_2016.AddFileInfoList(fc3_BuKtautau_mc_2016.GetList())
    t3_BuKtautau_mc_2017.AddFileInfoList(fc3_BuKtautau_mc_2017.GetList())
    t3_BuKtautau_mc_2018.AddFileInfoList(fc3_BuKtautau_mc_2018.GetList())

    t3_BuKtautau_mc_2016.Add(t3_BuKtautau_mc_2017)
    t3_BuKtautau_mc_2016.Add(t3_BuKtautau_mc_2018)

    if(name == "response_merged"):
        # sklearn merged
        fc4_BuKtautau_mc_2016 = ROOT.TFileCollection("fc4_BuKtautau_mc_2016", "fc4_BuKtautau_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_1/bdt_output.txt")
        fc4_BuKtautau_mc_2017 = ROOT.TFileCollection("fc4_BuKtautau_mc_2017", "fc4_BuKtautau_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_1/bdt_output.txt")
        fc4_BuKtautau_mc_2018 = ROOT.TFileCollection("fc4_BuKtautau_mc_2018", "fc4_BuKtautau_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_1/bdt_output.txt")

        t4_BuKtautau_mc_2016 = ROOT.TChain("DecayTree")
        t4_BuKtautau_mc_2017 = ROOT.TChain("DecayTree")
        t4_BuKtautau_mc_2018 = ROOT.TChain("DecayTree")

        t4_BuKtautau_mc_2016.AddFileInfoList(fc4_BuKtautau_mc_2016.GetList())
        t4_BuKtautau_mc_2017.AddFileInfoList(fc4_BuKtautau_mc_2017.GetList())
        t4_BuKtautau_mc_2018.AddFileInfoList(fc4_BuKtautau_mc_2018.GetList())

        t4_BuKtautau_mc_2016.Add(t4_BuKtautau_mc_2017)
        t4_BuKtautau_mc_2016.Add(t4_BuKtautau_mc_2018)

    t_BuKtautau_mc_2016.AddFriend(t1_BuKtautau_mc_2016)
    t_BuKtautau_mc_2016.AddFriend(t2_BuKtautau_mc_2016)
    t_BuKtautau_mc_2016.AddFriend(t3_BuKtautau_mc_2016)
    if(name == "response_merged"):
        t_BuKtautau_mc_2016.AddFriend(t4_BuKtautau_mc_2016)


    ######################################################################################################################################################
    
    ###################################################### Loose truth-matched cocktail MCs #####################################################################
    cut_cocktail_mc = (trigger+ktautau_mc_rectangular_cuts+pass_mass_fit+fit_range+mass_vetoes).GetTitle()

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
    fc2_BuDDKp_mc_2016 = ROOT.TFileCollection("fc2_BuDDKp_mc_2016", "fc2_BuDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_100/bdt_output.txt")
    fc2_BdDDKp_mc_2016 = ROOT.TFileCollection("fc2_BdDDKp_mc_2016", "fc2_BdDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_110/bdt_output.txt")
    fc2_BsDDKp_mc_2016 = ROOT.TFileCollection("fc2_BsDDKp_mc_2016", "fc2_BsDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_120/bdt_output.txt")
    fc2_BuDDK0_mc_2016 = ROOT.TFileCollection("fc2_BuDDK0_mc_2016", "fc2_BuDDK0_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_130/bdt_output.txt")
    fc2_BuDD_mc_2016 = ROOT.TFileCollection("fc2_BuDD_mc_2016", "fc2_BuDD_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_150/bdt_output.txt")

    fc2_BuDDKp_mc_2017 = ROOT.TFileCollection("fc2_BuDDKp_mc_2017", "fc2_BuDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_100/bdt_output.txt")
    fc2_BdDDKp_mc_2017 = ROOT.TFileCollection("fc2_BdDDKp_mc_2017", "fc2_BdDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_110/bdt_output.txt")
    fc2_BsDDKp_mc_2017 = ROOT.TFileCollection("fc2_BsDDKp_mc_2017", "fc2_BsDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_120/bdt_output.txt")
    fc2_BuDDK0_mc_2017 = ROOT.TFileCollection("fc2_BuDDK0_mc_2017", "fc2_BuDDK0_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_130/bdt_output.txt")
    fc2_BuDD_mc_2017 = ROOT.TFileCollection("fc2_BuDD_mc_2017", "fc2_BuDD_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_150/bdt_output.txt")

    fc2_BuDDKp_mc_2018 = ROOT.TFileCollection("fc2_BuDDKp_mc_2018", "fc2_BuDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_100/bdt_output.txt")
    fc2_BdDDKp_mc_2018 = ROOT.TFileCollection("fc2_BdDDKp_mc_2018", "fc2_BdDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_110/bdt_output.txt")
    fc2_BsDDKp_mc_2018 = ROOT.TFileCollection("fc2_BsDDKp_mc_2018", "fc2_BsDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_120/bdt_output.txt")
    fc2_BuDDK0_mc_2018 = ROOT.TFileCollection("fc2_BuDDK0_mc_2018", "fc2_BuDDK0_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_130/bdt_output.txt")
    fc2_BuDD_mc_2018 = ROOT.TFileCollection("fc2_BuDD_mc_2018", "fc2_BuDD_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_150/bdt_output.txt")

    t2_BuDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t2_BdDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t2_BsDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t2_BuDDK0_mc_2016 = ROOT.TChain("DecayTree")
    t2_BuDD_mc_2016 = ROOT.TChain("DecayTree")

    t2_BuDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t2_BdDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t2_BsDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t2_BuDDK0_mc_2017 = ROOT.TChain("DecayTree")
    t2_BuDD_mc_2017 = ROOT.TChain("DecayTree")

    t2_BuDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t2_BdDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t2_BsDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t2_BuDDK0_mc_2018 = ROOT.TChain("DecayTree")
    t2_BuDD_mc_2018 = ROOT.TChain("DecayTree")

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

    # mass tree
    fc3_BuDDKp_mc_2016 = ROOT.TFileCollection("fc3_BuDDKp_mc_2016", "fc3_BuDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_100/invariant_mass_tree.txt")
    fc3_BdDDKp_mc_2016 = ROOT.TFileCollection("fc3_BdDDKp_mc_2016", "fc3_BdDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_110/invariant_mass_tree.txt")
    fc3_BsDDKp_mc_2016 = ROOT.TFileCollection("fc3_BsDDKp_mc_2016", "fc3_BsDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_120/invariant_mass_tree.txt")
    fc3_BuDDK0_mc_2016 = ROOT.TFileCollection("fc3_BuDDK0_mc_2016", "fc3_BuDDK0_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_130/invariant_mass_tree.txt")
    fc3_BuDD_mc_2016 = ROOT.TFileCollection("fc3_BuDD_mc_2016", "fc3_BuDD_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_150/invariant_mass_tree.txt")

    fc3_BuDDKp_mc_2017 = ROOT.TFileCollection("fc3_BuDDKp_mc_2017", "fc3_BuDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_100/invariant_mass_tree.txt")
    fc3_BdDDKp_mc_2017 = ROOT.TFileCollection("fc3_BdDDKp_mc_2017", "fc3_BdDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_110/invariant_mass_tree.txt")
    fc3_BsDDKp_mc_2017 = ROOT.TFileCollection("fc3_BsDDKp_mc_2017", "fc3_BsDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_120/invariant_mass_tree.txt")
    fc3_BuDDK0_mc_2017 = ROOT.TFileCollection("fc3_BuDDK0_mc_2017", "fc3_BuDDK0_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_130/invariant_mass_tree.txt")
    fc3_BuDD_mc_2017 = ROOT.TFileCollection("fc3_BuDD_mc_2017", "fc3_BuDD_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_150/invariant_mass_tree.txt")

    fc3_BuDDKp_mc_2018 = ROOT.TFileCollection("fc3_BuDDKp_mc_2018", "fc3_BuDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_100/invariant_mass_tree.txt")
    fc3_BdDDKp_mc_2018 = ROOT.TFileCollection("fc3_BdDDKp_mc_2018", "fc3_BdDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_110/invariant_mass_tree.txt")
    fc3_BsDDKp_mc_2018 = ROOT.TFileCollection("fc3_BsDDKp_mc_2018", "fc3_BsDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_120/invariant_mass_tree.txt")
    fc3_BuDDK0_mc_2018 = ROOT.TFileCollection("fc3_BuDDK0_mc_2018", "fc3_BuDDK0_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_130/invariant_mass_tree.txt")
    fc3_BuDD_mc_2018 = ROOT.TFileCollection("fc3_BuDD_mc_2018", "3_BuDD_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_150/invariant_mass_tree.txt")

    t3_BuDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t3_BdDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t3_BsDDKp_mc_2016 = ROOT.TChain("DecayTree")
    t3_BuDDK0_mc_2016 = ROOT.TChain("DecayTree")
    t3_BuDD_mc_2016 = ROOT.TChain("DecayTree")

    t3_BuDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t3_BdDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t3_BsDDKp_mc_2017 = ROOT.TChain("DecayTree")
    t3_BuDDK0_mc_2017 = ROOT.TChain("DecayTree")
    t3_BuDD_mc_2017 = ROOT.TChain("DecayTree")

    t3_BuDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t3_BdDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t3_BsDDKp_mc_2018 = ROOT.TChain("DecayTree")
    t3_BuDDK0_mc_2018 = ROOT.TChain("DecayTree")
    t3_BuDD_mc_2018 = ROOT.TChain("DecayTree")

    t3_BuDDKp_mc_2016.AddFileInfoList(fc3_BuDDKp_mc_2016.GetList())
    t3_BdDDKp_mc_2016.AddFileInfoList(fc3_BdDDKp_mc_2016.GetList())
    t3_BsDDKp_mc_2016.AddFileInfoList(fc3_BsDDKp_mc_2016.GetList())
    t3_BuDDK0_mc_2016.AddFileInfoList(fc3_BuDDK0_mc_2016.GetList())
    t3_BuDD_mc_2016.AddFileInfoList(fc3_BuDD_mc_2016.GetList())

    t3_BuDDKp_mc_2017.AddFileInfoList(fc3_BuDDKp_mc_2017.GetList())
    t3_BdDDKp_mc_2017.AddFileInfoList(fc3_BdDDKp_mc_2017.GetList())
    t3_BsDDKp_mc_2017.AddFileInfoList(fc3_BsDDKp_mc_2017.GetList())
    t3_BuDDK0_mc_2017.AddFileInfoList(fc3_BuDDK0_mc_2017.GetList())
    t3_BuDD_mc_2017.AddFileInfoList(fc3_BuDD_mc_2017.GetList())

    t3_BuDDKp_mc_2018.AddFileInfoList(fc3_BuDDKp_mc_2018.GetList())
    t3_BdDDKp_mc_2018.AddFileInfoList(fc3_BdDDKp_mc_2018.GetList())
    t3_BsDDKp_mc_2018.AddFileInfoList(fc3_BsDDKp_mc_2018.GetList())
    t3_BuDDK0_mc_2018.AddFileInfoList(fc3_BuDDK0_mc_2018.GetList())
    t3_BuDD_mc_2018.AddFileInfoList(fc3_BuDD_mc_2018.GetList())

    t3_BuDDKp_mc_2016.Add(t3_BuDDKp_mc_2017)
    t3_BuDDKp_mc_2016.Add(t3_BuDDKp_mc_2018)

    t3_BdDDKp_mc_2016.Add(t3_BdDDKp_mc_2017)
    t3_BdDDKp_mc_2016.Add(t3_BdDDKp_mc_2018)

    t3_BsDDKp_mc_2016.Add(t3_BsDDKp_mc_2017)
    t3_BsDDKp_mc_2016.Add(t3_BsDDKp_mc_2018)

    t3_BuDDK0_mc_2016.Add(t3_BuDDK0_mc_2017)
    t3_BuDDK0_mc_2016.Add(t3_BuDDK0_mc_2018)

    t3_BuDD_mc_2016.Add(t3_BuDD_mc_2017)
    t3_BuDD_mc_2016.Add(t3_BuDD_mc_2018)

    if(name == "response_merged"):
        fc4_BuDDKp_mc_2016 = ROOT.TFileCollection("fc4_BuDDKp_mc_2016", "fc4_BuDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_100/bdt_output.txt")
        fc4_BdDDKp_mc_2016 = ROOT.TFileCollection("fc4_BdDDKp_mc_2016", "fc4_BdDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_110/bdt_output.txt")
        fc4_BsDDKp_mc_2016 = ROOT.TFileCollection("fc4_BsDDKp_mc_2016", "fc4_BsDDKp_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_120/bdt_output.txt")
        fc4_BuDDK0_mc_2016 = ROOT.TFileCollection("fc4_BuDDK0_mc_2016", "fc4_BuDDK0_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_130/bdt_output.txt")
        fc4_BuDD_mc_2016 = ROOT.TFileCollection("fc4_BuDD_mc_2016", "fc4_BuDD_mc_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_150/bdt_output.txt")

        fc4_BuDDKp_mc_2017 = ROOT.TFileCollection("fc4_BuDDKp_mc_2017", "fc4_BuDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_100/bdt_output.txt")
        fc4_BdDDKp_mc_2017 = ROOT.TFileCollection("fc4_BdDDKp_mc_2017", "fc4_BdDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_110/bdt_output.txt")
        fc4_BsDDKp_mc_2017 = ROOT.TFileCollection("fc4_BsDDKp_mc_2017", "fc4_BsDDKp_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_120/bdt_output.txt")
        fc4_BuDDK0_mc_2017 = ROOT.TFileCollection("fc4_BuDDK0_mc_2017", "fc4_BuDDK0_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_130/bdt_output.txt")
        fc4_BuDD_mc_2017 = ROOT.TFileCollection("fc4_BuDD_mc_2017", "fc4_BuDD_mc_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_150/bdt_output.txt")

        fc4_BuDDKp_mc_2018 = ROOT.TFileCollection("fc4_BuDDKp_mc_2018", "fc4_BuDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_100/bdt_output.txt")
        fc4_BdDDKp_mc_2018 = ROOT.TFileCollection("fc4_BdDDKp_mc_2018", "fc4_BdDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_110/bdt_output.txt")
        fc4_BsDDKp_mc_2018 = ROOT.TFileCollection("fc4_BsDDKp_mc_2018", "fc4_BsDDKp_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_120/bdt_output.txt")
        fc4_BuDDK0_mc_2018 = ROOT.TFileCollection("fc4_BuDDK0_mc_2018", "fc4_BuDDK0_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_130/bdt_output.txt")
        fc4_BuDD_mc_2018 = ROOT.TFileCollection("fc4_BuDD_mc_2018", "fc4_BuDD_mc_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_150/bdt_output.txt")

        t4_BuDDKp_mc_2016 = ROOT.TChain("DecayTree")
        t4_BdDDKp_mc_2016 = ROOT.TChain("DecayTree")
        t4_BsDDKp_mc_2016 = ROOT.TChain("DecayTree")
        t4_BuDDK0_mc_2016 = ROOT.TChain("DecayTree")
        t4_BuDD_mc_2016 = ROOT.TChain("DecayTree")

        t4_BuDDKp_mc_2017 = ROOT.TChain("DecayTree")
        t4_BdDDKp_mc_2017 = ROOT.TChain("DecayTree")
        t4_BsDDKp_mc_2017 = ROOT.TChain("DecayTree")
        t4_BuDDK0_mc_2017 = ROOT.TChain("DecayTree")
        t4_BuDD_mc_2017 = ROOT.TChain("DecayTree")

        t4_BuDDKp_mc_2018 = ROOT.TChain("DecayTree")
        t4_BdDDKp_mc_2018 = ROOT.TChain("DecayTree")
        t4_BsDDKp_mc_2018 = ROOT.TChain("DecayTree")
        t4_BuDDK0_mc_2018 = ROOT.TChain("DecayTree")
        t4_BuDD_mc_2018 = ROOT.TChain("DecayTree")

        t4_BuDDKp_mc_2016.AddFileInfoList(fc4_BuDDKp_mc_2016.GetList())
        t4_BdDDKp_mc_2016.AddFileInfoList(fc4_BdDDKp_mc_2016.GetList())
        t4_BsDDKp_mc_2016.AddFileInfoList(fc4_BsDDKp_mc_2016.GetList())
        t4_BuDDK0_mc_2016.AddFileInfoList(fc4_BuDDK0_mc_2016.GetList())
        t4_BuDD_mc_2016.AddFileInfoList(fc4_BuDD_mc_2016.GetList())

        t4_BuDDKp_mc_2017.AddFileInfoList(fc4_BuDDKp_mc_2017.GetList())
        t4_BdDDKp_mc_2017.AddFileInfoList(fc4_BdDDKp_mc_2017.GetList())
        t4_BsDDKp_mc_2017.AddFileInfoList(fc4_BsDDKp_mc_2017.GetList())
        t4_BuDDK0_mc_2017.AddFileInfoList(fc4_BuDDK0_mc_2017.GetList())
        t4_BuDD_mc_2017.AddFileInfoList(fc4_BuDD_mc_2017.GetList())

        t4_BuDDKp_mc_2018.AddFileInfoList(fc4_BuDDKp_mc_2018.GetList())
        t4_BdDDKp_mc_2018.AddFileInfoList(fc4_BdDDKp_mc_2018.GetList())
        t4_BsDDKp_mc_2018.AddFileInfoList(fc4_BsDDKp_mc_2018.GetList())
        t4_BuDDK0_mc_2018.AddFileInfoList(fc4_BuDDK0_mc_2018.GetList())
        t4_BuDD_mc_2018.AddFileInfoList(fc4_BuDD_mc_2018.GetList())

        t4_BuDDKp_mc_2016.Add(t4_BuDDKp_mc_2017)
        t4_BuDDKp_mc_2016.Add(t4_BuDDKp_mc_2018)

        t4_BdDDKp_mc_2016.Add(t4_BdDDKp_mc_2017)
        t4_BdDDKp_mc_2016.Add(t4_BdDDKp_mc_2018)

        t4_BsDDKp_mc_2016.Add(t4_BsDDKp_mc_2017)
        t4_BsDDKp_mc_2016.Add(t4_BsDDKp_mc_2018)

        t4_BuDDK0_mc_2016.Add(t4_BuDDK0_mc_2017)
        t4_BuDDK0_mc_2016.Add(t4_BuDDK0_mc_2018)

        t4_BuDD_mc_2016.Add(t4_BuDD_mc_2017)
        t4_BuDD_mc_2016.Add(t4_BuDD_mc_2018)

    t_cocktail_mc = ROOT.TChain("DecayTree")
    t1_cocktail_mc = ROOT.TChain("DecayTree")
    t2_cocktail_mc = ROOT.TChain("DecayTree")
    t3_cocktail_mc = ROOT.TChain("DecayTree")
    t4_cocktail_mc = ROOT.TChain("DecayTree")

    t_cocktail_mc.Add(t_BuDDKp_mc_2016)
    t_cocktail_mc.Add(t_BdDDKp_mc_2016)
    t_cocktail_mc.Add(t_BsDDKp_mc_2016)
    t_cocktail_mc.Add(t_BuDDK0_mc_2016)
    t_cocktail_mc.Add(t_BuDD_mc_2016)

    t1_cocktail_mc.Add(t1_BuDDKp_mc_2016)
    t1_cocktail_mc.Add(t1_BdDDKp_mc_2016)
    t1_cocktail_mc.Add(t1_BsDDKp_mc_2016)
    t1_cocktail_mc.Add(t1_BuDDK0_mc_2016)
    t1_cocktail_mc.Add(t1_BuDD_mc_2016)

    t2_cocktail_mc.Add(t2_BuDDKp_mc_2016)
    t2_cocktail_mc.Add(t2_BdDDKp_mc_2016)
    t2_cocktail_mc.Add(t2_BsDDKp_mc_2016)
    t2_cocktail_mc.Add(t2_BuDDK0_mc_2016)
    t2_cocktail_mc.Add(t2_BuDD_mc_2016)

    t3_cocktail_mc.Add(t3_BuDDKp_mc_2016)
    t3_cocktail_mc.Add(t3_BdDDKp_mc_2016)
    t3_cocktail_mc.Add(t3_BsDDKp_mc_2016)
    t3_cocktail_mc.Add(t3_BuDDK0_mc_2016)
    t3_cocktail_mc.Add(t3_BuDD_mc_2016)

    t_cocktail_mc.AddFriend(t1_cocktail_mc)
    t_cocktail_mc.AddFriend(t2_cocktail_mc)
    t_cocktail_mc.AddFriend(t3_cocktail_mc)
    
    if(name == "response_merged"):
        t4_cocktail_mc.Add(t4_BuDDKp_mc_2016)
        t4_cocktail_mc.Add(t4_BdDDKp_mc_2016)
        t4_cocktail_mc.Add(t4_BsDDKp_mc_2016)
        t4_cocktail_mc.Add(t4_BuDDK0_mc_2016)
        t4_cocktail_mc.Add(t4_BuDD_mc_2016)
    
        t_cocktail_mc.AddFriend(t4_cocktail_mc)

    t_BuDDKp_mc_2016.AddFriend(t1_BuDDKp_mc_2016)
    t_BuDDKp_mc_2016.AddFriend(t2_BuDDKp_mc_2016)
    t_BuDDKp_mc_2016.AddFriend(t3_BuDDKp_mc_2016)

    t_BdDDKp_mc_2016.AddFriend(t1_BdDDKp_mc_2016)
    t_BdDDKp_mc_2016.AddFriend(t2_BdDDKp_mc_2016)
    t_BdDDKp_mc_2016.AddFriend(t3_BdDDKp_mc_2016)

    t_BsDDKp_mc_2016.AddFriend(t1_BsDDKp_mc_2016)
    t_BsDDKp_mc_2016.AddFriend(t2_BsDDKp_mc_2016)
    t_BsDDKp_mc_2016.AddFriend(t3_BsDDKp_mc_2016)

    t_BuDDK0_mc_2016.AddFriend(t1_BuDDK0_mc_2016)
    t_BuDDK0_mc_2016.AddFriend(t2_BuDDK0_mc_2016)
    t_BuDDK0_mc_2016.AddFriend(t3_BuDDK0_mc_2016)

    t_BuDD_mc_2016.AddFriend(t1_BuDD_mc_2016)
    t_BuDD_mc_2016.AddFriend(t2_BuDD_mc_2016)
    t_BuDD_mc_2016.AddFriend(t3_BuDD_mc_2016)

    if(name == "response_merged"):
        t_BuDDKp_mc_2016.AddFriend(t4_BuDDKp_mc_2016)
        t_BdDDKp_mc_2016.AddFriend(t4_BdDDKp_mc_2016)
        t_BsDDKp_mc_2016.AddFriend(t4_BsDDKp_mc_2016)
        t_BuDDK0_mc_2016.AddFriend(t4_BuDDK0_mc_2016)
        t_BuDD_mc_2016.AddFriend(t4_BuDD_mc_2016)

    ######################################################################################################################################################

    ###################################################### RS and WS data #####################################################################    
    # pre-sel tree
    fc_rs_data_2016 = ROOT.TFileCollection("fc_rs_data_2016", "fc_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_2/pre_sel_tree.txt", 100)
    fc_rs_data_2017 = ROOT.TFileCollection("fc_rs_data_2017", "fc_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_2/pre_sel_tree.txt", 100)
    fc_rs_data_2018 = ROOT.TFileCollection("fc_rs_data_2018", "fc_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_2/pre_sel_tree.txt", 100)

    fc_ws_data_2016 = ROOT.TFileCollection("fc_ws_data_2016", "fc_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 100)
    fc_ws_data_2017 = ROOT.TFileCollection("fc_ws_data_2017", "fc_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 100)
    fc_ws_data_2018 = ROOT.TFileCollection("fc_ws_data_2018", "fc_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 100)

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
    fc1_rs_data_2016 = ROOT.TFileCollection("fc1_rs_data_2016", "fc1_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt", 100)
    fc1_rs_data_2017 = ROOT.TFileCollection("fc1_rs_data_2017", "fc1_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt", 100)
    fc1_rs_data_2018 = ROOT.TFileCollection("fc1_rs_data_2018", "fc1_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt", 100)

    fc1_ws_data_2016 = ROOT.TFileCollection("fc1_ws_data_2016", "fc1_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 100)
    fc1_ws_data_2017 = ROOT.TFileCollection("fc1_ws_data_2017", "fc1_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 100)
    fc1_ws_data_2018 = ROOT.TFileCollection("fc1_ws_data_2018", "fc1_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 100)

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
    fc2_rs_data_2016 = ROOT.TFileCollection("fc2_rs_data_2016", "fc2_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_2/bdt_output.txt", 100)
    fc2_rs_data_2017 = ROOT.TFileCollection("fc2_rs_data_2017", "fc2_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_2/bdt_output.txt", 100)
    fc2_rs_data_2018 = ROOT.TFileCollection("fc2_rs_data_2018", "fc2_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_2/bdt_output.txt", 100)

    fc2_ws_data_2016 = ROOT.TFileCollection("fc2_ws_data_2016", "fc2_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_3/bdt_output.txt", 100)
    fc2_ws_data_2017 = ROOT.TFileCollection("fc2_ws_data_2017", "fc2_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_3/bdt_output.txt", 100)
    fc2_ws_data_2018 = ROOT.TFileCollection("fc2_ws_data_2018", "fc2_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_3/bdt_output.txt", 100)

    t2_rs_data_2016 = ROOT.TChain("DecayTree")
    t2_rs_data_2017 = ROOT.TChain("DecayTree")
    t2_rs_data_2018 = ROOT.TChain("DecayTree")

    t2_ws_data_2016 = ROOT.TChain("DecayTree")
    t2_ws_data_2017 = ROOT.TChain("DecayTree")
    t2_ws_data_2018 = ROOT.TChain("DecayTree")

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

    # mass
    fc3_rs_data_2016 = ROOT.TFileCollection("fc3_rs_data_2016", "fc3_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_2/invariant_mass_tree.txt", 100)
    fc3_rs_data_2017 = ROOT.TFileCollection("fc3_rs_data_2017", "fc3_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_2/invariant_mass_tree.txt", 100)
    fc3_rs_data_2018 = ROOT.TFileCollection("fc3_rs_data_2018", "fc3_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_2/invariant_mass_tree.txt", 100)

    fc3_ws_data_2016 = ROOT.TFileCollection("fc3_ws_data_2016", "fc3_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_3/invariant_mass_tree.txt", 100)
    fc3_ws_data_2017 = ROOT.TFileCollection("fc3_ws_data_2017", "fc3_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_3/invariant_mass_tree.txt", 100)
    fc3_ws_data_2018 = ROOT.TFileCollection("fc3_ws_data_2018", "fc3_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_3/invariant_mass_tree.txt", 100)

    t3_rs_data_2016 = ROOT.TChain("DecayTree")
    t3_rs_data_2017 = ROOT.TChain("DecayTree")
    t3_rs_data_2018 = ROOT.TChain("DecayTree")

    t3_ws_data_2016 = ROOT.TChain("DecayTree")
    t3_ws_data_2017 = ROOT.TChain("DecayTree")
    t3_ws_data_2018 = ROOT.TChain("DecayTree")

    t3_rs_data_2016.AddFileInfoList(fc3_rs_data_2016.GetList())
    t3_rs_data_2017.AddFileInfoList(fc3_rs_data_2017.GetList())
    t3_rs_data_2018.AddFileInfoList(fc3_rs_data_2018.GetList())

    t3_ws_data_2016.AddFileInfoList(fc3_ws_data_2016.GetList())
    t3_ws_data_2017.AddFileInfoList(fc3_ws_data_2017.GetList())
    t3_ws_data_2018.AddFileInfoList(fc3_ws_data_2018.GetList())

    t3_rs_data_2016.Add(t3_rs_data_2017)
    t3_rs_data_2016.Add(t3_rs_data_2018)

    t3_ws_data_2016.Add(t3_ws_data_2017)
    t3_ws_data_2016.Add(t3_ws_data_2018)

    if(name == "response_merged"):
        # sklearn merged 
        fc4_rs_data_2016 = ROOT.TFileCollection("fc4_rs_data_2016", "fc4_rs_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_2/bdt_output.txt", 100)
        fc4_rs_data_2017 = ROOT.TFileCollection("fc4_rs_data_2017", "fc4_rs_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_2/bdt_output.txt", 100)
        fc4_rs_data_2018 = ROOT.TFileCollection("fc4_rs_data_2018", "fc4_rs_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_2/bdt_output.txt", 100)

        fc4_ws_data_2016 = ROOT.TFileCollection("fc4_ws_data_2016", "fc4_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_3/bdt_output.txt", 100)
        fc4_ws_data_2017 = ROOT.TFileCollection("fc4_ws_data_2017", "fc4_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_3/bdt_output.txt", 100)
        fc4_ws_data_2018 = ROOT.TFileCollection("fc4_ws_data_2018", "fc4_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_3/bdt_output.txt", 100)

        t4_rs_data_2016 = ROOT.TChain("DecayTree")
        t4_rs_data_2017 = ROOT.TChain("DecayTree")
        t4_rs_data_2018 = ROOT.TChain("DecayTree")

        t4_ws_data_2016 = ROOT.TChain("DecayTree")
        t4_ws_data_2017 = ROOT.TChain("DecayTree")
        t4_ws_data_2018 = ROOT.TChain("DecayTree")

        t4_rs_data_2016.AddFileInfoList(fc4_rs_data_2016.GetList())
        t4_rs_data_2017.AddFileInfoList(fc4_rs_data_2017.GetList())
        t4_rs_data_2018.AddFileInfoList(fc4_rs_data_2018.GetList())

        t4_ws_data_2016.AddFileInfoList(fc4_ws_data_2016.GetList())
        t4_ws_data_2017.AddFileInfoList(fc4_ws_data_2017.GetList())
        t4_ws_data_2018.AddFileInfoList(fc4_ws_data_2018.GetList())

        t4_rs_data_2016.Add(t4_rs_data_2017)
        t4_rs_data_2016.Add(t4_rs_data_2018)

        t4_ws_data_2016.Add(t4_ws_data_2017)
        t4_ws_data_2016.Add(t4_ws_data_2018)


    t_rs_data_2016.AddFriend(t1_rs_data_2016)
    t_rs_data_2016.AddFriend(t2_rs_data_2016)
    t_rs_data_2016.AddFriend(t3_rs_data_2016)

    t_ws_data_2016.AddFriend(t1_ws_data_2016)
    t_ws_data_2016.AddFriend(t2_ws_data_2016)
    t_ws_data_2016.AddFriend(t3_ws_data_2016)

    if(name == "response_merged"):
        t_rs_data_2016.AddFriend(t4_rs_data_2016)
        t_ws_data_2016.AddFriend(t4_ws_data_2016)


    ######################################################################################################################################################

    ###################################################### Normalisation mode MC (w/o rectangular cuts) #####################################################################    
    # pre-sel tree
    fc_norm_mc_no_rect_cuts_2016 = ROOT.TFileCollection("fc_norm_mc_no_rect_cuts_2016", "fc_norm_mc_no_rect_cuts_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_71/pre_sel_tree.txt")
    fc_norm_mc_no_rect_cuts_2017 = ROOT.TFileCollection("fc_norm_mc_no_rect_cuts_2017", "fc_norm_mc_no_rect_cuts_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_71/pre_sel_tree.txt")
    fc_norm_mc_no_rect_cuts_2018 = ROOT.TFileCollection("fc_norm_mc_no_rect_cuts_2018", "fc_norm_mc_no_rect_cuts_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_71/pre_sel_tree.txt")

    t_norm_mc_no_rect_cuts_2016 = ROOT.TChain("DecayTree")
    t_norm_mc_no_rect_cuts_2017 = ROOT.TChain("DecayTree")
    t_norm_mc_no_rect_cuts_2018 = ROOT.TChain("DecayTree")

    t_norm_mc_no_rect_cuts_2016.AddFileInfoList(fc_norm_mc_no_rect_cuts_2016.GetList())
    t_norm_mc_no_rect_cuts_2017.AddFileInfoList(fc_norm_mc_no_rect_cuts_2017.GetList())
    t_norm_mc_no_rect_cuts_2018.AddFileInfoList(fc_norm_mc_no_rect_cuts_2018.GetList())

    t_norm_mc_no_rect_cuts_2016.Add(t_norm_mc_no_rect_cuts_2017)
    t_norm_mc_no_rect_cuts_2016.Add(t_norm_mc_no_rect_cuts_2018)

    # sklearn
    fc2_norm_mc_no_rect_cuts_2016 = ROOT.TFileCollection("fc2_norm_mc_no_rect_cuts_2016", "fc2_norm_mc_no_rect_cuts_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_71/bdt_output.txt")
    fc2_norm_mc_no_rect_cuts_2017 = ROOT.TFileCollection("fc2_norm_mc_no_rect_cuts_2017", "fc2_norm_mc_no_rect_cuts_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_71/bdt_output.txt")
    fc2_norm_mc_no_rect_cuts_2018 = ROOT.TFileCollection("fc2_norm_mc_no_rect_cuts_2018", "fc2_norm_mc_no_rect_cuts_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_71/bdt_output.txt")

    t2_norm_mc_no_rect_cuts_2016 = ROOT.TChain("DecayTree")
    t2_norm_mc_no_rect_cuts_2017 = ROOT.TChain("DecayTree")
    t2_norm_mc_no_rect_cuts_2018 = ROOT.TChain("DecayTree")

    t2_norm_mc_no_rect_cuts_2016.AddFileInfoList(fc2_norm_mc_no_rect_cuts_2016.GetList())
    t2_norm_mc_no_rect_cuts_2017.AddFileInfoList(fc2_norm_mc_no_rect_cuts_2017.GetList())
    t2_norm_mc_no_rect_cuts_2018.AddFileInfoList(fc2_norm_mc_no_rect_cuts_2018.GetList())

    t2_norm_mc_no_rect_cuts_2016.Add(t2_norm_mc_no_rect_cuts_2017)
    t2_norm_mc_no_rect_cuts_2016.Add(t2_norm_mc_no_rect_cuts_2018)

    # sklearn merged
    if(name == "response_merged"):
        fc3_norm_mc_no_rect_cuts_2016 = ROOT.TFileCollection("fc3_norm_mc_no_rect_cuts_2016", "fc3_norm_mc_no_rect_cuts_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_71/bdt_output.txt")
        fc3_norm_mc_no_rect_cuts_2017 = ROOT.TFileCollection("fc3_norm_mc_no_rect_cuts_2017", "fc3_norm_mc_no_rect_cuts_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_71/bdt_output.txt")
        fc3_norm_mc_no_rect_cuts_2018 = ROOT.TFileCollection("fc3_norm_mc_no_rect_cuts_2018", "fc3_norm_mc_no_rect_cuts_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_71/bdt_output.txt")

        t3_norm_mc_no_rect_cuts_2016 = ROOT.TChain("DecayTree")
        t3_norm_mc_no_rect_cuts_2017 = ROOT.TChain("DecayTree")
        t3_norm_mc_no_rect_cuts_2018 = ROOT.TChain("DecayTree")

        t3_norm_mc_no_rect_cuts_2016.AddFileInfoList(fc3_norm_mc_no_rect_cuts_2016.GetList())
        t3_norm_mc_no_rect_cuts_2017.AddFileInfoList(fc3_norm_mc_no_rect_cuts_2017.GetList())
        t3_norm_mc_no_rect_cuts_2018.AddFileInfoList(fc3_norm_mc_no_rect_cuts_2018.GetList())

        t3_norm_mc_no_rect_cuts_2016.Add(t3_norm_mc_no_rect_cuts_2017)
        t3_norm_mc_no_rect_cuts_2016.Add(t3_norm_mc_no_rect_cuts_2018)

    t_norm_mc_no_rect_cuts_2016.AddFriend(t2_norm_mc_no_rect_cuts_2016)
    if(name == "response_merged"):
        t_norm_mc_no_rect_cuts_2016.AddFriend(t3_norm_mc_no_rect_cuts_2016)

    ######################################################################################################################################################
    
    ###################################################### Normalisation mode data (w/o rectangular cuts) #####################################################################    
    # pre-sel tree
    fc_norm_data_no_rect_cuts_2016 = ROOT.TFileCollection("fc_norm_data_no_rect_cuts_2016", "fc_norm_data_no_rect_cuts_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_81/pre_sel_tree.txt", 500)
    fc_norm_data_no_rect_cuts_2017 = ROOT.TFileCollection("fc_norm_data_no_rect_cuts_2017", "fc_norm_data_no_rect_cuts_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_81/pre_sel_tree.txt", 500)
    fc_norm_data_no_rect_cuts_2018 = ROOT.TFileCollection("fc_norm_data_no_rect_cuts_2018", "fc_norm_data_no_rect_cuts_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_81/pre_sel_tree.txt", 500)

    t_norm_data_no_rect_cuts_2016 = ROOT.TChain("DecayTree")
    t_norm_data_no_rect_cuts_2017 = ROOT.TChain("DecayTree")
    t_norm_data_no_rect_cuts_2018 = ROOT.TChain("DecayTree")

    t_norm_data_no_rect_cuts_2016.AddFileInfoList(fc_norm_data_no_rect_cuts_2016.GetList())
    t_norm_data_no_rect_cuts_2017.AddFileInfoList(fc_norm_data_no_rect_cuts_2017.GetList())
    t_norm_data_no_rect_cuts_2018.AddFileInfoList(fc_norm_data_no_rect_cuts_2018.GetList())

    t_norm_data_no_rect_cuts_2016.Add(t_norm_data_no_rect_cuts_2017)
    t_norm_data_no_rect_cuts_2016.Add(t_norm_data_no_rect_cuts_2018)

    # sklearn
    fc2_norm_data_no_rect_cuts_2016 = ROOT.TFileCollection("fc2_norm_data_no_rect_cuts_2016", "fc2_norm_data_no_rect_cuts_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2016/Species_81/bdt_output.txt", 500)
    fc2_norm_data_no_rect_cuts_2017 = ROOT.TFileCollection("fc2_norm_data_no_rect_cuts_2017", "fc2_norm_data_no_rect_cuts_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2017/Species_81/bdt_output.txt", 500)
    fc2_norm_data_no_rect_cuts_2018 = ROOT.TFileCollection("fc2_norm_data_no_rect_cuts_2018", "fc2_norm_data_no_rect_cuts_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response/2018/Species_81/bdt_output.txt", 500)

    t2_norm_data_no_rect_cuts_2016 = ROOT.TChain("DecayTree")
    t2_norm_data_no_rect_cuts_2017 = ROOT.TChain("DecayTree")
    t2_norm_data_no_rect_cuts_2018 = ROOT.TChain("DecayTree")

    t2_norm_data_no_rect_cuts_2016.AddFileInfoList(fc2_norm_data_no_rect_cuts_2016.GetList())
    t2_norm_data_no_rect_cuts_2017.AddFileInfoList(fc2_norm_data_no_rect_cuts_2017.GetList())
    t2_norm_data_no_rect_cuts_2018.AddFileInfoList(fc2_norm_data_no_rect_cuts_2018.GetList())

    t2_norm_data_no_rect_cuts_2016.Add(t2_norm_data_no_rect_cuts_2017)
    t2_norm_data_no_rect_cuts_2016.Add(t2_norm_data_no_rect_cuts_2018)

    # sklearn merged
    if(name == "response_merged"):
        fc3_norm_data_no_rect_cuts_2016 = ROOT.TFileCollection("fc3_norm_data_no_rect_cuts_2016", "fc3_norm_data_no_rect_cuts_2016", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2016/Species_81/bdt_output.txt", 500)
        fc3_norm_data_no_rect_cuts_2017 = ROOT.TFileCollection("fc3_norm_data_no_rect_cuts_2017", "fc3_norm_data_no_rect_cuts_2017", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2017/Species_81/bdt_output.txt", 500)
        fc3_norm_data_no_rect_cuts_2018 = ROOT.TFileCollection("fc3_norm_data_no_rect_cuts_2018", "fc3_norm_data_no_rect_cuts_2018", "/panfs/felician/B2Ktautau/workflow/sklearn_response_merged/2018/Species_81/bdt_output.txt", 500)

        t3_norm_data_no_rect_cuts_2016 = ROOT.TChain("DecayTree")
        t3_norm_data_no_rect_cuts_2017 = ROOT.TChain("DecayTree")
        t3_norm_data_no_rect_cuts_2018 = ROOT.TChain("DecayTree")

        t3_norm_data_no_rect_cuts_2016.AddFileInfoList(fc3_norm_data_no_rect_cuts_2016.GetList())
        t3_norm_data_no_rect_cuts_2017.AddFileInfoList(fc3_norm_data_no_rect_cuts_2017.GetList())
        t3_norm_data_no_rect_cuts_2018.AddFileInfoList(fc3_norm_data_no_rect_cuts_2018.GetList())

        t3_norm_data_no_rect_cuts_2016.Add(t3_norm_data_no_rect_cuts_2017)
        t3_norm_data_no_rect_cuts_2016.Add(t3_norm_data_no_rect_cuts_2018)

    t_norm_data_no_rect_cuts_2016.AddFriend(t2_norm_data_no_rect_cuts_2016)
    if(name == "response_merged"):
        t_norm_data_no_rect_cuts_2016.AddFriend(t3_norm_data_no_rect_cuts_2016)

    ######################################################################################################################################################


    ### PLOTS: 3 MC components
    # isolation BDT
    h_3pi3pi_mc_iso = ROOT.TH1D("h_3pi3pi_mc_iso", "h_3pi3pi_mc_iso", 30, 0, 1)
    h_3pi3pipi0_mc_iso = ROOT.TH1D("h_3pi3pipi0_mc_iso", "h_3pi3pipi0_mc_iso", 30, 0, 1)
    h_3pi3pi2pi0_mc_iso = ROOT.TH1D("h_3pi3pi2pi0_mc_iso", "h_3pi3pi2pi0_mc_iso", 30, 0, 1)

    t_BuKtautau_mc_2016.Draw("BDT1 >> h_3pi3pi_mc_iso", cut_Ktautau_mc+" && (component==0)")
    t_BuKtautau_mc_2016.Draw("BDT1 >> h_3pi3pipi0_mc_iso", cut_Ktautau_mc+" && (component==1)")
    t_BuKtautau_mc_2016.Draw("BDT1 >> h_3pi3pi2pi0_mc_iso", cut_Ktautau_mc+" && (component==2)")

    c = ROOT.TCanvas()
    c.cd()
    h_3pi3pi_mc_iso.GetXaxis().SetTitle("BDT1")
    h_3pi3pi_mc_iso.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_3pi3pi_mc_iso.SetTitle("Isolation BDT")
    h_3pi3pi_mc_iso.SetLineColor(4)
    h_3pi3pipi0_mc_iso.SetLineColor(8)
    h_3pi3pi2pi0_mc_iso.SetLineColor(2)
    h_3pi3pi_mc_iso.DrawNormalized()
    h_3pi3pipi0_mc_iso.DrawNormalized("same")
    h_3pi3pi2pi0_mc_iso.DrawNormalized("same")
    iso_leg = ROOT.TLegend(0.1, 0.7, 0.4, 0.89)
    iso_leg.AddEntry(h_3pi3pi_mc_iso, "3#pi3#pi MC", "lp")
    iso_leg.AddEntry(h_3pi3pipi0_mc_iso, "3#pi3#pi #pi^{0} MC", "lp")
    iso_leg.AddEntry(h_3pi3pi2pi0_mc_iso, "3#pi3#pi 2#pi^{0} MC", "lp")
    iso_leg.SetBorderSize(0)
    iso_leg.Draw("same")
    c.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/iso_ktautau_mc_components.pdf")

    # Topology BDT
    h_3pi3pi_mc_topo = ROOT.TH1D("h_3pi3pi_mc_topo", "h_3pi3pi_mc_topo", 30, 0, 1)
    h_3pi3pipi0_mc_topo = ROOT.TH1D("h_3pi3pipi0_mc_topo", "h_3pi3pipi0_mc_topo", 30, 0, 1)
    h_3pi3pi2pi0_mc_topo = ROOT.TH1D("h_3pi3pi2pi0_mc_topo", "h_3pi3pi2pi0_mc_topo", 30, 0, 1)

    t_BuKtautau_mc_2016.Draw("BDT2 >> h_3pi3pi_mc_topo", cut_Ktautau_mc+" && (component==0)")
    t_BuKtautau_mc_2016.Draw("BDT2 >> h_3pi3pipi0_mc_topo", cut_Ktautau_mc+" && (component==1)")
    t_BuKtautau_mc_2016.Draw("BDT2 >> h_3pi3pi2pi0_mc_topo", cut_Ktautau_mc+" && (component==2)")

    c1 = ROOT.TCanvas()
    c1.cd()
    h_3pi3pi_mc_topo.GetXaxis().SetTitle("BDT2")
    h_3pi3pi_mc_topo.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_3pi3pi_mc_topo.SetTitle("Topology BDT")
    h_3pi3pi_mc_topo.SetLineColor(4)
    h_3pi3pipi0_mc_topo.SetLineColor(8)
    h_3pi3pi2pi0_mc_topo.SetLineColor(2)
    h_3pi3pi_mc_topo.DrawNormalized()
    h_3pi3pipi0_mc_topo.DrawNormalized("same")
    h_3pi3pi2pi0_mc_topo.DrawNormalized("same")
    topo_leg = ROOT.TLegend(0.1, 0.7, 0.4, 0.89)
    topo_leg.AddEntry(h_3pi3pi_mc_topo, "3#pi3#pi MC", "lp")
    topo_leg.AddEntry(h_3pi3pipi0_mc_topo, "3#pi3#pi #pi^{0} MC", "lp")
    topo_leg.AddEntry(h_3pi3pi2pi0_mc_topo, "3#pi3#pi 2#pi^{0} MC", "lp")
    topo_leg.SetBorderSize(0)
    topo_leg.Draw("same")
    c1.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/topo_ktautau_mc_components.pdf")

    if(name == "response_merged"):
        # Merged BDT
        h_3pi3pi_mc_merged = ROOT.TH1D("h_3pi3pi_mc_merged", "h_3pi3pi_mc_merged", 30, 0, 1)
        h_3pi3pipi0_mc_merged = ROOT.TH1D("h_3pi3pipi0_mc_merged", "h_3pi3pipi0_mc_merged", 30, 0, 1)
        h_3pi3pi2pi0_mc_merged = ROOT.TH1D("h_3pi3pi2pi0_mc_merged", "h_3pi3pi2pi0_mc_merged", 30, 0, 1)

        t_BuKtautau_mc_2016.Draw("BDT >> h_3pi3pi_mc_merged", cut_Ktautau_mc+" && (component==0)")
        t_BuKtautau_mc_2016.Draw("BDT >> h_3pi3pipi0_mc_merged", cut_Ktautau_mc+" && (component==1)")
        t_BuKtautau_mc_2016.Draw("BDT >> h_3pi3pi2pi0_mc_merged", cut_Ktautau_mc+" && (component==2)")

        c11 = ROOT.TCanvas()
        c11.cd()
        h_3pi3pi_mc_merged.GetXaxis().SetTitle("BDT")
        h_3pi3pi_mc_merged.GetYaxis().SetTitle("Normalized entries / (0.03)")
        h_3pi3pi_mc_merged.SetTitle("Merged BDT")
        h_3pi3pi_mc_merged.SetLineColor(4)
        h_3pi3pipi0_mc_merged.SetLineColor(8)
        h_3pi3pi2pi0_mc_merged.SetLineColor(2)
        h_3pi3pi_mc_merged.DrawNormalized()
        h_3pi3pipi0_mc_merged.DrawNormalized("same")
        h_3pi3pi2pi0_mc_merged.DrawNormalized("same")
        merged_leg = ROOT.TLegend(0.1, 0.7, 0.4, 0.89)
        merged_leg.AddEntry(h_3pi3pi_mc_merged, "3#pi3#pi MC", "lp")
        merged_leg.AddEntry(h_3pi3pipi0_mc_merged, "3#pi3#pi #pi^{0} MC", "lp")
        merged_leg.AddEntry(h_3pi3pi2pi0_mc_merged, "3#pi3#pi 2#pi^{0} MC", "lp")
        merged_leg.SetBorderSize(0)
        merged_leg.Draw("same")
        c11.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/merged_ktautau_mc_components.pdf")

    ### PLOTS: cocktail MCs
    # Isolation BDT
    h_BuDDKp_mc_iso = ROOT.TH1D("h_BuDDKp_mc_iso", "h_BuDDKp_mc_iso", 30, 0, 1)
    h_BdDDKp_mc_iso = ROOT.TH1D("h_BdDDKp_mc_iso", "h_BdDDKp_mc_iso", 30, 0, 1)
    h_BsDDKp_mc_iso = ROOT.TH1D("h_BsDDKp_mc_iso", "h_BsDDKp_mc_iso", 30, 0, 1)
    h_BuDDK0_mc_iso = ROOT.TH1D("h_BuDDK0_mc_iso", "h_BuDDK0_mc_iso", 30, 0, 1)
    h_BuDD_mc_iso = ROOT.TH1D("h_BuDD_mc_iso", "h_BuDD_mc_iso", 30, 0, 1)

    t_BuDDKp_mc_2016.Draw("BDT1 >> h_BuDDKp_mc_iso", cut_cocktail_mc)
    t_BdDDKp_mc_2016.Draw("BDT1 >> h_BdDDKp_mc_iso", cut_cocktail_mc)
    t_BsDDKp_mc_2016.Draw("BDT1 >> h_BsDDKp_mc_iso", cut_cocktail_mc)
    t_BuDDK0_mc_2016.Draw("BDT1 >> h_BuDDK0_mc_iso", cut_cocktail_mc)
    t_BuDD_mc_2016.Draw("BDT1 >> h_BuDD_mc_iso", cut_cocktail_mc)

    c2 = ROOT.TCanvas()
    c2.cd()
    h_BsDDKp_mc_iso.GetXaxis().SetTitle("BDT1")
    h_BsDDKp_mc_iso.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_BsDDKp_mc_iso.SetTitle("Isolation BDT")
    h_BuDDKp_mc_iso.SetLineColor(4)
    h_BdDDKp_mc_iso.SetLineColor(8)
    h_BsDDKp_mc_iso.SetLineColor(2)
    h_BuDDK0_mc_iso.SetLineColor(51)
    h_BuDD_mc_iso.SetLineColor(94)
    h_BsDDKp_mc_iso.DrawNormalized()
    h_BdDDKp_mc_iso.DrawNormalized("same")
    h_BuDDKp_mc_iso.DrawNormalized("same")
    h_BuDDK0_mc_iso.DrawNormalized("same")
    h_BuDD_mc_iso.DrawNormalized("same")
    iso_leg_1 = ROOT.TLegend(0.4, 0.6, 0.7, 0.89)
    iso_leg_1.AddEntry(h_BuDDKp_mc_iso, "B^{+} #rightarrow D D K^{+}", "lp")
    iso_leg_1.AddEntry(h_BdDDKp_mc_iso, "B^{0} #rightarrow D D K^{+}", "lp")
    iso_leg_1.AddEntry(h_BsDDKp_mc_iso, "B^{0}_{s} #rightarrow D D K^{+}", "lp")
    iso_leg_1.AddEntry(h_BuDDK0_mc_iso, "B^{+} #rightarrow D D K^{0}", "lp")
    iso_leg_1.AddEntry(h_BuDD_mc_iso, "B^{+} #rightarrow D D ", "lp")
    iso_leg_1.SetBorderSize(0)
    iso_leg_1.Draw("same")
    c2.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/iso_cocktail_mcs.pdf")

    # Topology BDT
    h_BuDDKp_mc_topo = ROOT.TH1D("h_BuDDKp_mc_topo", "h_BuDDKp_mc_topo", 30, 0, 1)
    h_BdDDKp_mc_topo = ROOT.TH1D("h_BdDDKp_mc_topo", "h_BdDDKp_mc_topo", 30, 0, 1)
    h_BsDDKp_mc_topo = ROOT.TH1D("h_BsDDKp_mc_topo", "h_BsDDKp_mc_topo", 30, 0, 1)
    h_BuDDK0_mc_topo = ROOT.TH1D("h_BuDDK0_mc_topo", "h_BuDDK0_mc_topo", 30, 0, 1)
    h_BuDD_mc_topo = ROOT.TH1D("h_BuDD_mc_topo", "h_BuDD_mc_topo", 30, 0, 1)

    t_BuDDKp_mc_2016.Draw("BDT2 >> h_BuDDKp_mc_topo", cut_cocktail_mc)
    t_BdDDKp_mc_2016.Draw("BDT2 >> h_BdDDKp_mc_topo", cut_cocktail_mc)
    t_BsDDKp_mc_2016.Draw("BDT2 >> h_BsDDKp_mc_topo", cut_cocktail_mc)
    t_BuDDK0_mc_2016.Draw("BDT2 >> h_BuDDK0_mc_topo", cut_cocktail_mc)
    t_BuDD_mc_2016.Draw("BDT2 >> h_BuDD_mc_topo", cut_cocktail_mc)

    c3 = ROOT.TCanvas()
    c3.cd()
    h_BuDDK0_mc_topo.GetXaxis().SetTitle("BDT2")
    h_BuDDK0_mc_topo.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_BuDDK0_mc_topo.SetTitle("Topology BDT")
    h_BuDDKp_mc_topo.SetLineColor(4)
    h_BdDDKp_mc_topo.SetLineColor(8)
    h_BsDDKp_mc_topo.SetLineColor(2)
    h_BuDDK0_mc_topo.SetLineColor(51)
    h_BuDD_mc_topo.SetLineColor(94)
    h_BuDDK0_mc_topo.DrawNormalized()
    h_BdDDKp_mc_topo.DrawNormalized("same")
    h_BsDDKp_mc_topo.DrawNormalized("same")
    h_BuDD_mc_topo.DrawNormalized("same")
    h_BuDDKp_mc_topo.DrawNormalized("same")
    topo_leg_1 = ROOT.TLegend(0.4, 0.6, 0.7, 0.89)
    topo_leg_1.AddEntry(h_BuDDKp_mc_topo, "B^{+} #rightarrow D D K^{+}", "lp")
    topo_leg_1.AddEntry(h_BdDDKp_mc_topo, "B^{0} #rightarrow D D K^{+}", "lp")
    topo_leg_1.AddEntry(h_BsDDKp_mc_topo, "B^{0}_{s} #rightarrow D D K^{+}", "lp")
    topo_leg_1.AddEntry(h_BuDDK0_mc_topo, "B^{+} #rightarrow D D K^{0}", "lp")
    topo_leg_1.AddEntry(h_BuDD_mc_topo, "B^{+} #rightarrow D D ", "lp")
    topo_leg_1.SetBorderSize(0)
    topo_leg_1.Draw("same")
    c3.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/topo_cocktail_mcs.pdf")

    if(name == "response_merged"):
        # Merged BDT
        h_BuDDKp_mc_merged = ROOT.TH1D("h_BuDDKp_mc_merged", "h_BuDDKp_mc_merged", 30, 0, 1)
        h_BdDDKp_mc_merged = ROOT.TH1D("h_BdDDKp_mc_merged", "h_BdDDKp_mc_merged", 30, 0, 1)
        h_BsDDKp_mc_merged = ROOT.TH1D("h_BsDDKp_mc_merged", "h_BsDDKp_mc_merged", 30, 0, 1)
        h_BuDDK0_mc_merged = ROOT.TH1D("h_BuDDK0_mc_merged", "h_BuDDK0_mc_merged", 30, 0, 1)
        h_BuDD_mc_merged = ROOT.TH1D("h_BuDD_mc_merged", "h_BuDD_mc_merged", 30, 0, 1)

        t_BuDDKp_mc_2016.Draw("BDT >> h_BuDDKp_mc_merged", cut_cocktail_mc)
        t_BdDDKp_mc_2016.Draw("BDT >> h_BdDDKp_mc_merged", cut_cocktail_mc)
        t_BsDDKp_mc_2016.Draw("BDT >> h_BsDDKp_mc_merged", cut_cocktail_mc)
        t_BuDDK0_mc_2016.Draw("BDT >> h_BuDDK0_mc_merged", cut_cocktail_mc)
        t_BuDD_mc_2016.Draw("BDT >> h_BuDD_mc_merged", cut_cocktail_mc)

        c33 = ROOT.TCanvas()
        c33.cd()
        h_BuDD_mc_merged.GetXaxis().SetTitle("BDT")
        h_BuDD_mc_merged.GetYaxis().SetTitle("Normalized entries / (0.03)")
        h_BuDD_mc_merged.SetTitle("Merged BDT")
        h_BuDDKp_mc_merged.SetLineColor(4)
        h_BdDDKp_mc_merged.SetLineColor(8)
        h_BsDDKp_mc_merged.SetLineColor(2)
        h_BuDDK0_mc_merged.SetLineColor(51)
        h_BuDD_mc_merged.SetLineColor(94)
        h_BuDD_mc_merged.DrawNormalized()
        h_BdDDKp_mc_merged.DrawNormalized("same")
        h_BsDDKp_mc_merged.DrawNormalized("same")
        h_BuDDK0_mc_merged.DrawNormalized("same")
        h_BuDDKp_mc_merged.DrawNormalized("same")
        merged_leg_1 = ROOT.TLegend(0.4, 0.6, 0.7, 0.89)
        merged_leg_1.AddEntry(h_BuDDKp_mc_merged, "B^{+} #rightarrow D D K^{+}", "lp")
        merged_leg_1.AddEntry(h_BdDDKp_mc_merged, "B^{0} #rightarrow D D K^{+}", "lp")
        merged_leg_1.AddEntry(h_BsDDKp_mc_merged, "B^{0}_{s} #rightarrow D D K^{+}", "lp")
        merged_leg_1.AddEntry(h_BuDDK0_mc_merged, "B^{+} #rightarrow D D K^{0}", "lp")
        merged_leg_1.AddEntry(h_BuDD_mc_merged, "B^{+} #rightarrow D D ", "lp")
        merged_leg_1.SetBorderSize(0)
        merged_leg_1.Draw("same")
        c33.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/merged_cocktail_mcs.pdf")

    # PLOTS: RS and WS data
    # Isolation BDT
    h_rs_data_iso = ROOT.TH1D("h_rs_data_iso", "h_rs_data_iso", 30, 0, 1)
    h_ws_data_iso = ROOT.TH1D("h_ws_data_iso", "h_ws_data_iso", 30, 0, 1)

    t_rs_data_2016.Draw("BDT1 >> h_rs_data_iso", cut_Ktautau_data)
    t_ws_data_2016.Draw("BDT1 >> h_ws_data_iso", cut_Ktautau_data)

    c4 = ROOT.TCanvas()
    c4.cd()
    h_ws_data_iso.GetXaxis().SetTitle("BDT1")
    h_ws_data_iso.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_ws_data_iso.SetTitle("Isolation BDT")
    h_rs_data_iso.SetLineColor(4)
    h_ws_data_iso.SetLineColor(8)
    h_ws_data_iso.DrawNormalized()
    h_rs_data_iso.DrawNormalized("same")
    iso_leg_2 = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
    iso_leg_2.AddEntry(h_rs_data_iso, "RS data", "lp")
    iso_leg_2.AddEntry(h_ws_data_iso, "WS data", "lp")
    iso_leg_2.SetBorderSize(0)
    iso_leg_2.Draw("same")
    c4.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/iso_ktautau_data.pdf")

    # Topology BDT
    h_rs_data_topo = ROOT.TH1D("h_rs_data_topo", "h_rs_data_topo", 30, 0, 1)
    h_ws_data_topo = ROOT.TH1D("h_ws_data_topo", "h_ws_data_topo", 30, 0, 1)

    t_rs_data_2016.Draw("BDT2 >> h_rs_data_topo", cut_Ktautau_data)
    t_ws_data_2016.Draw("BDT2 >> h_ws_data_topo", cut_Ktautau_data)

    c5 = ROOT.TCanvas()
    c5.cd()
    h_rs_data_topo.GetXaxis().SetTitle("BDT2")
    h_rs_data_topo.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_rs_data_topo.SetTitle("Topology BDT")
    h_rs_data_topo.SetLineColor(4)
    h_ws_data_topo.SetLineColor(8)
    h_rs_data_topo.DrawNormalized()
    h_ws_data_topo.DrawNormalized("same")
    topo_leg_2 = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
    topo_leg_2.AddEntry(h_rs_data_topo, "RS data", "lp")
    topo_leg_2.AddEntry(h_ws_data_topo, "WS data", "lp")
    topo_leg_2.SetBorderSize(0)
    topo_leg_2.Draw("same")
    c5.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/topo_ktautau_data.pdf")

    if(name == "response_merged"):
        # Merged BDT
        h_rs_data_merged = ROOT.TH1D("h_rs_data_merged", "h_rs_data_merged", 30, 0, 1)
        h_ws_data_merged = ROOT.TH1D("h_ws_data_merged", "h_ws_data_merged", 30, 0, 1)

        t_rs_data_2016.Draw("BDT >> h_rs_data_merged", cut_Ktautau_data)
        t_ws_data_2016.Draw("BDT >> h_ws_data_merged", cut_Ktautau_data)

        c55 = ROOT.TCanvas()
        c55.cd()
        h_rs_data_merged.GetXaxis().SetTitle("BDT")
        h_rs_data_merged.GetYaxis().SetTitle("Normalized entries / (0.03)")
        h_rs_data_merged.SetTitle("Merged BDT")
        h_rs_data_merged.SetLineColor(4)
        h_ws_data_merged.SetLineColor(8)
        h_rs_data_merged.DrawNormalized()
        h_ws_data_merged.DrawNormalized("same")
        merged_leg_2 = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
        merged_leg_2.AddEntry(h_rs_data_merged, "RS data", "lp")
        merged_leg_2.AddEntry(h_ws_data_merged, "WS data", "lp")
        merged_leg_2.SetBorderSize(0)
        merged_leg_2.Draw("same")
        c55.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/merged_ktautau_data.pdf")

    # PLOTS: Signal vs background
    # Isolation BDT
    h_ktautau_mc_iso = ROOT.TH1D("h_ktautau_mc_iso", "h_ktautau_mc_iso", 30, 0, 1)
    t_BuKtautau_mc_2016.Draw("BDT1 >> h_ktautau_mc_iso", cut_Ktautau_mc)

    c5 = ROOT.TCanvas()
    c5.cd()
    h_ktautau_mc_iso.GetXaxis().SetTitle("BDT1")
    h_ktautau_mc_iso.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_ktautau_mc_iso.SetTitle("Isolation BDT")
    h_ktautau_mc_iso.SetLineColor(4)
    h_ws_data_iso.SetLineColor(2)
    h_ktautau_mc_iso.DrawNormalized("hist")
    h_ws_data_iso.DrawNormalized("same")
    iso_leg_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    iso_leg_3.AddEntry(h_ktautau_mc_iso, "Signal", "lp")
    iso_leg_3.AddEntry(h_ws_data_iso, "Background", "lp")
    iso_leg_3.SetBorderSize(0)
    iso_leg_3.SetBorderSize(0)
    iso_leg_3.Draw("same")
    c5.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/iso_sig_vs_bkg.pdf")

    # Topology BDT
    h_ktautau_mc_topo = ROOT.TH1D("h_ktautau_mc_topo", "h_ktautau_mc_topo", 30, 0, 1)
    t_BuKtautau_mc_2016.Draw("BDT2 >> h_ktautau_mc_topo", cut_Ktautau_mc)

    c6 = ROOT.TCanvas()
    c6.cd()
    h_ws_data_topo.GetXaxis().SetTitle("BDT2")
    h_ws_data_topo.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_ws_data_topo.SetTitle("Topology BDT")
    h_ktautau_mc_topo.SetLineColor(4)
    h_ws_data_topo.SetLineColor(2)
    h_ws_data_topo.DrawNormalized()
    h_ktautau_mc_topo.DrawNormalized("same")
    topo_leg_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_3.AddEntry(h_ktautau_mc_topo, "Signal", "lp")
    topo_leg_3.AddEntry(h_ws_data_topo, "Background", "lp")
    topo_leg_3.SetBorderSize(0)
    topo_leg_3.Draw("same")
    c6.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/topo_sig_vs_bkg.pdf")

    if(name == "response_merged"):
        # Merged BDT
        h_ktautau_mc_merged = ROOT.TH1D("h_ktautau_mc_merged", "h_ktautau_mc_merged", 30, 0, 1)
        t_BuKtautau_mc_2016.Draw("BDT >> h_ktautau_mc_merged", cut_Ktautau_mc)

        c66 = ROOT.TCanvas()
        c66.cd()
        h_ktautau_mc_merged.GetXaxis().SetTitle("BDT")
        h_ktautau_mc_merged.GetYaxis().SetTitle("Normalized entries / (0.03)")
        h_ktautau_mc_merged.SetTitle("Merged BDT")
        h_ktautau_mc_merged.SetLineColor(4)
        h_ws_data_merged.SetLineColor(2)
        h_ktautau_mc_merged.DrawNormalized()
        h_ws_data_merged.DrawNormalized("same")
        merged_leg_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_3.AddEntry(h_ktautau_mc_merged, "Signal", "lp")
        merged_leg_3.AddEntry(h_ws_data_merged, "Background", "lp")
        merged_leg_3.SetBorderSize(0)
        merged_leg_3.Draw("same")
        c66.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/merged_sig_vs_bkg.pdf")

    # PLOTS: BDT1 vs BDT2
    h_bdt1_bdt2_Ktautau_mc = ROOT.TH2D("h_bdt1_bdt2_Ktautau_mc", "h_bdt1_bdt2_Ktautau_mc", 30, 0, 1, 30, 0, 1)
    h_bdt1_bdt2_cocktail_mc = ROOT.TH2D("h_bdt1_bdt2_cocktail_mc", "h_bdt1_bdt2_cocktail_mc", 30, 0, 1, 30, 0, 1)
    h_bdt1_bdt2_ws_data = ROOT.TH2D("h_bdt1_bdt2_ws_data", "h_bdt1_bdt2_ws_data", 30, 0, 1, 30, 0, 1)
    h_bdt1_bdt2_rs_data = ROOT.TH2D("h_bdt1_bdt2_rs_data", "h_bdt1_bdt2_rs_data", 30, 0, 0.7, 30, 0, 0.7)

    t_BuKtautau_mc_2016.Draw("BDT1 : BDT2 >> h_bdt1_bdt2_Ktautau_mc", cut_Ktautau_mc)
    t_cocktail_mc.Draw("BDT1 : BDT2 >> h_bdt1_bdt2_cocktail_mc", cut_cocktail_mc)
    t_ws_data_2016.Draw("BDT1 : BDT2 >> h_bdt1_bdt2_ws_data", cut_Ktautau_data)
    t_rs_data_2016.Draw("BDT1 : BDT2 >> h_bdt1_bdt2_rs_data", cut_Ktautau_data)

    c7 = ROOT.TCanvas()
    c7.cd()
    c7.SetLogz()
    h_bdt1_bdt2_Ktautau_mc.GetXaxis().SetTitle("Topology (BDT2)")
    h_bdt1_bdt2_Ktautau_mc.GetYaxis().SetTitle("Isolation (BDT1)")
    h_bdt1_bdt2_Ktautau_mc.SetTitle(f"Ktautau MC: isolation vs topology BDT (#rho = {h_bdt1_bdt2_Ktautau_mc.GetCorrelationFactor():.2f})")
    h_bdt1_bdt2_Ktautau_mc.Draw("COLZ")
    c7.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bdt1_vs_bdt2_ktautau_mc.pdf")

    c8 = ROOT.TCanvas()
    c8.cd()
    c8.SetLogz()
    h_bdt1_bdt2_cocktail_mc.GetXaxis().SetTitle("Topology (BDT2)")
    h_bdt1_bdt2_cocktail_mc.GetYaxis().SetTitle("Isolation (BDT1)")
    h_bdt1_bdt2_cocktail_mc.SetTitle(f"Cocktail MCs: isolation vs topology BDT (#rho = {h_bdt1_bdt2_cocktail_mc.GetCorrelationFactor():.2f})")
    h_bdt1_bdt2_cocktail_mc.Draw("COLZ")
    c8.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bdt1_vs_bdt2_cocktail_mc.pdf")

    c9 = ROOT.TCanvas()
    c9.cd()
    c9.SetLogz()
    h_bdt1_bdt2_ws_data.GetXaxis().SetTitle("Topology (BDT2)")
    h_bdt1_bdt2_ws_data.GetYaxis().SetTitle("Isolation (BDT1)")
    h_bdt1_bdt2_ws_data.SetTitle(f"WS data: isolation vs topology BDT (#rho = {h_bdt1_bdt2_ws_data.GetCorrelationFactor():.2f})")
    h_bdt1_bdt2_ws_data.Draw("COLZ")
    c9.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bdt1_vs_bdt2_ws_data.pdf")

    c10 = ROOT.TCanvas()
    c10.cd()
    c10.SetLogz()
    h_bdt1_bdt2_rs_data.GetXaxis().SetTitle("Topology (BDT2)")
    h_bdt1_bdt2_rs_data.GetYaxis().SetTitle("Isolation (BDT1)")
    h_bdt1_bdt2_rs_data.SetTitle(f"RS data: isolation vs topology BDT (#rho = {h_bdt1_bdt2_rs_data.GetCorrelationFactor():.2f})")
    h_bdt1_bdt2_rs_data.Draw("COLZ")
    c10.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bdt1_vs_bdt2_rs_data.pdf")

    ### PLTOS: B+ mass evolution w/ BDT cut
    bdt_cuts = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9 ]
    colors = [1, 4, 6, 2, 8, 51]
    N = len(bdt_cuts)

    histos_ktautau_mc = []
    histos_cocktail_mc = []
    histos_rs_data = []
    histos_ws_data = []

    histos_bdt1_ktautau_mc = []
    histos_bdt1_cocktail_mc = []
    histos_bdt1_rs_data = []
    histos_bdt1_ws_data = []

    histos_bdt2_ktautau_mc = []
    histos_bdt2_cocktail_mc = []
    histos_bdt2_rs_data = []
    histos_bdt2_ws_data = []

    histos_bdt_ktautau_mc = []
    histos_bdt_cocktail_mc = []
    histos_bdt_rs_data = []
    histos_bdt_ws_data = []

    for i in range(N):
        h_mass_corr_ktautau_mc = ROOT.TH1D(f"h_mass_corr_ktautau_mc_{i}", f"h_mass_corr_ktautau_mc_{i}", 100, 4000, 9000)
        h_mass_corr_cocktail_mc = ROOT.TH1D(f"h_mass_corr_cocktail_mc_{i}", f"h_mass_corr_cocktail_mc_{i}", 100, 4000, 9000)
        h_mass_corr_rs_data = ROOT.TH1D(f"h_mass_corr_rs_data_{i}", f"h_mass_corr_rs_data_{i}", 100, 4000, 9000)
        h_mass_corr_ws_data = ROOT.TH1D(f"h_mass_corr_ws_data_{i}", f"h_mass_corr_ws_data_{i}", 100, 4000, 9000)

        h_mass_bdt1_corr_ktautau_mc = ROOT.TH1D(f"h_mass_bdt1_corr_ktautau_mc_{i}", f"h_mass_bdt1_corr_ktautau_mc_{i}", 100, 4000, 9000)
        h_mass_bdt1_corr_cocktail_mc = ROOT.TH1D(f"h_mass_bdt1_corr_cocktail_mc_{i}", f"h_mass_bdt1_corr_cocktail_mc_{i}", 100, 4000, 9000)
        h_mass_bdt1_corr_rs_data = ROOT.TH1D(f"h_mass_bdt1_corr_rs_data_{i}", f"h_mass_bdt1_corr_rs_data_{i}", 100, 4000, 9000)
        h_mass_bdt1_corr_ws_data = ROOT.TH1D(f"h_mass_bdt1_corr_ws_data_{i}", f"h_mass_bdt1_corr_ws_data_{i}", 100, 4000, 9000)

        h_mass_bdt2_corr_ktautau_mc = ROOT.TH1D(f"h_mass_bdt2_corr_ktautau_mc_{i}", f"h_mass_bdt2_corr_ktautau_mc_{i}", 100, 4000, 9000)
        h_mass_bdt2_corr_cocktail_mc = ROOT.TH1D(f"h_mass_bdt2_corr_cocktail_mc_{i}", f"h_mass_bdt2_corr_cocktail_mc_{i}", 100, 4000, 9000)
        h_mass_bdt2_corr_rs_data = ROOT.TH1D(f"h_mass_bdt2_corr_rs_data_{i}", f"h_mass_bdt2_corr_rs_data_{i}", 100, 4000, 9000)
        h_mass_bdt2_corr_ws_data = ROOT.TH1D(f"h_mass_bdt2_corr_ws_data_{i}", f"h_mass_bdt2_corr_ws_data_{i}", 100, 4000, 9000)

        if(name == "response_merged"):
            h_mass_bdt_corr_ktautau_mc = ROOT.TH1D(f"h_mass_bdt_corr_ktautau_mc_{i}", f"h_mass_bdt_corr_ktautau_mc_{i}", 100, 4000, 9000)
            h_mass_bdt_corr_cocktail_mc = ROOT.TH1D(f"h_mass_bdt_corr_cocktail_mc_{i}", f"h_mass_bdt_corr_cocktail_mc_{i}", 100, 4000, 9000)
            h_mass_bdt_corr_rs_data = ROOT.TH1D(f"h_mass_bdt_corr_rs_data_{i}", f"h_mass_bdt_corr_rs_data_{i}", 100, 4000, 9000)
            h_mass_bdt_corr_ws_data = ROOT.TH1D(f"h_mass_bdt_corr_ws_data_{i}", f"h_mass_bdt_corr_ws_data_{i}", 100, 4000, 9000)

        t_BuKtautau_mc_2016.Draw(f"df_Bp_M >> h_mass_corr_ktautau_mc_{i}", cut_Ktautau_mc+" && "+f"(BDT1 > {bdt_cuts[i]}) && (BDT2 > {bdt_cuts[i]})")
        t_cocktail_mc.Draw(f"df_Bp_M >> h_mass_corr_cocktail_mc_{i}", cut_cocktail_mc+" && "+f"(BDT1 > {bdt_cuts[i]}) && (BDT2 > {bdt_cuts[i]})")
        t_rs_data_2016.Draw(f"df_Bp_M >> h_mass_corr_rs_data_{i}", cut_Ktautau_data+" && "+f"(BDT1 > {bdt_cuts[i]}) && (BDT2 > {bdt_cuts[i]}) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_ws_data_2016.Draw(f"df_Bp_M >> h_mass_corr_ws_data_{i}", cut_Ktautau_data+" && "+f"(BDT1 > {bdt_cuts[i]}) && (BDT2 > {bdt_cuts[i]})")

        t_BuKtautau_mc_2016.Draw(f"df_Bp_M >> h_mass_bdt1_corr_ktautau_mc_{i}", cut_Ktautau_mc+" && "+f"(BDT1 > {bdt_cuts[i]})")
        t_cocktail_mc.Draw(f"df_Bp_M >> h_mass_bdt1_corr_cocktail_mc_{i}", cut_cocktail_mc+" && "+f"(BDT1 > {bdt_cuts[i]})")
        t_rs_data_2016.Draw(f"df_Bp_M >> h_mass_bdt1_corr_rs_data_{i}", cut_Ktautau_data+" && "+f"(BDT1 > {bdt_cuts[i]}) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_ws_data_2016.Draw(f"df_Bp_M >> h_mass_bdt1_corr_ws_data_{i}", cut_Ktautau_data+" && "+f"(BDT1 > {bdt_cuts[i]})")

        t_BuKtautau_mc_2016.Draw(f"df_Bp_M >> h_mass_bdt2_corr_ktautau_mc_{i}", cut_Ktautau_mc+" && "+f"(BDT2 > {bdt_cuts[i]})")
        t_cocktail_mc.Draw(f"df_Bp_M >> h_mass_bdt2_corr_cocktail_mc_{i}", cut_cocktail_mc+" && "+f"(BDT2 > {bdt_cuts[i]})")
        t_rs_data_2016.Draw(f"df_Bp_M >> h_mass_bdt2_corr_rs_data_{i}", cut_Ktautau_data+" && "+f"(BDT2 > {bdt_cuts[i]}) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_ws_data_2016.Draw(f"df_Bp_M >> h_mass_bdt2_corr_ws_data_{i}", cut_Ktautau_data+" && "+f"(BDT2 > {bdt_cuts[i]})")

        if(name == "response_merged"):
            t_BuKtautau_mc_2016.Draw(f"df_Bp_M >> h_mass_bdt_corr_ktautau_mc_{i}", cut_Ktautau_mc+" && "+f"(BDT > {bdt_cuts[i]})")
            t_cocktail_mc.Draw(f"df_Bp_M >> h_mass_bdt_corr_cocktail_mc_{i}", cut_cocktail_mc+" && "+f"(BDT > {bdt_cuts[i]})")
            t_rs_data_2016.Draw(f"df_Bp_M >> h_mass_bdt_corr_rs_data_{i}", cut_Ktautau_data+" && "+f"(BDT > {bdt_cuts[i]}) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
            t_ws_data_2016.Draw(f"df_Bp_M >> h_mass_bdt_corr_ws_data_{i}", cut_Ktautau_data+" && "+f"(BDT > {bdt_cuts[i]})")

        h_mass_corr_ktautau_mc.SetLineColor(colors[i])
        h_mass_corr_cocktail_mc.SetLineColor(colors[i])
        h_mass_corr_rs_data.SetLineColor(colors[i])
        h_mass_corr_ws_data.SetLineColor(colors[i])

        h_mass_bdt1_corr_ktautau_mc.SetLineColor(colors[i])
        h_mass_bdt1_corr_cocktail_mc.SetLineColor(colors[i])
        h_mass_bdt1_corr_rs_data.SetLineColor(colors[i])
        h_mass_bdt1_corr_ws_data.SetLineColor(colors[i])

        h_mass_bdt2_corr_ktautau_mc.SetLineColor(colors[i])
        h_mass_bdt2_corr_cocktail_mc.SetLineColor(colors[i])
        h_mass_bdt2_corr_rs_data.SetLineColor(colors[i])
        h_mass_bdt2_corr_ws_data.SetLineColor(colors[i])

        if(name == "response_merged"):
            h_mass_bdt_corr_ktautau_mc.SetLineColor(colors[i])
            h_mass_bdt_corr_cocktail_mc.SetLineColor(colors[i])
            h_mass_bdt_corr_rs_data.SetLineColor(colors[i])
            h_mass_bdt_corr_ws_data.SetLineColor(colors[i])

        histos_ktautau_mc.append(h_mass_corr_ktautau_mc)
        histos_cocktail_mc.append(h_mass_corr_cocktail_mc)
        histos_rs_data.append(h_mass_corr_rs_data)
        histos_ws_data.append(h_mass_corr_ws_data)

        histos_bdt1_ktautau_mc.append(h_mass_bdt1_corr_ktautau_mc)
        histos_bdt1_cocktail_mc.append(h_mass_bdt1_corr_cocktail_mc)
        histos_bdt1_rs_data.append(h_mass_bdt1_corr_rs_data)
        histos_bdt1_ws_data.append(h_mass_bdt1_corr_ws_data)

        histos_bdt2_ktautau_mc.append(h_mass_bdt2_corr_ktautau_mc)
        histos_bdt2_cocktail_mc.append(h_mass_bdt2_corr_cocktail_mc)
        histos_bdt2_rs_data.append(h_mass_bdt2_corr_rs_data)
        histos_bdt2_ws_data.append(h_mass_bdt2_corr_ws_data)

        if(name == "response_merged"):
            histos_bdt_ktautau_mc.append(h_mass_bdt_corr_ktautau_mc)
            histos_bdt_cocktail_mc.append(h_mass_bdt_corr_cocktail_mc)
            histos_bdt_rs_data.append(h_mass_bdt_corr_rs_data)
            histos_bdt_ws_data.append(h_mass_bdt_corr_ws_data)

    c11 = ROOT.TCanvas()
    c11.cd()
    leg11 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_ktautau_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_ktautau_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_ktautau_mc[j].SetTitle("Ktautau MC: B mass vs BDT cuts")
        leg11.AddEntry(histos_ktautau_mc[j], f"(BDT_{{iso}},BDT_{{topo}}) > ({bdt_cuts[j]},{bdt_cuts[j]})", "fp")
        if(i == 0):
            histos_ktautau_mc[j].DrawNormalized()
        else:
            histos_ktautau_mc[j].DrawNormalized("same")
    leg11.SetBorderSize(0)
    leg11.SetTextSize(0.03)
    leg11.Draw("same")
    c11.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_ktautau_mc.pdf")

    c11_bdt1 = ROOT.TCanvas()
    c11_bdt1.cd()
    leg11_bdt1 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt1_ktautau_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt1_ktautau_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt1_ktautau_mc[j].SetTitle("Ktautau MC: B mass vs isolation BDT")
        leg11_bdt1.AddEntry(histos_bdt1_ktautau_mc[j], f"BDT_{{iso}} > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt1_ktautau_mc[j].DrawNormalized()
        else:
            histos_bdt1_ktautau_mc[j].DrawNormalized("same")
    leg11_bdt1.SetBorderSize(0)
    leg11_bdt1.SetTextSize(0.03)
    leg11_bdt1.Draw("same")
    c11_bdt1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt1_corr_ktautau_mc.pdf")

    c11_bdt2 = ROOT.TCanvas()
    c11_bdt2.cd()
    leg11_bdt2 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt2_ktautau_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt2_ktautau_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt2_ktautau_mc[j].SetTitle("Ktautau MC: B mass vs topology BDT")
        leg11_bdt2.AddEntry(histos_bdt2_ktautau_mc[j], f"BDT_{{topo}} > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt2_ktautau_mc[j].DrawNormalized()
        else:
            histos_bdt2_ktautau_mc[j].DrawNormalized("same")
    leg11_bdt2.SetBorderSize(0)
    leg11_bdt2.SetTextSize(0.03)
    leg11_bdt2.Draw("same")
    c11_bdt2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt2_corr_ktautau_mc.pdf")

    if(name == "response_merged"):
        c11_bdt = ROOT.TCanvas()
        c11_bdt.cd()
        leg11_bdt = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
        for i in range(N):
            j = N-i-1
            histos_bdt_ktautau_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
            histos_bdt_ktautau_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
            histos_bdt_ktautau_mc[j].SetTitle("Ktautau MC: B mass vs merged BDT")
            leg11_bdt.AddEntry(histos_bdt_ktautau_mc[j], f"BDT > {bdt_cuts[j]}", "fp")
            if(i == 0):
                histos_bdt_ktautau_mc[j].DrawNormalized()
            else:
                histos_bdt_ktautau_mc[j].DrawNormalized("same")
        leg11_bdt.SetBorderSize(0)
        leg11_bdt.SetTextSize(0.03)
        leg11_bdt.Draw("same")
        c11_bdt.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_ktautau_mc.pdf")


    c12 = ROOT.TCanvas()
    c12.cd()
    leg12 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_cocktail_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_cocktail_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_cocktail_mc[j].SetTitle("Cocktail MCs: B mass vs BDT cuts")
        leg12.AddEntry(histos_cocktail_mc[j], f"(BDT_{{iso}},BDT_{{topo}}) > ({bdt_cuts[j]},{bdt_cuts[j]})", "fp")
        if(i == 0):
            histos_cocktail_mc[j].DrawNormalized()
        else:
            histos_cocktail_mc[j].DrawNormalized("same")
    leg12.SetBorderSize(0)
    leg12.SetTextSize(0.03)
    leg12.Draw("same")
    c12.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_cocktail_mc.pdf")

    c12_bdt1 = ROOT.TCanvas()
    c12_bdt1.cd()
    leg12_bdt1 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt1_cocktail_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt1_cocktail_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt1_cocktail_mc[j].SetTitle("Cocktail MCs: B mass vs isolation BDT")
        leg12_bdt1.AddEntry(histos_bdt1_cocktail_mc[j], f"BDT_{{iso}} > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt1_cocktail_mc[j].DrawNormalized()
        else:
            histos_bdt1_cocktail_mc[j].DrawNormalized("same")
    leg12_bdt1.SetBorderSize(0)
    leg12_bdt1.SetTextSize(0.03)
    leg12_bdt1.Draw("same")
    c12_bdt1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt1_corr_cocktail_mc.pdf")

    c12_bdt2 = ROOT.TCanvas()
    c12_bdt2.cd()
    leg12_bdt2 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt2_cocktail_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt2_cocktail_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt2_cocktail_mc[j].SetTitle("Cocktail MCs: B mass vs topology BDT")
        leg12_bdt2.AddEntry(histos_bdt2_cocktail_mc[j], f"BDT_{{topo}} > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt2_cocktail_mc[j].DrawNormalized()
        else:
            histos_bdt2_cocktail_mc[j].DrawNormalized("same")
    leg12_bdt2.SetBorderSize(0)
    leg12_bdt2.SetTextSize(0.03)
    leg12_bdt2.Draw("same")
    c12_bdt2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt2_corr_cocktail_mc.pdf")

    if(name == "response_merged"):
        c12_bdt = ROOT.TCanvas()
        c12_bdt.cd()
        leg12_bdt = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
        for i in range(N):
            j = N-i-1
            histos_bdt_cocktail_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
            histos_bdt_cocktail_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
            histos_bdt_cocktail_mc[j].SetTitle("Cocktail MCs: B mass vs merged BDT")
            leg12_bdt.AddEntry(histos_bdt_cocktail_mc[j], f"BDT > {bdt_cuts[j]}", "fp")
            if(i == 0):
                histos_bdt_cocktail_mc[j].DrawNormalized()
            else:
                histos_bdt_cocktail_mc[j].DrawNormalized("same")
        leg12_bdt.SetBorderSize(0)
        leg12_bdt.SetTextSize(0.03)
        leg12_bdt.Draw("same")
        c12_bdt.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_cocktail_mc.pdf")

    c13 = ROOT.TCanvas()
    c13.cd()
    leg13 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        histos_rs_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_rs_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_rs_data[i].SetTitle("RS data: B mass vs BDT cuts")
        leg13.AddEntry(histos_rs_data[i], f"(BDT_{{iso}},BDT_{{topo}}) > ({bdt_cuts[i]},{bdt_cuts[i]})", "fp")
        if(i == 0):
            histos_rs_data[i].DrawNormalized()
        else:
            histos_rs_data[i].DrawNormalized("same")
    leg13.SetBorderSize(0)
    leg13.SetTextSize(0.03)
    leg13.Draw("same")
    c13.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_rs_data.pdf")

    c13_bdt1 = ROOT.TCanvas()
    c13_bdt1.cd()
    leg13_bdt1 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        histos_bdt1_rs_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt1_rs_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt1_rs_data[i].SetTitle("RS data: B mass vs isolation BDT")
        leg13_bdt1.AddEntry(histos_bdt1_rs_data[i], f"BDT_{{iso}} > {bdt_cuts[i]}", "fp")
        if(i == 0):
            histos_bdt1_rs_data[i].DrawNormalized()
        else:
            histos_bdt1_rs_data[i].DrawNormalized("same")
    leg13_bdt1.SetBorderSize(0)
    leg13_bdt1.SetTextSize(0.03)
    leg13_bdt1.Draw("same")
    c13_bdt1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt1_corr_rs_data.pdf")

    c13_bdt2 = ROOT.TCanvas()
    c13_bdt2.cd()
    leg13_bdt2 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt2_rs_data[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt2_rs_data[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt2_rs_data[j].SetTitle("RS data: B mass vs topology BDT")
        leg13_bdt2.AddEntry(histos_bdt2_rs_data[i], f"BDT_{{topo}} > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt2_rs_data[j].DrawNormalized()
        else:
            histos_bdt2_rs_data[j].DrawNormalized("same")
    leg13_bdt2.SetBorderSize(0)
    leg13_bdt2.SetTextSize(0.03)
    leg13_bdt2.Draw("same")
    c13_bdt2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt2_corr_rs_data.pdf")

    if(name == "response_merged"):
        c13_bdt = ROOT.TCanvas()
        c13_bdt.cd()
        leg13_bdt = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
        for i in range(N):
            j = N-i-1
            histos_bdt_rs_data[j].GetXaxis().SetTitle("m_{B} (MeV)")
            histos_bdt_rs_data[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
            histos_bdt_rs_data[j].SetTitle("RS data: B mass vs merged BDT")
            leg13_bdt.AddEntry(histos_bdt_rs_data[i], f"BDT > {bdt_cuts[j]}", "fp")
            if(i == 0):
                histos_bdt_rs_data[j].DrawNormalized()
            else:
                histos_bdt_rs_data[j].DrawNormalized("same")
        leg13_bdt.SetBorderSize(0)
        leg13_bdt.SetTextSize(0.03)
        leg13_bdt.Draw("same")
        c13_bdt.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_rs_data.pdf")

    c14 = ROOT.TCanvas()
    c14.cd()
    leg14 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        histos_ws_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_ws_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_ws_data[i].SetTitle("WS data: B mass vs BDT cuts")
        leg14.AddEntry(histos_ws_data[i], f"(BDT_{{iso}},BDT_{{topo}}) > ({bdt_cuts[i]},{bdt_cuts[i]})", "fp")
        if(i == 0):
            histos_ws_data[i].DrawNormalized()
        else:
            histos_ws_data[i].DrawNormalized("same")
    leg14.SetBorderSize(0)
    leg14.SetTextSize(0.03)
    leg14.Draw("same")
    c14.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_ws_data.pdf")

    c14_bdt1 = ROOT.TCanvas()
    c14_bdt1.cd()
    leg14_bdt1 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        histos_bdt1_ws_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt1_ws_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt1_ws_data[i].SetTitle("WS data: B mass vs isolation BDT")
        leg14_bdt1.AddEntry(histos_bdt1_ws_data[i], f"BDT_{{iso}} > {bdt_cuts[i]}", "fp")
        if(i == 0):
            histos_bdt1_ws_data[i].DrawNormalized()
        else:
            histos_bdt1_ws_data[i].DrawNormalized("same")
    leg14_bdt1.SetBorderSize(0)
    leg14_bdt1.SetTextSize(0.03)
    leg14_bdt1.Draw("same")
    c14_bdt1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt1_corr_ws_data.pdf")

    c14_bdt2 = ROOT.TCanvas()
    c14_bdt2.cd()
    leg14_bdt2 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        histos_bdt2_ws_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt2_ws_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt2_ws_data[i].SetTitle("WS data: B mass vs topology BDT")
        leg14_bdt2.AddEntry(histos_bdt2_ws_data[i], f"BDT_{{topo}} > {bdt_cuts[i]}", "fp")
        if(i == 0):
            histos_bdt2_ws_data[i].DrawNormalized()
        else:
            histos_bdt2_ws_data[i].DrawNormalized("same")
    leg14_bdt2.SetBorderSize(0)
    leg14_bdt2.SetTextSize(0.03)
    leg14_bdt2.Draw("same")
    c14_bdt2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt2_corr_ws_data.pdf")

    if(name == "response_merged"):
        c14_bdt = ROOT.TCanvas()
        c14_bdt.cd()
        leg14_bdt = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
        for i in range(N):
            histos_bdt_ws_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
            histos_bdt_ws_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
            histos_bdt_ws_data[i].SetTitle("WS data: B mass vs merged BDT")
            leg14_bdt.AddEntry(histos_bdt_ws_data[i], f"BDT > {bdt_cuts[i]}", "fp")
            if(i == 0):
                histos_bdt_ws_data[i].DrawNormalized()
            else:
                histos_bdt_ws_data[i].DrawNormalized("same")
        leg14_bdt.SetBorderSize(0)
        leg14_bdt.SetTextSize(0.03)
        leg14_bdt.Draw("same")
        c14_bdt.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_ws_data.pdf")

    # RS vs WS mass at different BDT cuts (0, 0.5, 0.8, 0.9, 0.95)
    # Isolation BDT
    h_mass_iso_rs_data_1 = ROOT.TH1D("h_mass_iso_rs_data_1", "h_mass_iso_rs_data_1", 100, 4000, 9000)
    h_mass_iso_rs_data_2 = ROOT.TH1D("h_mass_iso_rs_data_2", "h_mass_rs_data_2", 100, 4000, 9000)
    h_mass_iso_rs_data_3 = ROOT.TH1D("h_mass_iso_rs_data_3", "h_mass_iso_rs_data_3", 100, 4000, 9000)
    h_mass_iso_rs_data_4 = ROOT.TH1D("h_mass_iso_rs_data_4", "h_mass_iso_rs_data_4", 100, 4000, 9000)
    h_mass_iso_rs_data_5 = ROOT.TH1D("h_mass_iso_rs_data_5", "h_mass_iso_rs_data_5", 100, 4000, 9000)

    h_mass_iso_ws_data_1 = ROOT.TH1D("h_mass_iso_ws_data_1", "h_mass_iso_ws_data_1", 100, 4000, 9000)
    h_mass_iso_ws_data_2 = ROOT.TH1D("h_mass_iso_ws_data_2", "h_mass_iso_ws_data_2", 100, 4000, 9000)
    h_mass_iso_ws_data_3 = ROOT.TH1D("h_mass_iso_ws_data_3", "h_mass_iso_ws_data_3", 100, 4000, 9000)
    h_mass_iso_ws_data_4 = ROOT.TH1D("h_mass_iso_ws_data_4", "h_mass_iso_ws_data_4", 100, 4000, 9000)
    h_mass_iso_ws_data_5 = ROOT.TH1D("h_mass_iso_ws_data_5", "h_mass_iso_ws_data_5", 100, 4000, 9000)

    t_rs_data_2016.Draw("df_Bp_M >> h_mass_iso_rs_data_1")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_iso_rs_data_2", "(BDT1 > 0.5)")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_iso_rs_data_3", "(BDT1 > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_iso_rs_data_4", "(BDT1 > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_iso_rs_data_5", "(BDT1 > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

    t_ws_data_2016.Draw("df_Bp_M >> h_mass_iso_ws_data_1")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_iso_ws_data_2", "(BDT1 > 0.5)")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_iso_ws_data_3", "(BDT1 > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_iso_ws_data_4", "(BDT1 > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_iso_ws_data_5", "(BDT1 > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

    c100_1 = ROOT.TCanvas()
    c100_1.cd()
    h_mass_iso_rs_data_1.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_iso_rs_data_1.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_iso_rs_data_1.SetTitle("Isolation: BDT1 > 0.0")
    h_mass_iso_rs_data_1.SetLineColor(1)
    h_mass_iso_ws_data_1.SetLineColor(2)
    h_mass_iso_rs_data_1.DrawNormalized()
    h_mass_iso_ws_data_1.DrawNormalized("same")
    topo_leg_100_1 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_100_1.AddEntry(h_mass_iso_rs_data_1, "RS data", "lp")
    topo_leg_100_1.AddEntry(h_mass_iso_ws_data_1, "WS data", "lp")
    topo_leg_100_1.SetBorderSize(0)
    topo_leg_100_1.Draw("same")
    c100_1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_isolation_BDT_0.0.pdf") 

    c100_2 = ROOT.TCanvas()
    c100_2.cd()
    h_mass_iso_ws_data_2.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_iso_ws_data_2.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_iso_ws_data_2.SetTitle("Isolation: BDT1 > 0.5")
    h_mass_iso_rs_data_2.SetLineColor(1)
    h_mass_iso_ws_data_2.SetLineColor(2)
    h_mass_iso_ws_data_2.DrawNormalized()
    h_mass_iso_rs_data_2.DrawNormalized("same")
    topo_leg_200_2 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_200_2.AddEntry(h_mass_iso_rs_data_2, "RS data", "lp")
    topo_leg_200_2.AddEntry(h_mass_iso_ws_data_2, "WS data", "lp")
    topo_leg_200_2.SetBorderSize(0)
    topo_leg_200_2.Draw("same")
    c100_2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_isolation_BDT_0.5.pdf") 

    c100_3 = ROOT.TCanvas()
    c100_3.cd()
    h_mass_iso_ws_data_3.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_iso_ws_data_3.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_iso_ws_data_3.SetTitle("Isolation: BDT1 > 0.8")
    h_mass_iso_rs_data_3.SetLineColor(1)
    h_mass_iso_ws_data_3.SetLineColor(2)
    h_mass_iso_ws_data_3.DrawNormalized()
    h_mass_iso_rs_data_3.DrawNormalized("same")
    topo_leg_300_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_300_3.AddEntry(h_mass_iso_rs_data_3, "RS data", "lp")
    topo_leg_300_3.AddEntry(h_mass_iso_ws_data_3, "WS data", "lp")
    topo_leg_300_3.SetBorderSize(0)
    topo_leg_300_3.Draw("same")
    c100_3.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_isolation_BDT_0.8.pdf") 

    c100_4 = ROOT.TCanvas()
    c100_4.cd()
    h_mass_iso_ws_data_4.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_iso_ws_data_4.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_iso_ws_data_4.SetTitle("Isolation: BDT1 > 0.9")
    h_mass_iso_rs_data_4.SetLineColor(1)
    h_mass_iso_ws_data_4.SetLineColor(2)
    h_mass_iso_ws_data_4.DrawNormalized()
    h_mass_iso_rs_data_4.DrawNormalized("same")
    topo_leg_400_4 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_400_4.AddEntry(h_mass_iso_rs_data_4, "RS data", "lp")
    topo_leg_400_4.AddEntry(h_mass_iso_ws_data_4, "WS data", "lp")
    topo_leg_400_4.SetBorderSize(0)
    topo_leg_400_4.Draw("same")
    c100_4.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_isolation_BDT_0.9.pdf") 

    c100_5 = ROOT.TCanvas()
    c100_5.cd()
    h_mass_iso_rs_data_5.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_iso_rs_data_5.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_iso_rs_data_5.SetTitle("Isolation: BDT1 > 0.95")
    h_mass_iso_rs_data_5.SetLineColor(1)
    h_mass_iso_ws_data_5.SetLineColor(2)
    h_mass_iso_rs_data_5.DrawNormalized()
    h_mass_iso_ws_data_5.DrawNormalized("same")
    topo_leg_500_5 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_500_5.AddEntry(h_mass_iso_rs_data_5, "RS data", "lp")
    topo_leg_500_5.AddEntry(h_mass_iso_ws_data_5, "WS data", "lp")
    topo_leg_500_5.SetBorderSize(0)
    topo_leg_500_5.Draw("same")
    c100_5.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_isolation_BDT_0.95.pdf") 

    # Topology BDT
    h_mass_topo_rs_data_1 = ROOT.TH1D("h_mass_topo_rs_data_1", "h_mass_topo_rs_data_1", 100, 4000, 9000)
    h_mass_topo_rs_data_2 = ROOT.TH1D("h_mass_topo_rs_data_2", "h_mass_rs_data_2", 100, 4000, 9000)
    h_mass_topo_rs_data_3 = ROOT.TH1D("h_mass_topo_rs_data_3", "h_mass_topo_rs_data_3", 100, 4000, 9000)
    h_mass_topo_rs_data_4 = ROOT.TH1D("h_mass_topo_rs_data_4", "h_mass_topo_rs_data_4", 100, 4000, 9000)
    h_mass_topo_rs_data_5 = ROOT.TH1D("h_mass_topo_rs_data_5", "h_mass_topo_rs_data_5", 100, 4000, 9000)

    h_mass_topo_ws_data_1 = ROOT.TH1D("h_mass_topo_ws_data_1", "h_mass_topo_ws_data_1", 100, 4000, 9000)
    h_mass_topo_ws_data_2 = ROOT.TH1D("h_mass_topo_ws_data_2", "h_mass_topo_ws_data_2", 100, 4000, 9000)
    h_mass_topo_ws_data_3 = ROOT.TH1D("h_mass_topo_ws_data_3", "h_mass_topo_ws_data_3", 100, 4000, 9000)
    h_mass_topo_ws_data_4 = ROOT.TH1D("h_mass_topo_ws_data_4", "h_mass_topo_ws_data_4", 100, 4000, 9000)
    h_mass_topo_ws_data_5 = ROOT.TH1D("h_mass_topo_ws_data_5", "h_mass_topo_ws_data_5", 100, 4000, 9000)

    t_rs_data_2016.Draw("df_Bp_M >> h_mass_topo_rs_data_1")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_topo_rs_data_2", "(BDT2 > 0.5)")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_topo_rs_data_3", "(BDT2 > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_topo_rs_data_4", "(BDT2 > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_topo_rs_data_5", "(BDT2 > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

    t_ws_data_2016.Draw("df_Bp_M >> h_mass_topo_ws_data_1")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_topo_ws_data_2", "(BDT2 > 0.5)")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_topo_ws_data_3", "(BDT2 > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_topo_ws_data_4", "(BDT2 > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_topo_ws_data_5", "(BDT2 > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

    c100_1 = ROOT.TCanvas()
    c100_1.cd()
    h_mass_topo_rs_data_1.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_topo_rs_data_1.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_topo_rs_data_1.SetTitle("Topology: BDT1 > 0.0")
    h_mass_topo_rs_data_1.SetLineColor(1)
    h_mass_topo_ws_data_1.SetLineColor(2)
    h_mass_topo_rs_data_1.DrawNormalized()
    h_mass_topo_ws_data_1.DrawNormalized("same")
    topo_leg_100_1 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_100_1.AddEntry(h_mass_topo_rs_data_1, "RS data", "lp")
    topo_leg_100_1.AddEntry(h_mass_topo_ws_data_1, "WS data", "lp")
    topo_leg_100_1.SetBorderSize(0)
    topo_leg_100_1.Draw("same")
    c100_1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_topology_BDT_0.0.pdf") 

    c100_2 = ROOT.TCanvas()
    c100_2.cd()
    h_mass_topo_ws_data_2.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_topo_ws_data_2.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_topo_ws_data_2.SetTitle("Topology: BDT1 > 0.5")
    h_mass_topo_rs_data_2.SetLineColor(1)
    h_mass_topo_ws_data_2.SetLineColor(2)
    h_mass_topo_ws_data_2.DrawNormalized()
    h_mass_topo_rs_data_2.DrawNormalized("same")
    topo_leg_200_2 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_200_2.AddEntry(h_mass_topo_rs_data_2, "RS data", "lp")
    topo_leg_200_2.AddEntry(h_mass_topo_ws_data_2, "WS data", "lp")
    topo_leg_200_2.SetBorderSize(0)
    topo_leg_200_2.Draw("same")
    c100_2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_topology_BDT_0.5.pdf") 

    c100_3 = ROOT.TCanvas()
    c100_3.cd()
    h_mass_topo_ws_data_3.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_topo_ws_data_3.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_topo_ws_data_3.SetTitle("Topology: BDT1 > 0.8")
    h_mass_topo_rs_data_3.SetLineColor(1)
    h_mass_topo_ws_data_3.SetLineColor(2)
    h_mass_topo_ws_data_3.DrawNormalized()
    h_mass_topo_rs_data_3.DrawNormalized("same")
    topo_leg_300_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_300_3.AddEntry(h_mass_topo_rs_data_3, "RS data", "lp")
    topo_leg_300_3.AddEntry(h_mass_topo_ws_data_3, "WS data", "lp")
    topo_leg_300_3.SetBorderSize(0)
    topo_leg_300_3.Draw("same")
    c100_3.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_topology_BDT_0.8.pdf") 

    c100_4 = ROOT.TCanvas()
    c100_4.cd()
    h_mass_topo_ws_data_4.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_topo_ws_data_4.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_topo_ws_data_4.SetTitle("Topology: BDT1 > 0.9")
    h_mass_topo_rs_data_4.SetLineColor(1)
    h_mass_topo_ws_data_4.SetLineColor(2)
    h_mass_topo_ws_data_4.DrawNormalized()
    h_mass_topo_rs_data_4.DrawNormalized("same")
    topo_leg_400_4 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_400_4.AddEntry(h_mass_topo_rs_data_4, "RS data", "lp")
    topo_leg_400_4.AddEntry(h_mass_topo_ws_data_4, "WS data", "lp")
    topo_leg_400_4.SetBorderSize(0)
    topo_leg_400_4.Draw("same")
    c100_4.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_topology_BDT_0.9.pdf") 

    c100_5 = ROOT.TCanvas()
    c100_5.cd()
    h_mass_topo_ws_data_5.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_topo_ws_data_5.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_topo_ws_data_5.SetTitle("Topology: BDT1 > 0.95")
    h_mass_topo_rs_data_5.SetLineColor(1)
    h_mass_topo_ws_data_5.SetLineColor(2)
    h_mass_topo_ws_data_5.DrawNormalized()
    h_mass_topo_rs_data_5.DrawNormalized("same")
    topo_leg_500_5 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_500_5.AddEntry(h_mass_topo_rs_data_5, "RS data", "lp")
    topo_leg_500_5.AddEntry(h_mass_topo_ws_data_5, "WS data", "lp")
    topo_leg_500_5.SetBorderSize(0)
    topo_leg_500_5.Draw("same")
    c100_5.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_topology_BDT_0.95.pdf") 

    if(name == "response_merged"):
        # Merged BDT
        h_mass_merged_rs_data_1 = ROOT.TH1D("h_mass_merged_rs_data_1", "h_mass_merged_rs_data_1", 100, 4000, 9000)
        h_mass_merged_rs_data_2 = ROOT.TH1D("h_mass_merged_rs_data_2", "h_mass_rs_data_2", 100, 4000, 9000)
        h_mass_merged_rs_data_3 = ROOT.TH1D("h_mass_merged_rs_data_3", "h_mass_merged_rs_data_3", 100, 4000, 9000)
        h_mass_merged_rs_data_4 = ROOT.TH1D("h_mass_merged_rs_data_4", "h_mass_merged_rs_data_4", 100, 4000, 9000)
        h_mass_merged_rs_data_5 = ROOT.TH1D("h_mass_merged_rs_data_5", "h_mass_merged_rs_data_5", 100, 4000, 9000)

        h_mass_merged_ws_data_1 = ROOT.TH1D("h_mass_merged_ws_data_1", "h_mass_merged_ws_data_1", 100, 4000, 9000)
        h_mass_merged_ws_data_2 = ROOT.TH1D("h_mass_merged_ws_data_2", "h_mass_merged_ws_data_2", 100, 4000, 9000)
        h_mass_merged_ws_data_3 = ROOT.TH1D("h_mass_merged_ws_data_3", "h_mass_merged_ws_data_3", 100, 4000, 9000)
        h_mass_merged_ws_data_4 = ROOT.TH1D("h_mass_merged_ws_data_4", "h_mass_merged_ws_data_4", 100, 4000, 9000)
        h_mass_merged_ws_data_5 = ROOT.TH1D("h_mass_merged_ws_data_5", "h_mass_merged_ws_data_5", 100, 4000, 9000)

        t_rs_data_2016.Draw("df_Bp_M >> h_mass_merged_rs_data_1")
        t_rs_data_2016.Draw("df_Bp_M >> h_mass_merged_rs_data_2", "(BDT > 0.5)")
        t_rs_data_2016.Draw("df_Bp_M >> h_mass_merged_rs_data_3", "(BDT > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_rs_data_2016.Draw("df_Bp_M >> h_mass_merged_rs_data_4", "(BDT > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_rs_data_2016.Draw("df_Bp_M >> h_mass_merged_rs_data_5", "(BDT > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

        t_ws_data_2016.Draw("df_Bp_M >> h_mass_merged_ws_data_1")
        t_ws_data_2016.Draw("df_Bp_M >> h_mass_merged_ws_data_2", "(BDT > 0.5)")
        t_ws_data_2016.Draw("df_Bp_M >> h_mass_merged_ws_data_3", "(BDT > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_ws_data_2016.Draw("df_Bp_M >> h_mass_merged_ws_data_4", "(BDT > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_ws_data_2016.Draw("df_Bp_M >> h_mass_merged_ws_data_5", "(BDT > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

        c100_1 = ROOT.TCanvas()
        c100_1.cd()
        h_mass_merged_rs_data_1.GetXaxis().SetTitle("m_{B} (MeV)")
        h_mass_merged_rs_data_1.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        h_mass_merged_rs_data_1.SetTitle("Merged: BDT1 > 0.0")
        h_mass_merged_rs_data_1.SetLineColor(1)
        h_mass_merged_ws_data_1.SetLineColor(2)
        h_mass_merged_rs_data_1.DrawNormalized()
        h_mass_merged_ws_data_1.DrawNormalized("same")
        merged_leg_100_1 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_100_1.AddEntry(h_mass_merged_rs_data_1, "RS data", "lp")
        merged_leg_100_1.AddEntry(h_mass_merged_ws_data_1, "WS data", "lp")
        merged_leg_100_1.SetBorderSize(0)
        merged_leg_100_1.Draw("same")
        c100_1.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_merged_BDT_0.0.pdf") 

        c100_2 = ROOT.TCanvas()
        c100_2.cd()
        h_mass_merged_ws_data_2.GetXaxis().SetTitle("m_{B} (MeV)")
        h_mass_merged_ws_data_2.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        h_mass_merged_ws_data_2.SetTitle("Merged: BDT1 > 0.5")
        h_mass_merged_rs_data_2.SetLineColor(1)
        h_mass_merged_ws_data_2.SetLineColor(2)
        h_mass_merged_ws_data_2.DrawNormalized()
        h_mass_merged_rs_data_2.DrawNormalized("same")
        merged_leg_200_2 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_200_2.AddEntry(h_mass_merged_rs_data_2, "RS data", "lp")
        merged_leg_200_2.AddEntry(h_mass_merged_ws_data_2, "WS data", "lp")
        merged_leg_200_2.SetBorderSize(0)
        merged_leg_200_2.Draw("same")
        c100_2.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_merged_BDT_0.5.pdf") 

        c100_3 = ROOT.TCanvas()
        c100_3.cd()
        h_mass_merged_ws_data_3.GetXaxis().SetTitle("m_{B} (MeV)")
        h_mass_merged_ws_data_3.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        h_mass_merged_ws_data_3.SetTitle("Merged: BDT1 > 0.8")
        h_mass_merged_rs_data_3.SetLineColor(1)
        h_mass_merged_ws_data_3.SetLineColor(2)
        h_mass_merged_ws_data_3.DrawNormalized()
        h_mass_merged_rs_data_3.DrawNormalized("same")
        merged_leg_300_3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_300_3.AddEntry(h_mass_merged_rs_data_3, "RS data", "lp")
        merged_leg_300_3.AddEntry(h_mass_merged_ws_data_3, "WS data", "lp")
        merged_leg_300_3.SetBorderSize(0)
        merged_leg_300_3.Draw("same")
        c100_3.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_merged_BDT_0.8.pdf") 

        c100_4 = ROOT.TCanvas()
        c100_4.cd()
        h_mass_merged_ws_data_4.GetXaxis().SetTitle("m_{B} (MeV)")
        h_mass_merged_ws_data_4.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        h_mass_merged_ws_data_4.SetTitle("Merged: BDT1 > 0.9")
        h_mass_merged_rs_data_4.SetLineColor(1)
        h_mass_merged_ws_data_4.SetLineColor(2)
        h_mass_merged_ws_data_4.DrawNormalized()
        h_mass_merged_rs_data_4.DrawNormalized("same")
        merged_leg_400_4 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_400_4.AddEntry(h_mass_merged_rs_data_4, "RS data", "lp")
        merged_leg_400_4.AddEntry(h_mass_merged_ws_data_4, "WS data", "lp")
        merged_leg_400_4.SetBorderSize(0)
        merged_leg_400_4.Draw("same")
        c100_4.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_merged_BDT_0.9.pdf") 

        c100_5 = ROOT.TCanvas()
        c100_5.cd()
        h_mass_merged_ws_data_5.GetXaxis().SetTitle("m_{B} (MeV)")
        h_mass_merged_ws_data_5.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        h_mass_merged_ws_data_5.SetTitle("Merged: BDT1 > 0.95")
        h_mass_merged_rs_data_5.SetLineColor(1)
        h_mass_merged_ws_data_5.SetLineColor(2)
        h_mass_merged_ws_data_5.DrawNormalized()
        h_mass_merged_rs_data_5.DrawNormalized("same")
        merged_leg_500_5 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_500_5.AddEntry(h_mass_merged_rs_data_5, "RS data", "lp")
        merged_leg_500_5.AddEntry(h_mass_merged_ws_data_5, "WS data", "lp")
        merged_leg_500_5.SetBorderSize(0)
        merged_leg_500_5.Draw("same")
        c100_5.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_merged_BDT_0.95.pdf") 


    # PLOTS: Signal vs background (norm mode)
    # Isolation BDT
    h_norm_mc_iso = ROOT.TH1D("h_norm_mc_iso", "h_norm_mc_iso", 30, 0, 1)
    t_norm_mc_no_rect_cuts_2016.Draw("BDT1 >> h_norm_mc_iso", "(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)")

    h_norm_data_iso = ROOT.TH1D("h_norm_data_iso", "h_norm_data_iso", 30, 0, 1)
    t_norm_data_no_rect_cuts_2016.Draw("BDT1 >> h_norm_data_iso", "(Bp_dtf_M[0] > 5320)")

    c15 = ROOT.TCanvas()
    c15.cd()
    h_norm_data_iso.GetXaxis().SetTitle("BDT1")
    h_norm_data_iso.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_norm_data_iso.SetTitle("Isolation BDT")
    h_norm_data_iso.SetLineColor(2)
    h_norm_mc_iso.SetLineColor(4)
    h_norm_data_iso.DrawNormalized("hist")
    h_norm_mc_iso.DrawNormalized("same")
    iso_leg_4 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    iso_leg_4.AddEntry(h_norm_mc_iso, "Signal", "lp")
    iso_leg_4.AddEntry(h_norm_data_iso, "Background", "lp")
    iso_leg_4.SetBorderSize(0)
    iso_leg_4.SetBorderSize(0)
    iso_leg_4.Draw("same")
    c15.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/iso_sig_vs_bkg_norm.pdf")

    # Topology BDT
    h_norm_mc_topo = ROOT.TH1D("h_norm_mc_topo", "h_norm_mc_topo", 30, 0, 1)
    t_norm_mc_no_rect_cuts_2016.Draw("BDT2 >> h_norm_mc_topo", "(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)")

    h_norm_data_topo = ROOT.TH1D("h_norm_data_topo", "h_norm_data_topo", 30, 0, 1)
    t_norm_data_no_rect_cuts_2016.Draw("BDT2 >> h_norm_data_topo", "(Bp_dtf_M[0] > 5320)")

    c16 = ROOT.TCanvas()
    c16.cd()
    h_norm_data_topo.GetXaxis().SetTitle("BDT2")
    h_norm_data_topo.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_norm_data_topo.SetTitle("Topology BDT")
    h_norm_data_topo.SetLineColor(2)
    h_norm_mc_topo.SetLineColor(4)
    h_norm_data_topo.DrawNormalized("hist")
    h_norm_mc_topo.DrawNormalized("same")
    topo_leg_4 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    topo_leg_4.AddEntry(h_norm_mc_topo, "Signal", "lp")
    topo_leg_4.AddEntry(h_norm_data_topo, "Background", "lp")
    topo_leg_4.SetBorderSize(0)
    topo_leg_4.SetBorderSize(0)
    topo_leg_4.Draw("same")
    c16.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/topo_sig_vs_bkg_norm.pdf")

    if(name == "response_merged"):
        # Merged
        h_norm_mc_merged = ROOT.TH1D("h_norm_mc_merged", "h_norm_mc_merged", 30, 0, 1)
        t_norm_mc_no_rect_cuts_2016.Draw("BDT >> h_norm_mc_merged", "(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)")

        h_norm_data_merged = ROOT.TH1D("h_norm_data_merged", "h_norm_data_merged", 30, 0, 1)
        t_norm_data_no_rect_cuts_2016.Draw("BDT >> h_norm_data_merged", "(Bp_dtf_M[0] > 5320)")

        c166 = ROOT.TCanvas()
        c166.cd()
        h_norm_data_merged.GetXaxis().SetTitle("BDT2")
        h_norm_data_merged.GetYaxis().SetTitle("Normalized entries / (0.03)")
        h_norm_data_merged.SetTitle("Merged BDT")
        h_norm_data_merged.SetLineColor(2)
        h_norm_mc_merged.SetLineColor(4)
        h_norm_data_merged.DrawNormalized("hist")
        h_norm_mc_merged.DrawNormalized("same")
        merged_leg_4 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
        merged_leg_4.AddEntry(h_norm_mc_merged, "Signal", "lp")
        merged_leg_4.AddEntry(h_norm_data_merged, "Background", "lp")
        merged_leg_4.SetBorderSize(0)
        merged_leg_4.SetBorderSize(0)
        merged_leg_4.Draw("same")
        c166.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/merged_sig_vs_bkg_norm.pdf")


    # PLOTS: BDT1 vs BDT2 (norm mode)
    h_bdt1_bdt2_norm_mc = ROOT.TH2D("h_bdt1_bdt2_norm_mc", "h_bdt1_bdt2_norm_mc", 30, 0, 1, 30, 0, 1)
    h_bdt1_bdt2_norm_data = ROOT.TH2D("h_bdt1_bdt2_norm_data", "h_bdt1_bdt2_norm_data", 30, 0, 1, 30, 0, 1)

    t_norm_mc_no_rect_cuts_2016.Draw("BDT1 : BDT2 >> h_bdt1_bdt2_norm_mc")
    t_norm_data_no_rect_cuts_2016.Draw("BDT1 : BDT2 >> h_bdt1_bdt2_norm_data")

    c17 = ROOT.TCanvas()
    c17.cd()
    c17.SetLogz()
    h_bdt1_bdt2_norm_mc.GetXaxis().SetTitle("Topology (BDT2)")
    h_bdt1_bdt2_norm_mc.GetYaxis().SetTitle("Isolation (BDT1)")
    h_bdt1_bdt2_norm_mc.SetTitle(f"Norm MC: isolation vs topology BDT (#rho = {h_bdt1_bdt2_norm_mc.GetCorrelationFactor():.2f})")
    h_bdt1_bdt2_norm_mc.Draw("COLZ")
    c17.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bdt1_vs_bdt2_norm_mc.pdf")

    c18 = ROOT.TCanvas()
    c18.cd()
    c18.SetLogz()
    h_bdt1_bdt2_norm_data.GetXaxis().SetTitle("Topology (BDT2)")
    h_bdt1_bdt2_norm_data.GetYaxis().SetTitle("Isolation (BDT1)")
    h_bdt1_bdt2_norm_data.SetTitle(f"Norm data: isolation vs topology BDT (#rho = {h_bdt1_bdt2_norm_data.GetCorrelationFactor():.2f})")
    h_bdt1_bdt2_norm_data.Draw("COLZ")
    c18.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bdt1_vs_bdt2_norm_data.pdf")

    ### PLTOS: B+ mass evolution w/ BDT cut (norm mode)
    histos_norm_mc = []
    histos_norm_data = []

    for i in range(N):
        h_mass_corr_norm_mc = ROOT.TH1D(f"h_mass_corr_norm_mc_{i}", f"h_mass_corr_norm_mc_{i}", 100, 5235, 5355)
        h_mass_corr_norm_data = ROOT.TH1D(f"h_mass_corr_norm_data_{i}", f"h_mass_corr_norm_data_{i}", 100, 5235, 5355)

        t_norm_mc_no_rect_cuts_2016.Draw(f"Bp_dtf_M[0] >> h_mass_corr_norm_mc_{i}", f"(BDT1 > {bdt_cuts[i]}) && (BDT2 > {bdt_cuts[i]})")
        t_norm_data_no_rect_cuts_2016.Draw(f"Bp_dtf_M[0] >> h_mass_corr_norm_data_{i}", f"(BDT1 > {bdt_cuts[i]}) && (BDT2 > {bdt_cuts[i]})")

        h_mass_corr_norm_mc.SetLineColor(colors[i])
        h_mass_corr_norm_data.SetLineColor(colors[i])

        histos_norm_mc.append(h_mass_corr_norm_mc)
        histos_norm_data.append(h_mass_corr_norm_data)

    c19 = ROOT.TCanvas()
    c19.cd()
    leg19 = ROOT.TLegend(0.7, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_norm_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_norm_mc[j].GetYaxis().SetTitle("Normalized entries / (40 MeV)")
        histos_norm_mc[j].SetTitle("Norm MC: B mass vs (BDT1,BDT2) cut")
        leg19.AddEntry(histos_norm_mc[j], f"(BDT1,BDT2) > ({bdt_cuts[j]},{bdt_cuts[j]})")
        if(i == 0):
            histos_norm_mc[j].DrawNormalized()
        else:
            histos_norm_mc[j].DrawNormalized("same")
    leg19.SetBorderSize(0)
    leg19.Draw("same")
    c19.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_norm_mc.pdf")

    c20 = ROOT.TCanvas()
    c20.cd()
    leg20 = ROOT.TLegend(0.7, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_norm_data[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_norm_data[j].GetYaxis().SetTitle("Normalized entries / (40 MeV)")
        histos_norm_data[j].SetTitle("Norm data: B mass vs (BDT1,BDT2) cut")
        leg20.AddEntry(histos_norm_data[j], f"(BDT1,BDT2) > ({bdt_cuts[j]},{bdt_cuts[j]})")
        if(i == 0):
            histos_norm_data[j].DrawNormalized()
        else:
            histos_norm_data[j].DrawNormalized("same")
    leg20.SetBorderSize(0)
    leg20.Draw("same")
    c20.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_norm_data.pdf")


if __name__ == "__main__":
    main(sys.argv)