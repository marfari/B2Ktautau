import sys
import ROOT


def main(argv):
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

    t_BuKtautau_mc_2016.AddFriend(t1_BuKtautau_mc_2016)
    t_BuKtautau_mc_2016.AddFriend(t2_BuKtautau_mc_2016)
    t_BuKtautau_mc_2016.AddFriend(t3_BuKtautau_mc_2016)

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

    t_rs_data_2016.AddFriend(t1_rs_data_2016)
    t_rs_data_2016.AddFriend(t2_rs_data_2016)
    t_rs_data_2016.AddFriend(t3_rs_data_2016)

    t_ws_data_2016.AddFriend(t1_ws_data_2016)
    t_ws_data_2016.AddFriend(t2_ws_data_2016)
    t_ws_data_2016.AddFriend(t3_ws_data_2016)

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
    t_norm_mc_no_rect_cuts_2016.AddFriend(t2_norm_mc_no_rect_cuts_2016)

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
    t_norm_data_no_rect_cuts_2016.AddFriend(t2_norm_data_no_rect_cuts_2016)

    ######################################################################################################################################################


    ### PLOTS: 3 MC components
    h_3pi3pi_mc = ROOT.TH1D("h_3pi3pi_mc", "h_3pi3pi_mc", 30, 0, 1)
    h_3pi3pipi0_mc = ROOT.TH1D("h_3pi3pipi0_mc", "h_3pi3pipi0_mc", 30, 0, 1)
    h_3pi3pi2pi0_mc = ROOT.TH1D("h_3pi3pi2pi0_mc", "h_3pi3pi2pi0_mc", 30, 0, 1)

    t_BuKtautau_mc_2016.Draw("BDT >> h_3pi3pi_mc", cut_Ktautau_mc+" && (component==0)")
    t_BuKtautau_mc_2016.Draw("BDT >> h_3pi3pipi0_mc", cut_Ktautau_mc+" && (component==1)")
    t_BuKtautau_mc_2016.Draw("BDT >> h_3pi3pi2pi0_mc", cut_Ktautau_mc+" && (component==2)")

    c = ROOT.TCanvas()
    c.cd()
    h_3pi3pi_mc.GetXaxis().SetTitle("BDT")
    h_3pi3pi_mc.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_3pi3pi_mc.SetTitle("")
    h_3pi3pi_mc.SetLineColor(4)
    h_3pi3pipi0_mc.SetLineColor(8)
    h_3pi3pi2pi0_mc.SetLineColor(2)
    h_3pi3pi_mc.DrawNormalized()
    h_3pi3pipi0_mc.DrawNormalized("same")
    h_3pi3pi2pi0_mc.DrawNormalized("same")
    leg = ROOT.TLegend(0.1, 0.7, 0.4, 0.89)
    leg.AddEntry(h_3pi3pi_mc, "3#pi3#pi MC", "lp")
    leg.AddEntry(h_3pi3pipi0_mc, "3#pi3#pi #pi^{0} MC", "lp")
    leg.AddEntry(h_3pi3pi2pi0_mc, "3#pi3#pi 2#pi^{0} MC", "lp")
    leg.SetBorderSize(0)
    leg.Draw("same")
    c.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/ktautau_mc_components.pdf")

    ### PLOTS: cocktail MCs
    h_BuDDKp_mc = ROOT.TH1D("h_BuDDKp_mc", "h_BuDDKp_mc", 30, 0, 1)
    h_BdDDKp_mc = ROOT.TH1D("h_BdDDKp_mc", "h_BdDDKp_mc", 30, 0, 1)
    h_BsDDKp_mc = ROOT.TH1D("h_BsDDKp_mc", "h_BsDDKp_mc", 30, 0, 1)
    h_BuDDK0_mc = ROOT.TH1D("h_BuDDK0_mc", "h_BuDDK0_mc", 30, 0, 1)
    h_BuDD_mc = ROOT.TH1D("h_BuDD_mc", "h_BuDD_mc", 30, 0, 1)

    t_BuDDKp_mc_2016.Draw("BDT >> h_BuDDKp_mc", cut_cocktail_mc)
    t_BdDDKp_mc_2016.Draw("BDT >> h_BdDDKp_mc", cut_cocktail_mc)
    t_BsDDKp_mc_2016.Draw("BDT >> h_BsDDKp_mc", cut_cocktail_mc)
    t_BuDDK0_mc_2016.Draw("BDT >> h_BuDDK0_mc", cut_cocktail_mc)
    t_BuDD_mc_2016.Draw("BDT >> h_BuDD_mc", cut_cocktail_mc)

    c1 = ROOT.TCanvas()
    c1.cd()
    h_BuDD_mc.GetXaxis().SetTitle("BDT")
    h_BuDD_mc.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_BuDD_mc.SetTitle("")
    h_BuDDKp_mc.SetLineColor(4)
    h_BdDDKp_mc.SetLineColor(8)
    h_BsDDKp_mc.SetLineColor(2)
    h_BuDDK0_mc.SetLineColor(51)
    h_BuDD_mc.SetLineColor(94)
    h_BuDD_mc.DrawNormalized()
    h_BdDDKp_mc.DrawNormalized("same")
    h_BsDDKp_mc.DrawNormalized("same")
    h_BuDDK0_mc.DrawNormalized("same")
    h_BuDDKp_mc.DrawNormalized("same")
    leg1 = ROOT.TLegend(0.4, 0.6, 0.7, 0.89)
    leg1.AddEntry(h_BuDDKp_mc, "B^{+} #rightarrow D D K^{+}", "lp")
    leg1.AddEntry(h_BdDDKp_mc, "B^{0} #rightarrow D D K^{+}", "lp")
    leg1.AddEntry(h_BsDDKp_mc, "B^{0}_{s} #rightarrow D D K^{+}", "lp")
    leg1.AddEntry(h_BuDDK0_mc, "B^{+} #rightarrow D D K^{0}", "lp")
    leg1.AddEntry(h_BuDD_mc, "B^{+} #rightarrow D D ", "lp")
    leg1.SetBorderSize(0)
    leg1.Draw("same")
    c1.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/cocktail_mcs.pdf")

    # PLOTS: RS and WS data
    h_rs_data = ROOT.TH1D("h_rs_data", "h_rs_data", 30, 0, 1)
    h_ws_data = ROOT.TH1D("h_ws_data", "h_ws_data", 30, 0, 1)

    t_rs_data_2016.Draw("BDT >> h_rs_data", cut_Ktautau_data)
    t_ws_data_2016.Draw("BDT >> h_ws_data", cut_Ktautau_data)

    c2 = ROOT.TCanvas()
    c2.cd()
    h_rs_data.GetXaxis().SetTitle("BDT")
    h_rs_data.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_rs_data.SetTitle("")
    h_rs_data.SetLineColor(4)
    h_ws_data.SetLineColor(8)
    h_rs_data.DrawNormalized()
    h_ws_data.DrawNormalized("same")
    leg2 = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
    leg2.AddEntry(h_rs_data, "RS data", "lp")
    leg2.AddEntry(h_ws_data, "WS data", "lp")
    leg2.SetBorderSize(0)
    leg2.Draw("same")
    c2.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/ktautau_data.pdf")

    # PLOTS: Signal vs background
    h_ktautau_mc = ROOT.TH1D("h_ktautau_mc", "h_ktautau_mc", 30, 0, 1)
    t_BuKtautau_mc_2016.Draw("BDT >> h_ktautau_mc", cut_Ktautau_mc)

    c3 = ROOT.TCanvas()
    c3.cd()
    h_ktautau_mc.GetXaxis().SetTitle("BDT")
    h_ktautau_mc.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_ktautau_mc.SetTitle("")
    h_ktautau_mc.SetLineColor(4)
    h_ws_data.SetLineColor(2)
    h_ktautau_mc.DrawNormalized()
    h_ws_data.DrawNormalized("same")
    leg3 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg3.AddEntry(h_ktautau_mc, "Signal", "lp")
    leg3.AddEntry(h_ws_data, "Background", "lp")
    leg3.SetBorderSize(0)
    leg3.Draw("same")
    c3.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/sig_vs_bkg.pdf")

    ### PLTOS: B+ mass evolution w/ BDT cut
    bdt_cuts = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9 ]
    colors = [1, 4, 6, 2, 8, 51]
    N = len(bdt_cuts)

    histos_bdt_ktautau_mc = []
    histos_bdt_cocktail_mc = []
    histos_bdt_rs_data = []
    histos_bdt_ws_data = []

    for i in range(N):
        h_mass_bdt_corr_ktautau_mc = ROOT.TH1D(f"h_mass_bdt_corr_ktautau_mc_{i}", f"h_mass_bdt_corr_ktautau_mc_{i}", 100, 4000, 9000)
        h_mass_bdt_corr_cocktail_mc = ROOT.TH1D(f"h_mass_bdt_corr_cocktail_mc_{i}", f"h_mass_bdt_corr_cocktail_mc_{i}", 100, 4000, 9000)
        h_mass_bdt_corr_rs_data = ROOT.TH1D(f"h_mass_bdt_corr_rs_data_{i}", f"h_mass_bdt_corr_rs_data_{i}", 100, 4000, 9000)
        h_mass_bdt_corr_ws_data = ROOT.TH1D(f"h_mass_bdt_corr_ws_data_{i}", f"h_mass_bdt_corr_ws_data_{i}", 100, 4000, 9000)

        t_BuKtautau_mc_2016.Draw(f"df_Bp_M >> h_mass_bdt_corr_ktautau_mc_{i}", cut_Ktautau_mc+" && "+f"(BDT > {bdt_cuts[i]})")
        t_cocktail_mc.Draw(f"df_Bp_M >> h_mass_bdt_corr_cocktail_mc_{i}", cut_cocktail_mc+" && "+f"(BDT > {bdt_cuts[i]})")
        t_rs_data_2016.Draw(f"df_Bp_M >> h_mass_bdt_corr_rs_data_{i}", cut_Ktautau_data+" && "+f"(BDT > {bdt_cuts[i]}) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
        t_ws_data_2016.Draw(f"df_Bp_M >> h_mass_bdt_corr_ws_data_{i}", cut_Ktautau_data+" && "+f"(BDT > {bdt_cuts[i]})")

        h_mass_bdt_corr_ktautau_mc.SetLineColor(colors[i])
        h_mass_bdt_corr_cocktail_mc.SetLineColor(colors[i])
        h_mass_bdt_corr_rs_data.SetLineColor(colors[i])
        h_mass_bdt_corr_ws_data.SetLineColor(colors[i])

        histos_bdt_ktautau_mc.append(h_mass_bdt_corr_ktautau_mc)
        histos_bdt_cocktail_mc.append(h_mass_bdt_corr_cocktail_mc)
        histos_bdt_rs_data.append(h_mass_bdt_corr_rs_data)
        histos_bdt_ws_data.append(h_mass_bdt_corr_ws_data)

    c4 = ROOT.TCanvas()
    c4.cd()
    leg4 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt_ktautau_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt_ktautau_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt_ktautau_mc[j].SetTitle("Ktautau MC: B mass vs BDT")
        leg4.AddEntry(histos_bdt_ktautau_mc[j], f"BDT > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt_ktautau_mc[j].DrawNormalized()
        else:
            histos_bdt_ktautau_mc[j].DrawNormalized("same")
    leg4.SetBorderSize(0)
    leg4.SetTextSize(0.03)
    leg4.Draw("same")
    c4.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_ktautau_mc.pdf")

    c5 = ROOT.TCanvas()
    c5.cd()
    leg5 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt_cocktail_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt_cocktail_mc[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt_cocktail_mc[j].SetTitle("Cocktail MCs: B mass vs BDT")
        leg5.AddEntry(histos_bdt_cocktail_mc[j], f"BDT > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt_cocktail_mc[j].DrawNormalized()
        else:
            histos_bdt_cocktail_mc[j].DrawNormalized("same")
    leg5.SetBorderSize(0)
    leg5.SetTextSize(0.03)
    leg5.Draw("same")
    c5.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_cocktail_mc.pdf")

    c6 = ROOT.TCanvas()
    c6.cd()
    leg6 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_bdt_rs_data[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt_rs_data[j].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt_rs_data[j].SetTitle("RS data: B mass vs BDT")
        leg6.AddEntry(histos_bdt_rs_data[j], f"BDT > {bdt_cuts[j]}", "fp")
        if(i == 0):
            histos_bdt_rs_data[j].DrawNormalized()
        else:
            histos_bdt_rs_data[j].DrawNormalized("same")
    leg6.SetBorderSize(0)
    leg6.SetTextSize(0.03)
    leg6.Draw("same")
    c6.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_rs_data.pdf")

    c7 = ROOT.TCanvas()
    c7.cd()
    leg7 = ROOT.TLegend(0.5, 0.5, 0.89, 0.89)
    for i in range(N):
        histos_bdt_ws_data[i].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_bdt_ws_data[i].GetYaxis().SetTitle("Normalized entries / (50 MeV)")
        histos_bdt_ws_data[i].SetTitle("WS data: B mass vs BDT")
        leg7.AddEntry(histos_bdt_ws_data[i], f"BDT > {bdt_cuts[i]}", "fp")
        if(i == 0):
            histos_bdt_ws_data[i].DrawNormalized()
        else:
            histos_bdt_ws_data[i].DrawNormalized("same")
    leg7.SetBorderSize(0)
    leg7.SetTextSize(0.03)
    leg7.Draw("same")
    c7.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_bdt_corr_ws_data.pdf")

    # RS vs WS mass at different BDT cuts (0, 0.5, 0.8, 0.9, 0.95)
    h_mass_rs_data_1 = ROOT.TH1D("h_mass_rs_data_1", "h_mass_rs_data_1", 100, 4000, 9000)
    h_mass_rs_data_2 = ROOT.TH1D("h_mass_rs_data_2", "h_mass_rs_data_2", 100, 4000, 9000)
    h_mass_rs_data_3 = ROOT.TH1D("h_mass_rs_data_3", "h_mass_rs_data_3", 100, 4000, 9000)
    h_mass_rs_data_4 = ROOT.TH1D("h_mass_rs_data_4", "h_mass_rs_data_4", 100, 4000, 9000)
    h_mass_rs_data_5 = ROOT.TH1D("h_mass_rs_data_5", "h_mass_rs_data_5", 100, 4000, 9000)

    h_mass_ws_data_1 = ROOT.TH1D("h_mass_ws_data_1", "h_mass_ws_data_1", 100, 4000, 9000)
    h_mass_ws_data_2 = ROOT.TH1D("h_mass_ws_data_2", "h_mass_ws_data_2", 100, 4000, 9000)
    h_mass_ws_data_3 = ROOT.TH1D("h_mass_ws_data_3", "h_mass_ws_data_3", 100, 4000, 9000)
    h_mass_ws_data_4 = ROOT.TH1D("h_mass_ws_data_4", "h_mass_ws_data_4", 100, 4000, 9000)
    h_mass_ws_data_5 = ROOT.TH1D("h_mass_ws_data_5", "h_mass_ws_data_5", 100, 4000, 9000)

    t_rs_data_2016.Draw("df_Bp_M >> h_mass_rs_data_1")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_rs_data_2", "(BDT > 0.5)")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_rs_data_3", "(BDT > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_rs_data_4", "(BDT > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_rs_data_2016.Draw("df_Bp_M >> h_mass_rs_data_5", "(BDT > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

    t_ws_data_2016.Draw("df_Bp_M >> h_mass_ws_data_1")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_ws_data_2", "(BDT > 0.5)")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_ws_data_3", "(BDT > 0.8) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_ws_data_4", "(BDT > 0.9) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")
    t_ws_data_2016.Draw("df_Bp_M >> h_mass_ws_data_5", "(BDT > 0.95) && ((df_Bp_M < 5000) || (df_Bp_M > 6500))")

    c8 = ROOT.TCanvas()
    c8.cd()
    h_mass_rs_data_1.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_rs_data_1.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_rs_data_1.SetTitle("BDT > 0.0")
    h_mass_rs_data_1.SetLineColor(1)
    h_mass_ws_data_1.SetLineColor(2)
    h_mass_rs_data_1.DrawNormalized()
    h_mass_ws_data_1.DrawNormalized("same")
    leg8 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg8.AddEntry(h_mass_rs_data_1, "RS data", "lp")
    leg8.AddEntry(h_mass_ws_data_1, "WS data", "lp")
    leg8.SetBorderSize(0)
    leg8.Draw("same")
    c8.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_BDT_0.0.pdf") 

    c9 = ROOT.TCanvas()
    c9.cd()
    h_mass_ws_data_2.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_ws_data_2.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_ws_data_2.SetTitle("BDT > 0.5")
    h_mass_rs_data_2.SetLineColor(1)
    h_mass_ws_data_2.SetLineColor(2)
    h_mass_ws_data_2.DrawNormalized()
    h_mass_rs_data_2.DrawNormalized("same")
    leg9 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg9.AddEntry(h_mass_rs_data_2, "RS data", "lp")
    leg9.AddEntry(h_mass_ws_data_2, "WS data", "lp")
    leg9.SetBorderSize(0)
    leg9.Draw("same")
    c9.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_BDT_0.5.pdf") 

    c10 = ROOT.TCanvas()
    c10.cd()
    h_mass_ws_data_3.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_ws_data_3.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_ws_data_3.SetTitle("BDT > 0.8")
    h_mass_rs_data_3.SetLineColor(1)
    h_mass_ws_data_3.SetLineColor(2)
    h_mass_ws_data_3.DrawNormalized()
    h_mass_rs_data_3.DrawNormalized("same")
    leg10 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg10.AddEntry(h_mass_rs_data_3, "RS data", "lp")
    leg10.AddEntry(h_mass_ws_data_3, "WS data", "lp")
    leg10.SetBorderSize(0)
    leg10.Draw("same")
    c10.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_BDT_0.8.pdf") 

    c11 = ROOT.TCanvas()
    c11.cd()
    h_mass_ws_data_4.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_ws_data_4.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_ws_data_4.SetTitle("BDT > 0.9")
    h_mass_rs_data_4.SetLineColor(1)
    h_mass_ws_data_4.SetLineColor(2)
    h_mass_ws_data_4.DrawNormalized()
    h_mass_rs_data_4.DrawNormalized("same")
    leg11 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg11.AddEntry(h_mass_rs_data_4, "RS data", "lp")
    leg11.AddEntry(h_mass_ws_data_4, "WS data", "lp")
    leg11.SetBorderSize(0)
    leg11.Draw("same")
    c11.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_BDT_0.9.pdf") 

    c12 = ROOT.TCanvas()
    c12.cd()
    h_mass_ws_data_5.GetXaxis().SetTitle("m_{B} (MeV)")
    h_mass_ws_data_5.GetYaxis().SetTitle("Normalized entries / (50 MeV)")
    h_mass_ws_data_5.SetTitle("BDT > 0.95")
    h_mass_rs_data_5.SetLineColor(1)
    h_mass_ws_data_5.SetLineColor(2)
    h_mass_ws_data_5.DrawNormalized()
    h_mass_rs_data_5.DrawNormalized("same")
    leg12 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg12.AddEntry(h_mass_rs_data_5, "RS data", "lp")
    leg12.AddEntry(h_mass_ws_data_5, "WS data", "lp")
    leg12.SetBorderSize(0)
    leg12.Draw("same")
    c12.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/rs_vs_ws_data_mass_BDT_0.95.pdf") 


    # PLOTS: Signal vs background (norm mode)
    h_norm_mc = ROOT.TH1D("h_norm_mc", "h_norm_mc", 30, 0, 1)
    t_norm_mc_no_rect_cuts_2016.Draw("BDT >> h_norm_mc", "(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)")

    h_norm_data = ROOT.TH1D("h_norm_data", "h_norm_data", 30, 0, 1)
    t_norm_data_no_rect_cuts_2016.Draw("BDT >> h_norm_data", "(Bp_dtf_M[0] > 5320)")

    c13 = ROOT.TCanvas()
    c13.cd()
    h_norm_data.GetXaxis().SetTitle("BDT2")
    h_norm_data.GetYaxis().SetTitle("Normalized entries / (0.03)")
    h_norm_data.SetTitle("")
    h_norm_data.SetLineColor(2)
    h_norm_mc.SetLineColor(4)
    h_norm_data.DrawNormalized("hist")
    h_norm_mc.DrawNormalized("same")
    leg13 = ROOT.TLegend(0.4, 0.8, 0.7, 0.89)
    leg13.AddEntry(h_norm_mc, "Signal", "lp")
    leg13.AddEntry(h_norm_data, "Background", "lp")
    leg13.SetBorderSize(0)
    leg13.SetBorderSize(0)
    leg13.Draw("same")
    c13.SaveAs(f"/panfs/felician/B2Ktautau/workflow/sklearn_plots/sig_vs_bkg_norm.pdf")

    ### PLTOS: B+ mass evolution w/ BDT cut (norm mode)
    histos_norm_mc = []
    histos_norm_data = []

    for i in range(N):
        h_mass_corr_norm_mc = ROOT.TH1D(f"h_mass_corr_norm_mc_{i}", f"h_mass_corr_norm_mc_{i}", 100, 5235, 5355)
        h_mass_corr_norm_data = ROOT.TH1D(f"h_mass_corr_norm_data_{i}", f"h_mass_corr_norm_data_{i}", 100, 5235, 5355)

        t_norm_mc_no_rect_cuts_2016.Draw(f"Bp_dtf_M[0] >> h_mass_corr_norm_mc_{i}", f"(BDT > {bdt_cuts[i]})")
        t_norm_data_no_rect_cuts_2016.Draw(f"Bp_dtf_M[0] >> h_mass_corr_norm_data_{i}", f"(BDT > {bdt_cuts[i]})")

        h_mass_corr_norm_mc.SetLineColor(colors[i])
        h_mass_corr_norm_data.SetLineColor(colors[i])

        histos_norm_mc.append(h_mass_corr_norm_mc)
        histos_norm_data.append(h_mass_corr_norm_data)

    c14 = ROOT.TCanvas()
    c14.cd()
    leg14 = ROOT.TLegend(0.7, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_norm_mc[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_norm_mc[j].GetYaxis().SetTitle("Normalized entries / (40 MeV)")
        histos_norm_mc[j].SetTitle("Norm MC: B mass vs BDT cut")
        leg14.AddEntry(histos_norm_mc[j], f"BDT > {bdt_cuts[j]}")
        if(i == 0):
            histos_norm_mc[j].DrawNormalized()
        else:
            histos_norm_mc[j].DrawNormalized("same")
    leg14.SetBorderSize(0)
    leg14.Draw("same")
    c14.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_norm_mc.pdf")

    c15 = ROOT.TCanvas()
    c15.cd()
    leg15 = ROOT.TLegend(0.7, 0.5, 0.89, 0.89)
    for i in range(N):
        j = N-i-1
        histos_norm_data[j].GetXaxis().SetTitle("m_{B} (MeV)")
        histos_norm_data[j].GetYaxis().SetTitle("Normalized entries / (40 MeV)")
        histos_norm_data[j].SetTitle("Norm data: B mass vs BDT cut")
        leg15.AddEntry(histos_norm_data[j], f"BDT > {bdt_cuts[j]}")
        if(i == 0):
            histos_norm_data[j].DrawNormalized()
        else:
            histos_norm_data[j].DrawNormalized("same")
    leg15.SetBorderSize(0)
    leg15.Draw("same")
    c15.SaveAs("/panfs/felician/B2Ktautau/workflow/sklearn_plots/bmass_corr_norm_data.pdf")


if __name__ == "__main__":
    main(sys.argv)