import sys
import ROOT
import numpy as np

def main(argv):

    # branching fraction term
    B_tau_3pi_nu = 9.31/100; # +/- 0.05/100
    B_tau_3pi_pi0_nu = 4.62/100; # +/- 0.05/100
    B_Bp_D0bar_Dsp = 9.0/1000; # +/- 0.9/1000
    B_D0bar_Kp_pi = 3.947/100; # +/- 0.030/100
    B_Dsp_Kp_Km_pi = 5.37/100; #  +/- 0.10/100

    branching_fraction_term = (B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)/((B_tau_3pi_nu+B_tau_3pi_pi0_nu)**2)

    # normalisation channel term
    ### Norm MC after all selections (for N_norm_num)
    f_norm_mc = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_72/post_sel_tree_bdt1_0_bdt2_0.root")
    t_norm_mc = f_norm_mc.Get("DecayTree")

    ### Gen norm MC (for N_norm_den)
    fc_gen_norm_2016 = ROOT.TFileCollection("fc_gen_norm_2016", "fc_gen_norm_2016", "Files_on_grid/MC_D0Dps_2016.txt")
    fc_gen_norm_2017 = ROOT.TFileCollection("fc_gen_norm_2017", "fc_gen_norm_2017", "Files_on_grid/MC_D0Dps_2017.txt")
    fc_gen_norm_2018 = ROOT.TFileCollection("fc_gen_norm_2018", "fc_gen_norm_2018", "Files_on_grid/MC_D0Dps_2018.txt")

    t_norm_gen = ROOT.TChain("mc_ntuple/MCDecayTree")
    t_norm_gen_2017 = ROOT.TChain("mc_ntuple/MCDecayTree")
    t_norm_gen_2018 = ROOT.TChain("mc_ntuple/MCDecayTree")

    t_norm_gen.AddFileInfoList(fc_gen_norm_2016.GetList())
    t_norm_gen_2017.AddFileInfoList(fc_gen_norm_2017.GetList())
    t_norm_gen_2018.AddFileInfoList(fc_gen_norm_2018.GetList())

    print("Norm")
    print(t_norm_gen.GetEntries())
    print(t_norm_gen_2017.GetEntries())
    print(t_norm_gen_2018.GetEntries())

    t_norm_gen.Add(t_norm_gen_2017)
    t_norm_gen.Add(t_norm_gen_2018)

    ### Norm acceptance
    eps_acc_norm_2016 = 14.50/100
    eps_acc_norm_2017 = 14.48/100
    eps_acc_norm_2018 = 14.58/100
    eps_acc_norm = (10.8/34.2)*eps_acc_norm_2016 + (12./34.2)*eps_acc_norm_2017 + (11.4/34.2)*eps_acc_norm_2018

    ### eps_norm
    N_norm_num = t_norm_mc.GetEntries()
    N_norm_den = t_norm_gen.GetEntries()
    eps_norm = eps_acc_norm*(N_norm_num/N_norm_den)

    print(eps_norm*100)

    ### N_norm
    N_norm = 32858

    normalisation_channel_term = eps_norm/N_norm

    print(normalisation_channel_term)

    # signal efficiency term (all MC)
    ### Ktautau gen MC
    fc_gen_2016 = ROOT.TFileCollection("fc_gen_2016", "fc_gen_2016", "Files_on_grid/MC_2016.txt")
    fc_gen_2017 = ROOT.TFileCollection("fc_gen_2017", "fc_gen_2017", "Files_on_grid/MC_2017.txt")
    fc_gen_2018 = ROOT.TFileCollection("fc_gen_2018", "fc_gen_2018", "Files_on_grid/MC_2018.txt")

    #### 2016
    t_gen_3pi3pi = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
    t_gen_3pi3pi_pi0_2016 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
    t_gen_3pipi0_3pi_2016 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
    t_gen_3pi3pi2pi0_2016 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")
    t_gen_all = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")

    t_gen_3pi3pi.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_3pi3pi_pi0_2016.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_3pipi0_3pi_2016.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_3pi3pi2pi0_2016.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_all.AddFileInfoList(fc_gen_2016.GetList())

    print("2016")
    print(t_gen_3pi3pi.GetEntries())
    print(t_gen_3pi3pi_pi0_2016.GetEntries())
    print(t_gen_3pipi0_3pi_2016.GetEntries())
    print(t_gen_3pi3pi2pi0_2016.GetEntries())
    print(t_gen_all.GetEntries())

    t_gen_all.Add(t_gen_3pi3pi_pi0_2016)
    t_gen_all.Add(t_gen_3pipi0_3pi_2016)
    t_gen_all.Add(t_gen_3pi3pi2pi0_2016)

    #### 2017
    t_gen_3pi3pi_2017 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
    t_gen_3pi3pi_pi0_2017 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
    t_gen_3pipi0_3pi_2017 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
    t_gen_3pi3pi2pi0_2017 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")
    t_gen_all_2017 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")

    t_gen_3pi3pi_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_3pi3pi_pi0_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_3pipi0_3pi_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_3pi3pi2pi0_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_all_2017.AddFileInfoList(fc_gen_2017.GetList())

    print("2017")
    print(t_gen_3pi3pi_2017.GetEntries()) 
    print(t_gen_3pi3pi_pi0_2017.GetEntries())
    print(t_gen_3pipi0_3pi_2017.GetEntries())
    print(t_gen_3pi3pi2pi0_2017.GetEntries())
    print(t_gen_all_2017.GetEntries())

    t_gen_all_2017.Add(t_gen_3pi3pi_pi0_2017)
    t_gen_all_2017.Add(t_gen_3pipi0_3pi_2017)
    t_gen_all_2017.Add(t_gen_3pi3pi2pi0_2017)

    #### 2018
    t_gen_3pi3pi_2018 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
    t_gen_3pi3pi_pi0_2018 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
    t_gen_3pipi0_3pi_2018 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
    t_gen_3pi3pi2pi0_2018 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")
    t_gen_all_2018 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")

    t_gen_3pi3pi_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_3pi3pi_pi0_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_3pipi0_3pi_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_3pi3pi2pi0_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_all_2018.AddFileInfoList(fc_gen_2018.GetList())

    print("2018")
    print(t_gen_3pi3pi_2018.GetEntries())
    print(t_gen_3pi3pi_pi0_2018.GetEntries())
    print(t_gen_3pipi0_3pi_2018.GetEntries())
    print(t_gen_3pi3pi2pi0_2018.GetEntries())
    print(t_gen_all_2018.GetEntries())

    t_gen_all_2018.Add(t_gen_3pi3pi_pi0_2018)
    t_gen_all_2018.Add(t_gen_3pipi0_3pi_2018)
    t_gen_all_2018.Add(t_gen_3pi3pi2pi0_2018)

    t_gen_3pi3pi.Add(t_gen_3pi3pi_2017)
    t_gen_3pi3pi.Add(t_gen_3pi3pi_2018)

    t_gen_all.Add(t_gen_all_2017)
    t_gen_all.Add(t_gen_all_2018)

    ### Ktautau acceptance
    eps_acc_2016 = (5.318/100)
    eps_acc_2017 = (5.321/100)
    eps_acc_2018 = (5.330/100)
    eps_acc = (28.5/165)*eps_acc_2016 + (60.0/165)*eps_acc_2017 + (76.5/165)*eps_acc_2018

    ### Ktautau stripping 
    eps_strip_2016 = (1.028/100)
    eps_strip_2017 = (1.056/100)
    eps_strip_2018 = (1.050/100)
    eps_strip = (28.5/165)*eps_strip_2016 + (60.0/165)*eps_strip_2017 + (76.5/165)*eps_strip_2018

    ### Ktautua reco1 = 100%
    N_den = t_gen_all.GetEntries()

    signal_efficiency_term = N_den/(eps_acc*eps_strip)

    print(signal_efficiency_term)

    fixed_term = branching_fraction_term*normalisation_channel_term*signal_efficiency_term

    np.save("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy", fixed_term)

if __name__ == "__main__":
    main(sys.argv)