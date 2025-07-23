import sys
import ROOT
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy

def main(argv):

    # branching fraction term
    B_tau_3pi_nu = ufloat(9.31/100, 0.05/100)
    B_tau_3pi_pi0_nu = ufloat(4.62/100, 0.05/100)
    B_Bp_D0bar_Dsp = ufloat(9.0/1000, 0.9/1000)
    B_D0bar_Kp_pi = ufloat(3.947/100, 0.030/100)
    B_Dsp_Kp_Km_pi = ufloat(5.37/100, 0.10/100)

    # normalisation channel term
    ### Norm MC after all selections (for N_norm_num)
    f_norm_mc = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_72/post_sel_tree_bdt_0.0.root")
    t_norm_mc = f_norm_mc.Get("DecayTree")

    ### Gen norm MC (for N_norm_den)
    # fc_gen_norm_2016 = ROOT.TFileCollection("fc_gen_norm_2016", "fc_gen_norm_2016", "Files_on_grid/MC_D0Dps_2016.txt")
    # fc_gen_norm_2017 = ROOT.TFileCollection("fc_gen_norm_2017", "fc_gen_norm_2017", "Files_on_grid/MC_D0Dps_2017.txt")
    # fc_gen_norm_2018 = ROOT.TFileCollection("fc_gen_norm_2018", "fc_gen_norm_2018", "Files_on_grid/MC_D0Dps_2018.txt")

    # t_norm_gen = ROOT.TChain("mc_ntuple/MCDecayTree")
    # t_norm_gen_2017 = ROOT.TChain("mc_ntuple/MCDecayTree")
    # t_norm_gen_2018 = ROOT.TChain("mc_ntuple/MCDecayTree")

    # t_norm_gen.AddFileInfoList(fc_gen_norm_2016.GetList())
    # t_norm_gen_2017.AddFileInfoList(fc_gen_norm_2017.GetList())
    # t_norm_gen_2018.AddFileInfoList(fc_gen_norm_2018.GetList())

    # print("Norm")
    # print(t_norm_gen.GetEntries())
    # print(t_norm_gen_2017.GetEntries())
    # print(t_norm_gen_2018.GetEntries())

    # t_norm_gen.Add(t_norm_gen_2017)
    # t_norm_gen.Add(t_norm_gen_2018)

    ### Norm acceptance
    eps_acc_norm_2016 = ufloat(14.50/100, 0.10/100)
    eps_acc_norm_2017 = ufloat(14.48/100, 0.10/100)
    eps_acc_norm_2018 = ufloat(14.58/100, 0.10/100)
    eps_acc_norm = (10.8/34.2)*eps_acc_norm_2016 + (12./34.2)*eps_acc_norm_2017 + (11.4/34.2)*eps_acc_norm_2018

    ### eps_norm
    N_norm_num = t_norm_mc.GetEntries()
    N_norm_den_2016 = 5350785
    N_norm_den_2017 = 5782016
    N_norm_den_2018 = 5763580
    N_norm_den = N_norm_den_2016 + N_norm_den_2017 + N_norm_den_2018
    # N_norm_den = t_norm_gen.GetEntries()

    upper_bound = ROOT.TEfficiency.Wilson(N_norm_den, N_norm_num, 0.68, True)
    lower_bound = ROOT.TEfficiency.Wilson(N_norm_den, N_norm_num, 0.68, False)
    eps_norm_err = 0.5*(upper_bound - lower_bound)

    eps_norm_post_acc = ufloat( N_norm_num/N_norm_den, eps_norm_err )

    eps_norm = eps_acc_norm*eps_norm_post_acc
    print("eps_norm = ", eps_norm)

    ### N_norm
    f_norm_fit = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_8/mass_fit_result.root")
    norm_fit_result = f_norm_fit.Get("fitresult_model_data")
    n_norm_fit_var = norm_fit_result.floatParsFinal().find("n_signal")
    N_norm_value = n_norm_fit_var.getVal()
    N_norm_error = n_norm_fit_var.getError()

    print("N_norm_value = ", N_norm_value)
    print("N_norm_error = ", N_norm_error)
    
    N_norm = ufloat(N_norm_value, N_norm_error)
    print("N_norm = ", N_norm)

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
    eps_acc_2016 = ufloat(5.318/100, 0.010/100)
    eps_acc_2017 = ufloat(5.321/100, 0.011/100)
    eps_acc_2018 = ufloat(5.330/100, 0.010/100)
    eps_acc = (28.5/165)*eps_acc_2016 + (60.0/165)*eps_acc_2017 + (76.5/165)*eps_acc_2018

    ### Ktautau stripping 
    eps_strip_2016 = ufloat(1.028/100, 0.003/100)
    eps_strip_2017 = ufloat(1.056/100, 0.003/100)
    eps_strip_2018 = ufloat(1.050/100, 0.002/100)
    eps_strip = (28.5/165)*eps_strip_2016 + (60.0/165)*eps_strip_2017 + (76.5/165)*eps_strip_2018

    ### Ktautua reco1 = 100%
    # N_den_2016 = 46774+19242+19252+8028
    # N_den_2017 = 61870+25333+25516+10549
    # N_den_2018 = 95913+39754+39423+16649
    # N_den = N_den_2016+N_den_2017+N_den_2018
    N_den = t_gen_all.GetEntries()
    # print(N_den)

    np.save('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_s_den.npy', N_den)

    # signal_efficiency_term = 1/(eps_acc*eps_strip)
    # print(signal_efficiency_term)

    A = (1/((B_tau_3pi_nu+B_tau_3pi_pi0_nu)**2))*(1/(eps_acc*eps_strip))
    B = (B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)*(eps_norm/N_norm)

    print("B = ", B.nominal_value, " +/- ", B.std_dev)

    np.save('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const.npy', A.nominal_value)
    np.save('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const_err.npy', A.std_dev)

    np.save('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy', B.nominal_value)
    np.save('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy', B.std_dev)


if __name__ == "__main__":
    main(sys.argv)