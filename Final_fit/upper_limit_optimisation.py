import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import ROOT
from uncertainties import ufloat

def main(argv):
    fit_type = argv[1]
    channel_type = argv[2]
    BF_sig = argv[3]

    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root")
    t_sig = f_sig.Get("DecayTree")
    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt_0.0.root")
    t_ws = f_ws.Get("DecayTree")
    t_BDDKp = ROOT.TChain("DecayTree")
    t_BDDKp.Add("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root")
    t_BDDKp.Add("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root")
    t_BDDKp.Add("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root")

    bdt_cuts = np.round( np.linspace(0.99,1,100), 4)
    N = len(bdt_cuts)
    branching_fraction_values = np.zeros(N)
    branching_fraction_values_up = np.zeros(N)
    branching_fraction_values_down = np.zeros(N)
    eps_sig = np.zeros(N)
    eps_sig_err = np.zeros(N)
    eps_ws = np.zeros(N)
    eps_ws_err = np.zeros(N)
    eps_BDDKp = np.zeros(N)
    eps_BDDKp_err = np.zeros(N)
    n_sig_mc = t_sig.GetEntries()
    n_ws = t_ws.GetEntries()
    n_BDDKp_mc = t_BDDKp.GetEntries()

    B = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B_err.npy')
    b = ufloat(B, B_err)

    N_sig = np.zeros(N)
    N_sig_err = np.zeros(N)
    N_comb = np.zeros(N)
    N_comb_err = np.zeros(N)
    N_BDDKp = np.zeros(N)
    N_BDDKp_err = np.zeros(N)

    for i in range(N):
        bdt = bdt_cuts[i]

        eps_sig_num = t_sig.GetEntries(f"BDT > {bdt}")
        eps_sig_up = ROOT.TEfficiency.Wilson(n_sig_mc, eps_sig_num, 0.68, True)
        eps_sig_down = ROOT.TEfficiency.Wilson(n_sig_mc, eps_sig_num, 0.68, False)
        eps_sig[i] = (eps_sig_num/n_sig_mc)*100
        eps_sig_err[i] = 0.5*(eps_sig_up - eps_sig_down)*100

        eps_ws_num = t_ws.GetEntries(f"BDT > {bdt}")
        eps_ws_up = ROOT.TEfficiency.Wilson(n_ws, eps_ws_num, 0.68, True)
        eps_ws_down = ROOT.TEfficiency.Wilson(n_ws, eps_ws_num, 0.68, False)
        eps_ws[i] = (eps_ws_num/n_ws)*100
        eps_ws_err[i] = 0.5*(eps_ws_up - eps_ws_down)*100

        eps_BDDKp_num = t_BDDKp.GetEntries(f"BDT > {bdt}")
        eps_BDDKp_up = ROOT.TEfficiency.Wilson(n_BDDKp_mc, eps_BDDKp_num, 0.68, True)
        eps_BDDKp_down = ROOT.TEfficiency.Wilson(n_BDDKp_mc, eps_BDDKp_num, 0.68, False)
        eps_BDDKp[i] = (eps_BDDKp_num/n_BDDKp_mc)*100
        eps_BDDKp_err[i] = 0.5*(eps_BDDKp_up - eps_BDDKp_down)*100

        A = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A.npy')[0]
        A_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_err.npy')[0]
        C = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C.npy')[0]
        C_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_err.npy')[0]

        a = ufloat(A, A_err)
        c = ufloat(C, C_err)

        try:
            branching_fraction_values[i] = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig}/BDT_{bdt}/CLs_limit/cls_limit.npy')
        except:
            branching_fraction_values[i] = np.inf

        try:
            branching_fraction_values_up[i] = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig}/BDT_{bdt}/CLs_limit/cls_limit_up.npy')
        except:
            branching_fraction_values_up[i] = branching_fraction_values[i]

        try:
            branching_fraction_values_down[i] = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig}/BDT_{bdt}/CLs_limit/cls_limit_down.npy')
        except:
            branching_fraction_values_down[i] = branching_fraction_values[i]

        BF = ufloat(branching_fraction_values[i], max(np.abs( branching_fraction_values[i] - branching_fraction_values_up[i] ), np.abs( branching_fraction_values[i] - branching_fraction_values_down[i] )))

        if(A != 0):
            N_sig[i] = (BF/(a*b)).nominal_value
            N_sig_err[i] = (BF/(a*b)).std_dev
        N_comb[i] = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb.npy')[0]
        N_comb_err[i] = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_err.npy')[0]
        N_BDDKp[i] = (c*b).nominal_value
        N_BDDKp_err[i] = (c*b).std_dev

        print(f"BDT = {bdt} | BF = {branching_fraction_values[i]} | eps_sig = {eps_sig[i]} | eps_comb = {eps_ws[i]} | eps_BDDKp = {eps_BDDKp[i]}")     
    
    min_index = np.argmin(branching_fraction_values)

    print(f"The BDT cut, BDT > {bdt_cuts[min_index]}, minimise the 90% C.L. upper limit to be \n B(B -> K tau tau ) < {branching_fraction_values[min_index]} ")
    print(f"At the minimum: N_sig = {N_sig[min_index]} +/- {N_sig_err[min_index]} | N_comb = {N_comb[min_index]} +/- {N_comb_err[min_index]} | N_BDDKp = {N_BDDKp[min_index]} +/- {N_BDDKp_err[min_index]}")
    print(f"At the minimum: eps_sig = {eps_sig[min_index]} +/- {eps_sig_err[min_index]} | eps_ws = {eps_ws[min_index]} +/- {eps_ws_err[min_index]} | eps_BDDKp = {eps_BDDKp[min_index]} +/- {eps_BDDKp_err[min_index]}")

    up_err_comb_yield_syst = np.abs( branching_fraction_values - branching_fraction_values_up )
    down_err_comb_yield_syst = np.abs( branching_fraction_values - branching_fraction_values_down )

    # Branching fraction
    plt.errorbar(bdt_cuts, branching_fraction_values, yerr=[up_err_comb_yield_syst, down_err_comb_yield_syst], fmt='.', linestyle='none')
    plt.plot(bdt_cuts[min_index], branching_fraction_values[min_index], marker="*", color="black", markersize=15)  
    plt.xlabel('BDT > x')
    plt.ylabel('90% C.L. upper limit on B($B^+ \\to K^+ \\tau^+ \\tau^-$)')
    if(channel_type == "AllEvents"):
        plt.title("Fit to all events")
    else:
        plt.title("Fit in error categories")
    plt.tight_layout()
    plt.grid()
    plt.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{channel_type}/BF_sig_{BF_sig}/BF_vs_bdt.pdf")
    plt.clf()

    # Signal and background efficiencies
    plt.errorbar(bdt_cuts, eps_sig, yerr=eps_sig_err, fmt='.', linestyle='none', c="blue", label="Signal MC")
    plt.plot(bdt_cuts[min_index], eps_sig[min_index], marker="*", color="blue", markersize=15)  
    plt.errorbar(bdt_cuts, eps_ws, yerr=eps_ws_err, fmt='.', linestyle='none', c="red", label="WS data")
    plt.plot(bdt_cuts[min_index], eps_ws[min_index], marker="*", color="red", markersize=15)  
    plt.errorbar(bdt_cuts, eps_BDDKp, yerr=eps_BDDKp_err, fmt='.', linestyle='none', c="green", label="$BDDK^{+}$ MC")
    plt.plot(bdt_cuts[min_index], eps_BDDKp[min_index], marker="*", color="green", markersize=15)  
    plt.xlabel('BDT > x')
    plt.ylabel('BDT efficiency (%)')
    if(channel_type == "AllEvents"):
        plt.title("Fit to all events")
    else:
        plt.title("Fit in error categories")
    plt.legend()
    plt.tight_layout()
    plt.grid()
    plt.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{channel_type}/BF_sig_{BF_sig}/sig_vs_bkg_eff.pdf")
    plt.clf()

    # Signal and background yields
    plt.errorbar(bdt_cuts, N_sig, yerr=N_sig_err, fmt='.', linestyle='none', c="blue", label="Signal MC")
    plt.plot(bdt_cuts[min_index], N_sig[min_index], marker="*", color="blue", markersize=15)  
    plt.errorbar(bdt_cuts, N_comb, yerr=N_comb_err, fmt='.', linestyle='none', c="red", label="WS data")
    plt.plot(bdt_cuts[min_index], N_comb[min_index], marker="*", color="red", markersize=15)  
    plt.errorbar(bdt_cuts, N_BDDKp, yerr=N_BDDKp_err, fmt='.', linestyle='none', c="green", label="$BDDK^{+}$ MC")
    plt.plot(bdt_cuts[min_index], N_BDDKp[min_index], marker="*", color="green", markersize=15)  
    plt.xlabel('BDT > x')
    plt.ylabel('Yield')
    if(channel_type == "AllEvents"):
        plt.title("Fit to all events")
    else:
        plt.title("Fit in error categories")
    plt.legend()
    plt.tight_layout()
    plt.grid()
    plt.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{channel_type}/BF_sig_{BF_sig}/sig_vs_bkg_yield.pdf")
    plt.clf()


if __name__ == "__main__":
    main(sys.argv)