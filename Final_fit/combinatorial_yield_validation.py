import sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from uncertainties import ufloat

def main(argv):

    bdt_cuts = np.round( np.linspace(0.99,1,200), 4)
    N = len(bdt_cuts)

    estimated_comb_yield = np.zeros(N)
    estimated_comb_yield_err = np.zeros(N)
    real_comb_yield = np.zeros(N)
    real_comb_yield_err = np.zeros(N)

    f_rs = ROOT.TFile('/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root')
    t_rs = f_rs.Get("DecayTree")

    i = 0
    for bdt in bdt_cuts:
        f = ROOT.TFile(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms.root')
        h_comb = f.Get(f"Channel_0/h_comb_0")
        the_bin = h_comb.FindBin(6500)
        N_ws_num = h_comb.Integral(the_bin, h_comb.GetNbinsX()+1) 
        N_ws_den = h_comb.Integral(0, h_comb.GetNbinsX()+1) 

        ws_sideband_eff = N_ws_num/N_ws_den
        ws_sideband_eff_up = ROOT.TEfficiency.Wilson(N_ws_den, N_ws_num, 0.68, True)
        ws_sideband_eff_down = ROOT.TEfficiency.Wilson(N_ws_den, N_ws_num, 0.68, False)
        ws_sideband_eff_err = 0.5*(ws_sideband_eff_up - ws_sideband_eff_down)
        ws_sideband_eff = ufloat(ws_sideband_eff, ws_sideband_eff_err)

        n_comb_estimate = ufloat( np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb.npy')[0], np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_err.npy')[0])
        n_comb_estimate *= ws_sideband_eff

        estimated_comb_yield[i] = n_comb_estimate.nominal_value
        estimated_comb_yield_err[i] = n_comb_estimate.std_dev

        n_comb_real = t_rs.GetEntries(f'(df_Bp_M > 6500) && (BDT > {bdt})')
        n_comb_real = ufloat( n_comb_real, np.sqrt(n_comb_real) )

        real_comb_yield[i] = n_comb_real.nominal_value
        real_comb_yield_err[i] = n_comb_real.std_dev

        print("BDT = ", bdt, " | N_comb_estimate = ", estimated_comb_yield[i], " | N_comb_real = ", real_comb_yield[i])
        i += 1

    plt.plot(bdt_cuts, estimated_comb_yield, marker='.', c='red', label="Estimated", linestyle='none')
    plt.plot(bdt_cuts, real_comb_yield, marker='.', c='blue', label="Real", linestyle='none')
    plt.xlabel('BDT cut')
    plt.ylabel('Combinatorial yield')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("/panfs/felician/B2Ktautau/workflow/combinatorial_yield_validation/ncomb_real_vs_expected.pdf")


if __name__ == "__main__":
    main(sys.argv)