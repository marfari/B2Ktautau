import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import ROOT

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    folder_name = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    nbins=50
    N = 1000

    fit_poi = np.zeros(N, dtype=float)
    fit_poi_error = np.zeros(N, dtype=float)

    nbkg = np.zeros(N, dtype=float)
    nbkg_error = np.zeros(N, dtype=float)

    for i in range(N):
        fit_result = np.load('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, i ))
        fit_poi[i] = fit_result[0]
        fit_poi_error[i] = fit_result[1]

        nbkg[i] = fit_result[2]
        nbkg_error[i] = fit_result[3]

    bias = fit_poi - 0
    (mu_bias, sigma_bias) = norm.fit(bias)

    n, bins, patches = plt.hist(bias, bins=nbins, density=True, label="$\mu = $ {0} \n $\sigma = $ {1}".format( round(mu_bias,5), round(sigma_bias,5) ))

    y_bias = norm.pdf( bins, mu_bias, sigma_bias)
    plt.plot(bins, y_bias, 'r-', linewidth=2)
    plt.xlabel("Bias: $BF_{sig}$ - 0")
    plt.ylabel("Normalised entries / {0} bins".format(nbins))
    plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g}")
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/bias_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))
    plt.clf()

    pull = bias/fit_poi_error
    (mu_pull, sigma_pull) = norm.fit(pull)

    n1, bins1, patches1 = plt.hist( pull, bins=nbins, density=True, label="$\mu = $ {0} \n $\sigma = $ {1}".format( round(mu_pull,5), round(sigma_pull,5) ))

    y_pull = norm.pdf( bins1, mu_pull, sigma_pull)
    plt.plot(bins1, y_pull, 'r-', linewidth=2)
    plt.xlabel("Pull: ($BF_{sig}$ - 0)/error")
    plt.ylabel("Normalised entries / {0} bins".format(nbins))
    plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g}")
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/pull_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))


    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_ws = f_ws.Get("DecayTree")
    eps_ws_den = t_ws.GetEntries()
    eps_ws_num = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws = eps_ws_num/eps_ws_den
    N_bkg = 2*19006409*eps_ws

    bias1 = nbkg - N_bkg
    (mu_bias1, sigma_bias1) = norm.fit(bias1)

    n2, bins2, patches2 = plt.hist(bias1, bins=nbins, density=True, label="$\mu = $ {0} \n $\sigma = $ {1}".format( round(mu_bias1,5), round(sigma_bias1,5) ))

    y_bias1 = norm.pdf( bins2, mu_bias1, sigma_bias1)
    plt.plot(bins2, y_bias1, 'r-', linewidth=2)
    plt.xlabel("Bias: $N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$")
    plt.ylabel("Normalised entries / {0} bins".format(nbins))
    plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g}")
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/bias_nbkg_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))
    plt.clf()


    pull1 = bias1/nbkg_error
    (mu_pull1, sigma_pull1) = norm.fit(pull1)

    n3, bins3, patches3 = plt.hist( pull1, bins=nbins, density=True, label="$\mu = $ {0} \n $\sigma = $ {1}".format( round(mu_pull1,5), round(sigma_pull1,5) ))

    y_pull1 = norm.pdf( bins3, mu_pull, sigma_pull)
    plt.plot(bins3, y_pull1, 'r-', linewidth=2)
    plt.xlabel("Pull: ($N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$)/error")
    plt.ylabel("Normalised entries / {0} bins".format(nbins))
    plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g}")
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/pull_nbkg_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))




if __name__ == "__main__":
    main(sys.argv)