import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    folder_name = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    nbins=50

    fit_poi = np.zeros(1000, dtype=float)
    fit_poi_error = np.zeros(1000, dtype=float)

    for i in range(1000):
        fit_result = np.load('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, i ))
        fit_poi[i] = fit_result[0]
        fit_poi_error[i] = fit_result[1]

    bias = fit_poi - 0
    (mu_bias, sigma_bias) = norm.fit(bias)

    n, bins, patches = plt.hist(bias, bins=nbins, density=True, label="$\mu = $ {0} \n $\sigma = $ {1}".format( round(mu_bias,5), round(sigma_bias,5) ))

    y_bias = norm.pdf( bins, mu_bias, sigma_bias)
    plt.plot(bins, y_bias, 'r-', linewidth=2)
    plt.xlabel("Bias: $BF_{sig}$ - 0")
    plt.ylabel("Normalised entries / {0} bins".format(nbins))
    plt.title(f"BDT1 = {bdt1:.3g} | BDT2 = {bdt2:.3g}")
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
    plt.title(f"BDT1 = {bdt1:.3g} | BDT2 = {bdt2:.3g}")
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/pull_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))



if __name__ == "__main__":
    main(sys.argv)