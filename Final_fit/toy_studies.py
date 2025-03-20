import sys
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    fit_poi = np.zeros(100, dtype=float)
    fit_poi_error = np.zeros(100, dtype=float)

    for i in range(100):
        fit_result = np.load('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/fit_result_bdt1_{0}_bdt2_{1}_seed_{2}.npy'.format( bdt1, bdt2, i ))
        fit_poi[i] = fit_result[0]
        fit_poi_error[i] = fit_result[1]

    pull = (fit_poi - 0)/fit_poi_error

    plt.hist( pull , label="$\mu = $ {0}".format( round(np.mean(pull),5) ))
    plt.xlabel("Pull: (Nsig - mean)/error")
    plt.ylabel("Entries")
    plt.title(f"BDT1 = {bdt1:.1g} | BDT2 = {bdt2:.1g}")
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/toy_studies/toy_nsig_pulls_bdt1_{0}_bdt2_{1}.pdf'.format( bdt1, bdt2 ))



if __name__ == "__main__":
    main(sys.argv)