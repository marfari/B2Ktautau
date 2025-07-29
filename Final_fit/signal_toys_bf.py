import sys
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    fit_type = argv[1]
    fit_name = argv[2]
    bdt = argv[3]

    branching_fractions = [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01]
    N = len(branching_fractions)
    # obs_upper_limit = np.zeros(N) 
    BF_toys_mean = np.zeros(N) 
    BF_toys_error = np.zeros(N)

    for i in range(N):
        if(branching_fractions[i] == 0.0000001):
            # obs_upper_limit[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-07/BDT_{bdt}/cls_obs_limit.npy")
            BF_toys_mean[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-07/BDT_{bdt}/fit_poi_mean.npy")
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-07/BDT_{bdt}/fit_poi_mean_error.npy")
        elif(branching_fractions[i] == 0.000001):
            # obs_upper_limit[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-06/BDT_{bdt}/cls_obs_limit.npy")
            BF_toys_mean[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-06/BDT_{bdt}/fit_poi_mean.npy")
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-06/BDT_{bdt}/fit_poi_mean_error.npy")
        elif(branching_fractions[i] == 0.00001):
            # obs_upper_limit[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-05/BDT_{bdt}/cls_obs_limit.npy")
            BF_toys_mean[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-05/BDT_{bdt}/fit_poi_mean.npy")
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_1e-05/BDT_{bdt}/fit_poi_mean_error.npy")
        else:
            # obs_upper_limit[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{branching_fractions[i]}/BDT_{bdt}/cls_obs_limit.npy")
            BF_toys_mean[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{branching_fractions[i]}/BDT_{bdt}/fit_poi_mean.npy")
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{branching_fractions[i]}/BDT_{bdt}/fit_poi_mean_error.npy")

    # plt.errorbar(branching_fractions, obs_upper_limit, marker='.', label="Observed limit")
    plt.errorbar(branching_fractions, BF_toys_mean, yerr=BF_toys_error, marker='.', label="Fitted mean from toys")
    plt.xlabel('Injected BF')
    plt.ylabel('Fitted BF')
    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.title(f'{fit_type} | {fit_name} | BDT > {bdt}')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(f"/panfs/felician/B2Ktautau/workflow/signal_toys_bf/{fit_type}_{fit_name}/BDT_{bdt}/BF_injected_vs_expected.pdf")

if __name__ == "__main__":
    main(sys.argv)