import sys
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    fit_type = argv[1]
    channel_type = argv[2]
    bdt = argv[3]

    branching_fractions = [0.00001, 0.00002, 0.00004, 0.00006, 0.00008, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01]
    N = len(branching_fractions)
    # obs_upper_limit = np.zeros(N) 
    BF_toys_mean = np.zeros(N) 
    BF_toys_error = np.zeros(N)

    for i in range(N):
        if(branching_fractions[i] == 0.00001):
            BF_toys_mean[i] = abs(np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_1e-05/BDT_{bdt}/fit_poi_mean.npy"))
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_1e-05/BDT_{bdt}/fit_poi_mean_error.npy")
        elif(branching_fractions[i] == 0.00002):
            BF_toys_mean[i] = abs(np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_2e-05/BDT_{bdt}/fit_poi_mean.npy"))
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_2e-05/BDT_{bdt}/fit_poi_mean_error.npy")
        elif(branching_fractions[i] == 0.00004):
            BF_toys_mean[i] = abs(np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_4e-05/BDT_{bdt}/fit_poi_mean.npy"))
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_4e-05/BDT_{bdt}/fit_poi_mean_error.npy")
        elif(branching_fractions[i] == 0.00006):
            BF_toys_mean[i] = abs(np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_6e-05/BDT_{bdt}/fit_poi_mean.npy"))
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_6e-05/BDT_{bdt}/fit_poi_mean_error.npy")
        elif(branching_fractions[i] == 0.00008):
            BF_toys_mean[i] = abs(np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_8e-05/BDT_{bdt}/fit_poi_mean.npy"))
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_8e-05/BDT_{bdt}/fit_poi_mean_error.npy")    
        else:
            BF_toys_mean[i] = abs(np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{branching_fractions[i]}/BDT_{bdt}/fit_poi_mean.npy"))
            BF_toys_error[i] = np.load(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{branching_fractions[i]}/BDT_{bdt}/fit_poi_mean_error.npy")

    # plt.errorbar(branching_fractions, obs_upper_limit, marker='.', label="Observed limit")
    plt.errorbar(branching_fractions, BF_toys_mean, yerr=BF_toys_error, marker='.')
    plt.xlabel('Injected BF')
    plt.ylabel('Absolute value of fitted BF mean')
    plt.grid()
    plt.title(f'{fit_type} | {channel_type} | BDT > {bdt}')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(f"/panfs/felician/B2Ktautau/workflow/signal_toys_bf/{fit_type}_{channel_type}/BDT_{bdt}/BF_injected_vs_expected.pdf")

if __name__ == "__main__":
    main(sys.argv)