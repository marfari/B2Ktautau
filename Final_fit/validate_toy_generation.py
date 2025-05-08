import sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from scipy.stats import poisson

def main(argv):
    folder_name = argv[1]
    bdt1 = argv[2]
    bdt2 = argv[3]

    n_toys = 100
    nbins = 20

    f_bkg = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_bkg = f_bkg.Get("DecayTree")
    h_data = ROOT.TH1D("h_data", "h_data", 40, 4000, 8000)
    t_bkg.Draw("df_Bp_M >> h_data", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    N_toy = np.zeros(n_toys, dtype=float)
    N_bin_10 = np.zeros(n_toys, dtype=float)

    for i in range(n_toys):
        f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{0}/toy_data_bdt1_{1}_bdt2_{2}_seed_{3}.root".format(folder_name,bdt1,bdt2,i))
        h = f.Get("h_toy_data")
        N_toy[i] = h.GetEntries()
        N_bin_10[i] = h.GetBinContent(10)

    lbd = np.mean(N_toy)
    lbd_err = np.sqrt(lbd/n_toys)

    plt.hist(N_toy, bins=nbins, label="$\lambda = $ {0} $\pm$ {1}".format(int(lbd),int(lbd_err)))
    plt.xlabel('$N_{toy}$')
    plt.ylabel('Entries / {0} bins'.format(nbins))
    plt.title('BDT1 = {0} | BDT2 = {1} \n All events'.format(bdt1,bdt2))
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/validate_toy_generation/{0}/N_toy_bdt1_{1}_bdt2_{2}.pdf'.format(folder_name, bdt1, bdt2))
    plt.clf()

    lbd_10_true = h_data.GetBinContent(10)*1376/h_data.Integral(0,40-1)
    print(lbd_10_true)

    lbd_10 = np.mean(N_bin_10)
    lbd_10_err = np.sqrt(lbd_10/n_toys)

    # bin 10
    plt.hist(N_bin_10, bins=nbins, label="$\lambda = $ {0} $\pm$ {1}".format(round(lbd_10,2),round(lbd_10_err,2)))
    plt.xlabel('$N_{10}$')
    plt.ylabel('Entries / {0} bins'.format(nbins))
    plt.title('BDT1 = {0} | BDT2 = {1} \n All events'.format(bdt1,bdt2))
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/validate_toy_generation/{0}/N_10_bdt1_{1}_bdt2_{2}.pdf'.format(folder_name, bdt1, bdt2))
    plt.clf()

if __name__ == "__main__":
    main(sys.argv)