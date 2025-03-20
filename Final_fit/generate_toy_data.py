import sys
import pyhf 
import ROOT
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    # Toy data: combinatorial background only for now (WS data)
    # Expected combinatorial background yield in data: 2*N_RS_prebdt*eps_ws (eps_rs=2*eps_ws for BDT1=BDT2=0.8)
    f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t = f.Get("DecayTree")

    h_data = ROOT.TH1D("h_data", "h_data", 100, 4000, 8000)
    t.Draw("df_Bp_M >> h_data", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    eps_ws_den = t.GetEntries()
    eps_ws_num = t.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws = eps_ws_num/eps_ws_den

    N_bkg = 2*19006409*eps_ws
    N_bkg_error = np.sqrt(N_bkg)

    np.random.seed(random_seed)
    N = np.random.poisson(N_bkg)

    h_toy_data = ROOT.TH1D("h_toy_data", "h_toy_data", 100, 4000, 8000)

    rnd = ROOT.TRandom()
    rnd.SetSeed(random_seed)
    h_toy_data.FillRandom(h_data, N, rnd)

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2 ,random_seed), "RECREATE")
    fout.cd()
    h_toy_data.Write()
    fout.Close()


if __name__ == "__main__":
    main(sys.argv)
