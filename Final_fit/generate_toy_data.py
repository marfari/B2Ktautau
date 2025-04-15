import sys
import pyhf 
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import array 

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]
    folder_name = argv[4]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    if((folder_name != "background_only") and (folder_name != "signal_injection")):
        print("wrong folder name, try 'background_only' or 'signal_injection'")
        quit()

    # Toy data: combinatorial background only for now (WS data)
    # Expected combinatorial background yield in data: 2*N_RS_prebdt*eps_ws (eps_rs=2*eps_ws for BDT1=BDT2=0.8)
    f_bkg = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_bkg = f_bkg.Get("DecayTree")

    # Different binning schemes
    # f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    # t_sig = f_sig.Get("DecayTree")

    # if( (t_sig.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2)) < 500) or (t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2)) < 500) ):
    #     nbins = 10
    # if( (t_sig.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2)) < 5000) or (t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2)) < 5000) ):
    #     nbins = 10

    # variable bin-width
    # bins = array.array('d', [4000, 4500, 4580, 4660, 4740, 4820, 4900, 4980, 5060, 5140, 5220, 5300, 5380, 5460, 5540, 5620, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6700, 6900, 7200, 7500, 8000 ])
    # bins = array.array('d', [4000, 4100, 4200, 4300, 4400, 4500, 4540, 4580, 4620, 4660, 4700, 4740, 4780, 4820, 4860, 4900, 4940, 4980, 5020, 5060, 5100, 5140, 5180, 5220, 5260, 5300, 5340, 5380, 5420, 5460, 5500, 5540, 5580, 5620, 5660, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000])
    # bins = array.array('d', [4000, 4250, 4500, 4540, 4580, 4620, 4660, 4700, 4740, 4780, 4820, 4860, 4900, 4940, 4980, 5020, 5060, 5100, 5140, 5180, 5220, 5260, 5300, 5340, 5380, 5420, 5460, 5500, 5540, 5580, 5620, 5660, 5700, 5740, 5780, 5820, 5860, 5900, 5940, 5980, 6020, 6240, 6460, 6680, 6900, 7120, 7340, 7560, 7780, 8000])
    # if((bdt1 > 0.96) or (bdt2 > 0.96)): # 80 MeV
    #     bins = array.array('d', [4000, 4250, 4500, 4580, 4660, 4740, 4820, 4900, 4980, 5060, 5140, 5220, 5300, 5380, 5460, 5540, 5620, 5700, 5780, 5860, 5940, 6020, 6100, 6180, 6260, 6340, 6420, 6500, 6580, 6660, 6740, 6820, 6900, 6980, 7060, 7295, 7530, 7765, 8000])
    # if((bdt1 > 0.98) or (bdt2 > 0.98)): # 200 MeV
    #     bins =  array.array('d', [4000, 4250, 4500, 4700, 4900, 5100, 5300, 5500, 5700, 5900, 6100, 6300, 6500, 6700, 6900, 7450, 8000 ])
    # if((bdt1 > 0.99) or (bdt2 > 0.99)): # 4000
    #     bins =  array.array('d', [4000, 4400, 4800, 5200, 5600, 6000, 6400, 7000, 8000])

    # bins = array.array('d', [4000, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 7250, 8000])
    
    # bins = array.array('d', [4000, 4500, 4650, 4800, 4950, 5100, 5250, 5400, 5550, 5700, 5850, 6000, 6500, 8000])
    # nbins = len(bins) - 1

    # Compute optimal number of bins using the Freedman-Diaconis rule
    # nbins = 40
    # if((bdt1 >= 0.5) and (bdt2 >= 0.5)):
    #     h_data_norm = ROOT.TH1D("h_data_norm", "h_data_norm", nbins, 4000, 8000)
    #     t_bkg.Draw("df_Bp_M >> h_data_norm", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    #     quantiles = np.array([0.0, 0.0])
    #     prob = np.array([0.25, 0.75])
    #     h_data_norm.GetQuantiles(2, quantiles, prob)

    #     q1, q3 = quantiles[0], quantiles[1]
    #     iqr = q3-q1

    #     N = h_data_norm.GetEntries()
    #     if(N == 0):
    #         bin_width = 4000
    #     else:
    #         bin_width = 2*iqr/(N**(1/3))

    #     nbins = int((8000-4000)/bin_width)
    #     print(bin_width)

    # 1) Create histogram (WS data B+ mass)
    nbins = 40 # 100 MeV bin width (~half of signal mass resolution)

    # Use optimal number of bins
    h_data = ROOT.TH1D("h_data", "h_data", nbins, 4000, 8000)
    t_bkg.Draw("df_Bp_M >> h_data", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    # 2) Sample from WS data histogram
    eps_ws_den = t_bkg.GetEntries()
    eps_ws_num = t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws = eps_ws_num/eps_ws_den
    N_bkg = int(2*19006409*eps_ws)
    
    np.random.seed(random_seed)
    N = np.random.poisson(N_bkg)

    h_toy_data = ROOT.TH1D("h_toy_data", "h_toy_data", nbins, 4000, 8000)

    rnd = ROOT.TRandom()
    rnd.SetSeed(random_seed)

    if(random_seed == 1000):
        print("Using B")
        h_toy_data.FillRandom(h_data, N_bkg, rnd)
    elif((random_seed < 1000) and (random_seed >= 0)):
        h_toy_data.FillRandom(h_data, N, rnd)
    else:
        print("random seed must be [0,1000]")
        quit()

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout.cd()
    h_toy_data.Write()
    fout.Close()


if __name__ == "__main__":
    main(sys.argv)
