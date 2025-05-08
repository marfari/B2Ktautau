import sys
import pyhf 
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import array 
from ROOT import TMath

def optimal_nbins(h, nbins):
    # Compute optimal number of bins using the Freedman-Diaconis rule
    quantiles = np.array([0.0, 0.0])
    prob = np.array([0.25, 0.75])
    h.GetQuantiles(2, quantiles, prob)

    q1, q3 = quantiles[0], quantiles[1]
    iqr = q3-q1

    N = h.GetEntries()
    if(N == 0):
        bin_width = 4000
    else:
        bin_width = 2*iqr/(N**(1/3))

    nbins_opt = int((8000-4000)/bin_width)
    return nbins_opt


def number_of_bins(t_sig, ch, bdt1, bdt2):
    h = ROOT.TH1D("h", "h", 100, 4000, 8000)

    if((bdt1 >= 0.9) and (bdt2 >= 0.9)):
        if(ch == 0):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1})".format(0.9,0.9))
        elif(ch == 1):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(0.9,0.9))
        elif(ch == 2):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(0.9,0.9))
        elif(ch == 3):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(0.9,0.9))
        else:
            print("Wrong channel value")
            quit()
    else:
        if(ch == 0):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
        elif(ch == 1):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(bdt1,bdt2))
        elif(ch == 2):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(bdt1,bdt2))
        elif(ch == 3):
            t_sig.Draw("df_Bp_M >> h", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(bdt1,bdt2))
        else:
            print("Wrong channel value")
            quit()

    bin1 = h.FindFirstBinAbove(h.GetMaximum()/2)
    bin2 = h.FindLastBinAbove(h.GetMaximum()/2)
    fwhm = h.GetBinCenter(bin2) - h.GetBinCenter(bin1)

    sigma = fwhm/2.4
    print("ch = ", ch, " | sigma = ", sigma)

    return int(4000/(sigma/2))


def merge_bins(t_sig, bdt1, bdt2, ch, x_low=4000, x_up=8000):
    
    nbins = number_of_bins(t_sig, ch, bdt1, bdt2)

    h = ROOT.TH1D("h", "h", nbins, 4000, 8000)

    initial_bins = np.zeros(nbins+1)

    for i in range(nbins+1):
        initial_bins[i] = h.GetBinLowEdge(i+1)

    bins = []
    for i in range(nbins+1):
        if( ((initial_bins[i] < x_low) or (initial_bins[i] > x_up)) and (i != 0) and (i != nbins) ):
            continue
        else:
            bins.append( float(initial_bins[i]) )

    return array.array('d', bins)


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
    f_bkg = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_bkg = f_bkg.Get("DecayTree")

    # at each BDT cut the number of bins is chosen s.t. bin_width ~ half of the signal resolution
    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")

    # 1) Create histogram (WS data B+ mass)
    # bin width ~half of signal mass resolution
    nbins = number_of_bins(t_sig, 0, bdt1, bdt2)
    nbins_1 = number_of_bins(t_sig, 1, bdt1, bdt2)
    nbins_2 = number_of_bins(t_sig, 2, bdt1, bdt2)
    nbins_3 = number_of_bins(t_sig, 3, bdt1, bdt2)

    bins = merge_bins(t_sig, bdt1, bdt2, 0, x_low=4000, x_up=8000)
    bins_1 = merge_bins(t_sig, bdt1, bdt2, 1, x_low=4000, x_up=5750)
    bins_2 = merge_bins(t_sig, bdt1, bdt2, 2, x_low=4000, x_up=6250)
    bins_3 = merge_bins(t_sig, bdt1, bdt2, 3, x_low=4000, x_up=7250)

    nbins = len(bins)-1
    nbins_1 = len(bins_1)-1
    nbins_2 = len(bins_2)-1
    nbins_3 = len(bins_3)-1

    h_data = ROOT.TH1D("h_data", "h_data", nbins, bins)
    h_data_1 = ROOT.TH1D("h_data_1", "h_data_1", nbins_1, bins_1)
    h_data_2 = ROOT.TH1D("h_data_2", "h_data_2", nbins_2, bins_2)
    h_data_3 = ROOT.TH1D("h_data_3", "h_data_3", nbins_3, bins_3)

    # I do this because the shape of the WS B+ mass does not change much at higher BDT cuts and above 0.97 expecially the background distributions has very low statistics
    if((bdt1 >= 0.9) and (bdt2 >= 0.9)):
        t_bkg.Draw("df_Bp_M >> h_data", "(BDT1 > {0}) && (BDT2 > {1})".format(0.9,0.9))
        t_bkg.Draw("df_Bp_M >> h_data_1", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(0.9,0.9))
        t_bkg.Draw("df_Bp_M >> h_data_2", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(0.9,0.9))
        t_bkg.Draw("df_Bp_M >> h_data_3", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(0.9,0.9))
    else:
        t_bkg.Draw("df_Bp_M >> h_data", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
        t_bkg.Draw("df_Bp_M >> h_data_1", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(bdt1,bdt2))
        t_bkg.Draw("df_Bp_M >> h_data_2", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(bdt1,bdt2))
        t_bkg.Draw("df_Bp_M >> h_data_3", "(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(bdt1,bdt2))

    # 2) Compute expected number of combinatorial background events in RS data
    eps_ws_den = t_bkg.GetEntries()
    eps_ws_num = t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws_num_1 = t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(bdt1,bdt2))
    eps_ws_num_2 = t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(bdt1,bdt2))
    eps_ws_num_3 = t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(bdt1,bdt2))

    eps_ws = eps_ws_num/eps_ws_den
    eps_ws_1 = eps_ws_num_1/eps_ws_den
    eps_ws_2 = eps_ws_num_2/eps_ws_den
    eps_ws_3 = eps_ws_num_3/eps_ws_den

    N_bkg = int(2*19006409*eps_ws)
    N_bkg_1 = int(2*19006409*eps_ws_1)
    N_bkg_2 = int(2*19006409*eps_ws_2)
    N_bkg_3 = int(2*19006409*eps_ws_3)

    # For validating the fit I let the number Poisson fluctuate
    if((random_seed < 1000) and (random_seed >= 0)):
        np.random.seed(random_seed)
        N_bkg = np.random.poisson(N_bkg)
        N_bkg_1 = np.random.poisson(N_bkg_1)
        N_bkg_2 = np.random.poisson(N_bkg_2)
        N_bkg_3 = np.random.poisson(N_bkg_3)

    # 3) Sample from WS data histogram
    h_toy_data = ROOT.TH1D("h_toy_data", "h_toy_data", nbins, 4000, 8000)
    h_toy_data_1 = ROOT.TH1D("h_toy_data_1", "h_toy_data_1", nbins_1, bins_1)
    h_toy_data_2 = ROOT.TH1D("h_toy_data_2", "h_toy_data_2", nbins_2, bins_2)
    h_toy_data_3 = ROOT.TH1D("h_toy_data_3", "h_toy_data_3", nbins_3, bins_3)

    rnd = ROOT.TRandom()
    rnd.SetSeed(random_seed)

    h_toy_data.FillRandom(h_data, N_bkg, rnd)
    h_toy_data_1.FillRandom(h_data_1, N_bkg_1, rnd)
    h_toy_data_2.FillRandom(h_data_2, N_bkg_2, rnd)
    h_toy_data_3.FillRandom(h_data_3, N_bkg_3, rnd)

    print("All events: ", h_toy_data.GetEntries())
    print("Channel 1: ", h_toy_data_1.GetEntries())
    print("Channel 2: ", h_toy_data_2.GetEntries())
    print("Channel 3: ", h_toy_data_3.GetEntries())

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout.cd()
    h_toy_data.Write()
    fout.Close()

    fout_1 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}_ch1.root".format( bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout_1.cd()
    h_toy_data_1.Write()
    fout_1.Close()

    fout_2 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}_ch2.root".format( bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout_2.cd()
    h_toy_data_2.Write()
    fout_2.Close()

    fout_3 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}_ch3.root".format( bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout_3.cd()
    h_toy_data_3.Write()
    fout_3.Close()

if __name__ == "__main__":
    main(sys.argv)


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




# # 2) Fit the WS data histogram
# def landau(x, par):
#     xx = x[0]
#     mu = par[0]
#     sigma = par[1]

#     return TMath.Landau(xx, mu, sigma)

# def chebyshev_pdf(x, par):
#     # Normalize x to [-1, 1]
#     xx = x[0]
#     xmin = 4000
#     xmax = 8000
#     xnorm = 2 * (xx - xmin) / (xmax - xmin) - 1

#     # Evaluate the Chebyshev sum
#     n = len(par)
#     T_prev = 1
#     T_curr = xnorm
#     result = par[0] * T_prev
#     if n > 1:
#         result += par[1] * T_curr

#     for i in range(2, n):
#         T_next = 2 * xnorm * T_curr - T_prev
#         result += par[i] * T_next
#         T_prev, T_curr = T_curr, T_next

#     return result

# def landau_plus_chebychev(x, par):
#     mu = par[0]
#     sigma = par[1]
#     a0 = par[2]
#     a1 = par[3]
#     a2 = par[4]
#     a3 = par[5]
#     f = par[6]
#     N = par[7]

#     p = [mu, sigma]
#     p1 = [a0, a1, a2, a3]

#     return N*( f*landau(x, p) + (1-f)*chebyshev_pdf(x, p1) )

# # Landau + chebychev
# pdf = ROOT.TF1("pdf", landau_plus_chebychev, 4000, 8000, 8)
# pdf.SetParameters(4800, 400, 0.5, 0.5, 0.5, 0.05, 0.5, eps_ws_num)
# pdf.SetParNames("mu", "sig", "a0", "a1", "a2", "a3", "f", "N")
# pdf.SetParLimits(0, 4500, 5500)
# pdf.SetParLimits(1, 50, 1000)
# pdf.SetParLimits(2, -1, 10)
# pdf.SetParLimits(3, -5, 5)
# pdf.SetParLimits(3, -5, 5)
# pdf.SetParLimits(4, -2, 2)
# pdf.SetParLimits(5, -1, 1)
# pdf.FixParameter(6, 0.95)
# if((bdt1 >= 0.98) and (bdt2 >= 0.98)):
#     pdf.FixParameter(6, 1)
# pdf.SetParLimits(7, 0, 2*eps_ws_num)

# fit_status = int(h_data.Fit("pdf", "LME"))
# fit_pdf = h_data.GetFunction("pdf")

# chi2_bkg = fit_pdf.GetChisquare() / fit_pdf.GetNDF()
# print("Background fit status = ", fit_status)
# print("Background fit chi2/ndf = ", chi2_bkg)

# if(eps_ws_num < 1000):
#     h_toy_data.FillRandom("pdf", N_bkg, rnd)
# else:



# # Uniform distribution:
# contents = np.array([h_data.GetBinContent(i+1) for i in range(nbins)])
# probs = contents/np.sum(contents) # probability per bin
# cdf = np.cumsum(probs)

# for i in range(N):
#     r = np.random.uniform()
#     bin_index = int(np.searchsorted(cdf, r) + 1)
#     h_toy_data.Fill( h_data.GetBinCenter(bin_index), 1)

# h_toy_data.FillRandom(h_data, N_bkg, rnd)