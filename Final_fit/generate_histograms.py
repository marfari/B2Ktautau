import sys
import ROOT
import numpy as np
from uncertainties import ufloat
import time

channel_cut = {0: "", 1: " && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)", 2: " && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)", 3: " && (df_Bp_MERR > 250)"}
error_threshold = 0
bdt_fix_sig = 0.985
bdt_fix = 0.9

def number_of_bins(t_sig, ch, bdt1, bdt2, channel_cut):
    h = ROOT.TH1D(f"h_{ch}", f"h_{ch}", 100, 4000, 8000)

    cut_string = channel_cut[ch]

    if((bdt1 >= bdt_fix_sig) and (bdt2 >= bdt_fix_sig)):
        t_sig.Draw(f"df_Bp_M >> h_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix_sig,bdt_fix_sig)+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    bin1 = h.FindFirstBinAbove(h.GetMaximum()/2)
    bin2 = h.FindLastBinAbove(h.GetMaximum()/2)
    fwhm = h.GetBinCenter(bin2) - h.GetBinCenter(bin1)

    sigma = fwhm/2.4

    if(ch == 1):
        sigma = min(100, sigma)
    if(ch == 2):
        if(sigma < 100):
            sigma = 100
        elif(sigma > 250):
            sigma = 250
    if(ch == 3):
        sigma = max(250, sigma)

    # h1 = ROOT.TH1D(f"h1_{ch}", f"h1_{ch}", 100, 0, 1000)
    # if((bdt1 >= bdt_fix_sig) and (bdt2 >= bdt_fix_sig)):
    #     t_sig.Draw(f"df_Bp_MERR >> h1_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix_sig,bdt_fix_sig)+cut_string)
    # else:
    #     t_sig.Draw(f"df_Bp_MERR >> h1_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    # mass_resolution = h1.GetMean()
    # print("ch = ", ch, " | sigma (MERR) = ", mass_resolution, " | sigma (FWHM) = ", sigma)

    # print("resolution from FWHM = ", fwhm/2.4)
    # print("resolution from MERR = ", mass_resolution)

    nbins = int(4000/(sigma/2))
    print(nbins)

    return nbins

# def merge_bins(t_sig, bdt1, bdt2, ch, x_low=4000, x_up=8000):
     
#     nbins = number_of_bins(t_sig, ch, bdt1, bdt2)

#     h = ROOT.TH1D("h", "h", nbins, 4000, 8000)

#     initial_bins = np.zeros(nbins+1)

#     for i in range(nbins+1):
#         initial_bins[i] = h.GetBinLowEdge(i+1)

#     bins = []
#     for i in range(nbins+1):
#         if( ((initial_bins[i] < x_low) or (initial_bins[i] > x_up)) and (i != 0) and (i != nbins) ):
#             continue
#         else:
#             bins.append( float(initial_bins[i]) )

#     return array.array('d', bins)

def create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, ch, channel_cut):
    # bin width ~half of signal mass resolution
    nbins = number_of_bins(t_sig, ch, bdt1, bdt2, channel_cut)

    # Save histograms after BDT cuts
    h_sig = ROOT.TH1D(f"h_sig_{ch}", f"h_sig_{ch}", nbins, 4000, 8000)
    h_comb = ROOT.TH1D(f"h_comb_{ch}", f"h_comb_{ch}", nbins, 4000, 8000)
    ## B -> D D K+
    h_BDDKp = ROOT.TH1D(f"h_BDDKp_{ch}", f"h_BDDKp_{ch}", nbins, 4000, 8000)
    h_BuD0D0Kp = ROOT.TH1D(f"h_BuD0D0Kp_{ch}", f"h_BuD0D0Kp_{ch}", nbins, 4000, 8000)
    h_BuDpDmKp = ROOT.TH1D(f"h_BuDpDmKp_{ch}", f"h_BuDpDmKp_{ch}", nbins, 4000, 8000)
    h_BuDsDsKp = ROOT.TH1D(f"h_BuDsDsKp_{ch}", f"h_BuDsDsKp_{ch}", nbins, 4000, 8000)
    h_BdD0DmKp = ROOT.TH1D(f"h_BdD0DmKp_{ch}", f"h_BdD0DmKp_{ch}", nbins, 4000, 8000)
    h_BsD0DsKp = ROOT.TH1D(f"h_BsD0DsKp_{ch}", f"h_BsD0DsKp_{ch}", nbins, 4000, 8000)

    cut_string = channel_cut[ch]

    if((bdt1 >= bdt_fix_sig) and (bdt2 >= bdt_fix_sig)):
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix_sig,bdt_fix_sig)+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    if((bdt1 >= bdt_fix) and (bdt2 >= bdt_fix)):
        t_comb.Draw(f"df_Bp_M >> h_comb_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuD0D0Kp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 100) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDpDmKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 101) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDsDsKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 102) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdD0DmKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsD0DsKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
    else:
        t_comb.Draw(f"df_Bp_M >> h_comb_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuD0D0Kp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 100) ".format(bdt1,bdt2)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDpDmKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 101) ".format(bdt1,bdt2)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDsDsKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 102) ".format(bdt1,bdt2)+cut_string)
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdD0DmKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsD0DsKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    h_BDDKp.Sumw2()

    # Add shapes
    h_BDDKp.Add(h_BuD0D0Kp)
    h_BDDKp.Add(h_BuDpDmKp)
    h_BDDKp.Add(h_BuDsDsKp)
    h_BDDKp.Add(h_BdD0DmKp)
    h_BDDKp.Add(h_BsD0DsKp)

    if(h_sig.Integral() != 0):
        h_sig.Scale(1/h_sig.Integral())
    if(h_comb.Integral() != 0):
        h_comb.Scale(1/h_comb.Integral())
    if(h_BDDKp.Integral() != 0):
        h_BDDKp.Scale(1/h_BDDKp.Integral())

    return [h_sig, h_comb, h_BDDKp]


def create_histogram_errors(h_sig, h_comb, h_BDDKp, ch):
    nbins = h_sig.GetNbinsX()

    h_sig_err = ROOT.TH1D(f"h_sig_err_{ch}", f"h_sig_err_{ch}", nbins, 4000, 8000)
    h_comb_err = ROOT.TH1D(f"h_comb_err_{ch}", f"h_comb_err_{ch}", nbins, 4000, 8000)
    h_BDDKp_err = ROOT.TH1D(f"h_BDDKp_err_{ch}", f"h_BDDKp_err_{ch}", nbins, 4000, 8000)

    # All events
    for i in range(nbins):
        # Signal
        if(h_sig.GetBinContent(i+1) == 0):
            h_sig_err.Fill(h_sig.GetBinCenter(i+1), 0.0) 
        elif(h_sig.GetBinError(i+1)/h_sig.GetBinContent(i+1) < error_threshold):
            h_sig_err.Fill(h_sig.GetBinCenter(i+1), 0.0) 
        else:
            h_sig_err.Fill(h_sig.GetBinCenter(i+1), h_sig.GetBinError(i+1)) 
        h_sig_err.SetBinError(i+1, 0)

        # Combinatorial background
        if(h_comb.GetBinContent(i+1) == 0):
            h_comb_err.Fill(h_comb.GetBinCenter(i+1), 0.0) 
        elif(h_comb.GetBinError(i+1)/h_comb.GetBinContent(i+1) < error_threshold):
            h_comb_err.Fill(h_comb.GetBinCenter(i+1), 0.0) 
        else:
            h_comb_err.Fill(h_comb.GetBinCenter(i+1), h_comb.GetBinError(i+1)) 
        h_comb_err.SetBinError(i+1, 0) 

        # B -> D D K+
        if(h_BDDKp.GetBinContent(i+1) == 0):
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), 0.0) 
        elif(h_BDDKp.GetBinError(i+1)/h_BDDKp.GetBinContent(i+1) < error_threshold):
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), 0.0) 
        else:
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), h_BDDKp.GetBinError(i+1)) 
        h_BDDKp_err.SetBinError(i+1, 0) 

        # print("SIGNAL : ", h_sig.GetBinError(i+1), " | BACKGROUND : ", h_comb.GetBinError(i+1))
    
    return [h_sig_err, h_comb_err, h_BDDKp_err]


def combinatorial_background_yield(t_rs_data, t_comb, bdt1, bdt2):
    # 2) Estimate and save all needed yields after the BDT cuts
    # Combinatorial background yield
    eps_ws_den = t_comb.GetEntries()
    eps_rs_den = t_rs_data.GetEntries()

    eps_num_bdt = t_comb.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2))

    if((bdt1 <= 0.8) and (bdt2 <= 0.8)):
        eps_ws_num = t_comb.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2))
        eps_rs_num = t_rs_data.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2))
    else:
        eps_ws_num = t_comb.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(0.8,0.8))
        eps_rs_num = t_rs_data.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(0.8,0.8))

    eps_ws = eps_ws_num/eps_ws_den
    eps_rs = eps_rs_num/eps_rs_den

    r = eps_rs/eps_ws

    eps_ws_bdt = eps_num_bdt/eps_ws_den

    N_comb = 19006409*eps_ws_bdt*r

    return N_comb


def error_category_efficiency(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, channel_cut):
    n = len(channel_cut)-2

    eps_sig = np.zeros(n)
    eps_comb = np.zeros(n)
    eps_BDDKp = np.zeros(n)
    eps_BuDDK0 = np.zeros(n)
    eps_BuDD = np.zeros(n)

    eps_sig_err = np.zeros(n)
    eps_comb_err = np.zeros(n)
    eps_BDDKp_err = np.zeros(n)
    eps_BuDDK0_err = np.zeros(n)
    eps_BuDD_err = np.zeros(n)

    for ch in range(n):
        cut_string = f"(BDT1 > {bdt1}) && (BDT2 > {bdt2})"

        # signal
        N_num_sig = t_sig.GetEntries(cut_string+channel_cut[ch+1])
        N_den_sig = t_sig.GetEntries(cut_string)
        up_sig = ROOT.TEfficiency.Wilson(N_den_sig, N_num_sig, 0.68, True)
        lo_sig = ROOT.TEfficiency.Wilson(N_den_sig, N_num_sig, 0.68, False)
        if(N_den_sig != 0):
            eps_sig[ch] = N_num_sig/N_den_sig
        else:
            eps_sig[ch] = 0
        eps_sig_err[ch] = 0.5*(up_sig - lo_sig)

        # combinatorial
        N_num_comb = t_comb.GetEntries(cut_string+channel_cut[ch+1])
        N_den_comb = t_comb.GetEntries(cut_string)
        up_comb = ROOT.TEfficiency.Wilson(N_den_comb, N_num_comb, 0.68, True)
        lo_comb = ROOT.TEfficiency.Wilson(N_den_comb, N_num_comb, 0.68, False)
        if(N_den_comb != 0):
            eps_comb[ch] = N_num_comb/N_den_comb
        else:
            eps_comb[ch] = 0
        eps_comb_err[ch] = 0.5*(up_comb - lo_comb)

        # B -> DD K+
        N_num_BDDKp = t_BuDDKp.GetEntries(cut_string+channel_cut[ch+1]) + t_BdDDKp.GetEntries(cut_string+channel_cut[ch+1]) + t_BsDDKp.GetEntries(cut_string+channel_cut[ch+1])
        N_den_BDDKp = t_BuDDKp.GetEntries(cut_string) + t_BdDDKp.GetEntries(cut_string) + t_BsDDKp.GetEntries(cut_string)
        up_BDDKp = ROOT.TEfficiency.Wilson(N_den_BDDKp, N_num_BDDKp, 0.68, True)
        lo_BDDKp = ROOT.TEfficiency.Wilson(N_den_BDDKp, N_num_BDDKp, 0.68, False)
        if(N_den_BDDKp != 0):
            eps_BDDKp[ch] = N_num_BDDKp/N_den_BDDKp
        else:
            eps_BDDKp[ch] = 0
        eps_BDDKp_err[ch] = 0.5*(up_BDDKp - lo_BDDKp)

        # B -> DD K0
        N_num_BuDDK0 = t_BuDDK0.GetEntries(cut_string+channel_cut[ch+1])
        N_den_BuDDK0 = t_BuDDK0.GetEntries(cut_string)
        up_BuDDK0 = ROOT.TEfficiency.Wilson(N_den_BuDDK0, N_num_BuDDK0, 0.68, True)
        lo_BuDDK0 = ROOT.TEfficiency.Wilson(N_den_BuDDK0, N_num_BuDDK0, 0.68, False)
        if(N_den_BuDDK0 != 0):
            eps_BuDDK0[ch] = N_num_BuDDK0/N_den_BuDDK0
        else:
            eps_BuDDK0[ch] = 0
        eps_BuDDK0_err[ch] = 0.5*(up_BuDDK0 - lo_BuDDK0)

        # B -> DD
        N_num_BuDD = t_BuDD.GetEntries(cut_string+channel_cut[ch+1])
        N_den_BuDD = t_BuDD.GetEntries(cut_string)
        up_BuDD = ROOT.TEfficiency.Wilson(N_den_BuDD, N_num_BuDD, 0.68, True)
        lo_BuDD = ROOT.TEfficiency.Wilson(N_den_BuDD, N_num_BuDD, 0.68, False)
        if(N_den_BuDD != 0):
            eps_BuDD[ch] = N_num_BuDD/N_den_BuDD
        else:
            eps_BuDD[ch] = 0
        eps_BuDD_err[ch] = 0.5*(up_BuDD - lo_BuDD)

    return eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD, eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err


def physics_backgrounds_yields(C_values, C_errors, bdt1, bdt2):
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]

    #### All events
    C_BuD0D0Kp = 0 # 100
    C_BuDpDmKp = 0 # 101
    C_BuDsDsKp = 0 # 102
    C_BdDDKp = 0 # 110
    C_BsDDKp = 0 # 120
    C_BuDDK0 = 0 # 130
    C_BuDD = 0 # 150, 151

    C_BuD0D0Kp_err = 0 # 100
    C_BuDpDmKp_err = 0 # 101
    C_BuDsDsKp_err = 0 # 102
    C_BdDDKp_err = 0 # 110
    C_BsDDKp_err = 0 # 120
    C_BuDDK0_err = 0 # 130
    C_BuDD_err = 0 # 150, 151

    s = 0
    for i in range(len(all_species)):
        species = all_species[i]

        for j in range(4):
            if(species == 100):
                C_BuD0D0Kp += C_values[s,j]
                C_BuD0D0Kp_err += C_errors[s,j]**2
            elif(species == 101):
                C_BuDpDmKp += C_values[s,j]
                C_BuDpDmKp_err += C_errors[s,j]**2
            elif(species == 102):
                C_BuDsDsKp += C_values[s,j]
                C_BuDsDsKp_err += C_errors[s,j]**2
            elif(species == 110):
                C_BdDDKp += C_values[s,j]
                C_BdDDKp_err += C_errors[s,j]**2
            elif(species == 120):
                C_BsDDKp += C_values[s,j]
                C_BsDDKp_err += C_errors[s,j]**2
            elif(species == 130):
                C_BuDDK0 += C_values[s,j]
                C_BuDDK0_err += C_errors[s,j]**2
            elif(species == 150):
                C_BuDD += C_values[s,j]
                C_BuDD_err += C_errors[s,j]**2

        s += 1

    C_BuD0D0Kp_err = np.sqrt(C_BuD0D0Kp_err)
    C_BuDpDmKp_err = np.sqrt(C_BuDpDmKp_err)
    C_BuDsDsKp_err = np.sqrt(C_BuDsDsKp_err)
    C_BdDDKp_err = np.sqrt(C_BdDDKp_err)
    C_BsDDKp_err = np.sqrt(C_BsDDKp_err)
    C_BuDDK0_err = np.sqrt(C_BuDDK0_err)
    C_BuDD_err = np.sqrt(C_BuDD_err)

    C_BDDKp = C_BuD0D0Kp + C_BuDpDmKp + C_BuDsDsKp + C_BdDDKp + C_BsDDKp
    C_BDDKp_err = np.sqrt( C_BuD0D0Kp_err**2 + C_BuDpDmKp_err**2 + C_BuDsDsKp_err**2 + C_BdDDKp_err**2 + C_BsDDKp_err**2 )

    C_BDDKp = ufloat(C_BDDKp, C_BDDKp_err)
    C_BuDDK0 = ufloat(C_BuDDK0, C_BuDDK0_err)
    C_BuDD = ufloat(C_BuDD, C_BuDD_err)

    c = [C_BDDKp.nominal_value, C_BuDDK0.nominal_value, C_BuDD.nominal_value]
    c_err = [C_BDDKp.std_dev, C_BuDDK0.std_dev, C_BuDD.std_dev]

    return c, c_err

def main(argv):
    start = time.time()

    bdt1 = argv[1]
    bdt2 = argv[2]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    # Signal MC
    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")

    # RS data (yields are blinded until 0.8,0.8)
    f_rs_data = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt1_0_bdt2_0.root")
    t_rs_data = f_rs_data.Get("DecayTree")

    # Combinatorial background (WS data)
    f_comb = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_comb = f_comb.Get("DecayTree")

    # B+ -> DD K+
    f_BuDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt1_0_bdt2_0.root")
    t_BuDDKp = f_BuDDKp.Get("DecayTree")

    # B0 -> DD K+
    f_BdDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt1_0_bdt2_0.root")
    t_BdDDKp = f_BdDDKp.Get("DecayTree")

    # B0s -> DD K+
    f_BsDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt1_0_bdt2_0.root")
    t_BsDDKp = f_BsDDKp.Get("DecayTree")

    # B+ -> DD K0
    f_BuDDK0 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt1_0_bdt2_0.root")
    t_BuDDK0 = f_BuDDK0.Get("DecayTree")

    # B+ -> DD 
    f_BuDD = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt1_0_bdt2_0.root")
    t_BuDD = f_BuDD.Get("DecayTree")
    
    ###################################################### 1) Retrieve histograms after BDT cuts (all histograms are normalised to unity)
    [h_sig, h_comb, h_BDDKp] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 0, channel_cut)
    [h_sig_1, h_comb_1, h_BDDKp_1] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 1, channel_cut)
    [h_sig_2, h_comb_2, h_BDDKp_2] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 2, channel_cut)
    [h_sig_3, h_comb_3, h_BDDKp_3] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 3, channel_cut)
    
    print("Saving histograms into file")
    f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/histograms_bdt1_{bdt1}_bdt2_{bdt2}.root", "RECREATE")
    f.cd()
    f.mkdir("Channel_0")
    f.mkdir("Channel_1")
    f.mkdir("Channel_2")
    f.mkdir("Channel_3")
    f.cd("Channel_0")
    h_sig.Write()
    h_comb.Write()
    h_BDDKp.Write()
    f.cd("Channel_1")
    h_sig_1.Write()
    h_comb_1.Write()
    h_BDDKp_1.Write()
    f.cd("Channel_2")
    h_sig_2.Write()
    h_comb_2.Write()
    h_BDDKp_2.Write()
    f.cd("Channel_3")
    h_sig_3.Write()
    h_comb_3.Write()
    h_BDDKp_3.Write()
    f.Close()
    #######################################################################################################################################################33

    ###################################################### 2) Retrieve histograms with errors
    [h_sig_err, h_comb_err, h_BDDKp_err] = create_histogram_errors(h_sig, h_comb, h_BDDKp, 0)
    [h_sig_err_1, h_comb_err_1, h_BDDKp_err_1] = create_histogram_errors(h_sig_1, h_comb_1, h_BDDKp_1, 1)
    [h_sig_err_2, h_comb_err_2, h_BDDKp_err_2] = create_histogram_errors(h_sig_2, h_comb_2, h_BDDKp_2, 2)
    [h_sig_err_3, h_comb_err_3, h_BDDKp_err_3] = create_histogram_errors(h_sig_3, h_comb_3, h_BDDKp_3, 3)

    print("Saving histograms errors into file")
    f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/histograms_errors_bdt1_{bdt1}_bdt2_{bdt2}.root", "RECREATE")
    f1.cd()
    f1.mkdir("Channel_0")
    f1.mkdir("Channel_1")
    f1.mkdir("Channel_2")
    f1.mkdir("Channel_3")
    f1.cd("Channel_0")
    h_sig_err.Write()
    h_comb_err.Write()
    h_BDDKp_err.Write()
    f1.cd("Channel_1")
    h_sig_err_1.Write()
    h_comb_err_1.Write()
    h_BDDKp_err_1.Write()
    f1.cd("Channel_2")
    h_sig_err_2.Write()
    h_comb_err_2.Write()
    h_BDDKp_err_2.Write()
    f1.cd("Channel_3")
    h_sig_err_3.Write()
    h_comb_err_3.Write()
    h_BDDKp_err_3.Write()
    f1.Close()
    #######################################################################################################################################################33

    ### Combinatorial background yield
    N_comb = combinatorial_background_yield(t_rs_data, t_comb, bdt1, bdt2)
    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/N_comb_bdt1_{bdt1}_bdt2_{bdt2}.npy', N_comb)

    ### Error category efficiency (eps1, eps2)
    eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD, eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err = error_category_efficiency(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, channel_cut)

    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/eff_category_value_bdt1_{bdt1}_bdt2_{bdt2}.npy', [eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD])
    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/eff_category_error_bdt1_{bdt1}_bdt2_{bdt2}.npy', [eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err])

    ### Physics backgrounds C values
    C_values = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/C_bdt1_{bdt1}_bdt2_{bdt2}.npy')
    C_errors = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/C_err_bdt1_{bdt1}_bdt2_{bdt2}.npy')

    c, c_err = physics_backgrounds_yields(C_values, C_errors, bdt1, bdt2)

    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/C_bdt1_{bdt1}_bdt2_{bdt2}.npy", c)
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/C_err_bdt1_{bdt1}_bdt2_{bdt2}.npy", c_err)

    end = time.time()
    print(f"Elapsed time: {end - start:.2f} seconds")

if __name__ == "__main__":
    main(sys.argv)