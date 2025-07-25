import sys
import ROOT
import numpy as np
from uncertainties import ufloat
import time

# Definition of the signal region
xmin = [5000, 4800, 5000, 5000]
xmax = [6500, 6000, 6500, 7500]

# Definition of the error categories
channel_cut = {0: "", 1: " && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)", 2: " && (df_Bp_MERR > 100) && (df_Bp_MERR <= 240)", 3: " && (df_Bp_MERR > 240)"}

error_threshold = 0
bdt_fix_sig = 0.999
bdt_fix_phys = 0.985
bdt_fix_comb = 0.995

def number_of_bins(t_sig, ch, bdt, channel_cut):
    h = ROOT.TH1D(f"h_{ch}", f"h_{ch}", 100, 4000, 9000)

    cut_string = channel_cut[ch]

    if(bdt >= bdt_fix_sig):
        t_sig.Draw(f"df_Bp_M >> h_{ch}", f"(BDT > {bdt})"+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_{ch}", f"(BDT > {bdt})"+cut_string)

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

    if(sigma == 0):
        nbins = 1
    else:
        nbins = int(4000/(sigma/2))
    print(nbins)

    return nbins

# def merge_bins(t_sig, bdt1, bdt2, ch, x_low=4000, x_up=9000):
     
#     nbins = number_of_bins(t_sig, ch, bdt1, bdt2)

#     h = ROOT.TH1D("h", "h", nbins, 4000, 9000)

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

def create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, ch, channel_cut):
    # bin width ~half of signal mass resolution
    nbins = number_of_bins(t_sig, ch, bdt, channel_cut)

    # Save histograms after BDT cuts
    h_sig = ROOT.TH1D(f"h_sig_{ch}", f"h_sig_{ch}", nbins, 4000, 9000)
    h_comb = ROOT.TH1D(f"h_comb_{ch}", f"h_comb_{ch}", nbins, 4000, 9000)
    ## B -> D D K+
    h_BDDKp = ROOT.TH1D(f"h_BDDKp_{ch}", f"h_BDDKp_{ch}", nbins, 4000, 9000)
    h_BuD0D0Kp = ROOT.TH1D(f"h_BuD0D0Kp_{ch}", f"h_BuD0D0Kp_{ch}", nbins, 4000, 9000)
    h_BuDpDmKp = ROOT.TH1D(f"h_BuDpDmKp_{ch}", f"h_BuDpDmKp_{ch}", nbins, 4000, 9000)
    h_BuDsDsKp = ROOT.TH1D(f"h_BuDsDsKp_{ch}", f"h_BuDsDsKp_{ch}", nbins, 4000, 9000)
    h_BdD0DmKp = ROOT.TH1D(f"h_BdD0DmKp_{ch}", f"h_BdD0DmKp_{ch}", nbins, 4000, 9000)
    h_BsD0DsKp = ROOT.TH1D(f"h_BsD0DsKp_{ch}", f"h_BsD0DsKp_{ch}", nbins, 4000, 9000)

    if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
        cut_string = channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )"
    elif((fit_type == "ToyData") or (fit_type == "RSData")):
        cut_string = channel_cut[ch]
    else:
        print("Wrong fit type. Try: 'ToyDataSidebands': validation of signal region definition \n 'RSDataSidebands': extraction of Ncomb \n 'ToyData': calculation of expected upper limit (blinded) \n 'RSData': calculation of observed upper limit (unblinded) ")
        quit()

    if(bdt >= bdt_fix_sig):
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", f"(BDT > {bdt_fix_sig})"+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", f"(BDT > {bdt})"+cut_string)

    if(bdt >= bdt_fix_comb):
        if(fit_type == "RSDataSidebands"):
            t_rs_data.Draw(f"df_Bp_M >> h_comb_{ch}", f"(BDT > {bdt_fix_comb})"+cut_string)
        else:
            t_comb.Draw(f"df_Bp_M >> h_comb_{ch}", f"(BDT > {bdt_fix_comb})"+cut_string)
    else:
        if(fit_type == "RSDataSidebands"):
            t_rs_data.Draw(f"df_Bp_M >> h_comb_{ch}", f"(BDT > {bdt})"+cut_string)
        else:
            t_comb.Draw(f"df_Bp_M >> h_comb_{ch}", f"(BDT > {bdt})"+cut_string)

    if(bdt >= bdt_fix_phys):
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuD0D0Kp_{ch}", f"(BDT > {bdt_fix_phys}) && (species == 100) "+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDpDmKp_{ch}", f"(BDT > {bdt_fix_phys}) && (species == 101) "+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDsDsKp_{ch}", f"(BDT > {bdt_fix_phys}) && (species == 102) "+cut_string)
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdD0DmKp_{ch}", f"(BDT > {bdt_fix_phys})"+cut_string)
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsD0DsKp_{ch}", f"(BDT > {bdt_fix_phys})"+cut_string)
    else:
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuD0D0Kp_{ch}", f"(BDT > {bdt}) && (species == 100) "+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDpDmKp_{ch}", f"(BDT > {bdt}) && (species == 101) "+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDsDsKp_{ch}", f"(BDT > {bdt}) && (species == 102) "+cut_string)
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdD0DmKp_{ch}", f"(BDT > {bdt}) "+cut_string)
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsD0DsKp_{ch}", f"(BDT > {bdt}) "+cut_string)

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

    if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
        excluded_bins = []
        for i in range(nbins):
            bin_left_edge = h_sig.GetBinLowEdge(i+1)
            bin_right_edge = bin_left_edge + h_sig.GetBinWidth(i+1)
            if((bin_right_edge >= xmin[ch]) and (bin_left_edge <= xmax[ch])):
                excluded_bins.append(i+1)

        i_min = excluded_bins[0]+1
        i_max = excluded_bins[len(excluded_bins)-1]-1
    else:
        i_min = -1
        i_max = -1

    return [h_sig, h_comb, h_BDDKp], i_min, i_max


def create_histogram_errors(h_sig, h_comb, h_BDDKp, ch):
    nbins = h_sig.GetNbinsX()

    h_sig_err = ROOT.TH1D(f"h_sig_err_{ch}", f"h_sig_err_{ch}", nbins, 4000, 9000)
    h_comb_err = ROOT.TH1D(f"h_comb_err_{ch}", f"h_comb_err_{ch}", nbins, 4000, 9000)
    h_BDDKp_err = ROOT.TH1D(f"h_BDDKp_err_{ch}", f"h_BDDKp_err_{ch}", nbins, 4000, 9000)

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


def combinatorial_background_yield(fit_type, t_rs_data, t_comb, bdt, N_BDDKp):
    if(fit_type == "ToyDataSidebands"):
        # Combinatorial background yield (WS data estimate)
        eps_ws_den = t_comb.GetEntries()
        eps_rs_den = t_rs_data.GetEntries()

        eps_num_bdt = t_comb.GetEntries(f"BDT > {bdt}")

        if(bdt <= 0.9):
            eps_ws_num = t_comb.GetEntries(f"BDT > {bdt}")
            eps_rs_num = t_rs_data.GetEntries(f"BDT > {bdt}")
        else:
            eps_ws_num = t_comb.GetEntries(f"BDT > {0.9}")
            eps_rs_num = t_rs_data.GetEntries(f"BDT > {0.9}")

        eps_ws = eps_ws_num/eps_ws_den
        eps_ws_up = ROOT.TEfficiency.Wilson(eps_ws_den, eps_ws_num, 0.68, True)
        eps_ws_down = ROOT.TEfficiency.Wilson(eps_ws_den, eps_ws_num, 0.68, False)
        eps_ws_err = 0.5*(eps_ws_up-eps_ws_down)
        eps_ws = ufloat(eps_ws, eps_ws_err)

        eps_rs = eps_rs_num/eps_rs_den
        eps_rs_up = ROOT.TEfficiency.Wilson(eps_rs_den, eps_rs_num, 0.68, True)
        eps_rs_down = ROOT.TEfficiency.Wilson(eps_rs_den, eps_rs_num, 0.68, False)
        eps_rs_err = 0.5*(eps_rs_up-eps_rs_down)
        eps_rs = ufloat(eps_rs, eps_rs_err)

        r = eps_rs/eps_ws
        print("eps_rs / eps_ws = ", r)

        eps_ws_bdt = eps_num_bdt/eps_ws_den
        eps_ws_bdt_up = ROOT.TEfficiency.Wilson(eps_ws_den, eps_num_bdt, 0.68, True)
        eps_ws_bdt_down = ROOT.TEfficiency.Wilson(eps_ws_den, eps_num_bdt, 0.68, False)
        eps_ws_bdt_err = 0.5*(eps_ws_bdt_up-eps_ws_bdt_down)
        eps_ws_bdt = ufloat(eps_ws_bdt, eps_ws_bdt_err)
        print("eps_ws = ", eps_ws_bdt)

        n_rs_prebdt_value = t_rs_data.GetEntries()

        n_rs_prebdt = ufloat(n_rs_prebdt_value, np.sqrt(n_rs_prebdt_value))
        print("N_rs_prebdt = ", n_rs_prebdt)

        N_comb = n_rs_prebdt*eps_ws_bdt*r

    elif(fit_type == "RSDataSidebands"):
        N_comb = t_rs_data.GetEntries(f"(BDT > {bdt})") - N_BDDKp
        N_comb = ufloat(N_comb, np.sqrt(N_comb))

    elif(fit_type == "ToyData"):
        N_comb_sidebands = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{bdt}/ncomb_value.npy')

        # I_sidebands = h_comb.Integral(h_comb.FindBin(4000), h_comb.FindBin(xmin[0])) + h_comb.Integral(h_comb.FindBin(xmax[0]), h_comb.FindBin(9000))
        # I_fit_region = h_comb.Integral(h_comb.FindBin(4000),h_comb.FindBin(9000))

        if(t_comb.GetEntries(f"BDT > {bdt}") != 0):
            eps_comb = (t_comb.GetEntries(f"(BDT > {bdt}) && ((df_Bp_M < {xmin[0]}) || (df_Bp_M > {xmax[0]}))"))/(t_comb.GetEntries(f"BDT > {bdt}"))
        else:
            eps_comb = 0

        if(eps_comb != 0):
            N_comb = N_comb_sidebands/eps_comb
        else:
            N_comb = 0
        N_comb = ufloat(N_comb, np.sqrt(N_comb))

    return N_comb.nominal_value, N_comb.std_dev


def physics_backgrounds_yields(C_values, C_errors):
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


def channel_efficiency(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, channel_cut):
    n = len(channel_cut)

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

    bdt_string = f"(BDT > {bdt})"

    for ch in range(n):
        cut_string = bdt_string + channel_cut[ch]
        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            cut_string += f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}))"

        # signal
        N_num_sig = t_sig.GetEntries(cut_string)
        N_den_sig = t_sig.GetEntries(bdt_string)
        up_sig = ROOT.TEfficiency.Wilson(N_den_sig, N_num_sig, 0.68, True)
        lo_sig = ROOT.TEfficiency.Wilson(N_den_sig, N_num_sig, 0.68, False)
        if(N_den_sig != 0):
            eps_sig[ch] = N_num_sig/N_den_sig
        else:
            eps_sig[ch] = 0
        eps_sig_err[ch] = 0.5*(up_sig - lo_sig)

        # combinatorial
        if(fit_type == "RSDataSidebands"):
            N_num_comb = t_rs_data.GetEntries(cut_string)
            N_den_comb = t_rs_data.GetEntries(bdt_string)
            up_comb = ROOT.TEfficiency.Wilson(N_den_comb, N_num_comb, 0.68, True)
            lo_comb = ROOT.TEfficiency.Wilson(N_den_comb, N_num_comb, 0.68, False)
            if(N_den_comb != 0):
                eps_comb[ch] = N_num_comb/N_den_comb
            else:
                eps_comb[ch] = 0
            eps_comb_err[ch] = 0.5*(up_comb - lo_comb)
        else:
            N_num_comb = t_comb.GetEntries(cut_string)
            N_den_comb = t_comb.GetEntries(bdt_string)
            up_comb = ROOT.TEfficiency.Wilson(N_den_comb, N_num_comb, 0.68, True)
            lo_comb = ROOT.TEfficiency.Wilson(N_den_comb, N_num_comb, 0.68, False)
            if(N_den_comb != 0):
                eps_comb[ch] = N_num_comb/N_den_comb
            else:
                eps_comb[ch] = 0
            eps_comb_err[ch] = 0.5*(up_comb - lo_comb)

        # B -> DD K+
        N_num_BDDKp = t_BuDDKp.GetEntries(cut_string) + t_BdDDKp.GetEntries(cut_string) + t_BsDDKp.GetEntries(cut_string)
        N_den_BDDKp = t_BuDDKp.GetEntries(bdt_string) + t_BdDDKp.GetEntries(bdt_string) + t_BsDDKp.GetEntries(bdt_string)
        up_BDDKp = ROOT.TEfficiency.Wilson(N_den_BDDKp, N_num_BDDKp, 0.68, True)
        lo_BDDKp = ROOT.TEfficiency.Wilson(N_den_BDDKp, N_num_BDDKp, 0.68, False)
        if(N_den_BDDKp != 0):
            eps_BDDKp[ch] = N_num_BDDKp/N_den_BDDKp
        else:
            eps_BDDKp[ch] = 0
        eps_BDDKp_err[ch] = 0.5*(up_BDDKp - lo_BDDKp)

        # B -> DD K0
        N_num_BuDDK0 = t_BuDDK0.GetEntries(cut_string)
        N_den_BuDDK0 = t_BuDDK0.GetEntries(bdt_string)
        up_BuDDK0 = ROOT.TEfficiency.Wilson(N_den_BuDDK0, N_num_BuDDK0, 0.68, True)
        lo_BuDDK0 = ROOT.TEfficiency.Wilson(N_den_BuDDK0, N_num_BuDDK0, 0.68, False)
        if(N_den_BuDDK0 != 0):
            eps_BuDDK0[ch] = N_num_BuDDK0/N_den_BuDDK0
        else:
            eps_BuDDK0[ch] = 0
        eps_BuDDK0_err[ch] = 0.5*(up_BuDDK0 - lo_BuDDK0)

        # B -> DD
        N_num_BuDD = t_BuDD.GetEntries(cut_string)
        N_den_BuDD = t_BuDD.GetEntries(bdt_string)
        up_BuDD = ROOT.TEfficiency.Wilson(N_den_BuDD, N_num_BuDD, 0.68, True)
        lo_BuDD = ROOT.TEfficiency.Wilson(N_den_BuDD, N_num_BuDD, 0.68, False)
        if(N_den_BuDD != 0):
            eps_BuDD[ch] = N_num_BuDD/N_den_BuDD
        else:
            eps_BuDD[ch] = 0
        eps_BuDD_err[ch] = 0.5*(up_BuDD - lo_BuDD)

    return eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD, eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err


def main(argv):
    start = time.time()

    fit_type = argv[1]
    bdt = argv[2]
    bdt = float(bdt)

    # Signal MC
    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root")
    t_sig = f_sig.Get("DecayTree")

    # RS data (yields are blinded until 0.8,0.8)
    f_rs_data = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root")
    t_rs_data = f_rs_data.Get("DecayTree")

    # Combinatorial background (WS data)
    f_comb = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt_0.0.root")
    t_comb = f_comb.Get("DecayTree")

    # B+ -> DD K+
    f_BuDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root")
    t_BuDDKp = f_BuDDKp.Get("DecayTree")

    # B0 -> DD K+
    f_BdDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root")
    t_BdDDKp = f_BdDDKp.Get("DecayTree")

    # B0s -> DD K+
    f_BsDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root")
    t_BsDDKp = f_BsDDKp.Get("DecayTree")

    # B+ -> DD K0
    f_BuDDK0 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt_0.0.root")
    t_BuDDK0 = f_BuDDK0.Get("DecayTree")

    # B+ -> DD 
    f_BuDD = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt_0.0.root")
    t_BuDD = f_BuDD.Get("DecayTree")
    
    ###################################################### 1) Retrieve histograms after BDT cuts (all histograms are normalised to unity)
    [h_sig, h_comb, h_BDDKp], i_min, i_max = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 0, channel_cut)
    [h_sig_1, h_comb_1, h_BDDKp_1], i_min_1, i_max_1 = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 1, channel_cut)
    [h_sig_2, h_comb_2, h_BDDKp_2], i_min_2, i_max_2 = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 2, channel_cut)
    [h_sig_3, h_comb_3, h_BDDKp_3], i_min_3, i_max_3 = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 3, channel_cut)
    
    print("Saving histograms into file")
    f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/histograms.root", "RECREATE")
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
    f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/histogram_errors.root", "RECREATE")
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

    ### Physics backgrounds C values
    C_values = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{bdt}/C.npy')
    C_errors = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{bdt}/C_err.npy')

    c, c_err = physics_backgrounds_yields(C_values, C_errors)

    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/C.npy", c)
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/C_err.npy", c_err)

    ### Combinatorial background yield
    B = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy')
    N_comb, N_comb_err = combinatorial_background_yield(fit_type, t_rs_data, t_comb, bdt, c[0]/B)
    if(fit_type != "RSDataSidebands"):
        print("N_comb = ", N_comb, " +/- ", N_comb_err)
    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/N_comb.npy', [N_comb, N_comb_err])

    ### Channel efficiency
    eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD, eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err = channel_efficiency(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, channel_cut)

    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/channel_eff_value.npy', [eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD])
    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/channel_eff_error.npy', [eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err])

    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/signal_region_indices.npy", [ [i_min, i_min_1, i_min_2, i_max_3], [i_max, i_max_1, i_max_2, i_max_3] ])

    end = time.time()
    print(f"Elapsed time: {end - start:.2f} seconds")

if __name__ == "__main__":
    main(sys.argv)