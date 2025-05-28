import sys
import ROOT
import numpy as np
from uncertainties import ufloat
import time

channel_cut = {0: "", 1: " && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)", 2: " && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)", 3: " && (df_Bp_MERR > 250)"}
error_threshold = 0.05

def eps_error(Num, Den):
    return (Num/Den)*np.sqrt( 1/Num + 1/Den )

def number_of_bins(t_sig, ch, bdt1, bdt2, channel_cut):
    h = ROOT.TH1D(f"h_{ch}", f"h_{ch}", 100, 4000, 8000)

    cut_string = channel_cut[ch]

    bdt_fix = 0.9
    if((bdt1 >= bdt_fix) and (bdt2 >= bdt_fix)):
        t_sig.Draw(f"df_Bp_M >> h_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    bin1 = h.FindFirstBinAbove(h.GetMaximum()/2)
    bin2 = h.FindLastBinAbove(h.GetMaximum()/2)
    fwhm = h.GetBinCenter(bin2) - h.GetBinCenter(bin1)

    sigma = fwhm/2.4

    # h1 = ROOT.TH1D(f"h1_{ch}", f"h1_{ch}", 100, 0, 1000)
    # if((bdt1 >= bdt_fix) and (bdt2 >= bdt_fix)):
    #     t_sig.Draw(f"df_Bp_MERR >> h1_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
    # else:
    #     t_sig.Draw(f"df_Bp_MERR >> h1_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    # mass_resolution = h1.GetMean()
    # print("ch = ", ch, " | sigma = ", mass_resolution)

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
    h_BDDKp = ROOT.TH1D(f"h_BDDKp_{ch}", f"h_BDDKp_{ch}", nbins, 4000, 8000) # the B -> DD K+ backgrounds are merged because they have a similar shape
    h_BuDDKp = ROOT.TH1D(f"h_BuDDKp_{ch}", f"h_BuDDKp_{ch}", nbins, 4000, 8000)
    h_BdDDKp = ROOT.TH1D(f"h_BdDDKp_{ch}", f"h_BdDDKp_{ch}", nbins, 4000, 8000)
    h_BsDDKp = ROOT.TH1D(f"h_BsDDKp_{ch}", f"h_BsDDKp_{ch}", nbins, 4000, 8000)
    h_BuDDK0 = ROOT.TH1D(f"h_BuDDK0_{ch}", f"h_BuDDK0_{ch}", nbins, 4000, 8000)
    h_BuDD = ROOT.TH1D(f"h_BuDD_{ch}", f"h_BuDD_{ch}", nbins, 4000, 8000)

    cut_string = channel_cut[ch]

    bdt_fix = 0.9
    if((bdt1 >= bdt_fix) and (bdt2 >= bdt_fix)):
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_comb.Draw(f"df_Bp_M >> h_comb_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDDKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdDDKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsDDKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BuDDK0.Draw(f"df_Bp_M >> h_BuDDK0_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
        t_BuDD.Draw(f"df_Bp_M >> h_BuDD_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt_fix,bdt_fix)+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_comb.Draw(f"df_Bp_M >> h_comb_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDDKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdDDKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsDDKp_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BuDDK0.Draw(f"df_Bp_M >> h_BuDDK0_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)
        t_BuDD.Draw(f"df_Bp_M >> h_BuDD_{ch}", "(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+cut_string)

    h_sig.Sumw2()
    h_comb.Sumw2()
    h_BDDKp.Sumw2()
    h_BuDDKp.Sumw2()
    h_BdDDKp.Sumw2()
    h_BsDDKp.Sumw2()
    h_BuDDK0.Sumw2()
    h_BuDD.Sumw2()

    # Add B -> DD K+ shapes
    h_BDDKp.Add(h_BuDDKp)
    h_BDDKp.Add(h_BdDDKp)
    h_BDDKp.Add(h_BsDDKp)

    if(h_sig.Integral() != 0):
        h_sig.Scale(1/h_sig.Integral())
    if(h_comb.Integral() != 0):
        h_comb.Scale(1/h_comb.Integral())
    if(h_BDDKp.Integral() != 0):
        h_BDDKp.Scale(1/h_BDDKp.Integral())
    if(h_BuDDK0.Integral() != 0):
        h_BuDDK0.Scale(1/h_BuDDK0.Integral())
    if(h_BuDD.Integral() != 0):
        h_BuDD.Scale(1/h_BuDD.Integral())

    return [h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD]

def create_histogram_errors(h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD, ch):
    nbins = h_sig.GetNbinsX()

    h_sig_err = ROOT.TH1D(f"h_sig_err_{ch}", f"h_sig_err_{ch}", nbins, 4000, 8000)
    h_comb_err = ROOT.TH1D(f"h_comb_err_{ch}", f"h_comb_err_{ch}", nbins, 4000, 8000)
    h_BDDKp_err = ROOT.TH1D(f"h_BDDKp_err_{ch}", f"h_BDDKp_err_{ch}", nbins, 4000, 8000)
    h_BuDDK0_err = ROOT.TH1D(f"h_BuDDK0_err_{ch}", f"h_BuDDK0_err_{ch}", nbins, 4000, 8000)
    h_BuDD_err = ROOT.TH1D(f"h_BuDD_err_{ch}", f"h_BuDD_err_{ch}", nbins, 4000, 8000)

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

        # BDDKp
        if(h_BDDKp.GetBinContent(i+1) == 0):
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), 0.0) 
        elif(h_BDDKp.GetBinError(i+1)/h_BDDKp.GetBinContent(i+1) < error_threshold):
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), 0.0) 
        else:
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), h_BDDKp.GetBinError(i+1)) 
        h_BDDKp_err.SetBinError(i+1, 0) 

        # BuDDK0
        if(h_BuDDK0.GetBinContent(i+1) == 0):
            h_BuDDK0_err.Fill(h_BuDDK0.GetBinCenter(i+1), 0.0)
        elif(h_BuDDK0.GetBinError(i+1)/h_BuDDK0.GetBinContent(i+1) < error_threshold):
            h_BuDDK0_err.Fill(h_BuDDK0.GetBinCenter(i+1), 0.0) 
        else:
            h_BuDDK0_err.Fill(h_BuDDK0.GetBinCenter(i+1), h_BuDDK0.GetBinError(i+1)) 
        h_BuDDK0_err.SetBinError(i+1, 0) 

        # BuDD
        if(h_BuDD.GetBinContent(i+1) == 0):
            h_BuDD_err.Fill(h_BuDD.GetBinCenter(i+1), 0.0) 
        elif(h_BuDD.GetBinError(i+1)/h_BuDD.GetBinContent(i+1) < error_threshold):
            h_BuDD_err.Fill(h_BuDD.GetBinCenter(i+1), 0.0) 
        else:
            h_BuDD_err.Fill(h_BuDD.GetBinCenter(i+1), h_BuDD.GetBinError(i+1)) 
        h_BuDD_err.SetBinError(i+1, 0) 

        h_sig_err.Sumw2()
        h_comb_err.Sumw2()
        h_BDDKp_err.Sumw2()
        h_BuDDK0_err.Sumw2()
        h_BuDD_err.Sumw2()

        # print("SIGNAL : ", h_sig.GetBinError(i+1), " | BACKGROUND : ", h_comb.GetBinError(i+1))
    
    return [h_sig_err, h_comb_err, h_BDDKp_err, h_BuDDK0_err, h_BuDD_err]


def retrieve_yields(t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, ch, channel_cut, N_BDDKp, N_BuDDK0, N_BuDD):
    # 2) Estimate and save all needed yields after the BDT cuts
    # Combinatorial background yield
    eps_ws_den = t_comb.GetEntries()
    eps_rs_den = t_rs_data.GetEntries()

    eps_num_bdt = t_comb.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+channel_cut[ch])

    if((bdt1 <= 0.8) and (bdt2 <= 0.8)):
        eps_ws_num = t_comb.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+channel_cut[ch])
        eps_rs_num = t_rs_data.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+channel_cut[ch])
    else:
        eps_ws_num = t_comb.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(0.8,0.8)+channel_cut[ch])
        eps_rs_num = t_rs_data.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(0.8,0.8)+channel_cut[ch])

    eps_ws = eps_ws_num/eps_ws_den
    eps_rs = eps_rs_num/eps_rs_den

    r = eps_rs/eps_ws

    eps_ws_bdt = eps_num_bdt/eps_ws_den

    N_comb = 19006409*eps_ws_bdt*r
    N_comb = ufloat(N_comb, np.sqrt(N_comb))

    # Physics backgrounds
    if(ch != 0):
        cut_string = channel_cut[ch].replace("&&", "", 1)

        N_den_BDDKp = t_BuDDKp.GetEntries() + t_BdDDKp.GetEntries() + t_BsDDKp.GetEntries()
        N_num_BDDKp = t_BuDDKp.GetEntries(cut_string) + t_BdDDKp.GetEntries(cut_string) + t_BsDDKp.GetEntries(cut_string)

        N_den_BuDDK0 = t_BuDDK0.GetEntries()
        N_num_BuDDK0 = t_BuDDK0.GetEntries(cut_string) 

        N_den_BuDD = t_BuDD.GetEntries()
        N_num_BuDD = t_BuDD.GetEntries(cut_string) 

        eps_BDDKp = ufloat( N_num_BDDKp/N_den_BDDKp, eps_error(N_num_BDDKp, N_den_BDDKp) )
        eps_BuDDK0 = ufloat( N_num_BuDDK0/N_den_BuDDK0, eps_error(N_num_BuDDK0, N_den_BuDDK0) )

        eps_BuDD = ufloat( N_num_BuDD/N_den_BuDD, eps_error(N_num_BuDD, N_den_BuDD) )

        N_BDDKp = N_BDDKp*eps_BDDKp
        N_BuDDK0 = N_BuDDK0*eps_BuDDK0
        N_BuDD = N_BuDD*eps_BuDD

    return [N_comb, N_BDDKp, N_BuDDK0, N_BuDD]


def retrieve_constant(t_sig, bdt1, bdt2, ch, channel_cut, bdt_independent_constant, N_den):

    N_num = t_sig.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) ".format(bdt1,bdt2)+channel_cut[ch])

    upper = ROOT.TEfficiency.Wilson(N_den, N_num, 0.68, True)
    lower = ROOT.TEfficiency.Wilson(N_den, N_num, 0.68, False)
    eps_s_err = 0.5*(upper - lower)

    eps_s = ufloat( N_num/N_den, eps_s_err )

    c = eps_s/bdt_independent_constant

    return c

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
    [h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 0, channel_cut)
    [h_sig_1, h_comb_1, h_BDDKp_1, h_BuDDK0_1, h_BuDD_1] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 1, channel_cut)
    [h_sig_2, h_comb_2, h_BDDKp_2, h_BuDDK0_2, h_BuDD_2] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 2, channel_cut)
    [h_sig_3, h_comb_3, h_BDDKp_3, h_BuDDK0_3, h_BuDD_3] = create_histograms(t_sig, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 3, channel_cut)
    
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
    h_BuDDK0.Write()
    h_BuDD.Write()
    f.cd("Channel_1")
    h_sig_1.Write()
    h_comb_1.Write()
    h_BDDKp_1.Write()
    h_BuDDK0_1.Write()
    h_BuDD_1.Write()
    f.cd("Channel_2")
    h_sig_2.Write()
    h_comb_2.Write()
    h_BDDKp_2.Write()
    h_BuDDK0_2.Write()
    h_BuDD_2.Write()
    f.cd("Channel_3")
    h_sig_3.Write()
    h_comb_3.Write()
    h_BDDKp_3.Write()
    h_BuDDK0_3.Write()
    h_BuDD_3.Write()
    f.Close()
    #######################################################################################################################################################33

    ###################################################### 2) Retrieve histograms with errors
    [h_sig_err, h_comb_err, h_BDDKp_err, h_BuDDK0_err, h_BuDD_err] = create_histogram_errors(h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD, 0)
    [h_sig_err_1, h_comb_err_1, h_BDDKp_err_1, h_BuDDK0_err_1, h_BuDD_err_1] = create_histogram_errors(h_sig_1, h_comb_1, h_BDDKp_1, h_BuDDK0_1, h_BuDD_1, 1)
    [h_sig_err_2, h_comb_err_2, h_BDDKp_err_2, h_BuDDK0_err_2, h_BuDD_err_2] = create_histogram_errors(h_sig_2, h_comb_2, h_BDDKp_2, h_BuDDK0_2, h_BuDD_2, 2)
    [h_sig_err_3, h_comb_err_3, h_BDDKp_err_3, h_BuDDK0_err_3, h_BuDD_err_3] = create_histogram_errors(h_sig_3, h_comb_3, h_BDDKp_3, h_BuDDK0_3, h_BuDD_3, 3)

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
    h_BuDDK0_err.Write()
    h_BuDD_err.Write()
    f1.cd("Channel_1")
    h_sig_err_1.Write()
    h_comb_err_1.Write()
    h_BDDKp_err_1.Write()
    h_BuDDK0_err_1.Write()
    h_BuDD_err_1.Write()
    f1.cd("Channel_2")
    h_sig_err_2.Write()
    h_comb_err_2.Write()
    h_BDDKp_err_2.Write()
    h_BuDDK0_err_2.Write()
    h_BuDD_err_2.Write()
    f1.cd("Channel_3")
    h_sig_err_3.Write()
    h_comb_err_3.Write()
    h_BDDKp_err_3.Write()
    h_BuDDK0_err_3.Write()
    h_BuDD_err_3.Write()
    f1.Close()
    #######################################################################################################################################################33

    ###################################################### 3) Retrieve yields (ufloats)
    # Physics backgrounds yields
    yield_estimates = np.load(f'/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{bdt1}_bdt2_{bdt2}.npy')
    yield_errors = np.load(f'/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_errors_bdt1_{bdt1}_bdt2_{bdt2}.npy')
    # print(yield_estimates)

    all_species = [100, 101, 102, 110, 120, 130, 150, 151]

    #### All events
    N_BuDDKp = 0
    N_BdDDKp = 0
    N_BsDDKp = 0
    N_BuDDK0 = 0
    N_BuDD = 0

    N_BuDDKp_err = 0
    N_BdDDKp_err = 0
    N_BsDDKp_err = 0
    N_BuDDK0_err = 0
    N_BuDD_err = 0

    s = 0
    for i in range(len(all_species)):
        species = all_species[i]

        for j in range(4):
            if((species == 100) or (species == 101) or (species == 102)):
                N_BuDDKp += yield_estimates[s,j]
                N_BuDDKp_err += yield_errors[s,j]**2
            elif(species == 110):
                N_BdDDKp += yield_estimates[s,j]
                N_BdDDKp_err += yield_errors[s,j]**2
            elif(species == 120):
                N_BsDDKp += yield_estimates[s,j]
                N_BsDDKp_err += yield_errors[s,j]**2
            elif(species == 130):
                N_BuDDK0 += yield_estimates[s,j]
                N_BuDDK0_err += yield_errors[s,j]**2
            elif(species == 150):
                N_BuDD += yield_estimates[s,j]
                N_BuDD_err += yield_errors[s,j]**2

        s += 1

    N_BuDDKp_err = np.sqrt(N_BuDDKp_err)
    N_BdDDKp_err = np.sqrt(N_BdDDKp_err)
    N_BsDDKp_err = np.sqrt(N_BsDDKp_err)
    N_BuDDK0_err = np.sqrt(N_BuDDK0_err)
    N_BuDD_err = np.sqrt(N_BuDD_err)

    N_BDDKp = N_BuDDKp+N_BdDDKp+N_BsDDKp
    N_BDDKp_err = np.sqrt( N_BuDDKp_err**2 + N_BdDDKp_err**2 + N_BsDDKp_err**2 )

    N_BDDKp = ufloat( N_BDDKp, N_BDDKp_err )
    N_BuDDK0 = ufloat( N_BuDDK0, N_BuDDK0_err )
    N_BuDD = ufloat( N_BuDD, N_BuDD_err )

    [N_comb, N_BDDKp, N_BuDDK0, N_BuDD] = retrieve_yields(t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 0, channel_cut, N_BDDKp, N_BuDDK0, N_BuDD)
    [N_comb_1, N_BDDKp_1, N_BuDDK0_1, N_BuDD_1] = retrieve_yields(t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 1, channel_cut, N_BDDKp, N_BuDDK0, N_BuDD)
    [N_comb_2, N_BDDKp_2, N_BuDDK0_2, N_BuDD_2] = retrieve_yields(t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 2, channel_cut, N_BDDKp, N_BuDDK0, N_BuDD)
    [N_comb_3, N_BDDKp_3, N_BuDDK0_3, N_BuDD_3] = retrieve_yields(t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt1, bdt2, 3, channel_cut, N_BDDKp, N_BuDDK0, N_BuDD)

    yields_all_events = [ [N_comb.nominal_value, N_BDDKp.nominal_value, N_BuDDK0.nominal_value, N_BuDD.nominal_value], [N_comb_1.nominal_value, N_BDDKp_1.nominal_value, N_BuDDK0_1.nominal_value, N_BuDD_1.nominal_value], [N_comb_2.nominal_value, N_BDDKp_2.nominal_value, N_BuDDK0_2.nominal_value, N_BuDD_2.nominal_value], [N_comb_3.nominal_value, N_BDDKp_3.nominal_value, N_BuDDK0_3.nominal_value, N_BuDD_3.nominal_value] ]
    yield_errors_all_events = [ [N_BDDKp.std_dev, N_BuDDK0.std_dev, N_BuDD.std_dev], [N_BDDKp_1.std_dev, N_BuDDK0_1.std_dev, N_BuDD_1.std_dev], [N_BDDKp_2.std_dev, N_BuDDK0_2.std_dev, N_BuDD_2.std_dev], [N_BDDKp_3.std_dev, N_BuDDK0_3.std_dev, N_BuDD_3.std_dev] ] 

    print("Saving yield values and errors into files")
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/yield_values_bdt1_{bdt1}_bdt2_{bdt2}.npy", yields_all_events)
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/yield_errors_bdt1_{bdt1}_bdt2_{bdt2}.npy", yield_errors_all_events)
    #######################################################################################################################################################33

    ###################################################### 4) Retrieve constant term in branching fraction (ufloat)
    fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy")
    fixed_value_error = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs_error.npy")
    bdt_independent_constant = ufloat( fixed_value, fixed_value_error )
    N_den = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_s_den.npy')

    c = retrieve_constant(t_sig, bdt1, bdt2, 0, channel_cut, bdt_independent_constant, N_den)
    c_1 = retrieve_constant(t_sig, bdt1, bdt2, 1, channel_cut, bdt_independent_constant, N_den)
    c_2 = retrieve_constant(t_sig, bdt1, bdt2, 2, channel_cut, bdt_independent_constant, N_den)
    c_3 = retrieve_constant(t_sig, bdt1, bdt2, 3, channel_cut, bdt_independent_constant, N_den)

    c_values = [c.nominal_value, c_1.nominal_value, c_2.nominal_value, c_3.nominal_value]
    c_errors = [c.std_dev, c_1.std_dev, c_2.std_dev, c_3.std_dev]

    print("Saving constant value in the branching fraction")
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/c_value_bdt1_{bdt1}_bdt2_{bdt2}.npy", c_values)
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/c_error_bdt1_{bdt1}_bdt2_{bdt2}.npy", c_errors)

    end = time.time()
    print(f"Elapsed time: {end - start:.2f} seconds")

if __name__ == "__main__":
    main(sys.argv)