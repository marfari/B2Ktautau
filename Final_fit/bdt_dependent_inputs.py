import sys
import ROOT
import numpy as np
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy

# Signal region definition (xmin,xmax)
xmin = [5000, 4800, 5000, 5000]
xmax = [6500, 6000, 6500, 7500]

# Error category definition
channel_cut = {0: "", 1: " && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)", 2: " && (df_Bp_MERR > 100) && (df_Bp_MERR <= 240)", 3: " && (df_Bp_MERR > 240)"}

# Points from which the shape is taken at higher BDT cuts
bdt_fix_sig = 0.999
bdt_fix_phys = 0.995
bdt_fix_comb = 0.995
bdt_fix_data = 0.995

# Error threshold (ShapeSys)
error_threshold = 0

unblind_rs_data = False

def number_of_bins(t_sig, ch, bdt):
    h = ROOT.TH1D(f"h_{ch}", f"h_{ch}", 100, 4000, 9000)

    if(bdt > bdt_fix_sig):
        t_sig.Draw(f"df_Bp_M >> h_{ch}", f"(BDT > {bdt_fix_sig})"+channel_cut[ch])
    else:
        t_sig.Draw(f"df_Bp_M >> h_{ch}", f"(BDT > {bdt})"+channel_cut[ch])

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

    print(f"Channel {ch}: nbins = {nbins}")

    return nbins


def create_histograms(t_sig, t_rs, t_ws, t_BDDKp, bdt, ch, A, B, C, A_sideband, C_sideband, unblind_rs_data):
    # bin width ~half of signal mass resolution
    nbins = number_of_bins(t_sig, ch, bdt)

    # Fit templates:
    # In the fit region:
    h_sig = ROOT.TH1D(f"h_sig_{ch}", f"h_sig_{ch}", nbins, 4000, 9000) # signal
    h_BDDKp = ROOT.TH1D(f"h_BDDKp_{ch}", f"h_BDDKp_{ch}", nbins, 4000, 9000) # B -> DD K+ physics backgrounds
    h_downward = ROOT.TH1D(f"h_downward_{ch}", f"h_downward_{ch}", nbins, 4000, 9000)  # downward shape = WS data
    h_data = ROOT.TH1D(f"h_data_{ch}", f"h_data_{ch}", nbins, 4000, 9000) # RS data

    # In the sidebands
    h_sig_sideband = ROOT.TH1D(f"h_sig_sideband_{ch}", f"h_sig_sideband_{ch}", nbins, 4000, 9000) # signal
    h_BDDKp_sideband = ROOT.TH1D(f"h_BDDKp_sideband_{ch}", f"h_BDDKp_sideband_{ch}", nbins, 4000, 9000) # B -> DD K+ physics backgrounds
    h_downward_sideband = ROOT.TH1D(f"h_downward_sideband_{ch}", f"h_downward_sideband_{ch}", nbins, 4000, 9000)  # downward shape = WS data
    h_data_sideband = ROOT.TH1D(f"h_data_sideband_{ch}", f"h_data_sideband_{ch}", nbins, 4000, 9000) # RS data

    cut_string = channel_cut[ch]
    cut_string_sidebands = channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}))"

    # Signal shape = Ktautau MC 
    if(bdt > bdt_fix_sig):
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", f"(BDT > {bdt_fix_sig})"+cut_string)
        t_sig.Draw(f"df_Bp_M >> h_sig_sideband_{ch}", f"(BDT > {bdt_fix_sig})"+cut_string_sidebands)
    else:
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", f"(BDT > {bdt})"+cut_string)
        t_sig.Draw(f"df_Bp_M >> h_sig_sideband_{ch}", f"(BDT > {bdt})"+cut_string_sidebands)

    # B -> DD K+ background shape = cocktail MC (no shape uncertainty)
    if(bdt > bdt_fix_phys):
        t_BDDKp.Draw(f"df_Bp_M >> h_BDDKp_{ch}", f"(BDT > {bdt_fix_phys}) "+cut_string)
        t_BDDKp.Draw(f"df_Bp_M >> h_BDDKp_sideband_{ch}", f"(BDT > {bdt_fix_phys}) "+cut_string_sidebands)
    else:
        t_BDDKp.Draw(f"df_Bp_M >> h_BDDKp_{ch}", f"(BDT > {bdt}) "+cut_string)
        t_BDDKp.Draw(f"df_Bp_M >> h_BDDKp_sideband_{ch}", f"(BDT > {bdt}) "+cut_string_sidebands)

    # Combinatorial downward shape = WS data (shape uncertainty)
    if(bdt > bdt_fix_comb):
        t_ws.Draw(f"df_Bp_M >> h_downward_{ch}", f"(BDT > {bdt_fix_comb})"+cut_string)
        t_ws.Draw(f"df_Bp_M >> h_downward_sideband_{ch}", f"(BDT > {bdt_fix_comb})"+cut_string_sidebands)
    else:
        t_ws.Draw(f"df_Bp_M >> h_downward_{ch}", f"(BDT > {bdt})"+cut_string)
        t_ws.Draw(f"df_Bp_M >> h_downward_sideband_{ch}", f"(BDT > {bdt})"+cut_string_sidebands)

    # RS data (when the shape is taken from a lower BDT, the histogram is scaled to the correct yield)
    if(bdt > bdt_fix_data):
        if(unblind_rs_data):
            t_rs.Draw(f"df_Bp_M >> h_data_{ch}", f"(BDT > {bdt_fix_data})"+cut_string)
            if(h_data.Integral() != 0):
                h_data.Scale(t_rs.GetEntries(f"(BDT > {bdt}) "+cut_string)/h_data.Integral())
        t_rs.Draw(f"df_Bp_M >> h_data_sideband_{ch}", f"(BDT > {bdt_fix_data})"+cut_string_sidebands)
        t_rs.Draw(f"df_Bp_M >> h_weight_sideband_{ch}", f"(BDT > {bdt_fix_data})"+cut_string_sidebands)
        if(h_data_sideband.Integral() != 0):
            h_data_sideband.Scale(t_rs.GetEntries(f"(BDT > {bdt}) "+cut_string_sidebands)/h_data_sideband.Integral())
    else:
        if(unblind_rs_data):
            t_rs.Draw(f"df_Bp_M >> h_data_{ch}", f"(BDT > {bdt})"+cut_string)
        t_rs.Draw(f"df_Bp_M >> h_data_sideband_{ch}", f"(BDT > {bdt})"+cut_string_sidebands)

    # Combinatorial background upward shape
    h_upward = h_downward.Clone(f"h_upward_{ch}")
    h_upward.SetTitle(f"h_upward_{ch}")
    h_upward_sideband = h_downward_sideband.Clone(f"h_upward_sideband_{ch}")
    h_upward_sideband.SetTitle(f"h_upward_sideband_{ch}")

    # Combinatorial background nominal shape
    h_comb = h_downward.Clone(f"h_comb_{ch}")
    h_comb.SetTitle(f"h_comb_{ch}")
    h_comb_sideband = h_downward_sideband.Clone(f"h_comb_sideband_{ch}")
    h_comb_sideband.SetTitle(f"h_comb_sideband_{ch}")

    if(h_sig.Integral() != 0):
        h_sig.Scale((A[ch]*B)/h_sig.Integral()) # Scale signal template to A*B 
    if(h_sig_sideband.Integral() != 0):
        h_sig_sideband.Scale((A_sideband[ch]*B)/h_sig_sideband.Integral())
    if(h_BDDKp.Integral() != 0):
        h_BDDKp.Scale((C[ch]*B)/h_BDDKp.Integral()) # Scale BDDK+ template to C*B 
    if(h_BDDKp_sideband.Integral() != 0):
        h_BDDKp_sideband.Scale((C_sideband[ch]*B)/h_BDDKp_sideband.Integral())
    if(h_downward.Integral() != 0):
        h_downward.Scale(1.0/h_downward.Integral())
    if(h_downward_sideband.Integral() != 0):
        h_downward_sideband.Scale(1.0/h_downward_sideband.Integral())

    ######################################### Ratio R = (RS-BDDK+)/WS in the sidebands ####################################################
    # 1) Draw the RS data shape in the sidebands
    h_weight_sideband = h_data_sideband.Clone(f"h_weight_sideband_{ch}")
    h_weight_sideband.SetTitle(f"h_weight_sideband_{ch}")

    # 2) Subtract BDDK+ histogram scaled to expected yield from RS data histogram
    h_weight_sideband.Add(h_BDDKp_sideband, -1)
    for i in range(nbins):
        if(h_weight_sideband.GetBinContent(i+1) < 0):
            h_weight_sideband.SetBinContent(i+1, 0)

    # 3) Normalise RS-BDDK+ and WS to unity
    if(h_weight_sideband.Integral() != 0):
        h_weight_sideband.Scale(1.0/h_weight_sideband.Integral())
    h_numerator_sideband = h_weight_sideband.Clone(f"h_numerator_sideband_{ch}")
    h_numerator_sideband.SetTitle(f"h_numerator_sideband_{ch}")

    # 4) Divide normalised histograms
    h_weight_sideband.Divide(h_downward_sideband)

    # Fit the RS/WS ratio in the sidebands and extrapolate fit function to the signal region
    if(h_weight_sideband.Integral() != 0):
        n = 3
        xm = 4000
        xM = 9000

        bernstein = ROOT.TF1("bernstein", 
            f"""[0]*TMath::Power((1-(x-{xm})/({xM}-{xm})),{n}) +
                [1]*3*((x-{xm})/({xM}-{xm}))*TMath::Power((1-(x-{xm})/({xM}-{xm})),{n-1}) +
                [2]*3*TMath::Power((x-{xm})/({xM}-{xm}),2)*((1-(x-{xm})/({xM}-{xm})) ) +
                [3]*TMath::Power((x-{xm})/({xM}-{xm}),3)
            """, xm, xM)

        for i in range(n+1):
            bernstein.SetParameter(i, 1.0)
            bernstein.SetParLimits(i, 0, 1e6)

        fit_result = h_weight_sideband.Fit("bernstein", "RS", "", 4000, 9000)

        fit_status = fit_result.Status()
        print("Fit status:", fit_status)
        if(fit_status != 0):
            print(f"Fit failed on channel {ch}")
            quit()
        else:
            if(fit_result.Ndf() == 0):
                chi2 = fit_result.Chi2()
            else:
                chi2 = fit_result.Chi2()/fit_result.Ndf()

            pdf = h_weight_sideband.GetFunction("bernstein")

            ROOT.gStyle.SetOptStat(0)
            c = ROOT.TCanvas("c", "c")
            c.cd()
            h_weight_sideband.GetXaxis().SetTitle("m_{B} (MeV)")
            h_weight_sideband.GetYaxis().SetTitle("RS/WS (normalised) ratio")
            h_weight_sideband.SetTitle(f"BDT > {bdt} | #chi^{{2}}/ndf = {chi2:.2f}")
            h_weight_sideband.Draw()
            pdf.Draw("same")
            c.SaveAs(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/fit_result_{ch}.pdf")

            h_weight_all = ROOT.TH1D(f"h_weight_all_{ch}", f"h_weight_all_{ch}", nbins, 4000, 9000)
            for i in range(nbins):
                bin_center = h_weight_all.GetXaxis().GetBinCenter(i+1)
                h_weight_all.SetBinContent(i+1, pdf.Eval(bin_center) )

            weights = [h_weight_all.GetBinContent(i+1) for i in range(nbins)]
            print("WEIGHTS")
            print(weights)

            h_upward.Multiply(h_weight_all)
            h_upward.Multiply(h_weight_all)
            h_comb.Multiply(h_weight_all)
            h_upward_sideband.Multiply(h_weight_all)
            h_upward_sideband.Multiply(h_weight_all)
            h_comb_sideband.Multiply(h_weight_all)

    c1 = ROOT.TCanvas("c1", "c1")
    c1.cd()
    h_numerator_sideband.GetXaxis().SetTitle("m_{B} (MeV)")
    h_numerator_sideband.GetYaxis().SetTitle(f"Normalised entries / ({5000/nbins:.1f} MeV)")
    h_numerator_sideband.SetTitle(f"BDT > {bdt}")
    h_numerator_sideband.SetLineColor(1)
    h_downward_sideband.SetLineColor(2)
    rp = ROOT.TRatioPlot(h_numerator_sideband, h_downward_sideband)
    rp.SetH1DrawOpt("E")
    rp.Draw()
    rp.GetLowerRefYaxis().SetTitle("(RS-BDDK+)/WS")
    leg = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    leg.AddEntry(h_numerator_sideband, "RS-BDDK", "fp")
    leg.AddEntry(h_downward_sideband, "WS", "fp")
    leg.SetBorderSize(0)
    leg.Draw("same")
    c1.SaveAs(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/rs_vs_ws_{ch}.pdf")

    if(h_upward.Integral() != 0):
        h_upward.Scale(1.0/h_upward.Integral())
    if(h_comb.Integral() != 0):
        h_comb.Scale(1.0/h_comb.Integral())
    if(h_upward_sideband.Integral() != 0):
        h_upward_sideband.Scale(1.0/h_upward_sideband.Integral())
    if(h_comb_sideband.Integral() != 0):
        h_comb_sideband.Scale(1.0/h_comb_sideband.Integral())

    excluded_bins = []
    for i in range(nbins):
        bin_left_edge = h_sig.GetBinLowEdge(i+1)
        bin_right_edge = bin_left_edge + h_sig.GetBinWidth(i+1)
        if((bin_right_edge >= xmin[ch]) and (bin_left_edge <= xmax[ch])):
            excluded_bins.append(i+1)

    i_min = excluded_bins[0]+1
    i_max = excluded_bins[len(excluded_bins)-1]-1

    print("Signal: ", h_sig.Integral())
    print("Combinatorial: ", h_comb.Integral())
    print("BDDK+: ", h_BDDKp.Integral())
    print("Upward: ", h_upward.Integral())
    print("Downward: ", h_downward.Integral())

    return [h_sig, h_comb, h_BDDKp, h_upward, h_downward, h_data], [h_sig_sideband, h_comb_sideband, h_BDDKp_sideband, h_upward_sideband, h_downward_sideband, h_data_sideband], i_min, i_max


def create_histogram_errors(h_sig, h_comb, h_BDDKp, h_sig_sideband, h_comb_sideband, h_BDDKp_sideband, ch):
    nbins = h_sig.GetNbinsX()

    h_sig_err = ROOT.TH1D(f"h_sig_err_{ch}", f"h_sig_err_{ch}", nbins, 4000, 9000)
    h_comb_err = ROOT.TH1D(f"h_comb_err_{ch}", f"h_comb_err_{ch}", nbins, 4000, 9000)
    h_BDDKp_err = ROOT.TH1D(f"h_BDDKp_err_{ch}", f"h_BDDKp_err_{ch}", nbins, 4000, 9000)
    h_sig_sideband_err = ROOT.TH1D(f"h_sig_sideband_err_{ch}", f"h_sig_sideband_err_{ch}", nbins, 4000, 9000)
    h_comb_sideband_err = ROOT.TH1D(f"h_comb_sideband_err_{ch}", f"h_comb_sideband_err_{ch}", nbins, 4000, 9000)
    h_BDDKp_sideband_err = ROOT.TH1D(f"h_BDDKp_sideband_err_{ch}", f"h_BDDKp_sideband_err_{ch}", nbins, 4000, 9000)

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

        if(h_sig_sideband.GetBinContent(i+1) == 0):
            h_sig_sideband_err.Fill(h_sig_sideband.GetBinCenter(i+1), 0.0) 
        elif(h_sig_sideband.GetBinError(i+1)/h_sig_sideband.GetBinContent(i+1) < error_threshold):
            h_sig_sideband_err.Fill(h_sig_sideband.GetBinCenter(i+1), 0.0) 
        else:
            h_sig_sideband_err.Fill(h_sig_sideband.GetBinCenter(i+1), h_sig_sideband.GetBinError(i+1)) 
        h_sig_sideband_err.SetBinError(i+1, 0)

        # Combinatorial background
        if(h_comb.GetBinContent(i+1) == 0):
            h_comb_err.Fill(h_comb.GetBinCenter(i+1), 0.0) 
        elif(h_comb.GetBinError(i+1)/h_comb.GetBinContent(i+1) < error_threshold):
            h_comb_err.Fill(h_comb.GetBinCenter(i+1), 0.0) 
        else:
            h_comb_err.Fill(h_comb.GetBinCenter(i+1), h_comb.GetBinError(i+1)) 
        h_comb_err.SetBinError(i+1, 0) 

        if(h_comb_sideband.GetBinContent(i+1) == 0):
            h_comb_sideband_err.Fill(h_comb_sideband.GetBinCenter(i+1), 0.0) 
        elif(h_comb_sideband.GetBinError(i+1)/h_comb_sideband.GetBinContent(i+1) < error_threshold):
            h_comb_sideband_err.Fill(h_comb_sideband.GetBinCenter(i+1), 0.0) 
        else:
            h_comb_sideband_err.Fill(h_comb_sideband.GetBinCenter(i+1), h_comb_sideband.GetBinError(i+1)) 
        h_comb_sideband_err.SetBinError(i+1, 0) 

        # B -> D D K+
        if(h_BDDKp.GetBinContent(i+1) == 0):
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), 0.0) 
        elif(h_BDDKp.GetBinError(i+1)/h_BDDKp.GetBinContent(i+1) < error_threshold):
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), 0.0) 
        else:
            h_BDDKp_err.Fill(h_BDDKp.GetBinCenter(i+1), h_BDDKp.GetBinError(i+1)) 
        h_BDDKp_err.SetBinError(i+1, 0) 

        if(h_BDDKp_sideband.GetBinContent(i+1) == 0):
            h_BDDKp_sideband_err.Fill(h_BDDKp_sideband.GetBinCenter(i+1), 0.0) 
        elif(h_BDDKp_sideband.GetBinError(i+1)/h_BDDKp_sideband.GetBinContent(i+1) < error_threshold):
            h_BDDKp_sideband_err.Fill(h_BDDKp_sideband.GetBinCenter(i+1), 0.0) 
        else:
            h_BDDKp_sideband_err.Fill(h_BDDKp_sideband.GetBinCenter(i+1), h_BDDKp_sideband.GetBinError(i+1)) 
        h_BDDKp_sideband_err.SetBinError(i+1, 0) 


    return [h_sig_err, h_comb_err, h_BDDKp_err], [h_sig_sideband_err, h_comb_sideband_err, h_BDDKp_sideband_err]


def main(argv):
    bdt = argv[1]
    bdt = float(bdt)
    n_channels = len(xmin)

    if(bdt == 0.0):
        global unblind_rs_data
        unblind_rs_data = True

    ########################################################### A = [A0, A1, A2, A3], A_sideband #####################################################
    A_fix_value = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/A_fix_value.npy')
    A_fix_error = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/A_fix_error.npy')
    A_fix = ufloat(A_fix_value, A_fix_error)

    Ns_gen = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/Ns_gen.npy')

    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root")
    t_sig = f_sig.Get("DecayTree")

    A = np.zeros(n_channels)
    A_err = np.zeros(n_channels)
    A_sideband = np.zeros(n_channels)
    A_sideband_err = np.zeros(n_channels)

    for ch in range(n_channels):
        Ns_num = t_sig.GetEntries(f"(BDT > {bdt}) "+channel_cut[ch])
        Ns_num_sideband = t_sig.GetEntries(f"(BDT > {bdt}) "+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}))")

        sig_eff = Ns_num/Ns_gen
        sig_eff_up = ROOT.TEfficiency.Wilson(Ns_gen, Ns_num, 0.68, True)
        sig_eff_down = ROOT.TEfficiency.Wilson(Ns_gen, Ns_num, 0.68, False)
        sig_eff_err = 0.5*(sig_eff_up - sig_eff_down)
        sig_eff = ufloat(sig_eff, sig_eff_err)
        
        sig_eff_sideband = Ns_num_sideband/Ns_gen
        sig_eff_sideband_up = ROOT.TEfficiency.Wilson(Ns_gen, Ns_num_sideband, 0.68, True)
        sig_eff_sideband_down = ROOT.TEfficiency.Wilson(Ns_gen, Ns_num_sideband, 0.68, False)
        sig_eff_sideband_err = 0.5*(sig_eff_sideband_up - sig_eff_sideband_down)
        sig_eff_sideband = ufloat(sig_eff_sideband, sig_eff_sideband_err)

        A[ch] = (A_fix*sig_eff).nominal_value
        A_err[ch] = (A_fix*sig_eff).std_dev
        A_sideband[ch] = (A_fix*sig_eff_sideband).nominal_value
        A_sideband_err[ch] = (A_fix*sig_eff_sideband).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A.npy', A)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_err.npy', A_err)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_sideband.npy', A_sideband)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_sideband_err.npy', A_sideband_err)
    ##############################################################################################################################################

    ######################################################### c = [C0, C1, C2, C3] ###############################################################
    ######################################################### C = [c_BuDDK+, c_BdDDK+, c_BsDDK+, c_BuDDK0, c_BuDD]  ##############################
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]
    N_species = len(all_species)
    N_components = 4

    f_BuDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root")
    f_BdDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root")
    f_BsDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root")
    f_BuDDK0 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt_0.0.root")
    f_BuDD = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt_0.0.root")

    Nb_gen = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/Nb_gen.npy')

    C_fix_values = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/C_fix_values.npy')
    C_fix_errors = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/C_fix_errors.npy')

    # C_matrix : C values for each species and component (n_species, n_components)
    C_matrix = np.zeros((N_species,N_components))
    C_matrix_errors = np.zeros((N_species,N_components))

    for i in range(N_species):
        species = all_species[i]
        n = 0

        if((species == 100) or (species == 101) or (species == 102)):
            t = f_BuDDKp.Get("DecayTree")
        elif(species == 110):
            t = f_BdDDKp.Get("DecayTree")
        elif(species == 120):
            t = f_BsDDKp.Get("DecayTree")
        elif(species == 130):
            t = f_BuDDK0.Get("DecayTree")
        elif((species == 150) or (species == 151)):
            t = f_BuDD.Get("DecayTree")

        for j in range(N_components):
            Nb_reco = t.GetEntries(f"(BDT > {bdt}) && (species == {species}) && (component == {j})")
            eps_b = Nb_reco/Nb_gen[i,j]

            eps_b_up = ROOT.TEfficiency.Wilson(Nb_gen[i,j], Nb_reco , 0.68, True)
            eps_b_down = ROOT.TEfficiency.Wilson(Nb_gen[i,j], Nb_reco, 0.68, False)
            eps_b_err = 0.5*(eps_b_up - eps_b_down)
            eps_b = ufloat(eps_b, eps_b_err)
            
            C_matrix[i,j] = (ufloat(C_fix_values[i,j], C_fix_errors[i,j])*eps_b).nominal_value
            C_matrix_errors[i,j] = (ufloat(C_fix_values[i,j], C_fix_errors[i,j])*eps_b).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_matrix_values.npy', C_matrix)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_matrix_errors.npy', C_matrix_errors)

    # C_vector : C value for each species (summing all components)
    C_vector = np.zeros(N_species)
    C_vector_errors = np.zeros(N_species)

    for i in range(N_species):
        for j in range(N_components):
            c_ij = ufloat(C_matrix[i,j], C_matrix_errors[i,j])

            C_vector[i] = (ufloat(C_vector[i], C_vector_errors[i]) + c_ij).nominal_value
            C_vector_errors[i] = (ufloat(C_vector[i], C_vector_errors[i]) + c_ij).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_vector_values.npy', C_vector)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_vector_errors.npy', C_vector_errors)

    # B -> DD K+ C values
    C_BDDKp = ufloat(0, 0)
    for i in range(N_species-3):
        C_BDDKp += ufloat(C_vector[i], C_vector_errors[i])

    t_BDDKp = ROOT.TChain("DecayTree")
    t_BDDKp.Add('/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root')
    t_BDDKp.Add('/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root')
    t_BDDKp.Add('/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root')

    C = np.zeros(n_channels)
    C_err = np.zeros(n_channels)
    C_sideband = np.zeros(n_channels)
    C_sideband_err = np.zeros(n_channels)

    for ch in range(n_channels):
        BDDKp_num = t_BDDKp.GetEntries(f"(BDT > {bdt}) "+channel_cut[ch])
        BDDKp_den = t_BDDKp.GetEntries(f"BDT > {bdt}")

        if(BDDKp_den != 0):
            BDDKp_eff = BDDKp_num/BDDKp_den
        else:
            BDDKp_eff = 0
        BDDKp_eff_up = ROOT.TEfficiency.Wilson(BDDKp_den, BDDKp_num, 0.68, True)
        BDDKp_eff_down = ROOT.TEfficiency.Wilson(BDDKp_den, BDDKp_num, 0.68, False)
        BDDKp_eff_err = 0.5*(BDDKp_eff_up - BDDKp_eff_down)
        BDDKp_eff = ufloat(BDDKp_eff, BDDKp_eff_err)

        C[ch] = (C_BDDKp*BDDKp_eff).nominal_value
        C_err[ch] = (C_BDDKp*BDDKp_eff).std_dev

        BDDKp_sideband_num = t_BDDKp.GetEntries(f"(BDT > {bdt}) "+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}))")

        if(BDDKp_den != 0):
            BDDKp_sideband_eff = BDDKp_sideband_num/BDDKp_den
        else:
            BDDKp_sideband_eff = 0
        BDDKp_sideband_eff_up = ROOT.TEfficiency.Wilson(BDDKp_den, BDDKp_sideband_num, 0.68, True)
        BDDKp_sideband_eff_down = ROOT.TEfficiency.Wilson(BDDKp_den, BDDKp_sideband_num, 0.68, False)
        BDDKp_sideband_eff_err = 0.5*(BDDKp_sideband_eff_up - BDDKp_sideband_eff_down)
        BDDKp_sideband_eff = ufloat(BDDKp_sideband_eff, BDDKp_sideband_eff_err)

        C_sideband[ch] = (C_BDDKp*BDDKp_sideband_eff).nominal_value
        C_sideband_err[ch] = (C_BDDKp*BDDKp_sideband_eff).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C.npy', C)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_err.npy', C_err)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_sideband.npy', C_sideband)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_sideband_err.npy', C_sideband_err)
    ###########################################################################################################################################

    B = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B_err.npy')
    B = ufloat(B, B_err)

    ################################################# Physics backgrounds yield estimates #####################################################
    N_phys_matrix = np.zeros((N_species,N_components))
    N_phys_matrix_err = np.zeros((N_species,N_components))

    for i in range(N_species):
        for j in range(N_components):
            N_phys_matrix[i,j] = (ufloat(C_matrix[i,j], C_matrix_errors[i,j])*B).nominal_value
            N_phys_matrix_err[i,j] = (ufloat(C_matrix[i,j], C_matrix_errors[i,j])*B).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_phys_matrix_values.npy', N_phys_matrix)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_phys_matrix_errors.npy', N_phys_matrix_err)

    N_phys_vector = np.zeros(N_species)
    N_phys_vector_err = np.zeros(N_species)

    for i in range(N_species):
        for j in range(N_components):
            n_ij = ufloat(N_phys_matrix[i,j], N_phys_matrix_err[i,j])

            N_phys_vector[i] = (ufloat(N_phys_vector[i], N_phys_vector_err[i]) + n_ij).nominal_value
            N_phys_vector_err[i] = (ufloat(N_phys_vector[i], N_phys_vector_err[i]) + n_ij).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_phys_vector_values.npy', N_phys_vector)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_phys_vector_errors.npy', N_phys_vector_err)

    N_BDDKp = C_BDDKp*B
    print("N_BDDKp = ", N_BDDKp.nominal_value, " +/- ", N_BDDKp.std_dev)

    N_phys = np.zeros(n_channels)
    N_phys_err = np.zeros(n_channels)
    N_phys_sideband = np.zeros(n_channels)
    N_phys_sideband_err = np.zeros(n_channels)

    for ch in range(n_channels):
        N_phys[ch] = (ufloat(C[ch], C_err[ch])*B).nominal_value
        N_phys_err[ch] = (ufloat(C[ch], C_err[ch])*B).std_dev
        N_phys_sideband[ch] = (ufloat(C_sideband[ch], C_sideband_err[ch])*B).nominal_value
        N_phys_sideband_err[ch] = (ufloat(C_sideband[ch], C_sideband_err[ch])*B).std_dev

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_BDDKp.npy', N_phys)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_BDDKp_err.npy', N_phys_err)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_BDDKp_sideband.npy', N_phys_sideband)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_BDDKp_sideband_err.npy', N_phys_sideband_err)
    ##########################################################################################################################################

    ################################################## Combinatorial yield estimate (= N_rs_pre_bdt*eps_ws) ##########################################################
    f_rs = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root")
    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt_0.0.root")

    t_rs = f_rs.Get("DecayTree")
    t_ws = f_ws.Get("DecayTree")

    N_rs_pre_bdt = t_rs.GetEntries()
    N_rs_pre_bdt = ufloat(N_rs_pre_bdt, np.sqrt(N_rs_pre_bdt))

    N_ws_den = t_ws.GetEntries()

    N_comb = np.zeros(n_channels)
    N_comb_err = np.zeros(n_channels)
    N_comb_sideband = np.zeros(n_channels)
    N_comb_sideband_err = np.zeros(n_channels)

    for ch in range(n_channels):
        N_ws_num = t_ws.GetEntries(f"(BDT > {bdt}) "+channel_cut[ch])

        ws_eff = N_ws_num/N_ws_den
        ws_eff_up = ROOT.TEfficiency.Wilson(N_ws_den, N_ws_num, 0.68, True)
        ws_eff_down = ROOT.TEfficiency.Wilson(N_ws_den, N_ws_num, 0.68, False)
        ws_eff_err = 0.5*(ws_eff_up - ws_eff_down)
        ws_eff = ufloat(ws_eff, ws_eff_err)

        N_ws_sideband_num = t_ws.GetEntries(f"(BDT > {bdt}) "+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}))")

        ws_sideband_eff = N_ws_sideband_num/N_ws_den
        ws_sideband_eff_up = ROOT.TEfficiency.Wilson(N_ws_den, N_ws_sideband_num, 0.68, True)
        ws_sideband_eff_down = ROOT.TEfficiency.Wilson(N_ws_den, N_ws_sideband_num, 0.68, False)
        ws_sideband_eff_err = 0.5*(ws_sideband_eff_up - ws_sideband_eff_down)
        ws_sideband_eff = ufloat(ws_sideband_eff, ws_sideband_eff_err)

        N_comb[ch] = (N_rs_pre_bdt*ws_eff).nominal_value
        N_comb_err[ch] = (N_rs_pre_bdt*ws_eff).std_dev
        N_comb_sideband[ch] = (N_rs_pre_bdt*ws_sideband_eff).nominal_value
        N_comb_sideband_err[ch] = (N_rs_pre_bdt*ws_sideband_eff).std_dev

    print("N_comb = ", N_comb[0], " +/- ", N_comb_err[0])

    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb.npy', N_comb)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_err.npy', N_comb_err)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_sideband.npy', N_comb_sideband)
    np.save(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_sideband_err.npy', N_comb_sideband_err)
    ###################################################################################################################################

    ################################################### Fit templates #################################################################

    [h_sig, h_comb, h_BDDKp, h_upward, h_downward, h_data], [h_sig_sideband, h_comb_sideband, h_BDDKp_sideband, h_upward_sideband, h_downward_sideband, h_data_sideband], i_min, i_max = create_histograms(t_sig, t_rs, t_ws, t_BDDKp, bdt, 0, A, B.nominal_value, C, A_sideband, C_sideband, unblind_rs_data)  
    [h_sig_1, h_comb_1, h_BDDKp_1, h_upward_1, h_downward_1, h_data_1], [h_sig_sideband_1, h_comb_sideband_1, h_BDDKp_sideband_1, h_upward_sideband_1, h_downward_sideband_1, h_data_sideband_1], i_min_1, i_max_1 = create_histograms(t_sig, t_rs, t_ws, t_BDDKp, bdt, 1, A, B.nominal_value, C, A_sideband, C_sideband, unblind_rs_data)  
    [h_sig_2, h_comb_2, h_BDDKp_2, h_upward_2, h_downward_2, h_data_2], [h_sig_sideband_2, h_comb_sideband_2, h_BDDKp_sideband_2, h_upward_sideband_2, h_downward_sideband_2, h_data_sideband_2], i_min_2, i_max_2 = create_histograms(t_sig, t_rs, t_ws, t_BDDKp, bdt, 2, A, B.nominal_value, C, A_sideband, C_sideband, unblind_rs_data)  
    [h_sig_3, h_comb_3, h_BDDKp_3, h_upward_3, h_downward_3, h_data_3], [h_sig_sideband_3, h_comb_sideband_3, h_BDDKp_sideband_3, h_upward_sideband_3, h_downward_sideband_3, h_data_sideband_3], i_min_3, i_max_3 = create_histograms(t_sig, t_rs, t_ws, t_BDDKp, bdt, 3, A, B.nominal_value, C, A_sideband, C_sideband, unblind_rs_data)  

    f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms.root", "RECREATE")
    f.cd()
    f.mkdir("Channel_0")
    f.mkdir("Channel_1")
    f.mkdir("Channel_2")
    f.mkdir("Channel_3")
    f.cd("Channel_0")
    h_sig.Write()
    h_comb.Write()
    h_BDDKp.Write()
    h_upward.Write()
    h_downward.Write()
    h_data.Write()
    f.cd("Channel_1")
    h_sig_1.Write()
    h_comb_1.Write()
    h_BDDKp_1.Write()
    h_upward_1.Write()
    h_downward_1.Write()
    h_data_1.Write()
    f.cd("Channel_2")
    h_sig_2.Write()
    h_comb_2.Write()
    h_BDDKp_2.Write()
    h_upward_2.Write()
    h_downward_2.Write()
    h_data_2.Write()
    f.cd("Channel_3")
    h_sig_3.Write()
    h_comb_3.Write()
    h_BDDKp_3.Write()
    h_upward_3.Write()
    h_downward_3.Write()
    h_data_3.Write()
    f.Close()

    f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms_sidebands.root", "RECREATE")
    f1.cd()
    f1.mkdir("Channel_0")
    f1.mkdir("Channel_1")
    f1.mkdir("Channel_2")
    f1.mkdir("Channel_3")
    f1.cd("Channel_0")
    h_sig_sideband.Write()
    h_comb_sideband.Write()
    h_BDDKp_sideband.Write()
    h_upward_sideband.Write()
    h_downward_sideband.Write()
    h_data_sideband.Write()
    f1.cd("Channel_1")
    h_sig_sideband_1.Write()
    h_comb_sideband_1.Write()
    h_BDDKp_sideband_1.Write()
    h_upward_sideband_1.Write()
    h_downward_sideband_1.Write()
    h_data_sideband_1.Write()
    f1.cd("Channel_2")
    h_sig_sideband_2.Write()
    h_comb_sideband_2.Write()
    h_BDDKp_sideband_2.Write()
    h_upward_sideband_2.Write()
    h_downward_sideband_2.Write()
    h_data_sideband_2.Write()
    f1.cd("Channel_3")
    h_sig_sideband_3.Write()
    h_comb_sideband_3.Write()
    h_BDDKp_sideband_3.Write()
    h_upward_sideband_3.Write()
    h_downward_sideband_3.Write()
    h_data_sideband_3.Write()
    f1.Close()

    np.save(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/signal_region_indices.npy", [ [i_min, i_min_1, i_min_2, i_max_3], [i_max, i_max_1, i_max_2, i_max_3] ])
    ##########################################################################################################################################################################################

    ################################################################## Fit templates errors ##################################################################################################
    [h_sig_err, h_comb_err, h_BDDKp_err], [h_sig_sideband_err, h_comb_sideband_err, h_BDDKp_sideband_err] = create_histogram_errors(h_sig, h_comb, h_BDDKp, h_sig_sideband, h_comb_sideband, h_BDDKp_sideband, 0)
    [h_sig_err_1, h_comb_err_1, h_BDDKp_err_1], [h_sig_sideband_err_1, h_comb_sideband_err_1, h_BDDKp_sideband_err_1] = create_histogram_errors(h_sig_1, h_comb_1, h_BDDKp_1, h_sig_sideband_1, h_comb_sideband_1, h_BDDKp_sideband_1, 1)
    [h_sig_err_2, h_comb_err_2, h_BDDKp_err_2], [h_sig_sideband_err_2, h_comb_sideband_err_2, h_BDDKp_sideband_err_2] = create_histogram_errors(h_sig_2, h_comb_2, h_BDDKp_2, h_sig_sideband_2, h_comb_sideband_2, h_BDDKp_sideband_2, 2)
    [h_sig_err_3, h_comb_err_3, h_BDDKp_err_3], [h_sig_sideband_err_3, h_comb_sideband_err_3, h_BDDKp_sideband_err_3] = create_histogram_errors(h_sig_3, h_comb_3, h_BDDKp_3, h_sig_sideband_3, h_comb_sideband_3, h_BDDKp_sideband_3, 3)

    f2 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms_errors.root", "RECREATE")
    f2.cd()
    f2.mkdir("Channel_0")
    f2.mkdir("Channel_1")
    f2.mkdir("Channel_2")
    f2.mkdir("Channel_3")
    f2.cd("Channel_0")
    h_sig_err.Write()
    h_comb_err.Write()
    h_BDDKp_err.Write()
    f2.cd("Channel_1")
    h_sig_err_1.Write()
    h_comb_err_1.Write()
    h_BDDKp_err_1.Write()
    f2.cd("Channel_2")
    h_sig_err_2.Write()
    h_comb_err_2.Write()
    h_BDDKp_err_2.Write()
    f2.cd("Channel_3")
    h_sig_err_3.Write()
    h_comb_err_3.Write()
    h_BDDKp_err_3.Write()
    f2.Close()

    f3 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms_sidebands_errors.root", "RECREATE")
    f3.cd()
    f3.mkdir("Channel_0")
    f3.mkdir("Channel_1")
    f3.mkdir("Channel_2")
    f3.mkdir("Channel_3")
    f3.cd("Channel_0")
    h_sig_sideband_err.Write()
    h_comb_sideband_err.Write()
    h_BDDKp_sideband_err.Write()
    f3.cd("Channel_1")
    h_sig_sideband_err_1.Write()
    h_comb_sideband_err_1.Write()
    h_BDDKp_sideband_err_1.Write()
    f3.cd("Channel_2")
    h_sig_sideband_err_2.Write()
    h_comb_sideband_err_2.Write()
    h_BDDKp_sideband_err_2.Write()
    f3.cd("Channel_3")
    h_sig_sideband_err_3.Write()
    h_comb_sideband_err_3.Write()
    h_BDDKp_sideband_err_3.Write()
    f3.Close()



if __name__ == "__main__":
    main(sys.argv)