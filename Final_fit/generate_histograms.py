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
bdt_fix_phys = 0.995
bdt_fix_comb = 0.995

def number_of_bins(t_sig, ch, bdt, channel_cut):
    h = ROOT.TH1D(f"h_{ch}", f"h_{ch}", 100, 4000, 9000)

    cut_string = channel_cut[ch]

    if(bdt >= bdt_fix_sig):
        t_sig.Draw(f"df_Bp_M >> h_{ch}", f"(BDT > {bdt_fix_sig})"+cut_string)
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

    print(f"Channel {ch}: nbins = {nbins}")

    return nbins


def create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, ch, channel_cut, N_BDDKp):
    # bin width ~half of signal mass resolution
    nbins = number_of_bins(t_sig, ch, bdt, channel_cut)

    # Save histograms after BDT cuts
    # Fit templates:
    # Signal
    h_sig = ROOT.TH1D(f"h_sig_{ch}", f"h_sig_{ch}", nbins, 4000, 9000)
    # B -> DD K+ physics backgrounds
    h_BDDKp = ROOT.TH1D(f"h_BDDKp_{ch}", f"h_BDDKp_{ch}", nbins, 4000, 9000)
    h_BuD0D0Kp = ROOT.TH1D(f"h_BuD0D0Kp_{ch}", f"h_BuD0D0Kp_{ch}", nbins, 4000, 9000)
    h_BuDpDmKp = ROOT.TH1D(f"h_BuDpDmKp_{ch}", f"h_BuDpDmKp_{ch}", nbins, 4000, 9000)
    h_BuDsDsKp = ROOT.TH1D(f"h_BuDsDsKp_{ch}", f"h_BuDsDsKp_{ch}", nbins, 4000, 9000)
    h_BdD0DmKp = ROOT.TH1D(f"h_BdD0DmKp_{ch}", f"h_BdD0DmKp_{ch}", nbins, 4000, 9000)
    h_BsD0DsKp = ROOT.TH1D(f"h_BsD0DsKp_{ch}", f"h_BsD0DsKp_{ch}", nbins, 4000, 9000)
    # Combinatorial background
    h_comb = ROOT.TH1D(f"h_comb_{ch}", f"h_comb_{ch}", nbins, 4000, 9000) # nominal shape = WS data * R
    h_upward = ROOT.TH1D(f"h_upward_{ch}", f"h_upward_{ch}", nbins, 4000, 9000)  # upward shape = WS data * R^2
    h_downward = ROOT.TH1D(f"h_downward_{ch}", f"h_downward_{ch}", nbins, 4000, 9000)  # downward shape = WS data
    # RS data
    h_data = ROOT.TH1D(f"h_data_{ch}", f"h_data_{ch}", nbins, 4000, 9000)
    
    if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
        cut_string = channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )"
    elif((fit_type == "ToyData") or (fit_type == "RSData")):
        cut_string = channel_cut[ch]
    else:
        print("Wrong fit type. Try: 'ToyDataSidebands': validation of signal region definition \n 'RSDataSidebands': extraction of Ncomb \n 'ToyData': calculation of expected upper limit (blinded) \n 'RSData': calculation of observed upper limit (unblinded) ")
        quit()

    # Signal shape = Ktautau MC (no shape uncertainty)
    if(bdt >= bdt_fix_sig):
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", f"(BDT > {bdt_fix_sig})"+cut_string)
    else:
        t_sig.Draw(f"df_Bp_M >> h_sig_{ch}", f"(BDT > {bdt})"+cut_string)

    # B -> DD K+ background shape = cocktail MC (no shape uncertainty)
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
    h_BDDKp.Add(h_BuD0D0Kp)
    h_BDDKp.Add(h_BuDpDmKp)
    h_BDDKp.Add(h_BuDsDsKp)
    h_BDDKp.Add(h_BdD0DmKp)
    h_BDDKp.Add(h_BsD0DsKp)

    # Combinatorial downward shape = WS data (shape uncertainty)
    if(bdt >= bdt_fix_comb):
        t_comb.Draw(f"df_Bp_M >> h_downward_{ch}", f"(BDT > {bdt_fix_comb})"+cut_string)
    else:
        t_comb.Draw(f"df_Bp_M >> h_downward_{ch}", f"(BDT > {bdt})"+cut_string)

    ############################## (RS-BDDK)/WS ratio R in the sidebands, fitted and extrapolated to the signal region ###################################
    h_weight_sideband = ROOT.TH1D(f"h_weight_sideband_{ch}", f"h_weight_sideband_{ch}", nbins, 4000, 9000)
    h_comb_sideband = ROOT.TH1D(f"h_comb_sideband_{ch}", f"h_comb_sideband_{ch}", nbins, 4000, 9000)
    h_BDDKp_sideband = ROOT.TH1D(f"h_BDDKp_sideband_{ch}", f"h_BDDKp_sideband_{ch}", nbins, 4000, 9000)
    h_BuDDKp_sideband = ROOT.TH1D(f"h_BuDDKp_sideband_{ch}", f"h_BuDDKp_sideband_{ch}", nbins, 4000, 9000)
    h_BdDDKp_sideband = ROOT.TH1D(f"h_BdDDKp_sideband_{ch}", f"h_BdDDKp_sideband_{ch}", nbins, 4000, 9000)
    h_BsDDKp_sideband = ROOT.TH1D(f"h_BsDDKp_sideband_{ch}", f"h_BsDDKp_sideband_{ch}", nbins, 4000, 9000)

    # h_weight = ROOT.TH1D(f"h_weight_{ch}", f"h_weight_{ch}", nbins, 4000, 9000)
    # h_comb_weight = ROOT.TH1D(f"h_comb_weight_{ch}", f"h_comb_weight_{ch}", nbins, 4000, 9000)
    # h_BDDKp_weight = ROOT.TH1D(f"h_BDDKp_weight_{ch}", f"h_BDDKp_weight_{ch}", nbins, 4000, 9000)
    # h_BuDDKp_weight = ROOT.TH1D(f"h_BuDDKp_weight_{ch}", f"h_BuDDKp_weight_{ch}", nbins, 4000, 9000)
    # h_BdDDKp_weight = ROOT.TH1D(f"h_BdDDKp_weight_{ch}", f"h_BdDDKp_weight_{ch}", nbins, 4000, 9000)
    # h_BsDDKp_weight = ROOT.TH1D(f"h_BsDDKp_weight_{ch}", f"h_BsDDKp_weight_{ch}", nbins, 4000, 9000)

    if(bdt >= bdt_fix_comb):
        t_rs_data.Draw(f"df_Bp_M >> h_weight_sideband_{ch}", f"(BDT > {bdt_fix_comb})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")
        t_comb.Draw(f"df_Bp_M >> h_comb_sideband_{ch}", f"(BDT > {bdt_fix_comb})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")

        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDDKp_sideband_{ch}", f"(BDT > {bdt_fix_comb})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdDDKp_sideband_{ch}", f"(BDT > {bdt_fix_comb})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsDDKp_sideband_{ch}", f"(BDT > {bdt_fix_comb})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")

    else:
        t_rs_data.Draw(f"df_Bp_M >> h_weight_sideband_{ch}", f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")
        t_comb.Draw(f"df_Bp_M >> h_comb_sideband_{ch}", f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")

        t_BuDDKp.Draw(f"df_Bp_M >> h_BuDDKp_sideband_{ch}", f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")
        t_BdDDKp.Draw(f"df_Bp_M >> h_BdDDKp_sideband_{ch}", f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")
        t_BsDDKp.Draw(f"df_Bp_M >> h_BsDDKp_sideband_{ch}", f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )")

        # t_rs_data.Draw(f"df_Bp_M >> h_weight_{ch}", f"(BDT > {bdt})"+channel_cut[ch])
        # t_comb.Draw(f"df_Bp_M >> h_comb_weight_{ch}", f"(BDT > {bdt})"+channel_cut[ch])
        # t_BuDDKp.Draw(f"df_Bp_M >> h_BuDDKp_weight_{ch}", f"(BDT > {bdt})"+channel_cut[ch])
        # t_BdDDKp.Draw(f"df_Bp_M >> h_BdDDKp_weight_{ch}", f"(BDT > {bdt})"+channel_cut[ch])
        # t_BsDDKp.Draw(f"df_Bp_M >> h_BsDDKp_weight_{ch}", f"(BDT > {bdt})"+channel_cut[ch])

    # h_BDDKp_weight.Sumw2()
    # h_BDDKp_weight.Add(h_BuDDKp_weight)
    # h_BDDKp_weight.Add(h_BdDDKp_weight)
    # h_BDDKp_weight.Add(h_BsDDKp_weight)
    # h_BDDKp_weight.Scale( (N_BDDKp*(h_BDDKp_weight.GetEntries()/( t_BuDDKp.GetEntries(f"(BDT > {bdt})") + t_BdDDKp.GetEntries(f"(BDT > {bdt})") + t_BsDDKp.GetEntries(f"(BDT > {bdt})") ))) /h_BDDKp_weight.Integral())
    # h_rs_data_all = h_weight.Clone("h_weight_all{ch}")
    # h_weight_all = h_weight.Clone("h_weight_all{ch}")
    # h_weight.Add(h_BDDKp_weight, -1)
    # h_weight.Scale(1.0/h_weight_sideband.Integral())
    # h_comb_weight.Scale(1.0/h_comb_sideband.Integral())
    # h_weight.Divide(h_comb_weight)
    # h_weight_all.Scale(1.0/h_weight_sideband.Integral())
    # h_weight_all.Divide(h_comb_weight)

    h_BDDKp_sideband.Add(h_BuDDKp_sideband)
    h_BDDKp_sideband.Add(h_BdDDKp_sideband)
    h_BDDKp_sideband.Add(h_BsDDKp_sideband)

    eps_BDDKp_num = t_BuDDKp.GetEntries(f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )") + t_BdDDKp.GetEntries(f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )") + t_BsDDKp.GetEntries(f"(BDT > {bdt})"+channel_cut[ch]+f" && ((df_Bp_M < {xmin[ch]}) || (df_Bp_M > {xmax[ch]}) )") 
    eps_BDDKp_den = t_BuDDKp.GetEntries(f"(BDT > {bdt})") + t_BdDDKp.GetEntries(f"(BDT > {bdt})") + t_BsDDKp.GetEntries(f"(BDT > {bdt})") 
    sigma = 0
    if(eps_BDDKp_den != 0):
        eps_BDDKp = eps_BDDKp_num/eps_BDDKp_den
        eps_BDDKp_up = ROOT.TEfficiency.Wilson(eps_BDDKp_den, eps_BDDKp_num, 0.68, True)
        eps_BDDKp_down = ROOT.TEfficiency.Wilson(eps_BDDKp_den, eps_BDDKp_num, 0.68, False)
        eps_BDDKp_err = 0.5*(eps_BDDKp_up-eps_BDDKp_down)
        eps_BDDKp = ufloat(eps_BDDKp, eps_BDDKp_err)

        n_BDDKp = N_BDDKp*eps_BDDKp
        sigma = n_BDDKp.nominal_value/n_BDDKp.std_dev
        print("Sigma = ", sigma)

        h_BDDKp_sideband.Scale(n_BDDKp.nominal_value/h_BDDKp_sideband.Integral()) # B -> DD K+ shape scaled to expected yield in RS data (yield in the sidebands)

    if(sigma > 3):
        h_weight_sideband.Add(h_BDDKp_sideband, -1)

    # Normalise histograms before dividing
    h_weight_sideband.Scale(1.0/h_weight_sideband.Integral())
    h_comb_sideband.Scale(1.0/h_comb_sideband.Integral())
    h_weight_sideband.Divide(h_comb_sideband) # R = RS / WS

    # Fit the RS/WS ratio in the sidebands and extrapolate fit function to the signal region
    # In the signal region the fitted ratio values are used to build the upward and downward variations

    fit_result = h_weight_sideband.Fit("pol3", "RS", "", 4000, 9000)

    fit_status = fit_result.Status()
    print("Fit status:", fit_status)
    if(fit_status != 0):
        print(f"Fit failed on channel {ch}")
        quit()
    if(fit_result.Ndf() == 0):
        chi2 = fit_result.Chi2()
    else:
        chi2 = fit_result.Chi2()/fit_result.Ndf()

    pdf = h_weight_sideband.GetFunction("pol3")
    pdf.SetRange(4000,9000)

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "c")
    c.cd()
    h_weight_sideband.GetXaxis().SetTitle("m_{B} (MeV)")
    h_weight_sideband.GetYaxis().SetTitle("RS/WS (normalised) ratio")
    h_weight_sideband.SetTitle(f"BDT > {bdt} | #chi^{{2}}/ndf = {chi2:.2f}")
    h_weight_sideband.Draw()
    # h_weight.GetXaxis().SetTitle("m_{B} (MeV)")
    # h_weight.GetYaxis().SetTitle("RS/WS (normalised) ratio")
    # h_weight.SetTitle(f"BDT > {bdt} | #chi^{{2}}/ndf = {chi2:.2f}")
    # h_weight.SetLineColor(4)
    # h_weight_all.SetLineColor(1)
    # h_weight.Draw()
    # h_weight_all.Draw("same")
    pdf.Draw("same")
    c.SaveAs(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/fit_result_{ch}.pdf")

    # c1 = ROOT.TCanvas("c1", "c1")
    # c1.cd()
    # h_rs_data_all.Scale(1.0/h_rs_data_all.Integral())
    # h_comb_weight.Scale(1.0/h_comb_weight.Integral())
    # h_rs_data_all.GetXaxis().SetTitle("m_{B} (MeV)")
    # h_rs_data_all.GetYaxis().SetTitle("Normalized entries")
    # h_rs_data_all.SetTitle(f"BDT > {bdt} | #chi^{{2}}/ndf = {chi2:.2f}")
    # h_rs_data_all.SetLineColor(1)
    # h_comb_weight.SetLineColor(2)
    # # h_rs_data_all.Draw()
    # # h_comb_weight.Draw("same")
    # rp = ROOT.TRatioPlot(h_rs_data_all, h_comb_weight)
    # rp.SetH1DrawOpt("E")
    # rp.Draw()
    # c1.SaveAs(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/rs_vs_ws_data_{ch}.pdf")

    h_weight_all = ROOT.TH1D(f"h_weight_all_{ch}", f"h_weight_all_{ch}", nbins, 4000, 9000)
    for i in range(nbins):
        bin_center = h_weight_all.GetXaxis().GetBinCenter(i+1)
        h_weight_all.SetBinContent(i+1, pdf.Eval(bin_center) )

    for i in range(nbins):
        if(h_weight_all.GetBinContent(i+1) < 0):
            print("Before: ", h_weight_all.GetBinContent(i+1), " | ", h_weight_all.GetBinError(i+1))
            h_weight_all.SetBinContent(i+1, h_weight_all.GetBinContent(i+1)+h_weight_all.GetBinError(i+1))
            print("After: ", h_weight_all.GetBinContent(i+1), " | ", h_weight_all.GetBinError(i+1))
            if(h_weight_all.GetBinContent(i+1) < 0):
                h_weight_all.SetBinContent(i+1, 0)

    weights = [h_weight_all.GetBinContent(i+1) for i in range(nbins)]
    # print("WEIGHTS")
    # print(weights)

    ###########################################################################################################################################
    
    h_upward = h_downward.Clone(f"h_upward_{ch}")
    h_upward.SetTitle(f"h_upward_{ch}")
    h_upward.Multiply(h_weight_all)
    h_upward.Multiply(h_weight_all)

    h_comb = h_downward.Clone(f"h_comb_{ch}")
    h_comb.SetTitle(f"h_comb_{ch}")
    h_comb.Multiply(h_weight_all)

    # RS data shape (used in the fit to the RS data sidebands and in the fit to the RS data)
    if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
        t_rs_data.Draw(f"df_Bp_M >> h_data_{ch}", f"(BDT > {bdt}) "+cut_string)

    upward_values = [h_upward.GetBinContent(i+1) for i in range(nbins)]
    downward_values = [h_downward.GetBinContent(i+1) for i in range(nbins)]
    nominal_values = [h_comb.GetBinContent(i+1) for i in range(nbins)]
    # print("NOMINAL")
    # print(nominal_values)
    # print("UPWARD")
    # print(upward_values)
    # print("DOWNWARD")
    # print(downward_values)

    # eps = 0.1
    # for i in range(nbins):
    #     if(h_sig.GetBinContent(i+1) == 0):
    #         h_sig.SetBinContent(i+1, eps)
    #         h_sig.SetBinError(i+1, np.sqrt(eps))
    #     if(h_comb.GetBinContent(i+1) == 0):
    #         h_comb.SetBinContent(i+1, eps)
    #         h_comb.SetBinError(i+1, np.sqrt(eps))
    #     if(h_BDDKp.GetBinContent(i+1) == 0):
    #         h_BDDKp.SetBinContent(i+1, eps)
    #         h_BDDKp.SetBinError(i+1, np.sqrt(eps))

    if(h_sig.Integral() != 0):
        h_sig.Scale(1/h_sig.Integral())
    if(h_comb.Integral() != 0):
        h_comb.Scale(1/h_comb.Integral())
    if(h_BDDKp.Integral() != 0):
        h_BDDKp.Scale(1/h_BDDKp.Integral())
    if(h_upward.Integral() != 0):
        h_upward.Scale(1/h_upward.Integral())
    if(h_downward.Integral() != 0):
        h_downward.Scale(1/h_downward.Integral())

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

    return [h_sig, h_comb, h_BDDKp, h_upward, h_downward, h_data], i_min, i_max


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


def combinatorial_background_yield(fit_type, t_rs_data, t_comb, bdt, B, B_err, c, c_err):

    if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands") or (fit_type == "ToyData")): # WS data scaling
        C_values_0 = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_0.0/C.npy')
        C_errors_0 = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_0.0/C_err.npy')
        c_0, c_err_0 = physics_backgrounds_yields(C_values_0, C_errors_0)
        c_0 = ufloat(c_0[0], c_err_0[0])
        B = ufloat(B, B_err)
        n_BDDKp_0 = c_0/B
    
        C_values_09 = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_0.9/C.npy')
        C_errors_09 = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_0.9/C_err.npy')
        c_09, c_err_09 = physics_backgrounds_yields(C_values_09, C_errors_09)
        c_09 = ufloat(c_09[0], c_err_09[0])
        n_BDDKp_09 = c_09/B

        eps_ws_den = t_comb.GetEntries()
        eps_rs_den = t_rs_data.GetEntries() - n_BDDKp_0.nominal_value
        eps_num_bdt = t_comb.GetEntries(f"BDT > {bdt}")

        if(bdt <= 0.9):
            eps_ws_num = t_comb.GetEntries(f"BDT > {bdt}")
            eps_rs_num = t_rs_data.GetEntries(f"BDT > {bdt}") - c/B.nominal_value
        else:
            eps_ws_num = t_comb.GetEntries(f"BDT > {0.9}")
            eps_rs_num = t_rs_data.GetEntries(f"BDT > {0.9}") - n_BDDKp_09.nominal_value

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

        n_rs_data_entries_pre_bdt = t_rs_data.GetEntries()
        n_rs_prebdt = ufloat( n_rs_data_entries_pre_bdt, np.sqrt(n_rs_data_entries_pre_bdt) ) - n_BDDKp_0

        N_comb = n_rs_prebdt*eps_ws_bdt*r

    # elif(fit_type == "ToyData"): 
    #     # from fit to RS data sidebands divided by the WS data sideband efficiency
    #     # N_comb_sidebands = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{bdt}/ncomb_value.npy')
    #     # N_comb_sidebands_err = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{bdt}/ncomb_error.npy')
    #     # N_comb_sidebands = ufloat(N_comb_sidebands, N_comb_sidebands_err)

    #     n_rs_data_sidebands = t_rs_data.GetEntries(f"(BDT > {bdt}) && ((df_Bp_M < {xmin[0]}) || (df_Bp_M > {xmax[0]}))")
    #     n_rs_data_sidebands = ufloat(n_rs_data_sidebands, np.sqrt(n_rs_data_sidebands))

    #     N_BDDKp = ufloat(c, c_err)/ufloat(B, B_err)

    #     eps_comb_num = t_comb.GetEntries(f"(BDT > {bdt}) && ((df_Bp_M < {xmin[0]}) || (df_Bp_M > {xmax[0]}))")
    #     eps_comb_den = t_comb.GetEntries(f"BDT > {bdt}")
    #     if(eps_comb_den != 0):
    #         eps_comb_num = t_comb.GetEntries(f"(BDT > {bdt}) && ((df_Bp_M < {xmin[0]}) || (df_Bp_M > {xmax[0]}))")
    #         eps_comb_den = t_comb.GetEntries(f"BDT > {bdt}")

    #         eps_comb = eps_comb_num/eps_comb_den
    #         eps_comb_up = ROOT.TEfficiency.Wilson(eps_comb_den, eps_comb_num, 0.68, True)
    #         eps_comb_down = ROOT.TEfficiency.Wilson(eps_comb_den, eps_comb_num, 0.68, False)
    #         eps_comb_err = 0.5*(eps_comb_up-eps_comb_down)
    #         eps_comb = ufloat(eps_comb, eps_comb_err)

    #         N_comb = n_rs_data_sidebands/eps_comb - N_BDDKp
    #         # N_comb = N_comb_sidebands/eps_comb
    #     else:
    #         N_comb = ufloat(0,0)

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
    
    ### Physics backgrounds C values
    C_values = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{bdt}/C.npy')
    C_errors = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{bdt}/C_err.npy')
    c, c_err = physics_backgrounds_yields(C_values, C_errors)

    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/C.npy", c)
    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/C_err.npy", c_err)

    ### Combinatorial background yield
    B = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy')
    N_comb, N_comb_err = combinatorial_background_yield(fit_type, t_rs_data, t_comb, bdt, B, B_err, c[0], c_err[0])
    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/N_comb.npy', [N_comb, N_comb_err])

    ### Channel efficiency
    eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD, eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err = channel_efficiency(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, channel_cut)

    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/channel_eff_value.npy', [eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD])
    np.save(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/channel_eff_error.npy', [eps_sig_err, eps_comb_err, eps_BDDKp_err, eps_BuDDK0_err, eps_BuDD_err])

    ###################################################### 1) Retrieve histograms after BDT cuts (all histograms are normalised to unity)
    # I want to use the upward variation as the nominal combinatorial background shape (h_comb <-> h_upward)
    N_phys = ufloat(c[0], c_err[0])/ufloat(B, B_err)
    [h_sig, h_comb, h_BDDKp, h_upward, h_downward, h_data], i_min, i_max = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 0, channel_cut, N_phys)
    [h_sig_1, h_comb_1, h_BDDKp_1, h_upward_1, h_downward_1, h_data_1], i_min_1, i_max_1 = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 1, channel_cut, N_phys)
    [h_sig_2, h_comb_2, h_BDDKp_2, h_upward_2, h_downward_2, h_data_2], i_min_2, i_max_2 = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 2, channel_cut, N_phys)
    [h_sig_3, h_comb_3, h_BDDKp_3, h_upward_3, h_downward_3, h_data_3], i_min_3, i_max_3 = create_histograms(fit_type, t_sig, t_rs_data, t_comb, t_BuDDKp, t_BdDDKp, t_BsDDKp, t_BuDDK0, t_BuDD, bdt, 3, channel_cut, N_phys)

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
    h_upward.Write()
    h_downward.Write()
    if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
        h_data.Write()
    f.cd("Channel_1")
    h_sig_1.Write()
    h_comb_1.Write()
    h_BDDKp_1.Write()
    h_upward_1.Write()
    h_downward_1.Write()
    if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
        h_data_1.Write()
    f.cd("Channel_2")
    h_sig_2.Write()
    h_comb_2.Write()
    h_BDDKp_2.Write()
    h_upward_2.Write()
    h_downward_2.Write()
    if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
        h_data_2.Write()
    f.cd("Channel_3")
    h_sig_3.Write()
    h_comb_3.Write()
    h_BDDKp_3.Write()
    h_upward_3.Write()
    h_downward_3.Write()
    if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
        h_data_3.Write()
    f.Close()

    np.save(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/signal_region_indices.npy", [ [i_min, i_min_1, i_min_2, i_max_3], [i_max, i_max_1, i_max_2, i_max_3] ])
    #######################################################################################################################################################33

    ###################################################### 2) Retrieve histograms with errors
    # I want to use the upward variation and the nominal combinatorial shape (h_comb <-> h_upward)
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

    n_BDDKp = ufloat(c[0], c_err[0])/B
    print("N_BDDKp = ", n_BDDKp.nominal_value, " +/- ", n_BDDKp.std_dev)
    if(fit_type != "RSDataEntries"):
        print("N_comb = ", N_comb, " +/- ", N_comb_err)

    end = time.time()
    print(f"Elapsed time: {end - start:.2f} seconds")

if __name__ == "__main__":
    main(sys.argv)