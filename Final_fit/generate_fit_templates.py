import sys
import ROOT
import numpy as np
import matplotlib.pyplot as plt

def eps_error(Num, Den):
    return (Num/Den)*np.sqrt( 1/Num + 1/Den )

def draw_rs_ws_bdt_eff():
    # RS data
    f_rs = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt1_0_bdt2_0.root")
    t_rs = f_rs.Get("DecayTree")

    # WS data
    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_ws = f_ws.Get("DecayTree")

    bdt_cuts = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    N = len(bdt_cuts)

    N_rs_den = t_rs.GetEntries()
    N_ws_den = t_ws.GetEntries()

    bdt1_eff_ratio = np.zeros(N)
    bdt2_eff_ratio = np.zeros(N)
    bdt_eff_ratio = np.zeros(N)

    bdt1_eff_ratio_err = np.zeros(N)
    bdt2_eff_ratio_err = np.zeros(N)
    bdt_eff_ratio_err = np.zeros(N)

    for i in range(N):
        print("BDT > {0}".format(bdt_cuts[i]))

        N1_rs_num = t_rs.GetEntries("BDT1 > {0}".format(bdt_cuts[i]))
        N2_rs_num = t_rs.GetEntries("BDT2 > {0}".format(bdt_cuts[i]))
        N_rs_num = t_rs.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt_cuts[i], bdt_cuts[i]))

        N1_ws_num = t_ws.GetEntries("BDT1 > {0}".format(bdt_cuts[i]))
        N2_ws_num = t_ws.GetEntries("BDT2 > {0}".format(bdt_cuts[i]))
        N_ws_num = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt_cuts[i], bdt_cuts[i]))

        eps1_rs = N1_rs_num/N_rs_den
        eps1_rs_err = eps_error(N1_rs_num,N_rs_den)

        eps1_ws = N1_ws_num/N_ws_den
        eps1_ws_err = eps_error(N1_ws_num,N_ws_den)

        eps2_rs = N2_rs_num/N_rs_den
        eps2_rs_err = eps_error(N2_rs_num,N_rs_den)

        eps2_ws = N2_ws_num/N_ws_den
        eps2_ws_err = eps_error(N2_ws_num,N_ws_den)

        eps_rs = N_rs_num/N_rs_den
        eps_rs_err = eps_error(N_rs_num,N_rs_den)

        eps_ws = N_ws_num/N_ws_den
        eps_ws_err = eps_error(N_ws_num,N_ws_den)

        bdt1_eff_ratio[i] = eps1_rs/eps1_ws
        bdt2_eff_ratio[i] = eps2_rs/eps2_ws
        bdt_eff_ratio[i] = eps_rs/eps_ws

        bdt1_eff_ratio_err[i] = np.sqrt( ((1/eps1_ws)**2)*(eps1_rs_err**2) + ((eps1_rs/(eps1_ws**2))**2)*(eps1_ws_err**2) )
        bdt2_eff_ratio_err[i] = np.sqrt( ((1/eps2_ws)**2)*(eps2_rs_err**2) + ((eps2_rs/(eps2_ws**2))**2)*(eps2_ws_err**2) )
        bdt_eff_ratio_err[i] = np.sqrt( ((1/eps_ws)**2)*(eps_rs_err**2) + ((eps_rs/(eps_ws**2))**2)*(eps_ws_err**2) )


    plt.errorbar(bdt_cuts, bdt1_eff_ratio, yerr=bdt1_eff_ratio_err, label="Isolation", marker="o")
    plt.errorbar(bdt_cuts, bdt2_eff_ratio, yerr=bdt2_eff_ratio_err, label="Topology", marker="o")
    plt.errorbar(bdt_cuts, bdt_eff_ratio, yerr=bdt_eff_ratio_err, label="Total", marker="o")
    plt.xlabel("BDT cut")
    plt.ylabel("RS/WS BDT efficiency ratio")
    plt.legend()
    plt.savefig("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/rs_vs_ws_bdt_eff.pdf")
    plt.clf()

def main(argv):
    nbins = 100
    # ROOT.gStyle.SetOptStat(0)

    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]
    folder_name = argv[4]
    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    # Signal 
    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_1/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")

    # Combinatorial background
    f_bkg = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_bkg = f_bkg.Get("DecayTree")

    # Signal + background templates (nominal)
    h_sig = ROOT.TH1D("h_sig", "h_sig", nbins, 4000, 8000)
    h_bkg = ROOT.TH1D("h_bkg", "h_bkg", nbins, 4000, 8000)
    
    t_sig.Draw("df_Bp_M >> h_sig", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    t_bkg.Draw("df_Bp_M >> h_bkg", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    # Scale histograms to unit area
    if((h_sig.Integral() != 0) and (h_bkg.Integral() != 0)):
        h_sig.Scale(1/h_sig.Integral())
        h_bkg.Scale(1/h_bkg.Integral())

        h_sig.Sumw2()
        h_bkg.Sumw2()

    # c = ROOT.TCanvas("c", "c")
    # c.cd()
    # h_sig.SetLineColor(4)
    # h_bkg.SetLineColor(2)
    # h_sig.SetMarkerColor(4)
    # h_bkg.SetMarkerColor(2)
    # h_sig.SetTitle("BDT1 > {0} and BDT2 > {1}".format( bdt1, bdt2 ))
    # h_sig.GetXaxis().SetTitle("m_{B} (MeV)")
    # h_sig.GetYaxis().SetTitle("Normalized entries / (40 MeV)")
    # h_sig.Draw()
    # h_bkg.Draw("same")
    # leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
    # leg.AddEntry(h_sig, "Signal")
    # leg.AddEntry(h_bkg, "Comb. background")
    # leg.Draw("same")
    # c.SaveAs("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/fit_template_histograms_bdt1_{0}_bdt2_{1}.pdf".format(bdt1, bdt2 ))

    # Signal + background uncertainties (error)
    h_sig_err = ROOT.TH1D("h_sig_err", "h_sig_err", nbins, 4000, 8000)
    h_bkg_err = ROOT.TH1D("h_bkg_err", "h_bkg_err", nbins, 4000, 8000)

    x = np.zeros(nbins)
    sig_yerr = np.zeros(nbins)
    bkg_yerr = np.zeros(nbins)

    for i in range(nbins):
        x[i] = h_sig.GetBinCenter(i)
        sig_yerr[i] = h_sig.GetBinError(i)
        bkg_yerr[i] = h_bkg.GetBinError(i)

    for i in range(nbins):
        h_sig_err.Fill(x[i], sig_yerr[i])
        h_sig_err.SetBinError(i, 0)

        h_bkg_err.Fill(x[i], bkg_yerr[i])
        h_bkg_err.SetBinError(i, 0)

    h_sig_err.SetMarkerColor(4)
    h_bkg_err.SetMarkerColor(2)

    # Data
    f_data = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data = f_data.Get("h_toy_data")

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{3}/fit_templates_bdt1_{0}_bdt2_{1}_seed_{2}.root".format(bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout.mkdir("Data")
    fout.mkdir("Signal")
    fout.mkdir("Background")

    fout.cd("Data")
    h_data.SetTitle("data")
    h_data.SetName("data")
    h_data.Write()

    fout.cd("Signal")
    h_sig.SetName("nominal")
    h_sig.SetTitle("nominal")
    h_sig_err.SetName("error")
    h_sig_err.SetTitle("error")
    h_sig.Write()
    h_sig_err.Write()

    fout.cd("Background")
    h_bkg.SetName("nominal")
    h_bkg.SetName("nominal")
    h_bkg_err.SetName("error")
    h_bkg_err.SetTitle("error")
    h_bkg.Write()
    h_bkg_err.Write()

    fout.Close()

    # draw_rs_ws_bdt_eff()


if __name__ == "__main__":
    main(sys.argv)