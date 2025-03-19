import sys
import ROOT
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy

def main(argv):
    ROOT.gStyle.SetOptStat(0)

    bdt1 = argv[1]
    bdt2 = argv[2]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    yields = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ))

    all_species = [100, 101, 102, 110, 120, 130, 150, 151]
    species_names = ['$B^+ \to \bar{D}^0 D^0 K^+$', '$B^+ \to D^+ D^- K^+$', '$B^+ \to D^+_s D^-_s K^+$', '$B^0 \to D^- D^0 K^+$', '$B^0_s \to D^-_s D^0 K^+$', '$B^+ \to \bar{D}^0 D^+ K^0$', '$B^+ \to \bar{D}^0 D^+_s$', '$B^+ \to \bar{D}^0 D^+$']

    N_species = len(species_names)
    N_components = 4

    N = { 100: 0, 101: 0, 102: 0, 110: 0, 120: 0, 130: 0, 150: 0, 151: 0 }

    for i in range(N_species):
        for j in range(N_components):
            N[all_species[i]] += float(yields[i,j])
    
    branching_fraction = 0.001
    fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy")
    fixed_value_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs_error.npy')
    fixed_value = ufloat( fixed_value, fixed_value_err )

    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")
    N_num = t_sig.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    N_num = ufloat( N_num, np.sqrt(N_num) )
    S = (branching_fraction*N_num)/fixed_value

    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_ws = f_ws.Get("DecayTree")
    N_bkg_num = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1, bdt2))
    N_bkg_den = t_ws.GetEntries()
    eps_bkg = 2*(N_bkg_num/N_bkg_den) # factor 2 comes from eps_rs = 2*eps_ws @ BDT1=BDT2=0.8
    B = 19006409*eps_bkg

    f_BuKtautau = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    f_BuDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt1_0_bdt2_0.root")
    f_BdDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt1_0_bdt2_0.root")
    f_BsDDKp = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt1_0_bdt2_0.root")
    f_BuDDK0 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt1_0_bdt2_0.root")
    f_BuDD = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt1_0_bdt2_0.root")

    t_BuKtautau = f_BuKtautau.Get("DecayTree")
    t_BuDDKp = f_BuDDKp.Get("DecayTree")
    t_BdDDKp = f_BdDDKp.Get("DecayTree")
    t_BsDDKp = f_BsDDKp.Get("DecayTree")
    t_BuDDK0 = f_BuDDK0.Get("DecayTree")
    t_BuDD = f_BuDD.Get("DecayTree")

    nbins = 50

    h_BuKtautau = ROOT.TH1D("h_BuKtautau", "h_BuKtautau", nbins, 4000, 8000)
    h_BuD0D0Kp = ROOT.TH1D("h_BuD0D0Kp", "h_BuD0D0Kp", nbins, 4000, 8000)
    h_BuDpDmKp = ROOT.TH1D("h_BuDpDmKp", "h_BuDpDmKp", nbins, 4000, 8000)
    h_BuDsDsKp = ROOT.TH1D("h_BuDsDsKp", "h_BuDsDsKp", nbins, 4000, 8000)
    h_BdDpD0Kp = ROOT.TH1D("h_BdDpD0Kp", "h_BdDpD0Kp", nbins, 4000, 8000)
    h_BsDsD0Kp = ROOT.TH1D("h_BsDsD0Kp", "h_BsDsD0Kp", nbins, 4000, 8000)
    h_BuD0DpK0 = ROOT.TH1D("h_BuD0DpK0", "h_BuD0DpK0", nbins, 4000, 8000)
    h_BuD0Ds = ROOT.TH1D("h_BuD0Ds", "h_BuD0Ds", nbins, 4000, 8000)
    h_BuD0Dp = ROOT.TH1D("h_BuD0Dp", "h_BuD0Dp", nbins, 4000, 8000)
    h_ws = ROOT.TH1D("h_ws", "h_ws", nbins, 4000, 8000)

    t_BuKtautau.Draw("df_Bp_M >> h_BuKtautau", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    t_BuDDKp.Draw("df_Bp_M >> h_BuD0D0Kp", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 100)".format(bdt1,bdt2))
    t_BuDDKp.Draw("df_Bp_M >> h_BuDpDmKp", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 101)".format(bdt1,bdt2))
    t_BuDDKp.Draw("df_Bp_M >> h_BuDsDsKp", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 102)".format(bdt1,bdt2))
    t_BdDDKp.Draw("df_Bp_M >> h_BdDpD0Kp", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 110)".format(bdt1,bdt2))
    t_BsDDKp.Draw("df_Bp_M >> h_BsDsD0Kp", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 120)".format(bdt1,bdt2))
    t_BuDDK0.Draw("df_Bp_M >> h_BuD0DpK0", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 130)".format(bdt1,bdt2))
    t_BuDD.Draw("df_Bp_M >> h_BuD0Ds", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 150)".format(bdt1,bdt2))
    t_BuDD.Draw("df_Bp_M >> h_BuD0Dp", "(BDT1 > {0}) && (BDT2 > {1}) && (species == 151)".format(bdt1,bdt2))
    t_ws.Draw("df_Bp_M >> h_ws", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    h_BuKtautau.Scale(unumpy.nominal_values(S)/h_BuKtautau.Integral())
    h_BuD0D0Kp.Scale(N[100]/h_BuD0D0Kp.Integral())
    h_BuDpDmKp.Scale(N[101]/h_BuDpDmKp.Integral())
    h_BuDsDsKp.Scale(N[102]/h_BuDsDsKp.Integral())
    h_BdDpD0Kp.Scale(N[110]/h_BdDpD0Kp.Integral())
    h_BsDsD0Kp.Scale(N[120]/h_BsDsD0Kp.Integral())
    h_BuD0DpK0.Scale(N[130]/h_BuD0DpK0.Integral())
    h_BuD0Ds.Scale(N[150]/h_BuD0Ds.Integral())
    h_BuD0Dp.Scale(N[151]/h_BuD0Dp.Integral())
    h_ws.Scale(B/h_ws.Integral())

    if((bdt1 == 0.979) and (bdt2 == 0.974)):
        h_BuKtautau.GetXaxis().SetTitle("m_{B} (MeV)")
        h_BuKtautau.GetYaxis().SetTitle("Entries / ({0} MeV)".format(4000/nbins))
        h_BuKtautau.SetTitle("N_{{s}}^{{K #tau #tau}} = {0}".format(S))
    else:
        h_BuD0Ds.GetXaxis().SetTitle("m_{B} (MeV)")
        h_BuD0Ds.GetYaxis().SetTitle("Entries / ({0} MeV)".format(4000/nbins))
        h_BuD0Ds.SetTitle("N_{{s}}^{{K #tau #tau}} = {0}".format(S))

    h_BuKtautau.SetLineColor(1)
    h_BuD0D0Kp.SetLineColor(600)
    h_BuDpDmKp.SetLineColor(432+1)
    h_BuDsDsKp.SetLineColor(416+1)
    h_BdDpD0Kp.SetLineColor(632)
    h_BsDsD0Kp.SetLineColor(880+1)
    h_BuD0DpK0.SetLineColor(800-3)
    h_BuD0Ds.SetLineColor(900-4)
    h_BuD0Dp.SetLineColor(616+2)
    h_ws.SetLineColor(1)
    h_ws.SetLineStyle(9)

    h_BuKtautau.SetFillColorAlpha(1, 0.25)
    h_BuD0D0Kp.SetFillColorAlpha(600, 0.25)
    h_BuDpDmKp.SetFillColorAlpha(432+1, 0.25)
    h_BuDsDsKp.SetFillColorAlpha(416+1, 0.25)
    h_BdDpD0Kp.SetFillColorAlpha(632, 0.25)
    h_BsDsD0Kp.SetFillColorAlpha(880+1, 0.25)
    h_BuD0DpK0.SetFillColorAlpha(800-3, 0.25)
    h_BuD0Ds.SetFillColorAlpha(900-4, 0.25)
    h_BuD0Dp.SetFillColorAlpha(616+2, 0.25)

    leg = ROOT.TLegend(0.65,0.3,0.9,0.9)
    leg.AddEntry(h_BuKtautau, "B^{+} #rightarrow K^{+} #tau^{+} #tau^{-}")
    leg.AddEntry(h_BuD0D0Kp, "B^{+} #rightarrow #bar{D}^{0} D^{0} K^{+}")
    leg.AddEntry(h_BuDpDmKp, "B^{+} #rightarrow D^{+} D^{-} K^{+}")
    leg.AddEntry(h_BuDsDsKp, "B^{+} #rightarrow D^{+}_{s} D^{+}_{s} K^{+}")
    leg.AddEntry(h_BdDpD0Kp, "B^{0} #rightarrow D^{-} D^{0} K^{+}")
    leg.AddEntry(h_BsDsD0Kp, "B^{0}_{s} #rightarrow D^{-}_{s} D^{0} K^{+}")
    leg.AddEntry(h_BuD0DpK0, "B^{+} #rightarrow #bar{D}^{0} D^{+} K^{0}")
    leg.AddEntry(h_BuD0Ds, "B^{+} #rightarrow  #bar{D}^{0} D^{+}_{s}")
    leg.AddEntry(h_BuD0Dp, "B^{+} #rightarrow  #bar{D}^{0} D^{+}")

    c = ROOT.TCanvas("c", "c")
    c.cd()
    if((bdt1 == 0.979) and (bdt2 == 0.974)):
        # c.SetLogy()
        h_BuKtautau.Draw("HIST")   
        h_BuD0Ds.Draw("HIST same")
    else:
        h_BuD0Ds.Draw("HIST")
        h_BuKtautau.Draw("HIST same")  

    h_BuD0D0Kp.Draw("HIST same")
    h_BuDpDmKp.Draw("HIST same")
    h_BuDsDsKp.Draw("HIST same")
    h_BdDpD0Kp.Draw("HIST same")
    h_BsDsD0Kp.Draw("HIST same")
    h_BuD0DpK0.Draw("HIST same")
    h_BuD0Dp.Draw("HIST same")
    leg.Draw("same")
    c.SaveAs("/panfs/felician/B2Ktautau/workflow/yield_estimates_plots/yield_plots_bdt1_{0}_bdt2_{1}.pdf".format( bdt1, bdt2 ))

    c1 = ROOT.TCanvas('c1', 'c1')
    c1.cd()
    c1.SetLogy()

    h_ws.GetXaxis().SetTitle("m_{B} (MeV)")
    h_ws.GetYaxis().SetTitle("Entries / ({0} MeV)".format(4000/nbins))
    h_ws.SetTitle("S / #sqrt{{(S+B)}} = {0}".format( round( unumpy.nominal_values(S) /np.sqrt( unumpy.nominal_values(S) + B), 1) ))

    h_ws.Draw("HIST")
    h_BuKtautau.Draw("HIST same")   
    h_BuD0Ds.Draw("HIST same")
    h_BuD0D0Kp.Draw("HIST same")
    h_BuDpDmKp.Draw("HIST same")
    h_BuDsDsKp.Draw("HIST same")
    h_BdDpD0Kp.Draw("HIST same")
    h_BsDsD0Kp.Draw("HIST same")
    h_BuD0DpK0.Draw("HIST same")
    h_BuD0Dp.Draw("HIST same")
    leg.AddEntry(h_ws, "Comb. bkg.")
    leg.Draw("same")
    c1.SaveAs("/panfs/felician/B2Ktautau/workflow/yield_estimates_plots/yield_plots_bdt1_{0}_bdt2_{1}_comb_bkg.pdf".format( bdt1, bdt2 ))


if __name__ == "__main__":
    main(sys.argv)
