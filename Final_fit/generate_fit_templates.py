import sys
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import array
from ROOT import TMath

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
    ROOT.gStyle.SetOptFit(1)

    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]
    folder_name = argv[4]
    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    # Data
    f_data = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data = f_data.Get("h_toy_data")

    nbins = h_data.GetNbinsX()
    # bins = [h_data.GetBinLowEdge(i+1) for i in range(nbins+1)]
    # bins = array.array('d', bins)

    # Signal (w/o TM)
    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")

    # Combinatorial background
    f_bkg = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_bkg = f_bkg.Get("DecayTree")

    # Signal + background templates (nominal)
    h_sig = ROOT.TH1D("h_sig", "h_sig", nbins, 4000, 8000)
    h_bkg = ROOT.TH1D("h_bkg", "h_bkg", nbins, 4000, 8000)

    if((bdt1 >= 0.9) and (bdt2 >= 0.9)):
        t_sig.Draw("df_Bp_M >> h_sig", "(BDT1 > {0}) && (BDT2 > {1})".format(0.9,0.9))
        t_bkg.Draw("df_Bp_M >> h_bkg", "(BDT1 > {0}) && (BDT2 > {1})".format(0.9,0.9))
    else:
        t_sig.Draw("df_Bp_M >> h_sig", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
        t_bkg.Draw("df_Bp_M >> h_bkg", "(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    # nsig = t_sig.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    # nbkg = t_bkg.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    # print("nbkg = ", nbkg)

    # if(nbkg < 1000):

    #     ############################### Fit to signal template #######################################################################################
    #     def crystalball(x, par):
    #         mu = par[0]  # Mean
    #         sigma = par[1]  # Standard deviation (sigma)
    #         alpha = par[2]  # Left tail parameter (alpha)
    #         n = par[3]  # Right tail exponent (n)

    #         A = TMath.Power(n / TMath.Abs(alpha), n) * TMath.Exp(-0.5 * alpha * alpha)  # Amplitude for the left tail

    #         t = (x[0] - mu) / sigma  # Standardized variable
    #         if t < -alpha:
    #             return A*TMath.Power(n / TMath.Abs(alpha), n) * TMath.Power(n - (t + alpha) / alpha, -n)
    #         else:
    #             return A*TMath.Exp(-0.5 * t * t)

    #     def double_sided_crystalball(x, par):
    #         xx = x[0]
    #         mu     = par[0]
    #         sigma  = par[1]
    #         alphaL = abs(par[2])
    #         nL     = par[3]
    #         alphaR = abs(par[4])
    #         nR     = par[5]

    #         t = (xx - mu)/sigma

    #         fact1TLessMinosAlphaL = alphaL/nL
    #         fact2TLessMinosAlphaL = (nL/alphaL) - alphaL -t
    #         fact1THihgerAlphaH = alphaR/nR;
    #         fact2THigherAlphaH = (nR/alphaR) - alphaR +t
            
    #         if( (-alphaL <= t) and (alphaR >= t)):
    #             result = TMath.Exp(-0.5*t*t)
            
    #         elif(t < -alphaL):
    #             result = TMath.Exp(-0.5*alphaL*alphaL)*TMath.Power(fact1TLessMinosAlphaL*fact2TLessMinosAlphaL, -nL)

    #         elif(t > alphaR):
    #             result = TMath.Exp(-0.5*alphaR*alphaR)*TMath.Power(fact1THihgerAlphaH*fact2THigherAlphaH, -nR)
                
    #         return result
        
    #     def double_sided_crystalball_extended(x, par):
    #         xx = x[0]
    #         mu     = par[0]
    #         sigma  = par[1]
    #         alphaL = abs(par[2])
    #         nL     = par[3]
    #         alphaR = abs(par[4])
    #         nR     = par[5]
    #         N = par[6]

    #         p = [mu, sigma, alphaL, nL, alphaR, nR]

    #         return N*double_sided_crystalball(x, p)

    #     def gaussian(x, par):
    #         """
    #         par[0] = mean (mu)
    #         par[1] = sigma (standard deviation)
    #         par[2] = norm (amplitude)
    #         """
    #         xx = x[0]
    #         mu = par[0]
    #         sigma = par[1]

    #         # return  TMath.Exp(-0.5 * ((xx - mu) / sigma) ** 2)
    #         return TMath.Gaus(xx, mu, sigma)
        
    #     def landau(x, par):
    #         xx = x[0]
    #         mu = par[0]
    #         sigma = par[1]

    #         return TMath.Landau(xx, mu, sigma)
        
    #     def exponential(x, par):
    #         xx = x[0]
    #         lbd = par[0]
    #         xp = par[1]

    #         return TMath.Exp(-lbd*(xx - xp))
        
    #     def landau_plus_exponential_plus_gaussian(x, par):
    #         mu = par[0]
    #         sigma = par[1]
    #         sigma1 = par[2]
    #         lbd = par[3]    
    #         f = par[4]
    #         f1 = par[5]
    #         N = par[6]

    #         p1 = [mu, sigma]
    #         p2 = [lbd]
    #         p3 = [mu, sigma1]

    #         return N*( f*landau(x, p1) + f1*exponential(x, p2) + (1-f-f1)*gaussian(x, p3) )

    #     def cb_plus_exponential(x, par):
    #         xx = x[0]
    #         mu = par[0]  
    #         sigma = par[1]  
    #         alpha = par[2]  
    #         n = par[3]  
    #         lbd = par[4]
    #         xp = par[5]
    #         f = par[6]
    #         N = par[7]

    #         p = [mu, sigma, alpha, n]
    #         p1 = [lbd, xp]

    #         return N*( f*crystalball(x, p) + (1-f)*exponential(x, p1) )

    #     def landau_extended(x, par):
    #         mu = par[0]
    #         sigma = par[1]
    #         N = par[2]

    #         p = [mu, sigma]

    #         return N*landau(x, p)
        
    #     def landau_plus_landau(x, par):
    #         mu = par[0]
    #         sigma1 = par[1]
    #         sigma2 = par[2]
    #         f = par[3]
    #         N = par[4]

    #         p1 = [mu, sigma1]
    #         p2 = [mu, sigma2]

    #         return N*( f*landau(x, p1) + (1-f)*landau(x, p2) )

    #     def gaussian_plus_double_sided_CB(x, par):
    #         xx = x[0]
    #         mu     = par[0]
    #         sigma  = par[1]
    #         alphaL = abs(par[2])
    #         nL     = par[3]
    #         alphaR = abs(par[4])
    #         nR     = par[5]
    #         parameters = [mu, sigma, alphaL, nL, alphaR, nR]

    #         sigma1 = par[6]
    #         parameters1 = [mu, sigma1]

    #         f = par[7]
    #         N   = par[8]

    #         return N*(f*gaussian(x, parameters1) + (1-f)*double_sided_crystalball(x, parameters))
            
    #     def landau_plus_double_sided_CB(x, par):
    #         xx = x[0]
    #         mu     = par[0]
    #         sigma  = par[1]
    #         alphaL = abs(par[2])
    #         nL     = par[3]
    #         alphaR = abs(par[4])
    #         nR     = par[5]
    #         parameters = [mu, sigma, alphaL, nL, alphaR, nR]

    #         sigma1 = par[6]
    #         parameters1 = [mu, sigma1]

    #         f = par[7]
    #         N   = par[8]

    #         return N*(f*landau(x, parameters1) + (1-f)*double_sided_crystalball(x, parameters))
        
    #     def gaussian_plus_landau(x, par):
    #         xx = x[0]
    #         mu     = par[0]
    #         mu1     = par[1]
    #         sigma  = par[2]
    #         sigma1 = par[3]
    #         f = par[4]
    #         N = par[5]

    #         parameters  = [mu, sigma]
    #         parameters1 = [mu1, sigma1]

    #         return N*( f*gaussian(x, parameters) + (1-f)*landau(x, parameters1) )
        
    #     def CBplus_landau(x, par):
    #         xx = x[0]
    #         mu = par[0]  # Mean
    #         sigma = par[1]  # Standard deviation (sigma)
    #         alpha = par[2]  # Left tail parameter (alpha)
    #         n = par[3]  # Right tail exponent (n)
    #         sigma1 = par[4]
    #         f = par[5]
    #         N = par[6]

    #         parameters = [mu, sigma, alpha, n]
    #         parameters1 = [mu, sigma1]

    #         return N*( f*crystalball(x, parameters) + (1-f)*landau(x, parameters1) )

    #     def poly(x, par):
    #         xx = x[0]
    #         a = par[0]  
    #         b = par[1]  
    #         mu = par[2]

    #         return a*((xx-mu)**2) + b*xx 

    #     def poly_plus_landau(x, par):
    #         xx = x[0]
    #         a = par[0] 
    #         b = par[1]  
    #         mu = par[2]
    #         sigma = par[3]
    #         f = par[4]
    #         N = par[5]

    #         p = [a, b, mu]
    #         p1 = [mu, sigma]

    #         return N*(f*poly(x,p) + (1-f)*landau(x, p1))
        
    #     def poly_plus_gaussian_plus_landau(x, par):
    #         a = par[0] 
    #         b = par[1]  
    #         mu = par[2]
    #         sigma = par[3]
    #         f = par[4]

    #         sigma1 = par[5]
    #         f1 = par[6]
    #         N = par[7]

    #         p = [a, b, mu, sigma, f]
    #         p1 = [mu, sigma1]

    #         return N*( f1*gaussian(x, p1) + (1-f1)*poly_plus_landau(x, p) )
        
    #     def poly_plus_double_CB(x, par):
    #         xx = x[0]
    #         a = par[0] 
    #         b = par[1]  
    #         mu = par[2]
    #         sigma  = par[3]
    #         alphaL = abs(par[4])
    #         nL     = par[5]
    #         alphaR = abs(par[6])
    #         nR     = par[7]
    #         f = par[8]
    #         N = par[9]

    #         p1 = [a, b, mu]
    #         p2 = [mu, sigma, alphaL, nL, alphaR, nR]

    #         return N*( f*poly(x, p1) + (1-f)*double_sided_crystalball(x, p2) )

    #     # pdf_sig = ROOT.TF1("pdf_sig", "[4]*([0]*TMath::Gaus(x, [1], [2]) +(1-[0])*TMath::Landau(x, [1], [3]))", 4000, 8000)
    #     # pdf_sig = ROOT.TF1("pdf_sig", "[6]*( [0]*TMath::Gaus(x, [2], [3]) + [1]*TMath::Landau(x, [2], [4]) + (1-[0]-[1])*TMath::Landau(-x, [2], [5]) )", 4000, 8000)
    #     # pdf_sig = ROOT.TF1("pdf_sig", "[6]*([0]*TMath::Gaus(x, [2], [3]) + [1]*TMath::Gaus(x, [2], [4]) + (1-[0]-[1])*TMath::Landau(x, [2], [5]))", 4000, 8000)
        
    #     # Single crystal ball function
    #     # pdf_sig = ROOT.TF1("pdf_sig", crystalball, 4000, 8000, 4)
    #     # pdf_sig.SetParameters(5279, 200, 1, 30)
    #     # pdf_sig.SetParLimits(0, 5200, 5400)
    #     # pdf_sig.SetParLimits(1, 50, 500)
    #     # pdf_sig.SetParLimits(2, -10, 10)
    #     # pdf_sig.SetParLimits(3, 0, 1000)

    #     # Crystal ball
    #     # pdf_sig = ROOT.TF1("pdf_sig", crystalball, 4000, 8000, 4)
    #     # pdf_sig.SetParameters(5279, 200, 1, 30)
    #     # pdf_sig.SetParLimits(0, 5200, 5400)
    #     # pdf_sig.SetParLimits(1, 50, 500)
    #     # pdf_sig.SetParLimits(2, -10, 10)
    #     # pdf_sig.SetParLimits(3, 0, 1000)

    #     # Double sided crystal ball
    #     # pdf_sig = ROOT.TF1("pdf_sig", double_sided_crystalball, 4000, 8000, 7)
    #     # pdf_sig.SetParameters(5279, 200, 1, 10, 2, 30, nsig)
    #     # pdf_sig.SetParLimits(0, 5200, 5400)
    #     # pdf_sig.SetParLimits(1, 50, 500)
    #     # pdf_sig.SetParLimits(2, -10, 10)
    #     # pdf_sig.SetParLimits(3, 0, 1000)
    #     # pdf_sig.SetParLimits(4, -10, 10)
    #     # pdf_sig.SetParLimits(5, 0, 1000)
    #     # pdf_sig.SetParLimits(6, 0, 2*nsig)

    #     # Gaussian + double sided crystal ball
    #     pdf_sig = ROOT.TF1("pdf_sig", gaussian_plus_double_sided_CB, 4000, 8000, 9)
    #     pdf_sig.SetParameters(5279, 200, 1, 10, 2, 30, 200, 0.5, nsig)
    #     pdf_sig.SetParNames("mu", "sig", "alphaL", "nL", "alphaR", "nR", "sig1", "f", "nsig")
    #     pdf_sig.SetParLimits(0, 5200, 5400)
    #     pdf_sig.SetParLimits(1, 50, 1000)
    #     pdf_sig.SetParLimits(2, -10, 10)
    #     pdf_sig.FixParameter(3, 40)
    #     pdf_sig.SetParLimits(4, -10, 10)
    #     pdf_sig.FixParameter(5, 40)
    #     pdf_sig.SetParLimits(6, 50, 1000)
    #     pdf_sig.SetParLimits(7, 0, 1)
    #     pdf_sig.SetParLimits(8, 0, 2*nsig)

    #     signal_fit_status = int(h_sig.Fit("pdf_sig", "LME"))
    #     signal_pdf = h_sig.GetFunction("pdf_sig")

    #     chi2_sig = signal_pdf.GetChisquare() / signal_pdf.GetNDF()
    #     print("Signal fit status = ", signal_fit_status)
    #     print("Signal fit chi2/ndf = ", chi2_sig)

    #     c_sig = ROOT.TCanvas("c_sig", "c_sig")
    #     c_sig.cd()
    #     h_sig.GetXaxis().SetTitle("m_{B} (MeV)")
    #     h_sig.GetYaxis().SetTitle("Entries / (100 MeV)")
    #     h_sig.SetTitle("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed))
    #     rp_sig = ROOT.TRatioPlot(h_sig)
    #     rp_sig.Draw()
    #     rp_sig.GetLowerRefYaxis().SetTitle("(fit-histo)/error")
    #     rp_sig.GetUpperRefYaxis().SetTitle("Entries / (100 MeV)")
    #     rp_sig.SetH1DrawOpt("E")
    #     rp_sig.Draw("same")
    #     leg_sig = ROOT.TLegend(0.5, 0.5, 0.7, 0.7)
    #     leg_sig.AddEntry(h_sig, "Data") 
    #     leg_sig.AddEntry(signal_pdf, "Fit") 
    #     leg_sig.SetBorderSize(0)
    #     leg_sig.SetTextSize(0.04)
    #     leg_sig.SetFillStyle(0)
    #     leg_sig.Draw("same")
    #     c_sig.SaveAs("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{0}/signal_fit_bdt1_{1}_bdt2_{2}_seed_{3}.pdf".format(folder_name,bdt1,bdt2,random_seed))

    #     ##########################################################################################################################################

    #     ############################################## Fit to background template ################################################################
    #     def expGausConv(x, par):
    #         mu    = par[0]
    #         sigma = par[1]
    #         lam   = par[2]
    #         N     = par[3]

    #         arg1 = lam / 2.0 * (2 * mu + lam * sigma**2 - 2 * x[0])
    #         arg2 = (mu + lam * sigma**2 - x[0]) / (TMath.Sqrt(2) * sigma)
    #         return N * (lam / 2.0) * TMath.Exp(arg1) * TMath.Erfc(arg2)

    #     def chebyshev_pdf(x, par):
    #         # Normalize x to [-1, 1]
    #         xx = x[0]
    #         xmin = 4000
    #         xmax = 8000
    #         xnorm = 2 * (xx - xmin) / (xmax - xmin) - 1

    #         # Evaluate the Chebyshev sum
    #         n = len(par)
    #         T_prev = 1
    #         T_curr = xnorm
    #         result = par[0] * T_prev
    #         if n > 1:
    #             result += par[1] * T_curr

    #         for i in range(2, n):
    #             T_next = 2 * xnorm * T_curr - T_prev
    #             result += par[i] * T_next
    #             T_prev, T_curr = T_curr, T_next

    #         return result
                
    #     def landau_plus_chebychev(x, par):
    #         mu = par[0]
    #         sigma = par[1]
    #         a0 = par[2]
    #         a1 = par[3]
    #         a2 = par[4]
    #         a3 = par[5]
    #         f = par[6]
    #         N = par[7]

    #         p = [mu, sigma]
    #         p1 = [a0, a1, a2, a3]

    #         return N*( f*landau(x, p) + (1-f)*chebyshev_pdf(x, p1) )
        
    #     # # Landau
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", landau_extended, 4000, 8000, 3)
    #     # pdf_bkg.SetParameters(4800, 300, nbkg)
    #     # pdf_bkg.SetParNames("mu", "sig", "nbkg")
    #     # pdf_bkg.SetParLimits(0, 4000, 6000)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, 0, 2*nbkg)

    #     # # Landau + landau
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", landau_plus_landau, 4000, 8000, 5)
    #     # pdf_bkg.SetParameters(4800,  300, 300, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4500, 5500)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, 50, 500)
    #     # pdf_bkg.SetParLimits(3, 0, 1)
    #     # pdf_bkg.SetParLimits(4, 0, 2*nbkg)

    #     # # Landau + Gaussian
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", gaussian_plus_landau, 4000, 8000, 6)
    #     # pdf_bkg.SetParameters(4800, 4800, 100, 300, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4500, 5500)
    #     # pdf_bkg.SetParLimits(1, 4500, 5500)
    #     # pdf_bkg.SetParLimits(2, 50, 500)
    #     # pdf_bkg.SetParLimits(3, 50, 500)
    #     # pdf_bkg.SetParLimits(4, 0, 1)
    #     # pdf_bkg.SetParLimits(5, 0, 2*nbkg)

    #     # # Gaussian convoluted with exponential
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", expGausConv, 4000, 8000, 4)
    #     # pdf_bkg.SetParameters(4800, 100, 0.002, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4500, 5500)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, -5, 5)
    #     # pdf_bkg.SetParLimits(3, 0, 2*nbkg)

    #     # # Landau + CB
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", CBplus_landau, 4000, 8000, 7)
    #     # pdf_bkg.SetParameters(4800, 100, 5, 3, 300, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4000, 6000)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, -10, 10)
    #     # pdf_bkg.FixParameter(3, 3)
    #     # pdf_bkg.SetParLimits(4, 50, 500)
    #     # pdf_bkg.SetParLimits(5, 0, 1)
    #     # pdf_bkg.SetParLimits(6, 0, 2*nbkg)

    #     # # Landau + double sided crystal ball
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", landau_plus_double_sided_CB, 4000, 8000, 9)
    #     # pdf_bkg.SetParameters(5000, 500, 1, 10, 2, 30, 1, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4000, 6000)
    #     # pdf_bkg.SetParLimits(1, 100, 1000)
    #     # pdf_bkg.SetParLimits(2, -10, 10)
    #     # pdf_bkg.FixParameter(3, 100)
    #     # pdf_bkg.SetParLimits(4, -10, 10)
    #     # pdf_bkg.FixParameter(5, 100)
    #     # pdf_bkg.FixParameter(6, 0.01)
    #     # pdf_bkg.SetParLimits(7, 0, 1)
    #     # pdf_bkg.SetParLimits(8, 0, 2*nbkg)

    #     # # Double sided crystal ball
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", double_sided_crystalball_extended, 4000, 8000, 7)
    #     # pdf_bkg.SetParameters(5000, 500, 1, 10, 2, 30, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4000, 6000)
    #     # pdf_bkg.SetParLimits(1, 100, 1000)
    #     # pdf_bkg.SetParLimits(2, -10, 10)
    #     # pdf_bkg.FixParameter(3, 100)
    #     # pdf_bkg.SetParLimits(4, -10, 10)
    #     # pdf_bkg.FixParameter(5, 100)
    #     # pdf_bkg.SetParLimits(6, 0, 2*nbkg)

    #     # # Landau + 2nd order polynomial
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", poly_plus_landau, 4000, 8000, 6)
    #     # pdf_bkg.SetParameters(-1, 300, 4800, 400, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, -100, 10)
    #     # pdf_bkg.SetParLimits(1, -10, 1000)
    #     # pdf_bkg.SetParLimits(2, 4500, 6000)
    #     # pdf_bkg.SetParLimits(3, 100, 1000)
    #     # pdf_bkg.FixParameter(4, 0.8)
    #     # pdf_bkg.SetParLimits(5, 0, 2*nbkg)

    #     # # CB + exponential
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", cb_plus_exponential, 4000, 8000, 8)
    #     # pdf_bkg.SetParameters(4800, 200, 1.5, 3, 0.01, 5500, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4000, 6000)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, -10, 10)
    #     # pdf_bkg.FixParameter(3, 3)
    #     # pdf_bkg.SetParLimits(4,  0.001, 0.1)
    #     # pdf_bkg.SetParLimits(5,  4000, 8000)
    #     # pdf_bkg.SetParLimits(6, 0, 1)
    #     # pdf_bkg.SetParLimits(7, 0, 2*nbkg)

    #     # def double_gaussian_plus_landau(x, par):
    #     #     mu     = par[0]
    #     #     sigma  = par[1]
    #     #     sigma1 = par[2]
    #     #     sigma2 = par[3]
    #     #     f      = par[4]
    #     #     f1     = par[5]
    #     #     N      = par[6]

    #     #     p1 = [mu, sigma]
    #     #     p2 = [mu, sigma1]
    #     #     p3 = [mu, sigma2]

    #     #     return N*( f1*(f*gaussian(x, p1) + (1-f)*gaussian(x, p2)) + (1-f1)*landau(x, p3)   )

    #     # pdf_bkg = ROOT.TF1("pdf_bkg", double_gaussian_plus_landau, 4000, 8000, 7)
    #     # pdf_bkg.SetParameters(4800, 100, 300, 300, 0.25, 0.25, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4500, 5500)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, 50, 500)
    #     # pdf_bkg.SetParLimits(3, 50, 500)
    #     # pdf_bkg.SetParLimits(4, 0, 1)
    #     # pdf_bkg.SetParLimits(5, 0, 1)
    #     # pdf_bkg.SetParLimits(6, 0, 2*nbkg)

    #     # pdf_bkg = ROOT.TF1("pdf_bkg", landau_plus_exponential_plus_gaussian, 4000, 8000, 7)
    #     # pdf_bkg.SetParameters(4800, 300, 300, 0.01, 0.5, 0.5, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4500, 5500)
    #     # pdf_bkg.SetParLimits(1, 50, 500)
    #     # pdf_bkg.SetParLimits(2, 50, 500)
    #     # pdf_bkg.SetParLimits(3, 0.001, 0.1)
    #     # pdf_bkg.SetParLimits(4, 0, 1)
    #     # pdf_bkg.SetParLimits(5, 0, 1)
    #     # pdf_bkg.SetParLimits(6, 0, 2*nbkg)

    #     # # Convolution of gaussian and landau
    #     # pdf_landau = ROOT.TF1("pdf_landau", landau, 4000, 8000, 2)
    #     # pdf_gaussian = ROOT.TF1("pdf_gaussian", gaussian, 4000, 8000, 2)
    #     # convolution = ROOT.TF1Convolution("landau", "gaus", 4000, 8000)
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", convolution, 4000, 8000, convolution.GetNpar())
    #     # print(convolution.GetNpar())
    #     # pdf_bkg.SetParameters(4800, 200, 100, nbkg/2, nbkg/2)
    #     # # pdf_bkg.SetParLimits(0, 4000, 60000)
    #     # # pdf_bkg.SetParLimits(1, 50, 5000)
    #     # # pdf_bkg.SetParLimits(2, 4000, 60000)
    #     # # pdf_bkg.SetParLimits(3, 50, 500)
    #     # # pdf_bkg.SetParLimits(4, 0, 2*nbkg)
        
    #     # # Convolution of gaussian and exponential
    #     # def emg(x, par):
    #     #     xx     = x[0]
    #     #     mu     = par[0]
    #     #     sigma  = par[1]
    #     #     lambd  = par[2]
    #     #     N = par[3]

    #     #     if sigma <= 0 or lambd <= 0:
    #     #         return 0.0  # Avoid division by zero or invalid domain

    #     #     # Argument for erfc
    #     #     arg = (mu + lambd * sigma**2 - xx) / (TMath.Sqrt(2) * sigma)

    #     #     # Exponential part
    #     #     exp_part = TMath.Exp(0.5 * lambd * (2 * mu + lambd * sigma**2 - 2 * xx))

    #     #     # Error function part
    #     #     erfc_part = TMath.Erfc(arg)

    #     #     # Final EMG formula
    #     #     return N*(exp_part*erfc_part)
        
    #     # pdf_bkg = ROOT.TF1("pdf_bkg", emg, 4000, 8000, 4)
    #     # pdf_bkg.SetParameters(4800, 400, 0.005, nbkg)
    #     # pdf_bkg.SetParLimits(0, 4500, 5200)
    #     # pdf_bkg.SetParLimits(1, 50, 1000)
    #     # pdf_bkg.SetParLimits(2, -0.1, 0.1)
    #     # pdf_bkg.SetParLimits(3, 0, 10*nbkg)

    #     # Landau + chebychev
    #     pdf_bkg = ROOT.TF1("pdf_bkg", landau_plus_chebychev, 4000, 8000, 8)
    #     pdf_bkg.SetParameters(4800, 400, 0.5, 0.5, 0.5, 0.05, 0.5, nbkg)
    #     pdf_bkg.SetParNames("mu", "sig", "a0", "a1", "a2", "a3", "f", "N")
    #     pdf_bkg.SetParLimits(0, 4500, 5500)
    #     pdf_bkg.SetParLimits(1, 50, 1000)
    #     pdf_bkg.SetParLimits(2, -1, 10)
    #     pdf_bkg.SetParLimits(3, -5, 5)
    #     pdf_bkg.SetParLimits(3, -5, 5)
    #     pdf_bkg.SetParLimits(4, -2, 2)
    #     pdf_bkg.SetParLimits(5, -1, 1)
    #     pdf_bkg.FixParameter(6, 0.95)
    #     if((bdt1 >= 0.98) and (bdt2 >= 0.98)):
    #         pdf_bkg.FixParameter(6, 1)
    #     pdf_bkg.SetParLimits(7, 0, 2*nbkg)
        
    #     background_fit_status = int(h_bkg.Fit("pdf_bkg", "LME"))
    #     background_pdf = h_bkg.GetFunction("pdf_bkg")

    #     chi2_bkg = background_pdf.GetChisquare() / background_pdf.GetNDF()
    #     print("Background fit status = ", background_fit_status)
    #     print("Background fit chi2/ndf = ", chi2_bkg)

    #     c_bkg = ROOT.TCanvas("c_bkg", "c_bkg")
    #     c_bkg.cd()
    #     h_bkg.GetXaxis().SetTitle("m_{B} (MeV)")
    #     h_bkg.GetYaxis().SetTitle("Entries / (100 MeV)")
    #     h_bkg.SetTitle("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed))
    #     rp_bkg = ROOT.TRatioPlot(h_bkg)
    #     rp_bkg.Draw()
    #     rp_bkg.GetLowerRefYaxis().SetTitle("(fit-histo)/error")
    #     rp_bkg.GetUpperRefYaxis().SetTitle("Entries / (100 MeV)")
    #     rp_bkg.SetH1DrawOpt("E")
    #     rp_bkg.Draw("same")
    #     leg_bkg = ROOT.TLegend(0.5, 0.5, 0.7, 0.7)
    #     leg_bkg.AddEntry(h_bkg, "Data")
    #     leg_bkg.AddEntry(background_pdf, "Fit")
    #     leg_bkg.SetBorderSize(0)
    #     leg_bkg.SetTextSize(0.04)
    #     leg_bkg.SetFillStyle(0)
    #     leg_bkg.Draw("same")
    #     c_bkg.SaveAs("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{0}/background_fit_bdt1_{1}_bdt2_{2}_seed_{3}.pdf".format(folder_name,bdt1,bdt2,random_seed))

    #     ######################################################################################################################################################

    #     sig_fit_values = np.zeros(nbins)
    #     bkg_fit_values = np.zeros(nbins)

    #     for i in range(nbins):
    #         sig_fit_values[i] = signal_pdf.Eval( h_sig.GetBinCenter(i+1) )
    #         bkg_fit_values[i] = background_pdf.Eval( h_bkg.GetBinCenter(i+1) )

    #     h_sig_fit = ROOT.TH1D("h_sig_fit", "h_sig_fit", nbins, 4000, 8000)
    #     h_bkg_fit = ROOT.TH1D("h_bkg_fit", "h_bkg_fit", nbins, 4000, 8000)

    #     for i in range(nbins):
    #         h_sig_fit.Fill( h_sig.GetBinCenter(i+1), sig_fit_values[i] )
    #         h_bkg_fit.Fill( h_bkg.GetBinCenter(i+1), bkg_fit_values[i] )

    #         h_sig_fit.SetBinError(i+1, np.sqrt( np.abs(sig_fit_values[i]) ))
    #         h_bkg_fit.SetBinError(i+1, np.sqrt( np.abs(bkg_fit_values[i]) ))

    #     if((h_sig_fit.Integral() != 0) and (h_bkg_fit.Integral() != 0)):
    #         h_sig_fit.Scale(1/h_sig_fit.Integral())
    #         h_bkg_fit.Scale(1/h_bkg_fit.Integral())

    #     h_sig_err = ROOT.TH1D("h_sig_err", "h_sig_err", nbins, 4000, 8000)
    #     h_bkg_err = ROOT.TH1D("h_bkg_err", "h_bkg_err", nbins, 4000, 8000)

    #     for i in range(nbins):
    #         h_sig_err.Fill(h_sig_fit.GetBinCenter(i+1), h_sig_fit.GetBinError(i+1) )
    #         h_sig_err.SetBinError(i+1, 0)

    #         h_bkg_err.Fill(h_bkg_err.GetBinCenter(i+1), h_bkg_fit.GetBinError(i+1) )
    #         h_bkg_err.SetBinError(i+1, 0)

    # N_tot = h_sig.GetEntries()
    # print("N_tot = ", N_tot)

    # else:
    # Scale histograms to unit area
    if((h_sig.Integral() != 0) and (h_bkg.Integral() != 0)):
        h_sig.Scale(1/h_sig.Integral())
        h_bkg.Scale(1/h_bkg.Integral())

    # Signal + background uncertainties (error)
    h_sig_err = ROOT.TH1D("h_sig_err", "h_sig_err", nbins, 4000, 8000)
    h_bkg_err = ROOT.TH1D("h_bkg_err", "h_bkg_err", nbins, 4000, 8000)

    for i in range(nbins):
        # if( h_sig.GetBinContent(i+1) == 0 ):
        #     h_sig.SetBinError(i+1, 10**(-6))
        # if( h_bkg.GetBinContent(i+1) == 0 ):
        #     h_bkg.SetBinError(i+1, 10**(-6))

        # h_sig_err.Fill(h_sig.GetBinCenter(i+1), h_sig.GetBinError(i+1)) # staterror (absolute)
        if(h_sig.GetBinContent(i+1) == 0):
            h_sig_err.Fill(h_sig.GetBinCenter(i+1), 0.0 ) # shapesys (relative)
        else:
            if(h_sig.GetBinError(i+1)/h_sig.GetBinContent(i+1) < 0.05):
                h_sig_err.Fill(h_sig.GetBinCenter(i+1), 0.0 ) # shapesys (relative)
            else:
                h_sig_err.Fill(h_sig.GetBinCenter(i+1), h_sig.GetBinError(i+1)/h_sig.GetBinContent(i+1)) # shapesys (relative)
        h_sig_err.SetBinError(i+1, 0)

        # h_bkg_err.Fill(h_bkg.GetBinCenter(i+1), h_bkg.GetBinError(i+1)) # staterror (absolute)
        if(h_bkg.GetBinContent(i+1) == 0):
            h_bkg_err.Fill(h_bkg.GetBinCenter(i+1), 0.0) # shapesys (relative)
        else:
            if(h_bkg.GetBinError(i+1)/h_bkg.GetBinContent(i+1) < 0.05):
                h_bkg_err.Fill(h_bkg.GetBinCenter(i+1), 0.0) # shapesys (relative)
            else:
                h_bkg_err.Fill(h_bkg.GetBinCenter(i+1), h_bkg.GetBinError(i+1)/h_bkg.GetBinContent(i+1)) # shapesys (relative)
        h_bkg_err.SetBinError(i+1, 0)    

        # print("SIGNAL: ", h_sig.GetBinContent(i+1), " bin error = ", h_sig.GetBinError(i+1), " sqrt(N/Ntot) = ", np.sqrt(h_sig.GetBinContent(i+1)), "  sqrt(N)/Ntot = ", np.sqrt(h_sig_raw.GetBinContent(i+1))/N_tot )

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{3}/fit_templates_bdt1_{0}_bdt2_{1}_seed_{2}.root".format(bdt1, bdt2, random_seed, folder_name), "RECREATE")
    fout.mkdir("Data")
    fout.mkdir("Signal")
    fout.mkdir("Background")

    fout.cd("Data")
    h_data.SetTitle("data")
    h_data.SetName("data")
    h_data.Write()

    # if(nbkg < 1000):
    #     fout.cd("Signal")
    #     h_sig_fit.SetName("nominal")
    #     h_sig_fit.SetTitle("nominal")
    #     h_sig_err.SetName("error")
    #     h_sig_err.SetTitle("error")
    #     h_sig_fit.Write()
    #     h_sig_err.Write()

    #     fout.cd("Background")
    #     h_bkg_fit.SetName("nominal")
    #     h_bkg_fit.SetTitle("nominal")
    #     h_bkg_err.SetName("error")
    #     h_bkg_err.SetTitle("error")
    #     h_bkg_fit.Write()
    #     h_bkg_err.Write()
    # else:

    # DEFAULT
    fout.cd("Signal")
    h_sig.SetName("nominal")
    h_sig.SetTitle("nominal")
    h_sig_err.SetName("error")
    h_sig_err.SetTitle("error")
    h_sig.Write()
    h_sig_err.Write()

    fout.cd("Background")
    h_bkg.SetName("nominal")
    h_bkg.SetTitle("nominal")
    h_bkg_err.SetName("error")
    h_bkg_err.SetTitle("error")
    h_bkg.Write()
    h_bkg_err.Write()

    fout.Close()

    # draw_rs_ws_bdt_eff()

if __name__ == "__main__":
    main(sys.argv)
