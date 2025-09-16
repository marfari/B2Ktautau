import sys
import ROOT
import numpy as np
import time
import pyhf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
from uncertainties import ufloat
from pyhf.contrib.viz import brazil
from scipy.interpolate import griddata
import re
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import math
from itertools import chain
import matplotlib.gridspec as gridspec
from scipy.stats import chi2

# Global flags
add_BDDX_decays = False
add_statistical_error = False
toy_based_limit = False
comb_shape_sys = 0
validate_fit = False
n_toys = 1000

def save_dummy(fit_type, channel_type, BF_sig, bdt, comb_yield_syst=0):
    print("Fit failed. Saving dummy results.")
    if(validate_fit):
        f, ax = plt.subplots()
        f.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_validation_plot.pdf')
        f1, ax1 = plt.subplots()
        f1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/pulls/BF_sig.pdf')
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_poi_mean.npy', 0)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_poi_mean_error.npy', 0)  
    else:
        if(comb_yield_syst == 0):
            f, ax = plt.subplots()
            f.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot.pdf')  

            f3, ax3 = plt.subplots()
            f3.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit.pdf')  
        elif(comb_yield_syst == 1):
            f1, ax1 = plt.subplots()
            f1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot_up.pdf')  

            f4, ax4 = plt.subplots()
            f4.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit_up.pdf')  
        elif(comb_yield_syst == -1):
            f2, ax2 = plt.subplots()
            f2.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot_down.pdf')  

            f5, ax5 = plt.subplots()
            f5.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit_down.pdf')  


def retrieve_histograms(fit_type, f, f1, ch):
    # Fit templates (signal, combinatorial, BDDK+, upward, downward, data)
    if((fit_type == "ToyData") or (fit_type == "RSData")):
        h_sig = f.Get(f"Channel_{ch}/h_sig_{ch}")
        h_comb = f.Get(f"Channel_{ch}/h_comb_{ch}")
        h_BDDKp = f.Get(f"Channel_{ch}/h_BDDKp_{ch}")
        h_BuDDK0 = f.Get(f"Channel_{ch}/h_BuDDK0_{ch}")
        h_BuDD = f.Get(f"Channel_{ch}/h_BuDD_{ch}")
        h_upward = f.Get(f"Channel_{ch}/h_upward_{ch}") 
        h_downward = f.Get(f"Channel_{ch}/h_downward_{ch}")
        h_data = f.Get(f"Channel_{ch}/h_data_{ch}")

        h_sig_err = f1.Get(f"Channel_{ch}/h_sig_err_{ch}")
        h_comb_err = f1.Get(f"Channel_{ch}/h_comb_err_{ch}")
        h_BDDKp_err = f1.Get(f"Channel_{ch}/h_BDDKp_err_{ch}")
        h_BuDDK0_err = f1.Get(f"Channel_{ch}/h_BuDDK0_err_{ch}")
        h_BuDD_err = f1.Get(f"Channel_{ch}/h_BuDD_err_{ch}")

    elif((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
        h_sig = f.Get(f"Channel_{ch}/h_sig_sideband_{ch}")
        h_comb = f.Get(f"Channel_{ch}/h_comb_sideband_{ch}")
        h_BDDKp = f.Get(f"Channel_{ch}/h_BDDKp_sideband_{ch}")
        h_BuDDK0 = f.Get(f"Channel_{ch}/h_BuDDK0_sideband_{ch}")
        h_BuDD = f.Get(f"Channel_{ch}/h_BuDD_sideband_{ch}")
        h_upward = f.Get(f"Channel_{ch}/h_upward_sideband_{ch}") 
        h_downward = f.Get(f"Channel_{ch}/h_downward_sideband_{ch}")
        h_data = f.Get(f"Channel_{ch}/h_data_sideband_{ch}")

        h_sig_err = f1.Get(f"Channel_{ch}/h_sig_sideband_err_{ch}")
        h_comb_err = f1.Get(f"Channel_{ch}/h_comb_sideband_err_{ch}")
        h_BDDKp_err = f1.Get(f"Channel_{ch}/h_BDDKp_sideband_err_{ch}")
        h_BuDDK0_err = f1.Get(f"Channel_{ch}/h_BuDDK0_sideband_err_{ch}")
        h_BuDD_err = f1.Get(f"Channel_{ch}/h_BuDD_sideband_err_{ch}")

    histograms = [h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD, h_upward, h_downward, h_data]
    histograms_errors = [h_sig_err, h_comb_err, h_BDDKp_err, h_BuDDK0_err, h_BuDD_err]

    return histograms, histograms_errors


def retrieve_normalisations(fit_type, bdt):
    # A, B, C, N_comb
    B = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B_err.npy')

    if((fit_type == "ToyData") or (fit_type == "RSData")):
        A = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A.npy')
        C_BDDKp = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BDDKp.npy')
        C_BuDDK0 = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDDK0.npy')
        C_BuDD = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDD.npy')
        N_comb = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb.npy')

        A_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_err.npy')
        C_BDDKp_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BDDKp_err.npy')
        C_BuDDK0_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDDK0_err.npy')
        C_BuDD_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDD_err.npy')
        N_comb_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_err.npy')
    
    elif((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
        A = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_sideband.npy')
        C_BDDKp = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BDDKp_sideband.npy')
        C_BuDDK0 = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDDK0_sideband.npy')
        C_BuDD = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDD_sideband.npy')
        N_comb = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_sideband.npy')

        A_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/A_sideband_err.npy')
        C_BDDKp_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BDDKp_sideband_err.npy')
        C_BuDDK0_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDDK0_sideband_err.npy')
        C_BuDD_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/C_BuDD_sideband_err.npy')
        N_comb_err = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/N_comb_sideband_err.npy')

    norm_values = [A, B, C_BDDKp, C_BuDDK0, C_BuDD, N_comb]
    norm_errors = [A_err, B_err, C_BDDKp_err, C_BuDDK0_err, C_BuDD_err, N_comb_err]
    return norm_values, norm_errors


def generate_nominal_data(fit_type, channel_type, BF_sig, histograms, norm_values, norm_errors, comb_yield_syst=0):
    if((fit_type == "RSData") or (fit_type == "RSDataSidebands")):
        print("Wrong fit type for generating CLs data")
        quit()

    A = norm_values[0]
    B = norm_values[1]
    C_BDDKp = norm_values[2]
    C_BuDDK0 = norm_values[3]
    C_BuDD = norm_values[4]
    N_comb = norm_values[5]
    N_comb_err = norm_errors[5]

    if(channel_type == "AllEvents"):
        n_channels = 1
    else:
        n_channels = 3
    
    histo_data = []
    for i in range(n_channels):
        if(n_channels == 1):
            ch = i
        else:
            ch = i+1

        h_sig = histograms[i][0]
        h_comb = histograms[i][1]
        h_BDDKp = histograms[i][2]
        h_BuDDK0 = histograms[i][3]
        h_BuDD = histograms[i][4]
        nbins = h_sig.GetNbinsX()

        N_sig = BF_sig*A[ch]*B
        N_BDDKp = C_BDDKp[ch]*B
        N_BuDDK0 = C_BuDDK0[ch]*B
        N_BuDD = C_BuDD[ch]*B

        if(comb_yield_syst == 0):
            n_comb = N_comb[ch]
        elif(comb_yield_syst == 1):
            n_comb = N_comb[ch]+N_comb_err[ch]
        elif(comb_yield_syst == -1):
            n_comb = N_comb[ch]-N_comb_err[ch]
        else:
            print("Wrong value for comb_yield_syst. Should be -1, 0 or 1.")
            quit()

        # data
        h_data_sig = h_sig.Clone(f"h_data_sig_{ch}")
        h_data_sig.SetTitle(f"h_data_sig_{ch}")
        h_data_sig.Scale(N_sig/h_data_sig.Integral())

        h_data_comb = h_comb.Clone(f"h_data_comb_{ch}")
        h_data_comb.SetTitle(f"h_data_comb_{ch}")
        h_data_comb.Scale(n_comb/h_data_comb.Integral())

        h_data_BDDKp = h_BDDKp.Clone(f"h_data_BDDKp_{ch}")
        h_data_BDDKp.SetTitle(f"h_data_BDDKp_{ch}")
        h_data_BDDKp.Scale(N_BDDKp/h_data_BDDKp.Integral())

        h_data_BuDDK0 = h_BuDDK0.Clone(f"h_data_BuDDK0_{ch}")
        h_data_BuDDK0.SetTitle(f"h_data_BuDDK0_{ch}")
        h_data_BuDDK0.Scale(N_BuDDK0/h_data_BuDDK0.Integral())

        h_data_BuDD = h_BuDD.Clone(f"h_data_BuDD_{ch}")
        h_data_BuDD.SetTitle(f"h_data_BuDD_{ch}")
        h_data_BuDD.Scale(N_BuDD/h_data_BuDD.Integral())

        h_data = ROOT.TH1D(f"h_data_{ch}", f"h_data_{ch}", nbins, 4000, 9000)
        h_data.Sumw2()
        h_data.Add(h_data_sig)
        h_data.Add(h_data_comb)
        h_data.Add(h_data_BDDKp)
        if(add_BDDX_decays):
            h_data.Add(h_data_BuDDK0)
            h_data.Add(h_data_BuDD)

        histo_data.append(h_data)
    
    return histo_data


def build_model(fit_type, channel_type, BF_sig, h_data, histograms, histograms_errors, norm_values, norm_errors, i_min, i_max):
    # Normalisations
    A = norm_values[0]
    B = norm_values[1]
    C_BDDKp = norm_values[2]
    C_BuDDK0 = norm_values[3]
    C_BuDD = norm_values[4]
    N_comb = norm_values[5]

    A_err = norm_errors[0]
    B_err = norm_errors[1]
    C_BDDKp_err = norm_errors[2]
    C_BuDDK0_err = norm_errors[3]
    C_BuDD_err = norm_errors[4]
    N_comb_err = norm_errors[5]

    if((B == 0) or (B_err > B)):
        upward_B = 1.999
        downward_B = 0.001
    else:
        upward_B = 1+(B_err/B)
        downward_B = 1-(B_err/B)

    if((channel_type == "AllEvents")):
        n_channels = 1
    else:
        n_channels = 3
    
    # Build samples and parameters
    channel_samples = []
    parameters = [{"name": "lumi", "auxdata": [1.0], "bounds": [[0.5,1.5]], "fixed": True, "inits": [1.0], "sigmas": [0.1]}, {"name": "BF_sig", "bounds": [[-1,1]], "inits": [BF_sig]}]
    
    for i in range(n_channels):
        if(n_channels == 1):
            ch = i
        else:
            ch = i+1

        # NormFactor and NormSys modifiers
        a = A[ch]
        c_BDDKp = C_BDDKp[ch]
        c_BuDDK0 = C_BuDDK0[ch]
        c_BuDD = C_BuDD[ch]
        n_comb = N_comb[ch]

        a_err = A_err[ch]
        c_BDDKp_err = C_BDDKp_err[ch]
        c_BuDDK0_err = C_BuDDK0_err[ch]
        c_BuDD_err = C_BuDD_err[ch]
        n_comb_err = N_comb_err[ch]

        if((a == 0) or (a_err > a)):
            upward_A = 1.999
            downward_A = 0.001
        else:
            upward_A = 1+(a_err/a)
            downward_A = 1-(a_err/a)

        if((c_BDDKp == 0) or (c_BDDKp_err > c_BDDKp)):
            upward_C_BDDKp = 1.999
            downward_C_BDDKp = 0.001
        else:
            upward_C_BDDKp = 1+(c_BDDKp_err/c_BDDKp)
            downward_C_BDDKp = 1-(c_BDDKp_err/c_BDDKp)

        if((c_BuDDK0 == 0) or (c_BuDDK0_err > c_BuDDK0)):
            upward_C_BuDDK0 = 1.999
            downward_C_BuDDK0 = 0.001
        else:
            upward_C_BuDDK0 = 1+(c_BuDDK0_err/c_BuDDK0)
            downward_C_BuDDK0 = 1-(c_BuDDK0_err/c_BuDDK0)

        if((c_BuDD == 0) or (c_BuDD_err > c_BuDD)):
            upward_C_BuDD = 1.999
            downward_C_BuDD = 0.001
        else:
            upward_C_BuDD = 1+(c_BuDD_err/c_BuDD)
            downward_C_BuDD = 1-(c_BuDD_err/c_BuDD)

        sig_modifiers = [{"data": None, "name": "BF_sig", "type": "normfactor"}, {"data": {"hi": upward_A, "lo": downward_A}, "name": f"alpha_A_{ch}", "type": "normsys"}, {"data": {"hi": upward_B, "lo": downward_B}, "name": "alpha_B", "type": "normsys"}]
        comb_modifiers = [{"data": None, "name": f"N_comb_{ch}", "type": "normfactor"}]   
        BDDKp_modifiers = [{"data": {"hi": upward_B, "lo": downward_B}, "name": "alpha_B", "type": "normsys"}, {"data": {"hi": upward_C_BDDKp, "lo": downward_C_BDDKp}, "name": f"alpha_C_BDDKp_{ch}", "type": "normsys"}]
        BuDDK0_modifiers = [{"data": {"hi": upward_B, "lo": downward_B}, "name": "alpha_B", "type": "normsys"}, {"data": {"hi": upward_C_BuDDK0, "lo": downward_C_BuDDK0}, "name": f"alpha_C_BuDDK0_{ch}", "type": "normsys"}]
        BuDD_modifiers = [{"data": {"hi": upward_B, "lo": downward_B}, "name": "alpha_B", "type": "normsys"}, {"data": {"hi": upward_C_BuDD, "lo": downward_C_BuDD}, "name": f"alpha_C_BuDD_{ch}", "type": "normsys"}]

        # ShapeSys and HistoSys modifiers
        h_sig = histograms[i][0]
        h_comb = histograms[i][1]
        h_BDDKp = histograms[i][2]
        h_BuDDK0 = histograms[i][3]
        h_BuDD = histograms[i][4]
        h_upward = histograms[i][5]
        h_downward = histograms[i][6]
        h_sig_err = histograms_errors[i][0]
        h_comb_err = histograms_errors[i][1]
        h_BDDKp_err = histograms_errors[i][2]
        h_BuDDK0_err = histograms_errors[i][3]
        h_BuDD_err = histograms_errors[i][4]

        nbins = h_sig.GetNbinsX()

        data_sig = [h_sig.GetBinContent(k+1) for k in range(nbins)]
        data_comb = [h_comb.GetBinContent(k+1) for k in range(nbins)]
        data_BDDKp = [h_BDDKp.GetBinContent(k+1) for k in range(nbins)]
        data_BuDDK0 = [h_BuDDK0.GetBinContent(k+1) for k in range(nbins)]
        data_BuDD = [h_BuDD.GetBinContent(k+1) for k in range(nbins)]
        data_upward = [h_upward.GetBinContent(k+1) for k in range(nbins)]
        data_downward = [h_downward.GetBinContent(k+1) for k in range(nbins)] 
        data_sig_err = [h_sig_err.GetBinContent(k+1) for k in range(nbins)]
        data_comb_err = [h_comb_err.GetBinContent(k+1) for k in range(nbins)]
        data_BDDKp_err = [h_BDDKp_err.GetBinContent(k+1) for k in range(nbins)]
        data_BuDDK0_err = [h_BuDDK0_err.GetBinContent(k+1) for k in range(nbins)]
        data_BuDD_err = [h_BuDD_err.GetBinContent(k+1) for k in range(nbins)]

        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            data_sig = [data_sig[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_sig)) )]
            data_comb = [data_comb[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_comb)) )]
            data_BDDKp = [data_BDDKp[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BDDKp)))]
            data_BuDDK0 = [data_BuDDK0[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BuDDK0)))]
            data_BuDD = [data_BuDD[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BuDD)))]
            data_upward = [data_upward[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_upward)))]
            data_downward = [data_downward[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_downward)))]
            data_sig_err = [data_sig_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_sig_err)) )]
            data_comb_err = [data_comb_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_comb_err)) )]
            data_BDDKp_err = [data_BDDKp_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BDDKp_err)))]
            data_BuDDK0_err = [data_BuDDK0_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BuDDK0_err)))]
            data_BuDD_err = [data_BuDD_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BuDD_err)))]

        if(add_statistical_error):
            sig_modifiers.append({"name": f"gamma_sig_{ch}", "data": data_sig_err, "type": "shapesys"})
            comb_modifiers.append({"name": f"gamma_comb_{ch}", "data": data_comb_err, "type": "shapesys"})
            BDDKp_modifiers.append({"name": f"gamma_BDDKp_{ch}", "data": data_BDDKp_err, "type": "shapesys"})
            BuDDK0_modifiers.append({"name": f"gamma_BuDDK0_{ch}", "data": data_BuDDK0_err, "type": "shapesys"})
            BuDD_modifiers.append({"name": f"gamma_BuDD_{ch}", "data": data_BuDD_err, "type": "shapesys"})

        if(comb_shape_sys == 0):
            comb_modifiers.append({"name": f"beta_{ch}", "data": {"hi_data": data_upward, "lo_data": data_downward}, "type": "histosys"})

        samples = [{"name": "Signal", "data": data_sig, "modifiers": sig_modifiers}]
        if(comb_shape_sys == 1):
            samples.append({"name": "Combinatorial", "data": data_upward, "modifiers": comb_modifiers})
        elif(comb_shape_sys == -1):
            samples.append({"name": "Combinatorial", "data": data_downward, "modifiers": comb_modifiers})
        else:
            samples.append({"name": "Combinatorial", "data": data_comb, "modifiers": comb_modifiers})
        samples.append({"name": "BDDKp", "data": data_BDDKp, "modifiers": BDDKp_modifiers})   
        if(add_BDDX_decays):
            samples.append({"name": "BuDDK0", "data": data_BuDDK0, "modifiers": BuDDK0_modifiers})   
            samples.append({"name": "BuDD", "data": data_BuDD, "modifiers": BuDD_modifiers})   
        channel_samples.append(samples)

        parameters.append({"name": f"N_comb_{ch}", "bounds": [[-2*(np.abs(N_comb[ch])+1), 2*(np.abs(N_comb[ch])+1)]], "inits": [N_comb[ch]]})

    # Build spec
    if(channel_type == "AllEvents"):
        nbins = h_data[0].GetNbinsX()

        data_values = [h_data[0].GetBinContent(i+1) for i in range(nbins)]
        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            data_values = [data_values[k] for k in chain(range(0, i_min[0]-1), range(i_max[0], nbins) )]

        channels = [{"name": "AllEvents", "samples": channel_samples[0]}]
        observations = [{"name": "AllEvents", "data": data_values}]        
    else:
        nbins_1 = h_data[0].GetNbinsX()
        nbins_2 = h_data[1].GetNbinsX()
        nbins_3 = h_data[2].GetNbinsX()

        data_values_1 = [h_data[0].GetBinContent(k+1) for k in range(nbins_1)]
        data_values_2 = [h_data[1].GetBinContent(k+1) for k in range(nbins_2)]
        data_values_3 = [h_data[2].GetBinContent(k+1) for k in range(nbins_3)]

        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            data_values_1 = [data_values_1[k] for k in chain(range(0, i_min[1]-1), range(i_max[1], nbins_1) )]
            data_values_2 = [data_values_2[k] for k in chain(range(0, i_min[2]-1), range(i_max[2], nbins_2) )]
            data_values_3 = [data_values_3[k] for k in chain(range(0, i_min[3]-1), range(i_max[3], nbins_3) )]

        channels = [{"name": "Channel_1", "samples": channel_samples[0]}, {"name": "Channel_2", "samples": channel_samples[1]}, {"name": "Channel_3", "samples": channel_samples[2]}]
        observations = [{"name": "Channel_1", "data": data_values_1}, {"name": "Channel_2", "data": data_values_2}, {"name": "Channel_3", "data": data_values_3}]

    spec = {
        "channels": channels,
        "measurements": [
            {
                "config": {
                    "parameters": parameters,
                    "poi": "BF_sig"
                },
                "name": "MinimalConfiguration"
            }
        ],
        "observations": observations,
        "version": "1.0.0"}

    return spec


def plot(fit_type, channel_type, BF_sig, bdt, spec, h_data, fit_pars, fit_errors, i_min, i_max, comb_yield_syst=0, log_scale=True, seed=-1):
    workspace = pyhf.Workspace(spec)
    model = workspace.model()

    if(channel_type == "AllEvents"):
        nbins = h_data[0].GetNbinsX()
        mass = [h_data[0].GetBinCenter(i+1) for i in range(nbins)] 
        edges = [h_data[0].GetBinLowEdge(i+1) for i in range(nbins+1)] 
    else:
        nbins_1 = h_data[0].GetNbinsX()
        nbins_2 = h_data[1].GetNbinsX()
        nbins_3 = h_data[2].GetNbinsX()

        mass_1 = [h_data[0].GetBinCenter(i+1) for i in range(nbins_1)] 
        mass_2 = [h_data[1].GetBinCenter(i+1) for i in range(nbins_2)] 
        mass_3 = [h_data[2].GetBinCenter(i+1) for i in range(nbins_3)] 

        edges_1 = [h_data[0].GetBinLowEdge(i+1) for i in range(nbins_1+1)] 
        edges_2 = [h_data[1].GetBinLowEdge(i+1) for i in range(nbins_2+1)] 
        edges_3 = [h_data[2].GetBinLowEdge(i+1) for i in range(nbins_3+1)] 

        nbins = [nbins_1, nbins_2, nbins_3]
        mass = [mass_1, mass_2, mass_3]
        edges = [edges_1, edges_2, edges_3]

    # Retrieve fit parameters     
    par_names = model.config.par_names
    table_data = []
    table_data_1 = []
    table_data_2 = []
    table_data_3 = []
    column_names = ["Parameter", "Value ± Error"]
    a = 0
    a1 = 0
    a2 = 0
    a3 = 0
    for name, val, err in zip(par_names, fit_pars, fit_errors):
        if err is not None:
            if(err == 0):
                exp = 0
            else:
                exp = math.floor(math.log10(abs(err)))
            rounded = round(err, -exp)
            decimals = max(-exp, 0)

            if(a < 20):
                table_data.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])

            if("BF_sig" == name):
                table_data_1.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])
                table_data_2.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])
                table_data_3.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])
            if("_1" in name):
                if(a1 < 10):
                    table_data_1.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])
            if("_2" in name):
                if(a2 < 10):
                    table_data_2.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])
            if("_3" in name):
                if(a3 < 10):
                    table_data_3.append([name, f"{val:.{decimals}f} ± {err:.{decimals}f}"])
        else:
            if(a < 20):
                table_data.append([name, f"{val:.1g}"])
            if("BF_sig" == name):
                table_data_1.append([name, f"{val:.1g}"])
                table_data_2.append([name, f"{val:.1g}"])
                table_data_3.append([name, f"{val:.1g}"])        
            if("_1" in name):
                if(a1 < 10):
                    table_data_1.append([name, f"{val:.1g}"])
            if("_2" in name):
                if(a2 < 10):
                    table_data_2.append([name, f"{val:.1g}"])
            if("_3" in name):
                if(a3 < 10):
                    table_data_3.append([name, f"{val:.1g}"])

        a += 1
        if(("_1" in name) or (name == "BF_sig")):
            a1 += 1
        if(("_2" in name) or (name == "BF_sig")):
            a2 += 1
        if(("_3" in name) or (name == "BF_sig")):
            a3 += 1

    if(a >= 20):
        table_data.append(["...", "..."])
    if(a1 >= 10):
        table_data_1.append(["...", "..."])
    if(a2 >= 10):
        table_data_2.append(["...", "..."])
    if(a3 >= 10):
        table_data_3.append(["...", "..."])

    expected = model.expected_data(fit_pars)
    data = workspace.data(model)
    expected = np.array(expected)
    data = np.array(data)
    mask = expected > 0
    chi2_value = np.sum((data[mask] - expected[mask]) ** 2 / expected[mask])
    ndf = len(data) - len(fit_pars) 
    chi2_reduced = chi2_value/ndf
    p_value = 1 - chi2.cdf(chi2_value, df=ndf)
    textstr = f"$\\chi^2/ndf = {chi2_reduced:.2f}$\n$p = {p_value:.5f}$"

    # Get counts from fit
    def get_mc_counts(pars):
        deltas, factors = model._modifications(pars)
        allsum = pyhf.tensorlib.concatenate(
            deltas + [pyhf.tensorlib.astensor(model.nominal_rates)]
        )
        nom_plus_delta = pyhf.tensorlib.sum(allsum, axis=0)
        nom_plus_delta = pyhf.tensorlib.reshape(
            nom_plus_delta, (1,) + pyhf.tensorlib.shape(nom_plus_delta)
        )
        allfac = pyhf.tensorlib.concatenate(factors + [nom_plus_delta])
        return pyhf.tensorlib.product(allfac, axis=0)

    par_name_dict = {k: v["slice"].start for k, v in model.config.par_map.items()}
    par_settings = {k: fit_pars[v] for k, v in par_name_dict.items()}

    pars = pyhf.tensorlib.astensor(model.config.suggested_init())
    for k, v in par_settings.items():
        pars[par_name_dict[k]] = v

    mc_counts = get_mc_counts(pars)
    data_sig = mc_counts[model.config.samples.index("Signal")][0]
    data_comb = mc_counts[model.config.samples.index("Combinatorial")][0]
    data_BDDKp = mc_counts[model.config.samples.index("BDDKp")][0]
    if(add_BDDX_decays):
        data_BuDDK0 = mc_counts[model.config.samples.index("BuDDK0")][0]
        data_BuDD = mc_counts[model.config.samples.index("BuDD")][0]

    if(channel_type == "AllEvents"):
        if(add_BDDX_decays):
            all_data = data_sig+data_comb+data_BDDKp+data_BuDDK0+data_BuDD
            all_data = np.append(all_data, all_data[-1])

            data_sig = np.append(data_sig, data_sig[-1])
            data_comb = np.append(data_comb, data_comb[-1])
            data_BDDKp = np.append(data_BDDKp, data_BDDKp[-1])
            data_BuDDK0 = np.append(data_BuDDK0, data_BuDDK0[-1])
            data_BuDD = np.append(data_BuDD, data_BuDD[-1])

            the_data = workspace.data(model, include_auxdata=False)
            the_data_err = np.sqrt(the_data)

            if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
                full_data = np.zeros(nbins+1)
                full_data_sig = np.zeros(nbins+1)
                full_data_comb = np.zeros(nbins+1)
                full_data_BDDKp = np.zeros(nbins+1)
                full_data_BuDDK0 = np.zeros(nbins+1)
                full_data_BuDD = np.zeros(nbins+1)

                the_full_data = np.zeros(nbins)
                the_full_data_err = np.zeros(nbins)

                k = 0
                for i in range(nbins+1):
                    if((i < i_min[0]-1) or (i > i_max[0])):
                        full_data[i] = all_data[k]
                        full_data_sig[i] = data_sig[k]
                        full_data_comb[i] = data_comb[k]
                        full_data_BDDKp[i] = data_BDDKp[k]
                        full_data_BuDDK0[i] = data_BuDDK0[k]
                        full_data_BuDD[i] = data_BuDD[k]
                        k += 1
                    else:
                        full_data[i] = 0.0
                        full_data_sig[i] = 0.0
                        full_data_comb[i] = 0.0
                        full_data_BDDKp[i] = 0.0
                        full_data_BuDDK0[i] = 0.0
                        full_data_BuDD[i] = 0.0
                l = 0
                for i in range(nbins):
                    if((i < i_min[0]-1) or (i > i_max[0])):
                        the_full_data[i] = the_data[l]
                        the_full_data_err[i] = the_data_err[l]
                        l += 1
                    else:
                        the_full_data[i] = 0.0
                        the_full_data_err[i] = 0.0

                all_data = full_data
                data_sig = full_data_sig
                data_comb = full_data_comb
                data_BDDKp = full_data_BDDKp
                data_BuDDK0 = full_data_BuDDK0
                data_BuDD = full_data_BuDD
                the_data = the_full_data
                the_data_err = the_full_data_err

        else:
            all_data = data_sig+data_comb+data_BDDKp
            all_data = np.append(all_data, all_data[-1])

            data_sig = np.append(data_sig, data_sig[-1])
            data_comb = np.append(data_comb, data_comb[-1])
            data_BDDKp = np.append(data_BDDKp, data_BDDKp[-1])

            the_data = workspace.data(model, include_auxdata=False)
            the_data_err = np.sqrt(the_data)

            if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
                full_data = np.zeros(nbins+1)
                full_data_sig = np.zeros(nbins+1)
                full_data_comb = np.zeros(nbins+1)
                full_data_BDDKp = np.zeros(nbins+1)
                the_full_data = np.zeros(nbins)
                the_full_data_err = np.zeros(nbins)

                k = 0
                for i in range(nbins+1):
                    if((i < i_min[0]-1) or (i > i_max[0])):
                        full_data[i] = all_data[k]
                        full_data_sig[i] = data_sig[k]
                        full_data_comb[i] = data_comb[k]
                        full_data_BDDKp[i] = data_BDDKp[k]
                        k += 1
                    else:
                        full_data[i] = 0.0
                        full_data_sig[i] = 0.0
                        full_data_comb[i] = 0.0
                        full_data_BDDKp[i] = 0.0
                
                l = 0
                for i in range(nbins):
                    if((i < i_min[0]-1) or (i > i_max[0])):
                        the_full_data[i] = the_data[l]
                        the_full_data_err[i] = the_data_err[l]
                        l += 1
                    else:
                        the_full_data[i] = 0.0
                        the_full_data_err[i] = 0.0

                all_data = full_data
                data_sig = full_data_sig
                data_comb = full_data_comb
                data_BDDKp = full_data_BDDKp
                the_data = the_full_data
                the_data_err = the_full_data_err

        pulls = (the_data - all_data[:-1]) / the_data_err

        fig = plt.figure(figsize=(13, 8))
        gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[3, 1], wspace=0.2, hspace=0.2)

        ax_main = fig.add_subplot(gs[0, 0])
        ax_pull = fig.add_subplot(gs[1, 0], sharex=ax_main)
        ax_table = fig.add_subplot(gs[:, 1])

        # --- Main fit panel ---
        ax_main.step(edges, all_data, where='post', color='black', label="Total fit")
        ax_main.step(edges, data_sig, where='post', color='blue', label='Signal')
        ax_main.step(edges, data_comb, where='post', color='red', label='Combinatorial')
        ax_main.step(edges, data_BDDKp, where='post', color='green', label='$B \\to D D K^+$')
        if(add_BDDX_decays):
            ax_main.step(edges, data_BuDDK0, where='post', color='magenta', label='$B^+ \\to D D K^0$')
            ax_main.step(edges, data_BuDD, where='post', color='cyan', label='$B^+ \\to D D$')
        if((fit_type == "ToyData") or (fit_type == "ToyDataSidebands")):    
            ax_main.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
        else:
            ax_main.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Data")
        ax_main.text(0.55, 0.95, textstr, transform=ax_main.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        if(validate_fit):
            ax_main.set_title(f'BDT = {bdt} | seed = {seed}')
        else:
            ax_main.set_title(f'BDT = {bdt}')
        ax_main.set_ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins, 1) ))
        if(log_scale):
            ax_main.set_yscale("symlog")
        ax_main.set_xlim(4000,9000)
        ax_main.legend()

        # --- Pulls panel ---
        ax_pull.axhline(0, color='black', linewidth=1)
        ax_pull.axhline(1, color='gray', linestyle='--', linewidth=1)
        ax_pull.axhline(-1, color='gray', linestyle='--', linewidth=1)
        ax_pull.axhline(2, color='gray', linestyle=':', linewidth=1)
        ax_pull.axhline(-2, color='gray', linestyle=':', linewidth=1)
        ax_pull.errorbar(mass, pulls, yerr=np.ones_like(pulls), fmt='o', color='k')
        ax_pull.set_ylabel('Pull')
        ax_pull.set_xlabel('$m_B$ (MeV)')
        ax_pull.set_ylim(-5, 5)

        # --- Table panel ---
        ax_table.axis("off")
        ax_table.table(cellText=table_data, colLabels=column_names, loc="center", cellLoc="left", colWidths=[0.5] * len(table_data), edges='closed')

        plt.tight_layout()

    else:
        channel_nbins = model.config.channel_nbins
        ch_nbins_1 = channel_nbins["Channel_1"]
        ch_nbins_2 = channel_nbins["Channel_2"]
        ch_nbins_3 = channel_nbins["Channel_3"]

        data_sig_1 = data_sig[0:ch_nbins_1]
        data_sig_2 = data_sig[ch_nbins_1:ch_nbins_1+ch_nbins_2]
        data_sig_3 = data_sig[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

        data_comb_1 = data_comb[0:ch_nbins_1]
        data_comb_2 = data_comb[ch_nbins_1:ch_nbins_1+ch_nbins_2]
        data_comb_3 = data_comb[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]
    
        data_BDDKp_1 = data_BDDKp[0:ch_nbins_1]
        data_BDDKp_2 = data_BDDKp[ch_nbins_1:ch_nbins_1+ch_nbins_2]
        data_BDDKp_3 = data_BDDKp[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

        if(add_BDDX_decays):
            data_BuDDK0_1 = data_BuDDK0[0:ch_nbins_1]
            data_BuDDK0_2 = data_BuDDK0[ch_nbins_1:ch_nbins_1+ch_nbins_2]
            data_BuDDK0_3 = data_BuDDK0[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

            data_BuDD_1 = data_BuDD[0:ch_nbins_1]
            data_BuDD_2 = data_BuDD[ch_nbins_1:ch_nbins_1+ch_nbins_2]
            data_BuDD_3 = data_BuDD[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

            all_data_1 = data_sig_1+data_comb_1+data_BDDKp_1+data_BuDDK0_1+data_BuDD_1
            all_data_2 = data_sig_2+data_comb_2+data_BDDKp_2+data_BuDDK0_2+data_BuDD_2
            all_data_3 = data_sig_3+data_comb_3+data_BDDKp_3+data_BuDDK0_3+data_BuDD_3

            the_data = workspace.data(model, include_auxdata=False)
            the_data_err = np.sqrt(the_data)

            the_data_1 = the_data[0:ch_nbins_1]
            the_data_2 = the_data[ch_nbins_1:ch_nbins_1+ch_nbins_2]
            the_data_3 = the_data[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

            the_data_err_1 = the_data_err[0:ch_nbins_1]
            the_data_err_2 = the_data_err[ch_nbins_1:ch_nbins_1+ch_nbins_2]
            the_data_err_3 = the_data_err[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

            all_data_1 = np.append(all_data_1, all_data_1[-1])
            data_sig_1 = np.append(data_sig_1, data_sig_1[-1])
            data_comb_1 = np.append(data_comb_1, data_comb_1[-1])
            data_BDDKp_1 = np.append(data_BDDKp_1, data_BDDKp_1[-1])
            data_BuDDK0_1 = np.append(data_BuDDK0_1, data_BuDDK0_1[-1])
            data_BuDD_1 = np.append(data_BuDD_1, data_BuDD_1[-1])

            all_data_2 = np.append(all_data_2, all_data_2[-1])
            data_sig_2 = np.append(data_sig_2, data_sig_2[-1])
            data_comb_2 = np.append(data_comb_2, data_comb_2[-1])
            data_BDDKp_2 = np.append(data_BDDKp_2, data_BDDKp_2[-1])
            data_BuDDK0_2 = np.append(data_BuDDK0_2, data_BuDDK0_2[-1])
            data_BuDD_2 = np.append(data_BuDD_2, data_BuDD_2[-1])

            all_data_3 = np.append(all_data_3, all_data_3[-1])
            data_sig_3 = np.append(data_sig_3, data_sig_3[-1])
            data_comb_3 = np.append(data_comb_3, data_comb_3[-1])
            data_BDDKp_3 = np.append(data_BDDKp_3, data_BDDKp_3[-1])
            data_BuDDK0_3 = np.append(data_BuDDK0_3, data_BuDDK0_3[-1])
            data_BuDD_3 = np.append(data_BuDD_3, data_BuDD_3[-1])

            if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
                full_data_1 = np.zeros(nbins_1+1)
                full_data_2 = np.zeros(nbins_2+1)
                full_data_3 = np.zeros(nbins_3+1)

                full_data_sig_1 = np.zeros(nbins_1+1)
                full_data_sig_2 = np.zeros(nbins_2+1)
                full_data_sig_3 = np.zeros(nbins_3+1)

                full_data_comb_1 = np.zeros(nbins_1+1)
                full_data_comb_2 = np.zeros(nbins_2+1)
                full_data_comb_3 = np.zeros(nbins_3+1)

                full_data_BDDKp_1 = np.zeros(nbins_1+1)
                full_data_BDDKp_2 = np.zeros(nbins_2+1)
                full_data_BDDKp_3 = np.zeros(nbins_3+1)

                full_data_BuDDK0_1 = np.zeros(nbins_1+1)
                full_data_BuDDK0_2 = np.zeros(nbins_2+1)
                full_data_BuDDK0_3 = np.zeros(nbins_3+1)

                full_data_BuDD_1 = np.zeros(nbins_1+1)
                full_data_BuDD_2 = np.zeros(nbins_2+1)
                full_data_BuDD_3 = np.zeros(nbins_3+1)

                the_full_data_1 = np.zeros(nbins_1)
                the_full_data_2 = np.zeros(nbins_2)
                the_full_data_3 = np.zeros(nbins_3)

                the_full_data_err_1 = np.zeros(nbins_1)
                the_full_data_err_2 = np.zeros(nbins_2)
                the_full_data_err_3 = np.zeros(nbins_3)

                k1 = 0
                for i in range(nbins_1+1):
                    if((i < i_min[1]-1) or (i > i_max[1])):
                        full_data_1[i] = all_data_1[k1]
                        full_data_sig_1[i] = data_sig_1[k1]
                        full_data_comb_1[i] = data_comb_1[k1]
                        full_data_BDDKp_1[i] = data_BDDKp_1[k1]
                        full_data_BuDDK0_1[i] = data_BuDDK0_1[k1]
                        full_data_BuDD_1[i] = data_BuDD_1[k1]
                        k1 += 1
                    else:
                        full_data_1[i] = 0.0
                        full_data_sig_1[i] = 0.0
                        full_data_comb_1[i] = 0.0
                        full_data_BDDKp_1[i] = 0.0
                        full_data_BuDDK0_1[i] = 0.0
                        full_data_BuDD_1[i] = 0.0
                k2 = 0
                for i in range(nbins_2+1):
                    if((i < i_min[2]-1) or (i > i_max[2])):
                        full_data_2[i] = all_data_2[k2]
                        full_data_sig_2[i] = data_sig_2[k2]
                        full_data_comb_2[i] = data_comb_2[k2]
                        full_data_BDDKp_2[i] = data_BDDKp_2[k2]
                        full_data_BuDDK0_2[i] = data_BuDDK0_2[k2]
                        full_data_BuDD_2[i] = data_BuDD_2[k2]
                        k2 += 1
                    else:
                        full_data_2[i] = 0.0
                        full_data_sig_2[i] = 0.0
                        full_data_comb_2[i] = 0.0
                        full_data_BDDKp_2[i] = 0.0
                        full_data_BuDDK0_2[i] = 0.0
                        full_data_BuDD_2[i] = 0.0
                k3 = 0
                for i in range(nbins_3+1):
                    if((i < i_min[3]-1) or (i > i_max[3])):
                        full_data_3[i] = all_data_3[k3]
                        full_data_sig_3[i] = data_sig_3[k3]
                        full_data_comb_3[i] = data_comb_3[k3]
                        full_data_BDDKp_3[i] = data_BDDKp_3[k3]
                        full_data_BuDDK0_3[i] = data_BuDDK0_3[k3]
                        full_data_BuDD_3[i] = data_BuDD_3[k3]
                        k3 += 1
                    else:
                        full_data_3[i] = 0.0
                        full_data_sig_3[i] = 0.0
                        full_data_comb_3[i] = 0.0
                        full_data_BDDKp_3[i] = 0.0
                        full_data_BuDDK0_3[i] = 0.0
                        full_data_BuDD_3[i] = 0.0

                l1 = 0
                for i in range(nbins_1):
                    if((i < i_min[1]-1) or (i > i_max[1])):
                        the_full_data_1[i] = the_data_1[l1]
                        the_full_data_err_1[i] = the_data_err_1[l1]
                        l1 += 1
                    else:
                        the_full_data_1[i] = 0.0
                        the_full_data_err_1[i] = 0.0

                l2 = 0
                for i in range(nbins_2):
                    if((i < i_min[2]-1) or (i > i_max[2])):
                        the_full_data_2[i] = the_data_2[l2]
                        the_full_data_err_2[i] = the_data_err_2[l2]
                        l2 += 1
                    else:
                        the_full_data_2[i] = 0.0
                        the_full_data_err_2[i] = 0.0

                l3 = 0
                for i in range(nbins_3):
                    if((i < i_min[3]-1) or (i > i_max[3])):
                        the_full_data_3[i] = the_data_3[l3]
                        the_full_data_err_3[i] = the_data_err_3[l3]
                        l3 += 1
                    else:
                        the_full_data_3[i] = 0.0
                        the_full_data_err_3[i] = 0.0

                all_data_1 = full_data_1
                all_data_2 = full_data_2
                all_data_3 = full_data_3
                data_sig_1 = full_data_sig_1
                data_sig_2 = full_data_sig_2
                data_sig_3 = full_data_sig_3
                data_comb_1 = full_data_comb_1
                data_comb_2 = full_data_comb_2
                data_comb_3 = full_data_comb_3
                data_BDDKp_1 = full_data_BDDKp_1
                data_BDDKp_2 = full_data_BDDKp_2
                data_BDDKp_3 = full_data_BDDKp_3
                data_BuDDK0_1 = full_data_BuDDK0_1
                data_BuDDK0_2 = full_data_BuDDK0_2
                data_BuDDK0_3 = full_data_BuDDK0_3
                data_BuDD_1 = full_data_BuDD_1
                data_BuDD_2 = full_data_BuDD_2
                data_BuDD_3 = full_data_BuDD_3
                the_data_1 = the_full_data_1
                the_data_2 = the_full_data_2
                the_data_3 = the_full_data_3
                the_data_err_1 = the_full_data_err_1
                the_data_err_2 = the_full_data_err_2
                the_data_err_3 = the_full_data_err_3

            all_data = [all_data_1, all_data_2, all_data_3]
            data_sig = [data_sig_1, data_sig_2, data_sig_3]
            data_comb = [data_comb_1, data_comb_2, data_comb_3]
            data_BDDKp = [data_BDDKp_1, data_BDDKp_2, data_BDDKp_3]
            data_BuDDK0 = [data_BuDDK0_1, data_BuDDK0_2, data_BuDDK0_3]
            data_BuDD = [data_BuDD_1, data_BuDD_2, data_BuDD_3]
            the_data = [the_data_1, the_data_2, the_data_3]
            the_data_err = [the_data_err_1, the_data_err_2, the_data_err_3]
        
        else:
            all_data_1 = data_sig_1+data_comb_1+data_BDDKp_1
            all_data_2 = data_sig_2+data_comb_2+data_BDDKp_2
            all_data_3 = data_sig_3+data_comb_3+data_BDDKp_3

            the_data = workspace.data(model, include_auxdata=False)
            the_data_err = np.sqrt(the_data)

            the_data_1 = the_data[0:ch_nbins_1]
            the_data_2 = the_data[ch_nbins_1:ch_nbins_1+ch_nbins_2]
            the_data_3 = the_data[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

            the_data_err_1 = the_data_err[0:ch_nbins_1]
            the_data_err_2 = the_data_err[ch_nbins_1:ch_nbins_1+ch_nbins_2]
            the_data_err_3 = the_data_err[ch_nbins_1+ch_nbins_2:ch_nbins_1+ch_nbins_2+ch_nbins_3]

            all_data_1 = np.append(all_data_1, all_data_1[-1])
            data_sig_1 = np.append(data_sig_1, data_sig_1[-1])
            data_comb_1 = np.append(data_comb_1, data_comb_1[-1])
            data_BDDKp_1 = np.append(data_BDDKp_1, data_BDDKp_1[-1])

            all_data_2 = np.append(all_data_2, all_data_2[-1])
            data_sig_2 = np.append(data_sig_2, data_sig_2[-1])
            data_comb_2 = np.append(data_comb_2, data_comb_2[-1])
            data_BDDKp_2 = np.append(data_BDDKp_2, data_BDDKp_2[-1])

            all_data_3 = np.append(all_data_3, all_data_3[-1])
            data_sig_3 = np.append(data_sig_3, data_sig_3[-1])
            data_comb_3 = np.append(data_comb_3, data_comb_3[-1])
            data_BDDKp_3 = np.append(data_BDDKp_3, data_BDDKp_3[-1])

            if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
                full_data_1 = np.zeros(nbins_1+1)
                full_data_2 = np.zeros(nbins_2+1)
                full_data_3 = np.zeros(nbins_3+1)

                full_data_sig_1 = np.zeros(nbins_1+1)
                full_data_sig_2 = np.zeros(nbins_2+1)
                full_data_sig_3 = np.zeros(nbins_3+1)

                full_data_comb_1 = np.zeros(nbins_1+1)
                full_data_comb_2 = np.zeros(nbins_2+1)
                full_data_comb_3 = np.zeros(nbins_3+1)

                full_data_BDDKp_1 = np.zeros(nbins_1+1)
                full_data_BDDKp_2 = np.zeros(nbins_2+1)
                full_data_BDDKp_3 = np.zeros(nbins_3+1)

                the_full_data_1 = np.zeros(nbins_1)
                the_full_data_2 = np.zeros(nbins_2)
                the_full_data_3 = np.zeros(nbins_3)

                the_full_data_err_1 = np.zeros(nbins_1)
                the_full_data_err_2 = np.zeros(nbins_2)
                the_full_data_err_3 = np.zeros(nbins_3)

                k1 = 0
                for i in range(nbins_1+1):
                    if((i < i_min[1]-1) or (i > i_max[1])):
                        full_data_1[i] = all_data_1[k1]
                        full_data_sig_1[i] = data_sig_1[k1]
                        full_data_comb_1[i] = data_comb_1[k1]
                        full_data_BDDKp_1[i] = data_BDDKp_1[k1]
                        k1 += 1
                    else:
                        full_data_1[i] = 0.0
                        full_data_sig_1[i] = 0.0
                        full_data_comb_1[i] = 0.0
                        full_data_BDDKp_1[i] = 0.0
                
                k2 = 0
                for i in range(nbins_2+1):
                    if((i < i_min[2]-1) or (i > i_max[2])):
                        full_data_2[i] = all_data_2[k2]
                        full_data_sig_2[i] = data_sig_2[k2]
                        full_data_comb_2[i] = data_comb_2[k2]
                        full_data_BDDKp_2[i] = data_BDDKp_2[k2]
                        k2 += 1
                    else:
                        full_data_2[i] = 0.0
                        full_data_sig_2[i] = 0.0
                        full_data_comb_2[i] = 0.0
                        full_data_BDDKp_2[i] = 0.0

                k3 = 0
                for i in range(nbins_3+1):
                    if((i < i_min[3]-1) or (i > i_max[3])):
                        full_data_3[i] = all_data_3[k3]
                        full_data_sig_3[i] = data_sig_3[k3]
                        full_data_comb_3[i] = data_comb_3[k3]
                        full_data_BDDKp_3[i] = data_BDDKp_3[k3]
                        k3 += 1
                    else:
                        full_data_3[i] = 0.0
                        full_data_sig_3[i] = 0.0
                        full_data_comb_3[i] = 0.0
                        full_data_BDDKp_3[i] = 0.0

                l1 = 0
                for i in range(nbins_1):
                    if((i < i_min[1]-1) or (i > i_max[1])):
                        the_full_data_1[i] = the_data_1[l1]
                        the_full_data_err_1[i] = the_data_err_1[l1]
                        l1 += 1
                    else:
                        the_full_data_1[i] = 0.0
                        the_full_data_err_1[i] = 0.0

                l2 = 0
                for i in range(nbins_2):
                    if((i < i_min[2]-1) or (i > i_max[2])):
                        the_full_data_2[i] = the_data_2[l2]
                        the_full_data_err_2[i] = the_data_err_2[l2]
                        l2 += 1
                    else:
                        the_full_data_2[i] = 0.0
                        the_full_data_err_2[i] = 0.0

                l3 = 0
                for i in range(nbins_3):
                    if((i < i_min[3]-1) or (i > i_max[3])):
                        the_full_data_3[i] = the_data_3[l3]
                        the_full_data_err_3[i] = the_data_err_3[l3]
                        l3 += 1
                    else:
                        the_full_data_3[i] = 0.0
                        the_full_data_err_3[i] = 0.0

                all_data_1 = full_data_1
                all_data_2 = full_data_2
                all_data_3 = full_data_3
                data_sig_1 = full_data_sig_1
                data_sig_2 = full_data_sig_2
                data_sig_3 = full_data_sig_3
                data_comb_1 = full_data_comb_1
                data_comb_2 = full_data_comb_2
                data_comb_3 = full_data_comb_3
                data_BDDKp_1 = full_data_BDDKp_1
                data_BDDKp_2 = full_data_BDDKp_2
                data_BDDKp_3 = full_data_BDDKp_3
                the_data_1 = the_full_data_1
                the_data_2 = the_full_data_2
                the_data_3 = the_full_data_3
                the_data_err_1 = the_full_data_err_1
                the_data_err_2 = the_full_data_err_2
                the_data_err_3 = the_full_data_err_3

            all_data = [all_data_1, all_data_2, all_data_3]
            data_sig = [data_sig_1, data_sig_2, data_sig_3]
            data_comb = [data_comb_1, data_comb_2, data_comb_3]
            data_BDDKp = [data_BDDKp_1, data_BDDKp_2, data_BDDKp_3]
            the_data = [the_data_1, the_data_2, the_data_3]
            the_data_err = [the_data_err_1, the_data_err_2, the_data_err_3]

        pulls_1 = (the_data_1 - all_data_1[:-1]) / the_data_err_1
        pulls_2 = (the_data_2 - all_data_2[:-1]) / the_data_err_2
        pulls_3 = (the_data_3 - all_data_3[:-1]) / the_data_err_3
        pulls = [pulls_1, pulls_2, pulls_3]

        fig = plt.figure(figsize=(20, 12))
        gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[3, 1, 3], wspace=0.2, hspace=0.15)

        axes = []
        for i in range(3):
            # --- Main fit panel ---
            ax = fig.add_subplot(gs[0, i])
            ax.step(edges[i], all_data[i], where='post', color='black', label="Total fit")
            ax.step(edges[i], data_sig[i], where='post', color='blue', label='Signal')
            ax.step(edges[i], data_comb[i], where='post', color='red', label='Combinatorial')
            ax.step(edges[i], data_BDDKp[i], where='post', color='green', label='$B \\to D D K^+$')
            if(add_BDDX_decays):
                ax.step(edges[i], data_BuDDK0[i], where='post', color='magenta', label='$B^+ \\to D D K^0$')
                ax.step(edges[i], data_BuDD[i], where='post', color='cyan', label='$B^+ \\to D D$')
            if((fit_type == "ToyData") or (fit_type == "ToyDataSidebands")):
                ax.errorbar(mass[i], the_data[i], the_data_err[i], c="k", marker='.', linestyle='', zorder=99, label="Toy data")
            else:
                ax.errorbar(mass[i], the_data[i], the_data_err[i], c="k", marker='.', linestyle='', zorder=99, label="Data")
            ax.text(0.55, 0.95, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
            if(validate_fit):
                ax.set_title(f'Channel {i+1} | BDT = {bdt:.4g} | seed = {seed}')
            else:
                ax.set_title(f'Channel {i+1} | BDT = {bdt:.4g}')
            ax.set_ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins[i], 1) ))
            if(log_scale):
                ax.set_yscale("symlog")
            ax.legend()
            axes.append(ax)

            # --- Pulls panel ---
            ax1 = fig.add_subplot(gs[1, i], sharex=ax)
            ax1.axhline(0, color='black', linewidth=1)
            ax1.axhline(1, color='gray', linestyle='--', linewidth=1)
            ax1.axhline(-1, color='gray', linestyle='--', linewidth=1)
            ax1.axhline(2, color='gray', linestyle=':', linewidth=1)
            ax1.axhline(-2, color='gray', linestyle=':', linewidth=1)
            ax1.errorbar(mass[i], pulls[i], yerr=np.ones_like(pulls[i]), fmt='o', color='k')
            ax1.set_ylabel('Pull')
            ax1.set_xlabel('$m_B$ (MeV)')
            ax1.set_ylim(-5, 5)
            axes.append(ax1)

        ax_table_1 = fig.add_subplot(gs[2, 0])
        ax_table_1.axis("off")
        table1 = ax_table_1.table(cellText=table_data_1, colLabels=["Parameter", "Value ± Error"], loc="center", cellLoc="left", colWidths=[0.3]*len(table_data_1))
        table1.scale(1.5, 1.5)

        ax_table_2 = fig.add_subplot(gs[2, 1])
        ax_table_2.axis("off")
        table2 = ax_table_2.table(cellText=table_data_2, colLabels=["Parameter", "Value ± Error"], loc="center", cellLoc="left", colWidths=[0.3]*len(table_data_2))
        table2.scale(1.5, 1.5)

        ax_table_3 = fig.add_subplot(gs[2, 2])
        ax_table_3.axis("off")
        table3 = ax_table_3.table(cellText=table_data_3, colLabels=["Parameter", "Value ± Error"], loc="center", cellLoc="left", colWidths=[0.3]*len(table_data_3))
        table3.scale(1.5, 1.5)


    plt.tight_layout()
    if(validate_fit):
        fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot_seed_{seed}.pdf')
    else:
        if(comb_yield_syst == 0):
            fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot.pdf')
        elif(comb_yield_syst == 1):
            fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot_up.pdf')
        elif(comb_yield_syst == -1):
            fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_plot_down.pdf')
    plt.clf()


def run_cls_limit(fit_type, channel_type, BF_sig, bdt, fit_poi, fit_poi_error, data, model, comb_yield_syst):
    if(toy_based_limit): # low stats
        poi_values = np.linspace(0, np.abs(fit_poi)+4*fit_poi_error, 10)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q", calctype="toybased", ntoys=100)
    else: # high stats
        poi_values = np.linspace(0, np.abs(fit_poi)+5*fit_poi_error, 50)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q")

    print(f"Upper limit (obs): μ = {obs_limit:.6f}")
    print(f"Upper limit (exp): μ = {exp_limits[2]:.6f}")

    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    if(channel_type == "AllEvents"):
        ax.set_title(f"Fit to all events | BDT = {bdt}")
    else:
        ax.set_title(f"Fit in error categories | BDT = {bdt}")
    ax.set_xlabel(r"$BF_{sig}$")
    ax.set_ylabel(r"$\mathrm{CL}_{s}$")
    brazil.plot_results(scan, results, test_size=0.1, ax=ax)
    textstr1 = f"Upper limit (exp): = {exp_limits[2]:.6f}"
    plt.text(0.6, 0.6, textstr1, fontsize=15, transform=ax.transAxes)
    if(comb_yield_syst == 0):
        fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit.pdf")
    elif(comb_yield_syst == 1):
        fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit_up.pdf")
    elif(comb_yield_syst == -1):
        fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit_down.pdf")
    plt.clf()

    return exp_limits[2], obs_limit


def do_fit(fit_type, channel_type, BF_sig, bdt, i_min, i_max, spec, h_data, save_plot, comb_yield_syst=0, seed=-1):
    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    data = workspace.data(model)

    # Fit
    optimizer = pyhf.optimize.minuit_optimizer(verbose=0, maxiter=500000)
    pyhf.set_backend('numpy', optimizer) # minuit returns errors and correlation matrix; we ned the errors to make pulls

    fit_result, res_obj = pyhf.infer.mle.fit(data, model, return_uncertainties=True, return_result_obj=True)

    fit_pars, fit_errors = fit_result.T

    if(save_plot):
        print(res_obj)
        plot(fit_type, channel_type, BF_sig, bdt, spec, h_data, fit_pars, fit_errors, i_min, i_max, comb_yield_syst, False, seed)

    if(validate_fit):
        return fit_pars, fit_errors
    else:
        sig_index = model.config.par_names.index("BF_sig")
        fit_poi = fit_pars[sig_index]
        fit_poi_error = fit_errors[sig_index]

        exp_limit, obs_limit = run_cls_limit(fit_type, channel_type, BF_sig, bdt, fit_poi, fit_poi_error, data, model, comb_yield_syst)
        return exp_limit, obs_limit


def toy_studies(fit_type, channel_type, BF_sig, bdt, spec, N_fail, toy_fit_values, toy_fit_errors, n_toys, nbins=30):
    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    expected_values =  model.config.suggested_init()
    toy_labels = model.config.par_names

    sig_index = model.config.par_names.index("BF_sig")
    if(channel_type == "AllEvents"):
        comb_index = model.config.par_names.index("N_comb_0")

        N_comb = expected_values[comb_index]
        N_comb = [N_comb]
    else:
        comb_index_1 = model.config.par_names.index("N_comb_1")
        N_comb_1 = expected_values[comb_index_1]

        comb_index_2 = model.config.par_names.index("N_comb_2")
        N_comb_2 = expected_values[comb_index_2] 

        comb_index_3 = model.config.par_names.index("N_comb_3")
        N_comb_3 = expected_values[comb_index_3]
        
        N_comb = [N_comb_1, N_comb_2, N_comb_3]
    
    N = len(toy_fit_values)
    BF_toys = np.zeros(N)
    BF_toys_err = np.zeros(N)
    Nbkg_toys = np.zeros(N)
    Nbkg_toys_err = np.zeros(N)
    Nbkg_toys_1 = np.zeros(N)
    Nbkg_toys_2 = np.zeros(N)
    Nbkg_toys_3 = np.zeros(N)
    Nbkg_toys_err_1 = np.zeros(N)
    Nbkg_toys_err_2 = np.zeros(N)
    Nbkg_toys_err_3 = np.zeros(N)

    for i in range(N):
        BF_toys[i] = toy_fit_values[i][sig_index]
        BF_toys_err[i] = toy_fit_errors[i][sig_index]
        if(channel_type == "AllEvents"):
            Nbkg_toys[i] = toy_fit_values[i][comb_index]
            Nbkg_toys_err[i] = toy_fit_errors[i][comb_index]
        else:
            Nbkg_toys_1[i] = toy_fit_values[i][comb_index_1]
            Nbkg_toys_2[i] = toy_fit_values[i][comb_index_2]
            Nbkg_toys_3[i] = toy_fit_values[i][comb_index_3]

            Nbkg_toys_err_1[i] = toy_fit_errors[i][comb_index_1]
            Nbkg_toys_err_2[i] = toy_fit_errors[i][comb_index_2]
            Nbkg_toys_err_3[i] = toy_fit_errors[i][comb_index_3]

    if(channel_type == "AllEvents"):
        Nbkg_toys = [Nbkg_toys]
        Nbkg_toys_err = [Nbkg_toys_err]
    else:
        Nbkg_toys = [Nbkg_toys_1, Nbkg_toys_2, Nbkg_toys_3]
        Nbkg_toys_err = [Nbkg_toys_err_1, Nbkg_toys_err_2, Nbkg_toys_err_3]

    np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_poi_mean.npy', np.mean(BF_toys))
    np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_poi_mean_error.npy', np.sqrt(np.sum(BF_toys_err**2)) / len(BF_toys))

    for i in range(len(expected_values)):
        param_values = [toy_fit_values[j][i] for j in range(N)]
        param_errors = [toy_fit_errors[j][i] for j in range(N)]

        if(0 in param_errors):
            param_pulls = [param_values[j] - expected_values[i] for j in range(N)]
        else:
            param_pulls = [(param_values[j] - expected_values[i])/param_errors[j] for j in range(N)]

        fig1, ax1 = plt.subplots()
        (mu, sigma) = norm.fit(param_pulls)
        nn, bb, ptches = ax1.hist(param_pulls, bins=nbins)
        x = (bb[:-1] + bb[1:]) / 2
        ax1.errorbar(x, nn, yerr=np.sqrt(nn), ecolor='black', fmt='k.')
        y = sum(nn)*(bb[1] - bb[0])*norm.pdf( x, mu, sigma)
        ax1.plot(x, y, 'r-', linewidth=2)
        if(0 in param_errors):
            ax1.set_xlabel(f"Bias of {toy_labels[i]}")
        else:
            ax1.set_xlabel(f"Pull of {toy_labels[i]}")
        ax1.set_ylabel(f"Entries / {nbins}")
        ax1.set_title(f"Fit to all events ({N_fail}/{n_toys} toy fits fail) \n BDT = {bdt}")
        ax1.axvline(0, color='black', linestyle='--', linewidth=2)
        ax1.text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu, sigma/np.sqrt(N), sigma, sigma/np.sqrt(2*N)),
            transform=ax1.transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
        fig1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/pulls/{toy_labels[i]}.pdf')
        plt.clf()
        
    n_vars = len(Nbkg_toys)+1

    fig, ax = plt.subplots(n_vars, 3, figsize=(15, 5*n_vars))

    # BF pull
    if(fit_type == "ToyDataSidebands"):
        pull = (BF_toys - 0)/BF_toys_err
    else:
        pull = (BF_toys - BF_sig)/BF_toys_err
    (mu_pull, sigma_pull) = norm.fit(pull)
    n, bins, patches = ax[0,0].hist(pull, bins=nbins)
    xcenters = (bins[:-1] + bins[1:]) / 2
    ax[0,0].errorbar(xcenters, n, yerr=np.sqrt(n), ecolor='black', fmt='k.')
    y_pull = sum(n)*(bins[1] - bins[0])*norm.pdf(xcenters, mu_pull, sigma_pull)
    ax[0,0].plot(xcenters, y_pull, 'r-', linewidth=2)
    if(fit_type == "ToyDataSidebands"):
        ax[0,0].set_xlabel("Pull: ($BF_{sig}$ - 0)/error")
    else:
        ax[0,0].set_xlabel("Pull: ($BF_{sig} - BF_{sig}^{expected}$)/error")
    ax[0,0].set_ylabel(f"Entries / {nbins} bins")
    chi2_value = np.sum( ((n - y_pull)**2) / y_pull) / (nbins - 2)
    ax[0,0].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[0,0].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull, sigma_pull/np.sqrt(N), sigma_pull, sigma_pull/np.sqrt(2*N)),
            transform=ax[0,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[0,0].set_title(f"Branching fraction pull: $\chi^2$/ndf = {chi2_value:.3g}")

    # BF
    ax[0,1].hist(BF_toys, bins=nbins)
    ax[0,1].set_xlabel("$BF_{sig}$")
    ax[0,1].set_ylabel(f"Entries / {nbins} bins")
    ax[0,1].set_title("Branching fraction")
    if(fit_type == "ToyDataSidebands"):
        ax[0,1].axvline(0, color='black', linestyle='--', linewidth=2)
    else:
        ax[0,1].axvline(BF_sig, color='black', linestyle='--', linewidth=2)
    ax[0,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(BF_toys), np.std(BF_toys) ),
            transform=ax[0,1].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # BF error
    ax[0,2].hist(BF_toys_err, bins=nbins)
    ax[0,2].set_xlabel("$BF_{sig}$ error")
    ax[0,2].set_ylabel(f"Entries / {nbins} bins")
    ax[0,2].set_title("Branching fraction error")
    ax[0,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(BF_toys_err), np.std(BF_toys_err) ),
            transform=ax[0,2].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    for ch in range(n_vars-1):
        if(channel_type == "AllEvents"):
            channel_number = ch
        else:
            channel_number = ch+1

        # N_comb pull
        if(0 in Nbkg_toys_err[ch]):
            pull_2 = Nbkg_toys[ch] - N_comb[ch]
        else:
            pull_2 = (Nbkg_toys[ch] - N_comb[ch])/Nbkg_toys_err[ch]
        (mu_pull_2, sigma_pull_2) = norm.fit(pull_2)
        n_2, bins_2, patches_2 = ax[ch+1,0].hist( pull_2, bins=nbins)
        xcenters_2 = (bins_2[:-1] + bins_2[1:]) / 2
        ax[ch+1,0].errorbar(xcenters_2, n_2, yerr=np.sqrt(n_2), ecolor='black', fmt='k.')
        y_pull_2 = sum(n_2)*(bins_2[1] - bins_2[0])*norm.pdf( xcenters_2, mu_pull_2, sigma_pull_2)
        ax[ch+1,0].plot(xcenters_2, y_pull_2, 'r-', linewidth=2)
        chi2_2 = np.sum( ((n_2 - y_pull_2)**2) / y_pull_2) / (nbins - 2)
        if(0 in Nbkg_toys_err[ch]):
            ax[ch+1,0].set_xlabel("Bias: $N_{bkg} - N_{bkg}^{expected}$")
            ax[ch+1,0].set_title(f"Background yield bias {channel_number}: $\chi^2$/ndf = {chi2_2:.3g}")
        else:
            ax[ch+1,0].set_xlabel("Pull: ($N_{bkg} - N_{bkg}^{expected}$)/error")
            ax[ch+1,0].set_title(f"Background yield pull {channel_number}: $\chi^2$/ndf = {chi2_2:.3g}")
        ax[ch+1,0].set_ylabel(f"Entries / {nbins} bins")
        ax[ch+1,0].axvline(0, color='black', linestyle='--', linewidth=2)
        ax[ch+1,0].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_2, sigma_pull_2/np.sqrt(N), sigma_pull_2, sigma_pull_2/np.sqrt(2*N)),
                transform=ax[ch+1,0].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        # N_comb
        ax[ch+1,1].hist(Nbkg_toys[ch], bins=nbins)
        ax[ch+1,1].set_xlabel("$N_{bkg}$")
        ax[ch+1,1].set_ylabel(f"Entries / {nbins} bins")
        ax[ch+1,1].set_title(f"Background yield {channel_number}")
        ax[ch+1,1].axvline(N_comb[ch], color='black', linestyle='--', linewidth=2)
        ax[ch+1,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys[ch]), np.std(Nbkg_toys[ch]) ),
                transform=ax[ch+1,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        # N_comb error
        ax[ch+1,2].hist(Nbkg_toys_err[ch], bins=nbins)
        ax[ch+1,2].set_xlabel("$N_{bkg}$ error")
        ax[ch+1,2].set_ylabel(f"Entries / {nbins} bins")
        ax[ch+1,2].set_title(f"Background yield error {channel_number}")
        ax[ch+1,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_err[ch]), np.std(Nbkg_toys_err[ch]) ),
                transform=ax[ch+1,2].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))        

    if(channel_type == "AllEvents"):
        fig.suptitle(f"Fit to all events ({N_fail}/{n_toys} toy fits fail) \n BDT = {bdt}", fontsize=24)
    else:
        fig.suptitle(f"Fit in error categories ({N_fail}/{n_toys} toy fits fail) \n BDT = {bdt}", fontsize=24)
 
    plt.tight_layout()
    fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/fit_validation_plot.pdf')
    plt.clf()


def main(argv):
    fit_type = argv[1]
    channel_type = argv[2]
    BF_sig = argv[3]
    bdt = argv[4]

    BF_sig = float(BF_sig)
    bdt = float(bdt)

    if((fit_type != "ToyData") and (fit_type != "ToyDataSidebands") and (fit_type != "RSData") and (fit_type != "RSDataSidebands")):
        print("Wrong fit type. Try: 'ToyData', 'ToyDataSidebands', 'RSData' or 'RSDataSidebands' ")
        quit()

    if((channel_type != "AllEvents") and (channel_type != "ErrorCategories")):
        print("Wrong channel type. Try: 'AllEvents' or 'ErrorCategories' ")
        quit()

    # Fit templates + fit templates errors
    if((fit_type == "ToyData") or (fit_type == "RSData")):
        f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms.root")
        f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms_errors.root")
    elif((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
        f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms_sidebands.root")
        f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/histograms_sidebands_errors.root")

    if(channel_type == "AllEvents"):
        histograms, histograms_errors = retrieve_histograms(fit_type, f, f1, 0)
        histograms = [histograms]
        histograms_errors = [histograms_errors]
    elif(channel_type == "ErrorCategories"):
        histograms_1, histograms_errors_1 = retrieve_histograms(fit_type, f, f1, 1)
        histograms_2, histograms_errors_2 = retrieve_histograms(fit_type, f, f1, 2)
        histograms_3, histograms_errors_3 = retrieve_histograms(fit_type, f, f1, 3)

        histograms = [histograms_1, histograms_2, histograms_3]
        histograms_errors = [histograms_errors_1, histograms_errors_2, histograms_errors_3]

    # Normalisation values + errors
    norm_values, norm_errors = retrieve_normalisations(fit_type, bdt)

    # Indices defining the signal region
    indices = np.load(f'/panfs/felician/B2Ktautau/workflow/bdt_dependent_inputs/BDT_{bdt}/signal_region_indices.npy')
    i_min = indices[0]
    i_max = indices[1]

    if((fit_type == "RSData") and (bdt > 0.5)):
        print("Cannot unblind the RS data")
        quit()

    if(validate_fit):
        toy_fit_values = []
        toy_fit_errors = []
        N_fail = 0

        if((fit_type == "RSData") or (fit_type == "RSDataSidebands")):
            if(channel_type == "AllEvents"):
                h_data = [histograms[0][7]]
            else:
                h_data = [histograms[0][7], histograms[1][7], histograms[2][7]]
        else:
            h_data = generate_nominal_data(fit_type, channel_type, BF_sig, histograms, norm_values, norm_errors)
        spec = build_model(fit_type, channel_type, BF_sig, h_data, histograms, histograms_errors, norm_values, norm_errors, i_min, i_max)

        np.random.seed(0)
        workspace = pyhf.Workspace(spec)
        model = workspace.model()
        pars = model.config.suggested_init()
        pdf = model.make_pdf(pyhf.tensorlib.astensor(pars))
        toy_data = pdf.sample((n_toys,))
        
        optimizer = pyhf.optimize.minuit_optimizer(verbose=0, maxiter=500000)
        pyhf.set_backend('numpy', optimizer) # minuit returns errors and correlation matrix; we ned the errors to make pulls

        for seed in range(n_toys):
            if(n_toys > 10):
                if seed % (n_toys // 10) == 0:
                    print(f"Progress: {100 * seed // n_toys}% (seed={seed})")

            try:
                fit_result, res_obj = pyhf.infer.mle.fit(toy_data[seed], model, return_uncertainties=True, return_result_obj=True)
                fit_pars, fit_errors = fit_result.T

                if(seed == 0):
                    print(res_obj)

                toy_fit_values.append(fit_pars)
                toy_fit_errors.append(fit_errors)
            except:
                N_fail += 1
                continue

        toy_studies(fit_type, channel_type, BF_sig, bdt, spec, N_fail, toy_fit_values, toy_fit_errors, n_toys)
        print(f"{N_fail}/{n_toys} toys fail")

    else:
        if((fit_type == "RSData") or (fit_type == "RSDataSidebands")):
            if(channel_type == "AllEvents"):
                h_data_nominal = [histograms[0][7]]
            else:
                h_data_nominal = [histograms[0][7], histograms[1][7], histograms[2][7]]
        else:
            h_data_nominal = generate_nominal_data(fit_type, channel_type, BF_sig, histograms, norm_values, norm_errors, 0)
            h_data_up = generate_nominal_data(fit_type, channel_type, BF_sig, histograms, norm_values, norm_errors, 1)
            h_data_down = generate_nominal_data(fit_type, channel_type, BF_sig, histograms, norm_values, norm_errors, -1)

        spec_nominal = build_model(fit_type, channel_type, BF_sig, h_data_nominal, histograms, histograms_errors, norm_values, norm_errors, i_min, i_max)
        # try:
        print("CLs calculation: nominal")
        exp_limit_nominal, obs_limit_nominal = do_fit(fit_type, channel_type, BF_sig, bdt, i_min, i_max, spec_nominal, h_data_nominal, True, 0)
    
        # except:
        #     save_dummy(fit_type, channel_type, BF_sig, bdt, 0)
        #     exp_limit_nominal = np.inf
        #     obs_limit_nominal = np.inf
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit.npy', exp_limit_nominal)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_obs_limit.npy', obs_limit_nominal)

        if((fit_type == "ToyData") or (fit_type == "ToyDataSidebands")):
            spec_up = build_model(fit_type, channel_type, BF_sig, h_data_up, histograms, histograms_errors, norm_values, norm_errors, i_min, i_max)
            try:
                print("CLs calculation: nominal + error")
                exp_limit_up, obs_limit_up = do_fit(fit_type, channel_type, BF_sig, bdt, i_min, i_max, spec_up, h_data_up, True, 1)
            except:
                save_dummy(fit_type, channel_type, BF_sig, bdt, 1)
                exp_limit_up = np.inf
                obs_limit_up = np.inf
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit_up.npy', exp_limit_up)
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_obs_limit_up.npy', obs_limit_up)

            spec_down = build_model(fit_type, channel_type, BF_sig, h_data_down, histograms, histograms_errors, norm_values, norm_errors, i_min, i_max)
            try:
                print("CLs calculation: nominal - error")
                exp_limit_down, obs_limit_down = do_fit(fit_type, channel_type, BF_sig, bdt, i_min, i_max, spec_down, h_data_down, True, -1)
            except:
                save_dummy(fit_type, channel_type, BF_sig, bdt, -1)
                exp_limit_down = np.inf
                obs_limit_down = np.inf
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_limit_down.npy', exp_limit_down)
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{channel_type}/BF_sig_{BF_sig:.6g}/BDT_{bdt}/cls_obs_limit_down.npy', obs_limit_down)
    
if __name__ == "__main__":
    main(sys.argv)