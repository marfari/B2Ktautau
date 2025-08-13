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
from itertools import chain

# Global flags
add_statistical_error = False
toy_based_limit = False
validate_fit = False
comb_variation_for_limit = 0
if(validate_fit):
    n_toys = 5000

# Functions
def save_dummy(validate_fit, fit_type, fit_name, BF_sig, bdt):
    print("Fit failed. Saving dummy results.")
    if(validate_fit and not fit_type == "RSDataSidebands"):
        f, ax = plt.subplots()
        f.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/pulls/BF_sig.pdf')
        f1, ax1 = plt.subplots()
        f1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_validation_plot.pdf')
    elif(fit_type == "RSDataSidebands"):
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/ncomb_value.npy', 0)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/ncomb_error.npy', 0)    
    else:
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit.npy', np.inf)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_obs_limit.npy', np.inf)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_up_comb_yield_syst.npy', np.inf)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_obs_limit_up_comb_yield_syst.npy', np.inf)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_down_comb_yield_syst.npy', np.inf)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_obs_limit_down_comb_yield_syst.npy', np.inf)
        f, ax = plt.subplots()
        f.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot.pdf')
        f1, ax1 = plt.subplots()
        f1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit.pdf')

        f2, a2 = plt.subplots()
        f2.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot_up_comb_yield_syst.pdf')
        f3, ax3 = plt.subplots()
        f3.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_up_comb_yield_syst.pdf')

        f4, a4 = plt.subplots()
        f4.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot_down_comb_yield_syst.pdf')
        f5, ax5 = plt.subplots()
        f5.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_down_comb_yield_syst.pdf')


def retrieve_histograms(fit_type, f, f1, ch):
    h_sig = f.Get(f"Channel_{ch}/h_sig_{ch}")
    h_comb = f.Get(f"Channel_{ch}/h_comb_{ch}") # I am using WS*R as the nominal (testing)
    h_BDDKp = f.Get(f"Channel_{ch}/h_BDDKp_{ch}")
    h_upward = f.Get(f"Channel_{ch}/h_upward_{ch}") # These variations are not used in the test
    h_downward = f.Get(f"Channel_{ch}/h_downward_{ch}")

    h_sig_err = f1.Get(f"Channel_{ch}/h_sig_err_{ch}")
    h_comb_err = f1.Get(f"Channel_{ch}/h_comb_err_{ch}")
    h_BDDKp_err = f1.Get(f"Channel_{ch}/h_BDDKp_err_{ch}")

    histograms = [h_sig, h_comb, h_BDDKp]
    histogram_errors = [h_sig_err, h_comb_err, h_BDDKp_err]
    
    if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
        h_data = f.Get(f"Channel_{ch}/h_data_{ch}")
    
        return histograms, histogram_errors, h_upward, h_downward, h_data
    else:
        return histograms, histogram_errors, h_upward, h_downward


def error_category_efficiency(epsilon, epsilon_error, ch):
    # epsilon = [eps1, eps2]

    if(ch == 1):
        eps = epsilon[0]
        eps_err = epsilon_error[0]  
    elif(ch == 2):
        eps = epsilon[1]
        eps_err = epsilon_error[1]
    elif(ch == 3):
        eps_1 = epsilon[0]
        eps_err_1 = epsilon_error[0]  

        eps_2 = epsilon[1]
        eps_err_2 = epsilon_error[1]

        eps_1 = ufloat(eps_1, eps_err_1)
        eps_2 = ufloat(eps_2, eps_err_2)

        eps_ufloat = ufloat(1,0) - eps_1 - eps_2

        eps = eps_ufloat.nominal_value
        eps_err = eps_ufloat.std_dev
    else:
        eps = 1
        eps_err = 0

    return eps, eps_err


def convert_histogram_to_array(histograms, histogram_errors):
    h_sig = histograms[0]  
    h_comb = histograms[1]  
    h_BDDKp = histograms[2]

    h_sig_err = histogram_errors[0]  
    h_comb_err = histogram_errors[1]  
    h_BDDKp_err = histogram_errors[2]

    nbins = h_sig.GetNbinsX()
    data_sig = [h_sig.GetBinContent(i+1) for i in range(nbins)]
    data_comb = [h_comb.GetBinContent(i+1) for i in range(nbins)]
    data_BDDKp = [h_BDDKp.GetBinContent(i+1) for i in range(nbins)]

    data_sig_err = [h_sig_err.GetBinContent(i+1) for i in range(nbins)]
    data_comb_err = [h_comb_err.GetBinContent(i+1) for i in range(nbins)]
    data_BDDKp_err = [h_BDDKp_err.GetBinContent(i+1) for i in range(nbins)]        

    data = [data_sig, data_comb, data_BDDKp]
    data_err = [data_sig_err, data_comb_err, data_BDDKp_err]

    return data, data_err


def generate_cls_data(fit_name, fit_type, BF_sig, nominal_templates, h_upward, h_downward, norm_parameters, norm_parameters_errors, add_efficiencies, comb_yield_syst):
    # This is the data we fit to when fit_type == ToyData or fit_type == ToyDataSidebands
    # It simply scales the fit templates to the expected yields in RS data and adds them to build h_data
    # This is used to compute the expected 90% CL upper limit on the branching fraction 

    if((fit_type != "ToyData") and (fit_type != "ToyDataSidebands")):
        print("Wrong fit type for generate_cls_data. Try 'ToyData' or 'ToyDataSidebands'")

    a = norm_parameters['a']
    b = norm_parameters['b']
    C_BDDKp = norm_parameters['C_BDDKp']
    eps_sig = norm_parameters['eps_sig']
    eps_BDDKp = norm_parameters['eps_BDDKp']

    N_comb_nominal = ufloat(norm_parameters['N_comb'], norm_parameters_errors['N_comb'])
    N_sig_nominal = BF_sig*a*b
    N_BDDKp_nominal = C_BDDKp*b

    if(fit_name == "AllEvents"):
        n_channels = 1
    else:
        n_channels = 3

    histo_data = []
    for i in range(n_channels):
        if(fit_name == "AllEvents"):
            ch = i
        else:
            ch = i+1

        h_sig = nominal_templates[i][0]  
        h_comb = nominal_templates[i][1]  
        h_BDDKp = nominal_templates[i][2]
        upward_variation = h_upward[i]
        downward_variation = h_downward[i]
        nbins = h_sig.GetNbinsX()

        # data
        h_data = ROOT.TH1D(f"h_data_{ch}", f"h_data_{ch}", nbins, 4000, 9000)
        h_data_comb = ROOT.TH1D(f"h_data_comb_{ch}", f"h_data_comb_{ch}", nbins, 4000, 9000)
        h_data_sig = ROOT.TH1D(f"h_data_sig_{ch}", f"h_data_sig_{ch}", nbins, 4000, 9000)
        h_data_BDDKp = ROOT.TH1D(f"h_data_BDDKp_{ch}", f"h_data_BDDKp_{ch}", nbins, 4000, 9000)

        eps_comb_ch = ufloat( norm_parameters['eps_comb'][ch], norm_parameters_errors['eps_comb'][ch] )

        if(add_efficiencies):
            N_comb_ch = N_comb_nominal*eps_comb
            N_sig_ch =  N_sig_nominal*eps_sig
            N_BDDKp_ch = N_BDDKp_nominal*eps_BDDKp
        else:
            N_comb_ch = N_comb_nominal
            N_sig_ch = N_sig_nominal
            N_BDDKp_ch = N_BDDKp_nominal

        if(comb_yield_syst == 0):
            N_comb_ch = N_comb_ch.nominal_value
        elif(comb_yield_syst == 1):
            N_comb_ch = N_comb_ch.nominal_value + N_comb_ch.std_dev
        elif(comb_yield_syst == -1):
            if(N_comb_ch.std_dev > N_comb_ch.nominal_value):
                N_comb_ch = 0
            else:
                N_comb_ch = N_comb_ch.nominal_value - N_comb_ch.std_dev

        print("Ncomb (fit) = ", N_comb_ch)

        if(comb_variation_for_limit == 0):
            h_data_comb = h_comb.Clone(f"h_data_comb_{ch}")
            h_data_comb.Scale(N_comb_ch)
        else:
            h_comb_interpolation = ROOT.TH1D(f"h_comb_interpolation_{ch}", f"h_comb_interpolation_{ch}", nbins, 4000, 9000)

            if(comb_variation_for_limit == 1):
                histosys_alpha = 1
            elif(comb_variation_for_limit == 2):
                histosys_alpha = -1

            if(histosys_alpha >= 0):
                h_comb_interpolation = upward_variation.Clone(f"h_comb_interpolation_{ch}")
                h_comb_interpolation.SetTitle(f"h_comb_interpolation_{ch}")
                h_comb_interpolation.Sumw2()

                h_comb_interpolation.Add(h_comb, -1)
            else:
                h_comb_interpolation = h_comb.Clone(f"h_comb_interpolation_{ch}")
                h_comb_interpolation.SetTitle(f"h_comb_interpolation_{ch}")
                h_comb_interpolation.Sumw2()

                h_comb_interpolation.Add(downward_variation, -1)

            h_data_comb = h_comb.Clone(f"h_data_comb_{ch}")
            h_data_comb.SetTitle(f"h_data_comb_{ch}")
            h_data_comb.Sumw2()
            h_data_comb.Add(h_comb_interpolation, histosys_alpha)

            h_data_comb.Scale(N_comb_ch)

        h_data_sig = h_sig.Clone(f"h_data_sig_{ch}")
        h_data_sig.Scale(N_sig_ch)

        h_data_BDDKp = h_BDDKp.Clone(f"h_data_BDDKp_{ch}")
        h_data_BDDKp.Scale(N_BDDKp_ch)

        h_data.Sumw2()
        h_data.Add(h_data_sig)
        h_data.Add(h_data_comb)
        h_data.Add(h_data_BDDKp)

        histo_data.append(h_data)

    return histo_data


def generate_toy_data(fit_name, BF_sig, seed, nominal_templates, h_upward, h_downward, norm_parameters, norm_parameters_errors, add_efficiencies):
    N_comb = norm_parameters['N_comb']
    a = norm_parameters['a']
    b = norm_parameters['b']
    C_BDDKp = norm_parameters['C_BDDKp']
    eps_sig = norm_parameters['eps_sig']
    eps_comb = norm_parameters['eps_comb']
    eps_BDDKp = norm_parameters['eps_BDDKp']

    a_err = norm_parameters_errors['a']
    b_err = norm_parameters_errors['b']
    C_BDDKp_err = norm_parameters_errors['C_BDDKp']
    eps_sig_err = norm_parameters_errors['eps_sig']
    eps_comb_err = norm_parameters_errors['eps_comb']
    eps_BDDKp_err = norm_parameters_errors['eps_BDDKp']

    np.random.seed(seed)
    a = max(0, np.random.normal(a, a_err))
    b = max(0, np.random.normal(b, b_err))
    C_BDDKp = max(0, np.random.normal(C_BDDKp, C_BDDKp_err))
    
    if(fit_name == "AllEvents"):
        n_channels = 1
    else:
        n_channels = 3

    N_comb_nominal = N_comb
    N_BDDKp_nominal = C_BDDKp*b
    N_sig_nominal = BF_sig*a*b

    histo_data = []
    for i in range(n_channels):
        if(fit_name == "AllEvents"):
            ch = i
        else:
            ch = i+1

        h_sig = nominal_templates[i][0]  
        h_comb = nominal_templates[i][1]  
        h_BDDKp = nominal_templates[i][2]
        nbins = h_sig.GetNbinsX()
        upward_variation = h_upward[i]
        downward_variation = h_downward[i]

        # varied templates
        # shapesys variation
        h_comb_varied = ROOT.TH1D(f"h_comb_varied_{seed}_{ch}", f"h_comb_varied_{seed}_{ch}", nbins, 4000, 9000)
        h_sig_varied = ROOT.TH1D(f"h_sig_varied_{seed}_{ch}", f"h_sig_varied_{seed}_{ch}", nbins, 4000, 9000)
        h_BDDKp_varied = ROOT.TH1D(f"h_BDDKp_varied_{seed}_{ch}", f"h_BDDKp_varied_{seed}_{ch}", nbins, 4000, 9000)

        # # histosys variation
        # h_comb_histosys_varied = ROOT.TH1D(f"h_comb_histosys_varied_{seed}_{ch}", f"h_comb_histosys_varied_{seed}_{ch}", nbins, 4000, 9000)
        # h_comb_interpolation = ROOT.TH1D(f"h_comb_interpolation_{seed}_{ch}", f"h_comb_interpolation_{seed}_{ch}", nbins, 4000, 9000)

        # histosys_alpha = np.random.normal(0, 1)
        # if(histosys_alpha >= 0):
        #     h_comb_interpolation = upward_variation.Clone(f"h_comb_interpolation_{seed}_{ch}")
        #     h_comb_interpolation.SetTitle(f"h_comb_interpolation_{seed}_{ch}")
        #     h_comb_interpolation.Sumw2()

        #     h_comb_interpolation.Add(h_comb, -1)
        # else:
        #     h_comb_interpolation = h_comb.Clone(f"h_comb_interpolation_{seed}_{ch}")
        #     h_comb_interpolation.SetTitle(f"h_comb_interpolation_{seed}_{ch}")
        #     h_comb_interpolation.Sumw2()

        #     h_comb_interpolation.Add(downward_variation, -1)

        # h_comb_histosys_varied = h_comb.Clone(f"h_comb_histosys_varied_{seed}_{ch}")
        # h_comb_histosys_varied.SetTitle(f"h_comb_histosys_varied_{seed}_{ch}")
        # h_comb_histosys_varied.Sumw2()
        # h_comb_histosys_varied.Add(h_comb_interpolation, histosys_alpha)

        # data
        h_data = ROOT.TH1D(f"h_data_{seed}_{ch}", f"h_data_{seed}_{ch}", nbins, 4000, 9000)
        h_data_comb = ROOT.TH1D(f"h_data_comb_{seed}_{ch}", f"h_data_comb_{seed}_{ch}", nbins, 4000, 9000)
        h_data_sig = ROOT.TH1D(f"h_data_sig_{seed}_{ch}", f"h_data_sig_{seed}_{ch}", nbins, 4000, 9000)
        h_data_BDDKp = ROOT.TH1D(f"h_data_BDDKp_{seed}_{ch}", f"h_data_BDDKp_{seed}_{ch}", nbins, 4000, 9000)

        # Generate toy data from the varied templates
        # Yields
        if(add_efficiencies):
            N_comb_ch = N_comb_nominal*eps_comb[ch]

            sig_eff = max(0, np.random.normal(eps_sig[ch], eps_sig_err[ch]))
            BDDKp_eff = max(0, np.random.normal(eps_BDDKp[ch], eps_BDDKp_err[ch]))

            N_sig_ch = N_sig_nominal*sig_eff
            N_BDDKp_ch = N_BDDKp_nominal*BDDKp_eff
        else:
            N_comb_ch = N_comb_nominal
            N_sig_ch = N_sig_nominal
            N_BDDKp_ch = N_BDDKp_nominal
                
        if(add_statistical_error):
            # Vary the templates within their statistical uncertainties:
            # Combinatorial
            toy_counts = np.array([np.random.poisson( max(0, h_comb.GetBinContent(i+1)*h_comb.GetEntries() ) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts):
                h_comb_varied.SetBinContent(i+1, count)
            h_comb_varied.Sumw2()

            # Signal
            toy_counts_sig = np.array([np.random.poisson( h_sig.GetBinContent(i+1)*h_sig.GetEntries() ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_sig):
                h_sig_varied.SetBinContent(i+1, count)
            h_sig_varied.Sumw2()

            # B -> DD K+
            toy_counts_BDDKp = np.array([np.random.poisson( h_BDDKp.GetBinContent(i+1)*h_BDDKp.GetEntries() ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BDDKp):
                h_BDDKp_varied.SetBinContent(i+1, count)
            h_BDDKp_varied.Sumw2()

            # Sample from the varied templates:
            # Combinatorial
            toy_counts_1 = np.array([np.random.poisson( h_comb_varied.GetBinContent(i+1)*(N_comb_ch/h_comb_varied.Integral()) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_1):
                h_data_comb.SetBinContent(i+1, count)
            h_data_comb.Sumw2()

            # Signal
            toy_counts_sig_1 = np.array([np.random.poisson( h_sig_varied.GetBinContent(i+1)*(N_sig_ch/h_sig_varied.Integral()) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_sig_1):
                h_data_sig.SetBinContent(i+1, count)
            h_data_sig.Sumw2()

            # B -> DD K+
            toy_counts_BDDKp_1 = np.array([np.random.poisson( h_BDDKp_varied.GetBinContent(i+1)*(N_BDDKp_ch/h_BDDKp_varied.Integral()) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BDDKp_1):
                h_data_BDDKp.SetBinContent(i+1, count)
            h_data_BDDKp.Sumw2()

        else:
            # Combinatorial
            toy_counts = np.array([np.random.poisson( max(0, h_comb.GetBinContent(i+1)*N_comb_ch ) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts):
                h_data_comb.SetBinContent(i+1, count)
            h_data_comb.Sumw2()

            # Signal
            toy_counts_sig = np.array([np.random.poisson( h_sig.GetBinContent(i+1)*N_sig_ch ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_sig):
                h_data_sig.SetBinContent(i+1, count)
            h_data_sig.Sumw2()

            # B -> DD K+
            toy_counts_BDDKp = np.array([np.random.poisson( h_BDDKp.GetBinContent(i+1)*N_BDDKp_ch ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BDDKp):
                h_data_BDDKp.SetBinContent(i+1, count)
            h_data_BDDKp.Sumw2()

        h_data.Sumw2()
        h_data.Add(h_data_sig)
        h_data.Add(h_data_comb)
        h_data.Add(h_data_BDDKp)

        histo_data.append(h_data)

    return histo_data


def build_model(fit_name, fit_type, BF_sig, nominal_templates, error_templates, norm_parameters, norm_parameters_errors, h_data, h_upward, h_downward, i_min, i_max, add_efficiencies):
    ### NormFactor expected values
    N_comb = norm_parameters['N_comb']
    a = norm_parameters['a']
    b = norm_parameters['b']
    C_BDDKp = norm_parameters['C_BDDKp']
    eps_sig = norm_parameters['eps_sig']
    eps_comb = norm_parameters['eps_comb']
    eps_BDDKp = norm_parameters['eps_BDDKp']

    # print("Channel 1: ", "eps_sig = ", eps_sig[1], " | eps_comb = ", eps_comb[1], " | eps_BDDKp = ", eps_BDDKp[1])
    # print("Channel 2: ", "eps_sig = ", eps_sig[2], " | eps_comb = ", eps_comb[2], " | eps_BDDKp = ", eps_BDDKp[2])
    # print("Channel 3: ", "eps_sig = ", eps_sig[3], " | eps_comb = ", eps_comb[3], " | eps_BDDKp = ", eps_BDDKp[3])

    a_err = norm_parameters_errors['a']
    b_err = norm_parameters_errors['b']
    C_BDDKp_err = norm_parameters_errors['C_BDDKp']
    eps_sig_err = norm_parameters_errors['eps_sig']
    eps_comb_err = norm_parameters_errors['eps_comb']
    eps_BDDKp_err = norm_parameters_errors['eps_BDDKp']

    ### Values for Gaussian constraints
    if((a == 0) or (a_err > a)):
        upward_a = 1.999
        downward_a = 0.001
    else:
        upward_a = 1+(a_err/a)
        downward_a = 1-(a_err/a)

    if((b == 0) or (b_err > b)):
        upward_b = 1.999
        downward_b = 0.001
    else:
        upward_b = 1+(b_err/b)
        downward_b = 1-(b_err/b)

    if((C_BDDKp == 0) or (C_BDDKp_err > C_BDDKp)):
        upward_BDDKp = 1.999
        downward_BDDKp = 0.001
    else:
        upward_BDDKp = 1+(C_BDDKp_err/C_BDDKp)
        downward_BDDKp = 1-(C_BDDKp_err/C_BDDKp)

    if(fit_name == "AllEvents"):
        n_channels = 1
    else:
        n_channels = 3

    # Build samples and parameters
    channel_samples = []
    parameters = [{"name": "lumi", "auxdata": [1.0], "bounds": [[0.5,1.5]], "fixed": True, "inits": [1.0], "sigmas": [0.1]}, {"name": "BF_sig", "bounds": [[-1,1]], "inits": [BF_sig]}, {"name": "a", "bounds": [[-2*(a+1),2*(a+1)]], "fixed": True, "inits": [a]}, {"name": "b", "bounds": [[-2*(b+1),2*(b+1)]], "fixed": True, "inits": [b]}, {"bounds": [[-2*(C_BDDKp+1),2*(C_BDDKp+1)]], "fixed": True, "inits": [C_BDDKp], "name": "C_BDDKp"}]
    
    for i in range(n_channels):
        if(fit_name == "AllEvents"):
            ch = i
        else:
            ch = i+1

        # Templates
        template_data, template_error = convert_histogram_to_array(nominal_templates[i], error_templates[i])

        data_sig = template_data[0]
        data_comb = template_data[1]
        data_BDDKp = template_data[2]

        data_sig_err = template_error[0]
        data_comb_err = template_error[1]
        data_BDDKp_err = template_error[2]

        if(add_efficiencies):
            if((eps_sig[ch] == 0) or (eps_sig_err[ch] > eps_sig[ch])):
                upward_eps_sig = 1.999
                downward_eps_sig = 0.001
            else:
                upward_eps_sig = 1+(eps_sig_err[ch]/eps_sig[ch])
                downward_eps_sig = 1-(eps_sig_err[ch]/eps_sig[ch])

            if((eps_BDDKp[ch] == 0) or (eps_BDDKp_err[ch] > eps_BDDKp[ch])):
                upward_eps_BDDKp = 1.999
                downward_eps_BDDKp = 0.001
            else:
                upward_eps_BDDKp = 1+(eps_BDDKp_err[ch]/eps_BDDKp[ch])
                downward_eps_BDDKp = 1-(eps_BDDKp_err[ch]/eps_BDDKp[ch])


        data_comb_upward = [h_upward[ch].GetBinContent(i+1) for i in range(len(data_comb))]
        data_comb_downward = [h_downward[ch].GetBinContent(i+1) for i in range(len(data_comb))]
        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            data_sig = [data_sig[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_sig)) )]
            data_comb = [data_comb[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_comb)) )]
            data_BDDKp = [data_BDDKp[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BDDKp)))]

            data_sig_err = [data_sig_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_sig_err)) )]
            data_comb_err = [data_comb_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_comb_err)) )]
            data_BDDKp_err = [data_BDDKp_err[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_BDDKp_err)))]

            data_comb_upward = [data_comb_upward[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_comb_upward)))]
            data_comb_downward = [data_comb_downward[k] for k in chain(range(0, i_min[ch]-1), range(i_max[ch], len(data_comb_downward)))]

        # Modifiers 
        sig_modifiers = [{"data": None, "name": "BF_sig", "type": "normfactor"}, {"data": None, "name": "a", "type": "normfactor"}, {"data": None, "name": "b", "type": "normfactor"}]
        comb_modifiers = [{"data": None, "name": f"N_comb_{ch}", "type": "normfactor"}]   
        BDDKp_modifiers = [{"data": None, "name": "C_BDDKp", "type": "normfactor"}, {"data": None, "name": "b", "type": "normfactor"}]
        
        sig_modifiers.append({"data": {"hi": upward_a, "lo": downward_a}, "name": "alpha_a", "type": "normsys"})    
        sig_modifiers.append({"data": {"hi": upward_b, "lo": downward_b}, "name": "alpha_b", "type": "normsys"})
        if(C_BDDKp != 0):
            BDDKp_modifiers.append({"data": {"hi": upward_b, "lo": downward_b}, "name": "alpha_b", "type": "normsys"})
            BDDKp_modifiers.append({"data": {"hi": upward_BDDKp, "lo": downward_BDDKp}, "name": "alpha_C", "type": "normsys"})

        if(add_efficiencies):
            sig_modifiers.append({"data": None, "name": f"eps_sig_{ch}", "type": "normfactor"}) 
            sig_modifiers.append({"data": {"hi": upward_eps_sig, "lo": downward_eps_sig}, "name": f"alpha_eps_sig_{ch}", "type": "normsys"})

            BDDKp_modifiers.append({"data": None, "name": f"eps_BDDKp_{ch}", "type": "normfactor"})
            if(eps_BDDKp[ch] != 0):
                BDDKp_modifiers.append({"data": {"hi": upward_eps_BDDKp, "lo": downward_eps_BDDKp}, "name": f"alpha_eps_BDDKp_{ch}", "type": "normsys"})

        if(add_statistical_error):
            sig_modifiers.append({"name": f"shapesys_sig_{ch}", "data": data_sig_err, "type": "shapesys"})
            comb_modifiers.append({"name": f"shapesys_comb_{ch}", "data": data_comb_err, "type": "shapesys"})
            BDDKp_modifiers.append({"name": f"shapesys_BDDKp_{ch}", "data": data_BDDKp_err, "type": "shapesys"})

        # histosys: shape differences between RS and WS data
        comb_modifiers.append({"name": f"histosys_comb_{ch}", "data": {"hi_data": data_comb_upward, "lo_data": data_comb_downward}, "type": "histosys"})

        samples = [{"name": "Signal", "data": data_sig, "modifiers": sig_modifiers}]
        samples.append({"name": "Combinatorial", "data": data_comb, "modifiers": comb_modifiers})
        samples.append({"name": "BDDKp", "data": data_BDDKp, "modifiers": BDDKp_modifiers})   

        channel_samples.append(samples)

        # Parameters
        if(add_efficiencies):
            parameters.append({"name": f"eps_sig_{ch}", "bounds": [[0.0,2.0]], "fixed": True, "inits": [eps_sig[ch]]})
            parameters.append({"name": f"eps_BDDKp_{ch}", "bounds": [[0.0,2.0]], "fixed": True, "inits": [eps_BDDKp[ch]]})

        parameters.append({"name": f"N_comb_{ch}", "bounds": [[-2*(N_comb*eps_comb[ch]+1), 2*(N_comb*eps_comb[ch]+1)]], "inits": [N_comb*eps_comb[ch]]})

    # Build spec
    if(fit_name == "AllEvents"):
        nbins = nominal_templates[0][0].GetNbinsX()
        channels = [{"name": "AllEvents", "samples": channel_samples[0]}]

        data_values = [h_data[0].GetBinContent(i+1) for i in range(nbins)]
        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            data_values = [data_values[k] for k in chain(range(0, i_min[0]-1), range(i_max[0], len(data_values)) )]

        observations = [{"name": "AllEvents", "data": data_values}]        
    else:
        nbins_1 = nominal_templates[0][0].GetNbinsX()
        nbins_2 = nominal_templates[1][0].GetNbinsX()
        nbins_3 = nominal_templates[2][0].GetNbinsX()

        data_values_1 = [h_data[0].GetBinContent(i+1) for i in range(nbins_1)]
        data_values_2 = [h_data[1].GetBinContent(i+1) for i in range(nbins_2)]
        data_values_3 = [h_data[2].GetBinContent(i+1) for i in range(nbins_3)]

        if((fit_type == "ToyDataSidebands") or (fit_type == "RSDataSidebands")):
            data_values_1 = [data_values_1[k] for k in chain(range(0, i_min[1]-1), range(i_max[1], len(data_values_1)) )]
            data_values_2 = [data_values_2[k] for k in chain(range(0, i_min[2]-1), range(i_max[2], len(data_values_2)) )]
            data_values_3 = [data_values_3[k] for k in chain(range(0, i_min[3]-1), range(i_max[3], len(data_values_3)) )]

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


def plot(fit_name, fit_type, h_data, workspace, fit_pars, fit_errors, BF_sig, bdt, i_min, i_max, add_efficiencies, comb_yield_syst, seed=-1):
    model = workspace.model()

    if(fit_name == "AllEvents"):
        nbins = h_data[0][0].GetNbinsX()
        mass = [h_data[0][0].GetBinCenter(i+1) for i in range(nbins)] 
        edges = [h_data[0][0].GetBinLowEdge(i+1) for i in range(nbins+1)] 
    else:
        nbins_1 = h_data[0][0].GetNbinsX()
        nbins_2 = h_data[1][0].GetNbinsX()
        nbins_3 = h_data[2][0].GetNbinsX()

        mass_1 = [h_data[0][0].GetBinCenter(i+1) for i in range(nbins_1)] 
        mass_2 = [h_data[1][0].GetBinCenter(i+1) for i in range(nbins_2)] 
        mass_3 = [h_data[2][0].GetBinCenter(i+1) for i in range(nbins_3)] 

        edges_1 = [h_data[0][0].GetBinLowEdge(i+1) for i in range(nbins_1+1)] 
        edges_2 = [h_data[1][0].GetBinLowEdge(i+1) for i in range(nbins_2+1)] 
        edges_3 = [h_data[2][0].GetBinLowEdge(i+1) for i in range(nbins_3+1)] 

    ################################################ Retrieve parameters and indices ###############################################################
    # Branching fraction
    sig_index = model.config.par_names.index("BF_sig")
    fit_poi = fit_pars[sig_index] 
    fit_poi_error = fit_errors[sig_index] 
    BF_ufloat = ufloat(fit_poi, fit_poi_error)
    print("BF_sig = ", BF_ufloat)

    # N_comb
    if(fit_name == "AllEvents"):
        comb_index = model.config.par_names.index("N_comb_0")
        ncomb = fit_pars[comb_index]
        ncomb_err = fit_errors[comb_index]
        ncomb_ufloat = ufloat(ncomb, ncomb_err)
        print("N_comb = ", ncomb_ufloat)
    else:
        comb_index_1 = model.config.par_names.index("N_comb_1")
        ncomb_1 = fit_pars[comb_index_1]
        ncomb_err_1 = fit_errors[comb_index_1]
        ncomb_ufloat_1 = ufloat(ncomb_1, ncomb_err_1)
        print("N_comb_1 = ", ncomb_ufloat_1)


        comb_index_2 = model.config.par_names.index("N_comb_2")
        ncomb_2 = fit_pars[comb_index_2]
        ncomb_err_2 = fit_errors[comb_index_2]
        ncomb_ufloat_2 = ufloat(ncomb_2, ncomb_err_2)
        print("N_comb_2 = ", ncomb_ufloat_2)

        comb_index_3 = model.config.par_names.index("N_comb_3")
        ncomb_3 = fit_pars[comb_index_3]
        ncomb_err_3 = fit_errors[comb_index_3]
        ncomb_ufloat_3 = ufloat(ncomb_3, ncomb_err_3)
        print("N_comb_3 = ", ncomb_ufloat_3)

    # a
    a_index = model.config.par_names.index("a")
    a = fit_pars[a_index] 

    # b
    b_index = model.config.par_names.index("b")
    b = fit_pars[b_index]

    # C
    C_BDDKp_index = model.config.par_names.index("C_BDDKp")
    C = fit_pars[C_BDDKp_index]

    # Channel efficiencies
    if(add_efficiencies):
        if(fit_name == "AllEvents"):
            eps_sig_index = model.config.par_names.index("eps_sig_0")
            eps_sig = fit_pars[eps_sig_index]

            eps_BDDKp_index = model.config.par_names.index("eps_BDDKp_0")
            eps_BDDKp = fit_pars[eps_BDDKp_index]
        else:
            eps_sig_1_index = model.config.par_names.index("eps_sig_1")
            eps_sig_1 = fit_pars[eps_sig_1_index]
            eps_BDDKp_1_index = model.config.par_names.index("eps_BDDKp_1")
            eps_BDDKp_1 = fit_pars[eps_BDDKp_1_index]
            
            eps_sig_2_index = model.config.par_names.index("eps_sig_2")
            eps_sig_2 = fit_pars[eps_sig_2_index]
            eps_BDDKp_2_index = model.config.par_names.index("eps_BDDKp_2")
            eps_BDDKp_2 = fit_pars[eps_BDDKp_2_index]

            eps_sig_3_index = model.config.par_names.index("eps_sig_3")
            eps_sig_3 = fit_pars[eps_sig_3_index]
            eps_BDDKp_3_index = model.config.par_names.index("eps_BDDKp_3")
            eps_BDDKp_3 = fit_pars[eps_BDDKp_3_index]

    nsig_ufloat = BF_ufloat*a*b
    n_BDDKp_ufloat = C*b
    if(add_efficiencies):
        if(fit_name == "AllEvents"):
            nsig_ufloat *= eps_sig
            n_BDDKp_ufloat *= eps_BDDKp
        else:
            nsig_ufloat_1 = nsig_ufloat*eps_sig_1
            n_BDDKp_ufloat_1 = n_BDDKp_ufloat*eps_BDDKp_1

            nsig_ufloat_2 = nsig_ufloat*eps_sig_2
            n_BDDKp_ufloat_2 = n_BDDKp_ufloat*eps_BDDKp_2

            nsig_ufloat_3 = nsig_ufloat*eps_sig_3
            n_BDDKp_ufloat_3 = n_BDDKp_ufloat*eps_BDDKp_3

    ###################################################################################################################################################

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

    if(fit_name == "AllEvents"):
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

            f, ax = plt.subplots()
            plt.step(edges, full_data, where='post', color='black', label="Total fit")
            plt.step(edges, full_data_sig, where='post', color='blue', label='Signal')
            plt.step(edges, full_data_comb, where='post', color='red', label='Combinatorial')
            plt.step(edges, full_data_BDDKp, where='post', color='green', label='$B \\to D D K^+$')

            plt.errorbar(mass, the_full_data, the_full_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
            plt.title(f'All events \n BDT = {bdt} | seed = {seed}')
            plt.xlabel('m_B (MeV)')
            plt.ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins, 1) ))
            # ax.set_yscale("symlog")
            plt.xlim(4000,9000)
        
        else:
            f, ax = plt.subplots()
            plt.step(edges, all_data, where='post', color='black', label="Total fit")
            plt.step(edges, data_sig, where='post', color='blue', label='Signal')
            plt.step(edges, data_comb, where='post', color='red', label='Combinatorial')
            plt.step(edges, data_BDDKp, where='post', color='green', label='$B \\to D D K^+$')

            plt.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
            plt.title(f'All events \n BDT = {bdt} | seed = {seed}')
            plt.xlabel('m_B (MeV)')
            plt.ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins, 1) ))
            ax.set_yscale("symlog")
            plt.xlim(4000,9000)

        textstr = f"$BF_{{sig}} = $ {BF_ufloat} \n $N_{{sig}} = $ {nsig_ufloat:.1f} "
        textstr += f"\n $N_{{comb}} = $ {ncomb_ufloat} "
        textstr += f"\n $N_{{BDDKp}} = $ {n_BDDKp_ufloat:.1f} "

        plt.text(0.3, 0.8, textstr, fontsize=10, transform=ax.transAxes)
        plt.tight_layout()
        plt.legend()

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

            figure, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
            ax1.step(edges_1, full_data_1, where='post', color='black', label="Total fit")
            ax1.step(edges_1, full_data_sig_1, where='post', color='blue', label='Signal')
            ax1.step(edges_1, full_data_comb_1, where='post', color='red', label='Combinatorial')
            ax1.step(edges_1, full_data_BDDKp_1, where='post', color='green', label='$B \\to D D K^+$')
            ax1.errorbar(mass_1, the_full_data_1, the_full_data_err_1, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

            ax2.step(edges_2, full_data_2, where='post', color='black', label="Total fit")
            ax2.step(edges_2, full_data_sig_2, where='post', color='blue', label='Signal')
            ax2.step(edges_2, full_data_comb_2, where='post', color='red', label='Combinatorial')
            ax2.step(edges_2, full_data_BDDKp_2, where='post', color='green', label='$B \\to D D K^+$')
            ax2.errorbar(mass_2, the_full_data_2, the_full_data_err_2, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

            ax3.step(edges_3, full_data_3, where='post', color='black', label="Total fit")
            ax3.step(edges_3, full_data_sig_3, where='post', color='blue', label='Signal')
            ax3.step(edges_3, full_data_comb_3, where='post', color='red', label='Combinatorial')
            ax3.step(edges_3, full_data_BDDKp_3, where='post', color='green', label='$B \\to D D K^+$')
            ax3.errorbar(mass_3, the_full_data_3, the_full_data_err_3, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
        else:
            figure, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
            ax1.step(edges_1, all_data_1, where='post', color='black', label="Total fit")
            ax1.step(edges_1, data_sig_1, where='post', color='blue', label='Signal')
            ax1.step(edges_1, data_comb_1, where='post', color='red', label='Combinatorial')
            ax1.step(edges_1, data_BDDKp_1, where='post', color='green', label='$B \\to D D K^+$')
            ax1.errorbar(mass_1, the_data_1, the_data_err_1, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

            ax2.step(edges_2, all_data_2, where='post', color='black', label="Total fit")
            ax2.step(edges_2, data_sig_2, where='post', color='blue', label='Signal')
            ax2.step(edges_2, data_comb_2, where='post', color='red', label='Combinatorial')
            ax2.step(edges_2, data_BDDKp_2, where='post', color='green', label='$B \\to D D K^+$')
            ax2.errorbar(mass_2, the_data_2, the_data_err_2, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

            ax3.step(edges_3, all_data_3, where='post', color='black', label="Total fit")
            ax3.step(edges_3, data_sig_3, where='post', color='blue', label='Signal')
            ax3.step(edges_3, data_comb_3, where='post', color='red', label='Combinatorial')
            ax3.step(edges_3, data_BDDKp_3, where='post', color='green', label='$B \\to D D K^+$')
            ax3.errorbar(mass_3, the_data_3, the_data_err_3, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax1.set_title(f'Channel 1 \n BDT = {bdt:.4g} | seed = {seed}')
        ax1.set_xlabel('m_B (MeV)')
        ax1.set_ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins_1, 1) ))
        ax1.set_yscale("symlog")
        textstr_1 = f"$BF_{{sig}} = $ {BF_ufloat}  \n $N_{{sig}} = $ {nsig_ufloat_1:.1f} "
        textstr_1 += f"\n $N_{{comb}} = $ {ncomb_ufloat_1} "
        textstr_1 += f"\n $N_{{BDDKp}} = $ {n_BDDKp_ufloat_1:.1f} "
        ax1.text(0.4, 0.8, textstr_1, fontsize=12, transform=ax1.transAxes)
        ax1.legend()

        ax2.set_title(f'Channel 2 \n BDT = {bdt:.4g} | seed = {seed}')
        ax2.set_xlabel('m_B (MeV)')
        ax2.set_ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins_2, 1) ))
        ax2.set_yscale("symlog")
        textstr_2 = f"$BF_{{sig}} = $ {BF_ufloat}  \n $N_{{sig}} = $ {nsig_ufloat_2:.1f} "
        textstr_2 += f"\n $N_{{comb}} = $ {ncomb_ufloat_2} "
        textstr_2 += f"\n $N_{{BDDKp}} = $ {n_BDDKp_ufloat_2:.1f} "
        ax2.text(0.4, 0.8, textstr_2, fontsize=12, transform=ax2.transAxes)
        ax2.legend()

        ax3.set_title(f'Channel 3 \n BDT = {bdt:.4g} | seed = {seed}')
        ax3.set_xlabel('m_B (MeV)')
        ax3.set_ylabel('Entries / ({0} MeV)'.format( round((9000-4000)/nbins_3, 1) ))
        ax3.set_yscale("symlog")
        textstr_3 = f"$BF_{{sig}} = $ {BF_ufloat}  \n $N_{{sig}} = $ {nsig_ufloat_3:.1f} "
        textstr_3 += f"\n $N_{{comb}} = $ {ncomb_ufloat_3} "
        textstr_3 += f"\n $N_{{BDDKp}} = $ {n_BDDKp_ufloat_3:.1f} "
        ax3.text(0.4, 0.8, textstr_3, fontsize=12, transform=ax3.transAxes)
        ax3.legend()

    plt.tight_layout()
    if(validate_fit and not fit_type == "RSDataSidebands"):
        plt.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot_seed_{seed}.pdf')
    else:
        if(comb_yield_syst == 0):
            plt.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot.pdf') 
        elif(comb_yield_syst == 1):
            plt.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot_up_comb_yield_syst.pdf') 
        elif(comb_yield_syst == -1):
            plt.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_plot_down_comb_yield_syst.pdf') 
    plt.clf()


def run_cls_limit(fit_poi, fit_poi_error, data, model, fit_name, fit_type, BF_sig, bdt, comb_yield_syst):
    print("### CLS LIMIT CALCULATION")
    # CLs limit: evaluate upper limit on parameter of interest
    if(toy_based_limit): # low stats
        poi_values = np.linspace(0, np.abs(fit_poi)+4*fit_poi_error, 10)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q", calctype="toybased", ntoys=100)
    else: # high stats
        poi_values = np.linspace(0, np.abs(fit_poi)+4*fit_poi_error, 50)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q")

    print(f"Upper limit (obs): μ = {obs_limit:.6f}")
    print(f"Upper limit (exp): μ = {exp_limits[2]:.6f}")

    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    ax.set_title(f"Fit {fit_name} \n BDT = {bdt}")
    ax.set_xlabel(r"$BF_{sig}$")
    ax.set_ylabel(r"$\mathrm{CL}_{s}$")
    brazil.plot_results(scan, results, test_size=0.1, ax=ax)
    textstr1 = f"Upper limit (exp): = {exp_limits[2]:.6f}"
    plt.text(0.6, 0.6, textstr1, fontsize=15, transform=ax.transAxes)
    if(comb_yield_syst == 0):
        fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit.pdf")
    elif(comb_yield_syst == 1):
        fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_up_comb_yield_syst.pdf")
    elif(comb_yield_syst == -1):
        fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_down_comb_yield_syst.pdf")
    plt.clf()

    return exp_limits[2], obs_limit


def do_fit(fit_name, fit_type, BF_sig, bdt, spec, h_templates, save_plot, i_min, i_max, add_efficiencies, comb_yield_syst, seed=-1):        
    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    data = workspace.data(model)

    # initial_parameters = model.config.suggested_init()
    # ll = pyhf.infer.mle.twice_nll(initial_parameters, data, model)
    # print(ll)

    # Fitting
    optimizer = pyhf.optimize.minuit_optimizer(verbose=0, maxiter=500000)
    pyhf.set_backend('numpy', optimizer) # minuit returns errors and correlation matrix; we ned the errors to make pulls

    fit_result, res_obj = pyhf.infer.mle.fit(data, model, return_uncertainties=True, return_result_obj=True)

    # cov_matrix = res_obj.minuit.covariance
    # num_rows, num_cols = cov_matrix.shape
    # corr_matrix = np.zeros((num_rows,num_cols))

    # for i in range(num_rows):
    #     for j in range(num_cols):
    #         corr_matrix[i][j] = cov_matrix[i][j]/np.sqrt( cov_matrix[i][i] *  cov_matrix[j][j] )

    # # df = pd.DataFrame(corr_matrix)
    # correlation_value = corr_matrix[0][1]
    # print("BF_sig and Nbkg correlation = ", correlation_value)

    fit_pars, fit_errors = fit_result.T

    if(save_plot and not fit_type == "RSDataSidebands"):
        print(res_obj)
        plot(fit_name, fit_type, h_templates, workspace, fit_pars, fit_errors, BF_sig, bdt, i_min, i_max, add_efficiencies, comb_yield_syst, seed)

    if(not validate_fit and not fit_type == "RSDataSidebands"):
        sig_index = model.config.par_names.index("BF_sig")
        fit_poi = fit_pars[sig_index]
        fit_poi_error = fit_errors[sig_index]

        exp_limit, obs_limit = run_cls_limit(fit_poi, fit_poi_error, data, model, fit_name, fit_type, BF_sig, bdt, comb_yield_syst)

    if(fit_type == "RSDataSidebands"):
        comb_index = model.config.par_names.index("N_comb_0")
        ncomb = fit_pars[comb_index]
        ncomb_err = fit_errors[comb_index]
        ncomb_ufloat = ufloat(ncomb, ncomb_err)
        print("N_comb = ", ncomb_ufloat)
        return ncomb_ufloat
    else:
        if(validate_fit):
            return fit_pars, fit_errors
        else:
            return exp_limit, obs_limit


def toy_studies(fit_name, fit_type, BF_sig, bdt, model, N_fail, toy_fit_values, toy_fit_errors, nbins=30):
    expected_values =  model.config.suggested_init()
    toy_labels = model.config.par_names

    sig_index = model.config.par_names.index("BF_sig")
    if(fit_name == "AllEvents"):
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
        if(fit_name == "AllEvents"):
            Nbkg_toys[i] = toy_fit_values[i][comb_index]
            Nbkg_toys_err[i] = toy_fit_errors[i][comb_index]
        else:
            Nbkg_toys_1[i] = toy_fit_values[i][comb_index_1]
            Nbkg_toys_2[i] = toy_fit_values[i][comb_index_2]
            Nbkg_toys_3[i] = toy_fit_values[i][comb_index_3]

            Nbkg_toys_err_1[i] = toy_fit_errors[i][comb_index_1]
            Nbkg_toys_err_2[i] = toy_fit_errors[i][comb_index_2]
            Nbkg_toys_err_3[i] = toy_fit_errors[i][comb_index_3]

    if(fit_name == "AllEvents"):
        Nbkg_toys = [Nbkg_toys]
        Nbkg_toys_err = [Nbkg_toys_err]
    else:
        Nbkg_toys = [Nbkg_toys_1, Nbkg_toys_2, Nbkg_toys_3]
        Nbkg_toys_err = [Nbkg_toys_err_1, Nbkg_toys_err_2, Nbkg_toys_err_3]

    np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_poi_mean.npy', np.mean(BF_toys))
    np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_poi_mean_error.npy', np.sqrt(np.sum(BF_toys_err**2)) / len(BF_toys))

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
        ax1.set_title(f"Fit to all events ({N_fail}/{n_toys} fits fail) \n BDT = {bdt}")
        ax1.axvline(0, color='black', linestyle='--', linewidth=2)
        ax1.text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu, sigma/np.sqrt(N), sigma, sigma/np.sqrt(2*N)),
            transform=ax1.transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
        fig1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/pulls/{toy_labels[i]}.pdf')
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
    chi2 = np.sum( ((n - y_pull)**2) / y_pull) / (nbins - 2)
    ax[0,0].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[0,0].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull, sigma_pull/np.sqrt(N), sigma_pull, sigma_pull/np.sqrt(2*N)),
            transform=ax[0,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[0,0].set_title(f"Branching fraction pull: $\chi^2$/ndf = {chi2:.3g}")

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
        if(fit_name == "AllEvents"):
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

    if(fit_name == "AllEvents"):
        fig.suptitle(f"Fit to all events ({N_fail}/{n_toys} fits fail) \n BDT = {bdt}", fontsize=24)
    else:
        fig.suptitle(f"Fit in error categories ({N_fail}/{n_toys} fits fail) \n BDT = {bdt}", fontsize=24)
 
    plt.tight_layout()
    fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/fit_validation_plot.pdf')
    plt.clf()


def main(argv):

    fit_name = argv[1]
    fit_type = argv[2]
    BF_sig = argv[3]
    bdt = argv[4]

    BF_sig = float(BF_sig)
    bdt = float(bdt)

    add_efficiencies = True
    if((fit_name == "AllEvents") and ((fit_type == "ToyData") or (fit_type == "RSData"))):
        add_efficiencies = False
    if((fit_name == "RSDataSidebands")):
        add_statistical_error = False
    
    start = time.time()

    run_fit = True

    if((fit_type != "ToyDataSidebands") and (fit_type != "RSDataSidebands") and (fit_type != "ToyData") and (fit_type != "RSData")):
        print("Wrong fit type. Try: 'ToyDataSidebands': validation of signal region definition \n 'RSDataSidebands': extraction of Ncomb \n 'ToyData': calculation of expected upper limit (blinded) \n 'RSData': calculation of observed upper limit (unblinded) ")
        quit()

    if((fit_name != "AllEvents") and (fit_name != "ErrorCategories")):
        print("Wrong fit name. Try 'AllEvents' or 'ErrorCategories'")
        quit()

    signal_region_range = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/signal_region_indices.npy')
    i_min = signal_region_range[0]
    i_max = signal_region_range[1]

    # 1) Retrieve histograms and normalisations
    f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/histograms.root")
    f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/histogram_errors.root")

    # A
    A = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{bdt}/A.npy')
    A_err = np.load(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{bdt}/A_err.npy')
    A_ufloat = ufloat(A, A_err)
    if(A_ufloat.nominal_value == 0):
        run_fit = False
        a_ufloat = ufloat(-1, 0)
    else:
        a_ufloat = 1/A_ufloat

    # B
    B = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy')
    B_ufloat = ufloat(B, B_err)
    b_ufloat = 1/B_ufloat

    # C = [C_BDDKp, C_BuDDK0, C_BuDD]
    C = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/C.npy')
    C_err = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/C_err.npy')
    C_BDDKp_ufloat = ufloat( C[0], C_err[0] )

    # N_comb (this changes depending on the fit_configuration) 
    combinatoria_yield = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/N_comb.npy')
    N_comb_ufloat = ufloat( combinatoria_yield[0], combinatoria_yield[1] )

    # channel_efficiency = [eps_sig, eps_comb, eps_BDDKp, eps_BuDDK0, eps_BuDD]
    # eps = [eps_0, eps_1, eps_2, eps_3]
    channel_eps = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/channel_eff_value.npy')
    channel_eps_err = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{bdt}/channel_eff_error.npy')
    
    norm_parameters = { "N_comb": N_comb_ufloat.nominal_value, "a": a_ufloat.nominal_value, "b": b_ufloat.nominal_value, "C_BDDKp": C_BDDKp_ufloat.nominal_value, "eps_sig": channel_eps[0], "eps_comb": channel_eps[1], "eps_BDDKp": channel_eps[2] }
    norm_parameters_errors = { "N_comb": N_comb_ufloat.std_dev, "a": a_ufloat.std_dev, "b": b_ufloat.std_dev, "C_BDDKp": C_BDDKp_ufloat.std_dev, "eps_sig": channel_eps_err[0], "eps_comb": channel_eps_err[1], "eps_BDDKp": channel_eps_err[2] }

    # Retrieve template histograms
    if(fit_name == "AllEvents"):
        if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
            histograms, histogram_errors, h_upward, h_downward, h_data = retrieve_histograms(fit_type, f,f1,0)
            h_data = [h_data]
        else:
            histograms, histogram_errors, h_upward, h_downward = retrieve_histograms(fit_type, f,f1,0)

        nominal_templates  = [histograms]
        error_templates    = [histogram_errors]
        h_upward = [h_upward]
        h_downward = [h_downward]
    else:
        if((fit_type == "RSDataSidebands") or (fit_type == "RSData")):
            histograms_1, histogram_errors_1, h_upward_1, h_downward_1, h_data_1 = retrieve_histograms(fit_type, f,f1,1)
            histograms_2, histogram_errors_2, h_upward_2, h_downward_2, h_data_2 = retrieve_histograms(fit_type, f,f1,2)
            histograms_3, histogram_errors_3, h_upward_3, h_downward_3, h_data_3 = retrieve_histograms(fit_type, f,f1,3)

            h_data = [h_data_1, h_data_2, h_data_3]
        else:
            histograms_1, histogram_errors_1, h_upward_1, h_downward_1 = retrieve_histograms(fit_type, f,f1,1)
            histograms_2, histogram_errors_2, h_upward_2, h_downward_2 = retrieve_histograms(fit_type, f,f1,2)
            histograms_3, histogram_errors_3, h_upward_3, h_downward_3 = retrieve_histograms(fit_type, f,f1,3)

        nominal_templates = [histograms_1, histograms_2, histograms_3]
        error_templates = [histogram_errors_1, histogram_errors_2, histogram_errors_3]
        h_upward = [h_upward_1, h_upward_2, h_upward_3]
        h_downward = [h_downward_1, h_downward_2, h_downward_3]

    end = time.time()
    print(f"Elapsed time (retrieving files): {end - start:.2f} seconds")

    ### Do not run the fit in the signal efficiency is 0 (no signal passes the BDT selections)
    if(fit_name == "AllEvents"):
        n_sig = BF_sig*a_ufloat*b_ufloat*ufloat( norm_parameters['eps_sig'][0], norm_parameters_errors['eps_sig'][0])
        n_BDDKp = C_BDDKp_ufloat*b_ufloat*ufloat(norm_parameters['eps_BDDKp'][0], norm_parameters_errors['eps_BDDKp'][0])
        n_comb = ufloat(norm_parameters['N_comb'], norm_parameters_errors['N_comb'])*ufloat(norm_parameters['eps_comb'][0], norm_parameters_errors['eps_comb'][0])
    else:
        n_sig = 0
        n_BDDKp = ufloat(0, 0)
        n_comb = ufloat(0, 0)

        for i in range(len(norm_parameters['eps_sig'])-1):
            n_sig += BF_sig*a_ufloat*b_ufloat*ufloat( norm_parameters['eps_sig'][i+1], norm_parameters_errors['eps_sig'][i+1])
            n_BDDKp +=  C_BDDKp_ufloat*b_ufloat*ufloat( norm_parameters['eps_BDDKp'][i+1], norm_parameters_errors['eps_BDDKp'][i+1] )
            n_comb += ufloat( norm_parameters['N_comb'], norm_parameters_errors['N_comb'])*ufloat( norm_parameters['eps_comb'][i+1], norm_parameters_errors['eps_comb'][i+1] )

    print("Yields generated (all events) : ")
    print("N_sig = ", n_sig.nominal_value, " +/- ", n_sig.std_dev)
    print("N_comb = ", n_comb.nominal_value, " +/- ", n_comb.std_dev)
    print("N_BDDKp = ", n_BDDKp.nominal_value, " +/- ", n_BDDKp.std_dev)
    n_total = n_sig+n_BDDKp+n_comb

    print("N_total = ", n_total.nominal_value, " +/- ", n_total.std_dev)
    # if(int(n_total.nominal_value) < 10):
    #     print("Not running the fit because there are less < 10 events in data")
    #     run_fit = False
    
    if(run_fit):
        ############################################################# Fit validation ########################################################################
        if(validate_fit and not fit_type == "RSDataSidebands"):
            toy_fit_values = []
            toy_fit_errors = []
            N_fail = 0

            for seed in range(n_toys):
                if(n_toys > 10):
                    if seed % (n_toys // 10) == 0:
                        print(f"Progress: {100 * seed // n_toys}% (seed={seed})")

                h_data = generate_toy_data(fit_name, BF_sig, seed, nominal_templates, h_upward, h_downward, norm_parameters, norm_parameters_errors, add_efficiencies)
                
                spec = build_model(fit_name, fit_type, BF_sig, nominal_templates, error_templates, norm_parameters, norm_parameters_errors, h_data, h_upward, h_downward, i_min, i_max, add_efficiencies)

                if(seed == 0):
                    save_plot = True
                else:
                    save_plot = False

                try:
                    fit_pars, fit_errors = do_fit(fit_name, fit_type, BF_sig, bdt, spec, nominal_templates, save_plot, i_min, i_max, add_efficiencies, 0, seed)

                    toy_fit_values.append(fit_pars)
                    toy_fit_errors.append(fit_errors)
                except:
                    N_fail += 1
                    continue

            workspace = pyhf.Workspace(spec)
            model = workspace.model()
            toy_studies(fit_name, fit_type, BF_sig, bdt, model, N_fail, toy_fit_values, toy_fit_errors)
            print(f"{N_fail}/{n_toys} toys fail")

        ############################################################# Extract Ncomb ########################################################################
        elif(fit_type == "RSDataSidebands"):
            spec = build_model(fit_name, fit_type, BF_sig, nominal_templates, error_templates, norm_parameters, norm_parameters_errors, h_data, h_upward, h_downward, i_min, i_max, add_efficiencies)
            try:
                ncomb = do_fit(fit_name, fit_type, BF_sig, bdt, spec, nominal_templates, True, i_min, i_max, add_efficiencies, 0)
                np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/ncomb_value.npy', ncomb.nominal_value)
                np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/ncomb_error.npy', ncomb.std_dev)
            except:
                print("CLs calculation failed")
                save_dummy(validate_fit, fit_type, fit_name, BF_sig, bdt)
        
        ############################################################# CLs calculation ########################################################################
        else:
            h_data_nominal = generate_cls_data(fit_name, fit_type, BF_sig, nominal_templates, h_upward, h_downward, norm_parameters, norm_parameters_errors, add_efficiencies, 0)
            h_data_upward = generate_cls_data(fit_name, fit_type, BF_sig, nominal_templates, h_upward, h_downward, norm_parameters, norm_parameters_errors, add_efficiencies, 1)
            h_data_downward = generate_cls_data(fit_name, fit_type, BF_sig, nominal_templates, h_upward, h_downward, norm_parameters, norm_parameters_errors, add_efficiencies, -1)

            spec_nominal = build_model(fit_name, fit_type, BF_sig, nominal_templates, error_templates, norm_parameters, norm_parameters_errors, h_data_nominal, h_upward, h_downward, i_min, i_max, add_efficiencies)
            spec_upward = build_model(fit_name, fit_type, BF_sig, nominal_templates, error_templates, norm_parameters, norm_parameters_errors, h_data_upward, h_upward, h_downward, i_min, i_max, add_efficiencies)
            spec_downward = build_model(fit_name, fit_type, BF_sig, nominal_templates, error_templates, norm_parameters, norm_parameters_errors, h_data_downward, h_upward, h_downward, i_min, i_max, add_efficiencies)

            # try:
            print("NOMINAL")
            exp_limit_nominal, obs_limit_nominal = do_fit(fit_name, fit_type, BF_sig, bdt, spec_nominal, nominal_templates, True, i_min, i_max, add_efficiencies, 0)
            print("NOMINAL + ERROR")
            exp_limit_upward, obs_limit_upward = do_fit(fit_name, fit_type, BF_sig, bdt, spec_upward, nominal_templates, True, i_min, i_max, add_efficiencies, 1)
            print("NOMINAL - ERROR")
            exp_limit_downward, obs_limit_downward = do_fit(fit_name, fit_type, BF_sig, bdt, spec_downward, nominal_templates, True, i_min, i_max, add_efficiencies, -1)

            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit.npy', exp_limit_nominal)
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_obs_limit.npy', obs_limit_nominal)
        
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_up_comb_yield_syst.npy', exp_limit_upward)
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_obs_limit_up_comb_yield_syst.npy', obs_limit_upward)

            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_limit_down_comb_yield_syst.npy', exp_limit_downward)
            np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig:.5g}/BDT_{bdt}/cls_obs_limit_down_comb_yield_syst.npy', obs_limit_downward)
        
            # except:
            #     print("CLs calculation failed")
            #     save_dummy(validate_fit, fit_type, fit_name, BF_sig, bdt)

    else:
        save_dummy(validate_fit, fit_type, fit_name, BF_sig, bdt)

    ############################################################################################################################################################

    end1 = time.time()
    print(f"Elapsed time : {end1 - end:.2f} seconds")


if __name__ == "__main__":
    main(sys.argv)




