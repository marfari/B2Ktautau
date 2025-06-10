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

add_physics_backgrounds = True
add_statistical_error = True
validate_fit = True

n_toys = 1000
if(n_toys > 1000):
    validate_fit = False

def retrieve_histograms(f, f1, ch):
    h_sig = f.Get(f"Channel_{ch}/h_sig_{ch}")
    h_comb = f.Get(f"Channel_{ch}/h_comb_{ch}")
    h_BDDKp = f.Get(f"Channel_{ch}/h_BDDKp_{ch}")

    h_sig_err = f1.Get(f"Channel_{ch}/h_sig_err_{ch}")
    h_comb_err = f1.Get(f"Channel_{ch}/h_comb_err_{ch}")
    h_BDDKp_err = f1.Get(f"Channel_{ch}/h_BDDKp_err_{ch}")

    histograms = [h_sig, h_comb, h_BDDKp]
    histogram_errors = [h_sig_err, h_comb_err, h_BDDKp_err]

    return histograms, histogram_errors


def generate_toy_data(histograms, yield_values, yield_errors, seed, nsig, add_physics_backgrounds, ch):
    h_sig = histograms[0]  
    h_comb = histograms[1]  
    h_BDDKp = histograms[2]

    nbins = h_sig.GetNbinsX()

    N_comb = yield_values[ch][0]
    if(add_physics_backgrounds):
        N_BDDKp = yield_values[0][1]
        N_BDDKp_err = yield_errors[0][1]

        epsilon = yield_values[ch][1]/N_BDDKp # currently this is fixed in the fit so I don't need to vary it across toys (but later need to sample it from a Gaussian)

    # varied templates
    h_comb_varied = ROOT.TH1D(f"h_comb_varied_{seed}_{ch}", f"h_comb_varied_{seed}_{ch}", nbins, 4000, 8000)
    if(nsig != 0):
        h_sig_varied = ROOT.TH1D(f"h_sig_varied_{seed}_{ch}", f"h_sig_varied_{seed}_{ch}", nbins, 4000, 8000)
    if(add_physics_backgrounds):
        h_BDDKp_varied = ROOT.TH1D(f"h_BDDKp_varied_{seed}_{ch}", f"h_BDDKp_varied_{seed}_{ch}", nbins, 4000, 8000)

    # data
    h_data = ROOT.TH1D(f"h_data_{seed}_{ch}", f"h_data_{seed}_{ch}", nbins, 4000, 8000)
    h_data_comb = ROOT.TH1D(f"h_data_comb_{seed}_{ch}", f"h_data_comb_{seed}_{ch}", nbins, 4000, 8000)
    if(nsig != 0):
        h_data_sig = ROOT.TH1D(f"h_data_sig_{seed}_{ch}", f"h_data_sig_{seed}_{ch}", nbins, 4000, 8000)
    if(add_physics_backgrounds):
        h_data_BDDKp = ROOT.TH1D(f"h_data_BDDKp_{seed}_{ch}", f"h_data_BDDKp_{seed}_{ch}", nbins, 4000, 8000)

    np.random.seed(seed)
    if(validate_fit):
        if(add_physics_backgrounds):
            N_BDDKp = max(0, np.random.normal(N_BDDKp, N_BDDKp_err) ) # this is the total N, it is sampled from a Gaussian distribution because a Gaussian constraint is applied in the fit
            N_BDDKp *= epsilon # multiply it by the error category efficiency (in the fit to all events this is 1); error category efficiency is fixed across toys
 
        # vary the templates
        toy_counts = np.array([np.random.poisson( h_comb.GetBinContent(i+1)*h_comb.GetEntries() ) for i in range(nbins)])
        for i, count in enumerate(toy_counts):
            h_comb_varied.SetBinContent(i+1, count)
        h_comb_varied.Sumw2()

        if(nsig != 0):
            toy_counts_sig = np.array([np.random.poisson( h_sig.GetBinContent(i+1)*h_sig.GetEntries() ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_sig):
                h_sig_varied.SetBinContent(i+1, count)
            h_sig_varied.Sumw2()

        if(add_physics_backgrounds):
            toy_counts_BDDKp = np.array([np.random.poisson( h_BDDKp.GetBinContent(i+1)*h_BDDKp.GetEntries() ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BDDKp):
                h_BDDKp_varied.SetBinContent(i+1, count)
            h_BDDKp_varied.Sumw2()

        # generate toy data from the varied templates
        toy_counts_1 = np.array([np.random.poisson( h_comb_varied.GetBinContent(i+1)*(N_comb/h_comb_varied.Integral()) ) for i in range(nbins)])
        for i, count in enumerate(toy_counts_1):
            h_data_comb.SetBinContent(i+1, count)
        h_data_comb.Sumw2()

        if(nsig != 0):
            toy_counts_sig_1 = np.array([np.random.poisson( h_sig_varied.GetBinContent(i+1)*(nsig/h_sig_varied.Integral()) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_sig_1):
                h_data_sig.SetBinContent(i+1, count)
            h_data_sig.Sumw2()

        if(add_physics_backgrounds):
            toy_counts_BDDKp_1 = np.array([np.random.poisson( h_BDDKp_varied.GetBinContent(i+1)*( N_BDDKp/h_BDDKp_varied.Integral()) ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BDDKp_1):
                h_data_BDDKp.SetBinContent(i+1, count)
            h_data_BDDKp.Sumw2()

    else:
        h_data_comb = h_comb.Clone("h_data_comb")
        h_data_comb.Scale(N_comb)

        if(nsig != 0):
            h_data_sig = h_sig.Clone("h_data_sig")
            h_data_sig.Scale(nsig)
        if(add_physics_backgrounds):
            h_data_BDDKp = h_BDDKp.Clone("h_data_BDDKp")
            h_data_BDDKp.Scale(N_BDDKp)

    h_data.Sumw2()
    h_data.Add(h_data_comb)
    if(nsig != 0):
        h_data.Add(h_data_sig)
    if(add_physics_backgrounds):
        h_data.Add(h_data_BDDKp)

    return h_data


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


def convert_toy_data_histogram_to_array(histo):
    nbins = histo.GetNbinsX()
    return [histo.GetBinContent(i+1) for i in range(nbins)]


def write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yield_values, yield_errors, c_value, c_error, add_physics_backgrounds):
    if(fit_name == "all_events"):
        n_channels = 1
    else:
        n_channels = 3

    N_BDDKp = yield_values[0][1]
    N_BDDKp_err = yield_errors[0][1]
    c = c_value[0]
    c_err = c_error[0]

    # Compute values for Gaussian constraints (these are common among the error categories)
    # Signal
    if(c == 0):
        upward = 1.999
        downward = 0.001
    else:
        upward = 1+(c_err/c)
        downward = 1-(c_err/c)

    # B -> D D K+
    if(N_BDDKp == 0):
        upward_BDDKp = 1.999
        downward_BDDKp = 0.001
    else:
        upward_BDDKp = 1+(N_BDDKp_err/N_BDDKp)
        downward_BDDKp = 1-(N_BDDKp_err/N_BDDKp)

    channel_samples = []
    parameters = [{"name": "lumi", "auxdata": [1.0], "bounds": [[0.5,1.5]], "fixed": True, "inits": [1.0], "sigmas": [0.1]}, {"name": "BF_sig", "bounds": [[-1.0,1.0]], "inits": [nsig/c]}, {"name": "c_bdt1_bdt2", "bounds": [[0.0,1000000000.0]], "fixed": True, "inits": [c]}]
    
    if(add_physics_backgrounds and (N_BDDKp != 0) and (upward_BDDKp < 2) and (downward_BDDKp > 0)):
        parameters.append({"bounds": [[-10*(N_BDDKp+1),10*(N_BDDKp+1)]], "fixed": True, "inits": [N_BDDKp], "name": "N_BDDKp"})

    for i in range(n_channels):
        data_sig = histo_data[i][0]
        data_comb = histo_data[i][1]
        data_BDDKp = histo_data[i][2]

        data_sig_err = histo_data_err[i][0]
        data_comb_err = histo_data_err[i][1]
        data_BDDKp_err = histo_data_err[i][2]

        if(fit_name == "all_events"):
            ch = i
        else:
            ch = i+1

        # Modifiers
        sig_modifiers = [{"data": None, "name": "c_bdt1_bdt2", "type": "normfactor"}, {"data": None, "name": "BF_sig", "type": "normfactor"}, {"data": {"hi": upward, "lo": downward}, "name": "c_bdt1_bdt2_syst", "type": "normsys"}]
        comb_modifiers = [{"data": None, "name": f"N_comb_{ch}", "type": "normfactor"}]   
        BDDKp_modifiers = [{"data": None, "name": "N_BDDKp", "type": "normfactor"}, {"data": {"hi": upward_BDDKp, "lo": downward_BDDKp}, "name": "N_BDDKp_syst", "type": "normsys"}]                 
        
        if(fit_name == "error_categories"):
            sig_modifiers.append({"data": None, "name": f"c_epsilon_{ch}", "type": "normfactor"})
            BDDKp_modifiers.append({"data": None, "name": f"N_BDDKp_epsilon_{ch}", "type": "normfactor"})
            parameters.append({"name": f"c_epsilon_{ch}", "bounds": [[0.0,2.0]], "fixed": True, "inits": [c_value[ch]/c]})

        if(add_statistical_error):
            sig_modifiers.append({"name": f"shapesys_sig_{ch}", "data": data_sig_err, "type": "shapesys"})
            comb_modifiers.append({"name": f"shapesys_comb_{ch}", "data": data_comb_err, "type": "shapesys"})
            BDDKp_modifiers.append({"name": f"shapesys_BDDKp_{ch}", "data": data_BDDKp_err, "type": "shapesys"})

        # Samples & parameters
        samples = [{"name": "Signal", "data": data_sig, "modifiers": sig_modifiers}, {"name": "Combinatorial", "data": data_comb, "modifiers": comb_modifiers}]
        N_comb = yield_values[ch][0]
        parameters.append({"name": f"N_comb_{ch}", "bounds": [[-10*(N_comb+1),10*(N_comb+1)]], "inits": [N_comb]})

        if(add_physics_backgrounds and (N_BDDKp != 0) and (upward_BDDKp < 2) and (downward_BDDKp > 0)):
            samples.append({"name": "BDDKp", "data": data_BDDKp, "modifiers": BDDKp_modifiers})
            if(fit_name == "error_categories"):
                parameters.append({"name": f"N_BDDKp_epsilon_{ch}", "bounds": [[0.0,2.0]], "fixed": True, "inits": [yield_values[ch][1]/N_BDDKp]})

        channel_samples.append(samples)

    if(fit_name == "all_events"):
        channels = [{"name": "AllEvents", "samples": channel_samples[0]}]
        observations = [{"name": "AllEvents", "data": histo_toy_data[0]}]
    else:
        channels = [{"name": "Channel_1", "samples": channel_samples[0]}, {"name": "Channel_2", "samples": channel_samples[1]}, {"name": "Channel_3", "samples": channel_samples[2]}]
        observations = [{"name": "Channel_1", "data": histo_toy_data[0]}, {"name": "Channel_2", "data": histo_toy_data[1]}, {"name": "Channel_3", "data": histo_toy_data[2]}]

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


def plot(fit_name, h_data, workspace, fit_pars, fit_errors, nsig, bdt1, bdt2, seed):
    model = workspace.model()

    channel_nbins = model.config.channel_nbins
    if(fit_name == "all_events"):
        nbins = channel_nbins["AllEvents"]
        mass = [h_data[0].GetBinCenter(i+1) for i in range(nbins)] 
        edges = [h_data[0].GetBinLowEdge(i+1) for i in range(nbins+1)] 
    else:
        nbins_1 = channel_nbins["Channel_1"]
        nbins_2 = channel_nbins["Channel_2"]
        nbins_3 = channel_nbins["Channel_3"]

        mass_1 = [h_data[0].GetBinCenter(i+1) for i in range(nbins_1)] 
        mass_2 = [h_data[1].GetBinCenter(i+1) for i in range(nbins_2)] 
        mass_3 = [h_data[2].GetBinCenter(i+1) for i in range(nbins_3)] 

        edges_1 = [h_data[0].GetBinLowEdge(i+1) for i in range(nbins_1+1)] 
        edges_2 = [h_data[1].GetBinLowEdge(i+1) for i in range(nbins_2+1)] 
        edges_3 = [h_data[2].GetBinLowEdge(i+1) for i in range(nbins_3+1)] 

    sig_index = model.config.par_names.index("BF_sig")
    fit_poi = fit_pars[sig_index] 
    fit_poi_error = fit_errors[sig_index] 
    print("BF_sig = ", fit_poi, " +/- ", fit_poi_error)

    c_index = model.config.par_names.index("c_bdt1_bdt2")
    
    try:
        BDDKp_index = model.config.par_names.index("N_BDDKp")
    except:
        BDDKp_index = -1

    if(fit_name == "all_events"):
        comb_index = model.config.par_names.index("N_comb_0")

        ncomb = fit_pars[comb_index]
        ncomb_err = fit_errors[comb_index]
        print("N_comb = ", ncomb, " +/- ", ncomb_err)

    else:
        comb_index_1 = model.config.par_names.index("N_comb_1")
        comb_index_2 = model.config.par_names.index("N_comb_2")
        comb_index_3 = model.config.par_names.index("N_comb_3")

        nbkg_1 = fit_pars[comb_index_1]
        nbkg_2 = fit_pars[comb_index_2]
        nbkg_3 = fit_pars[comb_index_3]

        nbkg_1_error = fit_errors[comb_index_1]
        nbkg_2_error = fit_errors[comb_index_2]
        nbkg_3_error = fit_errors[comb_index_3]

        print("N_bkg_1  = ", nbkg_1, " +/- ", nbkg_1_error)
        print("N_bkg_2  = ", nbkg_2, " +/- ", nbkg_2_error)
        print("N_bkg_3  = ", nbkg_3, " +/- ", nbkg_3_error)

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
    if(BDDKp_index >= 0):
        data_BDDKp = mc_counts[model.config.samples.index("BDDKp")][0]

    if(fit_name == "all_events"):
        all_data = data_sig+data_comb
        if(BDDKp_index >= 0):
            all_data += data_BDDKp

        f, ax = plt.subplots()
        plt.step(edges, np.append(all_data, all_data[-1]), where='post', color='black', label="Total fit")
        plt.step(edges, np.append(data_sig, data_sig[-1]), where='post', color='blue', label='Signal')
        plt.step(edges, np.append(data_comb, data_comb[-1]), where='post', color='red', label='Combinatorial')
        if(BDDKp_index >= 0):
            plt.step(edges, np.append(data_BDDKp, data_BDDKp[-1]), where='post', color='green', label='$B \\to D D K^+$')

        the_data = workspace.data(model, include_auxdata=False)
        the_data_err = np.sqrt(the_data)

        plt.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
        plt.title(f'All events \n BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed}')
        plt.xlabel('m_B (MeV)')
        plt.ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins, 1) ))
        ax.set_yscale("symlog")
        plt.xlim(4000,8000)

        param1 = ufloat(fit_poi, fit_poi_error)
        param2 = ufloat(ncomb, ncomb_err) 

        textstr = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1}".format(param1, param2)
        plt.text(0.3, 0.8, textstr, fontsize=10, transform=ax.transAxes)
        plt.tight_layout()
        plt.legend()

    else:
        data_sig_1 = data_sig[0:nbins_1]
        data_sig_2 = data_sig[nbins_1:nbins_1+nbins_2]
        data_sig_3 = data_sig[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

        data_bkg_1 = data_comb[0:nbins_1]
        data_bkg_2 = data_comb[nbins_1:nbins_1+nbins_2]
        data_bkg_3 = data_comb[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]
    
        if(BDDKp_index >= 0):
            data_BDDKp_1 = data_BDDKp[0:nbins_1]
            data_BDDKp_2 = data_BDDKp[nbins_1:nbins_1+nbins_2]
            data_BDDKp_3 = data_BDDKp[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

        all_data_1 = data_sig_1+data_bkg_1
        all_data_2 = data_sig_2+data_bkg_2
        all_data_3 = data_sig_3+data_bkg_3
        if(BDDKp_index >= 0):
            all_data_1 += data_BDDKp_1
            all_data_2 += data_BDDKp_2
            all_data_3 += data_BDDKp_3

        the_data = workspace.data(model, include_auxdata=False)
        the_data_err = np.sqrt(the_data)

        the_data_1 = the_data[0:nbins_1]
        the_data_2 = the_data[nbins_1:nbins_1+nbins_2]
        the_data_3 = the_data[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

        the_data_err_1 = the_data_err[0:nbins_1]
        the_data_err_2 = the_data_err[nbins_1:nbins_1+nbins_2]
        the_data_err_3 = the_data_err[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

        figure, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
        ax1.step(edges_1, np.append(all_data_1, all_data_1[-1]), where='post', color='black', label="Total fit")
        ax1.step(edges_1, np.append(data_sig_1, data_sig_1[-1]), where='post', color='blue', label='Signal')
        ax1.step(edges_1, np.append(data_bkg_1, data_bkg_1[-1]), where='post', color='red', label='Combinatorial')
        if(BDDKp_index >= 0):
            ax1.step(edges_1, np.append(data_BDDKp_1, data_BDDKp_1[-1]), where='post', color='green', label='$B \\to D D K^+$')
        ax1.errorbar(mass_1, the_data_1, the_data_err_1, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax2.step(edges_2, np.append(all_data_2, all_data_2[-1]), where='post', color='black', label="Total fit")
        ax2.step(edges_2, np.append(data_sig_2, data_sig_2[-1]), where='post', color='blue', label='Signal')
        ax2.step(edges_2, np.append(data_bkg_2, data_bkg_2[-1]), where='post', color='red', label='Combinatorial')
        if(BDDKp_index >= 0):
            ax2.step(edges_2, np.append(data_BDDKp_2, data_BDDKp_2[-1]), where='post', color='green', label='$B \\to D D K^+$')
        ax2.errorbar(mass_2, the_data_2, the_data_err_2, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax3.step(edges_3, np.append(all_data_3, all_data_3[-1]), where='post', color='black', label="Total fit")
        ax3.step(edges_3, np.append(data_sig_3, data_sig_3[-1]), where='post', color='blue', label='Signal')
        ax3.step(edges_3, np.append(data_bkg_3, data_bkg_3[-1]), where='post', color='red', label='Combinatorial')
        if(BDDKp_index >= 0):
            ax3.step(edges_3, np.append(data_BDDKp_3, data_BDDKp_3[-1]), where='post', color='green', label='$B \\to D D K^+$')
        ax3.errorbar(mass_3, the_data_3, the_data_err_3, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax1.set_title(f'Channel 1 \n BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed}')
        ax2.set_title(f'Channel 2 \n BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed}')
        ax3.set_title(f'Channel 3 \n BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed}')

        ax1.set_xlabel('m_B (MeV)')
        ax1.set_ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins_1, 1) ))

        ax2.set_xlabel('m_B (MeV)')
        ax2.set_ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins_2, 1) ))

        ax3.set_xlabel('m_B (MeV)')
        ax3.set_ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins_3, 1) ))

        ax1.set_yscale("symlog")
        ax2.set_yscale("symlog")
        ax3.set_yscale("symlog")

        the_poi = ufloat(fit_poi, fit_poi_error)
        the_nbkg_1 = ufloat(nbkg_1, nbkg_1_error) 
        the_nbkg_2 = ufloat(nbkg_2, nbkg_2_error) 
        the_nbkg_3 = ufloat(nbkg_3, nbkg_3_error) 

        textstr_1 = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1}".format(the_poi, the_nbkg_1)
        textstr_2 = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1}".format(the_poi, the_nbkg_2)
        textstr_3 = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1}".format(the_poi, the_nbkg_3)

        ax1.text(0.4, 0.8, textstr_1, fontsize=12, transform=ax1.transAxes)
        ax2.text(0.4, 0.8, textstr_2, fontsize=12, transform=ax2.transAxes)
        ax3.text(0.4, 0.8, textstr_3, fontsize=12, transform=ax3.transAxes)

        ax1.legend()
        ax2.legend()
        ax3.legend()

    plt.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_name}/N_sig_{nsig}/fit_plots/fit_plot_bdt1_{bdt1}_bdt2_{bdt2}_seed_{seed}.pdf')
    plt.clf()


def run_cls_limit(fit_poi, fit_poi_error, data, model, fit_name, nsig, bdt1, bdt2):
    print("### CLS LIMIT CALCULATION")
    # CLs limit: evaluate upper limit on parameter of interest
    poi_values = np.linspace(-np.abs(fit_poi), np.abs(fit_poi)+4*fit_poi_error, 50)
    obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q")
    print(f"Upper limit (obs): μ = {obs_limit:.6f}")
    print(f"Upper limit (exp): μ = {exp_limits[2]:.6f}")

    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    ax.set_title(f"Fit {fit_name} \n BDT1 = {bdt1} | BDT2 = {bdt2}")
    ax.set_xlabel(r"$BF_{sig}$")
    ax.set_ylabel(r"$\mathrm{CL}_{s}$")
    brazil.plot_results(scan, results, test_size=0.1, ax=ax)
    textstr1 = f"Upper limit (exp): = {exp_limits[2]:.6f}"
    plt.text(0.6, 0.6, textstr1, fontsize=15, transform=ax.transAxes)
    fig.savefig(f"/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_name}/N_sig_{nsig}/fit_plots/cls_limit_bdt1_{bdt1}_bdt2_{bdt2}.pdf")
    plt.clf()

    return exp_limits[2]


def fit(fit_name, nsig, bdt1, bdt2, h_data, spec, seed, save_plot, run_cls):

    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    data = workspace.data(model)

    run_fit = True    

    # Fitting
    optimizer = pyhf.optimize.minuit_optimizer(verbose=0)
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

    if(save_plot):
        print(res_obj)
        plot(fit_name, h_data, workspace, fit_pars, fit_errors, nsig, bdt1, bdt2, seed)
    if(run_cls):
        sig_index = model.config.par_names.index("BF_sig")
        fit_poi = fit_pars[sig_index]
        fit_poi_error = fit_errors[sig_index]

        exp_limit = run_cls_limit(fit_poi, fit_poi_error, data, model, fit_name, nsig, bdt1, bdt2)

    if(validate_fit):
        return fit_pars, fit_errors
    else:
        return exp_limit


def toy_studies(fit_name, nsig, bdt1, bdt2, model, N_fail, toy_fit_values, toy_fit_errors, nbins=30):
    sig_index = model.config.par_names.index("BF_sig")
    if(fit_name == "all_events"):
        comb_index = model.config.par_names.index("N_comb_0")
    else:
        comb_index_1 = model.config.par_names.index("N_comb_1")
        comb_index_2 = model.config.par_names.index("N_comb_2")
        comb_index_3 = model.config.par_names.index("N_comb_3")

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
        if(fit_name == "all_events"):
            Nbkg_toys[i] = toy_fit_values[i][comb_index]
            Nbkg_toys_err[i] = toy_fit_errors[i][comb_index]
        else:
            Nbkg_toys_1[i] = toy_fit_values[i][comb_index_1]
            Nbkg_toys_2[i] = toy_fit_values[i][comb_index_2]
            Nbkg_toys_3[i] = toy_fit_values[i][comb_index_3]
            Nbkg_toys_err_1[i] = toy_fit_errors[i][comb_index_1]
            Nbkg_toys_err_2[i] = toy_fit_errors[i][comb_index_2]
            Nbkg_toys_err_3[i] = toy_fit_errors[i][comb_index_3]

    expected_values =  model.config.suggested_init()
    toy_labels = model.config.par_names

    for i in range(len(expected_values)):
        param_values = [toy_fit_values[j][i] for j in range(N)]
        param_errors = [toy_fit_errors[j][i] for j in range(N)]

        if(toy_labels[i] != "c_bdt1_bdt2"):
            for j in range(N):
                if(param_errors[j] == 0):
                    param_errors[j] = 0.01

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
        ax1.set_title(f"Fit to all events ({N_fail}/{n_toys} fits fail) \n BDT1 = {bdt1} | BDT2 = {bdt2}")
        ax1.axvline(0, color='black', linestyle='--', linewidth=2)
        ax1.text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu, sigma/np.sqrt(N), sigma, sigma/np.sqrt(2*N)),
            transform=ax1.transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
        fig1.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_name}/N_sig_{nsig}/fit_validation_plots/pulls/BDT1_{bdt1}_bdt2_{bdt2}/{toy_labels[i]}.pdf')
        plt.clf()
        
    if(fit_name == "all_events"):
        fig, ax = plt.subplots(3, 2, figsize=(10, 15))
    else:
        fig, ax = plt.subplots(3, 4, figsize=(20, 15))

    BF_sig = expected_values[sig_index]
    if(fit_name == "all_events"):
        N_comb = expected_values[comb_index]
    else:
        N_comb_1 = expected_values[comb_index_1]
        N_comb_2 = expected_values[comb_index_2]
        N_comb_3 = expected_values[comb_index_3]

    # BF pull
    pull = (BF_toys - BF_sig)/BF_toys_err
    (mu_pull, sigma_pull) = norm.fit(pull)
    n, bins, patches = ax[0,0].hist(pull, bins=nbins)
    xcenters = (bins[:-1] + bins[1:]) / 2
    ax[0,0].errorbar(xcenters, n, yerr=np.sqrt(n), ecolor='black', fmt='k.')
    y_pull = sum(n)*(bins[1] - bins[0])*norm.pdf(xcenters, mu_pull, sigma_pull)
    ax[0,0].plot(xcenters, y_pull, 'r-', linewidth=2)
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
    ax[1,0].hist(BF_toys, bins=nbins)
    ax[1,0].set_xlabel("$BF_{sig}$")
    ax[1,0].set_ylabel(f"Entries / {nbins} bins")
    ax[1,0].set_title("Branching fraction")
    ax[1,0].axvline(BF_sig, color='black', linestyle='--', linewidth=2)
    ax[1,0].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(BF_toys), np.std(BF_toys) ),
            transform=ax[1,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # BF error
    ax[2,0].hist(BF_toys_err, bins=nbins)
    ax[2,0].set_xlabel("$BF_{sig}$ error")
    ax[2,0].set_ylabel(f"Entries / {nbins} bins")
    ax[2,0].set_title("Branching fraction error")
    ax[2,0].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(BF_toys_err), np.std(BF_toys_err) ),
            transform=ax[2,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    if(fit_name == "all_events"):
        pull_2 = (Nbkg_toys - N_comb)/Nbkg_toys_err
        (mu_pull_2, sigma_pull_2) = norm.fit(pull_2)
        n_2, bins_2, patches_2 = ax[0,1].hist( pull_2, bins=nbins)
        xcenters_2 = (bins_2[:-1] + bins_2[1:]) / 2
        ax[0,1].errorbar(xcenters_2, n_2, yerr=np.sqrt(n_2), ecolor='black', fmt='k.')
        y_pull_2 = sum(n_2)*(bins_2[1] - bins_2[0])*norm.pdf( xcenters_2, mu_pull_2, sigma_pull_2)
        ax[0,1].plot(xcenters_2, y_pull_2, 'r-', linewidth=2)
        ax[0,1].set_xlabel("Pull: ($N_{bkg} - N_{bkg}^{expected}$)/error")
        ax[0,1].set_ylabel(f"Entries / {nbins} bins")
        chi2_2 = np.sum( ((n_2 - y_pull_2)**2) / y_pull_2) / (nbins - 2)
        ax[0,1].set_title(f"Background yield pull: $\chi^2$/ndf = {chi2_2:.3g}")
        ax[0,1].axvline(0, color='black', linestyle='--', linewidth=2)
        ax[0,1].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_2, sigma_pull_2/np.sqrt(N), sigma_pull_2, sigma_pull_2/np.sqrt(2*N)),
                transform=ax[0,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[1,1].hist(Nbkg_toys, bins=nbins)
        ax[1,1].set_xlabel("$N_{bkg}$")
        ax[1,1].set_ylabel(f"Entries / {nbins} bins")
        ax[1,1].set_title("Background yield")
        ax[1,1].axvline(N_comb, color='black', linestyle='--', linewidth=2)
        ax[1,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys), np.std(Nbkg_toys) ),
                transform=ax[1,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,1].hist(Nbkg_toys_err, bins=nbins)
        ax[2,1].set_xlabel("$N_{bkg}$ error")
        ax[2,1].set_ylabel(f"Entries / {nbins} bins")
        ax[2,1].set_title("Background yield error")
        ax[2,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_err), np.std(Nbkg_toys_err) ),
                transform=ax[2,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    else:
        pull_2 = (Nbkg_toys_1 - N_comb_1)/Nbkg_toys_err_1
        (mu_pull_2, sigma_pull_2) = norm.fit(pull_2)
        n_2, bins_2, patches_2 = ax[0,1].hist( pull_2, bins=nbins)
        xcenters_2 = (bins_2[:-1] + bins_2[1:]) / 2
        ax[0,1].errorbar(xcenters_2, n_2, yerr=np.sqrt(n_2), ecolor='black', fmt='k.')
        y_pull_2 = sum(n_2)*(bins_2[1] - bins_2[0])*norm.pdf( xcenters_2, mu_pull_2, sigma_pull_2)
        ax[0,1].plot(xcenters_2, y_pull_2, 'r-', linewidth=2)
        ax[0,1].set_xlabel("Pull: ($N_{bkg, 1} - N_{bkg, 1}^{expected}$)/error")
        ax[0,1].set_ylabel(f"Entries / {nbins} bins")
        chi2_2 = np.sum( ((n_2 - y_pull_2)**2) / y_pull_2) / (nbins - 2)
        ax[0,1].set_title(f"Background yield 1 pull: $\chi^2$/ndf = {chi2_2:.3g}")
        ax[0,1].axvline(0, color='black', linestyle='--', linewidth=2)
        ax[0,1].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_2, sigma_pull_2/np.sqrt(N), sigma_pull_2, sigma_pull_2/np.sqrt(2*N)),
                transform=ax[0,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[1,1].hist(Nbkg_toys_1, bins=nbins)
        ax[1,1].set_xlabel("$N_{bkg, 1}$")
        ax[1,1].set_ylabel(f"Entries / {nbins} bins")
        ax[1,1].set_title("Background yield 1")
        ax[1,1].axvline(N_comb_1, color='black', linestyle='--', linewidth=2)
        ax[1,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_1), np.std(Nbkg_toys_1) ),
                transform=ax[1,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,1].hist(Nbkg_toys_err_1, bins=nbins)
        ax[2,1].set_xlabel("$N_{bkg, 1}$ error")
        ax[2,1].set_ylabel(f"Entries / {nbins} bins")
        ax[2,1].set_title("Background yield 1 error")
        ax[2,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_err_1), np.std(Nbkg_toys_err_1) ),
                transform=ax[2,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        pull_3 = (Nbkg_toys_2 - N_comb_2)/Nbkg_toys_err_2
        (mu_pull_3, sigma_pull_3) = norm.fit(pull_3)
        n_3, bins_3, patches_3 = ax[0,2].hist( pull_3, bins=nbins)
        xcenters_3 = (bins_3[:-1] + bins_3[1:]) / 2
        ax[0,2].errorbar(xcenters_3, n_3, yerr=np.sqrt(n_3), ecolor='black', fmt='k.')
        y_pull_3 = sum(n_3)*(bins_3[1] - bins_3[0])*norm.pdf( xcenters_3, mu_pull_3, sigma_pull_3)
        ax[0,2].plot(xcenters_3, y_pull_3, 'r-', linewidth=2)
        ax[0,2].set_xlabel("Pull: ($N_{bkg, 2} - N_{bkg, 2}^{expected}$)/error")
        ax[0,2].set_ylabel(f"Entries / {nbins} bins")
        chi2_3 = np.sum( ((n_3 - y_pull_3)**2) / y_pull_3) / (nbins - 2)
        ax[0,2].set_title(f"Background yield 2 pull: $\chi^2$/ndf = {chi2_3:.3g}")
        ax[0,2].axvline(0, color='black', linestyle='--', linewidth=2)
        ax[0,2].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_3, sigma_pull_3/np.sqrt(N), sigma_pull_3, sigma_pull_3/np.sqrt(2*N)),
                transform=ax[0,2].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[1,2].hist(Nbkg_toys_2, bins=nbins)
        ax[1,2].set_xlabel("$N_{bkg, 2}$")
        ax[1,2].set_ylabel(f"Entries / {nbins} bins")
        ax[1,2].set_title("Background yield 2")
        ax[1,2].axvline(N_comb_2, color='black', linestyle='--', linewidth=2)
        ax[1,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_2), np.std(Nbkg_toys_2) ),
                transform=ax[1,2].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,2].hist(Nbkg_toys_err_2, bins=nbins)
        ax[2,2].set_xlabel("$N_{bkg, 2}$ error")
        ax[2,2].set_ylabel(f"Entries / {nbins} bins")
        ax[2,2].set_title("Background yield 2 error")
        ax[2,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_err_2), np.std(Nbkg_toys_err_2) ),
                transform=ax[2,2].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        pull_4 = (Nbkg_toys_3 - N_comb_3)/Nbkg_toys_err_3
        (mu_pull_4, sigma_pull_4) = norm.fit(pull_4)
        n_4, bins_4, patches_4 = ax[0,3].hist( pull_4, bins=nbins)
        xcenters_4 = (bins_4[:-1] + bins_4[1:]) / 2
        ax[0,3].errorbar(xcenters_4, n_4, yerr=np.sqrt(n_4), ecolor='black', fmt='k.')
        y_pull_4 = sum(n_4)*(bins_4[1] - bins_4[0])*norm.pdf( xcenters_4, mu_pull_4, sigma_pull_4)
        ax[0,3].plot(xcenters_4, y_pull_4, 'r-', linewidth=2)
        ax[0,3].set_xlabel("Pull: ($N_{bkg, 3} - N_{bkg, 3}^{expected}$)/error")
        ax[0,3].set_ylabel("Entries / {0} bins".format(nbins))
        chi2_4 = np.sum( ((n_4 - y_pull_4)**2) / y_pull_4) / (nbins - 2)
        ax[0,3].set_title(f"Background yield 3 pull: $\chi^2$/ndf = {chi2_4:.3g}")
        ax[0,3].axvline(0, color='black', linestyle='--', linewidth=2)
        ax[0,3].text(0.05, 0.95, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_4, sigma_pull_4/np.sqrt(N), sigma_pull_4, sigma_pull_4/np.sqrt(2*N)),
                transform=ax[0,3].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[1,3].hist(Nbkg_toys_3, bins=nbins)
        ax[1,3].set_xlabel("$N_{bkg, 3}$")
        ax[1,3].set_ylabel("Entries / {0} bins".format(nbins))
        ax[1,3].set_title("Background yield 3")
        ax[1,3].axvline(N_comb_3, color='black', linestyle='--', linewidth=2)
        ax[1,3].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_3), np.std(Nbkg_toys_3) ),
                transform=ax[1,3].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,3].hist(Nbkg_toys_err_3, bins=nbins)
        ax[2,3].set_xlabel("$N_{bkg, 3}$ error")
        ax[2,3].set_ylabel("Entries / {0} bins".format(nbins))
        ax[2,3].set_title("Background yield 3 error")
        ax[2,3].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_toys_err_3), np.std(Nbkg_toys_err_3) ),
                transform=ax[2,3].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    if(fit_name == "all_events"):
        fig.suptitle(f"Fit to all events ({N_fail}/{n_toys} fits fail) \n BDT1 = {bdt1} | BDT2 = {bdt2}", fontsize=24)
    else:
        fig.suptitle(f"Fit in error categories ({N_fail}/{n_toys} fits fail) \n BDT1 = {bdt1} | BDT2 = {bdt2}", fontsize=24)
 
    fig.savefig(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_name}/N_sig_{nsig}/fit_validation_plots/bdt1_{bdt1}_bdt2_{bdt2}.pdf')
    plt.clf()

def main(argv):

    fit_name = argv[1]
    nsig = argv[2]
    bdt1 = argv[3]
    bdt2 = argv[4]

    nsig = int(nsig)
    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    start = time.time()

    # 1) Retrieve histograms, constant and yields
    f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/histograms_bdt1_{bdt1}_bdt2_{bdt2}.root")
    f1 = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/histograms_errors_bdt1_{bdt1}_bdt2_{bdt2}.root")

    A = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/A_values_bdt1_{bdt1}_bdt2_{bdt2}.npy')
    A_err = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/A_errors_bdt1_{bdt1}_bdt2_{bdt2}.npy')

    B = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy')

    C = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/C_values_bdt1_{bdt1}_bdt2_{bdt2}.npy', allow_pickle=True)
    C_err = np.load(f'/panfs/felician/B2Ktautau/workflow/generate_histograms/C_errors_bdt1_{bdt1}_bdt2_{bdt2}.npy', allow_pickle=True)

    normalisation_values = [A, B, C]
    normalisation_errors = [A_err, B_err, C_err]
    
    if(fit_name == "all_events"):
        histograms, histogram_errors = retrieve_histograms(f,f1,0)
        histo_data, histo_data_err = convert_histogram_to_array(histograms, histogram_errors)

        histo_data = [histo_data]
        histo_data_err = [histo_data_err]
    else:
        histograms_1, histogram_errors_1 = retrieve_histograms(f,f1,1)
        histograms_2, histogram_errors_2 = retrieve_histograms(f,f1,2)
        histograms_3, histogram_errors_3 = retrieve_histograms(f,f1,3)

        histo_data_1, histo_data_err_1 = convert_histogram_to_array(histograms_1, histogram_errors_1)
        histo_data_2, histo_data_err_2 = convert_histogram_to_array(histograms_2, histogram_errors_2)
        histo_data_3, histo_data_err_3 = convert_histogram_to_array(histograms_3, histogram_errors_3)

        histo_data = [histo_data_1, histo_data_2, histo_data_3]
        histo_data_err = [histo_data_err_1, histo_data_err_2, histo_data_err_3]

    end = time.time()
    print(f"Elapsed time (retrieving files): {end - start:.2f} seconds")

    toy_fit_values = []
    toy_fit_errors = []
    N_fail = 0

    if(validate_fit):
        for seed in range(n_toys):
            if(n_toys > 10):
                if seed % (n_toys // 10) == 0:
                    print(f"Progress: {100 * seed // n_toys}% (seed={seed})")

            if(seed == 0):
                save_plot = True
            else:
                save_plot = False

            # Generate toy data
            if(fit_name == "all_events"):
                h_data = generate_toy_data(histograms, normalisation_values, normalisation_errors, seed, nsig, add_physics_backgrounds, 0)
                histo_toy_data = convert_toy_data_histogram_to_array(h_data)

                h_data = [h_data]
                histo_toy_data = [histo_toy_data]
            else:
                h_data_1 = generate_toy_data(histograms_1, normalisation_values, normalisation_errors, seed, nsig, add_physics_backgrounds, 1)
                h_data_2 = generate_toy_data(histograms_2, normalisation_values, normalisation_errors, seed, nsig, add_physics_backgrounds, 2)
                h_data_3 = generate_toy_data(histograms_3, normalisation_values, normalisation_errors, seed, nsig, add_physics_backgrounds, 3)

                histo_toy_data_1 = convert_toy_data_histogram_to_array(h_data_1)
                histo_toy_data_2 = convert_toy_data_histogram_to_array(h_data_2)
                histo_toy_data_3 = convert_toy_data_histogram_to_array(h_data_3)

                h_data = [h_data_1, h_data_2, h_data_3]
                histo_toy_data = [histo_toy_data_1, histo_toy_data_2, histo_toy_data_3]
            
            spec = write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yield_values, yield_errors, c_value, c_error, add_physics_backgrounds)

            try:
                fit_pars, fit_errors = fit(fit_name, nsig, bdt1, bdt2, h_data, spec, seed, save_plot, False)

                toy_fit_values.append(fit_pars)
                toy_fit_errors.append(fit_errors)
            except:
                N_fail += 1
                continue
            
        workspace = pyhf.Workspace(spec)
        model = workspace.model()
        toy_studies(fit_name, nsig, bdt1, bdt2, model, N_fail, toy_fit_values, toy_fit_errors)
        print(f"{N_fail}/{n_toys} toys fail")

    else: # do not vary the data, run the CLs calculation
        # Generate toy data
        if(fit_name == "all_events"):
            h_data = generate_toy_data(histograms, yield_values, yield_errors, 1000, nsig, add_physics_backgrounds, 0)
            histo_toy_data = convert_toy_data_histogram_to_array(h_data)

            h_data = [h_data]
            histo_toy_data = [histo_toy_data]
        else:
            h_data_1 = generate_toy_data(histograms_1, yield_values, yield_errors, 1000, nsig, add_physics_backgrounds, 1)
            h_data_2 = generate_toy_data(histograms_2, yield_values, yield_errors, 1000, nsig, add_physics_backgrounds, 2)
            h_data_3 = generate_toy_data(histograms_3, yield_values, yield_errors, 1000, nsig, add_physics_backgrounds, 3)

            histo_toy_data_1 = convert_toy_data_histogram_to_array(h_data_1)
            histo_toy_data_2 = convert_toy_data_histogram_to_array(h_data_2)
            histo_toy_data_3 = convert_toy_data_histogram_to_array(h_data_3)

            h_data = [h_data_1, h_data_2, h_data_3]
            histo_toy_data = [histo_toy_data_1, histo_toy_data_2, histo_toy_data_3]
        
        spec = write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yield_values, yield_errors, c_value, c_error, add_physics_backgrounds)

        exp_limit = fit(fit_name, nsig, bdt1, bdt2, h_data, spec, 1000, True, True)
        np.save(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_name}/N_sig_{nsig}/fit_results/cls_limit_bdt1_{bdt1}_bdt2_{bdt2}.npy', exp_limit)

    end1 = time.time()
    print(f"Elapsed time (toy data generation): {end1 - end:.2f} seconds")


if __name__ == "__main__":
    main(sys.argv)