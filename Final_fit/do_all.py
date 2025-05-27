import sys
import ROOT
import numpy as np
import time
import pyhf
import matplotlib.pyplot as plt
from scipy.stats import norm
from uncertainties import ufloat

add_physics_backgrounds = False
add_statistical_error = True
vary_fit_templates = False
validate_fit = True
if(validate_fit):
    vary_fit_templates = True

n_toys = 1000

def retrieve_histograms(f, f1, ch):
    h_sig = f.Get(f"Channel_{ch}/h_sig_{ch}")
    h_comb = f.Get(f"Channel_{ch}/h_comb_{ch}")
    h_BDDKp = f.Get(f"Channel_{ch}/h_BDDKp_{ch}")
    h_BuDDK0 = f.Get(f"Channel_{ch}/h_BuDDK0_{ch}")
    h_BuDD = f.Get(f"Channel_{ch}/h_BuDD_{ch}")

    h_sig_err = f1.Get(f"Channel_{ch}/h_sig_err_{ch}")
    h_comb_err = f1.Get(f"Channel_{ch}/h_comb_err_{ch}")
    h_BDDKp_err = f1.Get(f"Channel_{ch}/h_BDDKp_err_{ch}")
    h_BuDDK0_err = f1.Get(f"Channel_{ch}/h_BuDDK0_err_{ch}")
    h_BuDD_err = f1.Get(f"Channel_{ch}/h_BuDD_err_{ch}")

    histograms = [h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD]
    histogram_errors = [h_sig_err, h_comb_err, h_BDDKp_err, h_BuDDK0_err, h_BuDD_err]

    return histograms, histogram_errors


def generate_toy_data(histograms, yields, seed, nsig, add_physics_backgrounds):
    h_sig = histograms[0]
    h_comb = histograms[1]
    h_BDDKp = histograms[2]
    h_BuDDK0 = histograms[3]
    h_BuDD = histograms[4]

    nbins = h_sig.GetNbinsX()

    N_comb = yields[0]
    if(add_physics_backgrounds):
        N_BDDKp = yields[1]
        N_BuDDK0 = yields[2]
        N_BuDD = yields[3]
    
    h_data = ROOT.TH1D(f"h_data_{seed}", f"h_data_{seed}", nbins, 4000, 8000)
    h_data_comb = ROOT.TH1D(f"h_data_comb_{seed}", f"h_data_comb_{seed}", nbins, 4000, 8000)
    if(nsig != 0):
        h_data_sig = ROOT.TH1D(f"h_data_sig_{seed}", f"h_data_sig_{seed}", nbins, 4000, 8000)
    if(add_physics_backgrounds):
        h_data_BDDKp = ROOT.TH1D(f"h_data_BDDKp_{seed}", f"h_data_BDDKp_{seed}", nbins, 4000, 8000)
        h_data_BuDDK0 = ROOT.TH1D(f"h_data_BuDDK0_{seed}", f"h_data_BuDDK0_{seed}", nbins, 4000, 8000)
        h_data_BuDD = ROOT.TH1D(f"h_data_BuDD_{seed}", f"h_data_BuDD_{seed}", nbins, 4000, 8000)

    np.random.seed(seed)
    # if(validate_fit):
    #     N_comb = np.random.poisson(N_comb)
    #     if(add_physics_backgrounds):
    #         N_BDDKp = np.random.poisson(N_BDDKp)
    #         N_BuDDK0 = np.random.poisson(N_BuDDK0)
    #         N_BuDD = np.random.poisson(N_BuDD)

    if(validate_fit):
        toy_counts = np.array([np.random.poisson( h_comb.GetBinContent(i+1)*N_comb ) for i in range(nbins)])
        for i, count in enumerate(toy_counts):
            h_data_comb.SetBinContent(i+1, count)

        if(nsig != 0):
            toy_counts_sig = np.array([np.random.poisson( h_sig.GetBinContent(i+1)*nsig ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_sig):
                h_data_sig.SetBinContent(i+1, count)

        if(add_physics_backgrounds):
            toy_counts_BDDKp = np.array([np.random.poisson( h_BDDKp.GetBinContent(i+1)*N_BDDKp ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BDDKp):
                h_data_BDDKp.SetBinContent(i+1, count)

            toy_counts_BuDDK0 = np.array([np.random.poisson( h_BuDDK0.GetBinContent(i+1)*N_BuDDK0 ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BuDDK0):
                h_data_BuDDK0.SetBinContent(i+1, count)

            toy_counts_BuDD = np.array([np.random.poisson( h_BuDD.GetBinContent(i+1)*N_BuDD ) for i in range(nbins)])
            for i, count in enumerate(toy_counts_BuDD):
                h_data_BuDD.SetBinContent(i+1, count)

        # ROOT.gRandom.SetSeed(0)
        # rnd = ROOT.TRandom()
        # rnd.SetSeed(seed+1)

        # h_data_comb.FillRandom(h_comb, int(N_comb), rnd)
        # if(add_physics_backgrounds):
        #     h_data_BDDKp.FillRandom(h_BDDKp, int(N_BDDKp), rnd)
        #     h_data_BuDDK0.FillRandom(h_BuDDK0, int(N_BuDDK0), rnd)
        #     h_data_BuDD.FillRandom(h_BuDD, int(N_BuDD), rnd)

        h_data.Sumw2()
        h_data_comb.Sumw2()
        if(nsig != 0):
            h_data_sig.Sumw2()
        if(add_physics_backgrounds):
            h_data_BDDKp.Sumw2()
            h_data_BuDDK0.Sumw2()
            h_data_BuDD.Sumw2()

        h_data.Add(h_data_comb)
        if(nsig != 0):
            h_data.Add(h_data_sig)
        if(add_physics_backgrounds):
            h_data.Add(h_data_BDDKp)
            h_data.Add(h_data_BuDDK0)
            h_data.Add(h_data_BuDD)

    else:
        h_data_comb = h_comb.Clone("h_data_comb")
        if(nsig != 0):
            h_data_sig = h_sig.Clone("h_data_sig")
        if(add_physics_backgrounds):
            h_data_BDDKp = h_BDDKp.Clone("h_data_BDDKp")
            h_data_BuDDK0 = h_BuDDK0.Clone("h_data_BuDDK0")
            h_data_BuDD = h_BuDD.Clone("h_data_BuDD")

        h_data.Sumw2()
        h_data_comb.Sumw2()
        if(nsig != 0):
            h_data_sig.Sumw2()
        if(add_physics_backgrounds):
            h_data_BDDKp.Sumw2()
            h_data_BuDDK0.Sumw2()
            h_data_BuDD.Sumw2()

        h_data_comb.Scale(N_comb)
        if(nsig != 0):
            h_data_sig.Scale(nsig)
        if(add_physics_backgrounds):
            h_data_BDDKp.Scale(N_BDDKp)
            h_data_BuDDK0.Scale(N_BuDDK0)
            h_data_BuDD.Scale(N_BuDD)

        h_data.Add(h_data_comb)
        if(nsig != 0):
            h_data.Add(h_data_sig)
        if(add_physics_backgrounds):
            h_data.Add(h_data_BDDKp)
            h_data.Add(h_data_BuDDK0)
            h_data.Add(h_data_BuDD)

    return h_data


def vary_templates(histograms, seed):
    h_sig = histograms[0]  
    h_comb = histograms[1]  
    h_BDDKp = histograms[2]
    h_BuDDK0 = histograms[3]
    h_BuDD = histograms[4]

    nbins = h_sig.GetNbinsX()

    h_sig_varied = ROOT.TH1D(f"h_sig_varied_{seed}", f"h_sig_varied_{seed}", nbins, 4000, 8000)
    h_comb_varied = ROOT.TH1D(f"h_comb_varied_{seed}", f"h_comb_varied_{seed}", nbins, 4000, 8000)
    h_BDDKp_varied = ROOT.TH1D(f"h_BDDKp_varied_{seed}", f"h_BDDKp_varied_{seed}", nbins, 4000, 8000)
    h_BuDDK0_varied = ROOT.TH1D(f"h_BuDDK0_varied_{seed}", f"h_BuDDK0_varied_{seed}", nbins, 4000, 8000)
    h_BuDD_varied = ROOT.TH1D(f"h_BuDD_varied_{seed}", f"h_BuDD_varied_{seed}", nbins, 4000, 8000)


    np.random.seed(seed)
    counts = np.array([np.random.poisson( h_sig.GetBinContent(i+1)*h_sig.GetEntries() ) for i in range(nbins)])
    for i, count in enumerate(counts):
        h_sig_varied.SetBinContent(i+1, count)

    counts_1 = np.array([np.random.poisson( h_comb.GetBinContent(i+1)*h_comb.GetEntries() ) for i in range(nbins)])
    for i, count in enumerate(counts_1):
        h_comb_varied.SetBinContent(i+1, count)
    
    counts_2 = np.array([np.random.poisson( h_BDDKp.GetBinContent(i+1)*h_BDDKp.GetEntries() ) for i in range(nbins)])
    for i, count in enumerate(counts_2):
        h_BDDKp_varied.SetBinContent(i+1, count)

    counts_3 = np.array([np.random.poisson( h_BuDDK0.GetBinContent(i+1)*h_BuDDK0.GetEntries() ) for i in range(nbins)])
    for i, count in enumerate(counts_3):
        h_BuDDK0_varied.SetBinContent(i+1, count)

    counts_4 = np.array([np.random.poisson( h_BuDD.GetBinContent(i+1)*h_BuDD.GetEntries() ) for i in range(nbins)])
    for i, count in enumerate(counts_4):
        h_BuDD_varied.SetBinContent(i+1, count)

    h_sig_varied.Sumw2()
    h_comb_varied.Sumw2()
    h_BDDKp_varied.Sumw2()
    h_BuDDK0_varied.Sumw2()
    h_BuDD_varied.Sumw2()
    
    h_sig_varied.Scale(1/h_sig_varied.Integral())
    h_comb_varied.Scale(1/h_comb_varied.Integral())
    h_BDDKp_varied.Scale(1/h_BDDKp_varied.Integral())
    h_BuDDK0_varied.Scale(1/h_BuDDK0_varied.Integral())
    h_BuDD_varied.Scale(1/h_BuDD_varied.Integral())

    data_sig = [h_sig_varied.GetBinContent(i+1) for i in range(nbins)]
    data_comb = [h_comb_varied.GetBinContent(i+1) for i in range(nbins)]
    data_BDDKp = [h_BDDKp_varied.GetBinContent(i+1) for i in range(nbins)]
    data_BuDDK0 = [h_BuDDK0_varied.GetBinContent(i+1) for i in range(nbins)]
    data_BuDD = [h_BuDD_varied.GetBinContent(i+1) for i in range(nbins)]

    data_sig_err = [h_sig_varied.GetBinError(i+1) for i in range(nbins)]
    data_comb_err = [h_comb_varied.GetBinError(i+1) for i in range(nbins)]
    data_BDDKp_err = [h_BDDKp_varied.GetBinError(i+1) for i in range(nbins)]
    data_BuDDK0_err = [h_BuDDK0_varied.GetBinError(i+1) for i in range(nbins)]
    data_BuDD_err = [h_BuDD_varied.GetBinError(i+1) for i in range(nbins)]

    varied_data = [data_sig, data_comb, data_BDDKp, data_BuDDK0, data_BuDD]
    varied_data_err = [data_sig_err, data_comb_err, data_BDDKp_err, data_BuDDK0_err, data_BuDD_err]

    return varied_data, varied_data_err

def convert_histogram_to_array(histograms, histogram_errors):
    h_sig = histograms[0]  
    h_comb = histograms[1]  
    h_BDDKp = histograms[2]
    h_BuDDK0 = histograms[3]
    h_BuDD = histograms[4]

    h_sig_err = histogram_errors[0]  
    h_comb_err = histogram_errors[1]  
    h_BDDKp_err = histogram_errors[2]
    h_BuDDK0_err = histogram_errors[3]
    h_BuDD_err = histogram_errors[4]

    nbins = h_sig.GetNbinsX()
    data_sig = [h_sig.GetBinContent(i+1) for i in range(nbins)]
    data_comb = [h_comb.GetBinContent(i+1) for i in range(nbins)]
    data_BDDKp = [h_BDDKp.GetBinContent(i+1) for i in range(nbins)]
    data_BuDDK0 = [h_BuDDK0.GetBinContent(i+1) for i in range(nbins)]
    data_BuDD = [h_BuDD.GetBinContent(i+1) for i in range(nbins)]

    data_sig_err = [h_sig_err.GetBinContent(i+1) for i in range(nbins)]
    data_comb_err = [h_comb_err.GetBinContent(i+1) for i in range(nbins)]
    data_BDDKp_err = [h_BDDKp_err.GetBinContent(i+1) for i in range(nbins)]
    data_BuDDK0_err = [h_BuDDK0_err.GetBinContent(i+1) for i in range(nbins)]
    data_BuDD_err = [h_BuDD_err.GetBinContent(i+1) for i in range(nbins)]

    data = [data_sig, data_comb, data_BDDKp, data_BuDDK0, data_BuDD]
    data_err = [data_sig_err, data_comb_err, data_BDDKp_err, data_BuDDK0_err, data_BuDD_err]

    return data, data_err


def convert_toy_data_histogram_to_array(histo):
    nbins = histo.GetNbinsX()
    return [histo.GetBinContent(i+1) for i in range(nbins)]


def write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yields, yields_err, c, c_err, add_physics_backgrounds):
    data_sig = histo_data[0]
    data_comb = histo_data[1]
    data_BDDKp = histo_data[2]
    data_BuDDK0 = histo_data[3]
    data_BuDD = histo_data[4]

    data_sig_err = histo_data_err[0]
    data_comb_err = histo_data_err[1]
    data_BDDKp_err = histo_data_err[2]
    data_BuDDK0_err = histo_data[3]
    data_BuDD_err = histo_data[4]

    N_comb = yields[0]
    N_BDDKp = yields[1]
    N_BuDDK0 = yields[2]
    N_BuDD = yields[3]

    N_BDDKp_err = yields_err[0]
    N_BuDDK0_err = yields_err[1]
    N_BuDD_err = yields_err[2]

    # Compute values for Gaussian constraints
    # Signal
    if(c == 0):
        upward = 1.999
        downward = 0.001
    else:
        upward = 1+(c_err/c)
        downward = 1-(c_err/c)

    # B -> DD K+
    if(N_BDDKp == 0):
        upward_BDDKp = 1.999
        downward_BDDKp = 0.001
    else:
        upward_BDDKp = 1+(N_BDDKp_err/N_BDDKp)
        downward_BDDKp = 1-(N_BDDKp_err/N_BDDKp)

    # B+ -> DD K0
    if(N_BuDDK0 == 0):
        upward_BuDDK0 = 1.999
        downward_BuDDK0 = 0.001
    else:
        upward_BuDDK0 = 1+(N_BuDDK0_err/N_BuDDK0)
        downward_BuDDK0 = 1-(N_BuDDK0_err/N_BuDDK0)

    # B+ -> DD
    if(N_BuDD == 0):
        upward_BuDD = 1.999
        downward_BuDD = 0.001
    else:
        upward_BuDD = 1+(N_BuDD_err/N_BuDD)
        downward_BuDD = 1-(N_BuDD_err/N_BuDD)

    # Modifiers
    sig_modifiers = [{"data": None, "name": "c_bdt1_bdt2", "type": "normfactor"}, {"data": None, "name": "BF_sig", "type": "normfactor"}, {"data": {"hi": upward, "lo": downward}, "name": "c_bdt1_bdt2_syst", "type": "normsys"}]
    comb_modifiers = [{"data": None, "name": "N_comb", "type": "normfactor"}]   
    BDDKp_modifiers = [{"data": None, "name": "N_BDDKp", "type": "normfactor"}, {"data": {"hi": upward_BDDKp, "lo": downward_BDDKp}, "name": "N_BDDKp_syst", "type": "normsys"}]                 
    BuDDK0_modifiers = [{"data": None, "name": "N_BuDDK0", "type": "normfactor"}, {"data": {"hi": upward_BuDDK0, "lo": downward_BuDDK0}, "name": "N_BuDDK0_syst", "type": "normsys"}]                 
    BuDD_modifiers = [{"data": None, "name": "N_BuDD", "type": "normfactor"}, {"data": {"hi": upward_BuDD, "lo": downward_BuDD}, "name": "N_BuDD_syst", "type": "normsys"}]                 

    if(add_statistical_error):
        # shapesys
        # sig_modifiers.append({"name": "shapesys_sig", "data": data_sig_err, "type": "shapesys"})
        # comb_modifiers.append({"name": "shapesys_comb", "data": data_comb_err, "type": "shapesys"})
        # BDDKp_modifiers.append({"name": "shapesys_BDDKp", "data": data_BDDKp_err, "type": "shapesys"})
        # BuDDK0_modifiers.append({"name": "shapesys_BuDDK0", "data": data_BuDDK0_err, "type": "shapesys"})
        # BuDD_modifiers.append({"name": "shapesys_BuDD", "data": data_BuDD_err, "type": "shapesys"})
        # staterror (light)
        # sig_modifiers.append({"name": "staterror", "data": data_sig_err, "type": "staterror"})
        # comb_modifiers.append({"name": "staterror", "data": data_comb_err, "type": "staterror"})
        # BDDKp_modifiers.append({"name": "staterror", "data": data_BDDKp_err, "type": "staterror"})
        # BuDDK0_modifiers.append({"name": "staterror", "data": data_BuDDK0_err, "type": "staterror"})
        # BuDD_modifiers.append({"name": "staterror", "data": data_BuDD_err, "type": "staterror"})
        # staterror (full)
        sig_modifiers.append({"name": "staterror_sig", "data": data_sig_err, "type": "staterror"})
        comb_modifiers.append({"name": "staterror_comb", "data": data_comb_err, "type": "staterror"})
        BDDKp_modifiers.append({"name": "staterror_BDDKp", "data": data_BDDKp_err, "type": "staterror"})
        BuDDK0_modifiers.append({"name": "staterror_BuDDK0", "data": data_BuDDK0_err, "type": "staterror"})
        BuDD_modifiers.append({"name": "staterror_BuDD", "data": data_BuDD_err, "type": "staterror"})

    # Samples
    samples = [{"name": "Signal", "data": data_sig, "modifiers": sig_modifiers}, {"name": "Combinatorial", "data": data_comb, "modifiers": comb_modifiers}]
    
    # Parameters
    parameters = [{"name": "lumi", "auxdata": [1.0], "bounds": [[0.5,1.5]], "fixed": True, "inits": [1.0], "sigmas": [0.1]}, {"name": "BF_sig", "bounds": [[-1.0,1.0]],  "inits": [nsig/c]}, {"name": "N_comb", "bounds": [[-10*(N_comb+1),10*(N_comb+1)]], "inits": [N_comb]}, {"name": "c_bdt1_bdt2", "bounds": [[0.0,1000000000.0]], "fixed": True, "inits": [c]}]

    if(add_physics_backgrounds):
        if((N_BDDKp != 0) and (upward_BDDKp < 2) and (downward_BDDKp > 0)):
            parameters.append({"bounds": [[-10*(N_BDDKp+1),10*(N_BDDKp+1)]], "fixed": True, "inits": [N_BDDKp], "name": "N_BDDKp"})
            samples.append({"name": "BDDKp", "data": data_BDDKp, "modifiers": BDDKp_modifiers})
        if((N_BuDDK0 != 0) and (upward_BuDDK0 < 2) and (downward_BuDDK0)):
            parameters.append({"bounds": [[-10*(N_BuDDK0+1),10*(N_BuDDK0+1)]], "fixed": True, "inits": [N_BuDDK0], "name": "N_BuDDK0"})
            samples.append({"name": "BuDDK0", "data": data_BuDDK0, "modifiers": BuDDK0_modifiers})
        if((N_BuDD != 0) and (upward_BuDD < 2) and (downward_BuDD > 0)):
            parameters.append({"bounds": [[-10*(N_BuDD+1),10*(N_BuDD+1)]], "fixed": True, "inits": [N_BuDD], "name": "N_BuDD"})
            samples.append({"name": "BuDD", "data": data_BuDD, "modifiers": BuDD_modifiers})

    spec = {
        "channels": [
            {
                "name": "AllEvents",
                "samples": samples
            }
        ],
        "measurements": [
            {
                "config": {
                    "parameters": parameters,
                    "poi": "BF_sig"
                },
                "name": "MinimalConfiguration"
            }
        ],
        "observations": [
            {
                "name": "AllEvents",
                "data": histo_toy_data,
            }
        ],
        "version": "1.0.0"}

    return spec


def plot(fit_name, h_data, workspace, fit_pars, fit_errors, nsig, bdt1, bdt2, seed):
    model = workspace.model()

    nbins = sum(model.config.channel_nbins.values())
    mass = [h_data.GetBinCenter(i+1) for i in range(nbins)] 
    edges = [h_data.GetBinLowEdge(i+1) for i in range(nbins+1)] 

    if(fit_name == "all_events"):
        sig_index = model.config.par_names.index("BF_sig")
        comb_index = model.config.par_names.index("N_comb")
        c_index = model.config.par_names.index("c_bdt1_bdt2")
        try:
            BDDKp_index = model.config.par_names.index("N_BDDKp")
        except:
            BDDKp_index = -1
        try:
            BuDDK0_index = model.config.par_names.index("N_BuDDK0")
        except:
            BuDDK0_index = -1
        try:
            BuDD_index = model.config.par_names.index("N_BuDD")
        except:
            BuDD_index = -1

        fit_poi = fit_pars[sig_index] 
        fit_poi_error = fit_errors[sig_index] 
        print("BF_sig = ", fit_poi, " +/- ", fit_poi_error)

        ncomb = fit_pars[comb_index]
        ncomb_err = fit_errors[comb_index]
        print("N_comb = ", ncomb, " +/- ", ncomb_err)

        if(BDDKp_index >= 0):
            n_BDDKp = fit_pars[BDDKp_index]
            n_BDDKp_err = fit_errors[BDDKp_index]
            print("N_BDDKp = ", n_BDDKp, " +/- ", n_BDDKp_err)

        if(BuDDK0_index >= 0):
            n_BuDDK0 = fit_pars[BuDDK0_index]
            n_BuDDK0_err = fit_errors[BuDDK0_index]
            print("N_BuDDK0 = ", n_BuDDK0, " +/- ", n_BuDDK0_err)

        if(BuDD_index >= 0):
            n_BuDD = fit_pars[BuDD_index]
            n_BuDD_err = fit_errors[BuDD_index]
            print("N_BuDD = ", n_BuDD, " +/- ", n_BuDD_err)

    else:
        sig_index = model.config.par_names.index("BF_sig")
        comb_index_1 = model.config.par_names.index("N_comb_1")
        comb_index_2 = model.config.par_names.index("N_comb_2")
        comb_index_3 = model.config.par_names.index("N_comb_3")
        c_index_1 = model.config.par_names.index("c_bdt1_bdt2_1")
        c_index_2 = model.config.par_names.index("c_bdt1_bdt2_2")
        c_index_3 = model.config.par_names.index("c_bdt1_bdt2_3")
        try:
            BDDKp_index_1 = model.config.par_names.index("N_BDDKp_1")
            BDDKp_index_2 = model.config.par_names.index("N_BDDKp_2")
            BDDKp_index_3 = model.config.par_names.index("N_BDDKp_3")
        except:
            BDDKp_index_1 = -1
            BDDKp_index_2 = -1
            BDDKp_index_3 = -1
        try:
            BuDDK0_index_1 = model.config.par_names.index("N_BuDDK0_1")
            BuDDK0_index_2 = model.config.par_names.index("N_BuDDK0_2")
            BuDDK0_index_3 = model.config.par_names.index("N_BuDDK0_3")
        except:
            BuDDK0_index_1 = -1
            BuDDK0_index_2 = -1
            BuDDK0_index_3 = -1
        try:
            BuDD_index_1 = model.config.par_names.index("N_BuDD_1")
            BuDD_index_2 = model.config.par_names.index("N_BuDD_2")
            BuDD_index_3 = model.config.par_names.index("N_BuDD_3")
        except:
            BuDD_index_1 = -1
            BuDD_index_2 = -1
            BuDD_index_3 = -1

        fit_poi = fit_pars[sig_index]
        nbkg_1 = fit_pars[comb_index_1]
        nbkg_2 = fit_pars[comb_index_2]
        nbkg_3 = fit_pars[comb_index_3]

        fit_poi_error = fit_errors[sig_index]
        nbkg_1_error = fit_errors[comb_index_1]
        nbkg_2_error = fit_errors[comb_index_2]
        nbkg_3_error = fit_errors[comb_index_3]

        print("BF_sig = ", fit_poi, " +/- ", fit_poi_error)
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
    
    if(fit_name == "all_events"):

        data_sig = mc_counts[sig_index][0]
        data_comb = mc_counts[comb_index][0]
        if(BDDKp_index >= 0):
            data_BDDKp = mc_counts[BDDKp_index][0]
        if(BuDDK0_index >= 0):
            data_BuDDK0 = mc_counts[BuDDK0_index][0]
        if(BuDD_index >= 0):
            data_BuDD = mc_counts[BuDD_index][0]

        all_data = data_sig+data_comb
        if(BDDKp_index >= 0):
            all_data += data_BDDKp
        if(BuDDK0_index >= 0):
            all_data += data_BuDDK0
        if(BuDD_index >= 0):
            all_data += data_BuDD

        f, ax = plt.subplots()

        plt.step(edges, np.append(all_data, all_data[-1]), where='post', color='black', label="Total fit")
        plt.step(edges, np.append(data_sig, data_sig[-1]), where='post', color='blue', label='Signal')
        plt.step(edges, np.append(data_comb, data_comb[-1]), where='post', color='red', label='Combinatorial')
        if(BDDKp_index >= 0):
            plt.step(edges, np.append(data_BDDKp, data_BDDKp[-1]), where='post', color='green', label='$B \\to DD K^+$')
        if(BuDDK0_index >= 0):
            plt.step(edges, np.append(data_BuDDK0, data_BuDDK0[-1]), where='post', color='darkviolet', label='$B^+ \\to DD K^0$')
        if(BuDD_index >= 0):
            plt.step(edges, np.append(data_BuDD, data_BuDD[-1]), where='post', color='orange', label='$B^+ \\to DD$')

        the_data = workspace.data(model, include_auxdata=False)
        the_data_err = np.sqrt(the_data)

        plt.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
        plt.title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed} \n All events')
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
        data_sig = mc_counts[sig_index][0]
        data_sig_1 = data_sig[0:nbins_1]
        data_sig_2 = data_sig[nbins_1:nbins_1+nbins_2]
        data_sig_3 = data_sig[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

        data_comb = mc_counts[comb_index_1][0]
        data_bkg_1 = data_comb[0:nbins_1]
        data_bkg_2 = data_comb[nbins_1:nbins_1+nbins_2]
        data_bkg_3 = data_comb[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

        all_data_1 = data_sig_1+data_bkg_1
        all_data_2 = data_sig_2+data_bkg_2
        all_data_3 = data_sig_3+data_bkg_3

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
        ax1.errorbar(mass_1, the_data_1, the_data_err_1, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax2.step(edges_2, np.append(all_data_2, all_data_2[-1]), where='post', color='black', label="Total fit")
        ax2.step(edges_2, np.append(data_sig_2, data_sig_2[-1]), where='post', color='blue', label='Signal')
        ax2.step(edges_2, np.append(data_bkg_2, data_bkg_2[-1]), where='post', color='red', label='Combinatorial')
        ax2.errorbar(mass_2, the_data_2, the_data_err_2, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax3.step(edges_3, np.append(all_data_3, all_data_3[-1]), where='post', color='black', label="Total fit")
        ax3.step(edges_3, np.append(data_sig_3, data_sig_3[-1]), where='post', color='blue', label='Signal')
        ax3.step(edges_3, np.append(data_bkg_3, data_bkg_3[-1]), where='post', color='red', label='Combinatorial')
        ax3.errorbar(mass_3, the_data_3, the_data_err_3, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        ax1.set_title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed} \n Channel 1')
        ax2.set_title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed} \n Channel 2')
        ax3.set_title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {seed} \n Channel 3')

        ax1.set_xlabel('m_B (MeV)')
        ax1.set_ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins_1, 1) ))

        ax2.set_xlabel('m_B (MeV)')
        ax2.set_ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins_2, 1) ))

        ax3.set_xlabel('m_B (MeV)')
        ax3.set_ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins_3, 1) ))

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
    fit_parnames = model.config.par_names

    if(save_plot):
        print(res_obj)
        plot(fit_name, h_data, workspace, fit_pars, fit_errors, nsig, bdt1, bdt2, seed)
    # if(run_cls):

    return fit_pars, fit_errors, fit_parnames


def toy_studies(fit_name, nsig, bdt1, bdt2, BF_toys, BF_toys_err, Nbkg_toys, Nbkg_toys_err, Nbkg_1_toys, Nbkg_1_toys_err, Nbkg_2_toys, Nbkg_2_toys_err, Nbkg_3_toys, Nbkg_3_toys_err, c_values, yield_values, nbins=30):
    N = len(BF_toys)

    if(fit_name == "all_events"):
        fig, ax = plt.subplots(3, 2, figsize=(10, 15))
    else:
        fig, ax = plt.subplots(3, 4, figsize=(16, 12))

    c = c_values[0]

    # BF pull
    pull = (BF_toys - nsig)/BF_toys_err
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
    ax[1,0].axvline(nsig/c, color='black', linestyle='--', linewidth=2)
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

    # N_comb
    N_comb = yield_values[0][0]
    N_comb_1 = yield_values[1][0]
    N_comb_2 = yield_values[2][0]
    N_comb_3 = yield_values[3][0]

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
        pull_2 = (Nbkg_1_toys - N_comb_1)/Nbkg_1_toys_err
        (mu_pull_2, sigma_pull_2) = norm.fit(pull_2)
        n_2, bins_2, patches_2 = ax[0,1].hist( pull_2, bins=nbins)
        xcenters_2 = (bins_2[:-1] + bins_2[1:]) / 2
        ax[0,1].errorbar(xcenters_2, n_2, yerr=np.sqrt(n_2), ecolor='black', fmt='k.')
        y_pull_2 = sum(n_2)*(bins_2[1] - bins_2[0])*norm.pdf( xcenters_2, mu_pull_2, sigma_pull_2)
        ax[0,1].plot(xcenters_2, y_pull_2, 'r-', linewidth=2)
        ax[0,1].set_xlabel(f"Pull: ($N_{{bkg, 1}} - {N_comb_1})/error")
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

        ax[1,1].hist(Nbkg_1_toys, bins=nbins)
        ax[1,1].set_xlabel("$N_{bkg, 1}$")
        ax[1,1].set_ylabel(f"Entries / {nbins} bins")
        ax[1,1].set_title("Background yield 1")
        ax[1,1].axvline(N_comb_1, color='black', linestyle='--', linewidth=2)
        ax[1,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_1_toys), np.std(Nbkg_1_toys) ),
                transform=ax[1,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,1].hist(Nbkg_1_toys_err, bins=nbins)
        ax[2,1].set_xlabel("$N_{bkg, 1}$ error")
        ax[2,1].set_ylabel(f"Entries / {nbins} bins")
        ax[2,1].set_title("Background yield 1 error")
        ax[2,1].hist(Nbkg_1_toys_err, bins=nbins)
        ax[2,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_1_toys_err), np.std(Nbkg_1_toys_err) ),
                transform=ax[2,1].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        pull_3 = (Nbkg_2_toys - N_comb_2)/Nbkg_2_toys_err
        (mu_pull_3, sigma_pull_3) = norm.fit(pull_3)
        n_3, bins_3, patches_3 = ax[0,2].hist( pull_3, bins=nbins)
        xcenters_3 = (bins_3[:-1] + bins_3[1:]) / 2
        ax[0,2].errorbar(xcenters_3, n_3, yerr=np.sqrt(n_3), ecolor='black', fmt='k.')
        y_pull_3 = sum(n_3)*(bins_3[1] - bins_3[0])*norm.pdf( xcenters_3, mu_pull_3, sigma_pull_3)
        ax[0,2].plot(xcenters_3, y_pull_3, 'r-', linewidth=2)
        ax[0,2].set_xlabel(f"Pull: ($N_{{bkg, 2}} - {N_comb_2})/error")
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

        ax[1,2].hist(Nbkg_2_toys, bins=nbins)
        ax[1,2].set_xlabel("$N_{bkg, 2}$")
        ax[1,2].set_ylabel(f"Entries / {nbins} bins")
        ax[1,2].set_title("Background yield 2")
        ax[1,2].axvline(N_comb_2, color='black', linestyle='--', linewidth=2)
        ax[1,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_2_toys), np.std(Nbkg_2_toys) ),
                transform=ax[1,2].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,2].hist(Nbkg_2_toys_err, bins=nbins)
        ax[2,2].set_xlabel("$N_{bkg, 2}$ error")
        ax[2,2].set_ylabel(f"Entries / {nbins} bins")
        ax[2,2].set_title("Background yield 2 error")
        ax[2,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_2_toys_err), np.std(Nbkg_2_toys_err) ),
                transform=ax[2,2].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        pull_4 = (Nbkg_3_toys - N_comb_3)/Nbkg_3_toys_err
        (mu_pull_4, sigma_pull_4) = norm.fit(pull_4)
        n_4, bins_4, patches_4 = ax[0,3].hist( pull_4, bins=nbins)
        xcenters_4 = (bins_4[:-1] + bins_4[1:]) / 2
        ax[0,3].errorbar(xcenters_4, n_4, yerr=np.sqrt(n_4), ecolor='black', fmt='k.')
        y_pull_4 = sum(n_4)*(bins_4[1] - bins_4[0])*norm.pdf( xcenters_4, mu_pull_4, sigma_pull_4)
        ax[0,3].plot(xcenters_4, y_pull_4, 'r-', linewidth=2)
        ax[0,3].set_xlabel(f"Pull: ($N_{{bkg, 3}} - {Nbkg_3_toys})/error")
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

        ax[1,3].hist(Nbkg_3_toys, bins=nbins)
        ax[1,3].set_xlabel("$N_{bkg, 3}$")
        ax[1,3].set_ylabel("Entries / {0} bins".format(nbins))
        ax[1,3].set_title("Background yield 3")
        ax[1,3].axvline(N_comb_3, color='black', linestyle='--', linewidth=2)
        ax[1,3].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_3_toys), np.std(Nbkg_3_toys) ),
                transform=ax[1,3].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        ax[2,3].hist(Nbkg_3_toys_err, bins=nbins)
        ax[2,3].set_xlabel("$N_{bkg, 3}$ error")
        ax[2,3].set_ylabel("Entries / {0} bins".format(nbins))
        ax[2,3].set_title("Background yield 3 error")
        ax[2,3].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(Nbkg_3_toys_err), np.std(Nbkg_3_toys_err) ),
                transform=ax[2,3].transAxes,
                fontsize=12,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    fig.suptitle(f"Fit to all events \n BDT1 = {bdt1} | BDT2 = {bdt2}", fontsize=24)

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
    yield_values = np.load(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/yield_values_bdt1_{bdt1}_bdt2_{bdt2}.npy", allow_pickle=True)
    yield_errors = np.load(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/yield_errors_bdt1_{bdt1}_bdt2_{bdt2}.npy", allow_pickle=True)

    c_value = np.load(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/c_value_bdt1_{bdt1}_bdt2_{bdt2}.npy")
    c_error = np.load(f"/panfs/felician/B2Ktautau/workflow/generate_histograms/c_error_bdt1_{bdt1}_bdt2_{bdt2}.npy")
    
    if(fit_name == "all_events"):
        histograms, histogram_errors = retrieve_histograms(f,f1,0)
        histo_data, histo_data_err = convert_histogram_to_array(histograms, histogram_errors)

        yields = yield_values[0]
        yields_err = yield_errors[0]

        c = c_value[0]
        c_err = c_error[0]
    else:
        histograms_1, histogram_errors_1 = retrieve_histograms(f,f1,1)
        histograms_2, histogram_errors_2 = retrieve_histograms(f,f2,2)
        histograms_3, histogram_errors_3 = retrieve_histograms(f,f2,3)
        histo_data_1, histo_data_err_1 = convert_histogram_to_array(histograms_1, histogram_errors_1)
        histo_data_2, histo_data_err_2 = convert_histogram_to_array(histograms_2, histogram_errors_2)
        histo_data_3, histo_data_err_3 = convert_histogram_to_array(histograms_3, histogram_errors_3)

        yields_1 = yield_values[1]
        yields_err_1 = yield_errors[1]
    
        yields_2 = yield_values[2]
        yields_err_2 = yield_errors[2]

        yields_3 = yield_values[3]
        yields_err_3 = yield_errors[3]

        c_1 = c_value[1]
        c_err_1 = c_error[1]

        c_2 = c_value[2]
        c_err_2 = c_error[2]

        c_3 = c_value[3]
        c_err_3 = c_error[3]

    end = time.time()
    print(f"Elapsed time (retrieving files): {end - start:.2f} seconds")

    BF_toys = np.zeros(n_toys)
    Nbkg_toys = np.zeros(n_toys)
    Nbkg_1_toys = np.zeros(n_toys)
    Nbkg_2_toys = np.zeros(n_toys)
    Nbkg_3_toys = np.zeros(n_toys)

    BF_toys_err = np.zeros(n_toys)
    Nbkg_toys_err = np.zeros(n_toys)
    Nbkg_1_toys_err = np.zeros(n_toys)
    Nbkg_2_toys_err = np.zeros(n_toys)
    Nbkg_3_toys_err = np.zeros(n_toys)

    N_fail = 0

    if(validate_fit):
        for seed in range(n_toys):
            if seed % (n_toys // 10) == 0:
                print(f"Progress: {100 * seed // n_toys}% (seed={seed})")

            if(seed == 0):
                save_plot = True
            else:
                save_plot = False
            
            if(fit_name == "all_events"):
                h_data = generate_toy_data(histograms, yields, seed, nsig, add_physics_backgrounds)
                histo_toy_data = convert_toy_data_histogram_to_array(h_data)

                if(vary_fit_templates):
                    histo_data, histo_data_err = vary_templates(histograms, seed)

                spec = write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yields, yields_err, c, c_err, add_physics_backgrounds)

                workspace = pyhf.Workspace(spec)
                model = workspace.model()
                sig_index = model.config.par_names.index("BF_sig")
                comb_index = model.config.par_names.index("N_comb")

                try:
                    fit_pars, fit_errors, fit_parnames = fit(fit_name, nsig, bdt1, bdt2, h_data, spec, seed, save_plot, False)

                    BF_toys[seed] = fit_pars[sig_index]
                    BF_toys_err[seed] = fit_errors[sig_index]
        
                    Nbkg_toys[seed] = fit_pars[comb_index]
                    Nbkg_toys_err[seed] = fit_errors[comb_index]
                except:
                    BF_toys[seed] = -999
                    BF_toys_err[seed] = 999
        
                    Nbkg_toys[seed] = -999
                    Nbkg_toys_err[seed] = 999

                    N_fail += 1

            else:
                h_data_1 = generate_toy_data(histograms_1, yields_1, seed, nsig, add_physics_backgrounds)
                h_data_2 = generate_toy_data(histograms_2, yields_2, seed, nsig, add_physics_backgrounds)
                h_data_3 = generate_toy_data(histograms_3, yields_3, seed, nsig, add_physics_backgrounds)

                histo_toy_data_1 = convert_toy_data_histogram_to_array(h_data_1)
                histo_toy_data_2 = convert_toy_data_histogram_to_array(h_data_2)
                histo_toy_data_3 = convert_toy_data_histogram_to_array(h_data_3)

                if((h_data_1.Integral() < 10) or (h_data_2.Integral() < 10) or (h_data_3.Integral() < 10)):
                    continue

        toy_studies(fit_name, nsig, bdt1, bdt2, BF_toys, BF_toys_err, Nbkg_toys, Nbkg_toys_err, Nbkg_1_toys, Nbkg_1_toys_err, Nbkg_2_toys, Nbkg_2_toys_err, Nbkg_3_toys, Nbkg_3_toys_err, c_value, yield_values)


    else:
        h_data = generate_toy_data(histograms, yields, 1000, nsig, add_physics_backgrounds)
        histo_toy_data = convert_toy_data_histogram_to_array(h_data)

        spec = write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yields, yields_err, c, c_err, add_physics_backgrounds)

        fit(fit_name, nsig, spec, True, True)


    end1 = time.time()
    print(f"{N_fail}/{n_toys} toys fail")
    print(f"Elapsed time (toy data generation): {end1 - end:.2f} seconds")


if __name__ == "__main__":
    main(sys.argv)