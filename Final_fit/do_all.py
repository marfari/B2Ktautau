import sys
import ROOT
import numpy as np
import time
import pyhf
import matplotlib.pyplot as plt
from scipy.stats import norm

add_physics_backgrounds = False
add_statistical_error = True
validate_fit = True

n_toys = 1000

def retrieve_histograms(f, ch):
    h_sig = f.Get(f"Channel_{ch}/h_sig_{ch}")
    h_comb = f.Get(f"Channel_{ch}/h_comb_{ch}")
    h_BDDKp = f.Get(f"Channel_{ch}/h_BDDKp_{ch}")
    h_BuDDK0 = f.Get(f"Channel_{ch}/h_BuDDK0_{ch}")
    h_BuDD = f.Get(f"Channel_{ch}/h_BuDD_{ch}")

    histograms = [h_sig, h_comb, h_BDDKp, h_BuDDK0, h_BuDD]

    return histograms


def retrieve_histogram_errors(f, ch):
    h_sig_err = f.Get(f"Channel_{ch}/h_sig_err_{ch}")
    h_comb_err = f.Get(f"Channel_{ch}/h_comb_err_{ch}")
    h_BDDKp_err = f.Get(f"Channel_{ch}/h_BDDKp_err_{ch}")
    h_BuDDK0_err = f.Get(f"Channel_{ch}/h_BuDDK0_err_{ch}")
    h_BuDD_err = f.Get(f"Channel_{ch}/h_BuDD_err_{ch}")

    histogram_errors = [h_sig_err, h_comb_err, h_BDDKp_err, h_BuDDK0_err, h_BuDD_err]

    return histogram_errors


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

    h_data.Sumw2()
    h_data_comb.Sumw2()
    if(nsig != 0):
        h_data_sig.Sumw2()
    if(add_physics_backgrounds):
        h_data_BDDKp.Sumw2()
        h_data_BuDDK0.Sumw2()
        h_data_BuDD.Sumw2()

    np.random.seed(seed+1)
    if((seed < 1000) and (seed >= 0)):
        N_comb = np.random.poisson(N_comb)
        if(add_physics_backgrounds):
            N_BDDKp = np.random.poisson(N_BDDKp)
            N_BuDDK0 = np.random.poisson(N_BuDDK0)
            N_BuDD = np.random.poisson(N_BuDD)


    if((seed < 1000) and (seed >= 0)):
        # toy_counts = np.array([np.random.poisson( h_comb.GetBinContent(i+1)*N_comb ) for i in range(nbins)])
        # for i, count in enumerate(toy_counts):
        #     h_data_comb.SetBinContent(i+1, count)

        # if(nsig != 0):
        #     toy_counts_sig = np.array([np.random.poisson( h_sig.GetBinContent(i+1)*nsig ) for i in range(nbins)])
        #     for i, count in enumerate(toy_counts_sig):
        #         h_data_sig.SetBinContent(i+1, count)

        # if(add_physics_backgrounds):
        #     toy_counts_BDDKp = np.array([np.random.poisson( h_BDDKp.GetBinContent(i+1)*N_BDDKp ) for i in range(nbins)])
        #     for i, count in enumerate(toy_counts_BDDKp):
        #         h_data_BDDKp.SetBinContent(i+1, count)

        #     toy_counts_BuDDK0 = np.array([np.random.poisson( h_BuDDK0.GetBinContent(i+1)*N_BuDDK0 ) for i in range(nbins)])
        #     for i, count in enumerate(toy_counts_BuDDK0):
        #         h_data_BuDDK0.SetBinContent(i+1, count)

        #     toy_counts_BuDD = np.array([np.random.poisson( h_BuDD.GetBinContent(i+1)*N_BuDD ) for i in range(nbins)])
        #     for i, count in enumerate(toy_counts_BuDD):
        #         h_data_BuDD.SetBinContent(i+1, count)

        # ROOT.gRandom.SetSeed(0)
        rnd = ROOT.TRandom()
        rnd.SetSeed(seed+1)

        h_data_comb.FillRandom(h_comb, int(N_comb), rnd)
        if(add_physics_backgrounds):
            h_data_BDDKp.FillRandom(h_BDDKp, int(N_BDDKp), rnd)
            h_data_BuDDK0.FillRandom(h_BuDDK0, int(N_BuDDK0), rnd)
            h_data_BuDD.FillRandom(h_BuDD, int(N_BuDD), rnd)

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


def convert_histogram_to_array(histograms):
    h_sig = histograms[0]  
    h_comb = histograms[1]  
    h_BDDKp = histograms[2]
    h_BuDDK0 = histograms[3]
    h_BuDD = histograms[4]

    nbins = h_sig.GetNbinsX()
    data_sig = [h_sig.GetBinContent(i+1) for i in range(nbins)]
    data_comb = [h_comb.GetBinContent(i+1) for i in range(nbins)]
    data_BDDKp = [h_BDDKp.GetBinContent(i+1) for i in range(nbins)]
    data_BuDDK0 = [h_BuDDK0.GetBinContent(i+1) for i in range(nbins)]
    data_BuDD = [h_BuDD.GetBinContent(i+1) for i in range(nbins)]

    return [data_sig, data_comb, data_BDDKp, data_BuDDK0, data_BuDD]


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

        if(upward > 2):
            upward = 1.999
        if(downward < 0):
            downward = 0.001

    # B -> DD K+
    if(N_BDDKp == 0):
        upward_BDDKp = 1.999
        downward_BDDKp = 0.001
    else:
        upward_BDDKp = 1+(N_BDDKp_err/N_BDDKp)
        downward_BDDKp = 1-(N_BDDKp_err/N_BDDKp)
        if(upward_BDDKp > 2):
            upward_BDDKp = 1.999
        if(downward_BDDKp < 0):
            downward_BDDKp = 0.001

    # B+ -> DD K0
    if(N_BuDDK0 == 0):
        upward_BuDDK0 = 1.999
        downward_BuDDK0 = 0.001
    else:
        upward_BuDDK0 = 1+(N_BuDDK0_err/N_BuDDK0)
        downward_BuDDK0 = 1-(N_BuDDK0_err/N_BuDDK0)
        if(upward_BuDDK0 > 2):
            upward_BuDDK0 = 1.999
        if(downward_BuDDK0 < 0):
            downward_BuDDK0 = 0.001

    # B+ -> DD
    if(N_BuDD_err == 0):
        upward_BuDD = 1.999
        downward_BuDD = 0.001
    else:
        upward_BuDD = 1+(N_BuDD_err/N_BuDD)
        downward_BuDD = 1-(N_BuDD_err/N_BuDD)
        if(upward_BuDD > 2):
            upward_BuDD = 1.999
        if(downward_BuDD < 0):
            downward_BuDD = 0.001

    # Modifiers
    sig_modifiers = [{"data": None, "name": "c_bdt1_bdt2", "type": "normfactor"}, {"data": None, "name": "BF_sig", "type": "normfactor"}, {"data": {"hi": upward, "lo": downward}, "name": "c_bdt1_bdt2_syst", "type": "normsys"}]
    comb_modifiers = [{"data": None, "name": "N_comb", "type": "normfactor"}]   
    BDDKp_modifiers = [{"data": None, "name": "N_BDDKp", "type": "normfactor"}, {"data": {"hi": upward_BDDKp, "lo": downward_BDDKp}, "name": "N_BDDKp_syst", "type": "normsys"}]                 
    BuDDK0_modifiers = [{"data": None, "name": "N_BuDDK0", "type": "normfactor"}, {"data": {"hi": upward_BuDDK0, "lo": downward_BuDDK0}, "name": "N_BuDDK0_syst", "type": "normsys"}]                 
    BuDD_modifiers = [{"data": None, "name": "N_BuDD", "type": "normfactor"}, {"data": {"hi": upward_BuDD, "lo": downward_BuDD}, "name": "N_BuDD_syst", "type": "normsys"}]                 

    if(add_statistical_error):
        # shapesys
        sig_modifiers.append({"name": "shapesys_sig", "data": data_sig_err, "type": "shapesys"})
        comb_modifiers.append({"name": "shapesys_comb", "data": data_comb_err, "type": "shapesys"})
        BDDKp_modifiers.append({"name": "shapesys_BDDKp", "data": data_BDDKp_err, "type": "shapesys"})
        BuDDK0_modifiers.append({"name": "shapesys_BuDDK0", "data": data_BuDDK0_err, "type": "shapesys"})
        BuDD_modifiers.append({"name": "shapesys_BuDD", "data": data_BuDD_err, "type": "shapesys"})
        # staterror
        # sig_modifiers.append({"name": "staterror_sig", "data": data_sig_err, "type": "staterror"})
        # comb_modifiers.append({"name": "staterror_comb", "data": data_comb_err, "type": "staterror"})
        # BDDKp_modifiers.append({"name": "staterror_BDDKp", "data": data_BDDKp_err, "type": "staterror"})
        # BuDDK0_modifiers.append({"name": "staterror_BuDDK0", "data": data_BuDDK0_err, "type": "staterror"})
        # BuDD_modifiers.append({"name": "staterror_BuDD", "data": data_BuDD_err, "type": "staterror"})

    # Samples
    samples = [{"name": "Signal", "data": data_sig, "modifiers": sig_modifiers}, {"name": "Combinatorial", "data": data_comb, "modifiers": comb_modifiers}]
    
    # Parameters
    parameters = [{"name": "lumi", "auxdata": [1.0], "bounds": [[0.5,1.5]], "fixed": True, "inits": [1.0], "sigmas": [0.1]}, {"name": "BF_sig", "bounds": [[-1.0,1.0]],  "inits": [nsig/c]}, {"name": "N_comb", "bounds": [[-10*(N_comb+1),10*(N_comb+1)]], "inits": [N_comb]}, {"name": "c_bdt1_bdt2", "bounds": [[0.0,1000000000.0]], "fixed": True, "inits": [c]}]


    if(add_physics_backgrounds):
        if(N_BDDKp != 0):
            parameters.append({"bounds": [[-10*(N_BDDKp+1),10*(N_BDDKp+1)]], "fixed": true, "inits": [N_BDDKp], "name": "N_BDDKp"})
            samples.append({"name": "BDDKp", "data": data_BDDKp, "modifiers": BDDKp_modifiers})
        if(N_BuDDK0 != 0):
            parameters.append({"bounds": [[-10*(N_BuDDK0+1),10*(N_BuDDK0+1)]], "fixed": true, "inits": [N_BuDDK0], "name": "N_BuDDK0"})
            samples.append({"name": "BuDDK0", "data": data_BuDDK0, "modifiers": BuDDK0_modifiers})
        if(N_BuDD != 0):
            parameters.append({"bounds": [[-10*(N_BuDD+1),10*(N_BuDD+1)]], "fixed": true, "inits": [N_BuDD], "name": "N_BuDD"})
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


def fit(fit_name, nsig, spec, save_plot, run_cls):

    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    data = workspace.data(model)

    run_fit = True    

    nbins = sum(model.config.channel_nbins.values())

    step_size = (8000-4000)/nbins
    mass = np.arange(4000,8000+step_size,step_size)
    observed_data = data[:nbins]

    # Fitting
    optimizer = pyhf.optimize.minuit_optimizer(verbose=0)
    pyhf.set_backend('numpy', optimizer) # minuit returns errors and correlation matrix; we ned the errors to make pulls

    fit_result, res_obj = pyhf.infer.mle.fit(data, model, return_uncertainties=True, return_result_obj=True)

    # print(res_obj)

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

    # if(save_plot):

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
        histograms = retrieve_histograms(f, 0)
        histo_data = convert_histogram_to_array(histograms)

        histogram_errors = retrieve_histogram_errors(f1, 0)
        histo_data_err = convert_histogram_to_array(histogram_errors)

        yields = yield_values[0]
        yields_err = yield_errors[0]

        c = c_value[0]
        c_err = c_error[0]
    else:
        histograms_1 = retrieve_histograms(f, 1)
        histograms_2 = retrieve_histograms(f, 2)
        histograms_3 = retrieve_histograms(f, 3)
        histo_data_1 = convert_histogram_to_array(histograms_1)
        histo_data_2 = convert_histogram_to_array(histograms_2)
        histo_data_3 = convert_histogram_to_array(histograms_3)

        histogram_errors_1 = retrieve_histogram_errors(f1, 1)
        histogram_errors_2 = retrieve_histogram_errors(f1, 2)
        histogram_errors_3 = retrieve_histogram_errors(f1, 3)
        histo_data_err_1 = convert_histogram_to_array(histogram_errors_1)
        histo_data_err_2 = convert_histogram_to_array(histogram_errors_2)
        histo_data_err_3 = convert_histogram_to_array(histogram_errors_3)

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

            if(fit_name == "all_events"):
                h_data = generate_toy_data(histograms, yields, seed, nsig, add_physics_backgrounds)
                histo_toy_data = convert_toy_data_histogram_to_array(h_data)

                spec = write_json_file(fit_name, nsig, histo_toy_data, histo_data, histo_data_err, yields, yields_err, c, c_err, add_physics_backgrounds)

                workspace = pyhf.Workspace(spec)
                model = workspace.model()
                sig_index = model.config.par_names.index("BF_sig")
                comb_index = model.config.par_names.index("N_comb")

                try:
                    fit_pars, fit_errors, fit_parnames = fit(fit_name, nsig, spec, True, True)

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