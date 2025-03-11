import sys
import pyhf
import json
import matplotlib.pyplot as plt
from IPython.display import display, HTML
import numpy as np
import matplotlib.ticker as mticker

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    with open("/panfs/felician/B2Ktautau/workflow/pyhf_fit/workspace/workspace_bdt1_{0}_bdt2_{1}_seed_{2}.json".format( round(bdt1,1), round(bdt2,1), random_seed )) as serialized:
        spec = json.load(serialized)

    workspace = pyhf.Workspace(spec)

    # Get model 
    model = workspace.model()

    # Get observations
    data = workspace.data(model)

    # Fitting
    ## maximum likelihood fit (produces the maximum likelihood = denominator of the likelihood ratio)
    pyhf.set_backend('numpy', 'minuit') # minuit returns errors and correlation matrix; we ned the errors to make pulls
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.minuit_optimizer(verbose=2))
    
    fit_result, twice_nll = pyhf.infer.mle.fit(data, model, return_uncertainties=True, return_fitted_val=True)
    fit_pars, fit_errors = fit_result.T
    # print(fit_pars)

    ## fitted POI value
    fit_poi = fit_pars[model.config.poi_index] # I need to save this for each toy (and make a pull plot later)
    fit_poi_error = fit_errors[model.config.poi_index] # I need to save this for each toy (and make a pull plot later)
    # print(fit_poi, "  ", fit_poi_error)

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

    def plot(order=[1,0], **par_settings):
        pars = pyhf.tensorlib.astensor(model.config.suggested_init())
        for k, v in par_settings.items():
            pars[par_name_dict[k]] = v

        mc_counts = get_mc_counts(pars)
        bottom = None
        names = ['Signal', 'Background']
        for i, sample_index in enumerate(order):
            data = mc_counts[sample_index][0]
            x = np.arange(len(data))
            plt.bar(x, data, 1, bottom=bottom, alpha=1.0, label=names[i])
            bottom = data if i == 0 else bottom + data
        plt.scatter(x, workspace.data(model, include_auxdata=False), c="k", alpha=1.0, zorder=99, label="Toy data (WS)")

    par_name_dict = {k: v["slice"].start for k, v in model.config.par_map.items()}
    plt.title(f'Nsig = {fit_poi:.3g} $\pm$ {fit_poi_error:.3g} \n BDT1 = {bdt1:.1g} | BDT2 = {bdt2:.1g} | seed = {random_seed}')
    plot(**{k: fit_pars[v] for k, v in par_name_dict.items()})
    plt.xlabel('m_B (MeV)')
    plt.ylabel('Entries / (40 MeV)')
    plt.legend()
    plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/fit_plot_bdt1_{0}_bdt2_{1}_seed_{2}.pdf'.format( round(bdt1,1), round(bdt2,1), random_seed))

    # Save fit result
    np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/fit_result_bdt1_{0}_bdt2_{1}_seed_{2}.npy'.format( round(bdt1,1), round(bdt2,1), random_seed ), [fit_poi,fit_poi_error])

    # Make pull plot (for every fit parameter)
    pulls = pyhf.tensorlib.concatenate(
        [
            (fit_pars[model.config.par_slice(k)] - model.config.param_set(k).suggested_init)
            / model.config.param_set(k).width()
            for k in model.config.par_order
            if model.config.param_set(k).constrained
        ]
    )

    pullerr = pyhf.tensorlib.concatenate(
        [
            fit_errors[model.config.par_slice(k)] / model.config.param_set(k).width()
            for k in model.config.par_order
            if model.config.param_set(k).constrained
        ]
    )

    labels = np.arange(len(pulls))
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 5)

    # set up axes labeling, ranges, etc...
    ax.xaxis.set_major_locator(mticker.FixedLocator(np.arange(labels.size).tolist()))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_xlim(-0.5, len(pulls) - 0.5)
    ax.set_title("Pull Plot", fontsize=18)
    ax.set_ylabel(r"$(\theta - \hat{\theta})\,/ \Delta \theta$", fontsize=18)

    # draw the +/- 2.0 horizontal lines
    ax.hlines([-2, 2], -0.5, len(pulls) - 0.5, colors="black", linestyles="dotted")
    # draw the +/- 1.0 horizontal lines
    ax.hlines([-1, 1], -0.5, len(pulls) - 0.5, colors="black", linestyles="dashdot")
    # draw the +/- 2.0 sigma band
    ax.fill_between([-0.5, len(pulls) - 0.5], [-2, -2], [2, 2], facecolor="yellow")
    # drawe the +/- 1.0 sigma band
    ax.fill_between([-0.5, len(pulls) - 0.5], [-1, -1], [1, 1], facecolor="green")
    # draw a horizontal line at pull=0.0
    ax.hlines([0], -0.5, len(pulls) - 0.5, colors="black", linestyles="dashed")
    # finally draw the pulls
    ax.scatter(range(len(pulls)), pulls, color="black")
    # and their uncertainties
    ax.errorbar(
        range(len(pulls)),
        pulls,
        color="black",
        xerr=0,
        yerr=pullerr,
        marker=".",
        fmt="none",
    )
    fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/fit_pull_plot_bdt1_{0}_bdt2_{1}_seed_{2}.pdf'.format( round(bdt1,1), round(bdt2,1), random_seed))

if __name__ == "__main__":
    main(sys.argv)