import sys
import pyhf
import json
import matplotlib.pyplot as plt
from IPython.display import display, HTML
import numpy as np
import matplotlib.ticker as mticker
from pyhf.contrib.viz import brazil
from uncertainties import ufloat
import ROOT
import array
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]
    folder_name = argv[4]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    f_data = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data = f_data.Get("h_toy_data")
    nbins = h_data.GetNbinsX()
    mass = [h_data.GetBinCenter(i+1) for i in range(nbins)] # for uniform bins
    # print(mass)
    # print(nbins)

    with open("/panfs/felician/B2Ktautau/workflow/generate_fit_workspaces/{0}/workspace_bdt1_{1}_bdt2_{2}_seed_{3}.json".format( folder_name, bdt1, bdt2, random_seed )) as serialized:
        spec = json.load(serialized)

    workspace = pyhf.Workspace(spec)

    # Get model 
    model = workspace.model()

    # Get observations
    data = workspace.data(model)

    # Check if signal or background data is empty (because nothing passes the BDT cuts)
    # If yes, return a large value for BF; otherwise perform the fit normally
    f_fit_templates = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{0}/fit_templates_bdt1_{1}_bdt2_{2}_seed_{3}.root".format(folder_name,bdt1,bdt2,random_seed))
    h_data = f_fit_templates.Get("Data/data")
    h_sig = f_fit_templates.Get("Signal/nominal")
    h_bkg = f_fit_templates.Get("Background/nominal")

    if(h_data.GetEntries() == 0):
        np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed ), [0,0])
        np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed ), [np.inf, np.inf])
    
        fig = plt.figure(figsize=(6,6))
        fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
    
        fig1 = plt.figure(figsize=(6,6))
        fig1.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
    else:
        # Fitting
        ## maximum likelihood fit (produces the maximum likelihood = denominator of the likelihood ratio)
        optimizer = pyhf.optimize.minuit_optimizer(verbose=2)
        pyhf.set_backend('numpy', optimizer) # minuit returns errors and correlation matrix; we ned the errors to make pulls
        
        fit_result, res_obj = pyhf.infer.mle.fit(data, model, return_uncertainties=True, return_result_obj=True)

        # run MINOS
        # res_obj.minuit.minos()
        # minos_errors = res_obj.minuit.merrors
        # print(minos_errors)

        cov_matrix = res_obj.minuit.covariance
        num_rows, num_cols = cov_matrix.shape
        corr_matrix = np.zeros((num_rows,num_cols))

        for i in range(num_rows):
            for j in range(num_cols):
                corr_matrix[i][j] = cov_matrix[i][j]/np.sqrt( cov_matrix[i][i] *  cov_matrix[j][j] )

        # df = pd.DataFrame(corr_matrix)
        correlation_value = corr_matrix[0][1]
        print("BF_sig and Nbkg correlation = ", correlation_value)

        fit_pars, fit_errors = fit_result.T
        # print(fit_pars)
        # print(fit_errors)

        ## fitted POI value
        fit_poi = fit_pars[model.config.poi_index] # I need to save this for each toy (and make a pull plot later)
        fit_poi_error = fit_errors[model.config.poi_index] # I need to save this for each toy (and make a pull plot later)
        print("POI = ", fit_poi, " +/- ", fit_poi_error)

        nbkg = fit_pars[0]
        nbkg_err = fit_errors[0]
        print("N_bkg = ", nbkg, " +/- ", nbkg_err, " sqrt(N_bkg) = ", np.sqrt(nbkg))

        # x1, x2, y = res_obj.minuit.contour(0, 1)
        # fig = plt.figure(figsize=(6,6))
        # ax = Axes3D(fig, auto_add_to_figure=False)
        # fig.add_axes(ax)
        # cm = plt.get_cmap("coolwarm")
        # ax.scatter(x1, x2, y, marker='o', c=y, cmap=cm)
        # ax.set_xlabel('Nbkg', rotation=-15)
        # ax.set_ylabel('BF_sig', rotation=45)
        # ax.set_zlabel('Function')
        # ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed))
        # ax.set_box_aspect(None, zoom=0.9)
        # plt.tight_layout()
        # plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/function_contour_2D_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        # plt.clf()

        # plt.plot(x1, y[0])
        # plt.xlabel("Nbkg")
        # plt.ylabel("Function")
        # plt.title("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed))
        # plt.axvspan(nbkg-nbkg_err, nbkg+nbkg_err, alpha=0.3, color='blue')
        # plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/function_contour_Nbkg_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        # plt.clf()

        # plt.plot(x2, y[1])
        # plt.xlabel("BF_sig")
        # plt.ylabel("Function")
        # plt.title("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed))
        # plt.axvspan(fit_poi-fit_poi_error, fit_poi+fit_poi_error, alpha=0.3, color='blue')
        # plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/function_contour_BF_sig_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        # plt.clf()

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

        def plot(**par_settings):
            pars = pyhf.tensorlib.astensor(model.config.suggested_init())
            for k, v in par_settings.items():
                pars[par_name_dict[k]] = v

            mc_counts = get_mc_counts(pars)

            data_sig = mc_counts[1][0]
            data_bkg = mc_counts[0][0]
            all_data = data_sig+data_bkg

            x = np.arange(len(data_sig))
            # mass = [ bins[i] for i in range(len(x)) ]

            plt.step(mass, all_data, where='mid', color='blue', label="Total fit")
            plt.step(mass, data_sig, where='mid', color='green', label='Signal')
            plt.step(mass, data_bkg, where='mid', color='red', label='Comb. background')

            the_data = workspace.data(model, include_auxdata=False)
            the_data_err = np.sqrt(the_data)

            plt.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

        par_name_dict = {k: v["slice"].start for k, v in model.config.par_map.items()}
        plt.title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {random_seed}')
        
        f, ax = plt.subplots()
        plot(**{k: fit_pars[v] for k, v in par_name_dict.items()})
        plt.xlabel('m_B (MeV)')
        plt.ylabel('Entries / ({0} MeV)'.format((8000-4000)/nbins))
        plt.title('BDT1 = {0} | BDT2 = {1} | seed = {2}'.format(bdt1,bdt2,random_seed))
        param1 = ufloat(fit_poi, fit_poi_error)
        param2 = ufloat(nbkg, nbkg_err) 
        textstr = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1}".format(param1, param2)
        plt.text(0.55, 0.6, textstr, fontsize=10, transform=ax.transAxes)
        plt.legend()
        plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))

        # Save fit result
        np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed ), [fit_poi, fit_poi_error, nbkg, nbkg_err, correlation_value])

        # # Make pull plot (for every fit parameter)
        # pulls = pyhf.tensorlib.concatenate(
        #     [
        #         (fit_pars[model.config.par_slice(k)] - model.config.param_set(k).suggested_init)
        #         / model.config.param_set(k).width()
        #         for k in model.config.par_order
        #         if model.config.param_set(k).constrained
        #     ]
        # )

        # pullerr = pyhf.tensorlib.concatenate(
        #     [
        #         fit_errors[model.config.par_slice(k)] / model.config.param_set(k).width()
        #         for k in model.config.par_order
        #         if model.config.param_set(k).constrained
        #     ]
        # )

        # labels = np.arange(len(pulls))
        # fig, ax = plt.subplots()
        # fig.set_size_inches(20, 5)

        # # set up axes labeling, ranges, etc...
        # ax.xaxis.set_major_locator(mticker.FixedLocator(np.arange(labels.size).tolist()))
        # ax.set_xticklabels(labels, rotation=30, ha="right")
        # ax.set_xlim(-0.5, len(pulls) - 0.5)
        # ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed), fontsize=18)
        # ax.set_ylabel(r"$(\theta - \hat{\theta})\,/ \Delta \theta$", fontsize=18)

        # # draw the +/- 2.0 horizontal lines
        # ax.hlines([-2, 2], -0.5, len(pulls) - 0.5, colors="black", linestyles="dotted")
        # # draw the +/- 1.0 horizontal lines
        # ax.hlines([-1, 1], -0.5, len(pulls) - 0.5, colors="black", linestyles="dashdot")
        # # draw the +/- 2.0 sigma band
        # ax.fill_between([-0.5, len(pulls) - 0.5], [-2, -2], [2, 2], facecolor="yellow")
        # # drawe the +/- 1.0 sigma band
        # ax.fill_between([-0.5, len(pulls) - 0.5], [-1, -1], [1, 1], facecolor="green")
        # # draw a horizontal line at pull=0.0
        # ax.hlines([0], -0.5, len(pulls) - 0.5, colors="black", linestyles="dashed")
        # # finally draw the pulls
        # ax.scatter(range(len(pulls)), pulls, color="black")
        # # and their uncertainties
        # ax.errorbar(
        #     range(len(pulls)),
        #     pulls,
        #     color="black",
        #     xerr=0,
        #     yerr=pullerr,
        #     marker=".",
        #     fmt="none",
        # )
        # fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_pull_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))

        if(random_seed == 1000):
            print("############################################# CLS LIMIT CALCULATION #######################################################")
            # CLs limit: evaluate upper limit on parameter of interest
            poi_values = np.linspace(-np.abs(fit_poi), np.abs(fit_poi)+4*fit_poi_error, 50)
            obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q")
            print(f"Upper limit (obs): μ = {obs_limit:.6f}")
            print(f"Upper limit (exp): μ = {exp_limits[2]:.6f}")
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed ), [obs_limit, exp_limits[2]])

            fig, ax = plt.subplots()
            fig.set_size_inches(10.5, 7)
            ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2}".format(bdt1,bdt2,random_seed))
            ax.set_xlabel(r"$BF_{sig}$")
            ax.set_ylabel(r"$\mathrm{CL}_{s}$")
            brazil.plot_results(scan, results, test_size=0.1, ax=ax)
            textstr1 = f"Upper limit (exp): = {exp_limits[2]:.6f}"
            plt.text(0.6, 0.6, textstr1, fontsize=15, transform=ax.transAxes)
            fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))


            # print("Using toys")
            # toy_calc = pyhf.infer.calculators.ToyCalculator(data, model, test_stat="q", ntoys=1000, track_progress=False)
            # teststat = toy_calc.teststatistic(poi_test=1.0)
            # print(f"q = {teststat}")

            # inits = model.config.suggested_init()
            # bounds = model.config.suggested_bounds()
            # fixeds = model.config.suggested_fixed()
            # pyhf.infer.test_statistics.qmu(1.0, data, model, inits, bounds, fixeds)

            # sb_dist, b_dist = toy_calc.distributions(poi_test=1.0)

            # p_exp_sb, p_exp_b, p_exp_s = toy_calc.expected_pvalues(sb_dist, b_dist)
            # print(f"exp. CL_sb = {p_exp_sb}")
            # print(f"exp. CL_b = {p_exp_b}")
            # print(f"exp. CL_s = CL_sb / CL_b = {p_exp_s}")



if __name__ == "__main__":
    main(sys.argv)