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

def do_fit(h_data, h_data_1, h_data_2, h_data_3, spec, ch, folder_name, bdt1, bdt2, random_seed):

    if(ch == 0):
        nbins = h_data.GetNbinsX()
        mass = [h_data.GetBinCenter(i+1) for i in range(nbins)] 
        edges = [h_data.GetBinLowEdge(i+1) for i in range(nbins+1)] 
    else:
        nbins_1 = h_data_1.GetNbinsX()
        nbins_2 = h_data_2.GetNbinsX()
        nbins_3 = h_data_3.GetNbinsX()

        mass_1 = [h_data_1.GetBinCenter(i+1) for i in range(nbins_1)] 
        mass_2 = [h_data_2.GetBinCenter(i+1) for i in range(nbins_2)] 
        mass_3 = [h_data_3.GetBinCenter(i+1) for i in range(nbins_3)] 

        edges_1 =  [h_data_1.GetBinLowEdge(i+1) for i in range(nbins_1+1)] 
        edges_2 =  [h_data_2.GetBinLowEdge(i+1) for i in range(nbins_2+1)] 
        edges_3 =  [h_data_3.GetBinLowEdge(i+1) for i in range(nbins_3+1)] 

    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    data = workspace.data(model)

    print("Data : ", h_data.GetEntries())
    print("Data 1 : ", h_data_1.GetEntries())
    print("Data 2 : ", h_data_2.GetEntries())
    print("Data 3 : ", h_data_3.GetEntries())

    # Check if signal or background data is empty (because nothing passes the BDT cuts)
    # If yes, return inifinity otherwise perform the fit normally
    if( ((ch == 0) and (h_data.GetEntries() < 10)) or (((h_data_1.GetEntries() < 10) or (h_data_2.GetEntries() < 10) or (h_data_3.GetEntries() < 10)) and (ch != 0)) ):
        if(ch == 0):
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed), [0, 0, 0, 0, 0, 0])
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed), [np.inf, np.inf])
        
            fig = plt.figure(figsize=(6,6))
            fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        
            fig1 = plt.figure(figsize=(6,6))
            fig1.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        
            fig2 = plt.figure(figsize=(6,6))
            fig2.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/pull_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        else:
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}_ch.npy'.format( folder_name, bdt1, bdt2, random_seed), [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}_ch.npy'.format( folder_name, bdt1, bdt2, random_seed), [np.inf, np.inf])
        
            fig = plt.figure(figsize=(6,6))
            fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_plot_bdt1_{1}_bdt2_{2}_seed_{3}_ch.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        
            fig1 = plt.figure(figsize=(6,6))
            fig1.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}_ch.pdf'.format( folder_name, bdt1, bdt2, random_seed))

            fig2 = plt.figure(figsize=(6,6))
            fig2.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/pull_plot_bdt1_{1}_bdt2_{2}_seed_{3}_ch.pdf'.format( folder_name, bdt1, bdt2, random_seed))
    else:
        # Fitting
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
        # print("fit parameters")
        # print(len(fit_pars))
        # print(fit_pars)
        # print(fit_errors)

        if(ch == 0):
            fit_poi = fit_pars[1] 
            fit_poi_error = fit_errors[1] 
            print("BF_sig = ", fit_poi, " +/- ", fit_poi_error)

            nbkg = fit_pars[0]
            nbkg_err = fit_errors[0]
            print("N_bkg = ", nbkg, " +/- ", nbkg_err)

            constant = fit_pars[2]
            constant_error = fit_errors[2]

            # constant_syst = fit_pars[3]
            # constant_syst_error = fit_errors[3]

            fit_poi_unc = ufloat(fit_poi, fit_poi_error)
            constant_unc = ufloat(constant, constant_error)
            nsig = constant_unc*fit_poi_unc

        else:
            nbkg_1 = fit_pars[0]
            nbkg_2 = fit_pars[3]
            nbkg_3 = fit_pars[5]
            fit_poi = fit_pars[1]

            nbkg_1_error = fit_errors[0]
            nbkg_2_error = fit_errors[3]
            nbkg_3_error = fit_errors[5]
            fit_poi_error = fit_errors[1]

            print("BF_sig = ", fit_poi, " +/- ", fit_poi_error)
            print("N_bkg_1  = ", nbkg_1, " +/- ", nbkg_1_error)
            print("N_bkg_2  = ", nbkg_2, " +/- ", nbkg_2_error)
            print("N_bkg_3  = ", nbkg_3, " +/- ", nbkg_3_error)

            constant_1 = fit_pars[2]
            constant_2 = fit_pars[4]
            constant_3 = fit_pars[6]

            constant_1_error = fit_errors[2]
            constant_2_error = fit_errors[4]
            constant_3_error = fit_errors[6]

            fit_poi_unc = ufloat(fit_poi, fit_poi_error)
            constant_1_unc = ufloat(constant_1, constant_1_error)
            constant_2_unc = ufloat(constant_2, constant_2_error)
            constant_3_unc = ufloat(constant_3, constant_3_error)

            nsig_1 = constant_1*fit_poi_unc
            nsig_2 = constant_2*fit_poi_unc
            nsig_3 = constant_3*fit_poi_unc

            # constant_syst_1 = fit_pars[7]
            # constant_syst_2 = fit_pars[8]
            # constant_syst_3 = fit_pars[9]

            # constant_syst_1_error = fit_errors[7]
            # constant_syst_2_error = fit_errors[8]
            # constant_syst_3_error = fit_errors[9]


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
            
            if(ch == 0):
                data_sig = mc_counts[1][0]
                data_bkg = mc_counts[0][0]
                all_data = data_sig+data_bkg

                x = np.arange(len(data_sig))

                f, ax = plt.subplots()

                plt.step(edges, np.append(all_data, all_data[-1]), where='post', color='blue', label="Total fit")
                plt.step(edges, np.append(data_sig, data_sig[-1]), where='post', color='green', label='Signal')
                plt.step(edges, np.append(data_bkg, data_bkg[-1]), where='post', color='red', label='Comb. background')

                the_data = workspace.data(model, include_auxdata=False)
                the_data_err = np.sqrt(the_data)

                plt.errorbar(mass, the_data, the_data_err, c="k", marker='.', linestyle='', zorder=99, label="Toy data")
                plt.title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {random_seed} \n All events')
                plt.xlabel('m_B (MeV)')
                plt.ylabel('Entries / ({0} MeV)'.format( round((8000-4000)/nbins, 1) ))
                plt.xlim(4000,8000)

                param1 = ufloat(fit_poi, fit_poi_error)
                param2 = ufloat(nbkg, nbkg_err) 
                textstr = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1} \n $N_{{sig}} = $ {2}".format(param1, param2, nsig)
                plt.text(0.55, 0.6, textstr, fontsize=10, transform=ax.transAxes)
                plt.legend()

            else:
                data_sig = mc_counts[1][0]
                data_sig_1 = data_sig[0:nbins_1]
                data_sig_2 = data_sig[nbins_1:nbins_1+nbins_2]
                data_sig_3 = data_sig[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

                data_bkg = mc_counts[0][0]
                data_bkg_1 = data_bkg[0:nbins_1]
                data_bkg_2 = data_bkg[nbins_1:nbins_1+nbins_2]
                data_bkg_3 = data_bkg[nbins_1+nbins_2:nbins_1+nbins_2+nbins_3]

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
                ax1.step(edges_1, np.append(all_data_1, all_data_1[-1]), where='post', color='blue', label="Total fit")
                ax1.step(edges_1, np.append(data_sig_1, data_sig_1[-1]), where='post', color='green', label='Signal')
                ax1.step(edges_1, np.append(data_bkg_1, data_bkg_1[-1]), where='post', color='red', label='Comb. background')
                ax1.errorbar(mass_1, the_data_1, the_data_err_1, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

                ax2.step(edges_2, np.append(all_data_2, all_data_2[-1]), where='post', color='blue', label="Total fit")
                ax2.step(edges_2, np.append(data_sig_2, data_sig_2[-1]), where='post', color='green', label='Signal')
                ax2.step(edges_2, np.append(data_bkg_2, data_bkg_2[-1]), where='post', color='red', label='Comb. background')
                ax2.errorbar(mass_2, the_data_2, the_data_err_2, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

                ax3.step(edges_3, np.append(all_data_3, all_data_3[-1]), where='post', color='blue', label="Total fit")
                ax3.step(edges_3, np.append(data_sig_3, data_sig_3[-1]), where='post', color='green', label='Signal')
                ax3.step(edges_3, np.append(data_bkg_3, data_bkg_3[-1]), where='post', color='red', label='Comb. background')
                ax3.errorbar(mass_3, the_data_3, the_data_err_3, c="k", marker='.', linestyle='', zorder=99, label="Toy data")

                ax1.set_title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {random_seed} \n Channel 1')
                ax2.set_title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {random_seed} \n Channel 2')
                ax3.set_title(f'BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | seed = {random_seed} \n Channel 3')

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

                textstr_1 = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1} \n $N_{{sig}} = $ {2}".format(the_poi, the_nbkg_1, nsig_1)
                textstr_2 = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1} \n $N_{{sig}} = $ {2}".format(the_poi, the_nbkg_2, nsig_2)
                textstr_3 = "$BF_{{sig}} = $ {0} \n $N_{{bkg}} = $ {1} \n $N_{{sig}} = $ {2}".format(the_poi, the_nbkg_3, nsig_3)

                ax1.text(0.55, 0.6, textstr_1, fontsize=12, transform=ax1.transAxes)
                ax2.text(0.55, 0.6, textstr_2, fontsize=12, transform=ax2.transAxes)
                ax3.text(0.55, 0.6, textstr_3, fontsize=12, transform=ax3.transAxes)

                ax1.legend()
                ax2.legend()
                ax3.legend()

        par_name_dict = {k: v["slice"].start for k, v in model.config.par_map.items()}

        # f, ax = plt.subplots()
        plot(**{k: fit_pars[v] for k, v in par_name_dict.items()})
        if(ch == 0):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        else:
            plt.tight_layout()
            plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/fit_plot_bdt1_{1}_bdt2_{2}_seed_{3}_ch.pdf'.format( folder_name, bdt1, bdt2, random_seed))

        # Save fit result
        if(ch == 0):
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed), [fit_poi, fit_poi_error, nbkg, nbkg_err, nsig.nominal_value, nsig.std_dev])
            # np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed), [fit_poi, fit_poi_error, nbkg, nbkg_err, nsig.nominal_value, nsig.std_dev, constant_syst, constant_syst_error])
        else:
            np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}_ch.npy'.format( folder_name, bdt1, bdt2, random_seed), [fit_poi, fit_poi_error, nbkg_1, nbkg_1_error, nbkg_2, nbkg_2_error, nbkg_3, nbkg_3_error, nsig_1.nominal_value, nsig_1.std_dev, nsig_2.nominal_value, nsig_2.std_dev, nsig_3.nominal_value, nsig_3.std_dev])
            # np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}_ch.npy'.format( folder_name, bdt1, bdt2, random_seed), [fit_poi, fit_poi_error, nbkg_1, nbkg_1_error, nbkg_2, nbkg_2_error, nbkg_3, nbkg_3_error, nsig_1.nominal_value, nsig_1.std_dev, nsig_2.nominal_value, nsig_2.std_dev, nsig_3.nominal_value, nsig_3.std_dev, constant_syst_1, constant_syst_1_error, constant_syst_2, constant_syst_2_error, constant_syst_3, constant_syst_3_error])

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
        # if(ch == 0):
        #     ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2} \n All events".format(bdt1,bdt2,random_seed), fontsize=18)
        # else:
        #     ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2} \n Error categories".format(bdt1,bdt2,random_seed), fontsize=18)
        # ax.set_ylabel(r"$(\theta - \hat{\theta})\,/ \Delta \theta$", fontsize=18)

        # ax.hlines([-2, 2], -0.5, len(pulls) - 0.5, colors="black", linestyles="dotted")
        # ax.hlines([-1, 1], -0.5, len(pulls) - 0.5, colors="black", linestyles="dashdot")
        # ax.fill_between([-0.5, len(pulls) - 0.5], [-2, -2], [2, 2], facecolor="yellow")
        # ax.fill_between([-0.5, len(pulls) - 0.5], [-1, -1], [1, 1], facecolor="green")
        # ax.hlines([0], -0.5, len(pulls) - 0.5, colors="black", linestyles="dashed")
        # ax.scatter(range(len(pulls)), pulls, color="black")
        # ax.errorbar(
        #     range(len(pulls)),
        #     pulls,
        #     color="black",
        #     xerr=0,
        #     yerr=pullerr,
        #     marker=".",
        #     fmt="none",
        # )
        # if(ch == 0):
        #     fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/pull_plot_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
        # else:
        #     fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/pull_plot_bdt1_{1}_bdt2_{2}_seed_{3}_ch.pdf'.format( folder_name, bdt1, bdt2, random_seed))

        if(random_seed == 1000):
            print("### CLS LIMIT CALCULATION")
            # CLs limit: evaluate upper limit on parameter of interest
            poi_values = np.linspace(-np.abs(fit_poi), np.abs(fit_poi)+4*fit_poi_error, 50)
            obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, poi_values, level=0.1, return_results=True, test_stat="q")
            print(f"Upper limit (obs) ch {ch}: μ = {obs_limit:.6f}")
            print(f"Upper limit (exp) ch {ch}: μ = {exp_limits[2]:.6f}")
            if(ch == 0):
                np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, random_seed), [obs_limit, exp_limits[2]])
            else:
                np.save('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}_ch.npy'.format( folder_name, bdt1, bdt2, random_seed), [obs_limit, exp_limits[2]])

            fig, ax = plt.subplots()
            fig.set_size_inches(10.5, 7)
            if(ch == 0):
                ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2} \n All events".format(bdt1,bdt2,random_seed))
            else:
                ax.set_title("BDT1 = {0} | BDT2 = {1} | seed = {2} \n Error categories".format(bdt1,bdt2,random_seed))
            ax.set_xlabel(r"$BF_{sig}$")
            ax.set_ylabel(r"$\mathrm{CL}_{s}$")
            brazil.plot_results(scan, results, test_size=0.1, ax=ax)
            textstr1 = f"Upper limit (exp): = {exp_limits[2]:.6f}"
            plt.text(0.6, 0.6, textstr1, fontsize=15, transform=ax.transAxes)
            if(ch == 0):
                fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}.pdf'.format( folder_name, bdt1, bdt2, random_seed))
            else:
                fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit/plots/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_{3}_ch.pdf'.format( folder_name, bdt1, bdt2, random_seed))
    

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]
    folder_name = argv[4]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    ##### Get number of bins to build the mass array
    f_data = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data = f_data.Get("h_toy_data")

    f_data_1 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}_ch1.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data_1 = f_data_1.Get("h_toy_data_1")

    f_data_2 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}_ch2.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data_2 = f_data_2.Get("h_toy_data_2")

    f_data_3 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{3}/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}_ch3.root".format( bdt1, bdt2 ,random_seed, folder_name))
    h_data_3 = f_data_3.Get("h_toy_data_3")
    #####

    with open("/panfs/felician/B2Ktautau/workflow/generate_fit_workspaces/{0}/workspace_bdt1_{1}_bdt2_{2}_seed_{3}.json".format( folder_name, bdt1, bdt2, random_seed )) as serialized:
        spec = json.load(serialized)

    with open("/panfs/felician/B2Ktautau/workflow/generate_fit_workspaces/{0}/workspace_bdt1_{1}_bdt2_{2}_seed_{3}_ch.json".format( folder_name, bdt1, bdt2, random_seed )) as serialized_ch:
        spec_ch = json.load(serialized_ch)

    print("############################################# Fitting all ################################################")
    do_fit(h_data, h_data_1, h_data_2, h_data_3, spec, 0, folder_name, bdt1, bdt2, random_seed)

    print("############################################# Fitting with error categories ################################################")
    do_fit(h_data, h_data_1, h_data_2, h_data_3, spec_ch, 1, folder_name, bdt1, bdt2, random_seed)

            
if __name__ == "__main__":
    main(sys.argv)


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
