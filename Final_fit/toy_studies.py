import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import ROOT
from uncertainties import ufloat
from uncertainties import unumpy

def do_toy_study(folder_name, bdt1, bdt2, ch):
    N     = 100
    nbins = 20

    fit_poi = np.zeros(N, dtype=float)
    fit_poi_error = np.zeros(N, dtype=float)

    nbkg = np.zeros(N, dtype=float)
    nbkg_error = np.zeros(N, dtype=float)

    nbkg_1 = np.zeros(N, dtype=float)
    nbkg_error_1 = np.zeros(N, dtype=float)

    nbkg_2 = np.zeros(N, dtype=float)
    nbkg_error_2 = np.zeros(N, dtype=float)

    nbkg_3 = np.zeros(N, dtype=float)
    nbkg_error_3 = np.zeros(N, dtype=float)

    nsig = np.zeros(N, dtype=float)
    nsig_error = np.zeros(N, dtype=float)

    nsig_1 = np.zeros(N, dtype=float)
    nsig_error_1 = np.zeros(N, dtype=float)

    nsig_2 = np.zeros(N, dtype=float)
    nsig_error_2 = np.zeros(N, dtype=float)

    nsig_3 = np.zeros(N, dtype=float)
    nsig_error_3 = np.zeros(N, dtype=float)

    # c_syst = np.zeros(N, dtype=float)
    # c_syst_err = np.zeros(N, dtype=float)

    # c_syst_1 = np.zeros(N, dtype=float)
    # c_syst_err_1 = np.zeros(N, dtype=float)

    # c_syst_2 = np.zeros(N, dtype=float)
    # c_syst_err_2 = np.zeros(N, dtype=float)

    # c_syst_3 = np.zeros(N, dtype=float)
    # c_syst_err_3 = np.zeros(N, dtype=float)

    for i in range(N):
        if(ch == 0):
            fit_result = np.load('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}.npy'.format( folder_name, bdt1, bdt2, i ))                
        else:
            fit_result = np.load('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/fit_result_bdt1_{1}_bdt2_{2}_seed_{3}_ch.npy'.format( folder_name, bdt1, bdt2, i ))

        fit_poi[i] = fit_result[0]
        fit_poi_error[i] = fit_result[1]

        if(fit_poi_error[i] == 0):
            if(ch == 0):
                fig = plt.figure(figsize=(6,6))
                fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/bias_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2))

                fig1 = plt.figure(figsize=(6,6))
                fig1.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/pull_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2))
                quit()
            else:
                fig = plt.figure(figsize=(6,6))
                fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/bias_bdt1_{1}_bdt2_{2}_ch.pdf'.format( folder_name, bdt1, bdt2))

                fig1 = plt.figure(figsize=(6,6))
                fig1.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/pull_bdt1_{1}_bdt2_{2}_ch.pdf'.format( folder_name, bdt1, bdt2))
                quit()

        if(ch == 0):
            nbkg[i] = fit_result[2]
            nbkg_error[i] = fit_result[3]

            nsig[i] = fit_result[4]
            nsig_error[i] = fit_result[5]

            # c_syst[i] = fit_result[6]
            # c_syst_err[i] = fit_result[7]
        else:
            nbkg_1[i] = fit_result[2]
            nbkg_error_1[i] = fit_result[3]

            nbkg_2[i] = fit_result[4]
            nbkg_error_2[i] = fit_result[5]

            nbkg_3[i] = fit_result[6]
            nbkg_error_3[i] = fit_result[7]

            nsig_1[i] = fit_result[8]
            nsig_error_1[i] = fit_result[9]

            nsig_2[i] = fit_result[10]
            nsig_error_2[i] = fit_result[11]

            nsig_3[i] = fit_result[12]
            nsig_error_3[i] = fit_result[13]

            # c_syst_1[i] = fit_result[14]
            # c_syst_err_1[i] = fit_result[15]

            # c_syst_2[i] = fit_result[16]
            # c_syst_err_2[i] = fit_result[17]

            # c_syst_3[i] = fit_result[18]
            # c_syst_err_3[i] = fit_result[19]

    fig, ax = plt.subplots(3, 3, figsize=(15, 15))

    # BF
    ax[0,0].hist(fit_poi, bins=nbins)
    ax[0,0].set_xlabel("$BF_{sig}$")
    ax[0,0].set_ylabel("Entries / {0} bins".format(nbins))
    ax[0,0].set_title("Branching fraction")
    ax[0,0].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[0,0].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(fit_poi), np.std(fit_poi) ),
            transform=ax[0,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # BF error
    ax[1,0].hist(fit_poi_error, bins=nbins)
    ax[1,0].set_xlabel("$BF_{sig}$ error")
    ax[1,0].set_ylabel("Entries / {0} bins".format(nbins))
    ax[1,0].set_title("Branching fraction error")
    ax[1,0].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(fit_poi_error), np.std(fit_poi_error) ),
            transform=ax[1,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # BF pull
    pull = (fit_poi - 0)/fit_poi_error
    (mu_pull, sigma_pull) = norm.fit(pull)
    n, bins, patches = ax[2,0].hist( pull, bins=nbins)
    xcenters = (bins[:-1] + bins[1:]) / 2
    ax[2,0].errorbar(xcenters, n, yerr=np.sqrt(n), ecolor='black', fmt='k.')
    y_pull = sum(n)*(bins[1] - bins[0])*norm.pdf(xcenters, mu_pull, sigma_pull)
    ax[2,0].plot(xcenters, y_pull, 'r-', linewidth=2)
    ax[2,0].set_xlabel("Pull: ($BF_{sig}$ - 0)/error")
    ax[2,0].set_ylabel("Entries / {0} bins".format(nbins))
    chi2 = np.sum( ((n - y_pull)**2) / y_pull) / (nbins - 2)
    ax[2,0].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[2,0].text(0.05, 0.95, "From histogram \n $\mu$ = {:.4f} \n $\sigma$ = {:.4f}".format( np.mean(pull), np.std(pull) ),
            transform=ax[2,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[2,0].text(0.05, 0.65, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull, sigma_pull/np.sqrt(N), sigma_pull, sigma_pull/np.sqrt(2*N)),
            transform=ax[2,0].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[2,0].set_title(f"Branching fraction pull: $\chi^2$/ndf = {chi2:.3g}")

    # N_sig
    ax[0,1].hist(nsig, bins=nbins)
    ax[0,1].set_xlabel("$N_{sig}$")
    ax[0,1].set_ylabel("Entries / {0} bins".format(nbins))
    ax[0,1].set_title("Signal yield")
    ax[0,1].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[0,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(nsig), np.std(nsig) ),
            transform=ax[0,1].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # N_sig error
    ax[1,1].hist(nsig_error, bins=nbins)
    ax[1,1].set_xlabel("$N_{sig}$ error")
    ax[1,1].set_ylabel("Entries / {0} bins".format(nbins))
    ax[1,1].set_title("Signal yield error")
    ax[1,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(nsig_error), np.std(nsig_error) ),
            transform=ax[1,1].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # N_sig pull
    pull_1 = (nsig - 0)/nsig_error
    (mu_pull_1, sigma_pull_1) = norm.fit(pull_1)
    n_1, bins_1, patches_1 = ax[2,1].hist( pull_1, bins=nbins)
    xcenters_1 = (bins_1[:-1] + bins_1[1:]) / 2
    ax[2,1].errorbar(xcenters_1, n_1, yerr=np.sqrt(n_1), ecolor='black', fmt='k.')
    y_pull_1 = sum(n_1)*(bins_1[1] - bins_1[0])*norm.pdf(xcenters_1, mu_pull_1, sigma_pull_1)
    ax[2,1].plot(xcenters_1, y_pull_1, 'r-', linewidth=2)
    ax[2,1].set_xlabel("Pull: ($N_{sig}$ - 0)/error")
    ax[2,1].set_ylabel("Entries / {0} bins".format(nbins))
    chi2_1 = np.sum( ((n_1 - y_pull_1)**2) / y_pull_1) / (nbins - 2)
    ax[2,1].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[2,1].text(0.05, 0.95, "From histogram \n $\mu$ = {:.4f} \n $\sigma$ = {:.4f}".format( np.mean(pull_1), np.std(pull_1) ),
            transform=ax[2,1].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[2,1].text(0.05, 0.65, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_1, sigma_pull_1/np.sqrt(N), sigma_pull_1, sigma_pull_1/np.sqrt(2*N)),
            transform=ax[2,1].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[2,1].set_title(f"Signal yield pull: $\chi^2$/ndf = {chi2_1:.3g}")

    # N_bkg
    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_ws = f_ws.Get("DecayTree")
    eps_ws_den = t_ws.GetEntries()
    eps_ws_num = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    eps_ws = eps_ws_num/eps_ws_den
    N_bkg = 19006409*eps_ws

    ax[0,2].hist(nbkg, bins=nbins)
    ax[0,2].set_xlabel("$N_{bkg}$")
    ax[0,2].set_ylabel("Entries / {0} bins".format(nbins))
    ax[0,2].set_title("Background yield")
    ax[0,2].axvline(N_bkg, color='black', linestyle='--', linewidth=2)
    ax[0,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(nbkg), np.std(nbkg) ),
            transform=ax[0,2].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    ax[1,2].hist(nbkg_error, bins=nbins)
    ax[1,2].set_xlabel("$N_{bkg}$ error")
    ax[1,2].set_ylabel("Entries / {0} bins".format(nbins))
    ax[1,2].set_title("Background yield error")
    ax[1,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.1e} \n $\sigma$ = {:.1e}".format( np.mean(nbkg_error), np.std(nbkg_error) ),
            transform=ax[1,2].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    pull_2 = (nbkg - N_bkg)/nbkg_error
    (mu_pull_2, sigma_pull_2) = norm.fit(pull_2)
    n_2, bins_2, patches_2 = ax[2,2].hist( pull_2, bins=nbins)
    xcenters_2 = (bins_2[:-1] + bins_2[1:]) / 2
    ax[2,2].errorbar(xcenters_2, n_2, yerr=np.sqrt(n_2), ecolor='black', fmt='k.')
    y_pull_2 = sum(n_2)*(bins_2[1] - bins_2[0])*norm.pdf( xcenters_2, mu_pull_2, sigma_pull_2)
    ax[2,2].plot(xcenters_2, y_pull_2, 'r-', linewidth=2)
    ax[2,2].set_xlabel("Pull: ($N_{bkg} - N_{bkg}^{expected})/error")
    ax[2,2].set_ylabel("Entries / {0} bins".format(nbins))
    chi2_2 = np.sum( ((n_2 - y_pull_2)**2) / y_pull_2) / (nbins - 2)
    ax[2,2].set_title(f"Background yield pull: $\chi^2$/ndf = {chi2_2:.3g}")
    ax[2,2].axvline(0, color='black', linestyle='--', linewidth=2)
    ax[2,2].text(0.05, 0.95, "From histogram \n $\mu$ = {:.4f} \n $\sigma$ = {:.4f}".format( np.mean(pull_2), np.std(pull_2) ),
            transform=ax[2,2].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax[2,2].text(0.05, 0.65, "From fit \n $\mu$ = {:.4f} $\pm$ {:.4f} \n $\sigma$ = {:.4f} $\pm$ {:.4f}".format( mu_pull_2, sigma_pull_2/np.sqrt(N), sigma_pull_2, sigma_pull_2/np.sqrt(2*N)),
            transform=ax[2,2].transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    fig.suptitle("Fit to all events \n BDT1 = {0} | BDT2 = {1}".format(bdt1,bdt2), fontsize=24)
    fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/fit_validation_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))
    plt.clf()

    # else:
    #     eps_ws_num_1 = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(bdt1,bdt2))
    #     eps_ws_num_2 = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(bdt1,bdt2))
    #     eps_ws_num_3 = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(bdt1,bdt2))

    #     eps_ws_1 = eps_ws_num_1/eps_ws_den
    #     eps_ws_2 = eps_ws_num_2/eps_ws_den
    #     eps_ws_3 = eps_ws_num_3/eps_ws_den

    #     N_bkg_1 = 2*19006409*eps_ws_1
    #     N_bkg_2 = 2*19006409*eps_ws_2
    #     N_bkg_3 = 2*19006409*eps_ws_3

    #     bias1_1 = nbkg_1 - N_bkg_1
    #     bias1_2 = nbkg_2 - N_bkg_2
    #     bias1_3 = nbkg_3 - N_bkg_3

    #     (mu_bias1_1, sigma_bias1_1) = norm.fit(bias1_1)
    #     (mu_bias1_2, sigma_bias1_2) = norm.fit(bias1_2)
    #     (mu_bias1_3, sigma_bias1_3) = norm.fit(bias1_3)

    #     fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    #     n2_1, bins2_1, patches2_1 = ax1.hist(bias1_1, bins=nbins, label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu_bias1_1,3), round(sigma_bias1_1/np.sqrt(N),3), round(sigma_bias1_1,3), round(sigma_bias1_1/np.sqrt(2*N),3) ))
    #     n2_2, bins2_2, patches2_2 = ax2.hist(bias1_2, bins=nbins, label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu_bias1_2,3), round(sigma_bias1_2/np.sqrt(N),3), round(sigma_bias1_2,3), round(sigma_bias1_2/np.sqrt(2*N),3) ))
    #     n2_3, bins2_3, patches2_3 = ax3.hist(bias1_3, bins=nbins, label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu_bias1_3,3), round(sigma_bias1_3/np.sqrt(N),3), round(sigma_bias1_3,3), round(sigma_bias1_3/np.sqrt(2*N),3) ))

    #     xcenters2_1 = (bins2_1[:-1] + bins2_1[1:]) / 2
    #     xcenters2_2 = (bins2_2[:-1] + bins2_2[1:]) / 2
    #     xcenters2_3 = (bins2_3[:-1] + bins2_3[1:]) / 2

    #     ax1.errorbar(xcenters2_1, n2_1, yerr=np.sqrt(n2_1), ecolor='black', fmt='k.')
    #     ax2.errorbar(xcenters2_2, n2_2, yerr=np.sqrt(n2_2), ecolor='black', fmt='k.')
    #     ax3.errorbar(xcenters2_3, n2_3, yerr=np.sqrt(n2_3), ecolor='black', fmt='k.')
    #     y_bias1_1 = sum(n2_1)*(bins2_1[1] - bins2_1[0])*norm.pdf( xcenters2_1, mu_bias1_1, sigma_bias1_1)
    #     y_bias1_2 = sum(n2_2)*(bins2_2[1] - bins2_2[0])*norm.pdf( xcenters2_2, mu_bias1_2, sigma_bias1_2)
    #     y_bias1_3 = sum(n2_3)*(bins2_3[1] - bins2_3[0])*norm.pdf( xcenters2_3, mu_bias1_3, sigma_bias1_3)
    #     ax1.plot(xcenters2_1, y_bias1_1, 'r-', linewidth=2)
    #     ax2.plot(xcenters2_2, y_bias1_2, 'r-', linewidth=2)
    #     ax3.plot(xcenters2_3, y_bias1_3, 'r-', linewidth=2)
    #     ax1.set_xlabel("Bias: $N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$")
    #     ax2.set_xlabel("Bias: $N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$")
    #     ax3.set_xlabel("Bias: $N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$")
    #     ax1.set_ylabel("Entries / {0} bins".format(nbins))
    #     ax2.set_ylabel("Entries / {0} bins".format(nbins))
    #     ax3.set_ylabel("Entries / {0} bins".format(nbins))
    #     chi2_2_1 = np.sum( ((n2_1 - y_bias1_1)**2) / y_bias1_1) / (nbins - 2)
    #     chi2_2_2 = np.sum( ((n2_2 - y_bias1_2)**2) / y_bias1_2) / (nbins - 2)
    #     chi2_2_3 = np.sum( ((n2_3 - y_bias1_3)**2) / y_bias1_3) / (nbins - 2)
    #     ax1.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_2_1:.3g} \n Channel 1")
    #     ax2.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_2_2:.3g} \n Channel 2")
    #     ax3.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_2_3:.3g} \n Channel 3")
    #     ax1.legend()
    #     ax2.legend()
    #     ax3.legend()
    #     plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/bias_nbkg_bdt1_{1}_bdt2_{2}_ch.pdf'.format( folder_name, bdt1, bdt2 ))
    #     plt.clf()

    #     pull1_1 = bias1_1/nbkg_error_1
    #     pull1_2 = bias1_2/nbkg_error_2
    #     pull1_3 = bias1_3/nbkg_error_3

    #     (mu_pull1_1, sigma_pull1_1) = norm.fit(pull1_1)
    #     (mu_pull1_2, sigma_pull1_2) = norm.fit(pull1_2)
    #     (mu_pull1_3, sigma_pull1_3) = norm.fit(pull1_3)

    #     fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    #     n3_1, bins3_1, patches3_1 = ax1.hist( pull1_1, bins=nbins, label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu_pull1_1,3), round(sigma_pull1_1/np.sqrt(N),3), round(sigma_pull1_1,3), round(sigma_pull1_1/np.sqrt(2*N),3) ))
    #     n3_2, bins3_2, patches3_2 = ax2.hist( pull1_2, bins=nbins, label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu_pull1_2,3), round(sigma_pull1_2/np.sqrt(N),3), round(sigma_pull1_2,3), round(sigma_pull1_2/np.sqrt(2*N),3) ))
    #     n3_3, bins3_3, patches3_3 = ax3.hist( pull1_3, bins=nbins, label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu_pull1_3,3), round(sigma_pull1_3/np.sqrt(N),3), round(sigma_pull1_3,3), round(sigma_pull1_3/np.sqrt(2*N),3) ))

    #     xcenters3_1 = (bins3_1[:-1] + bins3_1[1:]) / 2
    #     xcenters3_2 = (bins3_2[:-1] + bins3_2[1:]) / 2
    #     xcenters3_3 = (bins3_3[:-1] + bins3_3[1:]) / 2

    #     ax1.errorbar(xcenters3_1, n3_1, yerr=np.sqrt(n3_1), ecolor='black', fmt='k.')
    #     ax2.errorbar(xcenters3_2, n3_2, yerr=np.sqrt(n3_2), ecolor='black', fmt='k.')
    #     ax3.errorbar(xcenters3_3, n3_3, yerr=np.sqrt(n3_3), ecolor='black', fmt='k.')
    #     y_pull1_1 = sum(n3_1)*(bins3_1[1] - bins3_1[0])*norm.pdf( xcenters3_1, mu_pull1_1, sigma_pull1_1)
    #     y_pull1_2 = sum(n3_2)*(bins3_2[1] - bins3_2[0])*norm.pdf( xcenters3_2, mu_pull1_2, sigma_pull1_2)
    #     y_pull1_3 = sum(n3_3)*(bins3_3[1] - bins3_3[0])*norm.pdf( xcenters3_3, mu_pull1_3, sigma_pull1_3)
    #     ax1.plot(xcenters3_1, y_pull1_1, 'r-', linewidth=2)
    #     ax2.plot(xcenters3_2, y_pull1_2, 'r-', linewidth=2)
    #     ax3.plot(xcenters3_3, y_pull1_3, 'r-', linewidth=2)
    #     ax1.set_xlabel("Pull: ($N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$)/error")
    #     ax2.set_xlabel("Pull: ($N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$)/error")
    #     ax3.set_xlabel("Pull: ($N_{bkg} - 2*N_{RS}^{pre-BDT}*\epsilon_{WS}$)/error")
    #     ax1.set_ylabel("Entries / {0} bins".format(nbins))
    #     ax2.set_ylabel("Entries / {0} bins".format(nbins))
    #     ax3.set_ylabel("Entries / {0} bins".format(nbins))
    #     chi2_3_1 = np.sum( ((n3_1 - y_pull1_1)**2) / y_pull1_1) / (nbins - 2)
    #     chi2_3_2 = np.sum( ((n3_2 - y_pull1_2)**2) / y_pull1_2) / (nbins - 2)
    #     chi2_3_3 = np.sum( ((n3_3 - y_pull1_3)**2) / y_pull1_3) / (nbins - 2)
    #     ax1.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_3_1:.3g} \n Channel 1")
    #     ax2.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_3_2:.3g} \n Channel 2")
    #     ax3.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_3_3:.3g} \n Channel 3")
    #     ax1.legend()
    #     ax2.legend()
    #     ax3.legend()
    #     plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/pull_nbkg_bdt1_{1}_bdt2_{2}_ch.pdf'.format( folder_name, bdt1, bdt2 ))
    #     plt.clf()

    # ### Normsys ####
    # if(ch == 0):
    #     pull_4 = (c_syst-0)/c_syst_err
    #     (mu4, sigma4) = norm.fit(pull_4)

    #     n4, bins4, patches4 = plt.hist(pull_4, bins=nbins)
    #     xcenters4 = (bins4[:-1] + bins4[1:]) / 2

    #     y4 = sum(n4)*(bins4[1] - bins4[0])*norm.pdf( xcenters4, mu4, sigma4)
    #     chi2_4 = np.sum( ((n4 - y4)**2) / y4) / (nbins - 2)

    #     plt.errorbar(xcenters4, n4, yerr=np.sqrt(n4), ecolor='black', fmt='k.', label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu4,5), round(sigma4/np.sqrt(N),5), round(sigma4,5), round(sigma4/np.sqrt(2*N),5) ))
    #     plt.plot(xcenters4, y4, 'r-', linewidth=2)
    #     plt.xlabel("(c_syst - 0)/c_syst_err")
    #     plt.ylabel("Entries / {0} bins".format(nbins))
    #     plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_4:.3g} \n All events")
    #     plt.legend()
    #     plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/normsys_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))
    #     plt.clf()
    # else:
    #     pull_4_1 = c_syst_1
    #     pull_4_2 = c_syst_2
    #     pull_4_3 = c_syst_3

    #     # pull_4_1 = (c_syst_1-0)/c_syst_err_1
    #     # pull_4_2 = (c_syst_2-0)/c_syst_err_2
    #     # pull_4_3 = (c_syst_3-0)/c_syst_err_3

    #     (mu4_1, sigma4_1) = norm.fit(pull_4_1)
    #     (mu4_2, sigma4_2) = norm.fit(pull_4_2)
    #     (mu4_3, sigma4_3) = norm.fit(pull_4_3)

    #     fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    #     n4_1, bins4_1, patches4_1 = ax1.hist(pull_4_1, bins=nbins)
    #     n4_2, bins4_2, patches4_2 = ax2.hist(pull_4_2, bins=nbins)
    #     n4_3, bins4_3, patches4_3 = ax3.hist(pull_4_3, bins=nbins)

    #     xcenters4_1 = (bins4_1[:-1] + bins4_1[1:]) / 2
    #     xcenters4_2 = (bins4_2[:-1] + bins4_2[1:]) / 2
    #     xcenters4_3 = (bins4_3[:-1] + bins4_3[1:]) / 2

    #     y4_1 = sum(n4_1)*(bins4_1[1] - bins4_1[0])*norm.pdf( xcenters4_1, mu4_1, sigma4_1)
    #     y4_2 = sum(n4_2)*(bins4_2[1] - bins4_2[0])*norm.pdf( xcenters4_2, mu4_2, sigma4_2)
    #     y4_3 = sum(n4_3)*(bins4_3[1] - bins4_3[0])*norm.pdf( xcenters4_3, mu4_3, sigma4_3)

    #     chi2_4_1 = np.sum( ((n4_1 - y4_1)**2) / y4_1) / (nbins - 2)
    #     chi2_4_2 = np.sum( ((n4_2 - y4_2)**2) / y4_2) / (nbins - 2)
    #     chi2_4_3 = np.sum( ((n4_3 - y4_3)**2) / y4_3) / (nbins - 2)

    #     ax1.errorbar(xcenters4_1, n4_1, yerr=np.sqrt(n4_1), ecolor='black', fmt='k.', label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu4_1,5), round(sigma4_1/np.sqrt(N),5), round(sigma4_1,5), round(sigma4_1/np.sqrt(2*N),5) ))
    #     ax2.errorbar(xcenters4_2, n4_2, yerr=np.sqrt(n4_2), ecolor='black', fmt='k.', label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu4_2,5), round(sigma4_2/np.sqrt(N),5), round(sigma4_2,5), round(sigma4_2/np.sqrt(2*N),5) ))
    #     ax3.errorbar(xcenters4_3, n4_3, yerr=np.sqrt(n4_3), ecolor='black', fmt='k.', label="$\mu = $ {0} $\pm$ {1} \n $\sigma = $ {2} $\pm$ {3}".format( round(mu4_3,5), round(sigma4_3/np.sqrt(N),5), round(sigma4_3,5), round(sigma4_3/np.sqrt(2*N),5) ))

    #     ax1.plot(xcenters4_1, y4_1, 'r-', linewidth=2)
    #     ax2.plot(xcenters4_2, y4_2, 'r-', linewidth=2)
    #     ax3.plot(xcenters4_3, y4_3, 'r-', linewidth=2)

    #     ax1.set_xlabel("(c_syst - 0)/c_syst_err")
    #     ax2.set_xlabel("(c_syst - 0)/c_syst_err")
    #     ax3.set_xlabel("(c_syst - 0)/c_syst_err")

    #     ax1.set_ylabel("Entries / {0} bins".format(nbins))
    #     ax2.set_ylabel("Entries / {0} bins".format(nbins))
    #     ax3.set_ylabel("Entries / {0} bins".format(nbins))

    #     ax1.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_4_1:.3g} \n Channel 1")
    #     ax2.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_4_2:.3g} \n Channel 2")
    #     ax3.set_title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} | chi2/ndf = {chi2_4_3:.3g} \n Channel 3")

    #     ax1.legend()
    #     ax2.legend()
    #     ax3.legend()

    #     fig.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/normsys_bdt1_{1}_bdt2_{2}_ch.pdf'.format( folder_name, bdt1, bdt2 ))
    #     plt.clf()

    # #### Nsig + Nbkg == N ? ####    
    # if(ch == 0):
    #     Nsig = unumpy.uarray( nsig, nsig_error )
    #     Nbkg = unumpy.uarray( nbkg, nbkg_error )

    #     N_total = Nsig + Nbkg

    #     n_total = unumpy.nominal_values(N_total)
    #     n_total_error = unumpy.std_devs(N_total)

    #     n_total_real = np.zeros(N)
    #     for i in range(N):
    #         f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{0}/toy_data_bdt1_{1}_bdt2_{2}_seed_{3}.root".format(folder_name,bdt1,bdt2,i))
    #         h = f.Get("h_toy_data")
    #         n_total_real[i] = h.GetEntries()
    # else:
    #     Nsig_1 = unumpy.uarray( nsig_1, nsig_error_1 )
    #     Nsig_2 = unumpy.uarray( nsig_2, nsig_error_2 )
    #     Nsig_3 = unumpy.uarray( nsig_3, nsig_error_3 )

    #     Nbkg_1 = unumpy.uarray( nbkg_1, nbkg_error_1 )
    #     Nbkg_2 = unumpy.uarray( nbkg_2, nbkg_error_2 )
    #     Nbkg_3 = unumpy.uarray( nbkg_3, nbkg_error_3 )

    #     N_total = Nsig_1 + Nsig_2 + Nsig_3 + Nbkg_1 + Nbkg_2 + Nbkg_3

    #     n_total = unumpy.nominal_values(N_total)
    #     n_total_error = unumpy.std_devs(N_total)

    #     n_total_real = np.zeros(N)
    #     for i in range(N):
    #         f1 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{0}/toy_data_bdt1_{1}_bdt2_{2}_seed_{3}_ch1.root".format(folder_name,bdt1,bdt2,i))
    #         f2 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{0}/toy_data_bdt1_{1}_bdt2_{2}_seed_{3}_ch2.root".format(folder_name,bdt1,bdt2,i))
    #         f3 = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/generate_toy_data/{0}/toy_data_bdt1_{1}_bdt2_{2}_seed_{3}_ch3.root".format(folder_name,bdt1,bdt2,i))

    #         h1 = f1.Get("h_toy_data_1")
    #         h2 = f2.Get("h_toy_data_2")
    #         h3 = f3.Get("h_toy_data_3")

    #         n_total_real[i] = h1.GetEntries() + h2.GetEntries() + h3.GetEntries()

    # # plt.figure(figsize=(8,8))
    # plt.errorbar(n_total, [i for i in range(len(n_total))], xerr=n_total_error, ecolor='black', fmt='ko', label="Fitted")
    # plt.errorbar(n_total_real, [i for i in range(len(n_total))], fmt='b.', label="Generated")
    # plt.axvline(np.mean(n_total), color='red', linestyle='--', label="Fitted mean")
    # if(ch == 0):
    #     plt.axvline(N_bkg, color='red', linestyle='-', label="True")
    # else:
    #     plt.axvline(N_bkg_1+N_bkg_2+N_bkg_3, color='red', linestyle='-', label="True")

    # plt.xlabel("$N_{total} = N_{sig} + N_{bkg}$")
    # plt.ylabel("Toy data")
    # plt.legend()
    # if(ch == 0):
    #     plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} \n All events")
    #     plt.tight_layout()
    #     plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/ntotal_bdt1_{1}_bdt2_{2}.pdf'.format( folder_name, bdt1, bdt2 ))
    #     plt.clf()
    # else:
    #     plt.title(f"BDT1 = {bdt1:.4g} | BDT2 = {bdt2:.4g} \n Error categories")
    #     plt.tight_layout()
    #     plt.savefig('/panfs/felician/B2Ktautau/workflow/pyhf_fit_validation/{0}/ntotal_bdt1_{1}_bdt2_{2}_ch.pdf'.format( folder_name, bdt1, bdt2 ))
    #     plt.clf()

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    folder_name = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    print("############################## Toy study for all ##############################################")
    do_toy_study(folder_name, bdt1, bdt2, 0)

    # print("############################## Toy study for mass error categories ##############################################")
    # do_toy_study(folder_name, bdt1, bdt2, 1)

if __name__ == "__main__":
    main(sys.argv)