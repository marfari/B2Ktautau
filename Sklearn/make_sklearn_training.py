import random
import ROOT
import sys

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sb
from array import array

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier, GradientBoostingClassifier
from sklearn.gaussian_process.kernels import RBF

from sklearn.metrics import classification_report, roc_auc_score, roc_curve, auc
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.utils import compute_class_weight

import pickle
import xgboost as xgb

# Functions
def draw_variables(sig_data, bkg_data, mc_weights, columns, setup_name, step_name, species_name):
    nbins = 40

    hist_settings = {'bins': nbins, 'density': True, 'alpha': 0.4} # density = True means the histograms are normalised
    plt.figure(figsize=[15, 10])

    for id, column in enumerate(columns, 1):
        print(column)

        if(column == "VTXISODCHI2MASSONETRACK_B"):
            xlim = [3000,10000]
        elif(column == "VTXISODCHI2MASSTWOTRACK_B"):
            xlim = [0,20000]
        elif(column == "VTXISODCHI2ONETRACK_tau_min"):
            xlim = [0,25]
        elif(column == "VTXISODCHI2ONETRACK_B"):
            xlim = [0,1000]
        elif(column == "VTXISODCHI2TWOTRACK_B"):
            xlim = [0,5000]
        elif(column == "VTXISODCHI2TWOTRACK_tau_max"):
            xlim = [0,100]
        elif(column == "df_chi2"):
            xlim = [0,100]
        elif(column == "DV1_DV2_distance"):
            xlim = [0,20]
        elif(column == "tau_FD_BV_max"):
            xlim = [0,30]
        elif(column == "taum_pi2_IPCHI2_OWNPV"):
            xlim = [0,6000]
        elif(column == "DV1_DV2_distance_chi2"):
            xlim = [0,5]
        elif(column == "B_FD_PV_chi2"):
            xlim = [0,40000]
        elif(column == "D_FD_BV_min"):
            xlim = [0,150]
        elif(column == "Bp_VTXISODCHI2ONETRACK_B"):
            xlim = [0,5000]
        elif(column == "Bp_VTXISODCHI2ONETRACK_taup"):
            xlim = [0,100]
        elif(column == "Bp_VTXISODCHI2ONETRACK_taum"):
            xlim = [0,100]
        elif(column == "Bp_VTXISODCHI2TWOTRACK_taup"):
            xlim = [0,200]
        elif(column == "Bp_VTXISODCHI2TWOTRACK_taum"):
            xlim = [0,150]
        elif(column == "Bp_VTXISODCHI2ONETRACK_tau_max"):
            xlim = [0,150]
        elif(column == "Bp_VTXISODCHI2ONETRACK_tau_min"):
            xlim = [0,30]
        elif(column == "Bp_VTXISODCHI2TWOTRACK_tau_min"):
            xlim = [0,100]
        elif(column == "A_BV_DV1_DV2"):
            xlim = [0,2]
        elif(column == "A_PV_DV1_DV2"):
            xlim = [0,10]
        else:
            xlim = np.percentile(np.hstack([sig_data[column]]), [0.01, 99.99])
    
        # plt.hist(sig_data[column], color='b', weights=mc_weights, range=xlim, **hist_settings, label="Signal")
        plt.hist(sig_data[column], color='b', range=xlim, **hist_settings, label="Signal")
        plt.hist(bkg_data[column], color='r', range=xlim, **hist_settings, label="Background")
        plt.legend(loc="upper right", fontsize=20)
        plt.title(column, fontsize=25)
        plt.xlabel(column, fontsize=15)
        plt.ylabel("Normalised entries / ({0})".format(nbins), fontsize=15)

        plt.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/'+column+'.pdf')
        plt.clf()


def draw_scatter_plots(data, columns, setup_name, step_name, species_name):
    fig = sb.pairplot(data, hue="y", corner=True)
    fig._legend.set_title("")
    new_labels = ['Signal', 'Background']
    for t, l in zip(fig._legend.texts, new_labels):
        t.set_text(l)
    plt.setp(fig._legend.get_texts(), fontsize='32')

    fig.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/scatter_plot.pdf') 
    plt.clf()    


def correlation_matrix(data, name, setup_name, step_name, species_name, **kwds):
    """Calculate pairwise correlation between features.
    Extra arguments are passed on to DataFrame.corr()
    """
    # simply call df.corr() to get a table of
    # correlation values if you do not need
    # the fancy plotting
    corrmat = data.corr(**kwds)

    fig1, ax1 = plt.subplots(ncols=1, figsize=(6,5))
    
    opts = {'cmap': plt.get_cmap("RdBu"),
            'vmin': -1, 'vmax': +1}
    heatmap1 = ax1.pcolor(corrmat, **opts)
    plt.colorbar(heatmap1, ax=ax1)

    if(name == "sig"): 
        ax1.set_title("Signal")
    else:
        ax1.set_title("Background")

    labels = corrmat.columns.values
    for ax in (ax1,):
        # shift location of ticks to center of the bins
        ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
        ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
        ax.set_xticklabels(labels, minor=False, ha='right', rotation=70)
        ax.set_yticklabels(labels, minor=False)
        
    plt.tight_layout()
    fig1.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/correlation_matrix_'+name+'.pdf') 
    fig1.clf()


def do_cross_validation(name, clf, X_train, y_train, X_test, y_test):
    # Optimise hyper-parameters
    if(name == "AdaBDT"):
        param_grid = {'learning_rate': [0.05, 0.1, 0.5, 1.0], 'n_estimators': [50, 100, 500, 800]}
    elif(name == "GradBDT"):
        param_grid = {'learning_rate': [0.1, 0.5, 1.0], 'n_estimators': [100,500,800], 'max_depth': [1,3,6]}
    elif(name == "RForest"):
        param_grid = {'n_estimators': [50, 100, 500, 800], 'max_depth': [1,3,6,9]}
    elif(name == "XGBoost"):
        param_grid = {'learning_rate': [0.1, 0.5, 1.0], 'n_estimators': [100,500,800], 'max_depth': [1,3,6]}

    print("Performing grid search for hyper-paramter optimisation of "+name)
    # RandomizedSearchCV samples from distributions of the hyper-parameters and evaluates random points. It can be quicker than GridSearchCV
    grid_result = GridSearchCV(clf, param_grid, cv=3, scoring='roc_auc', n_jobs=6)
        
    if( (name=="AdaBDT") or (name=="GradBDT") or (name == "XGBoost") ):
        class_weights = compute_class_weight(class_weight='balanced', classes=np.unique(y_train), y=y_train)
        class_weight = np.zeros(len(X_train))

        # signal: y > 0.5 (class 1)
        # bakcground: y < 0.5 (class 0)
        for i in range(len(y_train)):
            if y_train[i] > 0.5:
                class_weight[i] = class_weights[1]
            else:
                class_weight[i] = class_weights[0]

        # print(class_weight)

        # grid_result.fit(X_train, y_train, sample_weight=X_dev_weights*class_weight)
        grid_result.fit(X_train, y_train, sample_weight=class_weight)
    else:
        # grid_result.fit(X_train, y_train, sample_weight=X_dev_weights)
        grid_result.fit(X_train, y_train)

    print("Best parameter set found on development set:")
    print(grid_result.best_estimator_)

    print("Grid scores on a subset of the development set:")
    cv_results = grid_result.cv_results_
    mean_score = cv_results['mean_test_score'] # mean f1 score from the (3) cross-validation splits
    std_score = cv_results['std_test_score'] # standard deviation of f1 score from the (3) cross validation splits
    params = cv_results['params'] # dictionary of the valur of the parameters for each combination in param_grid

    for i in range(len(mean_score)):
        print( "{0} +/- {1} for {2}".format(round(mean_score[i],4), round(std_score[i],4), params[i]) )

    if((name == "RForest") or (name == "XGBoost")):
        y_true, y_pred = y_test, grid_result.predict_proba(X_test)[:, 1]
    else:
        y_true, y_pred = y_test, grid_result.decision_function(X_test)
    print( "It scores {0} on the evaluation set".format(  round( roc_auc_score(y_true, y_pred), 4 ) ) )

    return grid_result.best_estimator_


def roc_curve_plot(classifiers, X_test_list, y_test_list, setup_name, step_name, species_name):
    clf_fpr = []
    clf_tpr = []
    clf_thresholds = []

    fpr_list = []
    tpr_list = []
    thresholds_list = []

    for i in range(len(classifiers)):
        decisions = classifiers[i].predict_proba(X_test_list[i])[:, 1]

        # Compute ROC curve and area under the curve
        fpr, tpr, thresholds = roc_curve(y_test_list[i], decisions)
        clf_fpr.append(fpr)
        clf_tpr.append(tpr)
        clf_thresholds.append(thresholds)

        roc_auc = auc(fpr, tpr)

        plt.plot(fpr, tpr, lw=1, label=f"fold = {i}"+' (area = %0.2f)'%(roc_auc))

        fpr_list.append(fpr)
        tpr_list.append(tpr)
        thresholds_list.append(thresholds)

    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Background efficiency')
    plt.ylabel('Signal efficiency')
    plt.title('ROC curve')
    plt.legend(loc="lower right")
    plt.grid()
    plt.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/roc_curve.pdf')
    plt.clf()

    return fpr_list, tpr_list, thresholds_list

# def bdt_2d_cut(X_test, y_test, thresholds_1st_step, thresholds_2nd_step, decisions_first_step, decisions_second_step):
#     # Add bdt columns to dataset
#     bdt_output_values_1 = pd.DataFrame({'bdt1': decisions_first_step})
#     bdt_output_values_2 = pd.DataFrame({'bdt2': decisions_second_step})

#     print(decisions_first_step)
#     print(decisions_second_step)

#     print(thresholds_1st_step)
#     print(thresholds_2nd_step)

#     metric = np.zeros((len(thresholds_1st_step),len(thresholds_2nd_step)))

#     print(X_test)
#     print(y_test)

#     # for bdt1 in thresholds_1st_step:
#     #     for bdt2 in thresholds_2nd_step:
#     #         eps_s = X_test[]


#     #     bdt_output_values = pd.DataFrame({'bdt': output})
#     #     bdt_output_values_signal =  bdt_output_values[0:signal.shape[0]]
#     #     bdt_output_values_background = bdt_output_values[signal.shape[0]: signal.shape[0] + background.shape[0]].reset_index()

#     #     print("before bdt cut")
#     #     print(signal)
#     #     print(background)

#     #     print(bdt_output_values_signal)
#     #     print(bdt_output_values_background)

#     #     signal = pd.concat([signal,bdt_output_values_signal], axis=1)
#     #     background = pd.concat([background,bdt_output_values_background], axis=1)


def draw_signal_significance_vs_bdt_cut(fold, fpr, tpr, thresholds, N_bkg, setup_name, step_name, species_name):
    eps_s = tpr
    B = fpr*N_bkg
    a = 5

    metric = eps_s/(a/2 + np.sqrt(B))

    plt.plot(thresholds, metric*20, color='g', label='Punzi (x20)')
    plt.plot(thresholds, tpr, color='b', label='Signal eff.')
    plt.plot(thresholds, fpr, color='r', label='Background eff.')
    plt.xlabel('BDT cut value')
    plt.title(f"fold = {fold}")
    plt.grid()
    plt.legend()

    plt.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/significance_vs_bdt_cut_fold_{fold}.pdf')
    plt.clf()

    optimal_index = np.argmax( np.ma.masked_invalid(metric) )
    optimal_metric = metric[optimal_index]
    optimal_cut = thresholds[optimal_index]
    print('The optimal cut value is {0} for a max significance of = {1}'.format(optimal_cut,optimal_metric))
    print('For this BDT cut, signal eff. = {0} and background eff. = {1}'.format(tpr[optimal_index], fpr[optimal_index]))

    return optimal_cut, optimal_metric


def compare_train_test(clf, X_train, y_train, X_test, y_test, setup_name, step_name, species_name, fold, bins=30):
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):
        d1 = clf.predict_proba(X[y>0.5])[:,1].ravel()
        d2 = clf.predict_proba(X[y<0.5])[:,1].ravel()
        decisions += [d1, d2]

    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    low_high = (low,high)
    
    plt.hist(decisions[0],
             color='b', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='S (train)')
    plt.hist(decisions[1],
             color='r', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='B (train)')

    hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='S (test)')
    
    hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='B (test)')

    # if(name == "RForest"): 
    #     twoclass_output = clf.predict_proba(X_train)
    # else:
    #     twoclass_output = clf.decision_function(X_train)
    # plot_range = (twoclass_output.min(), twoclass_output.max())
    # plt.hist(twoclass_output[y_train > 0.5], bins=bins, range=plot_range, facecolor='r', label="S (train)", alpha=0.5, edgecolor="k", density=True)
    # plt.hist(twoclass_output[y_train < 0.5], bins=bins, range=plot_range, facecolor='b', label="B (train)", alpha=0.5, edgecolor="k", density=True)
    # x1, x2, y1, y2 = plt.axis()

    h_sig_train = ROOT.TH1D("h_sig_train", "h_sig_train", 30, 0, 1)
    h_sig_test = ROOT.TH1D("h_sig_test", "h_sig_test", 30, 0, 1)
    h_bkg_train = ROOT.TH1D("h_bkg_train", "h_bkg_train", 30, 0, 1)
    h_bkg_test = ROOT.TH1D("h_bkg_test", "h_bkg_test", 30, 0, 1)

    for i in range(len(decisions[0])):
        h_sig_train.Fill(decisions[0][i])
    
    for i in range(len(decisions[1])):
        h_bkg_train.Fill(decisions[1][i])

    for i in range(len(decisions[2])):
        h_sig_test.Fill(decisions[2][i])

    for i in range(len(decisions[3])):
        h_bkg_test.Fill(decisions[3][i])
    
    h_sig_train.Scale(1.0/h_sig_train.Integral())
    h_sig_test.Scale(1.0/h_sig_test.Integral())
    h_bkg_train.Scale(1.0/h_bkg_train.Integral())
    h_bkg_test.Scale(1.0/h_bkg_test.Integral())

    ks_sig = h_sig_train.KolmogorovTest(h_sig_test)
    ks_bkg = h_bkg_train.KolmogorovTest(h_bkg_test)

    plt.xlabel("BDT output")
    plt.ylabel("Scaled entries")
    plt.title(f"fold = {fold} \n sig KS ({ks_sig:.4f}) | bkg KS ({ks_bkg:.4f})")
    plt.legend(loc='best')

    plt.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/classifier_output_fold_{fold}.pdf')
    plt.clf()


def confusion_matrix(name, clf, X_test, y_test, setup_name, step_name, species_name):
    fig = ConfusionMatrixDisplay.from_estimator(clf, X_test, y_test)
    lt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/{0}/{1}/{2}/confusion_matrix_'.format(species_name,setup_name,step_name)+name+'.pdf')
    plt.clf()


def draw_feature_importance(fold, clf, columns, setup_name, step_name, species_name):
    plt.figure(figsize=(20,10))
    importances = clf.feature_importances_
    sorted_idx = np.argsort(importances)
    sorted_columns = [columns[i] for i in sorted_idx]
    sorted_importances = importances[sorted_idx]
    plt.barh(sorted_columns, sorted_importances)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.title(f"Feature importance: fold = {fold}", fontsize=20)
    plt.tight_layout()

    plt.savefig(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/feature_importance_fold_{fold}.pdf')
    plt.clf()


def make_classification(sig_df, bkg_df, species_name, step_name, setup_name, cross_validation):
    signal = pd.DataFrame()
    background = pd.DataFrame()

    # Variables needed to build the input features for the 1st and 2nd steps
    if(species_name == "Ktautau"):
        branch_names = ['Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2ONETRACK_taup', 'Bp_VTXISODCHI2ONETRACK_taum', 'Bp_VTXISODCHI2TWOTRACK_taup', 'Bp_VTXISODCHI2TWOTRACK_taum', 'Bp_VTXISONUMVTX_taup', 'Bp_VTXISONUMVTX_taum', 
                        'Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup', 'Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum', 'Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup', 'Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum',
                        'Bp_TRKISOBDTFIRSTVALUE_taup_pi1', 'Bp_TRKISOBDTFIRSTVALUE_taup_pi2', 'Bp_TRKISOBDTFIRSTVALUE_taup_pi3', 
                        'Bp_TRKISOBDTSECONDVALUE_taup_pi1', 'Bp_TRKISOBDTSECONDVALUE_taup_pi2', 'Bp_TRKISOBDTSECONDVALUE_taup_pi3', 
                        'Bp_TRKISOBDTTHIRDVALUE_taup_pi1', 'Bp_TRKISOBDTTHIRDVALUE_taup_pi2', 'Bp_TRKISOBDTTHIRDVALUE_taup_pi3', 
                        'Bp_TRKISOBDTFIRSTVALUE_taum_pi1', 'Bp_TRKISOBDTFIRSTVALUE_taum_pi2', 'Bp_TRKISOBDTFIRSTVALUE_taum_pi3', 
                        'Bp_TRKISOBDTSECONDVALUE_taum_pi1', 'Bp_TRKISOBDTSECONDVALUE_taum_pi2', 'Bp_TRKISOBDTSECONDVALUE_taum_pi3', 
                        'Bp_TRKISOBDTTHIRDVALUE_taum_pi1', 'Bp_TRKISOBDTTHIRDVALUE_taum_pi2', 'Bp_TRKISOBDTTHIRDVALUE_taum_pi3',
                        'Bp_NC_05_PTASYM_taup', 'Bp_NC_05_PTASYM_taum', 'Bp_CCNC_05_IT_B', 'Bp_NC_05_PTASYM_taup', 'Bp_NC_05_PTASYM_taum', 'Bp_CC_05_MULT_B', 'Bp_NC_05_IT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_taup', 'Bp_CCNC_05_IT_taum',
                        'taup_M', 'taup_M12', 'taup_M23', 'taup_M13', 'taup_DIRA_ORIVX', 'taum_M', 'taum_M12', 'taum_M23', 'taum_M13', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                        'df_BVx', 'df_BVy', 'df_BVz', 'df_DV1x', 'df_DV1y', 'df_DV1z', 'df_DV2x', 'df_DV2y', 'df_DV2z', 'df_PVx', 'df_PVy', 'df_PVz',
                        'df_BVx_err', 'df_BVy_err', 'df_BVz_err', 'df_DV1x_err', 'df_DV1y_err', 'df_DV1z_err', 'df_DV2x_err', 'df_DV2y_err', 'df_DV2z_err', 'df_PVx_err', 'df_PVy_err', 'df_PVz_err',
                        'df_Kp_PX', 'df_Kp_PY', 'df_Kp_PZ', 'df_chi2', 
                        'Kp_IPCHI2_OWNPV', 'taup_pi1_IPCHI2_OWNPV', 'taup_pi2_IPCHI2_OWNPV', 'taup_pi3_IPCHI2_OWNPV', 'taum_pi1_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV', 'taum_pi3_IPCHI2_OWNPV',
                        'taup_AMAXDOCA', 'taup_AMINDOCA', 'taup_DOCACHI2MAX', 'taum_AMAXDOCA', 'taum_AMINDOCA', 'taum_DOCACHI2MAX',
                        'Bp_M02', 'Bp_M04', 'Bp_M06', 'Bp_M0456', 
                        'taup_ENDVERTEX_CHI2', 'taum_ENDVERTEX_CHI2', 'Bp_ENDVERTEX_CHI2', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_taum', 'Bp_VTXISODCHI2MASSTWOTRACK_taup',
                        'Bp_NC_05_PZASYM_taup', 'Bp_NC_05_PZASYM_taum', 'Bp_CC_05_PTASYM_taum', 'Bp_CC_05_PTASYM_taup', 
                        'eventNumber', 'runNumber', 'nCandidate', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_CC_05_PZASYM_B', 'Bp_NC_05_IT_taup', 'Bp_NC_05_IT_taum']

    else:
        branch_names = ['Bp_CC_05_IT_B',  'Bp_VTXISONUMVTX_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 
                        'Bp_VTXISONUMVTX_D0bar', 'Bp_VTXISODCHI2ONETRACK_D0bar', 'Bp_VTXISODCHI2MASSONETRACK_D0bar', 'Bp_VTXISODCHI2TWOTRACK_D0bar', 'Bp_VTXISODCHI2MASSTWOTRACK_D0bar', 
                        'Bp_VTXISONUMVTX_Dsp', 'Bp_VTXISODCHI2ONETRACK_Dsp', 'Bp_VTXISODCHI2MASSONETRACK_Dsp', 'Bp_VTXISODCHI2TWOTRACK_Dsp', 'Bp_VTXISODCHI2MASSTWOTRACK_Dsp',
                        'Bp_CC_05_DELTAPHI_B', 'Bp_NC_05_IT_B', 'Bp_CC_05_PTASYM_B', 'Bp_CC_05_PX_B', 'Bp_CC_05_PZASYM_B', 'Bp_CC_05_DELTAETA_B', 'Bp_NC_05_MULT_B', 'Bp_NC_05_VPT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_B', 'Bp_CC_05_MULT_B',
                        'Bp_CC_05_PASYM_D0bar', 'Bp_CC_05_DELTAETA_D0bar', 'Bp_NC_05_MULT_D0bar', 'Bp_NC_05_SPT_D0bar', 'Bp_CCNC_05_IT_D0bar', 'Bp_NC_05_PTASYM_D0bar', 'Bp_CC_05_PTASYM_D0bar', 'Bp_CC_05_IT_D0bar', 'Bp_TRKISOBDTSECONDVALUE_D0bar_K', 'Bp_TRKISOBDTSECONDVALUE_D0bar_pi', 'Bp_TRKISOBDTFIRSTVALUE_D0bar_K', 'Bp_TRKISOBDTFIRSTVALUE_D0bar_pi', 'Bp_TRKISOBDTTHIRDVALUE_D0bar_K', 'Bp_TRKISOBDTTHIRDVALUE_D0bar_pi',
                        'Bp_CC_05_PASYM_Dsp', 'Bp_CC_05_DELTAETA_Dsp', 'Bp_NC_05_MULT_Dsp', 'Bp_NC_05_SPT_Dsp', 'Bp_CCNC_05_IT_Dsp', 'Bp_NC_05_PTASYM_Dsp', 'Bp_CC_05_PTASYM_Dsp', 'Bp_CC_05_IT_Dsp', 'Bp_TRKISOBDTSECONDVALUE_Dsp_K1', 'Bp_TRKISOBDTSECONDVALUE_Dsp_K2', 'Bp_TRKISOBDTSECONDVALUE_Dsp_pi', 'Bp_TRKISOBDTFIRSTVALUE_Dsp_K1', 'Bp_TRKISOBDTFIRSTVALUE_Dsp_K2', 'Bp_TRKISOBDTFIRSTVALUE_Dsp_pi', 'Bp_TRKISOBDTTHIRDVALUE_Dsp_K1', 'Bp_TRKISOBDTTHIRDVALUE_Dsp_K2', 'Bp_TRKISOBDTTHIRDVALUE_Dsp_pi',
                        'Bp_Bstautau_ISOBDTFIRSTVALUE_taup', 'Bp_Bstautau_ISOBDTSECONDVALUE_taup', 'Bp_Bstautau_ISOBDTTHIRDVALUE_taup', 'Bp_Bstautau_ISOBDTFIRSTVALUE_taum', 
                        'Bp_Bstautau_ISOBDTSECONDVALUE_taum', 'Bp_Bstautau_ISOBDTTHIRDVALUE_taum',
                        'D0bar_M', 'D0bar_DIRA_ORIVX', 'Dsp_M', 'Dsp_M12', 'Dsp_M23', 'Dsp_M13', 'Dsp_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                        'Bp_ENDVERTEX_X', 'Bp_ENDVERTEX_Y', 'Bp_ENDVERTEX_Z', 'D0bar_ENDVERTEX_X', 'D0bar_ENDVERTEX_Y', 'D0bar_ENDVERTEX_Z', 'Dsp_ENDVERTEX_X', 'Dsp_ENDVERTEX_Y', 'Dsp_ENDVERTEX_Z', 'Bp_OWNPV_X', 'Bp_OWNPV_Y', 'Bp_OWNPV_Z',
                        'Bp_ENDVERTEX_XERR', 'Bp_ENDVERTEX_YERR', 'Bp_ENDVERTEX_ZERR', 'D0bar_ENDVERTEX_XERR', 'D0bar_ENDVERTEX_YERR', 'D0bar_ENDVERTEX_ZERR', 'Dsp_ENDVERTEX_XERR', 'Dsp_ENDVERTEX_YERR', 'Dsp_ENDVERTEX_ZERR', 'Bp_OWNPV_XERR', 'Bp_OWNPV_YERR', 'Bp_OWNPV_ZERR', 'Bp_dtf_chi2', 
                        'D0bar_AMAXDOCA', 'D0bar_AMINDOCA', 'D0bar_DOCACHI2MAX', 'Dsp_AMAXDOCA', 'Dsp_AMINDOCA', 'Dsp_DOCACHI2MAX', 'Bp_FDCHI2_OWNPV', 'D0bar_FD_ORIVX', 'Dsp_FD_ORIVX',
                        'D0bar_K_PX', 'D0bar_K_PY', 'D0bar_K_PZ', 'D0bar_K_PE', 'D0bar_pi_PX', 'D0bar_pi_PY', 'D0bar_pi_PZ', 'D0bar_pi_PE', 
                        'Dsp_K1_PX', 'Dsp_K1_PY', 'Dsp_K1_PZ', 'Dsp_K1_PE', 'Dsp_K2_PX', 'Dsp_K2_PY', 'Dsp_K2_PZ', 'Dsp_K2_PE', 'Dsp_pi_PX', 'Dsp_pi_PY', 'Dsp_pi_PZ', 'Dsp_pi_PE',
                        'eventNumber', 'runNumber', 'nCandidate', 'D0bar_ENDVERTEX_CHI2', 'Dsp_ENDVERTEX_CHI2']

    if(species_name == "Ktautau"):
        branch_names += ['df_Bp_M']
    else:
        branch_names += ['Bp_dtf_M']

    print("Loading signal")
    if(species_name == "Ktautau"):
        sig = sig_df.AsNumpy(branch_names+['Kp_ProbNNk_pidgen_default', 'taum_pi1_ProbNNpi_pidgen_default', 'taum_pi2_ProbNNpi_pidgen_default', 'taum_pi3_ProbNNpi_pidgen_default', 'taup_pi1_ProbNNpi_pidgen_default', 'taup_pi2_ProbNNpi_pidgen_default', 'taup_pi3_ProbNNpi_pidgen_default'])
        if((step_name == "combinatorial") or (step_name == "big_combinatorial")):
            bkg = bkg_df.AsNumpy(branch_names+['Kp_ProbNNk', 'taum_pi1_ProbNNpi', 'taum_pi2_ProbNNpi', 'taum_pi3_ProbNNpi', 'taup_pi1_ProbNNpi', 'taup_pi2_ProbNNpi', 'taup_pi3_ProbNNpi']) 
        else:
            bkg = bkg_df.AsNumpy(branch_names)
    else:
        sig = sig_df.AsNumpy(branch_names+['D0bar_K_ProbNNk_pidgen_default', 'Dsp_K1_ProbNNk_pidgen_default', 'Dsp_K2_ProbNNk_pidgen_default', 'D0bar_pi_ProbNNpi_pidgen_default', 'Dsp_pi_ProbNNpi_pidgen_default'])
        if((step_name == "combinatorial") or (step_name == "big_combinatorial")):
            bkg = bkg_df.AsNumpy(branch_names+['D0bar_K_ProbNNk', 'Dsp_K1_ProbNNk', 'Dsp_K2_ProbNNk', 'D0bar_pi_ProbNNpi', 'Dsp_pi_ProbNNpi']) 
        else:
            bkg = bkg_df.AsNumpy(branch_names)

    # Define the input feature
    if(species_name == "Ktautau"):
        if(step_name == "physics"):
            # SIGNAL
            signal['Bp_VTXISODCHI2ONETRACK_taup'] = sig['Bp_VTXISODCHI2ONETRACK_taup']
            signal['Bp_VTXISODCHI2ONETRACK_taum'] = sig['Bp_VTXISODCHI2ONETRACK_taum']
            signal['Bp_VTXISODCHI2TWOTRACK_taup'] = sig['Bp_VTXISODCHI2TWOTRACK_taup']
            signal['Bp_VTXISODCHI2TWOTRACK_taum'] = sig['Bp_VTXISODCHI2TWOTRACK_taum']
            # signal['Bp_VTXISODCHI2MASSTWOTRACK_taup'] = sig['Bp_VTXISODCHI2MASSTWOTRACK_taup']
            # signal['Bp_VTXISODCHI2MASSTWOTRACK_taum'] = sig['Bp_VTXISODCHI2MASSTWOTRACK_taum']
            signal['Bp_VTXISODCHI2ONETRACK_B'] = sig['Bp_VTXISODCHI2ONETRACK_B']
            # signal['Bp_VTXISODCHI2TWOTRACK_B'] = sig['Bp_VTXISODCHI2TWOTRACK_B']
            # signal['Bp_VTXISODCHI2MASSONETRACK_B'] = sig['Bp_VTXISODCHI2MASSONETRACK_B']
            signal['taup_iso_second_value'] = sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup']
            signal['taum_iso_second_value'] = sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum']
            signal['taup_iso_third_value'] = sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup']
            signal['taum_iso_third_value'] = sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum']
            # signal['Bp_CC_05_PZASYM_B'] = sig['Bp_CC_05_PZASYM_B']
            signal['Bp_NC_05_IT_B'] = sig['Bp_NC_05_IT_B']
            # signal['Bp_NC_05_IT_taup'] = sig['Bp_NC_05_IT_taup']
            # signal['Bp_NC_05_IT_taum'] = sig['Bp_NC_05_IT_taum']
            # signal['Bp_M02'] = sig['Bp_M02']
            # signal['Bp_M04'] = sig['Bp_M04']
            # signal['Bp_M06'] = sig['Bp_M06']

            a1_sig = np.sqrt( (sig['df_BVx'] - sig['df_DV1x'])**2 + (sig['df_BVy'] - sig['df_DV1y'])**2 + (sig['df_BVz'] - sig['df_DV1z'])**2 )
            b1_sig = np.sqrt( (sig['df_BVx'] - sig['df_DV2x'])**2 + (sig['df_BVy'] - sig['df_DV2y'])**2 + (sig['df_BVz'] - sig['df_DV2z'])**2 )
            c1_sig= np.sqrt( (sig['df_DV1x'] - sig['df_DV2x'])**2 + (sig['df_DV1y'] - sig['df_DV2y'])**2 + (sig['df_DV1z'] - sig['df_DV2z'])**2 )
            s1_sig = (a1_sig+b1_sig+c1_sig)*0.5
            signal['log10_A_BV_DV1_DV2'] = np.log10( np.sqrt(s1_sig*(s1_sig-a1_sig)*(s1_sig-b1_sig)*(s1_sig-c1_sig)) )

            # BACKGROUND
            background['Bp_VTXISODCHI2ONETRACK_taup'] = bkg['Bp_VTXISODCHI2ONETRACK_taup']
            background['Bp_VTXISODCHI2ONETRACK_taum'] = bkg['Bp_VTXISODCHI2ONETRACK_taum']
            background['Bp_VTXISODCHI2TWOTRACK_taup'] = bkg['Bp_VTXISODCHI2TWOTRACK_taup']
            background['Bp_VTXISODCHI2TWOTRACK_taum'] = bkg['Bp_VTXISODCHI2TWOTRACK_taum']
            # background['Bp_VTXISODCHI2MASSTWOTRACK_taup'] = bkg['Bp_VTXISODCHI2MASSTWOTRACK_taup']
            # background['Bp_VTXISODCHI2MASSTWOTRACK_taum'] = bkg['Bp_VTXISODCHI2MASSTWOTRACK_taum']
            background['Bp_VTXISODCHI2ONETRACK_B'] = bkg['Bp_VTXISODCHI2ONETRACK_B']
            # background['Bp_VTXISODCHI2TWOTRACK_B'] = bkg['Bp_VTXISODCHI2TWOTRACK_B']
            # background['Bp_VTXISODCHI2MASSONETRACK_B'] = bkg['Bp_VTXISODCHI2MASSONETRACK_B']
            background['taup_iso_second_value'] = bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup']
            background['taum_iso_second_value'] = bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum']
            background['taup_iso_third_value'] = bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup']
            background['taum_iso_third_value'] = bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum']
            # background['Bp_CC_05_PZASYM_B'] = bkg['Bp_CC_05_PZASYM_B']
            background['Bp_NC_05_IT_B'] = bkg['Bp_NC_05_IT_B']
            # background['Bp_NC_05_IT_taup'] = bkg['Bp_NC_05_IT_taup']
            # background['Bp_NC_05_IT_taum'] = bkg['Bp_NC_05_IT_taum']
            # background['Bp_M02'] = bkg['Bp_M02']
            # background['Bp_M04'] = bkg['Bp_M04']
            # background['Bp_M06'] = bkg['Bp_M06']

            a1_bkg = np.sqrt( (bkg['df_BVx'] - bkg['df_DV1x'])**2 + (bkg['df_BVy'] - bkg['df_DV1y'])**2 + (bkg['df_BVz'] - bkg['df_DV1z'])**2 )
            b1_bkg = np.sqrt( (bkg['df_BVx'] - bkg['df_DV2x'])**2 + (bkg['df_BVy'] - bkg['df_DV2y'])**2 + (bkg['df_BVz'] - bkg['df_DV2z'])**2 )
            c1_bkg= np.sqrt( (bkg['df_DV1x'] - bkg['df_DV2x'])**2 + (bkg['df_DV1y'] - bkg['df_DV2y'])**2 + (bkg['df_DV1z'] - bkg['df_DV2z'])**2 )
            s1_bkg = (a1_bkg+b1_bkg+c1_bkg)*0.5
            background['log10_A_BV_DV1_DV2'] = np.log10( np.sqrt(s1_bkg*(s1_bkg-a1_bkg)*(s1_bkg-b1_bkg)*(s1_bkg-c1_bkg)) )

        elif(step_name == "small_physics"):

            signal['tau_iso_second_value_min'] = np.minimum( sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] )
            signal['tau_iso_third_value_min'] = np.minimum( sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] ) 
            signal['Bp_NC_05_PTASYM_taup'] = sig['Bp_NC_05_PTASYM_taup']
            signal['Bp_NC_05_PTASYM_taum'] = sig['Bp_NC_05_PTASYM_taum']
            signal['Bp_NC_05_IT_B'] = sig['Bp_NC_05_IT_B']
            signal['Bp_M02'] = sig['Bp_M02']
            signal['Bp_M04'] = sig['Bp_M04']
            signal['Bp_M06'] = sig['Bp_M06']

            taup_FD_BV_sig = np.sqrt( (sig['df_DV1x'] - sig['df_BVx'])**2 + (sig['df_DV1y'] - sig['df_BVy'])**2 + (sig['df_DV1z'] - sig['df_BVz'])**2 )
            taum_FD_BV_sig = np.sqrt( (sig['df_DV2x'] - sig['df_BVx'])**2 + (sig['df_DV2y'] - sig['df_BVy'])**2 + (sig['df_DV2z'] - sig['df_BVz'])**2 )
            signal['tau_FD_BV_max'] = np.maximum(taup_FD_BV_sig, taum_FD_BV_sig)
      
            a_sig = np.sqrt( (sig['df_PVx'] - sig['df_DV1x'])**2 + (sig['df_PVy'] - sig['df_DV1y'])**2 + (sig['df_PVz'] - sig['df_DV1z'])**2 )
            b_sig = np.sqrt( (sig['df_PVx'] - sig['df_DV2x'])**2 + (sig['df_PVy'] - sig['df_DV2y'])**2 + (sig['df_PVz'] - sig['df_DV2z'])**2 )
            c_sig = np.sqrt( (sig['df_DV1x'] - sig['df_DV2x'])**2 + (sig['df_DV1y'] - sig['df_DV2y'])**2 + (sig['df_DV1z'] - sig['df_DV2z'])**2 )
            s_sig = 0.5*(a_sig+b_sig+c_sig)
            signal['A_DV1_DV2_PV_triangle'] = np.sqrt( s_sig*(s_sig-a_sig)*(s_sig-b_sig)*(s_sig-c_sig) )

            a1_sig = np.sqrt( (sig['df_BVx'] - sig['df_DV1x'])**2 + (sig['df_BVy'] - sig['df_DV1y'])**2 + (sig['df_BVz'] - sig['df_DV1z'])**2 )
            b1_sig = np.sqrt( (sig['df_BVx'] - sig['df_DV2x'])**2 + (sig['df_BVy'] - sig['df_DV2y'])**2 + (sig['df_BVz'] - sig['df_DV2z'])**2 )
            c1_sig = np.sqrt( (sig['df_DV1x'] - sig['df_DV2x'])**2 + (sig['df_DV1y'] - sig['df_DV2y'])**2 + (sig['df_DV1z'] - sig['df_DV2z'])**2 )
            s1_sig = 0.5*(a1_sig+b1_sig+c1_sig)
            signal['A_DV1_DV2_BV_triangle'] = np.sqrt( s1_sig*(s1_sig-a1_sig)*(s1_sig-b1_sig)*(s1_sig-c1_sig) )

            
            background['tau_iso_second_value_min'] = np.minimum( bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] )
            background['tau_iso_third_value_min'] = np.minimum( bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] ) 
            background['Bp_NC_05_PTASYM_taup'] = bkg['Bp_NC_05_PTASYM_taup']
            background['Bp_NC_05_PTASYM_taum'] = bkg['Bp_NC_05_PTASYM_taum']
            background['Bp_NC_05_IT_B'] = bkg['Bp_NC_05_IT_B']
            background['Bp_M02'] = bkg['Bp_M02']
            background['Bp_M04'] = bkg['Bp_M04']
            background['Bp_M06'] = bkg['Bp_M06']

            taup_FD_BV_bkg = np.sqrt( (bkg['df_DV1x'] - bkg['df_BVx'])**2 + (bkg['df_DV1y'] - bkg['df_BVy'])**2 + (bkg['df_DV1z'] - bkg['df_BVz'])**2 )
            taum_FD_BV_bkg = np.sqrt( (bkg['df_DV2x'] - bkg['df_BVx'])**2 + (bkg['df_DV2y'] - bkg['df_BVy'])**2 + (bkg['df_DV2z'] - bkg['df_BVz'])**2 )
            background['tau_FD_BV_max'] = np.maximum(taup_FD_BV_bkg, taum_FD_BV_bkg)
      
            a_bkg = np.sqrt( (bkg['df_PVx'] - bkg['df_DV1x'])**2 + (bkg['df_PVy'] - bkg['df_DV1y'])**2 + (bkg['df_PVz'] - bkg['df_DV1z'])**2 )
            b_bkg = np.sqrt( (bkg['df_PVx'] - bkg['df_DV2x'])**2 + (bkg['df_PVy'] - bkg['df_DV2y'])**2 + (bkg['df_PVz'] - bkg['df_DV2z'])**2 )
            c_bkg = np.sqrt( (bkg['df_DV1x'] - bkg['df_DV2x'])**2 + (bkg['df_DV1y'] - bkg['df_DV2y'])**2 + (bkg['df_DV1z'] - bkg['df_DV2z'])**2 )
            s_bkg = 0.5*(a_bkg+b_bkg+c_bkg)
            background['A_DV1_DV2_PV_triangle'] = np.sqrt( s_bkg*(s_bkg-a_bkg)*(s_bkg-b_bkg)*(s_bkg-c_bkg) )
    
            a1_bkg = np.sqrt( (bkg['df_BVx'] - bkg['df_DV1x'])**2 + (bkg['df_BVy'] - bkg['df_DV1y'])**2 + (bkg['df_BVz'] - bkg['df_DV1z'])**2 )
            b1_bkg = np.sqrt( (bkg['df_BVx'] - bkg['df_DV2x'])**2 + (bkg['df_BVy'] - bkg['df_DV2y'])**2 + (bkg['df_BVz'] - bkg['df_DV2z'])**2 )
            c1_bkg = np.sqrt( (bkg['df_DV1x'] - bkg['df_DV2x'])**2 + (bkg['df_DV1y'] - bkg['df_DV2y'])**2 + (bkg['df_DV1z'] - bkg['df_DV2z'])**2 )
            s1_bkg = 0.5*(a1_bkg+b1_bkg+c1_bkg)
            background['A_DV1_DV2_BV_triangle'] = np.sqrt( s1_bkg*(s1_bkg-a1_bkg)*(s1_bkg-b1_bkg)*(s1_bkg-c1_bkg) )

        elif(step_name == "combinatorial"):
            # SIGNAL:
            # signal['tau_M_max'] = np.maximum( sig['taup_M'], sig['taum_M'] )
            signal['tau_M12_M23_max_max'] = np.maximum( np.maximum( sig['taup_M12'], sig['taup_M23'] ), np.maximum( sig['taum_M12'], sig['taum_M23'] ) )
            # signal['tau_M12_M23_min_max'] = np.minimum( np.maximum( sig['taup_M12'], sig['taup_M23'] ), np.maximum( sig['taum_M12'], sig['taum_M23'] ) )
            # signal['tau_M12_M23_max_min'] = np.maximum( np.minimum( sig['taup_M12'], sig['taup_M23'] ), np.minimum( sig['taum_M12'], sig['taum_M23'] ) )
            signal['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taum_DIRA_ORIVX'] ) )
            signal['log10_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log10(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taum_DIRA_ORIVX'] ) )
            signal['log10_df_chi2'] = np.log10( sig['df_chi2'] )
            signal['Kp_ProbNNk'] = sig['Kp_ProbNNk_pidgen_default']
            signal['tau_prod_pi_min'] = np.minimum(sig['taum_pi1_ProbNNpi_pidgen_default']*sig['taum_pi2_ProbNNpi_pidgen_default']*sig['taum_pi3_ProbNNpi_pidgen_default'], sig['taup_pi1_ProbNNpi_pidgen_default']*sig['taup_pi2_ProbNNpi_pidgen_default']*sig['taup_pi3_ProbNNpi_pidgen_default'])
            signal['log10_Bp_FD_PV'] = np.log10( np.sqrt( (sig['df_BVx'] - sig['df_PVx'])**2 + (sig['df_BVy'] - sig['df_PVy'])**2 + (sig['df_BVz'] - sig['df_PVz'])**2 ) )
            signal['log10_tau_FD_BV_min'] = np.log10( np.minimum( np.sqrt( (sig['df_DV1x'] - sig['df_BVx'])**2 + (sig['df_DV1y'] - sig['df_BVy'])**2 + (sig['df_DV1z'] - sig['df_BVz'])**2), np.sqrt( (sig['df_DV2x'] - sig['df_BVx'])**2 + (sig['df_DV2y'] - sig['df_BVy'])**2 + (sig['df_DV2z'] - sig['df_BVz'])**2) ) )
            signal['log10_tau_ENDVERTEX_chi2_max'] = np.log10( np.maximum( sig['taup_ENDVERTEX_CHI2'], sig['taum_ENDVERTEX_CHI2'] ) )

            Cx_taup_sig =  (sig['df_DV1y'] - sig['df_BVy'])*sig['df_Kp_PZ']  - ( sig['df_DV1z'] - sig['df_BVz'])*sig['df_Kp_PY']
            Cy_taup_sig =  (sig['df_DV1z'] - sig['df_BVz'])*sig['df_Kp_PX']  - ( sig['df_DV1x'] - sig['df_BVx'])*sig['df_Kp_PZ']
            Cz_taup_sig =  (sig['df_DV1x'] - sig['df_BVx'])*sig['df_Kp_PY']  - ( sig['df_DV1y'] - sig['df_BVy'])*sig['df_Kp_PX']
            C_taup_sig = np.sqrt( Cx_taup_sig**2 + Cy_taup_sig**2 + Cz_taup_sig**2  )
            IP_taup_Kp_sig = (2*C_taup_sig)/( np.sqrt( sig['df_Kp_PX']**2 + sig['df_Kp_PY']**2 + sig['df_Kp_PZ']**2 ) )

            Cx_taum_sig =  (sig['df_DV2y'] - sig['df_BVy'])*sig['df_Kp_PZ']  - ( sig['df_DV2z'] - sig['df_BVz'])*sig['df_Kp_PY']
            Cy_taum_sig =  (sig['df_DV2z'] - sig['df_BVz'])*sig['df_Kp_PX']  - ( sig['df_DV2x'] - sig['df_BVx'])*sig['df_Kp_PZ']
            Cz_taum_sig =  (sig['df_DV2x'] - sig['df_BVx'])*sig['df_Kp_PY']  - ( sig['df_DV2y'] - sig['df_BVy'])*sig['df_Kp_PX']
            C_taum_sig = np.sqrt( Cx_taum_sig**2 + Cy_taum_sig**2 + Cz_taum_sig**2  )
            IP_taum_Kp_sig = (2*C_taum_sig)/( np.sqrt( sig['df_Kp_PX']**2 + sig['df_Kp_PY']**2 + sig['df_Kp_PZ']**2 ) )

            signal['log10_IP_tau_Kp_max'] = np.log10( np.maximum( IP_taup_Kp_sig, IP_taum_Kp_sig ) )
            signal['log10_IP_tau_Kp_min'] = np.log10( np.minimum( IP_taup_Kp_sig, IP_taum_Kp_sig ) )

            a_sig = np.sqrt( (sig['df_PVx'] - sig['df_DV1x'])**2 + (sig['df_PVy'] - sig['df_DV1y'])**2 + (sig['df_PVz'] - sig['df_DV1z'])**2 )
            b_sig = np.sqrt( (sig['df_PVx'] - sig['df_DV2x'])**2 + (sig['df_PVy'] - sig['df_DV2y'])**2 + (sig['df_PVz'] - sig['df_DV2z'])**2 )
            c_sig= np.sqrt( (sig['df_DV1x'] - sig['df_DV2x'])**2 + (sig['df_DV1y'] - sig['df_DV2y'])**2 + (sig['df_DV1z'] - sig['df_DV2z'])**2 )
            s_sig = (a_sig+b_sig+c_sig)*0.5
            signal['log10_A_PV_DV1_DV2'] = np.log10( np.sqrt(s_sig*(s_sig-a_sig)*(s_sig-b_sig)*(s_sig-c_sig)) )

            # BACKGROUND:
            # background['tau_M_max'] = np.maximum( bkg['taup_M'], bkg['taum_M'] )
            background['tau_M12_M23_max_max'] = np.maximum( np.maximum( bkg['taup_M12'], bkg['taup_M23'] ), np.maximum( bkg['taum_M12'], bkg['taum_M23'] ) )
            # background['tau_M12_M23_min_max'] = np.minimum( np.maximum( bkg['taup_M12'], bkg['taup_M23'] ), np.maximum( bkg['taum_M12'], bkg['taum_M23'] ) )
            # background['tau_M12_M23_max_min'] = np.maximum( np.minimum( bkg['taup_M12'], bkg['taup_M23'] ), np.minimum( bkg['taum_M12'], bkg['taum_M23'] ) )
            background['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taum_DIRA_ORIVX'] ) )
            background['log10_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log10(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taum_DIRA_ORIVX'] ) )
            background['log10_df_chi2'] = np.log10( bkg['df_chi2'] )
            background['Kp_ProbNNk'] = bkg['Kp_ProbNNk']
            background['tau_prod_pi_min'] = np.minimum(bkg['taum_pi1_ProbNNpi']*bkg['taum_pi2_ProbNNpi']*bkg['taum_pi3_ProbNNpi'], bkg['taup_pi1_ProbNNpi']*bkg['taup_pi2_ProbNNpi']*bkg['taup_pi3_ProbNNpi'])
            background['log10_Bp_FD_PV'] = np.log10( np.sqrt( (bkg['df_BVx'] - bkg['df_PVx'])**2 + (bkg['df_BVy'] - bkg['df_PVy'])**2 + (bkg['df_BVz'] - bkg['df_PVz'])**2 ) )
            background['log10_tau_FD_BV_min'] = np.log10( np.minimum( np.sqrt( (bkg['df_DV1x'] - bkg['df_BVx'])**2 + (bkg['df_DV1y'] - bkg['df_BVy'])**2 + (bkg['df_DV1z'] - bkg['df_BVz'])**2), np.sqrt( (bkg['df_DV2x'] - bkg['df_BVx'])**2 + (bkg['df_DV2y'] - bkg['df_BVy'])**2 + (bkg['df_DV2z'] - bkg['df_BVz'])**2) ) )
            background['log10_tau_ENDVERTEX_chi2_max'] = np.log10( np.maximum( bkg['taup_ENDVERTEX_CHI2'], bkg['taum_ENDVERTEX_CHI2'] ) )

            Cx_taup_bkg =  (bkg['df_DV1y'] - bkg['df_BVy'])*bkg['df_Kp_PZ']  - ( bkg['df_DV1z'] - bkg['df_BVz'])*bkg['df_Kp_PY']
            Cy_taup_bkg =  (bkg['df_DV1z'] - bkg['df_BVz'])*bkg['df_Kp_PX']  - ( bkg['df_DV1x'] - bkg['df_BVx'])*bkg['df_Kp_PZ']
            Cz_taup_bkg =  (bkg['df_DV1x'] - bkg['df_BVx'])*bkg['df_Kp_PY']  - ( bkg['df_DV1y'] - bkg['df_BVy'])*bkg['df_Kp_PX']
            C_taup_bkg = np.sqrt( Cx_taup_bkg**2 + Cy_taup_bkg**2 + Cz_taup_bkg**2  )
            IP_taup_Kp_bkg = (2*C_taup_bkg)/( np.sqrt( bkg['df_Kp_PX']**2 + bkg['df_Kp_PY']**2 + bkg['df_Kp_PZ']**2 ) )

            Cx_taum_bkg =  (bkg['df_DV2y'] - bkg['df_BVy'])*bkg['df_Kp_PZ']  - ( bkg['df_DV2z'] - bkg['df_BVz'])*bkg['df_Kp_PY']
            Cy_taum_bkg =  (bkg['df_DV2z'] - bkg['df_BVz'])*bkg['df_Kp_PX']  - ( bkg['df_DV2x'] - bkg['df_BVx'])*bkg['df_Kp_PZ']
            Cz_taum_bkg =  (bkg['df_DV2x'] - bkg['df_BVx'])*bkg['df_Kp_PY']  - ( bkg['df_DV2y'] - bkg['df_BVy'])*bkg['df_Kp_PX']
            C_taum_bkg = np.sqrt( Cx_taum_bkg**2 + Cy_taum_bkg**2 + Cz_taum_bkg**2  )
            IP_taum_Kp_bkg = (2*C_taum_bkg)/( np.sqrt( bkg['df_Kp_PX']**2 + bkg['df_Kp_PY']**2 + bkg['df_Kp_PZ']**2 ) )

            background['log10_IP_tau_Kp_max'] = np.log10( np.maximum( IP_taup_Kp_bkg, IP_taum_Kp_bkg ) )
            background['log10_IP_tau_Kp_min'] = np.log10( np.minimum( IP_taup_Kp_bkg, IP_taum_Kp_bkg ) )

            a_bkg = np.sqrt( (bkg['df_PVx'] - bkg['df_DV1x'])**2 + (bkg['df_PVy'] - bkg['df_DV1y'])**2 + (bkg['df_PVz'] - bkg['df_DV1z'])**2 )
            b_bkg = np.sqrt( (bkg['df_PVx'] - bkg['df_DV2x'])**2 + (bkg['df_PVy'] - bkg['df_DV2y'])**2 + (bkg['df_PVz'] - bkg['df_DV2z'])**2 )
            c_bkg= np.sqrt( (bkg['df_DV1x'] - bkg['df_DV2x'])**2 + (bkg['df_DV1y'] - bkg['df_DV2y'])**2 + (bkg['df_DV1z'] - bkg['df_DV2z'])**2 )
            s_bkg = (a_bkg+b_bkg+c_bkg)*0.5
            background['log10_A_PV_DV1_DV2'] = np.log10( np.sqrt(s_bkg*(s_bkg-a_bkg)*(s_bkg-b_bkg)*(s_bkg-c_bkg)) )

        elif(step_name == "big_combinatorial"):
            # SIGNAL:
            signal['Bp_FD_PV'] = np.sqrt( (sig['df_BVx'] - sig['df_PVx'])**2 + (sig['df_BVy'] - sig['df_PVy'])**2 + (sig['df_BVz'] - sig['df_PVz'])**2 )
            signal['Bp_VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( sig['Bp_VTXISODCHI2TWOTRACK_taup'], sig['Bp_VTXISODCHI2TWOTRACK_taum'] )
            signal['Bp_VTXISODCHI2ONETRACK_tau_min'] = np.minimum( sig['Bp_VTXISODCHI2ONETRACK_taup'], sig['Bp_VTXISODCHI2ONETRACK_taum'])
            signal['tau_iso_second_value_max'] = np.maximum( sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] )
            signal['tau_iso_third_value_max'] = np.maximum( sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] ) 
            signal['TRKISOBDTFIRSTVALUE_tau_pi_min_min'] = np.minimum( np.minimum( sig['Bp_TRKISOBDTFIRSTVALUE_taup_pi1'], sig['Bp_TRKISOBDTFIRSTVALUE_taup_pi2'], sig['Bp_TRKISOBDTFIRSTVALUE_taup_pi3'] ), np.minimum( sig['Bp_TRKISOBDTFIRSTVALUE_taum_pi1'], sig['Bp_TRKISOBDTFIRSTVALUE_taum_pi2'], sig['Bp_TRKISOBDTFIRSTVALUE_taum_pi3'] ) )
            signal['TRKISOBDTTHIRDVALUE_tau_pi_min_min'] = np.minimum( np.minimum( sig['Bp_TRKISOBDTTHIRDVALUE_taup_pi1'], sig['Bp_TRKISOBDTTHIRDVALUE_taup_pi2'], sig['Bp_TRKISOBDTTHIRDVALUE_taup_pi3'] ), np.minimum( sig['Bp_TRKISOBDTTHIRDVALUE_taum_pi1'], sig['Bp_TRKISOBDTTHIRDVALUE_taum_pi2'], sig['Bp_TRKISOBDTTHIRDVALUE_taum_pi3'] ) )
            signal['Bp_VTXISODCHI2ONETRACK_B'] = sig['Bp_VTXISODCHI2ONETRACK_B']
            signal['Bp_CC_05_MULT_B'] = sig['Bp_CC_05_MULT_B']
            signal['Bp_NC_05_PTASYM_B'] = sig['Bp_NC_05_PTASYM_B']
            signal['tau_M_max'] = np.maximum( sig['taup_M'], sig['taum_M'] )
            signal['tau_M12_M23_max_max'] = np.maximum( np.maximum( sig['taup_M12'], sig['taup_M23'] ), np.maximum( sig['taum_M12'], sig['taum_M23'] ) )
            signal['tau_M12_M23_min_max'] = np.minimum( np.maximum( sig['taup_M12'], sig['taup_M23'] ), np.maximum( sig['taum_M12'], sig['taum_M23'] ) )
            signal['tau_M12_M23_max_min'] = np.maximum( np.minimum( sig['taup_M12'], sig['taup_M23'] ), np.minimum( sig['taum_M12'], sig['taum_M23'] ) )
            signal['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taum_DIRA_ORIVX'] ) )
            signal['log10_df_chi2'] = np.log10( sig['df_chi2'] )
            signal['Kp_ProbNNk'] = sig['Kp_ProbNNk_pidgen_default']
            signal['tau_prod_pi_min'] = np.minimum(sig['taum_pi1_ProbNNpi_pidgen_default']*sig['taum_pi2_ProbNNpi_pidgen_default']*sig['taum_pi3_ProbNNpi_pidgen_default'], sig['taup_pi1_ProbNNpi_pidgen_default']*sig['taup_pi2_ProbNNpi_pidgen_default']*sig['taup_pi3_ProbNNpi_pidgen_default'])

            Cx_taup_sig =  (sig['df_DV1y'] - sig['df_BVy'])*sig['df_Kp_PZ']  - ( sig['df_DV1z'] - sig['df_BVz'])*sig['df_Kp_PY']
            Cy_taup_sig =  (sig['df_DV1z'] - sig['df_BVz'])*sig['df_Kp_PX']  - ( sig['df_DV1x'] - sig['df_BVx'])*sig['df_Kp_PZ']
            Cz_taup_sig =  (sig['df_DV1x'] - sig['df_BVx'])*sig['df_Kp_PY']  - ( sig['df_DV1y'] - sig['df_BVy'])*sig['df_Kp_PX']
            C_taup_sig = np.sqrt( Cx_taup_sig**2 + Cy_taup_sig**2 + Cz_taup_sig**2  )
            IP_taup_Kp_sig = (2*C_taup_sig)/( np.sqrt( sig['df_Kp_PX']**2 + sig['df_Kp_PY']**2 + sig['df_Kp_PZ']**2 ) )

            Cx_taum_sig =  (sig['df_DV2y'] - sig['df_BVy'])*sig['df_Kp_PZ']  - ( sig['df_DV2z'] - sig['df_BVz'])*sig['df_Kp_PY']
            Cy_taum_sig =  (sig['df_DV2z'] - sig['df_BVz'])*sig['df_Kp_PX']  - ( sig['df_DV2x'] - sig['df_BVx'])*sig['df_Kp_PZ']
            Cz_taum_sig =  (sig['df_DV2x'] - sig['df_BVx'])*sig['df_Kp_PY']  - ( sig['df_DV2y'] - sig['df_BVy'])*sig['df_Kp_PX']
            C_taum_sig = np.sqrt( Cx_taum_sig**2 + Cy_taum_sig**2 + Cz_taum_sig**2  )
            IP_taum_Kp_sig = (2*C_taum_sig)/( np.sqrt( sig['df_Kp_PX']**2 + sig['df_Kp_PY']**2 + sig['df_Kp_PZ']**2 ) )

            signal['IP_tau_Kp_max'] = np.maximum( IP_taup_Kp_sig, IP_taum_Kp_sig ) 
            signal['IP_tau_Kp_min'] = np.minimum( IP_taup_Kp_sig, IP_taum_Kp_sig ) 


            # BACKGROUND:
            background['Bp_FD_PV'] = np.sqrt( (bkg['df_BVx'] - bkg['df_PVx'])**2 + (bkg['df_BVy'] - bkg['df_PVy'])**2 + (bkg['df_BVz'] - bkg['df_PVz'])**2 )
            background['Bp_VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( bkg['Bp_VTXISODCHI2TWOTRACK_taup'], bkg['Bp_VTXISODCHI2TWOTRACK_taum'] )
            background['Bp_VTXISODCHI2ONETRACK_tau_min'] = np.minimum( bkg['Bp_VTXISODCHI2ONETRACK_taup'], bkg['Bp_VTXISODCHI2ONETRACK_taum'])
            background['tau_iso_second_value_max'] = np.maximum( bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] )
            background['tau_iso_third_value_max'] = np.maximum( bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] ) 
            background['TRKISOBDTFIRSTVALUE_tau_pi_min_min'] = np.minimum( np.minimum( bkg['Bp_TRKISOBDTFIRSTVALUE_taup_pi1'], bkg['Bp_TRKISOBDTFIRSTVALUE_taup_pi2'], bkg['Bp_TRKISOBDTFIRSTVALUE_taup_pi3'] ), np.minimum( bkg['Bp_TRKISOBDTFIRSTVALUE_taum_pi1'], bkg['Bp_TRKISOBDTFIRSTVALUE_taum_pi2'], bkg['Bp_TRKISOBDTFIRSTVALUE_taum_pi3'] ) )
            background['TRKISOBDTTHIRDVALUE_tau_pi_min_min'] = np.minimum( np.minimum( bkg['Bp_TRKISOBDTTHIRDVALUE_taup_pi1'], bkg['Bp_TRKISOBDTTHIRDVALUE_taup_pi2'], bkg['Bp_TRKISOBDTTHIRDVALUE_taup_pi3'] ), np.minimum( bkg['Bp_TRKISOBDTTHIRDVALUE_taum_pi1'], bkg['Bp_TRKISOBDTTHIRDVALUE_taum_pi2'], bkg['Bp_TRKISOBDTTHIRDVALUE_taum_pi3'] ) ) 
            background['Bp_VTXISODCHI2ONETRACK_B'] = bkg['Bp_VTXISODCHI2ONETRACK_B']
            background['Bp_CC_05_MULT_B'] = bkg['Bp_CC_05_MULT_B']
            background['Bp_NC_05_PTASYM_B'] = bkg['Bp_NC_05_PTASYM_B']
            background['Kp_ProbNNk'] = bkg['Kp_ProbNNk']
            background['tau_prod_pi_min'] = np.minimum(bkg['taum_pi1_ProbNNpi']*bkg['taum_pi2_ProbNNpi']*bkg['taum_pi3_ProbNNpi'], bkg['taup_pi1_ProbNNpi']*bkg['taup_pi2_ProbNNpi']*bkg['taup_pi3_ProbNNpi'])


            background['tau_M_max'] = np.maximum( bkg['taup_M'], bkg['taum_M'] )
            background['tau_M12_M23_max_max'] = np.maximum( np.maximum( bkg['taup_M12'], bkg['taup_M23'] ), np.maximum( bkg['taum_M12'], bkg['taum_M23'] ) )
            background['tau_M12_M23_min_max'] = np.minimum( np.maximum( bkg['taup_M12'], bkg['taup_M23'] ), np.maximum( bkg['taum_M12'], bkg['taum_M23'] ) )
            background['tau_M12_M23_max_min'] = np.maximum( np.minimum( bkg['taup_M12'], bkg['taup_M23'] ), np.minimum( bkg['taum_M12'], bkg['taum_M23'] ) )
            background['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taum_DIRA_ORIVX'] ) )
            background['log10_df_chi2'] = np.log10( bkg['df_chi2'] )

            Cx_taup_bkg =  (bkg['df_DV1y'] - bkg['df_BVy'])*bkg['df_Kp_PZ']  - ( bkg['df_DV1z'] - bkg['df_BVz'])*bkg['df_Kp_PY']
            Cy_taup_bkg =  (bkg['df_DV1z'] - bkg['df_BVz'])*bkg['df_Kp_PX']  - ( bkg['df_DV1x'] - bkg['df_BVx'])*bkg['df_Kp_PZ']
            Cz_taup_bkg =  (bkg['df_DV1x'] - bkg['df_BVx'])*bkg['df_Kp_PY']  - ( bkg['df_DV1y'] - bkg['df_BVy'])*bkg['df_Kp_PX']
            C_taup_bkg = np.sqrt( Cx_taup_bkg**2 + Cy_taup_bkg**2 + Cz_taup_bkg**2  )
            IP_taup_Kp_bkg = (2*C_taup_bkg)/( np.sqrt( bkg['df_Kp_PX']**2 + bkg['df_Kp_PY']**2 + bkg['df_Kp_PZ']**2 ) )

            Cx_taum_bkg =  (bkg['df_DV2y'] - bkg['df_BVy'])*bkg['df_Kp_PZ']  - ( bkg['df_DV2z'] - bkg['df_BVz'])*bkg['df_Kp_PY']
            Cy_taum_bkg =  (bkg['df_DV2z'] - bkg['df_BVz'])*bkg['df_Kp_PX']  - ( bkg['df_DV2x'] - bkg['df_BVx'])*bkg['df_Kp_PZ']
            Cz_taum_bkg =  (bkg['df_DV2x'] - bkg['df_BVx'])*bkg['df_Kp_PY']  - ( bkg['df_DV2y'] - bkg['df_BVy'])*bkg['df_Kp_PX']
            C_taum_bkg = np.sqrt( Cx_taum_bkg**2 + Cy_taum_bkg**2 + Cz_taum_bkg**2  )
            IP_taum_Kp_bkg = (2*C_taum_bkg)/( np.sqrt( bkg['df_Kp_PX']**2 + bkg['df_Kp_PY']**2 + bkg['df_Kp_PZ']**2 ) )

            background['IP_tau_Kp_max'] = np.maximum( IP_taup_Kp_bkg, IP_taum_Kp_bkg ) 
            background['IP_tau_Kp_min'] = np.minimum( IP_taup_Kp_bkg, IP_taum_Kp_bkg ) 

    else:
        mkaon = 493.677
        mpion = 139.57039
        if(step_name == "physics"):
            # SIGNAL
            signal['Bp_VTXISODCHI2ONETRACK_D0bar'] = sig['Bp_VTXISODCHI2ONETRACK_D0bar']
            signal['Bp_VTXISODCHI2ONETRACK_Dsp'] = sig['Bp_VTXISODCHI2ONETRACK_Dsp']
            signal['Bp_VTXISODCHI2TWOTRACK_D0bar'] = sig['Bp_VTXISODCHI2TWOTRACK_D0bar']
            signal['Bp_VTXISODCHI2TWOTRACK_Dsp'] = sig['Bp_VTXISODCHI2TWOTRACK_Dsp']
            signal['Bp_VTXISODCHI2MASSTWOTRACK_D0bar'] = sig['Bp_VTXISODCHI2MASSTWOTRACK_D0bar']
            signal['Bp_VTXISODCHI2MASSTWOTRACK_Dsp'] = sig['Bp_VTXISODCHI2MASSTWOTRACK_Dsp']
            signal['Bp_VTXISODCHI2ONETRACK_B'] = sig['Bp_VTXISODCHI2ONETRACK_B']
            signal['Bp_VTXISODCHI2MASSONETRACK_B'] = sig['Bp_VTXISODCHI2MASSONETRACK_B']
            signal['tau_iso_second_value_max'] = np.maximum( sig['Bp_Bstautau_ISOBDTSECONDVALUE_taup'], sig['Bp_Bstautau_ISOBDTSECONDVALUE_taum'] )
            signal['tau_iso_third_value_max'] = np.maximum( sig['Bp_Bstautau_ISOBDTTHIRDVALUE_taup'], sig['Bp_Bstautau_ISOBDTTHIRDVALUE_taum'] )
            signal['Bp_NC_05_IT_B'] = sig['Bp_NC_05_IT_B']

            M15_squared_sig = mkaon**2 + mpion**2 + 2*sig['D0bar_K_PE']*sig['Dsp_pi_PE'] - 2*( sig['D0bar_K_PX']*sig['Dsp_pi_PX'] + sig['D0bar_K_PY']*sig['Dsp_pi_PY'] + sig['D0bar_K_PZ']*sig['Dsp_pi_PZ'] )
            M23_squared_sig = mkaon**2 + mpion**2 + 2*sig['D0bar_pi_PE']*sig['Dsp_K1_PE'] - 2*( sig['D0bar_pi_PX']*sig['Dsp_K1_PX'] + sig['D0bar_pi_PY']*sig['Dsp_K1_PY'] + sig['D0bar_pi_PZ']*sig['Dsp_K1_PZ'] )
            M24_squared_sig = mkaon**2 + mpion**2 + 2*sig['D0bar_pi_PE']*sig['Dsp_K2_PE'] - 2*( sig['D0bar_pi_PX']*sig['Dsp_K2_PX'] + sig['D0bar_pi_PY']*sig['Dsp_K2_PY'] + sig['D0bar_pi_PZ']*sig['Dsp_K2_PZ'] )

            signal['Bp_M15'] = np.sqrt( np.abs(M15_squared_sig) )*np.sign(M15_squared_sig)
            signal['Bp_M23'] = np.sqrt( np.abs(M23_squared_sig) )*np.sign(M23_squared_sig)
            signal['Bp_M24'] = np.sqrt( np.abs(M24_squared_sig) )*np.sign(M24_squared_sig)

            a1_sig = np.sqrt( (sig['Bp_ENDVERTEX_X'] - sig['D0bar_ENDVERTEX_X'])**2 + (sig['Bp_ENDVERTEX_Y'] - sig['D0bar_ENDVERTEX_Y'])**2 + (sig['Bp_ENDVERTEX_Z'] - sig['D0bar_ENDVERTEX_Z'])**2 )
            b1_sig = np.sqrt( (sig['Bp_ENDVERTEX_X'] - sig['Dsp_ENDVERTEX_X'])**2 + (sig['Bp_ENDVERTEX_Y'] - sig['Dsp_ENDVERTEX_Y'])**2 + (sig['Bp_ENDVERTEX_Z'] - sig['Dsp_ENDVERTEX_Z'])**2 )
            c1_sig= np.sqrt( (sig['D0bar_ENDVERTEX_X'] - sig['Dsp_ENDVERTEX_X'])**2 + (sig['D0bar_ENDVERTEX_Y'] - sig['Dsp_ENDVERTEX_Y'])**2 + (sig['D0bar_ENDVERTEX_Z'] - sig['Dsp_ENDVERTEX_Z'])**2 )
            s1_sig = (a1_sig+b1_sig+c1_sig)*0.5
            signal['log10_A_BV_DV1_DV2'] = np.log10( np.sqrt(s1_sig*(s1_sig-a1_sig)*(s1_sig-b1_sig)*(s1_sig-c1_sig)) )

            # BACKGROUND
            background['Bp_VTXISODCHI2ONETRACK_D0bar'] = bkg['Bp_VTXISODCHI2ONETRACK_D0bar']
            background['Bp_VTXISODCHI2ONETRACK_Dsp'] = bkg['Bp_VTXISODCHI2ONETRACK_Dsp']
            background['Bp_VTXISODCHI2TWOTRACK_D0bar'] = bkg['Bp_VTXISODCHI2TWOTRACK_D0bar']
            background['Bp_VTXISODCHI2TWOTRACK_Dsp'] = bkg['Bp_VTXISODCHI2TWOTRACK_Dsp']
            background['Bp_VTXISODCHI2MASSTWOTRACK_D0bar'] = bkg['Bp_VTXISODCHI2MASSTWOTRACK_D0bar']
            background['Bp_VTXISODCHI2MASSTWOTRACK_Dsp'] = bkg['Bp_VTXISODCHI2MASSTWOTRACK_Dsp']
            background['Bp_VTXISODCHI2ONETRACK_B'] = bkg['Bp_VTXISODCHI2ONETRACK_B']
            background['Bp_VTXISODCHI2MASSONETRACK_B'] = bkg['Bp_VTXISODCHI2MASSONETRACK_B']
            background['tau_iso_second_value_max'] = np.maximum( bkg['Bp_Bstautau_ISOBDTSECONDVALUE_taup'], bkg['Bp_Bstautau_ISOBDTSECONDVALUE_taum'] )
            background['tau_iso_third_value_max'] = np.maximum( bkg['Bp_Bstautau_ISOBDTTHIRDVALUE_taup'], bkg['Bp_Bstautau_ISOBDTTHIRDVALUE_taum'] )
            background['Bp_NC_05_IT_B'] = bkg['Bp_NC_05_IT_B']

            M15_squared_bkg = mkaon**2 + mpion**2 + 2*bkg['D0bar_K_PE']*bkg['Dsp_pi_PE'] - 2*( bkg['D0bar_K_PX']*bkg['Dsp_pi_PX'] + bkg['D0bar_K_PY']*bkg['Dsp_pi_PY'] + bkg['D0bar_K_PZ']*bkg['Dsp_pi_PZ'] )
            M23_squared_bkg = mkaon**2 + mpion**2 + 2*bkg['D0bar_pi_PE']*bkg['Dsp_K1_PE'] - 2*( bkg['D0bar_pi_PX']*bkg['Dsp_K1_PX'] + bkg['D0bar_pi_PY']*bkg['Dsp_K1_PY'] + bkg['D0bar_pi_PZ']*bkg['Dsp_K1_PZ'] )
            M24_squared_bkg = mkaon**2 + mpion**2 + 2*bkg['D0bar_pi_PE']*bkg['Dsp_K2_PE'] - 2*( bkg['D0bar_pi_PX']*bkg['Dsp_K2_PX'] + bkg['D0bar_pi_PY']*bkg['Dsp_K2_PY'] + bkg['D0bar_pi_PZ']*bkg['Dsp_K2_PZ'] )

            background['Bp_M15'] = np.sqrt( np.abs(M15_squared_bkg) )*np.sign(M15_squared_bkg)
            background['Bp_M23'] = np.sqrt( np.abs(M23_squared_bkg) )*np.sign(M23_squared_bkg)
            background['Bp_M24'] = np.sqrt( np.abs(M24_squared_bkg) )*np.sign(M24_squared_bkg)

            a1_bkg = np.sqrt( (bkg['Bp_ENDVERTEX_X'] - bkg['D0bar_ENDVERTEX_X'])**2 + (bkg['Bp_ENDVERTEX_Y'] - bkg['D0bar_ENDVERTEX_Y'])**2 + (bkg['Bp_ENDVERTEX_Z'] - bkg['D0bar_ENDVERTEX_Z'])**2 )
            b1_bkg = np.sqrt( (bkg['Bp_ENDVERTEX_X'] - bkg['Dsp_ENDVERTEX_X'])**2 + (bkg['Bp_ENDVERTEX_Y'] - bkg['Dsp_ENDVERTEX_Y'])**2 + (bkg['Bp_ENDVERTEX_Z'] - bkg['Dsp_ENDVERTEX_Z'])**2 )
            c1_bkg= np.sqrt( (bkg['D0bar_ENDVERTEX_X'] - bkg['Dsp_ENDVERTEX_X'])**2 + (bkg['D0bar_ENDVERTEX_Y'] - bkg['Dsp_ENDVERTEX_Y'])**2 + (bkg['D0bar_ENDVERTEX_Z'] - bkg['Dsp_ENDVERTEX_Z'])**2 )
            s1_bkg = (a1_bkg+b1_bkg+c1_bkg)*0.5
            background['log10_A_BV_DV1_DV2'] = np.log10( np.sqrt(s1_bkg*(s1_bkg-a1_bkg)*(s1_bkg-b1_bkg)*(s1_bkg-c1_bkg)) )

        else:
            # SIGNAL:
            # signal['D_M_max'] = np.maximum( sig['D0bar_M'], sig['Dsp_M'] )
            # signal['log10_1_minus_D_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(sig['D0bar_DIRA_ORIVX'] ))*np.sign( sig['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['Dsp_DIRA_ORIVX'] ))*np.sign( sig['Dsp_DIRA_ORIVX'] ) )
            signal['log10_1_minus_D_DIRA_BV_max'] = np.maximum( np.log10(1 - np.abs(sig['D0bar_DIRA_ORIVX'] ))*np.sign( sig['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['Dsp_DIRA_ORIVX'] ))*np.sign( sig['Dsp_DIRA_ORIVX'] ) )

            dtf_chi2_sig = sig['Bp_dtf_chi2']
            for i in range(len(dtf_chi2_sig)):
                dtf_chi2_sig[i] = np.log10( dtf_chi2_sig[i][0] )
            # signal['log10_DTF_chi2'] = dtf_chi2_sig

            signal['D_prod_K_min'] = np.minimum(sig['D0bar_K_ProbNNk_pidgen_default'], sig['Dsp_K1_ProbNNk_pidgen_default']*sig['Dsp_K2_ProbNNk_pidgen_default'])
            signal['D_prod_pi_min'] = np.minimum(sig['D0bar_pi_ProbNNpi_pidgen_default'], sig['Dsp_pi_ProbNNpi_pidgen_default'])
            signal['log10_Bp_FD_PV'] = np.log10( np.sqrt( (sig['Bp_ENDVERTEX_X'] - sig['Bp_OWNPV_X'])**2 + (sig['Bp_ENDVERTEX_Y'] - sig['Bp_OWNPV_Y'])**2 + (sig['Bp_ENDVERTEX_Z'] - sig['Bp_OWNPV_Z'])**2 ) )
            signal['log10_D_FD_BV_min'] = np.log10( np.minimum( np.sqrt( (sig['D0bar_ENDVERTEX_X'] - sig['Bp_ENDVERTEX_X'])**2 + (sig['D0bar_ENDVERTEX_Y'] - sig['Bp_ENDVERTEX_Y'])**2 + (sig['D0bar_ENDVERTEX_Z'] - sig['Bp_ENDVERTEX_Z'])**2), np.sqrt( (sig['Dsp_ENDVERTEX_X'] - sig['Bp_ENDVERTEX_X'])**2 + (sig['Dsp_ENDVERTEX_Y'] - sig['Bp_ENDVERTEX_Y'])**2 + (sig['Dsp_ENDVERTEX_Z'] - sig['Bp_ENDVERTEX_Z'])**2) ) )
            signal['log10_tau_ENDVERTEX_chi2_max'] = np.log10( np.maximum( sig['D0bar_ENDVERTEX_CHI2'], sig['Dsp_ENDVERTEX_CHI2'] ) )

            a_sig = np.sqrt( (sig['Bp_OWNPV_X'] - sig['D0bar_ENDVERTEX_X'])**2 + (sig['Bp_OWNPV_Y'] - sig['D0bar_ENDVERTEX_Y'])**2 + (sig['Bp_OWNPV_Z'] - sig['D0bar_ENDVERTEX_Z'])**2 )
            b_sig = np.sqrt( (sig['Bp_OWNPV_X'] - sig['Dsp_ENDVERTEX_X'])**2 + (sig['Bp_OWNPV_Y'] - sig['Dsp_ENDVERTEX_Y'])**2 + (sig['Bp_OWNPV_Z'] - sig['Dsp_ENDVERTEX_Z'])**2 )
            c_sig= np.sqrt( (sig['D0bar_ENDVERTEX_X'] - sig['Dsp_ENDVERTEX_X'])**2 + (sig['D0bar_ENDVERTEX_Y'] - sig['Dsp_ENDVERTEX_Y'])**2 + (sig['D0bar_ENDVERTEX_Z'] - sig['Dsp_ENDVERTEX_Z'])**2 )
            s_sig = (a_sig+b_sig+c_sig)*0.5
            signal['log10_A_PV_DV1_DV2'] = np.log10( np.sqrt(s_sig*(s_sig-a_sig)*(s_sig-b_sig)*(s_sig-c_sig)) )

            # BACKGROUND
            # background['D_M_max'] = np.maximum( bkg['D0bar_M'], bkg['Dsp_M'] )
            # background['log10_1_minus_D_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(bkg['D0bar_DIRA_ORIVX'] ))*np.sign( bkg['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['Dsp_DIRA_ORIVX'] ))*np.sign( bkg['Dsp_DIRA_ORIVX'] ) )
            background['log10_1_minus_D_DIRA_BV_max'] = np.maximum( np.log10(1 - np.abs(bkg['D0bar_DIRA_ORIVX'] ))*np.sign( bkg['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['Dsp_DIRA_ORIVX'] ))*np.sign( bkg['Dsp_DIRA_ORIVX'] ) )

            dtf_chi2_bkg = bkg['Bp_dtf_chi2']
            for i in range(len(dtf_chi2_bkg)):
                dtf_chi2_bkg[i] = np.log10( dtf_chi2_bkg[i][0] )
            # background['log10_DTF_chi2'] = dtf_chi2_bkg

            background['D_prod_K_min'] = np.minimum(bkg['D0bar_K_ProbNNk'], bkg['Dsp_K1_ProbNNk']*bkg['Dsp_K2_ProbNNk'])
            background['D_prod_pi_min'] = np.minimum(bkg['D0bar_pi_ProbNNpi'], bkg['Dsp_pi_ProbNNpi'])
            background['log10_Bp_FD_PV'] = np.log10( np.sqrt( (bkg['Bp_ENDVERTEX_X'] - bkg['Bp_OWNPV_X'])**2 + (bkg['Bp_ENDVERTEX_Y'] - bkg['Bp_OWNPV_Y'])**2 + (bkg['Bp_ENDVERTEX_Z'] - bkg['Bp_OWNPV_Z'])**2 ) )
            background['log10_D_FD_BV_min'] = np.log10( np.minimum( np.sqrt( (bkg['D0bar_ENDVERTEX_X'] - bkg['Bp_ENDVERTEX_X'])**2 + (bkg['D0bar_ENDVERTEX_Y'] - bkg['Bp_ENDVERTEX_Y'])**2 + (bkg['D0bar_ENDVERTEX_Z'] - bkg['Bp_ENDVERTEX_Z'])**2), np.sqrt( (bkg['Dsp_ENDVERTEX_X'] - bkg['Bp_ENDVERTEX_X'])**2 + (bkg['Dsp_ENDVERTEX_Y'] - bkg['Bp_ENDVERTEX_Y'])**2 + (bkg['Dsp_ENDVERTEX_Z'] - bkg['Bp_ENDVERTEX_Z'])**2) ) )
            background['log10_tau_ENDVERTEX_chi2_max'] = np.log10( np.maximum( bkg['D0bar_ENDVERTEX_CHI2'], bkg['Dsp_ENDVERTEX_CHI2'] ) )

            a_bkg = np.sqrt( (bkg['Bp_OWNPV_X'] - bkg['D0bar_ENDVERTEX_X'])**2 + (bkg['Bp_OWNPV_Y'] - bkg['D0bar_ENDVERTEX_Y'])**2 + (bkg['Bp_OWNPV_Z'] - bkg['D0bar_ENDVERTEX_Z'])**2 )
            b_bkg = np.sqrt( (bkg['Bp_OWNPV_X'] - bkg['Dsp_ENDVERTEX_X'])**2 + (bkg['Bp_OWNPV_Y'] - bkg['Dsp_ENDVERTEX_Y'])**2 + (bkg['Bp_OWNPV_Z'] - bkg['Dsp_ENDVERTEX_Z'])**2 )
            c_bkg= np.sqrt( (bkg['D0bar_ENDVERTEX_X'] - bkg['Dsp_ENDVERTEX_X'])**2 + (bkg['D0bar_ENDVERTEX_Y'] - bkg['Dsp_ENDVERTEX_Y'])**2 + (bkg['D0bar_ENDVERTEX_Z'] - bkg['Dsp_ENDVERTEX_Z'])**2 )
            s_bkg = (a_bkg+b_bkg+c_bkg)*0.5
            background['log10_A_PV_DV1_DV2'] = np.log10( np.sqrt(s_bkg*(s_bkg-a_bkg)*(s_bkg-b_bkg)*(s_bkg-c_bkg)) )

    input_features = signal.columns.tolist()

    # k-fold
    k = 3
    signal['fold'] = (sig['eventNumber'] + sig['runNumber'] + sig['nCandidate']) % k
    background['fold'] = (bkg['eventNumber'] + bkg['runNumber'] + bkg['nCandidate']) % k

    # Add weights to the signal and background datasets
    mc_weights = []
    # mc_weights = sig_df.AsNumpy(['weight'])
    # mc_weights = pd.DataFrame(mc_weights)
    # data_weights = np.ones((bkg_df.Count()).GetValue())
    # data_weights = pd.DataFrame({'weight': data_weights})

    # print(mc_weights)
    # print(data_weights)

    # signal = pd.concat([signal,mc_weights], axis=1)
    # background = pd.concat([background,data_weights], axis=1) 

    # for sklearn data is usually organised into one 2D array of shape (n_samples x n_features) containing all the data and one array of categories of length n_samples
    # signal = 1 and background = 0
    X = np.concatenate((signal, background))
    y = np.concatenate((np.ones(signal.shape[0]), np.zeros(background.shape[0])))    

    # data = pd.DataFrame(np.hstack((X, y.reshape(y.shape[0], -1))),columns=input_features+['weight']+['y'])
    data = pd.DataFrame(np.hstack((X, y.reshape(y.shape[0], -1))), columns=input_features+['fold','y'])

    # Plotting variables and correlations
    if(setup_name == "input_features"):
        # Draw plots of the variables (signal vs background)
        print("Plotting input features")
        draw_variables(signal, background, mc_weights, input_features, setup_name, step_name, species_name)

        # data_scatter = data.drop('weight', axis=1)
        # draw_scatter_plots(data_scatter, input_features, step_name, species_name)
        # sig_corr = signal.drop('weight', axis=1)
        # bkg_corr = background.drop('weight', axis=1)

        sig_corr = signal
        bkg_corr = background
        if(species_name == "Ktautau"):
            sig_corr['Bp_M'] = sig['df_Bp_M']
            bkg_corr['Bp_M'] = bkg['df_Bp_M']
        else:
            dtf_M_sig = sig['Bp_dtf_M']
            for i in range(len(dtf_M_sig)):
                dtf_M_sig[i] = dtf_M_sig[i][0]

            dtf_M_bkg = bkg['Bp_dtf_M']
            for i in range(len(dtf_M_bkg)):
                dtf_M_bkg[i] = dtf_M_bkg[i][0]

            sig_corr['Bp_M'] = dtf_M_sig
            bkg_corr['Bp_M'] = dtf_M_bkg

        print("Plotting correlation matrix")
        correlation_matrix(signal, "sig", setup_name, step_name, species_name)
        correlation_matrix(background, "bkg", setup_name, step_name, species_name)

        quit()

    # Train and test split
    # X_dev,X_eval, y_dev,y_eval = train_test_split(X, y, test_size=0.1, random_state=42)
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    # X_train_weights = X_train[:, -1]
    # X_test_weights = X_test[:, -1]
    # X_train = np.delete(X_train, -1, 1)
    # X_test = np.delete(X_test, -1, 1)
    # X = np.delete(X, -1, 1)

    # print("Signal train size = ", X_train[y_train>0.5].shape[0] )
    # print("Signal test size = ", X_test[y_test>0.5].shape[0] )
    # print("Background train size = ", X_train[y_train<0.5].shape[0] )
    # print("Background test size = ", X_test[y_test<0.5].shape[0])

    ##############################################################################################################################################

    ###################################################### Train a classifier ####################################################################
    # I am using the optimised hyper-parameters here:
    name = "XGBoost"
    if(step_name == "physics"):
        clf = xgb.XGBClassifier(tree_method="hist", n_estimators=500, max_depth=2, learning_rate=0.1, random_state=42, importance_type='weight') # hist is the fastest tree method
    elif(step_name == "combinatorial"):
        clf = xgb.XGBClassifier(tree_method="hist", n_estimators=500, max_depth=2, learning_rate=0.1, random_state=42, importance_type='weight') # hist is the fastest tree method

    X = pd.DataFrame(X, columns=input_features+['fold'])
    y = pd.DataFrame(y, columns=['y'])
    df = pd.concat([X, y], axis=1)

    classifiers = []
    X_test_list = []
    y_test_list = []
    for fold in range(k):
        df_train = df[df['fold'] != fold]
        df_test = df[df['fold'] == fold]

        fold_train = df_train[['fold']]
        fold_test = df_test[['fold']]
        X_train = df_train[input_features]
        X_test = df_test[input_features]
        y_train = df_train[['y']]
        y_test = df_test[['y']]

        X_train = X_train.to_numpy()
        X_test = X_test.to_numpy()
        y_train = y_train.to_numpy().flatten()
        y_test = y_test.to_numpy().flatten()

        class_weights = compute_class_weight(class_weight='balanced', classes=np.unique(y_train), y=y_train)
        class_weight = np.zeros(len(X_train[:, -1]))

        # signal: y > 0.5 (class 1)
        # bakcground: y < 0.5 (class 0)
        for i in range(len(y_train)):
            if y_train[i] > 0.5:
                class_weight[i] = class_weights[1]
            else:
                class_weight[i] = class_weights[0]

        clf.fit(X_train, y_train, class_weight, eval_set=[(X_test, y_test)], verbose=False)

        if(cross_validation == "True"):
            # do_cross_validation(name, clf, X_train, y_train, X_test, y_test, X_train_weights) 
            do_cross_validation(name, clf, X_train, y_train, X_test, y_test) 
        else:
            print(f"Classification report: fold = {fold}")
            y_predicted = clf.predict(X_test)
            print(classification_report(y_test, y_predicted, target_names=["background", "signal"]))

            # print("Drawing confusion matrix")
            # confusion_matrix(name, clf, X_test, y_test, setup_name, step_name, species_name)

            print("Comparing train and test sets")
            compare_train_test(clf, X_train, y_train, X_test, y_test, setup_name, step_name, species_name, fold) 

        classifiers.append(clf)
        X_test_list.append(X_test)
        y_test_list.append(y_test)

    if(cross_validation == "True"):
        quit()


    # Compare different classifiers
    print("Drawing roc curve")
    fpr_list, tpr_list, thresholds_list = roc_curve_plot(classifiers, X_test_list, y_test_list, setup_name, step_name, species_name)

    print("Drawing signal significance vs BDT cut")
    for i in range(len(fpr_list)):
        x_test = X_test_list[i]
        N_bkg = x_test[y_test_list[i]<0.5].shape[0]
        bdt_cut, significance = draw_signal_significance_vs_bdt_cut(i, fpr_list[i], tpr_list[i], thresholds_list[i], N_bkg, setup_name, step_name, species_name)

    print("Drawing feature importance")
    for i in range(len(classifiers)):
        draw_feature_importance(i, classifiers[i], input_features, setup_name, step_name, species_name)

    # find which classifier gives the maximum value for the metric
    # idx = np.argmax( metrics )
    # if(names[idx] == "RForest"):
    #     decisions = classifiers[idx].predict_proba(X)[:, 1]
    # else:
    #     decisions = classifiers[idx].decision_function(X)

    # if(names[i] == "XGBoost"):
    #     classifiers[i].save_model('/panfs/felician/B2Ktautau/workflow/sklearn_training/{0}/clf_{1}.json'.format(species_name,step_name))
    # else:

    with open(f'/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/clf_{step_name}.pkl', 'wb') as f:
        for i in range(len(classifiers)):
            pickle.dump(classifiers[i], f)

def main(argv):

    species_name = argv[1]
    step_name = argv[2]
    setup_name = argv[3]
    makeCrossVal = argv[4] 

    ###################################################################### Retrieving files #####################################################################################
    if(species_name == "Ktautau"):
        # Prepare datasets
        # Signal proxy: All MC (fit region) (2016+2017+2018) (truth-matched - species 1)
        # pre-sel branches
        fc_sig_2016 = ROOT.TFileCollection("fc_sig_2016", "fc_sig_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1/pre_sel_tree.txt")
        fc_sig_2017 = ROOT.TFileCollection("fc_sig_2017", "fc_sig_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1/pre_sel_tree.txt")
        fc_sig_2018 = ROOT.TFileCollection("fc_sig_2018", "fc_sig_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1/pre_sel_tree.txt")

        # gsl branches
        fc_sig1_2016 = ROOT.TFileCollection("fc_sig1_2016", "fc_sig1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1/fit_results.txt")
        fc_sig1_2017 = ROOT.TFileCollection("fc_sig1_2017", "fc_sig1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1/fit_results.txt")
        fc_sig1_2018 = ROOT.TFileCollection("fc_sig1_2018", "fc_sig1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1/fit_results.txt")

        # invariant mass branches
        fc_sig2_2016 = ROOT.TFileCollection("fc_sig2_2016", "fc_sig2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1/invariant_mass_tree.txt")
        fc_sig2_2017 = ROOT.TFileCollection("fc_sig2_2017", "fc_sig2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1/invariant_mass_tree.txt")
        fc_sig2_2018 = ROOT.TFileCollection("fc_sig2_2018", "fc_sig2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1/invariant_mass_tree.txt")

        t_sig_2016 = ROOT.TChain("DecayTree")
        t_sig_2017 = ROOT.TChain("DecayTree")
        t_sig_2018 = ROOT.TChain("DecayTree")

        t_sig1_2016 = ROOT.TChain("DecayTree")
        t_sig1_2017 = ROOT.TChain("DecayTree")
        t_sig1_2018 = ROOT.TChain("DecayTree")

        t_sig2_2016 = ROOT.TChain("DecayTree")
        t_sig2_2017 = ROOT.TChain("DecayTree")
        t_sig2_2018 = ROOT.TChain("DecayTree")

        t_sig_2016.AddFileInfoList(fc_sig_2016.GetList())
        t_sig_2017.AddFileInfoList(fc_sig_2017.GetList())
        t_sig_2018.AddFileInfoList(fc_sig_2018.GetList())

        t_sig1_2016.AddFileInfoList(fc_sig1_2016.GetList())
        t_sig1_2017.AddFileInfoList(fc_sig1_2017.GetList())
        t_sig1_2018.AddFileInfoList(fc_sig1_2018.GetList())

        t_sig2_2016.AddFileInfoList(fc_sig2_2016.GetList())
        t_sig2_2017.AddFileInfoList(fc_sig2_2017.GetList())
        t_sig2_2018.AddFileInfoList(fc_sig2_2018.GetList())

        t_sig_2016.GetEntries()
        t_sig_2017.GetEntries()
        t_sig_2018.GetEntries()
        t_sig_2016.Add(t_sig_2017)
        t_sig_2016.Add(t_sig_2018)

        t_sig1_2016.GetEntries()
        t_sig1_2017.GetEntries()
        t_sig1_2018.GetEntries()
        t_sig1_2016.Add(t_sig1_2017)
        t_sig1_2016.Add(t_sig1_2018)

        t_sig2_2016.GetEntries()
        t_sig2_2017.GetEntries()
        t_sig2_2018.GetEntries()
        t_sig2_2016.Add(t_sig2_2017)
        t_sig2_2016.Add(t_sig2_2018)

        print("SIGNAL")
        print("pre sel entries = ", t_sig_2016.GetEntries())
        print("gsl entries = ", t_sig1_2016.GetEntries())
        print("inv mass entries = ", t_sig2_2016.GetEntries())

        t_sig_2016.AddFriend(t_sig1_2016)
        t_sig_2016.AddFriend(t_sig2_2016)

        if((step_name == "physics") or (step_name == "small_physics")):
            # Background proxy 1: B -> DD K+ (fit region) (2016+2017+2018) (B and final state particles truth-matched)
            # pre-sel branches
            fc_BuDDKp_2016 = ROOT.TFileCollection("fc_BuDDKp_2016", "fc_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_100/pre_sel_tree.txt")
            fc_BdDDKp_2016 = ROOT.TFileCollection("fc_BdDDKp_2016", "fc_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_110/pre_sel_tree.txt")
            fc_BsDDKp_2016 = ROOT.TFileCollection("fc_BsDDKp_2016", "fc_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_120/pre_sel_tree.txt")
            fc_BuDDK0_2016 = ROOT.TFileCollection("fc_BuDDK0_2016", "fc_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_130/pre_sel_tree.txt")
            fc_BuDD_2016 = ROOT.TFileCollection("fc_BuDD_2016", "fc_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_150/pre_sel_tree.txt")
            fc_BuDDKp_2017 = ROOT.TFileCollection("fc_BuDDKp_2017", "fc_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_100/pre_sel_tree.txt")
            fc_BdDDKp_2017 = ROOT.TFileCollection("fc_BdDDKp_2017", "fc_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_110/pre_sel_tree.txt")
            fc_BsDDKp_2017 = ROOT.TFileCollection("fc_BsDDKp_2017", "fc_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_120/pre_sel_tree.txt")
            fc_BuDDK0_2017 = ROOT.TFileCollection("fc_BuDDK0_2017", "fc_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_130/pre_sel_tree.txt")
            fc_BuDD_2017 = ROOT.TFileCollection("fc_BuDD_2017", "fc_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_150/pre_sel_tree.txt")
            fc_BuDDKp_2018 = ROOT.TFileCollection("fc_BuDDKp_2018", "fc_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_100/pre_sel_tree.txt")
            fc_BdDDKp_2018 = ROOT.TFileCollection("fc_BdDDKp_2018", "fc_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_110/pre_sel_tree.txt")
            fc_BsDDKp_2018 = ROOT.TFileCollection("fc_BsDDKp_2018", "fc_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_120/pre_sel_tree.txt")
            fc_BuDDK0_2018 = ROOT.TFileCollection("fc_BuDDK0_2018", "fc_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_130/pre_sel_tree.txt")
            fc_BuDD_2018 = ROOT.TFileCollection("fc_BuDD_2018", "fc_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_150/pre_sel_tree.txt")

            # gsl branches
            fc_BuDDKp1_2016 = ROOT.TFileCollection("fc_BuDDKp1_2016", "fc_BuDDKp1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_100/fit_results.txt")
            fc_BdDDKp1_2016 = ROOT.TFileCollection("fc_BdDDKp1_2016", "fc_BdDDKp1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_110/fit_results.txt")
            fc_BsDDKp1_2016 = ROOT.TFileCollection("fc_BsDDKp1_2016", "fc_BsDDKp1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_120/fit_results.txt")
            fc_BuDDK01_2016 = ROOT.TFileCollection("fc_BuDDK01_2016", "fc_BuDDK01_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_130/fit_results.txt")
            fc_BuDD1_2016 = ROOT.TFileCollection("fc_BuDD1_2016", "fc_BuDD1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_150/fit_results.txt")
            fc_BuDDKp1_2017 = ROOT.TFileCollection("fc_BuDDKp1_2017", "fc_BuDDKp1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_100/fit_results.txt")
            fc_BdDDKp1_2017 = ROOT.TFileCollection("fc_BdDDKp1_2017", "fc_BdDDKp1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_110/fit_results.txt")
            fc_BsDDKp1_2017 = ROOT.TFileCollection("fc_BsDDKp1_2017", "fc_BsDDKp1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_120/fit_results.txt")
            fc_BuDDK01_2017 = ROOT.TFileCollection("fc_BuDDK01_2017", "fc_BuDDK01_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_130/fit_results.txt")
            fc_BuDD1_2017 = ROOT.TFileCollection("fc_BuDD1_2017", "fc_BuDD1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_150/fit_results.txt")
            fc_BuDDKp1_2018 = ROOT.TFileCollection("fc_BuDDKp1_2018", "fc_BuDDKp1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_100/fit_results.txt")
            fc_BdDDKp1_2018 = ROOT.TFileCollection("fc_BdDDKp1_2018", "fc_BdDDKp1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_110/fit_results.txt")
            fc_BsDDKp1_2018 = ROOT.TFileCollection("fc_BsDDKp1_2018", "fc_BsDDKp1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_120/fit_results.txt")
            fc_BuDDK01_2018 = ROOT.TFileCollection("fc_BuDDK01_2018", "fc_BuDDK01_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_130/fit_results.txt")
            fc_BuDD1_2018 = ROOT.TFileCollection("fc_BuDD1_2018", "fc_BuDD1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_150/fit_results.txt")


            # invariant mass branches
            fc_BuDDKp2_2016 = ROOT.TFileCollection("fc_BuDDKp2_2016", "fc_BuDDKp2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_100/invariant_mass_tree.txt")
            fc_BdDDKp2_2016 = ROOT.TFileCollection("fc_BdDDKp2_2016", "fc_BdDDKp2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_110/invariant_mass_tree.txt")
            fc_BsDDKp2_2016 = ROOT.TFileCollection("fc_BsDDKp2_2016", "fc_BsDDKp2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_120/invariant_mass_tree.txt")
            fc_BuDDK02_2016 = ROOT.TFileCollection("fc_BuDDK02_2016", "fc_BuDDK02_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_130/invariant_mass_tree.txt")
            fc_BuDD2_2016 = ROOT.TFileCollection("fc_BuDD2_2016", "fc_BuDD2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_150/invariant_mass_tree.txt")
            fc_BuDDKp2_2017 = ROOT.TFileCollection("fc_BuDDKp2_2017", "fc_BuDDKp2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_100/invariant_mass_tree.txt")
            fc_BdDDKp2_2017 = ROOT.TFileCollection("fc_BdDDKp2_2017", "fc_BdDDKp2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_110/invariant_mass_tree.txt")
            fc_BsDDKp2_2017 = ROOT.TFileCollection("fc_BsDDKp2_2017", "fc_BsDDKp2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_120/invariant_mass_tree.txt")
            fc_BuDDK02_2017 = ROOT.TFileCollection("fc_BuDDK02_2017", "fc_BuDDK02_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_130/invariant_mass_tree.txt")
            fc_BuDD2_2017 = ROOT.TFileCollection("fc_BuDD2_2017", "fc_BuDD2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_150/invariant_mass_tree.txt")            
            fc_BuDDKp2_2018 = ROOT.TFileCollection("fc_BuDDKp2_2018", "fc_BuDDKp2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_100/invariant_mass_tree.txt")
            fc_BdDDKp2_2018 = ROOT.TFileCollection("fc_BdDDKp2_2018", "fc_BdDDKp2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_110/invariant_mass_tree.txt")
            fc_BsDDKp2_2018 = ROOT.TFileCollection("fc_BsDDKp2_2018", "fc_BsDDKp2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_120/invariant_mass_tree.txt")
            fc_BuDDK02_2018 = ROOT.TFileCollection("fc_BuDDK02_2018", "fc_BuDDK02_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_130/invariant_mass_tree.txt")
            fc_BuDD2_2018 = ROOT.TFileCollection("fc_BuDD2_2018", "fc_BuDD2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_150/invariant_mass_tree.txt")

            t_BuDDKp_2016 = ROOT.TChain("DecayTree")
            t_BdDDKp_2016 = ROOT.TChain("DecayTree")
            t_BsDDKp_2016 = ROOT.TChain("DecayTree")
            t_BuDDK0_2016 = ROOT.TChain("DecayTree")
            t_BuDD_2016 = ROOT.TChain("DecayTree")
            t_BuDDKp_2017 = ROOT.TChain("DecayTree")
            t_BdDDKp_2017 = ROOT.TChain("DecayTree")
            t_BsDDKp_2017 = ROOT.TChain("DecayTree")
            t_BuDDK0_2017 = ROOT.TChain("DecayTree")
            t_BuDD_2017 = ROOT.TChain("DecayTree")
            t_BuDDKp_2018 = ROOT.TChain("DecayTree")
            t_BdDDKp_2018 = ROOT.TChain("DecayTree")
            t_BsDDKp_2018 = ROOT.TChain("DecayTree")
            t_BuDDK0_2018 = ROOT.TChain("DecayTree")
            t_BuDD_2018 = ROOT.TChain("DecayTree")

            t_BuDDKp1_2016 = ROOT.TChain("DecayTree")
            t_BdDDKp1_2016 = ROOT.TChain("DecayTree")
            t_BsDDKp1_2016 = ROOT.TChain("DecayTree")
            t_BuDDK01_2016 = ROOT.TChain("DecayTree")
            t_BuDD1_2016 = ROOT.TChain("DecayTree")
            t_BuDDKp1_2017 = ROOT.TChain("DecayTree")
            t_BdDDKp1_2017 = ROOT.TChain("DecayTree")
            t_BsDDKp1_2017 = ROOT.TChain("DecayTree")
            t_BuDDK01_2017 = ROOT.TChain("DecayTree")
            t_BuDD1_2017 = ROOT.TChain("DecayTree")
            t_BuDDKp1_2018 = ROOT.TChain("DecayTree")
            t_BdDDKp1_2018 = ROOT.TChain("DecayTree")
            t_BsDDKp1_2018 = ROOT.TChain("DecayTree")
            t_BuDDK01_2018 = ROOT.TChain("DecayTree")
            t_BuDD1_2018 = ROOT.TChain("DecayTree")

            t_BuDDKp2_2016 = ROOT.TChain("DecayTree")
            t_BdDDKp2_2016 = ROOT.TChain("DecayTree")
            t_BsDDKp2_2016 = ROOT.TChain("DecayTree")
            t_BuDDK02_2016 = ROOT.TChain("DecayTree")
            t_BuDD2_2016 = ROOT.TChain("DecayTree")
            t_BuDDKp2_2017 = ROOT.TChain("DecayTree")
            t_BdDDKp2_2017 = ROOT.TChain("DecayTree")
            t_BsDDKp2_2017 = ROOT.TChain("DecayTree")
            t_BuDDK02_2017 = ROOT.TChain("DecayTree")
            t_BuDD2_2017 = ROOT.TChain("DecayTree")
            t_BuDDKp2_2018 = ROOT.TChain("DecayTree")
            t_BdDDKp2_2018 = ROOT.TChain("DecayTree")
            t_BsDDKp2_2018 = ROOT.TChain("DecayTree")
            t_BuDDK02_2018 = ROOT.TChain("DecayTree")
            t_BuDD2_2018 = ROOT.TChain("DecayTree")

            t_BuDDKp_2016.AddFileInfoList(fc_BuDDKp_2016.GetList())
            t_BdDDKp_2016.AddFileInfoList(fc_BdDDKp_2016.GetList())
            t_BsDDKp_2016.AddFileInfoList(fc_BsDDKp_2016.GetList())
            t_BuDDK0_2016.AddFileInfoList(fc_BuDDK0_2016.GetList())
            t_BuDD_2016.AddFileInfoList(fc_BuDD_2016.GetList())
            t_BuDDKp_2017.AddFileInfoList(fc_BuDDKp_2017.GetList())
            t_BdDDKp_2017.AddFileInfoList(fc_BdDDKp_2017.GetList())
            t_BsDDKp_2017.AddFileInfoList(fc_BsDDKp_2017.GetList())
            t_BuDDK0_2017.AddFileInfoList(fc_BuDDK0_2017.GetList())
            t_BuDD_2017.AddFileInfoList(fc_BuDD_2017.GetList())
            t_BuDDKp_2018.AddFileInfoList(fc_BuDDKp_2018.GetList())
            t_BdDDKp_2018.AddFileInfoList(fc_BdDDKp_2018.GetList())
            t_BsDDKp_2018.AddFileInfoList(fc_BsDDKp_2018.GetList())
            t_BuDDK0_2018.AddFileInfoList(fc_BuDDK0_2018.GetList())
            t_BuDD_2018.AddFileInfoList(fc_BuDD_2018.GetList())

            t_BuDDKp1_2016.AddFileInfoList(fc_BuDDKp1_2016.GetList())
            t_BdDDKp1_2016.AddFileInfoList(fc_BdDDKp1_2016.GetList())
            t_BsDDKp1_2016.AddFileInfoList(fc_BsDDKp1_2016.GetList())
            t_BuDDK01_2016.AddFileInfoList(fc_BuDDK01_2016.GetList())
            t_BuDD1_2016.AddFileInfoList(fc_BuDD1_2016.GetList())
            t_BuDDKp1_2017.AddFileInfoList(fc_BuDDKp1_2017.GetList())
            t_BdDDKp1_2017.AddFileInfoList(fc_BdDDKp1_2017.GetList())
            t_BsDDKp1_2017.AddFileInfoList(fc_BsDDKp1_2017.GetList())
            t_BuDDK01_2017.AddFileInfoList(fc_BuDDK01_2017.GetList())
            t_BuDD1_2017.AddFileInfoList(fc_BuDD1_2017.GetList())
            t_BuDDKp1_2018.AddFileInfoList(fc_BuDDKp1_2018.GetList())
            t_BdDDKp1_2018.AddFileInfoList(fc_BdDDKp1_2018.GetList())
            t_BsDDKp1_2018.AddFileInfoList(fc_BsDDKp1_2018.GetList())
            t_BuDDK01_2018.AddFileInfoList(fc_BuDDK01_2018.GetList())
            t_BuDD1_2018.AddFileInfoList(fc_BuDD1_2018.GetList())

            t_BuDDKp2_2016.AddFileInfoList(fc_BuDDKp2_2016.GetList())
            t_BdDDKp2_2016.AddFileInfoList(fc_BdDDKp2_2016.GetList())
            t_BsDDKp2_2016.AddFileInfoList(fc_BsDDKp2_2016.GetList())
            t_BuDDK02_2016.AddFileInfoList(fc_BuDDK02_2016.GetList())
            t_BuDD2_2016.AddFileInfoList(fc_BuDD2_2016.GetList())
            t_BuDDKp2_2017.AddFileInfoList(fc_BuDDKp2_2017.GetList())
            t_BdDDKp2_2017.AddFileInfoList(fc_BdDDKp2_2017.GetList())
            t_BsDDKp2_2017.AddFileInfoList(fc_BsDDKp2_2017.GetList())
            t_BuDDK02_2017.AddFileInfoList(fc_BuDDK02_2017.GetList())
            t_BuDD2_2017.AddFileInfoList(fc_BuDD2_2017.GetList())
            t_BuDDKp2_2018.AddFileInfoList(fc_BuDDKp2_2018.GetList())
            t_BdDDKp2_2018.AddFileInfoList(fc_BdDDKp2_2018.GetList())
            t_BsDDKp2_2018.AddFileInfoList(fc_BsDDKp2_2018.GetList())
            t_BuDDK02_2018.AddFileInfoList(fc_BuDDK02_2018.GetList())
            t_BuDD2_2018.AddFileInfoList(fc_BuDD2_2018.GetList())

            t_BuDDKp_2016.GetEntries()
            t_BuDDKp_2017.GetEntries()
            t_BuDDKp_2018.GetEntries()
            t_BuDDKp_2016.Add(t_BuDDKp_2017)
            t_BuDDKp_2016.Add(t_BuDDKp_2018)

            t_BuDDKp1_2016.GetEntries()
            t_BuDDKp1_2017.GetEntries()
            t_BuDDKp1_2018.GetEntries()
            t_BuDDKp1_2016.Add(t_BuDDKp1_2017)
            t_BuDDKp1_2016.Add(t_BuDDKp1_2018)

            t_BuDDKp2_2016.GetEntries()
            t_BuDDKp2_2017.GetEntries()
            t_BuDDKp2_2018.GetEntries()
            t_BuDDKp2_2016.Add(t_BuDDKp2_2017)
            t_BuDDKp2_2016.Add(t_BuDDKp2_2018)

            t_BdDDKp_2016.GetEntries()
            t_BdDDKp_2017.GetEntries()
            t_BdDDKp_2018.GetEntries()
            t_BdDDKp_2016.Add(t_BdDDKp_2017)
            t_BdDDKp_2016.Add(t_BdDDKp_2018)

            t_BdDDKp1_2016.GetEntries()
            t_BdDDKp1_2017.GetEntries()
            t_BdDDKp1_2018.GetEntries()
            t_BdDDKp1_2016.Add(t_BdDDKp1_2017)
            t_BdDDKp1_2016.Add(t_BdDDKp1_2018)

            t_BdDDKp2_2016.GetEntries()
            t_BdDDKp2_2017.GetEntries()
            t_BdDDKp2_2018.GetEntries()
            t_BdDDKp2_2016.Add(t_BdDDKp2_2017)
            t_BdDDKp2_2016.Add(t_BdDDKp2_2018)

            t_BsDDKp_2016.GetEntries()
            t_BsDDKp_2017.GetEntries()
            t_BsDDKp_2018.GetEntries()
            t_BsDDKp_2016.Add(t_BsDDKp_2017)
            t_BsDDKp_2016.Add(t_BsDDKp_2018)

            t_BsDDKp1_2016.GetEntries()
            t_BsDDKp1_2017.GetEntries()
            t_BsDDKp1_2018.GetEntries()
            t_BsDDKp1_2016.Add(t_BsDDKp1_2017)
            t_BsDDKp1_2016.Add(t_BsDDKp1_2018)

            t_BsDDKp2_2016.GetEntries()
            t_BsDDKp2_2017.GetEntries()
            t_BsDDKp2_2018.GetEntries()
            t_BsDDKp2_2016.Add(t_BsDDKp2_2017)
            t_BsDDKp2_2016.Add(t_BsDDKp2_2018)

            t_BuDDK0_2016.GetEntries()
            t_BuDDK0_2017.GetEntries()
            t_BuDDK0_2018.GetEntries()
            t_BuDDK0_2016.Add(t_BuDDK0_2017)
            t_BuDDK0_2016.Add(t_BuDDK0_2018)

            t_BuDDK01_2016.GetEntries()
            t_BuDDK01_2017.GetEntries()
            t_BuDDK01_2018.GetEntries()
            t_BuDDK01_2016.Add(t_BuDDK01_2017)
            t_BuDDK01_2016.Add(t_BuDDK01_2018)

            t_BuDDK02_2016.GetEntries()
            t_BuDDK02_2017.GetEntries()
            t_BuDDK02_2018.GetEntries()
            t_BuDDK02_2016.Add(t_BuDDK02_2017)
            t_BuDDK02_2016.Add(t_BuDDK02_2018)

            t_BuDD_2016.GetEntries()
            t_BuDD_2017.GetEntries()
            t_BuDD_2018.GetEntries()
            t_BuDD_2016.Add(t_BuDD_2017)
            t_BuDD_2016.Add(t_BuDD_2018)

            t_BuDD1_2016.GetEntries()
            t_BuDD1_2017.GetEntries()
            t_BuDD1_2018.GetEntries()
            t_BuDD1_2016.Add(t_BuDD1_2017)
            t_BuDD1_2016.Add(t_BuDD1_2018)

            t_BuDD2_2016.GetEntries()
            t_BuDD2_2017.GetEntries()
            t_BuDD2_2018.GetEntries()
            t_BuDD2_2016.Add(t_BuDD2_2017)
            t_BuDD2_2016.Add(t_BuDD2_2018)

            t_BuDDKp_2016.Add(t_BdDDKp_2016)
            t_BuDDKp_2016.Add(t_BsDDKp_2016)
            t_BuDDKp_2016.Add(t_BuDDK0_2016)
            t_BuDDKp_2016.Add(t_BuDD_2016)

            t_BuDDKp1_2016.Add(t_BdDDKp1_2016)
            t_BuDDKp1_2016.Add(t_BsDDKp1_2016)
            t_BuDDKp1_2016.Add(t_BuDDK01_2016)
            t_BuDDKp1_2016.Add(t_BuDD1_2016)

            t_BuDDKp2_2016.Add(t_BdDDKp2_2016)
            t_BuDDKp2_2016.Add(t_BsDDKp2_2016)
            t_BuDDKp2_2016.Add(t_BuDDK02_2016)
            t_BuDDKp2_2016.Add(t_BuDD2_2016)

            print("BACKGROUND : cocktail MCs")
            print("pre sel entries = ", t_BuDDKp_2016.GetEntries())
            print("gsl entries = ", t_BuDDKp1_2016.GetEntries())
            print("inv mass entries = ", t_BuDDKp2_2016.GetEntries())

            t_BuDDKp_2016.AddFriend(t_BuDDKp1_2016)
            t_BuDDKp_2016.AddFriend(t_BuDDKp2_2016)
        
        
        elif((step_name == "combinatorial") or (step_name == "big_combinatorial")):
            # Background proxy: WS data (K- tau+ tau+) (fit region) (100 files of 2016, 2017 and 2018)
            fc_ws_data_2016 = ROOT.TFileCollection("fc_ws_data_2016", "fc_ws_data_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 10) 
            fc_ws_data_2017 = ROOT.TFileCollection("fc_ws_data_2017", "fc_ws_data_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 10) 
            fc_ws_data_2018 = ROOT.TFileCollection("fc_ws_data_2018", "fc_ws_data_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 10) 
            
            fc_ws_data1_2016 = ROOT.TFileCollection("fc_ws_data1_2016", "fc_ws_data1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 10) 
            fc_ws_data1_2017 = ROOT.TFileCollection("fc_ws_data1_2017", "fc_ws_data1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 10) 
            fc_ws_data1_2018 = ROOT.TFileCollection("fc_ws_data1_2018", "fc_ws_data1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 10) 

            fc_ws_data2_2016 = ROOT.TFileCollection("fc_ws_data2_2016", "fc_ws_data2_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_3/invariant_mass_tree.txt", 10) 
            fc_ws_data2_2017 = ROOT.TFileCollection("fc_ws_data2_2017", "fc_ws_data2_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_3/invariant_mass_tree.txt", 10) 
            fc_ws_data2_2018 = ROOT.TFileCollection("fc_ws_data2_2018", "fc_ws_data2_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_3/invariant_mass_tree.txt", 10) 

            t_ws_data_2016 = ROOT.TChain("DecayTree")
            t_ws_data_2017 = ROOT.TChain("DecayTree")
            t_ws_data_2018 = ROOT.TChain("DecayTree")

            t_ws_data1_2016 = ROOT.TChain("DecayTree")
            t_ws_data1_2017 = ROOT.TChain("DecayTree")
            t_ws_data1_2018 = ROOT.TChain("DecayTree")

            t_ws_data2_2016 = ROOT.TChain("DecayTree")
            t_ws_data2_2017 = ROOT.TChain("DecayTree")
            t_ws_data2_2018 = ROOT.TChain("DecayTree")

            t_ws_data_2016.AddFileInfoList(fc_ws_data_2016.GetList())
            t_ws_data_2017.AddFileInfoList(fc_ws_data_2017.GetList())
            t_ws_data_2018.AddFileInfoList(fc_ws_data_2018.GetList())

            t_ws_data1_2016.AddFileInfoList(fc_ws_data1_2016.GetList())
            t_ws_data1_2017.AddFileInfoList(fc_ws_data1_2017.GetList())
            t_ws_data1_2018.AddFileInfoList(fc_ws_data1_2018.GetList())

            t_ws_data2_2016.AddFileInfoList(fc_ws_data2_2016.GetList())
            t_ws_data2_2017.AddFileInfoList(fc_ws_data2_2017.GetList())
            t_ws_data2_2018.AddFileInfoList(fc_ws_data2_2018.GetList())

            t_ws_data_2016.GetEntries()
            t_ws_data_2017.GetEntries()
            t_ws_data_2018.GetEntries()
            t_ws_data_2016.Add(t_ws_data_2017)
            t_ws_data_2016.Add(t_ws_data_2018)

            t_ws_data1_2016.GetEntries()
            t_ws_data1_2017.GetEntries()
            t_ws_data1_2018.GetEntries()
            t_ws_data1_2016.Add(t_ws_data1_2017)
            t_ws_data1_2016.Add(t_ws_data1_2018)

            t_ws_data2_2016.GetEntries()
            t_ws_data2_2017.GetEntries()
            t_ws_data2_2018.GetEntries()
            t_ws_data2_2016.Add(t_ws_data2_2017)
            t_ws_data2_2016.Add(t_ws_data2_2018)

            t_ws_data_2016.AddFriend(t_ws_data1_2016)
            t_ws_data_2016.AddFriend(t_ws_data2_2016)

            print("BACKGROUND : WS data")
            print("pre sel entries = ", t_ws_data_2016.GetEntries())
            print("gsl entries = ", t_ws_data1_2016.GetEntries())
            print("inv mass entries = ", t_ws_data1_2016.GetEntries())
        
        else:
            print("Wrong step name")
            quit()

        # Convert TChain in RDataframe and apply cuts
        sig_df = ROOT.RDataFrame(t_sig_2016)
        sig_df = sig_df.Filter('(df_status==0) && (df_Bp_M > 4000) && (df_Bp_M < 8000) && ( TMath::Min(TMath::Log10(1 - TMath::Abs(taup_DIRA_ORIVX)) * TMath::Sign(1.0, taup_DIRA_ORIVX), TMath::Log10(1 - TMath::Abs(taum_DIRA_ORIVX)) * TMath::Sign(1.0, taum_DIRA_ORIVX)) < -1 )')

        if((step_name == "physics")  or (step_name == "small_physics")):
            bkg_df = ROOT.RDataFrame(t_BuDDKp_2016) # B+ -> DD K+, B0 -> DD K+, B0s -> DD K+, B+ -> DD K0 and B+ -> DD
            bkg_df = bkg_df.Filter('(df_status==0) && (df_Bp_M > 4000) && (df_Bp_M < 8000) && (taup_BKGCAT <= 60) && (taum_BKGCAT <= 60) && (Bp_BKGCAT <= 60) && ( TMath::Min(TMath::Log10(1 - TMath::Abs(taup_DIRA_ORIVX)) * TMath::Sign(1.0, taup_DIRA_ORIVX), TMath::Log10(1 - TMath::Abs(taum_DIRA_ORIVX)) * TMath::Sign(1.0, taum_DIRA_ORIVX)) < -1 )')
        else:
            bkg_df = ROOT.RDataFrame(t_ws_data_2016) # WS data
            bkg_df = bkg_df.Filter('(df_status==0) && (df_Bp_M > 4000) && (df_Bp_M < 8000) && ( TMath::Min(TMath::Log10(1 - TMath::Abs(taup_DIRA_ORIVX)) * TMath::Sign(1.0, taup_DIRA_ORIVX), TMath::Log10(1 - TMath::Abs(taum_DIRA_ORIVX)) * TMath::Sign(1.0, taum_DIRA_ORIVX)) < -1 )')

    elif(species_name == "DDs"):
        # Prepare datasets
        # Signal proxy: MC without rectangular cuts (fit region) (2016+2017+2018)
        fc_sig_2016 = ROOT.TFileCollection("fc_sig_2016", "fc_sig_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_71/pre_sel_tree.txt")
        fc_sig_2017 = ROOT.TFileCollection("fc_sig_2017", "fc_sig_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_71/pre_sel_tree.txt")
        fc_sig_2018 = ROOT.TFileCollection("fc_sig_2018", "fc_sig_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_71/pre_sel_tree.txt")

        t_sig_2016 = ROOT.TChain("DecayTree")
        t_sig_2017 = ROOT.TChain("DecayTree")
        t_sig_2018 = ROOT.TChain("DecayTree")

        t_sig_2016.AddFileInfoList(fc_sig_2016.GetList())
        t_sig_2017.AddFileInfoList(fc_sig_2017.GetList())
        t_sig_2018.AddFileInfoList(fc_sig_2018.GetList())

        t_sig_2016.GetEntries()
        t_sig_2017.GetEntries()
        t_sig_2018.GetEntries()

        t_sig_2016.Add(t_sig_2017)
        t_sig_2016.Add(t_sig_2018)

        print("SIGNAL")
        print("pre sel entries = ", t_sig_2016.GetEntries())

        # Background proxy: data right sideband (fit region) (1000 files of 2016, 2017 and 2018)
        fc_bkg_2016 = ROOT.TFileCollection("fc_bkg_2016", "fc_bkg_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_81/pre_sel_tree.txt", 500) 
        fc_bkg_2017 = ROOT.TFileCollection("fc_bkg_2017", "fc_bkg_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_81/pre_sel_tree.txt", 500) 
        fc_bkg_2018 = ROOT.TFileCollection("fc_bkg_2018", "fc_bkg_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_81/pre_sel_tree.txt", 500) 

        t_bkg_2016 = ROOT.TChain("DecayTree")
        t_bkg_2017 = ROOT.TChain("DecayTree")
        t_bkg_2018 = ROOT.TChain("DecayTree")

        t_bkg_2016.AddFileInfoList(fc_bkg_2016.GetList())
        t_bkg_2017.AddFileInfoList(fc_bkg_2017.GetList())
        t_bkg_2018.AddFileInfoList(fc_bkg_2018.GetList())

        t_bkg_2016.GetEntries()
        t_bkg_2017.GetEntries()
        t_bkg_2018.GetEntries()
    
        t_bkg_2016.Add(t_bkg_2017)
        t_bkg_2016.Add(t_bkg_2018)
    
        print("BACKGROUND")
        print("pre sel entries = ", t_bkg_2016.GetEntries())

        # Convert TChain in RDataframe
        sig_df = ROOT.RDataFrame(t_sig_2016)
        bkg_df = ROOT.RDataFrame(t_bkg_2016)

        # Apply cuts to signal and background
        # mass (pass DTF is alread required in the create_pre_selection tree)
        sig_df = sig_df.Filter('(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)') # in the fit region
        bkg_df = bkg_df.Filter('(Bp_dtf_M[0] > 5320)') # in the right sideband
    
    else:
        print("Wrong species name. Try either Ktautau or DDs.")
        quit()
    ###################################################################### Retrieving files #####################################################################################

    # print((sig_df.Count()).GetValue())
    # print((bkg_df.Count()).GetValue())

    if((setup_name == "input_features") or (setup_name == "output_performance")):
        print("Species name: ", species_name)
        print("Step name: ", step_name)
        print("Setup name: ", setup_name)
    else:
        print("Wrong setup name. Try input_features or output_performance")
        quit()

    make_classification(sig_df, bkg_df, species_name, step_name, setup_name, cross_validation=makeCrossVal)

    # bdt_2d_cut(X_test_1st_step, y_test_1st_step, thresholds_1st_step, thresholds_2nd_step, decisions_first_step, decisions_second_step)

if __name__ == "__main__":
    main(sys.argv)
