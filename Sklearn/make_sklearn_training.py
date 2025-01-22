import random
import ROOT
import sys

import pandas as pd
import numpy as np
import matplotlib as mpl
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
def draw_variables(sig_data, bkg_data, mc_weights, columns, first_step, isKtautau):
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
        else:
            xlim = np.percentile(np.hstack([sig_data[column]]), [0.01, 99.99])
    
        # plt.hist(sig_data[column], color='b', weights=mc_weights, range=xlim, **hist_settings, label="Signal")
        plt.hist(sig_data[column], color='b', range=xlim, **hist_settings, label="Signal")
        plt.hist(bkg_data[column], color='r', range=xlim, **hist_settings, label="Background")
        plt.legend(loc="upper right", fontsize=20)
        plt.title(column, fontsize=25)
        plt.xlabel(column, fontsize=15)
        plt.ylabel("Normalised entries / ({0})".format(nbins), fontsize=15)

        if(isKtautau == "True"):
            if(first_step):
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Input_features_plots_first_step/'+column+'.pdf')
            else:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Input_features_plots_second_step/'+column+'.pdf')
        else:
            if(first_step):
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Input_features_plots_first_step/'+column+'.pdf')
            else:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Input_features_plots_second_step/'+column+'.pdf') 
        plt.clf()

def draw_profile_plots(data, columns, label, first_step, isKtautau):
    nbins = 30

    for id, column in enumerate(columns, 1):
        if column == 'Bp_M':
            continue

        x_values = np.array(data['Bp_M']) # B+ mass
        y_values = np.array(data[column]) # variable

        h, xe, ye, im = plt.hist2d(x_values, y_values, nbins, range=[[4000,8000],[np.min(data[column]),np.max(data[column])]])

        y_mean = np.zeros(nbins)
        y_error = np.zeros(nbins)
        x = np.zeros(nbins)

        xbinw = xe[1]-xe[0]

        for i in range(len(xe)-1):
            y_mean[i] = y_values[ (x_values > xe[i]) & (x_values <= xe[i+1]) ].mean()
            if(y_mean[i] == 0):
                y_error[i] = 0
            else:
                y_error[i] = (y_values[ (x_values > xe[i]) & (x_values <= xe[i+1]) ].std())/np.sqrt( np.abs(y_mean[i]) )
            x[i] = xe[i]+xbinw/2
    
        plt.clf()
        if(label == 'sig'):
            plt.errorbar(x, y_mean, y_error, fmt='_', ecolor='b', color='b')
        else:
            plt.errorbar(x, y_mean, y_error, fmt='_', ecolor='r', color='r')

        plt.title(column+'   (c = {0})'.format( round(np.corrcoef(x_values, y_values)[0,1],4) ), fontsize=15)
        plt.xlabel("m_B (MeV)", fontsize=15)
        plt.ylabel("Average {0} per m_B bin".format(column), fontsize=15)

        if(isKtautau == "True"):
            if(first_step):
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Input_features_plots_first_step/Profile_plots/'+column+'_'+label+'.pdf')
            else:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Input_features_plots_second_step/Profile_plots/'+column+'_'+label+'.pdf')
        else:
            if(first_step):
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Input_features_plots_first_step/Profile_plots/'+column+'_'+label+'.pdf')
            else:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Input_features_plots_second_step/Profile_plots/'+column+'_'+label+'.pdf') 
        plt.clf()

def draw_scatter_plots(data, columns, first_step):
    fig = sb.pairplot(data, hue="y", corner=True)
    fig._legend.set_title("")
    new_labels = ['Signal', 'Background']
    for t, l in zip(fig._legend.texts, new_labels):
        t.set_text(l)
    plt.setp(fig._legend.get_texts(), fontsize='32')
    if(first_step):
        fig.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Input_features_plots_first_step/scatter_plot.pdf') 
    else:
        fig.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Input_features_plots_second_step/scatter_plot.pdf') 
    plt.clf()    

def correlation_matrix(data, name, first_step, isKtautau, **kwds):
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
    if(isKtautau == "True"):
        if(first_step):
            fig1.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Input_features_plots_first_step/correlation_matrix_'+name+'.pdf') 
        else:
            fig1.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Input_features_plots_second_step/correlation_matrix_'+name+'.pdf') 
    else:
        if(first_step):
            fig1.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Input_features_plots_first_step/correlation_matrix_'+name+'.pdf') 
        else:
            fig1.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Input_features_plots_second_step/correlation_matrix_'+name+'.pdf') 
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

        print(class_weight)

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

def roc_curve_plot(names, classifiers, X_test, y_test, first_step, isKtautau):
    clf_fpr = []
    clf_tpr = []
    clf_thresholds = []

    for i in range(len(names)):
        if((names[i] == "RForest") or (names[i] == "XGBoost")):
            decisions = classifiers[i].predict_proba(X_test)[:, 1]
        else:
            decisions = classifiers[i].decision_function(X_test)

        # Compute ROC curve and area under the curve
        fpr, tpr, thresholds = roc_curve(y_test, decisions)
        clf_fpr.append(fpr)
        clf_tpr.append(tpr)
        clf_thresholds.append(thresholds)

        roc_auc = auc(fpr, tpr)

        plt.plot(fpr, tpr, lw=1, label=names[i]+' (area = %0.2f)'%(roc_auc))

    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Background efficiency')
    plt.ylabel('Signal efficiency')
    plt.title('ROC curve')
    plt.legend(loc="lower right")
    plt.grid()
    if(isKtautau == "True"):
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_first_step/roc_curve.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_second_step/roc_curve.pdf')
    else:
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_first_step/roc_curve.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_second_step/roc_curve.pdf')
    plt.clf()
    return clf_fpr, clf_tpr, clf_thresholds

def bdt_2d_cut(X_test, y_test, thresholds_1st_step, thresholds_2nd_step, decisions_first_step, decisions_second_step):
    # Add bdt columns to dataset
    bdt_output_values_1 = pd.DataFrame({'bdt1': decisions_first_step})
    bdt_output_values_2 = pd.DataFrame({'bdt2': decisions_second_step})

    print(decisions_first_step)
    print(decisions_second_step)

    print(thresholds_1st_step)
    print(thresholds_2nd_step)

    metric = np.zeros((len(thresholds_1st_step),len(thresholds_2nd_step)))

    print(X_test)
    print(y_test)

    # for bdt1 in thresholds_1st_step:
    #     for bdt2 in thresholds_2nd_step:
    #         eps_s = X_test[]


    #     bdt_output_values = pd.DataFrame({'bdt': output})
    #     bdt_output_values_signal =  bdt_output_values[0:signal.shape[0]]
    #     bdt_output_values_background = bdt_output_values[signal.shape[0]: signal.shape[0] + background.shape[0]].reset_index()

    #     print("before bdt cut")
    #     print(signal)
    #     print(background)

    #     print(bdt_output_values_signal)
    #     print(bdt_output_values_background)

    #     signal = pd.concat([signal,bdt_output_values_signal], axis=1)
    #     background = pd.concat([background,bdt_output_values_background], axis=1)

def draw_signal_significance_vs_bdt_cut(name, fpr, tpr, thresholds, N_bkg, first_step, isKtautau):
    eps_s = tpr
    B = fpr*N_bkg
    a = 5

    metric = eps_s/(a/2 + np.sqrt(B))

    plt.plot(thresholds, metric*20, color='g', label='Punzi (x20)')
    plt.plot(thresholds, tpr, color='b', label='Signal eff.')
    plt.plot(thresholds, fpr, color='r', label='Background eff.')
    plt.xlabel('BDT cut value')
    plt.grid()
    plt.legend()
    if(isKtautau == "True"):
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_first_step/significance_vs_bdt_cut_'+name+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_second_step/significance_vs_bdt_cut_'+name+'.pdf')
    else:
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_first_step/significance_vs_bdt_cut_'+name+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_second_step/significance_vs_bdt_cut_'+name+'.pdf')
    plt.clf()

    optimal_index = np.argmax( np.ma.masked_invalid(metric) )
    optimal_metric = metric[optimal_index]
    optimal_cut = thresholds[optimal_index]
    print(name)
    print('The optimal cut value is {0} for a max significance of = {1}'.format(optimal_cut,optimal_metric))
    print('For this BDT cut, signal eff. = {0} and background eff. = {1}'.format(tpr[optimal_index], fpr[optimal_index]))

    return optimal_cut, optimal_metric

def compare_train_test(name, clf, X_train, y_train, X_test, y_test, first_step, isKtautau, bins=30):
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):

        if((name == "RForest") or (name == "XGBoost")):
            d1 = clf.predict_proba(X[y>0.5])[:,1].ravel()
            d2 = clf.predict_proba(X[y<0.5])[:,1].ravel()
        else:
            d1 = clf.decision_function(X[y>0.5]).ravel()
            d2 = clf.decision_function(X[y<0.5]).ravel()
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

    plt.xlabel("Classifier output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='best')
    if(isKtautau == "True"):
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_first_step/classifier_output_'+name+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_second_step/classifier_output_'+name+'.pdf')
    else:
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_first_step/classifier_output_'+name+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_second_step/classifier_output_'+name+'.pdf')
    plt.clf()

def confusion_matrix(name, clf, X_test, y_test, first_step):
    fig = ConfusionMatrixDisplay.from_estimator(clf, X_test, y_test)
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Output_performance_first_step/confusion_matrix_'+name+'.pdf')
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Output_performance_second_step/confusion_matrix_'+name+'.pdf')
    plt.clf()

def draw_feature_importance(name, clf, columns, first_step, input_features, isKtautau):
    if(name == 'XGBoost'): 
        plt.figure(figsize=(20,10))
        plt.barh(columns, clf.feature_importances_)
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        plt.title("Feature importance for {0}".format(name), fontsize=20)
        plt.tight_layout()

        # feature_important = clf.get_booster().get_score(importance_type='weight')
        # keys = list(feature_important.keys())
        # values = list(feature_important.values())
        # data = pd.DataFrame(data=values, index=keys, columns=["score"]).sort_values(by = "score", ascending=False)
        # fig = data.nlargest(40, columns="score").plot(kind='barh', figsize = (12,10)).figure ## plot top 40 features

    else:
        importances = clf.feature_importances_
        feature_importances = pd.Series(importances, index=input_features)

        fig, ax = plt.subplots()
        std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
        feature_importances.plot.bar(yerr=std, ax=ax)
        ax.set_title("Feature importances using MDI")
        ax.set_ylabel("Mean decrease in impurity")
        fig.tight_layout()

    if(isKtautau == "True"):
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_first_step/feature_importance_'+name+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/Output_performance_second_step/feature_importance_'+name+'.pdf')
    else:
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_first_step/feature_importance_'+name+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/Output_performance_second_step/feature_importance_'+name+'.pdf')
    plt.clf()

def make_classification(sig_df, bkg_df, isKtautau, output, cut, first_step, draw_input_features, cross_validation):
    signal = pd.DataFrame()
    background = pd.DataFrame()

    # Variables needed to build the input features for the 1st and 2nd steps
    if(isKtautau == "True"):
        branch_names_1 = ['Bp_VTXISONUMVTX_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 
                        'Bp_VTXISONUMVTX_taup', 'Bp_VTXISONUMVTX_taum', 'Bp_VTXISODCHI2ONETRACK_taup', 'Bp_VTXISODCHI2ONETRACK_taum', 
                        'Bp_VTXISODCHI2MASSONETRACK_taup', 'Bp_VTXISODCHI2MASSONETRACK_taum', 'Bp_VTXISODCHI2TWOTRACK_taup', 'Bp_VTXISODCHI2TWOTRACK_taum',
                        'Bp_VTXISODCHI2MASSTWOTRACK_taup', 'Bp_VTXISODCHI2MASSTWOTRACK_taum', 'Bp_CC_05_DELTAPHI_B', 'Bp_NC_05_IT_B', 'Bp_CC_05_PTASYM_B', 'Bp_CC_05_PX_B',
                        'Bp_CC_05_PZASYM_B', 'Bp_CC_05_DELTAETA_B', 'Bp_NC_05_MULT_B', 'Bp_NC_05_VPT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_B', 'Bp_CC_05_MULT_B',
                        'Bp_CC_05_PASYM_taup', 'Bp_CC_05_DELTAETA_taup', 'Bp_NC_05_MULT_taup', 'Bp_NC_05_SPT_taup', 'Bp_CCNC_05_IT_taup', 'Bp_NC_05_PTASYM_taup', 'Bp_CC_05_PTASYM_taup', 'Bp_CC_05_IT_taup',
                        'Bp_CC_05_PASYM_taum', 'Bp_CC_05_DELTAETA_taum', 'Bp_NC_05_MULT_taum', 'Bp_NC_05_SPT_taum', 'Bp_CCNC_05_IT_taum', 'Bp_NC_05_PTASYM_taum', 'Bp_CC_05_PTASYM_taum', 'Bp_CC_05_IT_taum',
                        'Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taup', 'Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup', 'Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup', 'Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taum', 'Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum', 'Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum', 'df_Bp_M'] 

        branch_names_2 = ['taup_M12', 'taum_M12', 'taup_M23', 'taum_M23', 'taup_M13', 'taum_M13', 'taup_DIRA_ORIVX', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                        'df_BVx', 'df_BVy', 'df_BVz', 'df_DV1x', 'df_DV1y', 'df_DV1z', 'df_DV2x', 'df_DV2y', 'df_DV2z', 'df_PVx', 'df_PVy', 'df_PVz',
                        'df_BVx_err', 'df_BVy_err', 'df_BVz_err', 'df_DV1x_err', 'df_DV1y_err', 'df_DV1z_err', 'df_DV2x_err', 'df_DV2y_err', 'df_DV2z_err', 'df_PVx_err', 'df_PVy_err', 'df_PVz_err',
                        'df_Kp_PX', 'df_Kp_PY', 'df_Kp_PZ', 'taup_M', 'taum_M', 'df_chi2', 'Kp_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV', 
                        'taup_AMAXDOCA', 'taup_AMINDOCA', 'taup_DOCACHI2MAX', 'taum_AMAXDOCA', 'taum_AMINDOCA', 'taum_DOCACHI2MAX', 'df_Bp_M']
    else:
        branch_names_1 = ['Bp_VTXISONUMVTX_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 
                        'Bp_VTXISONUMVTX_D0bar', 'Bp_VTXISONUMVTX_Dsp', 'Bp_VTXISODCHI2ONETRACK_D0bar', 'Bp_VTXISODCHI2ONETRACK_Dsp', 
                        'Bp_VTXISODCHI2MASSONETRACK_D0bar', 'Bp_VTXISODCHI2MASSONETRACK_Dsp', 'Bp_VTXISODCHI2TWOTRACK_D0bar', 'Bp_VTXISODCHI2TWOTRACK_Dsp',
                        'Bp_VTXISODCHI2MASSTWOTRACK_D0bar', 'Bp_VTXISODCHI2MASSTWOTRACK_Dsp', 'Bp_CC_05_DELTAPHI_B', 'Bp_NC_05_IT_B', 'Bp_CC_05_PTASYM_B', 'Bp_CC_05_PX_B',
                        'Bp_CC_05_PZASYM_B', 'Bp_CC_05_DELTAETA_B', 'Bp_NC_05_MULT_B', 'Bp_NC_05_VPT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_B', 'Bp_CC_05_MULT_B',
                        'Bp_CC_05_PASYM_D0bar', 'Bp_CC_05_DELTAETA_D0bar', 'Bp_NC_05_MULT_D0bar', 'Bp_NC_05_SPT_D0bar', 'Bp_CCNC_05_IT_D0bar', 'Bp_NC_05_PTASYM_D0bar', 'Bp_CC_05_PTASYM_D0bar', 'Bp_CC_05_IT_D0bar',
                        'Bp_CC_05_PASYM_Dsp', 'Bp_CC_05_DELTAETA_Dsp', 'Bp_NC_05_MULT_Dsp', 'Bp_NC_05_SPT_Dsp', 'Bp_CCNC_05_IT_Dsp', 'Bp_NC_05_PTASYM_Dsp', 'Bp_CC_05_PTASYM_Dsp', 'Bp_CC_05_IT_Dsp',
                        'Bp_Bstautau_ISOBDTFIRSTVALUE_taup', 'Bp_Bstautau_ISOBDTSECONDVALUE_taup', 'Bp_Bstautau_ISOBDTTHIRDVALUE_taup', 'Bp_Bstautau_ISOBDTFIRSTVALUE_taum', 'Bp_Bstautau_ISOBDTSECONDVALUE_taum', 'Bp_Bstautau_ISOBDTTHIRDVALUE_taum', 'Bp_dtf_M'] 

        branch_names_2 = ['Dsp_M12', 'Dsp_M23', 'Dsp_M13', 'D0bar_DIRA_ORIVX', 'Dsp_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                'Bp_ENDVERTEX_X', 'Bp_ENDVERTEX_Y', 'Bp_ENDVERTEX_Z', 'D0bar_ENDVERTEX_X', 'D0bar_ENDVERTEX_Y', 'D0bar_ENDVERTEX_Z', 'Dsp_ENDVERTEX_X', 'Dsp_ENDVERTEX_Y', 'Dsp_ENDVERTEX_Z', 'Bp_OWNPV_X', 'Bp_OWNPV_Y', 'Bp_OWNPV_Z',
                'Bp_ENDVERTEX_XERR', 'Bp_ENDVERTEX_YERR', 'Bp_ENDVERTEX_ZERR', 'D0bar_ENDVERTEX_XERR', 'D0bar_ENDVERTEX_YERR', 'D0bar_ENDVERTEX_ZERR', 'Dsp_ENDVERTEX_XERR', 'Dsp_ENDVERTEX_YERR', 'Dsp_ENDVERTEX_ZERR', 'Bp_OWNPV_XERR', 'Bp_OWNPV_YERR', 'Bp_OWNPV_ZERR',
                'D0bar_M', 'Dsp_M', 'Bp_dtf_chi2', 'D0bar_AMAXDOCA', 'D0bar_AMINDOCA', 'D0bar_DOCACHI2MAX', 'Dsp_AMAXDOCA', 'Dsp_AMINDOCA', 'Dsp_DOCACHI2MAX',
                'Bp_FDCHI2_OWNPV', 'D0bar_FD_ORIVX', 'Dsp_FD_ORIVX', 'Bp_dtf_M']
    branch_names = branch_names_1 + branch_names_2

    sig = sig_df.AsNumpy(branch_names)
    bkg = bkg_df.AsNumpy(branch_names) 

    # Define the input features
    if(first_step):
        if(isKtautau == "True"):
            # SIGNAL
            signal['VTXISODCHI2ONETRACK_B'] = sig['Bp_VTXISODCHI2ONETRACK_B']
            signal['VTXISONUMVTX_tau_max'] = np.maximum( sig['Bp_VTXISONUMVTX_taup'], sig['Bp_VTXISONUMVTX_taum'] )
            signal['VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( sig['Bp_VTXISODCHI2TWOTRACK_taup'], sig['Bp_VTXISODCHI2TWOTRACK_taum'] ) 
            signal['CC_05_IT_tau_min'] = np.minimum( sig['Bp_CC_05_IT_taup'], sig['Bp_CC_05_IT_taum'] )
            signal['tau_iso_second_value_max'] = np.maximum( sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], sig['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] ) 
            signal['tau_iso_third_value_max'] = np.maximum( sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], sig['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] )

            # BACKGROUND
            background['VTXISODCHI2ONETRACK_B'] = bkg['Bp_VTXISODCHI2ONETRACK_B']
            background['VTXISONUMVTX_tau_max'] = np.maximum( bkg['Bp_VTXISONUMVTX_taup'], bkg['Bp_VTXISONUMVTX_taum'] )
            background['VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( bkg['Bp_VTXISODCHI2TWOTRACK_taup'], bkg['Bp_VTXISODCHI2TWOTRACK_taum'] ) 
            background['CC_05_IT_tau_min'] = np.minimum( bkg['Bp_CC_05_IT_taup'], bkg['Bp_CC_05_IT_taum'] )
            background['tau_iso_second_value_max'] = np.maximum( bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], bkg['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] ) 
            background['tau_iso_third_value_max'] = np.maximum( bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], bkg['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] )

        else:
            # SIGNAL
            signal['VTXISODCHI2ONETRACK_B'] = sig['Bp_VTXISODCHI2ONETRACK_B']
            signal['VTXISONUMVTX_D_max'] = np.maximum( sig['Bp_VTXISONUMVTX_D0bar'], sig['Bp_VTXISONUMVTX_Dsp'] )
            signal['VTXISODCHI2TWOTRACK_D_max'] = np.maximum( sig['Bp_VTXISODCHI2TWOTRACK_D0bar'], sig['Bp_VTXISODCHI2TWOTRACK_Dsp'] ) 
            signal['CC_05_IT_D_min'] = np.minimum( sig['Bp_CC_05_IT_D0bar'], sig['Bp_CC_05_IT_Dsp'] )
            signal['D_iso_second_value_max'] = np.maximum( sig['Bp_Bstautau_ISOBDTSECONDVALUE_taup'], sig['Bp_Bstautau_ISOBDTSECONDVALUE_taum'] ) 
            signal['D_iso_third_value_max'] = np.maximum( sig['Bp_Bstautau_ISOBDTTHIRDVALUE_taup'], sig['Bp_Bstautau_ISOBDTTHIRDVALUE_taum'] )

            # BACKGROUND
            background['VTXISODCHI2ONETRACK_B'] = bkg['Bp_VTXISODCHI2ONETRACK_B']
            background['VTXISONUMVTX_D_max'] = np.maximum( bkg['Bp_VTXISONUMVTX_D0bar'], bkg['Bp_VTXISONUMVTX_Dsp'] )
            background['VTXISODCHI2TWOTRACK_D_max'] = np.maximum( bkg['Bp_VTXISODCHI2TWOTRACK_D0bar'], bkg['Bp_VTXISODCHI2TWOTRACK_Dsp'] ) 
            background['CC_05_IT_D_min'] = np.minimum( bkg['Bp_CC_05_IT_D0bar'], bkg['Bp_CC_05_IT_Dsp'] )
            background['D_iso_second_value_max'] = np.maximum( bkg['Bp_Bstautau_ISOBDTSECONDVALUE_taup'], bkg['Bp_Bstautau_ISOBDTSECONDVALUE_taum'] ) 
            background['D_iso_third_value_max'] = np.maximum( bkg['Bp_Bstautau_ISOBDTTHIRDVALUE_taup'], bkg['Bp_Bstautau_ISOBDTTHIRDVALUE_taum'] )

        input_features = signal.columns.tolist()

    else:
        if(isKtautau == "True"):
            # SIGNAL:
            signal['tau_M_max'] = np.maximum( sig['taup_M'], sig['taum_M'] )
            signal['tau_M12_M23_max_max'] = np.maximum( np.maximum( sig['taup_M12'], sig['taup_M23'] ), np.maximum( sig['taum_M12'], sig['taum_M23'] ) )
            signal['tau_M12_M23_min_max'] = np.minimum( np.maximum( sig['taup_M12'], sig['taup_M23'] ), np.maximum( sig['taum_M12'], sig['taum_M23'] ) )
            signal['tau_M13_max'] = np.maximum( sig['taup_M13'], sig['taum_M13'] )
            signal['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX'] ) )
            signal['tau_FD_BV_min'] = np.minimum( np.sqrt( (sig['df_DV1x'] - sig['df_BVx'])**2 + (sig['df_DV1y'] - sig['df_BVy'])**2 + (sig['df_DV1z'] - sig['df_BVz'])**2 ), np.sqrt( (sig['df_DV2x'] - sig['df_BVx'])**2 + (sig['df_DV2y'] - sig['df_BVy'])**2 + (sig['df_DV2z'] - sig['df_BVz'])**2 ) )
            signal['DV1_DV2_distance'] = np.sqrt( (sig['df_DV1x'] - sig['df_DV2x'])**2 + (sig['df_DV1y'] - sig['df_DV2y'])**2 + (sig['df_DV1z'] - sig['df_DV2z'])**2 )
            signal['log10_df_chi2'] = np.log10( sig['df_chi2'] )
            signal['tau_AMAXDOCA_max'] = np.maximum( sig['taup_AMAXDOCA'], sig['taum_AMAXDOCA'] )

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

            # signal['Kp_IPCHI2_OWNPV'] = sig['Kp_IPCHI2_OWNPV']
            # signal['taum_pi2_IPCHI2_OWNPV'] = sig['taum_pi2_IPCHI2_OWNPV']
            # signal['tau_M_min'] = np.minimum( sig['taup_M'], sig['taum_M'] )
            # signal['tau_M12_M23_max_min'] = np.maximum( np.minimum( sig['taup_M12'], sig['taup_M23'] ), np.minimum( sig['taum_M12'], sig['taum_M23'] ) )
            # signal['tau_M12_M23_min_min'] = np.minimum( np.minimum( sig['taup_M12'], sig['taup_M23'] ), np.minimum( sig['taum_M12'], sig['taum_M23'] ) )
            # signal['tau_M13_max'] = np.maximum( sig['taup_M13'], sig['taum_M13'] )
            # signal['tau_M13_min'] = np.minimum( sig['taup_M13'], sig['taum_M13'] )
            # signal['log_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX'] ) )
            # signal['log_1_minus_B_DIRA_PV'] = np.log(1 - np.abs(sig['Bp_DIRA_OWNPV']) )
            # signal['B_FD_PV'] = np.sqrt( (sig['df_BVx'] - sig['df_PVx'])**2 + (sig['df_BVy'] - sig['df_PVy'])**2 + (sig['df_BVz'] - sig['df_PVz'])**2 )
            # signal['B_FD_PV_chi2'] = np.sqrt( (sig['df_BVx'] - sig['df_PVx'])**2 + (sig['df_BVy'] - sig['df_PVy'])**2 + (sig['df_BVz'] - sig['df_PVz'])**2 )/( np.abs( sig['df_BVx'] - sig['df_PVx'] )*np.sqrt( sig['df_BVx_err']**2 + sig['df_PVx_err']**2  ) + np.abs( sig['df_BVy'] - sig['df_PVy'] )*np.sqrt( sig['df_BVy_err']**2 + sig['df_PVy_err']**2 ) + np.abs( sig['df_BVz'] - sig['df_PVz'] )*np.sqrt( sig['df_BVz_err']**2 + sig['df_PVz_err']**2 ) )
            # signal['tau_FD_BV_max'] = np.maximum( np.sqrt( (sig['df_DV1x'] - sig['df_BVx'])**2 + (sig['df_DV1y'] - sig['df_BVy'])**2 + (sig['df_DV1z'] - sig['df_BVz'])**2 ), np.sqrt( (sig['df_DV2x'] - sig['df_BVx'])**2 + (sig['df_DV2y'] - sig['df_BVy'])**2 + (sig['df_DV2z'] - sig['df_BVz'])**2 ) )
            # signal['tau_FD_BV_chi2_max'] = np.maximum( np.sqrt( (sig['df_DV1x'] - sig['df_BVx'])**2 + (sig['df_DV1y'] - sig['df_BVy'])**2 + (sig['df_DV1z'] - sig['df_BVz'])**2 )/(  np.abs( sig['df_DV1x'] - sig['df_BVx'] )*np.sqrt( sig['df_DV1x_err']**2 + sig['df_BVx_err']**2  ) + np.abs( sig['df_DV1y'] - sig['df_BVy'] )*np.sqrt( sig['df_DV1y_err']**2 + sig['df_BVy_err']**2 ) + np.abs( sig['df_DV1z'] - sig['df_BVz'] )*np.sqrt( sig['df_DV1z_err']**2 + sig['df_BVz_err']**2 ) ) , np.sqrt( (sig['df_DV2x'] - sig['df_BVx'])**2 + (sig['df_DV2y'] - sig['df_BVy'])**2 + (sig['df_DV2z'] - sig['df_BVz'])**2 )/( np.abs( sig['df_DV2x'] - sig['df_BVx'] )*np.sqrt( sig['df_DV2x_err']**2 + sig['df_BVx_err']**2 ) + np.abs( sig['df_DV2y'] - sig['df_BVy'] )*np.sqrt( sig['df_DV2y_err']**2 + sig['df_BVy_err']**2 ) + np.abs( sig['df_DV2z'] - sig['df_BVz'] )*np.abs( sig['df_DV2z_err']**2 + sig['df_BVz_err']**2 ) ) )        
            # signal['tau_FD_BV_chi2_min'] = np.minimum( np.sqrt( (sig['df_DV1x'] - sig['df_BVx'])**2 + (sig['df_DV1y'] - sig['df_BVy'])**2 + (sig['df_DV1z'] - sig['df_BVz'])**2 )/(  np.abs( sig['df_DV1x'] - sig['df_BVx'] )*np.sqrt( sig['df_DV1x_err']**2 + sig['df_BVx_err']**2  ) + np.abs( sig['df_DV1y'] - sig['df_BVy'] )*np.sqrt( sig['df_DV1y_err']**2 + sig['df_BVy_err']**2 ) + np.abs( sig['df_DV1z'] - sig['df_BVz'] )*np.sqrt( sig['df_DV1z_err']**2 + sig['df_BVz_err']**2 ) ) , np.sqrt( (sig['df_DV2x'] - sig['df_BVx'])**2 + (sig['df_DV2y'] - sig['df_BVy'])**2 + (sig['df_DV2z'] - sig['df_BVz'])**2 )/( np.abs( sig['df_DV2x'] - sig['df_BVx'] )*np.sqrt( sig['df_DV2x_err']**2 + sig['df_BVx_err']**2 ) + np.abs( sig['df_DV2y'] - sig['df_BVy'] )*np.sqrt( sig['df_DV2y_err']**2 + sig['df_BVy_err']**2 ) + np.abs( sig['df_DV2z'] - sig['df_BVz'] )*np.abs( sig['df_DV2z_err']**2 + sig['df_BVz_err']**2 ) ) )
            # signal['DV1_DV2_distance_chi2'] = np.sqrt( (sig['df_DV1x'] - sig['df_DV2x'])**2 + (sig['df_DV1y'] - sig['df_DV2y'])**2 + (sig['df_DV1z'] - sig['df_DV2z'])**2 ) / (   np.abs( sig['df_DV1x'] - sig['df_DV2x'] )*np.sqrt( sig['df_DV1x_err']**2 + sig['df_DV2x_err']**2 ) + np.abs( sig['df_DV1y'] - sig['df_DV2y'] )*np.sqrt( sig['df_DV1y_err']**2 + sig['df_DV2y_err']**2 ) + np.abs( sig['df_DV1z'] - sig['df_DV2z'] )*np.sqrt( sig['df_DV1z_err']**2 + sig['df_DV2z_err']**2  ) )

            # a_sig = np.sqrt( ( sig['df_PVx'] - sig['df_DV1x'] )**2 + ( sig['df_PVy'] - sig['df_DV1y'] )**2 + ( sig['df_PVz'] - sig['df_DV1z'] )**2 )
            # b_sig = np.sqrt( ( sig['df_PVx'] - sig['df_DV2x'] )**2 + ( sig['df_PVy'] - sig['df_DV2y'] )**2 + ( sig['df_PVz'] - sig['df_DV2z'] )**2 )
            # c_sig = np.sqrt( ( sig['df_DV1x'] - sig['df_DV2x'] )**2 + ( sig['df_DV1y'] - sig['df_DV2y'] )**2 + ( sig['df_DV1z'] - sig['df_DV2z'] )**2 )
            # s_sig = (a_sig+b_sig+c_sig)/2
            # signal['DV1_DV2_PV_area'] = np.sqrt( s_sig + (s_sig-a_sig)*(s_sig-b_sig)*(s_sig-c_sig) )

            # BACKGROUND
            background['tau_M_max'] = np.maximum( bkg['taup_M'], bkg['taum_M'] )
            background['tau_M12_M23_max_max'] = np.maximum( np.maximum( bkg['taup_M12'], bkg['taup_M23'] ), np.maximum( bkg['taum_M12'], bkg['taum_M23'] ) )
            background['tau_M12_M23_min_max'] = np.minimum( np.maximum( bkg['taup_M12'], bkg['taup_M23'] ), np.maximum( bkg['taum_M12'], bkg['taum_M23'] ) )
            background['tau_M13_max'] = np.maximum( bkg['taup_M13'], bkg['taum_M13'] )
            background['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX'] ) )
            background['tau_FD_BV_min'] = np.minimum( np.sqrt( (bkg['df_DV1x'] - bkg['df_BVx'])**2 + (bkg['df_DV1y'] - bkg['df_BVy'])**2 + (bkg['df_DV1z'] - bkg['df_BVz'])**2 ), np.sqrt( (bkg['df_DV2x'] - bkg['df_BVx'])**2 + (bkg['df_DV2y'] - bkg['df_BVy'])**2 + (bkg['df_DV2z'] - bkg['df_BVz'])**2 ) )
            background['DV1_DV2_distance'] = np.sqrt( (bkg['df_DV1x'] - bkg['df_DV2x'])**2 + (bkg['df_DV1y'] - bkg['df_DV2y'])**2 + (bkg['df_DV1z'] - bkg['df_DV2z'])**2 )
            background['log10_df_chi2'] = np.log10( bkg['df_chi2'] )
            background['tau_AMAXDOCA_max'] = np.maximum( bkg['taup_AMAXDOCA'], bkg['taum_AMAXDOCA'] )

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

            # background['tau_M_min'] = np.minimum( bkg['taup_M'], bkg['taum_M'] )
            # background['tau_M12_M23_max_min'] = np.maximum( np.minimum( bkg['taup_M12'], bkg['taup_M23'] ), np.minimum( bkg['taum_M12'], bkg['taum_M23'] ) )
            # background['tau_M12_M23_min_min'] = np.minimum( np.minimum( bkg['taup_M12'], bkg['taup_M23'] ), np.minimum( bkg['taum_M12'], bkg['taum_M23'] ) )
            # background['tau_M13_min'] = np.minimum( bkg['taup_M13'], bkg['taum_M13'] )
            # background['log_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX'] ) )
            # background['log_1_minus_B_DIRA_PV'] = np.log(1 - np.abs(bkg['Bp_DIRA_OWNPV']) )
            # background['B_FD_PV'] = np.sqrt( (bkg['df_BVx'] - bkg['df_PVx'])**2 + (bkg['df_BVy'] - bkg['df_PVy'])**2 + (bkg['df_BVz'] - bkg['df_PVz'])**2 )
            # background['B_FD_PV_chi2'] = np.sqrt( (bkg['df_BVx'] - bkg['df_PVx'])**2 + (bkg['df_BVy'] - bkg['df_PVy'])**2 + (bkg['df_BVz'] - bkg['df_PVz'])**2 )/( np.abs( bkg['df_BVx'] - bkg['df_PVx'] )*np.sqrt( bkg['df_BVx_err']**2 + bkg['df_PVx_err']**2  ) + np.abs( bkg['df_BVy'] - bkg['df_PVy'] )*np.sqrt( bkg['df_BVy_err']**2 + bkg['df_PVy_err']**2 ) + np.abs( bkg['df_BVz'] - bkg['df_PVz'] )*np.sqrt( bkg['df_BVz_err']**2 + bkg['df_PVz_err']**2 ) )
            # background['tau_FD_BV_max'] = np.maximum( np.sqrt( (bkg['df_DV1x'] - bkg['df_BVx'])**2 + (bkg['df_DV1y'] - bkg['df_BVy'])**2 + (bkg['df_DV1z'] - bkg['df_BVz'])**2 ), np.sqrt( (bkg['df_DV2x'] - bkg['df_BVx'])**2 + (bkg['df_DV2y'] - bkg['df_BVy'])**2 + (bkg['df_DV2z'] - bkg['df_BVz'])**2 ) )
            # background['tau_FD_BV_chi2_max'] = np.maximum( np.sqrt( (bkg['df_DV1x'] - bkg['df_BVx'])**2 + (bkg['df_DV1y'] - bkg['df_BVy'])**2 + (bkg['df_DV1z'] - bkg['df_BVz'])**2 )/(  np.abs( bkg['df_DV1x'] - bkg['df_BVx'] )*np.sqrt( bkg['df_DV1x_err']**2 + bkg['df_BVx_err']**2  ) + np.abs( bkg['df_DV1y'] - bkg['df_BVy'] )*np.sqrt( bkg['df_DV1y_err']**2 + bkg['df_BVy_err']**2 ) + np.abs( bkg['df_DV1z'] - bkg['df_BVz'] )*np.sqrt( bkg['df_DV1z_err']**2 + bkg['df_BVz_err']**2 ) ) , np.sqrt( (bkg['df_DV2x'] - bkg['df_BVx'])**2 + (bkg['df_DV2y'] - bkg['df_BVy'])**2 + (bkg['df_DV2z'] - bkg['df_BVz'])**2 )/( np.abs( bkg['df_DV2x'] - bkg['df_BVx'] )*np.sqrt( bkg['df_DV2x_err']**2 + bkg['df_BVx_err']**2 ) + np.abs( bkg['df_DV2y'] - bkg['df_BVy'] )*np.sqrt( bkg['df_DV2y_err']**2 + bkg['df_BVy_err']**2 ) + np.abs( bkg['df_DV2z'] - bkg['df_BVz'] )*np.abs( bkg['df_DV2z_err']**2 + bkg['df_BVz_err']**2 ) ) )        
            # background['tau_FD_BV_chi2_min'] = np.minimum( np.sqrt( (bkg['df_DV1x'] - bkg['df_BVx'])**2 + (bkg['df_DV1y'] - bkg['df_BVy'])**2 + (bkg['df_DV1z'] - bkg['df_BVz'])**2 )/(  np.abs( bkg['df_DV1x'] - bkg['df_BVx'] )*np.sqrt( bkg['df_DV1x_err']**2 + bkg['df_BVx_err']**2  ) + np.abs( bkg['df_DV1y'] - bkg['df_BVy'] )*np.sqrt( bkg['df_DV1y_err']**2 + bkg['df_BVy_err']**2 ) + np.abs( bkg['df_DV1z'] - bkg['df_BVz'] )*np.sqrt( bkg['df_DV1z_err']**2 + bkg['df_BVz_err']**2 ) ) , np.sqrt( (bkg['df_DV2x'] - bkg['df_BVx'])**2 + (bkg['df_DV2y'] - bkg['df_BVy'])**2 + (bkg['df_DV2z'] - bkg['df_BVz'])**2 )/( np.abs( bkg['df_DV2x'] - bkg['df_BVx'] )*np.sqrt( bkg['df_DV2x_err']**2 + bkg['df_BVx_err']**2 ) + np.abs( bkg['df_DV2y'] - bkg['df_BVy'] )*np.sqrt( bkg['df_DV2y_err']**2 + bkg['df_BVy_err']**2 ) + np.abs( bkg['df_DV2z'] - bkg['df_BVz'] )*np.abs( bkg['df_DV2z_err']**2 + bkg['df_BVz_err']**2 ) ) )
            # background['DV1_DV2_distance_chi2'] = np.sqrt( (bkg['df_DV1x'] - bkg['df_DV2x'])**2 + (bkg['df_DV1y'] - bkg['df_DV2y'])**2 + (bkg['df_DV1z'] - bkg['df_DV2z'])**2 ) / (   np.abs( bkg['df_DV1x'] - bkg['df_DV2x'] )*np.sqrt( bkg['df_DV1x_err']**2 + bkg['df_DV2x_err']**2 ) + np.abs( bkg['df_DV1y'] - bkg['df_DV2y'] )*np.sqrt( bkg['df_DV1y_err']**2 + bkg['df_DV2y_err']**2 ) + np.abs( bkg['df_DV1z'] - bkg['df_DV2z'] )*np.sqrt( bkg['df_DV1z_err']**2 + bkg['df_DV2z_err']**2  ) )
            # background['Kp_IPCHI2_OWNPV'] = bkg['Kp_IPCHI2_OWNPV']
            # background['taum_pi2_IPCHI2_OWNPV'] = bkg['taum_pi2_IPCHI2_OWNPV']

            # a_bkg = np.sqrt( ( bkg['df_PVx'] - bkg['df_DV1x'] )**2 + ( bkg['df_PVy'] - bkg['df_DV1y'] )**2 + ( bkg['df_PVz'] - bkg['df_DV1z'] )**2 )
            # b_bkg = np.sqrt( ( bkg['df_PVx'] - bkg['df_DV2x'] )**2 + ( bkg['df_PVy'] - bkg['df_DV2y'] )**2 + ( bkg['df_PVz'] - bkg['df_DV2z'] )**2 )
            # c_bkg = np.sqrt( ( bkg['df_DV1x'] - bkg['df_DV2x'] )**2 + ( bkg['df_DV1y'] - bkg['df_DV2y'] )**2 + ( bkg['df_DV1z'] - bkg['df_DV2z'] )**2 )
            # s_bkg = (a_bkg+b_bkg+c_bkg)/2
            # background['DV1_DV2_PV_area'] = np.sqrt( s_bkg + (s_bkg-a_bkg)*(s_bkg-b_bkg)*(s_bkg-c_bkg) )
        else:
            # SIGNAL:
            signal['log10_1_minus_D_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(sig['D0bar_DIRA_ORIVX'] ))*np.sign( sig['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(sig['Dsp_DIRA_ORIVX'] ))*np.sign( sig['Dsp_DIRA_ORIVX'] ) )
            signal['D_FD_BV_min'] = np.minimum( sig['D0bar_FD_ORIVX'], sig['Dsp_FD_ORIVX'] )
            signal['DV1_DV2_distance'] = np.sqrt( (sig['D0bar_ENDVERTEX_X'] - sig['Dsp_ENDVERTEX_X'])**2 + (sig['D0bar_ENDVERTEX_Y'] - sig['Dsp_ENDVERTEX_Y'])**2 + (sig['D0bar_ENDVERTEX_Z'] - sig['Dsp_ENDVERTEX_Z'])**2 )
            
            dtf_chi2_sig = sig['Bp_dtf_chi2']
            for i in range(len(dtf_chi2_sig)):
                dtf_chi2_sig[i] = np.log10( dtf_chi2_sig[i][0] )
            signal['log10_DTF_chi2'] = dtf_chi2_sig

            signal['D_AMAXDOCA_max'] = np.maximum( sig['D0bar_AMAXDOCA'], sig['Dsp_AMAXDOCA'] )

            # BACKGROUND
            background['log10_1_minus_D_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(bkg['D0bar_DIRA_ORIVX'] ))*np.sign( bkg['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(bkg['Dsp_DIRA_ORIVX'] ))*np.sign( bkg['Dsp_DIRA_ORIVX'] ) )
            background['D_FD_BV_min'] = np.minimum( bkg['D0bar_FD_ORIVX'], bkg['Dsp_FD_ORIVX'] )
            background['DV1_DV2_distance'] = np.sqrt( (bkg['D0bar_ENDVERTEX_X'] - bkg['Dsp_ENDVERTEX_X'])**2 + (bkg['D0bar_ENDVERTEX_Y'] - bkg['Dsp_ENDVERTEX_Y'])**2 + (bkg['D0bar_ENDVERTEX_Z'] - bkg['Dsp_ENDVERTEX_Z'])**2 )

            dtf_chi2_bkg = bkg['Bp_dtf_chi2']
            for i in range(len(dtf_chi2_bkg)):
                dtf_chi2_bkg[i] = np.log10(  dtf_chi2_bkg[i][0] )
            background['log10_DTF_chi2'] = dtf_chi2_bkg

            background['D_AMAXDOCA_max'] = np.maximum( bkg['D0bar_AMAXDOCA'], bkg['Dsp_AMAXDOCA'] )      

        input_features = signal.columns.tolist()

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

    # for sklearn data is usually organised into one 2D array of shape (n_samples x n_features)
    # containing all the data and one array of categories of length n_samples
    X = np.concatenate((signal, background))
    y = np.concatenate((np.ones(signal.shape[0]), np.zeros(background.shape[0])))

    # data = pd.DataFrame(np.hstack((X, y.reshape(y.shape[0], -1))),columns=input_features+['weight']+['y'])
    data = pd.DataFrame(np.hstack((X, y.reshape(y.shape[0], -1))),columns=input_features+['y'])

    # Plotting variables and correlations
    if(draw_input_features):
        print("Drawing variables and correlations")

        draw_variables(signal, background, mc_weights, input_features, first_step, isKtautau)
        # data_scatter = data.drop('weight', axis=1)
        # draw_scatter_plots(data_scatter, input_features, first_step)
        # sig_corr = signal.drop('weight', axis=1)
        # bkg_corr = background.drop('weight', axis=1)
        # correlation_matrix(sig_corr, "sig", first_step)
        # correlation_matrix(bkg_corr, "bkg", first_step)

        sig_corr = signal
        bkg_corr = background
        if(isKtautau == "True"):
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

        # draw_profile_plots(signal, input_features, "sig", first_step, isKtautau)
        # draw_profile_plots(background, input_features, "bkg", first_step, isKtautau)

        correlation_matrix(signal, "sig", first_step, isKtautau)
        correlation_matrix(background, "bkg", first_step, isKtautau)

        quit()

    # Train and test split
    # X_dev,X_eval, y_dev,y_eval = train_test_split(X, y, test_size=0.1, random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    # X_train_weights = X_train[:, -1]
    # X_test_weights = X_test[:, -1]
    # X_train = np.delete(X_train, -1, 1)
    # X_test = np.delete(X_test, -1, 1)
    # X = np.delete(X, -1, 1)

    print("Signal train size = ", X_train[y_train>0.5].shape[0] )
    print("Signal test size = ", X_test[y_test>0.5].shape[0] )
    print("Background train size = ", X_train[y_train<0.5].shape[0] )
    print("Background test size = ", X_test[y_test<0.5].shape[0])

    ##############################################################################################################################################

    ###################################################### Train a classifier ####################################################################
    # I am using the optimised hyper-parameters here:
    # classifiers = best_classifiers
    # names = ["AdaBDT", "GradBDT", "RForest", "XGBoost"]
    names = ["XGBoost"]
    if(first_step):
        classifiers = [ # AdaBoostClassifier(DecisionTreeClassifier(max_depth=3), learning_rate=0.05, n_estimators=800, random_state=42),
                        # GradientBoostingClassifier(max_depth=3, learning_rate=0.1, n_estimators=800, random_state=42), 
                        # RandomForestClassifier(class_weight='balanced', max_depth=9, n_estimators=800, random_state=42),
                        xgb.XGBClassifier(tree_method="hist", early_stopping_rounds=2, n_estimators=100, max_depth=6, learning_rate=0.1, random_state=42, importance_type='weight') ] # hist is the fastest tree method
    else:
        classifiers = [ # AdaBoostClassifier(DecisionTreeClassifier(max_depth=3), learning_rate=0.1, n_estimators=800, random_state=42),
                        # GradientBoostingClassifier(max_depth=3, learning_rate=0.1, n_estimators=500, random_state=42),
                        # RandomForestClassifier(class_weight='balanced', max_depth=9, n_estimators=800, random_state=42),
                        xgb.XGBClassifier(tree_method="hist", early_stopping_rounds=2, n_estimators=100, max_depth=6, learning_rate=0.1, random_state=42, importance_type='weight') ] # hist is the fastest tree method

    clf_fpr = []
    clf_tpr = []
    clf_thresholds = []

    clf_output = []

    for i in range(len(names)):
        name = names[i]
        clf = classifiers[i]

        print("Fitting "+name)

        if( (name =="AdaBDT") or (name =="GradBDT") or(name == "XGBoost") ):
            class_weights = compute_class_weight(class_weight='balanced', classes=np.unique(y_train), y=y_train)
            class_weight = np.zeros(len(X_train[:, -1]))

            # signal: y > 0.5 (class 1)
            # bakcground: y < 0.5 (class 0)
            for i in range(len(y_train)):
                if y_train[i] > 0.5:
                    class_weight[i] = class_weights[1]
                else:
                    class_weight[i] = class_weights[0]

            # clf.fit(X_train, y_train, X_train_weights*class_weight)

            if(name == "XGBoost"):
                clf.fit(X_train, y_train, class_weight, eval_set=[(X_test, y_test)], verbose=False)
            else:
                clf.fit(X_train, y_train, class_weight)
        else:
            # clf.fit(X_train, y_train, X_train_weights)
            clf.fit(X_train, y_train)

        if(cross_validation):
            # do_cross_validation(name, clf, X_train, y_train, X_test, y_test, X_train_weights) 
            do_cross_validation(name, clf, X_train, y_train, X_test, y_test) 
        else:
            y_predicted = clf.predict(X_test)
            print(classification_report(y_test, y_predicted, target_names=["background", "signal"]))

            # print("Drawing confusion matrix")
            # confusion_matrix(name, clf, X_test, y_test, first_step)

            print("Comparing train and test sets")
            compare_train_test(name, clf, X_train, y_train, X_test, y_test, first_step, isKtautau) 

    if(cross_validation):
        quit()

    # Compare different classifiers
    print("Drawing roc curve")
    fpr, tpr, thresholds = roc_curve_plot(names, classifiers, X_test, y_test, first_step, isKtautau)

    cuts = []
    metrics = []

    print("Drawing signal significance vs BDT cut")
    N_bkg = X_test[y_test<0.5].shape[0]
    print("N_bkg = ", N_bkg)

    for i in range(len(fpr)):
        bdt_cut, significance = draw_signal_significance_vs_bdt_cut(names[i], fpr[i], tpr[i], thresholds[i], N_bkg, first_step, isKtautau)
        cuts.append(bdt_cut)
        metrics.append(significance)

    for i in range(len(names)):
        name = names[i]
        print("Drawing feature importances for "+ name)
        draw_feature_importance(names[i], classifiers[i], input_features, first_step, input_features, isKtautau)

    # find which classifier gives the maximum value for the metric
    # idx = np.argmax( metrics )
    # if(names[idx] == "RForest"):
    #     decisions = classifiers[idx].predict_proba(X)[:, 1]
    # else:
    #     decisions = classifiers[idx].decision_function(X)

    # Write trained classifiers to a file
    if(isKtautau == "True"):
        if(first_step):
            with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/clf_first_step.pkl', 'wb') as f:
                for i in range(len(classifiers)):
                    pickle.dump(classifiers[i], f)
        else:
            with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/clf_second_step.pkl', 'wb') as f:
                for i in range(len(classifiers)):
                    pickle.dump(classifiers[i], f)
    else:
        if(first_step):
            with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/clf_first_step.pkl', 'wb') as f:
                for i in range(len(classifiers)):
                    pickle.dump(classifiers[i], f)
        else:
            with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/clf_second_step.pkl', 'wb') as f:
                for i in range(len(classifiers)):
                    pickle.dump(classifiers[i], f)

    # Correlation coefficient between BDT response on signal and background (all dataset) vs B+ mass
    # decisions_sig = []
    # decisions_bkg = []
    # for i in range(len(names)):
    #     if( (names[i]=="RForest") or (names[i] == "XGBoost") ):
    #         decisions_sig.append( classifiers[i].predict_proba(np.array(signal))[:, 1] )
    #         decisions_bkg.append( classifiers[i].predict_proba(np.array(background))[:, 1] )
    #     else:
    #         decisions_sig.append( classifiers[i].decision_function(np.array(signal)) )
    #         decisions_bkg.append( classifiers[i].decision_function(np.array(background)) )

    # mass_values_sig = np.array(sig['df_Bp_M'])
    # mass_values_bkg = np.array(bkg['df_Bp_M'])

    # bdt_values_sig = decisions_sig[0]
    # bdt_values_bkg = decisions_bkg[0]

    # print("Correlation of BDT w/ B+ mass (signal) = {0}".format(round(np.corrcoef(bdt_values_sig, mass_values_sig)[0,1],4)))
    # print("Correlation of BDT w/ B+ mass (background) = {0}".format(round(np.corrcoef(bdt_values_bkg, mass_values_bkg)[0,1],4)))


def main(argv):

    isKtautau = argv[1]
    print(isKtautau)

    if(isKtautau == "True"):
        # Prepare datasets
        # Signal proxy: 3pi3pi MC (signal region) (2016+2017+2018)
        fc_sig_2016 = ROOT.TFileCollection("fc_sig_2016", "fc_sig_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_10/pre_sel_tree.txt")
        fc_sig_2017 = ROOT.TFileCollection("fc_sig_2017", "fc_sig_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_10/pre_sel_tree.txt")
        fc_sig_2018 = ROOT.TFileCollection("fc_sig_2018", "fc_sig_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_10/pre_sel_tree.txt")

        fc_sig1_2016 = ROOT.TFileCollection("fc_sig1_2016", "fc_sig1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_10/fit_results.txt")
        fc_sig1_2017 = ROOT.TFileCollection("fc_sig1_2017", "fc_sig1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_10/fit_results.txt")
        fc_sig1_2018 = ROOT.TFileCollection("fc_sig1_2018", "fc_sig1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_10/fit_results.txt")

        # t_sig = ROOT.TChain()
        # t_sig.Add("/panfs/felician/B2Ktautau/workflow/bdt_reweighter/2016/DDs_correction/Species_10/tree_with_weights.root/DecayTree")
        # t_sig.Add("/panfs/felician/B2Ktautau/workflow/bdt_reweighter/2017/DDs_correction/Species_10/tree_with_weights.root/DecayTree")
        # t_sig.Add("/panfs/felician/B2Ktautau/workflow/bdt_reweighter/2018/DDs_correction/Species_10/tree_with_weights.root/DecayTree")

        t_sig_2016 = ROOT.TChain("DecayTree")
        t_sig_2017 = ROOT.TChain("DecayTree")
        t_sig_2018 = ROOT.TChain("DecayTree")

        t_sig1_2016 = ROOT.TChain("DecayTree")
        t_sig1_2017 = ROOT.TChain("DecayTree")
        t_sig1_2018 = ROOT.TChain("DecayTree")

        t_sig_2016.AddFileInfoList(fc_sig_2016.GetList())
        t_sig_2017.AddFileInfoList(fc_sig_2017.GetList())
        t_sig_2018.AddFileInfoList(fc_sig_2018.GetList())

        t_sig1_2016.AddFileInfoList(fc_sig1_2016.GetList())
        t_sig1_2017.AddFileInfoList(fc_sig1_2017.GetList())
        t_sig1_2018.AddFileInfoList(fc_sig1_2018.GetList())

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

        t_sig_2016.AddFriend(t_sig1_2016)

        # Background proxy: WS data (signal region) 
        fc_bkg_2016 = ROOT.TFileCollection("fc_bkg_2016", "fc_bkg_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 100) 
        fc_bkg_2017 = ROOT.TFileCollection("fc_bkg_2017", "fc_bkg_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 100) 
        fc_bkg_2018 = ROOT.TFileCollection("fc_bkg_2018", "fc_bkg_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 100) 
        
        fc_bkg1_2016 = ROOT.TFileCollection("fc_bkg1_2016", "fc_bkg1_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 100) 
        fc_bkg1_2017 = ROOT.TFileCollection("fc_bkg1_2017", "fc_bkg1_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 100) 
        fc_bkg1_2018 = ROOT.TFileCollection("fc_bkg1_2018", "fc_bkg1_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 100) 

        t_bkg_2016 = ROOT.TChain("DecayTree")
        t_bkg_2017 = ROOT.TChain("DecayTree")
        t_bkg_2018 = ROOT.TChain("DecayTree")

        t_bkg1_2016 = ROOT.TChain("DecayTree")
        t_bkg1_2017 = ROOT.TChain("DecayTree")
        t_bkg1_2018 = ROOT.TChain("DecayTree")

        t_bkg_2016.AddFileInfoList(fc_bkg_2016.GetList())
        t_bkg_2017.AddFileInfoList(fc_bkg_2017.GetList())
        t_bkg_2018.AddFileInfoList(fc_bkg_2018.GetList())

        t_bkg1_2016.AddFileInfoList(fc_bkg1_2016.GetList())
        t_bkg1_2017.AddFileInfoList(fc_bkg1_2017.GetList())
        t_bkg1_2018.AddFileInfoList(fc_bkg1_2018.GetList())

        t_bkg_2016.GetEntries()
        t_bkg_2017.GetEntries()
        t_bkg_2018.GetEntries()
        t_bkg_2016.Add(t_bkg_2017)
        t_bkg_2016.Add(t_bkg_2018)

        t_bkg1_2016.GetEntries()
        t_bkg1_2017.GetEntries()
        t_bkg1_2018.GetEntries()
        t_bkg1_2016.Add(t_bkg1_2017)
        t_bkg1_2016.Add(t_bkg1_2018)

        t_bkg_2016.AddFriend(t_bkg1_2016)

        # Convert TChain in RDataframe
        sig_df = ROOT.RDataFrame(t_sig_2016)
        bkg_df = ROOT.RDataFrame(t_bkg_2016)

        # Apply cuts to signal and background
        # pass fitter
        sig_df = sig_df.Filter('(df_status==0)')
        bkg_df = bkg_df.Filter('(df_status==0)')

        sig_df = sig_df.Filter('(df_Bp_M > 4000) && (df_Bp_M < 8000)')
        bkg_df = bkg_df.Filter('(df_Bp_M > 4000) && (df_Bp_M < 8000)')
    else:
        # Prepare datasets
        # Signal proxy: 3pi3pi MC (signal region) (2016+2017+2018)
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

        # Background proxy: WS data (signal region) for both steps (only 2016 for now)
        fc_bkg_2016 = ROOT.TFileCollection("fc_bkg_2016", "fc_bkg_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_81/pre_sel_tree.txt", 1000) 
        fc_bkg_2017 = ROOT.TFileCollection("fc_bkg_2017", "fc_bkg_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_81/pre_sel_tree.txt", 1000) 
        fc_bkg_2018 = ROOT.TFileCollection("fc_bkg_2018", "fc_bkg_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_81/pre_sel_tree.txt", 1000) 

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

        # Convert TChain in RDataframe
        sig_df = ROOT.RDataFrame(t_sig_2016)
        bkg_df = ROOT.RDataFrame(t_bkg_2016)

        # Apply cuts to signal and background
        # pass fitter
        sig_df = sig_df.Filter('(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)')
        bkg_df = bkg_df.Filter('(Bp_dtf_M[0] > 5320)')

    # print((sig_df.Count()).GetValue())
    # print((bkg_df.Count()).GetValue())

    # print("1st step")
    # make_classification(sig_df, bkg_df, isKtautau, output=None, cut=None, first_step=True, draw_input_features=False, cross_validation=False)

    print("2nd step")
    make_classification(sig_df, bkg_df, isKtautau, output=None, cut=None, first_step=False, draw_input_features=False, cross_validation=False)

    # bdt_2d_cut(X_test_1st_step, y_test_1st_step, thresholds_1st_step, thresholds_2nd_step, decisions_first_step, decisions_second_step)

if __name__ == "__main__":
    main(sys.argv)