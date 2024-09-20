import random
import ROOT

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sb

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF

from sklearn.metrics import classification_report, roc_auc_score, roc_curve, auc
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate

first_step = True
draw_input_features = False
cross_validation = False
year = 6

# Functions
def FWHM(X,Y):
    half_max = max(Y) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    #plot(X[0:len(d)],d) #if you are interested
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0]
    right_idx = np.where(d < 0)[-1]
    return X[right_idx] - X[left_idx], left_idx, right_idx #return the difference (full width)

def draw_variables(sig_data, bkg_data, mc_weights, columns, first_step):
    nbins = 40

    hist_settings = {'bins': nbins, 'density': True, 'alpha': 0.4} # density = True means the histograms are normalised
    plt.figure(figsize=[15, 10])

    for id, column in enumerate(columns, 1):
        if(column == "VTXISODCHI2TWOTRACK_min"):
            xlim = [0,100]
        else:
            xlim = np.percentile(np.hstack([sig_data[column]]), [0.01, 99.99])
    
        plt.hist(sig_data[column], color='b', weights=mc_weights, range=xlim, **hist_settings, label="Signal")
        plt.hist(bkg_data[column], color='r', range=xlim, **hist_settings, label="Background")
        plt.legend(loc="upper right", fontsize=20)
        plt.title(column, fontsize=25)
        plt.xlabel(column, fontsize=15)
        plt.ylabel("Normalised entries / ({0})".format(nbins), fontsize=15)
        if(first_step):
            plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots_first_step/'.format(year)+column+'.pdf')
        else:
            plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots_second_step/'.format(year)+column+'.pdf')
        plt.clf()
    
def draw_scatter_plots(data, columns, first_step):
    fig = sb.pairplot(data, hue="y", corner=True)
    fig._legend.set_title("")
    new_labels = ['Signal', 'Background']
    for t, l in zip(fig._legend.texts, new_labels):
        t.set_text(l)
    plt.setp(fig._legend.get_texts(), fontsize='32')
    if(first_step):
        fig.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots_first_step/scatter_plot.pdf'.format(year)) 
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots_second_step/'.format(year)+column+'.pdf')
    plt.clf()

def correlation_matrix(data, name, first_step, **kwds):
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
    if(first_step):
        fig1.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots/correlation_matrix_'.format(year)+name+'_first_step.pdf'.format(year)) 
    else:
        fig1.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots/correlation_matrix_'.format(year)+name+'_second_step.pdf'.format(year)) 
    fig1.clf()

def do_cross_validation(name, clf, X_dev, y_dev, X_val, y_val, X_dev_weights):
    scoring = {'accuracy' : make_scorer(accuracy_score), 
               'precision' : make_scorer(precision_score),
               'recall' : make_scorer(recall_score), 
               'f1_score' : make_scorer(f1_score)}

    print(name+" initial cross validation scores")
    scores = cross_validate(clf, X_dev, y_dev, n_jobs=6, cv=3, scoring=scoring) 

    accuracy_val = scores['test_accuracy']
    precision_val = scores['test_precision']
    recall_val = scores['test_recall']
    f1_score_val = scores['test_f1_score']

    print("Accuracy: {0} +/- {1}".format(accuracy_val.mean(), accuracy_val.std()))
    print("Precision: {0} +/- {1}".format(precision_val.mean(), precision_val.std()))
    print("Recall: {0} +/- {1}".format(recall_val.mean(), recall_val.std()))
    print("F1 score: {0} +/- {1}".format(f1_score_val.mean(), f1_score_val.std()))

    # Optimise hyper-parameters
    if(name == "DTree"):
        param_grid = {"max_depth": [1,2,3,4,5,6,7,8,9,10]}
    elif(name == "RForest"):
        param_grid = {"max_features": [1,2,3,4,5,6,7,8,9,10]}
    elif(name == "AdaBoost"):
        param_grid = {'learning_rate': [0.01, 0.02, 0.03, 0.04, 0.05,0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.5, 1.8, 2]}

    print("Performing grid search for hyper-paramter optimisation of "+name)
    # RandomizedSearchCV samples from distributions of the hyper-parameters and evaluates random points. It can be quicker than GridSearchCV
    grid_result = GridSearchCV(clf, param_grid, cv=3, scoring='roc_auc', n_jobs=6)
    grid_result.fit(X_dev, y_dev, sample_weight=X_dev_weights)
    # grid_result.fit(X_dev, y_dev)

    print("Best parameter set found on development set:")
    print(grid_result.best_estimator_)

    print("Grid scores on a subset of the development set:")
    cv_results = grid_result.cv_results_
    mean_score = cv_results['mean_test_score'] # mean f1 score from the (3) cross-validation splits
    std_score = cv_results['std_test_score'] # standard deviation of f1 score from the (3) cross validation splits
    params = cv_results['params'] # dictionary of the valur of the parameters for each combination in param_grid

    for i in range(len(mean_score)):
        print( "{0} +/- {1} for {2}".format(mean_score[i], std_score[i], params[i]) )
    
    # y_true, y_pred = y_dev, grid_result.decision_function(X_dev) # Call decision_function on the estimator with the best found parameters.
    # print( "It scores {0} on the full development set".format( roc_auc_score(y_true, y_pred) ) )

    # y_true, y_pred = y_eval, grid_result.decision_function(X_eval)
    # print( "It scores {0} on the full evaluation set".format( roc_auc_score(y_true, y_pred) ) )

    print(name+" best cross validation scores")
    scores = cross_validate(grid_result.best_estimator_, X_dev, y_dev, n_jobs=6, cv=3, scoring=scoring) 

    accuracy_val = scores['test_accuracy']
    precision_val = scores['test_precision']
    recall_val = scores['test_recall']
    f1_score_val = scores['test_f1_score']

    print("Accuracy: {0} +/- {1}".format(accuracy_val.mean(), accuracy_val.std()))
    print("Precision: {0} +/- {1}".format(precision_val.mean(), precision_val.std()))
    print("Recall: {0} +/- {1}".format(recall_val.mean(), recall_val.std()))
    print("F1 score: {0} +/- {1}".format(f1_score_val.mean(), f1_score_val.std()))

    return grid_result.best_estimator_

def roc_curve_plot(names, classifiers, X_test, y_test, first_step):
    clf_fpr = []
    clf_tpr = []
    clf_thresholds = []

    for i in range(len(names)):
        if((names[i] == "DTree") or (names[i] == "RForest")):
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
    plt.xlabel('False Positive Rate = background efficiency')
    plt.ylabel('True Positive Rate = signal efficiency')
    plt.title('ROC curve')
    plt.legend(loc="lower right")
    plt.grid()
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_first_step/roc_curve.pdf'.format(year))
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_second_step/roc_curve.pdf'.format(year))
    plt.clf()
    return clf_fpr, clf_tpr, clf_thresholds

def draw_signal_significance_vs_bdt_cut(name, fpr, tpr, thresholds, N_sig, N_bkg, first_step):
    S = N_sig*tpr
    B = N_bkg*fpr
    metric = S/np.sqrt(S+B)

    plt.plot(thresholds, metric, color='g', label='$\\frac{S}{\\sqrt{S+B}}$')
    plt.plot(thresholds, tpr*70, color='b', label='Signal eff. (x70)')
    plt.plot(thresholds, fpr*70, color='r', label='Background eff. (x70)')
    plt.plot(thresholds, metric, color='g')
    plt.xlabel('BDT cut value')
    plt.grid()
    plt.legend()
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_first_step/significance_vs_bdt_cut_'.format(year)+name+'.pdf')
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_second_step/significance_vs_bdt_cut_'.format(year)+name+'.pdf')
    plt.clf()

    optimal_index = np.argmax( np.ma.masked_invalid(metric) )
    optimal_metric = metric[optimal_index]
    optimal_cut = thresholds[optimal_index]
    print(name)
    print('The optimal cut value is {0} for a max significance of S/sqrt(S+B) = {1}'.format(optimal_cut,optimal_metric))
    print('For this BDT cut, signal eff. = {0} and background eff. = {1}'.format(tpr[optimal_index], fpr[optimal_index]))

    return optimal_cut, optimal_metric

def compare_train_test(name, clf, X_train, y_train, X_test, y_test, first_step, bins=30):
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):

        if((name == "DTree") or (name == "RForest")):
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
             color='r', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='S (train)')
    plt.hist(decisions[1],
             color='b', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='B (train)')

    hist, bins = np.histogram(decisions[2],
                              bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')
    
    hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')

    plt.xlabel("Classifier output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='best')
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_first_step/classifier_output_'.format(year)+name+'.pdf')
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_second_step/classifier_output_'.format(year)+name+'.pdf')
    plt.clf()

def confusion_matrix(name, clf, X_test, y_test, first_step):
    fig = ConfusionMatrixDisplay.from_estimator(clf, X_test, y_test)
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_first_step/confusion_matrix_'.format(year)+name+'.pdf')
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_second_step/confusion_matrix_'.format(year)+name+'.pdf')
    plt.clf()

def draw_feature_importance(name, clf, first_step):
    importances = clf.feature_importances_
    feature_importances = pd.Series(importances, index=input_features)

    fig, ax = plt.subplots()
    std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
    feature_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()
    if(first_step):
        fig.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_first_step/feature_importance_'.format(year)+name+'.pdf')
    else:
        fig.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_second_step/feature_importance_'.format(year)+name+'.pdf')
    fig.clf()

######################################################## Prepare dataset ###################################################################
sig_filenames = open('/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_10/pre_sel_tree.txt'.format(year), 'r')
bkg_filenames = open('/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_3/pre_sel_tree.txt'.format(year), 'r')

if(first_step): # tau isolation
    branch_names = ['Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 'Bp_CCNC_20_IT_B',
                    'Bp_VTXISODCHI2MASSONETRACK_taup', 'Bp_VTXISODCHI2MASSONETRACK_taup', 'Bp_CC_20_MULT_taup', 'Bp_CC_20_PTASYM_taup', 'Bp_NC_20_PZASYM_taup',
                    'Bp_VTXISODCHI2TWOTRACK_taum', 'Bp_VTXISODCHI2MASSONETRACK_taum', 'Bp_CC_20_MULT_taum', 'Bp_CC_20_PTASYM_taum', 'Bp_NC_20_PZASYM_taum',
                    'df_Bp_M',
                    'taup_M12', 'taup_M23', 'taum_M12', 'taum_M23', 'Bp_M13',
                    'taup_DIRA_ORIVX', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                    'df_BVx', 'df_BVy', 'df_BVz',
                    'df_DV1x', 'df_DV1y', 'df_DV1z', 'df_DV2x', 'df_DV2y', 'df_DV2z',
                    'df_Kp_PX', 'df_Kp_PY', 'df_Kp_PZ',
                    'df_PVx', 'df_PVy', 'df_PVz'
                    ]

sig_root_files = ROOT.std.vector('string')()
for file in sig_filenames.readlines(): 
    file = file.replace("\n", "")
    sig_root_files.push_back(file)

bkg_root_files = ROOT.std.vector('string')()
for file in bkg_filenames.readlines(): 
    file = file.replace("\n", "")
    bkg_root_files.push_back(file)

sig_df = ROOT.RDataFrame("DecayTree", sig_root_files)
bkg_df = ROOT.RDataFrame("DecayTree", bkg_root_files)
# print((sig_df.Count()).GetValue())

sig = sig_df.AsNumpy(branch_names)
bkg = bkg_df.AsNumpy(branch_names)

# Transform the input features
def make_classification(sig, bkg, output, cut, first_step, draw_input_features, cross_validation):
    signal = pd.DataFrame()
    background = pd.DataFrame()

    if(first_step):
        input_features = ['VTXISODCHI2ONETRACK_B', 'VTXISODCHI2MASSTWOTRACK_B', 'CCNC_20_IT_B', 'VTXISODCHI2TWOTRACK_tau_max', 'VTXISODMASSONETRACK_tau_min', 'CC_20_PTASYM_tau_min', 'NC_20_PZASYM_tau_max']

        signal['VTXISODCHI2ONETRACK_B'] = sig['Bp_VTXISODCHI2ONETRACK_B'] 
        signal['VTXISODCHI2MASSTWOTRACK_B'] = sig['Bp_VTXISODCHI2MASSTWOTRACK_B']
        signal['CCNC_20_IT_B'] = sig['Bp_CCNC_20_IT_B']
        signal['VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( sig['Bp_VTXISODCHI2MASSONETRACK_taup'], sig['Bp_VTXISODCHI2TWOTRACK_taum'] )
        signal['VTXISODMASSONETRACK_tau_min'] = np.minimum( sig['Bp_VTXISODCHI2MASSONETRACK_taup'], sig['Bp_VTXISODCHI2MASSONETRACK_taum'] )
        signal['CC_20_PTASYM_tau_min'] = np.minimum( sig['Bp_CC_20_PTASYM_taup'], sig['Bp_CC_20_PTASYM_taum'] )
        signal['NC_20_PZASYM_tau_max'] = np.maximum( sig['Bp_NC_20_PZASYM_taup'], sig['Bp_NC_20_PZASYM_taum'] )

        background['VTXISODCHI2ONETRACK_B'] = bkg['Bp_VTXISODCHI2ONETRACK_B']
        background['VTXISODCHI2MASSTWOTRACK_B'] = bkg['Bp_VTXISODCHI2MASSTWOTRACK_B']
        background['CCNC_20_IT_B'] = bkg['Bp_CCNC_20_IT_B']
        background['VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( bkg['Bp_VTXISODCHI2MASSONETRACK_taup'], bkg['Bp_VTXISODCHI2TWOTRACK_taum'] )
        background['VTXISODMASSONETRACK_tau_min'] = np.minimum( bkg['Bp_VTXISODCHI2MASSONETRACK_taup'], bkg['Bp_VTXISODCHI2MASSONETRACK_taum'] )
        background['CC_20_PTASYM_tau_min'] = np.minimum( bkg['Bp_CC_20_PTASYM_taup'], bkg['Bp_CC_20_PTASYM_taum'] )
        background['NC_20_PZASYM_tau_max'] = np.maximum( bkg['Bp_NC_20_PZASYM_taup'], bkg['Bp_NC_20_PZASYM_taum'] )
    else:

        signal['tau_M12_max'] = np.maximum( sig['taup_M12'], sig['taum_M12'] )
        signal['tau_M23_max'] = np.minimum( sig['taup_M23'], sig['taum_M23'] )
        signal['B_M13'] = sig['Bp_M13']
        signal['log_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX'] ) )
        signal['log_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log(1 - np.abs(sig['taup_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX']),  np.log(1 - np.abs(sig['taum_DIRA_ORIVX'] ))*np.sign( sig['taup_DIRA_ORIVX'] ) )
        signal['log_1_minus_B_DIRA_PV'] = np.log(1 - np.abs(sig['Bp_DIRA_OWNPV'] ))*np.sign( sig['Bp_DIRA_OWNPV'])

        Cx_taup_sig =  (sig['df_DV1y'] - sig['df_BVy'])*sig['df_Kp_PZ']  - ( sig['df_DV1z'] - sig['df_BVz'])*sig['df_Kp_PY']
        Cy_taup_sig =  (sig['df_DV1z'] - sig['df_BVz'])*sig['df_Kp_PX']  - ( sig['df_DV1x'] - sig['df_BVx'])*sig['df_Kp_PZ']
        Cz_taup_sig =  (sig['df_DV1x'] - sig['df_BVx'])*sig['df_Kp_PY']  - ( sig['df_DV1y'] - sig['df_BVy'])*sig['df_Kp_PX']
        C_taup_sig = np.sqrt( Cx_taup_sig**2 + Cy_taup_sig**2 + Cz_taup_sig**2  )
        signal['IP_taup_Kp'] = (2*C_taup_sig)/( np.sqrt( sig['df_Kp_PX']**2 + sig['df_Kp_PY']**2 + sig['df_Kp_PZ']**2 ) )

        Cx_taum_sig =  (sig['df_DV2y'] - sig['df_BVy'])*sig['df_Kp_PZ']  - ( sig['df_DV2z'] - sig['df_BVz'])*sig['df_Kp_PY']
        Cy_taum_sig =  (sig['df_DV2z'] - sig['df_BVz'])*sig['df_Kp_PX']  - ( sig['df_DV2x'] - sig['df_BVx'])*sig['df_Kp_PZ']
        Cz_taum_sig =  (sig['df_DV2x'] - sig['df_BVx'])*sig['df_Kp_PY']  - ( sig['df_DV2y'] - sig['df_BVy'])*sig['df_Kp_PX']
        C_taum_sig = np.sqrt( Cx_taum_sig**2 + Cy_taum_sig**2 + Cz_taum_sig**2  )
        signal['IP_taum_Kp'] = (2*C_taum_sig)/( np.sqrt( sig['df_Kp_PX']**2 + sig['df_Kp_PY']**2 + sig['df_Kp_PZ']**2 ) )

        a_sig = np.sqrt( ( sig['df_PVx'] - sig['df_DV1x'] )**2 + ( sig['df_PVy'] - sig['df_DV1y'] )**2 + ( sig['df_PVz'] - sig['df_DV1z'] )**2 )
        b_sig = np.sqrt( ( sig['df_PVx'] - sig['df_DV2x'] )**2 + ( sig['df_PVy'] - sig['df_DV2y'] )**2 + ( sig['df_PVz'] - sig['df_DV2z'] )**2 )
        c_sig = np.sqrt( ( sig['df_DV1x'] - sig['df_DV2x'] )**2 + ( sig['df_DV1y'] - sig['df_DV2y'] )**2 + ( sig['df_DV1z'] - sig['df_DV2z'] )**2 )
        s_sig = (a_sig+b_sig+c_sig)/2
        signal['DV1_DV2_PV_area'] = np.sqrt( s_sig + (s_sig-a_sig)*(s_sig-b_sig)*(s_sig-c_sig) )

        background['tau_M12_max'] = np.maximum( bkg['taup_M12'], bkg['taum_M12'] )
        background['tau_M23_max'] = np.minimum( bkg['taup_M23'], bkg['taum_M23'] )
        background['B_M13'] = bkg['Bp_M13']
        background['log_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX'] ) )
        background['log_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log(1 - np.abs(bkg['taup_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX']),  np.log(1 - np.abs(bkg['taum_DIRA_ORIVX'] ))*np.sign( bkg['taup_DIRA_ORIVX'] ) )
        background['log_1_minus_B_DIRA_PV'] = np.log(1 - np.abs(bkg['Bp_DIRA_OWNPV'] ))*np.sign( bkg['Bp_DIRA_OWNPV'])

        Cx_taup_bkg =  (bkg['df_DV1y'] - bkg['df_BVy'])*bkg['df_Kp_PZ']  - ( bkg['df_DV1z'] - bkg['df_BVz'])*bkg['df_Kp_PY']
        Cy_taup_bkg =  (bkg['df_DV1z'] - bkg['df_BVz'])*bkg['df_Kp_PX']  - ( bkg['df_DV1x'] - bkg['df_BVx'])*bkg['df_Kp_PZ']
        Cz_taup_bkg =  (bkg['df_DV1x'] - bkg['df_BVx'])*bkg['df_Kp_PY']  - ( bkg['df_DV1y'] - bkg['df_BVy'])*bkg['df_Kp_PX']
        C_taup_bkg = np.sqrt( Cx_taup_bkg**2 + Cy_taup_bkg**2 + Cz_taup_bkg**2  )
        background['IP_taup_Kp'] = (2*C_taup_bkg)/( np.sqrt( bkg['df_Kp_PX']**2 + bkg['df_Kp_PY']**2 + bkg['df_Kp_PZ']**2 ) )

        Cx_taum_bkg =  (bkg['df_DV2y'] - bkg['df_BVy'])*bkg['df_Kp_PZ']  - ( bkg['df_DV2z'] - bkg['df_BVz'])*bkg['df_Kp_PY']
        Cy_taum_bkg =  (bkg['df_DV2z'] - bkg['df_BVz'])*bkg['df_Kp_PX']  - ( bkg['df_DV2x'] - bkg['df_BVx'])*bkg['df_Kp_PZ']
        Cz_taum_bkg =  (bkg['df_DV2x'] - bkg['df_BVx'])*bkg['df_Kp_PY']  - ( bkg['df_DV2y'] - bkg['df_BVy'])*bkg['df_Kp_PX']
        C_taum_bkg = np.sqrt( Cx_taum_bkg**2 + Cy_taum_bkg**2 + Cz_taum_bkg**2  )
        background['IP_taum_Kp'] = (2*C_taum_bkg)/( np.sqrt( bkg['df_Kp_PX']**2 + bkg['df_Kp_PY']**2 + bkg['df_Kp_PZ']**2 ) )

        a_bkg = np.sqrt( ( bkg['df_PVx'] - bkg['df_DV1x'] )**2 + ( bkg['df_PVy'] - bkg['df_DV1y'] )**2 + ( bkg['df_PVz'] - bkg['df_DV1z'] )**2 )
        b_bkg = np.sqrt( ( bkg['df_PVx'] - bkg['df_DV2x'] )**2 + ( bkg['df_PVy'] - bkg['df_DV2y'] )**2 + ( bkg['df_PVz'] - bkg['df_DV2z'] )**2 )
        c_bkg = np.sqrt( ( bkg['df_DV1x'] - bkg['df_DV2x'] )**2 + ( bkg['df_DV1y'] - bkg['df_DV2y'] )**2 + ( bkg['df_DV1z'] - bkg['df_DV2z'] )**2 )
        s_bkg = (a_bkg+b_bkg+c_bkg)/2
        background['DV1_DV2_PV_area'] = np.sqrt( s_bkg + (s_bkg-a_bkg)*(s_bkg-b_bkg)*(s_bkg-c_bkg) )

    # Add mc weights to signal dataframe (and ones to background)
    mc_weights_filename = '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_10/tree_with_weights.root'.format(year)
    mc_weights_df = ROOT.RDataFrame("DecayTree", mc_weights_filename)
    # print((mc_weights_df.Count()).GetValue())
    #print(mc_weights_df.GetColumnNames())
    mc_weights_dict = mc_weights_df.AsNumpy(['weight'])
    mc_weights = mc_weights_dict['weight']

    mc_weights_pandas = pd.DataFrame(mc_weights_dict)
    signal = pd.concat([signal,mc_weights_pandas], axis=1)

    data_weights = np.ones((bkg_df.Count()).GetValue())
    data_weights_pandas = pd.DataFrame({'weight': data_weights})
    background = pd.concat([background,data_weights_pandas], axis=1)

    # add BDT output from 1st step and apply weight
    if(not first_step):
        # bdt_output_values = pd.DataFrame({'bdt': output})

        print(len(data_weights)+len(mc_weights))
        print(len(output))

        print(output)

        quit()

    # for sklearn data is usually organised into one 2D array of shape (n_samples x n_features)
    # containing all the data and one array of categories of length n_samples
    X = np.concatenate((signal, background))
    y = np.concatenate((np.ones(signal.shape[0]), np.zeros(background.shape[0])))

    data = pd.DataFrame(np.hstack((X, y.reshape(y.shape[0], -1))),columns=input_features+['weight']+['y'])

    # Compute number of signal and background events under the signal peak (+/- 2 sigma)
    hist_settings = {'bins': 100, 'density': True, 'alpha': 0.4} 
    n, bins, patches = plt.hist(sig['df_Bp_M'], color='b', weights=mc_weights, range=[4000,8000], **hist_settings, label="Signal")
    plt.hist(bkg['df_Bp_M'], color='r', range=[4000,8000], **hist_settings, label="Background")
    plt.legend()
    plt.xlabel('Fitted mass (MeV)')
    plt.ylabel('Normalised entries / (40 MeV)')

    fwhm, left_idx, right_idx = FWHM(bins, n)
    fwhm = fwhm[0]
    sigma = fwhm/2.4
    mass_mean = bins[np.where( n == max(n) )[0][0]]

    r1 = mass_mean - 2*sigma
    r2 = mass_mean + 2*sigma

    plt.axvspan(r1, r2, facecolor='g', alpha=0.3)
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots_first_step/mass.pdf'.format(year))
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Input_features_plots_second_step/mass.pdf'.format(year))
    plt.clf()

    N_sig = ((sig_df.Filter("(df_Bp_M > {0}) && (df_Bp_M < {1})".format(r1, r2))).Count()).GetValue()
    N_bkg = ((bkg_df.Filter("(df_Bp_M > {0}) && (df_Bp_M < {1})".format(r1, r2))).Count()).GetValue()

    # Plotting variables and correlations
    if(draw_input_features):
        print("Drawing variables and correlations")
        draw_variables(signal, background, mc_weights, input_features, first_step)
        data_scatter = data.drop('weight', axis=1)
        draw_scatter_plots(data_scatter, input_features, first_step)
        sig_corr = signal.drop('weight', axis=1)
        bkg_corr = background.drop('weight', axis=1)
        correlation_matrix(sig_corr, "sig", first_step)
        correlation_matrix(bkg_corr, "bkg", first_step)
        quit()

    # Train and test split
    X_dev,X_eval, y_dev,y_eval = train_test_split(X, y, test_size=0.33, random_state=42)
    X_train,X_test, y_train,y_test = train_test_split(X_dev, y_dev, test_size=0.33, random_state=42)

    X_dev_weights = X_dev[:,-1]
    X_eval_weights = X_eval[:, -1]
    X_train_weights = X_train[:, -1]
    X_test_weights = X_test[:, -1]

    X_dev = np.delete(X_dev, -1, 1)
    X_eval = np.delete(X_eval, -1, 1)
    X_train = np.delete(X_train, -1, 1)
    X_test = np.delete(X_test, -1, 1)
    X = np.delete(X, -1, 1)

    ##############################################################################################################################################

    ###################################################### Train a classifier ####################################################################
    names = ["DTree", "RForest", "AdaBoost"]

    # I am using the optimised hyper-parameters here:
    # classifiers = best_classifiers
    classifiers = [ DecisionTreeClassifier(max_depth=5),
                    RandomForestClassifier(max_depth=6, max_features=3, n_estimators=1000, random_state=42),
                    AdaBoostClassifier(learning_rate=1.2, n_estimators=1000, algorithm='SAMME', random_state=42)] # by default the estimator is a decision tree with max_depth = 1

    best_classifiers = []

    clf_fpr = []
    clf_tpr = []
    clf_thresholds = []

    clf_output = []

    for i in range(len(names)):
        name = names[i]
        clf = classifiers[i]

        print("Fitting "+names[i])
        clf.fit(X_train, y_train, X_train_weights)

        if(cross_validation):
            do_cross_validation(name, clf, X_dev, y_dev, X_eval, y_eval, X_dev_weights) 
        else:
            y_predicted = clf.predict(X_test)
            print(classification_report(y_test, y_predicted, target_names=["background", "signal"]))

            print("Drawing confusion matrix")
            confusion_matrix(name, clf, X_test, y_test, first_step)

            print("Comparing train and test sets")
            compare_train_test(name, clf, X_train, y_train, X_test, y_test, first_step) 

            # print("Drawing feature importances")
            # draw_feature_importance(name, clf)

    if(cross_validation):
        quit()

    # Compare different classifiers
    print("Drawing roc curve")
    fpr, tpr, thresholds = roc_curve_plot(names, classifiers, X_test, y_test, first_step)

    cuts = []
    metrics = []

    print("Drawing signal significance vs BDT cut")
    for i in range(len(fpr)):
        bdt_cut, significance = draw_signal_significance_vs_bdt_cut(names[i], fpr[i], tpr[i], thresholds[i], N_sig, N_bkg, first_step)
        cuts.append(bdt_cut)
        metrics.append(significance)

    # find which classifier gives the maximum value for the metric
    idx = np.argmax( metrics )

    # compute classifier response on signal and background datasets
    if((names[idx] == "DTree") or (names[idx] == "RForest")):
        decisions = classifiers[idx].predict_proba(X)[:, 1]
        d1 = classifiers[idx].predict_proba(X[y>0.5])[:,1].ravel()
        d2 = classifiers[idx].predict_proba(X[y<0.5])[:,1].ravel()
    else:
        decisions = classifiers[idx].decision_function(X)
        d1 = classifiers[idx].decision_function(X[y>0.5]).ravel()
        d2 = classifiers[idx].decision_function(X[y<0.5]).ravel()

    # plot classifier response on signal and background datasets  
    dummy = []
    dummy += [d1,d2]      
    low = min(np.min(d) for d in dummy)
    high = max(np.max(d) for d in dummy)
    low_high = (low,high)
    
    plt.hist(dummy[0],
             color='r', alpha=0.5, range=low_high, bins=30,
             histtype='stepfilled', density=True,
             label='S')
    plt.hist(dummy[1],
             color='b', alpha=0.5, range=low_high, bins=30,
             histtype='stepfilled', density=True,
             label='B')
    plt.xlabel('Classifier output')
    plt.ylabel('Arbitrary units')
    plt.legend()
    if(first_step):
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_first_step/'.format(year)+'clf_response.pdf')
    else:
        plt.savefig('/panfs/felician/B2Ktautau/workflow/Sklearn/201{0}/Output_performance_second_step/'.format(year)+'clf_response.pdf')
    plt.clf()

    return decisions, cuts[idx]

print("FIRST STEP: ISOLATION - TRAINING")
first_step_output, first_step_cut = make_classification(sig, bkg, output=None, cut=None, first_step=True, draw_input_features=False, cross_validation=False)

# print("SECOND STEP: GEOMETRY + KINEMATICS")
make_classification(sig, bkg, output=first_step_output, cut=first_step_cut, first_step=False, draw_input_features=True, cross_validation=False)

# TO DO: APPLY CUT AND TRAIN 2ND STEP