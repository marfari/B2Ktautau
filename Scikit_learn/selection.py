import os, re, pickle
import sklearn.utils
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.metrics import RocCurveDisplay
from sklearn.inspection import permutation_importance
from sklearn.model_selection import cross_val_score

####### FUNCTIONS ##############################################################################
def plot_input_variable(sample, variable, bins, var_min, var_max, var_label, year, dir_name):
    signal = sample[sample.Class==1]
    background = sample[sample.Class==0]

    plt.hist(signal[variable], bins=bins, range=[var_min,var_max], histtype='stepfilled', lw=2, alpha=0.5, color='blue', label=[r'Signal'], density=True)
    plt.hist(background[variable], bins=bins, range=[var_min,var_max], histtype='stepfilled', lw=2, alpha=0.5, color='red', label=[r'Background'], density=True)

    plt.ylabel('Normalised Entries ('+str(bins)+' bins)')
    plt.xlabel(var_label)

    plt.legend(loc="upper right")

    plt.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Input_variables/'+dir_name+'/'+variable+'.png')
    plt.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Input_variables/'+dir_name+'/'+variable+'.pdf')
    plt.clf()

def plot_scatter_plot_matrix(sample, variables_to_model, year, dir_name):
    signal = sample[sample.Class==1]
    background = sample[sample.Class==0]

    signal_scatter_plot_matrix = sns.pairplot(signal[variables_to_model],corner=True, kind='reg', diag_kind='kde', plot_kws={'line_kws':{'color':'red'}})
    signal_scatter_plot_matrix.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/scatter_plot_matrix_signal.png') 
    signal_scatter_plot_matrix.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/scatter_plot_matrix_signal.pdf') 
    plt.clf()

    background_scatter_plot_matrix = sns.pairplot(background[variables_to_model],corner=True, kind='reg', diag_kind='kde', plot_kws={'line_kws':{'color':'red'}})
    background_scatter_plot_matrix.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/scatter_plot_matrix_background.png') 
    background_scatter_plot_matrix.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/scatter_plot_matrix_background.pdf') 
    plt.clf()

def plot_heatmap(sample, variables_to_model, year, dir_name):
    signal = sample[sample.Class==1]
    background = sample[sample.Class==0]

    signal_heatmap = sns.heatmap(signal[variables_to_model].corr(),annot=True, annot_kws={"size": 25})
    signal_heatmap.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/heatmap_signal.png')
    signal_heatmap.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/heatmap_signal.pdf')
    plt.clf()
    
    background_heatmap = sns.heatmap(background[variables_to_model].corr(),annot=True, annot_kws={"size": 25})
    background_heatmap.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/heatmap_background.png')
    background_heatmap.figure.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/Scatter_plots/'+dir_name+'/heatmap_background.pdf')
    plt.clf()

def rank_input_variables(names, classifiers, X_train, y_train, X_test, y_test, year, dir_name):
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    # Ranking of input variables (w/ permutation feature importance)
    r = permutation_importance(clf, X_test, y_test,n_repeats=30,random_state=0)
    perm_sorted_idx = r.importances_mean.argsort()
    ax.boxplot(r.importances[perm_sorted_idx].T,vert=False,labels=X_train.columns[perm_sorted_idx],)
    # Mean impurity decrease feature importance (only works for tree-based models)
    # importances = pd.DataFrame(data={'Attribute': X_train.columns,'Importance': clf.feature_importances_})
    # importances = importances.sort_values(by='Importance', ascending=False)
    # plt.bar(x=importances['Attribute'], height=importances['Importance'], color='#087E8B')
    # plt.title('Ranking of input variables (from permutation)', size=10)
    # plt.xticks(rotation='vertical')
    plt.xlabel('Drop in accuracy')
    fig.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/'+name+'_ranking_input_vars.png')
    fig.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/'+name+'_ranking_input_vars.pdf')
    plt.clf()

def compare_train_test(clf, X_train, y_train, X_test, y_test, bins=30):
    fig2, ax2 = plt.subplots()
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):
        d1 = clf.decision_function(X[y>0.5]).ravel()
        d2 = clf.decision_function(X[y<0.5]).ravel()
        decisions += [d1, d2]
        
    #low = min(np.min(d) for d in decisions)
    #high = max(np.max(d) for d in decisions)
    low = -0.25
    high = 0.25
    low_high = (low,high)
    
    ax2.hist(decisions[0],
             color='r', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='S (train)')
    ax2.hist(decisions[1],
             color='b', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='B (train)')

    hist, bins = np.histogram(decisions[2],
                              bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax2.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')
    
    hist, bins = np.histogram(decisions[3],
                              bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    ax2.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')

    plt.xlabel(name+" output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='best')
    fig2.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/'+name+'_overtraining_check.png')
    fig2.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/'+name+'_overtraining_check.pdf')
    plt.clf()

def plot_train_test(clf, name, sig_train, bkg_train, sig_test, bkg_test, year, dir_name):
    fig2, ax2 = plt.subplots()

    # Train outputs
    ax2.hist(sig_train[name+" output"],color='r', alpha=0.5, histtype='stepfilled', density=True, label='S (train)')
    ax2.hist(bkg_train[name+" output"],color='b', alpha=0.5, histtype='stepfilled', density=True, label='B (train)')
    

    # Test output
    hist, bins = np.histogram(sig_test[name+" output"], density=True)
    scale = len(sig_test[name+" output"]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax2.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')

    hist, bins = np.histogram(bkg_test[name+" output"], density=True)
    scale = len(bkg_test[name+" output"]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax2.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')

    plt.xlabel(name+" output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='best')
    fig2.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/'+name+'_overtraining_check.png')
    fig2.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/'+name+'_overtraining_check.pdf')
    plt.clf()


#############################################################################################################

year = '2018'
print("YEAR: "+year)
isWS = True
variables = 0
# 0 -> tau- kinematic
# 1 -> tau+ kinematic
# 2 -> tau- isolation
# 3 -> tau+ isolation
# 4 -> B+ kinematic
# 5 -> B+ isolation
plot = False

if isWS:
    data_sign = 'WS'
    print("BACKGROUND PROXY: WS data")
else:
    data_sign = 'RS'
    print("BACKGROUND PROXY: RS data, right sideband")

if variables == 0:
    selection_type = 'taum_kinematic'
    dir_name = 'Taum_kinematic'
    print("VARIABLES USED: tau- kinematic")
elif variables == 1:
    selection_type = 'taup_kinematic'
    dir_name = 'Taup_kinematic'
    print("VARIABLES USED: tau+ kinematic")
elif variables == 2:
    selection_type = 'taum_isolation'
    dir_name = 'Taum_isolation'
    print("VARIABLES USED: tau- isolation")
elif variables == 3:
    selection_type = 'taup_isolation'
    dir_name = 'Taup_isolation'
    print("VARIABLES USED: tau+ isolation")
elif variables == 4:
    selection_type = 'B_kinematic'
    dir_name = 'Bp_kinematic'
    print("VARIABLES USED: B+ kinematic")
elif variables == 5:
    selection_type = 'B_isolation'
    dir_name = 'Bp_isolation'
    print("VARIABLES USED: B+ isolation")

# Open data and MC dataframes
df_data_variables = pd.read_csv('Pandas/data_'+data_sign+'_'+year+'_'+selection_type+'.csv')
df_data_selections = pd.read_csv('Pandas/data_'+data_sign+'_'+year+'_forSelections.csv')
duplicates = list( set(df_data_variables.columns).intersection(df_data_selections.columns)  )
df_data_selections = df_data_selections.drop(duplicates,axis=1)
df_data = pd.concat([df_data_variables,df_data_selections],axis=1)

df_mc_variables = pd.read_csv('Pandas/mc_'+year+'_'+selection_type+'.csv')
df_mc_selections = pd.read_csv('Pandas/mc_'+year+'_forSelections.csv')
duplicates = list( set(df_mc_variables.columns).intersection(df_mc_selections.columns)  )
df_mc_selections = df_mc_selections.drop(duplicates,axis=1)
df_mc = pd.concat([df_mc_variables,df_mc_selections],axis=1)

# Manage Pandas DataFrames
## Do specific pre-selections
### Mass window

if not isWS:
    preselData = (df_data.Bp_ConsBp_0_M>6500)
    df_data = df_data[preselData]
preselMC = (df_mc.Bp_ConsBp_0_M>4500) & (df_mc.Bp_ConsBp_0_M<6000)
df_mc = df_mc[preselMC]

## Create column to assign to which category each dataset belongs
df_data['Class'] = 0 # background
df_mc['Class'] = 1 # signal

## Add signal and background datasets
### make sure the two have the same columns first

sample = pd.concat([df_data,df_mc])

## Do common pre-selections
passDTF = (sample.Bp_ConsBp_0_status==0)

L0_trigger = (sample.Bp_L0HadronDecision_TOS==1) | ((sample.Bp_L0HadronDecision_TIS==1) | (sample.Bp_L0MuonDecision_TIS==1) | (sample.Bp_L0ElectronDecision_TIS==1) | (sample.Bp_L0PhotonDecision_TIS==1))
HLT1_trigger = (sample.Bp_Hlt1TrackMVADecision_TOS==1) | (sample.Bp_Hlt1TwoTrackMVADecision_TOS==1)
HLT2_trigger = (sample.Bp_Hlt2Topo2BodyDecision_TOS==1) | (sample.Bp_Hlt2Topo3BodyDecision_TOS==1) | (sample.Bp_Hlt2Topo4BodyDecision_TOS==1)
trigger = L0_trigger & HLT1_trigger & HLT2_trigger

others = (sample.Bp_ConsBp_0_decayLength>5) & (sample.Bp_ConsBp_0_tauminus_0_decayLength>0.5) & (sample.Bp_ConsBp_0_tauminus_decayLength>0.5) & (sample.taup_DIRA_ORIVX > 0.99) & (sample.taum_DIRA_ORIVX > 0.99) & (sample.Bp_DIRA_OWNPV > 0.99) & (np.sqrt( (sample.taup_ENDVERTEX_X - sample.taum_ENDVERTEX_X)**2 + (sample.taup_ENDVERTEX_Y - sample.taum_ENDVERTEX_Y)**2 + (sample.taup_ENDVERTEX_Z - sample.taum_ENDVERTEX_Z)**2 ) > 0.5) & (sample.Bp_ConsBp_0_chi2 < 20)

preselCommon = passDTF & trigger & others

sample = sample[preselCommon]

##########################################################################################################################################
# Training

if variables == 0:
    sample['DTF_taum_decayLength_sig'] = (sample['Bp_ConsBp_0_tauminus_0_decayLengthErr']/sample['Bp_ConsBp_0_tauminus_0_decayLength'])
    sample['log(1-taum_DIRA_ORIVX)'] = np.log(1 - sample['taum_DIRA_ORIVX'])
    sample['log(1-taum_DIRA_OWNPV)'] = np.log(1 - sample['taum_DIRA_OWNPV'])
    sample['taum_daughters_min_IPCHI2_OWNPV'] = sample[['taum_pi1_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV']].min(axis=1)
    sample['Bp_ConsBp_0_tauminus_0_PT'] = np.sqrt( sample['Bp_ConsBp_0_tauminus_0_PX']**2 + sample['Bp_ConsBp_0_tauminus_0_PY']**2)

    variables_to_model = ['taum_M', 'taum_M12', 'taum_M13', 'taum_M23', 'taum_ENDVERTEX_CHI2', 'DTF_taum_decayLength_sig', 'log(1-taum_DIRA_ORIVX)', 'log(1-taum_DIRA_OWNPV)', 'taum_daughters_min_IPCHI2_OWNPV', 'Bp_ConsBp_0_tauminus_0_PE', 'Bp_ConsBp_0_tauminus_0_PT', 'Bp_ConsBp_0_tauminus_0_PZ']
    var_min = [800, 200, 200, 200, 0, 0, -25, -20, 0, 0, 0, 0]
    var_max = [1600, 1400, 1600, 1600, 20, 5, 0, 0, 1500, 600000, 16000, 600000]
    var_label = [r'Visible m($\tau^{-}$) (MeV)', r'Visible $\tau^{-}$ m($\pi^{+} \pi^{-}$) (MeV)', r'Visible $\tau^{-}$ m($\pi^{-} \pi^{-}$) (MeV)', r'Visible $\tau^{-}$ m($\pi^{+} \pi^{-}$) (MeV)', r'Visible $\tau^{-}$ DV $\chi^{2}$', r'DTF $\tau^{-}$ decay length significance', r'log(1 - $\tau^{-}$ DIRA to BV)', r'log(1 - $\tau^{-}$ DIRA to PV)', r'$\tau^{-}$ min($\pi_{1}$, $\pi_{2}$, $\pi_{3}$) IP $\chi^{2}$', r'DTF E($\tau^{-}$) (MeV)', r'DTF PT($\tau^{-}$) (MeV)', r'DTF PZ($\tau^{-}$) (MeV)']
elif variables == 1:
    variables_to_model = ['taup_M', 'taup_M12', 'taup_M13', 'taup_M23', 'taup_ENDVERTEX_CHI2', 'Bp_ConsBp_0_tauminus_decayLength', 'Bp_ConsBp_0_tauminus_decayLengthErr', 'taup_DIRA_ORIVX', 'taup_DIRA_OWNPV', 'taup_pi1_IPCHI2_OWNPV', 'taup_pi2_IPCHI2_OWNPV', 'taup_pi3_IPCHI2_OWNPV', 'Bp_ConsBp_0_tauminus_PE', 'Bp_ConsBp_0_tauminus_PX', 'Bp_ConsBp_0_tauminus_PY', 'Bp_ConsBp_0_tauminus_PZ']
elif variables == 2:
    variables_to_model = ['Bp_CC_IT_tau1', 'Bp_CC_MAXPT_PE_tau1', 'Bp_CC_MAXPT_PT_tau1', 'Bp_CC_MAXPT_Q_tau1', 'Bp_CC_MULT_tau1', 'Bp_CC_PZ_tau1', 'Bp_CC_PZASYM_tau1', 'Bp_CC_SPT_tau1', 'Bp_CC_VPT_tau1', 'Bp_NC_DELTAETA_tau1', 'Bp_NC_DELTAPHI_tau1', 'Bp_NC_IT_tau1', 'Bp_NC_MAXPT_PT_tau1', 'Bp_NC_MAXPT_PZ_tau1', 'Bp_NC_MULT_tau1', 'Bp_NC_PTASYM_tau1', 'Bp_NC_PZ_tau1', 'Bp_NC_SPT_tau1', 'Bp_NC_VPT_tau1', 'Bp_VTXISODCHI2ONETRACK_tau1', 'Bp_VTXISODCHI2TWOTRACK_tau1', 'Bp_VTXISONUMVTX_tau1']
elif variables == 3: 
    variables_to_model = ['Bp_CC_IT_tau2', 'Bp_CC_MAXPT_PE_tau2', 'Bp_CC_MAXPT_PT_tau2', 'Bp_CC_MAXPT_Q_tau2', 'Bp_CC_MULT_tau2', 'Bp_CC_PZ_tau2', 'Bp_CC_PZASYM_tau2', 'Bp_CC_SPT_tau2', 'Bp_CC_VPT_tau2', 'Bp_NC_DELTAETA_tau2', 'Bp_NC_DELTAPHI_tau2', 'Bp_NC_IT_tau2', 'Bp_NC_MAXPT_PT_tau2', 'Bp_NC_MAXPT_PZ_tau2', 'Bp_NC_MULT_tau2', 'Bp_NC_PTASYM_tau2', 'Bp_NC_PZ_tau2', 'Bp_NC_SPT_tau2', 'Bp_NC_VPT_tau2', 'Bp_VTXISODCHI2ONETRACK_tau2', 'Bp_VTXISODCHI2TWOTRACK_tau2', 'Bp_VTXISONUMVTX_tau2']
elif variables == 4:
    variables_to_model = ['Bp_ConsBp_0_PX', 'Bp_ConsBp_0_PY[09', 'Bp_ConsBp_0_PZ', 'Bp_ConsBp_0_decayLength', 'Bp_ConsBp_0_decayLengthErr', 'Bp_ETA', 'Bp_M01', 'Bp_M02', 'Bp_M03', 'Bp_M04', 'Bp_M05', 'Bp_M06', 'Bp_M012', 'Bp_M013', 'Bp_M014', 'Bp_M015', 'Bp_M016', 'Bp_M023', 'Bp_M024', 'Bp_M025', 'Bp_M026', 'Bp_M034', 'Bp_M035', 'Bp_M036', 'Bp_M045', 'Bp_M046', 'Bp_M056']
elif variables == 5:
    variables_to_model = ['Bp_CC_IT_B', 'Bp_CC_MAXPT_PE_B', 'Bp_CC_MAXPT_PT_B', 'Bp_CC_MAXPT_Q_B', 'Bp_CC_MULT_B', 'Bp_CC_PZ_B', 'Bp_CC_PZASYM_B', 'Bp_CC_SPT_B', 'Bp_CC_VPT_B', 'Bp_NC_DELTAETA_B', 'Bp_NC_DELTAPHI_B', 'Bp_NC_IT_B', 'Bp_NC_MAXPT_PT_B', 'Bp_NC_MAXPT_PZ_B', 'Bp_NC_MULT_B', 'Bp_NC_PTASYM_B', 'Bp_NC_PZ_B', 'Bp_NC_SPT_B', 'Bp_NC_VPT_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISONUMVTX_B']  

if plot:
    ## Plot input variables
    print("Plotting input variables")
    for i in range(len(variables_to_model)):
        plot_input_variable(sample, variables_to_model[i], 50, var_min[i], var_max[i], var_label[i], year, dir_name)

    ## Scatter plot matrix
    print("Plotting scatter plots")
    plot_scatter_plot_matrix(sample, variables_to_model, year, dir_name)

    ## Heatmap 
    plot_heatmap(sample, variables_to_model, year, dir_name)

## Shuffle the events
sample = sklearn.utils.shuffle(sample, random_state=123) #'123' is the random seed

## Split in test and train datasets
Ntrain_stop = int(round(sample.shape[0] * 0.7)) # 2/3 of the sample are used for training, 1/3 is used for testing

sample_Train = sample[:Ntrain_stop]
sample_Test = sample[Ntrain_stop:]

## Build X and y arrays for the training
X_train = sample_Train[variables_to_model]
y_train = sample_Train[['Class']]

X_test = sample_Test[variables_to_model]
y_test = sample_Test[['Class']]

## Classifiers
names = [
    #"Nearest Neighbors",
    #"Linear SVM",
    #"RBF SVM",
    #"Gaussian Process",
    "Decision Tree",
    "Random Forest",
    "Neural Network",
    "AdaBoost",
    #"Naive Bayes",
    #"QDA",
]

classifiers = [
    #KNeighborsClassifier(3),
    #SVC(kernel="linear", C=0.025),
    #SVC(gamma=2, C=1),
    #GaussianProcessClassifier(1.0 * RBF(1.0)),
    DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    MLPClassifier(alpha=1, max_iter=1000),
    AdaBoostClassifier(),
    #GaussianNB(),
    #QuadraticDiscriminantAnalysis(),
]

i = 0
fig1, ax1 = plt.subplots()
for clf, name in zip(classifiers,names):
    print("Classifier: ", name)
    clf.fit(X_train, y_train['Class'])
    score = clf.score(X_test, y_test) # mean accuracy on the given test data and labels
    print("Accuracy on test data: {:.2f}".format(score))

    # rank of input variables
    rank_input_variables(name, clf, X_train, y_train, X_test, y_test, year, dir_name)
    
    # classifier output
    clf_output_train_sig = clf.predict_proba(X_train[y_train.Class==1])[:,1]
    clf_output_train_bkg = clf.predict_proba(X_train[y_train.Class==0])[:,1]
    clf_output_test_sig = clf.predict_proba(X_test[y_test.Class==1])[:,1]
    clf_output_test_bkg = clf.predict_proba(X_test[y_test.Class==0])[:,1]

    sig_train = pd.DataFrame(clf_output_train_sig, columns=[name+' output'])
    bkg_train = pd.DataFrame(clf_output_train_bkg, columns=[name+' output'])
    sig_test = pd.DataFrame(clf_output_test_sig, columns=[name+' output'])
    bkg_test = pd.DataFrame(clf_output_test_bkg, columns=[name+' output'])

    plot_train_test(clf, name, sig_train, bkg_train, sig_test, bkg_test, year, dir_name)

    # if(name == 'AdaBoost'):
    #     compare_train_test(clf, X_train, y_train, X_test, y_test, bins=30)

    # ROC curve
    clf_roc = RocCurveDisplay.from_estimator(clf, X_test, y_test, ax=ax1, alpha=0.8)

ax1.grid()
ax1.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
fig1.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/ROC_curve.png')
fig1.savefig('/home/felician/B2Ktautau/Scikit_learn/Plots/'+year+'/ML_output/'+dir_name+'/ROC_curve.pdf')

