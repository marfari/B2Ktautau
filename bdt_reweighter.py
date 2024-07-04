import numpy
#import root_numpy
import ROOT
import pandas
from hep_ml import reweight
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from hep_ml.metrics_utils import ks_2samp_weighted
import sys
from array import array
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import RocCurveDisplay

def main(argv):

    year = argv[1]
    species = argv[2]

    mc_filename = '/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/201{0}/Species_7/tree.root'.format(year)
    data_filename = '/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/201{0}/Species_8/signal_tree.root'.format(year)

    # Open data and MC files and define set of variables that we want to use for re-weighting (columns)
    columns = ['Bp_PT', 'Bp_ETA', 'Bp_ENDVERTEX_Z', 'nTracks', 'nSPDHits', 'nPVs']

    mc_df = ROOT.RDataFrame("DatasetTree", mc_filename)
    data_df = ROOT.RDataFrame("DatasetTree", data_filename)

    # print((mc_df.Count()).GetValue())

    original = mc_df.AsNumpy(columns)
    target = data_df.AsNumpy(columns)

    original = pandas.DataFrame(original)
    target = pandas.DataFrame(target)

    # Splot weights (signal)
    target_weights_dict = data_df.AsNumpy(['n_signal_sw'])
    target_weights = target_weights_dict['n_signal_sw']

    original_weights = numpy.ones(len(original)) # stores the weights to be applied to MC

    # Apply sWeights to target dataframe
    # prepare train and test samples
    original_train, original_test = train_test_split(original) # divide original samples into training ant test parts
    target_train, target_test = train_test_split(target) # divide target samples into training ant test parts

    original_weights_train = numpy.ones(len(original_train))
    original_weights_test = numpy.ones(len(original_test))

    target_weights_train = numpy.zeros(len(target_train))
    target_weights_test = numpy.zeros(len(target_test))

    # Dive splot weights for train and test samples
    target_train_indices = target_train.index
    for i in range(len(target_train)):
        idx = target_train_indices[i]
        target_weights_train[i] = target_weights[idx]

    target_test_indices = target_test.index
    for i in range(len(target_test)):
        idx = target_test_indices[i]
        target_weights_test[i] = target_weights[idx]

    def draw_distributions(columns, original, target, new_original_weights, data_signal_weights, flag):
        hist_settings = {'bins': 50, 'density': True, 'alpha': 0.7} # dendity = True means the histograms are normalised
        plt.figure(figsize=[15, 10])
        for id, column in enumerate(columns, 1):
            xlim = numpy.percentile(numpy.hstack([target[column]]), [0.01, 99.99])
            # plt.subplot(3, 3, id)
            plt.hist(original[column], weights=new_original_weights, range=xlim, **hist_settings, label="Signal in MC")
            plt.hist(target[column], weights=data_signal_weights, range=xlim, **hist_settings, label="Signal in data")
            plt.legend(loc="upper left", fontsize=20)
            ks_val = ks_2samp_weighted(original[column], target[column], weights1=new_original_weights, weights2=numpy.ones(len(target), dtype=float))
            plt.title('KS over {0} = {1}'.format(column,ks_val), fontsize=15)
            print('KS over ', column, ' = ', ks_val)
            if flag == 0:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/MC_unweighted/{2}.pdf'.format(year,species,id-1))
            elif flag == 1:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/MC_bins_reweighter/{2}.pdf'.format(year,species,id-1))
            elif flag == 2:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/MC_GBReweighter/{2}.pdf'.format(year,species,id-1))
            elif flag == 3:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/MC_unweighted_all/{2}.pdf'.format(year,species,id-1))
            elif flag == 4:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/MC_GBReweighter_all/{2}.pdf'.format(year,species,id-1))
            elif flag == 5:
                plt.savefig('/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/MC_DDs_reweighted/{2}.pdf'.format(year,species,id-1))
            plt.clf()
                
    # 1) original distributions
    print("Plotting unweighted feature distributions")
    draw_distributions(columns, original, target, original_weights, target_weights, 0)

    # 2) Bins based re-weighting in N dimensions = histogram reweighting (w/ BDT output) -> only works for a small number of variables
    # bins_reweighter = reweight.BinsReweighter(n_bins=20, n_neighs=1.)
    # bins_reweighter.fit(original_train, target_train)

    # bins_weights_test = bins_reweighter.predict_weights(original_test)
    # draw_distributions(columns, original_test, target_test, bins_weights_test, 1) # validate reweighting rule on the test part comparing 1d projections

    # Gradient Boosted Reweighter (GBReweighter)
    #reweighter = reweight.GBReweighter(n_estimators=50, learning_rate=0.1, max_depth=3, min_samples_leaf=1000, gb_args={'subsample': 0.4})
    #reweighter.fit(original_train, target_train, original_weights_train, target_weights_train)

    #gb_weights_test = reweighter.predict_weights(original_test)
    # print("Plotting re-weighted feature distributions")
    # draw_distributions(columns, original_test, target_test, gb_weights_test, target_weights_test, 2) # validate reweighting rule on the test part comparing 1d projections

    # Folding reweighter
    # define base reweighter
    reweighter_base = reweight.GBReweighter(n_estimators=100, learning_rate=0.5, max_depth=3, min_samples_leaf=1000, gb_args={'subsample': 0.4})
    reweighter = reweight.FoldingReweighter(reweighter_base, n_folds=2)
    reweighter.fit(original, target, original_weights, target_weights) # it is not needed divide data into train/test parts; rewighter can be train on the whole samples
    # predict method provides unbiased weights prediction for the whole sample
    # folding reweighter contains two reweighters, each is trained on one half of samples
    # during predictions each reweighter predicts another half of samples not used in training
    folding_weights = reweighter.predict_weights(original)

    print("Plotting re-weighted feature distributions")
    draw_distributions(columns, original, target, folding_weights, target_weights, 2)

    # Testing quality of re-weighter: how well can a classifier distinguish between the re-weighted MC and the data (sWeight)?
    data = numpy.concatenate([original, target])
    labels = numpy.array([0] * len(original) + [1] * len(target))

    weights = {}
    weights['original'] = original_weights
    weights['2-folding'] = folding_weights

    for name, new_weights in weights.items():
        W = numpy.concatenate([new_weights / new_weights.sum() * len(target), [1] * len(target)])
        Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(data, labels, W, random_state=42, train_size=0.51)
        clf = GradientBoostingClassifier(subsample=0.3, n_estimators=30)
        clf_fit = clf.fit(Xtr, Ytr, sample_weight=Wtr)
    
    print(name, roc_auc_score(Yts, clf_fit.predict_proba(Xts)[:, 1], sample_weight=Wts))

    # Plot all variables before re-weighting the MC:
    original_all = mc_df.AsNumpy()
    target_all = data_df.AsNumpy()
    original_all = pandas.DataFrame(original_all)
    target_all = pandas.DataFrame(target_all)
    original_weights_all = numpy.ones(len(original_all))
    columns_all = original_all.columns
    print("Plotting all unweighted distributions")
    draw_distributions(columns_all, original_all, target_all, original_weights_all, target_weights, 3)

    print("Plotting re-weighted all distributions")
    draw_distributions(columns_all, original_all, target_all, folding_weights, target_weights, 4)

    # Correct MC (species 7)
    fc_DDs_mc = ROOT.TFileCollection("fc_DDs_mc", "fc_DDs_mc", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_{1}/pre_sel_tree.txt".format(year,species))
    t_DDs_mc = ROOT.TChain("DecayTree")
    t_DDs_mc.AddFileInfoList(fc_DDs_mc.GetList())
    df_DDs_mc = ROOT.RDataFrame(t_DDs_mc)
    original_DDs_mc = df_DDs_mc.AsNumpy(columns)
    original_DDs_mc = pandas.DataFrame(original_DDs_mc) # pandas dataframe: x (events in DDs MC) y (input features)

    folding_weights_DDs_mc = reweighter.predict_weights(original_DDs_mc)

    file_DDs = ROOT.TFile.Open("/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{0}/DDs_correction/Species_{1}/tree_with_weights.root".format(year,species), "RECREATE")
    tree_DDs = ROOT.TTree("DecayTree", "DecayTree")
    weights_var_DDs = array('d', [0])
    tree_DDs.Branch('weight', weights_var_DDs, 'weights/D')

    N_DDs = len(folding_weights_DDs_mc) # number of events

    for i in range(N_DDs):
        weights_var_DDs[0] = folding_weights_DDs_mc[i]
        tree_DDs.Fill()

    file_DDs.cd()
    tree_DDs.Write()
    file_DDs.Close()

    print(folding_weights_DDs_mc)

    print("Plotting re-weighted all distributions")
    draw_distributions(columns, original_DDs_mc, target, folding_weights_DDs_mc, target_weights, 5)

if __name__ == "__main__":
    main(sys.argv)