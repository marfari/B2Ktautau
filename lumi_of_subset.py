import numpy as np
import pandas as pd
from ROOT import TFile, TTree
import sys
from array import array

def main(argv):

    runs_in_subset = np.loadtxt('/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.txt'.format(year=argv[1]))

    df = pd.read_csv("run2_pp_lumi.csv")

    all_runs = df['run'].values
    all_lumi = df['lumi'].values

    file = TFile('/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root'.format(year=argv[1]), 'recreate')
    file.cd()

    tree = TTree("DecayTree", "DecayTree")
    lumi_in_subset  = 0

    for i in range(len(all_runs)):
        for j in range(len(runs_in_subset)):
            if(all_runs[i] == runs_in_subset[j]):
                lumi_in_subset += all_lumi[i] 

    print(lumi_in_subset)
    
    tree.Branch("lumi_in_subset",  lumi_in_subset,  'lumi_in_subset/D')
    tree.Fill()
    tree.Write()
    file.Close()

    return

if __name__ == "__main__":
    main(sys.argv)

    