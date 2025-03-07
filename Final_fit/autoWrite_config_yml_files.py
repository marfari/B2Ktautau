import ROOT
import numpy as np
import sys, os
sys.path.append(os.getcwd())

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t = f.Get("DecayTree")

    eps_ws_den = t.GetEntries()
    eps_ws_num = t.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws = eps_ws_num/eps_ws_den

    N_bkg = 2*19006409*eps_ws
    N_bkg_error = np.sqrt(N_bkg)

    string = '''
General:
Measurement: "minimal_example"
POI: "Nsig"
HistogramFolder: "/panfs/felician/B2Ktautau/workflow/final_fit/histograms/"

Regions:
- Name: "mass_fit_region"
    Variable: "df_Bp_M"
    # Filter: "lep_charge > 0"
    Binning: [4000, 4080, 4160, 4240, 4320, 4400, 4480, 4560, 4640, 4720, 4800, 4880, 4960, 5040, 5120, 5200, 5280, 5360, 5440, 5520, 5600, 5680, 5760, 5840, 5920, 6000, 6080, 6160, 6240, 6320, 6400, 6480, 6560, 6640, 6720, 6800, 6880, 6960, 7040, 7120, 7200, 7280, 7360, 7440, 7520, 7600, 7680, 7760, 7840, 7920, 8000]

Samples:
- Name: "h_toy_data"
    SamplePath: '/panfs/felician/B2Ktautau/workflow/generate_toy_data/toy_data_bdt1_{0}_bdt2_{1}_seed_{2}.root'
    Data: True
    
- Name: "h_sig"
    SamplePath: "/panfs/felician/B2Ktautau/workflow/create_fit_templates/fit_templates_bdt1_{0}_bdt2_{1}.root"
    # Weight: "weight"

- Name: "h_bkg"
    SamplePath: "/panfs/felician/B2Ktautau/workflow/create_fit_templates/fit_templates_bdt1_{0}_bdt2_{1}.root"
    # Weight: "weight"

NormFactors:
- Name: "Nsig"
    Samples: "h_sig"
    Nominal: 0
    Bounds: [-10, 10]

- Name: "Nbkg"
    Samples: "h_bkg"
    Nominal: {3}
    Bounds: [{4}, {5}]

Systematics:
- Name: "Luminosity"
    Up:
    Normalization: 0.05
    Down:
    Normalization: -0.05
    Type: "Normalization"

# - Name: "Modeling"
#   Up:
#     SamplePath: "prediction.root"
#     Tree: "background_varied"
#   Down:
#     Symmetrize: True
#   Samples: "Background"
#   Type: "NormPlusShape"

# - Name: "WeightBasedModeling"
#   Up:
#     Weight: "weight_up"
#   Down:
#     Weight: "0.7*weight"
#   Samples: "Background"
#   Type: "NormPlusShape"
    '''.format(bdt1, bdt2, random_seed, N_bkg, N_bkg-N_bkg_error, N_bkg+N_bkg_error)

    file = open('/panfs/felician/B2Ktautau/workflow/config_yaml_files/config_histograms_bdt1_{0}_bdt2_{1}_seed_{2}.yml'.format(bdt1, bdt2, random_seed), 'w')
    file.write(string)
    file.close()

if __name__ == "__main__":
    main(sys.argv)