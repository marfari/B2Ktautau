import sys
import ROOT
import numpy as np
import pandas as pd

def eps_error(Num, Den):
    return (Num/Den)*np.sqrt( 1/Num + 1/Den )

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    # Gex fixed yield inputs
    # The only term that is missing is the numerator of the background efficiency (only term that depends on the BDT cuts)
    fixed_input_values = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_values_inputs.npy')
    fixed_input_errors = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_errors_inputs.npy')

    species_name = ["$B^+ \\to \\bar{D}^0 D^0 K^+$", "$B^+ \\to D^+ D^- K^+$", "$B^+ \\to D^+_s D^-_s K^+$", "$B^0 \\to D^- D^0 K^+$", "$B^0_s \\to D^-_s D^0 K^+$", "$B^+ \\to \\bar{D}^0 D^+ K^0$", "$B^+ \\to \\bar{D}^0 D^+_s$", "$B^+ \\to \\bar{D}^0 D^+$"]
    component_name = ["$DD$", "$D^*D$", "$DD^*$", "$D^*D^*$"]
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]

    N_species = len(all_species)
    N_components = len(component_name)

    all_yields = pd.DataFrame(index=species_name, columns=component_name)
    yields = np.zeros((N_species,N_components))
    yields_errors = np.zeros((N_species,N_components))

    for i in range(N_species):
        species = all_species[i]

        if((species == 100) or (species == 101) or (species == 102)):
            f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt1_0_bdt2_0.root")
        elif(species == 110):
            f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt1_0_bdt2_0.root")
        elif(species == 120):
            f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt1_0_bdt2_0.root")
        elif(species == 130):
            f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt1_0_bdt2_0.root")
        elif((species == 150) or (species == 151)):
            f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt1_0_bdt2_0.root")

        t = f.Get("DecayTree")

        for j in range(N_components):
            N_num_bkg = t.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (species == {2}) && (component == {3})".format(bdt1, bdt2, species, j))
            N_num_bkg_error = np.sqrt(N_num_bkg)
            
            yields[i,j] = fixed_input_values[i,j]*N_num_bkg
            yields_errors[i,j] = np.sqrt( (N_num_bkg**2)*(fixed_input_errors[i,j]**2) + (fixed_input_values[i,j]**2)*(N_num_bkg_error**2) )

            # Relative yield
            all_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$".format( int(yields[i,j]), int(yields_errors[i,j]) )  

    print(yields)
    print(yields_errors)

    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{0}_bdt2_{1}.npy'.format( round(bdt1,1), round(bdt2,1) ), yields)
    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_errors_bdt1_{0}_bdt2_{1}.npy'.format( round(bdt1,1), round(bdt2,1) ), yields_errors)

    with open('/panfs/felician/B2Ktautau/workflow/yield_estimates/yields_table_bdt1_{0}_bdt2_{1}.tex'.format( round(bdt1,1), round(bdt2,1) ),'w') as fout:
        fout.write(all_yields.to_latex())
    
if __name__ == "__main__":
    main(sys.argv)