import sys
import ROOT
import numpy as np
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    bdt1 = float(bdt1)
    bdt2 = float(bdt2)

    # Gex fixed yield inputs
    # The only term that is missing is the numerator of the background efficiency (only term that depends on the BDT cuts)
    fixed_input_values = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_values_inputs.npy')
    fixed_input_errors = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_errors_inputs.npy')
    fixed_input_values = unumpy.uarray( fixed_input_values, fixed_input_errors )
    # print(fixed_input_values)

    # Get numerators of the background efficiency
    N_den_bkg = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_den.npy')

    species_name = ["$B^+ \\to \\bar{D}^0 D^0 K^+$", "$B^+ \\to D^+ D^- K^+$", "$B^+ \\to D^+_s D^-_s K^+$", "$B^0 \\to D^- D^0 K^+$", "$B^0_s \\to D^-_s D^0 K^+$", "$B^+ \\to \\bar{D}^0 D^+ K^0$", "$B^+ \\to \\bar{D}^0 D^+_s$", "$B^+ \\to \\bar{D}^0 D^+$"]
    component_name = ["$DD$", "$D^*D$", "$DD^*$", "$D^*D^*$", "Total"]
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]

    N_species = len(all_species)
    N_components = len(component_name)

    # expected number of signal events at BF = 10-3
    branching_fraction = 0.001
    fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy")
    fixed_value_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs_error.npy')
    fixed_value = ufloat( fixed_value, fixed_value_err )
    N_den = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_s_den.npy')
    # print(fixed_value)

    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")
    N_num = t_sig.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))

    upper = ROOT.TEfficiency.Wilson(N_den, N_num, 0.68, True)
    lower = ROOT.TEfficiency.Wilson(N_den, N_num, 0.68, False)
    eps_s_err = 0.5*(upper - lower)

    eps_s = ufloat( N_num/N_den, eps_s_err )

    S = branching_fraction/(fixed_value/eps_s)
    # print(S)

    all_yields = pd.DataFrame(index=species_name, columns=component_name)
    yields = np.zeros((N_species,N_components))
    yields_errors = np.zeros((N_species,N_components))
    # relative to the expected number of signal events at BF = 10-3
    all_relative_yields = pd.DataFrame(index=species_name, columns=component_name)
    relative_yields = np.zeros((N_species,N_components))
    relative_yields_errors = np.zeros((N_species,N_components))

    N_tot = np.zeros(N_species)
    N_tot_err = np.zeros(N_species)

    N_tot_relative = np.zeros(N_species)
    N_tot_err_relative = np.zeros(N_species)

    for i in range(N_species):
        species = all_species[i]
        n = 0
        n_rel = 0

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
            if(j < 4):
                N_num_bkg = t.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (species == {2}) && (component == {3})".format(bdt1, bdt2, species, j))
                
                eps_b = N_num_bkg/N_den_bkg[i,j]

                upper_bound = ROOT.TEfficiency.Wilson(N_den_bkg[i,j], N_num_bkg , 0.68, True)
                lower_bound = ROOT.TEfficiency.Wilson(N_den_bkg[i,j], N_num_bkg, 0.68, False)
                eps_b_error = 0.5*(upper_bound - lower_bound)

                eps_b_post_acc = ufloat(eps_b, eps_b_error)
                
                # Absolute yields
                ufloat_yields = fixed_input_values[i,j]*eps_b_post_acc

                yields[i,j] = unumpy.nominal_values(ufloat_yields)
                yields_errors[i,j] = unumpy.std_devs(ufloat_yields)
                all_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$".format( int(yields[i,j]), int(yields_errors[i,j]) )  

                # Relative yields
                if(S == 0):
                    ufloat_relative = ufloat(0,0)
                else:
                    ufloat_relative = ufloat_yields/S
                relative_yields[i,j] = unumpy.nominal_values(ufloat_relative)
                relative_yields_errors[i,j] = unumpy.std_devs(ufloat_relative)
                all_relative_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$ \\%".format( int(relative_yields[i,j]*100), int(relative_yields_errors[i,j]*100) )  

                n += ufloat_yields
                n_rel += ufloat_relative
            else:
                N_tot[i] = unumpy.nominal_values(n)
                N_tot_err[i] = unumpy.std_devs(n)

                N_tot_relative[i] = unumpy.nominal_values(n_rel)
                N_tot_err_relative[i] = unumpy.std_devs(n_rel)
                
                all_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$".format( int(N_tot[i]), int(N_tot_err[i]) )  

                all_relative_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$ \\%".format( int(N_tot_relative[i]*100), int(N_tot_err_relative[i]*100) )  

    # print(N_tot)
    # print(N_tot_relative)

    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), yields)
    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_errors_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), yields_errors)

    with open('/panfs/felician/B2Ktautau/workflow/yield_estimates/yields_table_bdt1_{0}_bdt2_{1}.tex'.format( bdt1, bdt2 ),'w') as fout:
        fout.write(all_yields.to_latex())

    with open('/panfs/felician/B2Ktautau/workflow/yield_estimates/relative_yields_table_bdt1_{0}_bdt2_{1}.tex'.format( bdt1, bdt2 ),'w') as fout1:
        fout1.write(all_relative_yields.to_latex())
    
if __name__ == "__main__":
    main(sys.argv)