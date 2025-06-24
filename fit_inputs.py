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
    # A is a value
    A_const = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const.npy')
    A_const_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const_err.npy')
    A_const = unumpy.uarray(A_const, A_const_err)

    # B is a value
    B = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy')
    B_err = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy')
    B = ufloat(B, B_err)

    # C is a matrix: row = decay; column = component (DD, D*D, DD*, D*D*)
    C_const = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/C_const.npy')
    C_const_err = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/C_const_err.npy')
    C_const = unumpy.uarray(C_const, C_const_err)

    # Efficiency denominators 
    eps_s_den = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_s_den.npy')
    eps_b_den = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_den.npy')

    ########################################################################### A value ########################################################################################
    f_sig = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_sig = f_sig.Get("DecayTree")
    eps_s_num = t_sig.GetEntries(f"(BDT1 > {bdt1}) && (BDT2 > {bdt2}) ")

    upper_sig = ROOT.TEfficiency.Wilson(eps_s_den, eps_s_num, 0.68, True)
    lower_sig = ROOT.TEfficiency.Wilson(eps_s_den, eps_s_num, 0.68, False)
    eps_s_err = 0.5*(upper_sig - lower_sig)

    eps_s = ufloat(eps_s_num/eps_s_den, eps_s_err)

    if(eps_s.nominal_value != 0):
        A = A_const/eps_s
    else:
        A = ufloat(0.0, 0.0)
    np.save(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/A_bdt1_{bdt1}_bdt2_{bdt2}.npy', A.nominal_value)
    np.save(f'/panfs/felician/B2Ktautau/workflow/fit_inputs/A_err_bdt1_{bdt1}_bdt2_{bdt2}.npy', A.std_dev)

    ########################################################################### Yield estimates | C values ########################################################################################
    species_name = ["$B^+ \\to \\bar{D}^0 D^0 K^+$", "$B^+ \\to D^+ D^- K^+$", "$B^+ \\to D^+_s D^-_s K^+$", "$B^0 \\to D^- D^0 K^+$", "$B^0_s \\to D^-_s D^0 K^+$", "$B^+ \\to \\bar{D}^0 D^+ K^0$", "$B^+ \\to \\bar{D}^0 D^+_s$", "$B^+ \\to \\bar{D}^0 D^+$"]
    component_name = ["$DD$", "$D^*D$", "$DD^*$", "$D^*D^*$", "Total"]
    component_name_1 = ["$DD$", "$D^*D$", "$DD^*$", "$D^*D^*$"]
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]

    N_species = len(all_species)
    N_components = 4

    C_values = np.zeros((N_species,N_components))
    C_errors = np.zeros((N_species,N_components))

    yields = np.zeros((N_species,N_components))
    yields_errors = np.zeros((N_species,N_components))

    eps_bkg_values = np.zeros((N_species,N_components))
    eps_bkg_errors = np.zeros((N_species,N_components))

    # Yields and total background efficiency tables
    all_yields = pd.DataFrame(index=species_name, columns=component_name)
    all_eps_b = pd.DataFrame(index=species_name, columns=component_name_1)

    # This is needed to build the table of the total background efficiencies (eps_acc*eps_strip*eps_reco)
    eps_bkg_acc_strip_reco = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_acc_strip_reco.npy')
    eps_bkg_acc_strip_reco_err = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_acc_strip_reco_err.npy')
    eps_bkg_acc_strip_reco = unumpy.uarray(eps_bkg_acc_strip_reco, eps_bkg_acc_strip_reco)

    N_tot = np.zeros(N_species)
    N_tot_err = np.zeros(N_species)

    for i in range(N_species):
        species = all_species[i]
        n = 0

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

        for j in range(len(component_name)):
            if(j < 4):
                eps_b_num = t.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (species == {2}) && (component == {3})".format(bdt1, bdt2, species, j))
                eps_b = eps_b_num/eps_b_den[i,j]

                upper_bound = ROOT.TEfficiency.Wilson(eps_b_den[i,j], eps_b_num , 0.68, True)
                lower_bound = ROOT.TEfficiency.Wilson(eps_b_den[i,j], eps_b_num, 0.68, False)
                eps_b_error = 0.5*(upper_bound - lower_bound)
                eps_b_post_acc = ufloat(eps_b, eps_b_error)
                ufloat_eps_bkg = eps_bkg_acc_strip_reco[i,j]*eps_b_post_acc
                eps_bkg_values[i,j] = unumpy.nominal_values(ufloat_eps_bkg)
                eps_bkg_errors[i,j] = unumpy.std_devs(ufloat_eps_bkg)
                
                # Absolute yields
                ufloat_yields = (C_const[i,j]*eps_b_post_acc)/B # yield = C / B
                yields[i,j] = unumpy.nominal_values(ufloat_yields)
                yields_errors[i,j] = unumpy.std_devs(ufloat_yields)

                ufloat_c = C_const[i,j]*eps_b_post_acc
                C_values[i,j] = unumpy.nominal_values(ufloat_c)
                C_errors[i,j] = unumpy.std_devs(ufloat_c)

                all_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$".format( round(yields[i,j],1), round(yields_errors[i,j],1) )  
                all_eps_b.loc[species_name[i], component_name[j]] = "{0}".format( ufloat_eps_bkg )  

                n += ufloat_yields
            else:
                N_tot[i] = unumpy.nominal_values(n)
                N_tot_err[i] = unumpy.std_devs(n)

                all_yields.loc[species_name[i], component_name[j]] = "${0} \\pm {1}$".format( round(N_tot[i],1), round(N_tot_err[i],1) )  

    # print(N_tot)

    np.save('/panfs/felician/B2Ktautau/workflow/fit_inputs/C_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), C_values)
    np.save('/panfs/felician/B2Ktautau/workflow/fit_inputs/C_err_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), C_errors)

    np.save('/panfs/felician/B2Ktautau/workflow/fit_inputs/yield_values_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), yields)
    np.save('/panfs/felician/B2Ktautau/workflow/fit_inputs/yield_errors_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), yields_errors)

    np.save('/panfs/felician/B2Ktautau/workflow/fit_inputs/eps_bkg_values_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), eps_bkg_values)
    np.save('/panfs/felician/B2Ktautau/workflow/fit_inputs/eps_bkg_errors_bdt1_{0}_bdt2_{1}.npy'.format( bdt1, bdt2 ), eps_bkg_errors)

    with open('/panfs/felician/B2Ktautau/workflow/fit_inputs/yields_table_bdt1_{0}_bdt2_{1}.tex'.format( bdt1, bdt2 ),'w') as fout:
        fout.write(all_yields.to_latex())   

    with open('/panfs/felician/B2Ktautau/workflow/fit_inputs/eps_b_table_bdt1_{0}_bdt2_{1}.tex'.format( bdt1, bdt2 ),'w') as fout:
        fout.write(all_eps_b.to_latex())    
    
if __name__ == "__main__":
    main(sys.argv)