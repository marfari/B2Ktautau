import sys
import ROOT
import numpy as np
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy

def main(argv):

    ############################################### B (normalisation mode parameter) ###################################################
    ### Normalisation model signal efficiency
    f_norm_mc = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_72/post_sel_tree_bdt_0.0.root")
    t_norm_mc = f_norm_mc.Get("DecayTree")

    # Gen norm MC (for N_norm_den)
    fc_gen_norm_2016 = ROOT.TFileCollection("fc_gen_norm_2016", "fc_gen_norm_2016", "Files_on_grid/MC_D0Dps_2016.txt")
    fc_gen_norm_2017 = ROOT.TFileCollection("fc_gen_norm_2017", "fc_gen_norm_2017", "Files_on_grid/MC_D0Dps_2017.txt")
    fc_gen_norm_2018 = ROOT.TFileCollection("fc_gen_norm_2018", "fc_gen_norm_2018", "Files_on_grid/MC_D0Dps_2018.txt")

    t_norm_gen = ROOT.TChain("mc_ntuple/MCDecayTree")
    t_norm_gen_2017 = ROOT.TChain("mc_ntuple/MCDecayTree")
    t_norm_gen_2018 = ROOT.TChain("mc_ntuple/MCDecayTree")

    t_norm_gen.AddFileInfoList(fc_gen_norm_2016.GetList())
    t_norm_gen_2017.AddFileInfoList(fc_gen_norm_2017.GetList())
    t_norm_gen_2018.AddFileInfoList(fc_gen_norm_2018.GetList())

    print("Norm")
    print(t_norm_gen.GetEntries())
    print(t_norm_gen_2017.GetEntries())
    print(t_norm_gen_2018.GetEntries())

    t_norm_gen.Add(t_norm_gen_2017)
    t_norm_gen.Add(t_norm_gen_2018)

    # Norm acceptance
    eps_acc_norm_2016 = ufloat(14.50/100, 0.10/100)
    eps_acc_norm_2017 = ufloat(14.48/100, 0.10/100)
    eps_acc_norm_2018 = ufloat(14.58/100, 0.10/100)
    eps_acc_norm = (10.8/34.2)*eps_acc_norm_2016 + (12./34.2)*eps_acc_norm_2017 + (11.4/34.2)*eps_acc_norm_2018

    N_norm_num = t_norm_mc.GetEntries()
    # N_norm_den_2016 = 5350785
    # N_norm_den_2017 = 5782016
    # N_norm_den_2018 = 5763580
    # N_norm_den = N_norm_den_2016 + N_norm_den_2017 + N_norm_den_2018
    N_norm_den = t_norm_gen.GetEntries()

    upper_bound = ROOT.TEfficiency.Wilson(N_norm_den, N_norm_num, 0.68, True)
    lower_bound = ROOT.TEfficiency.Wilson(N_norm_den, N_norm_num, 0.68, False)
    eps_norm_err = 0.5*(upper_bound - lower_bound)

    eps_norm_post_acc = ufloat( N_norm_num/N_norm_den, eps_norm_err )

    eps_norm = eps_acc_norm*eps_norm_post_acc
    print("eps_norm = ", eps_norm)

    ### Normalisation mode signal yield
    f_norm_fit = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/fit_mass/Species_8/mass_fit_result.root")
    norm_fit_result = f_norm_fit.Get("fitresult_model_data")
    n_norm_fit_var = norm_fit_result.floatParsFinal().find("n_signal")
    N_norm_value = n_norm_fit_var.getVal()
    N_norm_error = n_norm_fit_var.getError()
    
    N_norm = ufloat(N_norm_value, N_norm_error)
    print("N_norm = ", N_norm)

    ### PDG branching fractions
    B_Bp_D0bar_Dsp = ufloat(9.0/1000, 0.9/1000)
    B_D0bar_Kp_pi = ufloat(3.947/100, 0.030/100)
    B_Dsp_Kp_Km_pi = ufloat(5.37/100, 0.10/100)

    B = (N_norm/eps_norm)/(B_Bp_D0bar_Dsp*B_D0bar_Kp_pi*B_Dsp_Kp_Km_pi)
    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B.npy', B.nominal_value)
    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/B_err.npy', B.std_dev)
    ###########################################################################################################################################

    ############################################ Ktautau A_fixed #####################################################
    ### Branching fractions
    B_tau_3pi_nu = ufloat(9.31/100, 0.05/100)
    B_tau_3pi_pi0_nu = ufloat(4.62/100, 0.05/100)

    ### Signal efficiency (fixed)
    fc_gen_2016 = ROOT.TFileCollection("fc_gen_2016", "fc_gen_2016", "Files_on_grid/MC_2016.txt")
    fc_gen_2017 = ROOT.TFileCollection("fc_gen_2017", "fc_gen_2017", "Files_on_grid/MC_2017.txt")
    fc_gen_2018 = ROOT.TFileCollection("fc_gen_2018", "fc_gen_2018", "Files_on_grid/MC_2018.txt")

    # 2016
    t_gen_3pi3pi = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
    t_gen_3pi3pi_pi0_2016 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
    t_gen_3pipi0_3pi_2016 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
    t_gen_3pi3pi2pi0_2016 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")
    t_gen_all = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")

    t_gen_3pi3pi.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_3pi3pi_pi0_2016.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_3pipi0_3pi_2016.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_3pi3pi2pi0_2016.AddFileInfoList(fc_gen_2016.GetList())
    t_gen_all.AddFileInfoList(fc_gen_2016.GetList())

    print("2016")
    print(t_gen_3pi3pi.GetEntries())
    print(t_gen_3pi3pi_pi0_2016.GetEntries())
    print(t_gen_3pipi0_3pi_2016.GetEntries())
    print(t_gen_3pi3pi2pi0_2016.GetEntries())
    print(t_gen_all.GetEntries())

    t_gen_all.Add(t_gen_3pi3pi_pi0_2016)
    t_gen_all.Add(t_gen_3pipi0_3pi_2016)
    t_gen_all.Add(t_gen_3pi3pi2pi0_2016)

    # 2017
    t_gen_3pi3pi_2017 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
    t_gen_3pi3pi_pi0_2017 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
    t_gen_3pipi0_3pi_2017 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
    t_gen_3pi3pi2pi0_2017 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")
    t_gen_all_2017 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")

    t_gen_3pi3pi_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_3pi3pi_pi0_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_3pipi0_3pi_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_3pi3pi2pi0_2017.AddFileInfoList(fc_gen_2017.GetList())
    t_gen_all_2017.AddFileInfoList(fc_gen_2017.GetList())

    print("2017")
    print(t_gen_3pi3pi_2017.GetEntries()) 
    print(t_gen_3pi3pi_pi0_2017.GetEntries())
    print(t_gen_3pipi0_3pi_2017.GetEntries())
    print(t_gen_3pi3pi2pi0_2017.GetEntries())
    print(t_gen_all_2017.GetEntries())

    t_gen_all_2017.Add(t_gen_3pi3pi_pi0_2017)
    t_gen_all_2017.Add(t_gen_3pipi0_3pi_2017)
    t_gen_all_2017.Add(t_gen_3pi3pi2pi0_2017)

    # 2018
    t_gen_3pi3pi_2018 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
    t_gen_3pi3pi_pi0_2018 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
    t_gen_3pipi0_3pi_2018 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
    t_gen_3pi3pi2pi0_2018 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")
    t_gen_all_2018 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")

    t_gen_3pi3pi_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_3pi3pi_pi0_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_3pipi0_3pi_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_3pi3pi2pi0_2018.AddFileInfoList(fc_gen_2018.GetList())
    t_gen_all_2018.AddFileInfoList(fc_gen_2018.GetList())

    print("2018")
    print(t_gen_3pi3pi_2018.GetEntries())
    print(t_gen_3pi3pi_pi0_2018.GetEntries())
    print(t_gen_3pipi0_3pi_2018.GetEntries())
    print(t_gen_3pi3pi2pi0_2018.GetEntries())
    print(t_gen_all_2018.GetEntries())

    t_gen_all_2018.Add(t_gen_3pi3pi_pi0_2018)
    t_gen_all_2018.Add(t_gen_3pipi0_3pi_2018)
    t_gen_all_2018.Add(t_gen_3pi3pi2pi0_2018)

    t_gen_3pi3pi.Add(t_gen_3pi3pi_2017)
    t_gen_3pi3pi.Add(t_gen_3pi3pi_2018)

    t_gen_all.Add(t_gen_all_2017)
    t_gen_all.Add(t_gen_all_2018)

    # Ktautau acceptance
    eps_acc_2016 = ufloat(5.318/100, 0.010/100)
    eps_acc_2017 = ufloat(5.321/100, 0.011/100)
    eps_acc_2018 = ufloat(5.330/100, 0.010/100)
    eps_acc = (28.5/165)*eps_acc_2016 + (60.0/165)*eps_acc_2017 + (76.5/165)*eps_acc_2018

    # Ktautau stripping 
    eps_strip_2016 = ufloat(1.028/100, 0.003/100)
    eps_strip_2017 = ufloat(1.056/100, 0.003/100)
    eps_strip_2018 = ufloat(1.050/100, 0.002/100)
    eps_strip = (28.5/165)*eps_strip_2016 + (60.0/165)*eps_strip_2017 + (76.5/165)*eps_strip_2018

    ### Ktautua reco1 = 100%

    ### Ktautau eps_gen
    # N_gen_2016 = 46774+19242+19252+8028
    # N_gen_2017 = 61870+25333+25516+10549
    # N_gen_2018 = 95913+39754+39423+16649
    # N_gen = N_gen_2016+N_gen_2017+N_gen_2018
    N_gen = t_gen_all.GetEntries()

    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/Ns_gen.npy', N_gen)

    A_fix = ((B_tau_3pi_nu+B_tau_3pi_pi0_nu)**2)*(eps_acc*eps_strip)
    print("A_fix = ", A_fix)
    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/A_fix_value.npy', A_fix.nominal_value)
    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/A_fix_error.npy', A_fix.std_dev)
    ######################################################################################################################################

    ############################################### Physics backgrounds efficiency (fixed) ###############################################
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]
    all_components = [0, 1, 2, 3]

    N_species = len(all_species)
    N_components = len(all_components)

    # branching fractions
    # B(B -> Da Db X)
    Bvalues = pd.read_csv('/home/felician/B2Ktautau/Final_fit/physics_bkgs_branching_fractions/B_values.csv')
    Berrors = pd.read_csv('/home/felician/B2Ktautau/Final_fit/physics_bkgs_branching_fractions/B_errors.csv')

    Bvalues = Bvalues.transpose()
    Bvalues.reset_index(drop=True, inplace=True)
    Bvalues = Bvalues.rename(columns=Bvalues.iloc[0]).drop(Bvalues.index[0])
    Bvalues.columns = Bvalues.columns.astype('int')
    Bvalues = Bvalues.reset_index(drop=True)
    Bvalues = Bvalues.dropna()

    Berrors = Berrors.transpose()
    Berrors.reset_index(drop=True, inplace=True)
    Berrors = Berrors.rename(columns=Berrors.iloc[0]).drop(Berrors.index[0])
    Berrors.columns = Berrors.columns.astype('int')
    Berrors = Berrors.reset_index(drop=True)
    Berrors = Berrors.dropna()

    B_Dp_3piX = ufloat( 15.25/100, 0.20/100 )
    B_D0_3piX = ufloat( 17.60/100, 0.25/100 )
    B_Ds_3piX = ufloat( 32.81/100, 0.72/100 )

    B_Dstar0_D0X_1 = ufloat( 64.7/100, 0.9/100 ) # D0* -> D0 pi0
    B_Dstar0_D0X_2 = ufloat( 35.3/100, 0.9/100 ) # D*0 -> D0 gamma
    B_Dstar0_D0X = B_Dstar0_D0X_1 + B_Dstar0_D0X_2

    B_DstarP_D0pi = ufloat( 67.7/100, 0.5/100 ) # D*+ -> D0 pi+

    B_DstarP_DpX_1 = ufloat( 30.7/100, 0.5/100 ) # D*+ -> D+ pi0
    B_DstarP_DpX_2 = ufloat( 1.6/100, 0.4/100 ) # D*+ -> D+ gamma
    B_DstarP_DpX = B_DstarP_DpX_1 + B_DstarP_DpX_2

    B_DstarS_DsX_1 = ufloat( 93.6/100, 0.4/100 ) # Ds*+ -> Ds+ gamma
    B_DstarS_DsX_2 = ufloat( 5.77/100, 0.35/100 ) # Ds*+ -> Ds+ pi0
    B_DstarS_DsX = B_DstarS_DsX_1 + B_DstarS_DsX_2

    # DD -> 3piX branching fractions
    DD_values = np.zeros((N_species,N_components))
    DD_errors = np.zeros((N_species,N_components))

    s = 0
    for i in range(N_species):

        species = all_species[i]
        for j in range(N_components):
            component = all_components[j]

            if(species == 100): # B+ -> D0 D0 K+
                DD_product = B_D0_3piX*B_D0_3piX
                if((component == 1) or (component == 2)): # D*D
                    DD_product *= B_Dstar0_D0X
                elif(component == 2):
                    DD_product *= B_Dstar0_D0X*B_Dstar0_D0X

            elif(species == 101): # B+ -> D+ D- K+
                if(component == 0): # D+D-
                    DD_product = B_Dp_3piX*B_Dp_3piX
                elif((component == 1) or (component == 2)): # D*+ D-
                    DD_product = (B_DstarP_D0pi*B_D0_3piX + B_DstarP_DpX*B_Dp_3piX)*B_Dp_3piX
                elif(component == 3): # D*+ D*-
                    DD_product = (B_DstarP_D0pi*B_D0_3piX + B_DstarP_DpX*B_Dp_3piX)*(B_DstarP_D0pi*B_D0_3piX + B_DstarP_DpX*B_Dp_3piX)

            elif(species == 102): # B+ -> Ds+ Ds- K+
                DD_product = B_Ds_3piX*B_Ds_3piX
                if((component == 1) or (component == 2)): # D*D
                    DD_product *= B_DstarS_DsX
                elif(component == 3):
                    DD_product *= B_DstarS_DsX*B_DstarS_DsX

            elif((species == 110) or (species == 130) or (species == 151)): # B0 -> D- D0 K+ or B+ -> D0 D+ K0 or  B+ -> D0 D+
                if(component == 0):
                    DD_product = B_Dp_3piX*B_D0_3piX
                elif(component == 1): # D*- D0
                    DD_product = (B_DstarP_D0pi*B_D0_3piX + B_DstarP_DpX*B_Dp_3piX)*B_D0_3piX
                elif(component == 2): # D- D*0
                    DD_product = B_Dp_3piX*B_Dstar0_D0X*B_D0_3piX
                elif(component == 3): # D*- D*+     
                    DD_product = (B_DstarP_D0pi*B_D0_3piX + B_DstarP_DpX*B_Dp_3piX)*B_Dstar0_D0X*B_D0_3piX

            elif((species == 120) or (species == 150)): # B0s -> Ds- D0 K+ or B+ -> D0 Ds+
                DD_product = B_Ds_3piX*B_D0_3piX
                if(component == 1): # Ds*- D0
                    DD_product *= B_DstarS_DsX
                elif(component == 2): # Ds- D*0
                    DD_product *= B_Dstar0_D0X
                elif(component == 2):
                    DD_product *= B_DstarS_DsX*B_Dstar0_D0X

            DD_values[s,j] = DD_product.nominal_value
            DD_errors[s,j] = DD_product.std_dev

        s += 1

    variable_branching_fraction = np.zeros((N_species,N_components))
    variable_branching_fraction_err = np.zeros((N_species,N_components))

    s = 0
    for i in range(N_species):
        species = all_species[i]
        for j in range(N_components):
            component = all_components[j] 

            B_B_DaDbX = ufloat( Bvalues[species].iloc[component], Berrors[species].iloc[component])
            B_DaDb_3piX = ufloat( DD_values[s,j], DD_errors[s,j])

            variable_branching_fraction[i,j] = unumpy.nominal_values(B_B_DaDbX*B_DaDb_3piX)
            variable_branching_fraction_err[i,j] = unumpy.std_devs(B_B_DaDbX*B_DaDb_3piX)

        s += 1
    
    variable_branching_fraction = unumpy.uarray( variable_branching_fraction, variable_branching_fraction_err )

    # fragmentation fraction ratio
    fu_fu = 1
    fu_fu_error = 0

    fd_fu = 0.998
    fd_fu_error = 0.063

    fs_fu = 0.1319 
    fs_fu_error = 0.0024

    fragmentation_fraction_term = np.zeros((N_species,N_components))
    fragmentation_fraction_term_err = np.zeros((N_species,N_components))

    for i in range(N_species):
        species = all_species[i]

        for j in range(N_components):
            if(species == 110):
                frag_ratio = ufloat(fd_fu, fd_fu_error)
            elif(species == 120):
                frag_ratio = ufloat(fs_fu, fs_fu_error)
            else:
                frag_ratio = ufloat(fu_fu, fu_fu_error)
            
            fragmentation_fraction_term[i,j] = unumpy.nominal_values(frag_ratio)
            fragmentation_fraction_term_err[i,j] = unumpy.std_devs(frag_ratio)

    fragmentation_fraction_term = unumpy.uarray( fragmentation_fraction_term, fragmentation_fraction_term_err )

    # Background efficiency (fixed)
    bkg_eff = np.zeros((N_species,N_components))
    bkg_eff_err = np.zeros((N_species,N_components))
    Nb_gen = np.zeros((N_species,N_components))

    for i in range(N_species):
        species = all_species[i]

        if((species == 100) or (species == 101) or (species == 102)):
            the_species = 100
        elif((species == 150) or (species == 151)):
            the_species = 150
        else:
            the_species = species

        # Acceptance, stripping and reco 1 efficiencies
        if(the_species == 100):
            # Acceptance
            eps_acc_2016 = ufloat(1.056/100, 0.002/100)
            eps_acc_2017 = ufloat(1.055/100, 0.002/100)
            eps_acc_2018 = ufloat(1.051/100, 0.002/100)
            eps_acc = (7.2/32.2)*eps_acc_2016 + (10.0/32.2)*eps_acc_2017 + (15.0/32.2)*eps_acc_2018

            # Stripping
            eps_strip_2016 = ufloat(1.719/100, 0.006/100)
            eps_strip_2017 = ufloat(1.760/100, 0.005/100)
            eps_strip_2018 = ufloat(1.763/100, 0.004/100)
            eps_strip = (7.2/32.2)*eps_strip_2016 + (10.0/32.2)*eps_strip_2017 + (15.0/32.2)*eps_strip_2018

            # Reco ineff
            eps_reco_2016 = ufloat(99.8/100, 0.2/100)
            eps_reco_2017 = ufloat(99.7/100, 0.1/100)
            eps_reco_2018 = ufloat(99.7/100, 0.1/100)
            eps_reco = (7.2/32.2)*eps_reco_2016 + (10.0/32.2)*eps_reco_2017 + (15.0/32.2)*eps_reco_2018

        elif(the_species == 110):
            # Acceptance
            eps_acc_2016 = ufloat(1.364/100, 0.002/100)
            eps_acc_2017 = ufloat(1.371/100, 0.003/100)
            eps_acc_2018 = ufloat(1.369/100, 0.003/100)
            eps_acc = (7.2/32.2)*eps_acc_2016 + (10.0/32.2)*eps_acc_2017 + (15.0/32.2)*eps_acc_2018

            # Stripping
            eps_strip_2016 = ufloat(1.591/100, 0.008/100)
            eps_strip_2017 = ufloat(1.643/100, 0.007/100)
            eps_strip_2018 = ufloat(1.625/100, 0.006/100)
            eps_strip = (7.2/32.2)*eps_strip_2016 + (10.0/32.2)*eps_strip_2017 + (15.0/32.2)*eps_strip_2018

            # Reco ineff
            eps_reco_2016 = ufloat(99.8/100, 0.2/100)
            eps_reco_2017 = ufloat(99.7/100, 0.2/100)
            eps_reco_2018 = ufloat(99.7/100, 0.2/100)
            eps_reco = (7.2/32.2)*eps_reco_2016 + (10.0/32.2)*eps_reco_2017 + (15.0/32.2)*eps_reco_2018

        elif(the_species == 120):
            # Acceptance
            eps_acc_2016 = ufloat(1.522/100, 0.003/100)
            eps_acc_2017 = ufloat(1.517/100, 0.003/100)
            eps_acc_2018 = ufloat(1.524/100, 0.003/100)
            eps_acc = (2.4/10.8)*eps_acc_2016 + (3.4/10.8)*eps_acc_2017 + (5.0/10.8)*eps_acc_2018

            # Stripping
            eps_strip_2016 = ufloat(2.75/100, 0.01/100)
            eps_strip_2017 = ufloat(2.81/100, 0.01/100)
            eps_strip_2018 = ufloat(2.799/100, 0.009/100)
            eps_strip = (2.4/10.8)*eps_strip_2016 + (3.4/10.8)*eps_strip_2017 + (5.0/10.8)*eps_strip_2018

            # Reco ineff
            eps_reco_2016 = ufloat(99.7/100, 0.2/100)
            eps_reco_2017 = ufloat(99.7/100, 0.1/100)
            eps_reco_2018 = ufloat(99.7/100, 0.1/100)
            eps_reco = (2.4/10.8)*eps_reco_2016 + (3.4/10.8)*eps_reco_2017 + (5.0/10.8)*eps_reco_2018

        elif(the_species == 130):
            # Acceptance
            eps_acc_2016 = ufloat(2.293/100, 0.004/100)
            eps_acc_2017 = ufloat(2.291/100, 0.004/100)
            eps_acc_2018 = ufloat(2.289/100, 0.004/100)
            eps_acc = (14.3/64.3)*eps_acc_2016 + (20.0/64.3)*eps_acc_2017 + (30.0/64.3)*eps_acc_2018

            # Stripping
            eps_strip_2016 = ufloat(0.785/100, 0.004/100)
            eps_strip_2017 = ufloat(0.806/100, 0.004/100)
            eps_strip_2018 = ufloat(0.805/100, 0.004/100)
            eps_strip = (14.3/64.3)*eps_strip_2016 + (20.0/64.3)*eps_strip_2017 + (30.0/64.3)*eps_strip_2018

            # Reco ineff
            eps_reco_2016 = ufloat(99.8/100, 0.2/100)
            eps_reco_2017 = ufloat(99.7/100, 0.2/100)
            eps_reco_2018 = ufloat(99.8/100, 0.2/100)
            eps_reco = (14.3/64.3)*eps_reco_2016 + (20.0/64.3)*eps_reco_2017 + (30.0/64.3)*eps_reco_2018

        elif(the_species == 150):
            # Acceptance
            eps_acc_2016 = ufloat(2.922/100, 0.005/100)
            eps_acc_2017 = ufloat(2.926/100, 0.005/100)
            eps_acc_2018 = ufloat(2.919/100, 0.005/100)
            eps_acc = (14.3/64.3)*eps_acc_2016 + (20.0/64.3)*eps_acc_2017 + (30.0/64.3)*eps_acc_2018

            # Stripping
            eps_strip_2016 = ufloat(1.102/100, 0.005/100)
            eps_strip_2017 = ufloat(1.140/100, 0.005/100)
            eps_strip_2018 = ufloat(1.127/100, 0.005/100)
            eps_strip = (14.3/64.3)*eps_strip_2016 + (20.0/64.3)*eps_strip_2017 + (30.0/64.3)*eps_strip_2018

            # Reco ineff
            eps_reco_2016 = ufloat(99.9/100, 0.2/100)
            eps_reco_2017 = ufloat(99.8/100, 0.2/100)
            eps_reco_2018 = ufloat(99.8/100, 0.2/100)
            eps_reco = (14.3/64.3)*eps_reco_2016 + (20.0/64.3)*eps_reco_2017 + (30.0/64.3)*eps_reco_2018

        if(the_species == 100):
            fc_den_2016 = ROOT.TFileCollection("fc_den_2016", "fc_den_2016", "Files_on_grid/MC_2016_BuDDKp_cocktail.txt")
            fc_den_2017 = ROOT.TFileCollection("fc_den_2017", "fc_den_2017", "Files_on_grid/MC_2017_BuDDKp_cocktail.txt")
            fc_den_2018 = ROOT.TFileCollection("fc_den_2018", "fc_den_2018", "Files_on_grid/MC_2018_BuDDKp_cocktail.txt")
        elif(the_species == 110):
            fc_den_2016 = ROOT.TFileCollection("fc_den_2016", "fc_den_2016", "Files_on_grid/MC_2016_BdDDKp_cocktail.txt")
            fc_den_2017 = ROOT.TFileCollection("fc_den_2017", "fc_den_2017", "Files_on_grid/MC_2017_BdDDKp_cocktail.txt")
            fc_den_2018 = ROOT.TFileCollection("fc_den_2018", "fc_den_2018", "Files_on_grid/MC_2018_BdDDKp_cocktail.txt")
        elif(the_species == 120):
            fc_den_2016 = ROOT.TFileCollection("fc_den_2016", "fc_den_2016", "Files_on_grid/MC_2016_BsDDKp_cocktail.txt")
            fc_den_2017 = ROOT.TFileCollection("fc_den_2017", "fc_den_2017", "Files_on_grid/MC_2017_BsDDKp_cocktail.txt")
            fc_den_2018 = ROOT.TFileCollection("fc_den_2018", "fc_den_2018", "Files_on_grid/MC_2018_BsDDKp_cocktail.txt")
        elif(the_species == 130):
            fc_den_2016 = ROOT.TFileCollection("fc_den_2016", "fc_den_2016", "Files_on_grid/MC_2016_BuDDK0_cocktail.txt")
            fc_den_2017 = ROOT.TFileCollection("fc_den_2017", "fc_den_2017", "Files_on_grid/MC_2017_BuDDK0_cocktail.txt")
            fc_den_2018 = ROOT.TFileCollection("fc_den_2018", "fc_den_2018", "Files_on_grid/MC_2018_BuDDK0_cocktail.txt")
        elif(the_species == 150):
            fc_den_2016 = ROOT.TFileCollection("fc_den_2016", "fc_den_2016", "Files_on_grid/MC_2016_BuDD_cocktail.txt")
            fc_den_2017 = ROOT.TFileCollection("fc_den_2017", "fc_den_2017", "Files_on_grid/MC_2017_BuDD_cocktail.txt")
            fc_den_2018 = ROOT.TFileCollection("fc_den_2018", "fc_den_2018", "Files_on_grid/MC_2018_BuDD_cocktail.txt")

        if(species == 100):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BuD0D0Kp/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BuD0starD0Kp/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BuD0D0starKp/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BuD0starD0starKp/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BuD0D0Kp/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BuD0starD0Kp/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BuD0D0starKp/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BuD0starD0starKp/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BuD0D0Kp/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BuD0starD0Kp/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BuD0D0starKp/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BuD0starD0starKp/MCDecayTree")
        elif(species == 101):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BuDpDmKp/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BuDpstarDmKp/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BuDpDmstarKp/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BuDpstarDmstarKp/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BuDpDmKp/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BuDpstarDmKp/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BuDpDmstarKp/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BuDpstarDmstarKp/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BuDpDmKp/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BuDpstarDmKp/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BuDpDmstarKp/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BuDpstarDmstarKp/MCDecayTree")
        elif(species == 102):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BuDsDsKp/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BuDsstarDsKp/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BuDsDsstarKp/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BuDsstarDsstarKp/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BuDsDsKp/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BuDsstarDsKp/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BuDsDsstarKp/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BuDsstarDsstarKp/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BuDsDsKp/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BuDsstarDsKp/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BuDsDsstarKp/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BuDsstarDsstarKp/MCDecayTree")
        elif(species == 110):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BdDmD0Kp/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BdDmstarD0Kp/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BdDmD0starKp/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BdDmstarD0starKp/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BdDmD0Kp/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BdDmstarD0Kp/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BdDmD0starKp/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BdDmstarD0starKp/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BdDmD0Kp/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BdDmstarD0Kp/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BdDmD0starKp/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BdDmstarD0starKp/MCDecayTree")
        elif(species == 120):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BsDsD0Kp/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BsDsstarD0Kp/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BsDsD0starKp/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BsDsstarD0starKp/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BsDsD0Kp/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BsDsstarD0Kp/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BsDsD0starKp/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BsDsstarD0starKp/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BsDsD0Kp/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BsDsstarD0Kp/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BsDsD0starKp/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BsDsstarD0starKp/MCDecayTree")
        elif(species == 130):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BuD0DpK0/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BuD0starDpK0/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BuD0DpstarK0/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BuD0starDpstarK0/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BuD0DpK0/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BuD0starDpK0/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BuD0DpstarK0/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BuD0starDpstarK0/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BuD0DpK0/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BuD0starDpK0/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BuD0DpstarK0/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BuD0starDpstarK0/MCDecayTree")
        elif(species == 150):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BuD0Ds/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BuD0starDs/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BuD0Dsstar/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BuD0starDsstar/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BuD0Ds/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BuD0starDs/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BuD0Dsstar/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BuD0starDsstar/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BuD0Ds/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BuD0starDs/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BuD0Dsstar/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BuD0starDsstar/MCDecayTree")
        elif(species == 151):
            t_gen_DD_2016 = ROOT.TChain("mc_ntuple_BuD0Dp/MCDecayTree")
            t_gen_DstarD_2016 = ROOT.TChain("mc_ntuple_BuD0starDp/MCDecayTree")
            t_gen_DDstar_2016 = ROOT.TChain("mc_ntuple_BuD0Dpstar/MCDecayTree")
            t_gen_DstarDstar_2016 = ROOT.TChain("mc_ntuple_BuD0starDpstar/MCDecayTree")

            t_gen_DD_2017 = ROOT.TChain("mc_ntuple_BuD0Dp/MCDecayTree")
            t_gen_DstarD_2017 = ROOT.TChain("mc_ntuple_BuD0starDp/MCDecayTree")
            t_gen_DDstar_2017 = ROOT.TChain("mc_ntuple_BuD0Dpstar/MCDecayTree")
            t_gen_DstarDstar_2017 = ROOT.TChain("mc_ntuple_BuD0starDpstar/MCDecayTree")

            t_gen_DD_2018 = ROOT.TChain("mc_ntuple_BuD0Dp/MCDecayTree")
            t_gen_DstarD_2018 = ROOT.TChain("mc_ntuple_BuD0starDp/MCDecayTree")
            t_gen_DDstar_2018 = ROOT.TChain("mc_ntuple_BuD0Dpstar/MCDecayTree")
            t_gen_DstarDstar_2018 = ROOT.TChain("mc_ntuple_BuD0starDpstar/MCDecayTree")

        t_gen_DD_2016.AddFileInfoList(fc_den_2016.GetList())
        t_gen_DstarD_2016.AddFileInfoList(fc_den_2016.GetList())
        t_gen_DDstar_2016.AddFileInfoList(fc_den_2016.GetList())
        t_gen_DstarDstar_2016.AddFileInfoList(fc_den_2016.GetList())

        t_gen_DD_2017.AddFileInfoList(fc_den_2017.GetList())
        t_gen_DstarD_2017.AddFileInfoList(fc_den_2017.GetList())
        t_gen_DDstar_2017.AddFileInfoList(fc_den_2017.GetList())
        t_gen_DstarDstar_2017.AddFileInfoList(fc_den_2017.GetList())

        t_gen_DD_2018.AddFileInfoList(fc_den_2018.GetList())
        t_gen_DstarD_2018.AddFileInfoList(fc_den_2018.GetList())
        t_gen_DDstar_2018.AddFileInfoList(fc_den_2018.GetList())
        t_gen_DstarDstar_2018.AddFileInfoList(fc_den_2018.GetList())

        for j in range(N_components):
            component = all_components[j]

            if(component == 0):
                N_bkg_den = t_gen_DD_2016.GetEntries() + t_gen_DD_2017.GetEntries() + t_gen_DD_2018.GetEntries()   
            elif(component == 1):
                N_bkg_den = t_gen_DstarD_2016.GetEntries() + t_gen_DstarD_2017.GetEntries() + t_gen_DstarD_2018.GetEntries() 
            elif(component == 2):
                N_bkg_den = t_gen_DDstar_2016.GetEntries() + t_gen_DDstar_2017.GetEntries() + t_gen_DDstar_2018.GetEntries() 
            elif(component == 3):
                N_bkg_den = t_gen_DstarDstar_2016.GetEntries() + t_gen_DstarDstar_2017.GetEntries() + t_gen_DstarDstar_2018.GetEntries() 

            bkg_eff[i,j] = unumpy.nominal_values(eps_acc*eps_strip*eps_reco)
            bkg_eff_err[i,j] = unumpy.std_devs(eps_acc*eps_strip*eps_reco)
            Nb_gen[i,j] = N_bkg_den

    eps_bkg_term = unumpy.uarray( bkg_eff, bkg_eff_err )
    
    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/Nb_gen.npy', Nb_gen)

    C_fix = np.zeros((N_species,N_components))
    C_fix_err = np.zeros((N_species,N_components))

    for i in range(N_species):
        for j in range(N_components):
            fixed_value = variable_branching_fraction[i,j]*fragmentation_fraction_term[i,j]*eps_bkg_term[i,j]

            C_fix[i,j] = unumpy.nominal_values(fixed_value)
            C_fix_err[i,j] = unumpy.std_devs(fixed_value)

    C_fix = unumpy.uarray( C_fix, C_fix_err )
    print("C_fix = ", C_fix)

    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/C_fix_values.npy', unumpy.nominal_values(C_fix))
    np.save('/panfs/felician/B2Ktautau/workflow/bdt_independent_inputs/C_fix_errors.npy', unumpy.std_devs(C_fix))

if __name__ == "__main__":
    main(sys.argv)