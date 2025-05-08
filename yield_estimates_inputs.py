import sys
import numpy as np
import pandas as pd
import ROOT 
from uncertainties import ufloat
from uncertainties import unumpy

def eps_error(Num, Den):
    return (Num/Den)*np.sqrt( 1/Num + 1/Den )

def main(argv):

    ############################################################# Fixed terms, independing of background ###############################################################33
    # Branching fractions 
    B_D0Dsp = ufloat(9.0/1000, 0.9/1000)
    B_D0Kp = ufloat(3.947/100, 0.030/100)
    B_DsKKpi = ufloat(5.37/100,  0.10/100)
    fixed_branching_fraction = B_D0Dsp*B_D0Kp*B_DsKKpi
    print("Fixed BF term = ", fixed_branching_fraction)

    # Normalisation channel
    eps_norm_value = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_norm_value.npy')
    eps_norm_error = np.load('/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_norm_error.npy')
    eps_norm = ufloat(eps_norm_value, eps_norm_error)

    # Yield
    N_norm = ufloat(33187, 259)

    fixed_normalisation_channel = N_norm/eps_norm
    print("Fixed normalisation channel term = ", fixed_normalisation_channel)
    ##########################################################################################################################################################

    ############################################################# Fixed terms, depending on background ###############################################################33
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]
    all_components = [0, 1, 2, 3]

    N_species = len(all_species)
    N_components = len(all_components)

    # branching fractions
    # B(B -> Da Db X)
    Bvalues = pd.read_csv('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/B_values.csv')
    Berrors = pd.read_csv('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/B_errors.csv')

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
    print("Variable branching fraction term = ", variable_branching_fraction)

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
    print("Fragmentation fraction term = ", fragmentation_fraction_term)

    # eps_bkg term
    bkg_eff = np.zeros((N_species,N_components))
    bkg_eff_err = np.zeros((N_species,N_components))

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
            N_bkg_den_error = np.sqrt(N_bkg_den)
            N_bkg_den = ufloat(N_bkg_den, N_bkg_den_error)

            bkg_eff[i,j] = unumpy.nominal_values((eps_acc*eps_strip*eps_reco)/N_bkg_den)
            bkg_eff_err[i,j] = unumpy.std_devs((eps_acc*eps_strip*eps_reco)/N_bkg_den)

    eps_bkg_term = unumpy.uarray( bkg_eff, bkg_eff_err )
    print("Background efficiency term = ", eps_bkg_term)

    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_fixed_value.npy', unumpy.nominal_values(eps_bkg_term))
    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_fixed_error.npy', unumpy.std_devs(eps_bkg_term))

    fixed_term = np.zeros((N_species,N_components))
    fixed_term_err = np.zeros((N_species,N_components))

    for i in range(N_species):
        for j in range(N_components):
            fixed_value = (variable_branching_fraction[i,j]/fixed_branching_fraction)*fixed_normalisation_channel*fragmentation_fraction_term[i,j]*eps_bkg_term[i,j]

            fixed_term[i,j] = unumpy.nominal_values(fixed_value)
            fixed_term_err[i,j] = unumpy.std_devs(fixed_value)

    fixed_term = unumpy.uarray( fixed_term, fixed_term_err )
    print(fixed_term)

    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_values_inputs.npy', unumpy.nominal_values(fixed_term))
    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_errors_inputs.npy', unumpy.std_devs(fixed_term))

if __name__ == "__main__":
    main(sys.argv)