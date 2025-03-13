import sys
import numpy as np
import pandas as pd
import ROOT 

def eps_error(Num, Den):
    return (Num/Den)*np.sqrt( 1/Num + 1/Den )

def main(argv):

    ############################################################# Fixed terms, independing of background ###############################################################33
    # Branching fractions 
    B_D0Dsp = 9.0/1000
    B_D0Dps_error = 0.9/1000
    B_D0Kp = 3.947/100
    B_D0Kp_error = 0.030/100
    B_DsKKpi = 5.37/100
    B_DsKKpi_error = 0.10/100

    fixed_branching_fraction = B_D0Dsp*B_D0Kp*B_DsKKpi
    fixed_branching_fraction_error = np.sqrt( ((B_D0Kp*B_DsKKpi)**2)*(B_D0Dps_error**2) + ((B_D0Dsp*B_DsKKpi)**2)*(B_D0Kp_error**2) + ((B_D0Dsp*B_D0Kp)**2)*(B_DsKKpi_error**2) )
    # print(fixed_branching_fraction, " +/- ", fixed_branching_fraction_error)

    # Normalisation channel

    # Acceptance
    eps_norm_acc_2016 = 14.50/100
    eps_norm_acc_2017 = 14.48/100
    eps_norm_acc_2018 = 14.58/100
    eps_norm_acc_error = 0.10/100

    # Weighted average of the 3 years (based on how much MC was generated for the 3 years)
    eps_acc_norm = (10.8/34.2)*eps_norm_acc_2016 + (12./34.2)*eps_norm_acc_2017 + (11.4/34.2)*eps_norm_acc_2018
    eps_acc_norm_error = np.sqrt( ((10.8/34.2)**2)*(eps_norm_acc_error**2) + ((12./34.2)**2)*(eps_norm_acc_error**2) + ((11.4/34.2)**2)*(eps_norm_acc_error**2) )

    # Post-acc efficiency
    # Numerator of eps_norm: MC w/o TM cuts, after all selections and best candidate selection
    f_norm_num = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_72/post_sel_tree_bdt1_0_bdt2_0.root")
    t_norm_num = f_norm_num.Get("DecayTree")

    # Denominator of eps_norm: gen level MC
    fc_norm_den_2016 = ROOT.TFileCollection("fc_norm_den_2016", "fc_norm_den_2016", "Files_on_grid/MC_D0Dps_2016.txt")
    fc_norm_den_2017 = ROOT.TFileCollection("fc_norm_den_2017", "fc_norm_den_2017", "Files_on_grid/MC_D0Dps_2017.txt")
    fc_norm_den_2018 = ROOT.TFileCollection("fc_norm_den_2018", "fc_norm_den_2018", "Files_on_grid/MC_D0Dps_2018.txt")

    t_norm_den = ROOT.TChain("mc_ntuple/MCDecayTree")
    t_norm_den1 = ROOT.TChain("mc_ntuple/MCDecayTree")
    t_norm_den2 = ROOT.TChain("mc_ntuple/MCDecayTree")

    t_norm_den.AddFileInfoList(fc_norm_den_2016.GetList())
    t_norm_den1.AddFileInfoList(fc_norm_den_2017.GetList())
    t_norm_den2.AddFileInfoList(fc_norm_den_2018.GetList())

    N_norm_gen_2016 = t_norm_den.GetEntries()
    N_norm_gen_2017 = t_norm_den1.GetEntries()
    N_norm_gen_2018 = t_norm_den2.GetEntries()

    eps_norm_post_acc_num = t_norm_num.GetEntries()
    eps_norm_post_acc_den = N_norm_gen_2016+N_norm_gen_2017+N_norm_gen_2018
    
    eps_norm_post_acc = eps_norm_post_acc_num/eps_norm_post_acc_den
    eps_norm_post_acc_error = eps_error(eps_norm_post_acc_num, eps_norm_post_acc_den)

    eps_norm = eps_acc_norm*eps_norm_post_acc
    eps_norm_error = np.sqrt( (eps_acc_norm**2)*(eps_norm_post_acc_error**2) + (eps_acc_norm_error**2)*(eps_norm_post_acc**2) )
    # print(eps_norm, "+/-", eps_norm_error)

    # Yield
    N_norm = 32858
    N_norm_error = 195

    fixed_normalisation_channel = N_norm/eps_norm
    fixed_normalisation_channel_error = np.sqrt( ((1/eps_norm)**2)*(N_norm_error**2) + ((N_norm/(eps_norm**2))**2)*(eps_norm_error**2) )
    # print(fixed_normalisation_channel, " +/- ", fixed_normalisation_channel_error)

    fixed_independent_term = fixed_normalisation_channel/fixed_branching_fraction
    fixed_independent_term_error = np.sqrt( ((1/fixed_branching_fraction)**2)*(fixed_normalisation_channel_error**2) + (((fixed_normalisation_channel/(fixed_branching_fraction**2)))**2)*(fixed_branching_fraction_error**2) )
    ##########################################################################################################################################################

    ############################################################# Fixed terms, depending on background ###############################################################33
    all_species = [100, 101, 102, 110, 120, 130, 150, 151]
    all_components = [0, 1, 2, 3]

    N_species = len(all_species)
    N_components = len(all_components)

    # branching fractions
    Bvalues = pd.read_csv('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/B_values.csv')
    Berrors = pd.read_csv('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/B_errors.csv')
    DDvalues = pd.read_csv('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/DD_values.csv')
    DDerrors = pd.read_csv('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/DD_errors.csv')

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

    DDvalues = DDvalues.transpose()
    DDvalues.reset_index(drop=True, inplace=True)
    DDvalues = DDvalues.rename(columns=DDvalues.iloc[0]).drop(DDvalues.index[0])
    DDvalues.columns = DDvalues.columns.astype('int')
    DDvalues = DDvalues.reset_index(drop=True)
    DDvalues = DDvalues.dropna()

    DDerrors = DDerrors.transpose()
    DDerrors.reset_index(drop=True, inplace=True)
    DDerrors = DDerrors.rename(columns=DDerrors.iloc[0]).drop(DDerrors.index[0])
    DDerrors.columns = DDerrors.columns.astype('int')
    DDerrors = DDerrors.reset_index(drop=True)
    DDerrors = DDerrors.dropna()

    branching_fraction_term_1 = np.zeros((N_species,N_components))
    branching_fraction_term_1_error = np.zeros((N_species,N_components))
    for i in range(N_species):
        species = all_species[i]
        for j in range(N_components):
            component = all_components[j] 

            B_DaDbX = Bvalues[species].iloc[component]
            B_DaDbX_error = Berrors[species].iloc[component]
            B_DaDb_3piX = DDvalues[species].iloc[component]
            B_DaDb_3piX_error = DDerrors[species].iloc[component]

            branching_fraction_term_1[i,j] = B_DaDbX*B_DaDb_3piX
            branching_fraction_term_1_error[i,j] = np.sqrt( (B_DaDb_3piX**2)*(B_DaDbX_error**2) + (B_DaDbX**2)*(B_DaDb_3piX_error**2) )
    
    # print(branching_fraction_term_1)
    # print(branching_fraction_term_1_error)

    # fragmentation fraction ratio
    fu_fu = 1
    fu_fu_error = 0

    fd_fu = 0.998
    fd_fu_error = 0.063

    fs_fu = 0.1319 
    fs_fu_error = 0.0024

    fragmentation_fraction_term = np.zeros((N_species,N_components))
    fragmentation_fraction_term_error = np.zeros((N_species,N_components))

    for i in range(N_species):
        species = all_species[i]

        for j in range(N_components):
            if(species == 110):
                fragmentation_fraction_term[i,j] = fd_fu
                fragmentation_fraction_term_error[i,j] = fd_fu_error
            elif(species == 120):
                fragmentation_fraction_term[i,j] = fs_fu
                fragmentation_fraction_term_error[i,j] = fs_fu_error    
            else:
                fragmentation_fraction_term[i,j] = fu_fu
                fragmentation_fraction_term_error[i,j] = fu_fu_error

    # print(fragmentation_fraction_term)
    # print(fragmentation_fraction_term_error)

    # eps_bkg term
    eps_bkg_term = np.zeros((N_species,N_components))
    eps_bkg_term_error = np.zeros((N_species,N_components))

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
            eps_acc_2016 = 1.056/100
            eps_acc_2017 = 1.055/100
            eps_acc_2018 = 1.051/100

            eps_acc_2016_err = 0.002/100
            eps_acc_2017_err = 0.002/100
            eps_acc_2018_err = 0.002/100

            # 7.2 (2016), 10 (2017), 15 (2018) | tot = 32.2
            eps_acc = (7.2/32.2)*eps_acc_2016 + (10.0/32.2)*eps_acc_2017 + (15.0/32.2)*eps_acc_2018
            eps_acc_err = np.sqrt(  ((7.2/32.2)**2)*(eps_acc_2016_err**2) + ((10.0/32.2)**2)*(eps_acc_2017_err**2) + ((15.0/32.2)**2)*(eps_acc_2018_err**2) )

            # Stripping
            eps_strip_2016 = 1.719/100
            eps_strip_2017 = 1.760/100
            eps_strip_2018 = 1.763/100

            eps_strip_2016_err = 0.006/100
            eps_strip_2017_err = 0.005/100
            eps_strip_2018_err = 0.004/100

            eps_strip = (7.2/32.2)*eps_strip_2016 + (10.0/32.2)*eps_strip_2017 + (15.0/32.2)*eps_strip_2018
            eps_strip_err = np.sqrt(  ((7.2/32.2)**2)*(eps_strip_2016_err**2) + ((10.0/32.2)**2)*(eps_strip_2017_err**2) + ((15.0/32.2)**2)*(eps_strip_2018_err**2) )

            # Reco ineff
            eps_reco_2016 = 99.8/100
            eps_reco_2017 = 99.7/100
            eps_reco_2018 = 99.7/100

            eps_reco_2016_err = 0.2/100
            eps_reco_2017_err = 0.1/100
            eps_reco_2018_err = 0.1/100

            eps_reco = (7.2/32.2)*eps_reco_2016 + (10.0/32.2)*eps_reco_2017 + (15.0/32.2)*eps_reco_2018
            eps_reco_err = np.sqrt(  ((7.2/32.2)**2)*(eps_reco_2016_err**2) + ((10.0/32.2)**2)*(eps_reco_2017_err**2) + ((15.0/32.2)**2)*(eps_reco_2018_err**2) )

        elif(the_species == 110):
            # Acceptance
            eps_acc_2016 = 1.364/100
            eps_acc_2017 = 1.371/100
            eps_acc_2018 = 1.369/100

            eps_acc_2016_err = 0.002/100
            eps_acc_2017_err = 0.003/100
            eps_acc_2018_err = 0.003/100

            # 7.2 (2016), 10 (2017), 15 (2018) | tot = 32.2
            eps_acc = (7.2/32.2)*eps_acc_2016 + (10.0/32.2)*eps_acc_2017 + (15.0/32.2)*eps_acc_2018
            eps_acc_err = np.sqrt(  ((7.2/32.2)**2)*(eps_acc_2016_err**2) + ((10.0/32.2)**2)*(eps_acc_2017_err**2) + ((15.0/32.2)**2)*(eps_acc_2018_err**2) )

            # Stripping
            eps_strip_2016 = 1.591/100
            eps_strip_2017 = 1.643/100
            eps_strip_2018 = 1.625/100

            eps_strip_2016_err = 0.008/100
            eps_strip_2017_err = 0.007/100
            eps_strip_2018_err = 0.006/100

            eps_strip = (7.2/32.2)*eps_strip_2016 + (10.0/32.2)*eps_strip_2017 + (15.0/32.2)*eps_strip_2018
            eps_strip_err = np.sqrt(  ((7.2/32.2)**2)*(eps_strip_2016_err**2) + ((10.0/32.2)**2)*(eps_strip_2017_err**2) + ((15.0/32.2)**2)*(eps_strip_2018_err**2) )

            # Reco ineff
            eps_reco_2016 = 99.8/100
            eps_reco_2017 = 99.7/100
            eps_reco_2018 = 99.7/100

            eps_reco_2016_err = 0.2/100
            eps_reco_2017_err = 0.2/100
            eps_reco_2018_err = 0.2/100

            eps_reco = (7.2/32.2)*eps_reco_2016 + (10.0/32.2)*eps_reco_2017 + (15.0/32.2)*eps_reco_2018
            eps_reco_err = np.sqrt(  ((7.2/32.2)**2)*(eps_reco_2016_err**2) + ((10.0/32.2)**2)*(eps_reco_2017_err**2) + ((15.0/32.2)**2)*(eps_reco_2018_err**2) )

        elif(the_species == 120):
            # Acceptance
            eps_acc_2016 = 1.522/100
            eps_acc_2017 = 1.517/100
            eps_acc_2018 = 1.524/100

            eps_acc_2016_err = 0.003/100
            eps_acc_2017_err = 0.003/100
            eps_acc_2018_err = 0.003/100

            # 2.4 (2016), 3.4 (2017), 5 (2018) | tot = 10.8
            eps_acc = (2.4/10.8)*eps_acc_2016 + (3.4/10.8)*eps_acc_2017 + (5.0/10.8)*eps_acc_2018
            eps_acc_err = np.sqrt(  ((2.4/10.8)**2)*(eps_acc_2016_err**2) + ((3.4/10.8)**2)*(eps_acc_2017_err**2) + ((5.0/10.8)**2)*(eps_acc_2018_err**2) )

            # Stripping
            eps_strip_2016 = 2.75/100
            eps_strip_2017 = 2.81/100
            eps_strip_2018 = 2.799/100

            eps_strip_2016_err = 0.01/100
            eps_strip_2017_err = 0.01/100
            eps_strip_2018_err = 0.009/100

            eps_strip = (2.4/10.8)*eps_strip_2016 + (3.4/10.8)*eps_strip_2017 + (5.0/10.8)*eps_strip_2018
            eps_strip_err = np.sqrt(  ((2.4/10.8)**2)*(eps_strip_2016_err**2) + ((3.4/10.8)**2)*(eps_strip_2017_err**2) + ((5.0/10.8)**2)*(eps_strip_2018_err**2) )

            # Reco ineff
            eps_reco_2016 = 99.7/100
            eps_reco_2017 = 99.7/100
            eps_reco_2018 = 99.7/100

            eps_reco_2016_err = 0.2/100
            eps_reco_2017_err = 0.1/100
            eps_reco_2018_err = 0.1/100

            eps_reco = (2.4/10.8)*eps_reco_2016 + (3.4/10.8)*eps_reco_2017 + (5.0/10.8)*eps_reco_2018
            eps_reco_err = np.sqrt(  ((2.4/10.8)**2)*(eps_reco_2016_err**2) + ((3.4/10.8)**2)*(eps_reco_2017_err**2) + ((5.0/10.8)**2)*(eps_reco_2018_err**2) )

        elif(the_species == 130):
            # Acceptance
            eps_acc_2016 = 2.293/100
            eps_acc_2017 = 2.291/100
            eps_acc_2018 = 2.289/100

            eps_acc_2016_err = 0.004/100
            eps_acc_2017_err = 0.004/100
            eps_acc_2018_err = 0.004/100

            # 14.3 (2016), 20 (2017), 30 (2018) | tot = 64.3
            eps_acc = (14.3/64.3)*eps_acc_2016 + (20.0/64.3)*eps_acc_2017 + (30.0/64.3)*eps_acc_2018
            eps_acc_err = np.sqrt(  ((14.3/64.3)**2)*(eps_acc_2016_err**2) + ((20.0/64.3)**2)*(eps_acc_2017_err**2) + ((30.0/64.3)**2)*(eps_acc_2018_err**2) )

            # Stripping
            eps_strip_2016 = 0.785/100
            eps_strip_2017 = 0.806/100
            eps_strip_2018 = 0.805/100

            eps_strip_2016_err = 0.004/100
            eps_strip_2017_err = 0.004/100
            eps_strip_2018_err = 0.004/100

            eps_strip = (14.3/64.3)*eps_strip_2016 + (20.0/64.3)*eps_strip_2017 + (30.0/64.3)*eps_strip_2018
            eps_strip_err = np.sqrt(  ((14.3/64.3)**2)*(eps_strip_2016_err**2) + ((20.0/64.3)**2)*(eps_strip_2017_err**2) + ((30.0/64.3)**2)*(eps_strip_2018_err**2) )

            # Reco ineff
            eps_reco_2016 = 99.8/100
            eps_reco_2017 = 99.7/100
            eps_reco_2018 = 99.8/100

            eps_reco_2016_err = 0.2/100
            eps_reco_2017_err = 0.2/100
            eps_reco_2018_err = 0.2/100

            eps_reco = (14.3/64.3)*eps_reco_2016 + (20.0/64.3)*eps_reco_2017 + (30.0/64.3)*eps_reco_2018
            eps_reco_err = np.sqrt(  ((14.3/64.3)**2)*(eps_reco_2016_err**2) + ((20.0/64.3)**2)*(eps_reco_2017_err**2) + ((30.0/64.3)**2)*(eps_reco_2018_err**2) )

        elif(the_species == 150):
            # Acceptance
            eps_acc_2016 = 2.922/100
            eps_acc_2017 = 2.926/100
            eps_acc_2018 = 2.919/100

            eps_acc_2016_err = 0.005/100
            eps_acc_2017_err = 0.005/100
            eps_acc_2018_err = 0.005/100

            # 14.3 (2016), 20 (2017), 30 (2018) | tot = 64.3
            eps_acc = (14.3/64.3)*eps_acc_2016 + (20.0/64.3)*eps_acc_2017 + (30.0/64.3)*eps_acc_2018
            eps_acc_err = np.sqrt(  ((14.3/64.3)**2)*(eps_acc_2016_err**2) + ((20.0/64.3)**2)*(eps_acc_2017_err**2) + ((30.0/64.3)**2)*(eps_acc_2018_err**2) )

            # Stripping
            eps_strip_2016 = 1.102/100
            eps_strip_2017 = 1.140/100
            eps_strip_2018 = 1.127/100

            eps_strip_2016_err = 0.005/100
            eps_strip_2017_err = 0.005/100
            eps_strip_2018_err = 0.005/100

            eps_strip = (14.3/64.3)*eps_strip_2016 + (20.0/64.3)*eps_strip_2017 + (30.0/64.3)*eps_strip_2018
            eps_strip_err = np.sqrt(  ((14.3/64.3)**2)*(eps_strip_2016_err**2) + ((20.0/64.3)**2)*(eps_strip_2017_err**2) + ((30.0/64.3)**2)*(eps_strip_2018_err**2) )

            # Reco ineff
            eps_reco_2016 = 99.9/100
            eps_reco_2017 = 99.8/100
            eps_reco_2018 = 99.8/100

            eps_reco_2016_err = 0.2/100
            eps_reco_2017_err = 0.2/100
            eps_reco_2018_err = 0.2/100

            eps_reco = (14.3/64.3)*eps_reco_2016 + (20.0/64.3)*eps_reco_2017 + (30.0/64.3)*eps_reco_2018
            eps_reco_err = np.sqrt(  ((14.3/64.3)**2)*(eps_reco_2016_err**2) + ((20.0/64.3)**2)*(eps_reco_2017_err**2) + ((30.0/64.3)**2)*(eps_reco_2018_err**2) )


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

            eps_bkg_term[i,j] = (eps_acc*eps_strip*eps_reco)/N_bkg_den
            eps_bkg_term_error[i,j] = np.sqrt( ((eps_strip*eps_reco/N_bkg_den)**2)*(eps_acc_err**2) + ((eps_acc*eps_reco/N_bkg_den)**2)*(eps_strip_err**2) + ((eps_acc*eps_strip/N_bkg_den)**2)*(eps_reco_err**2) + (((eps_acc*eps_strip*eps_reco)/(N_bkg_den**2))**2)*(N_bkg_den_error**2) )

    # print(eps_bkg_term)
    # print(eps_bkg_term_error)

    fixed_dependent_term = np.zeros((N_species,N_components))
    fixed_dependent_term_error = np.zeros((N_species,N_components))

    for i in range(N_species):
        for j in range(N_components):
            fixed_dependent_term[i,j] = fixed_independent_term*branching_fraction_term_1[i,j]*fragmentation_fraction_term[i,j]*eps_bkg_term[i,j]
            fixed_dependent_term_error[i,j] = np.sqrt( ((branching_fraction_term_1[i,j]*fragmentation_fraction_term[i,j]*eps_bkg_term[i,j])**2)*(fixed_independent_term_error**2) + ((fixed_independent_term*fragmentation_fraction_term[i,j]*eps_bkg_term[i,j])**2)*(branching_fraction_term_1_error[i,j]**2) + ((fixed_independent_term*branching_fraction_term_1[i,j]*eps_bkg_term[i,j])**2)*(fragmentation_fraction_term_error[i,j]**2) + ((fixed_independent_term*branching_fraction_term_1[i,j]*fragmentation_fraction_term[i,j])**2)*(eps_bkg_term_error[i,j]**2)  )

    print(fixed_dependent_term)
    print(fixed_dependent_term_error)

    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_values_inputs.npy', fixed_dependent_term)
    np.save('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/yield_errors_inputs.npy', fixed_dependent_term_error)

if __name__ == "__main__":
    main(sys.argv)