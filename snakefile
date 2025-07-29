import numpy as np

localrules: exact_constraints, comparisons, fit_mass, make_sPlot_histos, compare_MC_sWeighted_data, addWeight, make_MLP_regression

# Species values
Ktautau_MC_species = [1, 10]
Ktautau_data_species = [2, 3]
Ktautau_species = Ktautau_MC_species+Ktautau_data_species
cocktail_MC_species = [100, 110, 120, 130, 150, 1000, 1100, 1200, 1300, 1500]
BuD0barDsp_MC_species = [7, 71, 72]
BuD0barDsp_data_species = [8, 81]
BuDDKp_species = [4, 5, 6]
BuD0D0Kp_species = [9, 0, -1]
all_species = Ktautau_species+cocktail_MC_species+BuD0barDsp_MC_species+BuD0barDsp_data_species+BuDDKp_species+BuD0D0Kp_species

# Years
years = [6, 7, 8]

# Lines
BuKtautau_MC_lines = {"6": 14, "7": 17, "8": 20}
BuKtautau_data_lines = {"6": 2368, "7": 2128, "8": 2669}
BuDDKp_cocktail_MC_lines = {"6": 10, "7": 10, "8": 12}
BdDDKp_cocktail_MC_lines = {"6": 8, "7": 11, "8": 8}
BsDDKp_cocktail_MC_lines = {"6": 7, "7": 5, "8": 9}
BuDDK0_cocktail_MC_lines = {"6": 8, "7": 8, "8": 7}
BuDD_cocktail_MC_lines = {"6": 7, "7": 7, "8": 7}
BuD0barDsp_MC_lines = {"6": 102, "7": 113, "8": 91}
BuD0barDsp_data_lines = {"6": 1880, "7": 1959, "8": 2484}

rule all:
    input:
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0, BDT=0.0),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0, BDT=0.5),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0, BDT=0.8),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0, BDT=0.95),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0, BDT=0.995),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0, BDT=0.9995),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.0),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.5),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.8),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.95),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.995),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.9995),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0, BDT=0.99915),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0.0000001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0.000001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0.00001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0.0001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0.001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="AllEvents", BF_sig=0.01, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0.0000001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0.000001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0.00001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0.0001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0.001, BDT=0.999),
        expand('/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf', fit_type="ToyData", fit_name="ErrorCategories", BF_sig=0.01, BDT=0.999)

        # expand('/panfs/felician/B2Ktautau/workflow/signal_toys_bf/{fit_type}_{fit_name}/BDT_{BDT}/BF_injected_vs_expected.pdf', fit_type="ToyData", fit_name="AllEvents", BDT=0.999),
        # expand('/panfs/felician/B2Ktautau/workflow/signal_toys_bf/{fit_type}_{fit_name}/BDT_{BDT}/BF_injected_vs_expected.pdf', fit_type="ToyData", fit_name="ErrorCategories", BDT=0.999)


# Every time grid files are updated need to:
# - run pidgen for ktautau / DDs over the pre-selection files (signal + normalisation channel)
# - seperate mc components in Ktautau
# - create pre-selection tree

rule pid_correction:
    ''' Runs PIDGen over MC (resampling) '''
    input:
        'PIDCalib/pid_resampling.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/PIDCalib/201{year}/Species_{species}/pid_corr_{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/PIDCalib/201{year}/Species_{species}/{line}.log'
    shell:
        'lb-conda pidgen/2024-05-08 python PIDCalib/pid_resampling.py {wildcards.year} {wildcards.species} {wildcards.line} &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/PIDCalib/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/PIDCalib/201{wildcards.year}/Species_{wildcards.species}/pid_corr.txt'

rule trigger_correction:
    ''' Runs trigger correction '''
    input:
        'TriggerCorr/L0Hadron_TOS_trigger_corr.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/TriggerCorrection/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/TriggerCorrection/201{year}/Species_{species}/{line}.log'
    shell:
        'root -l -b -q \'TriggerCorr/L0Hadron_TOS_trigger_corr.C( {wildcards.year}, {wildcards.species}, {wildcards.line} )\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/TriggerCorrection/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/TriggerCorrection/201{wildcards.year}/Species_{wildcards.species}/trig_corr.txt'

rule separate_reco_mc_components:
    ''' Creates a .root file with 1 branch (component) that allows to separate the 3pi3pi (component = 0), 3pi3pipi0 (component = 1) and 3pi3pi2pi0 (component = 2) components in the reco MC based on the true knowledge from gen MC '''
    input:
        'create_MC_component_branch.C',
        'Files_on_grid/MC_201{year}.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.log'
    shell:
        'root -l -b -q \'create_MC_component_branch.C( {wildcards.year}, {wildcards.line})\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/mc_components.txt'

rule separate_reco_cocktail_mc_components:
    ''' Creates a TTree with a component branch specifying wether the decay belongs to DD (0), DstarD (1), DDstar (2) or DstarDstar (3) '''
    input:
        'create_cocktail_MC_component_branch.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_{species}/{line}.log'
    shell:
        'root -l -b -q \'create_cocktail_MC_component_branch.C( {wildcards.year}, {wildcards.species}, {wildcards.line})\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{wildcards.year}/Species_{wildcards.species}/cocktail_mc_components.txt'

# rule abcd_method_plots:
#     ''' Creates plots for the ABCD method (before BDT selection) '''
#     input:
#         'abcd_method_plots.C'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/abcd_method_plots/2D_plot_MC.pdf'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/abcd_method_plots/out.log'
#     shell:
#         'root -l -b -q abcd_method_plots.C &> {log}'

# rule abcd_method_plots:
#     ''' Creates plots for the ABCD method (before BDT selection) '''
#     input:
#         'abcd_method_plots_post_selections.C'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/abcd_method_plots_post_selections/2D_plot_MC.pdf'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/abcd_method_plots_post_selections/out.log'
#     shell:
#         'root -l -b -q abcd_method_plots_post_selections.C &> {log}'

rule selections_efficiency_tables:
    ''' Saves all the selections of each species in a ROOT file. Creates the efficiency table (all years). Efficiencies are computed on TM MC.'''
    input:
        'selections_efficiency_tables.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/acceptance_table.tex'
        # '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/selections.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/out.log'
    shell:
        'python -u selections_efficiency_tables.py {wildcards.species} True &> {log}'

rule create_pre_selection_tree:
    ''' Creates a tree passing the pre-selections for 1 of the following species:
        1 - Ktautau MC (component = 0, 1, 2 for 3pi3pi, 3pi3pipi0 and 3pi3pi2pi0)
        2 - Ktautau RS data
        3 - Ktautau WS data 
        4 - DDK MC
        5 - DDK RS data
        6 - DDK WS data

        If the flag createTable is true, it computes the pre-selection efficiencies.
        CAREFUL: set to true only for 1 line / set to false by default.
    '''
    input:
        'create_pre_sel_tree.C',
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/mc_components.txt',
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_100/cocktail_mc_components.txt',
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_110/cocktail_mc_components.txt',
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_120/cocktail_mc_components.txt',
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_130/cocktail_mc_components.txt',
        '/panfs/felician/B2Ktautau/workflow/separate_reco_cocktail_mc_components/201{year}/Species_150/cocktail_mc_components.txt',
        '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/selections.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/{line}.root' 
    log:
        '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/{line}.log'
    resources:
        time = "99:00:00",
        mem_mb = 25000   
    shell:
        'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, {wildcards.species}, {wildcards.line})\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{wildcards.year}/Species_{wildcards.species}/pre_sel_tree.txt'

rule compare_vars:
    ''' Compares variables between 2 files (e.g. signal MC vs WS data) '''
    input:
        'compare_vars.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/compare_vars/2016/Ktautau_sig_vs_bkg_rect_cuts/var_0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/compare_vars/2016/Ktautau_sig_vs_bkg_rect_cuts/out.log'
    shell:
        'root -l -b -q compare_vars.C &> {log}'

rule make_MLP_regression:
    ''' TMVA regression for MLP initialisation '''
    input:
        'TMVARegression.C',
        RECO_files = '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_{species}/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/make_MLP_regression/TMVAReg_{species}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/make_MLP_regression/out_{species}.log'
    shell:
        'root -l -b -q \'TMVARegression.C( \"{input.RECO_files}\", {wildcards.species} ) \' &> {log};'

rule run_standalone_fitter:
    ''' Run offline minimisation '''
    # It needs lb-conda default/2024-06-08 (because of SymPy, pygsl work with new python version now, need to understand why sympy does not work)
    input:
        'decay_fit_gsl.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Species_{species}/{line}.log'
    resources:
        time = "99:00:00"
    shell:
        'python -u decay_fit_gsl.py {wildcards.year} {wildcards.species} {wildcards.line} &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/standalone_fitter/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/standalone_fitter/201{wildcards.year}/Species_{wildcards.species}/fit_results.txt'

rule exact_constraints:
    ''' Plots functions in system of equations '''
    input:
        'exact_constraints.C',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_2/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_2/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_2/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/exact_constraints/bmass_mc_components_2016.pdf',
        '/panfs/felician/B2Ktautau/workflow/exact_constraints/GSL_pulls/PVx.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/exact_constraints/out.log'
    shell:
        'root -l -b -q exact_constraints.C &> {log}'

rule comparisons:
    ''' Compares variables from 2 different files  '''
    input:
        'compare.C',
        FILES1 = 'Files_on_grid/MC_201{year}.txt',
        FILES2 = 'Files_on_grid/data_201{year}.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/comparisons/201{year}/TruthMatch_trigger_MC_vs_WS_data/var_0.gif'
    log: 
        '/panfs/felician/B2Ktautau/workflow/comparisons/201{year}/TruthMatch_trigger_MC_vs_WS_data/out.log'
    shell:
        'root -l -b -q \' compare.C( {wildcards.year}, \"{input.FILES1}\", \"{input.FILES2}\" ) \' &> {log}'

rule create_invariant_mass_tree:
    ''' Creates TTree with invariant mass combination branches for the mass vetoes '''  
    input:
        'create_invariant_mass_tree.C',
        '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201{year}/Species_{species}/{line}.log'
    shell:
       'root -l -b -q \'create_invariant_mass_tree.C( {wildcards.year}, {wildcards.species}, {wildcards.line})\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201{wildcards.year}/Species_{wildcards.species}/invariant_mass_tree.txt'

##############################################################################################################################################################################
##############################################################################   MC correction w/ normalisation channel   ####################################################################

rule create_dataset:
    ''' Creates the RooDataSet to make the mass fits for the normalisation channel '''
    ''' The RooDataSet is also used for sPlot / MC correction'''
    input:
        'createDataset.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_dataset/Species_{species}/mass_dataset.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_dataset/Species_{species}/out.log'
    shell:
        'root -l -b -q \'createDataset.C( {wildcards.species} )\' &> {log}'
    
rule fit_mass:
    ''' Makes mass fit with mass stored in RooDataSet data of normalisation channel after best candidate selection'''
    input:
        'fit_mass.C',
        '/panfs/felician/B2Ktautau/workflow/create_dataset/Species_{species}/mass_dataset.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/mass_fit_result.root',
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/mass_fit.pdf',
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/pulls_poisson.pdf',
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/params_poisson.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/out.log'
    shell:
        'root -l -b -q \'fit_mass.C( {wildcards.species} )\' &> {log}'
    
rule make_sPlot_histos:
    ''' Saves RooDataSet with sWeights. Species refers to data. '''
    input:
        'make_sPlot_histos.C',
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_8/mass_fit_result.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_8/splot_result.root',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_8/Sideband_plots/var_0.pdf',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_8/Splot_plots/var_0.pdf',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_8/Splot_vs_sideband_plots/var_0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_8/out.log'
    shell:
        'root -l -b -q \'make_sPlot_histos.C( 8 )\' &> {log}'

rule compare_MC_sWeighted_data:
    ''' Compares the MC with the sWeighted data (signal in MC vs signal in data). The species flag refers to the MC. False compares unweighted, True compares re-weighted. '''
    input:
        'compare_MC_sWeighted_data.C',
        '/panfs/felician/B2Ktautau/workflow/create_dataset/Species_72/mass_dataset.root',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/Species_8/splot_result.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/Species_72/MC_unweighted/var_0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/Species_72/MC_unweighted/out.log'
    shell:
        'root -l -b -q \'compare_MC_sWeighted_data.C( 72, false )\' &> {log}'

rule create_merged_root_from_dataset:
    ''' Creates a TTree from the merged datasets created with create_dataset '''
    input:
        'create_merged_root_from_dataset.C',
        '/panfs/felician/B2Ktautau/workflow/create_dataset/201{year}/Species_{species}/mass_dataset.root',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_8/splot_result.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/201{year}/Species_{species}/tree.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/201{year}/Species_{species}/out.log'
    shell:
        'root -l -b -q \'create_merged_root_from_dataset.C( {wildcards.year}, {wildcards.species} )\' &> {log}'

rule bdt_reweighter:
    ''' Saves re-weighter (splotted data vs MC) trained on DDs. The species flag refers to the MC. '''
    input:
        'bdt_reweighter.py',
        '/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/201{year}/Species_7/tree.root',
        '/panfs/felician/B2Ktautau/workflow/create_merged_root_from_dataset/201{year}/Species_8/tree.root',
        '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/tree_with_weights.root',
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/MC_unweighted/0.pdf',
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/MC_GBReweighter/0.pdf',
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/MC_unweighted_all/0.pdf',
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/MC_GBReweighter_all/0.pdf',
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/MC_DDs_reweighted/0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/out.log'
    shell:
        'python bdt_reweighter.py {wildcards.year} {wildcards.species} &> {log}'

rule addWeight:
    ''' Creates a TTree with data / MC weight that can be added as a friend. Species is MC. (Histogram method) '''
    input:
        'addWeight.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/addWeight/DDs_correction/201{year}/Species_{species}/tree_with_weight_1D.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/addWeight/DDs_correction/201{year}/Species_{species}/out.log'
    shell:
        'root -l -b -q \'addWeight.C( {wildcards.year}, {wildcards.species}, \"{input.MC_files}\", false )\' &> {log}'

######################################################### ML training ########################################################################################################
#########################################################    SkLearn    ######################################################################################################
rule sklearn_training:
    ''' Classifier to seprate signal from combinatorial background '''
    ''' species_name = (Ktautau,DDs) '''
    ''' step_name = (isolation,topology) '''
    ''' setup_name = (input_features,output_performance) '''
    input:
        'Sklearn/make_sklearn_training.py',
        # signal proxy - MC truth-matched
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_100/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_110/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_120/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_130/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_150/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_100/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_110/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_120/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_130/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_150/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_100/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_110/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_120/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_130/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_150/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_100/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_110/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_120/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_130/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_150/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_100/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_110/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_120/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_130/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_150/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_100/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_110/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_120/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_130/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_150/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_100/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_110/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_120/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_130/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_150/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_100/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_110/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_120/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_130/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_150/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_100/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_110/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_120/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_130/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_150/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_3/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_3/invariant_mass_tree.txt',
        # '/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_3/invariant_mass_tree.txt',
        '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_1/selections.root',
        '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_3/selections.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/out.log',
    log:
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/{species_name}/{setup_name}/{step_name}/out.log'
    shell:
        'python -u Sklearn/make_sklearn_training.py {wildcards.species_name} {wildcards.step_name} {wildcards.setup_name} False &> {log}'

rule sklearn_response:
    ''' Compute sklearn response '''
    input:
        'Sklearn/compute_sklearn_response.py',
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/clf_isolation.pkl',
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/clf_topology.pkl'
    output:
        '/panfs/felician/B2Ktautau/workflow/sklearn_{name}/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/sklearn_{name}/201{year}/Species_{species}/{line}.log'
    shell:
        'python -u Sklearn/compute_sklearn_response.py {wildcards.year} {wildcards.species} {wildcards.line} {wildcards.name} &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/sklearn_{wildcards.name}/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/sklearn_{wildcards.name}/201{wildcards.year}/Species_{wildcards.species}/bdt_output.txt'

rule sklearn_plots:
    ''' Plots to check if MVA output is correlated with B+ mass '''
    input:
        'Sklearn/sklearn_plots.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/sklearn_plots/iso_ktautau_mc_components.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/sklearn_plots/out.log'
    shell:
        'python -u Sklearn/sklearn_plots.py &> {log}'

##############################################################################################################################################################################

# rule get_run_numbers:
#     ''' Returns list with run numbers present in a particular .dst file of RS data (after pre-selections) '''
#     input:
#         'run_numbers_in_subset.C',
#         RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.txt'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.log'
#     shell:
#         'root -l -b -q \'run_numbers_in_subset.C( {wildcards.year}, \"{input.RS_DATA_files}\" )\' > {log}'

# rule get_lumi_of_subset:
#     ''' Computes luminosity of a .dst file of RS data (after pre-selections) '''
#     input:
#         'lumi_of_subset.py',
#         RUNS_files = '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.txt'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.log'
#     shell:
#         'python lumi_of_subset.py {wildcards.year} \"{input.RUNS_files}\" > {log};'

# rule sensitivity:
#     ''' Computes sensitivity for each BDT cut and each line of MC / RS data files'''
#     input:
#         'sensitivity.C',
#         MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_-1/pre_sel_tree.txt',
#         RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
#         MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_method1_type6_bkgProxy1/TMVA_response.txt',
#         RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method1_type6_bkgProxy1/TMVA_response.txt',
#         BKG_YIELD_file = '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.root'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.root'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.log'
#     shell:
#         'root -l -b -q \'sensitivity.C( {wildcards.year}, {wildcards.nprob}, {wildcards.bdt_it}, \"{input.MC_files}\", \"{input.RS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.BKG_YIELD_file}\" ) \' > {log};'
#         'ls /panfs/felician/B2Ktautau/workflow/sensitivity/201{wildcards.year}/N_prob_{wildcards.nprob}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/sensitivity/201{wildcards.year}/N_prob_{wildcards.nprob}/sensitivity.txt'

# rule sensitivity_plots:
#     ''' Estimates the 95% CL upper limit on B(B+ -> K+ tautau) '''
#     input:
#         'sensitivity_plots.C'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/bdt_efficiency.gif',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/bdt_efficiency.pdf',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/sensitivity.gif',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/sensitivity.pdf',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/signal_yield.gif',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/signal_yield.pdf',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/background_yield.gif',
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/background_yield.pdf'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/sensitivity.log'
#     shell:
#         'root -l -b -q \'sensitivity_plots.C( {wildcards.year}, {wildcards.nprob} ) \' > {log}'

#########################################################    Mass vetoes    ######################################################################################################

rule mass_vetoe_plots:
    ''' Creates plots for the mass vetoes '''   
    input:
        'mass_vetoes_plots.C'
        # '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_1/post_sel_tree_bdt1_0.6_bdt2_0.8.root',
        # '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt1_0.6_bdt2_0.8.root',
        # '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0.6_bdt2_0.8.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/2_particles/Bp_M01.pdf',
        '/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/3_particles/Bp_M012.pdf',
        '/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/4_particles/Bp_M0123.pdf',
        '/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/5_particles/Bp_M01234.pdf',
        '/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/6_particles/Bp_M012345.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/mass_vetoes_plots/out.log'
    shell:
        'root -l -b -q mass_vetoes_plots.C &> {log}'

#########################################################    Multiple events    ######################################################################################################

rule multiple_events:
    ''' Creates a TTree with a branch saying whether the candidate passes the best candidate selection or not '''
    ''' Before post-selection tree '''
    input:
        'multiple_events.py'
        # '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/multiple_events/201{year}/Species_{species}/{line}.root',
    log:
        '/panfs/felician/B2Ktautau/workflow/multiple_events/201{year}/Species_{species}/{line}.log'
    shell:
        'python -u multiple_events.py {wildcards.year} {wildcards.species} {wildcards.line} &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/multiple_events/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/multiple_events/201{wildcards.year}/Species_{wildcards.species}/multiple_events.txt'

rule best_candidate_selection_plots:
    ''' Makes plots after best candidate selection '''      
    input:
        'best_candidate_selection_plots.C',
        '/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_72/multiple_events.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/n_candidates_dds_mc.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/best_candidate_selection_plots/out.log'
    shell:
        'root -l -b -q best_candidate_selection_plots.C &> {log}'

#########################################################   Post selections Ktautau trees   ##########################################################################################
rule create_post_selection_tree:
    ''' Creates TTree passing GSL, mass vetoes and best candidate selection (if not TM) '''
    ''' All years and files are merged into a single output .root file '''
    input:
        'create_post_selection_tree.C',
        '/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/selections.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_{species}/post_sel_tree_bdt_{bdt}.root' 
        # '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_{species}/post_sel_tree_bdt1_{bdt1}_bdt2_{bdt2}.root' 
    log:
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_{species}/out_bdt_{bdt}.log'
        # '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_{species}/out_{bdt1}_{bdt2}.log'
    shell:  
        'root -l -b -q \'create_post_selection_tree.C( {wildcards.species}, {wildcards.bdt} )\' &> {log}' # to create the TTree

#########################################################   Fit inputs    ######################################################################################################

rule branching_fraction_inputs:
    ''' Saves the product of all fixed terms in the B->K tau tau branching fraction '''
    input:
        'Final_fit/branching_fraction_inputs.py',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_72/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_8/mass_fit_result.root',
        'Files_on_grid/MC_D0Dps_2016.txt',
        'Files_on_grid/MC_D0Dps_2017.txt',
        'Files_on_grid/MC_D0Dps_2018.txt',
        'Files_on_grid/MC_2016.txt',
        'Files_on_grid/MC_2017.txt',
        'Files_on_grid/MC_2018.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_s_den.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const_err.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy'
    log:
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/out.log'
    shell:
        'python -u Final_fit/branching_fraction_inputs.py &> {log}'


# rule rough_sensitivity_2D:
#     ''' Estimate 95% C.L. sensitivity using the normalisation channel (2D) '''
#     ''' the_case = "zero_bkg_3pi3pi", "zero_bkg_all_mc", "comb_bkg_3pi3pi", "comb_bkg_all_mc" '''
#     input: 
#         'rough_sensitivity_estimate_2D.py',
#         '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root',
#         # [ ['/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{0}_bdt2_{1}.npy'.format(BDT1,BDT2) for BDT1 in np.linspace(0,1,31)] for BDT2 in np.linspace(0,1,31) ]
#     output:
#         '/panfs/felician/B2Ktautau/workflow/rough_sensitivity_2D/BF_vs_bdt1_vs_bdt2_bkg_config_{bkg_config}.pdf'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/rough_sensitivity_2D/out_bkg_config_{bkg_config}.log'
#     shell:
#         'python -u rough_sensitivity_estimate_2D.py {wildcards.bkg_config} &> {log}'

rule yield_estimates_inputs:
    ''' Computes the inputs for the yield estimates that are fixed (need only to be computed once) '''
    input:
        'Final_fit/yield_estimates_inputs.py',
        'Files_on_grid/MC_2016_BuDDKp_cocktail.txt',
        'Files_on_grid/MC_2017_BuDDKp_cocktail.txt',
        'Files_on_grid/MC_2018_BuDDKp_cocktail.txt',
        'Files_on_grid/MC_2016_BdDDKp_cocktail.txt',
        'Files_on_grid/MC_2017_BdDDKp_cocktail.txt',
        'Files_on_grid/MC_2018_BdDDKp_cocktail.txt',
        'Files_on_grid/MC_2016_BsDDKp_cocktail.txt',
        'Files_on_grid/MC_2017_BsDDKp_cocktail.txt',
        'Files_on_grid/MC_2018_BsDDKp_cocktail.txt',
        'Files_on_grid/MC_2016_BuDDK0_cocktail.txt',
        'Files_on_grid/MC_2017_BuDDK0_cocktail.txt',
        'Files_on_grid/MC_2018_BuDDK0_cocktail.txt',
        'Files_on_grid/MC_2016_BuDD_cocktail.txt',
        'Files_on_grid/MC_2017_BuDD_cocktail.txt',
        'Files_on_grid/MC_2018_BuDD_cocktail.txt',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/B_values.csv',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/branching_fractions/B_errors.csv'
    output:
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_den.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_acc_strip_reco.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_acc_strip_reco_err.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/C_const.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/C_const_err.npy'
    log:
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/out.log'
    shell:
        'python Final_fit/yield_estimates_inputs.py &> {log}'

rule fit_inputs:
    ''' Estimate yields after all analysis selections using cocktail MC samples '''     
    input: 
        'Final_fit/fit_inputs.py',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/A_const_err.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/C_const.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/C_const_err.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_s_den.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_den.npy',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_acc_strip_reco.npy',
        '/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_acc_strip_reco_err.npy',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt_0.0.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/A.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/A_err.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/C.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/C_err.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/yield_values.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/yield_errors.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/eps_bkg_values.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/eps_bkg_errors.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/yields_table.tex',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/eps_b_table.tex'
    log:
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/out.log'
    shell:
        'python -u Final_fit/fit_inputs.py {wildcards.BDT} &> {log}'

# rule yield_estimates_plots:
#     ''' Makes plots normalised to the yield estimates '''
#     input:
#         'yield_estimates_plots.py',
#         '/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{BDT1}_bdt2_{BDT2}.npy'
#     output:
#         '/panfs/felician/B2Ktautau/workflow/yield_estimates_plots/yield_plots_bdt1_{BDT1}_bdt2_{BDT2}.pdf',
#         '/panfs/felician/B2Ktautau/workflow/yield_estimates_plots/yield_plots_bdt1_{BDT1}_bdt2_{BDT2}_comb_bkg.pdf'
#     log:
#         '/panfs/felician/B2Ktautau/workflow/yield_estimates_plots/out_bdt1_{BDT1}_bdt2_{BDT2}.log'
#     shell:
#         'python yield_estimates_plots.py {wildcards.BDT1} {wildcards.BDT2} &> {log}'


#########################################################   Final fit (binned, pyhf)  ######################################################################################################
rule fit_templates_choice:
    ''' Plots for choosing fit templates '''
    input:
        'Final_fit/template_shapes_choice.py',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt_0.0.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/templates_choice/BDDKp_shapes_bdt_0.0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/templates_choice/out.log'
    resources:
        mem_mb=64000
    shell:
        'python -u Final_fit/template_shapes_choice.py &> {log}'


def get_inputs(wildcards):
    inputs = [
        'Final_fit/generate_histograms.py',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt_0.0.root',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt_0.0.root',
        f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{wildcards.BDT}/C.npy',
        f'/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{wildcards.BDT}/C_err.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy',
    ]
    if wildcards.fit_type == "ToyData":
        inputs.append(
            f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{wildcards.BDT}/ncomb_value.npy'
        )
    return inputs

rule rs_data_pdf:
    ''' Extrapolates PDF of RS data sidebands to the signal region. Used to re-weight the WS data '''
    input:
        'Final_fit/rs_data_pdf.py',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/rs_data_pdf/BDT_{BDT}/fit_plot.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/rs_data_pdf/BDT_{BDT}/out.log'
    resources:
        mem_mb=64000
    shell:
        'python -u Final_fit/rs_data_pdf.py {wildcards.BDT} &> {log}'

rule generate_histograms:
    ''' Generates fit histograms after the BDT cuts and saves the expected normalisation factors 
        In this function the following parameters are defined:
            - error categories
            - signal region
            - relative error threshold (= 0): relative error below which we want to ignore the statistical error on a bin in the fit
    '''
    input:
        get_inputs
    output:
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/histograms.root',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/histogram_errors.root',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/N_comb.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/C.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/C_err.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/signal_region_indices.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/channel_eff_value.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/channel_eff_error.npy',
    log:
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/out.log'  
    resources:
        mem_mb=64000
    shell:
        'python -u Final_fit/generate_histograms.py {wildcards.fit_type} {wildcards.BDT} &> {log}'

rule pyhf_fit_rs_data_sidebands:
    input:
        'Final_fit/do_all.py',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/A.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/A_err.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/histograms.root',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/histogram_errors.root',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/N_comb.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/channel_eff_value.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/channel_eff_error.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/C.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/C_err.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/RSDataSidebands/BDT_{BDT}/signal_region_indices.npy'
    output:
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{BDT}/ncomb_value.npy',
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{BDT}/ncomb_error.npy'
    log:
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/RSDataSidebands_AllEvents/BF_sig_0/BDT_{BDT}/out.log'
    shell:
        'python -u Final_fit/do_all.py "AllEvents" "RSDataSidebands" 0 {wildcards.BDT} &> {log}'


rule pyhf_fit:
    ''' Does the fit, computes limit and does fit validation '''
    ''' fit_type = 
            "ToyDataSidebands": fit to toy data sidebands to validation signal region definition
            "ToyData": fit to toy data in the fit region, [4,8] GeV, for fit validation + expected upper limit calculation (blinded) 
            "RSDataSidebands": fit to RS data sidebands to extract estimate for N_comb
            "RSData": fit to RS data in the fit region, [4,8] GeV, for observed upper limit calculation (unblinded) 
        fit_name = 
            "AllEvents": fit to all events (1 channel)
            "ErrorCategories": fit in 3 mass error categories: [0,100] MeV (1), ]100,250] MeV (2) and > 250 MeV (3)
        BF_sig = any value between 0 and 1, this is the injected BF value (BF_sig = 0 for background only toy studies)
        BDT1 = physics BDT value, [0,1]
        BDT2 = combinatorial BDT value, [0,1]   
    '''
    input:
        'Final_fit/do_all.py',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/A.npy',
        '/panfs/felician/B2Ktautau/workflow/fit_inputs/BDT_{BDT}/A_err.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B.npy',
        '/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/B_err.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/histograms.root',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/histogram_errors.root',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/N_comb.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/channel_eff_value.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/channel_eff_error.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/C.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/C_err.npy',
        '/panfs/felician/B2Ktautau/workflow/generate_histograms/{fit_type}/BDT_{BDT}/signal_region_indices.npy'
    output:
        # Limit calculation outputs:
        # '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/cls_limit.npy',
        # '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot.pdf',
        # '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/cls_limit.pdf',
        # '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/cls_obs_limit.npy'
        # Fit validation outputs:
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_plot_seed_0.pdf',
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/pulls/BF_sig.pdf',
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_validation_plot.pdf',
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_poi_mean.npy',
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/fit_poi_mean_error.npy'
    log:
        '/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{BDT}/out.log'
    shell:
        'python -u Final_fit/do_all.py \"{wildcards.fit_name}\" \"{wildcards.fit_type}\" {wildcards.BF_sig} {wildcards.BDT} &> {log}'


rule upper_limit_optimisation:
    ''' Minimises the 90% C.L. upper limit on the branching fraction '''
    input:
        'Final_fit/upper_limit_optimisation.py',
        [f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{{fit_type}}_{{fit_name}}/BF_sig_{{BF_sig}}/BDT_{bdt}/cls_limit.npy' for bdt in np.round( np.linspace(0.99,1,200), 5)] 
    output:
        '/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BF_vs_bdt.pdf'
    log:    
        '/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{fit_name}/BF_sig_{BF_sig}/out.log'
    shell:
        'python -u Final_fit/upper_limit_optimisation.py {wildcards.fit_type} {wildcards.fit_name} {wildcards.BF_sig} &> {log}'

rule signal_toys_bf:
    ''' Scatter plot of BF fitted vs BF injected '''
    input:
        'Final_fit/signal_toys_bf.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/signal_toys_bf/{fit_type}_{fit_name}/BDT_{BDT}/BF_injected_vs_expected.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/signal_toys_bf/{fit_type}_{fit_name}/BDT_{BDT}/out.log'
    shell:
        'python -u Final_fit/signal_toys_bf.py {wildcards.fit_type} {wildcards.fit_name} {wildcards.BDT} &> {log}'



rule candidate_selection:
    ''' At each BDT cut counts the number of candidates per event and picks one randomly. Saves TTree with 2 branchnes: n_candidates and is_best_cand. '''
    input:
        'Final_fit/candidate_selection.py',
        '/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_{species}/post_sel_tree_bdt_0.0.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_{species}/candidate_bdt_{BDT}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_{species}/out_bdt_{BDT}.log'
    resources:
        mem_mb=64000
    shell:
        'python Final_fit/candidate_selection.py {wildcards.species} {wildcards.BDT} &> {log}'
#########################################################################################################################################################################################################