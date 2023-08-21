localrules: make_pre_sel_eff_tables, compare_TMVA, compute_bkg_yield, get_run_numbers, get_lumi_of_subset, sensitivity_plots, compare_DTF_fitter

rule all:
    input:
        expand('/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_{species}/{line}.root', year=8, species=2, line=range(1,14))

###################################################### Reconstruction ############################################################################
rule test_analytical_reconstruction:
    ''' Tests the analytical reconstruction calculations in MC (gen / reco) and in data '''
    input:
        'analytic_reconstruction.C',
        MC_files = 'Files_on_grid/full_MC_201{year}.txt',
        DATA_files = 'Files_on_grid/full_data_201{year}.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/analytical_reco/201{year}/ana_reco_vars_{type}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/analytical_reco/201{year}/ana_reco_{type}.log'
    shell:
        'root -l -b -q \' analytic_reconstruction.C( {wildcards.year}, \"{input.MC_files}\", \"{input.DATA_files}\", {wildcards.type} ) \' &> {log}'
###################################################################################################################################################

rule separate_reco_mc_components:
    ''' Creates a .root file with 1 branch (component) that allows to separate the 3pi3pi (component = 0), 3pi3pipi0 (component = 1) and 3pi3pi2pi0 (component = 2) components in the reco MC based on the true knowledge from gen MC '''
    input:
        'create_MC_component_branch.C',
        MC_files = 'Files_on_grid/full_MC_201{year}.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/MC_component_{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.log'
    shell:
        'root -l -b -q \'create_MC_component_branch.C( {wildcards.year}, \"{input.MC_files}\", {wildcards.line})\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/mc_components.txt'

rule create_MC_pre_selection_tree:
    ''' Creates a MC tree passing the pre-selections. Species must be set to 0. Component can be choosen: 0 (3pi3pi), 1 (3pi3pi pi0), 2 (3pi3pi 2pi0), -1 (all). '''
    input:
        'create_pre_sel_tree.C',
        MC_files = 'Files_on_grid/full_MC_201{year}.txt',
        MC_component_file = '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/MC_component_{line}.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_{component}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_{component}/{line}.log'
    shell:        
        'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, 0, {wildcards.component}, \"{input.MC_files}\", \"{input.MC_component_file}\", \" \", {wildcards.line}, false)\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{wildcards.year}/Component_{wildcards.component}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{wildcards.year}/Component_{wildcards.component}/pre_sel_tree.txt'

rule create_DATA_pre_selection_tree:
    ''' Creates a DATA tree passing the pre-selections. Species can be chosen: 1 (RS data) or 2 (WS data). '''
    input:
        'create_pre_sel_tree.C',
        DATA_files = 'Files_on_grid/full_data_201{year}.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_{species}/{line}.log'
    resources:
        time = "12:00:00"
    shell:
        'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, {wildcards.species}, 0, \" \", \" \", \"{input.DATA_files}\", {wildcards.line}, false)\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{wildcards.year}/Species_{wildcards.species}/pre_sel_tree.txt'

rule make_pre_sel_eff_tables:
    ''' Creates pre-selection efficiency tables for data and MC '''
    ''' If species==0, it produces the MC table. If species==1 it produces the data table. '''
    input:
        'create_pre_sel_tree.C',
        MC_files = 'Files_on_grid/full_MC_201{year}.txt',
        DATA_files = 'Files_on_grid/full_data_201{year}.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/pre_selection_efficiency_tables/201{year}/species{species}_pre_sel_efficiency.tex'
    log:
        '/panfs/felician/B2Ktautau/workflow/pre_selection_efficiency_tables/201{year}/species{species}_pre_sel_efficiency.log'
    shell:
        'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, {wildcards.species}, 0, \"{input.MC_files}\", \" ", \"{input.DATA_files}\", 0, true)\' &> {log}'

rule run_standalone_fitter:
    ''' Run offline minimisation '''
    input:
        'decay_fit.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_{component}/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Component_{component}/fit_result_{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Component_{component}/{line}.log'
    resources:
        time = "12:00:00",
        mem_mb=25000
    shell:
        'root -l -b -q \' decay_fit.C( {wildcards.year}, \"{input.MC_files}\", {wildcards.component}, {wildcards.line} ) \' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/standalone_fitter/201{wildcards.year}/Component_{wildcards.component}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/standalone_fitter/201{wildcards.year}/Component_{wildcards.component}/fit_results.txt'

rule compare_DTF_fitter:
    ''' Compares DTF and standalone fitter outputs '''
    input:
        'compare_fits.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_{component}/pre_sel_tree.txt',
        FIT_files = '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Component_{component}/fit_results.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/compare_DTF_fitter/201{year}/Component_{component}/Bmass.gif'
    log:
        '/panfs/felician/B2Ktautau/workflow/compare_DTF_fitter/201{year}/Component_{component}/out.log'
    shell:
        'root -l -b -q \' compare_fits.C( {wildcards.year}, {wildcards.component}, \"{input.MC_files}\", \"{input.FIT_files}\" ) \' &> {log}'

######################################################### ML training ########################################################################################################
#########################################################     TMVA    ######################################################################################################## 
rule make_TMVA_training:
    ''' Does the TMVA training. Uses only 1 data grid file for training. Uses all MC grid files for training. '''
    # tmva_method (TMVA method used in the training)
    # 0 -> BDT w/ gradient boost
    # 1 -> BDT
    # 2 -> artificial neural network
    # 3 -> deep neural network w/ GPU
    # tmva_type (input variables used in the training)
    # 0 -> tau geometry + kinematic
    # 1 -> 0 response cut (bdt1 > -0.5) + tau isolation
    # 2 -> 2 response + 
    # 6 -> Mix (throw sink at TMVA)
    # background_proxy
    # 0 -> use RS data, with M > 6.5 GeV
    # 1 -> use WS data
    input:
        'TMVA/TMVA_training.C',
        SIGNAL_MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_0/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_2/pre_sel_tree.txt',
        SIGNAL_MC_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_10/tmva_method0_type0_bkgProxy1.txt',
        RS_DATA_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method0_type0_bkgProxy1.txt',
        WS_DATA_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_method0_type0_bkgProxy1.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_output_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/log_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_training.C( {wildcards.year}, {wildcards.tmva_method}, {wildcards.tmva_type}, {wildcards.background_proxy}, \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.SIGNAL_MC_files}\", \"{input.SIGNAL_MC_BDT1_files}\", \"{input.RS_DATA_BDT1_files}\", \"{input.WS_DATA_BDT1_files}\" )\' &> {log};'
        'rm -rf /panfs/felician/B2Ktautau/workflow/TMVA_training/201{wildcards.year}/tmva_dataset_method{wildcards.tmva_method}_type{wildcards.tmva_type}_bkgProxy{wildcards.background_proxy};'
        'mv tmva_dataset_method{wildcards.tmva_method}_type{wildcards.tmva_type}_bkgProxy{wildcards.background_proxy} /panfs/felician/B2Ktautau/workflow/TMVA_training/201{wildcards.year}'

rule evaluate_TMVA_response:
    ''' Evaluates the BDT response.'''
    # species (in which files will the BDT response be evaluated)
    # 0 -> MC (all 3 components, component=-1)
    # 10 -> MC (3pi3pi)
    # 11 -> MC (3pi3pi pi0)
    # 12 -> MC (3pi3pi 2pi0) 
    # 1 -> RS data
    # 2 -> WS data
    # tmva_method (TMVA method used in the training)
    # 0 -> BDT w/ gradient boost
    # 1 -> BDT
    # 2 -> artificial neural network
    # 3 -> deep neural network w/ GPU
    # tmva_type (input variables used in the training)
    # 0 -> tau geometry + kinematic
    # 1 -> 0 response + tau isolation
    # 2 -> 2 response + 
    # 6 -> Mix (throw sink at TMVA)
    # background_proxy
    # 0 -> use RS data, with M > 6.5 GeV
    # 1 -> use WS data
    input:
        'TMVA/TMVA_response.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_-1/pre_sel_tree.txt',
        MC_0_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_0/pre_sel_tree.txt',
        MC_1_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_1/pre_sel_tree.txt',
        MC_2_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_2/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_2/pre_sel_tree.txt',
        weights = '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_dataset_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/weights/TMVAClassification_BDTG.weights.xml',
        MC_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_method0_type0_bkgProxy1.txt',
        MC_0_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_10/tmva_method0_type0_bkgProxy1.txt',
        MC_1_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_11/tmva_method0_type0_bkgProxy1.txt',
        MC_2_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_12/tmva_method0_type0_bkgProxy1.txt',
        RS_DATA_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method0_type0_bkgProxy1.txt',
        WS_DATA_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_method0_type0_bkgProxy1.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_{species}/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_{species}/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_{line}.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_response.C( {wildcards.year}, {wildcards.species}, {wildcards.tmva_method}, {wildcards.tmva_type}, {wildcards.background_proxy}, \"{input.MC_files}\", \"{input.MC_0_files}\", \"{input.MC_1_files}\", \"{input.MC_2_files}\", \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.weights}\", \"{input.MC_BDT1_files}\", \"{input.MC_0_BDT1_files}\", \"{input.MC_1_BDT1_files}\", \"{input.MC_2_BDT1_files}\", \"{input.RS_DATA_BDT1_files}\", \"{input.WS_DATA_BDT1_files}\", {wildcards.line} )\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/TMVA_response/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/TMVA_response/201{wildcards.year}/Species_{wildcards.species}/tmva_method{wildcards.tmva_method}_type{wildcards.tmva_type}_bkgProxy{wildcards.background_proxy}.txt;'

rule compare_TMVA:
    ''' Compares BDT distributions evaluated in MC, RS and WS data '''
    input:
        'TMVA/TMVA_compare.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_-1/pre_sel_tree.txt',
        MC_0_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_0/pre_sel_tree.txt',
        MC_1_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_1/pre_sel_tree.txt',
        MC_2_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_2/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_2/pre_sel_tree.txt',
        MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_0_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_10/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_1_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_11/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_2_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_12/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.txt',
        WS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_method0_type0_bkgProxy1.txt',
        MC_0_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_10/tmva_method0_type0_bkgProxy1.txt',
        MC_1_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_11/tmva_method0_type0_bkgProxy1.txt',
        MC_2_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_12/tmva_method0_type0_bkgProxy1.txt',
        RS_DATA_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method0_type0_bkgProxy1.txt',
        WS_DATA_BDT1_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_method0_type0_bkgProxy1.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/sig_bkg_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.gif',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/sig_bkg_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.pdf',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/mc_components_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.gif',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/mc_components_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.pdf',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/RS_WS_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.gif',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/RS_WS_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/compare_tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_compare.C( {wildcards.year}, {wildcards.tmva_method}, {wildcards.tmva_type}, {wildcards.background_proxy}, {wildcards.num_data_dst_files}, \"{input.MC_files}\", \"{input.MC_0_files}\", \"{input.MC_1_files}\", \"{input.MC_2_files}\", \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.MC_0_BDT_files}\", \"{input.MC_1_BDT_files}\", \"{input.MC_2_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.WS_DATA_BDT_files}\", \"{input.MC_BDT1_files}\", \"{input.MC_0_BDT1_files}\", \"{input.MC_1_BDT1_files}\", \"{input.MC_2_BDT1_files}\", \"{input.RS_DATA_BDT1_files}\", \"{input.WS_DATA_BDT1_files}\" )\' > {log}'

#########################################################    LUMIN    ######################################################################################################## 
rule convert_tree_to_pandas:
    ''' Converts pre-selection tree into a pandas dataframe that contains the variables that will be used in the training and the B+ DTF mass for defining the signal region '''
    # make_pandas_df(year, species, lumin_type)
    # species = 10 (signal 3pi3pi MC), 11 (3pi3pi pi0 MC), 12 (3pi3pi 2pi0 MC),  0 (MC), 1 (RS data), 2 (WS data)
    input:
        'lumin/make_pandas_df.py',
        MC_0_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_0/pre_sel_tree.txt', # signal: 3pi3pi MC
        MC_1_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_1/pre_sel_tree.txt', # signal: 3pi3pi pi0 MC
        MC_2_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_2/pre_sel_tree.txt', # signal: 3pi3pi 2pi0 MC
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_-1/pre_sel_tree.txt', # signal: 3pi3pi MC,
        RS_data_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt', # background: RS data
        WS_data_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_2/pre_sel_tree.txt' # background: WS data
    output:
        '/panfs/felician/B2Ktautau/workflow/LUMIN_pandas/201{year}/species{species}_lumin_type{lumin_type}.csv'
    log:
        '/panfs/felician/B2Ktautau/workflow/LUMIN_pandas/201{year}/species{species}_lumin_type{lumin_type}.log'
    shell:
        'python lumin/make_pandas_df.py {wildcards.year} {wildcards.species} {wildcards.lumin_type} \"{input.MC_0_files}\" \"{input.MC_1_files}\" \"{input.MC_2_files}\" \"{input.MC_files}\" \"{input.RS_data_files}\" \"{input.WS_data_files}\" > {log}'

# rule make_LUMIN_training:
#     ''' Does the lumin training '''
#     input:
#         'lumin/B2Ktautau_classification.py',
#         '/panfs/felician/B2Ktautau/workflow/LUMIN_pandas/201{year}/species{species}_lumin_type{lumin_type}.csv'
#     output:

##############################################################################################################################################################################

rule compute_bkg_yield:
    ''' Computes the background yield in the signal region '''
    input:
        'compute_B_yield.C',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_2/pre_sel_tree.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_type6_bkgProxy1.txt',
        WS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_type6_bkgProxy1.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.log'
    shell: 
        'root -l -b -q \'compute_B_yield.C( {wildcards.year}, \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.WS_DATA_BDT_files}\" )\' > {log}'

rule get_run_numbers:
    ''' Returns list with run numbers present in a particular .dst file of RS data (after pre-selections) '''
    input:
        'run_numbers_in_subset.C',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.txt'
    log:
        '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.log'
    shell:
        'root -l -b -q \'run_numbers_in_subset.C( {wildcards.year}, \"{input.RS_DATA_files}\" )\' > {log}'

rule get_lumi_of_subset:
    ''' Computes luminosity of a .dst file of RS data (after pre-selections) '''
    input:
        'lumi_of_subset.py',
        RUNS_files = '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.log'
    shell:
        'python lumi_of_subset.py {wildcards.year} \"{input.RUNS_files}\" > {log};'

rule sensitivity:
    ''' Computes sensitivity for each BDT cut and each line of MC / RS data files'''
    input:
        'sensitivity.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_MC_pre_selection_tree/201{year}/Component_-1/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
        MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_type6_bkgProxy1.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_type6_bkgProxy1.txt',
        BKG_YIELD_file = '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.root',
        LUMI_file = '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root',
        MC_PREBDT_EFF_file = '/panfs/felician/B2Ktautau/workflow/pre_selection_efficiency_tables/201{year}/species0_pre_sel_efficiency.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.log'
    shell:
        'root -l -b -q \'sensitivity.C( {wildcards.year}, {wildcards.nprob}, {wildcards.bdt_it}, \"{input.MC_files}\", \"{input.RS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.BKG_YIELD_file}\", \"{input.LUMI_file}\", \"{input.MC_PREBDT_EFF_file}\" ) \' > {log};'
        'ls /panfs/felician/B2Ktautau/workflow/sensitivity/201{wildcards.year}/N_prob_{wildcards.nprob}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/sensitivity/201{wildcards.year}/N_prob_{wildcards.nprob}/sensitivity.txt'

rule sensitivity_plots:
    ''' Estimates the 95% CL upper limit on B(B+ -> K+ tautau) '''
    input:
        'sensitivity_plots.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/bdt_efficiency.gif',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/bdt_efficiency.pdf',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/sensitivity.gif',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/sensitivity.pdf',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/signal_yield.gif',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/signal_yield.pdf',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/background_yield.gif',
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/background_yield.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/sensitivity_plots/201{year}/N_prob_{nprob}/sensitivity.log'
    shell:
        'root -l -b -q \'sensitivity_plots.C( {wildcards.year}, {wildcards.nprob} ) \' > {log}'