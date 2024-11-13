localrules: make_TMVA_training, compare_TMVA, compute_bkg_yield, get_run_numbers, get_lumi_of_subset, sensitivity_plots, exact_constraints, comparisons, compare_vars, create_dataset, fit_mass, make_sPlot_histos, compare_MC_sWeighted_data, addWeight, make_MLP_regression

rule all:
    input:
        expand('/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/{line}.root', year=7, species=8, line=range(1,1960))
        # expand('/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.root', year=7, line=8)
        # expand('/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.root', year=6, line=range(1,14)),
        # expand('/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.root', year=7, line=range(1,17)),
        # expand('/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.root', year=8, line=range(1,20))

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
        'lb-conda pidgen/2024-05-08 python PIDCalib/pid_resampling.py {wildcards.year} {wildcards.species} {wildcards.line} 2>&1 | tee {log};'
        'ls /panfs/felician/B2Ktautau/workflow/PIDCalib/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/PIDCalib/201{wildcards.year}/Species_{wildcards.species}/pid_corr.txt'

rule separate_reco_mc_components:
    ''' Creates a .root file with 1 branch (component) that allows to separate the 3pi3pi (component = 0), 3pi3pipi0 (component = 1) and 3pi3pi2pi0 (component = 2) components in the reco MC based on the true knowledge from gen MC '''
    input:
        'create_MC_component_branch.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/{line}.log'
    shell:
        'root -l -b -q \'create_MC_component_branch.C( {wildcards.year}, {wildcards.line})\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/mc_components.txt'

rule create_pre_selection_tree:
    ''' Creates a tree passing the pre-selections for 1 of the following species:
        10 - Ktautau MC (3pi3pi component)
        11 - Ktautau MC (3pi3pi pi0 component)
        12 - Ktautau MC (3pi3pi 2pi0 component)
        1 - Ktautau MC (all components)
        2 - Ktautau RS data
        3 - Ktautau WS data 
        4 - DDK MC
        5 - DDK RS data
        6 - DDK WS data

        If the flag createTable is true, it computes the pre-selection efficiencies.
        CAREFUL: set to true only for 1 line / set to false by default.
    '''
    input:
        'create_pre_sel_tree.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/{line}.log'
    resources:
        time = "99:00:00",
        mem_mb = 25000   
    shell:        
        'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, {wildcards.species}, {wildcards.line}, false)\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{wildcards.year}/Species_{wildcards.species}/pre_sel_tree.txt'

rule compare_vars:
    ''' Compares variables between 2 files (e.g. signal MC vs WS data) '''
    input:
        'compare_vars.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/compare_vars/201{year}/Ktautau_bkg_proxy/var_0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/compare_vars/201{year}/Ktautau_bkg_proxy/out.log'
    shell:
        'root -l -b -q \' compare_vars.C( {wildcards.year} ) \' &> {log}'

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
    input:
        'decay_fit_gsl.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Species_{species}/{line}.root'
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Species_{species}/Line_{line}/{i_first}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Species_{species}/{line}.log'
        # '/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{year}/Species_{species}/Line_{line}/{i_first}.log'
    resources:
        time = "99:00:00"
    shell:
        # 'root -l -b -q \'DECAY_FIT.C( {wildcards.year}, {wildcards.species}, {wildcards.line} ) \' &> {log};'
        'python -u decay_fit_gsl.py {wildcards.year} {wildcards.species} {wildcards.line} 2>&1 | tee {log};'
        'ls /panfs/felician/B2Ktautau/workflow/standalone_fitter/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/standalone_fitter/201{wildcards.year}/Species_{wildcards.species}/fit_results.txt'

rule exact_constraints:
    ''' Plots functions in system of equations '''
    input:
        'exact_constraints.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/exact_constraints/Species_{species}/F0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/exact_constraints/Species_{species}/out.log'
    shell:
        'root -l -b -q \' exact_constraints.C( {wildcards.species} ) \' &> {log}'

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

##############################################################################################################################################################################
############################################################################## MC correction ####################################################################

rule create_dataset:
    ''' Saves the mass variable and a colection of variables that are used in sideband subtraction and splot plots + TRatio plots in a RooDataSet called data '''
    input:
        'createDataset.C',
        FILES = '/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{year}/Species_{species}/pre_sel_tree.txt' 
    output:
        '/panfs/felician/B2Ktautau/workflow/create_dataset/201{year}/Species_{species}/mass_dataset.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/create_dataset/201{year}/Species_{species}/out.log'
    shell:
        'root -l -b -q \'createDataset.C( {wildcards.year}, {wildcards.species}, \"{input.FILES}\" )\' &> {log}'
    
rule fit_mass:
    ''' Makes mass fit with mass stored in RooDataSet data '''
    input:
        'fit_mass.C',
        '/panfs/felician/B2Ktautau/workflow/create_dataset/201{year}/Species_{species}/mass_dataset.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/fit_mass/201{year}/Species_{species}/mass_fit.pdf',
        '/panfs/felician/B2Ktautau/workflow/fit_mass/201{year}/Species_{species}/mass_fit_result.root'
        # '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/mass_fit.gif',
        # '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/mass_fit_result.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/fit_mass/201{year}/Species_{species}/out.log'
        # '/panfs/felician/B2Ktautau/workflow/fit_mass/Species_{species}/out.log'
    shell:
        'root -l -b -q \'fit_mass.C( {wildcards.year}, {wildcards.species} )\' &> {log}'
        # 'root -l -b -q \'fit_mass.C( 0, {wildcards.species} )\' &> {log}'
    
rule make_sPlot_histos:
    ''' Saves RooDataSet with sWeights. Species refers to data. '''
    input:
        'make_sPlot_histos.C',
        FIT_workspace = '/panfs/felician/B2Ktautau/workflow/fit_mass/201{year}/Species_{species}/mass_fit_result.root',
    output:
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_{species}/Sideband_plots/var_0.gif',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_{species}/splot_result.root',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_{species}/Splot_plots/Mass_var_plots/var_0.gif',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_{species}/Splot_plots/Sig_vs_bkg_plots/var_0.gif',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_{species}/Splot_vs_sideband_plots/var_0.gif'
    log:
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_{species}/out.log'
    shell:
        'root -l -b -q \'make_sPlot_histos.C( {wildcards.year}, {wildcards.species}, \"{input.FIT_workspace}\" )\' &> {log}'

rule compare_MC_sWeighted_data:
    ''' Compares the MC with the sWeighted data (signal in MC vs signal in data). The species flag refers to the MC. False compares unweighted, True compares re-weighted. '''
    input:
        'compare_MC_sWeighted_data.C',
        '/panfs/felician/B2Ktautau/workflow/create_dataset/201{year}/Species_{species}/mass_dataset.root',
        '/panfs/felician/B2Ktautau/workflow/make_sPlot_histos/201{year}/Species_8/splot_result.root',
        '/panfs/felician/B2Ktautau/workflow/bdt_reweighter/201{year}/DDs_correction/Species_{species}/tree_with_weights.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201{year}/DDs_correction/Species_{species}/MC_reweighted/var_0.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201{year}/DDs_correction/Species_{species}/MC_reweighted/out.log'
    shell:
        'root -l -b -q \'compare_MC_sWeighted_data.C( {wildcards.year}, {wildcards.species}, true )\' &> {log}'

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
        # WEIGHT_file  = '/panfs/felician/B2Ktautau/workflow/compare_MC_sWeighted_data/201{year}/DDs_correction/Species_{species}/weights.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/addWeight/DDs_correction/201{year}/Species_{species}/tree_with_weight.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/addWeight/DDs_correction/201{year}/Species_{species}/out.log'
    shell:
        'root -l -b -q \'addWeight.C( {wildcards.year}, {wildcards.species}, \"{input.MC_files}\", false )\' &> {log}'

######################################################### ML training ########################################################################################################
#########################################################     TMVA    ######################################################################################################## 

rule make_TMVA_training:
    ''' Does the TMVA training. '''
    # Uses TMVA to separate between to classes (FILES1 and FILES2)
    # year = 6,7,8
    # tmva_type = 0 (Ktautau sig vs bkg); 1 (DDK sig vs background)
    # tmva_method = 0 (BDT gradient boost) ; 1 (BDT adaptive boost); 2 (MLP); 3 (DNN)
    # tmva_vars = 0,1,2 (which set of input variables to use?)
    # background proxy = 0 (RS data, right sideband); (WS data, signal region)
    # FILES1 (signal), FILES2 (background) -> files containing the 2 species we want to separate (CHANGE ME for each tmva_type && background_proxy)
    input:
        'TMVA/TMVA_training.C'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_output_type{tmva_type}_method{tmva_method}_vars{tmva_vars}_bkgProxy{background_proxy}.root' 
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_output_type{tmva_type}_method{tmva_method}_vars{tmva_vars}_bkgProxy{background_proxy}.log' 
    shell:
        'root -l -b -q \'TMVA/TMVA_training.C( {wildcards.year}, {wildcards.tmva_type}, {wildcards.tmva_vars}, {wildcards.tmva_method}, {wildcards.background_proxy} )\' &> {log};'
        'rm -rf /panfs/felician/B2Ktautau/workflow/TMVA_training/201{wildcards.year}/tmva_output_type{wildcards.tmva_type}_method{wildcards.tmva_method}_vars{wildcards.tmva_vars}_bkgProxy{wildcards.background_proxy};' 
        'mv tmva_output_type{wildcards.tmva_type}_method{wildcards.tmva_method}_vars{wildcards.tmva_vars}_bkgProxy{wildcards.background_proxy} /panfs/felician/B2Ktautau/workflow/TMVA_training/201{wildcards.year}'

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
        weights = '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_output_type{tmva_type}_method{tmva_method}_vars{tmva_vars}_bkgProxy{background_proxy}/weights/TMVAClassification_BDT.weights.xml'
    resources:
        mem_mb = 25000
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_{species}/tmva_type{tmva_type}_method{tmva_method}_vars{tmva_vars}_bkgProxy{background_proxy}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_{species}/tmva_type{tmva_type}_method{tmva_method}_vars{tmva_vars}_bkgProxy{background_proxy}/{line}.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_response.C( {wildcards.year}, {wildcards.species}, {wildcards.tmva_method}, {wildcards.tmva_type}, {wildcards.tmva_vars}, {wildcards.background_proxy}, \"{input.weights}\", {wildcards.line} )\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/TMVA_response/201{wildcards.year}/Species_{wildcards.species}/tmva_type{wildcards.tmva_type}_method{wildcards.tmva_method}_vars{wildcards.tmva_vars}_bkgProxy{wildcards.background_proxy}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/TMVA_response/201{wildcards.year}/Species_{wildcards.species}/tmva_type{wildcards.tmva_type}_method{wildcards.tmva_method}_vars{wildcards.tmva_vars}_bkgProxy{wildcards.background_proxy}/TMVA_response.txt;'

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
        MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/TMVA_response.txt',
        MC_0_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_10/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/TMVA_response.txt',
        MC_1_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_11/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/TMVA_response.txt',
        MC_2_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_12/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/TMVA_response.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/TMVA_response.txt',
        WS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}/TMVA_response.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}/sig_bkg_comparison.gif'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/tmva_method{tmva_method}_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}/out.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_compare.C( {wildcards.year}, {wildcards.tmva_method}, {wildcards.tmva_type}, {wildcards.background_proxy}, {wildcards.num_data_dst_files}, \"{input.MC_files}\", \"{input.MC_0_files}\", \"{input.MC_1_files}\", \"{input.MC_2_files}\", \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.MC_0_BDT_files}\", \"{input.MC_1_BDT_files}\", \"{input.MC_2_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.WS_DATA_BDT_files}\" )\' > {log}'

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

#########################################################    SkLearn    ######################################################################################################## 

rule sklearn_training:
    ''' Classifier to seprate signal from combinatorial background '''
    input:
        'Sklearn/make_sklearn_training.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/clf_first_step.pkl',
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/clf_second_step.pkl'
    log:
        '/panfs/felician/B2Ktautau/workflow/sklearn_training/out.log'
    shell:
        'python -u Sklearn/make_sklearn_training.py 2>&1 | tee {log}'

rule sklearn_response:
    ''' Compute sklearn response '''
    input:
        'Sklearn/compute_sklearn_response.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/sklearn_response/201{year}/Species_{species}/{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/sklearn_response/201{year}/Species_{species}/{line}.log'
    shell:
        'python -u Sklearn/compute_sklearn_response.py {wildcards.year} {wildcards.species} {wildcards.line} 2>&1 | tee {log};'
        'ls /panfs/felician/B2Ktautau/workflow/sklearn_response/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/sklearn_response/201{wildcards.year}/Species_{wildcards.species}/bdt_output.txt'


##############################################################################################################################################################################
rule compute_bkg_yield:
    ''' Computes the background yield in the signal region '''
    input:
        'compute_B_yield.C',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_DATA_pre_selection_tree/201{year}/Species_2/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.log'
    shell: 
        'root -l -b -q \'compute_B_yield.C( {wildcards.year}, \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\" )\' > {log}'

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
        MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_method1_type6_bkgProxy1/TMVA_response.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_method1_type6_bkgProxy1/TMVA_response.txt',
        BKG_YIELD_file = '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.log'
    shell:
        'root -l -b -q \'sensitivity.C( {wildcards.year}, {wildcards.nprob}, {wildcards.bdt_it}, \"{input.MC_files}\", \"{input.RS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.BKG_YIELD_file}\" ) \' > {log};'
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
        