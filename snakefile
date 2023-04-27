localrules: make_TMVA_training, compare_TMVA, compute_bkg_yield, get_run_numbers, get_lumi_of_subset, sensitivity_plots, all

rule all:
    input:
        expand('/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root', year=6)

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
        'ls /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/*.root > /panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{wildcards.year}/mc_components.txt'

rule create_pre_sel_tree:
    ''' Creates a tree passing the pre-selections. One can choose the species: MC (0), RS data (1) or WS data (2). One can choose the component: 3pi3pi (0), 3pi3pipi0 (1), 3pi3pi2pi0 (2) or all (-1) for MC. '''
    input:
        'create_pre_sel_tree.C',
        MC_files = 'Files_on_grid/full_MC_201{year}.txt',
        MC_component_file = '/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/201{year}/MC_component_{line}.root',
        DATA_files = 'Files_on_grid/full_data_201{year}.txt'
    output:
        #'/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_{species}/Component_{component}/{line}.root' # MC
        '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_{species}/{line}.root' # data
    log:
        #'/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_{species}/Component_{component}/{line}.log' # MC
        '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_{species}/{line}.log' # data
    shell:
        #'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, {wildcards.species}, {wildcards.component}, \"{input.MC_files}\", \"{input.MC_component_file}\", \"{input.DATA_files}\", {wildcards.line})\' &> {log};' # MC
        'root -l -b -q \'create_pre_sel_tree.C( {wildcards.year}, {wildcards.species}, 0, \"{input.MC_files}\", \"{input.MC_component_file}\", \"{input.DATA_files}\", {wildcards.line})\' &> {log};' # data
        #'ls /panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{wildcards.year}/Species_{wildcards.species}/Component_{wildcards.component}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{wildcards.year}/Species_{wildcards.species}/Component_{wildcards.component}/pre_sel_tree.txt' # MC
        'ls /panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{wildcards.year}/Species_{wildcards.species}/pre_sel_tree.txt' # data

rule make_TMVA_training:
    ''' Does the TMVA training. Uses only 1 data grid file for training. Uses all MC grid files for training. '''
    # tmva_type (input variables used in the training)
    # 0 -> taus (kinematic)
    # 1 -> taus (isolation)
    # 2 -> B (kinematic)
    # 3 -> B (isolation)
    # 4 -> tau+ (isolation)
    # 5 -> tau- (isolation)
    # 6 -> Mix (initial BDT)
    # background_proxy
    # 0 -> use RS data, with M > 6.5 GeV
    # 1 -> use WS data
    input:
        'TMVA/TMVA_training.C',
        SIGNAL_MC_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_0/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_2/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_output_type{tmva_type}_bkgProxy{background_proxy}.root',
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/log_type{tmva_type}_bkgProxy{background_proxy}.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_training.C( {wildcards.year}, {wildcards.tmva_type}, {wildcards.background_proxy}, \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.SIGNAL_MC_files}\")\' &> {log};'
        'rm -rf /panfs/felician/B2Ktautau/workflow/TMVA_training/201{wildcards.year}/tmva_dataset_type{wildcards.tmva_type}_bkgProxy{wildcards.background_proxy};'
        'mv tmva_dataset_type{wildcards.tmva_type}_bkgProxy{wildcards.background_proxy} /panfs/felician/B2Ktautau/workflow/TMVA_training/201{wildcards.year}'

rule evaluate_TMVA_response:
    ''' Evaluates the BDT response.'''
    # species (in which files will the BDT response be evaluated)
    # 0 -> MC (all 3 components, component=-1)
    # 10 -> MC (3pi3pi)
    # 11 -> MC (3pi3pi pi0)
    # 12 -> MC (3pi3pi 2pi0) 
    # 1 -> RS data
    # 2 -> WS data
    input:
        'TMVA/TMVA_BDT_output.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_-1/pre_sel_tree.txt',
        MC_0_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_0/pre_sel_tree.txt',
        MC_1_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_1/pre_sel_tree.txt',
        MC_2_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_2/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_2/pre_sel_tree.txt',
        weights = '/panfs/felician/B2Ktautau/workflow/TMVA_training/201{year}/tmva_dataset_type{tmva_type}_bkgProxy{background_proxy}/weights/TMVAClassification_BDT.weights.xml'
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_{species}/BDT_output_tmva_type{tmva_type}_bkgProxy{background_proxy}_{line}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_{species}/tmva_type{tmva_type}_bkgProxy{background_proxy}_{line}.log'
    shell:
        'root -l -b -q \'TMVA/TMVA_BDT_output.C( {wildcards.year}, {wildcards.species}, {wildcards.tmva_type}, {wildcards.background_proxy}, \"{input.MC_files}\", \"{input.MC_0_files}\", \"{input.MC_1_files}\", \"{input.MC_2_files}\", \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.weights}\", {wildcards.line} )\' &> {log};'
        'ls /panfs/felician/B2Ktautau/workflow/TMVA_response/201{wildcards.year}/Species_{wildcards.species}/*.root | sort -V > /panfs/felician/B2Ktautau/workflow/TMVA_response/201{wildcards.year}/Species_{wildcards.species}/tmva_type{wildcards.tmva_type}_bkgProxy{wildcards.background_proxy}.txt;'

rule compare_TMVA:
    ''' Compares BDT distributions evaluated in MC, RS and WS data '''
    input:
        'TMVA/compare_BDTs.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_-1/pre_sel_tree.txt',
        MC_0_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_0/pre_sel_tree.txt',
        MC_1_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_1/pre_sel_tree.txt',
        MC_2_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_2/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_2/pre_sel_tree.txt',
        MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_0_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_10/tmva_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_1_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_11/tmva_type{tmva_type}_bkgProxy{background_proxy}.txt',
        MC_2_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_12/tmva_type{tmva_type}_bkgProxy{background_proxy}.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_type{tmva_type}_bkgProxy{background_proxy}.txt',
        WS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_2/tmva_type{tmva_type}_bkgProxy{background_proxy}.txt',
    output:
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/sig_bkg_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.gif',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/sig_bkg_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.pdf',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/mc_components_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.gif',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/mc_components_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.pdf',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/RS_WS_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.gif',
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/RS_WS_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}_comparison.pdf'
    log:
        '/panfs/felician/B2Ktautau/workflow/TMVA_plots_compare/201{year}/compare_BDTs_tmva_type{tmva_type}_bkgProxy{background_proxy}_numDST{num_data_dst_files}.log'
    shell:
        'root -l -b -q \'TMVA/compare_BDTs.C( {wildcards.year}, {wildcards.tmva_type}, {wildcards.background_proxy}, {wildcards.num_data_dst_files}, \"{input.MC_files}\", \"{input.MC_0_files}\", \"{input.MC_1_files}\", \"{input.MC_2_files}\", \"{input.RS_DATA_files}\", \"{input.WS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.MC_0_BDT_files}\", \"{input.MC_1_BDT_files}\", \"{input.MC_2_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.WS_DATA_BDT_files}\")\' > {log}'

rule compute_bkg_yield:
    ''' Computes the background yield in the signal region '''
    input:
        'compute_B_yield.C',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_1/pre_sel_tree.txt',
        WS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_2/pre_sel_tree.txt'
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
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_1/pre_sel_tree.txt'
    output:
        '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.txt'
    log:
        '/panfs/felician/B2Ktautau/workflow/get_run_numbers/201{year}/run_numbers_in_subset.log'
    shell:
        'root -l -b -q \'run_numbers_in_subset.C( {wildcards.year}, \"{input.RS_DATA_files}\" )\' > {log}'

rule get_lumi_of_subset:
    ''' Computes luminosity of a .dst file of RS data (after pre-selections) '''
    input:
        'lumi_of_subset.py'
    output:
        '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.log'
    shell:
        'python lumi_of_subset.py {wildcards.year} > {log};'

rule sensitivity:
    ''' Computes sensitivity for each BDT cut and each line of MC / RS data files'''
    input:
        'sensitivity.C',
        MC_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_0/Component_-1/pre_sel_tree.txt',
        RS_DATA_files = '/panfs/felician/B2Ktautau/workflow/create_pre_sel_tree/201{year}/Species_1/pre_sel_tree.txt',
        MC_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_0/tmva_type6_bkgProxy1.txt',
        RS_DATA_BDT_files = '/panfs/felician/B2Ktautau/workflow/TMVA_response/201{year}/Species_1/tmva_type6_bkgProxy1.txt',
        BKG_YIELD_file = '/panfs/felician/B2Ktautau/workflow/compute_bkg_yield/201{year}/bkg_yield.root',
        LUMI_file = '/panfs/felician/B2Ktautau/workflow/get_lumi_of_subset/201{year}/lumi_of_subset.root'
    output:
        '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.root'
    log:
        '/panfs/felician/B2Ktautau/workflow/sensitivity/201{year}/N_prob_{nprob}/BDT_cut_{bdt_it}.log'
    shell:
        'root -l -b -q \'sensitivity.C( {wildcards.year}, {wildcards.nprob}, {wildcards.bdt_it}, \"{input.MC_files}\", \"{input.RS_DATA_files}\", \"{input.MC_BDT_files}\", \"{input.RS_DATA_BDT_files}\", \"{input.BKG_YIELD_file}\", \"{input.LUMI_file}\" ) \' > {log};'
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