import ROOT
import numpy as np
import sys, os
sys.path.append(os.getcwd())

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]
    folder_name = argv[4]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    # WS data
    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_ws = f_ws.Get("DecayTree")

    eps_ws_den = t_ws.GetEntries()
    eps_ws_num = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws = eps_ws_num/eps_ws_den

    # Expected number of combinatorial backgorund events in RS data after BDT cuts
    # Factor 2 comes from the fact that eps_rs = 2*eps_ws at BDT1 = BDT2 = 0.8
    N_bkg = 2*19006409*eps_ws

    fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy")
    fixed_value_error = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs_error.npy")

    # MC (w/o TM cuts)
    f_mc = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_mc = f_mc.Get("DecayTree")
    N_num = t_mc.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    
    c_bdt1_bdt2 = N_num/fixed_value

    string1 = '''
<!DOCTYPE Channel  SYSTEM 'HistFactorySchema.dtd'>

<Channel Name="channel1" InputFile="/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{6}/fit_templates_bdt1_{0}_bdt2_{1}_seed_{2}.root" >
    <Data HistoName="data" HistoPath="Data" />
    
    <!-- Set the StatError type to Poisson.  Can also be Gaussian -->
    <!-- StatErrorConfig RelErrorThreshold="0.1" ConstraintType="Poisson" -->
    
    <Sample Name="Signal" HistoPath="Signal" HistoName="nominal">
        <ShapeSys Activate="True" Name="shapesys_sig" HistoPath="Signal" HistoName="error" />
        <!-- StatError Activate="True" -->
        <NormFactor Name="c_bdt1_bdt2" Val=\"{5}\" Low="0" High="1000000000." />
        <NormFactor Name="BF_sig" Val="0" Low="-100." High="100." />
    </Sample>
    <Sample Name="Background" HistoPath="Background" HistoName="nominal">
        <ShapeSys Activate="True" Name="shapesys_bkg" HistoPath="Background" HistoName="error" />
        <!-- StatError Activate="True" -->
        <NormFactor Name="Nbkg" Val=\"{3}\" Low="0" High=\"{4}\" />
    </Sample>
</Channel>
    '''.format( bdt1, bdt2, random_seed, N_bkg, 10*N_bkg, c_bdt1_bdt2, folder_name)

    file1 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}.xml'.format( folder_name, bdt1, bdt2, random_seed), 'w')
    file1.write(string1)
    file1.close()


    string2 = ''' 
<!DOCTYPE Combination  SYSTEM 'HistFactorySchema.dtd'>

<Combination OutputFilePrefix="/panfs/felician/B2Ktautau/workflow/generate_fit_workspaces/{0}/workspace_bdt1_{1}_bdt2_{2}_seed_{3}">
	<Input>/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}.xml</Input>
	<Measurement Name="MinimalConfiguration" Lumi="1." LumiRelErr="0.1">
		<POI>BF_sig</POI>
        <ParamSetting Const="True">Lumi c_bdt1_bdt2 </ParamSetting>
	</Measurement>
</Combination>
    '''.format( folder_name, bdt1, bdt2, random_seed )

    file2 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{3}/config_bdt1_{0}_bdt2_{1}_seed_{2}.xml'.format( bdt1, bdt2, random_seed, folder_name ), 'w')
    file2.write(string2)
    file2.close()

if __name__ == "__main__":
    main(sys.argv)