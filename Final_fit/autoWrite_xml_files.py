import ROOT
import numpy as np
import sys, os
sys.path.append(os.getcwd())
from uncertainties import ufloat
from uncertainties import unumpy

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
    eps_ws_num_1 = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(bdt1,bdt2))
    eps_ws_num_2 = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(bdt1,bdt2))
    eps_ws_num_3 = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(bdt1,bdt2))

    eps_ws = eps_ws_num/eps_ws_den
    eps_ws_1 = eps_ws_num_1/eps_ws_den
    eps_ws_2 = eps_ws_num_2/eps_ws_den
    eps_ws_3 = eps_ws_num_3/eps_ws_den

    N_bkg = 2*19006409*eps_ws
    N_bkg_1 = 2*19006409*eps_ws_1
    N_bkg_2 = 2*19006409*eps_ws_2
    N_bkg_3 = 2*19006409*eps_ws_3
    N_bkg_vector = [N_bkg, N_bkg_1, N_bkg_2, N_bkg_3]
    n_channels = len(N_bkg_vector)

    fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy")
    fixed_value_error = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs_error.npy")
    bdt_independent_constant = ufloat( fixed_value, fixed_value_error )

    # MC (w/o TM cuts)
    f_mc = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_mc = f_mc.Get("DecayTree")

    N_num = t_mc.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    N_num_1 = t_mc.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR >= 0) && (df_Bp_MERR <= 100)".format(bdt1,bdt2))
    N_num_2 = t_mc.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 100) && (df_Bp_MERR <= 250)".format(bdt1,bdt2))
    N_num_3 = t_mc.GetEntries("(BDT1 > {0}) && (BDT2 > {1}) && (df_Bp_MERR > 250)".format(bdt1,bdt2))

    numerator = ufloat( N_num, np.sqrt(N_num) )
    numerator_1 = ufloat( N_num_1, np.sqrt(N_num_1) )
    numerator_2 = ufloat( N_num_2, np.sqrt(N_num_2) )
    numerator_3 = ufloat( N_num_3, np.sqrt(N_num_3) )

    c = numerator/bdt_independent_constant
    c_1 = numerator_1/bdt_independent_constant
    c_2 = numerator_2/bdt_independent_constant
    c_3 = numerator_3/bdt_independent_constant
    c_vector = [c.nominal_value, c_1.nominal_value, c_2.nominal_value, c_3.nominal_value]
    c_vector_err = [c.std_dev, c_1.std_dev, c_2.std_dev, c_3.std_dev]

    upward = np.zeros(n_channels)
    downward = np.zeros(n_channels)

    for i in range(n_channels):
        if(c_vector[i] != 0):
            upward[i] = 1+(c_vector_err[i]/c_vector[i])
            downward[i] = 1-(c_vector_err[i]/c_vector[i])

    strings = []

    for i in range(n_channels):
        if(i == 0):
            filename = "/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{0}/fit_templates_bdt1_{1}_bdt2_{2}_seed_{3}.root".format(folder_name, bdt1, bdt2, random_seed)
            channel_name = "AllEvents"
        else:
            filename = "/panfs/felician/B2Ktautau/workflow/generate_fit_templates/{0}/fit_templates_bdt1_{1}_bdt2_{2}_seed_{3}_ch{4}.root".format(folder_name, bdt1, bdt2, random_seed, i)
            channel_name = "Mass_error_{0}".format(i)

        string = '''
<!DOCTYPE Channel  SYSTEM 'HistFactorySchema.dtd'>

<Channel Name="{1}" InputFile="{0}">
    <Data HistoName="data" HistoPath="Data" />

    <!-- StatErrorConfig RelErrorThreshold="0.05" ConstraintType="Poisson" -->

    <Sample Name="Signal" HistoPath="Signal" HistoName="nominal">
        <NormFactor Name="c_bdt1_bdt2_{2}" Val=\"{3}\" Low="0" High="1000000000." />
        <NormFactor Name="BF_sig" Val="0" Low="-100." High="100." />
        <OverallSys Name="c_bdt1_bdt2_syst_{2}" High=\"{6}\" Low=\"{7}\" ConstraintType="Gaussian" />
        <ShapeSys Activate="True" Name="shapesys_sig_{2}" HistoPath="Signal" HistoName="error" ConstraintType="Poisson" />
    </Sample>
    <Sample Name="Background" HistoPath="Background" HistoName="nominal">
        <NormFactor Name="Nbkg_{2}" Val=\"{4}\" Low="0" High=\"{5}\" />
        <ShapeSys Activate="True" Name="shapesys_bkg_{1}" HistoPath="Background" HistoName="error" ConstraintType="Poisson" />
    </Sample>
</Channel>
        '''.format(filename, channel_name, i, c_vector[i],  N_bkg_vector[i], 10*N_bkg_vector[i], upward[i], downward[i])

        strings.append(string)

    file = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}.xml'.format( folder_name, bdt1, bdt2, random_seed), 'w')
    file.write(strings[0])
    file.close()

    file1 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}_ch1.xml'.format( folder_name, bdt1, bdt2, random_seed), 'w')
    file1.write(strings[1])
    file1.close()

    file2 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}_ch2.xml'.format( folder_name, bdt1, bdt2, random_seed), 'w')
    file2.write(strings[2])
    file2.close()

    file3 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}_ch3.xml'.format( folder_name, bdt1, bdt2, random_seed), 'w')
    file3.write(strings[3])
    file3.close()

    string4 = ''' 
<!DOCTYPE Combination  SYSTEM 'HistFactorySchema.dtd'>

<Combination OutputFilePrefix="/panfs/felician/B2Ktautau/workflow/generate_fit_workspaces/{0}/workspace_bdt1_{1}_bdt2_{2}_seed_{3}">
	<Input>/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}.xml</Input>
	<Measurement Name="MinimalConfiguration" Lumi="1." LumiRelErr="0.1">
		<POI>BF_sig</POI>
        <ParamSetting Const="True">Lumi c_bdt1_bdt2_0 </ParamSetting>
	</Measurement>
</Combination>
    '''.format( folder_name, bdt1, bdt2, random_seed)

    string5 = ''' 
<!DOCTYPE Combination  SYSTEM 'HistFactorySchema.dtd'>

<Combination OutputFilePrefix="/panfs/felician/B2Ktautau/workflow/generate_fit_workspaces/{0}/workspace_bdt1_{1}_bdt2_{2}_seed_{3}_ch">
	<Input>/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}_ch1.xml</Input>
	<Input>/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}_ch2.xml</Input>
	<Input>/panfs/felician/B2Ktautau/workflow/write_xml_files/{0}/fit_region_bdt1_{1}_bdt2_{2}_seed_{3}_ch3.xml</Input>
	<Measurement Name="ErrorCategories" Lumi="1." LumiRelErr="0.1">
		<POI>BF_sig</POI>
        <ParamSetting Const="True">Lumi c_bdt1_bdt2_1 c_bdt1_bdt2_2 c_bdt1_bdt2_3 </ParamSetting>
	</Measurement>
</Combination>
    '''.format( folder_name, bdt1, bdt2, random_seed)

    file4 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{3}/config_bdt1_{0}_bdt2_{1}_seed_{2}.xml'.format( bdt1, bdt2, random_seed, folder_name ), 'w')
    file4.write(string4)
    file4.close()

    file5 = open('/panfs/felician/B2Ktautau/workflow/write_xml_files/{3}/config_bdt1_{0}_bdt2_{1}_seed_{2}_ch.xml'.format( bdt1, bdt2, random_seed, folder_name ), 'w')
    file5.write(string5)
    file5.close()


if __name__ == "__main__":
    main(sys.argv)