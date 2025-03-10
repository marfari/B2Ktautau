import ROOT
import numpy as np
import sys, os
sys.path.append(os.getcwd())

def main(argv):
    bdt1 = argv[1]
    bdt2 = argv[2]
    random_seed = argv[3]

    bdt1 = float(bdt1)
    bdt2 = float(bdt2)
    random_seed = int(random_seed)

    f = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t = f.Get("DecayTree")

    eps_ws_den = t.GetEntries()
    eps_ws_num = t.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    eps_ws = eps_ws_num/eps_ws_den

    # Expected number of combinatorial backgorund events in RS data after BDT cuts
    # Factor 2 comes from the fact that eps_rs = 2*eps_ws at BDT1 = BDT2 = 0.8
    N_bkg = 2*19006409*eps_ws

    string1 = '''
<!DOCTYPE Channel  SYSTEM 'HistFactorySchema.dtd'>

<Channel Name="channel1" InputFile="/panfs/felician/B2Ktautau/workflow/generate_fit_templates/fit_templates_bdt1_{0}_bdt2_{1}_seed_{2}.root" >
    <Data HistoName="data" HistoPath="Data" />
    <Sample Name="Signal" HistoPath="Signal" HistoName="nominal">
        <ShapeSys Name="staterror_sig" HistoPath="Signal" HistoName="error" ConstraintType="Poisson"/>
        <NormFactor Name="Nsig" Val="0" Low="-10" High="10."  />
    </Sample>
    <Sample Name="Background" HistoPath="Background" HistoName="nominal">
        <ShapeSys Name="staterror_bkg" HistoPath="Background" HistoName="error" ConstraintType="Poisson"/>
        <NormFactor Name="Nbkg" Val=\"{3}\" Low="0" High=\"{4}\"  />
    </Sample>
</Channel>
    '''.format( round(bdt1,1), round(bdt2,1), random_seed, N_bkg, 2*N_bkg)

    file1 = open('/panfs/felician/B2Ktautau/workflow/pyhf_fit/config/fit_region_bdt1_{0}_bdt2_{1}_seed_{2}.xml'.format( round(bdt1,1), round(bdt2,1), random_seed), 'w')
    file1.write(string1)
    file1.close()

    string2 = ''' 
<!DOCTYPE Combination  SYSTEM 'HistFactorySchema.dtd'>

<Combination OutputFilePrefix="/panfs/felician/B2Ktautau/workflow/pyhf_fit/workspace/workspace_bdt1_{0}_bdt2_{1}_seed_{2}">
	<Input>/panfs/felician/B2Ktautau/workflow/pyhf_fit/config/fit_region_bdt1_{0}_bdt2_{1}_seed_{2}.xml</Input>
	<Measurement Name="MinimalConfiguration" Lumi="1." LumiRelErr="0.1" ExportOnly="True">
		<POI>Nsig</POI>
	</Measurement>
</Combination>
    '''.format( round(bdt1,1), round(bdt2,1), random_seed )

    file2 = open('/panfs/felician/B2Ktautau/workflow/pyhf_fit/config/config_bdt1_{0}_bdt2_{1}_seed_{2}.xml'.format( round(bdt1,1), round(bdt2,1), random_seed), 'w')
    file2.write(string2)
    file2.close()

if __name__ == "__main__":
    main(sys.argv)