import ROOT
import sys
import pickle
import pandas as pd
import numpy as np
from array import array 

def main(argv):

    year = argv[1]
    species = argv[2]
    line = argv[3]

    year = int(year)
    species = int(species)
    line = int(line)

    fc = ROOT.TFileCollection("fc", "fc", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_{1}/pre_sel_tree.txt".format(year,species), 1, line)
    fc1 = ROOT.TFileCollection("fc1", "fc1", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{0}/Species_{1}/fit_results.txt".format(year,species), 1, line)

    t = ROOT.TChain("DecayTree")
    t1 = ROOT.TChain("DecayTree")

    t.AddFileInfoList(fc.GetList())
    t1.AddFileInfoList(fc1.GetList())
    t.AddFriend(t1)

    df = ROOT.RDataFrame(t)

    X1 = pd.DataFrame() # 1st step
    X2 = pd.DataFrame() # 2nd step

    branch_names_first_step = ['Bp_VTXISONUMVTX_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 
                    'Bp_VTXISONUMVTX_taup', 'Bp_VTXISONUMVTX_taum', 'Bp_VTXISODCHI2ONETRACK_taup', 'Bp_VTXISODCHI2ONETRACK_taum', 
                    'Bp_VTXISODCHI2MASSONETRACK_taup', 'Bp_VTXISODCHI2MASSONETRACK_taum', 'Bp_VTXISODCHI2TWOTRACK_taup', 'Bp_VTXISODCHI2TWOTRACK_taum',
                    'Bp_VTXISODCHI2MASSTWOTRACK_taup', 'Bp_VTXISODCHI2MASSTWOTRACK_taum', 'Bp_CC_05_DELTAPHI_B', 'Bp_NC_05_IT_B', 'Bp_CC_05_PTASYM_B', 'Bp_CC_05_PX_B',
                    'Bp_CC_05_PZASYM_B', 'Bp_CC_05_DELTAETA_B', 'Bp_NC_05_MULT_B', 'Bp_NC_05_VPT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_B', 'Bp_CC_05_MULT_B',
                    'Bp_CC_05_PASYM_taup', 'Bp_CC_05_DELTAETA_taup', 'Bp_NC_05_MULT_taup', 'Bp_NC_05_SPT_taup', 'Bp_CCNC_05_IT_taup', 'Bp_NC_05_PTASYM_taup', 'Bp_CC_05_PTASYM_taup',
                    'Bp_CC_05_PASYM_taum', 'Bp_CC_05_DELTAETA_taum', 'Bp_NC_05_MULT_taum', 'Bp_NC_05_SPT_taum', 'Bp_CCNC_05_IT_taum', 'Bp_NC_05_PTASYM_taum', 'Bp_CC_05_PTASYM_taum'] 

    branch_names_second_step = ['taup_M12', 'taum_M12', 'taup_M23', 'taum_M23', 'taup_DIRA_ORIVX', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                                'df_BVx', 'df_BVy', 'df_BVz', 'df_DV1x', 'df_DV1y', 'df_DV1z', 'df_DV2x', 'df_DV2y', 'df_DV2z',
                                'df_Kp_PX', 'df_Kp_PY', 'df_Kp_PZ', 'df_PVx', 'df_PVy', 'df_PVz', 'taup_M', 'taum_M', 'df_chi2']

    branch_names = branch_names_first_step + branch_names_second_step
    x  = df.AsNumpy(branch_names)

    ###################################################################### 1st step ##########################################################################
    X1['VTXISODCHI2ONETRACK_B'] = x['Bp_VTXISODCHI2ONETRACK_B']
    X1['VTXISODCHI2MASSONETRACK_B'] = x['Bp_VTXISODCHI2MASSONETRACK_B']
    X1['VTXISODCHI2TWOTRACK_B'] = x['Bp_VTXISODCHI2TWOTRACK_B']
    # X1['VTXISODCHI2MASSTWOTRACK_B'] = x['Bp_VTXISODCHI2MASSTWOTRACK_B']
    # X1['VTXISONUMVTX_tau_max'] = np.maximum( x['Bp_VTXISONUMVTX_taup'], x['Bp_VTXISONUMVTX_taum'] )
    # X1['VTXISODCHI2MASSONETRACK_tau_min'] = np.minimum( x['Bp_VTXISODCHI2MASSONETRACK_taup'], x['Bp_VTXISODCHI2MASSONETRACK_taum'] ) 
    X1['VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( x['Bp_VTXISODCHI2TWOTRACK_taup'], x['Bp_VTXISODCHI2TWOTRACK_taum'] ) 
    # X1['VTXISODCHI2MASSTWOTRACK_tau_min'] = np.minimum( x['Bp_VTXISODCHI2MASSTWOTRACK_taup'], x['Bp_VTXISODCHI2MASSTWOTRACK_taum'] )
    # X1['CC_05_MULT_B'] = x['Bp_CC_05_MULT_B']
    # X1['CC_05_DELTAETA_B'] = x['Bp_CC_05_DELTAETA_B']
    X1['NC_05_PTASYM_B'] = x['Bp_NC_05_PTASYM_B']
    X1['CCNC_05_IT_B'] = x['Bp_CCNC_05_IT_B']
    X1['CC_05_PTASYM_tau_max'] = np.maximum( x['Bp_CC_05_PTASYM_taup'], x['Bp_CC_05_PTASYM_taum'] )
    # X1['CCNC_05_IT_tau_min'] = np.minimum( x['Bp_CCNC_05_IT_taup'], x['Bp_CCNC_05_IT_taum'] )
    ############################################################################ 2nd step #####################################################################
    X2['taup_M'] = x['taup_M']
    X2['taum_M'] = x['taum_M']
    # X2['taup_M12'] = x['taup_M12']
    X2['taum_M23'] = x['taum_M23'] 
    # X2['log_1_minus_tau_DIRA_BV_max'] = np.maximum( np.log(1 - np.abs(x['taup_DIRA_ORIVX'] ))*np.sign( x['taup_DIRA_ORIVX']),  np.log(1 - np.abs(x['taum_DIRA_ORIVX'] ))*np.sign( x['taup_DIRA_ORIVX'] ) )
    X2['log_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log(1 - np.abs(x['taup_DIRA_ORIVX'] ))*np.sign( x['taup_DIRA_ORIVX']),  np.log(1 - np.abs(x['taum_DIRA_ORIVX'] ))*np.sign( x['taup_DIRA_ORIVX'] ) )
    X2['log_1_minus_B_DIRA_PV'] = np.log(1 - np.abs(x['Bp_DIRA_OWNPV']) )
    # X2['B_FD_PV'] = np.sqrt( (x['df_BVx'] - x['df_PVx'])**2 + (x['df_BVy'] - x['df_PVy'])**2 + (x['df_BVz'] - x['df_PVz'])**2 )
    X2['tau_FD_BV_max'] = np.maximum( np.sqrt( (x['df_DV1x'] - x['df_BVx'])**2 + (x['df_DV1y'] - x['df_BVy'])**2 + (x['df_DV1z'] - x['df_BVz'])**2 ), np.sqrt( (x['df_DV2x'] - x['df_BVx'])**2 + (x['df_DV2y'] - x['df_BVy'])**2 + (x['df_DV2z'] - x['df_BVz'])**2 ) )
    X2['DV1_DV2_distance'] = np.sqrt( (x['df_DV1x'] - x['df_DV2x'])**2 + (x['df_DV1y'] - x['df_DV2y'])**2 + (x['df_DV1z'] - x['df_DV2z'])**2 )
    # X2['K_PT'] = np.sqrt( x['df_Kp_PX']**2 + x['df_Kp_PY']**2 )
    X2['df_chi2'] = x['df_chi2']

    Cx_taup =  (x['df_DV1y'] - x['df_BVy'])*x['df_Kp_PZ']  - ( x['df_DV1z'] - x['df_BVz'])*x['df_Kp_PY']
    Cy_taup =  (x['df_DV1z'] - x['df_BVz'])*x['df_Kp_PX']  - ( x['df_DV1x'] - x['df_BVx'])*x['df_Kp_PZ']
    Cz_taup =  (x['df_DV1x'] - x['df_BVx'])*x['df_Kp_PY']  - ( x['df_DV1y'] - x['df_BVy'])*x['df_Kp_PX']
    C_taup = np.sqrt( Cx_taup**2 + Cy_taup**2 + Cz_taup**2  )
    IP_taup_Kp = (2*C_taup)/( np.sqrt( x['df_Kp_PX']**2 + x['df_Kp_PY']**2 + x['df_Kp_PZ']**2 ) )

    Cx_taum =  (x['df_DV2y'] - x['df_BVy'])*x['df_Kp_PZ']  - ( x['df_DV2z'] - x['df_BVz'])*x['df_Kp_PY']
    Cy_taum =  (x['df_DV2z'] - x['df_BVz'])*x['df_Kp_PX']  - ( x['df_DV2x'] - x['df_BVx'])*x['df_Kp_PZ']
    Cz_taum =  (x['df_DV2x'] - x['df_BVx'])*x['df_Kp_PY']  - ( x['df_DV2y'] - x['df_BVy'])*x['df_Kp_PX']
    C_taum = np.sqrt( Cx_taum**2 + Cy_taum**2 + Cz_taum**2  )
    IP_taum_Kp = (2*C_taum)/( np.sqrt( x['df_Kp_PX']**2 + x['df_Kp_PY']**2 + x['df_Kp_PZ']**2 ) )

    X2['IP_tau_Kp_max'] = np.maximum( IP_taup_Kp, IP_taum_Kp ) 

    a = np.sqrt( ( x['df_PVx'] - x['df_DV1x'] )**2 + ( x['df_PVy'] - x['df_DV1y'] )**2 + ( x['df_PVz'] - x['df_DV1z'] )**2 )
    b = np.sqrt( ( x['df_PVx'] - x['df_DV2x'] )**2 + ( x['df_PVy'] - x['df_DV2y'] )**2 + ( x['df_PVz'] - x['df_DV2z'] )**2 )
    c = np.sqrt( ( x['df_DV1x'] - x['df_DV2x'] )**2 + ( x['df_DV1y'] - x['df_DV2y'] )**2 + ( x['df_DV1z'] - x['df_DV2z'] )**2 )
    s = (a+b+c)/2
    # X2['DV1_DV2_PV_area'] = np.sqrt( s + (s-a)*(s-b)*(s-c) )
    ##################################################################################################################################################################
    
    names = ["AdaBDT", "GradBDT", "RForest"]

    classifiers_first_step = []
    classifiers_second_step = []
    with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/clf_first_step.pkl', 'rb') as f:
        while True:
            try:
                classifiers_first_step.append(pickle.load(f))
            except EOFError:
                break
    with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/clf_second_step.pkl', 'rb') as f:
        while True:
            try:
                classifiers_second_step.append(pickle.load(f))
            except EOFError:
                break

    decisions_first_step = []
    decisions_second_step = []
    for i in range(len(classifiers_first_step)):
        if(names[i]=="RForest"):
            decisions_first_step.append( classifiers_first_step[i].predict_proba(X1)[:, 1] )
        else:
            decisions_first_step.append( classifiers_first_step[i].decision_function(X1) )
    for i in range(len(classifiers_second_step)):
        if(names[i]=="RForest"):
            decisions_second_step.append( classifiers_second_step[i].predict_proba(X2)[:, 1] )
        else:
            decisions_second_step.append( classifiers_second_step[i].decision_function(X2) )

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/sklearn_response/201{0}/Species_{1}/{2}.root".format(year,species,line), "UPDATE")
    fout.cd()
    trees = []
    clf_outputs_first_step = []
    clf_outputs_second_step = []

    for i in range(len(names)):
        fout.mkdir(names[i])
        trees.append(ROOT.TTree("DecayTree", "Decaytree"))
        
        clf_outputs_first_step.append( array('d', [0]) )
        clf_outputs_second_step.append( array('d', [0]) )

        trees[i].Branch("BDT1", clf_outputs_first_step[i], "BDT1/D")
        trees[i].Branch("BDT2", clf_outputs_second_step[i], "BDT2/D")

        for j in range(len(decisions_first_step[i])):
            clf_outputs_first_step[i][0] = decisions_first_step[i][j]
            clf_outputs_second_step[i][0] = decisions_second_step[i][j]
            trees[i].Fill()

        fout.cd(names[i])
        trees[i].Write()
    fout.Close()

if __name__ == "__main__":
    main(sys.argv)