import ROOT
import sys
import pickle
import pandas as pd
import numpy as np
from array import array 
from numpy import inf

import xgboost as xgb

def compute_response(year, species, line, df, fout, isKtautau, is_cocktailMC, name):
    X1 = pd.DataFrame() # 1st step
    X2 = pd.DataFrame() # 2nd step
    
    if((isKtautau == True) or (is_cocktailMC == True)):
        branch_names_isolation = ['Bp_VTXISONUMVTX_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 
                        'Bp_VTXISONUMVTX_taup', 'Bp_VTXISODCHI2ONETRACK_taup', 'Bp_VTXISODCHI2MASSONETRACK_taup', 'Bp_VTXISODCHI2TWOTRACK_taup', 'Bp_VTXISODCHI2MASSTWOTRACK_taup',
                        'Bp_VTXISONUMVTX_taum', 'Bp_VTXISODCHI2ONETRACK_taum', 'Bp_VTXISODCHI2MASSONETRACK_taum', 'Bp_VTXISODCHI2TWOTRACK_taum', 'Bp_VTXISODCHI2MASSTWOTRACK_taum',
                        'Bp_CC_05_DELTAPHI_B', 'Bp_NC_05_IT_B', 'Bp_CC_05_PTASYM_B', 'Bp_CC_05_PX_B', 'Bp_CC_05_IT_B', 'Bp_CC_05_PZASYM_B', 'Bp_CC_05_DELTAETA_B', 'Bp_NC_05_MULT_B', 'Bp_NC_05_VPT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_B', 'Bp_CC_05_MULT_B',
                        'Bp_CC_05_PASYM_taup', 'Bp_CC_05_DELTAETA_taup', 'Bp_NC_05_MULT_taup', 'Bp_NC_05_SPT_taup', 'Bp_CCNC_05_IT_taup', 'Bp_NC_05_PTASYM_taup', 'Bp_CC_05_PTASYM_taup', 'Bp_CC_05_IT_taup', 'Bp_TRKISOBDTFIRSTVALUE_taup_pi1', 'Bp_TRKISOBDTFIRSTVALUE_taup_pi2', 'Bp_TRKISOBDTFIRSTVALUE_taup_pi3', 'Bp_TRKISOBDTSECONDVALUE_taup_pi1', 'Bp_TRKISOBDTSECONDVALUE_taup_pi2', 'Bp_TRKISOBDTSECONDVALUE_taup_pi3', 'Bp_TRKISOBDTTHIRDVALUE_taup_pi1', 'Bp_TRKISOBDTTHIRDVALUE_taup_pi2', 'Bp_TRKISOBDTTHIRDVALUE_taup_pi3',
                        'Bp_CC_05_PASYM_taum', 'Bp_CC_05_DELTAETA_taum', 'Bp_NC_05_MULT_taum', 'Bp_NC_05_SPT_taum', 'Bp_CCNC_05_IT_taum', 'Bp_NC_05_PTASYM_taum', 'Bp_CC_05_PTASYM_taum', 'Bp_CC_05_IT_taum', 'Bp_TRKISOBDTFIRSTVALUE_taum_pi1', 'Bp_TRKISOBDTFIRSTVALUE_taum_pi2', 'Bp_TRKISOBDTFIRSTVALUE_taum_pi3', 'Bp_TRKISOBDTSECONDVALUE_taum_pi1', 'Bp_TRKISOBDTSECONDVALUE_taum_pi2', 'Bp_TRKISOBDTSECONDVALUE_taum_pi3', 'Bp_TRKISOBDTTHIRDVALUE_taum_pi1', 'Bp_TRKISOBDTTHIRDVALUE_taum_pi2', 'Bp_TRKISOBDTTHIRDVALUE_taum_pi3',
                        'Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taup', 'Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup', 'Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup', 
                        'Bp_B2Ksttautau_ISOBDTFIRSTVALUE_taum', 'Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum', 'Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] 

        branch_names_topology = ['taup_M', 'taup_M12', 'taup_M23', 'taup_M13', 'taup_DIRA_ORIVX', 'taum_M', 'taum_M12', 'taum_M23', 'taum_M13', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                        'df_BVx', 'df_BVy', 'df_BVz', 'df_DV1x', 'df_DV1y', 'df_DV1z', 'df_DV2x', 'df_DV2y', 'df_DV2z', 'df_PVx', 'df_PVy', 'df_PVz',
                        'df_BVx_err', 'df_BVy_err', 'df_BVz_err', 'df_DV1x_err', 'df_DV1y_err', 'df_DV1z_err', 'df_DV2x_err', 'df_DV2y_err', 'df_DV2z_err', 'df_PVx_err', 'df_PVy_err', 'df_PVz_err',
                        'df_Kp_PX', 'df_Kp_PY', 'df_Kp_PZ', 'df_chi2', 
                        'Kp_IPCHI2_OWNPV', 'taup_pi1_IPCHI2_OWNPV', 'taup_pi2_IPCHI2_OWNPV', 'taup_pi3_IPCHI2_OWNPV', 'taum_pi1_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV', 'taum_pi3_IPCHI2_OWNPV',
                        'taup_AMAXDOCA', 'taup_AMINDOCA', 'taup_DOCACHI2MAX', 'taum_AMAXDOCA', 'taum_AMINDOCA', 'taum_DOCACHI2MAX',
                        'Bp_M01', 'Bp_M02', 'Bp_M03', 'Bp_M04', 'Bp_M05', 'Bp_M06']
    else:
        branch_names_isolation = ['Bp_CC_05_IT_B',  'Bp_VTXISONUMVTX_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2MASSONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISODCHI2MASSTWOTRACK_B', 
                        'Bp_VTXISONUMVTX_D0bar', 'Bp_VTXISODCHI2ONETRACK_D0bar', 'Bp_VTXISODCHI2MASSONETRACK_D0bar', 'Bp_VTXISODCHI2TWOTRACK_D0bar', 'Bp_VTXISODCHI2MASSTWOTRACK_D0bar', 
                        'Bp_VTXISONUMVTX_Dsp', 'Bp_VTXISODCHI2ONETRACK_Dsp', 'Bp_VTXISODCHI2MASSONETRACK_Dsp', 'Bp_VTXISODCHI2TWOTRACK_Dsp', 'Bp_VTXISODCHI2MASSTWOTRACK_Dsp',
                        'Bp_CC_05_DELTAPHI_B', 'Bp_NC_05_IT_B', 'Bp_CC_05_PTASYM_B', 'Bp_CC_05_PX_B', 'Bp_CC_05_PZASYM_B', 'Bp_CC_05_DELTAETA_B', 'Bp_NC_05_MULT_B', 'Bp_NC_05_VPT_B', 'Bp_NC_05_PTASYM_B', 'Bp_CCNC_05_IT_B', 'Bp_CC_05_MULT_B',
                        'Bp_CC_05_PASYM_D0bar', 'Bp_CC_05_DELTAETA_D0bar', 'Bp_NC_05_MULT_D0bar', 'Bp_NC_05_SPT_D0bar', 'Bp_CCNC_05_IT_D0bar', 'Bp_NC_05_PTASYM_D0bar', 'Bp_CC_05_PTASYM_D0bar', 'Bp_CC_05_IT_D0bar', 'Bp_TRKISOBDTSECONDVALUE_D0bar_K', 'Bp_TRKISOBDTSECONDVALUE_D0bar_pi', 'Bp_TRKISOBDTFIRSTVALUE_D0bar_K', 'Bp_TRKISOBDTFIRSTVALUE_D0bar_pi', 'Bp_TRKISOBDTTHIRDVALUE_D0bar_K', 'Bp_TRKISOBDTTHIRDVALUE_D0bar_pi',
                        'Bp_CC_05_PASYM_Dsp', 'Bp_CC_05_DELTAETA_Dsp', 'Bp_NC_05_MULT_Dsp', 'Bp_NC_05_SPT_Dsp', 'Bp_CCNC_05_IT_Dsp', 'Bp_NC_05_PTASYM_Dsp', 'Bp_CC_05_PTASYM_Dsp', 'Bp_CC_05_IT_Dsp', 'Bp_TRKISOBDTSECONDVALUE_Dsp_K1', 'Bp_TRKISOBDTSECONDVALUE_Dsp_K2', 'Bp_TRKISOBDTSECONDVALUE_Dsp_pi', 'Bp_TRKISOBDTFIRSTVALUE_Dsp_K1', 'Bp_TRKISOBDTFIRSTVALUE_Dsp_K2', 'Bp_TRKISOBDTFIRSTVALUE_Dsp_pi', 'Bp_TRKISOBDTTHIRDVALUE_Dsp_K1', 'Bp_TRKISOBDTTHIRDVALUE_Dsp_K2', 'Bp_TRKISOBDTTHIRDVALUE_Dsp_pi',
                        'Bp_Bstautau_ISOBDTFIRSTVALUE_taup', 'Bp_Bstautau_ISOBDTSECONDVALUE_taup', 'Bp_Bstautau_ISOBDTTHIRDVALUE_taup', 'Bp_Bstautau_ISOBDTFIRSTVALUE_taum', 
                        'Bp_Bstautau_ISOBDTSECONDVALUE_taum', 'Bp_Bstautau_ISOBDTTHIRDVALUE_taum'] 

        branch_names_topology = ['D0bar_M', 'D0bar_DIRA_ORIVX', 'Dsp_M', 'Dsp_M12', 'Dsp_M23', 'Dsp_M13', 'Dsp_DIRA_ORIVX', 'Bp_DIRA_OWNPV',
                'Bp_ENDVERTEX_X', 'Bp_ENDVERTEX_Y', 'Bp_ENDVERTEX_Z', 'D0bar_ENDVERTEX_X', 'D0bar_ENDVERTEX_Y', 'D0bar_ENDVERTEX_Z', 'Dsp_ENDVERTEX_X', 'Dsp_ENDVERTEX_Y', 'Dsp_ENDVERTEX_Z', 'Bp_OWNPV_X', 'Bp_OWNPV_Y', 'Bp_OWNPV_Z',
                'Bp_ENDVERTEX_XERR', 'Bp_ENDVERTEX_YERR', 'Bp_ENDVERTEX_ZERR', 'D0bar_ENDVERTEX_XERR', 'D0bar_ENDVERTEX_YERR', 'D0bar_ENDVERTEX_ZERR', 'Dsp_ENDVERTEX_XERR', 'Dsp_ENDVERTEX_YERR', 'Dsp_ENDVERTEX_ZERR', 'Bp_OWNPV_XERR', 'Bp_OWNPV_YERR', 'Bp_OWNPV_ZERR', 'Bp_dtf_chi2', 
                'D0bar_AMAXDOCA', 'D0bar_AMINDOCA', 'D0bar_DOCACHI2MAX', 'Dsp_AMAXDOCA', 'Dsp_AMINDOCA', 'Dsp_DOCACHI2MAX', 'Bp_FDCHI2_OWNPV', 'D0bar_FD_ORIVX', 'Dsp_FD_ORIVX']

    branch_names = branch_names_isolation + branch_names_topology

    x  = df.AsNumpy(branch_names)

    ###################################################################### Isolation ##########################################################################
    if((isKtautau == True) or (is_cocktailMC == True)):
        X1['VTXISODCHI2ONETRACK_tau_min'] = np.minimum( x['Bp_VTXISODCHI2ONETRACK_taup'], x['Bp_VTXISODCHI2ONETRACK_taum'] )
        X1['VTXISODCHI2TWOTRACK_tau_max'] = np.maximum( x['Bp_VTXISODCHI2TWOTRACK_taup'], x['Bp_VTXISODCHI2TWOTRACK_taum'] ) 
        X1['VTXISODCHI2ONETRACK_B'] = x['Bp_VTXISODCHI2ONETRACK_B']
        X1['tau_iso_second_value_max'] = np.maximum( x['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taup'], x['Bp_B2Ksttautau_ISOBDTSECONDVALUE_taum'] ) 
        X1['tau_iso_third_value_max'] = np.maximum( x['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup'], x['Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum'] )
        X1['VTXISONUMVTX_tau_max'] = np.maximum( x['Bp_VTXISONUMVTX_taup'], x['Bp_VTXISONUMVTX_taum'] )
        X1['TRKISOBDTFIRSTVALUE_tau_pi_min_min'] = np.minimum( np.minimum( x['Bp_TRKISOBDTFIRSTVALUE_taup_pi1'], x['Bp_TRKISOBDTFIRSTVALUE_taup_pi2'], x['Bp_TRKISOBDTFIRSTVALUE_taup_pi3'] ), np.minimum( x['Bp_TRKISOBDTFIRSTVALUE_taum_pi1'], x['Bp_TRKISOBDTFIRSTVALUE_taum_pi2'], x['Bp_TRKISOBDTFIRSTVALUE_taum_pi3'] ) )
        X1['TRKISOBDTTHIRDVALUE_tau_pi_min_min'] = np.minimum( np.minimum( x['Bp_TRKISOBDTTHIRDVALUE_taup_pi1'], x['Bp_TRKISOBDTTHIRDVALUE_taup_pi2'], x['Bp_TRKISOBDTTHIRDVALUE_taup_pi3'] ), np.minimum( x['Bp_TRKISOBDTTHIRDVALUE_taum_pi1'], x['Bp_TRKISOBDTTHIRDVALUE_taum_pi2'], x['Bp_TRKISOBDTTHIRDVALUE_taum_pi3'] ) )
    else:
        X1['VTXISODCHI2ONETRACK_D_min'] = np.minimum( x['Bp_VTXISODCHI2ONETRACK_D0bar'], x['Bp_VTXISODCHI2ONETRACK_Dsp'] )
        X1['VTXISODCHI2TWOTRACK_D_max'] = np.maximum( x['Bp_VTXISODCHI2TWOTRACK_D0bar'], x['Bp_VTXISODCHI2TWOTRACK_Dsp'] ) 
        X1['VTXISODCHI2ONETRACK_B'] = x['Bp_VTXISODCHI2ONETRACK_B']
        X1['D_iso_second_value_max'] = np.maximum( x['Bp_Bstautau_ISOBDTSECONDVALUE_taup'], x['Bp_Bstautau_ISOBDTSECONDVALUE_taum'] ) 
        X1['D_iso_third_value_max'] = np.maximum( x['Bp_Bstautau_ISOBDTTHIRDVALUE_taup'], x['Bp_Bstautau_ISOBDTTHIRDVALUE_taum'] )
        X1['VTXISONUMVTX_D_max'] = np.maximum( x['Bp_VTXISONUMVTX_D0bar'], x['Bp_VTXISONUMVTX_Dsp'] )
        X1['TRKISOBDTFIRSTVALUE_D_pi_min_min'] = np.minimum( np.minimum( x['Bp_TRKISOBDTFIRSTVALUE_D0bar_K'], x['Bp_TRKISOBDTFIRSTVALUE_D0bar_pi'] ), np.minimum( x['Bp_TRKISOBDTFIRSTVALUE_Dsp_K1'], x['Bp_TRKISOBDTFIRSTVALUE_Dsp_K2'], x['Bp_TRKISOBDTFIRSTVALUE_Dsp_pi'] ) )
        X1['TRKISOBDTTHIRDVALUE_D_pi_min_min'] = np.minimum( np.minimum( x['Bp_TRKISOBDTTHIRDVALUE_D0bar_K'], x['Bp_TRKISOBDTTHIRDVALUE_D0bar_pi'] ), np.minimum( x['Bp_TRKISOBDTTHIRDVALUE_Dsp_K1'], x['Bp_TRKISOBDTTHIRDVALUE_Dsp_K2'], x['Bp_TRKISOBDTTHIRDVALUE_Dsp_pi'] ) )
            
    ############################################################################ Topology #####################################################################
    if((isKtautau == True) or (is_cocktailMC == True)):
        X2['tau_M_max'] = np.maximum( x['taup_M'], x['taum_M'] )
        X2['tau_M12_M23_max_max'] = np.maximum( np.maximum( x['taup_M12'], x['taup_M23'] ), np.maximum( x['taum_M12'], x['taum_M23'] ) )
        X2['tau_M12_M23_min_max'] = np.minimum( np.maximum( x['taup_M12'], x['taup_M23'] ), np.maximum( x['taum_M12'], x['taum_M23'] ) )
        X2['tau_M12_M23_max_min'] = np.maximum( np.minimum( x['taup_M12'], x['taup_M23'] ), np.minimum( x['taum_M12'], x['taum_M23'] ) )
        X2['log10_1_minus_tau_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(x['taup_DIRA_ORIVX'] ))*np.sign( x['taup_DIRA_ORIVX']),  np.log10(1 - np.abs(x['taum_DIRA_ORIVX'] ))*np.sign( x['taum_DIRA_ORIVX'] ) )
        X2['log10_df_chi2'] = np.log10( x['df_chi2'] )

        Cx_taup_x =  (x['df_DV1y'] - x['df_BVy'])*x['df_Kp_PZ']  - ( x['df_DV1z'] - x['df_BVz'])*x['df_Kp_PY']
        Cy_taup_x =  (x['df_DV1z'] - x['df_BVz'])*x['df_Kp_PX']  - ( x['df_DV1x'] - x['df_BVx'])*x['df_Kp_PZ']
        Cz_taup_x =  (x['df_DV1x'] - x['df_BVx'])*x['df_Kp_PY']  - ( x['df_DV1y'] - x['df_BVy'])*x['df_Kp_PX']
        C_taup_x = np.sqrt( Cx_taup_x**2 + Cy_taup_x**2 + Cz_taup_x**2  )
        IP_taup_Kp_x = (2*C_taup_x)/( np.sqrt( x['df_Kp_PX']**2 + x['df_Kp_PY']**2 + x['df_Kp_PZ']**2 ) )

        Cx_taum_x =  (x['df_DV2y'] - x['df_BVy'])*x['df_Kp_PZ']  - ( x['df_DV2z'] - x['df_BVz'])*x['df_Kp_PY']
        Cy_taum_x =  (x['df_DV2z'] - x['df_BVz'])*x['df_Kp_PX']  - ( x['df_DV2x'] - x['df_BVx'])*x['df_Kp_PZ']
        Cz_taum_x =  (x['df_DV2x'] - x['df_BVx'])*x['df_Kp_PY']  - ( x['df_DV2y'] - x['df_BVy'])*x['df_Kp_PX']
        C_taum_x = np.sqrt( Cx_taum_x**2 + Cy_taum_x**2 + Cz_taum_x**2  )
        IP_taum_Kp_x = (2*C_taum_x)/( np.sqrt( x['df_Kp_PX']**2 + x['df_Kp_PY']**2 + x['df_Kp_PZ']**2 ) )

        X2['IP_tau_Kp_max'] = np.maximum( IP_taup_Kp_x, IP_taum_Kp_x ) 
        X2['IP_tau_Kp_min'] = np.minimum( IP_taup_Kp_x, IP_taum_Kp_x ) 

        X2['Bp_M03_M04_max'] = np.maximum( x['Bp_M03'], x['Bp_M04'] )
        X2['Bp_M05_M06_max'] = np.maximum( x['Bp_M05'], x['Bp_M06'] )

        for i in range(len(X2['log10_df_chi2'])):
            if np.isnan(X2['log10_df_chi2'][i]):
                X2['log10_df_chi2'][i] = -99999
    else:
        X2['log10_1_minus_D_DIRA_BV_min'] = np.minimum( np.log10(1 - np.abs(x['D0bar_DIRA_ORIVX'] ))*np.sign( x['D0bar_DIRA_ORIVX']),  np.log10(1 - np.abs(x['Dsp_DIRA_ORIVX'] ))*np.sign( x['Dsp_DIRA_ORIVX'] ) )

        dtf_chi2_x = x['Bp_dtf_chi2']
        for i in range(len(dtf_chi2_x)):
            if np.isnan( np.log10( dtf_chi2_x[i][0])  ):
                dtf_chi2_x[i] = -99999
            else:
                dtf_chi2_x[i] = np.log10( dtf_chi2_x[i][0]) 
        X2['log10_DTF_chi2'] = dtf_chi2_x.astype(float)

    ##################################################################################################################################################################
    names = ["XGBoost"]

    classifiers_first_step = []
    classifiers_second_step = []

    if((isKtautau == True) or (is_cocktailMC == True)):
        with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/clf_isolation.pkl', 'rb') as f:
            while True:
                try:
                    classifiers_first_step.append(pickle.load(f))
                except EOFError:
                    break
        with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/Ktautau/clf_topology.pkl', 'rb') as f:
            while True:
                try:
                    classifiers_second_step.append(pickle.load(f))
                except EOFError:
                    break
    else:
        with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/clf_isolation.pkl', 'rb') as f:
            while True:
                try:
                    classifiers_first_step.append(pickle.load(f))
                except EOFError:
                    break
        with open('/panfs/felician/B2Ktautau/workflow/sklearn_training/DDs/clf_topology.pkl', 'rb') as f:
            while True:
                try:
                    classifiers_second_step.append(pickle.load(f))
                except EOFError:
                    break

    decisions_first_step = []
    for i in range(len(classifiers_first_step)):
        if( (names[i]=="RForest") or (names[i] == "XGBoost") ):
            decisions_first_step.append( classifiers_first_step[i].predict_proba(X1)[:, 1] )
        else:
            decisions_first_step.append( classifiers_first_step[i].decision_function(X1) )

    decisions_second_step = []
    for i in range(len(classifiers_second_step)):
        if( (names[i]=="RForest") or (names[i] == "XGBoost") ):
            decisions_second_step.append( classifiers_second_step[i].predict_proba(X2)[:, 1] )
        else:
            decisions_second_step.append( classifiers_second_step[i].decision_function(X2) )

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

def main(argv):

    year = argv[1]
    species = argv[2]
    line = argv[3]

    year = int(year)
    species = int(species)
    line = int(line)
    
    if((species == 1) or (species == 10) or (species == 2) or (species == 3)):
        isKtautau = True
    else:
        isKtautau = False

    if(species == 100):
        isBuDDKp_cocktail = True
    else:
        isBuDDKp_cocktail = False
    if(species == 110):
        isBdDDKp_cocktail = True
    else:
        isBdDDKp_cocktail = False
    if(species == 120):
        isBsDDKp_cocktail = True
    else:
        isBsDDKp_cocktail = False
    if(species == 130):
        isBuDDK0_cocktail = True
    else:
        isBuDDK0_cocktail = False
    if(species == 140):
        isBdDDK0_cocktail = True
    else:
        isBdDDK0_cocktail = False
    if(species == 150):
        isBuDD_cocktail = True
    else: 
        isBuDD_cocktail = False
    if(species == 160):
        isBdDD_cocktail = True
    else:
        isBdDD_cocktail = False
    if(species == 170):
        isBsDD_cocktail = True
    else: 
        isBsDD_cocktail = False

    if(isBuDDKp_cocktail or isBdDDKp_cocktail or isBsDDKp_cocktail or isBuDDK0_cocktail or isBdDDK0_cocktail or isBuDD_cocktail or isBdDD_cocktail or isBsDD_cocktail):
        is_cocktailMC = True
    else:
        is_cocktailMC = False

    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/sklearn_response/201{0}/Species_{1}/{2}.root".format(year,species,line), "RECREATE")

    fc = ROOT.TFileCollection("fc", "fc", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_{1}/pre_sel_tree.txt".format(year,species), 1, line)
    t = ROOT.TChain("DecayTree")
    t.AddFileInfoList(fc.GetList())

    if((isKtautau == True) or (is_cocktailMC == True)):
        fc1 = ROOT.TFileCollection("fc1", "fc1", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{0}/Species_{1}/fit_results.txt".format(year,species), 1, line)
        t1 = ROOT.TChain("DecayTree")
        t1.AddFileInfoList(fc1.GetList())

        fc2 = ROOT.TFileCollection("fc2", "fc2", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/201{0}/Species_{1}/invariant_mass_tree.txt".format(year,species), 1, line)
        t2 = ROOT.TChain("DecayTree")
        t2.AddFileInfoList(fc2.GetList())

        N = t.GetEntries()
        N1 = t1.GetEntries()
        N2 = t2.GetEntries()

        print(N)
        print(N1)
        print(N2)

        if((N != N1) or (N != N2)):
            print("Wrong number of entries")
            quit()
        t.AddFriend(t1)
        t.AddFriend(t2)

    df = ROOT.RDataFrame(t)

    compute_response(year, species, line, df, fout, isKtautau, is_cocktailMC, "")

    fout.Close()

if __name__ == "__main__":
    main(sys.argv)