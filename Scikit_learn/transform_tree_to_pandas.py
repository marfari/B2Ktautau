import uproot
import pandas as pd
import numpy as np

year = '2018'
isMC = True
isWS = False
doSelection = True # variables necessary to make selections
variables = 0 # only necessary for data (for now) because too many events
# 0 -> tau- kinematic
# 1 -> tau+ kinematic
# 2 -> tau- isolation
# 3 -> tau+ isolation
# 4 -> B+ kinematic
# 5 -> B+ isolation
if doSelection:
    variables = -1

if year == '2018':
    if isMC:
        if doSelection:
            branches = ['Bp_ConsBp_seq_M', 'Bp_ConsBp_seq_status', 'Bp_L0HadronDecision_TOS', 'Bp_L0HadronDecision_TIS', 'Bp_L0MuonDecision_TIS', 'Bp_L0ElectronDecision_TIS', 'Bp_L0PhotonDecision_TIS', 'Bp_Hlt1TrackMVADecision_TOS', 'Bp_Hlt1TwoTrackMVADecision_TOS', 'Bp_Hlt2Topo2BodyDecision_TOS', 'Bp_Hlt2Topo3BodyDecision_TOS', 'Bp_Hlt2Topo4BodyDecision_TOS', 'Bp_ConsBp_seq_decayLength', 'Bp_ConsBp_seq_tauminus_0_decayLength', 'Bp_ConsBp_seq_tauminus_decayLength', 'taup_DIRA_ORIVX', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV', 'taup_ENDVERTEX_X', 'taum_ENDVERTEX_X', 'taup_ENDVERTEX_Y', 'taum_ENDVERTEX_Y', 'taup_ENDVERTEX_Z', 'taum_ENDVERTEX_Z', 'Bp_ConsBp_seq_chi2']
        elif variables == 0:
            branches = ['taum_M', 'taum_M12', 'taum_M13', 'taum_M23', 'taum_ENDVERTEX_CHI2', 'Bp_ConsBp_seq_tauminus_0_decayLength', 'Bp_ConsBp_seq_tauminus_0_decayLengthErr', 'taum_DIRA_ORIVX', 'taum_DIRA_OWNPV', 'taum_pi1_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV', 'taum_pi3_IPCHI2_OWNPV', 'Bp_ConsBp_seq_tauminus_0_PE', 'Bp_ConsBp_seq_tauminus_0_PX', 'Bp_ConsBp_seq_tauminus_0_PY', 'Bp_ConsBp_seq_tauminus_0_PZ']
        elif variables == 1:
            branches = ['taup_M', 'taup_M12', 'taup_M13', 'taup_M23', 'taup_ENDVERTEX_CHI2', 'Bp_ConsBp_seq_tauminus_decayLength', 'Bp_ConsBp_seq_tauminus_decayLengthErr', 'taup_DIRA_ORIVX', 'taup_DIRA_OWNPV', 'taup_pi1_IPCHI2_OWNPV', 'taup_pi2_IPCHI2_OWNPV', 'taup_pi3_IPCHI2_OWNPV', 'Bp_ConsBp_seq_tauminus_PE', 'Bp_ConsBp_seq_tauminus_PX', 'Bp_ConsBp_seq_tauminus_PY', 'Bp_ConsBp_seq_tauminus_PZ']
        elif variables == 2:
            branches = ['Bp_CC_IT_tau1', 'Bp_CC_MAXPT_PE_tau1', 'Bp_CC_MAXPT_PT_tau1', 'Bp_CC_MAXPT_Q_tau1', 'Bp_CC_MULT_tau1', 'Bp_CC_PZ_tau1', 'Bp_CC_PZASYM_tau1', 'Bp_CC_SPT_tau1', 'Bp_CC_VPT_tau1', 'Bp_NC_DELTAETA_tau1', 'Bp_NC_DELTAPHI_tau1', 'Bp_NC_IT_tau1', 'Bp_NC_MAXPT_PT_tau1', 'Bp_NC_MAXPT_PZ_tau1', 'Bp_NC_MULT_tau1', 'Bp_NC_PTASYM_tau1', 'Bp_NC_PZ_tau1', 'Bp_NC_SPT_tau1', 'Bp_NC_VPT_tau1', 'Bp_VTXISODCHI2ONETRACK_tau1', 'Bp_VTXISODCHI2TWOTRACK_tau1', 'Bp_VTXISONUMVTX_tau1']
        elif variables == 3: 
            branches = ['Bp_CC_IT_tau2', 'Bp_CC_MAXPT_PE_tau2', 'Bp_CC_MAXPT_PT_tau2', 'Bp_CC_MAXPT_Q_tau2', 'Bp_CC_MULT_tau2', 'Bp_CC_PZ_tau2', 'Bp_CC_PZASYM_tau2', 'Bp_CC_SPT_tau2', 'Bp_CC_VPT_tau2', 'Bp_NC_DELTAETA_tau2', 'Bp_NC_DELTAPHI_tau2', 'Bp_NC_IT_tau2', 'Bp_NC_MAXPT_PT_tau2', 'Bp_NC_MAXPT_PZ_tau2', 'Bp_NC_MULT_tau2', 'Bp_NC_PTASYM_tau2', 'Bp_NC_PZ_tau2', 'Bp_NC_SPT_tau2', 'Bp_NC_VPT_tau2', 'Bp_VTXISODCHI2ONETRACK_tau2', 'Bp_VTXISODCHI2TWOTRACK_tau2', 'Bp_VTXISONUMVTX_tau2']
        elif variables == 4:
            branches = ['Bp_ConsBp_seq_PX', 'Bp_ConsBp_seq_PY[09', 'Bp_ConsBp_seq_PZ', 'Bp_ConsBp_seq_decayLength', 'Bp_ConsBp_seq_decayLengthErr', 'Bp_ETA', 'Bp_M01', 'Bp_M02', 'Bp_M03', 'Bp_M04', 'Bp_M05', 'Bp_M06', 'Bp_M012', 'Bp_M013', 'Bp_M014', 'Bp_M015', 'Bp_M016', 'Bp_M023', 'Bp_M024', 'Bp_M025', 'Bp_M026', 'Bp_M034', 'Bp_M035', 'Bp_M036', 'Bp_M045', 'Bp_M046', 'Bp_M056']
        elif variables == 5:
            branches = ['Bp_CC_IT_B', 'Bp_CC_MAXPT_PE_B', 'Bp_CC_MAXPT_PT_B', 'Bp_CC_MAXPT_Q_B', 'Bp_CC_MULT_B', 'Bp_CC_PZ_B', 'Bp_CC_PZASYM_B', 'Bp_CC_SPT_B', 'Bp_CC_VPT_B', 'Bp_NC_DELTAETA_B', 'Bp_NC_DELTAPHI_B', 'Bp_NC_IT_B', 'Bp_NC_MAXPT_PT_B', 'Bp_NC_MAXPT_PZ_B', 'Bp_NC_MULT_B', 'Bp_NC_PTASYM_B', 'Bp_NC_PZ_B', 'Bp_NC_SPT_B', 'Bp_NC_VPT_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISONUMVTX_B']  
    else:
        if doSelection:
            branches = ['Bp_ConsBp_seq_M', 'Bp_ConsBp_seq_status', 'Bp_L0HadronDecision_TOS', 'Bp_L0HadronDecision_TIS', 'Bp_L0MuonDecision_TIS', 'Bp_L0ElectronDecision_TIS', 'Bp_L0PhotonDecision_TIS', 'Bp_Hlt1TrackMVADecision_TOS', 'Bp_Hlt1TwoTrackMVADecision_TOS', 'Bp_Hlt2Topo2BodyDecision_TOS', 'Bp_Hlt2Topo3BodyDecision_TOS', 'Bp_Hlt2Topo4BodyDecision_TOS', 'Bp_ConsBp_seq_decayLength', 'Bp_ConsBp_seq_tauminus_0_decayLength', 'Bp_ConsBp_seq_tauminus_decayLength', 'taup_DIRA_ORIVX', 'taum_DIRA_ORIVX', 'Bp_DIRA_OWNPV', 'taup_ENDVERTEX_X', 'taum_ENDVERTEX_X', 'taup_ENDVERTEX_Y', 'taum_ENDVERTEX_Y', 'taup_ENDVERTEX_Z', 'taum_ENDVERTEX_Z', 'Bp_ConsBp_seq_chi2']
        elif variables == 0:
            branches = ['taum_M', 'taum_M12', 'taum_M13', 'taum_M23', 'taum_ENDVERTEX_CHI2', 'Bp_ConsBp_seq_tauminus_0_decayLength', 'Bp_ConsBp_seq_tauminus_0_decayLengthErr', 'taum_DIRA_ORIVX', 'taum_DIRA_OWNPV', 'taum_pi1_IPCHI2_OWNPV', 'taum_pi2_IPCHI2_OWNPV', 'taum_pi3_IPCHI2_OWNPV', 'Bp_ConsBp_seq_tauminus_0_PE', 'Bp_ConsBp_seq_tauminus_0_PX', 'Bp_ConsBp_seq_tauminus_0_PY', 'Bp_ConsBp_seq_tauminus_0_PZ']
        elif variables == 1:
            branches = ['taup_M', 'taup_M12', 'taup_M13', 'taup_M23', 'taup_ENDVERTEX_CHI2', 'Bp_ConsBp_seq_tauminus_decayLength', 'Bp_ConsBp_seq_tauminus_decayLengthErr', 'taup_DIRA_ORIVX', 'taup_DIRA_OWNPV', 'taup_pi1_IPCHI2_OWNPV', 'taup_pi2_IPCHI2_OWNPV', 'taup_pi3_IPCHI2_OWNPV', 'Bp_ConsBp_seq_tauminus_PE', 'Bp_ConsBp_seq_tauminus_PX', 'Bp_ConsBp_seq_tauminus_PY', 'Bp_ConsBp_seq_tauminus_PZ']
        elif variables == 2:
            branches = ['Bp_CC_IT_tau1', 'Bp_CC_MAXPT_PE_tau1', 'Bp_CC_MAXPT_PT_tau1', 'Bp_CC_MAXPT_Q_tau1', 'Bp_CC_MULT_tau1', 'Bp_CC_PZ_tau1', 'Bp_CC_PZASYM_tau1', 'Bp_CC_SPT_tau1', 'Bp_CC_VPT_tau1', 'Bp_NC_DELTAETA_tau1', 'Bp_NC_DELTAPHI_tau1', 'Bp_NC_IT_tau1', 'Bp_NC_MAXPT_PT_tau1', 'Bp_NC_MAXPT_PZ_tau1', 'Bp_NC_MULT_tau1', 'Bp_NC_PTASYM_tau1', 'Bp_NC_PZ_tau1', 'Bp_NC_SPT_tau1', 'Bp_NC_VPT_tau1', 'Bp_VTXISODCHI2ONETRACK_tau1', 'Bp_VTXISODCHI2TWOTRACK_tau1', 'Bp_VTXISONUMVTX_tau1']
        elif variables == 3: 
            branches = ['Bp_CC_IT_tau2', 'Bp_CC_MAXPT_PE_tau2', 'Bp_CC_MAXPT_PT_tau2', 'Bp_CC_MAXPT_Q_tau2', 'Bp_CC_MULT_tau2', 'Bp_CC_PZ_tau2', 'Bp_CC_PZASYM_tau2', 'Bp_CC_SPT_tau2', 'Bp_CC_VPT_tau2', 'Bp_NC_DELTAETA_tau2', 'Bp_NC_DELTAPHI_tau2', 'Bp_NC_IT_tau2', 'Bp_NC_MAXPT_PT_tau2', 'Bp_NC_MAXPT_PZ_tau2', 'Bp_NC_MULT_tau2', 'Bp_NC_PTASYM_tau2', 'Bp_NC_PZ_tau2', 'Bp_NC_SPT_tau2', 'Bp_NC_VPT_tau2', 'Bp_VTXISODCHI2ONETRACK_tau2', 'Bp_VTXISODCHI2TWOTRACK_tau2', 'Bp_VTXISONUMVTX_tau2']
        elif variables == 4:
            branches = ['Bp_ConsBp_seq_PX', 'Bp_ConsBp_seq_PY[09', 'Bp_ConsBp_seq_PZ', 'Bp_ConsBp_seq_decayLength', 'Bp_ConsBp_seq_decayLengthErr', 'Bp_ETA', 'Bp_M01', 'Bp_M02', 'Bp_M03', 'Bp_M04', 'Bp_M05', 'Bp_M06', 'Bp_M012', 'Bp_M013', 'Bp_M014', 'Bp_M015', 'Bp_M016', 'Bp_M023', 'Bp_M024', 'Bp_M025', 'Bp_M026', 'Bp_M034', 'Bp_M035', 'Bp_M036', 'Bp_M045', 'Bp_M046', 'Bp_M056']
        elif variables == 5:
            branches = ['Bp_CC_IT_B', 'Bp_CC_MAXPT_PE_B', 'Bp_CC_MAXPT_PT_B', 'Bp_CC_MAXPT_Q_B', 'Bp_CC_MULT_B', 'Bp_CC_PZ_B', 'Bp_CC_PZASYM_B', 'Bp_CC_SPT_B', 'Bp_CC_VPT_B', 'Bp_NC_DELTAETA_B', 'Bp_NC_DELTAPHI_B', 'Bp_NC_IT_B', 'Bp_NC_MAXPT_PT_B', 'Bp_NC_MAXPT_PZ_B', 'Bp_NC_MULT_B', 'Bp_NC_PTASYM_B', 'Bp_NC_PZ_B', 'Bp_NC_SPT_B', 'Bp_NC_VPT_B', 'Bp_VTXISODCHI2ONETRACK_B', 'Bp_VTXISODCHI2TWOTRACK_B', 'Bp_VTXISONUMVTX_B']

print('Starting the transformation')
if isMC:
    tree = uproot.open('/panfs/felician/B2Ktautau/ROOT_Sim/'+year+'/mc_'+year+'_truth_matched.root:DecayTree')
    array = tree.arrays(branches, library="np")
    df = pd.DataFrame(data=array)

else:
    with open('/home/felician/B2Ktautau/Files_on_grid/data_'+year+'_seq_MagUp.txt') as file:
        files_MagUp = file.read().splitlines(True) 
    with open('/home/felician/B2Ktautau/Files_on_grid/data_'+year+'_seq_MagDown.txt') as file:
        files_MagDown = file.read().splitlines(True) 
    files = files_MagUp + files_MagDown

    num_files = len(files) # [0,len(files)]
    for i in range(num_files):
        if isWS:
            files[i] += ":ntuple_SS/DecayTree"
        else:
            files[i] += ":ntuple/DecayTree"

    i = 0
    for i in range(num_files): 
        file = files[i]
        tree = uproot.open(file)

        array = tree.arrays(branches,library="np")

        if(i==0):
            df = pd.DataFrame(data=array)
        else:
            #array | tree.arrays(branches,library="np")
            df.append(pd.DataFrame(data=array))
        i+=1
    # tree = uproot.open(files)
    # array = tree.arrays(branches,library="np")

    #df = pd.DataFrame(data=array)
print("Finished the transformation")

# Remove vectors from DTF variables
for i in range(len(df.columns)):
    if 'Bp_ConsBp' in df.columns[i]:
        df[df.columns[i]] = df[df.columns[i]].str[0]

print(df)

if isMC:    
    if doSelection:
        df.to_csv('Pandas/mc_'+year+'_forSelections.csv', index=False)
        print('Created Pandas/mc_'+year+'_forSelections.csv dataframe')
    elif variables == 0:
        df.to_csv('Pandas/mc_'+year+'_taum_kinematic.csv', index=False)
        print('Created Pandas/mc_'+year+'_taum_kinematic.csv dataframe')
    elif variables == 1:
        df.to_csv('Pandas/mc_'+year+'_taup_kinematic.csv', index=False)
        print('Created Pandas/mc_'+year+'_taup_kinematic.csv dataframe')
    elif variables == 2:
        df.to_csv('Pandas/mc_'+year+'_taum_isolation.csv', index=False)
        print('Created Pandas/mc_'+year+'_taum_isolation.csv dataframe')
    elif variables == 3:
        df.to_csv('Pandas/mc_'+year+'_taup_isolation.csv', index=False)
        print('Created Pandas/mc_'+year+'_taup_isolation.csv dataframe')
    elif variables == 4:
        df.to_csv('Pandas/mc_'+year+'_B_kinematic.csv', index=False)
        print('Created Pandas/mc_'+year+'_B_kinematic.csv dataframe')
    elif variables == 5:
        df.to_csv('Pandas/mc_'+year+'_B_isolation.csv', index=False)
        print('Created Pandas/mc_'+year+'_B_isolation.csv dataframe')
else:
    if isWS:
        if doSelection:
            df.to_csv('Pandas/data_WS_'+year+'_forSelections.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_forSelections.csv dataframe')
        elif variables == 0:
            df.to_csv('Pandas/data_WS_'+year+'_taum_kinematic.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_taum_kinematic.csv dataframe')
        elif variables == 1:
            df.to_csv('Pandas/data_WS_'+year+'_taup_kinematic.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_taup_kinematic.csv dataframe')
        elif variables == 2:
            df.to_csv('Pandas/data_WS_'+year+'_taum_isolation.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_taum_isolation.csv dataframe')
        elif variables == 3:
            df.to_csv('Pandas/data_WS_'+year+'_taup_isolation.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_taup_isolation.csv dataframe')
        elif variables == 4:
            df.to_csv('Pandas/data_WS_'+year+'_B_kinematic.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_B_kinematic.csv dataframe')
        elif variables == 5:
            df.to_csv('Pandas/data_WS_'+year+'_B_isolation.csv', index=False)
            print('Created Pandas/data_WS_'+year+'_B_isolation.csv dataframe')
    else:
        if doSelection:
            df.to_csv('Pandas/data_RS_'+year+'_forSelections.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_forSelections.csv dataframe')
        elif variables == 0:
            df.to_csv('Pandas/data_RS_'+year+'_taum_kinematic.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_taum_kinematic.csv dataframe')
        elif variables == 1:
            df.to_csv('Pandas/data_RS_'+year+'_taup_kinematic.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_taup_kinematic.csv dataframe')
        elif variables == 2:
            df.to_csv('Pandas/data_RS_'+year+'_taum_isolation.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_taum_isolation.csv dataframe')
        elif variables == 3:
            df.to_csv('Pandas/data_RS_'+year+'_taup_isolation.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_taup_isolation.csv dataframe')
        elif variables == 4:
            df.to_csv('Pandas/data_RS_'+year+'_B_kinematic.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_B_kinematic.csv dataframe')
        elif variables == 5:
            df.to_csv('Pandas/data_RS_'+year+'_B_isolation.csv', index=False)
            print('Created Pandas/data_RS_'+year+'_B_isolation.csv dataframe')