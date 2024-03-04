
import uproot
import numpy as np
import pandas as pd
import uncertainties

def main():

    i_list = [4,5,6, 7,8,9,10, 11,12,13, 14,15,16,17]
    N = len(i_list)

    branches_m = np.array(['aaaaaaa' for i in i_list])
    branches_V = np.array([ ['aaaaaaaaaaaaa' for i in i_list] for j in i_list])

    branches_m_new = np.array(['aaaaaaaaaaaaaaaaaaaaaaa' for i in i_list])
    branches_V_new = np.array([ ['aaaaaaaaaaaaaaaaaaaaaaaaaaaaa' for i in i_list] for j in i_list])
    
    for i in range(N):
        branches_m[i] = 'df_m_{0}'.format(i_list[i])
        branches_m_new[i] = 'df_m_prime_{0}'.format(i_list[i])
        for j in range(N):
            branches_V[i,j] = 'df_V_{0}_{1}'.format(i_list[i],i_list[j])
            branches_V_new[i,j] = 'df_V_prime_{0}_{1}'.format(i_list[i],i_list[j])

    tree = uproot.open('/panfs/felician/B2Ktautau/ROOT_Sim/2018/DVntuple_MC_2018_MagUp_visible.root:ntuple/DecayTree')
    df_m = tree.arrays(branches_m, library="pd")

    V_dfs = []
    for i in range(N):
        V_dfs.append( tree.arrays(branches_V[i], library="pd") )

    m = np.zeros(N)
    V = np.zeros((N,N))

    branches_dict = dict.fromkeys(branches_m_new, "float64")
    for i in range(N):
        branches_V_dict = dict.fromkeys(branches_V_new[i], "float64")
        branches_dict.update(branches_V_dict)

    fout = uproot.recreate("/panfs/felician/B2Ktautau/ROOT_Sim/2018/transformed_m_V.root", compression=uproot.LZ4(12))
    fout.mktree("DecayTree", branches_dict)
    
    num = tree.num_entries

    for evt in range(1):
        m = df_m.iloc[evt].to_numpy()

        for i in range(N):
            V[i,:] = V_dfs[i].iloc[evt].to_numpy()
        
        # make change of variables in m
        DV1x = m[0]
        DV1y = m[1]
        DV1z = m[2]
        p3pi1x = m[3]
        p3pi1y = m[4]
        p3pi1z = m[5]
        E3pi1 = m[6]

        DV2x = m[7]
        DV2y = m[8]
        DV2z = m[9]
        p3pi2x = m[10]
        p3pi2y = m[11]
        p3pi2z = m[12]
        E3pi2 = m[13]

        # p6piKx = m[14]
        # p6piKy = m[15]
        # p6piKz = m[16]
        # E6piK = m[17]

        m3pi1_squared = E3pi1**2 - p3pi1x**2 - p3pi1y**2 - p3pi1z**2
        m3pi2_squared = E3pi2**2 - p3pi2x**2 - p3pi2y**2 - p3pi2z**2
        # m6piK_squared = E6piK**2 - p6piKx**2 - p6piKy**2 - p6piKz**2

        m[6] = m3pi1_squared
        m[13] = m3pi2_squared
        # m[17] = m6piK_squared

        # make change of variables in V
        taup_cov = np.zeros((7,7))
        taum_cov = np.zeros((7,7))
        # Bp_cov = np.zeros((4,4))

        for i in range(7):
            for j in range(7):
                taup_cov[i,j] = V[i,j]
                taum_cov[i,j] = V[i+7,j+7]

        # for i in range(4):
        #     for j in range(4):
        #         Bp_cov[i,j] = V[i+14,j+14]

        (dv1x,dv1y,dv1z,p1x,p1y,p1z,E1) = uncertainties.correlated_values([DV1x,DV1y,DV1z,p3pi1x,p3pi1y,p3pi1z,E3pi1], taup_cov)
        (dv2x,dv2y,dv2z,p2x,p2y,p2z,E2) = uncertainties.correlated_values([DV2x,DV2y,DV2z,p3pi1x,p3pi2y,p3pi2z,E3pi2], taum_cov)
        # (px,py,pz,E) = uncertainties.correlated_values([p6piKx,p6piKy,p6piKz,E6piK], Bp_cov)

        m1_squared = E1**2 - p1x**2 - p1y**2 - p1z**2
        m2_squared = E2**2 - p2x**2 - p2y**2 - p2z**2
        # m_squared = E**2 - px**2 - py**2 - pz**2

        taup_cov_prime = uncertainties.covariance_matrix([dv1x,dv1y,dv1z,p1x,p1y,p1z,m1_squared])
        taum_cov_prime = uncertainties.covariance_matrix([dv2x,dv2y,dv2z,p2x,p2y,p2z,m2_squared])
        # Bp_cov_prime = uncertainties.covariance_matrix([px,py,pz,m_squared])

        # print("3pi1 covariance")
        # print(taup_cov_prime)

        # print("3pi2 covariance")
        # print(taum_cov_prime)

        # print("6piK covariance")
        # print(Bp_cov_prime)

        np_taup_cov_prime = np.array(taup_cov_prime)
        np_taum_cov_prime = np.array(taum_cov_prime)
        # np_Bp_cov_prime = np.array(Bp_cov_prime)

        for i in range(7):
            for j in range(7):
                V[i,j] = np_taup_cov_prime[i,j]
                V[i+7,j+7] = np_taum_cov_prime[i,j]

        # for i in range(4):
        #     for j in range(4):
        #         V[i+14,j+14] = np_Bp_cov_prime[i,j]
            
        # Save TTree with updated m and V
        for i in range(N):
            branches_dict['df_m_prime_{0}'.format(i_list[i])] = [m[i]]
            for j in range(N):
                branches_dict['df_V_prime_{0}_{1}'.format(i_list[i],i_list[j])] = [V[i,j]]
        
        fout["DecayTree"].extend(branches_dict) # this fills the TTree

    return


if __name__ == "__main__":
    main()

