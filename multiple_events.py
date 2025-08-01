import sys
import ROOT
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import math as m
from array import array
import random

def best_candidate(species, isKtautau, is_cokctailMC, t, evt, N, indices):   
    # bdt_sum = np.zeros(N)
    dtf_chi2 = np.zeros(N)

    gsl_chi2 = np.zeros(N)

    for cand in range(N): # loop over candidates and choose the one with the largest BDT value
        t.GetEntry(indices[cand])

        if(isKtautau or is_cokctailMC):
            # bdt1 = getattr(t, "BDT1")
            # bdt2 = getattr(t, "BDT2")
            # bdt_sum[cand] = bdt1 + bdt2 

            gsl_chi2[cand] = getattr(t, "df_chi2")

        elif((species == 72) or (species == 8)):
            chi2 = getattr(t, "Bp_dtf_chi2")
            dtf_chi2[cand] = chi2[0]

    if(isKtautau or is_cokctailMC):
        # return indices[np.argmax( bdt_sum )]
        return indices[np.argmin( gsl_chi2 )]

    elif((species == 72) or (species == 8)):
        return indices[np.argmin( dtf_chi2 )]


def candidate_category(t, evt, indices, nCand):
    N = nCand[evt]

    cand_cat = np.zeros(N, dtype=int)

    for cand in range(N):
        t.GetEntry(indices[evt][cand])

        taup_pi1_PX = getattr(t, "taup_pi1_PX")
        taup_pi1_PY = getattr(t, "taup_pi1_PY")
        taup_pi1_PZ = getattr(t, "taup_pi1_PZ")
        taup_pi1_PE = getattr(t, "taup_pi1_PE")

        taup_pi2_PX = getattr(t, "taup_pi2_PX")
        taup_pi2_PY = getattr(t, "taup_pi2_PY")
        taup_pi2_PZ = getattr(t, "taup_pi2_PZ")
        taup_pi2_PE = getattr(t, "taup_pi2_PE")

        taup_pi3_PX = getattr(t, "taup_pi3_PX")
        taup_pi3_PY = getattr(t, "taup_pi3_PY")
        taup_pi3_PZ = getattr(t, "taup_pi3_PZ")
        taup_pi3_PE = getattr(t, "taup_pi3_PE")

        taum_pi1_PX = getattr(t, "taum_pi1_PX")
        taum_pi1_PY = getattr(t, "taum_pi1_PY")
        taum_pi1_PZ = getattr(t, "taum_pi1_PZ")
        taum_pi1_PE = getattr(t, "taum_pi1_PE")

        taum_pi2_PX = getattr(t, "taum_pi2_PX")
        taum_pi2_PY = getattr(t, "taum_pi2_PY")
        taum_pi2_PZ = getattr(t, "taum_pi2_PZ")
        taum_pi2_PE = getattr(t, "taum_pi2_PE")

        taum_pi3_PX = getattr(t, "taum_pi3_PX")
        taum_pi3_PY = getattr(t, "taum_pi3_PY")
        taum_pi3_PZ = getattr(t, "taum_pi3_PZ")
        taum_pi3_PE = getattr(t, "taum_pi3_PE")

        Kp_PX = getattr(t, "Kp_PX")
        Kp_PY = getattr(t, "Kp_PY")
        Kp_PZ = getattr(t, "Kp_PZ")
        Kp_PE = getattr(t, "Kp_PE")

        P1 = ROOT.Math.PxPyPzEVector( taup_pi1_PX, taup_pi1_PY, taup_pi1_PZ, taup_pi1_PE )
        P2 = ROOT.Math.PxPyPzEVector( taup_pi2_PX, taup_pi2_PY, taup_pi2_PZ, taup_pi2_PE )
        P3 = ROOT.Math.PxPyPzEVector( taup_pi3_PX, taup_pi3_PY, taup_pi3_PZ, taup_pi3_PE )

        P4 = ROOT.Math.PxPyPzEVector( taum_pi1_PX, taum_pi1_PY, taum_pi1_PZ, taum_pi1_PE )
        P5 = ROOT.Math.PxPyPzEVector( taum_pi2_PX, taum_pi2_PY, taum_pi2_PZ, taum_pi2_PE )
        P6 = ROOT.Math.PxPyPzEVector( taum_pi3_PX, taum_pi3_PY, taum_pi3_PZ, taum_pi3_PE )

        PK = ROOT.Math.PxPyPzEVector( Kp_PX, Kp_PY, Kp_PZ, Kp_PE )

    return 0


def angle_between_two_tracks(t, i, j, isKtautau, is_cokctailMC, species):
    if(isKtautau or is_cokctailMC):
        taup_pi1_PX = getattr(t, "taup_pi1_PX")
        taup_pi1_PY = getattr(t, "taup_pi1_PY")
        taup_pi1_PZ = getattr(t, "taup_pi1_PZ")

        taup_pi2_PX = getattr(t, "taup_pi2_PX")
        taup_pi2_PY = getattr(t, "taup_pi2_PY")
        taup_pi2_PZ = getattr(t, "taup_pi2_PZ")

        taup_pi3_PX = getattr(t, "taup_pi3_PX")
        taup_pi3_PY = getattr(t, "taup_pi3_PY")
        taup_pi3_PZ = getattr(t, "taup_pi3_PZ")

        taum_pi1_PX = getattr(t, "taum_pi1_PX")
        taum_pi1_PY = getattr(t, "taum_pi1_PY")
        taum_pi1_PZ = getattr(t, "taum_pi1_PZ")

        taum_pi2_PX = getattr(t, "taum_pi2_PX")
        taum_pi2_PY = getattr(t, "taum_pi2_PY")
        taum_pi2_PZ = getattr(t, "taum_pi2_PZ")

        taum_pi3_PX = getattr(t, "taum_pi3_PX")
        taum_pi3_PY = getattr(t, "taum_pi3_PY")
        taum_pi3_PZ = getattr(t, "taum_pi3_PZ")

        Kp_PX = getattr(t, "Kp_PX")
        Kp_PY = getattr(t, "Kp_PY")
        Kp_PZ = getattr(t, "Kp_PZ")

        P1 = ROOT.Math.XYZVector( taup_pi1_PX, taup_pi1_PY, taup_pi1_PZ )
        P2 = ROOT.Math.XYZVector( taup_pi2_PX, taup_pi2_PY, taup_pi2_PZ )
        P3 = ROOT.Math.XYZVector( taup_pi3_PX, taup_pi3_PY, taup_pi3_PZ )

        P4 = ROOT.Math.XYZVector( taum_pi1_PX, taum_pi1_PY, taum_pi1_PZ )
        P5 = ROOT.Math.XYZVector( taum_pi2_PX, taum_pi2_PY, taum_pi2_PZ )
        P6 = ROOT.Math.XYZVector( taum_pi3_PX, taum_pi3_PY, taum_pi3_PZ )

        PK = ROOT.Math.XYZVector( Kp_PX, Kp_PY, Kp_PZ )

        mom = [P1, P2, P3, P4, P5, P6, PK]
    elif((species == 72) or (species == 8)):
        D0bar_K_PX = getattr(t, "D0bar_K_PX")
        D0bar_K_PY = getattr(t, "D0bar_K_PY")
        D0bar_K_PZ = getattr(t, "D0bar_K_PZ")

        D0bar_pi_PX = getattr(t, "D0bar_pi_PX")
        D0bar_pi_PY = getattr(t, "D0bar_pi_PY")
        D0bar_pi_PZ = getattr(t, "D0bar_pi_PZ")

        Dsp_K1_PX = getattr(t, "Dsp_K1_PX")
        Dsp_K1_PY = getattr(t, "Dsp_K1_PY")
        Dsp_K1_PZ = getattr(t, "Dsp_K1_PZ")

        Dsp_K2_PX = getattr(t, "Dsp_K2_PX")
        Dsp_K2_PY = getattr(t, "Dsp_K2_PY")
        Dsp_K2_PZ = getattr(t, "Dsp_K2_PZ")

        Dsp_pi_PX = getattr(t, "Dsp_pi_PX")
        Dsp_pi_PY = getattr(t, "Dsp_pi_PY")
        Dsp_pi_PZ = getattr(t, "Dsp_pi_PZ")

        P1 = ROOT.Math.XYZVector( D0bar_K_PX, D0bar_K_PY, D0bar_K_PZ )
        P2 = ROOT.Math.XYZVector( D0bar_pi_PX, D0bar_pi_PY, D0bar_pi_PZ )
        P3 = ROOT.Math.XYZVector( Dsp_K1_PX, Dsp_K1_PX, Dsp_K1_PX )
        P4 = ROOT.Math.XYZVector( Dsp_K2_PX, Dsp_K2_PX, Dsp_K2_PX )
        P5 = ROOT.Math.XYZVector( Dsp_pi_PX, Dsp_pi_PY, Dsp_pi_PZ )

        mom = [P1, P2, P3, P4, P5]

    theta = np.arccos( mom[i].Dot(mom[j])/( np.sqrt(mom[i].Mag2())*np.sqrt(mom[j].Mag2()) )  )
    
    return theta

def main(argv):
    year = argv[1]
    species = argv[2]
    line = argv[3]

    year = int(year)
    species = int(species) 
    line = int(line)

    is_cokctailMC = False
    if((species == 100) or (species == 110) or (species == 120) or (species == 130) or (species == 150)):
        is_cokctailMC = True
    
    isKtautau = False
    if((species == 1) or (species == 10) or (species == 2) or (species == 3)):
        isKtautau = True

    fc = ROOT.TFileCollection("fc", "fc", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_{1}/pre_sel_tree.txt".format(year,species), 1, line)
    t = ROOT.TChain("DecayTree")
    t.AddFileInfoList(fc.GetList())
    n_entries = t.GetEntries()

    if(isKtautau or is_cokctailMC):
        fc1 = ROOT.TFileCollection("fc1", "fc1", "/panfs/felician/B2Ktautau/workflow/sklearn_response/201{0}/Species_{1}/bdt_output.txt".format(year,species), 1, line)
        t1 = ROOT.TChain("DecayTree")
        t1.AddFileInfoList(fc1.GetList())
        n1_entries = t1.GetEntries()

        fc2 = ROOT.TFileCollection("fc1", "fc1", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{0}/Species_{1}/fit_results.txt".format(year,species), 1, line)
        t2 = ROOT.TChain("DecayTree")
        t2.AddFileInfoList(fc2.GetList())
        n2_entries = t2.GetEntries()

        if((n_entries == n1_entries) and (n_entries == n2_entries)):
            t.AddFriend(t1)
            t.AddFriend(t2)
        else:
            print("Wrong number of entries")
            quit()

    df = ROOT.RDataFrame(t)

    branch_names = ['runNumber', 'eventNumber']
    x = df.AsNumpy(branch_names)
    x_cand = [ (x['eventNumber'][i], x['runNumber'][i]) for i in range(n_entries)]

    nCand_per_event = Counter(x_cand)
    nCand = [nCand_per_event[evt] for evt in nCand_per_event] # number of multiple candidates per event
    n_events = len(nCand)

    # Loop over events
    best_candidates_index = -np.ones(n_events, dtype=int)
    # second_best_candidate_index = -np.ones(n_events, dtype=int)

    n_current = 0
    for evt in range(n_events):
        N = nCand[evt]

        indices = np.zeros(N, dtype=int)
        # bdt_sum = np.zeros(N)
        # gsl_chi2 = np.zeros(N)
        # tau_chi2_sum = np.zeros(N)
        # df_merr = np.zeros(N)
        # dtf_chi2 = np.zeros(N)

        # best candidate selection
        for i in range(N):
            indices[i] = i + n_current
        
        #     t.GetEntry(i + n_current)
        #     if(isKtautau or is_cokctailMC):
        #         bdt_sum[i] = getattr(t, "BDT1") + getattr(t, "BDT2") 
        #         gsl_chi2[i] = getattr(t, "df_chi2")
        #         tau_chi2_sum[i] = getattr(t, "taup_ENDVERTEX_CHI2") + getattr(t, "taum_ENDVERTEX_CHI2")
        #         df_merr[i] = getattr(t, "df_Bp_MERR")
        #     if((species == 72) or (species == 8)):
        #         chi2 = getattr(t, "Bp_dtf_chi2")
        #         dtf_chi2[i] = chi2[0]
                
        # if(isKtautau or is_cokctailMC):
            # bdt_sum_sorted_indices = bdt_sum.argsort()
            # best_candidates_index[evt] = indices[ bdt_sum_sorted_indices[-1] ]
            # if(N > 1):
            #     second_best_candidate_index[evt] = indices[ bdt_sum_sorted_indices[-2] ]

            # gsl_chi2_sorted_indices = gsl_chi2.argsort()
            # best_candidates_index[evt] = indices[ gsl_chi2_sorted_indices[0] ]
            # if(N > 1):
            #     second_best_candidate_index[evt] = indices[ gsl_chi2_sorted_indices[1] ]

            # tau_chi2_sorted_indices = tau_chi2_sum.argsort()
            # best_candidates_index[evt] = indices[ tau_chi2_sorted_indices[0] ]
            # if(N > 1):
            #     second_best_candidate_index[evt] = indices[ tau_chi2_sorted_indices[1] ]

            # df_merr_sorted_indices = df_merr.argsort()
            # best_candidates_index[evt] = indices[ df_merr_sorted_indices[0] ]
            # if(N > 1):
            #     second_best_candidate_index[evt] = indices[ df_merr_sorted_indices[1] ]
        # elif((species == 72) or (species == 8)):
        #     best_candidates_index[evt] = indices[np.argmin( dtf_chi2 )]

        random.seed(42)
        idx = random.randint(0, N-1)
        best_candidates_index[evt] = indices[ idx ]

        # print("N = ", N, " | best index = ", idx, " | event index = ", indices[ idx ])
        
        n_current += N

    print("Creating TTree")
    fout = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/multiple_events/201{0}/Species_{1}/{2}.root".format(year,species,line), "RECREATE")
    fout.cd()
    tout = ROOT.TTree("DecayTree", "DecayTree")
    isBest = array('i', [0])
    isSecondBest = array('i', [0])
    nCandidates = array('d', [0])
    tout.Branch("is_best_cand", isBest, "is_best_cand/O")
    # tout.Branch("is_second_best_cand", isSecondBest, "is_second_best_cand/O")
    tout.Branch("n_candidates", nCandidates, "n_candidates/D")

    evt = 0
    for i in range(n_entries):
        if(i in best_candidates_index):
            isBest[0] = True
            nCandidates[0] = nCand[evt]
            evt += 1
        else:
            isBest[0] = False
            nCandidates[0] = -1

        # if(i in second_best_candidate_index):
        #     isSecondBest[0] = True
        # else:
        #     isSecondBest[0] = False

        tout.Fill()

    tout.Write()
    fout.Close()
    print("Finished creating TTree")

if __name__ == "__main__":
    main(sys.argv)