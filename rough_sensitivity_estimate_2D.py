import sys
import numpy as np
import ROOT
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from uncertainties import ufloat 

def rough_sensitivity(bdt1, bdt2, t_mc, t_ws, fixed_value, bkg_config):
    N_num = t_mc.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1,bdt2))
    
    if(bkg_config == 0):
        B = 0
    else:
        N_bkg_num = t_ws.GetEntries("(BDT1 > {0}) && (BDT2 > {1})".format(bdt1, bdt2))
        N_bkg_den = t_ws.GetEntries()
        eps_bkg = 2*(N_bkg_num/N_bkg_den) # factor 2 comes from eps_rs = 2*eps_ws @ BDT1=BDT2=0.8
        B = 19006409*eps_bkg

        if(bkg_config == 2):
            yield_values = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates/yield_values_bdt1_{0}_bdt2_{1}.npy'.format(bdt1,bdt2))

            N_species = len(yield_values)
            N_components = len(yield_values[0])

            for i in range(N_species):
                for j in range(N_components):
                    B += yield_values[i,j]

    # 95% C.L. (1.64)
    # 90% C.L. (1.28)
    # S = 0.5*( 1.64**2 + np.sqrt(1.64**4 + 4*(1.64**2)*B) )
    S = 0.5*( 1.28**2 + np.sqrt(1.28**4 + 4*(1.28**2)*B) )

    branching_fraction = fixed_value*(S/N_num)
    return S, B, branching_fraction


def main(argv):
    bkg_config = argv[1] # 0 (B=0), 1 (combinatorial background only), 2 (combinatorial background + yield estimates)
    bkg_config = int(bkg_config)

    # Get fixed inputs to the branching fraction
    fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/BF_inputs.npy")
    print(fixed_value)

    # Ktautau MC post-selection tree (eps_sig)
    f_mc = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt1_0_bdt2_0.root")
    t_mc = f_mc.Get("DecayTree")

    # Efficiency ratio
    eps_norm_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_norm_value.npy")
    eps_norm_error = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_norm_error.npy")
    eps_norm = ufloat(eps_norm_value, eps_norm_error)

    eps_ktautau_fixed_value = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_ktautau_fixed_value.npy")
    eps_ktautau_fixed_error = np.load("/panfs/felician/B2Ktautau/workflow/branching_fraction_inputs/eps_ktautau_fixed_error.npy")

    N_num = t_mc.GetEntries()
    N_num = ufloat( N_num, np.sqrt(N_num) )
    eps_ktautau = N_num/eps_ktautau_fixed_value

    eff_ratio = eps_norm/eps_ktautau

    print("eps_norm = {0}",format(eps_norm*100))
    print("eps_ktautau (no BDT) = {0}".format(eps_ktautau*100))
    print("eff_ratio = {0}".format(eff_ratio))

    # Ktautau WS data tree (eps_bkg)
    f_ws = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt1_0_bdt2_0.root")
    t_ws = f_ws.Get("DecayTree")

    N = 20
    bdt_cuts = np.linspace(0.9,1,N)

    S_values = np.zeros(N*N)
    B_values = np.zeros(N*N)
    branching_fraction_values = np.zeros(N*N)

    BDT1_values = np.zeros(N*N)
    BDT2_values = np.zeros(N*N)

    a = 0
    for i in range(N):
        bdt1 = bdt_cuts[i]
        for j in range(N):
            bdt2 = bdt_cuts[j]

            results = rough_sensitivity(bdt1, bdt2, t_mc, t_ws, fixed_value, bkg_config)

            S_values[a] = results[0]
            B_values[a] = results[1]
            branching_fraction_values[a] = results[2]

            BDT1_values[a] = bdt1
            BDT2_values[a] = bdt2
            a += 1

            print("BDT1 = {0} | BDT2 = {1} | S = {2} | B = {3} | BF = {4}".format(bdt1, bdt2, results[0], results[1], results[2]))

    min_index = np.argmin(branching_fraction_values)

    print("The BDT cuts: BDT1 > {0} and BDT2 > {1}; minimise the 90% C.L. upper limit to be \n B(B -> K tau tau ) < {2} \n S = {3} and B = {4}".format(BDT1_values[min_index], BDT2_values[min_index], branching_fraction_values[min_index], S_values[min_index], B_values[min_index]))

    fig = plt.figure(figsize=(6,6))
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    cm = plt.get_cmap("coolwarm")
    ax.scatter(BDT1_values, BDT2_values, branching_fraction_values, marker='o', c=branching_fraction_values, cmap=cm)
    ax.set_xlabel('BDT1 (isolation)', rotation=-15)
    ax.set_ylabel('BDT2 (topology)', rotation=45)
    ax.set_zlabel('B($B^+ \\to K^+ \\tau^+ \\tau^-$)')
    ax.set_box_aspect(None, zoom=0.85)
    plt.tight_layout()
    plt.savefig("/panfs/felician/B2Ktautau/workflow/rough_sensitivity_2D/BF_vs_bdt1_vs_bdt2_bkg_config_{0}.pdf".format(bkg_config))

if __name__ == "__main__":
    main(sys.argv)


