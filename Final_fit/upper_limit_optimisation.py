import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def branching_fraction(bdt1,bdt2, folder_name):
    results = np.load('/panfs/felician/B2Ktautau/workflow/pyhf_fit/results/{0}/cls_limit_bdt1_{1}_bdt2_{2}_seed_1000.npy'.format(folder_name,bdt1,bdt2))
    exp_limit = results[1]
    return exp_limit

def main(argv):
    folder_name = argv[1]

    bdt_cuts = np.round( np.linspace(0.95,1,21), 6 )
    N = len(bdt_cuts)

    BDT1_values = np.zeros(N*N)
    BDT2_values = np.zeros(N*N)
    branching_fraction_values = np.zeros(N*N)

    a = 0
    for i in range(N):
        for j in range(N):
            bdt1 = bdt_cuts[i]
            bdt2 = bdt_cuts[j]
                
            branching_fraction_values[a] = branching_fraction(bdt1, bdt2, folder_name)
            print("BDT1 = {0} | BDT2 = {1} | BF = {2}".format(bdt1, bdt2, branching_fraction_values[a]))
            
            BDT1_values[a] = bdt1
            BDT2_values[a] = bdt2
            a += 1
    
    min_index = np.argmin(branching_fraction_values)

    print("The BDT cuts: BDT1 > {0} and BDT2 > {1}; minimise the 90% C.L. upper limit to be \n B(B -> K tau tau ) < {2} ".format(BDT1_values[min_index], BDT2_values[min_index], branching_fraction_values[min_index]))

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
    plt.savefig("/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{0}/BF_vs_bdt1_vs_bdt2.pdf".format(folder_name))


if __name__ == "__main__":
    main(sys.argv)