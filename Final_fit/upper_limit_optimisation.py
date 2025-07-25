import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

def main(argv):
    fit_type = argv[1]
    fit_name = argv[2]
    BF_sig = argv[3]

    bdt_cuts = np.round( np.linspace(0.99,1,200), 5)
    N = len(bdt_cuts)
    branching_fraction_values = np.zeros(N)

    for i in range(N):
        bdt = bdt_cuts[i]
        try:
            branching_fraction_values[i] = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BDT_{bdt}/cls_limit.npy')
        except:
            branching_fraction_values[i] = np.inf
        print(f"BDT = {bdt} | BF = {branching_fraction_values[i]}")
            
    min_index = np.argmin(branching_fraction_values)

    print(f"The BDT cut, BDT > {bdt_cuts[min_index]}, minimise the 90% C.L. upper limit to be \n B(B -> K tau tau ) < {branching_fraction_values[min_index]} ")

    plt.plot(bdt_cuts, branching_fraction_values, marker='.')
    plt.xlabel('BDT > x')
    plt.ylabel('90% C.L. upper limit on B($B^+ \\to K^+ \\tau^+ \\tau^-$)')
    plt.tight_layout()
    plt.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BF_vs_bdt.pdf")

    # fig = plt.figure(figsize=(6,6))
    # ax = Axes3D(fig, auto_add_to_figure=False)
    # fig.add_axes(ax)
    # cm = plt.get_cmap("coolwarm")
    # ax.scatter(BDT1_values, BDT2_values, branching_fraction_values, marker='o', c=branching_fraction_values, cmap=cm)
    # ax.set_xlabel('BDT1 (physics)', rotation=-15)
    # ax.set_ylabel('BDT2 (combinatorial)', rotation=45)
    # ax.set_zlabel('B($B^+ \\to K^+ \\tau^+ \\tau^-$)')
    # ax.set_box_aspect(None, zoom=0.85)
    # plt.tight_layout()
    # fig.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BF_vs_bdt1_vs_bdt2.pdf")
    # plt.clf()

    # BDT1_grid = BDT1_values.reshape(N,N)
    # BDT2_grid = BDT2_values.reshape(N,N)
    # BF_grid = branching_fraction_values.reshape(N, N)
    # fig1 = plt.figure(figsize=(8,6))
    # ax1 = fig1.add_subplot(1, 1, 1) 
    # p = ax1.pcolor(BDT1_grid, BDT2_grid, BF_grid, cmap=cm, shading='auto')
    # cbar = plt.colorbar(p, ax=ax1)
    # cbar.set_label('Branching fraction')  # optional: add a label to the legend
    # ax1.plot(BDT1_values[min_index], BDT2_values[min_index], marker='*', color='black', markersize=15, linestyle='None', label='Minimum')
    # plt.xlabel('physics BDT')
    # plt.ylabel('combinatorial BDT')
    # plt.title('Branching fraction vs BDT cuts')
    # plt.legend()
    # plt.tight_layout()
    # fig1.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_type}_{fit_name}/BF_sig_{BF_sig}/BF_vs_bdt1_vs_bdt2_1.pdf")

if __name__ == "__main__":
    main(sys.argv)