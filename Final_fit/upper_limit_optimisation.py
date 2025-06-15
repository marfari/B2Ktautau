import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

def main(argv):
    fit_name = argv[1]
    nsig = argv[2]

    bdt_cuts = np.round( np.linspace(0.965,0.985,21), 3 )
    N = len(bdt_cuts)

    BDT1_values = np.zeros(N*N)
    BDT2_values = np.zeros(N*N)
    branching_fraction_values = np.zeros(N*N)

    a = 0
    for i in range(N):
        for j in range(N):
            bdt1 = bdt_cuts[i]
            bdt2 = bdt_cuts[j]
                
            branching_fraction_values[a] = np.load(f'/panfs/felician/B2Ktautau/workflow/pyhf_fit/{fit_name}/N_sig_{nsig}/fit_results/cls_limit_bdt1_{bdt1}_bdt2_{bdt2}.npy')

            print(f"BDT1 = {bdt1} | BDT2 = {bdt2} | BF = {branching_fraction_values[a]}")
            
            BDT1_values[a] = bdt1
            BDT2_values[a] = bdt2
            a += 1
    
    min_index = np.argmin(branching_fraction_values)

    print(f"The BDT cuts: BDT1 > {BDT1_values[min_index]} and BDT2 > {BDT2_values[min_index]}; minimise the 90% C.L. upper limit to be \n B(B -> K tau tau ) < {branching_fraction_values[min_index]} ")

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
    fig.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_name}/N_sig_{nsig}/BF_vs_bdt1_vs_bdt2.pdf")
    plt.clf()

    BDT1_grid = BDT1_values.reshape(N,N)
    BDT2_grid = BDT2_values.reshape(N,N)
    BF_grid = branching_fraction_values.reshape(N, N)
    fig1 = plt.figure(figsize=(8,6))
    ax1 = fig1.add_subplot(1, 1, 1) 
    p = ax1.pcolor(BDT1_grid, BDT2_grid, BF_grid, cmap=cm, shading='auto')
    cbar = plt.colorbar(p, ax=ax1)
    cbar.set_label('Branching fraction')  # optional: add a label to the legend
    ax1.plot(BDT1_values[min_index], BDT2_values[min_index], marker='*', color='black', markersize=15, linestyle='None', label='Minimum')
    plt.xlabel('Isolation BDT')
    plt.ylabel('Topology BDT')
    plt.title('Branching fraction vs BDT cuts')
    plt.legend()
    plt.tight_layout()
    fig1.savefig(f"/panfs/felician/B2Ktautau/workflow/upper_limit_optimisation/{fit_name}/N_sig_{nsig}/BF_vs_bdt1_vs_bdt2_1.pdf")

if __name__ == "__main__":
    main(sys.argv)