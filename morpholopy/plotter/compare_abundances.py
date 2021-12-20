import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from .stellar_abundances import load_MW_data, load_GALAH_data, plot_GALAH_data, load_MW_data_with_Mg_Fe

def compare_stellar_abundances(sims_data, output_name_list, output_path):

    O_Fe_all = sims_data['O_Fe']
    Fe_H_all = sims_data['Fe_H']
    Mg_Fe_all = sims_data['Mg_Fe']
    counter = sims_data['counter']

    # Load MW data:
    FeH_MW, OFe_MW = load_MW_data()

    # Load MW data:
    GALAHdata = load_GALAH_data()
    galah_edges = np.array(GALAHdata["abundance_bin_edges"])

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)
    fig = plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(FeH_MW, OFe_MW, '+', color='tab:blue', ms=4, label='MW')
    plot_GALAH_data('O', galah_edges, GALAHdata)

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count+counter[i]]
        ym = O_Fe_all[count:count+counter[i]]
        count += counter[i]

        if i==0 :plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/O_Fe_comparison.png", dpi=200)

    ########################
    # Load MW data:
    FeH_MW, MgFe_MW = load_MW_data_with_Mg_Fe()

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(FeH_MW, MgFe_MW, '+', color='orange', ms=4, label='MW')
    plot_GALAH_data('Mg', galah_edges, GALAHdata)

    count = 0
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = Mg_Fe_all[count:count + counter[i]]
        count += counter[i]

        if i == 0: plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/Mg_Fe_comparison.png", dpi=200)
