import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from simulation_data import simulation_data
from .stellar_abundances import calculate_abundaces_from_MW_type_galaxies, \
    load_MW_data, load_GALAH_data, plot_GALAH_data, load_MW_data_with_Mg_Fe

def compare_stellar_abundances(config):

    bins = np.arange(-7.2, 1, 0.2)
    O_Fe_all = []
    Fe_H_all = []
    Mg_Fe_all = []
    counter = []
    output_name_list = []

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        # Fetch relevant input parameters from lists
        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]

        sim_info = simulation_data.SimInfo(directory=directory,
                                           snapshot=snapshot,
                                           catalogue=catalogue,
                                           name=sim_name,)

        output_name_list.append(sim_info.simulation_name)

        # Look for abundance ratios from COLIBRE snaps:
        ratios_MW, halo_stars = calculate_abundaces_from_MW_type_galaxies(sim_info)
        O_Fe = ratios_MW['O_Fe']
        Mg_Fe = ratios_MW['Mg_Fe']
        Fe_H = ratios_MW['Fe_H']

        ind = np.digitize(Fe_H, bins)
        xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
        ym = [np.median(O_Fe[ind == i]) for i in range(1, len(bins)) if len(O_Fe[ind == i]) > 10]
        zm = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]

        O_Fe_all = np.append(O_Fe_all, ym)
        Fe_H_all = np.append(Fe_H_all, xm)
        Mg_Fe_all = np.append(Mg_Fe_all, zm)
        counter = np.append(counter, len(xm))

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
        "figure.subplot.left": 0.1,
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
    for i in range(config.number_of_inputs):
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
    plt.savefig(f"{config.output_path}/O_Fe_comparison.png", dpi=200)

    ########################
    # Load MW data:
    FeH_MW, MgFe_MW = load_MW_data_with_Mg_Fe()

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(FeH_MW, MgFe_MW, '+', color='orange', ms=4, label='MW')
    plot_GALAH_data('Mg', galah_edges, GALAHdata)

    count = 0
    for i in range(config.number_of_inputs):
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
    plt.savefig(f"{config.output_path}/Mg_Fe_comparison.png", dpi=200)
