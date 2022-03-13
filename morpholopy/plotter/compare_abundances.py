import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from .stellar_abundances import load_MW_data, load_GALAH_data, \
    plot_GALAH_data, load_MW_data_with_Mg_Fe, \
    load_strontium_data_Roeder, load_strontium_data_Spite
from .plot_mass_metallicity import plot_Kirby_data, plot_Kirby_analysed, \
    plot_gallazzi, plot_gallazzi_2005, plot_Kudritzki_2016, plot_Zahid_2017


def compare_stellar_abundances(sims_data, output_name_list, output_path):

    O_Fe_all = sims_data['O_Fe']
    Fe_H_all = sims_data['Fe_H']
    Mg_Fe_all = sims_data['Mg_Fe']
    FeSNIa_Fe = sims_data['FeSNIa_Fe']
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

    count = 0
    color = ['tab:blue', 'tab:green', 'tab:orange', 'crimson', 'tab:purple']
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = FeSNIa_Fe[count:count + counter[i]]
        count += counter[i]
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("Fe(SNIa)/Fe", labelpad=2)
    plt.yscale('log')
    plt.axis([-4, 1, 1e-2,2])
    ax.tick_params(direction='in', axis='both', which='both')
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/FeSNIa_Fe_comparison.png", dpi=200)


    ##################
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

    ########################
    # Load data:
    C_Fe_all = sims_data['C_Fe']
    Ba_Fe_all = sims_data['Ba_Fe']
    Sr_Fe_all = sims_data['Sr_Fe']
    Eu_Fe_all = sims_data['Eu_Fe']
    Si_Fe_all = sims_data['Si_Fe']
    N_Fe_all = sims_data['N_Fe']
    Ne_Fe_all = sims_data['Ne_Fe']

    # make remaining plots with just GALAH data (Buder+21)
    for el in ['C', 'Si', 'Eu', 'Ba']:
        fig = plt.figure(figsize=(3.8, 3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        count = 0
        for i in range(len(output_name_list)):
            xm = Fe_H_all[count:count + counter[i]]
            if el == 'C': ym = C_Fe_all[count:count + counter[i]]
            if el == 'Ba': ym = Ba_Fe_all[count:count + counter[i]]
            if el == 'Eu': ym = Eu_Fe_all[count:count + counter[i]]
            if el == 'Si': ym = Si_Fe_all[count:count + counter[i]]
            count += counter[i]

            if i == 0: plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
            plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

        plot_GALAH_data(el, galah_edges, GALAHdata)
        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_comparison.png", dpi=200)

    for el in ['N', 'Ne']:
        fig = plt.figure(figsize=(3.8, 3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        count = 0
        for i in range(len(output_name_list)):
            xm = Fe_H_all[count:count + counter[i]]
            if el == 'N': ym = N_Fe_all[count:count + counter[i]]
            if el == 'Ne': ym = Ne_Fe_all[count:count + counter[i]]
            count += counter[i]
            plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

        #plot_GALAH_data(el, galah_edges, GALAHdata)
        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_comparison.png", dpi=200)

    FeH_Ro, SrFe_Ro = load_strontium_data_Roeder()
    FeH_Sp, SrFe_Sp = load_strontium_data_Spite()

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    count = 0
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = Sr_Fe_all[count:count + counter[i]]
        count += counter[i]
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.plot(FeH_Ro, SrFe_Ro, '+', color='crimson', ms=4, label='Roederer et al. (2014)')
    plt.plot(FeH_Sp, SrFe_Sp, 'o', color='blue', ms=4, label='Spite et al. (2018)')

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Sr/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/Sr_Fe_comparison.png", dpi=200)

def plot_Fe_H_mass_relation(Mstellar, Fe_H, counter, ylabel, output_name_list, output_file):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_Kirby_data()
    plot_Kirby_analysed()
    plot_gallazzi_2005()
    plot_Kudritzki_2016()
    plot_Zahid_2017()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = Fe_H[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(10**ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(10**ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(10**ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([1e5, 1e12, 1e-3, 1e2])
    handles, labels = plt.gca().get_legend_handles_labels()
    order = np.arange(len(handles)-2)+2
    order = np.append(order,0)
    order = np.append(order,1)

    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
               loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)


def plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe, counter, ylabel, output_name_list, output_file):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_gallazzi()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = Mg_Fe[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.axis([1e8, 1e12, -0.2, 0.6])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)

def plot_O_Fe_mass_relation(Mstellar, O_Fe, counter, ylabel, output_name_list, output_file):

    plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_gallazzi()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = O_Fe[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.axis([1e8, 1e12, 0.0, 0.6])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)

def plot_metallicity_mass_relation(Mstellar, Z, counter, ylabel, output_name_list, output_file):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_Kirby_data()
    plot_Kirby_analysed()
    plot_gallazzi_2005()
    plot_Kudritzki_2016()
    plot_Zahid_2017()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = Z[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([1e5, 1e12, 1e-3, 1e2])

    handles, labels = plt.gca().get_legend_handles_labels()
    order = np.arange(len(handles)-2)+2
    order = np.append(order,0)
    order = np.append(order,1)

    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
               loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)

def compare_mass_metallicity_relations(sim_data, output_name_list, output_path):

    Mstellar = sim_data['Mstellar']
    counter = sim_data['counter']

    ##############  Fe/H  #################

    Fe_H_median = sim_data['Fe_H_median']
    Fe_H_mass_weighted = sim_data['Fe_H_mass_weighted']
    Fe_H_ratio_weighted = sim_data['Fe_H_ratio_weighted']

    ylabel = "Stellar (median) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_median_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_median, counter, ylabel, output_name_list, output_file)

    ylabel = "Stellar (mass-weighted) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_mass_weighted_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_mass_weighted, counter, ylabel, output_name_list, output_file)

    ylabel = "Stellar (mass-weight of ratio) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_ratio_weighted_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_ratio_weighted, counter, ylabel, output_name_list, output_file)

    Fe_H_light_weighted = sim_data['Fe_H_light_weighted_r_band_with_mass']
    ylabel = "Stellar (light-weight of mass) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_light_weighted_r_band_with_mass_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_light_weighted, counter, ylabel, output_name_list, output_file)

    Fe_H_light_weighted_r_band = sim_data['Fe_H_light_weighted_r_band']
    Fe_H_light_weighted_i_band = sim_data['Fe_H_light_weighted_i_band']
    Fe_H_light_weighted_Z_band = sim_data['Fe_H_light_weighted_Z_band']

    ylabel = "Stellar (light-weight of ratio) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_light_weighted_r_band_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_light_weighted_r_band, counter, ylabel, output_name_list, output_file)

    ylabel = "Stellar (light-weighted i-band of ratio) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_light_weighted_i_band_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_light_weighted_i_band, counter, ylabel, output_name_list, output_file)

    ylabel = "Stellar (light-weighted Z-band of ratio) $10^{[\mathrm{Fe/H}]}$"
    output_file = f"{output_path}/Mstellar_Fe_H_light_weighted_Z_band_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, Fe_H_light_weighted_Z_band, counter, ylabel, output_name_list, output_file)

    ###################
    ###################

    Mg_Fe_median = sim_data['Mg_Fe_median']
    ylabel = "(median) [Mg/Fe]"
    output_file = f"{output_path}/Mstellar_Mg_Fe_median_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe_median, counter, ylabel, output_name_list, output_file)

    Mg_Fe_mass_weighted = sim_data['Mg_Fe_mass_weighted']
    ylabel = "(mass-weighted) [Mg/Fe]"
    output_file = f"{output_path}/Mstellar_Mg_Fe_mass_weighted_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe_mass_weighted, counter, ylabel, output_name_list, output_file)

    Mg_Fe_light_weighted = sim_data['Mg_Fe_light_weighted_r_band']
    ylabel = "(light-weighted r-band) [Mg/Fe]"
    output_file = f"{output_path}/Mstellar_Mg_Fe_light_weighted_r_band_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe_light_weighted, counter, ylabel, output_name_list, output_file)

    Mg_Fe_ratio_weighted = sim_data['Mg_Fe_ratio_weighted']
    ylabel = "(mass-ratio-weighted) [Mg/Fe]"
    output_file = f"{output_path}/Mstellar_Mg_Fe_ratio_weighted_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe_ratio_weighted, counter, ylabel, output_name_list, output_file)

    Mg_Fe_ratio_weighted = sim_data['Mg_Fe_light_weighted_r_band_with_mass']
    ylabel = "(light-weight r-band of mass) [Mg/Fe]"
    output_file = f"{output_path}/Mstellar_Mg_Fe_light_weighted_r_band_with_mass_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe_ratio_weighted, counter, ylabel, output_name_list, output_file)

    O_Fe_median = sim_data['O_Fe_median']
    ylabel = "(median) [O/Fe]"
    output_file = f"{output_path}/Mstellar_O_Fe_median_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, O_Fe_median, counter, ylabel, output_name_list, output_file)

    O_Fe_mass_weighted = sim_data['O_Fe_mass_weighted']
    ylabel = "(mass-weighted) [O/Fe]"
    output_file = f"{output_path}/Mstellar_O_Fe_mass_weighted_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, O_Fe_mass_weighted, counter, ylabel, output_name_list, output_file)

    O_Fe_ratio_weighted = sim_data['O_Fe_ratio_weighted']
    ylabel = "(mass-weight of ratio) [O/Fe]"
    output_file = f"{output_path}/Mstellar_O_Fe_ratio_weighted_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, O_Fe_ratio_weighted, counter, ylabel, output_name_list, output_file)

    O_Fe_light_weighted = sim_data['O_Fe_light_weighted_r_band']
    ylabel = "(light-weight r-band of ratio) [O/Fe]"
    output_file = f"{output_path}/Mstellar_O_Fe_light_weighted_r_band_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, O_Fe_light_weighted, counter, ylabel, output_name_list, output_file)

    O_Fe_ratio_weighted = sim_data['O_Fe_light_weighted_r_band_with_mass']
    ylabel = "(light-weight r-band of mass) [O/Fe]"
    output_file = f"{output_path}/Mstellar_O_Fe_light_weighted_r_band_with_mass_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, O_Fe_ratio_weighted, counter, ylabel, output_name_list, output_file)

    ###################
    ###################

    Z_mass_weighted = sim_data['Z_mass_weighted']
    ylabel = "Stellar (mass-weighted) $Z/Z_{\odot}$"
    output_file = f"{output_path}/Mstellar_Z_mass_weighted_comparison.png"
    plot_metallicity_mass_relation(Mstellar, Z_mass_weighted, counter, ylabel, output_name_list, output_file)

    Z_light_weighted_r_band = sim_data['Z_light_weighted_r_band']
    Z_light_weighted_i_band = sim_data['Z_light_weighted_i_band']
    Z_light_weighted_Z_band = sim_data['Z_light_weighted_Z_band']

    ylabel = "Stellar (light-weighted r-band) $Z/Z_{\odot}$"
    output_file = f"{output_path}/Mstellar_Z_light_weighted_r_band_comparison.png"
    plot_metallicity_mass_relation(Mstellar, Z_light_weighted_r_band, counter, ylabel, output_name_list, output_file)

    ylabel = "Stellar (light-weighted i-band) $Z/Z_{\odot}$"
    output_file = f"{output_path}/Mstellar_Z_light_weighted_i_band_comparison.png"
    plot_metallicity_mass_relation(Mstellar, Z_light_weighted_i_band, counter, ylabel, output_name_list, output_file)

    ylabel = "Stellar (light-weighted Z-band) $Z/Z_{\odot}$"
    output_file = f"{output_path}/Mstellar_Z_light_weighted_Z_band_comparison.png"
    plot_metallicity_mass_relation(Mstellar, Z_light_weighted_Z_band, counter, ylabel, output_name_list, output_file)


    # # Plot parameters
    # params = {
    #     "font.size": 12,
    #     "font.family": "Times",
    #     "text.usetex": True,
    #     "figure.figsize": (4, 3),
    #     "figure.subplot.left": 0.18,
    #     "figure.subplot.right": 0.95,
    #     "figure.subplot.bottom": 0.15,
    #     "figure.subplot.top": 0.95,
    #     "figure.subplot.wspace": 0.3,
    #     "figure.subplot.hspace": 0.3,
    #     "lines.markersize": 0.5,
    #     "lines.linewidth": 0.2,
    # }
    # rcParams.update(params)
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Fe_H_all[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(10**ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(10**ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(10**ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar $10^{[\mathrm{Fe/H}]}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Fe_H_comparison.png", dpi=200)


    #############

    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_gallazzi()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Mg_Fe_all[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("[Mg/Fe]", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.axis([1e8, 1e12, -0.2, 0.6])
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Mg_Fe_comparison.png", dpi=200)

    #############

    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_gallazzi()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = O_Fe_all[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("[O/Fe]", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.axis([1e8, 1e12, 0.0, 0.6])
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    #
    # plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_O_Fe_comparison.png", dpi=200)

    ##############
    # Mass-weighted relations
    ##############

    # O_Fe_all = sim_data['O_Fe_mass_weighted']
    # Fe_H_all = sim_data['Fe_H_mass_weighted']
    # Mg_Fe_all = sim_data['Mg_Fe_mass_weighted']
    # Z_mass_weighted = sim_data['Z_mass_weighted']

    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Fe_H_all[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(10**ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(10**ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(10**ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (mass-weighted) $10^{[\mathrm{Fe/H}]}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Fe_H_total_comparison.png", dpi=200)


    #############

    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_gallazzi()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Mg_Fe_all[count:count + counter[i]]
    #
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("[Mg/Fe] (mass-weighted)", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.axis([1e8, 1e12, -0.2, 0.6])
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Mg_Fe_total_comparison.png", dpi=200)

    #############

    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_gallazzi()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = O_Fe_all[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("[O/Fe] (mass-weighted)", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.axis([1e8, 1e12, 0.0, 0.6])
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_O_Fe_total_comparison.png", dpi=200)

    # ##############
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Z_mass_weighted[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (mass-weighted) $Z/Z_{\odot}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    #
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Z_mass_weighted_comparison.png", dpi=200)

    ##############
    # Light-weighted relations
    ##############

    # Z_light_weighted_r_band = sim_data['Z_light_weighted_r_band']
    # Z_light_weighted_i_band = sim_data['Z_light_weighted_i_band']
    # Z_light_weighted_Z_band = sim_data['Z_light_weighted_Z_band']
    # Fe_H_light_weighted_r_band = sim_data['Fe_H_light_weighted_r_band']
    # Fe_H_light_weighted_i_band = sim_data['Fe_H_light_weighted_i_band']
    # Fe_H_light_weighted_Z_band = sim_data['Fe_H_light_weighted_Z_band']
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Fe_H_light_weighted_r_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(10 ** ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(10 ** ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(10 ** ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (light-weighted r-band) $10^{[\mathrm{Fe/H}]}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Fe_H_light_weighted_r_band_comparison.png", dpi=200)
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Fe_H_light_weighted_i_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(10 ** ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(10 ** ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(10 ** ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (light-weighted i-band) $10^{[\mathrm{Fe/H}]}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles) - 2) + 2
    # order = np.append(order, 0)
    # order = np.append(order, 1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Fe_H_light_weighted_i_band_comparison.png", dpi=200)
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Fe_H_light_weighted_Z_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(10 ** ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(10 ** ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(10 ** ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (light-weighted z-band) $10^{[\mathrm{Fe/H}]}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles) - 2) + 2
    # order = np.append(order, 0)
    # order = np.append(order, 1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Fe_H_light_weighted_z_band_comparison.png", dpi=200)
    #
    # plt.figure()

    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Z_light_weighted_r_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (light-weighted r-band) $Z/Z_{\odot}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Z_light_weighted_r_band_comparison.png", dpi=200)
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Z_light_weighted_i_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (light-weighted i-band) $Z/Z_{\odot}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Z_light_weighted_i_band_comparison.png", dpi=200)
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_Kirby_data()
    # plot_Kirby_analysed()
    # plot_gallazzi_2005()
    # plot_Kudritzki_2016()
    # plot_Zahid_2017()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Z_light_weighted_Z_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("Stellar (light-weighted z-band) $Z/Z_{\odot}$", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e5, 1e12, 1e-3, 1e2])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)
    #
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Z_light_weighted_z_band_comparison.png", dpi=200)
    #

    # #########
    # O_Fe_light_weighted_r_band = sim_data['O_Fe_light_weighted_r_band']
    # Mg_Fe_light_weighted_r_band = sim_data['Mg_Fe_light_weighted_r_band']
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_gallazzi()
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = O_Fe_light_weighted_r_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("[O/Fe] (light-weighted r-band)", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e8, 1e12, 0, 0.6])
    # plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_O_Fe_light_weighted_r_band_comparison.png", dpi=200)
    #
    # ###
    #
    # plt.figure()
    #
    # # Box stellar abundance --------------------------------
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plot_gallazzi()
    #
    # count = 0
    # for i in range(len(output_name_list)):
    #     xm = Mstellar_all[count:count + counter[i]]
    #     ym = Mg_Fe_light_weighted_r_band[count:count + counter[i]]
    #     count += counter[i]
    #
    #     bins = np.arange(6, 12, 0.2)
    #     ind = np.digitize(xm, bins)
    #     ylo = [np.percentile(ym[ind == i], 16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     yhi = [np.percentile(ym[ind == i], 84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     xm = [np.median(10 ** xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
    #     plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
    #     plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])
    #
    # plt.ylabel("[Mg/Fe] (light-weighted r-band)", labelpad=2)
    # plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    # plt.xscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # plt.axis([1e8, 1e12, -0.2, 0.6])
    # plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    #
    # plt.savefig(f"{output_path}/Mstellar_Mg_Fe_light_weighted_r_band_comparison.png", dpi=200)
