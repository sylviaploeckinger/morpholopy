from pylab import *
import numpy as np
from .loadObservationalData import read_obs_data
from plotter.KS_relation import median_relations, Krumholz_eq39

def KS_relation_plots(siminfo,index,name_list):

    # read the observational data for the KS relations
    observational_data = read_obs_data("./plotter/obs_data")

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 2,
        "lines.linewidth": 1.0,
    }
    rcParams.update(params)
    methods = ['grid','radii']
    for method in methods:

        for mode in range(2):
            figure()
            ax = plt.subplot(1, 1, 1)

            Sigma_g = np.logspace(-1, 4, 1000)
            Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)
            plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), color="red",
                     label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$",
                     linestyle="--")

            Sigma_g = np.logspace(-1, 4, 1000)
            Sigma_star = KS(Sigma_g, 1.06, 2.511e-4)
            plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), lw=1, color="green",
                     label=r"2.51e-4 $\times$ $\Sigma_{g}^{1.06}$ (Pessa+ 2021)",
                     linestyle="-")

            # load the observational data
            if mode == 0:

                for ind, observation in enumerate(observational_data):
                    if observation.gas_surface_density is not None:
                        if (observation.description == "Bigiel et al. (2008) inner"):
                            data = observation.bin_data_KS(np.arange(-1, 3, .25), 0.4)
                            plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="v", ms=6,
                                         label=observation.description, color='darkred')
                        elif (observation.description == "Bigiel et al. (2010) outer"):
                            data2 = observation.bin_data_KS(np.arange(-1, 3, .25), 0.4)
                            plt.errorbar(data2[0], data2[1], yerr=[data2[2], data2[3]], fmt="o", ms=6,
                                         label=observation.description, color='darkred')
                plt.xlabel("log $\\Sigma_{HI}+ \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")

            elif mode == 1:

                for ind, observation in enumerate(observational_data):
                    if observation.gas_surface_density is not None:
                        if (observation.description == "Bigiel et al. (2008) inner"):
                            data = observation.bin_data_KS_molecular(np.arange(-1, 3, .25), 0.4)
                            plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="<", ms=6,
                                         label=observation.description, color='darkred')
                plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")

            color = ['tab:blue','tab:orange']
            i = 0
            for name in name_list:
                if mode == 0:
                    data = np.loadtxt(f"{siminfo.output_path}/KS_relation_best_"+method+"_%i_"%(index)+name+".txt")
                elif mode == 1:
                    data = np.loadtxt(f"{siminfo.output_path}/KS_molecular_relation_"+method+"_%i_" % (index)+name+".txt")

                surface_density = data[:,0]
                SFR_surface_density = data[:,1]

                # Get median lines
                median_surface_density, median_SFR_surface_density, \
                SFR_surface_density_err_down, SFR_surface_density_err_up = \
                    median_relations(surface_density, SFR_surface_density)

                plt.plot(surface_density, SFR_surface_density, 'o', color=color[i])
                plt.plot(median_surface_density, median_SFR_surface_density, '-', lw=2, color='white')
                plt.plot(median_surface_density, median_SFR_surface_density, '-',lw=1.2,
                         color=color[i],label=name)
                plt.fill_between(median_surface_density, SFR_surface_density_err_down,
                                 SFR_surface_density_err_up, alpha=0.2, color=color[i])


                i += 1

            plt.legend(loc='upper left',labelspacing=0.2, handlelength=1, handletextpad=0.2, frameon=False)
            ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
            plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
            plt.xlim(-1.0, 3.0)
            plt.ylim(-6.5, 1.0)
            if mode == 0:
                plt.savefig(f"{siminfo.output_path}/KS_relation_best_"+method+"_%i.png"%(index), dpi=200)
                plt.close()
            elif mode == 1:
                plt.savefig(f"{siminfo.output_path}/KS_molecular_relation_"+method+"_%i.png"%(index), dpi=200)
                plt.close()


def depletion_time_plots(siminfo,index,name_list):

    # read the observational data for the KS relations
    observational_data = read_obs_data("./plotter/obs_data")

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 2,
        "lines.linewidth": 1.0,
    }
    rcParams.update(params)
    methods = ['grid','radii']
    for method in methods:

        for mode in range(2):
            figure()
            ax = plt.subplot(1, 1, 1)

            Sigma_g = np.logspace(-1, 4, 1000)
            Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)
            plt.plot(np.log10(Sigma_g), np.log10(Sigma_g) - np.log10(Sigma_star) + 6., color="red",
                     label="KS law (Kennicutt 98)", linestyle="--")

            # load the observational data
            if mode == 0:
                for ind, observation in enumerate(observational_data):
                    if observation.gas_surface_density is not None:
                        if (observation.description == "Bigiel et al. (2008) inner"):
                            data = observation.bin_data_gas_depletion(np.arange(-1, 3, .25), 0.4)
                            plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt=">", ms=6,
                                         label=observation.description, color='darkred')
                        elif (observation.description == "Bigiel et al. (2010) outer"):
                            data2 = observation.bin_data_gas_depletion(np.arange(-1, 3, .25), 0.4)
                            plt.errorbar(data2[0], data2[1], yerr=[data2[2], data2[3]], fmt="o", ms=6,
                                         label=observation.description, color='darkred')
                plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
            elif mode == 1:
                for ind, observation in enumerate(observational_data):
                    if observation.gas_surface_density is not None:
                        if (observation.description == "Bigiel et al. (2008) inner"):
                            data = observation.bin_data_gas_depletion_molecular(np.arange(-1, 3, .25), 0.4)
                            plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="o", ms=6,
                                         label=observation.description, color='darkred')
                plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")

            color = ['tab:blue','tab:orange']
            i = 0
            for name in name_list:
                if mode == 0:
                    data = np.loadtxt(f"{siminfo.output_path}/gas_depletion_timescale_best_"+method+"_%i_" % (index)+name+".txt")
                elif mode == 1:
                    data = np.loadtxt(f"{siminfo.output_path}/molecular_gas_depletion_timescale_"+method+"_%i_"%(index)+name+".txt")

                surface_density = data[:,0]
                SFR_surface_density = data[:,1]

                # Get median lines
                median_surface_density, median_SFR_surface_density, \
                SFR_surface_density_err_down, SFR_surface_density_err_up = \
                    median_relations(surface_density, SFR_surface_density)

                plt.plot(surface_density, SFR_surface_density, 'o', color=color[i], alpha=0.6)
                plt.plot(median_surface_density, median_SFR_surface_density, '-', lw=2, color='white')
                plt.plot(median_surface_density, median_SFR_surface_density, '-',
                         color=color[i],label=name)
                plt.fill_between(median_surface_density, SFR_surface_density_err_down,
                                 SFR_surface_density_err_up, alpha=0.2, color=color[i])
                i += 1

            plt.legend(labelspacing=0.2, handlelength=2, handletextpad=0.4, frameon=False)
            ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
            plt.ylabel("log $\\rm t_{gas} = \\Sigma_{H_2} / \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$")
            plt.xlim(-1, 3.0)
            plt.ylim(7, 12)
            if mode == 0:
                plt.savefig(f"{siminfo.output_path}/gas_depletion_timescale_best_" + method + "_%i.png"%(index), dpi=200)
                plt.close()
            elif mode == 1:
                plt.savefig(f"{siminfo.output_path}/molecular_gas_depletion_timescale_"+method+"_%i.png"%(index), dpi=200)
                plt.close()

def surface_ratios_plots(siminfo, index, name_list):

    # Load data from Schruba +2021
    SchrubaData = np.loadtxt("./plotter/obs_data/Schruba2011_data.txt", usecols=(4, 5, 6))
    nonan = np.logical_and(np.isnan(SchrubaData[:, 0]) == False, np.isnan(SchrubaData[:, 1]) == False)
    Schruba_H1 = SchrubaData[nonan, 0]  # HI surface density [Msol / pc-2]
    Schruba_H2 = SchrubaData[nonan, 1]  # H2 surface density [Msol / pc-2]
    flags = SchrubaData[nonan, 2]
    select_flags = np.logical_or(flags == 1, flags == 0)  # Disregarding upper limits
    Schruba_H1 = Schruba_H1[select_flags]
    Schruba_H2 = Schruba_H2[select_flags]

    x_Schruba = np.log10(Schruba_H1 + Schruba_H2)
    y_Schruba = np.log10(Schruba_H2 / (Schruba_H1 + Schruba_H2))

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 2,
        "lines.linewidth": 1.0,
    }
    rcParams.update(params)

    methods = ['grid','radii']
    for method in methods:
        figure()
        ax = plt.subplot(1, 1, 1)

        # Krumholz 2009 lines
        Sigma_neutral = np.arange(-1, 3, 0.2)
        RH2 = 1. / Krumholz_eq39(10 ** Sigma_neutral, 0.5)
        FH2 = np.log10(1. / (1. + RH2))
        plt.plot(Sigma_neutral, FH2, '--', color='darkred', label="Krumholz+ (2009): f = 0.5")
        RH2 = 1. / Krumholz_eq39(10 ** Sigma_neutral, 0.1)
        FH2 = np.log10(1. / (1. + RH2))
        plt.plot(Sigma_neutral, FH2, ':', color='tab:red', label="Krumholz+ (2009): f = 0.1")
        plt.plot(x_Schruba, y_Schruba, 'o', color='darkred', label="Schruba+ (2011)")

        color = ['tab:blue', 'tab:orange']
        i = 0

        for name in name_list:
            if method == 'grid':
                data = np.loadtxt(f"{siminfo.output_path}/Surface_density_ratio_"+method+"_%i_"%(index)+name+".txt")
                Sigma_gas = data[:,0]
                Sigma_ratio = data[:,1]
                Median_Sigma_gas, Median_Sigma_ratio, Sigma_ratio_err_down, \
                Sigma_ratio_err_up = median_relations(Sigma_gas, Sigma_ratio)
                plt.plot(Sigma_gas, Sigma_ratio, 'o', color=color[i], alpha=0.6)
                plt.plot(Median_Sigma_gas, Median_Sigma_ratio, '-', lw=2, color='white')
                plt.plot(Median_Sigma_gas, Median_Sigma_ratio, '-', lw=1.2, color=color[i],label=name)
                plt.fill_between(Median_Sigma_gas, Sigma_ratio_err_down, Sigma_ratio_err_up, alpha=0.2, color=color[i])
            if method == 'radii':
                data1 = np.loadtxt(f"{siminfo.output_path}/Surface_density_ratio_radii_250pc_%i_"%(index)+name+".txt")
                data2 = np.loadtxt(f"{siminfo.output_path}/Surface_density_ratio_radii_800pc_%i_"%(index)+name+".txt")
                plt.plot(data1[:,0], data1[:,1], 'o', ms=4, color=color[i],alpha=0.2)
                plt.plot(data2[:,0], data2[:,1], 'o', ms=4, color=color[i],label=name)
            i+=1


        plt.xlim(-1.0, 3.0)
        plt.ylim(-8.0, 0.5)
        plt.ylabel(r"log $\Sigma_{\mathrm{H2}} / (\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}})$")
        plt.xlabel(r"log $\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}}$  [M$_{\odot}$ pc$^{-2}$]")
        plt.legend(loc='lower right', labelspacing=0.2, handlelength=2, handletextpad=0.4, frameon=False)
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.savefig(f"{siminfo.output_path}/Surface_density_ratio_"+method+"_%i.png"%(index), dpi=200)
        plt.close()




def make_comparison_plots(siminfo, name_list, num_galaxies):

    if num_galaxies > 10: num_galaxies = 10 #let's plot 10 most massive only

    for index in range(num_galaxies):

        KS_relation_plots(siminfo, index, name_list)
        depletion_time_plots(siminfo, index, name_list)
        surface_ratios_plots(siminfo, index, name_list)