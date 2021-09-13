from pylab import *
from .KS_relation import median_relations
import numpy as np


def plot_combined_surface_densities(
    sigma_SFR, sigma_gas, sigma_H2, stellar_mass, output_path, simulation_name
):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-2, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5.5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.85,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 4,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(
        np.log10(Sigma_g),
        np.log10(Sigma_star),
        "--",
        color="grey",
        label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$",
    )
    plt.scatter(
        sigma_H2,
        sigma_SFR,
        c=stellar_mass,
        alpha=0.9,
        s=25,
        vmin=6,
        vmax=10,
        cmap="CMRmap_r",
        edgecolors="none",
        zorder=2,
    )

    plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel(
        "log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$"
    )
    plt.xlim(-1.0, 4.0)
    plt.ylim(-6.0, 1.0)

    plt.legend()
    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[6, 7, 8, 9, 10, 11, 12], cax=cbar_ax)
    cb.set_label(label="$\log_{10}$ M$_{*}$/M$_{\odot}$", labelpad=0.5)
    ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
    plt.savefig(
        f"{output_path}/surface_density_gas_" + simulation_name + ".png", dpi=200
    )
    plt.close()

    #######
    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(
        np.log10(Sigma_g),
        np.log10(Sigma_star),
        "--",
        color="grey",
        label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$",
    )
    plt.scatter(
        sigma_gas,
        sigma_SFR,
        c=stellar_mass,
        alpha=0.9,
        s=25,
        vmin=6,
        vmax=10,
        cmap="CMRmap_r",
        edgecolors="none",
        zorder=2,
    )

    plt.xlabel("log $\\Sigma_{HI}+ \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel(
        "log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$"
    )
    plt.xlim(-1.0, 4.0)
    plt.ylim(-6.0, 1.0)
    plt.legend()
    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[6, 7, 8, 9, 10, 11, 12], cax=cbar_ax)
    cb.set_label(label="$\log_{10}$ M$_{*}$/M$_{\odot}$", labelpad=0.5)
    ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
    plt.savefig(
        f"{output_path}/surface_density_H2_" + simulation_name + ".png", dpi=200
    )
    plt.close()


def plot_combined_density_ratios(galaxy_data, output_path, simulation_name):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5.5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.85,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 2,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    Sigma_g = np.logspace(-2, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)
    plt.plot(
        np.log10(Sigma_g),
        np.log10(Sigma_star),
        "--",
        color="grey",
        label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$",
    )
    Sigma_g = np.logspace(-1, 4, 1000)
    Sigma_star = KS(Sigma_g, 1.06, 2.511e-4)
    plt.plot(
        np.log10(Sigma_g),
        np.log10(Sigma_star),
        lw=1,
        color="green",
        label=r"2.51e-4 $\times$ $\Sigma_{g}^{1.06}$ (Pessa+ 2021)",
        linestyle="-",
    )

    sigma_gas = galaxy_data.surface_density[1:]
    sigma_SFR = galaxy_data.SFR_density[1:]
    metals = galaxy_data.metallicity[1:]
    metals[metals == 0] = 1e-6
    metals = np.log10(np.array(metals, dtype="float"))
    arg_sort = np.argsort(metals)
    metals = metals[arg_sort[::-1]]
    sigma_gas = sigma_gas[arg_sort[::-1]]
    sigma_SFR = sigma_SFR[arg_sort[::-1]]

    plt.scatter(
        sigma_gas,
        sigma_SFR,
        c=metals,
        alpha=0.9,
        s=10,
        vmin=-3,
        vmax=1,
        cmap="CMRmap_r",
        edgecolors="none",
        zorder=2,
    )

    select = np.where((metals > -1.2) & (metals < -0.8))[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_SFR[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(
        x,
        y,
        "-",
        lw=1.5,
        color="crimson",
        label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$=-1",
    )

    select = np.where((metals > -0.2) & (metals < 0.2))[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_SFR[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(
        x,
        y,
        "-",
        lw=1.5,
        color="mediumpurple",
        label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$=0",
    )

    select = np.where(metals > 0.6)[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_SFR[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(
        x,
        y,
        "-",
        lw=1.5,
        color="lightblue",
        label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$=1",
    )

    x, y, y_down, y_up = median_relations(sigma_gas, sigma_SFR)
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(x, y, "-", lw=1.5, color="grey", label="All")

    select = np.where(sigma_SFR > -5.5)[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_SFR[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(x, y, "-", lw=1.5, color="black")

    plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel(
        "log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$"
    )
    plt.xlim(-1.0, 4.0)
    plt.ylim(-6.0, 1.0)
    plt.legend(
        loc=[0.0, 0.5],
        labelspacing=0.2,
        handlelength=1,
        handletextpad=0.2,
        frameon=False,
    )

    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[-3, -2, -1, 0, 1], cax=cbar_ax)
    cb.set_label(label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$", labelpad=0.5)
    ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
    plt.savefig(
        f"{output_path}/combined_surface_density_gas_" + simulation_name + ".png",
        dpi=200,
    )
    plt.close()

    #######
    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    sigma_gas = galaxy_data.surface_density[1:]
    sigma_ratio = galaxy_data.ratio_densities[1:]
    metals = galaxy_data.metallicity[1:]
    metals[metals == 0] = 1e-6
    metals = np.log10(np.array(metals, dtype="float"))
    arg_sort = np.argsort(metals)
    metals = metals[arg_sort[::-1]]
    sigma_gas = sigma_gas[arg_sort[::-1]]
    sigma_ratio = sigma_ratio[arg_sort[::-1]]

    plt.scatter(
        sigma_gas,
        sigma_ratio,
        c=metals,
        alpha=0.9,
        s=5,
        vmin=-3,
        vmax=1,
        cmap="CMRmap_r",
        edgecolors="none",
        zorder=2,
        label="method:grid",
    )

    select = np.where((metals > -1.2) & (metals < -0.8))[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_ratio[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(
        x,
        y,
        "-",
        lw=1.5,
        color="crimson",
        label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$=-1",
    )

    select = np.where((metals > -0.2) & (metals < 0.2))[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_ratio[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(
        x,
        y,
        "-",
        lw=1.5,
        color="mediumpurple",
        label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$=0",
    )

    select = np.where(metals > 0.6)[0]
    x, y, y_down, y_up = median_relations(sigma_gas[select], sigma_ratio[select])
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(
        x,
        y,
        "-",
        lw=1.5,
        color="lightblue",
        label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$=1",
    )

    x, y, y_down, y_up = median_relations(sigma_gas, sigma_ratio)
    plt.plot(x, y, "-", lw=2, color="white")
    plt.plot(x, y, "-", lw=1.5, color="grey", label="All")

    sigma_gas = galaxy_data.radii_surface_density[1:]
    sigma_ratio = galaxy_data.radii_surface_ratio[1:]
    plt.plot(
        sigma_gas,
        sigma_ratio,
        "o",
        ms=4,
        alpha=0.5,
        color="tab:blue",
        label="method:annuli",
    )

    plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel(
        r"log $\Sigma_{\mathrm{H2}} / (\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}})$"
    )
    plt.xlim(-1.0, 4.0)
    plt.ylim(-8.0, 0.5)
    plt.legend(
        loc=[0.65, 0.0],
        labelspacing=0.2,
        handlelength=1,
        handletextpad=0.2,
        frameon=False,
    )

    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[-3, -2, -1, 0, 1], cax=cbar_ax)
    cb.set_label(label="log Z$_{\mathrm{gas}}$/Z$_{\odot}$", labelpad=0.5)
    ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
    plt.savefig(
        f"{output_path}/combined_surface_density_ratios_" + simulation_name + ".png",
        dpi=200,
    )
    plt.close()


def plot_surface_densities(galaxy_data, output_path, simulation_name):

    # plot surface densities
    plot_combined_surface_densities(
        galaxy_data.sigma_SFR,
        galaxy_data.sigma_gas,
        galaxy_data.sigma_H2,
        galaxy_data.log10_stellar_mass,
        output_path,
        simulation_name,
    )

    # plot surface densities
    plot_combined_density_ratios(galaxy_data, output_path, simulation_name)
