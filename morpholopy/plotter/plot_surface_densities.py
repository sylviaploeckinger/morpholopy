import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
from .KS_relation import median_relations
import numpy as np


def plot_integrated_surface_densities(
    sigma_SFR, sigma_gas, sigma_H2, stellar_mass, output_path, simulation_name
):
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

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-2, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    rcParams.update(params)

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(
        np.log10(Sigma_g),
        np.log10(Sigma_star),
        "--",
        color="grey",
        label=r"1.51e-4 $\times$ $\Sigma_{\rm g}^{1.4}$",
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

    plt.xlabel("log $\\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
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
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(
        np.log10(Sigma_g),
        np.log10(Sigma_star),
        "--",
        color="grey",
        label=r"1.51e-4 $\times$ $\Sigma_{\rm g}^{1.4}$",
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

    plt.xlabel(
        "log $\\Sigma_{\\rm HI}+ \\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$"
    )
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


def plot_combined_surface_densities(combined_data, output_path, simulation_name):

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

    sigma_SFR = combined_data.SFR_surface_density
    metals = combined_data.gas_metallicity
    metals[metals == 0] = 1e-6
    metals = np.log10(np.array(metals, dtype="float"))
    arg_sort = np.argsort(metals)
    metals = metals[arg_sort[::-1]]
    sigma_SFR = sigma_SFR[arg_sort[::-1]]

    for plot in ["gas", "H2"]:

        # star formation rate surgface density vs surface density
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        Sigma_g = np.logspace(-2, 3, 1000)
        Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)
        plt.plot(
            np.log10(Sigma_g),
            np.log10(Sigma_star),
            "--",
            color="grey",
            label=r"1.51e-4 $\times$ $\Sigma_{\rm g}^{1.4}$",
        )
        Sigma_g = np.logspace(-1, 4, 1000)
        Sigma_star = KS(Sigma_g, 1.06, 2.511e-4)
        plt.plot(
            np.log10(Sigma_g),
            np.log10(Sigma_star),
            lw=1,
            color="green",
            label=r"2.51e-4 $\times$ $\Sigma_{\rm g}^{1.06}$ (Pessa+ 2021)",
            linestyle="-",
        )

        if plot == "gas":
            sigma_gas = combined_data.neutral_gas_surface_density
            t_gas = combined_data.depletion_time_neutral_gas
        else:
            sigma_gas = combined_data.molecular_gas_surface_density
            t_gas = combined_data.depletion_time_molecular_gas

        sigma_gas = sigma_gas[arg_sort[::-1]]
        t_gas = t_gas[arg_sort[::-1]]

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
        plt.plot(x, y, "-", lw=1.5, color="black", label="star-forming")

        if plot == "gas":
            plt.xlabel(
                "log $\\Sigma_{\\rm HI} + \\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$"
            )
        else:
            plt.xlabel("log $\\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")

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
            f"{output_path}/combined_surface_density_{plot}_"
            + simulation_name
            + ".png",
            dpi=200,
        )
        plt.close()

        # Depletion time vs surface density
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        Sigma_g = np.logspace(-1, 4, 1000)
        Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)
        plt.plot(
            np.log10(Sigma_g),
            np.log10(Sigma_g) - np.log10(Sigma_star) + 6.0,
            color="red",
            label="KS law (Kennicutt 98)",
            linestyle="--",
        )

        plt.scatter(
            sigma_gas,
            t_gas,
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
        x, y, y_down, y_up = median_relations(sigma_gas[select], t_gas[select])
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
        x, y, y_down, y_up = median_relations(sigma_gas[select], t_gas[select])
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

        x, y, y_down, y_up = median_relations(sigma_gas, t_gas)
        plt.plot(x, y, "-", lw=2, color="white")
        plt.plot(x, y, "-", lw=1.5, color="grey", label="All")

        select = np.where(sigma_SFR > -5.5)[0]
        x, y, y_down, y_up = median_relations(sigma_gas[select], t_gas[select])
        plt.plot(x, y, "-", lw=2, color="white")
        plt.plot(x, y, "-", lw=1.5, color="black", label="star-forming")

        if plot == "gas":
            plt.xlabel(
                "log $\\Sigma_{\\rm HI} + \\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$"
            )
            plt.ylabel(
                "log $\\rm t_{gas} = (\\Sigma_{\\rm HI} + \\Sigma_{\\rm H_2})/ \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$"
            )
        else:
            plt.xlabel("log $\\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
            plt.ylabel(
                "log $\\rm t_{H_2} = \\Sigma_{\\rm H_2} / \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$"
            )

        plt.xlim(-1, 4.0)
        plt.ylim(7, 12)
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
            f"{output_path}/depletion_time_combined_surface_density_{plot}_"
            + simulation_name
            + ".png",
            dpi=200,
        )
        plt.close()

    #######
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    sigma_gas = combined_data.neutral_gas_surface_density
    sigma_ratio = combined_data.H2_to_neutral_surface_density_ratio
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

    sigma_gas = combined_data.radii_neutral_gas_surface_density
    sigma_ratio = combined_data.radii_H2_to_neutral_surface_density_ratio
    plt.plot(
        sigma_gas,
        sigma_ratio,
        "o",
        ms=4,
        alpha=0.5,
        color="tab:blue",
        label="method:annuli",
    )

    plt.xlabel(
        "log $\\Sigma_{\\rm HI} + \\Sigma_{\\rm H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$"
    )
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


def plot_surface_densities(
    halo_catalogue_data, combined_data, output_path, simulation_name
):

    plot_integrated_surface_densities(
        halo_catalogue_data.sigma_SFR,
        halo_catalogue_data.sigma_gas,
        halo_catalogue_data.sigma_H2,
        halo_catalogue_data.log10_stellar_mass,
        output_path,
        simulation_name,
    )

    plot_combined_surface_densities(combined_data, output_path, simulation_name)
