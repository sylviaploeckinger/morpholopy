import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np

color = ["tab:blue", "tab:orange"]


def output_momentum(stellar_mass, momentum, parttype, output_path, simulation_name):
    index = [i for i, v in enumerate(momentum) if v is not None]
    x = [stellar_mass[i] for i in index]
    y = [momentum[i] for i in index]
    np.savetxt(
        f"{output_path}/momentum_parttype_%i_" % (parttype) + simulation_name + ".txt",
        np.transpose([x, y]),
    )


def output_kappa(stellar_mass, kappa, parttype, output_path, simulation_name):
    index = [i for i, v in enumerate(kappa) if v is not None]
    x = [stellar_mass[i] for i in index]
    y = [kappa[i] for i in index]
    np.savetxt(
        f"{output_path}/Kappa_co_parttype_%i_" % (parttype) + simulation_name + ".txt",
        np.transpose([x, y]),
    )


def output_axis_ratios(
    stellar_mass, axis_ratios, parttype, output_path, simulation_name
):
    y = axis_ratios[:, 0]
    y1 = axis_ratios[:, 1]
    y2 = axis_ratios[:, 2]
    index = [i for i, v in enumerate(y) if v is not None]
    x = [stellar_mass[i] for i in index]
    y0 = [y[i] for i in index]
    y1 = [y1[i] for i in index]
    y2 = [y2[i] for i in index]

    np.savetxt(
        f"{output_path}/Axis_ratios_parttype_%i_" % (parttype)
        + simulation_name
        + ".txt",
        np.transpose([x, y0, y1, y2]),
    )


def output_accumulative_densities(combined_data, output_path, simulation_name):

    np.savetxt(
        f"{output_path}/accumulative_surface_density_gas_grid"
        + simulation_name
        + ".txt",
        np.transpose(
            [
                combined_data.neutral_gas_surface_density,
                combined_data.SFR_surface_density,
                combined_data.gas_metallicity,
            ]
        ),
    )

    np.savetxt(
        f"{output_path}/accumulative_surface_density_H2_grid"
        + simulation_name
        + ".txt",
        np.transpose(
            [
                combined_data.molecular_gas_surface_density,
                combined_data.SFR_surface_density,
                combined_data.gas_metallicity,
            ]
        ),
    )

    np.savetxt(
        f"{output_path}/accumulative_surface_density_ratios_grid"
        + simulation_name
        + ".txt",
        np.transpose(
            [
                combined_data.SFR_surface_density,
                combined_data.H2_to_neutral_surface_density_ratio,
                combined_data.gas_metallicity,
            ]
        ),
    )

    np.savetxt(
        f"{output_path}/accumulative_surface_density_ratios_radii"
        + simulation_name
        + ".txt",
        np.transpose(
            [
                combined_data.radii_neutral_gas_surface_density,
                combined_data.radii_H2_to_neutral_surface_density_ratio,
            ]
        ),
    )


def output_HI_size_mass(HI_size, HI_mass, output_path, simulation_name):
    index = [i for i, v in enumerate(HI_size) if v is not None]
    x = [HI_size[i] for i in index]
    y = [HI_mass[i] for i in index]
    np.savetxt(
        f"{output_path}/HI_size_mass_{simulation_name}.txt", np.transpose([x, y])
    )


def write_morphology_data_to_file(
    galaxy_data, combined_data, output_path, simulation_name
):
    """
    Writes morphology data to a file
    """

    output_momentum(
        10 ** galaxy_data.log10_stellar_mass,
        galaxy_data.momentum,
        4,
        output_path,
        simulation_name,
    )
    output_momentum(
        10 ** galaxy_data.log10_stellar_mass,
        galaxy_data.gas_momentum,
        0,
        output_path,
        simulation_name,
    )

    # plot kappa for stars and gas :
    output_kappa(
        10 ** galaxy_data.log10_stellar_mass,
        galaxy_data.kappa_co,
        4,
        output_path,
        simulation_name,
    )
    output_kappa(
        10 ** galaxy_data.log10_stellar_mass,
        galaxy_data.gas_kappa_co,
        0,
        output_path,
        simulation_name,
    )

    # Axis ratios
    axis_ratios = np.zeros((len(galaxy_data.axis_ca), 3))
    axis_ratios[:, 0] = galaxy_data.axis_ca
    axis_ratios[:, 1] = galaxy_data.axis_cb
    axis_ratios[:, 2] = galaxy_data.axis_ba
    output_axis_ratios(
        10 ** galaxy_data.log10_stellar_mass,
        axis_ratios,
        4,
        output_path,
        simulation_name,
    )

    axis_ratios = np.zeros((len(galaxy_data.gas_axis_ca), 3))
    axis_ratios[:, 0] = galaxy_data.gas_axis_ca
    axis_ratios[:, 1] = galaxy_data.gas_axis_cb
    axis_ratios[:, 2] = galaxy_data.gas_axis_ba
    output_axis_ratios(
        10 ** galaxy_data.log10_stellar_mass,
        axis_ratios,
        0,
        output_path,
        simulation_name,
    )

    # plot surface densities
    output_accumulative_densities(combined_data, output_path, simulation_name)

    # HI size and mass
    output_HI_size_mass(
        galaxy_data.HI_size, galaxy_data.HI_mass, output_path, simulation_name
    )

    return


def plot_surface_densities(
    sigma_SFR, sigma_gas, sigma_H2, stellar_mass, MorphologyPlotsInWeb, output_path
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

    fig = plt.figure()
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
    plt.xlim(-1.0, 3.0)
    plt.ylim(-6.0, 0.0)

    plt.legend()
    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[6, 7, 8, 9, 10, 11, 12], cax=cbar_ax)
    cb.set_label(label="$\log_{10}$ M$_{*}$/M$_{\odot}$", labelpad=0.5)
    ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
    plt.savefig(f"{output_path}/surface_density_gas.png", dpi=200)
    plt.close()

    title = "$\\Sigma_{H_2}$ vs $\\Sigma_{\\rm SFR}$"
    caption = "Integrated surface densities of H2 gas and star-forming gas for each individual galaxy. "
    caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
    filename = "surface_density_gas.png"
    id = abs(hash("surface_density_gas"))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)

    #######
    fig = plt.figure()
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
    plt.xlim(-1.0, 3.0)
    plt.ylim(-6.0, 0.0)
    plt.legend()
    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[6, 7, 8, 9, 10, 11, 12], cax=cbar_ax)
    cb.set_label(label="$\log_{10}$ M$_{*}$/M$_{\odot}$", labelpad=0.5)
    ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
    plt.savefig(f"{output_path}/surface_density_H2.png", dpi=200)
    plt.close()

    caption = "Integrated surface densities of H2+HI gas and star-forming gas for each individual galaxy. "
    caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
    filename = "surface_density_H2.png"
    id = abs(hash("surface_density_H2"))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_morphology(output_path, name_list):

    plot_momentum(output_path, name_list)
    plot_kappa(output_path, name_list)
    plot_axis_ratios(output_path, name_list)


def plot_momentum(output_path, name_list):

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
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    parttype_list = [0, 4]

    for parttype in parttype_list:

        plt.figure()
        ax = plt.subplot(1, 1, 1)
        if parttype == 4:
            ax.set_title("Stellar component")
            ylabel = "$j_{\mathrm{stars}}$ [kpc km/s]"

        if parttype == 0:
            ax.set_title("HI+H2 gas")
            ylabel = "$j_{\mathrm{gas}}$ [kpc km/s]"

        plt.grid("True")

        for i, name in enumerate(name_list):
            data = np.loadtxt(
                f"{output_path}/momentum_parttype_%i_" % (parttype) + name + ".txt"
            )

            # In case of no objects len(shape) = 1
            if len(np.shape(data)) > 1:
                plt.plot(data[:, 0], data[:, 1], "o", color=color[i], label=name)

        plt.xlabel("Stellar Mass [M$_{\odot}$]")
        plt.ylabel(ylabel)
        plt.yscale("log")
        plt.xscale("log")
        plt.xlim(1e6, 1e12)
        plt.ylim(1e-1, 1e4)
        ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
        plt.legend(labelspacing=0.2, handlelength=2, handletextpad=0.4, frameon=False)
        plt.savefig(f"{output_path}/momentum_parttype_%i.png" % parttype, dpi=200)
        plt.close()


def plot_kappa(output_path, name_list):

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
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    parttype_list = [0, 4]
    for parttype in parttype_list:

        plt.figure()
        ax = plt.subplot(1, 1, 1)
        if parttype == 4:
            ax.set_title("Stellar component")
        if parttype == 0:
            ax.set_title("HI+H2 gas")

        plt.grid("True")

        for i, name in enumerate(name_list):
            data = np.loadtxt(
                f"{output_path}/Kappa_co_parttype_%i_" % (parttype) + name + ".txt"
            )

            # In case of no objects len(shape) = 1
            if len(np.shape(data)) > 1:
                plt.plot(data[:, 0], data[:, 1], "o", color=color[i], label=name)

        plt.xscale("log")
        plt.xlabel("Stellar Mass [M$_{\odot}$]")
        plt.ylabel(r"$\kappa_{\mathrm{co}}$")
        plt.xlim(1e6, 1e12)
        plt.ylim(0, 1)
        ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
        plt.legend(labelspacing=0.2, handlelength=2, handletextpad=0.4, frameon=False)
        plt.savefig(f"{output_path}/Kappa_co_parttype_%i.png" % parttype, dpi=200)
        plt.close()


def plot_axis_ratios(output_path, name_list):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (11, 3.5),
        "figure.subplot.left": 0.05,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 2,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    for title, parttype in zip(["HI+H2 gas", "Stellar component"], [0, 4]):

        ########
        plt.figure()
        ax = plt.subplot(1, 3, 1)
        ax.set_title(title)
        plt.grid("True")

        for i, name in enumerate(name_list):
            data = np.loadtxt(
                f"{output_path}/Axis_ratios_parttype_%i_" % (parttype) + name + ".txt"
            )

            # In case of no objects len(shape) = 1
            if len(np.shape(data)) > 1:
                plt.plot(data[:, 0], data[:, 1], "o", color=color[i], label=name)

        plt.xscale("log")
        plt.xlabel("Stellar Mass [M$_{\odot}$]")
        plt.ylabel("c/a")
        plt.xlim(1e6, 1e12)
        plt.ylim(0.0, 1.0)
        ax.tick_params(direction="in", axis="both", which="both", pad=4.5)

        ########
        ax = plt.subplot(1, 3, 2)
        plt.grid("True")

        for i, name in enumerate(name_list):
            data = np.loadtxt(
                f"{output_path}/Axis_ratios_parttype_%i_" % (parttype) + name + ".txt"
            )

            # In case of no objects len(shape) = 1
            if len(np.shape(data)) > 1:
                plt.plot(data[:, 0], data[:, 2], "o", color=color[i], label=name)

        plt.xscale("log")
        plt.xlabel("Stellar Mass [M$_{\odot}$]")
        plt.ylabel("c/b")
        plt.xlim(1e6, 1e12)
        plt.ylim(0.0, 1.0)
        ax.tick_params(direction="in", axis="both", which="both", pad=4.5)

        ########
        ax = plt.subplot(1, 3, 3)
        plt.grid("True")

        for i, name in enumerate(name_list):
            data = np.loadtxt(
                f"{output_path}/Axis_ratios_parttype_%i_" % (parttype) + name + ".txt"
            )

            # In case of no objects len(shape) = 1
            if len(np.shape(data)) > 1:
                plt.plot(data[:, 0], data[:, 3], "o", color=color[i], label=name)

        plt.xscale("log")
        plt.xlabel("Stellar Mass [M$_{\odot}$]")
        plt.ylabel("b/a")
        plt.xlim(1e6, 1e12)
        plt.ylim(0.0, 1.0)
        ax.tick_params(direction="in", axis="both", which="both", pad=4.5)
        plt.legend(labelspacing=0.2, handlelength=2, handletextpad=0.4, frameon=False)
        plt.savefig(f"{output_path}/Axis_ratios_parttype_%i.png" % parttype, dpi=200)
        plt.close()
