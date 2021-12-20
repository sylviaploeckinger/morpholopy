from pylab import *
import numpy as np
import os

from morpholopy.plotter.html import (
    make_web,
    add_metadata_to_web,
    PlotsInPipeline,
    render_population_web,
    add_web_section,
)


def plot_population_data(output_path, name_list):

    # Color list
    color = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

    # Mass bin range
    bins = np.arange(6, 11, 0.25)
    mass_bins = 0.5 * (bins[1:] + bins[:-1])

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 3.5),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "lines.markersize": 4,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        histogram_all, _ = np.histogram(stellar_mass, bins=bins)
        plt.plot(mass_bins, histogram_all, color=color[i], label=name)
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Number of galaxies per bin")
    plt.xlim(6, 12)
    # plt.ylim(-6.0, 0.0)

    plt.legend()
    plt.savefig(f"{output_path}/number_galaxies.png", dpi=200)
    plt.close()

    #############

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        sfr = data[:, 0]
        ssfr = sfr / 10 ** stellar_mass
        passive_galaxies = np.where(ssfr < 1e-11)[0]  # yr^-1
        histogram_passive, _ = np.histogram(stellar_mass[passive_galaxies], bins=bins)
        plt.plot(mass_bins, histogram_passive, color=color[i], label=name)
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Number of passive galaxies per bin")
    plt.xlim(6, 12)
    # plt.ylim(-6.0, 0.0)

    plt.legend()
    plt.savefig(f"{output_path}/number_passive_galaxies.png", dpi=200)
    plt.close()

    #############

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        sfr = data[:, 0]
        ssfr = sfr / 10 ** stellar_mass
        active_galaxies = np.where(ssfr >= 1e-11)[0]  # yr^-1
        histogram_active, _ = np.histogram(stellar_mass[active_galaxies], bins=bins)
        plt.plot(mass_bins, histogram_active, color=color[i], label=name)
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Number of active galaxies per bin")
    plt.xlim(6, 12)
    # plt.ylim(-6.0, 0.0)

    plt.legend()
    plt.savefig(f"{output_path}/number_active_galaxies.png", dpi=200)
    plt.close()

    #############

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        sfr = data[:, 0]
        ssfr = sfr / 10 ** stellar_mass
        passive_galaxies = np.where(ssfr < 1e-11)[0]  # yr^-1
        histogram_passive, _ = np.histogram(stellar_mass[passive_galaxies], bins=bins)
        histogram_all, _ = np.histogram(stellar_mass, bins=bins)
        zero_num = histogram_all == 0
        histogram_all[zero_num] = 1
        fraction_passive = histogram_passive / histogram_all
        minimum_range = histogram_all > 2
        plt.plot(
            mass_bins[minimum_range],
            fraction_passive[minimum_range],
            color=color[i],
            label=name,
        )
        plt.plot(mass_bins, fraction_passive, "-", lw=0.5, color=color[i])
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Fraction of passive galaxies per bin")
    plt.xlim(6, 12)
    plt.ylim(0, 1)

    plt.legend()
    plt.savefig(f"{output_path}/fraction_passive_galaxies.png", dpi=200)
    plt.close()

    #############

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        sfr = data[:, 0]
        ssfr = sfr / 10 ** stellar_mass
        active_galaxies = np.where(ssfr >= 1e-11)[0]  # yr^-1
        histogram_active, _ = np.histogram(stellar_mass[active_galaxies], bins=bins)
        histogram_all, _ = np.histogram(stellar_mass, bins=bins)
        zero_num = histogram_all == 0
        histogram_all[zero_num] = 1
        fraction_active = histogram_active / histogram_all
        minimum_range = histogram_all > 2
        plt.plot(
            mass_bins[minimum_range],
            fraction_active[minimum_range],
            color=color[i],
            label=name,
        )
        plt.plot(mass_bins, fraction_active, "-", lw=0.5, color=color[i])
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Fraction of active galaxies per bin")
    plt.xlim(6, 12)
    plt.ylim(0, 1)

    plt.legend()
    plt.savefig(f"{output_path}/fraction_active_galaxies.png", dpi=200)
    plt.close()

    #########

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        histogram_all, _ = np.histogram(stellar_mass, bins=bins)
        cumulative_histogram_all = np.cumsum(histogram_all)
        plt.plot(mass_bins, cumulative_histogram_all, color=color[i], label=name)
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Cumulative number of galaxies")
    plt.xlim(6, 12)
    # plt.ylim(-6.0, 0.0)

    plt.legend()
    plt.savefig(f"{output_path}/cumulative_number_galaxies.png", dpi=200)
    plt.close()

    #############

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        sfr = data[:, 0]
        ssfr = sfr / 10 ** stellar_mass
        passive_galaxies = np.where(ssfr < 1e-11)[0]  # yr^-1
        histogram_passive, _ = np.histogram(stellar_mass[passive_galaxies], bins=bins)
        cumulative_histogram_passive = np.cumsum(histogram_passive)
        plt.plot(mass_bins, cumulative_histogram_passive, color=color[i], label=name)
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Cumulative number of passive galaxies")
    plt.xlim(6, 12)
    # plt.ylim(-6.0, 0.0)

    plt.legend()
    plt.savefig(f"{output_path}/cumulative_number_passive_galaxies.png", dpi=200)
    plt.close()

    #############

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    i = 0
    for name in name_list:
        data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")
        stellar_mass = data[:, 1]
        sfr = data[:, 0]
        ssfr = sfr / 10 ** stellar_mass
        active_galaxies = np.where(ssfr >= 1e-11)[0]  # yr^-1
        histogram_active, _ = np.histogram(stellar_mass[active_galaxies], bins=bins)
        cumulative_histogram_active = np.cumsum(histogram_active)
        plt.plot(mass_bins, cumulative_histogram_active, color=color[i], label=name)
        i += 1

    plt.xlabel("log$_{10}$ M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Cumulative number of active galaxies")
    plt.xlim(6, 12)
    # plt.ylim(-6.0, 0.0)

    plt.legend()
    plt.savefig(f"{output_path}/cumulative_number_active_galaxies.png", dpi=200)
    plt.close()

    #############


def load_population_plots(web, name_list, output_path):

    PlotsInWeb = PlotsInPipeline()

    title = "Number of galaxies"
    id = abs(hash("Number of galaxies per stellar mass bin"))
    outfile = "number_galaxies.png"
    caption = "Number of galaxies per stellar mass bin."
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Number of active galaxies"
    id = abs(hash("Number of active galaxies per stellar mass bin"))
    outfile = "number_active_galaxies.png"
    caption = "Number of active galaxies per stellar mass bin."
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Number of passive galaxies"
    id = abs(hash("Number of passive galaxies per stellar mass bin"))
    outfile = "number_passive_galaxies.png"
    caption = "Number of passive galaxies per stellar mass bin."
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Fraction of active galaxies"
    id = abs(hash("Fraction of active galaxies per stellar mass bin"))
    outfile = "fraction_active_galaxies.png"
    caption = (
        "Fraction of active galaxies per stellar mass bin. Solid lines "
        "indicate the fraction of galaxies in those bins where the number of "
        "galaxies is larger than 2. Thin solid lines just indicate the fraction."
    )
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Fraction of passive galaxies"
    id = abs(hash("Fraction of passive galaxies per stellar mass bin"))
    outfile = "fraction_passive_galaxies.png"
    caption = (
        "Fraction of passive galaxies per stellar mass bin. Solid lines "
        "indicate the fraction of galaxies in those bins where the number of "
        "galaxies is larger than 2. Thin solid lines just indicate the fraction."
    )
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Cumulative number of galaxies"
    id = abs(hash("Cumulative number of galaxies per stellar mass bin"))
    outfile = "cumulative_number_galaxies.png"
    caption = "Cumulative number of galaxies per stellar mass bin."
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Cumulative number of active galaxies"
    id = abs(hash("Cumulative number of active galaxies per stellar mass bin"))
    outfile = "cumulative_number_active_galaxies.png"
    caption = "Cumulative number of active galaxies per stellar mass bin."
    PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Cumulative number of passive galaxies"
    id = abs(hash("Cumulative number of passive galaxies per stellar mass bin"))
    outfile = "cumulative_number_passive_galaxies.png"
    caption = "Cumulative number of passive galaxies per stellar mass bin."
    PlotsInWeb.load_plots(title, caption, outfile, id)

    id = abs(hash("Population plots"))
    plots = PlotsInWeb.plots_details
    title = "Galaxy population plots"
    caption = (
        "Stellar mass bins are [6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75 "
        "9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75] [log10 Msun]"
    )
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()


class SimInfo:
    def __init__(self, folder, snap, output_path, name):
        self.name = name
        self.output_path = output_path
        self.snapshot = os.path.join(folder, "colibre_%04i.hdf5" % snap)


if __name__ == "__main__":

    from morpholopy.argumentparser import ArgumentParser

    config_parameters = ArgumentParser()

    # Load MorpholoPy production details
    output_path = config_parameters.output_directory
    number_of_inputs = len(config_parameters.snapshot_list)
    directory_list = config_parameters.directory_list
    snapshot_list = config_parameters.snapshot_list
    name_list = config_parameters.name_list

    # Loop over simulation list
    for sims in range(number_of_inputs):
        directory = directory_list[sims]
        snap_number = int(snapshot_list[sims])
        sim_name = name_list[sims]
        siminfo = SimInfo(directory, snap_number, output_path, sim_name)

        # Make initial website
        if sims == 0:
            web = make_web(siminfo)
        if sims > 0:
            add_metadata_to_web(web, siminfo)

    # Make plots
    plot_population_data(output_path, name_list)

    # Load galaxy plots
    load_population_plots(web, name_list, siminfo.output_path)

    # Finish and output html file
    render_population_web(web, siminfo.output_path)
