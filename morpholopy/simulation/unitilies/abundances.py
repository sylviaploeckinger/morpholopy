from pylab import *
import numpy as np
import h5py
import os
import warnings

warnings.filterwarnings("ignore")
from morpholopy.plotter.html import (
    make_web,
    add_metadata_to_web,
    PlotsInPipeline,
    render_abundance_web,
    add_web_section,
)


def read_data(siminfo):

    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mFe_in_cgs = 55.845 * mp_in_cgs
    mO_in_cgs = 15.999 * mp_in_cgs
    mMg_in_cgs = 24.305 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    O_H_Sun = 8.69
    Mg_H_Sun = 7.6

    O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
    Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
    Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)

    # Read the simulation data
    sim = h5py.File(siminfo.snapshot_name, "r")
    star_abundances = sim["/PartType4/ElementMassFractions"][:][:]
    Fe_H = np.log10(star_abundances[:, 8] / star_abundances[:, 0]) - Fe_H_Sun
    O_Fe = np.log10(star_abundances[:, 4] / star_abundances[:, 8]) - O_Fe_Sun
    Mg_Fe = np.log10(star_abundances[:, 6] / star_abundances[:, 8]) - Mg_Fe_Sun
    redshift = sim["/Header"].attrs["Redshift"][0]
    Fe = star_abundances[:, 8]
    Mg = star_abundances[:, 6]
    O = star_abundances[:, 4]
    Fe_H[Fe == 0] = -7  # set lower limit
    Fe_H[Fe_H < -7] = -7  # set lower limit
    Mg_Fe[Fe == 0] = -2  # set lower limit
    Mg_Fe[Mg == 0] = -2  # set lower limit
    Mg_Fe[Mg_Fe < -2] = -2  # set lower limit
    O_Fe[Fe == 0] = -2  # set lower limit
    O_Fe[O == 0] = -2  # set lower limit
    O_Fe[O_Fe < -2] = -2  # set lower limit
    return Fe_H, O_Fe, Mg_Fe, redshift


def read_obs_data_OFe():

    # compute COLIBRE standard ratios
    Fe_over_H = 12.0 - 4.5
    Mg_over_H = 12.0 - 4.4
    O_over_H = 12.0 - 3.31
    Mg_over_Fe = Mg_over_H - Fe_over_H
    O_over_Fe = O_over_H - Fe_over_H

    # tabulate/compute the same ratios from Anders & Grevesse (1989)
    Fe_over_H_AG89 = 7.67
    Mg_over_H_AG89 = 7.58
    O_over_H_AG89 = 8.93

    # --
    Mg_over_Fe_AG89 = Mg_over_H_AG89 - Fe_over_H_AG89
    O_over_Fe_AG89 = O_over_H_AG89 - Fe_over_H_AG89

    ## I assume these works use Grevesser & Anders solar metallicity

    file = "../../plotter/obs_data/Letarte_2007.txt"
    data = np.loadtxt(file, skiprows=1)
    FeH_fornax = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_fornax = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = "../../plotter/obs_data/Sbordone_2007.txt"
    data = np.loadtxt(file, skiprows=1)
    FeH_sg = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_sg = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = "../../plotter/obs_data/Koch_2008.txt"
    data = np.loadtxt(file, skiprows=1)
    FeH_ca = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_ca = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = "../../plotter/obs_data/Geisler_2005.txt"
    data = np.loadtxt(file, skiprows=3)
    FeH_scu = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_scu = data[:, 4] - data[:, 0] + O_over_Fe_AG89 - O_over_Fe

    # MW data
    FeH_MW = []
    OFe_MW = []

    file = "../../plotter/obs_data/Koch_2008.txt"
    data = np.loadtxt(file, skiprows=3)
    FeH_koch = data[:, 1] + Fe_over_H_AG89 - Fe_over_H
    OH_koch = data[:, 2]
    OFe_koch = OH_koch - FeH_koch + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_koch)
    OFe_MW = np.append(OFe_MW, OFe_koch)

    file = "../../plotter/obs_data/Bai_2004.txt"
    data = np.loadtxt(file, skiprows=3, usecols=[1, 2])
    FeH_bai = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_bai = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_bai)
    OFe_MW = np.append(OFe_MW, OFe_bai)

    file = "../../plotter/obs_data/Cayrel_2004.txt"
    data = np.loadtxt(file, skiprows=18, usecols=[2, 6])
    FeH_cayrel = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_cayrel = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_cayrel)
    OFe_MW = np.append(OFe_MW, OFe_cayrel)

    file = "../../plotter/obs_data/Israelian_1998.txt"
    data = np.loadtxt(file, skiprows=3, usecols=[1, 3])
    FeH_isra = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_isra = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_isra)
    OFe_MW = np.append(OFe_MW, OFe_isra)

    file = "../../plotter/obs_data/Mishenina_1999.txt"
    data = np.loadtxt(file, skiprows=3, usecols=[1, 3])
    FeH_mish = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_mish = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_mish)
    OFe_MW = np.append(OFe_MW, OFe_mish)

    file = "../../plotter/obs_data/Zhang_Zhao_2005.txt"
    data = np.loadtxt(file, skiprows=3)
    FeH_zhang = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_zhang = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_zhang)
    OFe_MW = np.append(OFe_MW, OFe_zhang)
    return (
        FeH_MW,
        OFe_MW,
        FeH_ca,
        OFe_ca,
        FeH_scu,
        OFe_scu,
        FeH_fornax,
        OFe_fornax,
        FeH_sg,
        OFe_sg,
    )


def read_obs_data_MgFe():
    # -------------------------------------------------------------------------------
    # alpha-enhancement (Mg/Fe), extracted manually from Tolstoy, Hill & Tosi (2009)
    # -------------------------------------------------------------------------------
    file = "../../plotter/obs_data/Fornax.txt"
    data = np.loadtxt(file)
    FeH_fornax = data[:, 0]
    MgFe_fornax = data[:, 1]

    file = "../../plotter/obs_data/Sculptor.txt"
    data = np.loadtxt(file)
    FeH_sculptor = data[:, 0]
    MgFe_sculptor = data[:, 1]

    file = "../../plotter/obs_data/Sagittarius.txt"
    data = np.loadtxt(file)
    FeH_sagittarius = data[:, 0]
    MgFe_sagittarius = data[:, 1]

    file = "../../plotter/obs_data/Carina.txt"
    data = np.loadtxt(file)
    FeH_carina = data[:, 0]
    MgFe_carina = data[:, 1]

    file = "../../plotter/obs_data/MW.txt"
    data = np.loadtxt(file)
    FeH_mw = data[:, 0]
    MgFe_mw = data[:, 1]

    # compute COLIBRE standard ratios
    Fe_over_H = 12.0 - 4.5
    Mg_over_H = 12.0 - 4.4
    O_over_H = 12.0 - 3.31
    Mg_over_Fe = Mg_over_H - Fe_over_H
    O_over_Fe = O_over_H - Fe_over_H

    # tabulate/compute the same ratios from Anders & Grevesse (1989)
    Fe_over_H_AG89 = 7.67
    Mg_over_H_AG89 = 7.58
    O_over_H_AG89 = 8.93
    Mg_over_Fe_AG89 = Mg_over_H_AG89 - Fe_over_H_AG89
    O_over_Fe_AG89 = O_over_H_AG89 - Fe_over_H_AG89

    # correct normalisations for COLIBRE standard
    FeH_fornax += Fe_over_H_AG89 - Fe_over_H
    FeH_sculptor += Fe_over_H_AG89 - Fe_over_H
    FeH_sagittarius += Fe_over_H_AG89 - Fe_over_H
    FeH_carina += Fe_over_H_AG89 - Fe_over_H
    FeH_mw += Fe_over_H_AG89 - Fe_over_H

    MgFe_fornax += Mg_over_Fe_AG89 - Mg_over_Fe
    MgFe_sculptor += Mg_over_Fe_AG89 - Mg_over_Fe
    MgFe_sagittarius += Mg_over_Fe_AG89 - Mg_over_Fe
    MgFe_carina += Mg_over_Fe_AG89 - Mg_over_Fe
    MgFe_mw += Mg_over_Fe_AG89 - Mg_over_Fe
    return (
        FeH_mw,
        MgFe_mw,
        FeH_carina,
        MgFe_carina,
        FeH_sculptor,
        MgFe_sculptor,
        FeH_fornax,
        MgFe_fornax,
        FeH_sagittarius,
        MgFe_sagittarius,
    )


def calculate_relative_abundances(siminfo):

    Fe_H, O_Fe, Mg_Fe, redshift = read_data(siminfo)

    (
        FeH_MW,
        OFe_MW,
        FeH_ca,
        OFe_ca,
        FeH_scu,
        OFe_scu,
        FeH_fornax,
        OFe_fornax,
        FeH_sg,
        OFe_sg,
    ) = read_obs_data_OFe()

    # Plot parameters
    params = {
        "font.size": 13,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.12,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "lines.markersize": 4,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)
    # Plot the interesting quantities
    figure()

    # Box stellar abundance --------------------------------
    subplot(111)
    grid(True)

    plt.plot(Fe_H, O_Fe, "o", ms=0.5, color="grey")

    plt.plot(FeH_MW, OFe_MW, "+", color="orange", ms=4, label="MW")
    plt.plot(FeH_ca, OFe_ca, "o", color="crimson", ms=4, label="Carina")
    plt.plot(FeH_scu, OFe_scu, ">", color="khaki", ms=4, label="Sculptor")
    plt.plot(FeH_fornax, OFe_fornax, "<", color="royalblue", ms=4, label="Fornax")
    plt.plot(FeH_sg, OFe_sg, "*", ms=4, color="lightblue", label="Sagittarius")

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [
        np.median(Fe_H[ind == i])
        for i in range(1, len(bins))
        if len(Fe_H[ind == i]) > 10
    ]
    ym = [
        np.median(O_Fe[ind == i])
        for i in range(1, len(bins))
        if len(O_Fe[ind == i]) > 10
    ]
    plt.plot(xm, ym, "-", lw=1.5, color="black")

    plt.text(-1.9, 2.6, siminfo.simulation_name + " $z$=%0.2f" % redshift)

    xlabel("[Fe/H]", labelpad=2)
    ylabel("[O/Fe]", labelpad=2)
    axis([-7.2, 2, -2, 3])
    plt.legend(
        loc=[0.05, 0.02],
        labelspacing=0.1,
        handlelength=1.5,
        handletextpad=0.1,
        frameon=False,
        ncol=3,
        columnspacing=0.02,
    )
    plt.savefig(
        f"{siminfo.output_path}/FeH_OFe_" + siminfo.simulation_name + ".png", dpi=200
    )

    ###########################################################################

    (
        FeH_mw,
        MgFe_mw,
        FeH_carina,
        MgFe_carina,
        FeH_sculptor,
        MgFe_sculptor,
        FeH_fornax,
        MgFe_fornax,
        FeH_sagittarius,
        MgFe_sagittarius,
    ) = read_obs_data_MgFe()

    figure()

    # Box stellar abundance --------------------------------
    subplot(111)
    grid(True)

    plt.plot(Fe_H, Mg_Fe, "o", ms=0.5, color="grey")

    plt.plot(FeH_mw, MgFe_mw, "+", color="orange", ms=4, label="MW")
    plt.plot(FeH_carina, MgFe_carina, "o", color="crimson", ms=4, label="Carina")
    plt.plot(FeH_sculptor, MgFe_sculptor, ">", color="khaki", ms=4, label="Sculptor")
    plt.plot(FeH_fornax, MgFe_fornax, "<", color="royalblue", ms=4, label="Fornax")
    plt.plot(
        FeH_sagittarius,
        MgFe_sagittarius,
        "*",
        ms=4,
        color="lightblue",
        label="Sagittarius",
    )

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [
        np.median(Fe_H[ind == i])
        for i in range(1, len(bins))
        if len(Fe_H[ind == i]) > 10
    ]
    ym = [
        np.median(Mg_Fe[ind == i])
        for i in range(1, len(bins))
        if len(Mg_Fe[ind == i]) > 10
    ]
    plt.plot(xm, ym, "-", lw=1.5, color="black")

    plt.text(-1.9, 2.6, siminfo.simulation_name + " $z$=%0.2f" % redshift)

    xlabel("[Fe/H]", labelpad=2)
    ylabel("[Mg/Fe]", labelpad=2)
    axis([-7.2, 2, -2, 3])
    plt.legend(
        loc=[0.05, 0.02],
        labelspacing=0.1,
        handlelength=1.5,
        handletextpad=0.1,
        frameon=False,
        ncol=3,
        columnspacing=0.02,
    )
    plt.savefig(
        f"{siminfo.output_path}/FeH_MgFe_" + siminfo.simulation_name + ".png", dpi=200
    )


def load_plots(web, name_list, output_path):

    PlotsInWeb = PlotsInPipeline()

    for name in name_list:
        title = "[O/Fe] vs [Fe/H]"
        id = abs(hash("OFe" + name))
        outfile = "FeH_OFe_" + name + ".png"
        caption = (
            "[Fe/H] vs [O/Fe] using Asplund et al. (2009) values for [Fe/H]Sun = 7.5 and [O/H]Sun = 8.69. "
            "The observational data for MW compiles the works of Mishenina+99, Israelian+98, Cayrel+04, Bai+04, Zhang+05, "
            "Koch+08. Most of these works assume Grevesser & Anders (1989) values for solar metallicity, their were corrected "
            "to Asplund+09. Additional data includes Fornax (Letarte+07), Carina (Kock+05), Sculptor (Geisler+05) and "
            "Sagittarious (Sbordone+07)."
        )
        PlotsInWeb.load_plots(title, caption, outfile, id)

    id = abs(hash("O/Fe_Fe/H"))
    plots = PlotsInWeb.plots_details
    title = "[O/Fe] vs [Fe/H]"
    caption = ""
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for name in name_list:
        title = "[Mg/Fe] vs [Fe/H]"
        id = abs(hash("MgFe" + name))
        outfile = "FeH_MgFe_" + name + ".png"
        caption = (
            "[Fe/H] vs [Mg/Fe] using Asplund et al. (2009) values for [Fe/H]Sun = 7.5 and [Mg/H]Sun = 7.6. "
            "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious corresponds to a data compilation "
            "presented by Tolstoy, Hill & Tosi (2009) and extracted by Rob Crain. Note solar metallicity of "
            "Grevesser & Anders (1989) was corrected to Asplund+09."
        )
        PlotsInWeb.load_plots(title, caption, outfile, id)

    id = abs(hash("Mg/Fe_Fe/H"))
    plots = PlotsInWeb.plots_details
    title = "[Mg/Fe] vs [Fe/H]"
    caption = ""
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

        # Run morpholoPy
        calculate_relative_abundances(siminfo)

    # Load galaxy plots
    load_plots(web, name_list, siminfo.output_path)

    # Finish and output html file
    render_abundance_web(web, siminfo.output_path)
