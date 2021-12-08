import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm

def load_satellites_data_Mg_Fe():
    # -------------------------------------------------------------------------------
    # alpha-enhancement (Mg/Fe), extracted manually from Tolstoy, Hill & Tosi (2009)
    # -------------------------------------------------------------------------------
    file = './plotter/obs_data/Fornax.txt'
    data = np.loadtxt(file)
    FeH_fornax = data[:, 0]
    MgFe_fornax = data[:, 1]

    file = './plotter/obs_data/Sculptor.txt'
    data = np.loadtxt(file)
    FeH_sculptor = data[:, 0]
    MgFe_sculptor = data[:, 1]

    file = './plotter/obs_data/Sagittarius.txt'
    data = np.loadtxt(file)
    FeH_sagittarius = data[:, 0]
    MgFe_sagittarius = data[:, 1]

    file = './plotter/obs_data/Carina.txt'
    data = np.loadtxt(file)
    FeH_carina = data[:, 0]
    MgFe_carina = data[:, 1]

    return FeH_fornax, MgFe_fornax, FeH_sculptor, MgFe_sculptor, \
           FeH_sagittarius, MgFe_sagittarius, FeH_carina, MgFe_carina


def load_satellites_data():
    # compute COLIBRE standard ratios
    Fe_over_H = 12. - 4.5
    Mg_over_H = 12. - 4.4
    O_over_H = 12. - 3.31
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

    file = './plotter/obs_data/Letarte_2007.txt'
    data = np.loadtxt(file, skiprows=1)
    FeH_fornax = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_fornax = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = './plotter/obs_data/Sbordone_2007.txt'
    data = np.loadtxt(file, skiprows=1)
    FeH_sg = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_sg = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = './plotter/obs_data/Koch_2008.txt'
    data = np.loadtxt(file, skiprows=1)
    FeH_ca = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_ca = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = './plotter/obs_data/Geisler_2005.txt'
    data = np.loadtxt(file, skiprows=3)
    FeH_scu = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_scu = data[:, 4] - data[:, 0] + O_over_Fe_AG89 - O_over_Fe

    return FeH_fornax, OFe_fornax, FeH_sg, OFe_sg, FeH_ca, OFe_ca, FeH_scu, OFe_scu

def load_MW_data():
    #compute COLIBRE standard ratios
    Fe_over_H = 12. - 4.5
    Mg_over_H = 12. - 4.4
    O_over_H = 12. - 3.31
    Mg_over_Fe = Mg_over_H - Fe_over_H
    O_over_Fe = O_over_H - Fe_over_H

    # tabulate/compute the same ratios from Anders & Grevesse (1989)
    Fe_over_H_AG89 = 7.67
    Mg_over_H_AG89 = 7.58
    O_over_H_AG89 = 8.93

    # --
    Mg_over_Fe_AG89 = Mg_over_H_AG89 - Fe_over_H_AG89
    O_over_Fe_AG89 = O_over_H_AG89 - Fe_over_H_AG89

    # MW data
    FeH_MW = []
    OFe_MW = []

    file = './plotter/obs_data/Koch_2008.txt'
    data = np.loadtxt(file, skiprows=3)
    FeH_koch = data[:, 1] + Fe_over_H_AG89 - Fe_over_H
    OH_koch = data[:, 2]
    OFe_koch = OH_koch - FeH_koch + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_koch)
    OFe_MW = np.append(OFe_MW, OFe_koch)

    file = './plotter/obs_data/Bai_2004.txt'
    data = np.loadtxt(file, skiprows=3, usecols=[1, 2])
    FeH_bai = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_bai = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_bai)
    OFe_MW = np.append(OFe_MW, OFe_bai)

    file = './plotter/obs_data/Cayrel_2004.txt'
    data = np.loadtxt(file, skiprows=18, usecols=[2, 6])
    FeH_cayrel = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_cayrel = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_cayrel)
    OFe_MW = np.append(OFe_MW, OFe_cayrel)

    file = './plotter/obs_data/Israelian_1998.txt'
    data = np.loadtxt(file, skiprows=3, usecols=[1, 3])
    FeH_isra = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_isra = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_isra)
    OFe_MW = np.append(OFe_MW, OFe_isra)

    file = './plotter/obs_data/Mishenina_1999.txt'
    data = np.loadtxt(file, skiprows=3, usecols=[1, 3])
    FeH_mish = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_mish = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_mish)
    OFe_MW = np.append(OFe_MW, OFe_mish)

    file = './plotter/obs_data/Zhang_Zhao_2005.txt'
    data = np.loadtxt(file, skiprows=3)
    FeH_zhang = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_zhang = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_zhang)
    OFe_MW = np.append(OFe_MW, OFe_zhang)

    return FeH_MW, OFe_MW

def load_MW_data_with_Mg_Fe():
    file = './plotter/obs_data/MW.txt'
    data = np.loadtxt(file)
    FeH_mw = data[:, 0]
    MgFe_mw = data[:, 1]
    return FeH_mw, MgFe_mw

def compute_ratios(Hydrogen_fraction, Magnesium_fraction, Oxygen_fraction, Iron_fraction):

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

    Fe_H = np.log10(Iron_fraction / Hydrogen_fraction) - Fe_H_Sun
    O_Fe = np.log10(Oxygen_fraction / Iron_fraction) - O_Fe_Sun
    Mg_Fe = np.log10(Magnesium_fraction / Iron_fraction) - Mg_Fe_Sun

    # Let's set lower and upper limits:
    Fe_H[Iron_fraction == 0] = -7  # set lower limit
    Fe_H[Fe_H < -7] = -7  # set lower limit
    Mg_Fe[Iron_fraction == 0] = -2  # set lower limit
    Mg_Fe[Magnesium_fraction == 0] = -2  # set lower limit
    Mg_Fe[Mg_Fe < -2] = -2  # set lower limit
    O_Fe[Iron_fraction == 0] = -2  # set lower limit
    O_Fe[Oxygen_fraction == 0] = -2  # set lower limit
    O_Fe[O_Fe < -2] = -2  # set lower limit
    return Fe_H, O_Fe, Mg_Fe

def calculate_abundaces_from_MW_type_galaxies(sim_info):

    select_MW_mass = np.where((sim_info.halo_data.log10_stellar_mass >= 10.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 12.2))[0]
    select_centrals = np.where(sim_info.halo_data.type[select_MW_mass] == 10)[0]

    sample = select_MW_mass[select_centrals]

    Oxygen_fraction = []
    Iron_fraction = []
    Magnesium_fraction = []
    Hydrogen_fraction = []

    for i in tqdm(range(len(sample))):

        gas_data, stars_data = sim_info.make_particle_data(
            halo_id=sim_info.halo_data.halo_ids[sample[i]]
        )

        O_stars = stars_data[:, 12]
        Fe_stars = stars_data[:, 13]
        Mg_stars = stars_data[:, 14]
        H_stars = stars_data[:, 15]

        Oxygen_fraction = np.append(Oxygen_fraction, O_stars)
        Iron_fraction = np.append(Iron_fraction, Fe_stars)
        Magnesium_fraction = np.append(Magnesium_fraction, Mg_stars)
        Hydrogen_fraction = np.append(Hydrogen_fraction, H_stars)

    Fe_H, O_Fe, Mg_Fe = compute_ratios(Hydrogen_fraction, Magnesium_fraction,
                                       Oxygen_fraction, Iron_fraction)

    return Fe_H, O_Fe, Mg_Fe


def calculate_abundaces_from_satellite_galaxies(sim_info):

    select_mass = np.where((sim_info.halo_data.log10_stellar_mass >= 6.) &
                              (sim_info.halo_data.log10_stellar_mass <= 10.5))[0]

    select_satellites = np.where(sim_info.halo_data.type[select_mass] > 10)[0]

    sample = select_mass[select_satellites]

    Oxygen_fraction = []
    Iron_fraction = []
    Magnesium_fraction = []
    Hydrogen_fraction = []

    for i in tqdm(range(len(sample))):
        gas_data, stars_data = sim_info.make_particle_data(
            halo_id=sim_info.halo_data.halo_ids[sample[i]]
        )

        O_stars = stars_data[:, 12]
        Fe_stars = stars_data[:, 13]
        Mg_stars = stars_data[:, 14]
        H_stars = stars_data[:, 15]

        Oxygen_fraction = np.append(Oxygen_fraction, O_stars)
        Iron_fraction = np.append(Iron_fraction, Fe_stars)
        Magnesium_fraction = np.append(Magnesium_fraction, Mg_stars)
        Hydrogen_fraction = np.append(Hydrogen_fraction, H_stars)

    Fe_H, O_Fe, Mg_Fe = compute_ratios(Hydrogen_fraction, Magnesium_fraction,
                                       Oxygen_fraction, Iron_fraction)

    return Fe_H, O_Fe, Mg_Fe


def plot_stellar_abundances(sim_info, output_path):

    # Look for abundance ratios from COLIBRE snaps:
    Fe_H, O_Fe, Mg_Fe = calculate_abundaces_from_MW_type_galaxies(sim_info)

    # Load MW data:
    FeH_MW, OFe_MW = load_MW_data()

    # Load MW data:
    Fe_H_sat, O_Fe_sat, Mg_Fe_sat = calculate_abundaces_from_satellite_galaxies(sim_info)

    # Load Satellite data:
    FeH_fornax, OFe_fornax, FeH_sg, OFe_sg, FeH_ca, OFe_ca, FeH_scu, OFe_scu = load_satellites_data()

    # Plot the interesting quantities

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (6, 3),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)
    fig = plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, O_Fe, 'o', ms=0.2, color='grey')

    plt.plot(FeH_MW, OFe_MW, '+', color='tab:blue', ms=4, label='MW')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(O_Fe[ind == i]) for i in range(1, len(bins)) if len(O_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    ########################
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    plt.plot(Fe_H_sat, O_Fe_sat, 'o', ms=0.2, color='grey')

    plt.plot(FeH_ca, OFe_ca, '>', color='tab:purple', ms=4, label='Carina')
    plt.plot(FeH_scu, OFe_scu, '*', ms=4, color='tab:green', label='Sculptor')
    plt.plot(FeH_fornax, OFe_fornax, 'o', color='tab:orange', ms=4, label='Fornax')
    plt.plot(FeH_sg, OFe_sg, 'v', ms=4, color='crimson', label='Sagittarius')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H_sat, bins)
    xm = [np.median(Fe_H_sat[ind == i]) for i in range(1, len(bins)) if len(Fe_H_sat[ind == i]) > 10]
    ym = [np.median(O_Fe_sat[ind == i]) for i in range(1, len(bins)) if len(O_Fe_sat[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    plt.text(-3.8, 1.3, "Satellite galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/O_Fe" + sim_info.simulation_name + ".png", dpi=200)

    # New plot. Now turn Magnesium / Fe --------------------------------

    FeH_fornax, MgFe_fornax, FeH_sculptor, MgFe_sculptor, \
    FeH_sagittarius, MgFe_sagittarius, FeH_carina, MgFe_carina = load_satellites_data_Mg_Fe()

    # Load MW data:
    FeH_MW, MgFe_MW = load_MW_data_with_Mg_Fe()

    # Box stellar abundance --------------------------------
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, Mg_Fe, 'o', ms=0.5, color='grey')

    plt.plot(FeH_MW, MgFe_MW, '+', color='orange', ms=4, label='MW')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black',label=sim_info.simulation_name)

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    ###############
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    plt.plot(Fe_H_sat, Mg_Fe_sat, 'o', ms=0.5, color='grey')

    plt.plot(FeH_carina, MgFe_carina, 'o', color='crimson', ms=4, label='Carina')
    plt.plot(FeH_sculptor, MgFe_sculptor, '>', color='khaki', ms=4, label='Sculptor')
    plt.plot(FeH_fornax, MgFe_fornax, '<', color='royalblue', ms=4, label='Fornax')
    plt.plot(FeH_sagittarius, MgFe_sagittarius, '*', ms=4, color='lightblue', label='Sagittarius')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    plt.text(-3.8,1.2,"Satellite galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/Mg_Fe" + sim_info.simulation_name + ".png", dpi=200)