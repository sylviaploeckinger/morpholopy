import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
import h5py as h5
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

def load_GALAH_data():
    file = "./plotter/obs_data/Buder21_data.hdf5"
    GALAH_dataobj = h5.File(file, 'r')
    return GALAH_dataobj

def plot_GALAH_data(element, galah_edges, galah_data):
    obs_plane = np.array(galah_data[f"{element}_enrichment_vs_Fe_abundance"]).T
    obs_plane[obs_plane < 10] = None
    
    contour = plt.contour(np.log10(obs_plane), origin='lower',
                         extent=[galah_edges[0], galah_edges[-1],
                                 galah_edges[0], galah_edges[-1]],
                          zorder=100,
                          cmap='winter',
                          linewidths=1)
    #contour.collections[0].set_label(['GALAH DR3'])


def load_strontium_data_Roeder():
    file = './plotter/obs_data/Roeder_2014.txt'
    # Roeder are calculated wrt to Asplund+ solar metallicity
    data = np.loadtxt(file)
    FeH= data[:, 0]
    SrFe = data[:, 1]
    return FeH, SrFe

def load_strontium_data_Spite():
    file = './plotter/obs_data/Spite_2018.txt'
    # Roeder are calculated wrt to Lodder+ solar metallicity
    Fe_H_Sun = 7.50
    Eu_H_Sun = 0.53
    Ba_H_Sun = 2.18
    Sr_H_Sun = 2.9
    Sr_Fe_Sun = Sr_H_Sun - Fe_H_Sun

    # Asplund et al. (2009)
    Fe_H_Sun_As09 = 7.50
    Eu_H_Sun_As09 = 0.52
    Ba_H_Sun_As09 = 2.18
    Sr_H_Sun_As09 = 2.87

    Eu_Fe_Sun_As09 = Eu_H_Sun_As09 - Fe_H_Sun_As09
    Ba_Fe_Sun_As09 = Ba_H_Sun_As09 - Fe_H_Sun_As09
    Sr_Fe_Sun_As09 = Sr_H_Sun_As09 - Fe_H_Sun_As09

    data = np.loadtxt(file)
    FeH = data[:, 0]
    SrFe = data[:, 1] + Sr_Fe_Sun - Sr_Fe_Sun_As09

    return FeH, SrFe

def load_MW_data_with_Mg_Fe():
    file = './plotter/obs_data/MW.txt'
    data = np.loadtxt(file)
    FeH_mw = data[:, 0]
    MgFe_mw = data[:, 1]
    return FeH_mw, MgFe_mw

def compute_ratios(Hydrogen_fraction, Magnesium_fraction, Oxygen_fraction, Iron_fraction,
                   Carbon_fraction, Silicon_fraction, Europium_fraction, Barium_fraction, Strontium_fraction):
    
    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mFe_in_cgs = 55.845 * mp_in_cgs
    mO_in_cgs = 15.999 * mp_in_cgs
    mMg_in_cgs = 24.305 * mp_in_cgs
    mC_in_cgs =  12.0107 * mp_in_cgs
    mSi_in_cgs = 28.0855 * mp_in_cgs
    mEu_in_cgs = 151.964 * mp_in_cgs
    mBa_in_cgs = 137.327 * mp_in_cgs
    mSr_in_cgs = 87.62 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    O_H_Sun = 8.69
    Mg_H_Sun = 7.6
    C_H_Sun = 8.43
    Si_H_Sun = 7.51
    Eu_H_Sun = 0.52
    Ba_H_Sun = 2.18
    Sr_H_Sun = 2.87

    O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
    Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
    C_Fe_Sun = C_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mC_in_cgs)
    Si_Fe_Sun = Si_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mSi_in_cgs)
    Eu_Fe_Sun = Eu_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mEu_in_cgs)
    Ba_Fe_Sun = Ba_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mBa_in_cgs)
    Sr_Fe_Sun = Sr_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mSr_in_cgs)
    Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)

    Fe_H = np.log10(Iron_fraction / Hydrogen_fraction) - Fe_H_Sun
    O_Fe = np.log10(Oxygen_fraction / Iron_fraction) - O_Fe_Sun
    Mg_Fe = np.log10(Magnesium_fraction / Iron_fraction) - Mg_Fe_Sun
    C_Fe = np.log10(Carbon_fraction / Iron_fraction) - C_Fe_Sun
    Si_Fe = np.log10(Silicon_fraction / Iron_fraction) - Si_Fe_Sun
    Eu_Fe = np.log10(Europium_fraction / Iron_fraction) - Eu_Fe_Sun
    Ba_Fe = np.log10(Barium_fraction / Iron_fraction) - Ba_Fe_Sun
    Sr_Fe = np.log10(Strontium_fraction / Iron_fraction) - Sr_Fe_Sun

    # Let's set lower and upper limits:
    Fe_H[Iron_fraction == 0] = -7  # set lower limit
    Fe_H[Fe_H < -7] = -7  # set lower limit
    Mg_Fe[Iron_fraction == 0] = -2  # set lower limit
    Mg_Fe[Magnesium_fraction == 0] = -2  # set lower limit
    Mg_Fe[Mg_Fe < -2] = -2  # set lower limit
    O_Fe[Iron_fraction == 0] = -2  # set lower limit
    O_Fe[Oxygen_fraction == 0] = -2  # set lower limit
    O_Fe[O_Fe < -2] = -2  # set lower limit
    C_Fe[Iron_fraction == 0] = -2  # set lower limit
    C_Fe[Carbon_fraction == 0] = -2  # set lower limit
    C_Fe[C_Fe < -2] = -2  # set lower limit
    Si_Fe[Iron_fraction == 0] = -2  # set lower limit
    Si_Fe[Silicon_fraction == 0] = -2  # set lower limit
    Si_Fe[Si_Fe < -2] = -2  # set lower limit
    Eu_Fe[Iron_fraction == 0] = -2  # set lower limit
    Eu_Fe[Europium_fraction == 0] = -2  # set lower limit
    Eu_Fe[Eu_Fe < -2] = -2  # set lower limit
    Ba_Fe[Iron_fraction == 0] = -2  # set lower limit
    Ba_Fe[Barium_fraction == 0] = -2  # set lower limit
    Ba_Fe[Ba_Fe < -2] = -2  # set lower limit
    Sr_Fe[Iron_fraction == 0] = -2  # set lower limit
    Sr_Fe[Strontium_fraction == 0] = -2  # set lower limit
    Sr_Fe[Sr_Fe < -2] = -2  # set lower limit

    return {'Fe_H':Fe_H, 'O_Fe':O_Fe, 'Mg_Fe':Mg_Fe,
            'C_Fe':C_Fe, 'Si_Fe':Si_Fe, 'Eu_Fe':Eu_Fe,
            'Ba_Fe':Ba_Fe, 'Sr_Fe':Sr_Fe}

def calculate_abundaces_from_MW_type_galaxies(sim_info):

    select_MW_mass = np.where((sim_info.halo_data.log10_stellar_mass >= 9.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 10.5))[0]
    select_centrals = np.where(sim_info.halo_data.type[select_MW_mass] == 10)[0]

    sample = select_MW_mass[select_centrals]

    Oxygen_fraction = []
    Iron_fraction = []
    Magnesium_fraction = []
    Hydrogen_fraction = []
    Carbon_fraction = []
    Silicon_fraction = []
    Europium_fraction = []
    Strontium_fraction = []
    Barium_fraction = []

    for i in tqdm(range(len(sample))):

        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        O_stars = sim_info.stars.oxygen
        Fe_stars = sim_info.stars.iron
        Mg_stars = sim_info.stars.magnesium
        H_stars = sim_info.stars.hydrogen
        C_stars = sim_info.stars.carbon
        Si_stars = sim_info.stars.silicon
        Eu_stars = sim_info.stars.europium
        Ba_stars = sim_info.stars.barium
        Sr_stars = sim_info.stars.strontium
        halo_stars = np.where(sim_info.stars.in_halo == 1)[0]

        Oxygen_fraction = np.append(Oxygen_fraction, O_stars)
        Iron_fraction = np.append(Iron_fraction, Fe_stars)
        Magnesium_fraction = np.append(Magnesium_fraction, Mg_stars)
        Hydrogen_fraction = np.append(Hydrogen_fraction, H_stars)
        Carbon_fraction = np.append(Carbon_fraction, C_stars)
        Silicon_fraction = np.append(Silicon_fraction, Si_stars)
        Europium_fraction = np.append(Europium_fraction, Eu_stars)
        Barium_fraction = np.append(Barium_fraction, Ba_stars)
        Strontium_fraction = np.append(Strontium_fraction, Sr_stars)

    ratios = compute_ratios(Hydrogen_fraction, Magnesium_fraction, 
                            Oxygen_fraction, Iron_fraction,
                            Carbon_fraction, Silicon_fraction,
                            Europium_fraction, Barium_fraction, Strontium_fraction)
    
    return ratios, halo_stars


def calculate_abundaces_from_satellite_galaxies(sim_info):

    select_mass = np.where((sim_info.halo_data.log10_stellar_mass >= 6.) &
                           (sim_info.halo_data.log10_stellar_mass <= 10.))[0]

    select_satellites = np.where(sim_info.halo_data.type[select_mass] > 10)[0]

    sample = select_mass[select_satellites]

    Oxygen_fraction = []
    Iron_fraction = []
    Magnesium_fraction = []
    Hydrogen_fraction = []
    Carbon_fraction = []
    Silicon_fraction = []
    Europium_fraction = []
    Strontium_fraction = []
    Barium_fraction = []

    for i in tqdm(range(len(sample))):
        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        O_stars = sim_info.stars.oxygen
        Fe_stars = sim_info.stars.iron
        Mg_stars = sim_info.stars.magnesium
        H_stars = sim_info.stars.hydrogen
        C_stars = sim_info.stars.carbon
        Si_stars = sim_info.stars.silicon
        Eu_stars = sim_info.stars.europium
        Ba_stars = sim_info.stars.barium
        Sr_stars = sim_info.stars.strontium

        Oxygen_fraction = np.append(Oxygen_fraction, O_stars)
        Iron_fraction = np.append(Iron_fraction, Fe_stars)
        Magnesium_fraction = np.append(Magnesium_fraction, Mg_stars)
        Hydrogen_fraction = np.append(Hydrogen_fraction, H_stars)
        Carbon_fraction = np.append(Carbon_fraction, C_stars)
        Silicon_fraction = np.append(Silicon_fraction, Si_stars)
        Europium_fraction = np.append(Europium_fraction, Eu_stars)
        Barium_fraction = np.append(Barium_fraction, Ba_stars)
        Strontium_fraction = np.append(Strontium_fraction, Sr_stars)

    ratios = compute_ratios(Hydrogen_fraction, Magnesium_fraction, 
                            Oxygen_fraction, Iron_fraction,
                            Carbon_fraction, Silicon_fraction,
                            Europium_fraction, Barium_fraction, Strontium_fraction)
    return ratios


def plot_stellar_abundances(sim_info, output_path, abundance_data):

    # Look for abundance ratios from COLIBRE snaps:
    ratios_MW, halo_stars = calculate_abundaces_from_MW_type_galaxies(sim_info)
    O_Fe = ratios_MW['O_Fe']
    Mg_Fe = ratios_MW['Mg_Fe']
    Fe_H = ratios_MW['Fe_H']

    # Load MW data:
    FeH_MW, OFe_MW = load_MW_data()

    # Load MW data:
    GALAHdata =  load_GALAH_data()
    galah_edges = np.array(GALAHdata["abundance_bin_edges"])
    
    # Load Satellite data:
    ratios_sat = calculate_abundaces_from_satellite_galaxies(sim_info)
    Fe_H_sat = ratios_sat['Fe_H']
    O_Fe_sat = ratios_sat['O_Fe']
    Mg_Fe_sat = ratios_sat['Mg_Fe']

    # Load Satellite data:
    FeH_fornax, OFe_fornax, FeH_sg, OFe_sg, FeH_ca, OFe_ca, FeH_scu, OFe_scu = load_satellites_data()

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 3),
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

    plt.plot(Fe_H, O_Fe, 'o', ms=0.2, color='grey', alpha=0.2)
    plt.plot(Fe_H[halo_stars], O_Fe[halo_stars], 'o', ms=0.5, color='black')

    plt.plot(FeH_MW, OFe_MW, '+', color='tab:blue', ms=4, label='MW')

    plot_GALAH_data('O', galah_edges, GALAHdata)
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(O_Fe[ind == i]) for i in range(1, len(bins)) if len(O_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)
    Fe_H_median = xm.copy()
    O_Fe_median = ym.copy()
    counter = len(xm)

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

    plot_GALAH_data('O', galah_edges, GALAHdata)
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H_sat, bins)
    xm = [np.median(Fe_H_sat[ind == i]) for i in range(1, len(bins)) if len(Fe_H_sat[ind == i]) > 10]
    ym = [np.median(O_Fe_sat[ind == i]) for i in range(1, len(bins)) if len(O_Fe_sat[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue',label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    plt.text(-3.8, 1.3, "Satellite galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/O_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    # New plot. Now turn Magnesium / Fe --------------------------------

    FeH_fornax, MgFe_fornax, FeH_sculptor, MgFe_sculptor, \
    FeH_sagittarius, MgFe_sagittarius, FeH_carina, MgFe_carina = load_satellites_data_Mg_Fe()

    # Load MW data:
    FeH_MW, MgFe_MW = load_MW_data_with_Mg_Fe()

    # Box stellar abundance --------------------------------
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, Mg_Fe, 'o', ms=0.5, color='grey',alpha=0.5)

    plt.plot(FeH_MW, MgFe_MW, '+', color='orange', ms=4, label='MW')
    plot_GALAH_data('Mg', galah_edges, GALAHdata)
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black',label=sim_info.simulation_name)
    Mg_Fe_median = ym.copy()

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

    plot_GALAH_data('Mg', galah_edges, GALAHdata)
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H_sat, bins)
    xm = [np.median(Fe_H_sat[ind == i]) for i in range(1, len(bins)) if len(Fe_H_sat[ind == i]) > 10]
    ym = [np.median(Mg_Fe_sat[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe_sat[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    plt.text(-3.8,1.2,"Satellite galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/Mg_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    # make remaining plots with just GALAH data (Buder+21)
    for el in ['C','Si','Eu','Ba']:
        fig = plt.figure(figsize=(3.8,3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        plt.plot(Fe_H, ratios_MW[f'{el}_Fe'], 'o', ms=0.5, color='grey',alpha=0.5)
        plot_GALAH_data(el, galah_edges, GALAHdata)

        bins = np.arange(-7.2, 1, 0.2)
        ind = np.digitize(Fe_H, bins)
        xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
        ym = [np.median(ratios_MW[f'{el}_Fe'][ind == i]) for i in range(1, len(bins)) if len(ratios_MW[f'{el}_Fe'][ind == i]) > 10]
        plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
        plt.plot(xm, ym, '-', lw=1.5, color='black',label=sim_info.simulation_name)
        if el == 'C': C_Fe_median = ym.copy()
        if el == 'Si': Si_Fe_median = ym.copy()
        if el == 'Eu': Eu_Fe_median = ym.copy()
        if el == 'Ba': Ba_Fe_median = ym.copy()

        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    ######
    FeH_Ro, SrFe_Ro = load_strontium_data_Roeder()
    FeH_Sp, SrFe_Sp = load_strontium_data_Spite()

    fig = plt.figure(figsize=(3.8, 3))
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(Fe_H, ratios_MW['Sr_Fe'], 'o', ms=0.5, color='grey')
    plt.plot(FeH_Ro, SrFe_Ro, '+', color='crimson', ms=4, label='Roederer et al. (2014)')
    plt.plot(FeH_Sp, SrFe_Sp, 'o', color='blue', ms=4, label='Spite et al. (2018)')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(ratios_MW['Sr_Fe'][ind == i]) for i in range(1, len(bins)) if
          len(ratios_MW['Sr_Fe'][ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)
    Sr_Fe_median = ym.copy()

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Sr/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.tight_layout()
    plt.savefig(f"{output_path}/Sr_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    if abundance_data == None:
        abundance_data = {'Fe_H': Fe_H_median,
                          'O_Fe': O_Fe_median,
                          'Mg_Fe': Mg_Fe_median,
                          'C_Fe': C_Fe_median,
                          'Ba_Fe': Ba_Fe_median,
                          'Sr_Fe': Sr_Fe_median,
                          'Si_Fe': Si_Fe_median,
                          'Eu_Fe': Eu_Fe_median,
                          'counter':counter}
    else:
        Fe_H = abundance_data['Fe_H']
        Fe_H = np.append(Fe_H,Fe_H_median)
        O_Fe = abundance_data['O_Fe']
        O_Fe = np.append(O_Fe,O_Fe_median)
        Mg_Fe = abundance_data['Mg_Fe']
        Mg_Fe = np.append(Mg_Fe,Mg_Fe_median)
        C_Fe = abundance_data['C_Fe']
        C_Fe = np.append(C_Fe, C_Fe_median)
        Ba_Fe = abundance_data['Ba_Fe']
        Ba_Fe = np.append(Ba_Fe, Ba_Fe_median)
        Sr_Fe = abundance_data['Sr_Fe']
        Sr_Fe = np.append(Sr_Fe, Sr_Fe_median)
        Si_Fe = abundance_data['Si_Fe']
        Si_Fe = np.append(Si_Fe, Si_Fe_median)
        Eu_Fe = abundance_data['Eu_Fe']
        Eu_Fe = np.append(Eu_Fe, Eu_Fe_median)
        counter_sim = abundance_data['counter']
        counter_sim = np.append(counter_sim, counter)
        abundance_data = {'Fe_H': Fe_H,
                          'O_Fe': O_Fe,
                          'Mg_Fe': Mg_Fe,
                          'C_Fe': C_Fe,
                          'Ba_Fe': Ba_Fe,
                          'Sr_Fe': Sr_Fe,
                          'Eu_Fe': Eu_Fe,
                          'Si_Fe': Si_Fe,
                          'counter': counter_sim}

    return abundance_data

