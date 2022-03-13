import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm
import scipy.stats as stat


def compute_median_ratios(Oxygen_fraction, Iron_fraction,
                   Magnesium_fraction,Hydrogen_fraction):

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


    select = np.where(Fe_H >= -3)[0]
    Fe_H = Fe_H[select]
    O_Fe = O_Fe[select]
    Mg_Fe = Mg_Fe[select]

    Fe_H = np.median(Fe_H)
    O_Fe = np.median(O_Fe)
    Mg_Fe = np.median(Mg_Fe)

    return {'Fe_H': Fe_H, 'O_Fe': O_Fe, 'Mg_Fe': Mg_Fe}

def compute_mass_weighted_ratios(Oxygen_fraction, Iron_fraction, Magnesium_fraction,
                                 Hydrogen_fraction, stellar_mass):

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

    Fe_H = Iron_fraction / Hydrogen_fraction
    O_Fe = Oxygen_fraction / Iron_fraction
    Mg_Fe = Magnesium_fraction / Iron_fraction

    # Let's set lower and upper limits:
    Mg_Fe[Iron_fraction == 0] = 0.  # set lower limit
    O_Fe[Iron_fraction == 0] = 0.   # set lower limit

    Fe_H_mass_weighted = np.sum( stellar_mass * Iron_fraction) / np.sum( Hydrogen_fraction * stellar_mass)
    O_Fe_mass_weighted = np.sum( stellar_mass * Oxygen_fraction) / np.sum( Iron_fraction * stellar_mass)
    Mg_Fe_mass_weighted = np.sum( stellar_mass * Magnesium_fraction) / np.sum( Iron_fraction * stellar_mass)
    Fe_H_mass_weighted = np.log10(Fe_H_mass_weighted) - Fe_H_Sun
    O_Fe_mass_weighted = np.log10(O_Fe_mass_weighted) - O_Fe_Sun
    Mg_Fe_mass_weighted = np.log10(Mg_Fe_mass_weighted) - Mg_Fe_Sun

    Fe_H_mass_weighted_ratio = np.sum( stellar_mass * Fe_H) / np.sum( stellar_mass)
    O_Fe_mass_weighted_ratio = np.sum( stellar_mass * O_Fe) / np.sum( stellar_mass)
    Mg_Fe_mass_weighted_ratio = np.sum( stellar_mass * Mg_Fe) / np.sum( stellar_mass)
    Fe_H_mass_weighted_ratio = np.log10(Fe_H_mass_weighted_ratio) - Fe_H_Sun
    O_Fe_mass_weighted_ratio = np.log10(O_Fe_mass_weighted_ratio) - O_Fe_Sun
    Mg_Fe_mass_weighted_ratio = np.log10(Mg_Fe_mass_weighted_ratio) - Mg_Fe_Sun

    return {'Fe_H': Fe_H_mass_weighted,
            'O_Fe': O_Fe_mass_weighted,
            'Mg_Fe': Mg_Fe_mass_weighted,
            'Fe_H_weighted_ratio': Fe_H_mass_weighted_ratio,
            'O_Fe_weighted_ratio': O_Fe_mass_weighted_ratio,
            'Mg_Fe_weighted_ratio': Mg_Fe_mass_weighted_ratio}

def compute_metallicity_light_weighted_ratios(Z, Fe, H, mass, L_r_band, L_i_band, L_Z_band):

    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mFe_in_cgs = 55.845 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    Z_Sun = 0.0134
    Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)

    # Ratio
    Fe_H = Fe / H

    Z_weighted_r_band = np.sum(Z * L_r_band) / np.sum(L_r_band)
    Z_weighted_r_band /= Z_Sun

    Z_weighted_i_band = np.sum(Z * L_i_band) / np.sum(L_i_band)
    Z_weighted_i_band /= Z_Sun

    Z_weighted_Z_band = np.sum(Z * L_Z_band) / np.sum(L_Z_band)
    Z_weighted_Z_band /= Z_Sun

    Fe_H_weighted_r_band = np.sum(Fe_H * L_r_band) / np.sum(L_r_band)
    Fe_H_weighted_i_band = np.sum(Fe_H * L_i_band) / np.sum(L_i_band)
    Fe_H_weighted_Z_band = np.sum(Fe_H * L_Z_band) / np.sum(L_Z_band)
    Fe_H_weighted_r_band = np.log10(Fe_H_weighted_r_band) - Fe_H_Sun
    Fe_H_weighted_i_band = np.log10(Fe_H_weighted_i_band) - Fe_H_Sun
    Fe_H_weighted_Z_band = np.log10(Fe_H_weighted_Z_band) - Fe_H_Sun


    Fe_weighted = np.sum(Fe * mass * L_r_band) / np.sum(L_r_band)
    H_weighted = np.sum(H * mass * L_r_band) / np.sum(L_r_band)
    Fe_H_weighted_r_band_with_mass = np.log10(Fe_weighted / H_weighted) - Fe_H_Sun
    Fe_weighted = np.sum(Fe * mass * L_i_band) / np.sum(L_i_band)
    H_weighted = np.sum(H * mass * L_i_band) / np.sum(L_i_band)
    Fe_H_weighted_i_band_with_mass = np.log10(Fe_weighted / H_weighted) - Fe_H_Sun
    Fe_weighted = np.sum(Fe * mass * L_Z_band) / np.sum(L_Z_band)
    H_weighted = np.sum(H * mass * L_Z_band) / np.sum(L_Z_band)
    Fe_H_weighted_Z_band_with_mass = np.log10(Fe_weighted / H_weighted) - Fe_H_Sun

    return {'Fe_H_weighted_r_band': Fe_H_weighted_r_band,
            'Fe_H_weighted_i_band': Fe_H_weighted_i_band,
            'Fe_H_weighted_Z_band': Fe_H_weighted_Z_band,
            'Z_weighted_r_band': Z_weighted_r_band,
            'Z_weighted_i_band': Z_weighted_i_band,
            'Z_weighted_Z_band': Z_weighted_Z_band,
            'Fe_H_weighted_r_band_with_mass': Fe_H_weighted_r_band_with_mass,
            'Fe_H_weighted_i_band_with_mass': Fe_H_weighted_i_band_with_mass,
            'Fe_H_weighted_Z_band_with_mass': Fe_H_weighted_Z_band_with_mass}

def compute_alpha_light_weighted_ratios(Fe, O, Mg, mass, L_r_band):

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

    # Ratios
    O_Fe = O / Fe
    Mg_Fe = Mg / Fe
    Mg_Fe[Fe == 0] = 0.  # set lower limit
    O_Fe[Fe == 0] = 0.  # set lower limit

    O_Fe_weighted_r_band = np.sum(O_Fe * L_r_band) / np.sum(L_r_band)
    Mg_Fe_weighted_r_band = np.sum(Mg_Fe * L_r_band) / np.sum(L_r_band)
    O_Fe_weighted_r_band = np.log10(O_Fe_weighted_r_band) - O_Fe_Sun
    Mg_Fe_weighted_r_band = np.log10(Mg_Fe_weighted_r_band) - Mg_Fe_Sun

    Fe_weighted = np.sum(Fe * mass * L_r_band) / np.sum(L_r_band)
    O_weighted = np.sum(O * mass * L_r_band) / np.sum(L_r_band)
    O_Fe_weighted_r_band_with_mass = np.log10(O_weighted / Fe_weighted) - O_Fe_Sun
    Mg_weighted = np.sum(Mg * mass * L_r_band) / np.sum(L_r_band)
    Mg_Fe_weighted_r_band_with_mass = np.log10(Mg_weighted / Fe_weighted) - Mg_Fe_Sun

    return {'O_Fe_weighted_r_band': O_Fe_weighted_r_band,
            'Mg_Fe_weighted_r_band': Mg_Fe_weighted_r_band,
            'O_Fe_weighted_r_band_with_mass': O_Fe_weighted_r_band_with_mass,
            'Mg_Fe_weighted_r_band_with_mass': Mg_Fe_weighted_r_band_with_mass}

def calculate_galaxies_metallicity(sim_info):

    select_sample = np.where((sim_info.halo_data.log10_stellar_mass >= 7.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 12))[0]

    select_centrals = np.where(sim_info.halo_data.type[select_sample] == 10)[0]

    sample = select_sample[select_centrals]
    num_sample = len(sample)

    FeH_median = np.zeros(num_sample)
    MgFe_median = np.zeros(num_sample)
    OFe_median = np.zeros(num_sample)
    FeH_total_mass_weighted = np.zeros(num_sample)
    MgFe_total_mass_weighted = np.zeros(num_sample)
    OFe_total_mass_weighted = np.zeros(num_sample)
    FeH_ratio_weighted = np.zeros(num_sample)
    MgFe_ratio_weighted = np.zeros(num_sample)
    OFe_ratio_weighted = np.zeros(num_sample)

    Z_mass_weighted = np.zeros(num_sample)
    Z_light_weighted_r_band = np.zeros(num_sample)
    Z_light_weighted_i_band = np.zeros(num_sample)
    Z_light_weighted_Z_band = np.zeros(num_sample)

    FeH_light_weighted_r_band = np.zeros(num_sample)
    FeH_light_weighted_i_band = np.zeros(num_sample)
    FeH_light_weighted_Z_band = np.zeros(num_sample)

    FeH_light_weighted_r_band_with_mass = np.zeros(num_sample)
    FeH_light_weighted_i_band_with_mass = np.zeros(num_sample)
    FeH_light_weighted_Z_band_with_mass = np.zeros(num_sample)

    OFe_light_weighted_r_band = np.zeros(num_sample)
    MgFe_light_weighted_r_band = np.zeros(num_sample)
    OFe_light_weighted_r_band_with_mass = np.zeros(num_sample)
    MgFe_light_weighted_r_band_with_mass = np.zeros(num_sample)

    stellar_mass = sim_info.halo_data.log10_stellar_mass[sample]
    metallicity = sim_info.halo_data.metallicity[sample]

    for i in tqdm(range(len(sample))):
        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        O_stars = sim_info.stars.oxygen.copy()
        Fe_stars = sim_info.stars.iron.copy()
        Mg_stars = sim_info.stars.magnesium.copy()
        H_stars = sim_info.stars.hydrogen.copy()
        M_stars = sim_info.stars.mass.copy()
        Z_stars = sim_info.stars.metal_mass_fractions.copy()
        stars_mag = sim_info.calculate_luminosities(sim_info.stars.age, Z_stars, sim_info.stars.initmass)
        L_stars_r_band = pow(10.0, -0.4 * stars_mag["r"]) * 3631
        L_stars_i_band = pow(10.0, -0.4 * stars_mag["i"]) * 3631
        L_stars_z_band = pow(10.0, -0.4 * stars_mag["z"]) * 3631

        ratios = compute_alpha_light_weighted_ratios(Fe_stars, O_stars, Mg_stars, M_stars, L_stars_r_band)

        OFe_light_weighted_r_band[i] = ratios['O_Fe_weighted_r_band']
        MgFe_light_weighted_r_band[i] = ratios['Mg_Fe_weighted_r_band']

        OFe_light_weighted_r_band_with_mass[i] = ratios['O_Fe_weighted_r_band_with_mass']
        MgFe_light_weighted_r_band_with_mass[i] = ratios['Mg_Fe_weighted_r_band_with_mass']

        ratios = compute_metallicity_light_weighted_ratios(Z_stars, Fe_stars, H_stars, M_stars,
                                                           L_stars_r_band, L_stars_i_band, L_stars_z_band)
        Z_light_weighted_r_band[i] = ratios['Z_weighted_r_band']
        Z_light_weighted_i_band[i] = ratios['Z_weighted_i_band']
        Z_light_weighted_Z_band[i] = ratios['Z_weighted_Z_band']

        FeH_light_weighted_r_band[i] = ratios['Fe_H_weighted_r_band']
        FeH_light_weighted_i_band[i] = ratios['Fe_H_weighted_i_band']
        FeH_light_weighted_Z_band[i] = ratios['Fe_H_weighted_Z_band']

        FeH_light_weighted_r_band_with_mass[i] = ratios['Fe_H_weighted_r_band_with_mass']
        FeH_light_weighted_i_band_with_mass[i] = ratios['Fe_H_weighted_i_band_with_mass']
        FeH_light_weighted_Z_band_with_mass[i] = ratios['Fe_H_weighted_Z_band_with_mass']

        Zsun = 0.0134
        Z_mass_weighted[i] = ( np.sum( Z_stars * M_stars ) / np.sum(M_stars) ) / Zsun

        ratios = compute_median_ratios(O_stars, Fe_stars, Mg_stars, H_stars)
        FeH_median[i] = ratios['Fe_H']
        MgFe_median[i] = ratios['Mg_Fe']
        OFe_median[i] = ratios['O_Fe']


        ratios = compute_mass_weighted_ratios(O_stars, Fe_stars, Mg_stars, H_stars, M_stars)
        FeH_total_mass_weighted[i] = ratios['Fe_H']
        MgFe_total_mass_weighted[i] = ratios['Mg_Fe']
        OFe_total_mass_weighted[i] = ratios['O_Fe']
        FeH_ratio_weighted[i] = ratios['Fe_H_weighted_ratio']
        MgFe_ratio_weighted[i] = ratios['Mg_Fe_weighted_ratio']
        OFe_ratio_weighted[i] = ratios['O_Fe_weighted_ratio']


    return {'Fe_H_median': FeH_median,
            'O_Fe_median': OFe_median,
            'Mg_Fe_median': MgFe_median,

            'Fe_H_mass_weighted': FeH_total_mass_weighted,
            'O_Fe_mass_weighted': OFe_total_mass_weighted,
            'Mg_Fe_mass_weighted': MgFe_total_mass_weighted,

            'Fe_H_ratio_weighted': FeH_ratio_weighted,
            'O_Fe_ratio_weighted': OFe_ratio_weighted,
            'Mg_Fe_ratio_weighted': MgFe_ratio_weighted,

            'Z_mass_weighted': Z_mass_weighted,
            'Z_light_weighted_r_band': Z_light_weighted_r_band,
            'Z_light_weighted_i_band': Z_light_weighted_i_band,
            'Z_light_weighted_Z_band': Z_light_weighted_Z_band,
            'Fe_H_light_weighted_r_band': FeH_light_weighted_r_band,
            'Fe_H_light_weighted_i_band': FeH_light_weighted_i_band,
            'Fe_H_light_weighted_Z_band': FeH_light_weighted_Z_band,

            'Fe_H_light_weighted_r_band_with_mass': FeH_light_weighted_r_band_with_mass,
            'Fe_H_light_weighted_i_band_with_mass': FeH_light_weighted_i_band_with_mass,
            'Fe_H_light_weighted_Z_band_with_mass': FeH_light_weighted_Z_band_with_mass,

            'O_Fe_light_weighted_r_band': OFe_light_weighted_r_band,
            'Mg_Fe_light_weighted_r_band': MgFe_light_weighted_r_band,

            'O_Fe_light_weighted_r_band_with_mass': OFe_light_weighted_r_band_with_mass,
            'Mg_Fe_light_weighted_r_band_with_mass': MgFe_light_weighted_r_band_with_mass,

            'Mstellar':stellar_mass}

def plot_Kirby_distributions(output_path):

    kroupa_to_chabrier_mass = 0.912
    delimiter = "\t"

    # Read the data
    input_filename = "./plotter/obs_data/Kirby_2013_data.txt"
    galaxy_name = np.loadtxt(input_filename, usecols=[0], dtype=str)
    galaxy_mass = np.loadtxt(input_filename, usecols=[1])

    input_filename = "./plotter/obs_data/Kirby_2013_individual_stars.txt"
    raw_galaxy_list = np.loadtxt(input_filename, usecols=[0], dtype=str)
    raw = np.loadtxt(input_filename, usecols=[2, 3])

    num_galaxies = len(galaxy_mass)

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
    plt.grid(True)

    for i in range(num_galaxies):

        pick = np.where(raw_galaxy_list == galaxy_name[i])[0]
        if len(pick) == 0: continue

        Fe_H_data = 10 ** raw[pick, 0]

        bins = np.arange(-7.2,3,0.2)
        bins = 10 ** bins
        hist, _, _ = stat.binned_statistic(Fe_H_data , values=np.ones(len(Fe_H_data)), statistic="sum", bins=bins)
        bin_centers = 0.5*(bins[1:] + bins[:-1])
        hist = hist/np.sum(hist)
        plt.plot(bin_centers, hist, '-', lw=1, label=galaxy_name[i])

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("PDF", labelpad=2)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([2e-4, 2e0, 1e-3, 1e0])
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5,
               handletextpad=0.1, frameon=False, ncol=2,
               fontsize=9, columnspacing=0.02)

    plt.savefig(f"{output_path}/Kirby_distribution.png", dpi=200)


def plot_Kirby_analysed():

    kroupa_to_chabrier_mass = 0.912
    delimiter = "\t"

    # Read the data
    input_filename = "./plotter/obs_data/Kirby_2013_data.txt"
    galaxy_name = np.loadtxt(input_filename, usecols=[0], dtype=str)
    galaxy_mass = np.loadtxt(input_filename, usecols=[1])

    input_filename = "./plotter/obs_data/Kirby_2013_individual_stars.txt"
    raw_galaxy_list = np.loadtxt(input_filename, usecols=[0], dtype=str)
    raw = np.loadtxt(input_filename, usecols=[2, 3])

    num_galaxies = len(galaxy_mass)
    M_star = []
    Z_star = []
    Z_star_lo = []
    Z_star_hi = []

    for i in range(num_galaxies):

        pick = np.where(raw_galaxy_list == galaxy_name[i])[0]
        if len(pick) == 0: continue

        M_star = np.append(M_star, 10 ** galaxy_mass[i] * kroupa_to_chabrier_mass)

        Fe_H_data = 10 ** raw[pick, 0]

        Fe_H = np.median(Fe_H_data)
        Fe_H_lo = np.percentile(Fe_H_data, 16)
        Fe_H_hi = np.percentile(Fe_H_data, 84)

        Z_star = np.append(Z_star, Fe_H)
        Z_star_lo = np.append(Z_star_lo, Fe_H_lo)
        Z_star_hi = np.append(Z_star_hi, Fe_H_hi)

    # Define the scatter as offset from the mean value
    y_scatter = np.array((Z_star - Z_star_lo, Z_star_hi - Z_star))
    plt.errorbar(M_star, Z_star, yerr= y_scatter, fmt='*',
                 ls='none',color='black',lw=1,ms=2.5,label='Kirby et al. (2013)/Analysed')


def plot_Kirby_data():

    kroupa_to_chabrier_mass = 0.912

    input_filename = "./plotter/obs_data/Kirby_2013_ascii.dat"
    delimiter = "\t"

    # Read the data
    raw = np.loadtxt(input_filename, delimiter=delimiter)
    M_star = 10 ** raw[:, 0] * kroupa_to_chabrier_mass
    M_star_lo = 10 ** (raw[:, 0] - raw[:, 1]) * kroupa_to_chabrier_mass
    M_star_hi = 10 ** (raw[:, 0] + raw[:, 2]) * kroupa_to_chabrier_mass

    Z_star = 10 ** raw[:, 3]  # Z/Z_sun
    Z_star_lo = 10 ** (raw[:, 3] - raw[:, 4]) # Z/Z_sun
    Z_star_hi = 10 ** (raw[:, 3] + raw[:, 5]) # Z/Z_sun

    # Define the scatter as offset from the mean value
    x_scatter = np.array((M_star - M_star_lo, M_star_hi - M_star))
    y_scatter = np.array((Z_star - Z_star_lo, Z_star_hi - Z_star))
    plt.errorbar(M_star, Z_star, yerr= y_scatter, xerr=x_scatter, fmt='',
                 ls='none',color='lightgreen',lw=1,ms=1,label='Kirby et al. (2013)')

def plot_gallazzi():
    # Cosmology
    h_sim = 0.6777
    h_obs = 0.704  # WMAP7
    Z_solar_obs = 0.02

    input_filename = "./plotter/obs_data/gallazzi_2021_ascii.txt"
    delimiter = "\t"

    # Read the data
    raw = np.loadtxt(input_filename, delimiter=delimiter)
    M_star = (
            10 ** raw[:, 0] * (h_sim / h_obs) ** -2
    )

    # Correction factor due to the difference in (X_O/X_Fe)_Sun
    # from Grevesse & Sauval (1993) to Asplund+ (2009)

    O_over_H_Grevesse93 = 8.83  # Grevesse & Sauval
    Fe_over_H_Grevesse93 = 7.5  # Grevesse & Sauval

    O_over_H_Andres89 = 8.93
    Fe_over_H_Andres89 = 7.51

    O_over_H_Asplund09 = 8.69
    Fe_over_H_Asplund09 = 7.50

    O_over_Fe_solar_Grevesse93 = O_over_H_Grevesse93 - Fe_over_H_Grevesse93
    O_over_Fe_solar_Andres89 = O_over_H_Andres89 - Fe_over_H_Andres89
    O_over_Fe_solar_Asplund09 = O_over_H_Asplund09 - Fe_over_H_Asplund09

    correction_Sun_O_over_Fe = O_over_Fe_solar_Grevesse93 - O_over_Fe_solar_Asplund09

    Z_median = (raw[:, 1] + correction_Sun_O_over_Fe)
    Z_lo = (raw[:, 2] + correction_Sun_O_over_Fe)
    Z_hi = (raw[:, 3] + correction_Sun_O_over_Fe)

    # Define the scatter as offset from the mean value
    y_scatter = np.array((Z_median - Z_lo, Z_hi - Z_median))
    plt.errorbar(M_star, Z_median, yerr= y_scatter, fmt='o',
                 ls='none',color='purple',lw=1,ms=1,label='Gallazzi et al. (2021)')


def plot_gallazzi_2005():
    # Cosmology
    h_sim = 0.6777
    h_obs = 0.7
    Z_solar_obs = 0.02
    kroupa_to_chabrier_mass = 0.912
    solar_metallicity = 0.0126

    input_filename = "./plotter/obs_data/Gallazzi_2005_ascii.txt"
    delimiter = "\t"

    # Read the data
    raw = np.loadtxt(input_filename, delimiter=delimiter)
    M_star = (
            10 ** raw[:, 0] * (h_sim / h_obs) ** -2 * kroupa_to_chabrier_mass
    )
    Z_median = 10 ** raw[:, 1] * Z_solar_obs / solar_metallicity
    Z_lo = 10 ** raw[:, 2] * Z_solar_obs / solar_metallicity
    Z_hi = 10 ** raw[:, 3] * Z_solar_obs / solar_metallicity

    # Define the scatter as offset from the mean value
    y_scatter = np.array((Z_median - Z_lo, Z_hi - Z_median))

    plt.errorbar(M_star, Z_median, yerr= y_scatter, fmt='o',
                 ls='none',color='grey',lw=1,ms=1,label='Gallazzi et al. (2005)', zorder=10)


def plot_Zahid_2017():

    mass_bins = np.array([8.55, 8.65, 8.75, 8.85, 8.95, 9.05, 9.15, 9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85,
                          9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65, 10.75, 10.85])

    input_filename = "./plotter/obs_data/Zahid_2017.txt"

    # Read the data
    raw = np.loadtxt(input_filename)
    M_star = raw[:, 1]
    Z_star = 10 ** raw[:, 3]

    Z_median = np.zeros(len(mass_bins))
    for i in range(len(mass_bins)):
        select = np.where(M_star == mass_bins[i])[0]
        Z_median[i] = np.median(Z_star[select])

    plt.plot(10**mass_bins, Z_median, '+',color='purple',ms=5,label='Zahid et al. (2017)')


def plot_Kudritzki_2016():

    input_filename = "./plotter/obs_data/Kudritzki_2016.txt"

    # Read the data
    raw = np.loadtxt(input_filename)
    M_star = 10 ** raw[:, 0]
    Z_median = 10 ** raw[:, 1]

    plt.plot(M_star, Z_median, 'v',color='black',ms=3,label='Kudritzki et al. (2016)')

def compute_metallicity_relation(sim_info, metallicity_data):

    # Look for abundance ratios from COLIBRE snaps:
    data = calculate_galaxies_metallicity(sim_info)

    Fe_H_median = data['Fe_H_median']
    O_Fe_median = data['O_Fe_median']
    Mg_Fe_median = data['Mg_Fe_median']

    Mstellar_median = data['Mstellar']

    Fe_H_mass_weighted = data['Fe_H_mass_weighted']
    O_Fe_mass_weighted = data['O_Fe_mass_weighted']
    Mg_Fe_mass_weighted = data['Mg_Fe_mass_weighted']

    Fe_H_ratio_weighted = data['Fe_H_ratio_weighted']
    O_Fe_ratio_weighted = data['O_Fe_ratio_weighted']
    Mg_Fe_ratio_weighted = data['Mg_Fe_ratio_weighted']

    Fe_H_light_weighted_r_band = data['Fe_H_light_weighted_r_band']
    Fe_H_light_weighted_i_band = data['Fe_H_light_weighted_i_band']
    Fe_H_light_weighted_Z_band = data['Fe_H_light_weighted_Z_band']

    Fe_H_light_weighted_r_band_with_mass = data['Fe_H_light_weighted_r_band_with_mass']
    # Fe_H_light_weighted_i_band_with_mass = data['Fe_H_light_weighted_i_band_with_mass']
    # Fe_H_light_weighted_Z_band_with_mass = data['Fe_H_light_weighted_Z_band_with_mass']

    Z_mass_weighted = data['Z_mass_weighted']
    Z_light_weighted_r_band = data['Z_light_weighted_r_band']
    Z_light_weighted_i_band = data['Z_light_weighted_i_band']
    Z_light_weighted_Z_band = data['Z_light_weighted_Z_band']

    O_Fe_light_weighted_r_band = data['O_Fe_light_weighted_r_band']
    Mg_Fe_light_weighted_r_band = data['Mg_Fe_light_weighted_r_band']
    O_Fe_light_weighted_r_band_with_mass = data['O_Fe_light_weighted_r_band_with_mass']
    Mg_Fe_light_weighted_r_band_with_mass = data['Mg_Fe_light_weighted_r_band_with_mass']

    counter = np.array([len(Mstellar_median)])

    if metallicity_data == None:

        metallicity_data = {
            'Fe_H_median': Fe_H_median,
            'O_Fe_median': O_Fe_median,
            'Mg_Fe_median': Mg_Fe_median,
            'Fe_H_mass_weighted': Fe_H_mass_weighted,
            'O_Fe_mass_weighted': O_Fe_mass_weighted,
            'Mg_Fe_mass_weighted': Mg_Fe_mass_weighted,

            'Fe_H_ratio_weighted': Fe_H_ratio_weighted,
            'O_Fe_ratio_weighted': O_Fe_ratio_weighted,
            'Mg_Fe_ratio_weighted': Mg_Fe_ratio_weighted,

            'Mstellar': Mstellar_median,
            'Z_mass_weighted': Z_mass_weighted,
            'Z_light_weighted_r_band': Z_light_weighted_r_band,
            'Z_light_weighted_i_band': Z_light_weighted_i_band,
            'Z_light_weighted_Z_band': Z_light_weighted_Z_band,

            'Fe_H_light_weighted_r_band': Fe_H_light_weighted_r_band,
            'Fe_H_light_weighted_i_band': Fe_H_light_weighted_i_band,
            'Fe_H_light_weighted_Z_band': Fe_H_light_weighted_Z_band,

            'O_Fe_light_weighted_r_band': O_Fe_light_weighted_r_band,
            'Mg_Fe_light_weighted_r_band': Mg_Fe_light_weighted_r_band,

            'Fe_H_light_weighted_r_band_with_mass': Fe_H_light_weighted_r_band_with_mass,
            'O_Fe_light_weighted_r_band_with_mass': O_Fe_light_weighted_r_band_with_mass,
            'Mg_Fe_light_weighted_r_band_with_mass': Mg_Fe_light_weighted_r_band_with_mass,

            'counter':counter}

    else:

        Fe_H = metallicity_data['Fe_H_median']
        Fe_H = np.append(Fe_H,Fe_H_median)
        O_Fe = metallicity_data['O_Fe_median']
        O_Fe = np.append(O_Fe,O_Fe_median)
        Mg_Fe = metallicity_data['Mg_Fe_median']
        Mg_Fe = np.append(Mg_Fe,Mg_Fe_median)

        Fe_H_t = metallicity_data['Fe_H_mass_weighted']
        Fe_H_t = np.append(Fe_H_t,Fe_H_mass_weighted)
        O_Fe_t = metallicity_data['O_Fe_mass_weighted']
        O_Fe_t = np.append(O_Fe_t,O_Fe_mass_weighted)
        Mg_Fe_t = metallicity_data['Mg_Fe_mass_weighted']
        Mg_Fe_t = np.append(Mg_Fe_t,Mg_Fe_mass_weighted)

        Fe_H_r = metallicity_data['Fe_H_ratio_weighted']
        Fe_H_r = np.append(Fe_H_r,Fe_H_ratio_weighted)
        O_Fe_r = metallicity_data['O_Fe_ratio_weighted']
        O_Fe_r = np.append(O_Fe_r,O_Fe_ratio_weighted)
        Mg_Fe_r = metallicity_data['Mg_Fe_ratio_weighted']
        Mg_Fe_r = np.append(Mg_Fe_r,Mg_Fe_ratio_weighted)

        Mstellar = metallicity_data['Mstellar']
        Mstellar = np.append(Mstellar, Mstellar_median)
        counter_sim = metallicity_data['counter']
        counter_sim = np.append(counter_sim, counter)

        Fe_H_lw_r = metallicity_data['Fe_H_light_weighted_r_band']
        Fe_H_lw_r = np.append(Fe_H_lw_r, Fe_H_light_weighted_r_band)
        Fe_H_lw_i = metallicity_data['Fe_H_light_weighted_i_band']
        Fe_H_lw_i = np.append(Fe_H_lw_i, Fe_H_light_weighted_i_band)
        Fe_H_lw_Z = metallicity_data['Fe_H_light_weighted_Z_band']
        Fe_H_lw_Z = np.append(Fe_H_lw_Z, Fe_H_light_weighted_Z_band)

        O_Fe_lw_r = metallicity_data['O_Fe_light_weighted_r_band']
        O_Fe_lw_r = np.append(O_Fe_lw_r, O_Fe_light_weighted_r_band)
        Mg_Fe_lw_r = metallicity_data['Mg_Fe_light_weighted_r_band']
        Mg_Fe_lw_r = np.append(Mg_Fe_lw_r, Mg_Fe_light_weighted_r_band)

        Fe_H_lw_r_with_mass = metallicity_data['Fe_H_light_weighted_r_band_with_mass']
        Fe_H_lw_r_with_mass = np.append(Fe_H_lw_r_with_mass, Fe_H_light_weighted_r_band_with_mass)
        O_Fe_lw_r_with_mass = metallicity_data['O_Fe_light_weighted_r_band_with_mass']
        O_Fe_lw_r_with_mass = np.append(O_Fe_lw_r_with_mass, O_Fe_light_weighted_r_band_with_mass)
        Mg_Fe_lw_r_with_mass= metallicity_data['Mg_Fe_light_weighted_r_band_with_mass']
        Mg_Fe_lw_r_with_mass = np.append(Mg_Fe_lw_r_with_mass, Mg_Fe_light_weighted_r_band_with_mass)

        Z_mw = metallicity_data['Z_mass_weighted']
        Z_mw = np.append(Z_mw, Z_mass_weighted)

        Z_lw_r = metallicity_data['Z_light_weighted_r_band']
        Z_lw_r = np.append(Z_lw_r, Z_light_weighted_r_band)
        Z_lw_i = metallicity_data['Z_light_weighted_i_band']
        Z_lw_i = np.append(Z_lw_i, Z_light_weighted_i_band)
        Z_lw_Z = metallicity_data['Z_light_weighted_Z_band']
        Z_lw_Z = np.append(Z_lw_Z, Z_light_weighted_Z_band)

        metallicity_data = {
            'Fe_H_median': Fe_H,
            'O_Fe_median': O_Fe,
            'Mg_Fe_median': Mg_Fe,
            'Fe_H_mass_weighted': Fe_H_t,
            'O_Fe_mass_weighted': O_Fe_t,
            'Mg_Fe_mass_weighted': Mg_Fe_t,

            'Fe_H_ratio_weighted': Fe_H_r,
            'O_Fe_ratio_weighted': O_Fe_r,
            'Mg_Fe_ratio_weighted': Mg_Fe_r,

            'Mstellar': Mstellar,
            'Z_mass_weighted': Z_mw,
            'Z_light_weighted_r_band': Z_lw_r,
            'Z_light_weighted_i_band': Z_lw_i,
            'Z_light_weighted_Z_band': Z_lw_Z,

            'Fe_H_light_weighted_r_band': Fe_H_lw_r,
            'Fe_H_light_weighted_i_band': Fe_H_lw_i,
            'Fe_H_light_weighted_Z_band': Fe_H_lw_Z,
            'O_Fe_light_weighted_r_band': O_Fe_lw_r,
            'Mg_Fe_light_weighted_r_band': Mg_Fe_lw_r,

            'Fe_H_light_weighted_r_band_with_mass': Fe_H_lw_r_with_mass,
            'O_Fe_light_weighted_r_band_with_mass': O_Fe_lw_r_with_mass,
            'Mg_Fe_light_weighted_r_band_with_mass': Mg_Fe_lw_r_with_mass,

            'counter': counter_sim}

    return metallicity_data
