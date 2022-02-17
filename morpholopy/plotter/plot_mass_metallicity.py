import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm
import scipy.stats as stat


def compute_ratios(Oxygen_fraction, Iron_fraction,
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

def compute_total_ratios(Oxygen_fraction, Iron_fraction,
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

    return {'Fe_H': Fe_H, 'O_Fe': O_Fe, 'Mg_Fe': Mg_Fe}


def calculate_galaxies_metallicity(sim_info):

    select_sample = np.where((sim_info.halo_data.log10_stellar_mass >= 7.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 12))[0]

    select_centrals = np.where(sim_info.halo_data.type[select_sample] == 10)[0]

    sample = select_sample[select_centrals]
    num_sample = len(sample)

    FeH = np.zeros(num_sample)
    MgFe = np.zeros(num_sample)
    OFe = np.zeros(num_sample)
    totalFeH = np.zeros(num_sample)
    totalMgFe = np.zeros(num_sample)
    totalOFe = np.zeros(num_sample)

    stellar_mass = sim_info.halo_data.log10_stellar_mass[sample]
    metallicity = sim_info.halo_data.metallicity[sample]

    for i in tqdm(range(len(sample))):
        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        O_stars = sim_info.stars.oxygen.copy()
        Fe_stars = sim_info.stars.iron.copy()
        Mg_stars = sim_info.stars.magnesium.copy()
        H_stars = sim_info.stars.hydrogen.copy()

        ratios = compute_ratios(O_stars, Fe_stars, Mg_stars, H_stars)
        FeH[i] = ratios['Fe_H']
        MgFe[i] = ratios['Mg_Fe']
        OFe[i] = ratios['O_Fe']

        O_stars = np.sum( sim_info.stars.oxygen * sim_info.stars.mass )
        Fe_stars = np.sum( sim_info.stars.iron * sim_info.stars.mass )
        Mg_stars = np.sum( sim_info.stars.magnesium * sim_info.stars.mass )
        H_stars = np.sum( sim_info.stars.hydrogen * sim_info.stars.mass )

        ratios = compute_total_ratios(O_stars, Fe_stars, Mg_stars, H_stars)
        totalFeH[i] = ratios['Fe_H']
        totalMgFe[i] = ratios['Mg_Fe']
        totalOFe[i] = ratios['O_Fe']


    return {'Fe_H': FeH, 'O_Fe': OFe, 'Mg_Fe': MgFe,
            'Fe_H_total': totalFeH, 'O_Fe_total': totalOFe, 'Mg_Fe_total': totalMgFe,
            'Z':metallicity, 'Mstellar':stellar_mass}

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
    plt.axis([2e-4, 2e0, 1e-3, 1e0])
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5,
               handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

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

    input_filename = "./plotter/obs_data/Gallazzi_2021_ascii.txt"
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


def compute_metallicity_relation(sim_info, metallicity_data):

    # Look for abundance ratios from COLIBRE snaps:
    data = calculate_galaxies_metallicity(sim_info)
    Fe_H_median = data['Fe_H']
    O_Fe_median = data['O_Fe']
    Mg_Fe_median = data['Mg_Fe']
    Mstellar_median = data['Mstellar']
    Fe_H_total = data['Fe_H_total']
    O_Fe_total = data['O_Fe_total']
    Mg_Fe_total = data['Mg_Fe_total']
    counter = np.array([len(Mstellar_median)])

    if metallicity_data == None:

        metallicity_data = {
            'Fe_H': Fe_H_median,
            'O_Fe': O_Fe_median,
            'Mg_Fe': Mg_Fe_median,
            'Fe_H_total': Fe_H_total,
            'O_Fe_total': O_Fe_total,
            'Mg_Fe_total': Mg_Fe_total,
            'Mstellar': Mstellar_median,
            'counter':counter}

    else:

        Fe_H = metallicity_data['Fe_H']
        Fe_H = np.append(Fe_H,Fe_H_median)
        O_Fe = metallicity_data['O_Fe']
        O_Fe = np.append(O_Fe,O_Fe_median)
        Mg_Fe = metallicity_data['Mg_Fe']
        Mg_Fe = np.append(Mg_Fe,Mg_Fe_median)

        Fe_H_t = metallicity_data['Fe_H_total']
        Fe_H_t = np.append(Fe_H_t,Fe_H_total)
        O_Fe_t = metallicity_data['O_Fe_total']
        O_Fe_t = np.append(O_Fe_t,O_Fe_total)
        Mg_Fe_t = metallicity_data['Mg_Fe_total']
        Mg_Fe_t = np.append(Mg_Fe_t,Mg_Fe_total)

        Mstellar = metallicity_data['Mstellar']
        Mstellar = np.append(Mstellar, Mstellar_median)
        counter_sim = metallicity_data['counter']
        counter_sim = np.append(counter_sim, counter)

        metallicity_data = {
            'Fe_H': Fe_H,
            'O_Fe': O_Fe,
            'Mg_Fe': Mg_Fe,
            'Fe_H_total': Fe_H_t,
            'O_Fe_total': O_Fe_t,
            'Mg_Fe_total': Mg_Fe_t,
            'Mstellar': Mstellar,
            'counter': counter_sim}

    return metallicity_data