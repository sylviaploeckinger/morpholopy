"""
Acknowledgements:
The routines making and plotting the KS relation have been written by Folkert Nobels,
while the specific scatter and project pixel routines were developed by Josh Borrow.
"""

from pylab import *
import numpy as np
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import scipy.stats as stat
from .loadObservationalData import read_obs_data


def project_pixel_grid(data, mode, res, region, rotation_matrix):

    x_min, x_max = region
    y_min, y_max = region
    x_range = x_max - x_min
    y_range = y_max - y_min

    if mode == 0: m = data[:,9] #H2
    if mode == 1: m = data[:,9]+data[:,8] #H2+HI
    if mode == 2: m = data[:,10] #SFR
    if mode == 3: m = data[:,8] #HI

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:,:3].T)

    x = (x - x_min) / x_range
    y = (y - y_min) / y_range

    image = np.zeros((res, res))
    maximal_array_index = res
    inverse_cell_area = res * res

    for x_pos, y_pos, mass in zip(x, y, m):
        particle_cell_x = int(res * x_pos)
        particle_cell_y = int(res * y_pos)

        if not (
                particle_cell_x < 0
                or particle_cell_x >= maximal_array_index
                or particle_cell_y < 0
                or particle_cell_y >= maximal_array_index
        ):
            image[particle_cell_x, particle_cell_y] += mass * inverse_cell_area
    return image


def integrate_metallicity(data, res, region, rotation_matrix):

    x_min, x_max = region
    y_min, y_max = region
    x_range = x_max - x_min
    y_range = y_max - y_min

    m = data[:, 12]

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:, :3].T)

    x = (x - x_min) / x_range
    y = (y - y_min) / y_range

    image = np.zeros((res, res))
    num_parts = np.zeros((res, res))
    maximal_array_index = res

    for x_pos, y_pos, mass in zip(x, y, m):
        particle_cell_x = int(res * x_pos)
        particle_cell_y = int(res * y_pos)

        if not (
                particle_cell_x < 0
                or particle_cell_x >= maximal_array_index
                or particle_cell_y < 0
                or particle_cell_y >= maximal_array_index
        ):
            image[particle_cell_x, particle_cell_y] += mass
            num_parts[particle_cell_x, particle_cell_y] += 1

    num_parts[num_parts==0] = 1 #lower value to avoid error
    image /= num_parts # Mean metallicity
    return image

def project_gas(data, mode, resolution, region, rotation_matrix):


    image = project_pixel_grid(data, mode, resolution, region, rotation_matrix)

    x_range = region[1] - region[0]
    y_range = region[1] - region[0]
    area = 1.0 / (x_range * y_range)
    image *= area
    return image

def bin_surface(radial_bins):
    """Returns the surface of the bins. """

    single_surface = lambda x: np.pi * x ** 2
    outer = single_surface(radial_bins[1:])
    inner = single_surface(radial_bins[:-1])
    return outer - inner


def project_gas_with_azimutal_average(data, mode, rotation_matrix):

    if mode == 0: m = data[:,9] #H2
    if mode == 1: m = data[:,9]+data[:,8] #HI+H2
    if mode == 2: m = data[:,10] #SFR
    if mode == 3: m = data[:,8] #HI

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:,:3].T)
    r = np.sqrt( x**2 + y**2)

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0, 4, 0.01)
    radial_bins = 10 ** radial_bins

    SumMode, _, _ = stat.binned_statistic(x=r, values=np.ones(len(r)) * m, statistic="sum", bins=radial_bins, )
    surface_density = (SumMode / bin_surface(radial_bins))  # Msun/kpc^2

    return surface_density


def KS_relation(data, ang_momentum, mode):

    size = 0.25 #kpc
    image_diameter = 60
    extent = [-30, 30]  #kpc
    number_of_pixels = int(image_diameter / size + 1)

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    # Calculate the SPH smoothed maps
    map_mass = project_gas(data, mode, number_of_pixels, extent, face_on_rotation_matrix)

    star_formation_rate_mask = data[:,10] > 0.0
    partsDATA = data[star_formation_rate_mask,:].copy()
    map_SFR = project_gas(partsDATA, 2, number_of_pixels, extent, face_on_rotation_matrix)

    map_metals = integrate_metallicity(data, number_of_pixels, extent, face_on_rotation_matrix)
    #map_mass = project_gas_with_azimutal_average(data, mode, face_on_rotation_matrix)
    #map_SFR = project_gas_with_azimutal_average(partsDATA, 2, face_on_rotation_matrix)

    # Bounds
    map_SFR[map_SFR <= 0] = 1e-6
    map_mass[map_mass <= 0] = 1e-6

    surface_density = np.log10(map_mass.flatten()) #Msun / kpc^2
    surface_density -= 6  #Msun / pc^2
    SFR_surface_density = np.log10(map_SFR.flatten()) #Msun / yr / kpc^2
    tgas = surface_density - SFR_surface_density + 6.

    return surface_density, SFR_surface_density, tgas, map_metals

def median_relations(x, y):

    xrange = np.arange(-1, 3, 0.1)

    xvalues = np.ones(len(xrange) - 1) * (-10)
    yvalues = np.zeros(len(xrange) - 1)
    yvalues_err_down = np.zeros(len(xrange) - 1)
    yvalues_err_up = np.zeros(len(xrange) - 1)

    perc = [16, 86]

    for i in range(0, len(xrange) - 2):
        mask = (x > xrange[i]) & (x < xrange[i + 1])
        if len(x[mask]) > 4:
            xvalues[i] = np.median(x[mask])
            yvalues[i] = np.median(y[mask])
            yvalues_err_down[i], yvalues_err_up[i] = np.transpose(np.percentile(y[mask], perc))

    mask = xvalues>-10
    xvalues = xvalues[mask]
    yvalues = yvalues[mask]
    yvalues_err_down = yvalues_err_down[mask]
    yvalues_err_up = yvalues_err_up[mask]

    return xvalues, yvalues, yvalues_err_down, yvalues_err_up


def KS_plots(data, ang_momentum, mode, galaxy_data, index, output_path):

    # read the observational data for the KS relations
    observational_data = read_obs_data("./plotter/obs_data")

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-1, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Get the surface densities
    surface_density, SFR_surface_density, tgas, metals = KS_relation(data, ang_momentum, mode)
    # Get median lines
    median_surface_density, median_SFR_surface_density, \
    SFR_surface_density_err_down, SFR_surface_density_err_up = median_relations(surface_density, SFR_surface_density)

    # Let's append data points to haloes for final plot at the end
    if mode ==1:
        galaxy_data.surface_density = np.append(galaxy_data.surface_density, surface_density)
        galaxy_data.SFR_density = np.append(galaxy_data.SFR_density, SFR_surface_density)
        galaxy_data.metallicity = np.append(galaxy_data.metallicity, metals)

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.8,
        "lines.markersize": 2,
        "lines.linewidth": 1.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    sfr_galaxy = galaxy_data.star_formation_rate[index]
    mass_galaxy = galaxy_data.stellar_mass[index]
    gas_mass_galaxy = galaxy_data.gas_mass[index]
    mass_halo = galaxy_data.halo_mass[index]
    galaxy_metallicity_gas_sfr = galaxy_data.metallicity_gas_sfr[index]
    galaxy_metallicity_gas = galaxy_data.metallicity_gas[index]

    title = r"$\log_{10}$ M$_{200}$/M$_{\odot} = $%0.2f," % (mass_halo)
    title += " SFR = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
    title += "\n Z$_{SFR>0}$ = %0.1f," % (galaxy_metallicity_gas_sfr)
    title += " Z$_{gas}$ = %0.1f," % (galaxy_metallicity_gas)
    title += "\n $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
    title += " $\&$ $\log_{10}$ M$_{gas}$/M$_{\odot} = $%0.2f" % (gas_mass_galaxy)

    ax.set_title(title)
    plt.plot(surface_density, SFR_surface_density, 'o', color='tab:blue')
    plt.plot(median_surface_density, median_SFR_surface_density, '-', color='black')
    plt.fill_between(median_surface_density, SFR_surface_density_err_down,
                     SFR_surface_density_err_up, alpha=0.2)
    plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), color="red", label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$", linestyle="--")
    plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
    plt.xlim(-1.0, 3.0)
    plt.ylim(-6.5, 1.0)

    if mode == 0:
        # load the observational data
        for ind, observation in enumerate(observational_data):
            if observation.gas_surface_density is not None:
                if (observation.description == "Bigiel et al. (2008) inner"):
                    data = observation.bin_data_KS_molecular(np.arange(-1, 3, .25), 0.4)
                    plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="o",ms=6,
                                 label=observation.description, color='tab:orange')

        plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend(labelspacing=0.2,handlelength=2,handletextpad=0.4,frameon=False)
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.savefig(f"{output_path}/KS_molecular_relation_%i.png" % (index))
        #np.savetxt(f"{output_path}/KS_molecular_relation_file_{snapshot_number:04d}.txt", np.transpose(
        #    [surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up, tgas,
        #     tgas_err_down, tgas_err_up]))

    elif mode == 1:
        # load the observational data
        for ind, observation in enumerate(observational_data):
            if observation.gas_surface_density is not None:
                if (observation.description == "Bigiel et al. (2008) inner"):
                    data = observation.bin_data_KS(np.arange(-1, 3, .25), 0.4)
                    plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="o",ms=6,
                                 label=observation.description, color='tab:green')
                elif (observation.description == "Bigiel et al. (2010) outer"):
                    data2 = observation.bin_data_KS(np.arange(-1, 3, .25), 0.4)
                    plt.errorbar(data2[0], data2[1], yerr=[data2[2], data2[3]], fmt="o",ms=6,
                                 label=observation.description, color='tab:orange')

        plt.xlabel("log $\\Sigma_{HI}+ \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend(labelspacing=0.2,handlelength=2,handletextpad=0.4,frameon=False)
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.savefig(f"{output_path}/KS_relation_best_%i.png" % (index))
        #np.savetxt(f"{output_path}/KS_relation_best_file_{snapshot_number:04d}.txt", np.transpose(
        #    [surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up, tgas,
        #     tgas_err_down, tgas_err_up]))
    plt.close()

    median_surface_density, median_tgas, tgas_err_down, tgas_err_up = median_relations(surface_density, tgas)

    figure()
    ax = plt.subplot(1, 1, 1)
    sfr_galaxy = galaxy_data.star_formation_rate[index]
    mass_galaxy = galaxy_data.stellar_mass[index]
    gas_mass_galaxy = galaxy_data.gas_mass[index]
    mass_halo = galaxy_data.halo_mass[index]
    galaxy_metallicity_gas_sfr = galaxy_data.metallicity_gas_sfr[index]
    galaxy_metallicity_gas = galaxy_data.metallicity_gas[index]

    title = r"$\log_{10}$ M$_{200}$/M$_{\odot} = $%0.2f," % (mass_halo)
    title += " SFR = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
    title += "\n Z$_{SFR>0}$ = %0.1f," % (galaxy_metallicity_gas_sfr)
    title += " Z$_{gas}$ = %0.1f," % (galaxy_metallicity_gas)
    title += "\n $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
    title += " $\&$ $\log_{10}$ M$_{gas}$/M$_{\odot} = $%0.2f" % (gas_mass_galaxy)
    ax.set_title(title)

    plt.plot(surface_density, tgas, 'o', color='tab:blue')
    plt.plot(median_surface_density, median_tgas, '-', color='black')
    plt.fill_between(median_surface_density, tgas_err_down,  tgas_err_up, alpha=0.2)
    plt.plot(np.log10(Sigma_g),np.log10(Sigma_g)- np.log10(Sigma_star)+6.,color="red",
             label="KS law (Kennicutt 98)",linestyle="--")
    plt.xlim(-1,3.0)
    plt.ylim(7, 12)

    if mode == 0:
        # load the observational data
        for ind, observation in enumerate(observational_data):
            if observation.gas_surface_density is not None:
                if (observation.description == "Bigiel et al. (2008) inner"):
                    data = observation.bin_data_gas_depletion_molecular(np.arange(-1, 3, .25), 0.4)
                    plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="o",ms=6,
                                 label=observation.description, color='tab:orange')

        plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        #plt.legend(labelspacing=0.2,handlelength=2,handletextpad=0.4,frameon=False)
        plt.ylabel("log $\\rm t_{gas} = \\Sigma_{H_2} / \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.savefig(f"{output_path}/molecular_gas_depletion_timescale_%i.png" % (index))

    elif mode == 1:
        # load the observational data
        for ind, observation in enumerate(observational_data):
            if observation.gas_surface_density is not None:
                if (observation.description == "Bigiel et al. (2008) inner"):
                    data = observation.bin_data_gas_depletion(np.arange(-1, 3, .25), 0.4)
                    plt.errorbar(data[0], data[1], yerr=[data[2], data[3]], fmt="o",ms=6,
                                 label=observation.description, color='tab:green')
                elif (observation.description == "Bigiel et al. (2010) outer"):
                    data2 = observation.bin_data_gas_depletion(np.arange(-1, 3, .25), 0.4)
                    plt.errorbar(data2[0], data2[1], yerr=[data2[2], data2[3]], fmt="o",ms=6,
                                 label=observation.description, color='tab:orange')

        plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        #plt.legend(labelspacing=0.2,handlelength=2,handletextpad=0.4,frameon=False)
        plt.ylabel("log $\\rm t_{gas} = (\\Sigma_{HI} + \\Sigma_{H_2} )/ \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.savefig(f"{output_path}/gas_depletion_timescale_best_%i.png" % (index))

def surface_ratios(data, ang_momentum):

    size = 0.25 #kpc
    image_diameter = 60
    extent = [-30, 30]  #kpc
    number_of_pixels = int(image_diameter / size + 1)

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    # Calculate the SPH smoothed maps
    map_H2 = project_gas(data, 0, number_of_pixels, extent, face_on_rotation_matrix)
    map_HI = project_gas(data, 3, number_of_pixels, extent, face_on_rotation_matrix)
    map_gas = map_H2 + map_HI

    # Bounds
    map_H2[map_H2 <= 0] = 1e-6
    map_gas[map_gas <= 0] = 1e-6
    ratio = map_H2 / map_gas

    surface_density = np.log10(map_gas.flatten()) #HI+H2 Msun / kpc^2
    #H2_surface_density = np.log10(map_H2.flatten()) #Msun / kpc^2
    #ratio_density = H2_surface_density - surface_density # no units
    ratio_density = np.log10(ratio.flatten()) #no units
    surface_density -= 6  #HI+H2 Msun / pc^2
    return surface_density, ratio_density

def Krumholz_eq39(Sigma_neutral, f):
    Zprime = 1.0
    psi = 1.6  # fiducial from Krumholz
    s = 1. / f * Sigma_neutral * Zprime / psi
    RH2 = np.power((1. + np.power(s / 11., 3) * np.power((125. + s) / (96. + s), 3)), 1. / 3.) - 1.
    return RH2

def make_surface_density_ratios(data, ang_momentum, galaxy_data, index, output_path):

    # Get the surface densities
    Sigma_gas, Sigma_ratio = surface_ratios(data, ang_momentum)
    Median_Sigma_gas, Median_Sigma_ratio, Sigma_ratio_err_down, \
    Sigma_ratio_err_up = median_relations(Sigma_gas, Sigma_ratio)

    galaxy_data.ratio_densities = np.append(galaxy_data.ratio_densities, Sigma_ratio)

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.8,
        "lines.markersize": 2,
        "lines.linewidth": 1.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)

    sfr_galaxy = galaxy_data.star_formation_rate[index]
    mass_galaxy = galaxy_data.stellar_mass[index]
    gas_mass_galaxy = galaxy_data.gas_mass[index]
    mass_halo = galaxy_data.halo_mass[index]
    galaxy_metallicity_gas_sfr = galaxy_data.metallicity_gas_sfr[index]
    galaxy_metallicity_gas = galaxy_data.metallicity_gas[index]

    title = r"$\log_{10}$ M$_{200}$/M$_{\odot} = $%0.2f," % (mass_halo)
    title += " SFR = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
    title += "\n Z$_{SFR>0}$ = %0.1f," % (galaxy_metallicity_gas_sfr)
    title += " Z$_{gas}$ = %0.1f" % (galaxy_metallicity_gas)
    title += "\n $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
    title += " $\&$ $\log_{10}$ M$_{gas}$/M$_{\odot} = $%0.2f" % (gas_mass_galaxy)
    ax.set_title(title)

    plt.plot(Sigma_gas, Sigma_ratio, 'o', color='tab:blue')
    plt.plot(Median_Sigma_gas, Median_Sigma_ratio, '-', color='black')
    plt.fill_between(Median_Sigma_gas, Sigma_ratio_err_down, Sigma_ratio_err_up, alpha=0.2)
    plt.ylabel(r"log $\Sigma_{\mathrm{H2}} / (\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}})$")
    plt.xlabel(r"log $\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}}$  [M$_{\odot}$ pc$^{-2}$]")

    # Krumholz 2009 lines
    Sigma_neutral = np.arange(-1, 3, 0.2)
    RH2 = 1. / Krumholz_eq39(10**Sigma_neutral, 0.5)
    FH2 = np.log10(1. / (1. + RH2))
    plt.plot(Sigma_neutral, FH2, '--', color='tab:red', label="Krumholz+ (2009): f = 0.5")

    plt.xlim(-1.0, 3.0)
    plt.ylim(-8.0, 0.5)
    plt.legend(labelspacing=0.2, handlelength=2, handletextpad=0.4, frameon=False)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/Surface_density_ratio_%i.png" % (index))
    plt.close()

def calculate_integrated_quantities(data, ang_momentum, radius, mode):

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    x, y, _ = np.matmul(face_on_rotation_matrix, data[:, :3].T)
    r = np.sqrt(x ** 2 + y ** 2)
    select = r <= radius

    surface = np.pi * radius**2
    if mode == 0: m = data[select,9]
    if mode == 1: m = data[select,9]+data[select,8]

    # If we have gas within rhalfMs
    if len(m)>0:
        Sigma_gas = np.log10( np.sum(m)  / surface ) - 6. #Msun / pc^2

        sfr = data[select,10]
        sfr = sfr[sfr>0]
        Sigma_SFR = np.log10( np.sum(sfr) / surface ) #Msun / yr / kpc^2

    else:
        Sigma_gas = -6
        Sigma_SFR = -6

    return Sigma_gas, Sigma_SFR

def make_KS_plots(data, ang_momentum, galaxy_data, index, KSPlotsInWeb, output_path):


    for mode, project in enumerate(["molecular_hydrogen_masses", "not_ionized_hydrogen_masses"]):

        KS_plots(data, ang_momentum, mode, galaxy_data, index, output_path)

        if mode == 0:
            outfile = "KS_molecular_relation_%i.png" % (index)
            title = "KS relation (H2 mass)"
            id = abs(hash("galaxy KS relation H2 %i" % (index)))
        if mode == 1:
            outfile = "KS_relation_best_%i.png" % (index)
            title = "KS relation (H2+HI mass)"
            id = abs(hash("galaxy KS relation H2+HI %i" % (index)))

        caption = "KS relation."
        KSPlotsInWeb.load_plots(title, caption, outfile, id)

        if mode == 0:
            title = "Depletion time (H2 mass)"
            id = abs(hash("galaxy depletion H2 %i" % (index)))
            outfile = "molecular_gas_depletion_timescale_%i.png" % (index)
            caption = "Gas depletion times. The observational data-points correspond to"
            caption += " Bigiel et al. (2008) inner, same as KS relation (H2 mass) figure."
        if mode == 1:
            title = "Depletion time (H2+HI mass)"
            id = abs(hash("galaxy depletion H2+HI %i" % (index)))
            outfile = "gas_depletion_timescale_best_%i.png" % (index)
            caption = "Gas depletion times. The observational data-points correspond to"
            caption += " Bigiel et al. (2008, 2010), same as KS relation (H2+HI mass) figure."
        KSPlotsInWeb.load_plots(title, caption, outfile, id)


    make_surface_density_ratios(data, ang_momentum, galaxy_data, index, output_path)

    title = "Surface density ratios"
    id = abs(hash("density ratio H2+HI %i" % (index)))
    outfile = "Surface_density_ratio_%i.png" % (index)
    caption = "Surface density ratios."
    KSPlotsInWeb.load_plots(title, caption, outfile, id)


def calculate_surface_densities(data, ang_momentum, galaxy_data, index):

    # If we have gas, calculate ..
    radius = galaxy_data.halfmass_radius_star[index]

    # Mode ==0 : "molecular_hydrogen_masses"
    # Mode ==1 : "not_ionized_hydrogen_masses"
    Sigma_H2, Sigma_SFR_H2 = calculate_integrated_quantities(data, ang_momentum, radius, 0)
    Sigma_gas, Sigma_SFR = calculate_integrated_quantities(data, ang_momentum, radius, 1)
    Sigma = np.array([Sigma_H2, Sigma_gas, Sigma_SFR])
    galaxy_data.add_surface_density(Sigma, index)
