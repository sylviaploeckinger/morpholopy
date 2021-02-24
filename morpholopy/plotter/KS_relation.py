from pylab import *
import numpy as np
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import scipy.stats as stat


def scatter(x, y, m, res):

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

def project_pixel_grid(data, mode, res, region, rotation_matrix):

    x_min, x_max = region
    y_min, y_max = region
    x_range = x_max - x_min
    y_range = y_max - y_min

    if mode == 0: m = data[:,9]
    if mode == 1: m = data[:,9]+data[:,8]
    if mode == 2: m = data[:,10]

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:,:3].T)

    #binsize = 0.25  # kpc
    #edges = (np.arange(161) - (160 / 2.)) * binsize
    #Hmass = np.histogram2d(x, y, bins=(edges, edges), normed=False, weights=m)

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
    #return Hmass[0]

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

    if mode == 0: m = data[:,9]
    if mode == 1: m = data[:,9]+data[:,8]
    if mode == 2: m = data[:,10]

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

    size = 0.65 #kpc
    image_diameter = 60
    extent = [-30, 30]  #kpc
    number_of_pixels = int(image_diameter / size + 1)

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    # Calculate the SPH smoothed maps
    map_mass = project_gas(data, mode, number_of_pixels, extent, face_on_rotation_matrix)

    star_formation_rate_mask = data[:,10] > 0.0
    partsDATA = data[star_formation_rate_mask,:].copy()
    map_SFR = project_gas(partsDATA, 2, number_of_pixels, extent, face_on_rotation_matrix)

    #map_mass = project_gas_with_azimutal_average(data, mode, face_on_rotation_matrix)
    #map_SFR = project_gas_with_azimutal_average(partsDATA, 2, face_on_rotation_matrix)

    # Bounds
    map_SFR[map_SFR <= 0] = 1e-6
    map_mass[map_mass <= 0] = 1e-6

    surface_density = np.log10(map_mass.flatten()) #Msun / kpc^2
    surface_density -= 6  #Msun / pc^2
    SFR_surface_density = np.log10(map_SFR.flatten()) #Msun / yr / kpc^2
    tgas = surface_density - SFR_surface_density + 6.

    surface_range = np.arange(-1, 7, .25)

    SFR_values = np.zeros(len(surface_range)-1)
    SFR_values_err_down = np.zeros(len(surface_range)-1)
    SFR_values_err_up = np.zeros(len(surface_range)-1)
    tgas_values = np.zeros(len(surface_range) - 1)
    tgas_values_err_down = np.zeros(len(surface_range) - 1)
    tgas_values_err_up = np.zeros(len(surface_range) - 1)

    perc = [16, 86]

    for i in range(0, len(surface_range) - 2):
        mask = (surface_density > surface_range[i]) & (surface_density < surface_range[i + 1])
        surface = SFR_surface_density[mask]
        if len(surface)>0:
            SFR_values[i] = np.median(SFR_surface_density[mask])
            tgas_values[i] = np.median(tgas[mask])
        try:
            SFR_values_err_down[i], SFR_values_err_up[i] = np.transpose(np.percentile(SFR_surface_density[mask], perc))
            tgas_values_err_down[i], tgas_values_err_up[i] = np.transpose(np.percentile(tgas[mask], perc))
        except:
            SFR_values_err_down[i], SFR_values_err_up[i] = [0., 0.]
            tgas_values_err_down[i], tgas_values_err_up[i] = [0., 0.]

    plot_surface_range = (surface_range[1:] + surface_range[:-1]) / 2.

    SFR_values_err_down = np.abs(SFR_values_err_down - SFR_values)
    SFR_values_err_up = np.abs(SFR_values_err_up - SFR_values)
    tgas_values_err_down = np.abs(tgas_values_err_down - tgas_values)
    tgas_values_err_up = np.abs(tgas_values_err_up - tgas_values)

    return plot_surface_range, SFR_values, SFR_values_err_down, \
           SFR_values_err_up, tgas_values, tgas_values_err_down, \
           tgas_values_err_up


def make_KS_plots(data, ang_momentum, mode, galaxy_data, index, output_path):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-1, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Get the surface densities
    surface_density, SFR_surface_density, \
    SFR_surface_density_err_down, SFR_surface_density_err_up, \
    tgas, tgas_err_down, tgas_err_up = KS_relation(data, ang_momentum, mode)


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
        "lines.markersize": 6,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    sfr_galaxy = galaxy_data.star_formation_rate[index]
    mass_galaxy = galaxy_data.stellar_mass[index]
    gas_mass_galaxy = galaxy_data.gas_mass[index]
    title = r"Star formation rate = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
    title += "\n $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
    title += " $\&$ $\log_{10}$ M$_{gas}$/M$_{\odot} = $%0.2f" % (gas_mass_galaxy)
    ax.set_title(title)
    plt.plot(surface_density, SFR_surface_density)
    plt.fill_between(surface_density, SFR_surface_density - SFR_surface_density_err_down,
                     SFR_surface_density + SFR_surface_density_err_up, alpha=0.2)
    plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), color="red", label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$", linestyle="--")
    plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
    plt.xlim(-1.0, 3.0)
    plt.ylim(-7.0, 1.0)

    if mode == 0:
        plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend()
        plt.savefig(f"{output_path}/KS_molecular_relation_%i.png" % (index))
        #np.savetxt(f"{output_path}/KS_molecular_relation_file_{snapshot_number:04d}.txt", np.transpose(
        #    [surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up, tgas,
        #     tgas_err_down, tgas_err_up]))

    elif mode == 1:
        plt.xlabel("log $\\Sigma_{HI}+ \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend()
        plt.savefig(f"{output_path}/KS_relation_best_%i.png" % (index))
        #np.savetxt(f"{output_path}/KS_relation_best_file_{snapshot_number:04d}.txt", np.transpose(
        #    [surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up, tgas,
        #     tgas_err_down, tgas_err_up]))
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    title = r"Star formation rate = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
    title += "\n $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
    title += " $\&$ $\log_{10}$ M$_{gas}$/M$_{\odot} = $%0.2f" % (gas_mass_galaxy)
    ax.set_title(title)

    plt.plot(surface_density, tgas)
    plt.fill_between(surface_density, tgas-tgas_err_down,  tgas+tgas_err_up,alpha=0.2)
    plt.plot(np.log10(Sigma_g),np.log10(Sigma_g)- np.log10(Sigma_star)+6.,color="red",
             label="KS law (Kennicutt 98)",linestyle="--")
    plt.xlim(0,3.0)
    plt.ylim(2, 14)

    if mode == 0:
        plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend()
        plt.ylabel("log $\\rm t_{gas} = \\Sigma_{H_2} / \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$")
        plt.savefig(f"{output_path}/molecular_gas_depletion_timescale_%i.png" % (index))

    elif mode == 1:
        plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend()
        plt.ylabel("log $\\rm t_{gas} = (\\Sigma_{HI} + \\Sigma_{H_2} )/ \\Sigma_{\\rm SFR}$ $[{\\rm yr }]$")
        plt.savefig(f"{output_path}/gas_depletion_timescale_best_%i.png" % (index))

def surface_ratios(data, ang_momentum):

    size = 0.65 #kpc
    image_diameter = 60
    extent = [-30, 30]  #kpc
    number_of_pixels = int(image_diameter / size + 1)

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    # Calculate the SPH smoothed maps
    map_H2 = project_gas(data, 0, number_of_pixels, extent, face_on_rotation_matrix)
    map_gas = project_gas(data, 1, number_of_pixels, extent, face_on_rotation_matrix)

    # Bounds
    map_H2[map_H2 <= 0] = 1e-6
    map_gas[map_gas <= 0] = 1e-6

    surface_density = np.log10(map_gas.flatten()) #Msun / kpc^2
    H2_surface_density = np.log10(map_H2.flatten()) #Msun / kpc^2
    ratio_density = H2_surface_density - surface_density # no units
    surface_density -= 6  #Msun / pc^2

    surface_range = np.arange(-1, 7, .25)

    Ratio_values = np.zeros(len(surface_range)-1)
    Ratio_err_down = np.zeros(len(surface_range)-1)
    Ratio_err_up = np.zeros(len(surface_range)-1)

    perc = [16, 86]

    for i in range(0, len(surface_range) - 2):
        mask = (surface_density > surface_range[i]) & (surface_density < surface_range[i + 1])
        surface = ratio_density[mask]
        if len(surface)>0:
            Ratio_values[i] = np.median(ratio_density[mask])
        try:
            Ratio_err_down[i], Ratio_err_up[i] = np.transpose(np.percentile(ratio_density[mask], perc))
        except:
            Ratio_err_down[i], Ratio_err_up[i] = [0., 0.]

    plot_surface_range = (surface_range[1:] + surface_range[:-1]) / 2.
    Ratio_err_down = np.abs(Ratio_err_down - Ratio_values)
    Ratio_err_up = np.abs(Ratio_err_up - Ratio_values)

    return plot_surface_range, Ratio_values, Ratio_err_down, Ratio_err_up

def make_surface_density_ratios(data, ang_momentum, galaxy_data, index, output_path):

    # Get the surface densities
    Sigma_gas, Sigma_ratio, \
    Sigma_ratio_err_down, Sigma_ratio_err_up = surface_ratios(data, ang_momentum)


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
        "lines.markersize": 6,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    sfr_galaxy = galaxy_data.star_formation_rate[index]
    mass_galaxy = galaxy_data.stellar_mass[index]
    gas_mass_galaxy = galaxy_data.gas_mass[index]
    title = r"Star formation rate = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
    title += "\n $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
    title += " $\&$ $\log_{10}$ M$_{gas}$/M$_{\odot} = $%0.2f" % (gas_mass_galaxy)
    ax.set_title(title)

    plt.plot(Sigma_gas, Sigma_ratio)
    plt.fill_between(Sigma_gas, Sigma_ratio - Sigma_ratio_err_down,
                     Sigma_ratio + Sigma_ratio_err_up, alpha=0.2)
    plt.ylabel(r"log $\Sigma_{\mathrm{H2}} / (\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}})$")
    plt.xlabel(r"log $\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}}$  [M$_{\odot}$ pc$^{-2}$]")

    plt.xlim(-1.0, 3.0)
    plt.ylim(-7.0, 1.0)

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
    Sigma_gas = np.log10( np.sum(m)  / surface ) - 6. #Msun / pc^2

    sfr = data[select,10]
    sfr = sfr[sfr>0]
    Sigma_SFR = np.log10( np.sum(sfr) / surface ) #Msun / yr / kpc^2

    return Sigma_gas, Sigma_SFR

def KS_plots(data, ang_momentum, galaxy_data, index, KSPlotsInWeb, output_path):

    Sigma_gas = np.zeros(2)
    Sigma_SFR = np.zeros(2)
    radius = galaxy_data.halfmass_radius_star[index]

    for mode, project in enumerate(["molecular_hydrogen_masses", "not_ionized_hydrogen_masses"]):

        make_KS_plots(data, ang_momentum, mode, galaxy_data, index, output_path)
        Sigma_gas[mode], Sigma_SFR[mode] = calculate_integrated_quantities(data, ang_momentum, radius, mode)

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
        if mode == 1:
            title = "Depletion time (H2+HI mass)"
            id = abs(hash("galaxy depletion H2+HI %i" % (index)))
            outfile = "gas_depletion_timescale_best_%i.png" % (index)

        caption = "Gas depletion times."
        KSPlotsInWeb.load_plots(title, caption, outfile, id)


    make_surface_density_ratios(data, ang_momentum, galaxy_data, index, output_path)

    title = "Surface density ratios"
    id = abs(hash("density ratio H2+HI %i" % (index)))
    outfile = "Surface_density_ratio_%i.png" % (index)
    caption = "Surface density ratios."
    KSPlotsInWeb.load_plots(title, caption, outfile, id)

    return Sigma_gas, Sigma_SFR