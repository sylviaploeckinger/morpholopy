from pylab import *
import numpy as np
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector


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
    #binsize = 0.5
    #area = 1.0 / binsize**2
    #image *= area
    return image

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

    # Bounds
    map_SFR[map_SFR <= 0] = 1e-6
    map_mass[map_mass == 0] = 1e-6

    surface_density = np.log10(map_mass.flatten()) #Msun / kpc^2
    surface_density -= 6  #Msun / pc^2
    SFR_surface_density = np.log10(map_SFR.flatten()) #Msun / yr / kpc^2

    surface_range = np.arange(-1, 7, .25)

    SFR_values = np.zeros(len(surface_range)-1)
    SFR_values_err_down = np.zeros(len(surface_range)-1)
    SFR_values_err_up = np.zeros(len(surface_range)-1)

    perc = [16, 86]

    for i in range(0, len(surface_range) - 2):
        mask = (surface_density > surface_range[i]) & (surface_density < surface_range[i + 1])
        SFR_values[i] = np.median(SFR_surface_density[mask])
        try:
            SFR_values_err_down[i], SFR_values_err_up[i] = np.transpose(np.percentile(SFR_surface_density[mask], perc))
        except:
            SFR_values_err_down[i], SFR_values_err_up[i] = [0., 0.]

    plot_surface_range = (surface_range[1:] + surface_range[:-1]) / 2.

    SFR_values_err_down = np.abs(SFR_values_err_down - SFR_values)
    SFR_values_err_up = np.abs(SFR_values_err_up - SFR_values)

    return plot_surface_range, SFR_values, SFR_values_err_down, SFR_values_err_up


def make_KS_plots(data, ang_momentum, mode, index, output_path):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-1, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Get the surface densities
    surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up = KS_relation(
        data, ang_momentum, mode)

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
        "lines.markersize": 6,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)

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



def KS_plots(data, ang_momentum, index, KSPlotsInWeb, output_path):

    for mode, project in enumerate(["molecular_hydrogen_masses", "not_ionized_hydrogen_masses"]):

            make_KS_plots(data, ang_momentum, mode, index, output_path)

            if mode == 0:
                outfile = f"{output_path}/KS_molecular_relation_%i.png" % (index)
                title = "Galaxy %i / KS relation (H2 mass)" % (index)
                id = abs(hash("galaxy KS relation H2 %i" % (index)))
            if mode == 1:
                title = "Galaxy %i / KS relation (H2+HI mass)" % (index)
                id = abs(hash("galaxy KS relation H2+HI %i" % (index)))
                outfile = f"{output_path}/KS_relation_best_%i.png" % (index)

            caption = "KS relation."
            KSPlotsInWeb.load_plots(title, caption, outfile, id)
