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

    if mode == 0:
        m = data[:, 9]  # H2
    if mode == 1:
        m = data[:, 9] + data[:, 8]  # H2+HI
    if mode == 2:
        m = data[:, 10]  # SFR
    if mode == 3:
        m = data[:, 8]  # HI

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:, :3].T)

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


def integrate_metallicity_using_grid(data, res, region, rotation_matrix):

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

    num_parts[num_parts == 0] = 1  # lower value to avoid error
    image /= num_parts  # Mean metallicity
    return image


def project_gas(data, mode, resolution, region, rotation_matrix):

    image = project_pixel_grid(data, mode, resolution, region, rotation_matrix)

    x_range = region[1] - region[0]
    y_range = region[1] - region[0]
    area = 1.0 / (x_range * y_range)
    image *= area
    return image


def bin_surface(radial_bins):
    """Returns the surface of the bins."""

    single_surface = lambda x: np.pi * x ** 2
    outer = single_surface(radial_bins[1:])
    inner = single_surface(radial_bins[:-1])
    return outer - inner


def project_metals_with_azimuthal_average(data, rotation_matrix, bin_size):
    """Returns the mean gas metallicity from each concentric shell"""

    m = data[:, 12]

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:, :3].T)
    r = np.sqrt(x ** 2 + y ** 2)

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0, 40, bin_size)
    MeanMetals, _, _ = stat.binned_statistic(
        x=r,
        values=m,
        statistic="mean",
        bins=radial_bins,
    )

    return MeanMetals


def project_gas_with_azimuthal_average(data, mode, rotation_matrix, bin_size):

    if mode == 0:
        m = data[:, 9]  # H2
    if mode == 1:
        m = data[:, 9] + data[:, 8]  # HI+H2
    if mode == 2:
        m = data[:, 10]  # SFR
    if mode == 3:
        m = data[:, 8]  # HI

    # Rotate co-ordinates as required
    x, y, _ = np.matmul(rotation_matrix, data[:, :3].T)
    r = np.sqrt(x ** 2 + y ** 2)

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0, 30, bin_size)
    SumMode, _, _ = stat.binned_statistic(
        x=r,
        values=m,
        statistic="sum",
        bins=radial_bins,
    )
    surface_density = SumMode / bin_surface(radial_bins)  # Msun/kpc^2

    return surface_density


def KS_relation(data, ang_momentum, mode, method, size):

    image_diameter = 60
    extent = [-30, 30]  # kpc
    number_of_pixels = int(image_diameter / size + 1)

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    if method == "grid":
        # Calculate the surface density maps using grid of pixel size

        partsDATA = data.copy()
        map_mass = project_gas(
            partsDATA, mode, number_of_pixels, extent, face_on_rotation_matrix
        )
        map_metals = integrate_metallicity_using_grid(
            partsDATA, number_of_pixels, extent, face_on_rotation_matrix
        )

        star_formation_rate_mask = partsDATA[:, 10] > 0.0
        partsDATA = partsDATA[star_formation_rate_mask, :]
        map_SFR = project_gas(
            partsDATA, 2, number_of_pixels, extent, face_on_rotation_matrix
        )

    else:
        partsDATA = data.copy()
        map_mass = project_gas_with_azimuthal_average(
            partsDATA, mode, face_on_rotation_matrix, size
        )
        map_metals = project_metals_with_azimuthal_average(
            partsDATA, face_on_rotation_matrix, size
        )

        star_formation_rate_mask = np.where(partsDATA[:, 10] > 0.0)[0]
        partsDATA = partsDATA[star_formation_rate_mask, :]

        if len(star_formation_rate_mask) > 0:
            map_SFR = project_gas_with_azimuthal_average(
                partsDATA, 2, face_on_rotation_matrix, size
            )
        else:
            map_SFR = np.zeros(len(map_mass))

    # Bounds
    map_SFR[map_SFR <= 0] = 3e-7
    map_mass[map_mass <= 0] = 3e-7

    surface_density = np.log10(map_mass.flatten())  # Msun / kpc^2
    surface_density -= 6  # Msun / pc^2
    SFR_surface_density = np.log10(map_SFR.flatten())  # Msun / yr / kpc^2
    tgas = surface_density - SFR_surface_density + 6.0

    return surface_density, SFR_surface_density, tgas, map_metals


def median_relations(x, y):

    xrange = np.arange(-1, 3, 0.1)

    xvalues = np.ones(len(xrange) - 1) * (-10)
    yvalues = np.zeros(len(xrange) - 1)
    yvalues_err_down = np.zeros(len(xrange) - 1)
    yvalues_err_up = np.zeros(len(xrange) - 1)

    perc = [16, 84]

    for i in range(0, len(xrange) - 2):
        mask = (x > xrange[i]) & (x < xrange[i + 1])
        if len(x[mask]) > 4:
            xvalues[i] = np.median(x[mask])
            yvalues[i] = np.median(y[mask])
            yvalues_err_down[i], yvalues_err_up[i] = np.transpose(
                np.percentile(y[mask], perc)
            )

    mask = xvalues > -10
    xvalues = xvalues[mask]
    yvalues = yvalues[mask]
    yvalues_err_down = yvalues_err_down[mask]
    yvalues_err_up = yvalues_err_up[mask]

    return xvalues, yvalues, yvalues_err_down, yvalues_err_up


def make_KS_data(
    particles_data,
    ang_momentum,
    mode,
    index,
    output_path,
    simulation_name,
    combined_data,
):

    # Plotting KS relations with size
    method = "grid"
    size = 0.75  # kpc

    # Get the surface densities
    surface_density, SFR_surface_density, tgas, metals = KS_relation(
        particles_data, ang_momentum, mode, method, size
    )

    if mode == 0:

        combined_data.molecular_gas_surface_density = np.append(
            combined_data.molecular_gas_surface_density, surface_density
        )
        combined_data.depletion_time_molecular_gas = np.append(
            combined_data.depletion_time_molecular_gas, tgas
        )

        np.savetxt(
            f"{output_path}/KS_molecular_relation_grid_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, SFR_surface_density]),
        )

        np.savetxt(
            f"{output_path}/molecular_gas_depletion_timescale_grid_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, tgas]),
        )

    elif mode == 1:

        combined_data.neutral_gas_surface_density = np.append(
            combined_data.neutral_gas_surface_density, surface_density
        )
        combined_data.depletion_time_neutral_gas = np.append(
            combined_data.depletion_time_neutral_gas, tgas
        )
        combined_data.SFR_surface_density = np.append(
            combined_data.SFR_surface_density, SFR_surface_density
        )
        combined_data.gas_metallicity = np.append(combined_data.gas_metallicity, metals)

        np.savetxt(
            f"{output_path}/KS_relation_best_grid_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, SFR_surface_density]),
        )

        np.savetxt(
            f"{output_path}/gas_depletion_timescale_best_grid_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, tgas]),
        )

    ###### Making KS plots with azimuthally averaged method #################

    # Plotting KS relations with size
    method = "radii"
    size = 0.75  # kpc

    surface_density, SFR_surface_density, tgas, metals = KS_relation(
        particles_data, ang_momentum, mode, method, size
    )

    if mode == 0:
        np.savetxt(
            f"{output_path}/KS_molecular_relation_radii_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, SFR_surface_density]),
        )

        np.savetxt(
            f"{output_path}/molecular_gas_depletion_timescale_radii_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, tgas]),
        )

    elif mode == 1:
        np.savetxt(
            f"{output_path}/KS_relation_best_radii_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, SFR_surface_density]),
        )

        np.savetxt(
            f"{output_path}/gas_depletion_timescale_best_radii_%i_" % (index)
            + simulation_name
            + ".txt",
            np.transpose([surface_density, tgas]),
        )

    return


def surface_ratios(data, ang_momentum, method, size):

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    if method == "grid":
        image_diameter = 60
        extent = [-30, 30]  # kpc
        number_of_pixels = int(image_diameter / size + 1)

        # Calculate the maps using grid
        map_H2 = project_gas(data, 0, number_of_pixels, extent, face_on_rotation_matrix)
        map_HI = project_gas(data, 3, number_of_pixels, extent, face_on_rotation_matrix)

    else:
        # Calculate the maps using azimuthally-average shells
        map_H2 = project_gas_with_azimuthal_average(
            data, 0, face_on_rotation_matrix, size
        )
        map_HI = project_gas_with_azimuthal_average(
            data, 3, face_on_rotation_matrix, size
        )

    map_gas = map_H2 + map_HI

    # Bounds
    map_H2[map_H2 <= 0] = 3e-7
    map_gas[map_gas <= 0] = 3e-7
    H2_to_neutral_ratio = map_H2 / map_gas

    neutral_gas_surface_density = np.log10(map_gas.flatten())  # HI+H2 Msun / kpc^2
    molecular_gas_surface_density = np.log10(map_H2.flatten())  # H2 Msun / kpc^2

    neutral_gas_surface_density -= 6  # HI+H2 Msun / pc^2
    molecular_gas_surface_density -= 6  # H2 Msun / pc^2

    H2_to_neutral_ratio_density = np.log10(H2_to_neutral_ratio.flatten())  # no units

    return (
        neutral_gas_surface_density,
        molecular_gas_surface_density,
        H2_to_neutral_ratio_density,
    )


def Krumholz_eq39(Sigma_neutral, f):
    Zprime = 1.0
    psi = 1.6  # fiducial from Krumholz
    s = 1.0 / f * Sigma_neutral * Zprime / psi
    RH2 = (
        np.power(
            (1.0 + np.power(s / 11.0, 3) * np.power((125.0 + s) / (96.0 + s), 3)),
            1.0 / 3.0,
        )
        - 1.0
    )
    return RH2


def make_surface_density_ratios(
    data, ang_momentum, index, output_path, simulation_name, combined_data
):

    # Get the surface densities
    method = "grid"
    binsize = 0.75  # kpc
    Sigma_gas_neutral, Sigma_gas_molecular, Sigma_ratio = surface_ratios(
        data, ang_momentum, method, binsize
    )

    combined_data.H2_to_neutral_surface_density_ratio = np.append(
        combined_data.H2_to_neutral_surface_density_ratio, Sigma_ratio
    )

    np.savetxt(
        f"{output_path}/Surface_density_ratio_grid_%i_" % (index)
        + simulation_name
        + ".txt",
        np.transpose([Sigma_gas_neutral, Sigma_ratio]),
    )

    ########################################################################
    # Get the surface densities
    method = "radii"

    binsize = 0.25  # kpc
    (
        Sigma_gas_250pc_neutral,
        Sigma_gas_250pc_molecular,
        Sigma_ratio_250pc,
    ) = surface_ratios(data, ang_momentum, method, binsize)

    binsize = 0.75  # kpc
    (
        Sigma_gas_800pc_neutral,
        Sigma_gas_800pc_molecular,
        Sigma_ratio_800pc,
    ) = surface_ratios(data, ang_momentum, method, binsize)

    combined_data.radii_neutral_gas_surface_density = np.append(
        combined_data.radii_neutral_gas_surface_density, Sigma_gas_800pc_neutral
    )
    combined_data.radii_H2_to_neutral_surface_density_ratio = np.append(
        combined_data.radii_H2_to_neutral_surface_density_ratio, Sigma_ratio_800pc
    )

    np.savetxt(
        f"{output_path}/Surface_density_ratio_radii_250pc_%i_" % (index)
        + simulation_name
        + ".txt",
        np.transpose([Sigma_gas_250pc_neutral, Sigma_ratio_250pc]),
    )
    np.savetxt(
        f"{output_path}/Surface_density_ratio_radii_800pc_%i_" % (index)
        + simulation_name
        + ".txt",
        np.transpose([Sigma_gas_800pc_neutral, Sigma_ratio_800pc]),
    )


def calculate_integrated_quantities(data, ang_momentum, radius, mode):

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)

    x, y, _ = np.matmul(face_on_rotation_matrix, data[:, :3].T)
    r = np.sqrt(x ** 2 + y ** 2)
    select = r <= radius

    surface = np.pi * radius ** 2
    if mode == 0:
        m = data[select, 9]
    if mode == 1:
        m = data[select, 9] + data[select, 8]

    # If we have gas within rhalfMs
    if len(m) > 0:
        Sigma_gas = np.log10(np.sum(m) / surface) - 6.0  # Msun / pc^2

        sfr = data[select, 10]
        sfr = sfr[sfr > 0]
        Sigma_SFR = np.log10(np.sum(sfr) / surface)  # Msun / yr / kpc^2

    else:
        Sigma_gas = -6.5
        Sigma_SFR = -6.5

    return Sigma_gas, Sigma_SFR


def make_KS_plots(
    data, ang_momentum, index, output_path, simulation_name, combined_data
):

    for mode, project in enumerate(
        ["molecular_hydrogen_masses", "not_ionized_hydrogen_masses"]
    ):
        make_KS_data(
            data, ang_momentum, mode, index, output_path, simulation_name, combined_data
        )

    make_surface_density_ratios(
        data, ang_momentum, index, output_path, simulation_name, combined_data
    )


def calculate_surface_densities(data, ang_momentum, galaxy_data, index):

    # If we have gas, calculate ..
    radius = galaxy_data.half_mass_radius_star[index]

    # Mode ==0 : "molecular_hydrogen_masses"
    # Mode ==1 : "not_ionized_hydrogen_masses"
    Sigma_H2, Sigma_SFR_H2 = calculate_integrated_quantities(
        data, ang_momentum, radius, 0
    )
    Sigma_gas, Sigma_SFR = calculate_integrated_quantities(
        data, ang_momentum, radius, 1
    )
    Sigma = np.array([Sigma_H2, Sigma_gas, Sigma_SFR])
    galaxy_data.add_surface_density(Sigma, index)
