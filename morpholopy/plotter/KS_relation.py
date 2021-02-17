import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from swiftsimio.visualisation.projection import project_gas
from unyt import unyt_quantity


def KS_relation(data, project, size=1., spatial_size=50.):
    """
    Parameters
    ----------

    data: SWIFTData
        SWIFTsimIO Data set

    project: str
        Quantity to project for the gas mass.

    size: float
        Size of the pixels in kpc

    spatial_size:
        Region to visualise around the galaxy, in kpc (total diameter).
    """

    image_diameter = unyt_quantity(spatial_size, units="kpc")
    bs = 0.5 * image_diameter
    size = unyt_quantity(size, units="kpc")
    extent = [0.5 * data.metadata.boxsize[0] - bs, 0.5 * data.metadata.boxsize[0] + bs] * 2

    number_of_pixels = int(image_diameter / size + 1)

    # Calculate the SPH smoothed maps

    common_parameters = dict(
        data=data,
        resolution=number_of_pixels,
        region=extent,
        parallel=True,
        backend="subsampled_extreme",
    )

    map_mass = project_gas(
        project=project,
        **common_parameters
    )

    star_formation_rate_mask = data.gas.star_formation_rates > 0.0

    map_SFR = project_gas(
        project="star_formation_rates",
        mask=star_formation_rate_mask,
        **common_parameters
    )
    # Bounds for SFR
    map_SFR[map_SFR <= 0] = 1e-7

    surface_density = np.log10(map_mass.to("Solar_Mass / pc**2").flatten())
    SFR_surface_density = np.log10(map_SFR.to("Solar_Mass / kpc**2 / year").flatten())

    surface_range = np.arange(-1, 3, .25)

    SFR_values = np.zeros(len(surface_range) - 1)
    SFR_values_err_down = np.zeros(len(surface_range) - 1)
    SFR_values_err_up = np.zeros(len(surface_range) - 1)

    perc = [16, 86]

    for i in range(0, len(surface_range) - 1):
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


def make_KS_plots(data, project, mode, igal, output_path, size=0.65):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(1, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Get the surface densities
    surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up = KS_relation(
        data, project, size=size)

    plt.plot(surface_density, SFR_surface_density, label="Our run")
    plt.fill_between(surface_density, SFR_surface_density - SFR_surface_density_err_down,
                     SFR_surface_density + SFR_surface_density_err_up, alpha=0.2)
    plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), color="red", label="KS law (Kennicutt 98)", linestyle="--")
    plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
    plt.xlim(0, 3.0)
    plt.ylim(-6.0, 0.0)

    if mode == 0:
        # fix the x-axis and same image
        plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend()
        plt.savefig(f"{output_path}/KS_molecular_relation_%i.png" % (igal))
        #np.savetxt(f"{output_path}/KS_molecular_relation_file_{snapshot_number:04d}.txt", np.transpose(
        #    [surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up, tgas,
        #     tgas_err_down, tgas_err_up]))

    elif mode == 1:
        # fix the x-axis and same image
        plt.xlabel("log $\\Sigma_{HI}+ \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
        plt.legend()
        plt.savefig(f"{output_path}/KS_relation_best_%i.png" % (igal))
        #np.savetxt(f"{output_path}/KS_relation_best_file_{snapshot_number:04d}.txt", np.transpose(
        #    [surface_density, SFR_surface_density, SFR_surface_density_err_down, SFR_surface_density_err_up, tgas,
        #     tgas_err_down, tgas_err_up]))
    plt.close()



def KS_plots(data,igalaxy,output_path):

    for mode, project in enumerate(["molecular_hydrogen_masses", "not_ionized_hydrogen_masses"]):
            make_KS_plots(data, project, mode, igalaxy, output_path, size=0.65)
