"""
Description here
"""

from ArgumentParser import ArgumentParser
import glob

from plotter.plot_galaxy import visualize_galaxy
from object.unitilies.luminosities import MakeGrid
from plotter.KS_relation import make_KS_plots, calculate_surface_densities
from plotter.KS_comparison import make_comparison_plots
from plotter.plot_morphology import output_morphology, plot_morphology
from plotter.plot_surface_densities import plot_surface_densities
from object import simulation_metadata
from plotter.loadplots import loadGalaxyPlots
from plotter import html
from tqdm import tqdm
import unyt


def compute_morpholopy(
    siminfo: simulation_metadata.SimInfo, num_galaxies: int, output_path: str
):

    # Loading photometry grids for interpolation
    system = "GAMA"  # hard-coded for now
    pgrids = {}
    for pht in glob.glob(f"./photometry/{system}/*"):
        pgrids[pht[-1]] = MakeGrid(pht)

    # Loop over the sample to calculate morphological parameters
    for count, i in enumerate(tqdm(siminfo.halo_data.halo_ids)):

        print(count, i)

        # Read particle data
        gas_data, stars_data = siminfo.make_particle_data(halo_id=i)

        if len(gas_data) == 0:
            continue

        # Calculate morphology estimators: kappa, axial ratios for stars ..
        stars_ang_momentum, stars_data = siminfo.calculate_morphology(
            stars_data, count, 4
        )

        # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas ..
        gas_ang_momentum, gas_data = siminfo.calculate_morphology(gas_data, count, 0)

        # Calculate stellar luminosities
        star_abmags = siminfo.calculate_luminosities(stars_data, pgrids)

        # Calculate surface densities for HI+H2 gas ..
        calculate_surface_densities(
            gas_data, gas_ang_momentum, siminfo.halo_data, count
        )

        # Make plots for individual galaxies, perhaps.. only first 10
        if count < num_galaxies:
            visualize_galaxy(
                stars_data,
                gas_data,
                star_abmags,
                stars_ang_momentum,
                gas_ang_momentum,
                siminfo.halo_data,
                count,
                output_path,
                siminfo.simulation_name,
            )

            make_KS_plots(
                gas_data,
                stars_ang_momentum,
                siminfo.halo_data,
                count,
                output_path,
                siminfo.simulation_name,
            )

    # Finish plotting and output webpage
    output_morphology(siminfo.halo_data, output_path, siminfo.simulation_name)
    plot_surface_densities(siminfo.halo_data, output_path, siminfo.simulation_name)
    siminfo.output_galaxy_data(output_path=output_path)


def main(config: ArgumentParser):

    min_stellar_mass = unyt.unyt_quantity(1e8, "Msun")
    web = None

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]

        sim_info = simulation_metadata.SimInfo(
            directory=directory,
            snapshot=snapshot,
            catalogue=catalogue,
            name=sim_name,
            galaxy_min_stellar_mass=min_stellar_mass,
        )

        # Make initial website

        if sim == 0:
            web = html.make_web(sim_info.snapshot)
        elif web is not None:
            html.add_metadata_to_web(web, sim_info.snapshot)

        # Run morpholoPy
        compute_morpholopy(
            sim_info,
            config.number_of_galaxies,
            output_path=config_parameters.output_directory,
        )

    make_comparison_plots(
        output_path=config_parameters.output_directory,
        name_list=config_parameters.name_list,
        num_of_galaxies=config_parameters.number_of_galaxies,
    )
    plot_morphology(
        output_path=config_parameters.output_directory,
        name_list=config_parameters.name_list,
    )

    # After making individual plots finish up the website
    # Load galaxy plots
    loadGalaxyPlots(
        web,
        config_parameters.output_directory,
        config_parameters.number_of_galaxies,
        config_parameters.name_list,
    )

    # Finish and output html file
    html.render_web(web, config.output_directory)


if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)
