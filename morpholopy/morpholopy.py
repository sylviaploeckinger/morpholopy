"""
Description here
"""

from argumentparser import ArgumentParser

from plotter.plot_galaxy import visualize_galaxy
from plotter.KS_relation import make_KS_plots, calculate_surface_densities
from plotter.KS_comparison import make_comparison_plots
from plotter.plot_morphology import write_morphology_data_to_file, plot_morphology
from plotter.plot_surface_densities import plot_surface_densities
from object import simulation_data
from plotter.loadplots import loadGalaxyPlots
from plotter import html
from time import time
from tqdm import tqdm

from plotter.surface_maps_face_edge import surface_densities_overview
from plotter.species_transitions import species_transitions_combined

def compute_galaxy_morpholopy(
    sim_info: simulation_data.SimInfo,
    halo_counter: int,
    num_galaxies: int,
    output_path: str,
) -> None:
    """
    Computes morphological properties of galaxies from halo catalogue

    Parameters
    ----------
    sim_info: simulation_data.SimInfo
    Container with all simulation data

    halo_counter: int

    num_galaxies: int
    Number of galaxies to visualise

    output_path: str
    Path to the output directory
    """

    # Read particle data for a specific halo
    gas_data, stars_data = sim_info.make_particle_data(
        halo_id=sim_info.halo_data.halo_ids[halo_counter]
    )

    if len(gas_data) == 0:
        return

    # Calculate morphology estimators: kappa, axial ratios for stars ..
    stars_ang_momentum, stars_data = sim_info.calculate_morphology(
        stars_data, halo_counter, 4
    )

    # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas ..
    gas_ang_momentum, gas_data = sim_info.calculate_morphology(
        gas_data, halo_counter, 0
    )

    # Calculate stellar luminosities
    star_abmags = sim_info.calculate_luminosities(stars_data)

    # Calculate surface densities for HI+H2 gas ..
    calculate_surface_densities(
        gas_data, gas_ang_momentum, sim_info.halo_data, halo_counter
    )

    # Make plots for individual galaxies, perhaps.. only first 10
    if halo_counter < num_galaxies:
        '''
        visualize_galaxy(
            stars_data,
            gas_data,
            star_abmags,
            stars_ang_momentum,
            gas_ang_momentum,
            sim_info.halo_data,
            halo_counter,
            output_path,
            sim_info.simulation_name,
        )
        '''

        make_KS_plots(
            gas_data,
            stars_ang_momentum,
            halo_counter,
            output_path,
            sim_info.simulation_name,
            sim_info.combined_data,
        )

    return


def main(config: ArgumentParser):

    time_start = time()
    output_name_list = []
    output_number_of_galaxies_list = []
    web = None

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        # Fetch relevant input parameters from lists
        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]

        # Load all data and save it in SimInfo class
        sim_info = simulation_data.SimInfo(
            directory=directory,
            snapshot=snapshot,
            catalogue=catalogue,
            name=sim_name,
            galaxy_min_stellar_mass=config.min_stellar_mass,
        )

        output_name_list.append(sim_info.simulation_name)

        # Make initial part of the webpage
        if sim == 0:
            web = html.make_web(sim_info.snapshot)
        elif web is not None:
            html.add_metadata_to_web(web, sim_info.snapshot)

        # Load luminosity tables
        simulation_data.SimInfo.load_photometry_grid()

        # The actual number of galaxies to visualise for this run
        output_number_of_galaxies_list.append(
            min(sim_info.halo_data.number_of_haloes, config.number_of_galaxies)
        )

        print(
            f"Total number of haloes to analyse: {sim_info.halo_data.number_of_haloes}"
        )

        # Compute morphological properties (loop over haloes)
        print("Computing morphological properties...")

        surface_densities_overview(
                sim_name = sim_name,
                directory = directory,
                snapshot = snapshot,
                catalogue_file = catalogue,
                output_path=config.output_directory,
                nhalos = sim_info.halo_data.number_of_haloes,
                halo_min_stellar_mass = config.min_stellar_mass,
        )

        species_transitions_combined(
                sim_name = sim_name,
                output_path=config.output_directory,
                halo_min_stellar_mass = config.min_stellar_mass,
        )

        for i in tqdm(range(sim_info.halo_data.number_of_haloes)):
            compute_galaxy_morpholopy(
                sim_info=sim_info,
                num_galaxies=config.number_of_galaxies,
                output_path=config.output_directory,
                halo_counter=i,
            )

        write_morphology_data_to_file(
            sim_info.halo_data,
            sim_info.combined_data,
            config.output_directory,
            sim_info.simulation_name,
        )
        plot_surface_densities(
            sim_info.halo_data,
            sim_info.combined_data,
            config.output_directory,
            sim_info.simulation_name,
        )
        sim_info.write_galaxy_data_to_file(output_path=config.output_directory)

    num_galaxies_to_show = min(output_number_of_galaxies_list)

    make_comparison_plots(
        output_path=config.output_directory,
        name_list=output_name_list,
        num_of_galaxies_to_show=num_galaxies_to_show,
    )
    plot_morphology(
        output_path=config.output_directory,
        name_list=output_name_list,
    )

    # Load galaxy plots
    loadGalaxyPlots(
        web,
        config.output_directory,
        num_galaxies_to_show,
        output_name_list,
        config.min_stellar_mass,
    )

    # Finish and output html file
    html.render_web(web, config.output_directory)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

    return


if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)
