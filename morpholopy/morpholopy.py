"""
Description here
"""

import os
import h5py
import glob

from catalogue import HaloCatalogue, output_galaxy_data
from particles import calculate_morphology, make_particle_data, calculate_luminosities
from plotter.html import make_web, add_web_section, render_web, PlotsInPipeline, add_metadata_to_web
from plotter.plot_galaxy import visualize_galaxy
from luminosities import MakeGrid
from plotter.loadplots import loadGalaxyPlots
from plotter.KS_comparison import make_comparison_plots
from plotter.plot_morphology import plot_morphology, output_morphology
from plotter.KS_relation import make_KS_plots, calculate_surface_densities
from plotter.plot_surface_densities import plot_surface_densities
import unyt


class SimInfo:
    def __init__(self, folder, snap, output_path, name):
        self.name = name
        self.output_path = output_path
        self.snapshot = os.path.join(folder,"colibre_%04i.hdf5"%snap)

        properties = os.path.join(folder,"halo_%04i.properties.0"%snap)
        if os.path.exists(properties):
            self.subhalo_properties = os.path.join(folder,"halo_%04i.properties.0"%snap)
        else:
            self.subhalo_properties = os.path.join(folder,"halo_%04i.properties"%snap)

        catalog = os.path.join(folder,"halo_%04i.catalog_groups.0"%snap)
        if os.path.exists(catalog):
            self.catalog_groups = os.path.join(folder,"halo_%04i.catalog_groups.0"%snap)
        else :
            self.catalog_groups = os.path.join(folder,"halo_%04i.catalog_groups"%snap)

        catalog_particles = os.path.join(folder, "halo_%04i.catalog_particles.0" % snap)
        if os.path.exists(catalog_particles):
            self.catalog_particles = os.path.join(folder, "halo_%04i.catalog_particles.0" % snap)
        else :
            self.catalog_particles = os.path.join(folder, "halo_%04i.catalog_particles" % snap)

        snapshot_file = h5py.File(self.snapshot, "r")
        self.boxSize = snapshot_file["/Header"].attrs["BoxSize"][0] * 1e3 #kpc
        self.a = snapshot_file["/Header"].attrs["Scale-factor"]
        self.baryon_maxsoft = snapshot_file["/GravityScheme"].attrs['Maximal physical baryon softening length  [internal units]'] * 1e3 #kpc


def morpholopy(siminfo, web):

    # Loading photometry grids for interpolation
    system = 'GAMA'  # hard-coded for now
    pgrids = {}
    for pht in glob.glob(f'./photometry/{system}/*'):
        pgrids[pht[-1]] = MakeGrid(pht)

    # Loading halo catalogue and selecting galaxies more massive than lower limit
    lower_mass = 1e6 * unyt.msun  # ! Option of lower limit for gas mass
    halo_data = HaloCatalogue(siminfo, lower_mass)

    # Loop over the sample to calculate morphological parameters

    for i in range(halo_data.num):

        # Read particle data
        gas_data, stars_data = make_particle_data(siminfo, halo_data.halo_index[i])

        if len(gas_data) == 0: continue

        # Calculate morphology estimators: kappa, axial ratios for stars ..
        stars_ang_momentum, stars_data = calculate_morphology(halo_data, stars_data, siminfo, i, 4)

        # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas ..
        gas_ang_momentum, gas_data = calculate_morphology(halo_data, gas_data, siminfo, i, 0)

        # Calculate stellar luminosities
        star_abmags = calculate_luminosities(halo_data, stars_data, siminfo, i, 4, pgrids)

        # Calculate surface densities for HI+H2 gas ..
        calculate_surface_densities(gas_data, gas_ang_momentum, halo_data, i)

        # Make plots for individual galaxies, perhaps.. only first 10
        if i < 10:
            visualize_galaxy(stars_data, gas_data, star_abmags, stars_ang_momentum,
                             gas_ang_momentum, halo_data, i, siminfo)

            make_KS_plots(gas_data, stars_ang_momentum, halo_data, i, siminfo)

    # Finish plotting and output webpage
    output_morphology(halo_data, siminfo)
    plot_surface_densities(halo_data, siminfo)
    output_galaxy_data(halo_data,siminfo)

    return web


if __name__ == '__main__':
    from utils import *

    # Load MorpholoPy production details
    output_path = args.output
    number_of_inputs = len(args.snapshot)
    directory_list = args.directory
    snapshot_list = args.snapshot

    name_list = (
        args.run_names
        if args.run_names is not None
        else [None] * number_of_inputs
    )

    # Loop over simulation list
    for sims in range(number_of_inputs):
        directory = directory_list[sims]
        snap_number = int(snapshot_list[sims])
        sim_name = name_list[sims]
        siminfo = SimInfo(directory, snap_number,output_path, sim_name)

        # Make initial website
        if sims == 0: web = make_web(siminfo)
        if sims > 0: add_metadata_to_web(web, siminfo)

        # Run morpholoPy
        web = morpholopy(siminfo, web)

    make_comparison_plots(siminfo, name_list)
    plot_morphology(siminfo, name_list)

    # After making individual plots finish up the website
    # Load galaxy plots
    loadGalaxyPlots(web, name_list, siminfo.output_path)

    # Finish and output html file
    render_web(web, siminfo.output_path)



