"""
Description here
"""

import os
import sys
import h5py
import numpy as np
from velociraptor import load as load_catalogue
#from velociraptor.particles import load_groups

from catalogue import HaloCatalogue
from particles import calculate_morphology, make_particle_data
from plotter.html import make_web, add_web_section, render_web, PlotsInPipeline
from plotter.plot import plot_morphology
from plotter.plot_galaxy import plot_galaxy, plot_galaxy_parts

from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from plotter.KS_relation import KS_plots, KS_relation, project_gas

import unyt

class SimInfo:
    def __init__(self,folder,snap):
        self.snapshot = os.path.join(folder,"colibre_0%03i.hdf5"%snap)
        self.subhalo_properties = os.path.join(folder,"halo_0%03i.properties.0"%snap)
        self.catalog_groups = os.path.join(folder,"halo_0%03i.catalog_groups.0"%snap)
        self.catalog_particles = os.path.join(folder, "halo_0%03i.catalog_particles.0" % snap)
        snapshot_file = h5py.File(self.snapshot, "r")
        self.boxSize = snapshot_file["/Header"].attrs["BoxSize"][0] * 1e3 #kpc
        self.a = snapshot_file["/Header"].attrs["Scale-factor"]


if __name__ == '__main__':
    from utils import *

    output_path = args.output
    siminfo = SimInfo(args.directory, args.number)

    # Loading simulation data in website table
    web = make_web(siminfo)
    PartPlotsInWeb = PlotsInPipeline()
    GalPlotsInWeb = PlotsInPipeline()
    MorphologyPlotsInWeb = PlotsInPipeline()
    KSPlotsInWeb = PlotsInPipeline()

    # Loading halo catalogue and selecting galaxies more massive than lower limit
    lower_mass = 1e9 * unyt.msun  # ! Option of lower limit
    halo_data = HaloCatalogue(siminfo,lower_mass)

    # Loop over the sample to calculate morphological parameters
    for i in range(0,halo_data.num):

        # Read particle data
        gas_data, stars_data = make_particle_data(siminfo, halo_data.halo_index[i])

        # Calculate morphology estimators: kappa, axial ratios for stars ..
        stars_ang_momentum, stars_data = calculate_morphology(halo_data, stars_data, siminfo, i, 4)

        # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas ..
        gas_ang_momentum, gas_data = calculate_morphology(halo_data, gas_data, siminfo, i, 0)

        # Make galaxy plot perhaps.. only first 10.
        if i < 10:
            #plot_galaxy_parts(stars_data, 4, stars_ang_momentum, halo_data, i, PartPlotsInWeb, output_path)
            #plot_galaxy_parts(gas_data, 0, gas_ang_momentum, halo_data, i, PartPlotsInWeb, output_path)

            plot_galaxy(stars_data, 4, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
            plot_galaxy(gas_data, 0, gas_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
            KS_plots(gas_data, gas_ang_momentum, i, KSPlotsInWeb, output_path)

        last = halo_data.num-1
        if halo_data.num > 10 : last = 9
        if i == last :
            #title = 'Visualizations (Particles)'
            #id = abs(hash("galaxy particles"))
            #plots = PartPlotsInWeb.plots_details
            #add_web_section(web,title,id,plots)

            title = 'Visualizations (SPH-viewer)'
            id = abs(hash("galaxy sph"))
            plots = GalPlotsInWeb.plots_details
            add_web_section(web,title,id,plots)

            title = 'KS relation'
            id = abs(hash("galaxy KS relation"))
            plots = KSPlotsInWeb.plots_details
            add_web_section(web,title,id,plots)

    # Finish plotting and output hdf5 file
    plot_morphology(halo_data, web, MorphologyPlotsInWeb, output_path )

    render_web(web, output_path)




