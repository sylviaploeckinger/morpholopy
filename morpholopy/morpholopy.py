"""
Description here
"""

import os
import h5py

from catalogue import HaloCatalogue
from particles import calculate_morphology, make_particle_data
from plotter.html import make_web, add_web_section, render_web, PlotsInPipeline
from plotter.plot import plot_morphology
from plotter.plot_galaxy import visualize_galaxy

from plotter.KS_relation import make_KS_plots, calculate_surface_densities
import unyt


class SimInfo:
    def __init__(self,folder,snap):
        self.snapshot = os.path.join(folder,"colibre_%04i.hdf5"%snap)
        self.subhalo_properties = os.path.join(folder,"halo_%04i.properties.0"%snap)
        self.catalog_groups = os.path.join(folder,"halo_%04i.catalog_groups.0"%snap)
        self.catalog_particles = os.path.join(folder, "halo_%04i.catalog_particles.0" % snap)
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
    lower_mass = 1e9 * unyt.msun  # ! Option of lower limit for gas mass
    halo_data = HaloCatalogue(siminfo,lower_mass)

    # Loop over the sample to calculate morphological parameters
    for i in range(halo_data.num):

        # Read particle data
        gas_data, stars_data = make_particle_data(siminfo, halo_data.halo_index[i])

        if len(gas_data) ==0: continue

        # Calculate morphology estimators: kappa, axial ratios for stars ..
        stars_ang_momentum, stars_data = calculate_morphology(halo_data, stars_data, siminfo, i, 4)
        print('stars momentum', stars_ang_momentum, i)

        # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas ..
        gas_ang_momentum, gas_data = calculate_morphology(halo_data, gas_data, siminfo, i, 0)
        print('gas momentum', gas_ang_momentum, i)

        # Calculate surface densities for HI+H2 gas ..
        calculate_surface_densities(gas_data, gas_ang_momentum, halo_data, i)

        # Make plots for individual galaxies, perhaps.. only first 10
        if i < 10:

            visualize_galaxy(stars_data, gas_data, stars_ang_momentum, gas_ang_momentum,
                             halo_data, i, GalPlotsInWeb, output_path)

            make_KS_plots(gas_data, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)

            title = '%i Galaxy ' % (i+1)
            id = abs(hash("galaxy and ks relation %i" %i))
            plots = GalPlotsInWeb.plots_details
            add_web_section(web,title,id,plots)
            GalPlotsInWeb.reset_plots_list()


    # Finish plotting and output webpage
    plot_morphology(halo_data, web, MorphologyPlotsInWeb, output_path )
    render_web(web, output_path)




