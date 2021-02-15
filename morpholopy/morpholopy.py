"""
Description here
"""

import os
import h5py
import numpy as np
from velociraptor import load
from velociraptor.particles import load_groups

from catalogue import Galaxy_data
from particles import calculate_morphology, make_particle_data
from plotter.html import make_web, add_web_section, render_web, PlotsInPipeline
from plotter.plot import plot_galaxy, plot_morphology, plot_galaxy_sparts, plot_galaxy_gas_parts
import unyt

class Sim:
    def __init__(self,folder,snap):
        self.snapshot = os.path.join(folder,"colibre_0%03i.hdf5"%snap)
        self.subhalo_properties = os.path.join(folder,"halo_0%03i.properties.0"%snap)
        self.catalog_groups = os.path.join(folder,"halo_0%03i.catalog_groups.0"%snap)
        snapshot_file = h5py.File(self.snapshot, "r")
        self.boxSize = snapshot_file["/Header"].attrs["BoxSize"][0] * 1e3 #kpc
        self.a = snapshot_file["/Header"].attrs["Scale-factor"]


if __name__ == '__main__':

    file = '/snap7/scratch/dp004/dc-chai1/my_cosmological_box/XMAS2020_L006N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY'
    snapshot = 623
    #file = sys.argv[1]
    #snapshot = int(sys.argv[2])

    #file = '/Users/Camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data'
    #snapshot = 34
    siminfo = Sim(file, snapshot)

    # Loading simulation data in website table
    web = make_web(siminfo)
    PartPlotsInWeb = PlotsInPipeline()
    GalPlotsInWeb = PlotsInPipeline()
    MorphologyPlotsInWeb = PlotsInPipeline()

    # Loading halo catalogue
    properties = load(siminfo.subhalo_properties)
    groups = load_groups(siminfo.catalog_groups, catalogue=properties)
    stellar_mass = properties.masses.m_star_30kpc
    stellar_mass.convert_to_units("msun")
    gas_mass = properties.masses.m_gas_30kpc
    gas_mass.convert_to_units("msun")

    # Selecting galaxies more massive than lower limit
    lower_mass = 1e9 * unyt.msun #! Option of lower limit
    halo_catalogue = np.where(stellar_mass >= lower_mass)[0]
    
    # Selecting centrals only
    structure_type = properties.structure_type.structuretype.value
    centrals = np.where(structure_type[halo_catalogue] == 10)[0]
    halo_catalogue = halo_catalogue[centrals]

    # Sample :
    num_halos = len(halo_catalogue)
    stellar_mass = np.log10(stellar_mass[halo_catalogue])
    gas_mass = np.log10(gas_mass[halo_catalogue])
    galaxy_data = Galaxy_data(stellar_mass,num_halos)


    # Loop over sample to calculate morphological parameters
    for halo, i in zip(halo_catalogue, range(0,num_halos)):
        
        # Create subhalo data :
        # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz]
        subhalo_data = np.zeros(6)
        subhalo_data[0] = float(properties.positions.xcminpot[halo]) * 1e3
        subhalo_data[1] = float(properties.positions.ycminpot[halo]) * 1e3
        subhalo_data[2] = float(properties.positions.zcminpot[halo]) * 1e3
        subhalo_data[3] = float(properties.velocities.vxcminpot[halo])
        subhalo_data[4] = float(properties.velocities.vycminpot[halo])
        subhalo_data[5] = float(properties.velocities.vzcminpot[halo])

        # Calculate morphology estimators: kappa, axial ratios for stars
        stars_data = make_particle_data(siminfo,groups,halo,4)
        morphology, stars_data = calculate_morphology(subhalo_data, stars_data, siminfo)

        # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas
        gas_data = make_particle_data(siminfo,groups,halo,0)
        gas_morphology, gas_data = calculate_morphology(subhalo_data, gas_data, siminfo)

        # Make galaxy plot perhaps.. only first 10.
        if i < 10:
            plot_galaxy_sparts(stars_data,morphology[0],stellar_mass[i], i,PartPlotsInWeb)
            plot_galaxy_gas_parts(gas_data,gas_morphology[0],gas_mass[i],i,PartPlotsInWeb)
            #plot_galaxy(stars_data,morphology[0],stellar_mass[i],i,4,GalPlotsInWeb)
            #plot_galaxy(gas_data,gas_morphology[0],gas_mass[i],i,0,GalPlotsInWeb)

        last = num_halos-1
        if num_halos > 10 : last = 9
        if i == last :
            title = 'Visualizations (Particles)'
            id = abs(hash("galaxy particles"))
            plots = PartPlotsInWeb.plots_details
            add_web_section(web,title,id,plots)

            #title = 'Visualizations (SPH-viewer)'
            #id = abs(hash("galaxy sph"))
            #plots = GalPlotsInWeb.plots_details
            #add_web_section(web,title,id,plots)

        # Store info in galaxy class and continue
        morphology = np.append(morphology, gas_morphology)
        galaxy_data.add_morphology(morphology,i)

    # Finish plotting and output hdf5 file
    plot_morphology(galaxy_data,web,MorphologyPlotsInWeb)

    render_web(web)




