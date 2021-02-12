"""
Description here
"""
import sys
import os
import h5py
import numpy as np
from velociraptor import load
from catalogue import Galaxy_data
from particles import calculate_morphology, make_particle_data, make_gas_particle_data
from plotter.plot import plot_galaxy, plot_morphology, plot_galaxy_sparts, plot_galaxy_gas_parts
import unyt

class Sim:
    def __init__(self,folder,snap):
        self.snapshot = os.path.join(folder,"colibre_0%03i.hdf5"%snap)
        self.subhalo_properties = os.path.join(folder,"halo_0%03i.properties.0"%snap)
        snapshot_file = h5py.File(self.snapshot, "r")
        self.boxSize = snapshot_file["/Header"].attrs["BoxSize"][0] * 1e3 #kpc
        self.a = snapshot_file["/Header"].attrs["Scale-factor"]


if __name__ == '__main__':

    file = '/snap7/scratch/dp004/dc-chai1/my_cosmological_box/XMAS2020_L006N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY'
    snapshot = 623
    #file = sys.argv[1]
    #snapshot = int(sys.argv[2])

    siminfo = Sim(file, snapshot)
    properties = load(siminfo.subhalo_properties)
    stellar_mass = properties.masses.m_star_30kpc
    stellar_mass.convert_to_units("msun")
    
    # Selecting galaxies more massive than lower limit
    lower_mass = 1e6 * unyt.msun #! Option of lower limit
    halo_catalogue = np.where(stellar_mass >= lower_mass)[0]
    
    # Selecting centrals only
    structure_type = properties.structure_type.structuretype.value
    centrals = np.where(structure_type[halo_catalogue] == 10)[0]
    halo_catalogue = halo_catalogue[centrals]

    # Sample :
    num_halos = len(halo_catalogue)
    stellar_mass = np.log10(stellar_mass[halo_catalogue])
    galaxy_data = Galaxy_data(stellar_mass,num_halos)

    # Create data for stars and gas
    stars_data = make_particle_data(siminfo,4)
    gas_data = make_gas_particle_data(siminfo)

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
        morphology, particle_data = calculate_morphology(subhalo_data, stars_data, 4, siminfo)

        # Make galaxy plot perhaps.. only first 10.
        if i < 10: 
            plot_galaxy(particle_data,morphology[0],i,4)
            plot_galaxy_sparts(particle_data,morphology[0],i)
            
        # Calculate morphology estimators: kappa, axial ratios for HI+H2 gas
        gas_morphology, particle_data = calculate_morphology(subhalo_data, gas_data, 0, siminfo)

        if i < 10: # Only 10 most massive..
            plot_galaxy_gas_parts(particle_data,gas_morphology[0],i)
            plot_galaxy(particle_data,gas_morphology[0],i,0)

        # Add gas data
        morphology = np.append(morphology, gas_morphology)
        
        # Store info in galaxy class and continue
        galaxy_data.add_morphology(morphology,i)

    # Finish plotting and output hdf5 file
    plot_morphology(galaxy_data)
    #output_to_hdf5(galaxy_data)






