import numpy as np
import h5py
from functions import calculate_kappa_co, AxialRatios
#from velociraptor.swift.swift import to_swiftsimio_dataset
from swiftsimio import load

def make_masks(siminfo,halo):
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")
    snapshot_file = h5py.File(siminfo.snapshot, "r")

    star_ids = snapshot_file["/PartType4/ParticleIDs"][:]
    gas_ids = snapshot_file["/PartType0/ParticleIDs"][:]

    halo_start_position = group_file["Offset"][halo]
    halo_end_position = group_file["Offset"][halo + 1]

    particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]

    _, _, mask_stars = np.intersect1d(particle_ids_in_halo, star_ids, assume_unique=True,
                                             return_indices=True, )

    _, _, mask_gas = np.intersect1d(particle_ids_in_halo, gas_ids, assume_unique=True,
                                             return_indices=True, )

    return mask_gas, mask_stars

def make_particle_data(siminfo,halo_id):
    # Create particle data :
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz
    # | (7) hsml ]

    #particles, unbound_particles = groups.extract_halo(halo_id=int(halo_id))
    #data, mask = to_swiftsimio_dataset(particles,
    #    siminfo.snapshot,
    #    generate_extra_mask=True
    #)

    mask_gas, mask_stars = make_masks(siminfo,halo_id)
    data = load(siminfo.snapshot)

    gas_mass = data.gas.masses[mask_gas].value * 1e10 #Msun
    gas_n_parts = len(gas_mass)
    gas_data = np.zeros((gas_n_parts,9))
    gas_data[:,0:3] = data.gas.coordinates[mask_gas].value * siminfo.a * 1e3 #kpc
    gas_data[:, 8] = gas_mass # Msun
    gas_data[:,4:7] = data.gas.velocities[mask_gas].value #km/s
    gas_data[:,7] = data.gas.smoothing_lengths[mask_gas].value * siminfo.a * 1e3 #kpc
    
    XH = data.gas.element_mass_fractions.hydrogen[mask_gas].value
    gas_HI = data.gas.species_fractions.HI[mask_gas].value
    gas_H2 = data.gas.species_fractions.H2[mask_gas].value * 2.
    gas_HI_H2 = gas_HI * XH * gas_mass + gas_H2 * XH * gas_mass
    gas_data[:, 3] = gas_HI_H2 # Msun

    stars_mass = data.stars.masses[mask_stars].value * 1e10
    stars_n_parts = len(stars_mass)
    stars_data = np.zeros((stars_n_parts,8))
    stars_data[:,0:3] = data.stars.coordinates[mask_stars].value * siminfo.a * 1e3 #kpc
    stars_data[:,3] = stars_mass #Msun
    stars_data[:,4:7] = data.stars.velocities[mask_stars].value #km/s
    stars_data[:,7] = data.stars.smoothing_lengths[mask_stars].value * siminfo.a * 1e3 #kpc
    
    return gas_data, stars_data


def calculate_morphology(subhalo_data, part_data, siminfo):
    # Subhalo data :
    # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz]

    # Calculate kappa and specific angular momentum
    kappa, momentum, part_data = calculate_kappa_co(subhalo_data,part_data,siminfo)
    
    # Calculate axis ratios
    axis_1, axis_2, axis_3 = AxialRatios(part_data[:,:3], part_data[:,3])
    morphology = np.array([kappa, momentum, axis_1, axis_2, axis_3])

    return morphology, part_data
