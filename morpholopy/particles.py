import numpy as np
import h5py
from functions import calculate_kappa_co, AxialRatios
#from velociraptor.swift.swift import to_swiftsimio_dataset
from swiftsimio import load
from astropy.cosmology import WMAP9 as cosmo
import luminosities as lum

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
    z_snap = data.metadata.redshift
    
    gas_mass = data.gas.masses[mask_gas].value * 1e10 #Msun
    gas_n_parts = len(gas_mass)
    gas_data = np.zeros((gas_n_parts,13))
    gas_data[:,0:3] = data.gas.coordinates[mask_gas].value * siminfo.a * 1e3 #kpc
    gas_data[:, 3] = gas_mass # Msun
    gas_data[:,4:7] = data.gas.velocities[mask_gas].value #km/s
    gas_data[:,7] = data.gas.smoothing_lengths[mask_gas].value * siminfo.a * 1e3 #kpc
    
    XH = data.gas.element_mass_fractions.hydrogen[mask_gas].value
    gas_HI = data.gas.species_fractions.HI[mask_gas].value
    gas_H2 = data.gas.species_fractions.H2[mask_gas].value * 2.
    gas_data[:, 8] = gas_HI * XH * gas_mass # Msun
    gas_data[:, 9] = gas_H2 * XH * gas_mass # Msun
    gas_data[:, 10] = data.gas.star_formation_rates[mask_gas].value * 10227144.8879616 / 1e9 #Msun/yr
    gas_data[:, 11] = data.gas.densities[mask_gas].value * (1e10 / (siminfo.a * 1e3)**3) #Msun / kpc^3
    gas_data[:, 12] = data.gas.metal_mass_fractions[mask_gas].value/0.0134 # !assuming Solar metallicty

    stars_mass = data.stars.masses[mask_stars].value * 1e10
    stars_birthz = (pow(data.stars.birth_scale_factors,-1)-1.)[mask_stars].value
    stars_age = (cosmo.lookback_time(stars_birthz)-cosmo.lookback_time(z_snap)).value #in Gyr
    stars_Z = data.stars.metal_mass_fractions[mask_stars].value #absolute
    stars_initmass = data.stars.initial_masses[mask_stars].value * 1e10

    stars_n_parts = len(stars_mass)
    stars_data = np.zeros((stars_n_parts,12))
    stars_data[:,0:3] = data.stars.coordinates[mask_stars].value * siminfo.a * 1e3 #kpc
    stars_data[:,3] = stars_mass #Msun
    stars_data[:,4:7] = data.stars.velocities[mask_stars].value #km/s
    #stars_data[:,7] = data.stars.smoothing_lengths[mask_stars].value * siminfo.a * 1e3 #kpc
    stars_data[:, 7] = 0.5 * siminfo.baryon_maxsoft * np.ones(stars_mass.size)
    stars_data[:, 8] = stars_mass * (1.2348 / stars_data[:,7] )**3 #Msun/kpc^3
    stars_data[:,9] = stars_age
    stars_data[:,10] = stars_Z
    stars_data[:,11] = stars_initmass
    return gas_data, stars_data


def calculate_morphology(halo_data, part_data, siminfo, halo_index, parttype):

    # Calculate kappa and specific angular momentum
    kappa, specific_momentum, momentum, part_data = calculate_kappa_co(halo_data, part_data, siminfo, halo_index)

    # Calculate axis ratios
    axis_1, axis_2, axis_3 = AxialRatios(part_data[:,:3], part_data[:,3])
    morphology = np.array([kappa, specific_momentum, axis_1, axis_2, axis_3])

    #Store morphology parameters in halo data and continue
    if parttype == 4: halo_data.add_stellar_morphology(morphology, halo_index)
    if parttype == 0: halo_data.add_gas_morphology(morphology, halo_index)

    return momentum, part_data

def calculate_luminosities(halo_data, part_data, siminfo, halo_index, parttype, pgrids):
    star_abmags = {}
    
    for filt in pgrids.keys():
        grid, Z_p, t_p 	= pgrids[filt]
        fluxes = lum.BiPowInterp(t_p, Z_p, grid, part_data[:,9], part_data[:,10]) * part_data[:,11]
        star_abmags[filt]	= -2.5*np.log10(fluxes*lum.magfac)
    return star_abmags
    

