import numpy as np
import h5py
from functions import calculate_kappa_co, AxialRatios
from velociraptor.swift.swift import to_swiftsimio_dataset


def H2fraction(Temp,Den,Metal,z):
    
    filename = 'UV_dust1_CR1_G1_shield1.hdf5'
    
    with h5py.File(filename, "r") as f:
        RedshiftBins       = f['TableBins/RedshiftBins'][:]
        MetallicityBins    = f['TableBins/MetallicityBins'][:]
        TemperatureBins    = f['TableBins/TemperatureBins'][:]
        DensityBins        = f['TableBins/DensityBins'][:]
        HydrogenFractionsVol = f['Tdep/HydrogenFractionsVol'][:]
    
    H2frac = np.zeros(len(Temp))
    HIfrac = np.zeros(len(Temp))

    Ri = np.where(RedshiftBins>=z)[0]

    for j in range(0,len(Temp)):
        Ti = np.where(TemperatureBins>=np.log10(Temp[j]))[0]
        Di = np.where(DensityBins>=np.log10(Den[j]))[0]
        
        if len(Di)==0:Di=[len(DensityBins)]
        if Metal[j]==0.0:Metal[j]=1e-16

        Zi = np.where(MetallicityBins>=np.log10(Metal[j]))[0]
        
        if len(Ti)==0 or len(Di)==0 or len(Zi)==0:print(np.log10(Den[j]))
        if len(Ti)==0 or len(Di)==0 or len(Zi)==0:continue
        
        H2frac[j] = HydrogenFractionsVol[Ri[0]-1,Ti[0]-1,Zi[0]-1,Di[0]-1,2]
        HIfrac[j] = HydrogenFractionsVol[Ri[0]-1,Ti[0]-1,Zi[0]-1,Di[0]-1,0]
    
    return 10**H2frac, 10**HIfrac


def make_gas_particle_data(siminfo):
    
    # Create particle data :
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz]
    snapshot_file = h5py.File(siminfo.snapshot, "r")
    num_parttype = snapshot_file["/Header"].attrs["NumPart_Total"][0]
    snapshot_data = np.zeros((num_parttype,7))
    snapshot_data[:,0:3] = snapshot_file['PartType0/Coordinates'][:][:] * siminfo.a * 1e3 #kpc
    snapshot_data[:,4:7] = snapshot_file['PartType0/Velocities'][:][:] #km/s

    # Read units
    unit_length_in_cgs = snapshot_file["/Units"].attrs["Unit length in cgs (U_L)"]
    unit_mass_in_cgs = snapshot_file["/Units"].attrs["Unit mass in cgs (U_M)"]
    unit_density_in_cgs = unit_mass_in_cgs / unit_length_in_cgs**3

    masses = snapshot_file['PartType0/Masses'][:] * 1e10       # Msun
    density = snapshot_file['PartType0/Densities'][:] * unit_density_in_cgs  # g/cm^3
    density /= siminfo.a**3 # to physical units
    temperature = snapshot_file['PartType0/Temperatures'][:]   # K
    Z = snapshot_file['PartType0/MetalMassFractions'][:]
    XH = snapshot_file['PartType0/ElementMassFractions'][:,0]  # Hydrogen Mass fraction

    redshift = 1./siminfo.a - 1.
    # Calculate H2 and HI masses
    mu = 1.6726219e-24 # gr
    H2frac, HIfrac = H2fraction(temperature,XH * density/mu, Z, redshift)
    mH2 = H2frac * XH * masses
    mHI = HIfrac * XH * masses
    snapshot_data[:,3] = mH2 + mHI
    
    return snapshot_data

def make_particle_data(siminfo,groups,halo_id,parttype):
    #def make_particle_data(siminfo, parttype):
    # Create particle data :
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz
    # | (7) hsml ]

    particles, unbound_particles = groups.extract_halo(halo_id=int(halo_id))
    data, mask = to_swiftsimio_dataset(particles,
        siminfo.snapshot,
        generate_extra_mask=True
    )

    if parttype == 0:
        gas_mass = data.gas.masses[mask.gas].value * 1e10 #Msun
        gas_n_parts = len(gas_mass)
        gas_data = np.zeros((gas_n_parts,9))
        gas_data[:,0:3] = data.gas.coordinates[mask.gas].value * siminfo.a * 1e3 #kpc
        gas_data[:, 8] = gas_mass # Msun
        gas_data[:,4:7] = data.gas.velocities[mask.gas].value #km/s
        gas_data[:,7] = data.gas.smoothing_lengths[mask.gas].value * siminfo.a * 1e3 #kpc

        XH = data.gas.element_mass_fractions.hydrogen[mask.gas].value
        gas_HI = data.gas.species_fractions.HI[mask.gas].value
        gas_H2 = data.gas.species_fractions.H2[mask.gas].value * 2.
        gas_HI_H2 = gas_HI * XH * gas_mass + gas_H2 * XH * gas_mass
        gas_data[:, 3] = gas_HI_H2 # Msun
        return gas_data

    if parttype == 4:
        stars_mass = data.stars.masses[mask.stars].value * 1e10
        stars_n_parts = len(stars_mass)
        stars_data = np.zeros((stars_n_parts,8))
        stars_data[:,0:3] = data.stars.coordinates[mask.stars].value * siminfo.a * 1e3 #kpc
        stars_data[:,3] = stars_mass #Msun
        stars_data[:,4:7] = data.stars.velocities[mask.stars].value #km/s
        stars_data[:,7] = data.stars.smoothing_lengths[mask.stars].value * siminfo.a * 1e3 #kpc
        return stars_data


def calculate_morphology(subhalo_data, part_data, siminfo):
    # Subhalo data :
    # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz]

    # Calculate kappa and specific angular momentum
    kappa, momentum, part_data = calculate_kappa_co(subhalo_data,part_data,siminfo)
    
    # Calculate axis ratios
    axis_1, axis_2, axis_3 = AxialRatios(part_data[:,:3], part_data[:,3])
    morphology = np.array([kappa, momentum, axis_1, axis_2, axis_3])

    return morphology, part_data
