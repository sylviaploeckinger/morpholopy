import numpy as np
import h5py
from functions import calculate_kappa_co, AxialRatios

def H2fraction(Temp,Den,Metal,z):
    
    filename = '/Users/Camila/swift-sim-tests/UV_dust1_CR1_G1_shield1.hdf5'
    
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

def make_particle_data(siminfo,parttype):
    
    # Create particle data :
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz]
    snapshot_file = h5py.File(siminfo.snapshot, "r")
    num_parttype = snapshot_file["/Header"].attrs["NumPart_Total"][parttype]
    snapshot_data = np.zeros((num_parttype,7))
    snapshot_data[:,0:3] = snapshot_file['PartType%i/Coordinates'%parttype][:][:] * siminfo.a * 1e3 #kpc
    snapshot_data[:,3] = snapshot_file['PartType%i/Masses'%parttype][:] * 1e10 #Msun
    snapshot_data[:,4:7] = snapshot_file['PartType%i/Velocities'%parttype][:][:] #km/s
    return snapshot_data

def select_parts(subhalo,siminfo,parttype):
    # Selecting particles within 30kpc aperture
    # subhalo contain subhalo data and is strutured as follow
    # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz  | (6)R200c[kpc]]
    snapshot_file = h5py.File(siminfo.snapshot, "r")
    num_parttype = snapshot_file["/Header"].attrs["NumPart_Total"][parttype]
    parts_data = snapshot_file['PartType%i/Coordinates'%parttype][:][:] * siminfo.a * 1e3 #kpc
    mask = np.arange(0,num_parttype,1)
    
    # Center positions & unwrap space
    parts_data[:,:3]-=subhalo[0:3].astype('float')-siminfo.boxSize/2   # centering onto subhalo CoP, and unwrap the box
    parts_data[:,:3]%=(siminfo.boxSize)
    parts_data[:,:3]-=siminfo.boxSize/2                                # end the unwrap
    
    # Compute distances
    distancesDATA = np.linalg.norm(parts_data[:,:3],axis=1)
    
    # Restrict particles
    mask = distancesDATA < 30
    return mask


def calculate_morphology(subhalo_data, part_data, parttype, siminfo):

    # Subhalo data :
    # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz]

    # Select particles from this particular halo
    mask = select_parts(subhalo_data, siminfo, parttype)

    # Create particle data
    part_data = part_data[mask,:]

    # Calculate kappa and specific angular momentum
    kappa, momentum, part_data = calculate_kappa_co(subhalo_data,part_data,siminfo)
    
    # Calculate axis ratios
    axis_1, axis_2, axis_3 = AxialRatios(part_data[:,:3], part_data[:,3])
    morphology = np.array([kappa, momentum, axis_1, axis_2, axis_3])

    return morphology, part_data
