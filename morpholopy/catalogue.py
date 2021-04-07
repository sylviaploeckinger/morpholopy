import numpy as np
import h5py
#from velociraptor import load as load_catalogue

class HaloCatalogue:
    def __init__(self,siminfo, lower_mass):

        properties = h5py.File(siminfo.subhalo_properties, "r")
        #properties = load_catalogue(siminfo.subhalo_properties)
        stellar_mass = properties["Aperture_mass_star_30_kpc"][:] * 1e10 #msun
        #stellar_mass = properties.masses.m_star_30kpc
        #stellar_mass.convert_to_units("msun")
        gas_mass = properties["Aperture_mass_gas_30_kpc"][:] * 1e10 #msun
        #gas_mass = properties.masses.m_gas_30kpc
        #gas_mass.convert_to_units("msun")
        half_mass_radius_star = properties["R_HalfMass_star"][:] * 1e3 #kpc
        #half_mass_radius_star = properties.radii.r_halfmass_star
        #half_mass_radius_star.convert_to_units("kpc")
        half_mass_radius_gas = properties["R_HalfMass_gas"][:] * 1e3 #kpc
        #half_mass_radius_gas = properties.radii.r_halfmass_gas
        #half_mass_radius_gas.convert_to_units("kpc")
        sfr = properties["Aperture_SFR_gas_30_kpc"][:] * 10227144.8879616 / 1e9 #Msun/yr
        #sfr = properties.apertures.sfr_gas_30_kpc
        halo_mass = properties["Mass_200crit"][:] * 1e10 #msun
        #halo_mass = properties.masses.mass_200crit
        #halo_mass.convert_to_units("msun")
        metallicity_gas_sfr = properties["Aperture_Zmet_gas_sf_30_kpc"][:]
        #metallicity_gas_sfr = properties.metallicity.zmet_gas_sf
        metallicity_gas = properties["Aperture_Zmet_gas_30_kpc"][:]
        #metallicity_gas = properties.metallicity.zmet_gas

        # Selecting galaxies more massive than lower limit
        catalogue = np.where(gas_mass >= lower_mass)[0]
        select = np.where(stellar_mass[catalogue] >= lower_mass)[0]
        catalogue = catalogue[select]

        # Selecting centrals only
        structure_type = properties["Structuretype"][:]
        #structure_type = properties.structure_type.structuretype.value
        centrals = np.where(structure_type[catalogue] == 10)[0]
        catalogue = catalogue[centrals]

        self.num = len(catalogue)
        self.halo_index = [ catalogue[i] for i in range(self.num) ]

        # Sample :
        self.stellar_mass = np.log10(stellar_mass[catalogue])
        self.gas_mass = np.log10(gas_mass[catalogue])
        self.halfmass_radius_star = half_mass_radius_star[catalogue]
        self.halfmass_radius_gas = half_mass_radius_gas[catalogue]
        self.star_formation_rate = sfr[catalogue]
        self.halo_mass = np.log10(halo_mass[catalogue])
        self.metallicity_gas_sfr = metallicity_gas_sfr[catalogue]
        self.metallicity_gas = metallicity_gas[catalogue]

        self.kappa_co = [ None for i in range(self.num) ]
        self.momentum = [ None for i in range(self.num) ]
        self.axis_ca = [ None for i in range(self.num) ]
        self.axis_cb = [ None for i in range(self.num) ]
        self.axis_ba = [ None for i in range(self.num) ]
        
        self.gas_kappa_co = [ None for i in range(self.num) ]
        self.gas_momentum = [ None for i in range(self.num) ]
        self.gas_axis_ca = [ None for i in range(self.num) ]
        self.gas_axis_cb = [ None for i in range(self.num) ]
        self.gas_axis_ba = [ None for i in range(self.num) ]

        self.sigma_H2 = [ None for i in range(self.num) ]
        self.sigma_gas = [ None for i in range(self.num) ]
        self.sigma_SFR = [ None for i in range(self.num) ]

        self.surface_density =  None
        self.SFR_density =  None
        self.ratio_densities = None
        self.metallicity = None

        self.radii_surface_density = None
        self.radii_surface_ratio = None
        self.radii_nbins = None


        # Subhalo data :
        self.xminpot = properties["Xcminpot"][catalogue] * 1e3 #kpc
        self.yminpot = properties["Ycminpot"][catalogue] * 1e3 #kpc
        self.zminpot = properties["Zcminpot"][catalogue] * 1e3 #kpc
        #self.xminpot = properties.positions.xcminpot[catalogue].value * 1e3
        #self.yminpot = properties.positions.ycminpot[catalogue].value * 1e3
        #self.zminpot = properties.positions.zcminpot[catalogue].value * 1e3

        self.vxminpot = properties["VXcminpot"][catalogue] #km/s
        self.vyminpot = properties["VYcminpot"][catalogue] #km/s
        self.vzminpot = properties["VZcminpot"][catalogue] #km/s
        #self.vxminpot = properties.velocities.vxcminpot[catalogue].value
        #self.vyminpot = properties.velocities.vycminpot[catalogue].value
        #self.vzminpot = properties.velocities.vzcminpot[catalogue].value

    def add_stellar_morphology(self,data,index):
        self.kappa_co[index] = data[0]
        self.momentum[index] = data[1]
        self.axis_ca[index] = data[2]
        self.axis_cb[index] = data[3]
        self.axis_ba[index] = data[4]

    def add_gas_morphology(self, data, index):
        self.gas_kappa_co[index] = data[0]
        self.gas_momentum[index] = data[1]
        self.gas_axis_ca[index] = data[2]
        self.gas_axis_cb[index] = data[3]
        self.gas_axis_ba[index] = data[4]

    def add_surface_density(self, data, index):
        self.sigma_H2[index] = data[0]
        self.sigma_gas[index] = data[1]
        self.sigma_SFR[index] = data[2]



def output_galaxy_data(galaxy_data,siminfo):

    sfr_galaxy = galaxy_data.star_formation_rate
    mass_galaxy = galaxy_data.stellar_mass
    gas_mass_galaxy = galaxy_data.gas_mass
    mass_halo = galaxy_data.halo_mass
    galaxy_metallicity_gas_sfr = galaxy_data.metallicity_gas_sfr
    galaxy_metallicity_gas = galaxy_data.metallicity_gas

    np.savetxt(f"{siminfo.output_path}/galaxy_data_"+siminfo.name+".txt",
               np.transpose([sfr_galaxy, mass_galaxy, gas_mass_galaxy,
                             mass_halo, galaxy_metallicity_gas_sfr, galaxy_metallicity_gas]))