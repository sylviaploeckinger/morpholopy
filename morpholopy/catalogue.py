import numpy as np
from velociraptor import load as load_catalogue

class HaloCatalogue:
    def __init__(self,siminfo, lower_mass):

        properties = load_catalogue(siminfo.subhalo_properties)
        stellar_mass = properties.masses.m_star_30kpc
        stellar_mass.convert_to_units("msun")
        gas_mass = properties.masses.m_gas_30kpc
        gas_mass.convert_to_units("msun")

        # Selecting galaxies more massive than lower limit
        catalogue = np.where(stellar_mass >= lower_mass)[0]

        # Selecting centrals only
        structure_type = properties.structure_type.structuretype.value
        centrals = np.where(structure_type[catalogue] == 10)[0]
        catalogue = catalogue[centrals]

        self.num = len(catalogue)
        self.halo_index = catalogue

        # Sample :
        self.stellar_mass = np.log10(stellar_mass[catalogue])
        self.gas_mass = np.log10(gas_mass[catalogue])

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

        # Subhalo data :
        self.xminpot = float(properties.positions.xcminpot[catalogue]) * 1e3
        self.yminpot = float(properties.positions.ycminpot[catalogue]) * 1e3
        self.zminpot = float(properties.positions.zcminpot[catalogue]) * 1e3

        self.vxminplot = float(properties.velocities.vxcminpot[catalogue])
        self.vyminplot = float(properties.velocities.vycminpot[catalogue])
        self.vzminplot = float(properties.velocities.vzcminpot[catalogue])

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
