import numpy as np
from unyt import unit_object, unyt_quantity
from velociraptor import load


class HaloCatalogue:
    """
    General class containing halo properties
    """

    number_of_haloes: int
    path_to_catalogue: str

    def __init__(self, path_to_catalogue: str, galaxy_min_stellar_mass: unit_object):

        self.path_to_catalogue = path_to_catalogue

        catalogue = load(self.path_to_catalogue)

        # Selecting galaxies more massive than lower limit
        mask = np.logical_and(
            catalogue.apertures.mass_star_30_kpc >= galaxy_min_stellar_mass,
            catalogue.structure_type.structuretype == 10,
        )
        mask = np.logical_and(
            mask, catalogue.apertures.mass_gas_30_kpc > unyt_quantity(0.0, "Msun")
        )

        self.number_of_haloes = mask.sum()

        # Masses
        self.stellar_mass = np.log10(
            catalogue.apertures.mass_star_30_kpc.to("Msun").value[mask]
        )
        self.gas_mass = np.log10(
            catalogue.apertures.mass_gas_30_kpc.to("Msun").value[mask]
        )
        self.halo_mass = np.log10(catalogue.masses.mass_200crit.to("Msun").value[mask])

        # Radii
        self.half_mass_radius_star = catalogue.radii.r_halfmass_star.to("kpc").value[
            mask
        ]
        self.half_mass_radius_gas = catalogue.radii.r_halfmass_gas.to("kpc").value[mask]

        # Metallicity and star formation rate
        self.sfr = (
            catalogue.apertures.sfr_gas_30_kpc.value[mask] * 10227144.8879616 / 1e9
        )
        self.metallicity_gas_sfr = catalogue.apertures.zmet_gas_sf_30_kpc.value[mask]

        self.metallicity_gas = catalogue.apertures.zmet_gas_30_kpc.value[mask]

        self.halo_ids = np.array([i for i in range(len(mask)) if mask[i] == True])

        self.kappa_co = np.zeros(self.number_of_haloes)
        self.momentum = np.zeros(self.number_of_haloes)
        self.axis_ca = np.zeros(self.number_of_haloes)
        self.axis_cb = np.zeros(self.number_of_haloes)
        self.axis_ba = np.zeros(self.number_of_haloes)

        self.gas_kappa_co = np.zeros(self.number_of_haloes)
        self.gas_momentum = np.zeros(self.number_of_haloes)
        self.gas_axis_ca = np.zeros(self.number_of_haloes)
        self.gas_axis_cb = np.zeros(self.number_of_haloes)
        self.gas_axis_ba = np.zeros(self.number_of_haloes)

        self.sigma_H2 = np.zeros(self.number_of_haloes)
        self.sigma_gas = np.zeros(self.number_of_haloes)
        self.sigma_SFR = np.zeros(self.number_of_haloes)

        self.surface_density = np.zeros(self.number_of_haloes)
        self.SFR_density = np.zeros(self.number_of_haloes)
        self.ratio_densities = np.zeros(self.number_of_haloes)
        self.metallicity = np.zeros(self.number_of_haloes)

        self.radii_surface_density = np.zeros(self.number_of_haloes)
        self.radii_surface_ratio = np.zeros(self.number_of_haloes)
        self.radii_nbins = np.zeros(self.number_of_haloes)

        self.xminpot = catalogue.positions.xcminpot.to("kpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("kpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("kpc").value[mask]

        self.vxminpot = catalogue.velocities.vxcminpot.to("km/s").value[mask]
        self.vyminpot = catalogue.velocities.vycminpot.to("km/s").value[mask]
        self.vzminpot = catalogue.velocities.vzcminpot.to("km/s").value[mask]

    def add_stellar_morphology(self, data, index):
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
