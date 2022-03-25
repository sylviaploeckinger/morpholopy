import numpy as np
import unyt
from velociraptor import load
import glob


class HaloCatalogue:
    """
    General class containing halo properties
    """

    def __init__(
        self, path_to_catalogue: str, galaxy_min_stellar_mass: unyt.array.unyt_quantity
    ):
        """
        Parameters
        ----------
        path_to_catalogue: str
        Path to the catalogue with halo properties

        galaxy_min_stellar_mass: unyt.array.unyt_quantity
        Minimum stellar mass in units of Msun. Objects whose stellar mass is lower than this
        threshold are disregarded
        """

        self.path_to_catalogue = path_to_catalogue

        # Load catalogue using velociraptor python library
        # sometimes VR produces files with trailing 0
        fcat = glob.glob(f"{self.path_to_catalogue}*")[0]
        catalogue = load(fcat)

        # Selecting central galaxies whose stellar mass is larger than
        # 'galaxy_min_stellar_mass'
        mask = np.logical_and(
            catalogue.apertures.mass_star_30_kpc >= galaxy_min_stellar_mass,
            catalogue.structure_type.structuretype == 10,
        )
        # They also need to contain at least one gas particle
        mask = np.logical_and(
            mask, catalogue.apertures.mass_gas_30_kpc > unyt.unyt_quantity(0.0, "Msun")
        )
        # Compute the number of haloes following the selection mask
        self.number_of_haloes = mask.sum()

        # Log10 stellar mass in units of Msun
        self.log10_stellar_mass = np.log10(
            catalogue.apertures.mass_star_30_kpc.to("Msun").value[mask]
        )
        # Log10 gas mass in units of Msun
        self.log10_gas_mass = np.log10(
            catalogue.apertures.mass_gas_30_kpc.to("Msun").value[mask]
        )
        # Log10 halo mass in units of Msun
        self.log10_halo_mass = np.log10(
            catalogue.masses.mass_200crit.to("Msun").value[mask]
        )

        # Half mass radius in units of kpc (stars)
        self.half_mass_radius_star = catalogue.radii.r_halfmass_star.to("kpc").value[
            mask
        ]
        # Half mass radius in units of kpc (gas)
        self.half_mass_radius_gas = catalogue.radii.r_halfmass_gas.to("kpc").value[mask]

        # Star formation rate in units of Msun/yr
        self.sfr = (
            catalogue.apertures.sfr_gas_30_kpc.value[mask] * 10227144.8879616 / 1e9
        )

        # Metallicity of star-forming gas
        self.metallicity_gas_sfr = catalogue.apertures.zmet_gas_sf_30_kpc.value[mask]

        # Metallicity of all gas
        self.metallicity_gas = catalogue.apertures.zmet_gas_30_kpc.value[mask]

        # Ids of haloes satisfying the selection criterion
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

        self.HI_size = np.zeros(self.number_of_haloes)
        self.HI_mass = np.zeros(self.number_of_haloes)

        self.xminpot = catalogue.positions.xcminpot.to("kpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("kpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("kpc").value[mask]

        self.vxminpot = catalogue.velocities.vxcminpot.to("km/s").value[mask]
        self.vyminpot = catalogue.velocities.vycminpot.to("km/s").value[mask]
        self.vzminpot = catalogue.velocities.vzcminpot.to("km/s").value[mask]

    def add_stellar_morphology(self, data, index):
        """
        Add stellar morphology data

        @TODO rewrite this function such that it is obvious what arguments are passed in.
        """
        self.kappa_co[index] = data[0]
        self.momentum[index] = data[1]
        self.axis_ca[index] = data[2]
        self.axis_cb[index] = data[3]
        self.axis_ba[index] = data[4]

        return

    def add_gas_morphology(self, data, index):
        """
        Add gas morphology data

        @TODO rewrite this function such that it is obvious what arguments are passed in.
        """
        self.gas_kappa_co[index] = data[0]
        self.gas_momentum[index] = data[1]
        self.gas_axis_ca[index] = data[2]
        self.gas_axis_cb[index] = data[3]
        self.gas_axis_ba[index] = data[4]

        return

    def add_surface_density(self, data, index):
        """
        Add surface density data

        @TODO rewrite this function such that it is obvious what arguments are passed in.
        """
        self.sigma_H2[index] = data[0]
        self.sigma_gas[index] = data[1]
        self.sigma_SFR[index] = data[2]

        return

    def add_HI_size_mass(self, HI_size, HI_mass, index):
        self.HI_size[index] = HI_size
        self.HI_mass[index] = HI_mass
