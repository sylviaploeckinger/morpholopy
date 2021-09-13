from typing import List, Union, Tuple
from unyt import unit_object
import numpy as np
import glob

from .unitilies import constants
from .unitilies.functions import calculate_kappa_co, AxialRatios
from .unitilies import luminosities as lum

from .halo_catalogue import HaloCatalogue
from .particles import Particles

from swiftsimio import load, SWIFTDataset
from astropy.cosmology import WMAP9 as cosmo


class SimInfo(Particles):

    directory: str
    snapshot_name: str
    catalogue_name: str
    simulation_name: Union[str, None]
    catalogue_groups: str
    catalogue_particles: str

    a: float
    baryon_max_soft: float
    box_size: float

    snapshot: SWIFTDataset

    to_kpc_units: float
    to_Msun_units: float

    Zsolar = 0.0134

    halo_data: HaloCatalogue

    def __init__(
        self,
        directory: str,
        snapshot: str,
        catalogue: str,
        name: Union[str, None],
        galaxy_min_stellar_mass: unit_object,
    ):

        self.directory = directory
        self.snapshot_name = snapshot
        self.catalogue_name = catalogue
        self.simulation_name = name

        self.__find_groups_and_particles_catalogues()

        self.snapshot = load(f"{self.directory}/{self.snapshot_name}")

        self.to_kpc_units = (
            self.snapshot.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
            / constants.kpc
        )

        self.to_Msun_units = (
            self.snapshot.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
            / constants.Msun
        )

        self.to_Myr_units = (
            self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
            / constants.Myr
        )

        self.to_yr_units = (
            self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
            / constants.yr
        )

        self.boxSize = self.snapshot.metadata.boxsize.to("kpc").value[0]
        self.a = self.snapshot.metadata.scale_factor
        self.baryon_max_soft = (
            self.snapshot.metadata.gravity_scheme[
                "Maximal physical baryon softening length  [internal units]"
            ][0]
            * self.to_kpc_units
        )

        self.halo_data = HaloCatalogue(
            path_to_catalogue=f"{self.directory}/{self.catalogue_name}",
            galaxy_min_stellar_mass=galaxy_min_stellar_mass,
        )

        super().__init__(
            path_to_groups_file=f"{self.directory}/{self.catalogue_groups}",
            path_to_particles_file=f"{self.directory}/{self.catalogue_particles}",
            path_to_snapshot_file=f"{self.directory}/{self.snapshot_name}",
        )

    def __find_groups_and_particles_catalogues(self) -> None:
        """
        Finds paths to the fiels with particles catalogue and groups catalogue
        """

        catalogue_num = "".join([s for s in self.catalogue_name if s.isdigit()])
        catalogue_groups_paths: List[str] = glob.glob(
            f"{self.directory}/*{catalogue_num}.catalog_groups*"
        )
        catalogue_particles_paths: List[str] = glob.glob(
            f"{self.directory}/*{catalogue_num}.catalog_particles*"
        )

        # We expect one file for particle groups
        if len(catalogue_groups_paths) == 1:
            self.catalogue_groups = catalogue_groups_paths[0].split("/")[-1]
        else:
            raise IOError("Couldn't find catalogue_groups file")

        # We expect two files: one for bound and the other for unbound particles
        if len(catalogue_particles_paths) == 2:
            for path in catalogue_particles_paths:
                if path.find("unbound") == -1:
                    self.catalogue_particles = path.split("/")[-1]
        else:
            raise IOError("Couldn't find catalogue_particles file")

    def make_particle_data(self, halo_id: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Create particle data :
        [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz
        | (7) hsml ]

        Parameters
        ----------
        halo_id: int
        Halo id from the halo catalogue

        Returns
        -------
        Output: Tuple[np.ndarray, np.ndarray]
        Arrays with gas and stellar properties
        """

        mask_gas, mask_stars = self.make_masks_gas_and_stars(halo_idx=halo_id)

        gas_mass = self.snapshot.gas.masses[mask_gas].value * self.to_Msun_units
        gas_n_parts = len(gas_mass)
        gas_data = np.zeros((gas_n_parts, 13))

        gas_data[:, 0:3] = (
            self.snapshot.gas.coordinates[mask_gas].value * self.a * self.to_kpc_units
        )

        gas_data[:, 3] = gas_mass
        gas_data[:, 4:7] = self.snapshot.gas.velocities[mask_gas].value  # km/s
        gas_data[:, 7] = (
            self.snapshot.gas.smoothing_lengths[mask_gas].value
            * self.a
            * self.to_kpc_units
        )

        XH = self.snapshot.gas.element_mass_fractions.hydrogen[mask_gas].value
        gas_HI = self.snapshot.gas.species_fractions.HI[mask_gas].value
        gas_H2 = self.snapshot.gas.species_fractions.H2[mask_gas].value * 2.0

        gas_data[:, 8] = gas_HI * XH * gas_mass
        gas_data[:, 9] = gas_H2 * XH * gas_mass

        gas_data[:, 10] = (
            self.snapshot.gas.star_formation_rates[mask_gas].value
            * self.to_Msun_units
            / self.to_yr_units
        )
        gas_data[:, 11] = (
            self.snapshot.gas.densities[mask_gas].value
            * (self.a * self.to_Msun_units / self.to_kpc_units) ** 3
        )
        gas_data[:, 12] = (
            self.snapshot.gas.metal_mass_fractions[mask_gas].value / self.Zsolar
        )

        stars_mass = self.snapshot.stars.masses[mask_stars].value * self.to_Msun_units
        stars_birthz = 1.0 / self.snapshot.stars.birth_scale_factors[mask_stars].value - 1.0

        if len(stars_birthz) > 1:
            stars_age = cosmo.age(0.0).value - cosmo.age(stars_birthz).value  # in Gyr
        else:
            stars_age = 0.0

        stars_Z = self.snapshot.stars.metal_mass_fractions[mask_stars].value
        stars_initmass = (
            self.snapshot.stars.initial_masses[mask_stars].value * self.to_Msun_units
        )

        stars_n_parts = len(stars_mass)
        stars_data = np.zeros((stars_n_parts, 12))
        stars_data[:, 0:3] = (
            self.snapshot.stars.coordinates[mask_stars].value
            * self.a
            * self.to_kpc_units
        )
        stars_data[:, 3] = stars_mass
        stars_data[:, 4:7] = self.snapshot.stars.velocities[
            mask_stars
        ].value  # km/s  # km/s
        stars_data[:, 7] = 0.5 * self.baryon_max_soft * np.ones(stars_mass.size)
        stars_data[:, 8] = stars_mass * (1.2348 / stars_data[:, 7]) ** 3
        stars_data[:, 9] = stars_age
        stars_data[:, 10] = stars_Z
        stars_data[:, 11] = stars_initmass

        return gas_data, stars_data

    def calculate_morphology(
        self, part_data: np.ndarray, halo_index: int, parttype: int
    ):
        """
        Computes morphological properties of a given halo

        Parameters
        ----------
        part_data: ndarray
        Array with gas or stellar properties

        halo_index: int
        Halo index in the halo catalogue

        parttype: int
        Particle type: gas (0) or stars (4)
        """

        # Calculate kappa and specific angular momentum
        kappa, specific_momentum, momentum, part_data = calculate_kappa_co(
            halo_data=self.halo_data,
            partsDATA=part_data,
            box_size=self.boxSize,
            halo_index=halo_index,
        )

        # Calculate axis ratios
        axis_1, axis_2, axis_3 = AxialRatios(part_data[:, :3], part_data[:, 3])
        morphology = np.array([kappa, specific_momentum, axis_1, axis_2, axis_3])

        # Store morphology parameters in halo data and continue
        if parttype == 4:
            self.halo_data.add_stellar_morphology(morphology, halo_index)
        elif parttype == 0:
            self.halo_data.add_gas_morphology(morphology, halo_index)
        return momentum, part_data

    def output_galaxy_data(self, output_path: str) -> None:
        """
        Parameters
        ----------

        output_path: str
        Path to output directory
        """

        sfr_galaxy = self.halo_data.sfr
        mass_galaxy = self.halo_data.stellar_mass
        gas_mass_galaxy = self.halo_data.gas_mass
        mass_halo = self.halo_data.halo_mass
        galaxy_metallicity_gas_sfr = self.halo_data.metallicity_gas_sfr
        galaxy_metallicity_gas = self.halo_data.metallicity_gas

        np.savetxt(
            f"{output_path}/galaxy_data_" + self.simulation_name + ".txt",
            np.transpose(
                [
                    sfr_galaxy,
                    mass_galaxy,
                    gas_mass_galaxy,
                    mass_halo,
                    galaxy_metallicity_gas_sfr,
                    galaxy_metallicity_gas,
                ]
            ),
        )

        return

    @staticmethod
    def calculate_luminosities(spart_data: np.ndarray, pgrids):
        """
        Computes stellar luminosities
        Parameters
        ----------
        spart_data: np.ndarray

        An array with spart properties
        """

        star_abmags = {}

        for filt in pgrids.keys():
            grid, Z_p, t_p = pgrids[filt]
            fluxes = (
                lum.BiPowInterp(t_p, Z_p, grid, spart_data[:, 9], spart_data[:, 10])
                * spart_data[:, 11]
            )
            star_abmags[filt] = -2.5 * np.log10(fluxes * lum.magfac)
        return star_abmags
