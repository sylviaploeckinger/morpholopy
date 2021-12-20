from swiftsimio import load as swift_load
from velociraptor import load as velociraptor_load
from typing import List, Union, Tuple, Dict
import unyt
import numpy as np
import glob
import h5py

# Constants
kpc = 3.08567758e21
Msun = 1.9891e33
yr = 3.1556926e7
Myr = yr * 1e6


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
        catalogue = velociraptor_load(self.path_to_catalogue)

        mask = np.logical_and(
            catalogue.apertures.mass_star_30_kpc >= galaxy_min_stellar_mass,
            catalogue.apertures.mass_gas_30_kpc > unyt.unyt_quantity(0.0, "Msun")
        )

        # Compute the number of haloes following the selection mask
        self.number_of_haloes = mask.sum()

        # Log10 stellar mass in units of Msun
        self.log10_stellar_mass = np.log10(
            catalogue.apertures.mass_star_30_kpc.to("Msun").value[mask]
        )

        # Log10 halo mass in units of Msun
        self.log10_halo_mass = np.log10(
            catalogue.masses.mass_200crit.to("Msun").value[mask]
        )

        # Galaxy type, either central (=10) or satellite (>10)
        self.type = catalogue.structure_type.structuretype.value[mask]

        # Ids of haloes satisfying the selection criterion
        self.halo_ids = np.array([i for i in range(len(mask)) if mask[i] == True])

        self.xminpot = catalogue.positions.xcminpot.to("kpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("kpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("kpc").value[mask]

        self.vxminpot = catalogue.velocities.vxcminpot.to("km/s").value[mask]
        self.vyminpot = catalogue.velocities.vycminpot.to("km/s").value[mask]
        self.vzminpot = catalogue.velocities.vzcminpot.to("km/s").value[mask]

class ParticleIds:
    """
    A class providing the mapping between
    particle ids and halo ids
    """

    def __init__(
        self,
        path_to_groups_file: str,
        path_to_particles_file: str,
        path_to_snapshot_file: str,
    ):

        # Fetch ids
        group_file = h5py.File(path_to_groups_file, "r")
        particles_file = h5py.File(path_to_particles_file, "r")
        snapshot_file = h5py.File(path_to_snapshot_file, "r")

        # Ids of stellar particles from snapshot
        self.star_ids = snapshot_file["/PartType4/ParticleIDs"][:]
        # Ids of gas particles from snapshot
        self.gas_ids = snapshot_file["/PartType0/ParticleIDs"][:]
        # Particle ids from halo catalogue
        self.particle_ids_in_haloes = particles_file["Particle_IDs"][:]
        # Halo ids from group catalogue
        self.halo_ids = group_file["Offset"][:]

        group_file.close()
        particles_file.close()
        snapshot_file.close()

    def make_masks_gas_and_stars(self, halo_id: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Find gas and stellar particle ids that belong to a halo with the provided id

        Parameters
        ----------
        halo_id: int
        Halo id from the catalogue

        Returns
        -------
        Output: Tuple[np.ndarray, np.ndarray]
        A tuple containing ids of the stellar particles and gas particles
        """

        halo_start_position = self.halo_ids[halo_id]
        halo_end_position = self.halo_ids[halo_id + 1]
        particle_ids_in_halo = self.particle_ids_in_haloes[
            halo_start_position:halo_end_position
        ]

        _, _, mask_stars = np.intersect1d(
            particle_ids_in_halo,
            self.star_ids,
            assume_unique=True,
            return_indices=True,
        )

        _, _, mask_gas = np.intersect1d(
            particle_ids_in_halo,
            self.gas_ids,
            assume_unique=True,
            return_indices=True,
        )

        # Ensure that there are no negative indices
        mask_gas = mask_gas[mask_gas > 0]
        mask_stars = mask_stars[mask_stars > 0]

        return mask_gas, mask_stars

class SimInfo(ParticleIds):

    # Solar metallicity
    Zsolar = 0.0134

    # Dict with photometry tables to compute stars' luminosities
    pgrids: Dict = {}

    def __init__(
        self,
        directory: str,
        snapshot: str,
        catalogue: str,
        name: Union[str, None],
        galaxy_min_stellar_mass: unyt.array.unyt_quantity,
    ):
        """
        Parameters
        ----------

        directory: str
        Run directory

        snapshot: str
        Name of the snapshot file

        catalogue: str
        Name of the catalogue file

        name:
        Name of the run

        galaxy_min_stellar_mass: unyt.array.unyt_quantity
        """

        self.directory = directory
        self.snapshot_name = snapshot
        self.catalogue_name = catalogue

        # Find the group and particle catalogue files
        self.__find_groups_and_particles_catalogues()

        # Load snapshot via swiftsimio
        self.snapshot = swift_load(f"{self.directory}/{self.snapshot_name}")

        # Fetch the run name if not provided
        if name is not None:
            self.simulation_name = name
        else:
            self.simulation_name = self.snapshot.metadata.run_name

        # Conversion from internal units to kpc
        self.to_kpc_units = (
            self.snapshot.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
            / kpc
        )

        # Conversion from internal units to Msun
        self.to_Msun_units = (
            self.snapshot.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
            / Msun
        )

        # Conversion from internal units to Myr
        self.to_Myr_units = (
            self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
            / Myr
        )

        # Conversion from internal units to yr
        self.to_yr_units = (
            self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
            / yr
        )

        # Box size of the simulation in kpc
        self.boxSize = self.snapshot.metadata.boxsize.to("kpc").value[0]

        # Cosmic scale factor
        self.a = self.snapshot.metadata.scale_factor

        self.hubble_time_Gyr = self.snapshot.metadata.cosmology.hubble_time.value

        self.Omega_m = self.snapshot.metadata.cosmology.Om0

        # No curvature
        self.Omega_l = self.Omega_m

        # Maximum softening for baryons
        self.baryon_max_soft = (
            self.snapshot.metadata.gravity_scheme[
                "Maximal physical baryon softening length  [internal units]"
            ][0]
            * self.to_kpc_units
        )

        # Object containing halo properties (from halo catalogue)
        self.halo_data = HaloCatalogue(
            path_to_catalogue=f"{self.directory}/{self.catalogue_name}",
            galaxy_min_stellar_mass=1e6,
        )

        # Init parent class with particle ids
        super().__init__(
            path_to_groups_file=f"{self.directory}/{self.catalogue_groups}",
            path_to_particles_file=f"{self.directory}/{self.catalogue_particles}",
            path_to_snapshot_file=f"{self.directory}/{self.snapshot_name}",
        )

        return

    def __find_groups_and_particles_catalogues(self) -> None:
        """
        Finds paths to the fields with particles catalogue and groups catalogue
        """

        catalogue_num = "".join([s for s in self.catalogue_name.split('.')[0] if s.isdigit()])
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

        return

    def make_particle_data(self, halo_id: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes and saves gas and stellar particle data into numpy arrays.

        gas data:
        0-2 coordinates
        3 masses
        4-6 velocities
        7 smoothing lengths
        8 HI masses
        9 H2 masses
        10 SFR
        11 densities
        12 metallicities

        star data:
        0-2 coordinates
        3 masses
        4-7 velocities
        7 softenings
        8 Mass over softenning cubed
        9 stellar ages
        10 stellar metallicities
        11 stellar initial masses

        @TODO Use classes for the gas and stellar data containers, not numpy arrays

        Parameters
        ----------
        halo_id: int
        Halo id from the halo catalogue

        Returns
        -------
        Output: Tuple[np.ndarray, np.ndarray]
        Arrays with gas and stellar properties
        """

        mask_gas, mask_stars = self.make_masks_gas_and_stars(halo_id=halo_id)

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
        stars_birthz = (
            1.0 / self.snapshot.stars.birth_scale_factors[mask_stars].value - 1.0
        )

        stars_age = 0.0

        stars_Z = self.snapshot.stars.metal_mass_fractions[mask_stars].value
        stars_initmass = (
            self.snapshot.stars.initial_masses[mask_stars].value * self.to_Msun_units
        )

        stars_n_parts = len(stars_mass)
        stars_data = np.zeros((stars_n_parts, 20))
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

        stars_data[:, 12] = self.snapshot.stars.element_mass_fractions.oxygen[mask_stars].value
        stars_data[:, 13] = self.snapshot.stars.element_mass_fractions.iron[mask_stars].value
        stars_data[:, 14] = self.snapshot.stars.element_mass_fractions.magnesium[mask_stars].value
        stars_data[:, 15] = self.snapshot.stars.element_mass_fractions.hydrogen[mask_stars].value
        stars_data[:, 16] = self.snapshot.stars.element_mass_fractions.carbon[mask_stars].value
        stars_data[:, 17] = self.snapshot.stars.element_mass_fractions.silicon[mask_stars].value
        stars_data[:, 18] = self.snapshot.stars.element_mass_fractions.europium[mask_stars].value
        stars_data[:, 19] = np.zeros(stars_n_parts)

        indx = self.halo_data.halo_ids == halo_id
        x = stars_data[:, 0] - self.halo_data.xminpot[indx]
        y = stars_data[:, 1] - self.halo_data.yminpot[indx]
        z = stars_data[:, 2] - self.halo_data.zminpot[indx]
        r = x**2 + y**2 + z**2
        halo_stars = np.where(r > 10**2)[0] #further than 8kpc?
        stars_data[halo_stars, 19] = np.ones(len(halo_stars))

        return gas_data, stars_data