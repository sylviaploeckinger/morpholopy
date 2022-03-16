from typing import List, Union, Tuple, Dict
import unyt
import numpy as np
import glob

from .unitilies import constants
from .unitilies.helper_functions import (
    calculate_kappa_co,
    AxialRatios,
    cosmic_time_approx_Gyr,
)
from .unitilies import luminosities as lum

from .halo_catalogue import HaloCatalogue
from .particle_ids import ParticleIds

from swiftsimio import load

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
        self.snapshot = load(f"{self.directory}/{self.snapshot_name}")

        # Fetch the run name if not provided
        if name is not None:
            self.simulation_name = name
        else:
            self.simulation_name = self.snapshot.metadata.run_name

        # Conversion from internal units to kpc
        self.to_kpc_units = (
            self.snapshot.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
            / constants.kpc
        )

        # Conversion from internal units to Msun
        self.to_Msun_units = (
            self.snapshot.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
            / constants.Msun
        )

        # Conversion from internal units to Myr
        self.to_Myr_units = (
            self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
            / constants.Myr
        )

        # Conversion from internal units to yr
        self.to_yr_units = (
            self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
            / constants.yr
        )

        # Box size of the simulation in kpc
        self.boxSize = self.snapshot.metadata.boxsize.to("kpc").value[0]

        # Cosmic scale factor
        self.a = self.snapshot.metadata.scale_factor

        self.hubble_time_Gyr = self.snapshot.metadata.cosmology.hubble_time.value

        self.Omega_m = self.snapshot.metadata.cosmology.Om0

        # No curvature
        self.Omega_l = 1.0 - self.Omega_m

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
            galaxy_min_stellar_mass=galaxy_min_stellar_mass,
        )

        # Init parent class with particle ids
        super().__init__(
            path_to_groups_file=f"{self.directory}/{self.catalogue_groups}",
            path_to_particles_file=f"{self.directory}/{self.catalogue_particles}",
            path_to_snapshot_file=f"{self.directory}/{self.snapshot_name}",
        )

        # Contained with spatially resolved data for combined plots
        self.combined_data = CombinedData()

        print(f"Data from run '{self.simulation_name}' has been loaded! \n")

        return

    def __find_groups_and_particles_catalogues(self) -> None:
        """
        Finds paths to the fields with particles catalogue and groups catalogue
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

        if len(stars_birthz) > 1:
            stars_age = cosmic_time_approx_Gyr(
                z=0.0, Omega_L=self.Omega_l, Hubble_time=self.hubble_time_Gyr
            ) - cosmic_time_approx_Gyr(
                z=stars_birthz, Omega_L=self.Omega_l, Hubble_time=self.hubble_time_Gyr
            )
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

    def write_galaxy_data_to_file(self, output_path: str) -> None:
        """
        Writes data with halo main properties to a file

        Parameters
        ----------
        output_path: str
        Path to output directory where the file will be saved
        """

        np.savetxt(
            f"{output_path}/galaxy_data_{self.simulation_name}.txt",
            np.transpose(
                [
                    self.halo_data.sfr,
                    self.halo_data.log10_stellar_mass,
                    self.halo_data.log10_gas_mass,
                    self.halo_data.log10_halo_mass,
                    self.halo_data.metallicity_gas_sfr,
                    self.halo_data.metallicity_gas,
                ]
            ),
        )

        return

    @classmethod
    def load_photometry_grid(cls):
        """
        Loads photometry grids for interpolation
        """

        if not cls.pgrids:
            system = "GAMA"  # hard-coded for now
            for pht in glob.glob(f"./photometry/{system}/*"):
                cls.pgrids[pht[-1]] = lum.MakeGrid(pht)
            print("Photometry tables have been loaded! \n")

    @classmethod
    def calculate_luminosities(cls, spart_data: np.ndarray):
        """
        Computes stellar luminosities
        Parameters
        ----------
        spart_data: np.ndarray
        Array with spart data

        Returns
        -------
        Output: np.ndarray
        An array with spart luminosities
        """

        assert cls.pgrids is not {}, (
            "'pgrids' is empty! Load photometry tables before calculating"
            " stellar luminosities!"
        )

        star_abmags = {}

        for filt in cls.pgrids.keys():
            grid, Z_p, t_p = cls.pgrids[filt]
            fluxes = (
                lum.BiPowInterp(t_p, Z_p, grid, spart_data[:, 9], spart_data[:, 10])
                * spart_data[:, 11]
            )
            star_abmags[filt] = -2.5 * np.log10(fluxes * lum.magfac)
        return star_abmags


class CombinedData:
    """
    Container holding spatially resolved data for combined plots
    """

    def __init__(self):

        # Surface densities (spatially-resolved)
        self.neutral_gas_surface_density: np.ndarray = np.array([])
        self.molecular_gas_surface_density: np.ndarray = np.array([])
        self.atomic_gas_surface_density: np.ndarray = np.array([])

        # SFR surface density (spatially-resolved)
        self.SFR_surface_density: np.ndarray = np.array([])

        # Depletion time scales corrected for corresponding surface densities (spatially-resolved)
        self.depletion_time_molecular_gas: np.ndarray = np.array([])
        self.depletion_time_neutral_gas: np.ndarray = np.array([])
        self.depletion_time_atomic_gas: np.ndarray = np.array([])

        # H_2 to H_2 + HI ratio (spatially-resolved) 
        self.H2_to_neutral_surface_density_ratio: np.ndarray = np.array([])

        # Metallicity map (spatially-resolved)
        self.gas_metallicity: np.ndarray = np.array([])

        # Surface densities (azimuthally-averaged) 
        self.radii_neutral_gas_surface_density: np.ndarray = np.array([])
        self.radii_atomic_gas_surface_density: np.ndarray = np.array([])
        self.radii_molecular_gas_surface_density: np.ndarray = np.array([])

        # SFR surface density (azimuthally-averaged)
        self.radii_SFR_surface_density: np.ndarray = np.array([])

        # Depletion times (azimuthally-averaged)
        self.radii_depletion_time_molecular_gas: np.ndarray = np.array([])
        self.radii_depletion_time_neutral_gas: np.ndarray = np.array([])
        self.radii_depletion_time_atomic_gas: np.ndarray = np.array([])

        # H_2 to H_2 + HI ratio (azimuthally-averaged) 
        self.radii_H2_to_neutral_surface_density_ratio: np.ndarray = np.array([])
