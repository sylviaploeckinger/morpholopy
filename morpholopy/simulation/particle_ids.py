import numpy as np
import h5py
from typing import Tuple


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

    def make_mask_gas(self, halo_id: int) -> Tuple[np.ndarray, np.ndarray]:
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

        _, _, mask_gas = np.intersect1d(
            particle_ids_in_halo,
            self.gas_ids,
            assume_unique=True,
            return_indices=True,
        )

        # Ensure that there are no negative indices
        mask_gas = mask_gas[mask_gas > 0]

        return mask_gas

    def make_mask_stars(self, halo_id: int) -> Tuple[np.ndarray, np.ndarray]:
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

        # Ensure that there are no negative indices
        mask_stars = mask_stars[mask_stars > 0]

        return mask_stars

