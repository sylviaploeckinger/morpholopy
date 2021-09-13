import numpy as np
import h5py
from typing import Tuple


class Particles:

    gas_ids: np.ndarray
    star_ids: np.ndarray
    particle_ids_in_haloes: np.ndarray
    halo_ids: np.array

    def __init__(
        self,
        path_to_groups_file: str,
        path_to_particles_file: str,
        path_to_snapshot_file: str,
    ):

        group_file = h5py.File(path_to_groups_file, "r")
        particles_file = h5py.File(path_to_particles_file, "r")
        snapshot_file = h5py.File(path_to_snapshot_file, "r")

        self.star_ids = snapshot_file["/PartType4/ParticleIDs"][:]
        self.gas_ids = snapshot_file["/PartType0/ParticleIDs"][:]
        self.particle_ids_in_haloes = particles_file["Particle_IDs"][:]
        self.halo_ids = group_file["Offset"][:]

        group_file.close()
        particles_file.close()
        snapshot_file.close()

    def make_masks_gas_and_stars(self, halo_idx: int) -> Tuple[np.ndarray, np.ndarray]:

        halo_start_position = self.halo_ids[halo_idx]
        halo_end_position = self.halo_ids[halo_idx + 1]
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

        mask_gas = mask_gas[mask_gas > 0]
        mask_stars = mask_stars[mask_stars > 0]

        return mask_gas, mask_stars
