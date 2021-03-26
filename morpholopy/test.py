"""
Description here
"""

import os
import h5py
from plotter.html import make_web, add_web_section, render_web, PlotsInPipeline, add_metadata_to_web
import unyt


class SimInfo:
    def __init__(self, folder, snap, output_path, comparison, name):
        self.comparison = comparison
        self.name = name
        self.output_path = output_path
        self.snapshot = os.path.join(folder,"colibre_%04i.hdf5"%snap)
        self.subhalo_properties = os.path.join(folder,"halo_%04i.properties.0"%snap)
        self.catalog_groups = os.path.join(folder,"halo_%04i.catalog_groups.0"%snap)
        self.catalog_particles = os.path.join(folder, "halo_%04i.catalog_particles.0" % snap)
        snapshot_file = h5py.File(self.snapshot, "r")
        self.boxSize = snapshot_file["/Header"].attrs["BoxSize"][0] * 1e3 #kpc
        self.a = snapshot_file["/Header"].attrs["Scale-factor"]
        self.baryon_maxsoft = snapshot_file["/GravityScheme"].attrs['Maximal physical baryon softening length  [internal units]'] * 1e3 #kpc


if __name__ == '__main__':
    from utils import *

    # Load MorpholoPy production details
    output_path = args.output
    number_of_inputs = len(args.snapshot)
    directory_list = args.directory
    snapshot_list = args.snapshot
    name_list = (
        args.run_names
        if args.run_names is not None
        else [None] * number_of_inputs
    )

    # Are we comparing?
    if number_of_inputs > 1:
        comparison = True
    else: comparison = False

    # Loop over simulation list
    for sims in range(number_of_inputs):

        directory = directory_list[sims]
        snap_number = int(snapshot_list[sims])
        sim_name = name_list[sims]
        siminfo = SimInfo(directory, snap_number,
                          output_path, comparison, sim_name)

        # Loading simulation data in website table
        if sims==0:web = make_web(siminfo)
        if sims==1:
            add_metadata_to_web(web, siminfo)
            render_web(web, siminfo.output_path)
