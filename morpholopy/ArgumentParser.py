import argparse
from typing import List, Optional


class ArgumentParser(object):

    # List of snapshots that are to be processed.
    snapshot_list: List[str]
    # List of catalogues (in the same order as snapshots) to be processed.
    catalogue_list: List[str]
    # List of directories that contain the above snapshots and catalogues
    directory_list: List[str]
    # List of representative names for the snapshots; may be a list of Nones
    name_list: List[Optional[str]]
    # Directory to output the figure to
    output_directory: str
    # Number of inputs to the script.
    number_of_inputs: int
    # Number of galaxies to make individual plots for
    number_of_galaxies: int

    def __init__(self):

        parser = argparse.ArgumentParser(
            description="""General argument parser for morphology pipeline."""
        )

        parser.add_argument(
            "-s",
            "--snapshots",
            help="Snapshot list. Do not include directory. Example: snapshot_0000.hdf5",
            type=str,
            required=True,
            nargs="*",
        )

        parser.add_argument(
            "-c",
            "--catalogues",
            help=(
                "Catalogue list. Do not include directory. Example: "
                "catalogue_0000.properties"
            ),
            type=str,
            required=True,
            nargs="*",
        )

        parser.add_argument(
            "-d",
            "--input-directories",
            help="Input directory list. Catalogue and snapshot are in this directory.",
            type=str,
            required=True,
            nargs="*",
        )

        parser.add_argument(
            "-n",
            "--run-names",
            help="Names of the runs for placement in legends.",
            type=str,
            required=False,
            nargs="*",
        )

        parser.add_argument(
            "-o",
            "--output-directory",
            help="Output directory for the produced figure.",
            type=str,
            required=True,
        )

        parser.add_argument(
            "-g",
            "--galaxy_number",
            help="Option for outputting morphology plots of g number of galaxies",
            required=False,
            nargs="?",
            type=int,
            default=10,
        )

        args = parser.parse_args()

        self.snapshot_list = args.snapshots
        self.catalogue_list = args.catalogues
        self.directory_list = args.input_directories
        self.name_list = (
            args.run_names
            if args.run_names is not None
            else [None] * self.number_of_inputs
        )
        self.output_directory = args.output_directory

        self.number_of_inputs = len(args.snapshots)
        self.number_of_galaxies = args.galaxy_number
