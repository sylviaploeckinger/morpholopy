"""
Some helper functions .

Contains the argument parser and default parsing results.

See the README for available argument parsers.
"""

import argparse as ap
import glob
from typing import Optional, List



parser = ap.ArgumentParser(
    description="""General argument parser for isolated galaxy scripts."""
)

parser.add_argument(
    "-d",
    "--directory",
    help="Directory containing snapshots. Required.",
    type=str,
    required=True,
    nargs="*",
)

parser.add_argument(
    "-s",
    "--snapshot",
    help="Snapshot number to visualise. Required.",
    type=str,
    required=True,
    nargs="*",
)

parser.add_argument(
    "-n",
    "--run-names",
    help="Names of the runs for placement in legends.",
    type=str,
    required=True,
    nargs="*",
)

parser.add_argument(
    "-o",
    "--output",
    help="Output directory for the figures. If not present, the same directory as the snapshots are in is used.",
    required=True,
    type=str,
    default=None,
)

parser.add_argument(
    "-z",
    "--zoom",
    help="Option for whether the snapshot corresponds to a zoom simulation",
    required=False,
    nargs="*",
    type=str,
    default="no",
)

args = parser.parse_args()

if args.output is None:
    args.output = args.directory[0]

