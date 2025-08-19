"""Give a random subset of stars a name

Takes an AMUSE file as input and gives a random subset of particles a name
"""
import argparse

import numpy as np

from amuse.units import units
from amuse.io import write_set_to_file, read_set_from_file


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument("--name", default="Sun", help="name of the stars")
    parser.add_argument(
        "--nstars",
        type=int,
        default=1,
        help="number of stars with name",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=-1,
        help="random number seed",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()

    if args.seed > 0:
        np.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(args.filename, close_file=True)
    random_stars = bodies.random_sample(args.nstars)
    random_stars.name = args.name

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "added_names_to_subset_of_stars.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
