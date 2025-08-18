import argparse

import numpy as np

from amuse.io import read_set_from_file, write_set_to_file
from amuse.units import units


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "--pos",
        "--position",
        type=units.pc,
        nargs=3,
        default=[-8340.0, 0.0, 27.0] | units.pc,
        help="move com position here",
    )
    parser.add_argument(
        "--vel",
        dest="velocity",
        type=units.kms,
        nargs=3,
        default=[11.1, 240.0, 7.25] | units.kms,
        help="move com-velocity here",
    )
    parser.add_argument(
        "-t",
        "--timestamp",
        type=units.Myr,
        default=0 | units.Myr,
        help="set timestamp",
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

    bodies.move_to_center()
    bodies.position += args.position | units.parsec
    bodies.velocity += args.velocity | units.kms

    print(bodies)
    time = args.timestamp
    if args.outfile is None:
        filename = "moved_snapshot.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
