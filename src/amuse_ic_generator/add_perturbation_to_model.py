import argparse

import numpy as np

from amuse.units import units
from amuse.io import write_set_to_file, read_set_from_file


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "--perturbation",
        type=float,
        default=-8.0,
        help="log of the perturbation introduced",
    )
    parser.add_argument(
        "--parameter",
        default="x",
        help="parameter to be perturbed",
    )
    parser.add_argument(
        "--key",
        type=int,
        default=-1,
        help="key of the star to be perturbed",
    )
    parser.add_argument("--name", default="", help="name the perturbed star")
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
    if args.key > 0:
        perturber = bodies[bodies.key == args.key]
    else:
        perturber = bodies.random_sample(1)

    if args.parameter == "x":
        unit = perturber.x.unit
        perturber.x += 10 ** (args.perturbation) | unit
    elif args.parameter == "y":
        unit = perturber.y.unit
        perturber.y += 10 ** (args.perturbation) | unit
    elif args.parameter == "z":
        unit = perturber.z.unit
        perturber.z += 10 ** (args.perturbation) | unit
    elif args.parameter == "vx":
        unit = perturber.vx.unit
        perturber.vx += 10 ** (args.perturbation) | unit
    elif args.parameter == "vy":
        unit = perturber.vy.unit
        perturber.vy += 10 ** (args.perturbation) | unit
    elif args.parameter == "vz":
        unit = perturber.vz.unit
        perturber.vz += 10 ** (args.perturbation) | unit
    elif args.parameter == "mass":
        unit = perturber.mass.unit
        perturber.mass += 10 ** (args.perturbation) | unit

    if len(args.name) > 0:
        perturber.name = args.name

    time = 0 | units.Myr
    if args.outfile is None:
        namelist = args.filename.split(".")
        namelist[-2] = f"{namelist[-2]}.pert_{args.parameter}{int(args.perturbation)}"
        filename = namelist[0]
        for ni in namelist[1:]:
            filename += "." + ni
    else:
        filename = args.outfile

    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )
