import argparse

import numpy as np

from amuse.units import units, nbody_system
from amuse.datamodel import Particles
from amuse.io import write_set_to_file
from amuse.ic.salpeter import new_salpeter_mass_distribution


def make_cube_of_stars(masses, names, converter):
    N = len(masses)
    stars = Particles(N, mass=masses)
    # L = converter.to_si(1|nbody_system.length)
    # dv = converter.to_si(40|nbody_system.speed)
    L = 500 | units.au
    dv = 2.5 | units.kms
    np.random.seed(7654304)

    stars.x = L * np.random.uniform(-1.0, 1.0, N)
    stars.y = L * np.random.uniform(-1.0, 1.0, N)
    stars.z = L * 0.0
    stars.vx = dv * np.random.uniform(-1.0, 1.0, N)
    stars.vy = dv * np.random.uniform(-1.0, 1.0, N)
    stars.vz = dv * 0.0

    stars.radius = (1.0 | units.RSun) * (stars.mass / (1.0 | units.MSun)) ** (1.0 / 3.0)

    stars.type = "star"
    stars.name = names
    stars.mass = masses
    stars.move_to_center()
    return stars


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "--nstars",
        type=int,
        default=10,
        help="number of stars",
    )
    parser.add_argument(
        "--mass",
        type=units.MSun,
        default=1 | units.MSun,
        help="total cluster mass",
    )
    parser.add_argument(
        "--mmin",
        type=units.MSun,
        default=-0.1 | units.MSun,
        help="minimum stellar mass",
    )
    parser.add_argument(
        "--mmax",
        type=units.MSun,
        default=100 | units.MSun,
        help="maximum stellar mass",
    )
    parser.add_argument("-Q", "--Qvir", type="float", default=0.5, help="virial ratio")
    parser.add_argument(
        "--radius",
        type=units.parsec,
        default=1.0 | units.parsec,
        help="cluster radius",
    )
    parser.add_argument("--name", default="star", help="stellar name")
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
    # np.random.seed(7654304)

    if args.mmin > 0 | units.MSun:
        masses = new_salpeter_mass_distribution(args.nstars, args.mmin, args.mmax)
    else:
        masses = np.zeros(args.nstars) | units.MSun
        masses += args.mass

    converter = nbody_system.nbody_to_si(masses.sum(), args.radius)
    bodies = make_cube_of_stars(masses, args.name, converter)
    print(bodies)

    # bodies.scale_to_standard(convert_nbody=converter, virial_ratio=args.Qvir)
    time = 0 | units.Myr
    if args.outfile is None:
        filename = "cube.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
