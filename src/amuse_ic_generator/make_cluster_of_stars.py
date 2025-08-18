import argparse

import numpy as np

from amuse.units import units, nbody_system
from amuse.io import write_set_to_file, read_set_from_file
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.community.fractalcluster import new_fractal_cluster_model


def make_plummer_sphere(nstars, masses, name, converter):
    stars = new_plummer_model(nstars, converter)
    stars.type = "star"
    stars.name = name
    stars.mass = masses
    stars.move_to_center()
    return stars


def make_fractal_cluster(masses, name, converter, Fd=1.6):
    stars = new_fractal_cluster_model(
        N=len(masses), convert_nbody=converter, do_scale=False, fractal_dimension=Fd
    )
    stars.type = "star"
    stars.name = name
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
        default=100,
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
    parser.add_argument("-Q", "--Qvir", type=float, default=0.5, help="total mass")
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

    if args.mmin > 0 | units.MSun:
        masses = new_salpeter_mass_distribution(args.nstars, args.mmin, args.mmax)
    else:
        masses = np.array(args.nstars) | units.MSun
        masses += args.mass
    converter = nbody_system.nbody_to_si(masses.sum(), args.radius)
    bodies = make_plummer_sphere(args.nstars, masses, args.name, converter)
    bodies.scale_to_standard(convert_nbody=converter, virial_ratio=args.Qvir)
    time = 0 | units.Myr
    if args.outfile is None:
        filename = "plummer.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
