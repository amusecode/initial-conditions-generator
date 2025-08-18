import argparse

import numpy as np

from amuse.units import units, nbody_system
from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.plummer import new_plummer_sphere
from amuse.io import write_set_to_file
from amuse.ic.fractalcluster import new_fractal_cluster_model


def make_fractal_sphere(nstars, masses, name, Fd, converter):
    stars = new_fractal_cluster_model(
        nstars, convert_nbody=converter, fractal_dimension=Fd
    )
    stars.type = "star"
    stars.name = name
    stars.mass = masses
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
        default=0.08 | units.MSun,
        help="minimum stellar mass",
    )
    parser.add_argument(
        "--mmax",
        type=units.MSun,
        default=100 | units.MSun,
        help="maximum stellar mass",
    )
    parser.add_argument("-Q", "--Qvir", type=float, default=0.5, help="virial ratio")
    parser.add_argument(
        "--Fd",
        type=float,
        default=1.6,
        help="Fractal dimension",
    )
    parser.add_argument(
        "--radius",
        type=units.parsec,
        default=1.0 | units.parsec,
        help="cluster radius",
    )
    parser.add_argument(
        "--nsuns",
        "--nsun_like_stars",
        type=int,
        default=-1,
        help="number of Sun-like stars",
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
        masses = new_kroupa_mass_distribution(args.nstars, args.mmin, args.mmax)
    else:
        masses = np.array(args.nstars) | units.MSun
        masses += args.mass

    converter = nbody_system.nbody_to_si(masses.sum(), args.radius)
    bodies = make_fractal_sphere(args.nstars, masses, args.name, args.Fd, converter)

    if args.nsun_like_stars > 0:
        bodies = bodies.sorted_by_attributes("mass")[::-1]
        print(bodies.mass.in_(units.MSun))
        for bi in range(len(bodies)):
            if bodies[bi].mass < 1 | units.MSun:
                im = bi
                break
        bi = 1
        suns = bodies[im - bi : im + bi]
        while len(suns) < args.nsun_like_stars:
            bi += 1
            suns = bodies[im - bi : im + bi]

        print(suns.mass.in_(units.MSun))
        suns.mass = 1 | units.MSun
        suns.name = "Sun"

    bodies.move_to_center()
    bodies.scale_to_standard(convert_nbody=converter, virial_ratio=args.Qvir)
    time = 0 | units.Myr
    if args.outfile is None:
        filename = "fractal.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
