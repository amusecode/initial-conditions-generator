import argparse

import numpy as np

from amuse.units import units
from amuse.datamodel import Particles
from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.io import write_set_to_file


def ZAMS_radius(mass):
    log_mass = np.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun)) ** 2
    alpha = 0.08353 + 0.0565 * log_mass
    beta = 0.01291 + 0.2226 * log_mass
    gamma = 0.1151 + 0.06267 * log_mass
    r_zams = (
        pow(mass.value_in(units.MSun), 1.25)
        * (0.1148 + 0.8604 * mass_sq)
        / (0.04651 + mass_sq)
    )

    return r_zams | units.RSun


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-F", dest="outfile", default="star.amuse", help="output filename"
    )
    parser.add_argument(
        "-N", dest="n_star", type="int", default=10, help="number of stars"
    )
    parser.add_argument(
        "-m",
        unit=units.MSun,
        dest="m_min",
        type="float",
        default=0.1 | units.MSun,
        help="number of stars",
    )
    parser.add_argument(
        "--Salpeter",
        action="store_true",
        default="False",
        help="choose between Salpeter or Kroupa",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        type="int",
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

    if args.Salpeter:
        mass = new_kroupa_mass_distribution(args.n_star, mass_min=args.m_min)
    else:
        mass = new_salpeter_mass_distribution(args.n_star, mass_min=args.m_min)
    stars = Particles(len(mass))
    stars.mass = mass
    stars.position = (0, 0, 0) | units.pc
    stars.velocity = (0, 0, 0) | units.kms
    stars.type = "star"
    stars.name = "star"
    stars.radius = ZAMS_radius(stars.mass)

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "stars.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        stars, filename, "amuse", timestamp=time, overwrite_file=True, version="2.0"
    )


if __name__ == "__main__":
    main()
