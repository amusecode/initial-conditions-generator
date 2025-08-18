import argparse
import numpy

from amuse.units import units, nbody_system
from amuse.ic import new_salpeter_mass_distribution, new_plummer_model
from amuse.io import write_set_to_file


def ZAMS_radius(mass):
    log_mass = numpy.log10(mass.value_in(units.MSun))
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
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    result.add_argument("-f", "--filename", default=None, help="input filename")
    result.add_argument(
        "-F", "--outfile", default="Plummer.amuse", help="output filename"
    )
    result.add_argument(
        "-N", dest="--number_of_stars", type=int, default=10, help="number of stars"
    )
    result.add_argument(
        "-m",
        "--mass_min",
        type=units.MSun,
        default=0.1 | units.MSun,
        help="number of stars",
    )
    result.add_argument(
        "-R",
        "--radius_virial",
        type=units.parsec,
        default=1.0 | units.parsec,
        help="cluster Plummer radius",
    )
    result.add_argument(
        "--seed",
        type=int,
        default=-1,
        help="random number seed",
    )
    return result


def main():
    args = new_argument_parser().parse_args()
    if args.seed > 0:
        numpy.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    mass = new_salpeter_mass_distribution(args.number_of_stars, mass_min=args.mass_min)
    converter = nbody_system.nbody_to_si(mass.sum(), args.radius_virial)
    stars = new_plummer_model(len(mass), convert_nbody=converter)
    stars.type = "star"
    stars.name = "star"
    stars.radius = ZAMS_radius(stars.mass)

    time = 0 | units.Myr
    write_set_to_file(
        stars,
        args.outfile,
        "amuse",
        timestamp=time,
        append_to_file=False,
        version="2.0",
    )


if __name__ == "__main__":
    main()
