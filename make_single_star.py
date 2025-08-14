import numpy
import csv
import argparse

from amuse.units import units
from amuse.datamodel import Particle, Particles
from amuse.community.seba import Seba
from amuse.io import read_set_from_file, write_set_to_file


def read_Gaia_database(filename, index=0, age=0 | units.Myr):
    star = Particle()
    with open(filename, "r") as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            if i == index:
                star.name = row[1]
                star.type = "star"
                star.mass = float(row[2]) | units.MSun
                star.position = (
                    float(row[4]),
                    float(row[5]),
                    float(row[6]),
                ) | units.parsec
                star.velocity = (
                    float(row[7]),
                    float(row[8]),
                    float(row[9]),
                ) | units.kms
                break
            i += 1
    s = Seba()
    s.particles.add_particle(star)
    s.evolve_model(age)
    star.mass = s.particles[0].mass
    star.radius = s.particles[0].radius
    s.stop()
    print("star=", star)
    return star


def read_particle(filename, index=0):
    stars = read_set_from_file(filename, "hdf5")
    star = Particles(0)
    if index <= len(stars):
        star.add_particle(stars[index])
    return star


def Hill_radius(Mstar, a, Mplanet):
    return a * (Mplanet / (3.0 * Mstar)) ** (1.0 / 3.0)


def make_single_star(mass_star, radius, age, name):
    star = Particle()
    star.ZAMS_mass = mass_star
    star.mass = mass_star
    star.name = name
    star.type = "star"
    if radius < 0 | units.RSun:
        s = Seba()
        s.particles.add_particle(star)
        s.evolve_model(age)
        star.mass = s.particles[0].mass
        star.radius = s.particles[0].radius
        s.stop()
    else:
        star.radius = radius
    star.position = (0, 0, 0) | units.au
    star.velocity = (0, 0, 0) | units.kms
    return star


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default="star.amuse", help="output filename")
    parser.add_argument(
        "-M",
        "--mass_star",
        type=units.MSun,
        default=1 | units.MSun,
        help="stellar mass",
    )
    parser.add_argument(
        "-R",
        "--radius_star",
        type=units.RSun,
        default=-1.0 | units.RSun,
        help="stellar radius",
    )
    parser.add_argument(
        "-t",
        "--age",
        type=units.Myr,
        default=0.0 | units.Myr,
        help="stellar age",
    )
    parser.add_argument(
        "--gaia_database_file",
        default="gaia_nearby_stars.amuse",
        help="Gaia database dir and filename",
    )
    parser.add_argument(
        "--gaia_database_index",
        type=int,
        default=-1,
        help="Gaia database index",
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
        numpy.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    bodies = Particles(0)
    if args.gaia_database_index > 0:
        if ".amuse" in args.gaia_database_file:
            new_bodies = read_particle(
                args.gaia_database_file, args.gaia_database_index
            )
            print(new_bodies)
        elif ".csv" in args.gaia_database_file:
            new_bodies = read_Gaia_database(
                args.gaia_database_file, args.gaia_database_index
            )
        else:
            print("Unrecognized input filename: ", args.gaia_database_file)

    else:
        new_bodies = make_single_star(
            args.mass_star, args.radius_star, args.age, args.name
        )
    bodies.add_particle(new_bodies)

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "star.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies,
        filename,
        "amuse",
        timestamp=time,
        overwrite_file=True,
        append_to_file=False,
        version="2.0",
    )


if __name__ == "__main__":
    main()
