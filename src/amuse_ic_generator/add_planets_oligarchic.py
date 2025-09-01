"""Add planets to a subsection of stars in an AMUSE file.

Takes an AMUSE file as input, selects a subset of stars (specified by name and mass),
and adds a random number of planets to each star, following an oligarchic distribution.
The radius and mass of the disk can be specified, as well as the number of planets and their mass relative to the star mass.
"""

import argparse

import numpy as np

from amuse.units import units
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.io import write_set_to_file, read_set_from_file

from make_planets_oligarch import make_planets_oligarch
from move_bodies_to_stellar_position import move_bodies_to_stellar_position


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("-f", "--filename", default="star.amuse", help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument("--name", default=None, help="Add planets to named star")
    parser.add_argument(
        "--rmin_disk",
        type=units.au,
        default=1.0 | units.au,
        help="inner disk radius",
    )
    parser.add_argument(
        "--rmax_disk",
        type=units.au,
        default=100.0 | units.au,
        help="outer disk radius",
    )
    parser.add_argument(
        "--relative_diskmass",
        type=float,
        default=0.01,
        help="relative disk mass",
    )
    parser.add_argument(
        "--Nplanets",
        type=int,
        default=4,
        help="maximum number of planets",
    )
    parser.add_argument(
        "--fplanets",
        type=float,
        default=0.5,
        help="fraction of stars with planets",
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
    parser.add_argument(
        "--mscale",
        type=float,
        default=1,
        help="scaling of the planet mass",
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
    accurancy = 0.0001
    planet_density = 3 | units.g / units.cm**3

    if args.seed > 0:
        np.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies)
    stars = bodies[bodies.type == "star"]
    stars = stars[stars.name == args.name]
    stars = stars[stars.mass >= args.mmin]
    stars = stars[stars.mass <= args.mmax]
    if args.fplanets <= 1:
        nplanets = int(args.fplanets * len(stars))
    else:
        nplanets = min(len(stars), int(args.fplanets))
    stars_with_planets = stars.random_sample(nplanets)
    print("Number of stars with planets=", len(stars_with_planets))
    for star in stars_with_planets:
        disk_mass = args.relative_diskmass * star.mass
        planets = make_planets_oligarch(
            star.mass, star.radius, args.rmin_disk, args.rmax_disk, disk_mass
        )
        planets = planets[planets.type == "planet"]
        while len(planets) > args.Nplanets:
            mmin = planets.mass.min()
            lm_planet = planets[planets.mass <= mmin]
            planets.remove_particle(lm_planet)

        planets.type = "planet"
        planets.name = "planet"
        planets.host = star.key
        print(planets.position.in_(units.parsec))
        planets.position += star.position
        planets.velocity += star.velocity
        planets.mass *= args.mscale
        bodies.add_particles(planets)
    move_bodies_to_stellar_position(bodies)

    planets = bodies[bodies.type == "planet"]
    for i, planet in enumerate(planets):
        print(
            f"Planet: {i: 3d}, "
            f"mass: {planet.mass.value_in(units.MEarth): 8.3f}  MEarth, "
            f"a: {planet.semimajor_axis.value_in(units.au): 8.2f} au"
        )

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "added_oligarch_planets.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies,
        filename,
        "amuse",
        timestamp=time,
        append_to_file=False,
        overwrite_file=True,
        version="2.0",
    )


if __name__ == "__main__":
    main()
