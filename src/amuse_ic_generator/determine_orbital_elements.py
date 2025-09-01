import argparse

import matplotlib.pyplot as plt

from amuse.units import constants, units
from amuse.datamodel import Particles
from amuse.io import write_set_to_file, read_set_from_file
from amuse.ext.orbital_elements import orbital_elements_from_binary


def calculate_orbital_elements_for_single_planet(star, planet):
    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(
        p, G=constants.G
    )
    return a, e, inc


def determine_orbits(parent, planet):
    for pi in planet:
        a, e, inc = calculate_orbital_elements_for_single_planet(parent, pi)
        pi.semimajor_axis = a
        pi.eccentricity = e
        pi.inclination = inc


def determine_orbital_elements(bodies):
    # print bodies
    star = bodies[bodies.type == "star"]
    print(star)
    planets = bodies - star
    if len(planets) > 0:
        determine_orbits(star, planets)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--filename",
        default="starplanetplets_i0000.amuse",
        help="input filename",
    )
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    return parser


def main():
    args = new_argument_parser().parse_args()
    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies)
    determine_orbital_elements(bodies)
    if args.outfile:
        write_set_to_file(bodies, args.outfile, append_to_file=False)

    q = bodies[1:].semimajor_axis.value_in(units.au) * (1 - bodies[1:].eccentricity)
    # pyplot.scatter(q, bodies[1:].eccentricity)
    plt.scatter(bodies.x.value_in(units.au), bodies.y.value_in(units.au))
    plt.show()


if __name__ == "__main__":
    main()
