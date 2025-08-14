from amuse.lab import *
import sys
import numpy
from amuse.ext.solarsystem import new_solar_system
from amuse.ext.orbital_elements import (
    new_binary_from_orbital_elements,
    orbital_elements_from_binary,
)


def Hill_radius(Mstar, a, Mplanet):
    return a * (Mplanet / (3.0 * Mstar)) ** (1.0 / 3.0)


def determine_orbital_elements(bodies):
    from amuse.ext.orbital_elements import new_binary_from_orbital_elements

    star = bodies[bodies.type == "star"][0]
    planets = bodies[bodies.type == "planet"]
    for pi in planets:
        orbital_elements = orbital_elements_from_binary(star + pi)
        pi.semimajor_axis = orbital_elements[2]
        pi.eccentricity = orbital_elements[3]
        pi.inclination = orbital_elements[5]


def make_solar_system(iplanet_min, iplanet_max):
    ss = new_solar_system()
    ss.type = "planet"
    ss[0].type = "star"
    solar_system = ss
    """
    star = ss[ss.name=="SUN"]
    star.type = "star"
    star.id = 0
    planets = ss[iplanet_min:iplanet_max]
    planets.type = "planet"
    for i in range(len(planets)):
        planets[i].id = i
    solar_system = Particles(0)
    solar_system.add_particles(star)
    solar_system.add_particles(planets)
    solar_system.move_to_center()
    """

    determine_orbital_elements(solar_system)
    return solar_system


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option(
        "-f", dest="filename", default=None, help="input filename [%default]"
    )
    result.add_option(
        "-F",
        dest="outputfilename",
        default="star.amuse",
        help="outpute filename [%default]",
    )
    result.add_option(
        "--iplanet",
        dest="iplanet_min",
        type="int",
        default=1,
        help="minimal first planet (5=Jupiter) [%default]",
    )
    result.add_option(
        "--jplanet",
        dest="iplanet_max",
        type="int",
        default=9,
        help="maximal first planet (9=Neptune) [%default]",
    )
    result.add_option(
        "--seed",
        dest="seed",
        type="int",
        default=-1,
        help="random number seed [%default]",
    )
    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    if o.seed > 0:
        numpy.random.seed(o.seed)
    else:
        print("random number seed from clock.")

    if o.filename:
        bodies = read_set_from_file(o.filename, "hdf5")
    else:
        bodies = Particles(0)
    new_bodies = make_solar_system(o.iplanet_min, o.iplanet_max)
    bodies.add_particles(new_bodies)
    print(bodies)

    time = 0 | units.Myr
    write_set_to_file(
        bodies,
        o.outputfilename,
        "amuse",
        timestamp=time,
        append_to_file=False,
        version="2.0",
    )
