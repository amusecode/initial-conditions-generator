import numpy
import argparse
import logging

from amuse.units import units, constants
from amuse.units.trigo import sin
from amuse.datamodel import Particle
from amuse.io import read_set_from_file, write_set_to_file
from amuse.ext.orbital_elements import new_binary_from_orbital_elements


def true_anomaly_from_mean_anomaly(Ma, e):
    Ta = (
        Ma
        + (2 * e - e**3 / 4.0) * sin(Ma)
        + (5 * e**2 / 4.0) * sin(2 * Ma)
        + (13.0 / 12.0) * e**3 * sin(3 * Ma)
    )
    return Ta


def semimajor_axis_to_orbital_period(a, Mtot):
    return 2 * numpy.pi * (a**3 / (constants.G * Mtot)).sqrt()


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


def add_secondaries(
    parent_stars,
    companion_name,
    masses,
    semimajor_axis,
    eccentricity,
    inclination,
    mean_anomaly,
    LoAn,
    Aop,
    ctype="star",
):

    print("m=", masses)
    for i in range(len(parent_stars)):
        bi = parent_stars[i]
        binary_particle = Particle()
        binary_particle.position = bi.position
        binary_particle.velocity = bi.velocity
        binary_particle.type = ctype
        binary_particle.name = bi.name
        binary_particle.host = None

        binary_particle.type = "center_of_mass"

        mp = bi.mass
        ms = masses[i]

        Ta = true_anomaly_from_mean_anomaly(numpy.deg2rad(mean_anomaly), eccentricity)
        nb = new_binary_from_orbital_elements(
            mp,
            ms,
            semimajor_axis,
            eccentricity,
            Ta,
            inclination,
            LoAn,
            Aop,
            G=constants.G,
        )

        nb.position += binary_particle.position
        nb.velocity += binary_particle.velocity
        nb[0].type = bi.type
        # nb[0].host = binary_particle
        nb[0].name = bi.name
        nb[1].type = ctype
        # nb[1].host = binary_particle
        # nb[1].host = nb[0].name
        nb[1].host = bi.key

        nb[1].name = companion_name
        nb[1].radius = ZAMS_radius(nb[1].mass)
    return nb

    #    binary_particle.child1 = nb[0]
    #    binary_particle.child2 = nb[1]
    #    binary_particle.semi_major_axis = semimajor_axis
    #    binary_particle.eccentricity = eccentricity
    # return binary_particle


def calculate_orbital_elementss(bi, converter):
    kep = new_kepler(converter)
    comp1 = bi.child1
    comp2 = bi.child2
    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a, e = kep.get_elements()
    kep.stop()
    return a, e


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    result.add_argument(
        "-f", "--filename", default="Plummer.amuse", help="input filename"
    )
    result.add_argument("-F", "--outfile", default=None, help="output filename")
    result.add_argument(
        "--name",
        default="ABAuriga",
        help="name of the star that recieves companion",
    )
    result.add_argument(
        "--cname",
        "--companion_name",
        default="ABAur_b",
        help="name of the companion",
    )
    result.add_argument("--ctype", default="star", help="type of the companion")
    result.add_argument(
        "-m",
        "--mass",
        type=units.MSun,
        default=0 | units.MSun,
        help="binary companion mass",
    )
    result.add_argument(
        "-q",
        "--q_companion",
        type=float,
        default=-1,
        help="binary companion mass ratio",
    )
    result.add_argument(
        "-a",
        "--semimajor_axis",
        type=units.au,
        default=1.0 | units.au,
        help="binary separation",
    )
    result.add_argument(
        "-e",
        dest="eccentricity",
        type=float,
        default=0.0,
        help="binary eccenticity",
    )
    result.add_argument(
        "-i",
        "--inclination",
        type=float,
        default=0.0,
        help="inclination",
    )
    result.add_argument(
        "--ma",
        "--mean_anomaly",
        type=float,
        default=0.0,
        help="mean anomaly",
    )
    result.add_argument(
        "--LoAn",
        type=float,
        default=0.0,
        help="Longitude of the ascending node",
    )
    result.add_argument(
        "--Aop",
        type=float,
        default=0.0,
        help="Argument of pericenter",
    )
    result.add_argument(
        "--seed",
        type=int,
        default=None,
        help="random number seed",
    )
    return result


if __name__ in ("__main__", "__plot__"):
    args = new_argument_parser().parse_args()
    outfile = args.outfile

    if outfile is None:
        outfile = "stars_with_companion.amuse"

    bodies = read_set_from_file(args.filename, close_file=True)
    selected_stars = bodies[bodies.name == args.name]

    if args.q_companion > 0:
        mass = selected_stars.mass * args.q_companion
    else:
        mass = numpy.zeros(len(selected_stars)) | units.MSun
        mass += args.mass
    stars = add_secondaries(
        selected_stars,
        args.companion_name,
        mass,
        args.semimajor_axis,
        args.eccentricity,
        args.inclination,
        args.mean_anomaly,
        args.LoAn,
        args.Aop,
        args.ctype,
    )
    print(stars[1])
    bodies.add_particle(stars[1].as_set())
    print(bodies)
    time = 0 | units.Myr
    write_set_to_file(
        bodies,
        outfile,
        "amuse",
        timestamp=time,
        append_to_file=False,
        overwrite_file=True,
        version="2",
    )
