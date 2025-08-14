from amuse.lab import *
import sys
import numpy


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


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option(
        "-F", dest="outfile", default="star.amuse", help="output filename [%default]"
    )
    result.add_option(
        "-N", dest="n_star", type="int", default=10, help="number of stars [%default]"
    )
    result.add_option(
        "-m",
        unit=units.MSun,
        dest="m_min",
        type="float",
        default=0.1 | units.MSun,
        help="number of stars [%default]",
    )
    result.add_option(
        "--Salpeter",
        action="store_true",
        default="False",
        help="choose between Salpeter or Kroupa [%default]",
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

    if o.Salpeter:
        mass = new_kroupa_mass_distribution(o.n_star, mass_min=o.m_min)
    else:
        mass = new_salpeter_mass_distribution(o.n_star, mass_min=o.m_min)
    stars = Particles(len(mass))
    stars.mass = mass
    stars.position = (0, 0, 0) | units.pc
    stars.velocity = (0, 0, 0) | units.kms
    stars.type = "star"
    stars.name = "star"
    stars.radius = ZAMS_radius(stars.mass)

    time = 0 | units.Myr
    if o.outfile == None:
        filename = "stars.amuse".format(index)
    else:
        filename = o.outfile
    write_set_to_file(
        stars, filename, "amuse", timestamp=time, overwrite_file=True, version="2.0"
    )
