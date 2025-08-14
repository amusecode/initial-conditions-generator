import numpy
from amuse.lab import *
from amuse.ext.solarsystem import new_solar_system


def write_brutus_file(ss, time=0, cpu_time=0):
    mass = ss.mass.value_in(units.MSun)
    pos = ss.position.value_in(units.au)
    brutus_velocity_units = (2.0 * numpy.pi * 1 | units.au) / (1 | units.yr)
    vel = ss.velocity / brutus_velocity_units

    time = 0
    cpu_time = 0

    print(time, len(ss), cpu_time)
    for si in range(len(ss)):
        print(
            "{:.16f}".format(mass[si]),
            "{:.16f}".format(pos[si][0]),
            "{:.16f}".format(pos[si][1]),
            "{:.16f}".format(pos[si][2]),
            "{:.16f}".format(vel[si][0]),
            "{:.16f}".format(vel[si][1]),
            "{:.16f}".format(vel[si][2]),
        )


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option(
        "-f",
        dest="filename",
        default="solar_system.amuse",
        help="Output filename [%default]",
    )
    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()

    ss = read_set_from_file(o.filename, "amuse", close_file=True)
    write_brutus_file(ss)
