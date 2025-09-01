import argparse

import numpy as np

from amuse.units import units
from amuse.io import read_set_from_file
from amuse.ext.solarsystem import new_solar_system


def write_brutus_file(ss, time=0, cpu_time=0):
    mass = ss.mass.value_in(units.MSun)
    pos = ss.position.value_in(units.au)
    brutus_velocity_units = (2.0 * np.pi * 1 | units.au) / (1 | units.yr)
    vel = ss.velocity / brutus_velocity_units

    time = 0
    cpu_time = 0

    print(time, len(ss), cpu_time)
    for si in range(len(ss)):
        print(
            f"{mass[si]:.16f}",
            f"{pos[si][0]:.16f}",
            f"{pos[si][1]:.16f}",
            f"{pos[si][2]:.16f}",
            f"{vel[si][0]:.16f}",
            f"{vel[si][1]:.16f}",
            f"{vel[si][2]:.16f}",
        )


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--filename",
        default="solar_system.amuse",
        help="Output filename",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    ss = read_set_from_file(args.filename, "amuse", close_file=True)
    write_brutus_file(ss)


if __name__ == "__main__":
    main()
