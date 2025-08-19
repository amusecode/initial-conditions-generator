"""Add a debris disk to stars read from an AMUSE file.

This script reads an AMUSE file containing stellar particles and adds debris
disks to the stars.
"""

import argparse

import numpy as np

from amuse.units import units, nbody_system
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.couple import bridge
from amuse.io import write_set_to_file, read_set_from_file


def add_debris_disk(
    star, ndisk_per_MSun, fdisk, masteroid, Rmin, Rmax, alpha=1.5, Q=1.0
):

    Ndisk = int(ndisk_per_MSun * star.mass.value_in(units.MSun))
    print("Run with Ndisk = ", Ndisk)

    converter = nbody_system.nbody_to_si(star.mass, Rmin)
    disk = ProtoPlanetaryDisk(
        Ndisk,
        convert_nbody=converter,
        densitypower=alpha,
        Rmin=Rmin / Rmin,
        Rmax=Rmax / Rmin,
        q_out=Q,
        discfraction=fdisk,
    ).result

    disk.mass = masteroid
    disk.name = "asteroid"
    disk.type = "debris"
    disk.host = star.key
    disk.position += star.position
    disk.velocity += star.velocity
    return disk


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument("--name", default="star", help="name of star type")
    parser.add_argument("--dname", default="asteroid", help="name of disk type")
    parser.add_argument("--dtype", default="debris", help="type of disk")
    parser.add_argument(
        "--fdisk",
        type=float,
        default=-1,
        help="disk mass fraction",
    )
    parser.add_argument(
        "--Mdisk",
        type=units.MSun,
        default=0.01 | units.MSun,
        help="disk mass",
    )
    parser.add_argument(
        "--mdisk",
        type=units.kg,
        default=0.0 | units.kg,
        help="mass of individual disk particles",
    )
    parser.add_argument(
        "--ndisk",
        type=int,
        default=100,
        help="number of disk particles per MSun",
    )
    parser.add_argument(
        "--rmin",
        type=units.au,
        default=10 | units.au,
        help="minimal disk radius",
    )
    parser.add_argument(
        "--rmax",
        type=units.au,
        default=100 | units.au,
        help="maximal disk radius",
    )
    parser.add_argument("-Q", type=float, default=1.0, help="Toomre Q parameter")
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=1.5,
        help="Disk density profile",
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

    if args.seed > 0:
        np.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(args.filename, close_file=True)
    stars = bodies[bodies.name == args.name]
    if len(stars) == 0:
        stars = bodies[bodies.type.lower() in ["star", "sun"]]

    for star in stars:
        if args.fdisk < 0:
            fdisk = args.Mdisk / star.mass
        else:
            fdisk = args.fdisk
        debris = add_debris_disk(
            star,
            args.ndisk,
            fdisk,
            args.mdisk,
            args.rmin,
            args.rmax,
            args.alpha,
            args.Q,
        )
        debris.name = args.dname
        debris.type = args.dtype
        bodies.add_particles(debris)
        bodies.move_to_center()

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "added_debris_disk.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
