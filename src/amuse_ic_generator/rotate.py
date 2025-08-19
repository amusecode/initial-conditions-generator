"""Rotate (sub)sets of particles from an AMUSE file.
"""
import sys
import argparse
import numpy as np

from amuse.datamodel.rotation import rotate
from amuse.units.trigo import pi, arccos
from amuse.units import nbody_system, units
from amuse.io import write_set_to_file, read_set_from_file


def random_Euler_angles():
    phi = 2 * pi * np.random.random()
    theta = arccos(1 - 2 * np.random.random())
    chi = 2 * pi * np.random.random()
    return phi, theta, chi


def test_rotate(phi, theta, chi):
    from amuse.ext.protodisk import ProtoPlanetaryDisk

    converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.au)
    disk = ProtoPlanetaryDisk(
        1000, convert_nbody=converter, Rmin=1, Rmax=100, q_out=1.0, discfraction=0.01
    ).result

    rotate(disk, phi, theta, chi)

    from matplotlib import pyplot

    pyplot.scatter(disk.x.value_in(units.au), disk.y.value_in(units.au))
    pyplot.xlabel("x-axis [au]")
    pyplot.ylabel("y-axis [au]")
    pyplot.show()
    pyplot.scatter(disk.x.value_in(units.au), disk.z.value_in(units.au))
    pyplot.xlabel("x-axis [au]")
    pyplot.ylabel("z-axis [au]")
    pyplot.show()


def rotate_bodies_isotropically(particles):
    phi, theta, chi = random_Euler_angles()
    phi = np.rad2deg(phi)
    theta = np.rad2deg(theta)
    chi = np.rad2deg(chi)
    rotate(particles, phi, theta, chi)
    return particles


def rotate_minor_bodies_around_star(star, bodies):
    disk = bodies[bodies.host == star.key]
    disk.position -= star.position
    disk.velocity -= star.velocity
    rotate_bodies_isotropically(disk)
    disk.position += star.position
    disk.velocity += star.velocity


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "-p",
        "--phi",
        type=units.deg,
        default=0.0 | units.deg,
        help="rotate under x-axis",
    )
    parser.add_argument(
        "-t",
        "--theta",
        type=units.deg,
        default=-62.5 | units.deg,
        help="rotate under y-axis",
    )
    parser.add_argument(
        "--type", "--rotate_type", default="", help="rotate particle type"
    )
    parser.add_argument(
        "-c",
        "--chi",
        type=units.deg,
        default=0.0 | units.deg,
        help="rotate under z-axis",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    if args.filename is None:
        test_rotate(args.phi, args.theta, args.chi)
        sys.exit()

    bodies = read_set_from_file(args.filename, close_file=True)
    if len(args.rotate_type) > 0:
        stars = bodies[bodies.type == args.rotate_type]
        for star in stars:
            rotate_minor_bodies_around_star(star, bodies)
    else:
        com = bodies.center_of_mass()
        comv = bodies.center_of_mass_velocity()
        bodies.move_to_center()

        if args.phi is None:
            phi, theta, chi = random_Euler_angles()
            print(f"rotate randomly to (phi, theta, chi): {phi} {theta} {chi}")
        else:
            phi = args.phi
            theta = args.theta
            chi = args.chi

        rotate(bodies, phi, theta, chi)  # takes angles in degrees
        bodies.position += com
        bodies.velocity += comv

    if args.outfile is None:
        outfile = "particles_rotated.amuse"
    else:
        outfile = args.outfile
    time = 0 | units.Myr
    write_set_to_file(
        bodies, outfile, "amuse", timestamp=time, overwrite_file=True, version="2.0"
    )


if __name__ == "__main__":
    main()
