import argparse

from amuse.io import read_set_from_file, write_set_to_file


def move_bodies_to_stellar_position(bodies):
    star = bodies[bodies.type == "star"]
    # print "pos=", star.position.in_(units.parsec)
    com = star.center_of_mass()
    comv = star.center_of_mass_velocity()
    bodies.move_to_center()
    bodies.position += com
    bodies.velocity += comv
    # print "pos=", bodies.position.in_(units.parsec)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f", "--filename", default="star.amuse", help="input filename"
    )
    parser.add_argument(
        "-F", "--outfile", default=None, help="output filename"
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    bodies = read_set_from_file(args.filename, close_file=True)

    move_bodies_to_stellar_position(bodies)
    # print bodies.velocity

    if args.outfile is None:
        filename = "moved_bodies.amuse"
    else:
        filename = args.outfile
    write_set_to_file(bodies, filename, "amuse", append_to_file=False, version="2.0")


if __name__ == "__main__":
    main()
