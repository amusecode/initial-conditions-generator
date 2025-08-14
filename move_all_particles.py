from amuse.lab import *

# import sys
# import numpy


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option(
        "-f", dest="filename", default=None, help="input filename [%default]"
    )
    result.add_option(
        "-F", dest="outfile", default=None, help="output filename [%default]"
    )
    result.add_option(
        "--pos",
        dest="position",
        type="float",
        nargs=3,
        default=[-8340.0, 0.0, 27.0],  # parsec
        help="move com position to [in pc]",
    )
    result.add_option(
        "--vel",
        dest="velocity",
        type="float",
        nargs=3,
        default=[11.1, 240.0, 7.25],
        help="move com-velocity to [in km/s]",
    )
    result.add_option(
        "-t",
        unit=units.Myr,
        dest="timestamp",
        type="float",
        default=0 | units.Myr,
        help="set timestamp [%default]",
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

    bodies = read_set_from_file(o.filename, "hdf5", close_file=True)

    bodies.move_to_center()
    bodies.position += o.position | units.parsec
    bodies.velocity += o.velocity | units.kms

    print(bodies)
    time = o.timestamp
    if o.outfile == None:
        filename = "moved_snapshot.amuse"
    else:
        filename = o.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )
