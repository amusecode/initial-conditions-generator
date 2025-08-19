"""Plot the particles in an AMUSE file.

This script plots all particles in one or more files, and highlights the stars,
planets and debris particles (they need to have the "type" property set to
"star", "planet" or "debris").
The user can optionally specify the length unit for the axes.
"""

import argparse

import matplotlib.pyplot as plt

from amuse.units import units
from amuse.io import write_set_to_file, read_set_from_file


def plot_bodies(bodies, length_unit=units.au):
    figure = plt.figure()
    ax = figure.add_subplot(111)
    stars = bodies[bodies.type == "star"]
    planets = bodies[bodies.type == "planet"]
    debris = bodies[bodies.type == "debris"]
    ax.scatter(
        bodies.x.value_in(length_unit),
        bodies.y.value_in(length_unit),
        s=0.1,
        c="k",
    )
    ax.scatter(
        stars.x.value_in(length_unit),
        stars.y.value_in(length_unit),
        s=30,
        c="y",
    )
    ax.scatter(
        planets.x.value_in(length_unit),
        planets.y.value_in(length_unit),
        s=10,
        c="b",
    )
    ax.scatter(
        debris.x.value_in(length_unit),
        debris.y.value_in(length_unit),
        s=1,
        c="k",
    )
    ax.set_xlabel(f"x [{length_unit}]")
    ax.set_ylabel(f"y [{length_unit}]")
    return figure


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "-f", "--filename", nargs="+", default=None, help="input filename(s)"
    )
    parser.add_argument("filenames", nargs="*")
    parser.add_argument(
        "-o", "--output", nargs="+", default=None, help="output filename"
    )
    parser.add_argument("-l", "--length_unit", default="au", help="length unit")
    return parser


def main():
    args = new_argument_parser().parse_args()
    if args.filename is not None:
        filenames = args.filename
    elif args.filenames:
        filenames = args.filenames
    else:
        raise ValueError("Must specify a filename to read")

    if args.output is not None:
        if len(args.output) != len(filenames):
            raise ValueError(
                "Number of output filenames must match number of input filenames"
            )

    try:
        length_unit = getattr(units, args.length_unit.lower())
    except AttributeError:
        raise ValueError("Unknown length unit %s" % args.length_unit)
    for i, filename in enumerate(filenames):
        print("Reading %s" % filename)
        bodies = read_set_from_file(filename, close_file=True)
        figure = plot_bodies(bodies, length_unit=length_unit)
        if args.output is not None:
            plt.savefig(args.output[i])
        else:
            plt.show()
        figure.clear()


if __name__ == "__main__":
    main()
