import argparse

import matplotlib.pyplot as plt

from amuse.units import units
from amuse.io import write_set_to_file, read_set_from_file


def plot_bodies(bodies):

    stars = bodies[bodies.type == "star"]
    planets = bodies[bodies.type == "planet"]
    debris = bodies[bodies.type == "debris"]
    plt.scatter(bodies.x.value_in(units.au), bodies.y.value_in(units.au), s=0.1, c="k")
    plt.scatter(stars.x.value_in(units.au), stars.y.value_in(units.au), s=30, c="y")
    plt.scatter(planets.x.value_in(units.au), planets.y.value_in(units.au), s=10, c="b")
    plt.scatter(debris.x.value_in(units.au), debris.y.value_in(units.au), s=1, c="k")
    plt.xlabel("x [au]")
    plt.ylabel("y [au]")
    plt.show()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    return parser


def main():
    args = new_argument_parser().parse_args()

    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies)
    plot_bodies(bodies)


if __name__ == "__main__":
    main()
