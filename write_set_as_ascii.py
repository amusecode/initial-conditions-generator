import argparse

from amuse.io import read_set_from_file
from amuse.units import units


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-f", dest="filename", default=None, help="input filename"
    )
    parser.add_argument(
        "-F", dest="outfile", default="star.csv", help="output filename"
    )
    return parser


def main():
    args = new_argument_parser().parse_args()

    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies)

    with open(args.outfile, "w", encoding="utf-8") as f:
        print("# body key, name, type, mass (Msun), position (au), velocity (kms)", file=f)

        for bi in bodies:
            print(
                bi.key,
                bi.name,
                bi.type,
                bi.mass.value_in(units.MSun),
                bi.x.value_in(units.au),
                bi.y.value_in(units.au),
                bi.z.value_in(units.au),
                bi.vx.value_in(units.kms),
                bi.vy.value_in(units.kms),
                bi.vz.value_in(units.kms),
                sep=" ",
                file=f,
            )


if __name__ == "__main__":
    main()
