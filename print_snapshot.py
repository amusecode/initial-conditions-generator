import argparse

from amuse.io import read_set_from_file


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f", "--filename", default=None, help="input filename"
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies)


if __name__ == "__main__":
    main()
