from amuse.lab import *
import sys
import numpy

from amuse.lab import new_kroupa_mass_distribution
from amuse.lab import new_plummer_sphere
from amuse.community.fractalcluster.interface import new_fractal_cluster_model


def make_fractal_sphere(nstars, masses, name, Fd, converter):
    stars = new_fractal_cluster_model(
        nstars, convert_nbody=converter, fractal_dimension=Fd
    )
    stars.type = "star"
    stars.name = name
    stars.mass = masses
    return stars


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
        "--nstars",
        dest="nstars",
        type="int",
        default=100,
        help="number of stars [%default]",
    )
    result.add_option(
        "--mass",
        unit=units.MSun,
        dest="mass",
        type="float",
        default=1 | units.MSun,
        help="total cluster mass [%default]",
    )
    result.add_option(
        "--mmin",
        unit=units.MSun,
        dest="mmin",
        type="float",
        default=0.08 | units.MSun,
        help="minimum stellar mass [%default]",
    )
    result.add_option(
        "--mmax",
        unit=units.MSun,
        dest="mmax",
        type="float",
        default=100 | units.MSun,
        help="maximum stellar mass [%default]",
    )
    result.add_option(
        "-Q", dest="Qvir", type="float", default=0.5, help="total mass [%default]"
    )
    result.add_option(
        "--Fd",
        dest="Fd",
        type="float",
        default=1.6,
        help="Fractal dimension [%default]",
    )
    result.add_option(
        "--radius",
        unit=units.parsec,
        dest="radius",
        type="float",
        default=1.0 | units.parsec,
        help="cluster radius [%default]",
    )
    result.add_option(
        "--nsuns",
        dest="nsun_like_stars",
        type="int",
        default=-1,
        help="number of Sun-like stars [%default]",
    )
    result.add_option(
        "--name", dest="name", default="star", help="stellar name [%default]"
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

    if o.mmin > 0 | units.MSun:
        masses = new_kroupa_mass_distribution(o.nstars, o.mmin, o.mmax)
    else:
        masses = numpy.array(o.nstars) | units.MSun
        masses += o.mass

    converter = nbody_system.nbody_to_si(masses.sum(), o.radius)
    bodies = make_fractal_sphere(o.nstars, masses, o.name, o.Fd, converter)

    if o.nsun_like_stars > 0:
        bodies = bodies.sorted_by_attributes("mass")[::-1]
        print(bodies.mass.in_(units.MSun))
        imlo = 0
        imlo = 0
        for bi in range(len(bodies)):
            if bodies[bi].mass < 1 | units.MSun:
                im = bi
                break
        bi = 1
        suns = bodies[im - bi : im + bi]
        while len(suns) < o.nsun_like_stars:
            bi += 1
            suns = bodies[im - bi : im + bi]

        print(suns.mass.in_(units.MSun))
        suns.mass = 1 | units.MSun
        suns.name = "Sun"

    bodies.move_to_center()
    bodies.scale_to_standard(convert_nbody=converter, virial_ratio=o.Qvir)
    index = 0
    time = 0 | units.Myr
    if o.outfile == None:
        filename = "fractal.amuse".format(index)
    else:
        filename = o.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )
