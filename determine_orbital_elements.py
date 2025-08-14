from amuse.lab import *
from amuse.ext.orbital_elements import orbital_elements_from_binary


def calculate_orbital_elements_for_single_planet(star, planet):
    from amuse.ext.orbital_elements import orbital_elements_from_binary

    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(
        p, G=constants.G
    )
    return a, e, inc


def determine_orbits(parent, planet):
    for pi in planet:
        a, e, inc = calculate_orbital_elements_for_single_planet(parent, pi)
        pi.semimajor_axis = a
        pi.eccentricity = e
        pi.inclination = inc


def determine_orbital_elements(bodies):
    # print bodies
    star = bodies[bodies.type == "star"]
    print(star)
    planets = bodies - star
    if len(planets) > 0:
        determine_orbits(star, planets)


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option(
        "-f",
        dest="filename",
        default="starplanetplets_i0000.amuse",
        help="input filename [%default]",
    )
    result.add_option(
        "-F", dest="outfile", default=None, help="output filename [%default]"
    )
    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    bodies = read_set_from_file(o.filename, "hdf5", close_file=True)
    print(bodies)
    determine_orbital_elements(bodies)
    if o.outfile:
        write_set_to_file(bodies, o.outfile, "hdf5", append_to_file=False)

    from matplotlib import pyplot

    q = bodies[1:].semimajor_axis.value_in(units.au) * (1 - bodies[1:].eccentricity)
    # pyplot.scatter(q, bodies[1:].eccentricity)
    pyplot.scatter(bodies.x.value_in(units.au), bodies.y.value_in(units.au))
    pyplot.show()
