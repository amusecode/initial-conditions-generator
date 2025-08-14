from amuse.lab import *
import numpy
import logging


def semimajor_axis_to_orbital_period(a, Mtot):
    return 2 * numpy.pi * (a**3 / (constants.G * Mtot)).sqrt()


def new_kepler(converter):
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.set_longitudinal_unit_vector(1.0, 0.0, 0.0)
    kepler.set_transverse_unit_vector(0.0, 1.0, 0)
    return kepler


def new_binary_orbit(mass1, mass2, semi_major_axis, eccentricity=0, keyoffset=1):
    total_mass = mass1 + mass2
    mass_fraction_particle_1 = mass1 / (total_mass)

    binary = Particles(2)
    binary[0].mass = mass1
    binary[1].mass = mass2

    mu = constants.G * total_mass

    velocity_perihelion = numpy.sqrt(
        mu / semi_major_axis * ((1.0 + eccentricity) / (1.0 - eccentricity))
    )
    radius_perihelion = semi_major_axis * (1.0 - eccentricity)

    binary[0].position = (
        (1.0 - mass_fraction_particle_1) * radius_perihelion * [1.0, 0.0, 0.0]
    )
    binary[1].position = -(
        mass_fraction_particle_1 * radius_perihelion * [1.0, 0.0, 0.0]
    )

    binary[0].velocity = (
        (1.0 - mass_fraction_particle_1) * velocity_perihelion * [0.0, 1.0, 0.0]
    )
    binary[1].velocity = -(
        mass_fraction_particle_1 * velocity_perihelion * [0.0, 1.0, 0.0]
    )

    return binary


def random_semimajor_axis(amin, amax):
    lamin = numpy.log10(amin.value_in(units.au))
    lamax = numpy.log10(amax.value_in(units.au))
    rnd_min = lamin
    rnd_max = lamax
    rnd = numpy.random.uniform(rnd_min, rnd_max, 1)
    a = 10**rnd
    return a | units.au


# see Eggleton 2006 Equation 1.6.3 (2006epbm.book.....E)
def random_semimajor_axis_PPE(Mprim, Msec, amin, amax):

    Pmax = semimajor_axis_to_orbital_period(amax, Mprim + Msec).value_in(units.day)
    Pmin = semimajor_axis_to_orbital_period(amin, Mprim + Msec).value_in(units.day)
    mpf = (Mprim.value_in(units.MSun) ** 2.5) / 5.0e4
    rnd_max = (Pmax * mpf) ** (1.0 / 3.3) / (1 + (Pmin * mpf) ** (1.0 / 3.3))
    rnd_min = (Pmin * mpf) ** (1.0 / 3.3) / (1 + (Pmax * mpf) ** (1.0 / 3.3))
    rnd_max = min(rnd_max, 1)
    rnd = numpy.random.uniform(rnd_min, rnd_max, 1)
    Porb = ((rnd / (1.0 - rnd)) ** 3.3) / mpf | units.day
    print(Pmin, Pmax, Porb)
    xx

    Mtot = Mprim + Msec
    a = ((constants.G * Mtot) * (Porb / (2 * numpy.pi)) ** 2) ** (1.0 / 3.0)
    return a


def ZAMS_radius(mass):
    log_mass = numpy.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun)) ** 2
    alpha = 0.08353 + 0.0565 * log_mass
    beta = 0.01291 + 0.2226 * log_mass
    gamma = 0.1151 + 0.06267 * log_mass
    r_zams = (
        pow(mass.value_in(units.MSun), 1.25)
        * (0.1148 + 0.8604 * mass_sq)
        / (0.04651 + mass_sq)
    )

    return r_zams | units.RSun


def make_secondaries(center_of_masses, companion_mass, semimajor_axis, eccentricit):

    resulting_binaries = Particles()
    singles_in_binaries = Particles()
    binaries = center_of_masses.random_sample(Nbin)
    mmin = center_of_masses.mass.min()
    for bi in binaries:
        mp = bi.mass
        ms = (
            numpy.random.uniform(mmin.value_in(units.MSun), mp.value_in(units.MSun))
            | units.MSun
        )
        a = 0.0 | units.au
        e = 0.0
        while bi.radius > a * (1.0 - e):
            a = random_semimajor_axis(amin, amax)
            e = numpy.sqrt(numpy.random.random())
        print("Orbit a=", a.in_(units.au), "e=", e)

        nb = new_binary_orbit(mp, ms, a, e)
        nb.position += bi.position
        nb.velocity += bi.velocity
        nb = singles_in_binaries.add_particles(nb)
        nb.type = "star"
        nb[0].name = "primary"
        nb[1].name = "secondary"
        nb[1].radius = ZAMS_radius(nb[1].mass)

        nb.radius = 0.01 * a

        bi.radius = 3 * a
        binary_particle = bi.copy()
        binary_particle.child1 = nb[0]
        binary_particle.child2 = nb[1]
        binary_particle.semi_major_axis = a
        binary_particle.eccentricity = e
        resulting_binaries.add_particle(binary_particle)

    single_stars = center_of_masses - binaries
    # return single_stars, resulting_binaries, singles_in_binaries
    single_stars.add_particles(singles_in_binaries)
    return single_stars


def calculate_orbital_elementss(bi, converter):
    kep = new_kepler(converter)
    comp1 = bi.child1
    comp2 = bi.child2
    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a, e = kep.get_elements()
    kep.stop()
    return a, e


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option(
        "-f", dest="infile", default="Plummer.amuse", help="input filename [%default]"
    )
    result.add_option(
        "-F", dest="outfile", default=None, help="output filename [%default]"
    )
    result.add_option(
        "-m",
        unit=units.MSun,
        type="float",
        dest="mass",
        default=0 | units.MSun,
        help="mass for binary component [%default]",
    )
    result.add_option(
        "--name",
        dest="name",
        default="star",
        help="name the newly added companion star [%default]",
    )
    result.add_option(
        "--host_name",
        dest="host_name",
        default="star",
        help="name the newly added companion star [%default]",
    )
    result.add_option(
        "-a",
        unit=units.au,
        type="float",
        dest="semimajor_axis",
        default=1.0 | units.au,
        help="binary separation [%default]",
    )
    result.add_option(
        "-e",
        type="float",
        dest="eccentricity",
        default=0.0,
        help="binary eccentricity [%default]",
    )
    result.add_option(
        "--seed",
        dest="seed",
        type="int",
        default=None,
        help="random number seed [%default]",
    )
    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    outfile = o.outfile

    if o.outfile == None:
        outfile = o.infile

    bodies = read_set_from_file(o.filename, "hdf5", close_file=True)
    selected_stars = bodies[bodies.name == o.host_name]

    stars = make_secondaries(selected_stars, Nbin, o.amin, o.amax)
    time = 0 | units.Myr
    write_set_to_file(
        stars, outfile, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )
