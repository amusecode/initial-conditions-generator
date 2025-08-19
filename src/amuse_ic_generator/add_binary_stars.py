"""Makes binary stars.

Takes an AMUSE file as input, selects the stars from it within a specified mass
range, and adds binary companions to a specified fraction (randomly selected)
of these stars.
"""

import argparse

import numpy as np

from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.community.kepler import Kepler
from amuse.io import read_set_from_file, write_set_to_file


def semimajor_axis_to_orbital_period(a, mass_total):
    """Returns the orbital period for a given semimajor axis and total mass.
    """
    return 2 * np.pi * (a**3 / (constants.G * mass_total)).sqrt()


def new_kepler(converter):
    """Returns a new Kepler worker.
    """
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.set_longitudinal_unit_vector(1.0, 0.0, 0.0)
    kepler.set_transverse_unit_vector(0.0, 1.0, 0)
    return kepler


def new_binary_orbit(mass1, mass2, semi_major_axis, eccentricity=0):
    """Returns a binary orbit for the given parameters.
    """
    total_mass = mass1 + mass2
    mass_fraction_particle_1 = mass1 / (total_mass)

    binary = Particles(2)
    binary[0].mass = mass1
    binary[1].mass = mass2

    mu = constants.G * total_mass

    velocity_perihelion = np.sqrt(
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
    """Generates a random semimajor axis for a binary star system."""
    lamin = np.log10(amin.value_in(units.au))
    lamax = np.log10(amax.value_in(units.au))
    rnd_min = lamin
    rnd_max = lamax
    rnd = np.random.uniform(rnd_min, rnd_max, 1)
    a = 10**rnd
    return a | units.au


def random_semimajor_axis_ppe(mass_primary, mass_secondary, a_min, a_max):
    """Generates a random semimajor axis for a binary star system.

    See Eggleton 2006 Equation 1.6.3 (2006epbm.book.....E).
    """
    orbital_period_max = semimajor_axis_to_orbital_period(
        a_max, mass_primary + mass_secondary
    ).value_in(units.day)
    orbital_period_min = semimajor_axis_to_orbital_period(
        a_min, mass_primary + mass_secondary
    ).value_in(units.day)
    mpf = (mass_primary.value_in(units.MSun) ** 2.5) / 5.0e4
    rnd_max = (orbital_period_max * mpf) ** (1.0 / 3.3) / (
        1 + (orbital_period_min * mpf) ** (1.0 / 3.3)
    )
    rnd_min = (orbital_period_min * mpf) ** (1.0 / 3.3) / (
        1 + (orbital_period_max * mpf) ** (1.0 / 3.3)
    )
    rnd_max = min(rnd_max, 1)
    rnd = np.random.uniform(rnd_min, rnd_max, 1)
    orbital_period = ((rnd / (1.0 - rnd)) ** 3.3) / mpf | units.day
    print(orbital_period_min, orbital_period_max, orbital_period)

    mass_total = mass_primary + mass_secondary
    semi_major_axis = (
        (constants.G * mass_total) * (orbital_period / (2 * np.pi)) ** 2
    ) ** (1.0 / 3.0)
    return semi_major_axis


def zams_radius(mass):
    """Returns the ZAMS radius of a star with the given mass.
    """
    # log_mass = np.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun)) ** 2
    # alpha = 0.08353 + 0.0565 * log_mass
    # beta = 0.01291 + 0.2226 * log_mass
    # gamma = 0.1151 + 0.06267 * log_mass
    r_zams = (
        pow(mass.value_in(units.MSun), 1.25)
        * (0.1148 + 0.8604 * mass_sq)
        / (0.04651 + mass_sq)
    )

    return r_zams | units.RSun


def make_secondaries(center_of_masses, number_of_binaries, a_min, a_max):
    """Generates binary pairs from a list of stars.
    Takes a random subset of the given stars, and generates secondary
    companions for these. The semi-major axis is chosen randomly within the
    given range, and the eccentricity is chosen randomly.
    Returns a particle set containing the systems (type: star), the primaries
    in the systems (type: primary) and the secondaries in the systems (type:
    secondary).
    """
    resulting_binaries = Particles()
    singles_in_binaries = Particles()
    binaries = center_of_masses.random_sample(number_of_binaries)
    mmin = center_of_masses.mass.min()
    for bi in binaries:
        mp = bi.mass
        ms = (
            np.random.uniform(mmin.value_in(units.MSun), mp.value_in(units.MSun))
            | units.MSun
        )
        a = 0.0 | units.au
        e = 0.0
        # FIXME: this doesn't work since bi.radius is unset / zero!
        while bi.radius > a * (1.0 - e):
            a = random_semimajor_axis(a_min, a_max)
            e = np.sqrt(np.random.random())
        print("Orbit a=", a.in_(units.au), "e=", e)

        nb = new_binary_orbit(mp, ms, a, e)
        nb.position += bi.position
        nb.velocity += bi.velocity
        nb = singles_in_binaries.add_particles(nb)
        nb.type = "star"
        nb[0].name = "primary"
        nb[1].name = "secondary"
        nb[1].radius = zams_radius(nb[1].mass)

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


def calculate_orbital_elements(bi, converter):
    """Calculates the orbital elements of a binary system.
    """
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


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "-f", "--filename", default="Plummer.amuse", help="input filename"
    )
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "--f_binaries",
        type=float,
        default=0.5,
        help="binary fraction",
    )
    parser.add_argument(
        "-m",
        "--mmin",
        type=units.MSun,
        default=0 | units.MSun,
        help="minimum stellar mass for binary (inclusive)",
    )
    parser.add_argument(
        "-M",
        "--mmax",
        type=units.MSun,
        default=1000 | units.MSun,
        help="maximum stellar mass for binary (exclusive)",
    )
    parser.add_argument(
        "-a",
        "--amin",
        type=units.au,
        default=1.0 | units.au,
        help="minimum binary separation",
    )
    parser.add_argument(
        "-A",
        "--amax",
        type=units.au,
        default=1.0e8 | units.au,
        help="maximum binary separation",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="random number seed",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    outfile = args.outfile

    if outfile is None:
        outfile = "single_stars_and_binaries.amuse"

    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies[0])
    stars = bodies[bodies.type == "star"]
    print(f"star: {stars[0]}")
    print(f"Mass limits: {args.mmin} {args.mmax}")
    selected_stars = stars[stars.mass >= args.mmin]
    selected_stars = selected_stars[selected_stars.mass < args.mmax]
    number_of_binaries = int(args.f_binaries * len(selected_stars))
    print(f"Number of stars: {len(stars)} {len(selected_stars)} {number_of_binaries}")
    if number_of_binaries == 0:
        print("No binaries added (try increasing --f_binaries or relaxing mass limits)")
        write_set_to_file(
            selected_stars, outfile, "amuse", append_to_file=False, version="2.0"
        )
        return

    binary_stars = make_secondaries(selected_stars, number_of_binaries, args.amin, args.amax)
    time = 0 | units.Myr

    all_stars = stars
    all_stars.remove_particles(selected_stars)
    all_stars.add_particles(binary_stars)

    write_set_to_file(
        all_stars, outfile, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
