import argparse

import numpy as np
import matplotlib.pyplot as plt

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles, Particle
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.io import write_set_to_file

from add_asteroids_to_solar_system import True_anomaly_from_mean_anomaly
from make_cluster_of_stars import make_fractal_cluster
from determine_orbital_elements import calculate_orbital_elements_for_single_planet
from make_planets_oligarch import make_planets_oligarch
from rotate import rotate_bodies_isotropically


def ZAMS_radius(mass):
    log_mass = np.log10(mass.value_in(units.MSun))
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


def make_secondary(
    parent,
    companion_name,
    mass,
    semimajor_axis,
    eccentricity,
    inclination,
    mean_anomaly,
    LoAn,
    Aop,
    ctype="star",
):

    Ta = True_anomaly_from_mean_anomaly(np.deg2rad(mean_anomaly), eccentricity)
    bs = new_binary_from_orbital_elements(
        parent.mass,
        mass,
        semimajor_axis,
        eccentricity,
        Ta,
        inclination,
        LoAn,
        Aop,
        G=constants.G,
    )
    companion = bs[1]
    companion.position -= bs[0].position
    companion.velocity -= bs[0].velocity
    companion.position += parent.position
    companion.velocity += parent.velocity
    companion.type = "star"
    companion.name = companion_name
    # companion.host = parent.key
    companion.radius = ZAMS_radius(companion.mass)
    return companion


def random_semimajor_axis(amin, amax):
    lamin = np.log10(amin.value_in(units.au))
    lamax = np.log10(amax.value_in(units.au))
    rnd_min = lamin
    rnd_max = lamax
    rnd = np.random.uniform(rnd_min, rnd_max, 1)
    a = 10**rnd
    return a | units.au


def make_companions_for_primaries(primaries, amin=1 | units.au, amax=1000 | units.au):
    companion_name = "secondary"
    companions = Particles()

    for primary in primaries:
        c = make_companion_for_primary(primary, amin, amax)
        companions.add_particles(c.as_set())
    return companions


def make_companion_for_primary(primary, amin, amax):

    companion_name = "secondary"
    companion_mass = np.random.random() * primary.mass
    companion_mass = max(companion_mass, 0.1 | units.MSun)
    companion_radius = ZAMS_radius(companion_mass)
    rsum = primary.radius + companion_radius
    a = random_semimajor_axis(amin, amax)[0]
    e = np.sqrt(np.random.random())
    while rsum > a * (1.0 - e):
        a = random_semimajor_axis(amin, amax)[0]
        e = np.sqrt(np.random.random())
    if a * (1 - e) < 3 * rsum:
        a = a * (1 - e)
        e = 0
    print("Orbit a=", a.in_(units.au), "e=", e)
    inc = np.random.uniform(0, np.pi)  # what units..
    mean_anomaly = np.random.uniform(0, 2 * np.pi)
    LoAn = np.random.uniform(0, 2 * np.pi)
    Aop = np.random.uniform(0, 2 * np.pi)

    companion = make_secondary(
        primary,
        companion_name,
        companion_mass,
        a,
        e,
        inc,
        mean_anomaly,
        LoAn,
        Aop,
        ctype="star",
    )
    return companion


def Hill_radius(Mstar, a, Mplanet):
    return a * (Mplanet / (3.0 * Mstar)) ** (1.0 / 3.0)


def make_tertiary_companions(primaries, companions):
    tertiaries = Particles()
    planets = Particles()
    for primary, companion in zip(primaries, companions):
        a, e, i = calculate_orbital_elements_for_single_planet(primary, companion)
        print("Orbit:", a.in_(units.au), e)
        r = np.random.random()
        if r <= 0.5:  # make triple
            rr = np.random.random()
            if rr < 0.5:  # cirum-secondary
                com = companion
                amin = max(0.01 * a * (1 - e), 1 | units.au)  # 10*companion.radius)
                amax = min(amin, 0.5 * a * (1 - e))
                # t.host = com.key
            else:  # cirum-binary
                bs = Particles()
                bs.add_particle(primary)
                bs.add_particle(companion)
                com = Particle()
                com.mass = bs.mass.sum()
                com.position = bs.center_of_mass()
                com.velocity = bs.center_of_mass_velocity()
                com.radius = a * (1 - e)
                amin = 2 * a * (1 + e)
                amax = max(100 * a * (1 + e), 1000 | units.au)
            t = make_companion_for_primary(com, amin, amax)
            t.name = "tertiary"
            t.radius = ZAMS_radius(t.mass)
            # t.host = None
            tertiaries.add_particle(t)
        else:  # make circum-stellar planets
            if r < 0.25:
                com = primary
            else:
                com = companion
            amin = min(0.5 * a * (1 - e), 1 | units.au)  # 10*companion.radius)
            amax = max(10 * amin, 0.5 * a * (1 - e))
            p = make_planets_for_planetary(com.as_set(), amin, amax)
            planets.add_particle(p)
    return tertiaries, planets


def make_planets_for_planetaries(planetaries):
    planets = Particles()
    for star in planetaries:
        rmin = max(10 | units.au, 300 * star.radius)
        rmax = max(300 | units.au, 30000 * star.radius)
        p = make_planets_for_planetary(star, rmin, rmax)
        planets.add_particles(p)
    return planets


def make_planets_for_planetary(star, rmin, rmax):
    print("star=", star, rmin, rmax)

    disk_mass = 0.01 * star.mass
    p = make_planets_oligarch(star.mass, star.radius, rmin, rmax, disk_mass)
    p.type = "planet"
    p.name = "planet"
    # p.host=star.key

    p = rotate_bodies_isotropically(p)
    p.position += star.position
    p.velocity += star.velocity
    return p


def plot_bodies(bodies):
    length_unit = units.au
    stars = bodies[bodies.type == "star"]
    singletons = stars[stars.name == "singleton"]
    primaries = stars[stars.name == "primary"]
    secondaries = stars[stars.name == "secondary"]
    tertiaries = stars[stars.name == "tertiary"]

    planetaries = stars[stars.name == "planetary"]
    planets = bodies[bodies.name == "planet"]

    m = 30 * singletons.mass.value_in(units.MSun)
    plt.scatter(
        singletons.x.value_in(length_unit),
        singletons.y.value_in(length_unit),
        s=m,
        c="orange",
        alpha=0.5,
    )
    m = 30 * primaries.mass.value_in(units.MSun)
    plt.scatter(
        primaries.x.value_in(length_unit),
        primaries.y.value_in(length_unit),
        s=m,
        c="b",
        alpha=0.5,
    )
    m = 30 * secondaries.mass.value_in(units.MSun)
    plt.scatter(
        secondaries.x.value_in(length_unit),
        secondaries.y.value_in(length_unit),
        s=m,
        c="g",
    )
    m = 30 * tertiaries.mass.value_in(units.MSun)
    plt.scatter(
        tertiaries.x.value_in(length_unit),
        tertiaries.y.value_in(length_unit),
        s=m,
        c="r",
    )

    m = 30 * planetaries.mass.value_in(units.MSun)
    plt.scatter(
        planetaries.x.value_in(length_unit),
        planetaries.y.value_in(length_unit),
        s=m,
        c="y",
    )
    m = 1
    plt.scatter(
        planets.x.value_in(length_unit), planets.y.value_in(length_unit), s=m, c="k"
    )
    # planets = bodies[bodies.type == "planet"]
    # debris = bodies[bodies.type == "debris"]
    # plt.scatter(bodies.x.value_in(units.au), bodies.y.value_in(units.au), s=0.1, c="k")
    # plt.scatter(stars.x.value_in(units.au), stars.y.value_in(units.au), s=30, c="y")
    # plt.scatter(planets.x.value_in(units.au), planets.y.value_in(units.au), s=10, c="b")
    # plt.scatter(debris.x.value_in(units.au), debris.y.value_in(units.au), s=1, c="k")
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "--nstars",
        type=int,
        default=3,
        help="number of stars",
    )
    parser.add_argument(
        "--mass",
        type=units.MSun,
        default=1 | units.MSun,
        help="total cluster mass",
    )
    parser.add_argument(
        "--mmin",
        type=units.MSun,
        default=-0.1 | units.MSun,
        help="minimum stellar mass",
    )
    parser.add_argument(
        "--mmax",
        type=units.MSun,
        default=100 | units.MSun,
        help="maximum stellar mass",
    )
    parser.add_argument("-Q", "--Qvir", type=float, default=0.5, help="virial ratio")
    parser.add_argument(
        "--radius",
        type=units.parsec,
        default=1.0 | units.parsec,
        help="cluster radius",
    )
    parser.add_argument("--name", default="star", help="stellar name")
    parser.add_argument(
        "--seed",
        type=int,
        default=-1,
        help="random number seed",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    if args.seed > 0:
        np.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    n_star = args.nstars
    mmin = args.mmin
    cluster_radius = args.radius

    names = "star"
    masses = new_kroupa_mass_distribution(n_star, mass_min=mmin)
    converter = nbody_system.nbody_to_si(masses.sum(), cluster_radius)
    bodies = make_fractal_cluster(masses, names, converter)
    bodies.radius = ZAMS_radius(bodies.mass)
    bodies.type = "star"
    bodies.name = "star"
    # bodies.host = None
    bodies.scale_to_standard(convert_nbody=converter, virial_ratio=args.Qvir)
    print(bodies)

    # Select stars for single, binary companion or planetary host
    singletons = (
        bodies[bodies.mass < (0.6 | units.MSun)]
        + bodies[bodies.mass > (3 | units.MSun)]
    )
    singletons.name = "singleton"
    n = len(bodies) - len(singletons)
    nbinaries = int(n / 2)
    nplanetaries = n - nbinaries
    print("N=", len(bodies), len(singletons), nbinaries, nplanetaries)
    primaries = (bodies - singletons).random_sample(nbinaries)
    primaries.name = "primary"
    planetaries = bodies - singletons - primaries
    planetaries.name = "planetary"

    companion_stars = make_companions_for_primaries(primaries)
    bodies.add_particles(companion_stars)

    tertiaries, planets = make_tertiary_companions(primaries, companion_stars)
    bodies.add_particles(tertiaries)
    bodies.add_particles(planets)

    planets = make_planets_for_planetaries(planetaries)
    bodies.add_particles(planets)

    bodies.move_to_center()

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "star_cluster.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, overwrite_file=True, version="2.0"
    )
    plot_bodies(bodies)
    print("N=", len(bodies))


if __name__ == "__main__":
    main()
