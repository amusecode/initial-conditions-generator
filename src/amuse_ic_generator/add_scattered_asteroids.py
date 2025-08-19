import argparse

import numpy as np

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles
from amuse.io import write_set_to_file, read_set_from_file
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk

from move_bodies_to_stellar_position import move_bodies_to_stellar_position


def test_distribution(bodies):
    asteroids = bodies[bodies.type == b"debris"]
    star = bodies[bodies.type == b"star"]
    print(star)
    print(asteroids)
    p = Particles()
    p.add_particle(star)

    sma = [] | units.au
    ecc = []
    for ai in asteroids:
        p.add_particle(ai)
        M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(
            p, G=constants.G
        )
        sma.append(a)
        ecc.append(e)
        p.remove_particle(ai)

    # import matplotlib.pyplot as plt
    # plt.scatter(sma.value_in(units.au), ecc)
    # plt.semilogx()
    # plt.show()


def add_scattered_asteroids(star, Ndisk_per_MSun, mdisk, qmin, qmax, emin, emax):
    Ndisk = int(Ndisk_per_MSun * star.mass.value_in(units.MSun))
    print(f"Run with Ndisk = {Ndisk}")
    # converter = nbody_system.nbody_to_si(star.mass, qmin)
    e = emin + (emax - emin) * np.random.random(Ndisk)
    lqmin = np.log10(qmin.value_in(units.au))
    lqmax = np.log10(qmax.value_in(units.au))
    q = 10 ** np.random.uniform(lqmin, lqmax, Ndisk) | units.au
    a = q / (1 - e)
    asteroids = Particles(0)
    for i in range(Ndisk):
        ang_1, ang_2, ang_3 = np.random.uniform(0, 2 * np.pi, 3)
        b = new_binary_from_orbital_elements(
            mass1=star.mass,
            mass2=mdisk,
            semimajor_axis=a[i],
            eccentricity=e[i],
            true_anomaly=ang_1 | units.rad,
            inclination=0 | units.rad,
            longitude_of_the_ascending_node=ang_2 | units.rad,
            argument_of_periapsis=ang_3 | units.rad,
            G=constants.G,
        )
        b[1].position -= b[0].position
        b[1].velocity -= b[0].velocity
        b[1].semimajor_axis = a[i]
        b[1].eccentricity = e[i]
        b[1].inclination = 0
        asteroids.add_particle(b[1])
    asteroids.mass = mdisk
    asteroids.name = "comet"
    asteroids.type = "debris"
    asteroids.position += star.position
    asteroids.velocity += star.velocity
    return asteroids


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--filename", default=None, help="input filename")
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    parser.add_argument(
        "--mdisk",
        type=units.kg,
        default=0.0 | units.kg,
        help="mass of individual disk particles",
    )
    parser.add_argument(
        "--Ndisk",
        type=int,
        default=100,
        help="number of disk particles",
    )
    parser.add_argument(
        "--qmin",
        type=units.au,
        default=-5.2044 | units.au,
        help="minimal closest appraoch",
    )
    parser.add_argument(
        "--qmax",
        type=units.au,
        default=-30.33 | units.au,
        help="maximal closest appraoch",
    )
    parser.add_argument(
        "--emin",
        type=float,
        default=0.10,
        help="minimum eccentricity",
    )
    parser.add_argument(
        "--emax",
        type=float,
        default=0.90,
        help="maximum eccentricity",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=-1,
        help="random number seed",
    )
    return parser


def get_inner_planet(planets):
    amin = planets.semimajor_axis.min()
    inner_planet = planets[planets.semimajor_axis <= amin]
    return inner_planet[0]


def get_outer_planet(planets):
    amax = planets.semimajor_axis.max()
    outer_planet = planets[planets.semimajor_axis >= amax]
    return outer_planet[0]


def main():
    args = new_argument_parser().parse_args()

    if args.seed > 0:
        np.random.seed(args.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(args.filename, close_file=True)
    print(bodies)
    stars = bodies[bodies.type == "star"]
    planets = bodies[bodies.type == "planet"]
    if args.qmin < 0 | units.au:
        inner_planet = get_inner_planet(planets)
        qmin = 0.5 * inner_planet.semimajor_axis * (1.0 - inner_planet.eccentricity)
    else:
        qmin = args.qmin
    if args.qmax < 0 | units.au:
        outer_planet = get_outer_planet(planets)
        qmax = 2 * outer_planet.semimajor_axis * (1.0 - outer_planet.eccentricity)
    else:
        qmax = args.qmax
    if len(stars) == 0:
        stars = bodies[bodies.name == "star"]

    debris = add_scattered_asteroids(
        stars[0], args.Ndisk, args.mdisk, qmin, qmax, args.emin, args.emax
    )
    bodies.add_particles(debris)
    # bodies.move_to_center()
    move_bodies_to_stellar_position(bodies)

    # test_distribution(bodies)
    # import matplotlib.pyplot as plt
    # plt.scatter(
    #     bodies.x.value_in(units.parsec),
    #     bodies.y.value_in(units.parsec))
    # plt.show()

    time = 0 | units.Myr
    if args.outfile is None:
        filename = "added_scattered_asteroids.amuse"
    else:
        filename = args.outfile
    write_set_to_file(
        bodies, filename, "amuse", timestamp=time, append_to_file=False, version="2.0"
    )


if __name__ == "__main__":
    main()
