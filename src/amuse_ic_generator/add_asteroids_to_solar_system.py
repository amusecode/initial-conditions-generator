import os
import urllib
import argparse

import numpy as np
import matplotlib.pyplot as plt

from amuse.units import units, nbody_system, constants
from amuse.units.trigo import pi, sin, cos
from amuse.datamodel import Particles, Particle
from amuse.community.kepler.interface import Kepler
from amuse.io import read_set_from_file, write_set_to_file
from amuse.ic.solar_system_moons import new_lunar_system_in_time
from amuse.ic.solar_system_moons import solar_system_in_time
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.plot import scatter


def new_kepler():
    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.au)
    kepler = Kepler(converter)
    kepler.initialize_code()
    return kepler


def get_position(
    mass_sun,
    mass_planet,
    ecc,
    semi,
    mean_anomaly,
    incl,
    argument,
    longitude,
    delta_t=0.0 | units.day,
):
    """
    cartesian position and velocity from orbital elements,
    where the orbit is evolved from given mean_anomaly
    by time delta_t
    argument -- argument of perihelion
    longitude -- longitude of ascending node
    """
    kepler = new_kepler()
    kepler.initialize_from_elements(
        mass=(mass_sun + mass_planet), semi=semi, ecc=ecc, mean_anomaly=mean_anomaly
    )
    kepler.transform_to_time(time=delta_t)
    r = kepler.get_separation_vector()
    v = kepler.get_velocity_vector()

    kepler.stop()

    a1 = (
        [cos(longitude), -sin(longitude), 0.0],
        [sin(longitude), cos(longitude), 0.0],
        [0.0, 0.0, 1.0],
    )
    a2 = (
        [1.0, 0.0, 0.0],
        [0.0, cos(incl), -sin(incl)],
        [0.0, sin(incl), cos(incl)],
    )
    a3 = (
        [cos(argument), -sin(argument), 0.0],
        [sin(argument), cos(argument), 0.0],
        [0.0, 0.0, 1.0],
    )
    A = np.dot(np.dot(a1, a2), a3)
    print(A, r)

    r_vec = np.dot(A, np.reshape(r, 3, "F"))
    v_vec = np.dot(A, np.reshape(v, 3, "F"))

    # for relative vectors
    r[0] = r_vec[0]
    r[1] = r_vec[1]
    r[2] = r_vec[2]
    v[0] = v_vec[0]
    v[1] = v_vec[1]
    v[2] = v_vec[2]

    return r, v


def true_anomaly_from_mean_anomaly(Ma, e):
    Ta = (
        Ma
        + (2 * e - e**3 / 4.0) * sin(Ma)
        + (5 * e**2 / 4.0) * sin(2 * Ma)
        + (13.0 / 12.0) * e**3 * sin(3 * Ma)
    )
    return Ta


def construct_particle_set_from_orbital_elements(
    name, mass, a, ecc, inc, Ma, Aop, LoAn, Mparent=1 | units.MSun
):
    print(
        "length:", len(a), len(ecc), len(inc), len(Ma), len(Aop), len(LoAn), len(name)
    )

    p = Particles(len(a))
    print(name)
    p.name = name
    p.type = "asteroid"
    p.host = "Sun"
    p.mass = mass

    for i in range(len(p)):
        Ta = true_anomaly_from_mean_anomaly(Ma[i].in_(units.deg), ecc[i])
        b = new_binary_from_orbital_elements(
            Mparent, p[i].mass, a[i], ecc[i], Ta, inc[i], LoAn[i], Aop[i], G=constants.G
        )
        p[i].position = b[1].position - b[0].position
        p[i].velocity = b[1].velocity - b[0].velocity

        rho = 2.0 | units.g / units.cm**3
        p.radius = (p.mass / rho) ** (1.0 / 3.0)

    return p


def read_orbital_elements_from_MPCORB(filename="MPCORB.DAT", n=-1):
    a = [] | units.au
    ecc = []
    inc = [] | units.deg
    Ma = [] | units.deg
    Aop = [] | units.deg
    LoAn = [] | units.deg
    classify = []
    name = []
    Nobs = []
    U = []
    MPC_data = False
    for line in open(filename):
        if "00001" in line[0:6]:
            MPC_data = True
        if MPC_data and len(line.split()) > 10:
            No = line[117:122].strip()
            if len(No) > 0 and (
                "E" not in line[105:107]
                or "D" not in line[105:107]
                or "D" not in line[105:107]
            ):
                Nobs.append(int(No))
                Ma.append(units.deg(line[26:35]))
                Aop.append(units.deg(line[37:46]))
                LoAn.append(units.deg(line[48:57]))
                inc.append(units.deg(line[59:68]))
                ecc.append(float(line[70:79]))
                a.append(units.au(line[92:103]))
                U.append(line[105:107])
                classify.append(line[161:165])
                name.append(line[166:194])
                if n > 0 and len(a) >= n:
                    break

    return name, a, ecc, inc, Ma, Aop, LoAn


def read_orbital_elements_from_CometEls(filename="CometEls", n=-1):
    a = [] | units.au
    ecc = []
    inc = [] | units.deg
    Ma = [] | units.deg
    Aop = [] | units.deg
    LoAn = [] | units.deg
    classify = []
    name = []
    Nobs = []
    U = []
    MPC_data = False
    for line in open(filename):
        print(line)
        if True:
            p = float(line[31:40]) | units.au
            # ecc.append(float(line[50:57]))
            e = float(line[41:50])
            if True:  # e>0.9:
                ecc.append(e)
                name.append(line[0:12])
                print(name[-1])
                if ecc[-1] == 1:
                    ecc[-1] -= 1.0e-5
                a.append(p / (1 - ecc[-1]))
                Aop.append(units.deg(line[51:60]))
                LoAn.append(units.deg(line[61:70]))
                inc.append(units.deg(line[71:80]))
                yop = float(line[14:18])
                mop = float(line[19:21])
                dop = float(line[22:29])
                Ma.append((np.random.random() * 360 - 180) | units.deg)
                print(a[-1], ecc[-1], inc[-1], Aop[-1], LoAn[-1])
                if n > 0 and len(a) >= n:
                    break

    return name, a, ecc, inc, Ma, Aop, LoAn


def read_orbital_elements_from_MinorPlanetCenter(filename="MPCORB.DAT", n=-1):

    if filename == "MPCORB.DAT":
        name, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_from_MPCORB(
            filename, n
        )
        return name, a, ecc, inc, Ma, Aop, LoAn
    elif filename == "CometEls.txt":
        name, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_from_CometEls(
            filename, n
        )
        return name, a, ecc, inc, Ma, Aop, LoAn

    a = [] | units.au
    ecc = []
    inc = [] | units.deg
    Ma = [] | units.deg
    Aop = [] | units.deg
    LoAn = [] | units.deg
    name = []
    for line in open(filename):
        print(line)
        Ma.append(units.deg(line[26:35]))
        Aop.append(units.deg(line[37:46]))
        LoAn.append(units.deg(line[48:57]))
        inc.append(units.deg(line[59:68]))
        ecc.append(float(line[70:79]))
        a.append(units.au(line[92:103]))
        name.append(line[166:194])
        if n > 0 and len(a) >= n:
            break
    return name, a, ecc, inc, Ma, Aop, LoAn


def add_asteroids_to_solar_system(
    solar_system, MPC_filename, time_JD=2457099.5 | units.day, n=-1
):

    # download the asteroid ephemerids file
    if not os.path.isfile(MPC_filename):
        print("download ", MPC_filename)
        url = f"https://www.minorplanetcenter.net/iau/MPCORB/{MPC_filename}"
        datafile = urllib.request.urlopen(url)
        with open(MPC_filename, "wb") as f:
            f.write(datafile.content)

    sun = solar_system[solar_system.name == "Sun"][0]

    name, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_from_MinorPlanetCenter(
        filename=MPC_filename, n=n
    )
    mass = 1000 | units.kg
    p = construct_particle_set_from_orbital_elements(
        name, mass, a, ecc, inc, Ma, Aop, LoAn, sun.mass
    )

    print(p)
    print("number of particles N=", len(p))
    solar_system.add_particles(p)
    solar_system.move_to_center()

    return solar_system


def new_lunar_system(Julian_date=-1 | units.day):
    return new_lunar_system_in_time(Julian_date)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-n", "--number_of_asteroids", type=int, default=10, help="number of asteroids"
    )
    parser.add_argument(
        "--MPC_filename",
        default="NEA.txt",
        help="MPC data filename possibilities include: MPCORB.DAT, NEA.txt, PHA.txt, Distant.txt, Unusual.txt",
    )
    parser.add_argument(
        "-f",
        "--filename",
        default="solar_system.amuse",
        help="input filename",
    )
    parser.add_argument("-F", "--outfile", default=None, help="output filename")
    return parser


def main():
    args = new_argument_parser().parse_args()

    solar_system = read_set_from_file(args.filename, close_file=True)
    Julian_date = solar_system.Julian_date

    solar_system = add_asteroids_to_solar_system(
        solar_system, args.MPC_filename, Julian_date, args.number_of_asteroids
    )
    print(solar_system)

    star = solar_system[solar_system.type == "star"]
    planet = solar_system[solar_system.type == "planet"]
    moon = solar_system[solar_system.type == "moon"]
    asteroid = solar_system[solar_system.type == "asteroid"]
    scatter(star.x, star.y, s=100, c="y")
    scatter(planet.x, planet.y, s=30, c="b")
    scatter(moon.x, moon.y, s=10, c="r")
    scatter(asteroid.x, asteroid.y, s=1, c="k")
    plt.show()
    print(planet)

    if args.outfile:
        write_set_to_file(
            solar_system,
            args.outfile,
            "amuse",
            append_to_file=False,
            overwite_file=True,
            version="2.0",
        )


if __name__ == "__main__":
    main()
