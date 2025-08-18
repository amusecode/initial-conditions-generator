import os
import argparse
import urllib

import numpy as np
import matplotlib.pyplot as plt

from amuse.units.trigo import pi, sin, cos, atan2, sqrt
from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles, Particle
from amuse.plot import scatter
from amuse.io import write_set_to_file
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

from orbital_elements_to_cartesian import orbital_elements_to_pos_and_vel
from orbital_elements_to_cartesian import orbital_period


def eccentric_anomaly(mean_anomaly, e):
    ecc_anomaly = mean_anomaly + 2 * pi * (
        e * sin(mean_anomaly) + 0.5 * e * e * sin(2 * mean_anomaly)
    )
    m = ecc_anomaly - e * sin(ecc_anomaly)
    de = (mean_anomaly - m) / (1 - e * cos(ecc_anomaly))
    ecc_anomaly += de
    while de >= 0.001:
        m = ecc_anomaly - e * sin(ecc_anomaly)
        de = (mean_anomaly - m) / (1 - e * cos(ecc_anomaly))
        ecc_anomaly += de
    return ecc_anomaly


def True_anomaly_from_mean_anomaly(Ma, e):
    Ta = (
        Ma
        + (2 * e - e**3 / 4.0) * sin(Ma)
        + (5 * e**2 / 4.0) * sin(2 * Ma)
        + (13.0 / 12.0) * e**3 * sin(3 * Ma)
    )
    return Ta


def construct_particle_set_from_orbital_elements(
    name, mass, a, ecc, inc, Ma, Aop, LoAn, parent
):
    # print("length:", len(a), len(ecc), len(inc), len(Ma), len(Aop), len(LoAn), len(name))

    p = Particles(len(a))
    p.name = name
    p.type = "planet"
    p.host = parent.name
    Earth = p[p.name == "EM"]
    Earth.name = "EarthMoon"

    for i in range(len(p)):
        p[i].mass = mass[i]
        Ta = True_anomaly_from_mean_anomaly(np.deg2rad(Ma[i]), ecc[i])
        # print("True Anomaly:", Ta)
        b = new_binary_from_orbital_elements(
            parent.mass,
            p[i].mass,
            a[i],
            ecc[i],
            Ta,
            inc[i],
            LoAn[i],
            Aop[i],
            G=constants.G,
        )
        p[i].position = b[1].position - b[0].position
        p[i].velocity = b[1].velocity - b[0].velocity
        rho = 1.0 | units.g / units.cm**3
        p.radius = (p.mass / rho) ** (1.0 / 3.0)

    return p


def read_orbital_elements_for_planets(
    filename="MPCORB.DAT", Julian_date=2451545.0 | units.day
):
    # f = open(fdir+"MPCORB.DAT", "r")
    planet_names = [
        "Mercury",
        "Venus",
        "EM Bary",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",
    ]
    planet_mass = [
        3.302e23,
        48.685e23,
        5.97219e24,
        6.4185e23,
        1898.13e24,
        5.68319e26,
        86.8103e24,
        102.41e24,
        1.309e22,
    ] | units.kg

    sma = [] | units.au
    ecc = []
    inc = []  # degrees
    Ma = []  # degrees
    Aop = []  # degrees
    LoAn = []  # degrees
    name = []
    for line in open(filename):
        # see: https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
        for pi in planet_names:
            if pi in line:
                l = line.split(" ")
                name.append(l[0])
                sma.append(float(line[8:20]) | units.au)
                ecc.append(float(line[21:36]))
                inc.append(float(line[37:52]))
                L = float(line[53:70])
                Lperi = float(line[71:86])
                Lnode = float(line[87:102])
                # if n>0 and len(sma)>=n:
                #    break
                Aop.append(Lperi)
                LoAn.append(Lnode)
                T = (Julian_date - (2451545.0 | units.day)) / (36525 | units.day)
                omega = Lperi
                M = L - omega
                Ma.append(M)

    return name, planet_mass, sma, ecc, inc, Ma, Aop, LoAn


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", dest="filename", default="MPCORB.h5", help="read to file")
    parser.add_argument("-F", dest="outfilename", default=None, help="write to file")
    parser.add_argument(
        "-d",
        unit=units.day,
        type="float",
        dest="Julian_date",
        default=2451545.0 | units.day,
        help="Julan date",
    )
    parser.add_argument(
        "-n",
        dest="n_first",
        type="int",
        default=0,
        help="first planet (Earth==4)",
    )
    parser.add_argument(
        "-N",
        dest="n_last",
        type="int",
        default=10,
        help="last planet (Neptune = 8)[%default]",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()

    sun = Particles(1)
    sun.Julian_date = args.Julian_date
    sun.mass = 1 | units.MSun
    sun.radius = 1 | units.RSun
    sun.position = (0, 0, 0) | units.km
    sun.velocity = (0, 0, 0) | units.kms
    sun.type = "star"
    sun.name = "Sun"
    sun.host = "None"

    if not os.path.isfile("p_elem_t1.txt"):
        print("download p_elem_t1.txt")
        url = "https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt"
        datafile = urllib.request.urlopen(url)
        with open("p_elem_t1.txt", "wb") as f:
            f.write(datafile.content)

    name, mass, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_for_planets(
        filename="p_elem_t1.txt", Julian_date=args.Julian_date
    )
    p = construct_particle_set_from_orbital_elements(
        name, mass, a, ecc, inc, Ma, Aop, LoAn, sun[0]
    )

    sun.add_particles(p[args.n_first : args.n_last])
    sun.move_to_center()

    solar_system = sun

    print("number of particles N=", len(sun))

    if args.outfilename:
        write_set_to_file(
            sun,
            args.outfilename,
            "amuse",
            append_to_file=False,
            overwrite_file=True,
            version="2.0",
        )
    else:
        star = solar_system[solar_system.type == "star"]
        planet = solar_system[solar_system.type == "planet"]
        plt.scatter(star.x.value_in(units.au), star.y.value_in(units.au), s=100, c="y")
        plt.scatter(
            planet.x.value_in(units.au), planet.y.value_in(units.au), s=30, c="b"
        )
        plt.show()
        print(planet)


if __name__ == "__main__":
    main()
