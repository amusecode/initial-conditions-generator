from amuse.lab import *
from math import pi, sin, cos, sqrt, atan2
import numpy


def orbital_period(a, Mtot):
    return 2 * numpy.pi * (a**3 / (constants.G * Mtot)).sqrt()


def get_component_binary_elements(comp1, comp2, kepler, conv):
    mass = conv.to_nbody(comp1.mass + comp2.mass)
    pos = conv.to_nbody(comp2.position - comp1.position)
    vel = conv.to_nbody(comp2.velocity - comp1.velocity)
    kepler.initialize_from_dyn(mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a, e = kepler.get_elements()
    r = kepler.get_separation()
    E, J = kepler.get_integrals()  # per unit reduced mass, note

    return mass, a, e, r, E


def calculate_orbital_elements(primary, secondary, kepler, converter):
    m, a, e, r, E = get_component_binary_elements(primary, secondary, kepler, converter)
    m = converter.to_si(m).as_quantity_in(units.MSun)
    a = converter.to_si(a).as_quantity_in(units.AU)
    r = converter.to_si(r).as_quantity_in(units.AU)
    E = converter.to_si(r).as_quantity_in(units.AU)

    m0 = primary.mass
    m1 = secondary.mass
    return a, e, m0, m1


def orbital_parameters_for_the_planets(bodies, verbose=True):
    from amuse.community.kepler.interface import Kepler

    kepler = Kepler(redirection="none")
    kepler.initialize_code()
    #    kep_converter=nbody_system.nbody_to_si(1|units.MSun, 10|units.AU)
    converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1 | units.AU)
    a = [] | units.AU
    e = []
    m = [] | units.MSun
    name = []
    for bi in bodies[1:]:
        ai, ei, M, ms = calculate_orbital_elements(bodies[0], bi, kepler, converter)
        name.append(bi.name)
        a.append(ai)
        e.append(ei)
        m.append(ms)
    kepler.stop()
    if verbose:
        for i in range(len(a)):
            print("Planet: ", name[i], a[i], e[i], m[i])
    return a, e


# Solve Kepler equation by iterating: M = E - e sin E
# Lit.: Sterne, T.E., 1960, An introduction to Celestial Mechanics, p. 13-14
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


def mean_anomaly(time, tp, P):
    MA = 2.0 * pi * (time - tp) / P
    while (MA < 0) and (MA > 2.0 * pi):
        if MA < 0:
            time = time + P
        else:
            time = time - P
        MA = 2.0 * pi * (time - tp) / P
    return MA


def orbital_elements_to_pos_and_vel(
    time, a, ecc, inc, omra, omega, tp, Mbh, mstar, MA=-1
):

    P = orbital_period(a, Mbh)
    mu = constants.G * (Mbh + mstar)
    if MA < 0:
        MA = mean_anomaly(time, tp, P)

    EA = eccentric_anomaly(MA, ecc)  # eccentric anomaly
    # true anomaly in the correct quadrant
    ta = 2.0 * atan2(sqrt(1.0 + ecc) * sin(EA / 2.0), sqrt(1.0 - ecc) * cos(EA / 2.0))
    radius = a * (1.0 - ecc * cos(EA))  # radius from EA and ecc

    r = [] | units.AU  # Cartesian position
    r.append(
        radius * (cos(omra) * cos(omega + ta) - sin(omra) * sin(omega + ta) * cos(inc))
    )
    r.append(
        radius * (sin(omra) * cos(omega + ta) + cos(omra) * sin(omega + ta) * cos(inc))
    )
    r.append(radius * (sin(inc) * sin(omega + ta)))

    h = (mu * a * (1.0 - ecc * ecc)).sqrt()
    pp = a * (1 - ecc * ecc)

    v = [] | units.kms  # Cartesian velocity
    v.append(
        r.x * h * ecc / radius / pp * sin(ta)
        - h
        / radius
        * (cos(omra) * sin(omega + ta) + sin(omra) * cos(omega + ta) * cos(inc))
    )
    v.append(
        r.y * h * ecc / radius / pp * sin(ta)
        - h
        / radius
        * (sin(omra) * sin(omega + ta) - cos(omra) * cos(omega + ta) * cos(inc))
    )
    v.append(
        r.z * h * ecc / radius / pp * sin(ta) + h / radius * sin(inc) * cos(omega + ta)
    )

    return r, v


def main(T, a, e, i, o, O, t, P, M, m):
    T = T | units.yr
    a = a | units.AU
    t = t | units.yr
    M = M | units.MSun
    m = m | units.MSun
    i *= pi / 180.0
    o *= pi / 180.0
    O *= pi / 180.0
    r, v = orbital_elements_to_pos_and_vel(T, a, e, i, o, O, t, M, m)
    print("r=", r.in_(units.AU), "v=", v.in_(units.kms))


def new_option_parser():
    from optparse import OptionParser

    result = OptionParser()
    # data for S2 from 2009ApJ...692.1075G
    result.add_option("-T", dest="T", type="float", default=0)
    result.add_option("-a", dest="a", type="float", default=1042.5)
    result.add_option("-e", dest="e", type="float", default=0.88)
    result.add_option("-i", dest="i", type="float", default=135.25)
    result.add_option("-o", dest="o", type="float", default=225.39)
    result.add_option("-O", dest="O", type="float", default=63.56)
    result.add_option("-t", dest="t", type="float", default=2002.32)
    result.add_option("-M", dest="M", type="float", default=4.45e6)
    result.add_option("-m", dest="m", type="float", default=19.5)
    return result


if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
