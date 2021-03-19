from amuse.lab import *
import sys
import numpy
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.protodisk import ProtoPlanetaryDisk
from move_bodies_to_stellar_position import move_bodies_to_stellar_position


def test_distribution(bodies):
    from determine_orbital_elements import calculate_orbital_elements_for_single_planet
    asteroids = bodies[bodies.type==b"debris"]
    star = bodies[bodies.type==b"star"]
    print(star)
    print(asteroids)
    p = Particles()
    p.add_particle(star)
    from amuse.ext.orbital_elements import orbital_elements_from_binary
    sma = [] | units.au
    ecc = []
    for ai in asteroids:
        p.add_particle(ai)
        M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(p, G=constants.G)
        sma.append(a)
        ecc.append(e)
        p.remove_particle(ai)

    """
    from matplotlib import pyplot
    pyplot.scatter(sma.value_in(units.au), ecc)
    pyplot.semilogx()
    pyplot.show()
    """

def add_scattered_asteroids(star, Ndisk_per_MSun, mdisk,
                            qmin, qmax, emin, emax):

    Ndisk = int(Ndisk_per_MSun * star.mass.value_in(units.MSun))
    print("Run with Ndisk = ", Ndisk)
    converter = nbody_system.nbody_to_si(star.mass, qmin)
    e = emin + (emax-emin)*numpy.random.random(Ndisk)
    lqmin = numpy.log10(qmin.value_in(units.au))
    lqmax = numpy.log10(qmax.value_in(units.au))
    q = 10**numpy.random.uniform(lqmin, lqmax, Ndisk)  |units.au
    a = q/(1-e)
    asteroids = Particles(0)
    for i in range(Ndisk):
        ang_1, ang_2, ang_3  = numpy.random.uniform(0, 2*numpy.pi, 3)
        b = new_binary_from_orbital_elements(
            mass1=star.mass,
            mass2=mdisk,
            semimajor_axis=a[i],
            eccentricity=e[i],
            true_anomaly=ang_1 | units.rad,
            inclination=0 | units.rad,
            longitude_of_the_ascending_node=ang_2 | units.rad,
            argument_of_periapsis=ang_3 | units.rad,
            G=constants.G)
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

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = None,
                      help="output filename [%default]")
    result.add_option("--mdisk", unit=units.kg,
                      dest="mdisk", type="float",
                      default = 0.0 | units.kg,
                      help="mass of individual disk particles [%default]")
    result.add_option("--Ndisk", 
                      dest="Ndisk", type="int",
                      default = 100,
                      help="number of disk particles [%default]")
    result.add_option("--qmin", unit=units.au,
                      dest="qmin", type="float",
                      default = -5.2044 | units.au,
                      help="minimal closest appraoch [%default]")
    result.add_option("--qmax", unit=units.au,
                      dest="qmax", type="float",
                      default =  -30.33| units.au,
                      help="maximal closest appraoch [%default]")
    result.add_option("--emin", 
                      dest="emin", type="float",
                      default = 0.10,
                      help="minimum eccentricity [%default]")
    result.add_option("--emax", 
                      dest="emax", type="float",
                      default = 0.90,
                      help="maximum eccentricity [%default]")
    result.add_option("--seed", 
                      dest="seed", type="int",default = -1,
                      help="random number seed [%default]")
    return result

def get_inner_planet(planets):
    amin = planets.semimajor_axis.min()
    inner_planet = planets[planets.semimajor_axis<=amin]
    return inner_planet[0]

def get_outer_planet(planets):
    amax = planets.semimajor_axis.max()
    outer_planet = planets[planets.semimajor_axis>=amax]
    return outer_planet[0]

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    if o.seed>0:
        numpy.random.seed(o.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    print(bodies)
    stars = bodies[bodies.type==b'star']
    planets = bodies[bodies.type==b'planet']
    if o.qmin<0|units.au:
        inner_planet = get_inner_planet(planets)
        qmin = 0.5*inner_planet.semimajor_axis*(1.0-inner_planet.eccentricity)
    else:
        qmin = o.qmin
    if o.qmax<0|units.au:
        outer_planet = get_outer_planet(planets)
        qmax = 2*outer_planet.semimajor_axis*(1.0-outer_planet.eccentricity)
    else:
        qmax = o.qmax
    if len(stars)==0:
        stars = bodies[bodies.name==b"star"]

    debris = add_scattered_asteroids(stars[0], o.Ndisk, o.mdisk,
                                     qmin, qmax, o.emin, o.emax)
    bodies.add_particles(debris)
    #bodies.move_to_center()
    move_bodies_to_stellar_position(bodies)

    """
    test_distribution(bodies)
    from matplotlib import pyplot
    pyplot.scatter(bodies.x.value_in(units.parsec),
                   bodies.y.value_in(units.parsec))
    pyplot.show()
    """

    time = 0|units.Myr
    if o.outfile==None:
        filename = "added_scattered_asteroids.amuse"
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")
