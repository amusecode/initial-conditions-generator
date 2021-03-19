import os
import sys
import time
import numpy
from amuse.lab import * 

from amuse.ic.solar_system_moons import new_lunar_system
from amuse.ic.solar_system_moons import new_lunar_system_in_time
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.plot import scatter
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.solarsystem import solar_system_in_time

def semimajor_axis(Mtot, P):
    return ((constants.G * Mtot*(P**2)/(4 * numpy.pi**2)))**(1./3.)
def orbital_period(Mtot, a):
    return (((4 * numpy.pi**2) * a**3)/(constants.G * Mtot)).sqrt()

def rescale_planetary_orbit(sun, planets, name, sma):

    planet = planets[planets.name==name]
    pos = planet.position
    vel = planet.velocity
    ang_1, ang_2, ang_3 = numpy.random.uniform(0, 2*numpy.pi, 3)
    ecc = 0.01*numpy.random.random()
    inc = numpy.random.random()
    binary = new_binary_from_orbital_elements(sun.mass,
                                              planet.mass,
                                              sma,
                                              ecc,
                                              inclination=inc,
                                              true_anomaly=ang_1,
                                              argument_of_periapsis=ang_2,
                                              longitude_of_the_ascending_node=ang_3,
                                              G=constants.G)
    binary.position -= binary[0].position
    binary.velocity -= binary[0].velocity
    planet.position = binary[1].position
    planet.velocity = binary[1].velocity
    planet.semimajor_axis = sma
    planet.eccentricity = ecc
    planet.inclination = inc
    planet.type = "planet"
    planet.name = name
    planet.hostname = sun.name
    return planet

names = ["Uranus", "Neptune", "Saturn", "Jupiter"]
def get_nice_model_conditions():

    converter = nbody_system.nbody_to_si(1|units.MSun, 1 | units.AU)

    time_0 = 2457099.5 | units.day
    time_JD=2457099.5|units.day
    delta_JD = time_JD-time_0
    solar_system = solar_system_in_time(time_JD)
    solar_system.hosttype = "star"
    solar_system[0].type = "star"
    solar_system[0].name = "sun"
    solar_system[1:].type = "planet"
    solar_system[1:].hostname = "sun"

    sun = solar_system[solar_system.name=="sun"]
    solar_system.position -= sun.position
    solar_system.velocity -= sun.velocity
    
    nice_system = Particles()
    nice_system.add_particle(sun)
    nice_system.name = "Sun"
    nice_system.type = "star"
    sma = 15 | units.AU
    for name in names:
        planet = rescale_planetary_orbit(sun, solar_system, name, sma)
        nice_system.add_particles(planet)
        porb = 1./2. * orbital_period(sun.mass, sma)
        sma = semimajor_axis(sun.mass, porb)

    nice_system.move_to_center()
    return nice_system

def get_sun_jupiter_and_moons(Ndisk):
    converter = nbody_system.nbody_to_si(1|units.MSun, 1 | units.AU)
    solar_system = new_solar_system()
    pebbels = Particles(0)
    ss = Particles(0)
    ss.add_particle(solar_system[0])
    ss.add_particles(solar_system[solar_system.name=="Jupiter"])
    ss.add_particles(solar_system[solar_system.hostname=="Jupiter"])
    return ss, pebbels, converter


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outputfilename", default = "nice.amuse",
                      help="outpute filename [%default]")
    return result


if __name__ in ('__main__','__plot__'):
    o, r = new_option_parser().parse_args()

    nice_model = get_nice_model_conditions()
    print(nice_model)

    time = 0|units.Myr
    write_set_to_file(nice_model,
                      o.outputfilename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")

    

