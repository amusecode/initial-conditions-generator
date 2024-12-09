"""
Make planets based on the an assumption of oligarchical planetary growth.
(the formalism of Kokubo & Ida (1998))

After:

@article{0004-637X-807-2-157,
 author={Scott Tremaine},
 title={The Statistical Mechanics of Planet Orbits},
 journal={The Astrophysical Journal},
 volume={807},
 number={2},
 pages={157},
 url={http://stacks.iop.org/0004-637X/807/i=2/a=157},
 year={2015},
 abstract={The final phase of terrestrial planet formation is believed to begin with a large number of planetary on nearly circular, coplanar orbits. Mutual gravitational interactions gradually excite their eccentricities until their orbits cross and they collide and merge; through this process the number of surviving bodies declines until the system contains a small number of planets on well-separated, stable orbits. In this paper we explore a simple statistical model for the orbit distribution of planets formed by this process, based on the sheared-sheet approximation and the ansatz that the planets explore uniformly all of the stable region of phase space. The model provides analytic predictions for the distribution of eccentricities and semimajor axis differences, correlations between orbital elements of nearby planets, and the complete N -planet distribution function, in terms of a single parameter, the that is determined by the planetary masses. The predicted properties are generally consistent with N -body simulations of the giant-impact phase and with the distribution of semimajor axis differences in the Kepler catalog of extrasolar planets. A similar model may apply to the orbits of giant planets if these orbits are determined mainly by dynamical evolution after the planets have formed and the gas disk has disappeared.}
}

Refers to (with the actual method as implemented here):


@article{0004-637X-775-1-53,
  author={Brad M. S. Hansen and Norm Murray},
  title={Testing in Situ Assembly with the Kepler Planet Candidate Sample},
  journal={The Astrophysical Journal},
  volume={775},
  number={1},
  pages={53},
  url={http://stacks.iop.org/0004-637X/775/i=1/a=53},
  year={2013}
 }

Actual model:

@article{0004-637X-581-1-666,
  author={Eiichiro Kokubo and Shigeru Ida},
  title={Formation of Protoplanet Systems and Diversity of Planetary Systems},
  journal={The Astrophysical Journal},
  volume={581},
  number={1},
  pages={666},
  url={http://stacks.iop.org/0004-637X/581/i=1/a=666},
  year={2002},
}



"""
import numpy
from numpy import random
import time as pytime

import os.path


from amuse.units import quantities
from amuse.units import units
from amuse.units import nbody_system
from amuse.units.core import named_unit
from amuse.units.optparse import OptionParser
from amuse.datamodel import Particle, Particles
from amuse.community.kepler.interface import Kepler
from amuse.io import write_set_to_file
from amuse.lab import new_kroupa_mass_distribution


MEarth = named_unit("MEarth", "MEarth", 5.97219e+24 * units.kg)


def get_mass(a, surface_density_factor = 1.0, mstar = 1 | units.MSun):
    s = surface_density_factor * (10.0 | (units.g / units.cm**2))
    z = 2 * numpy.pi * a * s * ((a / (5.0|units.AU) ) ** -(3.0/2.0)) * 4.0 * a
    return (2.0 * z ** 3 / (3.0 * mstar)).sqrt()
    
def get_orbital_separation(mass, a, mstar = 1 | units.MSun):
    return 10.0 * a * (2*mass/(3.0 * mstar)) ** (1.0/3.0)


def get_mass_and_orbital_separation(a, surface_density_factor = 1.0, mstar = 1 | units.MSun):
    mass = get_mass(a, surface_density_factor = surface_density_factor, mstar = mstar)
    orbital_separation = get_orbital_separation(mass, a, mstar = mstar)
    return mass, orbital_separation

def new_distribution(x0, x1, surface_density_factor = 1.0, mstar = 1 | units.MSun):
    x = x0
    masses = []
    positions = []
    while x < x1:
        mass, orbital_separation = get_mass_and_orbital_separation(x, surface_density_factor = surface_density_factor, mstar = mstar)
        masses.append(mass)
        positions.append(x)
        x += orbital_separation
    
    return (
        quantities.as_vector_quantity(masses),
        quantities.as_vector_quantity(positions)
    )
    


def is_hit(mass, target_total_mass, accurancy):
    return abs((mass - target_total_mass) / target_total_mass) < accurancy

    

def new_planet_distribution(x0, x1,  target_total_mass, accurancy, max_iterations = 1000, mstar = 1 | units.MSun):
    min_surface_density_factor = 0.1
    max_surface_density_factor = 10000.0
    
    current_surface_density_factor = min_surface_density_factor
    mass, pos = new_distribution(x0,x1,current_surface_density_factor, mstar = mstar)
    total_mass = mass.sum()
    iteration = 0
    while not is_hit(total_mass,  target_total_mass, accurancy):
        #print iteration, total_mass / (1 | MEarth), current_surface_density_factor
        iteration += 1
        if iteration >= max_iterations:
            break
        
        if total_mass < target_total_mass:
            min_surface_density_factor = current_surface_density_factor
        elif total_mass > target_total_mass:
            max_surface_density_factor  = current_surface_density_factor
            
        current_surface_density_factor = (min_surface_density_factor + max_surface_density_factor) / 2.0
        
        mass, pos = new_distribution(x0,x1,current_surface_density_factor, mstar = mstar)
        total_mass = mass.sum()
    
    return mass, pos, current_surface_density_factor
        

    
    


def new_rotation_matrix_from_euler_angles(phi, theta, chi):
    cosp=numpy.cos(phi)
    sinp=numpy.sin(phi)
    cost=numpy.cos(theta)
    sint=numpy.sin(theta)
    cosc=numpy.cos(chi)
    sinc=numpy.sin(chi)
    #see wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix
    return numpy.array(
        [[cost*cosc, -cosp*sinc + sinp*sint*cosc, sinp*sinc + cosp*sint*cosc], 
         [cost*sinc, cosp*cosc + sinp*sint*sinc, -sinp*cosc + cosp*sint*sinc],
         [-sint,  sinp*cost,  cosp*cost]])



def rotate(position, velocity, phi, theta, psi): # theta and phi in radians
    Runit = position.unit
    Vunit = velocity.unit
    matrix = new_rotation_matrix_from_euler_angles(phi, theta, psi)
    return (numpy.dot(matrix, position.value_in(Runit)) | Runit,
           numpy.dot(matrix, velocity.value_in(Vunit)) | Vunit)

# select Eurler angles randomly. 

def random_Euler_angles():
    phi   = 2*numpy.pi*random()
    theta = numpy.acos(1-2*random())
    chi   = 2*numpy.pi*random()
    return phi, theta, chi


def posvel_from_orbital_elements(Mstar, semimajor_axis, eccentricity, kepler, rng = None):
    if rng is None:
        rng = random
        
    mean_anomaly = rng.uniform(0, 2*numpy.pi, 1)
    kepler.initialize_from_elements(
        Mstar, 
        semimajor_axis, 
        eccentricity, 
        mean_anomaly=mean_anomaly)
    position = quantities.as_vector_quantity(kepler.get_separation_vector())
    velocity = quantities.as_vector_quantity(kepler.get_velocity_vector())
    return position, velocity
    

def make_planets(central_particle, masses, radii, density = 3 | units.g/units.cm**3, phi=None, theta=None, eccentricity = 0.0, kepler = None, rng = None):
    volumes = masses / density
    planet_radii = (3.0 * volumes /  (4.0 * numpy.pi))**(1.0/3.0)
    n = len(masses)
    planet_particles = Particles(n)
    planet_particles.semimajor_axis = radii
    if eccentricity is None:
        eccentricity = numpy.abs(rng.normal(-0.00001,0.00001,n))
    planet_particles.eccentricity = eccentricity
    planet_particles.mass = masses
    planet_particles.radius = planet_radii


    if phi is None:
        phi = numpy.radians(rng.uniform(0.0, 90.0, 1)[0])#rotate under x
    if theta is None:
        theta0 = numpy.radians((rng.normal(-90.0,90.0,1)[0]))#rotate under y
        theta0 = 0
        theta_inclination = numpy.radians(rng.normal(0, 1.0, n )) 
        theta_inclination[0] = 0
        theta = theta0 + theta_inclination

    #psi = numpy.radians(rng.uniform(0, 180, 1))[0] #0 # numpy.radians(90) # numpy.radians(rng.uniform(0, 180, 1))[0]
    psi = numpy.radians(rng.uniform(0.0, 180.0, 1))[0] #0 # numpy.radians(90) # numpy.radians(rng.uniform(0, 180, 1))[0]
    com_particle = central_particle.copy()
    for x, t in zip(iter(planet_particles), theta):
        pos,vel = posvel_from_orbital_elements(com_particle.mass + x.mass, x.semimajor_axis, x.eccentricity, kepler, rng)
        pos,vel = rotate(pos, vel, 0, 0, psi) # theta and phi in radians            
        pos,vel = rotate(pos, vel, 0, t, 0) # theta and phi in radians            
        pos,vel = rotate(pos, vel, phi, 0, 0) # theta and phi in radians            
        x.position = pos + com_particle.position
        x.velocity = vel + com_particle.velocity
        if False:
            two_body = Particles(particles=[com_particle, x])
            print("dp:", (com_particle.position - two_body.center_of_mass()).as_quantity_in(units.AU))
            com_particle.mass = two_body.mass.sum()
            com_particle.position = two_body.center_of_mass()
            com_particle.velocity = two_body.center_of_mass_velocity()

    #planet_particles.position += central_particle.position
    #planet_particles.velocity += central_particle.velocity

    return planet_particles

def new_system(
        star_mass = 1|units.MSun, 
        star_radius = 1|units.RSun, 
        disk_minimum_radius = 0.05 | units.AU,
        disk_maximum_radius = 10 | units.AU,
        disk_mass = 20 | MEarth,
        accurancy = 0.0001, 
        planet_density =  3 | units.g/units.cm**3,
        rng = None,
        kepler = None):
            
    central_particle = Particle()
    central_particle.mass =  star_mass
    central_particle.position = (0,0,0) | units.AU
    central_particle.velocity = (0,0,0) | units.kms
    central_particle.radius = star_radius
    central_particle.name = "star"
    central_particle.type = "star"
    central_particle.id = 0
    
    if rng is None:
        rng = numpy.random
        
    converter = nbody_system.nbody_to_si(1|units.MSun, 1 | units.AU)
    
    if kepler is None:
        kepler = Kepler(converter)
        kepler.initialize_code()
        
    m, r, f = new_planet_distribution(
        disk_minimum_radius, disk_maximum_radius, 
        disk_mass,
        accurancy
    )
    
    planets = make_planets(
        central_particle, 
        m, r, 
        density = planet_density, 
        phi = 0, theta = None, 
        kepler = kepler, 
        rng = rng
    )
    planets.name = "planet"
    planets.type = "planet"
    for i in range(len(planets)):
        planets[i].id = i
    
    central_particle.planets = planets
    kepler.stop()
    p = Particles()
    p.add_particle(central_particle)
    return p
    

def new_option_parser():
    result = OptionParser()
    result.add_option("--seed", 
                      dest="seed", type="int", default = -1,
                      help="random number seed [%default]")
    result.add_option("-F", "--output",
                      dest="output_filename", default = "star.h5",
                      help="name of the output filename [%default]")
    result.add_option("-m", "--Mstar", unit=units.MSun, type="float", default = 1|units.MSun,
                      dest="star_mass", help='mass of the star [%default]')
    result.add_option("-r", "--Rstar", unit=units.RSun, type="float", default = -1|units.RSun,
                      dest="star_radius", help='radius of the star [%default]')
    result.add_option("--amin", unit=units.AU, type="float", default = 0.3|units.AU,
                      dest="disk_minimum_radius", help='minimum radius of the disk [%default]')
    result.add_option("--amax", unit=units.AU, type="float", default = 300|units.AU,
                      dest="disk_maximum_radius", help='maximum radius of the disk [%default]')
    result.add_option("--Nstars", type="int", default = 1,
                      dest="Nstars", help='number of stars [%default]')
    result.add_option("--Mmin", unit=units.MSun, type="float", default = 0.3|units.MSun,
                      dest="minimum_stellar_mass", help='minimum mass of primary star [%default]')
    result.add_option("--Mmax", unit=units.MSun, type="float", default = 3.0|units.MSun,
                      dest="maximum_stellar_mass", help='maximum mass of primary star [%default]')
    result.add_option("--Mdisk", unit=MEarth, type="float", default = -100|MEarth,
                      dest="disk_mass", help='mass of the disk [%default]')
    result.add_option("--fdisk", type="float", default = 0.01,
                      dest="disk_mass_fraction", help='fraction of the mass of the disk [%default]')
    result.add_option("--planet-density", unit=units.g/units.cm**3, type="float", default = 3 | units.g/units.cm**3,
                      dest="planet_density", help='density of the planets  (used for planet radii calculation) [%default]')
    result.add_option("--accuracy", type="float", default = 0.0001,
                      dest="accuracy", help='how accurate should the disk mass be integrated to[%default]')
    #result.add_option("--plot",
    #              action="store_true", dest="make_plot", default=False,
    #              help="make a plot of the planets")
    
    return result
    
def make_planets_oligarch(
        star_mass = 1|units.MSun, 
        star_radius = 1|units.RSun, 
        disk_minimum_radius = 0.05 | units.AU,
        disk_maximum_radius = 10 | units.AU,
        disk_mass = 200 | MEarth,
        accuracy = 0.0001, 
        planet_density =  3 | units.g/units.cm**3,
        seed = -1):
            
    if seed < 0:
        rng = random.RandomState()
    else:
        rng = random.RandomState(seed)
    
    
    output = new_system(
        star_mass, 
        star_radius, 
        disk_minimum_radius,
        disk_maximum_radius,
        disk_mass,
        accuracy, 
        planet_density,
        rng = rng)
    
    star = output[0].as_set()
    print("Stellar mass:", star_mass.in_(units.MSun))
    print("Number of planets generated:", len(star.planets))
    #print("Total mass:", output.planets.mass.sum().as_quantity_in(MEarth))
    for i, planet  in enumerate(star.planets[0]):
        print("Planet: {0: 3d} , mass: {1: 8.3f}  MEarth, a: {2: 8.2f} au".format(i, planet.mass.value_in(MEarth), planet.semimajor_axis.value_in(units.au)))

    star.add_particles(output.planets[0])
    return star

def ZAMS_radius(mass):
    log_mass = numpy.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun))**2
    alpha = 0.08353 + 0.0565*log_mass
    beta  = 0.01291 + 0.2226*log_mass
    gamma = 0.1151 + 0.06267*log_mass
    r_zams = pow(mass.value_in(units.MSun), 1.25) * (0.1148 + 0.8604*mass_sq) / (0.04651 + mass_sq)

    return r_zams | units.RSun

    
if __name__ == "__main__":
    o, arguments  = new_option_parser().parse_args()

    if o.seed>0:
        numpy.random.seed(o.seed)
    else:
        print("random number seed from clock.")

    if o.Nstars==1:
        if o.disk_mass>0|units.MSun:
            disk_mass = o.disk_mass
        else:
            disk_mass = o.disk_mass_fraction*o.star_mass
            
        if o.star_radius<=0|units.RSun:
            star_radius = ZAMS_radius(o.star_mass)
        else:
            star_radius = o.star_radius
        star_with_planets = make_planets_oligarch(star_mass = o.star_mass,
                                                  star_radius = star_radius,
                                                  disk_minimum_radius = o.disk_minimum_radius,
                                                  disk_maximum_radius = o.disk_maximum_radius,
                                                  disk_mass = disk_mass,
                                                  accuracy = o.accuracy,
                                                  planet_density =  o.planet_density,
                                                  seed = o.seed)
        #output = make_planets_oligarch(**new_option_parser().parse_args()[0].__dict__)
        write_set_to_file(star_with_planets,
                          o.output_filename,
                          'hdf5',
                          version="2.0",
                          append_to_file = False)
    else:
        masses = new_kroupa_mass_distribution(o.Nstars,
                                              mass_min=o.minimum_stellar_mass,
                                              mass_max=o.maximum_stellar_mass)

        for mi in range(len(masses)):
            output_filename = 'StarWithPlanets_i%6.6i.amuse'%mi
            if o.disk_mass>0|units.MSun:
                disk_mass = o.disk_mass
            else:
                disk_mass = o.disk_mass_fraction*masses[mi]
            
            star_radius = ZAMS_radius(masses[mi])
            star_with_planets = make_planets_oligarch(star_mass = masses[mi],
                                                      star_radius = star_radius,
                                                      disk_minimum_radius = o.disk_minimum_radius,
                                                      disk_maximum_radius = o.disk_maximum_radius,
                                                      disk_mass = disk_mass,
                                                      accuracy = o.accuracy,
                                                      planet_density =  o.planet_density,
                                                      seed = o.seed)
            write_set_to_file(star_with_planets,
                              output_filename,
                              'hdf5',
                              version="2.0",
                              append_to_file = False)

