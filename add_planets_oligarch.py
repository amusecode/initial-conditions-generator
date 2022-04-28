from amuse.lab import *
import sys
import numpy
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.protodisk import ProtoPlanetaryDisk
from make_planets_oligarch import make_planets_oligarch
from move_bodies_to_stellar_position import move_bodies_to_stellar_position

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "star.amuse",
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = None,
                      help="output filename [%default]")
    result.add_option("--rmin_disk", unit=units.au,
                      dest="rmin_disk", type="float",
                      default = 1.0 | units.au,
                      help="inner disk radius [%default]")
    result.add_option("--rmax_disk", unit=units.au,
                      dest="rmax_disk", type="float",
                      default = 100.0 | units.au,
                      help="outer disk radius [%default]")
    result.add_option("--relative_diskmass", 
                      dest="relative_diskmass", type="float",
                      default = 0.01, 
                      help="relative disk mass [%default]")
    result.add_option("--Nplanets", 
                      dest="Nplanets", type="int",
                      default = 4,
                      help="maximum number of planets [%default]")
    result.add_option("--fplanets", 
                      dest="fplanets", type="float",
                      default = 0.5,
                      help="fraction of stars with planets [%default]")
    result.add_option("--mmin", unit=units.MSun,
                      dest="mmin", type="float",
                      default = -0.1|units.MSun,
                      help="minimum stellar mass [%default]")
    result.add_option("--mmax", unit=units.MSun,
                      dest="mmax", type="float",
                      default = 100|units.MSun,
                      help="maximum stellar mass [%default]")
    
    result.add_option("--seed", 
                      dest="seed", type="int",default = -1,
                      help="random number seed [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    accurancy = 0.0001 
    planet_density =  3 | units.g/units.cm**3

    if o.seed>0:
        numpy.random.seed(o.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    stars = bodies[bodies.type=="star"]
    stars = stars[stars.mass>=o.mmin]
    stars = stars[stars.mass<=o.mmax]
    if o.fplanets<=1:
        nplanets = int(o.fplanets*len(stars))
    else:
        nplanets = min(len(stars), int(o.fplanets))
    stars_with_planets = stars.random_sample(nplanets)
    print("Number of stars with planets=", len(stars_with_planets))
    for star in stars_with_planets:
        disk_mass = o.relative_diskmass * star.mass
        planets = make_planets_oligarch(star.mass, star.radius,
                                        o.rmin_disk, o.rmax_disk, disk_mass)
        while len(planets)>o.Nplanets:
            mmin = planets.mass.min()
            lm_planet = planets[planets.mass<=mmin]
            planets.remove_particle(lm_planet)
            
        planets.type = "planet"
        planets.name = "planet"
        planets.host = star.key
        print(planets.position.in_(units.parsec))
        planets.position += star.position
        planets.velocity += star.velocity
        bodies.add_particles(planets)
    move_bodies_to_stellar_position(bodies)

    time = 0|units.Myr
    if o.outfile==None:
        filename = "added_oligarch_planets.amuse"
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")

    
