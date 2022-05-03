from amuse.lab import *
import sys
import numpy
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.couple import bridge

def add_debris_disk(star, ndisk_per_MSun, fdisk, masteroid,
                    Rmin, Rmax, alpha=1.5, Q=1.0):

    Ndisk = int(ndisk_per_MSun * star.mass.value_in(units.MSun))
    print("Run with Ndisk = ", Ndisk)

    converter = nbody_system.nbody_to_si(star.mass, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk,
                              convert_nbody=converter,
                              densitypower=alpha,
                              Rmin=Rmin/Rmin,
                              Rmax=Rmax/Rmin,
                              q_out=Q,
                              discfraction=fdisk).result
    
    disk.mass = masteroid
    disk.name = "asteroid"
    disk.type = "debris"
    disk.host = star.key
    disk.position += star.position
    disk.velocity += star.velocity
    return disk

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = None,
                      help="output filename [%default]")
    result.add_option("--name", 
                      dest="name", 
                      default = "Sun",
                      help="disk mass [%default]")
    result.add_option("--fdisk", 
                      dest="fdisk", type="float",
                      default = -1,
                      help="disk mass fraction [%default]")
    result.add_option("--Mdisk", 
                      dest="Mdisk", type="float", unit=units.MSun,
                      default = 0.01|units.MSun,
                      help="disk mass [%default]")
    result.add_option("--mdisk", unit=units.kg,
                      dest="mdisk", type="float",
                      default = 0.0 | units.kg,
                      help="mass of individual disk particles [%default]")
    result.add_option("--ndisk", 
                      dest="ndisk", type="int",
                      default = 100,
                      help="number of disk particles per MSun [%default]")
    result.add_option("--rmin", unit=units.au,
                      dest="rmin", type="float",
                      default = 10 | units.au,
                      help="minimal disk radius [%default]")
    result.add_option("--rmax", unit=units.au,
                      dest="rmax", type="float",
                      default = 100 | units.au,
                      help="maximal disk radius [%default]")
    result.add_option("-Q", 
                      dest="Q", type="float",
                      default = 1.0,
                      help="Toomre Q parameter [%default]")
    result.add_option("-a", 
                      dest="alpha", type="float",
                      default = 1.5,
                      help="Disk density profile [%default]")
    result.add_option("--seed", 
                      dest="seed", type="int",default = -1,
                      help="random number seed [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    if o.seed>0:
        numpy.random.seed(o.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    stars = bodies[bodies.name==o.name]
    if len(stars)==0:
        stars = bodies[bodies.name=="star"]
    if len(stars)==0:
        stars = bodies[bodies.name=="SUN"]
    if len(stars)==0:
        stars = bodies[bodies.type=="star"]

        
    for star in stars:
        if o.fdisk<0:
            fdisk = o.Mdisk/star.mass
        else:
            fdisk = o.fdisk
        debris = add_debris_disk(star, o.ndisk, fdisk, o.mdisk,
                                 o.rmin, o.rmax, o.alpha, o.Q)
        bodies.add_particles(debris)
        bodies.move_to_center()

    time = 0|units.Myr
    if o.outfile==None:
        filename = "added_debris_disk.amuse"
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")
