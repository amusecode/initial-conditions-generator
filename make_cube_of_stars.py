from amuse.lab import *
import sys
import numpy

def make_cube_of_stars(masses, names, converter):
  
    N = len(masses)
    stars=Particles(N, mass=masses)
    L = converter.to_si(1|nbody_system.length)
    dv = converter.to_si(40|nbody_system.speed)
  
    stars.x=L*numpy.random.uniform(-1.,1.,N)
    stars.y=L*numpy.random.uniform(-1.,1.,N)
    stars.z=L*0.
    stars.vx=dv*numpy.random.uniform(-1.,1.,N)
    stars.vy=dv*numpy.random.uniform(-1.,1.,N)
    stars.vz=dv*0.

    stars.radius=(1.|units.RSun)*(stars.mass/(1.|units.MSun))**(1./3.)

    stars.name = names
    stars.mass = masses
    stars.move_to_center()
    return stars

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile", default = None,
                      help="output filename [%default]")
    result.add_option("--nstars", 
                      dest="nstars", type="int",
                      default = 100,
                      help="number of stars [%default]")
    result.add_option("--mass", unit=units.MSun,
                      dest="mass", type="float",
                      default = 1|units.MSun,
                      help="total cluster mass [%default]")
    result.add_option("--mmin", unit=units.MSun,
                      dest="mmin", type="float",
                      default = -0.1|units.MSun,
                      help="minimum stellar mass [%default]")
    result.add_option("--mmax", unit=units.MSun,
                      dest="mmax", type="float",
                      default = 100|units.MSun,
                      help="maximum stellar mass [%default]")
    result.add_option("-Q", 
                      dest="Qvir", type="float",
                      default = 0.5,
                      help="total mass [%default]")
    result.add_option("--radius", unit=units.parsec,
                      dest="radius", type="float",
                      default = 1.0|units.parsec,
                      help="cluster radius [%default]")
    result.add_option("--name", 
                      dest="name", 
                      default = "star",
                      help="stellar name [%default]")
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

    if o.mmin>0|units.MSun:
        masses = new_salpeter_mass_distribution(o.nstars, o.mmin, o.mmax)
    else:
        masses = numpy.zeros(o.nstars) | units.MSun
        masses += o.mass

    converter = nbody_system.nbody_to_si(masses.sum(), o.radius)
    bodies = make_cube_of_stars(masses, o.name, converter)
    print(bodies)

    bodies.scale_to_standard(convert_nbody=converter, virial_ratio=o.Qvir)
    index = 0
    time = 0|units.Myr
    if o.outfile==None:
        filename = "cube.amuse".format(index)
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")
