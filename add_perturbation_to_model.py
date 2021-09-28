from amuse.lab import *
import sys
import numpy

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile", default = None,
                      help="output filename [%default]")
    result.add_option("--perturbation", 
                      dest="perturbation", type="float",
                      default = -8.0,
                      help="log of the perturbation introduced [%default]")
    result.add_option("--parameter", 
                      dest="parameter", 
                      default = "x",
                      help="parameter to be perturbed [%default]")
    result.add_option("--key", 
                      dest="key", type="int",
                      default = -0.1,
                      help="key of the star to be perturbed [%default]")
    result.add_option("--name", 
                      dest="name", 
                      default = "",
                      help="name the perturbed star [%default]")
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

    bodies = read_set_from_file(o.filename, "amuse", close_file=True)
    if o.key>0:
        perturber = bodies[bodies.key==o.key]
    else:
        perturber = bodies.random_sample(1)

    if o.parameter=="x":
        unit = perturber.x.unit
        perturber.x += 10**(o.perturbation) | unit
    elif o.parameter=="y":
        unit = perturber.y.unit
        perturber.y += 10**(o.perturbation) | unit
    elif o.parameter=="z":
        unit = perturber.z.unit
        perturber.z += 10**(o.perturbation) | unit
    elif o.parameter=="vx":
        unit = perturber.vx.unit
        perturber.vx += 10**(o.perturbation) | unit
    elif o.parameter=="vy":
        unit = perturber.vy.unit
        perturber.vy += 10**(o.perturbation) | unit
    elif o.parameter=="vz":
        unit = perturber.vz.unit
        perturber.vz += 10**(o.perturbation) | unit
    elif o.parameter=="mass":
        unit = perturber.mass.unit
        perturber.mass += 10**(o.perturbation) | unit
        
    if len(o.name)>0:
        perturber.name = o.name
    
    time = 0|units.Myr
    if o.outfile==None:
        namelist = o.filename.split(".")
        namelist[-2] = namelist[-2] + ".pert_"+o.parameter+str(int(o.perturbation))
        filename = namelist[0]
        for ni in namelist[1:]:
            filename += "." + ni
    else:
        filename = o.outfile
        
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")
