from amuse.lab import *
import sys

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = "star.csv",
                      help="output filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    print(bodies)
    
    sys.stdout = open(o.outfile, 'w')
    print("# body key, name, type, mass (Msun), position (au), velocity (kms)")

    for bi in bodies:
        print(bi.key, bi.name, bi.type, bi.mass.value_in(units.MSun),\
                   bi.x.value_in(units.au), bi.y.value_in(units.au),\
                   bi.z.value_in(units.au), bi.vx.value_in(units.kms),\
                   bi.vy.value_in(units.kms), bi.vz.value_in(units.kms))
        
