from amuse.lab import *

def move_bodies_to_stellar_position(bodies):
    star = bodies[bodies.type==b'star']
    #print "pos=", star.position.in_(units.parsec)
    com = star.center_of_mass()
    comv = star.center_of_mass_velocity()
    bodies.move_to_center()
    bodies.position +=  com
    bodies.velocity +=  comv
    #print "pos=", bodies.position.in_(units.parsec)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "star.amuse",
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = None,
                      help="output filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)

    move_bodies_to_stellar_position(bodies)
    #print bodies.velocity
    
    if o.outfile==None:
        filename = "moved_bodies.amuse"
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      append_to_file=False,
                      version="2.0")

    
