from amuse.lab import *
from matplotlib import pyplot

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    bodies = read_set_from_file(o.filename, "amuse", close_file=True)
    pyplot.scatter(bodies.x.value_in(units.au),
                   bodies.y.value_in(units.au))
    pyplot.xlabel("x [au]")
    pyplot.ylabel("y [au]")
    pyplot.show()
