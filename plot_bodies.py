from amuse.lab import *
from matplotlib import pyplot

def  plot_bodies(bodies):
    
    stars = bodies[bodies.type=="star"]
    planets = bodies[bodies.type=="planet"]
    debris = bodies[bodies.type=="debris"]
    pyplot.scatter(bodies.x.value_in(units.au),
                   bodies.y.value_in(units.au), s=0.1, c='k')
    pyplot.scatter(stars.x.value_in(units.au),
                   stars.y.value_in(units.au), s=30, c='y')
    pyplot.scatter(planets.x.value_in(units.au),
                   planets.y.value_in(units.au), s=10, c='b')
    pyplot.scatter(debris.x.value_in(units.au),
                   debris.y.value_in(units.au), s=1, c='k')
    pyplot.xlabel("x [au]")
    pyplot.ylabel("y [au]")
    pyplot.show()

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
    print(bodies)
    plot_bodies(bodies)
