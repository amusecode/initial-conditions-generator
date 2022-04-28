from amuse.lab import *
import numpy
import logging

def True_anomaly_from_mean_anomaly(Ma, e):
    from math import sin
    Ta = Ma + (2*e - e**3/4.)*sin(Ma)\
        + (5*e**2/4.)*sin(2*Ma)\
        + (13./12.)*e**3*sin(3*Ma)
    return Ta

def semimajor_axis_to_orbital_period(a, Mtot) :
    return 2*numpy.pi * (a**3/(constants.G*Mtot)).sqrt()

def ZAMS_radius(mass):
    log_mass = numpy.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun))**2
    alpha = 0.08353 + 0.0565*log_mass
    beta  = 0.01291 + 0.2226*log_mass
    gamma = 0.1151 + 0.06267*log_mass
    r_zams = pow(mass.value_in(units.MSun), 1.25) * (0.1148 + 0.8604*mass_sq) / (0.04651 + mass_sq)

    return r_zams | units.RSun

def add_secondary(parent_stars, companion_name, masses,
                   semimajor_axis, eccentricity,
                  inclination, mean_anomaly,
                  LoAn,
                  Aop, ctype="star"):

    print("m=", masses)
    for i in range(len(parent_stars)):
        bi = parent_stars[i]
        binary_particle = Particle()
        binary_particle.position = bi.position
        binary_particle.velocity = bi.velocity
        binary_particle.type = ctype
        binary_particle.name = bi.name
        binary_particle.host = None
        binary_particle.type = "center_of_mass"

        mp = bi.mass
        ms = masses[i]

        from amuse.ext.orbital_elements import new_binary_from_orbital_elements
        Ta = True_anomaly_from_mean_anomaly(numpy.deg2rad(mean_anomaly),
                                            eccentricity)
        nb = new_binary_from_orbital_elements(mp, ms,
                                              semimajor_axis,
                                              eccentricity,
                                              Ta, inclination,
                                              LoAn, Aop,
                                              G=constants.G)
        

        nb.position += binary_particle.position
        nb.velocity += binary_particle.velocity
        nb[0].type = bi.type
        #nb[0].host = binary_particle
        nb[0].name = bi.name
        nb[1].type = ctype
        #nb[1].host = binary_particle
        nb[1].host = nb[0].name
        nb[1].name = companion_name
        nb[1].radius = ZAMS_radius(nb[1].mass)
    return nb

    #    binary_particle.child1 = nb[0]
    #    binary_particle.child2 = nb[1]
    #    binary_particle.semi_major_axis = semimajor_axis
    #    binary_particle.eccentricity = eccentricity
    #return binary_particle

def calculate_orbital_elementss(bi, converter):
    kep = new_kepler(converter)
    comp1 = bi.child1
    comp2 = bi.child2
    mass = (comp1.mass + comp2.mass)
    pos = (comp2.position - comp1.position)
    vel = (comp2.velocity - comp1.velocity)
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()
    return a, e

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "Plummer.amuse",
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile", default = None,
                      help="output filename [%default]")
    result.add_option("--name", 
                      dest="name", default = "ABAuriga",
                      help="name of the star that recieves companion [%default]")
    result.add_option("--cname", 
                      dest="companion_name", default = "ABAur_b",
                      help="name of the companion [%default]")
    result.add_option("--ctype", 
                      dest="ctype", default = "star",
                      help="type of the companion [%default]")
    result.add_option("-m", unit=units.MSun,
                      type = "float",
                      dest="mass", default = 0|units.MSun,
                      help="binary companion mass [%default]")
    result.add_option("-q", 
                      type = "float",
                      dest="q_companion", default = -1,
                      help="binary companion mass ratio [%default]")
    result.add_option("-a", unit=units.au,
                      type = "float",
                      dest="semimajor_axis", default = 1.0|units.au,
                      help="binary separation [%default]")
    result.add_option("-e", 
                      type = "float",
                      dest="eccentricity", default = 0.0,
                      help="binary eccenticity [%default]")
    result.add_option("-i", 
                      type = "float",
                      dest="inclination", default = 0.0,
                      help="inclination [%default]")
    result.add_option("--ma", 
                      type = "float",
                      dest="mean_anomaly", default = 0.0,
                      help="mean anomaly[%default]")
    result.add_option("--LoAn", 
                      type = "float",
                      dest="LoAn", default = 0.0,
                      help="Longitude of the ascending node [%default]")
    result.add_option("--Aop", 
                      type = "float",
                      dest="Aop", default = 0.0,
                      help="Argument of pericenter, [%default]")
    result.add_option("--seed", 
                      dest="seed", type="int", default = None,
                      help="random number seed [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    outfile = o.outfile
    
    if o.outfile == None:
        outfile = "stars_with_companion.amuse"

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    selected_stars = bodies[bodies.name==o.name]

    if o.q_companion>0:
        mass = selected_stars.mass * o.q_companion
    else:
        mass = numpy.zeros(len(selected_stars)) | units.MSun
        mass += o.mass
    stars = add_secondary(selected_stars,
                           o.companion_name,
                           mass,
                           o.semimajor_axis, o.eccentricity,
                          o.inclination, o.mean_anomaly,
                          o.LoAn, o.Aop, o.ctype)
    print(stars[1])
    bodies.add_particle(stars[1].as_set())
    print(bodies)
    time = 0 | units.Myr
    write_set_to_file(bodies,
                      outfile, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      overwrite_file=True,
                      version="2")
