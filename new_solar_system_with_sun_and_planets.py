import numpy
import os
import requests
from amuse.lab import *
from math import pi, sin, cos, atan2, sqrt

def eccentric_anomaly(mean_anomaly, e) :
    ecc_anomaly = mean_anomaly + 2*pi*(e * sin(mean_anomaly) + 0.5*e*e*sin(2*mean_anomaly))
    m = ecc_anomaly - e*sin(ecc_anomaly)
    de = (mean_anomaly-m) / (1 - e*cos(ecc_anomaly))
    ecc_anomaly += de;
    while de >= 0.001 :
        m = ecc_anomaly - e*sin(ecc_anomaly)
        de = (mean_anomaly-m) / (1 - e*cos(ecc_anomaly))
        ecc_anomaly += de
    return ecc_anomaly

def True_anomaly_from_mean_anomaly(Ma, e):
    from math import sin
    Ta = Ma + (2*e - e**3/4.)*sin(Ma)\
        + (5*e**2/4.)*sin(2*Ma)\
        + (13./12.)*e**3*sin(3*Ma)
    return Ta

def construct_particle_set_from_orbital_elements(name,
                                                 mass,
                                                 a, ecc, inc, Ma, Aop, LoAn,
                                                 parent):
        
    #print("length:", len(a), len(ecc), len(inc), len(Ma), len(Aop), len(LoAn), len(name))

    p = Particles(len(a))
    p.name = name
    p.type = "planet"
    p.host = parent.name
    Earth = p[p.name=="EM"]
    Earth.name = "EarthMoon"
    
    from orbital_elements_to_cartesian import orbital_elements_to_pos_and_vel
    from orbital_elements_to_cartesian import orbital_period

    from amuse.ext.orbital_elements import new_binary_from_orbital_elements
    for i in range(len(p)):
        p[i].mass = mass[i]
        Ta = True_anomaly_from_mean_anomaly(numpy.deg2rad(Ma[i]), ecc[i])
        #print("True Anomaly:", Ta)
        b = new_binary_from_orbital_elements(parent.mass, p[i].mass, a[i], ecc[i], Ta, inc[i], LoAn[i], Aop[i], G=constants.G)
        p[i].position = b[1].position - b[0].position
        p[i].velocity = b[1].velocity - b[0].velocity
        rho = 1.0 | units.g/units.cm**3
        p.radius = (p.mass/rho)**(1./3.) 

    return p

def read_orbital_elements_for_planets(n, filename="MPCORB.DAT",
                                      Julian_date=2451545.0|units.day):
    #f = open(fdir+"MPCORB.DAT", "r")
    planet_names = ["Mercury", "Venus", "EM Bary", "Mars",
                    "Jupiter", "Saturn", "Uranus", "Neptune",
                    "Pluto"]
    planet_mass = [3.302e23,
                   48.685e23,
                   5.97219e24,
                   6.4185e23,
                   1898.13e24,
                   5.68319e26,
                   86.8103e24,
                   102.41e24,
                   1.309e+22] | units.kg

    sma = [] | units.AU
    ecc = []
    inc = [] # degrees
    Ma = [] # degrees
    Aop = [] # degrees
    LoAn = [] # degrees
    name = []
    for line in open(filename):
        #see: https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
        for pi in planet_names:
            if pi in line:
                l = line.split(" ")
                name.append(l[0])
                sma.append(float(line[8:20]) | units.au)
                ecc.append(float(line[21:36]))
                inc.append(float(line[37:52]))
                L = float(line[53:70])
                Lperi = float(line[71:86])
                Lnode = float(line[87:102])
                if n>0 and len(sma)>=n:
                    break
                Aop.append(Lperi)
                LoAn.append(Lnode)
                T = (Julian_date-(2451545.0|units.day))/(36525|units.day)
                omega = Lperi
                M = L - omega
                Ma.append(M)
            
    return name, planet_mass, sma, ecc, inc, Ma, Aop, LoAn

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "MPCORB.h5",
                      help="read to file [%default]")
    result.add_option("-F", 
                      dest="outfilename", default = None,
                      help="write to file [%default]")
    result.add_option("-d",  unit=units.day,
                      type = "float",
                      dest="Julian_date", default = 2451545.0|units.day,
                      help="Julan date [%default]")
    result.add_option("-n", 
                      dest="n", type="int", default = 0,
                      help="number of planetesimals [%default]")
    return result

if __name__ in ('__main__', '__plot__'):

    o, arguments  = new_option_parser().parse_args()

    sun = Particles(1)
    sun.Julian_date = o.Julian_date
    sun.mass = 1 | units.MSun
    sun.radius = 1 |units.RSun
    sun.position = (0,0,0) | units.km
    sun.velocity = (0,0,0) | units.kms
    sun.type = "star"
    sun.name = "Sun"
    sun.host = "None"
    
    if not os.path.isfile("p_elem_t1.txt"):
        print("download p_elem_t1.txt")
        url = 'https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt'
        datafile = requests.get(url)
        open('p_elem_t1.txt', 'wb').write(datafile.content)

    name, mass, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_for_planets(o.n, filename="p_elem_t1.txt", Julian_date=o.Julian_date)
    p = construct_particle_set_from_orbital_elements(name, mass,
                                                     a, ecc, inc, Ma, Aop, LoAn, sun[0])

    sun.add_particles(p)
    sun.move_to_center()

    solar_system = sun
    
    print("number of particles N=", len(sun))
    
    if o.outfilename:
        write_set_to_file(sun, o.outfilename, "amuse",
                          append_to_file=False,
                          overwrite_file=True,
                          version="2.0")
    else:
        from amuse.plot import scatter
        from matplotlib import pyplot
        star = solar_system[solar_system.type=="star"]
        planet = solar_system[solar_system.type=="planet"]
        pyplot.scatter(star.x.value_in(units.au), star.y.value_in(units.au), s=100, c='y')
        pyplot.scatter(planet.x.value_in(units.au), planet.y.value_in(units.au), s=30, c='b')
        pyplot.show()
        print(planet)
        
