import numpy
import os
import requests
from math import pi, sin, cos, atan2, sqrt

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles, Particle
from amuse.community.kepler.interface import Kepler
from amuse.lab import read_set_from_file, write_set_to_file
from amuse.ext.solarsystem import solar_system_in_time

def new_kepler():
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  kepler = Kepler(converter)
  kepler.initialize_code()
  return kepler

def get_position(mass_sun, mass_planet, ecc, semi, mean_anomaly, incl, argument, longitude, delta_t=0.|units.day):
  """
  cartesian position and velocity from orbital elements,
  where the orbit is evolved from given mean_anomaly 
  by time delta_t
  argument -- argument of perihelion
  longitude -- longitude of ascending node
  """
  kepler = new_kepler()
  kepler.initialize_from_elements(mass=(mass_sun+mass_planet),
                                  semi=semi,
                                  ecc=ecc,
                                  mean_anomaly=mean_anomaly)
  kepler.transform_to_time(time=delta_t)
  r = kepler.get_separation_vector()
  v = kepler.get_velocity_vector()
  
  kepler.stop()
  
  a1 = ([numpy.cos(longitude), -numpy.sin(longitude), 0.0], [numpy.sin(longitude), numpy.cos(longitude), 0.0], [0.0, 0.0, 1.0])
  a2 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
  a3 = ([numpy.cos(argument), -numpy.sin(argument), 0.0], [numpy.sin(argument), numpy.cos(argument), 0.0], [0.0, 0.0, 1.0])
  A = numpy.dot(numpy.dot(a1,a2),a3)
  print(A, r)
  
  # old version from P2.7
  #  r_vec = numpy.dot(A,numpy.reshape(r,3,1))
  #  v_vec = numpy.dot(A,numpy.reshape(v,3,1))
  r_vec = numpy.dot(A,numpy.reshape(r,3,'F'))
  v_vec = numpy.dot(A,numpy.reshape(v,3,'F'))
  
  # for relative vectors
  r[0] = r_vec[0]
  r[1] = r_vec[1]
  r[2] = r_vec[2]
  v[0] = v_vec[0]
  v[1] = v_vec[1]
  v[2] = v_vec[2]
  
  return r,v

def True_anomaly_from_mean_anomaly(Ma, e):
    from math import sin
    Ta = Ma + (2*e - e**3/4.)*sin(Ma)\
        + (5*e**2/4.)*sin(2*Ma)\
        + (13./12.)*e**3*sin(3*Ma)
    return Ta

def construct_particle_set_from_orbital_elements(name, mass,
                                                 a, ecc, inc, Ma, Aop, LoAn,
                                                 Mparent = 1|units.MSun):

    print("length:", len(a), len(ecc), len(inc), len(Ma), len(Aop), len(LoAn), len(name))

    p = Particles(len(a))
    print(name)
    p.name = name
    p.type = "asteroid"
    p.host = "Sun"
    p.mass = mass
    from orbital_elements_to_cartesian import orbital_elements_to_pos_and_vel
    from orbital_elements_to_cartesian import orbital_period

    from amuse.ext.orbital_elements import new_binary_from_orbital_elements
    for i in range(len(p)):
        Ta = True_anomaly_from_mean_anomaly(numpy.rad2deg(Ma[i]), ecc[i])
        b = new_binary_from_orbital_elements(Mparent, p[i].mass, a[i], ecc[i], Ta, inc[i], LoAn[i], Aop[i], G=constants.G)
        p[i].position = b[1].position - b[0].position
        p[i].velocity = b[1].velocity - b[0].velocity
        
        rho = 2.0 | units.g/units.cm**3
        p.radius = (p.mass/rho)**(1./3.) 

    return p

def read_orbital_elements_from_MPCORB(filename="MPCORB.DAT", n=-1):
    #f = open(fdir+"MPCORB.DAT", "r")
    a = [] | units.AU
    ecc = []
    inc = [] 
    Ma = []
    Aop = []
    LoAn = []
    classify = []
    name = []
    Nobs = []
    U = []
    MPC_data = False
    for line in open(filename):
        if "00001" in line[0:6]:
            MPC_data=True
        if MPC_data and len(line.split())>10:
          No = line[117:122].strip()
          if len(No)>0 and \
                  ('E' not in line[105:107] or \
                   'D' not in line[105:107] or \
                   'D' not in line[105:107]):
              Nobs.append(int(No))
              Ma.append(float(line[26:35]))
              Aop.append(float(line[37:46]))
              LoAn.append(float(line[48:57]))
              inc.append(float(line[59:68]))
              ecc.append(float(line[70:79]))
              a.append(float(line[92:103]) | units.au)
              U.append(line[105:107])
              classify.append(line[161:165])
              name.append(line[166:194])
              if n>0 and len(a)>=n:
                  break

    return name, a, ecc, inc, Ma, Aop, LoAn

def read_orbital_elements_from_MinorPlanetCenter(filename="MPCORB.DAT", n=-1):
    
    if filename == "MPCORB.DAT":
        name, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_from_MPCORB(filename, n)
        return name, a, ecc, inc, Ma, Aop, LoAn
    
    #f = open(fdir+"MPCORB.DAT", "r")
    a = [] | units.AU
    ecc = []
    inc = [] 
    Ma = []
    Aop = []
    LoAn = []
    name = []
    for line in open(filename):
        Ma.append(float(line[26:35]))
        Aop.append(float(line[37:46]))
        LoAn.append(float(line[48:57]))
        inc.append(float(line[59:68]))
        ecc.append(float(line[70:79]))
        a.append(float(line[92:103]) | units.au)
        name.append(line[166:194])
        if n>0 and len(a)>=n:
            break
    return name, a, ecc, inc, Ma, Aop, LoAn

def add_asteroids_to_solar_system(solar_system, MPC_filename,
                             time_JD=2457099.5|units.day, n=-1):

    # download the asteroid ephemerids file
    if not os.path.isfile(MPC_filename):
        print("download ", MPC_filename)
        url = 'https://www.minorplanetcenter.net/iau/MPCORB/'+MPC_filename
        datafile = requests.get(url)
        open(MPC_filename, 'wb').write(datafile.content)

    sun = solar_system[solar_system.name=="Sun"][0]

    name, a, ecc, inc, Ma, Aop, LoAn = read_orbital_elements_from_MinorPlanetCenter(filename=MPC_filename,
                                                                      n=n)
    mass = 1000 | units.kg
    p = construct_particle_set_from_orbital_elements(name, mass, a, ecc, inc, Ma, Aop, LoAn, sun.mass)

    print(p)
    print("number of particles N=", len(p))
    solar_system.add_particles(p)
    solar_system.move_to_center()
  
    return solar_system

def new_lunar_system(Julian_date=-1|units.day):
    return new_lunar_system_in_time(Julian_date)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n", 
                      type="int", default = 10,
                      help="number of asteroids [%default]")
    result.add_option("--MPC_filename", 
                      dest="MPC_filename",default = "NEA.txt",
                      help="MPC data filename possibilities include: MPCORB.DAT, NEA.txt, PHA.txt, Distant.txt, Unusual.txt, [%default]")
    result.add_option("-f", 
                      dest="filename", default = "solar_system.amuse",
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = None,
                      help="output filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    solar_system = read_set_from_file(o.filename, 'hdf5', close_file=True)
    Julian_date = solar_system.Julian_date

    solar_system = add_asteroids_to_solar_system(solar_system,
                                                 o.MPC_filename,
                                                 Julian_date, o.n)
    print(solar_system)

    from amuse.plot import scatter
    from matplotlib import pyplot
    star = solar_system[solar_system.type=="star"]
    planet = solar_system[solar_system.type=="planet"]
    moon = solar_system[solar_system.type=="moon"]
    asteroid = solar_system[solar_system.type=="asteroid"]
    scatter(star.x, star.y, s=100, c='y')
    scatter(planet.x, planet.y, s=30, c='b')
    scatter(moon.x, moon.y, s=10, c='r')
    scatter(asteroid.x, asteroid.y, s=1, c='k')
    pyplot.show()
    print(planet)
    
    if o.outfile:
        write_set_to_file(solar_system,
                          o.outfile, 'amuse',
                          append_to_file=False,
                          overwite_file=True,
                          version="2.0")

    
