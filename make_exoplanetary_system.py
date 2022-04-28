import numpy
from amuse.lab import *
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
import xml.etree.ElementTree as ET, urllib.request, gzip, io

def read_planetary_system_from_database(name):
    url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
    oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.request.urlopen(url).read())))

    mass = [] | units.MSun
    sma = [] | units.au
    ecc = []
    inc = []
    for system in oec.findall(".//star"):
        if len(system.findall(".//planet"))>1:
            print(system.findtext("name"), len(system.findall(".//planet")))
            if name == system.findtext("name"):
                Mstar = float(system.findtext("./mass")) | units.MSun
                for planet in system.findall(".//planet"):
                    try:
                        mass.append(float(planet.findtext("./mass"))
                                    | units.MJupiter)
                    except:
                        mass.append(-1.0| units.MJupiter)
                    try:
                        sma.append(float(planet.findtext("./semimajoraxis"))
                                   | units.au)
                    except:
                        sma.append(-1| units.au)
                    try:
                        ecc.append(float(planet.findtext("./eccentricity")))
                    except:
                        ecc.append(-1)
                    try:
                        inc.append(float(planet.findtext("./inclination")))
                    except:
                        inc.append(-1)
                break
    print("m=", Mstar, mass, sma, ecc, inc)
    return Mstar, mass, sma, ecc, inc

def read_names(nmin):
    url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
    oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.request.urlopen(url).read())))

    names = []
    nplanets = []
    for system in oec.findall(".//star"):
        if len(system.findall(".//planet"))>nmin:
            print(system.findtext("name"), len(system.findall(".//planet")))
            names.append(system.findtext("name"))
            nplanets.append(len(system.findall(".//planet")))
    return names, nplanets

def minimal_positive_value(x, unit):
    xmin = max(x)
    print("x=", x)
    for xi in x:
        if xi>0| unit:
            xmin = min(xmin, xi)
    if xmin<0 | unit:
        xmin = 0.1 | unit
    return xmin

def maximal_positive_value(x, unit):
    xmax = max(x)
    if xmax<0 | unit:
        xmax = 3 | unit
    return xmax

def fill_in_missing_planets(mass, sma, ecc, inc):

    mmin = 0.1*min(mass)
    mmax = min(mmin, max(mass))
    amin = min(sma)
    amax = max(sma)
    emin = min(ecc)
    emax = max(ecc)
    imin = min(inc)
    imax = max(inc)
    a_extreme = 45 | units.au
    if amax<a_extreme:
        dm = (200|units.MEarth) - mass.sum()
        print("dm=", dm.in_(units.MEarth))
        while dm>0|units.MEarth:
            mass.append(new_salpeter_mass_distribution(1,
                                               mass_min=0.1*mmin,
                                               mass_max=mmax,
                                               alpha=-1.2))
            dm -= mass[-1]
            sma.append(10**numpy.random.uniform(
                numpy.log10(amax.value_in(units.au)),
                numpy.log10(1.05*a_extreme.value_in(units.au)),
                size=1) | units.au)
            print("a=", sma.in_(units.au), amin.in_(units.au), amax.in_(units.au))

            ecc.append(numpy.random.uniform(emin, emax, size=1)**2)
            inc.append(numpy.random.uniform(imin, imax, size=1))
            amax = max(sma)
            print("amax=", amax.value_in(units.au), dm.in_(units.MEarth))
            if amax>=a_extreme:
                break
            
    return mass, sma, ecc, inc

def make_consistent_planetary_system(mass, sma, ecc, inc):

    mmin = minimal_positive_value(mass, units.MEarth)
    mmax = maximal_positive_value(mass, units.MJupiter)
    amin = minimal_positive_value(sma, units.au)
    amax = maximal_positive_value(sma, units.au)
    emin = 0.0
    emax = max(ecc)
    if emax<0:
        emax = 0.1
    imin = min(inc)
    if imin>0:
        for i in range(len(inc)):
            inc[i] = inc[i] - imin
    else:
        imin = 0.0
    imax = max(inc)
    if imax<0:
        imax = 1.0 # deg.
    #imax = numpy.deg2rad(imax)
    #inc = numpy.deg2rad(inc)
    print("extremes:", mmin, mmax, amin, amax, emin, emax, imin, imax)

    for i in range(len(mass)):
        if mass[i]<=0|units.MSun:
            mass[i] = new_salpeter_mass_distribution(1,
                                               mass_min=mmin,
                                               mass_max=mmax,
                                               alpha=-1.2)
        if sma[i]<0|units.au:
            print("a=", amin.value_in(units.au), amax.value_in(units.au))
            sma[i] = 10**numpy.random.uniform(
                numpy.log10(amin.value_in(units.au)),
                numpy.log10(amax.value_in(units.au)), size=1) | units.au
        if ecc[i]<0:
            ecc[i] = numpy.random.uniform(emin, emax)**2
        if inc[i]<0:
            inc[i] = numpy.random.uniform(imin, imax)

    return mass, sma, ecc, inc

def add_planet(star, np_obs, mass, sma, ecc, inc): 

    print(mass, sma, ecc, inc) 
    #converter = nbody_system.nbody_to_si(star.mass, max(sma))
    planets = Particles()
    for i in range(len(mass)):
        ang_1, ang_2, ang_3  = numpy.random.uniform(0, 2*numpy.pi, 3)
        print(i, star[0].mass, mass[i], sma[i], ecc[i], inc[i])
        print(star[0].mass)
        print("<<=", type(mass[i].number))
        b = new_binary_from_orbital_elements(
            mass1=star[0].mass,
            mass2=mass[i],
            semimajor_axis=sma[i],
            eccentricity=ecc[i],
            true_anomaly=ang_1 | units.rad,
            inclination= inc[i] | units.rad,
            longitude_of_the_ascending_node=ang_2 | units.rad,
            argument_of_periapsis=ang_3 | units.rad,
            G=constants.G)
        b[1].position -= b[0].position
        b[1].velocity -= b[0].velocity
        b[1].mass = mass[i]
        b[1].semimajor_axis = sma[i]
        b[1].eccentricity = ecc[i]
        b[1].inclination = inc[i]
        b[1].id = i
        density =  3 | units.g/units.cm**3
        b[1].radius =  (b[1].mass/density)**(1./3.)
        if i<np_obs:
            b[1].name = "planet"
        else:
            b[1].name = "virtual_planet"
        planets.add_particle(b[1])
    planets.type = "planet"
    planets.position += star.position
    planets.velocity += star.velocity
    star.add_particles(planets)
    return star

def make_planetary_system(name, add_missing=False):
    Mstar, mass, sma, ecc, inc = read_planetary_system_from_database(name)
    print(Mstar)
    print(mass, sma, ecc, inc)
    star = Particles(1)
    star.mass = Mstar
    star.position = [0,0,0] | units.au
    star.velocity = [0,0,0] | units.kms
    star.name = name
    star.type = "star"
    star.radius = 1|units.RSun
    star.id = 0
    mass, sma, ecc, inc = make_consistent_planetary_system(mass, sma, ecc, inc)
    np_obs = len(mass)

    if add_missing:
        mass, sma, ecc, inc = fill_in_missing_planets(mass, sma, ecc, inc)
    
    star = add_planet(star, np_obs, mass, sma, ecc, inc)
    return star
   
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--name", 
                      dest="name", default = "",
                      help="system name [%default]")
    result.add_option("-F", "--filename", 
                      dest="filename",
                      default = "planetary_system.amuse",
                      help="output filename [%default]")
    result.add_option("--nmin", 
                      dest="nmin", type="int", default = "1",
                      help="minimum number of planets [%default]")
    result.add_option("--add_missing", 
                      dest="add_missing", default = False,
                      action = 'store_true',
                      help="add potentially missing planets [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    if len(o.name) == 0:
        names, np = read_names(o.nmin)
        for name, number in zip(names, np):
            print("star:", name, "has: ", number, "planets.")
    else:
        system = make_planetary_system(o.name, o.add_missing)
        print(system)
        write_set_to_file(system, o.filename, "amuse")
                      
