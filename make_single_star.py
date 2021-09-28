from amuse.lab import *
import sys
import numpy
import csv

def read_Gaia_database(filename, index=0, age=0|units.Myr):
    star = Particle()
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            if i==index:
                star.name = row[1]
                star.type = "star"
                star.mass = float(row[2]) | units.MSun
                star.position = (float(row[4]),
                                 float(row[5]),
                                 float(row[6])) | units.parsec
                star.velocity = (float(row[7]),
                                 float(row[8]),
                                 float(row[9])) | units.kms
                break
            i+=1
    s = SeBa()
    s.particles.add_particle(star)
    s.evolve_model(age)
    star.mass = s.particles[0].mass
    star.radius = s.particles[0].radius
    s.stop()
    print("star=", star)
    return star

def read_particle(filename, index=0):
    stars = read_set_from_file(filename, "hdf5")
    star = Particles(0)
    if index<=len(stars):
        star.add_particle(stars[index])
    return star

def Hill_radius(Mstar, a, Mplanet):
    return a * (Mplanet/(3.0*Mstar))**(1./3.)

def determine_orbital_elements(solar_system):
        from amuse.ext.orbital_elements import new_binary_from_orbital_elements


def make_single_star(m_star, radius, age, name):
    star = Particle()
    star.ZAMS_mass = m_star
    star.mass = m_star
    star.name = name
    star.type = "star"
    if radius<0|units.RSun:
        s = SeBa()
        s.particles.add_particle(star)
        s.evolve_model(age)
        star.mass = s.particles[0].mass
        star.radius = s.particles[0].radius
        s.stop()
    else:
        star.radius = radius
    star.position = (0,0,0) | units.AU
    star.velocity = (0,0,0) | units.kms
    return star

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile", default = "star.amuse",
                      help="output filename [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="m_star", type="float",
                      default = 1|units.MSun,
                      help="stellar mass [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="r_star", type="float",
                      default = -1.0|units.RSun,
                      help="stellar radius [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="age", type="float",
                      default = 0.0|units.Myr,
                      help="stellar age [%default]")
    result.add_option("--gaia_database_file", 
                      dest="gaia_database_file", 
                      #default = "/home/spz/Lib/Catalog/Gaia_NearbyStars/pos_vel_stars_10pc.csv",
                      default = "/home/spz/Lib/Catalog/Gaia_NearbyStars/gaia_nearby_stars.amuse",
                      help="Gaia database dir and filename [%default]")
    result.add_option("--gaia_database_index", type="int",
                      dest="gaia_database_index", 
                      default = -1,
                      help="Gaia database index [%default]")
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

    bodies = Particles(0)
    if o.gaia_database_index>0:
        if ".amuse" in o.gaia_database_file:
            new_bodies = read_particle(o.gaia_database_file, o.gaia_database_index)
            print(new_bodies)
        elif ".csv" in o.gaia_database_file:
            new_bodies = read_Gaia_database(o.gaia_database_file, o.gaia_database_index)
        else:
            print("Unrecognized input filename: ", o.gaia_database_file)

    else:
        new_bodies = make_single_star(o.m_star, o.r_star, o.age, o.name)
    bodies.add_particle(new_bodies)
                   
    index = 0
    time = 0|units.Myr
    if o.outfile==None:
        filename = "star.amuse".format(index)
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      overwrite_file=True,
                      append_to_file=False,
                      version="2.0")
