#!/bin/sh
amuse_ic_make_plummer_sphere --nstars 16 --radius 0.1 -Q 0.3 -F plummer_N16R01pcQ03.amuse 
amuse_ic_name_stars -f plummer_N16R01pcQ03.amuse --nstars 6 --name system -F plummer_N16R01pcQ03N.amuse
amuse_ic_add_planets_oligarchic -f plummer_N16R01pcQ03N.amuse --fplanets 1 --name system -F plummer_N16R01pcQ03NOligarch.amuse --rmin_disk 3 --rmax_disk 300
amuse_ic_add_debris_disk -f plummer_N16R01pcQ03NOligarch.amuse --name system --rmin 3 --rmax 300 --ndisk 300 -F plummer_N16R01pcQ03NOligarch_A300.amuse
amuse_ic_rotate -f plummer_N16R01pcQ03NOligarch_A300.amuse -F plummer_N16R01pcQ03NOligarch_A1000R.amuse --type star
