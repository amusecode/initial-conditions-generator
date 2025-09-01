#!/bin/sh
amuse_ic_new_plummer --nstars 16 --radius 0.1 -Q 0.5 -F plummer_N16R01pcQ05.amuse
amuse_ic_name_stars -f plummer_N16R01pcQ05.amuse --nstars 8 --name system -F plummer_N16R01pcQ05N8.amuse
amuse_ic_add_planets_oligarchic -f plummer_N16R01pcQ05N8.amuse --fplanets 1 --name system -F plummer_N16R01pcQ05N8Oligarch.amuse --rmin_disk 30 --rmax_disk 3000
amuse_ic_rotate -f plummer_N16R01pcQ05N8Oligarch.amuse -F plummer_N16R01pcQ05N8OligarchR.amuse --type star

