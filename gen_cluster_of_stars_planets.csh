python make_plummer_sphere.py --nstars 16 --radius 0.1 -Q 0.5 -F plummer_N16R01pcQ05.amuse 
python name_stars.py -f plummer_N16R01pcQ05.amuse --nstars 8 --name system -F plummer_N16R01pcQ05N8.amuse
python add_planets_oligarch.py -f plummer_N16R01pcQ05N8.amuse --fplanets 1 --name system -F plummer_N16R01pcQ05N8Oligarch.amuse --rmin_disk 30 --rmax_disk 3000
python rotate.py -f plummer_N16R01pcQ05N8Oligarch.amuse -F plummer_N16R01pcQ05N8OligarchR.amuse --type star

