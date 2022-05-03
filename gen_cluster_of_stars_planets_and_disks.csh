python make_plummer_sphere.py --nstars 16 --radius 0.1 -Q 0.3 -F plummer_N16R01pcQ03.amuse 
python name_stars.py -f plummer_N16R01pcQ03.amuse --nstars 6 --name system -F plummer_N16R01pcQ03N.amuse
python add_planets_oligarch.py -f plummer_N16R01pcQ03N.amuse --fplanets 1 --name system -F plummer_N16R01pcQ03NOligarch.amuse --rmin_disk 3 --rmax_disk 300
python add_debris_disk.py -f plummer_N16R01pcQ03NOligarch.amuse --name system --rmin 3 --rmax 300 --ndisk 300 -F plummer_N16R01pcQ03NOligarch_A300.amuse
python rotate.py -f plummer_N16R01pcQ03NOligarch_A300.amuse -F plummer_N16R01pcQ03NOligarch_A1000R.amuse --type star

