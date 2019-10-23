grep -A3 "crystal axes: (cart. coord. in units of alat)" Mono_TiSe2.scf.out_reference | grep -v "crystal" | cut -b25-57 > cell_alat.dat

grep -A3 "reciprocal axes: (cart. coord. in units 2 pi/alat)" Mono_TiSe2.scf.out_reference | grep -v "reciprocal" | cut -b25-57 > cell_2pi_over_alat.dat
