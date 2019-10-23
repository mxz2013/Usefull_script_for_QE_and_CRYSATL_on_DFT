The scripts here are usefull when you want to transform atomic positions in alat units
to fractional units. I write this script because I have thousands of atomic configurations 
in alat units which can be a good input for Quantum Espresso, however, they are not suitable
to CRYSTAL code that only reads atomic positions in fractional units.

This script can be generalized to any basis transformation from one basis set to another.
"Mono_TiSe2.scf.out_reference" is one output from quantum espresso code where atomic positions
in both units are printed so that it can be our reference for testing our scipt.
In addition, it produces necessary inputs for our script, i.e., the two basis sets.

"grep-cell.sh" extracts the two basis sets from the output of quantum espresso code.


"alat2crystal.py" is a python script that takes a 3D vector (e.g., [x,y,z]) and transform
it from one basis to another.






