Here I write a script to find the atomic displacement by diagonalizing the Hessian matrix 
printed either from Quantum espresso "ph.x" or CRYSTAL code. The example I show here is using
the phonon calcualtion from CRYSTAL code where "HESSFREQ" is printed. 


The reason why I need this script is as follows:
1. I want to know the charge density wave distorsion of TiTe2. This is the aim of this small project.

2. I know in advance the high temperature phase of TiTe2 (space group 164, 3 atoms in unit cell), 
and I konw in low temperature phase it is a 2x2 distorsion, meaning that the stable strucutre is a 2x2 cell. But I do not know to which space group this low-T phase belongs.

3. I calculate the phonon dispersion by making a 2x2 cell in its high temperature phase. I should obtain imaginary (negative) phonon frequencies because the 2x2 cell is not stable in space group 164.

4. Each negative phonon frequency corresponds to one distorsion partten, i.e., displace of atoms according to the vector calculated by diagonalizing Hessian matrix.


5. The script "find_distortion.py" will do this work! It diagonalizes the correst Hessian reading from HESSFREQ.dat file. Then it finds the vectors corresponding to the unstable phonon, then transform these vectors from Cartesian coordinates to fractional coordinates.

6. "grep-hessian.sh" prepare inputs Hessian blocks.


With these vectors, we need to calcualte the total energies by displacing the atoms along 
each distorsion step by step. After that we will have a plot "energy_gain VS displacement" 
where we can identify approximately how big distorsion will make the 2x2 cell stable. 
Once we have the atomic position that has the largest energy gain, we can use "https://stokes.byu.edu/iso/findsym.php" to find the space group of the low T phase (i.e., space group 150). Finally we can relax using space group 150, and the atomic positions giving the largest energy gain,
 to relax the structure. 

I will add the process for finding the low T structure later.




