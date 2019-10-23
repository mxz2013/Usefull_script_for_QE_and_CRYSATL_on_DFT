import numpy as np
import sys
import csv
from numpy import linalg as LA
from numpy.linalg import inv

#######
Ry2cmm1 = 219474.63/2. # Units conversion 
Hartree2cmm1 = 1./0.4556335E-05   # Units conversion
nat = 12                          # number of atoms
dim = nat*3                       # dimension of hession matrix 
mass = np.zeros(nat) 
amu2Ry = 43628.1015381455/47.867    # Units conversion
m1_cry = 47.90                    # Mass of Ti in amu
m2_cry = 129.9067                # Mass of Te in amu
m1 = m1_cry*amu2Ry  #(43628.1015381455/47.867)
m2 = m2_cry*amu2Ry   #(71967.6373587642/78.96)

atom_order = [1,2,2,1,1,2,2,2,1,2,2,2]

for i in range(nat):         
   if atom_order[i]==1: 
      mass[i] = m1
   else:
      mass[i] = m2

def find_eigen_hessian(N, M):
## diagonalization of 1/sqrt(mA*mB) * \frac{\partial^2 E}{\partial u_1u_2}
## if mA=mB = 1, the eigenvectors are the displacement of atoms
## is mA and mB are correst masses, the eigenvalues are \omega_q^2
## where \omega_q is the phonon frequency
    ReMat =   np.zeros((N,N))  # in principle at Gamma all elemets are real,
    ImMat =   np.zeros((N,N))  #however, we still keep for general situations.
    for na in range(nat):
        for nb in range(nat):
            for jc in range(3): ## x, y, and z coordinates.
                jmode = 3*nb+jc
                m = 2*jc  ## real part
                n = 2*jc+1  ## imaginary part
                x = []
                y = []
                lines_1 = [line.rstrip('\n') for line in
                           open('hessian_'+str(na+1)+'-'+str(nb+1)+'.dat')]
                for plotPair in lines_1:
                    if not plotPair.startswith("#"):
                        xANDy = plotPair.split()
                        x.append(float(xANDy[m].rstrip('\r'))/np.sqrt(M[na]*M[nb])) ##*alat)
                        y.append(float(xANDy[n].rstrip('\r'))/np.sqrt(M[na]*M[nb])) ##*alat)
                for ic in range(3):
                    imode = 3*na+ic
                    ReMat[imode, jmode] = x[ic]*2.0 # Hartree/Bohr**2 to Ry/Bohr**2
                    ImMat[imode, jmode] = y[ic]*2.0 #


    Mat = ReMat + 1.j*ImMat 
    e_vals, e_vecs = LA.eig(Mat)
    re_evals = e_vals.real # sorting it according to the real part
    idx = re_evals.argsort()
    e_vals = e_vals[idx]
    e_vecs = e_vecs[:,idx]


    return e_vals, e_vecs
E_with_mass, V_with_mass = find_eigen_hessian(dim, mass)
E_no_mass, V_no_mass = find_eigen_hessian(dim, M=np.ones(nat))
####check if diagonalization is correct or not
####by comparing the eigenvalues with CRYSTAL output
print "the eigenvalue with mass (references are -115.25   -115.25   -115.25)"
print (1.j*(np.sqrt(E_with_mass[0:3])*Ry2cmm1)).real
#print "the vector for 1",
#print V_no_mass[:,0].real
#print "the vector for 2",
#print V_no_mass[:,1].real
#print "the vector for 3",
#print V_no_mass[:,2].real
####check if diagonalization is correct or not
####by comparing the eigenvalues with CRYSTAL output

##transform the vector in cartesian to fraction
b_cart = np.identity(3)
#       DtIRECT LATTICE VECTORS COMPONENTS (ANGSTROM)
#B1     6.542    -3.777     0.000
#B2     0.000     7.554     0.000
#B3     0.000     0.000    12.990
b_dirt = ([[ 6.542,    -3.777,     0.000],
           [ 0.000,     7.554,     0.000],
           [ 0.000,     0.000,    12.990]])
V_real = V_no_mass.real
#sys.exit()
# if b_dirt = P*b_cart then P = b_dirt*inv(b_cart)
# we say that the P transforms  
mat_trans_c2d = np.dot(b_dirt, inv(b_cart)) 
print mat_trans_c2d
for i in range(len(E_with_mass)): 
    if E_with_mass[i].real < 0: # nonstable phonon mode
        V_cart = np.zeros((nat, 3))
        V_frac = np.zeros((nat, 3))
        d = []
        x = []
        y = []
        z = []
        for j in range(dim):
            if (j % 3 == 0):
                mod_xy = np.sqrt((V_no_mass[j,i].real)**2+ (V_no_mass[j+1,i].real)**2)
                d.append(mod_xy)
                x.append(V_no_mass[j,i].real)
                y.append(V_no_mass[j+1,i].real)
                z.append(V_no_mass[j+1,i].real)
        V_cart [:, 0] = x
        V_cart [:, 1] = y
        V_cart [:, 2] = z
        atom_name = []
        for l in range(nat):
            #print "invmat", inv(mat_trans_c2d)
            #print "V", V_cart[l,:]
            V_frac[l,:] = np.dot(V_cart[l,:], inv(mat_trans_c2d))  # in
            if atom_order[l] == 1:
                atom_name.append(22)
            else:
                atom_name.append(252)
            #fractional units 
#        with open('disp_eigenvalue_'+str(i+1)+'.dat', 'w') as f:      
#           writer = csv.writer(f, delimiter = '\t')
#           writer.writerows (zip(atom_order, d, x ,y, np.array(x)/x[0], np.array(y)/x[0]))
        with open('eigenvector_in_fractional_units_for_eigenvalue-'+str(i+1)+'.dat', 'w') as f:      
           writer = csv.writer(f, delimiter = '\t')
           writer.writerows (zip(atom_name, V_frac[:,0] ,V_frac[:,1] ,V_frac[:,2] ))  









