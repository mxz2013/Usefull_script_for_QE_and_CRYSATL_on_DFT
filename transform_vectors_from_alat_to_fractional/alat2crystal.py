import numpy as np
import csv
from numpy.linalg import inv
import sys

def nonblank_lines(f):
   for l in f:
      line = l.rstrip()
      if line:
         yield line
pop=1
nf=1
cell_alat = np.zeros((3,3))
for i in range(3):
    x=[]
    lines_1 = [line.rstrip('\n') for line in
               open('cell_alat.dat')]
    for plotPair in nonblank_lines(lines_1):
        if not plotPair.startswith("#"):
            xANDy = plotPair.split()
            x.append(float(xANDy[i].rstrip('\r'))) ##*alat) 
    cell_alat[i,:] = np.array(x)

cell_2pi_over_alat = np.zeros((3,3))

for i in range(3):
    x=[]
    lines_1 = [line.rstrip('\n') for line in
               open('cell_2pi_over_alat.dat')]
    for plotPair in nonblank_lines(lines_1):
        if not plotPair.startswith("#"):
            xANDy = plotPair.split()
            x.append(float(xANDy[i].rstrip('\r'))) #*2*np.pi/alat) 
    cell_2pi_over_alat[i,:] = np.array(x)

mat_real2reci = np.dot(cell_alat,inv(cell_2pi_over_alat))  ## rc = mat_trans cc

mat_reci2real = np.dot(cell_2pi_over_alat,inv(cell_alat))  ## cc = mat_trans rc

############################# check transform alat to 2p over alat #######################
#vector_alat = np.array([-2.000000 ,  3.464102 ,  0.000000 ])*alat
#vector_2p_over_alat = np.array([0.000000,  0.288675,  0.000000])*2*np.pi/alat
#print np.dot(mat_reci2real,vector_alat)
#print vector_2p_over_alat
#print cell_2pi_over_alat[:,1] #/(2*np.pi/alat)
############################# check transform #######################


##############################check transform 2p over alat to frac ###
#vector_2p_over_alat = np.array([0.0000000 , -0.1443376 ,  0.0000000])  #*2*np.pi/alat
#vector_frac =  np.array([0.0000000 , -0.5000000,   0.0000000])
#print np.dot(vector_2p_over_alat, inv(cell_2pi_over_alat))
#print np.dot( inv(cell_2pi_over_alat), vector_2p_over_alat)
#print vector_frac
##############################check transform 2p over alat to frac ###

#print np.dot(mat_reci2real, atp_alat)
#print atp_frac
#print np.dot(inv(cell_2pi_over_alat), atp_2pi_over_alat)

#sys.exit(1)
#####################################
for j in range(nf):  ## xrange(1000):
    w=[]
    x=[]
    y=[]
    z=[]
    lines_1 = [line.rstrip('\n') for line in
               open('atp-QE-pop'+str(pop)+'/scf_population'+str(pop)+'_'+str(j+1)+'.dat')]
    for plotPair in nonblank_lines(lines_1):
        if not plotPair.startswith("#"):
            xANDy =  plotPair.split()
            w.append(str(xANDy[0].rstrip('\r')))
            x.append(float(xANDy[1].rstrip('\r')))
            y.append(float(xANDy[2].rstrip('\r')))
            z.append(float(xANDy[3].rstrip('\r')))
    x2=[]
    y2=[]
    z2=[]
    w2=[]
    print len(x)
    for i in range(len(x)):
        vector_alat = np.zeros(3)
        vector_alat[0] = x[i]
        vector_alat[1] = y[i]
        vector_alat[2] = z[i]
        vector_2p_over_alat =  np.dot(mat_reci2real , vector_alat)

        vector_frac = np.dot(inv(cell_2pi_over_alat), vector_2p_over_alat)

        if w[i] == 'Ti':
            w2.append(22)
        if w[i] == 'Se':
            w2.append(252)
        x2.append(vector_frac[0]) ## because it is 2x2 unit cell
        y2.append(vector_frac[1])
        z2.append(vector_frac[2])

    with open ('atp-crystal-pop'+str(pop)+'/scf_frac_population'+str(pop)+'_'+str(j+1)+'.dat', 'w') as f:
            writer=csv.writer(f, delimiter='\t')
            writer.writerows( zip (w2, x2,y2,z2))

