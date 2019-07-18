import os
from matplotlib import pyplot as plt
import numpy as np
import subprocess
import math


#lattice constants in cartesian co-ordinates
a = [  3.2058296329 , -5.5526598046,   0.0000000000 ]
b = [  3.2058296329 ,  5.5526598046  , 0.0000000000  ]
c = [  -2.1372197553 ,  0.0000000000 ,  6.0449703273]


# value of mass replaced using run_raman.py
mass_ge = value_from_py_2 # mass of Ge
mass_cu = value_from_py_1 # mass of Cu

   

# Variable information
# r_pstns stores positions of atoms in fractional co-ordinates, see 'bonds.txt' to understand file structure
# r stores unit vectors of the bonds from an atom to its nearest neighbours, the vector is directed outwards from the atom
# dis stores the lengths of these bonds
# See block 1 and 2 for better understanding

r_pstns = np.zeros((1,3))
dis = np.zeros(1)
r= np.zeros((1,3))

# block 1
# Data extraction from 'bonds.txt' file, check the file to see structure (extraction method largely dependent on file structure)
# The positions of atoms are in fractional co-ordinates in terms of lattice parameters

f=open("bonds.txt",mode='r')
i=0

while True :
	line = f.readline()
	if (line.find('Bond') != -1):
		while True :		
			line = f.readline()
			line = line[16:]	
			if (len(line) == 0) :
				line = f.readline()
				break
			r_pstns=np.append(r_pstns,np.fromstring(line,dtype=float, sep=' ').reshape((1,3)),axis=0)
	if (len(line) == 0) :
		break

#print(r_pstns) 

def dist(a,b): #To calculate distance 
    return np.linalg.norm(a-b)

def position(x): #To calculate position of atoms in cartesian co-ordinates given relative co-ordinates in terms of lattice vectors
    x_coord = x[0]*a[0] + x[1]*b[0] + x[2]*c[0]
    y_coord = x[0]*a[1] + x[1]*b[1] + x[2]*c[1]
    z_coord = x[0]*a[2] + x[1]*b[2] + x[2]*c[2]
    r = np.array([x_coord, y_coord, z_coord])
    return r

# block 2
# To get R(unit vector) corresponding to each bond; i_max is 1 less than double the total number of bonds (95)
# i runs from 1 to 95 in steps of 2

i=1
j=1

while (i<97) :
	dis=np.append(dis,dist(position(r_pstns[i]),position(r_pstns[i+1])))
	r=np.append(r,((position(r_pstns[i+1])-position(r_pstns[i]))/dis[j]).reshape((1,3)),axis=0)
	j=j+1
	i=i+2
#print(r)
#print(dis)


def delta(a,b) : # defining the kronecker delta function

	if (a==b):
		return 1
	else :
		return 0

#print(delta(61,60))

# Differentiation of polarizibility in i,j direction for bond with respect to displacement of atom gamma in k-direction 
# Given polarizibility constants a_p,a_l,a_p_d,a_l_d for the bonds 
# defining the contribution of a particular bond to the B tensor
# Variable information: a_p is perpendicular polarizability, a_l is longitudinal polarizability, and '_d' signify their derivatives

def diff_polar(bond,dis,a_p,a_l,a_p_d,a_l_d,i,j,k) : 

	return delta(i,j)*a_p_d*bond[k] + (a_l_d-a_p_d)*bond[i]*bond[j]*bond[k] + (a_l-a_p)*((bond[j]*((-1/dis)*(delta(i,k)-(bond[i]*bond[k])))) + (bond[i]*((-1/dis)*(delta(j,k)-(bond[j]*bond[k])))) ) 


#Finding elements of B by summing over all bonds of a particular atom using diff_polar
def B_elements(start,end,a_p,a_l,a_p_d,a_l_d,i,j,k) :
	res = 0.0
	while(start<=end) :
		res = res + diff_polar(r[start,:],dis[start],a_p,a_l,a_p_d,a_l_d,i,j,k)
		start = start + 1
	return res
#print(B_elements(1,(alpha+4),20,2,20,2,1,0,2))


# Getting eigenmodes from 'cuges.modes', check the file to see structure (extraction method largely dependent on file structure)
f_1=open("value_from_py_3",mode='r')
i=0
k=0
em = np.zeros((1,6))
e_m=np.zeros((36,3))
omega = np.zeros(36)
while (k<36) :
	line = f_1.readline()
	if (line.find('freq') != -1):
		omega[k] = np.fromstring(line[21:36],dtype=float, sep=' ')
		k=k+1
		i=0
		while (i<12) :		
			line = f_1.readline()
			line = line[3:73]		
			#print(((np.fromstring(line,dtype=float, sep=' ').reshape((1,3))).shape))
			em=np.append(em,np.fromstring(line,dtype=float, sep=' ').reshape((1,6)),axis=0)
			i=i+1
e_m = em[:, [0, 2,4]]

#Getting eigenvectors in order Cu0, Cu1, Cu2, Cu3, Ge4, Ge5, S6, S7, S8, S9, S10, S11, this is the order in which information for bonds is stored
m = 0
while(m<36) :
	e_m[[(m*12)+2, (m*12)+3]] = e_m[[(m*12)+3,(m*12)+ 2]]
	e_m[[(m*12)+7, (m*12)+3]] = e_m[[(m*12)+3, (m*12)+7]]
	e_m[[(m*12)+11, (m*12)+4]] = e_m[[(m*12)+4, (m*12)+11]]
	e_m[[(m*12)+9, (m*12)+6]] = e_m[[(m*12)+6, (m*12)+9]]
	e_m[[(m*12)+11, (m*12)+8]] = e_m[[(m*12)+8, (m*12)+11]]
	e_m[[(m*12)+11, (m*12)+10]] = e_m[[(m*12)+10, (m*12)+11]]
	m = m + 1
#print(np.linalg.norm(e_m[36,:]))

#declaring all the masses of all the atoms in the same order as given above in eigenvectors block

mass = [mass_cu,mass_cu,mass_cu,mass_cu,mass_ge,mass_ge,32.06,32.06,32.06,32.06,32.06,32.06] 



#Finding Raman constant elements for eigenmode mu in direction i,j by summing in directions of vibration k and atoms gamma over product of B_elements with kth cartesian component of corresponding eigenvectors for a mode mu of atom gamma(gm_1) divided by square root of mass of atom gamma(gm_1)
#Start mu,i,j with zero
#Change gm upper limit to change number of bonds(atoms) to be considered
def A_elements(mu,i,j,a_p,a_l,a_p_d,a_l_d) :
	A_res = 0.0
	k=0	
	while(k<3):	
		gm = 1 
		gm_1 = 0
		while(gm<49):
			A_res = A_res + ((B_elements(gm,(gm+3),a_p,a_l,a_p_d,a_l_d,i,j,k))* (e_m[((12*mu)+1+gm_1)][k]) / math.sqrt(mass[gm_1]))
			gm = gm+4
			gm_1 = gm_1 +1
		k=k+1
	return A_res


#n_bonds is 1 more than the total no. of bonds in the system
def B_matrix(n_bonds,a_p,a_l,a_p_d,a_l_d,k):
	B_m = np.zeros((12,3,3))
	atm = 1
	atm_1 = 0
	while(atm<n_bonds):
		i=0
 		while(i<3):
			j=0
			while(j<3):
				nb = 3
				B_m[atm_1][i][j] = rescale(B_elements(atm,(atm+nb),a_p,a_l,a_p_d,a_l_d,i,j,k))
				j=j+1
			i=i+1	
		atm = atm+nb+1
		atm_1 = atm_1 +1
	return B_m
# the line below helps print the B_matrix in one cartesian direction for all the elements, the last argument chooses the direction with 0 for X, 1 for Y and 2 for Z

#print(B_matrix(49, 2, 5, -4, -8, 2))

def rescale(a):	# To remove numerical noise and scale small numbers to zero
	if(abs(a)<0.00001):
		return 0
	else :
		return a
# Store values for mode mu taking all possible combinations i,j where i=0, j=0 is x,x mode 
# p_i and p_j are used to print the right element from A_matrix for all modes, by setting p_i=2, p_j=2 for zz direction
# A_max is to find the maximum value in A_matrix from all the modes and all directions
def A_matrix(mu_max,a_p,a_l,a_p_d,a_l_d,A_max,p_i,p_j) :
	A_m = np.zeros((mu_max,3,3))
	mu=0
	while(mu<mu_max):
		i=0
 		while(i<3):
			j=0
			while(j<3):
				A_m[mu][i][j] = pow(A_elements(mu,i,j,a_p,a_l,a_p_d,a_l_d),2) #Squaring and storing the values 
				#A_m[mu][i][j] = A_elements(mu,i,j,a_p,a_l,a_p_d,a_l_d) 
				#if(i==p_i and j==p_j) :
				#	print(A_m[mu][i][j]) # To print the A_matrix elements [p_i, p_j] for all modes 
				#if(A_m[mu][i][j] >= A_max) :
				#	A_max = A_m[mu][i][j] # To find and store the maximum value in A_matrix from all the modes and all directions
				j=j+1
			i=i+1
		#print(np.sum(A_m[mu][0, : ]) + np.sum(A_m[mu][0, : ]) - A_m[mu][0, 0 ] )
		print(np.sum(A_m[mu]))
		mu=mu+1
	#print(A_max)
	return A_m
print("frequencies")
print(omega)
print("End") 

A_matrix_1 = A_matrix(36,2,5,-4,-8,0.113990493,1,1) # Evaluating A_matrix for given values of polarizibilities and their derivatives

#print(A_matrix(36,-2,-5,4,8,0.113990493))

