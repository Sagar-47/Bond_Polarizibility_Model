import os
import numpy as np
import subprocess
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from numpy import linalg as LA
import math
#dataexp = np.genfromtxt('abinit_phondis.dat').T


a = [  0.498127022 ,  -0.867104052   , 0.000180692]
b = [        0.498127022   , 0.867104052   , 0.000180692]
c = [     -0.331751590  ,  0.000000000  ,  0.945842103 ]
	
def angle(a,b) :
	#print(LA.norm(a))
	#return np.dot(a,b)	
	return	np.arccos(np.dot(a,b)/(LA.norm(a)*LA.norm(b)))*180/math.pi


#Differentian of polarizibility in i,j direction for bond with respect to displacement of atom gamma in k-direction 
#Given polarizibility constants a_p,a_l,a_p_d,a_l_d for the bonds 
def diff_polar(bond,dis,a_p,a_l,a_p_d,a_l_d,i,j,k) : 
	if (i==j and i==k) :
		result = a_p_d*bond[k] + (a_l_d-a_p_d)*bond[i]*bond[j]*bond[k] + (a_l-a_p)*( (bond[j]*((-1/dis)*(1-(bond[i]*bond[k])))) + 			(bond[i]*((-1/dis)*(1-(bond[j]*bond[k])))) ) 
		#print(1)
	elif (i==j and i!=k) :
		result = a_p_d*bond[k] + (a_l_d-a_p_d)*bond[i]*bond[j]*bond[k] + (a_l-a_p)*( (bond[j]*(bond[i]*bond[k]/dis)) + 				(bond[i]*(bond[j]*bond[k]/dis)) ) 
		#print(2)
	elif (i!=j and i==k) :
		result = (a_l_d-a_p_d)*bond[i]*bond[j]*bond[k] + (a_l-a_p)*( (bond[j]*((-1/dis)*(1-(bond[i]*bond[k])))) + 				(bond[i]*(bond[j]*bond[k]/dis)) ) 
		#print(3)
	elif (i!=j and i!=k) :
		if(j==k) :
			result = (a_l_d-a_p_d)*bond[i]*bond[j]*bond[k] + (a_l-a_p)*( (bond[j]*(bond[i]*bond[k]/dis)) + 					(bond[i]*((-1/dis)*(1-(bond[j]*bond[k])))) ) 
			#print(4)
		else :
			result = (a_l_d-a_p_d)*bond[i]*bond[j]*bond[k] + (a_l-a_p)*( (bond[j]*(bond[i]*bond[k]/dis)) + 					(bond[i]*(bond[j]*bond[k]/dis)) )
			#print(5)
	return result

def dist(a,b): #To calculate distance 
    return np.linalg.norm(a-b)

def position(x): #To calculate position of atoms in cartesian co-ordinates given relative co-ordinates in terms of lattice vectors
    x_coord = x[0]*a[0] + x[1]*b[0] + x[2]*c[0]
    y_coord = x[0]*a[1] + x[1]*b[1] + x[2]*c[1]
    z_coord = x[0]*a[2] + x[1]*b[2] + x[2]*c[2]
    r = np.array([x_coord, y_coord, z_coord])
    return r
	
print(diff_polar(bond,dis,-2,-5,4,8,0,0,k))






