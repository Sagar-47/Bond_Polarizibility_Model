import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import time


# reading frequencies and A_xx from text files of the form '1_CuGeS'
# this file was used to plot graphs using a particular component of A_matrix and list of frequencies 
def fileread(file):
	f = open(file, mode="r")
	freq = np.zeros(1)
	a_xx = np.zeros(1)
	line = f.readline()

	if (line.find("freq") != -1):
	    while (line.find("nd") == -1):
	        line = f.readline()
	        line = line[1:]
	        freq = np.append(freq, np.fromstring( line, sep = ' '))
	    while (len(line) != 0):
	        line = f.readline()
	        a_xx = np.append(a_xx, np.fromstring( line, sep = ' '))
	freq *= 10 ** 12        
	return np.append(freq, a_xx)

# call fileread function


def x_y(file, offset, a):
 abc = fileread(file)
 freq = abc[0:37]
 a_xx = abc[37:74]

# calculating the frequency factor corresponding to optical modes

 freq_fact = np.zeros(37)
 i = 4
 while(i<37):
  freq_fact[i] = boltz(freq[i])/freq[i] 
  i += 1

# calculating raman intensity corresponding to optical modes 

 raman_int = np.zeros(37)
 i = 4
 while(i<37):
  raman_int[i] =  10**15 *  a_xx[i]**2 * freq_fact[i]
  i += 1
 #print(raman_int)

# plotting the raman intensities 

 freq= freq /(29979245800)
 sigma = 1
 i = 4
 x = np.zeros((37, 100))
 y = np.zeros((37, 100))
 while(i<37):
  x[i] = np.linspace( freq[i] - 5*sigma, freq[i] + 5*sigma, 100 )
  y[i] = raman_int[i] * mlab.normpdf( x[i], freq[i], sigma) * math.sqrt(2*math.pi) 
  plt.plot( x[i], y[i]+offset, color=a)
  if(i==36):
   plt.plot( x[i], y[i]+offset, color=a, label=file)
  i += 1
 return 

# defining the bose factor (n_mu + 1)

def boltz(x):
    return ( 1 + (np.exp(1.61 * 10**(-13) * x) - 1)**(-1))

# plotting all files of format '1_CuGeS' together, comment out if plotting just one
f=0.0
f_max = 10.0
while(f<=f_max) : 
	file_name = str(int(f)+1) + "_CuGeS"
	j = 0.0# - f/(f_max+6)  
	k = 0.2 + f/(f_max+6)
	l = 0.2 + f/(f_max+6) 
	a=(j,k,l) 
	print(a)
	graph1 = x_y(file_name, f*20 ,  a)
	#time.sleep(3)
	#graph2 = x_y("2_CuGeS", 'g')
	#graph3 = x_y("3_CuGeS", 'r')
	#graph4 = x_y("4_CuGeS", 'c')
	#graph5 = x_y("5_CuGeS", 'm')
	#graph6 = x_y("6_CuGeS", 'y')
	f=f+1.0
plt.xlabel('Frequency [cm^(-1)]')
plt.ylabel('Raman intensity')
plt.title("Raman spectra XX")
#plt.legend()
plt.show()

