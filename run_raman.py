import os
from matplotlib import pyplot as plt
import numpy as np
import subprocess
import math

i=0
i_max = 10
while (i<=i_max) :

	# To change the mass of atoms involved
	value_from_py1 = 65.2227 + (((63.7323 - 65.2227)/i_max)*i) # Changing mass of Cu
	value_from_py2 = 66.1311 + (((71.9079 - 66.1311)/i_max)*i) # Changing mass of Ge
	os.system('sed -e s/value_from_py_1/%f/g raman_1.py > raman_2.py' % value_from_py1) # Writing values in raman_1.py
	os.system('sed -i.bak -e s/value_from_py_2/%f/g raman_2.py' % value_from_py2)

	# To change the name of modes file from which eigenvectors have to be extracted
	value_from_py3 = "mix.modes_0." + str(i)  # Check file names of .modes files and edit accordingly 
	os.system('sed -i s/value_from_py_3/%s/g raman_2.py' % value_from_py3) 

	# Creating files with Frequencies and Raman tensor values for each mode
	file_name = str(i+1) + "_CuGeS"		
	#print(file_name)
	os.system('python raman_2.py > %s' % file_name)
	i=i+1

