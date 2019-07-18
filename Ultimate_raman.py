import os
from matplotlib import pyplot as plt
import numpy as np
import subprocess
import math

i=0
i_max = 5
while (i<=i_max) :
	value_from_py1 = 57918.6358105374 + (((59616.6564336298 - 57918.6358105374)/i_max)*i)
	value_from_py2 = 66198.1953060669 - (((66198.1953060669 - 59616.6564336298)/i_max)*i)
	
	os.system('sed -e s/value_from_py_1/%f/g raman_1.py > raman_2.py' % value_from_py1)
	os.system('sed -i.bak -e s/value_from_py_2/%f/g raman_2.py' % value_from_py2)
	if(i==0): 
		value_from_py3 = "mix.modes_0.1"
		os.system('sed -i s/value_from_py_3/%s/g raman_2.py' % value_from_py3)
		os.system('python raman_2.py > 1_CuGeS')
	elif(i==1):
		value_from_py3 = "mix.modes_0.2"
		os.system('sed -i s/value_from_py_3/"mix.modes_0.2"/g raman_2.py')
		os.system('python raman_2.py > 2_CuGeS')
	elif(i==2):
		value_from_py3 = "mix.modes_0.3"
		os.system('sed -i s/value_from_py_3/"mix.modes_0.3"/g raman_2.py')
		os.system('python raman_2.py > 3_CuGeS')

	elif(i==3):
		value_from_py3 = "mix.modes_0.4"
		os.system('sed -i s/value_from_py_3/"mix.modes_0.4"/g raman_2.py')
		os.system('python raman_2.py > 4_CuGeS')

	elif(i==4):
		value_from_py3 = "mix.modes_0.5"
		os.system('sed -i s/value_from_py_3/"mix.modes_0.5"/g raman_2.py')
		os.system('python raman_2.py > 5_CuGeS')
	else : 
		value_from_py3 = "mix.modes_0.6"
		os.system('sed -i s/value_from_py_3/"mix.modes_0.6"/g raman_2.py')
		os.system('python raman_2.py > 6_CuGeS')
	i=i+1

