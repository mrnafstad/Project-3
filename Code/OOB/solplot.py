from matplotlib.pylab import *
from numpy import *
'''
f = open("clusterVV_2_0.004.txt")

posx = []
posy = []

for line in f:
	pos = line.split()
	x = pos[1]
	y = pos[2]
	posx.append(float(x))
	posy.append(float(y))

posx = array(posx)
posy = array(posy)

plot(posx, posy)
show()
'''

f = open("EulerTest.txt")

posx = []
posy = []

for line in f:
	pos = line.split()
	x = pos[1]
	y = pos[2]
	posx.append(float(x))
	posy.append(float(y))

posx = array(posx)
posy = array(posy)

plot(posx, posy)
show()
