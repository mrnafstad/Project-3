from matplotlib.pylab import *
from numpy import *

f = open("clusterVV_2_0.000.txt")

posx = []
posy = []

for line in f:
	pos = line.split()
	x = pos[3]
	y = pos[4]
	posx.append(float(x))
	posy.append(float(y))

posx = array(posx)
posy = array(posy)

plot(posx, posy)
show()
