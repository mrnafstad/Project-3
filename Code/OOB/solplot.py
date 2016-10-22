from matplotlib.pylab import *
from numpy import *

f = open("VerletTest.txt")

num =  len(f.readline().split()) # Numbers per line
planets = (num -1)/3       # Three coordinates per planet

read_lines = []

for line in f:
	read_values = [float(i) for i in line.split()]
	read_lines.append(read_values[1:])

Pos = array(read_lines)

plot(0, 0, "yo")
for i in range(planets):
	plot(Pos[:, 0 +3*i], Pos[:, 1 + 3*i], "--")
	hold('on')
xlabel("x[AU]")
ylabel("y[AU]")
#legend(["Sun", "Orbit of earth"])
show()
