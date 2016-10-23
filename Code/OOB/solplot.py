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
color=iter(cm.rainbow(np.linspace(0.2, 1, planets)))
for i in range(planets):
	c = next(color)
	plot(Pos[:, 0 +3*i], Pos[:, 1 + 3*i], c = c)
	hold('on')
xlabel("x[AU]")
ylabel("y[AU]")
legend(["Sun", "Earth"])#, "Jupiter"])
show()
