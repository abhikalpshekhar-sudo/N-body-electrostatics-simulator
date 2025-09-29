import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Calculation of the electric field by a single charge
#
#	We know that E(r) \propto q * \vec{r}/r^3

def EFieldSingleCharge(q, r0, x, y):
	distance = np.hypot(x - r0[0], y - r0[1])
	return q * (x - r0[0]) / (distance**3), q * (y - r0[1]) / (distance**3)


nx = 128
ny = 128
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)

# And then we use meshgrid that creates an array of coordinates for X and Y

X, Y = np.meshgrid(x,y)

# Our charges array contains the charge carriers

charges = []

# We simply add all charges that we want ...

charges.append((1, (-0.5, 0)))
charges.append((-1, (0.5, 0)))

charges.append((1, (0.5, -0.5)))
charges.append((-1, (-0.5, -0.5)))

charges.append((1, (0.5, 0.5)))
charges.append((-1, (-0.5, 0.5)))


# Initialize our field components to be zero

Ex = np.zeros((ny, nx))
Ey = np.zeros((ny, nx))

# And iterate over charges. Since we apply superposition principle we can
# calculate the field created by each charge separatly and add up all fields in
# the field vectors

for charge in charges:
	ex, ey = EFieldSingleCharge(*charge, x=X, y=Y)
	Ex += ex
	Ey += ey

# Create a subplot

fig = plt.figure()
splot = fig.add_subplot(111)

# Color is determined by the magnitude of the field

color = np.log(np.hypot(Ex, Ey))

# Perform a plot of the vector arrows using streamplot

splot.streamplot(x,y,Ex, Ey, color=color, linewidth=0.5, cmap=plt.cm.inferno, density = 2, arrowstyle='->', arrowsize=1)

# Add circles for positive and negative charges

qColors = {
	True : '#FF0000',
	False : '#0000FF'
}
for q, pos in charges:
	splot.add_artist(Circle(pos, 0.05, color=qColors[q>0]))

# Set labels and areas

splot.set_xlabel('x')
splot.set_ylabel('y')
splot.set_xlim(-1,1)
splot.set_ylim(-1,1)
splot.set_aspect('equal')

# Done

plt.show()
