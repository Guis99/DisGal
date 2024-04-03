import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import json
import matplotlib.patches as patches

import subprocess

# discretization parameters
deg = 8
div = 2
Lx = 2
Ly = 2
meshInfo = [str(deg), str(div), str(div), str(Lx), str(Ly)] # pack into list of strings

# SIPG penalty param
penalty = 90

# forcing terms in x and y directions
forces = ["0","0"]

# setting up boundary conditions
numBoundaries = 4
velDir = ["1","1","1","1"]
pressDir = ["0","0","0","0"]
nat = ["0","0","0","0"]

# dirichlet boundary conditions for x and y velocity and pressure
u = ["1","0","0","0"]
v = ["0","0","0","0"]
p = []

# natural (outlet and free) boundary conditions
natBC = []

exeSelect = 2
toRun = ".\mainStokes.exe"

subprocess.run([toRun, *meshInfo, str(penalty), *forces, str(numBoundaries), *velDir, *pressDir, *nat, *u, *v, *p, *natBC]) 

def draw_cell_nr(cell, ax):
    # Adjusting for the center coordinates and level-dependent size
    for child in cell['children']:
        half_side_length = child['width'] / 2
        topLeftX = child['x'] - half_side_length
        topLeftY = child['y'] - half_side_length

        rect = patches.Rectangle((topLeftX, topLeftY), child['width'], child['width'], linewidth=.2, edgecolor='black', facecolor='none')

        ax.add_patch(rect)

# Step 1: Load your data
print("Loading results...")
data = np.loadtxt('output.txt', delimiter=',')
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]
v = data[:, 3]
p = data[:, 4]

# Step 2: Interpolate your data onto a regular grid
# Create grid coordinates
print("Interpolating data...")
xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate z values on the grid
ui = griddata((x, y), u, (xi, yi), method='cubic')
vi = griddata((x, y), v, (xi, yi), method='cubic')
pi = griddata((x, y), p, (xi, yi), method='cubic')

xq = np.linspace(min(x), max(x), 20)
yq = np.linspace(min(y), max(y), 20)
xq, yq = np.meshgrid(xq, yq)

uq = griddata((x, y), u, (xq, yq), method='cubic')
vq = griddata((x, y), v, (xq, yq), method='cubic')

print("Plotting data...")
# # Step 3: Plot the contour
fig,ax=plt.subplots()
c1=plt.contourf(xi, yi, pi, levels=15, cmap=plt.cm.jet)
cb1 = fig.colorbar(c1)
cb1.set_label('Pressure')
c2=plt.quiver(xq, yq, uq, vq, angles='xy', scale_units='xy', scale=5)
cb2 = fig.colorbar(c2)  # Show color scale
cb2.set_label('Velocity magnitude')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour plot')

print("Plotting quadtree...")
with open("quadtree.json", 'r') as f:
    quadtree = json.load(f)

    draw_cell_nr(quadtree, ax)
    plt.axis('equal')  # Ensures the plot is square in shape


# fig2 = plt.figure(figsize =(14, 9))
# ax2 = plt.axes(projection ='3d')
 
# # Creating plot
# ax2.plot_surface(xi, yi, pi)
# ax2.set_xlabel('x')


plt.show()


