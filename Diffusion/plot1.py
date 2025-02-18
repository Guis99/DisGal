import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import json
import matplotlib.patches as patches

import subprocess
import sys

sys.path.insert(0,'../Utils')
import Utils
toRun = Utils.getExecutableName("diffDG")

# discretization parameters
deg = 5
div = 8
Lx = 1
Ly = 1
meshInfo = [str(deg), str(div), str(div), str(Lx), str(Ly)] # pack into list of strings

penalty = 50

force = "2*pi^2*sin(pi*x/1)*sin(pi*y/1)"
# force = "0"

numBoundaries = 4

btm = {"1":"0", "0":"1"} # if bc is not one, it has to be the other


ess = ["1","1","1","1"]
nat = [btm[bc] for bc in ess]

dirichletBC = ["sin(5*pi*x)", "-sin(pi*y)", "sin(pi*x)", "-sin(2*pi*y)"]
neumannBC = ["0", "0", "0", "0"]

dirichletBC_trimmed = [dirichletBC[i] for i in range(numBoundaries) if ess[i] == "1"]
neumannBC_trimmed = [neumannBC[i] for i in range(numBoundaries) if nat[i] == "1"]
subprocess.run([toRun, *meshInfo, str(penalty), force, str(numBoundaries), *ess, *nat, *dirichletBC_trimmed, *neumannBC_trimmed]) 
# subprocess.run([toRun, "3", str(div2), str(div2), "1", "1", force, "0", "0", "0", "0", "50"])  

def draw_cell_nr(cell, ax):
    # Adjusting for the center coordinates and level-dependent size
    for child in cell['children']:
        half_side_length = child['width'] / 2
        topLeftX = child['x'] - half_side_length
        topLeftY = child['y'] - half_side_length

        rect = patches.Rectangle((topLeftX, topLeftY), child['width'], child['width'], linewidth=.2, edgecolor='black', facecolor='none')
        # ax.text(child['x'], child['y'], str(child['CID']), ha='center', va='center', fontsize=8, color='black')

        ax.add_patch(rect)

# Step 1: Load your data
print("Loading results")
data = np.loadtxt('outputdG.txt', delimiter=',')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Step 2: Interpolate your data onto a regular grid
# Create grid coordinates
print("Interpolating data")
xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate z values on the grid
zi = griddata((x, y), z, (xi, yi), method='cubic')

print("Plotting data")
# Step 3: Plot the contour
fig,ax=plt.subplots()
plt.contourf(xi, yi, zi, levels=15, cmap=plt.cm.jet)
plt.colorbar()  # Show color scale
# plt.scatter(x, y, c=z, cmap=plt.cm.jet)  # Optionally, plot your original data points on top
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour plot')

print("Plotting quadtree")
with open("quadtreedG.json", 'r') as f:
    quadtree = json.load(f)

    draw_cell_nr(quadtree, ax)
    plt.axis('equal')  # Ensures the plot is square in shape


fig2 = plt.figure(figsize =(14, 9))
ax2 = plt.axes(projection ='3d')
 
# Creating plot
ax2.plot_surface(xi, yi, zi)
ax2.set_xlabel('x')


plt.show()


