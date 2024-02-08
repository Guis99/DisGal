import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import subprocess

subprocess.run([".\main.exe", "1", "150", "150", "5", "5", "0", "sin(4*pi*x/5)", "0", "0", "0", "50"]) 

def draw_cell(cell, ax):
    # Adjusting for the center coordinates and level-dependent size
    half_side_length = cell['width'] / 2
    topLeftX = cell['x'] - half_side_length
    topLeftY = cell['y'] - half_side_length

    rect = patches.Rectangle((topLeftX, topLeftY), cell['width'], cell['width'], linewidth=.25, edgecolor='black', facecolor='none')

    ax.add_patch(rect)

    if not cell['isLeaf']:
        for child in cell['children']:
            draw_cell(child, ax)

# Step 1: Load your data
# Replace 'your_data_file.txt' with the path to your actual data file
print("Loading results")
data = np.loadtxt('output.txt', delimiter=',')
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
plt.title('Contour plot of irregularly spaced data')

print("Plotting quadtree")
with open("quadtree.json", 'r') as f:
    quadtree = json.load(f)

    draw_cell(quadtree, ax)
    plt.axis('equal')  # Ensures the plot is square in shape

fig2 = plt.figure(figsize =(14, 9))
ax2 = plt.axes(projection ='3d')
 
# Creating plot
ax2.plot_surface(xi, yi, zi)
ax2.set_xlabel('x')


plt.show()


