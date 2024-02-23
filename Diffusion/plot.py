import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import json
import matplotlib.patches as patches

import subprocess

div = 4
force = "2*pi^2*sin(pi*x/1)*sin(pi*y/1)"
zero = "10"

exeSelect = 2
toRun = ".\main.exe" if exeSelect == 1 else ".\mainSplit.exe"

div2 = 4
subprocess.run([toRun, "2", str(div), str(div), "1", "1", "100*x*y", "sin(pi*x)", "-sin(3*pi*y)", "sin(pi*x)", "-sin(3*pi*y)", "90"]) 
# subprocess.run([toRun, "3", str(div2), str(div2), "1", "1", force, "0", "0", "0", "0", "50"])  

def draw_cell_nr(cell, ax):
    # Adjusting for the center coordinates and level-dependent size
    for child in cell['children']:
        half_side_length = child['width'] / 2
        topLeftX = child['x'] - half_side_length
        topLeftY = child['y'] - half_side_length

        rect = patches.Rectangle((topLeftX, topLeftY), child['width'], child['width'], linewidth=.2, edgecolor='black', facecolor='none')
        ax.text(child['x'], child['y'], str(child['CID']), ha='center', va='center', fontsize=8, color='black')

        ax.add_patch(rect)

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
# # Step 3: Plot the contour
fig,ax=plt.subplots()
plt.contourf(xi, yi, zi, levels=15, cmap=plt.cm.jet)
plt.colorbar()  # Show color scale
# plt.scatter(x, y, c=z, cmap=plt.cm.jet)  # Optionally, plot your original data points on top
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour plot')

print("Plotting quadtree")
with open("quadtree.json", 'r') as f:
    quadtree = json.load(f)

    draw_cell_nr(quadtree, ax)
    plt.axis('equal')  # Ensures the plot is square in shape


# fig2 = plt.figure(figsize =(14, 9))
# ax2 = plt.axes(projection ='3d')
 
# # Creating plot
# ax2.plot_surface(xi, yi, zi)
# ax2.set_xlabel('x')


plt.show()


