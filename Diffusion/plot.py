import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_cell(cell, ax):
    # Adjusting for the center coordinates and level-dependent size
    half_side_length = cell['width'] / 2
    topLeftX = cell['x'] - half_side_length
    topLeftY = cell['y'] - half_side_length

    rect = patches.Rectangle((topLeftX, topLeftY), cell['width'], cell['width'], linewidth=1, edgecolor='black', facecolor='none')
    # if cell['isLeaf']:
        # ax.text(cell['x'], cell['y'], str(cell['CID']), ha='center', va='center', fontsize=8, color='blue')

    ax.add_patch(rect)

    if not cell['isLeaf']:
        for child in cell['children']:
            draw_cell(child, ax)

# Step 1: Load your data
# Replace 'your_data_file.txt' with the path to your actual data file
data = np.loadtxt('output.txt', delimiter=',')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Step 2: Interpolate your data onto a regular grid
# Create grid coordinates
xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate z values on the grid
zi = griddata((x, y), z, (xi, yi), method='cubic')

# Step 3: Plot the contour
fig,ax=plt.subplots()
plt.contourf(xi, yi, zi, levels=15, cmap=plt.cm.jet)
plt.colorbar()  # Show color scale
# plt.scatter(x, y, c=z, cmap=plt.cm.jet)  # Optionally, plot your original data points on top
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour plot of irregularly spaced data')

with open("quadtree.json", 'r') as f:
    quadtree = json.load(f)

    draw_cell(quadtree, ax)
    plt.axis('equal')  # Ensures the plot is square in shape


plt.show()


