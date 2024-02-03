import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

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
plt.figure()
plt.contourf(xi, yi, zi, levels=15, cmap=plt.cm.jet)
plt.colorbar()  # Show color scale
plt.scatter(x, y, c=z, cmap=plt.cm.jet)  # Optionally, plot your original data points on top
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour plot of irregularly spaced data')
plt.show()
