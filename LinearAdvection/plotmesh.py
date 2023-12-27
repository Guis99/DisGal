import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

 
 
# Creating dataset
x = np.loadtxt("xt.txt")

nx = x.size

 
# Creating figure
fig, ax = plt.subplots()

# ax.set_zlim([0,1])
# Creating plot
line, = ax.plot(x, np.zeros(x.size),"X")

# Show the animation
plt.show()