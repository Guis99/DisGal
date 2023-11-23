import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

 
 
# Creating dataset
x = np.loadtxt("xt.txt")
Z = np.loadtxt("zt.txt")
print(Z.shape)

nx = x.size

Z = Z.reshape(-1,nx)

print(Z.shape)
 
# Creating figure
fig, ax = plt.subplots()

# ax.set_zlim([0,1])
# Creating plot
line, = ax.plot(x, Z[0, :])
ax.axhline(1)
ax.set_ylim([-.5,2])

# Function to update the plot in each frame
def update(frame):
    line.set_ydata(Z[frame, :])
    return line,

# Set up the animation
num_frames = Z.shape[0]
stepsize = 1000/num_frames
ani = FuncAnimation(fig, update, frames=num_frames, interval=stepsize, blit=True)

# Show the animation
plt.show()