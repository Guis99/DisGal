import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import subprocess

# discretization parameters
deg = 5
div = 10
Lx = 2

timeLength = 1.75
timeSteps = 500

meshInfo = [str(deg), str(div), str(Lx)] # pack into list of strings

initialCondition = "-((4 * x - 1) ^ 20 - 1)* .75"
initialCondition = "2*x - .5"
# initialCondition = "x"

toRun = "t1DAdv.exe" 

subprocess.run([toRun, *meshInfo, initialCondition, str(timeLength), str(timeSteps)])
 
# Creating dataset
print("Loading results")
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
ax.axhline(-1)
ax.set_ylim([-1.5,1.5])

# Function to update the plot in each frame
def update(frame):
    line.set_ydata(Z[frame, :])
    return line,

# Set up the animation
num_frames = Z.shape[0]
stepsize = 1000/num_frames
ani = FuncAnimation(fig, update, frames=num_frames, interval=stepsize/1, blit=True)

fig2, ax2 = plt.subplots()
ax2.axhline(1)
ax2.axhline(-1)
ax2.set_ylim([-1.5,1.5])
ax2.plot(x, Z[0,:])

# Show the animation
plt.show()