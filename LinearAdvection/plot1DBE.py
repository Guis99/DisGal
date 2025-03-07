import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import subprocess
import platform

# discretization parameters
deg = 4
div = 50
Lx = 2

cfl = 1 / ((deg + 1)**2) / 4
# cfl = 1 / (2 * deg + 1) / 4

print(cfl)

timeLength = 4.
timeStepSize = cfl * Lx / div
timeSteps = timeLength / timeStepSize

print(timeLength, timeSteps)

integrators = {"Forward Euler": 0, "Crank-Nicholson": 1, "RK4": 2, "GL1": 3, "GL2": 4}
integratorType = "RK4"

systemForm = 0 # 0 = fully non-linear state vector, 1 = pseudo-linear matrix form

meshInfo = [str(deg), str(div), str(Lx)] # pack into list of strings

initialCondition = "-((4 * x - 1) ^ 20 - 1)* .75"
initialCondition = "2*x - .5"
initialCondition = ".5*sin(3*2*x*pi)"

baseline = 1.
cutoff = 10.
# initialCondition = "1"
# initialCondition = "x"

osName = platform.system()
if osName == "darwin" or osName.startswith("linux"):
    toRun = "./build/Adv1DBE"
elif osName == "windows":
    toRun = "./build/Adv1DBE.exe"


subprocess.run([toRun, *meshInfo, initialCondition, str(timeLength), 
                str(timeSteps), str(integrators[integratorType]), 
                str(systemForm), str(baseline), str(cutoff)])
 
# Creating dataset
print("Loading results")
x = np.loadtxt("x_be.txt")
Z = np.loadtxt("z_be.txt")

nx = x.size

Z = Z.reshape(-1,nx)
 
# Creating figure
fig, ax = plt.subplots()

# ax.set_zlim([0,1])
# Creating plot
line, = ax.plot(x, Z[0, :])
ax.axhline(1+baseline)
ax.axhline(-1+baseline)
ax.set_ylim([-1.5+baseline,1.5+baseline])

# Function to update the plot in each frame
def update(frame):
    line.set_ydata(Z[frame, :])
    return line,

# Set up the animation
num_frames = Z.shape[0]
stepsize = 1000/num_frames
ani = FuncAnimation(fig, update, frames=num_frames, interval=stepsize/1, blit=True)

# fig2, ax2 = plt.subplots()
# ax2.axhline(1)
# ax2.axhline(-1)
# ax2.set_ylim([-1.5,1.5])
# ax2.plot(x, Z[0,:])

# Show the animation
plt.show()