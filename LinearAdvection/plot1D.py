import numpy as np
import matplotlib.pyplot as plt
 
 
# Creating dataset
x = np.loadtxt("xt.txt")
Z = np.loadtxt("zt.txt")
print(Z.shape)

nx = x.size

Z = Z.reshape(-1,nx)

print(Z.shape)
 
# Creating figure
fig = plt.figure()
ax = plt.axes()
# ax.set_zlim([0,1])
# Creating plot
ax.plot(x, Z[2])

 
# show plot
plt.show()