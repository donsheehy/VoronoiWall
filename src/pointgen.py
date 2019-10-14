import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

# Create helix:
def make_helix(n):
    theta_max = 8 * np.pi
    theta = np.linspace(0, theta_max, n)
    x, y, z = theta, np.sin(theta), np.cos(theta)
    helix = np.vstack((x, y, z))

    return helix

# Attach 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define no. data points and create helix:
n = 100
data = [make_helix(n)]

for i in range(n):
    ax.scatter(data[0][0][i], data[0][1][i], data[0][2][i])

ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.show()