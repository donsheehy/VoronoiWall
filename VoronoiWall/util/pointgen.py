import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import pickle

# Create helix:
def make_helix(n):
    theta_max = 8 * np.pi
    theta = np.linspace(0, theta_max, n)
    x, y, z = theta/(8 * np.pi) , 0.5 * np.sin(theta)+ 0.5, 0.5*np.cos(theta) + 0.5

    helix = np.stack((x, y, z))

    return helix.T

# Attach 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define no. data points and create helix:
n = 100
data = make_helix(n)

print(data.shape)

with open('helix.p', 'wb') as f:
    # data = numpy.array(list(zip(*data)))
    pickle.dump(data, f)

for i in range(n):
    # ax.scatter(data[0][i], data[1][i], data[2][i])
    ax.scatter(data[i][0], data[i][1], data[i][2])
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.show()
