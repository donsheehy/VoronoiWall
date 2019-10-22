import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Define the 6 vertices of an octahedron.
vertices = np.array([\
    [0, 0, 0.7],
    [0, 1.1, 0],
    [1, 0, 0],
    [0, 0, -1.5],
    [0, -1.2, 0.5],
    [-1, 0, 0]
    ])
# Define the 8 triangular faces.
faces = np.array([\
    [0,1,2],
    [0,5,1],
    [1,5,3],
    [1,3,2],
    [2,3,4],
    [0,2,4],
    [0,4,5],
    [3,5,4],
    ])

# Create the mesh
oct = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j, v in enumerate(f):
        oct.vectors[i][j] = vertices[v,:]

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)
axes.set_axis_off()

triangles = mplot3d.art3d.Poly3DCollection(oct.vectors)
triangles.set_alpha(0.5)
triangles.set_edgecolor('k')
axes.add_collection3d(triangles)

# Auto scale to the mesh size
scale = oct.points.flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()
