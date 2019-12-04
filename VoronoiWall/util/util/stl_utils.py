import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot


def stl_plot(thething):
    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    axes.set_axis_off()

    triangles = mplot3d.art3d.Poly3DCollection(thething.vectors)
    triangles.set_alpha(0.5)
    triangles.set_edgecolor('k')
    axes.add_collection3d(triangles)

    # Auto scale to the mesh size
    scale = thething.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def polytope(vertices, faces):
    triangles = np.array([tri for face in faces for tri in triangulate(face)])
    poly = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
    for i, t in enumerate(triangles):
        for j, v in enumerate(t):
            poly.vectors[i][j] = vertices[v,:]
    return poly

def triangulate(f):
    for i in range(1, len(f)-1):
        yield (f[0], f[i], f[i + 1])
