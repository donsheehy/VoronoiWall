import numpy
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)
axes.set_axis_off()

# Or creating a new mesh (make sure not to overwrite the `mesh` import by
# naming it `mesh`):
# VERTICE_COUNT = 100
# data = numpy.zeros(VERTICE_COUNT, dtype=mesh.Mesh.dtype)
# your_mesh = mesh.Mesh(data, remove_empty_areas=False)

# Load the STL files and add the vectors to the plot
your_mesh = mesh.Mesh.from_file('dons_sandbox/halfdonut.stl')
polygons = mplot3d.art3d.Poly3DCollection(your_mesh.vectors)
polygons.set_alpha(0.5)
polygons.set_edgecolor('k')

axes.add_collection3d(polygons)


# Auto scale to the mesh size
scale = your_mesh.points.flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)


# Show the plot to the screen
pyplot.show()
