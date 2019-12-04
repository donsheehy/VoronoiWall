import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import scipy as sp
import pickle

from scipy.spatial import Voronoi

def reflect(points, x = 1, y = 1 , z = 1, dx = 0, dy = 0, dz = 0):
    newpoints = [[x * px + dx, y * py + dy, z * pz + dz] for (px, py, pz) in points]
    return newpoints

def containerize(points, box):
    left = reflect(points, x = -1, dx = box[0][0])
    right = reflect(points, x = -1, dx = 2*box[0][1])
    bottom = reflect(points, y = -1, dy = box[1][0])
    top = reflect(points, y = -1, dy = 2*box[1][1])
    front = reflect(points, z = -1, dz = box[2][0])
    back = reflect(points, z = -1, dz = 2 * box[2][1])
    points.extend(left + right + top + bottom+ front + back)


def main():
    # Prepare the points.
    n = 1200
    # data = [[i/n, 0.01 * i**2 / (n ** 2), i**3/ (n**3)] for i in range(n)]
    # data = [sp.rand(3) for i in range(n)]

    data = list(pickle.load(open("./helix.p", "rb")))
    containerize(data, [[0, 1], [0, 1], [0, 1]])
    input = np.array(data)

    # Build the Voronoi Diagram.
    vor = Voronoi(input)

    # Prepare to draw.
    ax = a3.Axes3D(plt.figure())
    ax.set_axis_off()
    facets = []

    # Generate the facets of some of the Voronoi regions.
    k = set(range(0,100,2)) # indices of Voronoi regions to plot.
    for i, pair in enumerate(vor.ridge_points):
        if k & set(pair):
            facets.append([list(vor.vertices[i]) for i in vor.ridge_vertices[i]])

    # Plot the points.
    # for x,y,z in data[:n]:
        # ax.scatter(x,y,z, color='k')

    # Plot the facets.
    polygons = a3.art3d.Poly3DCollection(facets)
    polygons.set_alpha(0.7)
    polygons.set_edgecolor('k')
    ax.add_collection3d(polygons)

    plt.show()


if __name__ == '__main__':
    main()
