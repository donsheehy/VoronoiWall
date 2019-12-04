from stl_utils import polytope, stl_plot

import numpy as np
import scipy as sp

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
    n = 10
    # data = [[i/n, 0.01 * i**2 / (n ** 2), i**3/ (n**3)] for i in range(n)]
    data = [sp.rand(3) for i in range(n)]
    containerize(data, [[0, 1], [0, 1], [0, 1]])
    input = np.array(data)

    # Build the Voronoi Diagram.
    vor = Voronoi(input)

    # Generate the facets of some of the Voronoi regions.
    k = 5 # index of Voronoi cell to draw.
    facets = [vor.ridge_vertices[i] for i, pair in enumerate(vor.ridge_points) if k in pair]

    P = polytope(vor.vertices, facets)

    stl_plot(P)


if __name__ == '__main__':
    main()
