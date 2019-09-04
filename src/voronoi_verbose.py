"""
Proof of concept: plotting a 2D Voronoi diagram computed from input points.
This implementation is verbose to reveal the Voronoi data exposed by scipy in case it ends up being useful.

See: https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.Voronoi.html
See: https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi


def main():
    # this is the input array of points
    points = np.array([[-1, 0], [1, 0], [0, 1], [0, -1], [0, 0]])
    # use scipy to compute the Voronoi
    vor = Voronoi(points)

    # plot the input points
    plt.plot(points[:, 0], points[:, 1], 'o')

    '''
    vor.vertices is the array of coordinates of Voronoi diagram vertices.
    For this example:
        0: [-0.5  0.5]
        1: [0.5 0.5]
        2: [-0.5 -0.5]
        3: [ 0.5 -0.5]
    '''

    # plot the Voronoi vertices
    plt.plot(vor.vertices[:, 0], vor.vertices[:, 1], '*')

    # set the axis limits for the matplotlib plot
    plt.xlim(-5, 5)  # -5 < x < 5
    plt.ylim(-5, 5)  # -5 < y < 5

    '''
    vor.ridge_vertices is an array of index pairs.
    The pair identifies the two Voronoi vertex ends of a Voronoi edge.
    In other words, [0, 1] is an edge between vertex 0 and vertex 1: [-0.5, 0.5] <-> [0.5, 0.5]
    Note: a negative value for the starting vertex indicates it extends from infinity to the ending vertex
    For this example:
        0: [0, 1]
        1: [0, 2]
        2: [1, 3]
        3: [2, 3]
        4: [-1, 0]
        5: [-1, 1]
        6: [-1, 3]
        7: [-1, 2]
    '''

    # plot the lines that don't extend out to infinity
    for edge in vor.ridge_vertices:
        edge = np.asarray(edge)
        if np.all(edge >= 0):
            plt.plot(vor.vertices[edge, 0], vor.vertices[edge, 1], 'k-')

    '''
    vor.ridge_points is an array of index pairs. Each array element corresponds to a Voronoi edge/ridge.
    The pair identifies the two input points which the edge falls between.
    For this example: (8 total Voronoi edges)
        0: [4 2]
        1: [4 0]
        2: [4 1]
        3: [4 3]
        4: [2 0]
        5: [2 1]
        6: [1 3]
        7: [0 3]
    '''

    # plot the lines which do extend out to infinity
    center = points.mean(axis=0)
    for pointIdx, vertex in zip(vor.ridge_points, vor.ridge_vertices):
        vertex = np.asarray(vertex)
        if np.any(vertex < 0):
            i = vertex[vertex >= 0][0]  # finite end Voronoi vertex
            t = points[pointIdx[1]] - points[pointIdx[0]]  # tangent
            t = t / np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal
            midpoint = points[pointIdx].mean(axis=0)
            far_point = vor.vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * 100
            plt.plot([vor.vertices[i, 0], far_point[0]], [vor.vertices[i, 1], far_point[1]], 'k--')

    # view the finalized plot
    plt.show()


main()
