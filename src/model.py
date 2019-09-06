import numpy as np
import matplotlib.pyplot as plt

"""
Generates an STL file for 3D printing a Voronoi diagram. (eventually)
"""


def prepare(voronoi, precision):
    """
    Prepares a 3D printable model of the given Voronoi diagram.

    :param voronoi: computed Voronoi attributes
    :param precision: some minimum precision value to produce an accurate 3D model
    :return:
    """
    # TODO: the rest


def __subdivideFace(points, maxLength):
    """
    Given the vertices of a 2D shape, subdivides the inner area into triangular shapes (necessary for 3D printing) using
    the Delaunay triangulation recursively until no triangle has an edge length larger than the specified maximum.

    Adapted from: https://stackoverflow.com/a/24952758

    :param points: a Python array of input points; this array is modified in-place
    :param maxLength: keep shrinking regions until all edges are shorter than this
    :return: unused
    """

    from scipy.spatial import Delaunay
    from math import sqrt, ceil

    print("Delaunay triangulation with " + str(len(points)) + " points.")

    triangulation = Delaunay(points)

    maxLengthSquared = maxLength ** 2
    # get set of edges from the simplices
    edges = set()
    for simplex in triangulation.simplices:
        # print(simplex)
        # simplex is one triangle: [4 5 17]
        edges.add((simplex[0], simplex[1]))
        edges.add((simplex[1], simplex[2]))
        edges.add((simplex[0], simplex[2]))
    # check if all edges are small enough, subdividing recursively if not
    isFinished = True
    for edge in edges:
        p1, p2 = edge
        [x1, y1] = points[p1]
        [x2, y2] = points[p2]
        lengthSquared = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)  # delay the costly sqrt until necessary
        if lengthSquared > maxLengthSquared:
            isFinished = False
            # split into how many pieces?
            nPieces = ceil(sqrt(lengthSquared / maxLengthSquared))
            for piece in range(1, int(nPieces)):
                points.append([x1 + piece / float(nPieces) * (x2 - x1), y1 + piece / float(nPieces) * (y2 - y1)])
    if not isFinished:
        __subdivideFace(points, maxLength)


def __plotDelaunayTriangles(points):
    """
    Display the subdivided face in 2D with matplotlib.

    Adapted from: https://stackoverflow.com/a/24952758
    :param points: points to plot, connecting them by simultaneously visualizing the Delaunary triangulation
    :return: unused
    """
    from scipy.spatial import Delaunay

    npPoints = np.array(points)
    triangulation = Delaunay(npPoints)
    plt.triplot(npPoints[:, 0], npPoints[:, 1], triangulation.simplices)
    plt.plot(npPoints[:, 0], npPoints[:, 1], 'o')
    plt.show()


def main():
    # here's an example of what's going on.
    points = [[0, 0], [10, 3], [9.5, 4]]
    __subdivideFace(points, 0.5)
    __plotDelaunayTriangles(points)
    # apparently that's how 3D models have to look.
    # the granularity could be tweaked to match a real world constraint (0.1mm faces, e.g.)


main()