from typing import List, Dict, Set

import numpy as np
from numpy.core._multiarray_umath import ndarray
from scipy.spatial import Voronoi, ConvexHull
import matplotlib.pyplot as plt


def barycenter(pts):
    # split each of x, y, z, w, etc. into its own array
    splitAxes = np.squeeze(np.hsplit(pts, len(pts[0])))
    # take the average of x, y, etc. to find a midpoint
    return list(map(np.average, splitAxes))  # type: ndarray


class CircularList(list):
    """
    Extends the built-in list class.
    """

    def cycle(self):
        """
        Generator that loops from the end of the list back to the beginning.
        Use next() to increment the generator.
        :return: essentially an infinite iterator
        """
        pos = 0
        while True:
            yield self[pos]
            pos = (pos + 1) % len(self)

    def __hash__(self):
        tot = 0
        for item in self:
            tot += hash(item)
        return tot


class Vertex:
    """
    Point > Region > Faces > Simplices > [Vertices]
    """

    def __init__(self, coordinate):
        self.coordinate = coordinate  # type: ndarray

    def __eq__(self, other):
        return np.array_equal(self.coordinate, other.coordinate)

    def __hash__(self):
        return hash(str(self.coordinate))

    def __str__(self):
        return str(self.coordinate)

    def __repr__(self):
        return self.__str__()


class Simplex:
    """
    Point > Region > Faces > [Simplices] > Vertices
    """

    def __init__(self, v1, v2, v3, equation):
        self.regions = CircularList()  # type: CircularList
        self.vertices = CircularList()  # type: CircularList
        for vertex in (v1, v2, v3):
            self.vertices.append(vertex)
        self.equation = equation  # type: ndarray
        self.neighbors = CircularList()  # type: CircularList

    def __eq__(self, other):
        """
        A Simplex is considered equivalent to another Simplex if they are composed of the same three vertices.
        :param other: another Simplex
        :return: True if equivalent
        """
        if len(self.vertices) != len(other.vertices):
            return False
        for thisVertex, otherVertex in zip(self.vertices, other.vertices):
            if thisVertex not in other.vertices or otherVertex not in self.vertices:
                return False
        return True

    def __hash__(self):
        """
        Computes a hash related to the average of the three vertices.
        This way the ordering of the vertices doesn't matter.
        :return: a unique but repeatable hash
        """
        midpoint = barycenter(
            np.array([self.vertices[0].coordinate, self.vertices[1].coordinate, self.vertices[2].coordinate]))
        return hash(str(midpoint))

    def __str__(self):
        return "Simplex composed of vertices: {}".format(str(self.vertices))

    def __repr__(self):
        return self.__str__()


class Region:
    """
    Point > [Region] > Faces > Simplices > Vertices
    """

    def __init__(self, point):
        self.point = point  # type: Point
        self.simplices = CircularList()  # type:CircularList

    def __eq__(self, other):
        """
        A Region is considered equivalent to another Region if the defining input point of each is the same.
        :param other: another Region
        :return: True if equivalent
        """
        return self.point == other.point

    def __hash__(self):
        return hash(self.point)

    def __str__(self):
        return "Region about point {}.".format(str(self.point.coordinate))

    def __repr__(self):
        return self.__str__()


class Point:
    """
    [Point] > Region > Faces > Simplices > Vertices
    """

    def __init__(self, coordinate):
        self.coordinate = coordinate  # type: ndarray
        self.region = Region(self)  # type: Region

    def __eq__(self, other):
        """
        A Point is considered equivalent to another Point if they have the same cartesian coordinate.
        :param other: another Point
        :return: True if equivalent
        """
        return self.coordinate == other.coordinate

    def __hash__(self):
        return hash(str(self.coordinate))

    def __str__(self):
        return "Point at {}.".format(str(self.coordinate))

    def __repr__(self):
        return self.__str__()


class Diagram:
    def __init__(self):
        class AutoDict(dict):
            """
            Attempting to get an item from this dictionary will either get the item as normal, or,
            in the case the the item did not exist, add the item to the dictionary and return it.

            Be aware that overridden __eq__ methods for these types means that equivalency requirements
            may be less strict than expected. For example, two simplices with antiparallel normal vectors but
            the same set of 3 vetices are considered equivalent.
            """

            def __getitem__(self, item):
                return super().setdefault(item, item)

        self.points = AutoDict()  # type: AutoDict
        self.regions = AutoDict()  # type: AutoDict
        self.faces = AutoDict()  # type: AutoDict
        self.simplices = AutoDict()  # type: AutoDict
        self.vertices = AutoDict()  # type: AutoDict


# class Face:
#     """
#     Point > Region > [Faces] > Simplices > Vertices
#     """

#     def __init__(self, region):
#         self.region = region
#         self.simplices = CircularList()

#     def __hash__(self):
#         # bad hash function, but inherently error detecting
#         return hash(self.region)

#     def __eq__(self, other):
#         if len(self.simplices) != len(other.simplices):
#             return False
#         for thisSimplex, otherSimplex in zip(self.simplices, other.simplices):
#             if thisSimplex not in other.simplices or otherSimplex not in self.simplices:
#                 return False
#         return True

#     def __str__(self):
#         return "Face composed of simplices: {}".format(str(self.simplices))

#     def __repr__(self):
#         return self.__str__()

def makeVoronoiDiagram(inputPoints: ndarray):
    # compute Voronoi using SciPy
    vor = Voronoi(inputPoints)
    # start tracking individual elements of the diagram
    diagram = Diagram()
    # process each region
    for inputIdx, regionIdx in enumerate(vor.point_region):
        # find the region that corresponds to the current input point
        regionByIdx = vor.regions[regionIdx]

        # ignore an unbounded or 2D region
        if -1 in regionByIdx or len(regionByIdx) < 4:
            continue

        # make a new point object from the input point for this region
        pointObject = diagram.points[Point(vor.points[inputIdx])]

        # convert from index-based vertices to an array of vertices defining this region
        currentRegion = np.array(list(map(lambda idx: vor.vertices[idx], regionByIdx)))

        # obtain and register the region created by constructing the point object
        regionObject = diagram.regions[pointObject.region]

        # since we've just got the set of vertices which bound the region, the convex hull should be the region itself
        hull = ConvexHull(currentRegion)

        for simplexIdx, simplexByIdx in enumerate(hull.simplices):
            # convert from index-based vertices to 3 coordinates defining this simplex
            p1, p2, p3 = list(map(lambda idx: hull.points[hull.vertices[idx]], simplexByIdx))

            vertices = []
            for coordinate in (p1, p2, p3):
                # if we've made this vertex before, the same object will be reused
                vertexObject = diagram.vertices[Vertex(coordinate)]
                vertices.append(vertexObject)

            v1, v2, v3 = vertices

            # if we've made this simplex before, the same object will be reused
            simplexObject = diagram.simplices[Simplex(v1, v2, v3, hull.equations[simplexIdx])]

            # make the simplex, old or new, a part of this region
            regionObject.simplices.append(simplexObject)
            simplexObject.regions.append(regionObject)

            # discover neighboring simplices
            for neighborIdx in hull.neighbors[simplexIdx]:
                neighborSimplexByIdx = hull.simplices[neighborIdx]
                p1, p2, p3 = list(map(lambda idx: hull.points[hull.vertices[idx]], neighborSimplexByIdx))

                vertices = []
                for coordinate in (p1, p2, p3):
                    # if we've made this vertex before, the same object will be reused
                    vertexObject = diagram.vertices[Vertex(coordinate)]
                    vertices.append(vertexObject)

                v1, v2, v3 = vertices

                # if we've made this simplex before, the same object will be reused
                neighborSimplexObject = diagram.simplices[Simplex(v1, v2, v3, hull.equations[simplexIdx])]

                # it's a neighboring simplex, so add it to our neighbors list
                simplexObject.neighbors.append(neighborSimplexObject)
                '''
                it could also be implied that we are a neighbor of that simplex, but the list
                of neighbor relationships will include that anyway, and I don't see a realistic way
                to skip that particular iteration when we get there, so I'm just letting it happen naturally
                by not adding simplexObj to neighborSimplexObj's neighbor list at this time
                '''
    return diagram


if __name__ == '__main__':
    import pickle
    data = None
    with open('helix.p', 'rb') as f:
        data = pickle.load(f)

    points = np.array(data)
    # points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
    #                    [7, 7, 4],
    #                    [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
    #                    [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
    #                    [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])

    diagram = makeVoronoiDiagram(points)
