from typing import List, Dict, Set

import numpy as np
from numpy.core._multiarray_umath import ndarray
from scipy.spatial import Voronoi, ConvexHull
import matplotlib.pyplot as plt


def barycenter(pts: ndarray) -> List[float]:
    # split each of x, y, z, w, etc. into its own array
    splitAxes = np.squeeze(np.hsplit(pts, len(pts[0])))
    # take the average of x, y, etc. to find a midpoint
    return list(map(np.average, splitAxes))


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

    def __init__(self, coordinate: ndarray) -> None:
        self.coordinate = coordinate

    def __eq__(self, other):
        return np.array_equal(self.coordinate, other.coordinate)

    def __hash__(self):
        return hash(str(self.coordinate))

    def __str__(self) -> str:
        return str(self.coordinate)

    def __repr__(self):
        return self.__str__()


class Simplex:
    """
    Point > Region > Faces > [Simplices] > Vertices
    """

    def __init__(self, v1: Vertex, v2: Vertex, v3: Vertex, equation: ndarray) -> None:
        self.regions = CircularList()
        self.vertices = CircularList()
        for v in (v1, v2, v3):
            self.vertices.append(v)
        self.equation = equation
        self.neighbors = CircularList()

    def __eq__(self, other):
        if len(self.vertices) != len(other.vertices):
            return False
        for thisVertex, otherVertex in zip(self.vertices, other.vertices):
            if thisVertex not in other.vertices or otherVertex not in self.vertices:
                return False
        return True

    def __hash__(self):
        return hash(str(barycenter(np.array([self.vertices[0].coordinate, self.vertices[1].coordinate, self.vertices[2].coordinate]))))

    def __str__(self) -> str:
        return "Simplex composed of vertices: {}".format(str(self.vertices))

    def __repr__(self):
        return self.__str__()


class Region:
    """
    Point > [Region] > Faces > Simplices > Vertices
    """

    def __init__(self, point) -> None:
        self.point: Point = point
        self.simplices = CircularList()

    def __eq__(self, other):
        return self.point == other.point

    def __hash__(self):
        return hash(self.point)

    def __str__(self) -> str:
        return "Region about point {}.".format(str(self.point.coordinate))

    def __repr__(self):
        return self.__str__()


class Point:
    """
    [Point] > Region > Faces > Simplices > Vertices
    """

    def __init__(self, coordinate: ndarray) -> None:
        self.coordinate = coordinate
        self.region: Region = Region(self)

    def __eq__(self, other):
        return self.coordinate == other.coordinate

    def __hash__(self):
        return hash(str(self.coordinate))

    def __str__(self) -> str:
        return "Point at {}.".format(str(self.coordinate))

    def __repr__(self):
        return self.__str__()


class Diagram:
    def __init__(self):
        self.points = dict()
        self.regions = dict()
        self.faces = dict()
        self.simplices = dict()
        self.vertices = dict()


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
    vor = Voronoi(inputPoints)
    diagram = Diagram()
    for inputIdx, regionIdx in enumerate(vor.point_region):
        regionByIdx = vor.regions[regionIdx]
        pointObj = Point(vor.points[inputIdx])
        if -1 in regionByIdx or len(regionByIdx) < 4: # unbounded or 2D
            continue
        currentRegion = np.array(list(map(lambda idx: vor.vertices[idx], regionByIdx)))

        diagram.points[pointObj] = pointObj
        regionObj = pointObj.region
        diagram.regions[regionObj] = regionObj
        hull = ConvexHull(currentRegion)

        for i, simplexByIdx in enumerate(hull.simplices):
            p1, p2, p3 = list(map(lambda idx: hull.points[hull.vertices[idx]], simplexByIdx))
            vertices = []
            for coordinate in (p1, p2, p3):
                # make a new vertex object
                vertexObj = Vertex(coordinate)
                # if we've made this vertex before, reuse it
                if vertexObj not in diagram.vertices:
                    diagram.vertices[vertexObj] = vertexObj
                else:
                    vertexObj = diagram.vertices[vertexObj]
                vertices.append(vertexObj)
            v1, v2, v3 = vertices
            # make a new simplex object
            simplexObj = Simplex(v1, v2, v3, hull.equations[i])
            # if we've made this simplex before, reuse it
            if simplexObj not in diagram.simplices:
                diagram.simplices[simplexObj] = simplexObj
            else:
                simplexObj = diagram.simplices[simplexObj]
            # make the simplex, old or new, a part of this region
            regionObj.simplices.append(simplexObj)
            simplexObj.regions.append(regionObj)
            # discover neighboring simplices
            for neighborIdx in hull.neighbors[i]:
                neighborSimplexByIdx = hull.simplices[neighborIdx]
                p1, p2, p3 = list(map(lambda idx: hull.points[hull.vertices[idx]], neighborSimplexByIdx))
                vertices = []
                for coordinate in (p1, p2, p3):
                    # make a new vertex object
                    vertexObj = Vertex(coordinate)
                    # if we've made this vertex before, reuse it
                    if vertexObj not in diagram.vertices:
                        diagram.vertices[vertexObj] = vertexObj
                    else:
                        vertexObj = diagram.vertices[vertexObj]
                    vertices.append(vertexObj)
                v1, v2, v3 = vertices
                # make a new simplex object
                neighborSimplexObj = Simplex(v1, v2, v3, hull.equations[i])
                # if we've made this simplex before, reuse it
                if neighborSimplexObj not in diagram.simplices:
                    diagram.simplices[neighborSimplexObj] = neighborSimplexObj
                else:
                    neighborSimplexObj = diagram.simplices[neighborSimplexObj]
                simplexObj.neighbors.append(neighborSimplexObj)
    return diagram

points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
                [7, 7, 4],
                [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
                [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
                [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])
diagram = makeVoronoiDiagram(points)