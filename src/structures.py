from typing import List, Dict, Set

import numpy as np
from numpy.core._multiarray_umath import ndarray
from scipy.spatial import Voronoi
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
            pos = (pos + 1) % len(self.items)

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
        return self.coordinate == other.coordinate

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

    def __init__(self, face: Face, v1: Vertex, v2: Vertex, v3: Vertex) -> None:
        self.face = face
        self.vertices = CircularList()
        for v in (v1, v2, v3):
            self.vertices.append(v)
        self.neighbors = CircularList()

    def __eq__(self, other):
        if len(self.vertices) != len(other.vertices):
            return False
        for thisVertex, otherVertex in zip(self.vertices, other.vertices):
            if thisVertex not in other.vertices or otherVertex not in self.vertices:
                return False
        return True

    def __hash__(self):
        # bad hash function, but inherently error detecting
        return hash(self.face)

    def __str__(self) -> str:
        return "Simplex composed of vertices: {}".format(str(self.vertices))

    def __repr__(self):
        return self.__str__()


class Face:
    """
    Point > Region > [Faces] > Simplices > Vertices
    """

    def __init__(self, region: Region):
        self.region = region
        self.simplices = CircularList()

    def __hash__(self):
        # bad hash function, but inherently error detecting
        return hash(self.region)

    def __eq__(self, other):
        if len(self.simplices) != len(other.simplices):
            return False
        for thisSimplex, otherSimplex in zip(self.simplices, other.simplices):
            if thisSimplex not in other.simplices or otherSimplex not in self.simplices:
                return False
        return True

    def __str__(self):
        return "Face composed of simplices: {}".format(str(self.simplices))

    def __repr__(self):
        return self.__str__()


class Region:
    """
    Point > [Region] > Faces > Simplices > Vertices
    """

    def __init__(self, point: Point) -> Region:
        self.point: Point = point
        self.faces = CircularList()

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

    def __init__(self, coordinate: ndarray) -> Point:
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
    class AutoDict(dict):
        def __getitem__(self, key):
            if key not in self:
                return self.setdefault(key, key)
            else:
                return super().__getitem__(key)

    def __init__(self):
        self.points = set()
        self.regions = set()
        self.faces = set()
        self.simplices = set()
        self.vertices = set()

"""
def parse(vor: Voronoi, pts: ndarray) -> CircularList:
    foundRegions: CircularList = CircularList()
    # edges might repeat; keep track of them
    foundEdges: Dict = dict()
    # make all of the Vertex objects, maintaining their original ordering
    vertices: List[Vertex] = list(map(Vertex, vor.vertices))
    for regionIdx, inputPoint in zip(vor.point_region, pts):
        # parse the region vor.regions[regionIdx]
        # ignoring unbounded regions and <2D regions
        vorRegion = vor.regions[regionIdx]
        if -1 not in vorRegion and len(vorRegion) > 2:
            # a Region has a Point and a Point has a Region
            point: Point = Point(inputPoint)
            region: Region = point.region
            # form edges between pairs of points, including last -> first
            for i in range(0, len(vorRegion)):
                v1: Vertex = vertices[vorRegion[i - 1]]
                v2: Vertex = vertices[vorRegion[i]]
                # already found this edge? use it. otherwise, make one and "find" it
                edge: Edge = foundEdges.setdefault(Edge(v1, v2), Edge(v1, v2))
                # tell the edge it's getting a new region
                edge.regions.append(region)
                # tell the region it's getting a new edge
                region.edges.append(edge)
            foundRegions.append(region)
    return foundRegions

# just a square, but reflected in every cardinal direction about its own edges
points: ndarray = np.array(
    [[0, 0], [6, 0], [0, 6], [6, 6], [-6, 0], [-6, 6], [12, 0], [12, 6], [0, -6], [6, -6], [0, 12], [6, 12]])
vor = Voronoi(points)

# debugging(vor, points)

regions: CircularList = parse(vor, points)
list(map(print, regions))

"""