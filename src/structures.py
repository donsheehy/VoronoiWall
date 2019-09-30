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


class Node:
    def __init__(self):
        self.next = None
        self.prev = None

    def __repr__(self) -> str:
        return self.__str__()


class CircularList:
    def __init__(self):
        self.data = []
        self.size = 0

    def append(self, node):
        self.data.append(node)
        self.size += 1

    def __iter__(self):
        return iter(self.data)

    def __len__(self) -> int:
        return self.size

    def __str__(self) -> str:
        return str(self.data)

    def __repr__(self) -> str:
        return self.__str__()


class Vertex(Node):
    """
    A vertex belongs to two or more Edges. It is defined by a Cartesian coordinate.
    """
    def __init__(self, coordinate: ndarray) -> None:
        super().__init__()
        self.coordinate = coordinate
        self.edges = CircularList()

    def __eq__(self, other):
        return self.coordinate == other.coordinate

    def __hash__(self):
        return hash(str(self.coordinate))

    def __str__(self) -> str:
        return str(self.coordinate)


class Edge(Node):
    """
    An Edge belongs to one or more Regions. It has two vertices which bound the Edge.
    """
    def __init__(self, v1: Vertex, v2: Vertex) -> None:
        super().__init__()
        self.regions = CircularList()
        self.v1 = v1
        self.v2 = v2
        v1.edges.append(self)
        v2.edges.append(self)

    def __eq__(self, other):
        return np.array_equal(self.v1, other.v1) and np.array_equal(self.v2, other.v2)

    def __hash__(self):
        # hash the vertices as a tuple
        return hash((self.v1, self.v2))

    def __str__(self) -> str:
        return str(self.v1) + " -> " + str(self.v2)


class Region(Node):
    """
    A Region has a Point around which it exists and a list of Edges which bound the Region.
    """
    def __init__(self, point) -> None:
        super().__init__()
        self.point = point
        self.edges = CircularList()

    def __eq__(self, other):
        return self.point == other.point

    def __hash__(self):
        return hash(self.point)

    def __str__(self) -> str:
        return "Region about " + str(self.point) + " with edges: " + str(self.edges)


class Point(Node):
    def __init__(self, coordinate: ndarray) -> None:
        super().__init__()
        self.coordinate = coordinate
        self.region = Region(self)

    def __eq__(self, other):
        return self.coordinate == other.coordinate

    def __hash__(self):
        return hash(str(self.coordinate))

    def __str__(self) -> str:
        return "Point at " + str(self.coordinate)


def debugging(vor, points):
    print("Voronoi Vertices:")
    [print("{}".format(v), end=" ") for v in vor.vertices]
    print()
    print("Edges/Ridges:")
    for v in vor.ridge_vertices:
        if -1 not in v:
            print(" - {} -> {}".format(vor.vertices[v[0]], vor.vertices[v[1]]))
    print("Bounded Regions:")
    for i, r in enumerate(vor.regions):
        if -1 not in r:
            print("Region {}:".format(i))
            for idx in r:
                print(" - {}".format(vor.vertices[idx]))
    print(vor.point_region)
    plt.scatter(*zip(*points), c="blue")
    plt.scatter(*zip(*vor.vertices), c="red")
    plt.show()


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
