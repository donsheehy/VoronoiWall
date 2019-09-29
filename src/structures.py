from typing import List

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
        self.tail = None
        self.size = 0

    def append(self, node):
        if self.size == 0:
            self.tail = node
            self.tail.next = self.tail
        else:
            node.next = self.tail.next
            self.tail.next = node
        self.tail = self.tail.next
        self.size += 1

    def __len__(self) -> int:
        return self.size

    def __str__(self) -> str:
        if self.size == 0:
            return "[]"
        s = []
        current = self.tail
        while True:
            current = current.next
            s.append(current)
            if current is self.tail:
                break
        return str(s)

    def __repr__(self) -> str:
        return self.__str__()


class Vertex(Node):
    def __init__(self, coordinate: List[float]) -> None:
        super().__init__()
        self.coordinate = coordinate
        self.edges = CircularList()

    def __str__(self) -> str:
        return str(self.coordinate)


class Edge(Node):
    def __init__(self, v1: Vertex, v2: Vertex, region):
        super().__init__()
        v1.edges.append(self)
        v2.edges.append(self)
        self.v1 = v1
        self.v2 = v2
        self.region = region

    def __str__(self) -> str:
        return str(self.v1) + " -> " + str(self.v2)


class Region(Node):
    def __init__(self, point) -> None:
        super().__init__()
        self.point = point
        self.edges = CircularList()

    def __str__(self) -> str:
        return "Region about " + str(self.point) + ":\n" + str(self.edges)


class Point(Node):
    def __init__(self, coordinate):
        super().__init__()
        self.coordinate = coordinate
        self.region = Region(self)

    def __str__(self) -> str:
        return "Point at " + str(self.coordinate) + "."


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


def parse(vor: Voronoi, pts: ndarray) -> List[Region]:
    regions: CircularList = CircularList()
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
                # make a new Edge and append it to this Region's edge list
                print("v1: {}".format(vertices[vorRegion[i - 1]]))
                print("v2: {}".format(vertices[vorRegion[i]]))
                region.edges.append(Edge(vertices[vorRegion[i - 1]], vertices[vorRegion[i]], region))
                print(region.edges)
            regions.append(region)
            print()
    return regions


points: ndarray = np.array(
    [[0, 0], [6, 0], [0, 6], [6, 6], [-6, 0], [-6, 6], [12, 0], [12, 6], [0, -6], [6, -6], [0, 12], [6, 12]])
vor = Voronoi(points)

#debugging(vor, points)
#
regions: CircularList = parse(vor, points)
#
# print(regions)
# print(vor.regions)
e = CircularList()
e.append(Point(1))
e.append(Point(2))
e.append(Point(3))
e.append(Point(4))
e.append(Point(5))
print(regions)
print(e)