import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from math import sqrt, ceil

"""
Generates an STL file for 3D printing a Voronoi diagram. (eventually)
"""

class DelaunayTris:
    def __init__(self, points=[]):
        self.points = points
        self.numPointsOriginal = len(points)

        self.target_point = -1



        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.cid_press = self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.triangulate()

    def triangulate(self):
        self.__subdivideFace(self.points, 0.5)
        self.__plotDelaunayTriangles(self.points, self.numPointsOriginal)

    def prepare(self, voronoi, precision):
        """
        Prepares a 3D printable model of the given Voronoi diagram.

        :param voronoi: computed Voronoi attributes
        :param precision: some minimum precision value to produce an accurate 3D model
        :return:
        """
        # TODO: the rest


    def __subdivideFace(self, points, maxLength, depth=0):
        """
        Given the vertices of a 2D shape, subdivides the inner area into triangular shapes (necessary for 3D printing) using
        the Delaunay triangulation recursively until no triangle has an edge length larger than the specified maximum.

        Adapted from: https://stackoverflow.com/a/24952758

        :param points: a Python array of input points; this array is modified in-place
        :param maxLength: keep shrinking regions until all edges are shorter than this
        :return: unused
        """

        #Randomly run into error where subdivision enters an infinite loop
        #have a depth counter to just exit if subdivision goes past the arbitrary level
        if depth > 10:
            return


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
            self.__subdivideFace(points, maxLength, depth=depth+1)


    def __plotDelaunayTriangles(self, points, numOriginalPoints):
        """
        Display the subdivided face in 2D with matplotlib.

        Adapted from: https://stackoverflow.com/a/24952758
        :param points: points to plot, connecting them by simultaneously visualizing the Delaunary triangulation
        :return: unused
        """




        self.ax.clear()
        npPoints = np.array(points)
        triangulation = Delaunay(npPoints)
        self.ax.triplot(npPoints[:, 0], npPoints[:, 1], triangulation.simplices)
        self.vertices = self.ax.plot(npPoints[0:numOriginalPoints, 0], npPoints[0:numOriginalPoints, 1], 'or', markersize=6)
        self.ax.plot(npPoints[numOriginalPoints:-1, 0], npPoints[numOriginalPoints:-1, 1], 'o', markersize=4)

        if self.target_point >= 0 and not self.is_inner(self.target_point):
            self.ax.plot(npPoints[self.target_point, 0], npPoints[self.target_point, 1], 'ok', markersize=6)

        #Not sure why this fixes a stack overflow error but it does
        for vert in self.vertices:
            vert.figure.canvas.draw()

    def on_press(self, event):
        print("press: " + str(event.xdata) + " " + str(event.ydata))
        print(self.target_point)
        min_dist_squared = 1e10
        close_point = -1;
        for i in range(self.numPointsOriginal):
            dist=(event.xdata - self.points[i][0])*(event.xdata - self.points[i][0]) + (event.ydata - self.points[i][1])*(event.ydata - self.points[i][1])
            if dist < min_dist_squared:
                min_dist_squared = dist
                close_point = i
        if sqrt(min_dist_squared)<0.1:
            self.target_point = close_point
            print(self.target_point)
        else:
            self.target_point = -1
            print(self.target_point)


    def on_release(self, event):
        if self.target_point == -1:
            self.add_point(event)
        else:
            self.points = self.points[0:self.numPointsOriginal]
            self.points[self.target_point] = [event.xdata, event.ydata]
            self.target_point = -1

            self.remove_inner_points()
            self.triangulate()

    def on_motion(self, event):
        if self.target_point >= 0:
            self.points[self.target_point] = [event.xdata, event.ydata]
            self.points = self.points[0:self.numPointsOriginal]
            self.triangulate()


    def add_point(self, event):
        """
        Function to add a point on a mouse click
        :param event:
        :return:
        """
        self.points = self.points[0:self.numPointsOriginal]
        self.points.append([event.xdata, event.ydata])
        self.numPointsOriginal += 1

        self.remove_inner_points()
        self.triangulate()


    def remove_inner_points(self):
        """
        Remove all original points that do not define the convex hull
        :return:
        """
        hull = ConvexHull(self.points[0:self.numPointsOriginal])
        vertices = hull.vertices

        #Tracker to make sure subsequent points are deleted correctly
        deleted_points = 0
        for i in range(self.numPointsOriginal):
            if not i in vertices:
                #Adjust i by the number of deleted points already
                self.points.pop(i - deleted_points)
                self.numPointsOriginal -= 1
                deleted_points += 1


    def is_inner(self, point_index):
        hull = ConvexHull(self.points[0:self.numPointsOriginal])
        vertices = hull.vertices

        return point_index in vertices


if __name__=="__main__":
    points = [[0, 0], [0.1, 0.3], [9.5, 4], [6, 0.5]]
    tris = DelaunayTris(points=points)
    plt.show()
