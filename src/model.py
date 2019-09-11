import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from math import sqrt, ceil

"""
Generates an STL file for 3D printing a Voronoi diagram.
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

    def prepare(voronoi, filePath):
        """
        Prepares a 3D printable model of the given Voronoi diagram.

        :param voronoi: computed Voronoi attributes
        :return: unused
        """

        output = open(filePath, "w")
        # a name is optional, but "solid " is not. this solid is named "Voronoi"
        output.write("solid Voronoi\n")
        # collect the faces which don't extend out to infinity
        faces = []
        for indexList in voronoi.ridge_vertices:
            if -1 not in indexList:
                face = []
                for index in indexList:
                    face.append(voronoi.vertices[index])
                faces.append(np.asarray(face))
        # split every surface into triangles and write those triangles to file in STL format
        for face in faces:
            triangles = triangulate(face)
            # compute a normal vector for this face
            normal = np.cross(face[1] - face[0], face[2] - face[1])
            # process points in batches of 3 (a triangle)
            for i in range(0, len(triangles), 3):
                # begin a new STL triangle
                output.write("facet normal {} {} {}\n".format(normal[0], normal[1], normal[2]))
                output.write("outer loop\n")
                trianglePoints = triangles[i:i + 3]
                for j in range(0, 3):
                    output.write(
                        "vertex {} {} {}\n".format(trianglePoints[j][0], trianglePoints[j][1], trianglePoints[j][2]))
            output.write("endloop\nendfacet\n")
        # end the STL file
        output.write("endsolid Voronoi\n")

    def triangulate(points):
        """
        Splits a 3D planar facet into triangles.
        :param points: vertex coordinates for a planar face in 3D
        :return: vertices of the divided plane
        """
        # move all points by this much so the shape has to be touching the origin
        offset = points[0]
        # get a normal vector to the plane for the rotation matrix
        normalVector = np.cross(points[1] - points[0], points[2] - points[1])
        # translate to the origin, then rotate from the current plane onto XY-plane with normal <0,0,1>
        points = __rotateToPlane(points, normalVector, np.array([0, 0, 1]), False, offset)
        # reduce the (N, 3) array to an (N, 2) array because QHull will complain about operating on planes in 3D
        points = __chopOffThirdDimension(points)
        # use the Delaunay triangulation to divide this plane into triangles (for 3D printing ability)
        points = __subdivideOnce(points)
        # expand the (N, 2) array back to an (N, 3) array by adding a zeroed column
        points = __addEmptyThirdDimension(points)
        # rotate back to the plane we were in originally, then translate back to the original location
        points = __rotateToPlane(points, np.array([0, 0, 1]), normalVector, True, offset)
        return points

    def __subdivideOnce(points):
        """
        Given the vertices of a 2D shape located in the XY coordinate plane, subdivides the inner area into triangular
        shapes (necessary for 3D printing) using the Delaunay triangulation.

        :param points: a numpy array of input points; this array is modified in-place
        :return: unused
        """

        from scipy.spatial import Delaunay

        triangulation = Delaunay(points)
        trianglePoints = []
        for indexList in triangulation.simplices:
            for index in indexList:
                trianglePoints.append(points[index])
        return trianglePoints

    def __chopOffThirdDimension(npArrayOf3DPoints):
        return np.delete(npArrayOf3DPoints, 2, 1)

    def __addEmptyThirdDimension(npArrayOf2DPoints):
        return np.insert(npArrayOf2DPoints, 2, values=0, axis=1)

    def __rotateToPlane(points, normalVectorOriginal, normalVectorNew, isAtOrigin=True, offset=np.array([0, 0, 0])):
        """
        Rotates a shape defined by its vertices about a defined axis. Useful for putting a planar shape located in 3D into
        a coordinate plane or restoring it to its original location in 3D space.

        :param points:                  list of points to rotate about an axis
        :param normalVectorOriginal:    vector (as numpy array) which is normal to the original plane
        :param normalVectorNew:         vector (as numpy array) which is normal to the desired plane
        :param isAtOrigin:              True if the shape defined by the given points is located at the origin
        :param offset:                  a vector (as numpy array) offset which is either subtracted from the given points or
                                            added to the resulting points if isAtOrigin is False or True
        :return: new numpy array of points rotated about the defined axis
        """
        from math import sqrt
        if not isAtOrigin:
            # translate points by the offset, typically moving the shape to the origin
            points = points - offset
        M = normalVectorOriginal
        N = normalVectorNew
        # compute costheta using the geometric dot product
        costheta = np.dot(M, N) / (np.linalg.norm(M) * np.linalg.norm(N))
        # cross the two axis vectors, make the result a unit vector
        mncross = np.cross(M, N)
        axis = mncross / np.linalg.norm(mncross)
        # shorten variable names (s = sintheta)
        c = costheta
        s = sqrt(1 - c * c)
        C = 1 - c
        [x, y, z] = axis

        # rotation matrix via https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
        rmat = np.array([[x * x * C + c, x * y * C - z * s, x * z * C + y * s],
                         [y * x * C + z * s, y * y * C + c, y * z * C - x * s],
                         [z * x * C - y * s, z * y * C + x * s, z * z * C + c]])

        if isAtOrigin:
            # rotate all of the points and then move the shape back to its original location
            return list(map(lambda point: np.dot(rmat, point) + offset, points))
        else:
            # rotate all of the points; will only work correctly if the shape is at the origin
            return list(map(lambda point: np.dot(rmat, point), points))

    def triangulate_vis(self):
        self.__subdivideFace(self.points, 0.5)
        self.__plotDelaunayTriangles(self.points, self.numPointsOriginal)

    def __subdivideFace(self, points, maxLength, depth=0):
        """
        Given the vertices of a 2D shape, subdivides the inner area into triangular shapes (necessary for 3D printing) using
        the Delaunay triangulation recursively until no triangle has an edge length larger than the specified maximum.

        Adapted from: https://stackoverflow.com/a/24952758

        :param points: a Python array of input points; this array is modified in-place
        :param maxLength: keep shrinking regions until all edges are shorter than this
        :return: unused
        """

        # Randomly run into error where subdivision enters an infinite loop
        # have a depth counter to just exit if subdivision goes past the arbitrary level
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
            self.__subdivideFace(points, maxLength, depth=depth + 1)

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
        self.vertices = self.ax.plot(npPoints[0:numOriginalPoints, 0], npPoints[0:numOriginalPoints, 1], 'or',
                                     markersize=6)
        self.ax.plot(npPoints[numOriginalPoints:-1, 0], npPoints[numOriginalPoints:-1, 1], 'o', markersize=4)

        if self.target_point >= 0 and not self.is_inner(self.target_point):
            self.ax.plot(npPoints[self.target_point, 0], npPoints[self.target_point, 1], 'ok', markersize=6)

        # Not sure why this fixes a stack overflow error but it does
        for vert in self.vertices:
            vert.figure.canvas.draw()

    def on_press(self, event):
        print("press: " + str(event.xdata) + " " + str(event.ydata))
        print(self.target_point)
        min_dist_squared = 1e10
        close_point = -1;
        for i in range(self.numPointsOriginal):
            dist = (event.xdata - self.points[i][0]) * (event.xdata - self.points[i][0]) + (
                    event.ydata - self.points[i][1]) * (event.ydata - self.points[i][1])
            if dist < min_dist_squared:
                min_dist_squared = dist
                close_point = i
        if sqrt(min_dist_squared) < 0.1:
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

        # Tracker to make sure subsequent points are deleted correctly
        deleted_points = 0
        for i in range(self.numPointsOriginal):
            if not i in vertices:
                # Adjust i by the number of deleted points already
                self.points.pop(i - deleted_points)
                self.numPointsOriginal -= 1
                deleted_points += 1

    def is_inner(self, point_index):
        hull = ConvexHull(self.points[0:self.numPointsOriginal])
        vertices = hull.vertices

        return point_index in vertices


if __name__ == "__main__":
    points = [[0, 0], [0.1, 0.3], [9.5, 4], [6, 0.5]]
    tris = DelaunayTris(points=points)
    plt.show()
