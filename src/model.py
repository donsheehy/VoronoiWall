import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial import Voronoi
from math import sqrt, ceil
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
import pickle

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
        plt.subplots_adjust(bottom=0.2)
        self.cid_press = self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        axprepare = plt.axes([0.7, 0.05, 0.15, 0.075])
        axsave = plt.axes([0.4, 0.05, 0.15, 0.075])
        axopen = plt.axes([0.2, 0.05, 0.15, 0.075])
        axfile_name = plt.axes([0.05, 0.05, 0.15, 0.075])

        self.bprepare = Button(axprepare, 'Prepare stl')
        self.bprepare.on_clicked(self.prepare)

        self.bsave = Button(axsave, 'Save Points to file')
        self.bsave.on_clicked(self.save_points)

        self.bopen = Button(axopen, 'Open Points File')
        self.bopen.on_clicked(self.open_points)

        self.points_file = "points.p"
        self.textbox_file_name = TextBox(axfile_name, "", initial="points.p")
        self.textbox_file_name.on_text_change(self.update_file_name)

        self.triangulate_vis()

    def load_points(self, points):
        """Function to load and display an array of points"""
        if len(points) == 0:
            print("No points provided!")
            return
        print("Loading " + str(len(points)) + " " + str(len(points[0])) + " dimensional points")
        self.points = points
        self.triangulate_vis()

    def open_points(self, event):
        file_points = pickle.load(open(self.points_file, "rb"))
        self.load_points(file_points)

    def save_points(self, event):
        pickle.dump(self.points, open(self.points_file, "wb+"))


    def update_file_name(self, text):
        self.points_file = text

    def prepare(self, event, filePath="out.stl"):
        """
        Prepares a 3D printable model of the given Voronoi diagram.

        :param voronoi: computed Voronoi attributes
        :return: unused
        """
        print("Preparing")
        points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
              [7, 7, 4],
              [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
              [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
              [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])
        voronoi = Voronoi(points)


        output = open(filePath, "w")
        output.write("solid Voronoi\n")
        faces = []
        for indexList in voronoi.ridge_vertices:
            if -1 not in indexList:
                face = []
                for index in indexList:
                    face.append(voronoi.vertices[index])
                faces.append(np.asarray(face))
        # I'm thinking order could be important for the triangle vertices and is being lost?
        for face in faces:
            triangles = self.triangulate(face)
            # compute a normal vector for this face
            normal = np.cross(face[1] - face[0], face[2] - face[1])
            # process points in batches of 3 (points of a triangle)
            for i in range(0, len(triangles), 3):
                # begin a new STL triangle
                output.write("facet normal {} {} {}\n".format(normal[0], normal[1], normal[2]))
                output.write("outer loop\n")
                trianglePoints = triangles[i:i + 3]
                for j in range(0, 3):
                    output.write("vertex {} {} {}\n".format(trianglePoints[j][0], trianglePoints[j][1], trianglePoints[j][2]))
            output.write("endloop\nendfacet\n")

        output.write("endsolid Voronoi\n")

    def triangulate(self, points):
        """
        Splits a 3D planar facet into triangles.
        :param points: vertex coordinates for a planar face in 3D
        :return: vertices of the divided plane
        """
        # move all points by this much so the shape has to be touching the origin
        offset = points[0]
        # get a normal vector to the plane for the rotation matrix
        normalVector = np.cross(points[1] - points[0], points[2] - points[1])
        # translate to the origin, then rotate from the current plane onto XY-plane with normal <0,0,-1>
        points = self.rotateToPlane(points, normalVector, np.array([0, 0, 1]), False, offset)
        # reduce the (N, 3) array to an (N, 2) array because QHull will complain about operating on planes in 3D
        points = self.chopOffThirdDimension(points)
        # use the Delaunay triangulation to divide this plane into triangles (for 3D printing ability)
        points = self.subdivideOnce(points)
        # expand the (N, 2) array back to an (N, 3) array by adding a zeroed column
        points = self.addEmptyThirdDimension(points)
        # rotate back to the plane we were in originally, then translate back to the original location
        points = self.rotateToPlane(points, np.array([0, 0, 1]), normalVector, True, offset)
        return points

    def subdivideOnce(self, points):
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

    def chopOffThirdDimension(self, npArrayOf3DPoints):
        return np.delete(npArrayOf3DPoints, 2, 1)


    def addEmptyThirdDimension(self, npArrayOf2DPoints):
        return np.insert(npArrayOf2DPoints, 2, values=0, axis=1)


    def rotateToPlane(self, points, normalVectorOriginal, normalVectorNew, isAtOrigin=True, offset=np.array([0, 0, 0])):
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
        self.subdivideFace(self.points)
        self.plotDelaunayTriangles(self.points, self.numPointsOriginal)


    def subdivideFace(self, points):
        """
        Given the vertices of a 2D shape located in the XY coordinate plane, subdivides the inner area into triangular
        shapes (necessary for 3D printing) using the Delaunay triangulation.
        :param points: a numpy array of input points; this array is modified in-place
        :return: unused
        """

        triangulation = Delaunay(points)

        trianglePoints = []
        for indexList in triangulation.simplices:
            for index in indexList:
                trianglePoints.append(points[index])

        points = trianglePoints
        return trianglePoints



    def plotDelaunayTriangles(self, points, numOriginalPoints):
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

        #Not sure why this fixes a stack overflow error but it does
        for vert in self.vertices:
            vert.figure.canvas.draw()

    def on_press(self, event):
        print("press: " + str(event.xdata) + " " + str(event.ydata))
        if event.inaxes == None:
            print("Not in axis!")
            return
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
        if event.inaxes == None:
            print("Not in axis!")
            return
        if self.target_point == -1:
            self.add_point(event)
        else:
            self.points = self.points[0:self.numPointsOriginal]
            self.points[self.target_point] = [event.xdata, event.ydata]
            self.target_point = -1

            self.triangulate_vis()

    def on_motion(self, event):
        if event.inaxes == None:
            return
        if self.target_point >= 0:
            self.points[self.target_point] = [event.xdata, event.ydata]
            self.points = self.points[0:self.numPointsOriginal]
            self.triangulate_vis()


    def add_point(self, event):
        """
        Function to add a point on a mouse click
        :param event:
        :return:
        """
        if event.inaxes == None:
            print("Not in axes!")
            return
        self.points = self.points[0:self.numPointsOriginal]
        self.points.append([event.xdata, event.ydata])
        self.numPointsOriginal += 1

        self.triangulate_vis()


if __name__=="__main__":
    points = [[0, 0], [0.1, 0.3], [9.5, 4], [6, 0.5]]
    tris = DelaunayTris(points=points)
    plt.show()
