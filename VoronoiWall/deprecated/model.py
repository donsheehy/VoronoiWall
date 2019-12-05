import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull
from math import sqrt, ceil
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
import pickle
from util.math_utils import barycenter
from mpl_toolkits.mplot3d import Axes3D

import halfedge
import diagram

"""
Generates an STL file for 3D printing a Voronoi diagram. (eventually)
"""
class DelaunayTris:
    def __init__(self, points=[]):
        self.points = points
        self.cells = []  #Faces for halfedge data structure

        self.target_point = -1

        self.fig = plt.figure()
        self.ax = self.fig.gca(projection="3d")
        plt.subplots_adjust(bottom=0.2)
        # self.cid_press = self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        # self.cid_release = self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        # self.cid_motion = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

        axprepare = plt.axes([0.7, 0.05, 0.15, 0.075])
        axsave = plt.axes([0.5, 0.05, 0.15, 0.075])
        axopen = plt.axes([0.3, 0.05, 0.15, 0.075])
        axfile_name = plt.axes([0.05, 0.05, 0.15, 0.075])

        self.bprepare = Button(axprepare, 'Prepare stl')
        self.bprepare.on_clicked(self.prepare)

        self.bsave = Button(axsave, 'Save Points')
        self.bsave.on_clicked(self.save_points)

        self.bopen = Button(axopen, 'Open Points')
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
        self.subdivideFace(0)
        print("Preparing")
        # points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
        #       [7, 7, 4],
        #       [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
        #       [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
        #       [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])


        output = open(filePath, "w")
        output.write("solid Voronoi\n")
        faces = []
        for indexList in self.voronoi.ridge_vertices:
            if -1 not in indexList:
                face = []
                for index in indexList:
                    face.append(self.voronoi.vertices[index])
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
                print(trianglePoints)
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
        average_point = np.zeros((1,3))
        for point in points:
            average_point += point
        average_point /= len(points)
        return np.append(points, average_point, axis=0)

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
        #self.subdivideFace(self.points)
        self.plotVoronoi(self.points)


    def subdivideFace(self, face_index):
        """
        Given the index of a 3D face located in self.faces, subdivides the inner area into triangular
        shapes (necessary for 3D printing)
        :param points: a numpy array of input points; this array is modified in-place
        :return: array of tris
        """

        face = self.faces[face_index]
        verts = face.getVertices()

        print(verts)
        center = barycenter(verts)

        tris = []
        for i in range(len(face.halfedges)):
            cur_point = face.halfedges[i].vertex
            next_point = face.halfedges[i].next.vertex
            tris.append([cur_point.location, next_point.location, center])

        print(tris)
        return
        triangulation = Delaunay(points)

        trianglePoints = []
        for indexList in triangulation.simplices:
            for index in indexList:
                trianglePoints.append(points[index])

        points = trianglePoints
        return trianglePoints



    def plotVoronoi(self, points):
        """
        Display the subdivided face in 2D with matplotlib.

        Adapted from: https://stackoverflow.com/a/24952758
        :param points: points to plot, connecting them by simultaneously visualizing the Delaunary triangulation
        :return: unused
        """

        self.cells = []
        self.voronoi = Voronoi(points)

        vertices = []
        for i in range(len(self.voronoi.vertices)):
            location = [self.voronoi.vertices[i, 0],self.voronoi.vertices[i, 1],self.voronoi.vertices[i, 2]]
            vertices.append(halfedge.vertex(location=location))

        for r in range(len(self.voronoi.regions)):
            #self.regions.append(halfedge.face())
            cell = halfedge.cell()
            faces = []
            region = self.voronoi.regions[r]
            region_points = []
            region_point_indices = []

            for index in region:
                if index == -1 or index >= len(vertices):
                    break
                region_points.append(vertices[index].location)
                region_point_indices.append(vertices[index])

            if len(region_points) != len(region) or len(region) < 3:
                continue

            hull = ConvexHull(region_points)
            for simplex in hull.simplices:
                face = halfedge.face()

                edges = []
                for i in range(len(simplex)):
                    edges.append(halfedge.halfedge(vertex=vertices[simplex[i]], face=face))
                    if i > 0:
                        edges[i].previous = edges[i-1]  #Previous edge is edge before this one in the list
                        edges[i-1].next = edges[i]      #This edge is the next edge for the one before
                        edges[i-1].vertex.halfedge = edges[i]   #This edge is the outgoing edge for the last vertex
                    if i == len(simplex)-1:
                        edges[0].previous = edges[len(simplex)-1]
                        edges[len(simplex)-1].next = edges[0]
                        edges[len(simplex)-1].vertex.halfedge = edges[0]
                face.halfedges = edges
                faces.append(face)
            cell.faces = faces
            self.cells.append(cell)
            print(len(self.cells[0].faces))

        #     face = halfedge.face()
        #
        #     edges = []
        #     for i in range(len(region)):
        #         vertex_index = region[i]
        #         if vertex_index == -1:
        #             break
        #
        #         edges.append(halfedge.halfedge(vertex=vertices[vertex_index], face=face))
        #         if i > 0:
        #             edges[i].previous = edges[i-1]  #Previous edge is edge before this one in the list
        #             edges[i-1].next = edges[i]      #This edge is the next edge for the one before
        #             edges[i-1].vertex.halfedge = edges[i]   #This edge is the outgoing edge for the last vertex
        #         if i == len(region)-1:
        #             edges[0].previous = edges[len(region)-1]
        #             edges[len(region)-1].next = edges[0]
        #             edges[len(region)-1].vertex.halfedge = edges[0]
        #
        #     if len(edges) < len(region):
        #         continue
        #     else:
        #         face.halfedges = edges
        #     self.cells.append(face)
        #
        #
        # #Algorithm to fill in the opposite field for halfedges
        # for f in range(len(self.cells)):        #Loop through every face
        #     halfedges = self.faces[f].halfedges         #Get the halfedges that make up the face
        #     for edge in halfedges:              #Do this for every halfedge
        #         if edge.opposite == None:       #only if the edge doesn't have an opposite
        #             next = edge.next            #Get the next edge
        #             vertex_next_edges = edge.vertex.halfedges       #List of outgoing edges from next vertex
        #             for vertex_next_edge in vertex_next_edges:      #loop through these outgoing edges
        #                 if not(vertex_next_edge == next) and vertex_next_edge.vertex == edge.previous.vertex:       #if the outgoing edge isn't the next halfedge and it goes
        #                                                                                                             #into the same vertex that the current edge originates from
        #                     edge.opposite = vertex_next_edge        #Set the opposite of the current edge
        #                     vertex_next_edge.opposite = edge        #to the outgoing edge




        self.ax.clear()

        #PLOTTING FROM SCIPY VORONOI DATA STRUCTURE
        self.voronoi_points = self.ax.scatter(self.voronoi.points[:,0],self.voronoi.points[:,1],self.voronoi.points[:,2])
        self.ax.scatter(self.voronoi.vertices[:,0],self.voronoi.vertices[:,1],self.voronoi.vertices[:,2], 'r')

        for i in range(len(self.voronoi.ridge_vertices)):
            points = np.zeros((0,3))
            for j in range(len(self.voronoi.ridge_vertices[i])):
                index = self.voronoi.ridge_vertices[i][j]
                if index >= 0:
                    points = np.append(points, np.array([[self.voronoi.vertices[index, 0],self.voronoi.vertices[index, 1],self.voronoi.vertices[index, 2]]]), axis=0)
            if len(points) > 1:
                self.ax.plot(points[:, 0], points[:, 1], points[:, 2], 'b')

        #PLOTTING FROM HALFEDGE DATA STRUCTURE
        for cell in self.cells:
            print(len(cell.faces))
            for face in cell.faces:
                if len(face.halfedges)<1:
                    continue
                start_edge = face.halfedges[-1]
                cur_edge = start_edge.next
                while True:
                    locations = np.array([cur_edge.previous.vertex.location, cur_edge.vertex.location])
                    self.ax.plot(locations[:,0], locations[:,1], locations[:,2], 'go-')

                    if cur_edge == start_edge:
                        break
                    cur_edge = cur_edge.next

        plt.show()

    def on_press(self, event):
        print(event.inaxes)
        print(event.xdata)
        print(event.ydata)
        print(self.voronoi_points.axes)
        if event.inaxes != self.voronoi_points.axes:
            print("Not in axis!")
            return
        print(self.target_point)
        min_dist_squared = 1e10
        close_point = -1;
        for i in range(len(self.points)):
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
        if event.inaxes != self.voronoi_points.axes:
            print("Not in axis!")
            return
        if self.target_point == -1:
            self.add_point(event)
        else:
            self.points[self.target_point] = [event.xdata, event.ydata, 0]
            self.target_point = -1

            self.triangulate_vis()

    def on_motion(self, event):
        print(self.target_point)
        if event.inaxes != self.voronoi_points.axes:
            return
        if self.target_point >= 0:
            self.points[self.target_point] = [event.xdata, event.ydata, 0]
            self.triangulate_vis()


    def add_point(self, event):
        """
        Function to add a point on a mouse click
        :param event:
        :return:
        """

        self.points = np.append(points, [[event.xdata, event.ydata, 0]], axis=0)
        print(self.points)

        self.triangulate_vis()


if __name__=="__main__":
    points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
        [7, 7, 4],
        [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
        [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
        [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])
    tris = DelaunayTris(points=points)
    plt.show()
