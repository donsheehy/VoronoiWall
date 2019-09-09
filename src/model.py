import vpython as v
import numpy as np

"""
Generates an STL file for 3D printing a Voronoi diagram. 
"""


# def displayShape(numpyArray3D):
#     # display the input points in white
#     vpyPoints = []
#     for point in numpyArray3D:
#         [x, y, z] = point
#         vpyVector = v.vec(x, y, z)
#         v.sphere(pos=vpyVector, radius=.2)
#         vpyPoints.append(vpyVector)
#
#     c = v.curve(vpyPoints)
#     c.append(vpyPoints[0])


def prepare(voronoi, filePath):
    """
    Prepares a 3D printable model of the given Voronoi diagram.

    :param voronoi: computed Voronoi attributes
    :return: unused
    """

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
        triangles = triangulate(face)
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

def triangulate(points):
    # move all points by this much so the shape has to be touching the origin
    offset = points[0]
    # get a normal vector to the plane for the rotation matrix
    normalVector = np.cross(points[1] - points[0], points[2] - points[1])
    # translate to the origin, then rotate from the current plane onto XY-plane with normal <0,0,-1>
    points = __rotateToPlane(points, normalVector, np.array([0, 0, 1]), False, offset)
    # reduce the (N, 3) array to an (N, 2) array because QHull will complain about operating on planes in 3D
    points = __chopOffThirdDimension(points)
    # use the Delaunay triangulation to divide this plane into triangles (for 3D printing ability)
    points = __subdivideFace(points)
    # expand the (N, 2) array back to an (N, 3) array by adding a zeroed column
    points = __addEmptyThirdDimension(points)
    # rotate back to the plane we were in originally, then translate back to the original location
    points = __rotateToPlane(points, np.array([0, 0, 1]), normalVector, True, offset)
    return points


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


def __subdivideFace(points):
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


def __plotDelaunayTriangles(points):
    """
    Display a subdivided face in 2D with matplotlib.

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
    points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
              [7, 7, 4],
              [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
              [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
              [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])

    from scipy.spatial import Voronoi

    vor = Voronoi(points)

    prepare(vor, "out.stl")


main()
