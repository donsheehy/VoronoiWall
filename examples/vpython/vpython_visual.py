"""
Visualizes 3D Voronoi diagrams using VPython.

Currently limited to visualization of the closed surfaces of the diagram, not those extending to infinity. This might
be ideal for the purpose of 3D printing.
"""

from vpython import *
from scipy.spatial import Voronoi

# some random 3D points
points = [[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9], [7, 7, 4],
          [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
          [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
          [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]]

# make the canvas bigger
scene.width = 800
scene.height = 600
scene.resizable = True

# display the input points in white
for point in points:
    [x, y, z] = point
    sphere(pos=vector(x, y, z), radius=.1)

# compute the Voronoi
vor = Voronoi(points)

# display the Voronoi vertices in red
for point in vor.vertices:
    [x, y, z] = point
    sphere(pos=vector(x, y, z), radius=.1, color=color.red)

"""print("vor.vertices: " + str(vor.vertices))"""
"""print("vor.ridge_vertices: " + str(vor.ridge_vertices))"""
"""print("vor.regions: " + str(vor.regions))"""

# collect the indices of the vertices for all of the closed faces (those which don't extend to infinity)
completeShapeVertices = []
for indexList in vor.ridge_vertices:
    if -1 not in indexList:# and indexList != []:
        completeShapeVertices.append(indexList)

"""print("completeShapeVertices: " + str(completeShapeVertices))"""

# display the outlines of Voronoi faces (arbitrary shapes is hard in vpython, so just the outline)
for vertexIndexList in completeShapeVertices:
    # collect the Voronoi vertices that make up this closed shape
    pts = []
    for index in vertexIndexList:
        pts.append(vor.vertices[index])

    """print("shape: " + str(pts))"""

    # convert point-lists to vpython vector objects
    pts = list(map((lambda x: vector(x[0], x[1], x[2])), pts))
    # draw the outline of the closed face
    c = curve(pts)
    # loop back to close the face
    c.color = color.blue
    c.append(pts[0])

    """print("shape: " + str(pts))"""
