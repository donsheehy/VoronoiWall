from vpython import *
from scipy.spatial import Voronoi

points = [[-2, -1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1], [1, 1, -1], [-1, 1, 1], [1, -1, 1], [1, 1, 1]]

scene.width = 1280
scene.height = 720

# display the 8 points in 3D

for point in points:
    [x, y, z] = point
    sphere(pos=vector(x, y, z), radius=.1)

vor = Voronoi(points)
# print(vor.vertices)

for point in vor.vertices:
    [x, y, z] = point
    sphere(pos=vector(x, y, z), radius=.1, color=color.red)

print(vor.ridge_vertices)
voronoiVertices = vor.vertices
completeFace = vor.ridge_vertices[2]
vertices = []
for idx in completeFace:
    vertices.append(vec(voronoiVertices[idx][0], voronoiVertices[idx][1], voronoiVertices[idx][2]))

t = triangle(
    v0=vertex(
        pos=vertices[0]),
    v1=
    vertex(
        pos=vertices[1]),
    v2=
    vertex(
        pos=vertices[2]))
