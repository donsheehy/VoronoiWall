from stl import mesh, stl
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from VoronoiWall.structures import Diagram


def readPointsFile(fileName):
    with open(fileName, 'r') as input_file:
        input_file.readline()
        # split each line at tabs, casting to 3 floats, storing as [[x, y, z], ...]
        return np.array([list(map(float, line.split('\t', 3))) for line in input_file])


def plotSTLMesh(mesh):
    # Create a new plot.
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)

    # Load the STL files and add the vectors to the plot.
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh.vectors))

    # Auto scale to the mesh size.
    scale = mesh.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen.
    pyplot.show()


def diagramToSTLMesh(diagram):
    # Flip any faces facing inwards to now face outwards.
    fixOrdering(diagram)

    vertices = diagram.vertices

    # Collect the facets of bounded regions.
    bounded_region_facets = [y for x in [region.facets for region in diagram.bounded_regions] for y in x]

    # Collect the vertices of triangulated facets.
    tri_facets_i = []
    for facet in bounded_region_facets:
        for tri in facet.triangles():
            tri_facets_i.append(tri.vertices_i)
    tri_facets_i = np.array(tri_facets_i)

    # Generate the Mesh object.
    m = mesh.Mesh(np.zeros(tri_facets_i.shape[0], dtype=mesh.Mesh.dtype))
    for i, facet_i in enumerate(tri_facets_i):
        for j in range(len(facet_i)):
            m.vectors[i][j] = vertices[facet_i[j], :]

    return m


def fixOrdering(diagram):
    # for every bounded facet
    for facet in diagram.bounded_facets:
        '''
        Every bounded facet is going to be in 2+ regions since even the outer facets are part of both a bounded
        and unbounded region. That means we can't just find the bounded facets that are part of 1 region and only check
        those. So we have to consider every bounded facet, and therefore there is extra work being done to 
        flip inner facets here.
        '''
        # find the bounded region it's a part of (or the last one if it's an internal face. how could those be skipped?)
        input_point = None
        for region in facet.regions:
            if -1 not in region.vertices_i:
                input_point = region.point

        # if this bounded facet was part of an unbounded region
        if input_point is None:
            continue

        normal_vector = np.cross(facet.vertices[1] - facet.vertices[0], facet.vertices[2] - facet.vertices[1])
        point_direction_vector = facet.vertices[0] - input_point.point
        costheta = np.dot(normal_vector, point_direction_vector) / (
                np.linalg.norm(normal_vector) * np.linalg.norm(point_direction_vector))

        '''
        Check if the angle between the (a facet vertex -> input point) vector and the (facet normal) vector is 
        greater than 90 degrees. In other words, the facet is like a plane, and in this case, the input point is on
        the opposite side of the plane than the normal vector points. I would have thought this was the correct,
        no-change-needed case, but in fact, reversing those facets seems to fix STL orientation.
        '''
        if (costheta < 0.0):  # any angle pi/2 < theta < 3pi/2, or 90deg < theta < 270deg, has a negative cos(theta)
            # will this check work for extreme, 270-360 degree angles? I'm not sure.
            facet.vertices_i = facet.vertices_i[::-1]
            facet.vertices = facet.vertices[::-1]


if __name__ == '__main__':
    import numpy as np
    from structures import Diagram

    points = np.array([[0, 0, 3], [0, 3, 3], [3, 0, 3], [3, 3, 3], [0, 0, 0], [0, 3, 0], [3, 0, 0], [3, 3, 0],
                       [1, 2, 2], [2, 2, 2], [1, 2, 1], [2, 2, 1], [1, 1, 2], [2, 1, 2], [1, 1, 1], [2, 1, 1]])
    points = np.array(
        [[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
         [7, 7, 4],
         [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
         [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
         [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])
    # Generate the Diagram data structure
    diagram = Diagram(points)

    # Generate a Mesh object with correct normals
    m = diagramToSTLMesh(diagram)

    # Write the Mesh to an STL file
    m.save('test.stl', mode=stl.Mode.ASCII)
