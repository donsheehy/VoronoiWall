from stl import mesh, stl
import numpy as np


def diagramToSTLMesh(diagram):
    vertices = diagram.vertices

    # consider making this act on the facets of bounded regions (some facets might be bounded, but not their region)
    tris = [[tri.vertices_i for tri in facet.triangles()] for facet in diagram.bounded_facets]
    tri_facets_i = np.array(tris).squeeze()

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

    # Generate the Diagram data structure
    diagram = Diagram(points)

    # Reorient any incorrectly ordered facets (normal vectors would have been the wrong way)
    fixOrdering(diagram)

    # Generate a Mesh object
    m = diagramToSTLMesh(diagram)

    # Write the Mesh to an STL file
    m.save('test.stl', mode=stl.Mode.ASCII)
