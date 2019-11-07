from stl import mesh, stl


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


if __name__ == '__main__':
    import numpy as np
    from structures import Diagram

    points = np.array([[0, 0, 3], [0, 3, 3], [3, 0, 3], [3, 3, 3], [0, 0, 0], [0, 3, 0], [3, 0, 0], [3, 3, 0],
                       [1, 2, 2], [2, 2, 2], [1, 2, 1], [2, 2, 1], [1, 1, 2], [2, 1, 2], [1, 1, 1], [2, 1, 1]])

    diagram = Diagram(points)

    vertices = diagram.vertices

    tris = [[tri.vertices_i for tri in facet.triangles()] for facet in diagram.bounded_facets]
    tri_facets_i = np.array(tris).squeeze()

    m = mesh.Mesh(np.zeros(tri_facets_i.shape[0], dtype=mesh.Mesh.dtype))
    for i, facet_i in enumerate(tri_facets_i):
        for j in range(len(facet_i)):
            m.vectors[i][j] = vertices[facet_i[j], :]
    m.update_normals()
    print(m.normals[3])
    print(m.vectors[3])
    m.vectors[3] = m.vectors[3][::-1]
    m.update_normals()
    print(m.normals[3])
    print(m.vectors[3])

    # TODO: fix normals... (I have an idea)

    m.save('test.stl', update_normals=False, mode=stl.Mode.ASCII)
