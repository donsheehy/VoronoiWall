from scipy.spatial import Voronoi


class Point:
    def __init__(self, point, region_obj):
        self.point = point
        self.region = region_obj


class Region:
    def __init__(self, voronoi, vertices_i):
        self.vertices_i = vertices_i
        self.vertices = [voronoi.vertices[i] for i in vertices_i if i != -1]  # -1 means unbounded; don't index to -1
        self.facets = []


class Facet:
    def __init__(self, voronoi, vertices_i):
        self.voronoi = voronoi
        self.vertices_i = vertices_i
        self.vertices = [voronoi.vertices[i] for i in vertices_i if i != -1]

    def triangles(self):
        for i in range(1, len(self.vertices) - 1):
            yield Facet(self.voronoi, [self.vertices_i[0], self.vertices_i[i], self.vertices_i[i + 1]])


class Diagram:
    def __init__(self, points):
        self.voronoi = Voronoi(points)
        self.vertices = self.voronoi.vertices

        self.facets = [Facet(self.voronoi, facet) for facet in self.voronoi.ridge_vertices]
        self.regions = [Region(self.voronoi, region) for region in self.voronoi.regions]

        # match facets to the region containing them
        for facet_obj in self.facets:
            facet_vert_set = set(facet_obj.vertices_i)
            for region_obj in self.regions:
                region_vert_set = set(region_obj.vertices_i)
                if facet_vert_set <= region_vert_set:
                    region_obj.facets.append(facet_obj)

        self.points = [Point(self.voronoi.points[i], self.regions[region_i]) for i, region_i in
                       enumerate(self.voronoi.point_region)]

        self.bounded_facets = [facet for facet in self.facets if -1 not in facet.vertices_i]
        self.bounded_regions = [region for region in self.regions if
                                -1 not in region.vertices_i and len(region.vertices_i) > 3]


if __name__ == '__main__':
    import numpy as np

    points = np.array([[0, 0, 3], [0, 3, 3], [3, 0, 3], [3, 3, 3], [0, 0, 0], [0, 3, 0], [3, 0, 0], [3, 3, 0],
                       [1, 2, 2], [2, 2, 2], [1, 2, 1], [2, 2, 1], [1, 1, 2], [2, 1, 2], [1, 1, 1], [2, 1, 1]])
    diagram = Diagram(points)

    # A list of Point objects (the input points)
    point_list = diagram.points

    # A Point object
    point_obj = point_list[0]

    # The point's location as [x, y, z]
    point_coord = point_obj.point

    # A Region object produced from the Point above
    region_obj = point_obj.region

    # A list of the Region's Facet objects
    facet_list = region_obj.facets

    # A Facet object
    facet_obj = facet_list[0]

    # A list of the Facet's vertices
    vertex_list = facet_obj.vertices

    # A vertex as [x, y, z]
    vertex = vertex_list[0]

    print(vertex)
