from matplotlib import pyplot
from mpl_toolkits import mplot3d

import numpy as np
from scipy.spatial import Voronoi


# class CircularList(list):
#     """
#     Extends the built-in list class.
#     """
#
#     def cycle(self):
#         """
#         Generator that loops from the end of the list back to the beginning.
#         Use next() to increment the generator.
#         :return: essentially an infinite iterator
#         """
#         pos = 0
#         while True:
#             yield self[pos]
#             pos = (pos + 1) % len(self)
#
#     def __hash__(self):
#         tot = 0
#         for item in self:
#             tot += hash(item)
#         return tot


# class Region:
#     def __init__(self, vertices_i):
#         self.vertices_i = vertices_i
#         self.facets_i = []
#
#     def __str__(self):
#         return str(self.vertices_i)
#
#     def __repr__(self):
#         return self.__str__()


class Diagram:
    def __init__(self, points):
        self.vor = Voronoi(points)
        self.bounded_facets_i = [face for face in self.vor.ridge_vertices if -1 not in face]
        self.bounded_regions_i = [region for region in self.vor.regions if -1 not in region and len(region) > 3]
        # regions = []
        # for region_i in bounded_regions_i:
        #     region = Region(region_i)
        #     region_set = set(region_i)
        #     for facet_i in bounded_facets_i:
        #         facet_set = set(facet_i)
        #         if facet_set <= region_set:  # if the facet is a subset of the region
        #             region.facets_i.append(facet_i)
        #     regions.append(region)

        class Plot:
            def __init__(self):
                self.selected_facet_idx = 0
                self.selected_facet = None

        self.plot = Plot()

    def view(self):
        facets_verts = []
        for facet_i in self.bounded_facets_i:
            verts = np.empty((0, 3), int)
            for coord in (self.vor.vertices[i] for i in facet_i):
                verts = np.concatenate((verts, [coord]), axis=0)
            facets_verts.append(verts)

        self.__plot(facets_verts)

    def __plot(self, facets):
        # Select the plot
        figure = pyplot.figure()
        axes = mplot3d.Axes3D(figure)
        axes.set_axis_off()

        facets_collection = mplot3d.art3d.Poly3DCollection(facets)  # type: mplot3d.art3d.Poly3DCollection
        facets_collection.set_alpha(0.5)
        facets_collection.set_edgecolor('k')
        axes.add_collection3d(facets_collection)

        def key_pressed(event):
            # highlight the next facet
            self.plot.selected_facet_idx = (self.plot.selected_facet_idx + 1) % len(facets)
            if self.plot.selected_facet is not None:
                # remove the old one so it's not highlighted
                axes.collections.remove(self.plot.selected_facet)

            self.plot.selected_facet = mplot3d.art3d.Poly3DCollection([facets[self.plot.selected_facet_idx]])

            # display the new facet
            self.plot.selected_facet.set_sort_zpos(20)
            self.plot.selected_facet.set_edgecolor('k')
            self.plot.selected_facet.set_facecolor('red')
            axes.add_collection3d(self.plot.selected_facet)

        # press any key to advance the iterator
        figure.canvas.mpl_connect('key_press_event', key_pressed)

        # Show the plot
        pyplot.show()
        return figure


if __name__ == '__main__':
    points = np.array([[0, 0, 3], [0, 3, 3], [3, 0, 3], [3, 3, 3], [0, 0, 0], [0, 3, 0], [3, 0, 0], [3, 3, 0],
                       [1, 2, 2], [2, 2, 2], [1, 2, 1], [2, 2, 1], [1, 1, 2], [2, 1, 2], [1, 1, 1], [2, 1, 1]])
    diagram = Diagram(points)
    diagram.view()
