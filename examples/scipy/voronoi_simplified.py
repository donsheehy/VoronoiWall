"""
If you've looked at vornoi_verbose.py, this accomplishes the same task just with a lot less code.
I found the verbose version helpful, though, since I bet the data it exposes will be useful in creating 3D model files.

See: https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.Voronoi.html
See: https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d


def main():
    # this is the input array of points
    points = np.array([[-1, 0], [1, 0], [0, 1], [0, -1], [0, 0]])

    # use scipy to compute the Voronoi
    vor = Voronoi(points)

    # use the built-in scipy plotter to generate a matplotlib plot
    voronoi_plot_2d(vor)

    # view the finalized plot
    plt.show()


main()
