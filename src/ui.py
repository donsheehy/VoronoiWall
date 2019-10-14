import matplotlib.pyplot as plt
import numpy as np
''' type annotations... '''
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Patch3DCollection
''' ... type annotations '''

points = np.array([[6, 4, 2], [9, 5, 8], [9, 1, 9], [8, 9, 1], [3, 8, 8], [2, 6, 2], [8, 2, 10], [3, 6, 1], [9, 8, 9],
                   [7, 7, 4],
                   [2, 10, 5], [4, 3, 10], [5, 3, 9], [4, 7, 4], [3, 6, 7], [7, 4, 3], [6, 4, 9], [5, 8, 4], [2, 9, 10],
                   [7, 8, 6], [9, 2, 7], [6, 10, 7], [9, 9, 3], [2, 9, 4], [5, 9, 6], [4, 8, 9], [9, 1, 2], [6, 9, 1],
                   [10, 6, 5], [1, 9, 9], [2, 1, 3], [10, 1, 5], [4, 10, 2]])


def displayRegion(region):
    pass


def displaySimplex(simplex):
    pass


def displayVertex(vertex):
    pass


def setSelectorMode(fig, mode):
    if mode == 'point':
        def onPickPoint(event):
            print("Clicked point: {}".format(points[event.ind]))
        eventHandler = onPickPoint
    elif mode == 'region':
        pass
    elif mode == 'simplex':
        pass
    else:
        raise ValueError('Valid selector modes: point, region, simplex')
    fig.canvas.mpl_connect('pick_event', eventHandler)


def update():
    plt.show()


def main():
    # begin a new mpl figure, set its title
    fig = plt.figure('Voronoi Editor')  # type: Figure
    # set the title at the top of the viewer
    fig.suptitle('Voronoi Wall')
    # make a new 3D plot
    ax = fig.add_subplot(111, projection='3d')  # type: Axes3D
    # set the title of the plot
    ax.set_title('Voronoi Wall')

    # plot the input points
    pts = ax.scatter(*zip(*points), picker=0)  # type: Patch3DCollection

    # an idea for toggling between region, point, simplex, etc. selection
    # (selection being clicking on one and being able to manipulate it)
    setSelectorMode(fig, 'point')

    # redraw the plot
    update()


main()
