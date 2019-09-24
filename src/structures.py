import numpy as np

A = [2, 0]
B = [12, 0]
C = [10, 10]
D = [0, 10]
points2d = np.array([A, B, C, D])
points3d = np.array([[2, 0, 0], [12, 0, 0], [10, 10, 0], [0, 10, 0]])


def barycenter(pts):
    # split each of x, y, z, w, etc. into its own array
    splitAxes = np.squeeze(np.hsplit(pts, len(pts[0])))
    # take the average of x, y, etc. to find a midpoint
    return list(map(np.average, splitAxes))


print(barycenter(points2d))
print(barycenter(points3d))
