import numpy as np

def barycenter(pts):
    # split each of x, y, z, w, etc. into its own array
    splitAxes = np.squeeze(np.hsplit(pts, len(pts[0])))
    # take the average of x, y, etc. to find a midpoint
    return list(map(np.average, splitAxes))  # type: ndarray