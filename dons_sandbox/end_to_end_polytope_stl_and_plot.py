import numpy as np
from stl_utils import polytope, stl_plot

# Define the 6 vertices of an octahedron.
vertices = np.array([\
    [0, 0, 0.7],
    [0, 1.1, 0],
    [1, 0, 0],
    [0, 0, -1.5],
    [0, -1.2, 0.5],
    [-1, 0, 0]
    ])
# Define the 8 triangular faces.
faces = np.array([\
    [0,1,2],
    [0,5,1],
    [1,5,3],
    [1,3,2],
    [2,3,4],
    [0,2,4],
    [0,4,5],
    [3,5,4],
    ])

p = polytope(vertices, faces)

stl_plot(p)
