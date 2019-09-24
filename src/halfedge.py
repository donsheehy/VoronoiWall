#Each vertex references one outgoing halfedge.
# Each face references one of the halfedges bounding it.
# Each halfedge provides a handle to
    # the vertex it points to
    # the face it belongs to
    # the next halfedge inside the face (ordered counter-clockwise)
    # the opposite halfedge
    # the previous halfedge in the face

class face:
    halfedges = []

    def __init__(self, halfedges=[]):
        self.halfedges = halfedges


class halfedge:
    vertex = None
    face = None
    next = None
    opposite = None
    previous = None

    def __init__(self, vertex=None, face=None, next=None, opposite=None, previous=None):
        self.vertex = vertex
        self.face = face
        self.next = next
        self.opposite = opposite
        self.previous = previous

class vertex:
    location = []
    halfedges = []

    def __init__(self, location, halfedges=[]):
        self.location = location
        self.halfedges = halfedges



