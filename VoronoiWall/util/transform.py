import sys
import os.path

class Facet:
    vertexList = []
    vectorList = []

    def __init__(self, vector, Ftype, v1, v2, v3):
        self.vector = vector
        self.Ftype = Ftype
        v1 = v1
        v2 = v2
        v3 = v3
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        Facet.vectorList.append(vector)
        Facet.vertexList.append(v1)
        Facet.vertexList.append(v2)
        Facet.vertexList.append(v3)

    def getVertices(self):
        listV = []
        listV.append(self.v1)
        listV.append(self.v2)
        listV.append(self.v3)
        return listV
    


if len(sys.argv) < 2:
    print("Set .stl file argumet")
    exit(0)

if not os.path.isfile(sys.argv[1]):
    print("File cannot be found")
    exit(0)

filename = sys.argv[1]

FacetList = []

with open(file=filename) as f:
    modelName = f.readline().split()[1]
    for line in f:
        line = line.split()
        Ftype = None
        vector = None
        if line[0] == 'facet':
            Ftype = line[1]
            vector = [line[2], line[3], line[4]]
        elif line[0] == 'endloop':
            continue
        elif line[0] == 'endfacet':
            continue
        elif line[0] == 'endsolid':
            break
        next(f)
        vLine = f.readline().split()
        v1 = [vLine[1], vLine[2], vLine[3]]
        vLine = f.readline().split()
        v2 = [vLine[1], vLine[2], vLine[3]]
        vLine = f.readline().split()
        v3 = [vLine[1], vLine[2], vLine[3]]            
        newFacet = Facet(vector, Ftype, v1, v2, v3)
        FacetList.append(newFacet)


print()
print(Facet.vectorList)




