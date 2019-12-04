from sys import argv, stderr
import numpy as np
from VoronoiWall.structures import Diagram
from VoronoiWall.printable import fixOrdering, diagramToSTLMesh

def main():
    if len(argv) != 2:
        stderr.write("Usage: python {} input_file\n".format(argv[0]))
        stderr.write("\tinput_file: a set of coordinates defined in three dimensions\n")
        exit(1)
    try:
        input_file = open(argv[1], 'r')
    except FileNotFoundError:
        stderr.write("ERROR: File not found, {}\n".format(argv[1]))
        exit(1)

    with input_file:
        # skip the first line
        input_file.readline()
        # split each line at tabs, casting to 3 ints, storing as [[x, y, z], ...]
        points = [list(map(int, line.split('\t', 3))) for line in input_file]

    points_array = np.array(points)

    diagram = Diagram(points_array)

if __name__ == '__main__':
    main()