from sys import argv, stderr
import numpy as np
from VoronoiWall.structures import Diagram
from VoronoiWall.printable import diagramToSTLMesh, plotSTLMesh, readPointsFile


def usage():
    if len(argv) < 2 or len(argv) > 3:
        stderr.write("Usage: python {} input_file [output_file]\n".format(argv[0]))
        stderr.write("\tinput_file: a set of coordinates defined in three dimensions\n")
        stderr.write("\toutput_file: optionally, provide a name for the output STL file\n")
        exit(1)


def main():
    usage()

    # Parse the coordinate points from the input file.
    try:
        points_array = readPointsFile(argv[1])
    except FileNotFoundError:
        stderr.write("ERROR: Input file not found, {}\n".format(argv[1]))
        exit(1)

    # Construct the Voronoi Diagram data structure.
    diagram = Diagram(points_array)

    # Warn the user if all regions are unbounded (can't be plotted).
    if len(diagram.bounded_regions) == 0:
        print("WARNING: The Voronoi diagram of the given input has no bounded regions.")

    # Generate an STL representation of the Diagram.
    mesh = diagramToSTLMesh(diagram)

    # Save it if an output file name was provided.
    if len(argv) == 3:
        try:
            mesh.save(argv[2])
        except FileNotFoundError:
            stderr.write("ERROR: Output file not found, {}\n".format(argv[1]))
            exit(1)

    # Display the STL mesh on the screen.
    plotSTLMesh(mesh)


if __name__ == '__main__':
    main()
