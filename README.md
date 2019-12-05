# VoronoiWall
It's what you think it is.

## Description



## Requirements

The following Python packages must be installed to use this utility.

- numpy
- numpy-stl
- scipy
- matplotlib

## Installation

1. Clone the repository.

   ```bash
   git clone https://github.com/donsheehy/VoronoiWall.git
   ```

2. Move into the project directory.

   ```bash
   cd VoronoiWall
   ```

3. Build and install the package.

   ```bash
   python setup.py install
   ```

4. Import and use the project like any Python package.

   ```python
   import VoronoiWall as vw
   ```

## Examples

The example driver program allows testing various functionality of the library.

Given a set of input points, the example driver will:

- Generate a Voronoi diagram data structure
- Generate an STL mesh object encoding the Voronoi vertices
- Display the Voronoi diagram using matplotlib
- Optionally, output the mesh to an STL file

To execute the example driver, follow the instructions below.

1. Navigate to the examples directory.

   ```bash
   cd VoronoiWall/examples
   ```

2. Run the driver program, providing an *input file* and optionally an *output file*.

   ```bash
   python example_driver.py input_file [output_file]
   ```

   Some sample input files are provided:

   ```bash
   python example_driver.py input-random-50.txt random-50.stl
   python example_driver.py input-random-10.txt
   ```