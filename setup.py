import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="VoronoiWall", # Replace with your own username
    version="1.0.0",
    author="CSC@NCSU",
    author_email="drsheehy@ncsu.edu",
    description="Pipeline for processing the Voronoi diagram of a point set, producing STL formatted output for the purpose of 3D printing.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/donsheehy/VoronoiWall",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)