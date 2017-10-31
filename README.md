# immunodbscan
A CUDA-based program for doing fast clone collapsing of B and T cell immune repertoires using DBSCAN. 

# Installation
This package is created using CMAKE. 
Within the main folder of this package, just type

 cmake .

and then 

 make

# Starting example run

Within the main folder, start the program with

./psdbscan

This will process the file "data/ExampleData.csv" using the configuration file "data/psdbscan.conf". 
The output file and a log will be found in folder "data".

If you need some help for program options, just type

./psdbscan --help

# Input data
This program uses a typical output of IMGT or other systems containing information about V, D and J genes as well as CDR3 sequences 
to do a clone collapsing using the DBSCAN clustering algorithm.

# Configuration file
The configuration file is written in the INFO file format described for the BOOST package "property_tree"-
It's nearly self-explaining.
