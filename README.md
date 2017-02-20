# immunodbscan
A CUDA-based program for doing fast clone collapsing of B and T cell immune repertoires using DBSCAN

# Principle
This program uses a typical output of IMGT or other systems containing information about V, D and J genes as well as CDR3 sequences 
to do a clone collapsing using the DBSCAN clustering algorithm.

# Extension for increased speed
The FM-DBSCAN variant of the DBSCAN algorithm looks promising to increase the speed to human-toleratable and reasonable levels.