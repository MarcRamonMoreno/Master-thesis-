# Master-thesis
Master Thesis Softmatter group ICMAB

Geometrical_carving.tcl was the custom tcl code used to carve out the square crystal, exposing the desired crystallographic planes (1 1 0) and (1 -1 0).
Carving_each_chain.tcl was the custom tcl code used to carve out the section of interest from the selection made with the Geometrical_carving.tcl code. This code selects the desired selection one cellulose chain at a time.

The file "glucose_orientation.tcl" creates a two-dimensional vector between carbon 2 and carbon 5 of each sugar molecule in the three spatial subdivisions of the nanocellulose crystal. It then calculates the orientation angle with respect to the vector (-1 1) at each time step.
The file "average_glucose_orientation.tcl" performs the same operation but averages the set of angles, obtaining a central tendency value for each spatial subdivision of the crystal.

The vector representation of the cellulose orientation angle was automated in the vector_representation_glucose_orientation.ipynb code, written in Python.

CNC_6_20.pdb and CNC_6_20.psf are the coordinates and the topology files of the generated structure of the cellulose nanocrystal model employed in this Master thesis.
