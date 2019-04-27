# mass2mol
A Python implementation of the algorithm developed by Sebastian Bocker and Zsuzsanna Liptak (A Fast and Simple Algorithm for the Money Changing Problem, Algorithmica, 2007, 48, 413–432)

The mass2mol object generated the extended residue table for carbon, hydrogen, oxygen, phosphorus and sulfur with a blowup factor of 10^5 upon initialisation.

The combinations of CHNOPS atoms which have the same monoisotopic masses as a given mass can then be calculated using the find_formula method.
