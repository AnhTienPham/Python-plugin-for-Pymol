# Python-plugin-for-Pymol
Create Python plugin for Pymol to assist the protein surface masking project
This plugin encodes the function "identifySA", which allows the user to determine which surface-accessible residues within a chain of a PDB is being masked by a binding partner within a co-crystal structure.

To do this, the function uses two objects: one object where the chain in question is encoded alongside the binding partner, and another object where the chain is isolated without any binding partners. The original object is specified as the first argument passed to the function (eg. 6lr7), and the query chain is passed as the second argument (chain A). The starting and ending residues within the query chain are specified by the following two arguments (eg. 1,236). The function will perform this differential analysis for each residue, and return all residue positions that exhibit a difference greater than 0.5 angstroms.

Example:

identifySA 6lr7,chain A,1,236

[17, 90, 111, 113, 115, 122] #returned result
