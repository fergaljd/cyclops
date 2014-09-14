PREPARING ALONE CYCLOPS EXECUTABLES
CycloPs can be compiled as a stand-alone executable using the setup.py scripts in py2exe (Windows),
py2app (Mac OS X), and cx_freeze (Linux).


CUSTOM AMINO-ACID STRUCTURES
CycloPs is capable of incorporating custom amino structures, which must be defined in a file called 
‘zinc_output’ in the current working directory.

To add custom amino-acid structures which will be available to peptide generation functions, add the peptides
to the zinc_output file in the format <SMILES>\t<NAME> where \t is the TAB character, <SMILES> is the
peptide SMILES string with the N-terminus at the start, and the C-terminus at the end, in
the form C(=O)O, and <NAME> is the amino-acid displayed name. 
The aa_converter.py script supplied with CycloPs can re-arrange an amino-acid SMILES string to the format usable by CycloPs.
These peptides will then be available to add to the library in the 'ZINC amino acids' 
option in the 'Modify Amino-acid library' menu.


GENERATING PEPTIDE LIBRARIES FROM TEXT FILES
CycloPs can generate peptide libraries from text files of peptide sequences and desired constraints.
See the included 'example_peptide_file.txt' file for the an example and description of syntax.
