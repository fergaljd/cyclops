CycloPs is a GUI python programme and accompanying library for building constrained
peptide structures from peptide sequences, as described in the paper:
  CycloPs: generating virtual libraries of cyclized and constrained peptides including nonnatural amino acids
  FJ Duffy, M Verniere, M Devocelle, E Bernard, DC Shields, AJ Chubb
  Journal of chemical information and modeling 51 (4), 829-836

The cyclops repository contains all the code necessary to run the GUI CycloPs app.

PepLibGen is a python package that contains a variety of libraries with code 
to convert peptide sequences into 2D and 3D molecular file formats.
Most of the conversion code is in the PepLibGen.StructGen.StructGen.py library.

To run the CycloPs.py GUI app, make sure the PepLibGen package is available on your PYTHONPATH,
and you have the RDKit (http://www.rdkit.org/) python extensions, 
and the Python Imaging Library (http://www.pythonware.com/products/pil/) installed. 
You can then call CycloPs.py as a normal python script.

All the CycloPs code in this repository is released under the terms of the
GNU General Public license V2, described in the LICENSE.txt file.
