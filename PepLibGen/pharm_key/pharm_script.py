#!/usr/bin/env python
'''
Script to read in a file of Smiles strings, and a tergat smiles, and 
do some pharmacophore related stuff
'''

from PepLibGen.pharm_key import pharm_key
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
import timeit

#1b11 ligand
target = '[C@@H](CO)(Cc1ccccc1)NC(=O)[C@H](C(C)C)NC(=O)[C@H](C)NC(=O)OCc1ccccc1'
cut_off = 0.9

def main():
	smiles_lib=pharm_key.read_delimited_file('trial_pharmacophore_library.txt',
	2)
	#target_mol = Chem.MolFromSmiles(target)
	#target_pharm = Generate.Gen2DFingerprint(target_mol, Gobbi_Pharm2D.factory)
	
	print 'Matching...'
	count = 0
	for match in pharm_key.matching_compounds(target, smiles_lib, 0.9):
		count += 1
		print 1;
	else:
		print count
		
if __name__ == '__main__':
	main()