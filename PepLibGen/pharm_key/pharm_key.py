#!/usr/bin/env python
'''
Pharmacophore key functions for large combinatorial libraries.

Uses rdkit pharmacophore functions, with Gobbi_Pharm2D features
'''

import os
import operator
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

class pharm_supplier(object):
	'''
	Reads in SMILES csv file, produces pharmacophores.
	
	'''
	def __init__(self, file, name_col, smi_col, delimiter):
		self.file_name = file
		self.name_col = name_col
		self.smi_col = smi_col
		self.delimiter=delimiter
		
	def supply_smiles(self):
		'''
		Returns a generator of smiles taken from the input file
		'''
		with open(self.file_name, 'r') as infile:
			for line in infile:
				yield line.split(self.delimiter)[self.smi_col-1].strip()
		
	def supply_mol(self):
		'''
		Returns a generator of rdkit molecule objects based on the input file
		'''
		with open(self.file_name, 'r') as infile:
			for line in infile:
				line = line.split(self.delimiter).strip()
				name = line[self.name_col-1]
				smiles = line[self.smiles_col-1]
				mol = Chem.MolFromSmiles(smiles)
				mol.SetProp('_Name', name)
				yield mol
		
	def supply(self):
		'''
		Returns a generator of pharmacaphore keys based on the input file
		'''
		for mol in self.supply_mol():
			yield Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
	
	


def read_delimited_file(name,column, delimiter=':'):
	'''
	Yields data from a delimited text file.
	
	If column is specified, returns only data from that column
	
	Otherwise splits line around delimiters, and returns each line

	'''
	with open(name, 'r') as smiles_file:
		for line in smiles_file:
			if column:
				yield line.split(delimiter)[column-1].strip()
			else:
				yield line.strip().split(delimiter)

def sub_struct_sim(target_pharm, query_pharm):
	'''
	Calculates similarity of two rdkit pharmacophore keys - 
	the tanimoto coefficient
	'''
	target_on_bits = list(target_pharm.GetOnBits())
	shared_on_bits = \
	[bit for bit in list(query_pharm.GetOnBits()) if bit in target_on_bits]
	query_on_bits = list(query_pharm.GetOnBits())
	return float(len(shared_on_bits)) / (float(len(target_on_bits)) \
	+ float(len(query_on_bits)) - float(len(shared_on_bits)))

def matching_compounds(target, library_list, cut_off, 
					   factory=Gobbi_Pharm2D.factory):
	'''
	Calculates which compounds in library_list are similar to a target compound.
	
	Peptides with a similarity above the cut_off value are returned. 
	Cut off is sub-structure similarity: (bits in common)/(bits in target)
	All compounds must be input as a SMILES string, or list of SMILES strings
	
	Returns a generator object
	'''
	
	target_mol = Chem.MolFromSmiles(target)
	target_pharm = Generate.Gen2DFingerprint(target_mol, factory)
	
	for smiles in library_list:
		mol = Chem.MolFromSmiles(smiles)
		query_pharm = Generate.Gen2DFingerprint(mol, factory)
		if sub_struct_sim(target_pharm, query_pharm) >= cut_off:
			yield smiles
		
def rank_compounds(target, library_list):
	'''
	Ranks compounds in terms of highest number of shared features with a target
	'''
	target_mol = Chem.MolFromSmiles(target)
	target_pharm = Generate.Gen2DFingerprint(target_mol, Gobbi_Pharm2D.factory)
	out_list = []

	for smiles in library_list:
		mol = Chem.MolFromSmiles(smiles)
		query_pharm = Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
		out_list.append((smiles, sub_struct_sim(target_pharm, query_pharm)))
		
	return out_list.sort(key=operator.itemgetter(1))
		