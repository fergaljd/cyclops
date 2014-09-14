#!/usr/bin/env python
'''
A set of tools to help prepare docking results for analysis by Ligplot, and also
to analyse the generated Ligplot output

from PepLibGen.Docking import LigplotTools as lig
'''

import os
import glob
import cPickle
import csv
import operator
import shutil
import math
import pprint
from os import path

import numpy as np

from PepLibGen.tools import tools

class VinaLogParseError(Exception): pass
class ParseResultsFoldersError(Exception): pass
class PDBFormatException(Exception):pass
class LigplotResultsError(Exception):pass
class NormalisationError(Exception):pass



def basic_prepare_ligplot(lig, rec, outdir, lig_resi_num = '500',\
 lig_chain_id ='Z'):
	'''
	Cut back version of prepare ligplot. Will combine a receptor and ligand
	given their paths, and write the result to outdir
	'''
	ligfile = open(lig, 'r')
	recfile = open(rec, 'r')
	ligfile_lines = ligfile.readlines()
	recfile_lines = recfile.readlines()
	ligfile.close()
	recfile.close()
	lig_name = path.basename(lig).split('.')[0]
	rec_name = path.basename(rec).split('.')[0]
	
	#Filter unwanted models and torsion tree data from ligand.
	#Only takes the first model
	add_line = False
	new_ligfile_lines = []
	
	for num,line in enumerate(ligfile_lines):
		#Either take first model, or if only one model, MODEL will not appear
		if 'MODEL 1' in line or ('MODEL' not in line and num==1):
			add_line = True
		elif 'ENDMDL' in line:
			break
		elif add_line == True: 
			if 'HETATM' in line or 'ATOM' in line:
				#Modify line to set chain ID and residue number
				#CHANGE ATOM recors to HETATM
				spacing = 5 -len(str(lig_chain_id)) - len(str(lig_resi_num))
				if len(str(lig_chain_id)) > 1:
					raise PDBFormatException,\
					'Chain identifier must be a single letter: not %s' \
					% lig_chain_id
				if len(str(lig_resi_num)) > 4:
					raise PDBFormatException,\
					'Resi number must be a max 4-digit number: not %s' \
					% lig_resi_num
				newline = \
				'HETATM'+line[6:17]+'INT ' + str(lig_chain_id) + ' '*spacing\
				+str(lig_resi_num)+line[26:]
						
				new_ligfile_lines.append(newline)
	
	#Now, write a Ligplot directory to the results
	ligplot_dir = path.join(outdir,'ligplot')
	if not path.exists(ligplot_dir):
		os.mkdir(ligplot_dir)
	
	#Add ligand to receptor lines...
	ligplotlines = list(recfile_lines) #take a copy
	ligplotlines.extend(new_ligfile_lines)

	
	#write file for ligplot to read...
	fname = rec_name+'_'+lig_name+'_ligplot.pdbqt'
	ligplotfile = open(path.join(ligplot_dir,fname), 'w')
	ligplotfile.writelines(ligplotlines)
	ligplotfile.close()
	
def prepare_ligplot(ligdir, recdir, lig_resi_num = '500', lig_chain_id ='Z'):
	'''
	Generates a structure for Ligplot by combining the receptor and ligand in
	a new file 'ligplot_recname_ligname.pdbqt' in a ligplot folder in the
	results directory, which ligplot can understand
	'''
	if path.isdir(ligdir) == True:
		try:
			receptor_name, ligand_name = path.basename(ligdir).split('_',1)
		except ValueError:
			print 'Error on the following input: '
			print 'ligdir', ligdir
			print 'recdir', recdir
			raise
	else:
		#Ligdir is not a directory - assume it points to the ligand file.
		try:
			receptor_name, ligand_name \
			= path.basename(path.dirname(ligdir)).split('_',1)
		except ValueError:
			print 'Error on the following input: '
			print 'ligdir', ligdir
			print 'recdir', recdir
			raise
	
	#Get contents of receptor and ligand PBDQT files
	try:
		#If ligdir is the path to the ligand, this will work, otherwise...
		ligfile = open(ligdir, 'r')
	except IOError:
		#IOError: ligdir is a directory
		try:
			ligfile = open(glob.glob(path.join(ligdir, '*.pdbqt')))
		except TypeError:
			#TypeError - more than one *.pdbqt found by glob - list returned
			ligfile = open(path.join(ligdir, 'out.pdbqt'))
			
	try:
		recfile = open(recdir, 'r')
	except IOError:
		try:
			#find receptor.pdbqt or receptor.pdb
			recfile = open(glob.glob(path.join(recdir,receptor_name+'.pdb*')))
		except TypeError: 
			#If receptor.pdb and receptor.pdbqt exist...
			recfile = open(path.join(recdir,receptor_name+'.pdbqt'))
			
	ligfile_lines = ligfile.readlines()
	ligfile.close()
	recfile_lines = recfile.readlines()
	recfile.close()
	
	#Filter unwanted models and torsion tree data from ligand.
	#Only takes the first model
	add_line = False
	new_ligfile_lines = []
	
	for num,line in enumerate(ligfile_lines):
		if 'MODEL 1' in line or ('MODEL' not in line and num==1):
			add_line = True
		elif 'ENDMDL' in line:
			break
		elif add_line == True and 'HETATM' in line or 'ATOM' in line:
			#Modify line to set chain ID and residue number
			#CHANGE ATOM recors to HETATM
			spacing = 5 -len(str(lig_chain_id)) - len(str(lig_resi_num))
			if len(str(lig_chain_id)) > 1:
				raise PDBFormatException,\
				'Chain identifier must be a single letter: not %s' \
				% lig_chain_id
			if len(str(lig_resi_num)) > 4:
				raise PDBFormatException,\
				'Resi number must be a max 4-digit number: not %s' \
				% lig_resi_num
			newline = \
			'HETATM'+line[6:17]+'INT ' + str(lig_chain_id) + ' '*spacing\
			+str(lig_resi_num)+line[26:]
					
			new_ligfile_lines.append(newline)
	
	#Now, write a Ligplot directory to the results
	
	try:
		if path.isdir(ligdir) == True:
			os.mkdir(path.join(ligdir,'ligplot'))
			ligplot_dir = path.join(ligdir,'ligplot')
		else:
			#Ligdir is not a directory - assume it points to the ligand file.
			os.mkdir(path.join(path.dirname(ligdir), 'ligplot'))
			ligplot_dir = path.join(path.dirname(ligdir), 'ligplot')
	except OSError:
		if path.exists(path.join(ligdir,'ligplot')) or \
		path.exists(path.join(path.dirname(ligdir), 'ligplot')):
			if path.exists(path.join(ligdir,'ligplot')):
				ligplot_dir = path.join(ligdir,'ligplot')
			else:
				ligplot_dir = path.join(path.dirname(ligdir), 'ligplot')
		else:
			raise
			
	#Add ligand to receptor lines...
	ligplotlines = list(recfile_lines) #take a copy
	ligplotlines.extend(new_ligfile_lines)
	
	#write file for ligplot to read...
	fname = receptor_name+'_'+ligand_name+'_ligplot.pdbqt'
	if not path.exists(path.join(ligplot_dir,fname)):
		ligplotfile = open(path.join(ligplot_dir,fname), 'w')
		ligplotfile.writelines(ligplotlines)
		ligplotfile.close()
	
def prepare_all_ligplot(top_dir, receptor_dir):
	'''
	Prepares a docking directory for ligploting by running through all results 
	directories, and creating the input structure for ligplot to analyse
	'''
	for dir, subdirs, files in os.walk(top_dir):
		if glob.glob(path.join(dir, '*.pdbqt'))\
		and 'ligplot' not in path.basename(dir): 
		#Will choke if tries to recreate directory
			prepare_ligplot(dir, receptor_dir)

def get_ligplot_bonds(res_dir, chain_id='INT Z  500 '):
	'''
	Looks in <res_dir> for ligplot.hhb (hydrogen bonds )and ligplot.nnb 
	(hydrophibic contacts) and returns the receptor atoms involved in hydrogen bonds
	and the residues involved in hydrophobic contacts.

	chain_id corresponds to the resi name, chain, and resi number of the ligand in
	the pdb sent to ligplot. Spacing is important:
	3 spaces for resi name, blank, 1 letter chain id, blank, 4 spaces for resi 
	number
	'''
	#Spaces are important when delaing with these files!
	chain_id = chain_id.strip()
	try:
		hbondf = open(path.join(res_dir, 'ligplot.hhb'))
		nbondf = open(path.join(res_dir, 'ligplot.nnb'))
	except IOError:
		print 'Could not find one of %s or %s' % \
		(path.join(res_dir, 'ligplot.hhb'),path.join(res_dir, 'ligplot.nnb'))
			

	hbonds = []
	nbonds = []

	for line_num,line in enumerate(hbondf):
		if line_num > 2 and chain_id in line:
		#Must have ligand only once in line - as donor, or acceptor, not both
		#Take detail of receptor - it can be donor or acceptor
			if line[:11].strip() == chain_id \
			and line[21:32].strip() != chain_id:
				acceptor_id = line[21:32].strip()
				acceptor_atom = line[34:37].strip()
				acceptor_distance = float(line[41:45].strip())
				hbonds.append(('acceptor', acceptor_id, \
				acceptor_atom, acceptor_distance))
			elif line[21:32].strip() == chain_id \
			and line[:11].strip() != chain_id:
				donor_id = line[:11].strip()
				donor_atom = line[13:16].strip()
				donor_distance = float(line[41:45].strip())
				hbonds.append(('donor',donor_id,donor_atom, donor_distance))
	hbondf.close()
	
	for line_num,line in enumerate(nbondf):
		if line_num > 2 and chain_id in line:
		#As above, find ligand once in line
			if line[:11].strip() == chain_id\
			and line[21:32].strip() != chain_id:
				atom_id = line[21:32].strip()
				atom = line[34:37].strip()
				atom_distance = float(line[41:45].strip())
				nbonds.append(('atom', atom_id, atom, atom_distance))
			elif line[21:32].strip() == chain_id\
			and line[:11].strip() != chain_id:
				atom_id = line[:11].strip()
				atom = line[13:16].strip()
				atom_distance = float(line[41:45].strip())
				nbonds.append(('atom',atom_id,atom, atom_distance))
	nbondf.close()
				
	return hbonds, nbonds
		
def compare_ligplot_bonds(main_ligdir, comp_ligdir):
	'''
	Compare ligands binding to a receptor, this function returns a score
	x = a/b where b is the total number of bonds recorded in the ligplot output
	files in main_ligdir and a is the number of bonds shared between the ligplot
	output in main_ligdir and comp_ligdir
	
	Should only be used to compare structures with 2 different ligands binding 
	to the same receptor, at the same location.
	'''
	
	main_hbonds, main_nbonds = get_ligplot_bonds(main_ligdir)
	comp_hbonds, comp_nbonds = get_ligplot_bonds(comp_ligdir)
		
	shared_hbonds = []
	shared_nbonds = []
	
	#Get rid of atom distances when making comparisons
	#Distance is last entry
	main_hbonds = [bond[:-1] for bond in main_hbonds]
	main_nbonds = [bond[:-1] for bond in main_nbonds]
	comp_hbonds = [bond[:-1] for bond in comp_hbonds]
	comp_nbonds = [bond[:-1] for bond in comp_nbonds]

	for hbond in main_hbonds:
		if hbond in comp_hbonds:
			#handle duplicate bonds in comp_bonds- only add twice if 
			#present twice in main_bonds - so pop out of comp_bonds list
			#during iteration
			#It is possible for two bonds to be made to the same 
			#place in a receptor
			loc = comp_hbonds.index(hbond)
			shared_hbonds.append(comp_hbonds.pop(loc))
			
	for nbond in main_nbonds:
		if nbond in comp_nbonds:
			loc = comp_nbonds.index(nbond)
			shared_nbonds.append(comp_nbonds.pop(loc))
			
	result = float(len(shared_hbonds)+len(shared_nbonds))\
	/float(len(main_hbonds)+len(main_nbonds))
	
	if result < 0:
		raise LigplotResultsError, '%s and %s give negative ligplot score!'\
		% (main_ligdir, comp_ligdir)
	elif result > 1:
		print 'Result', result
		print 'Num Shared hbonds, nbonds', shared_hbonds, shared_nbonds
		print 'Main hbonds, nbonds', main_hbonds, main_nbonds
		raise LigplotResultsError, '%s and %s give a ligplot score over 1'\
		% (main_ligdir, comp_ligdir)
	else:
		return result
		
def get_ligplot_scores(true_lig_dir, dock_results_top_dir,
obfile='ligplot_scores.obj', force_recalc = False):
	'''
	Looks for ligplot output of receptors and true positive ligands in 
	true_lig_dir. Results should be found in true_lig_dir/<receptorName>/ligplot
	Will take each receptorName/ligplot output data, and analyse the bonds in 
	ligplot.nnb and ligplot.hhb. It will take these results, and walk through 
	the dock_results_top_dir directory, looking for all docking outputs 
	corresponding to the particular receptor. It will the calculate the ligplot 
	interaction score for each docking.
	
	Returns a list of lists of the form [receptorname, ligandname, score]
	'''
	true_lig_dirs = {} #Will be of form 'name':'ligplot dir path'
	ligplot_scores = []
	
	if not force_recalc:
		try:
			print 'Attempting to read from %s...' % \
			path.join(dock_results_top_dir,obfile)
			ligplot_scores = tools.get_saved_data(obfile, dock_results_top_dir)
		except IOError:
			print 'Saved results not found, extracting from directories...'
			force_recalc=True
	
	if force_recalc:
		for dir, subdirs, files in os.walk(true_lig_dir):
			if path.basename(dir) == 'ligplot':
				#Name of parent folder is receptor name
				receptor_name = path.basename(path.dirname(dir))
				true_lig_dirs[receptor_name] = dir
				
		print 'Found %s receptor ligplot outputs' % len(true_lig_dirs)
			
		for dir, subdirs, files in os.walk(dock_results_top_dir):
			if path.basename(dir) == 'ligplot':
				#Get names of receptor and ligand
				#Get containing directory, and splitfields
				#Folder heirarchy should be 
				#of form /rec_lig/ligplot/rec_lig_ligplot.pdbqt
				updir = path.dirname(dir)
				receptor_name, ligand_name = path.basename(updir).split('_', 1)
				comp_ligdir = dir
				try:
					main_ligdir = true_lig_dirs[receptor_name]
				except KeyError:
					print 'Could not find receptor %s in true_lig_dirs %s' \
					% receptor_name, str(true_lig_dirs)
					raise
				#try get scores
				try:
					score = compare_ligplot_bonds(main_ligdir, comp_ligdir)
					ligplot_scores.append({'receptor':receptor_name, \
					'ligand':ligand_name, 'ligplot_score':score})
				except IOError:
					print 'Could not get scores for ',ligand_name, receptor_name
					print 'main_ligdir', main_ligdir
					print 'comp_ligdir', comp_ligdir
					raise
				
				if len(ligplot_scores) % 100 == 0 and len(ligplot_scores) != 0:
							print 'Processed %s scores' % len(ligplot_scores)
		print 'Processed %s scores' % len(ligplot_scores)	
			
	print 'Writing results to %s' % path.join(dock_results_top_dir, obfile)
	tools.write_object(ligplot_scores, obfile, dock_results_top_dir)
		
	return ligplot_scores

	
