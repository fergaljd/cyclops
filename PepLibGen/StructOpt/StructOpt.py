#!/usr/bin/env python
'''Module to allow easy minimisation of PDB 3d structures of small molecules. 
Wraps around openbabel tools like obminimize and obconformer
Openbabel tools must be in your PATH for this module to be of any use

Author:Fergal Duffy. email:fergal.duffy@ucd.ie

Recommended to import with:
from PepLibGen.StructOpt import StructOpt as so'''

import os
import subprocess
import glob

import numpy as np

from PepLibGen.Analysis import analysis as a

class ExternalProgrammeError(Exception): pass


def get_saved_data(file, folder=''):
	'''Opens a python object picked to a file, and returns that object'''
	if folder != '': file = os.path.join(folder, file)
	filehandle = open(file, 'r')
	data = cPickle.load(filehandle)
	filehandle.close()
	
	return data

def write_object(obj, file, folder=''):
	'''Takes a python object and writes it to a file with cPickle '''
	if folder != '': file = os.path.join(folder, file)
	filehandle = open(file, 'w')
	cPickle.dump(obj, filehandle)
	filehandle.close()
	return True

def minimise(pdbfile, steps=250, outdir='', overwrite=True):
	'''Wraps around obminimize'''
	maindir = os.path.split(os.path.dirname(pdbfile))[0]
	
	if outdir == '':
		outdir = os.path.join(maindir, '3D_files_minimised')
		if os.path.split(pdbfile)[0] == '':
			raise Exception, 'Output folder needs to be \
			specified for conformergen pdbs'
	
	if os.path.exists(outdir) == False:
		os.makedirs(outdir)
	outpdbfile = os.path.join(outdir, os.path.split(pdbfile)[1])
	
	if overwrite == False:
		if os.path.exists(outpdbfile) == True:
			return '%s exists - skipping' % outpdbfile
			
	programme = 'obminimize'
	commandstring = '%s -n %s "%s" > "%s"' % (programme, str(steps), pdbfile,\
	outpdbfile)
	try:
		proc = subprocess.Popen(commandstring, stderr=subprocess.PIPE,\
		stdout=subprocess.PIPE, shell=True)
		
		returncode = proc.wait()
		output_msg = proc.communicate()
		if returncode != 1:
			print output_msg[1] #stderror
			print 'return code', returncode
			raise ExternalProgrammeError
		elif returncode == 1:
			return output_msg[0]
			
	except OSError:
		print 'Could not run program - check current working directory and path'
		raise

		
def conformer_minimise(pdbfile, numconformers=10, steps=100, outdir='', 
overwrite=True):
	'''Wraps around obconformer'''
	maindir = os.path.split(os.path.dirname(pdbfile))[0]
	
	if outdir == '':
		outdir = os.path.join(maindir, '3D_files_conformergen')
		if os.path.split(pdbfile)[0] == '':
			raise Exception, 'Output folder needs to be\
			specified for conformergen pdbs'
			
	if os.path.exists(outdir) == False:
		os.makedirs(outdir)
	outpdbfile = os.path.join(outdir, os.path.split(pdbfile)[1])
	
	if overwrite == False:
		if os.path.exists(outpdbfile) == True:
			return '%s exists - skipping' % outpdbfile
		
	programme = 'obconformer'
	commandstring = '%s %s %s "%s" > "%s"' % (programme, \
	str(numconformers), str(steps), pdbfile, outpdbfile)
	try:
		proc = subprocess.Popen(commandstring, stderr=subprocess.PIPE,\
		stdout=subprocess.PIPE, shell=True)
		
		returncode = proc.wait()
		output_msg = proc.communicate()
		if returncode != 0:
			print output_msg[1] #stderror
			print 'return code', returncode
			raise ExternalProgrammeError
		elif returncode == 0:
			return output_msg[0]
			
	except OSError:
		print 'Could not run program - check current working directory and path'
		raise

		
def min_library(folder, type='.pdb', steps=100,genconformers=False, confs=10, 
															overwrite=False):
	'''Runs either obminimize (default) or obconformer on a folder of structure 
	files (default is pdb)'''
	outfile = os.path.join(folder, 'minimise_output.txt')
	f = open(outfile, 'w')
	count = 0
	print 'Processing %s files...' % len(os.listdir(folder))
	for pdbfile in os.listdir(folder):
		pdbfile = os.path.join(folder, pdbfile)
		if genconformers == False:
			out = minimise(pdbfile, steps, overwrite)
		elif genconformers == True:
			out = conformer_minimise(pdbfile, confs, steps, overwrite)
		else:
			raise Error, 'genconformers must have a value of True or False'
		count += 1
		if count % 100 == 0: print '%s files processed' % count
		f.write(out)
	f.close()
		
def smart_min_library(folder, type='.pdb', steps=100, confs=10, 
do_max=False, overwrite=False):
	''' Like min_library function, but chooses whether to generate conformers 
	based on the bond energy for the library calculated by the analysis module. 
	Molecules with more than 5 times the median energy get random conformers 
	generated for them'''
	
	if do_max == False:
		numenergies = len(glob.glob(os.path.join(folder,'*.pdb')))
	if do_max != False:
		numenergies = min(len(glob.glob(os.path.join(folder,'*.pdb'))), do_max)
	
	print 'Getting %s energies...' % numenergies

		
		
	energies = a.calc_lib_energies(folder, do_max)
	median_energy = np.median([dict['bond_energy'] for dict in energies])
	
	outfile = os.path.join(folder, 'minimise_output.txt')
	f = open(outfile, 'w')
	count = 0
	
	#Decide if pdbs need to be run through obconformer
	input = []
	for energy_dict in energies:
		if energy_dict['bond_energy'] > median_energy * 2:
			input.append((energy_dict['name'], True))
		elif energy_dict['bond_energy'] <= median_energy * 2:
			input.append((energy_dict['name'], False))
		else:
			raise Exception, 'Impossible error'
	print 'Minimising %s files...' % len(input)		
	#Do the minimisation
	for pdbtuple in input:
		pdbfile = os.path.join(folder, pdbtuple[0])

		if pdbtuple[1] == False:
			out = minimise(pdbfile, steps*10, overwrite=overwrite)
		elif pdbtuple[1] == True:
			out = conformer_minimise(pdbfile,confs,steps,overwrite=overwrite)
		else:
			raise Error, 'genconformers must have a value of True or False'
		count += 1
		if count % 100 == 0: print '%s files processed' % count
		f.write(out)
	f.close()