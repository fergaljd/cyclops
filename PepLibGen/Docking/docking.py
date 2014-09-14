#!/usr/bin/env python
'''Module containing functions to enable calling of Vina docking routines from 
python 
Must have autodock vina and prepare_ligands4.py script available in your PATH

Recommended to import as: 
from PepLibGen.Docking import docking as d'''

import os
import glob
import subprocess

max_unix_dirs = 30000 #Do not create more than this no. of dirs in one folder

def run_pdbqt_script(pdbfile, output_dir):
	'''runs the basic prepare_ligand4.py script, capturing output and error'''
	filebasename = os.path.split(pdbfile)[1].split('.')[0]
	output = os.path.join(output_dir,filebasename +'.pdbqt')
	
	try:
		proc = subprocess.Popen(['prepare_ligand4.py','-l', pdbfile,\
		'-o',output], stderr=subprocess.PIPE,stdout=subprocess.PIPE)
		
		prog_output = proc.communicate()
		proc.wait()
		if proc.returncode != 0: 
			print 'Stdout: ', prog_output[0]
			print 'Stderr: ', prog_output[1]
	except OSError:
		print 'Could not run program - check current working directory and path'
		raise
		
	return output

def convert_ligand_library(ligand_folder):
	output_dir =os.path.join(os.path.split(os.path.abspath(ligand_folder))[0]\
																	,'PDBQTs')
	if os.path.exists(output_dir) != True: os.makedirs(output_dir)	
	files = glob.glob(os.path.join(ligand_folder,'*.pdb'))
	if os.path.exists(ligand_folder) == False:
		print 'Ligand folder %s not found' % ligand_folder
	elif len(files) == 0: 
		print 'No input files found in %s' % ligand_folder
	counter = 0
	for f in files:
		pdbqtfile =run_pdbqt_script(f, output_dir)
		counter +=1
		if counter % 100 == 0: print '%s files converted' % counter
		
	print '%s files converted' % counter


def do_docking(receptor, ligand, config):
	'''Python wrapper to use vina to dock a ligand to a protein using a
	specified config file'''
	docking =subprocess.Popen(runstring,
				stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	docking.wait()
	prog_output = docking.communicate()
	if docking.returncode != 0: 
		return prog_output[1]
	else: return prog_output[0]
	
	


def generate_tasks(liganddir,receptordir, tasksdir='',makedirs=True):
	'''Prepare list of vina dockings for use by taskfarm utility on a cluster 
	Assumes config files are in same folder as receptors, with the same name, 
	but the extension changed to .conf'''
	if tasksdir =='': tasksdir = receptordir
	
	tasks = os.path.join(tasksdir, 'tasks')
	
	numwritten = 0
	#Use batches to break folders into max 30,000 docking subdirectories
	batch_num = 1 
	print 'Processing batch:  1'
	f = open(tasks, 'w')
	
	for ligand in glob.glob(os.path.join(liganddir,'*.pdbqt')):

		if numwritten % max_unix_dirs == 0 and numwritten != 0:
			batch_num += 1
			print 'Processing batch: ', batch_num
			
			
			
		outdir = os.path.join(receptordir, 'docking_tasks_batch_%s' % batch_num)
		ligandname = os.path.split(ligand.split('.')[0])[1]
	
		for receptor in glob.glob(os.path.join(receptordir,'*.pdbqt')):
			receptorname = os.path.split(receptor.split('.')[0])[1]
			config = receptor.split('.')[0] +'.conf'
			if os.path.exists(config) == False:
				print 'Warning: conf file for %s not found in receptor \
				directory' % receptorname
			dockingoutdirname = receptorname +'_' + ligandname
			dockingoutputdir = os.path.join(outdir, dockingoutdirname)
			dockingoutfile = os.path.join(dockingoutputdir, 'out.pdbqt')
			logfile = os.path.join(dockingoutputdir, 'log.txt')
		
			if os.path.exists(dockingoutputdir) != True:
				if makedirs == True:
					os.makedirs(dockingoutputdir)
			
			task = "cd '%s'; ulimit -s 300000; vina --cpu 1 --receptor '%s'\
			--ligand '%s' --config '%s' --out '%s' --log '%s'\n" \
			% (dockingoutputdir,receptor, ligand, config,dockingoutfile, 
			logfile)

			f.write(task)
			numwritten += 1
	f.close()
	return numwritten