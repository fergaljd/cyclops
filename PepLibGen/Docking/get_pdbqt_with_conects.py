#!/usr/bin/env python
#A script to take a PDB ligand, convert it to an autodock4 format ligand,
#and add conect records from the original pdb

import os
import sys
import subprocess
import glob
import random

def run_pdbqt_script(pdbfile, home_dir, output_dir):
	'''runs the basic prepare_ligand4.py script, capturing output and error'''
	filebasename = os.path.split(pdbfile)[1].split('.')[0]
	output = output_dir + '/' + filebasename +'.pdbqt'
	os.chdir(home_dir)
	try:
		proc = subprocess.Popen(['python2.5','prepare_ligand4.py','-l', pdbfile,\
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

def pdb_parser(pdbfile):
	'''Processes a pdb file (or pdb-style file) and returns dictionaries of atom
	serials, types and x,y,z locations, and conect records'''
	atoms = {}
	conects = []
	f = open(pdbfile, 'r')
	for line in f:
		if 'ATOM' in line:
			l = line.strip().split()
			d ={'element':l[2][0],'xyz':(float(l[6]),float(l[7]),float(l[8]))}
			atoms[int(l[1])] = d
		elif 'HETATM' in line:
			l = line.strip().split()
			d ={'element':l[2][0],'xyz':(float(l[5]),float(l[6]),float(l[7]))}												
			atoms[int(l[1])] = d
		elif 'CONECT' in line:
			l = line.strip().split()
			d = {'serial':int(l[1]),'connections':[ int(ele) for ele in l[2:] ]}
			conects.append(d)
		else: pass
	f.close()
	
	return atoms, conects		
		
		
def add_conects(pdbfile, pdbqtfile, pdb_parser):
	'''Takes conect records from PDB file and adds them to PDBQT file. 
	Identifies atoms from there x,y,z coordinates in each file'''
	print 'PDBFILE ', pdbfile
	print 'PDBQTFILE', pdbqtfile
	
	pdb_atom_locs, pdb_conects = pdb_parser(pdbfile)
	pdbqt_atom_locs, pdbqt_conects = pdb_parser(pdbqtfile) #pdbqt_conects is empty
	
	#print 'pdb atoms', pdb_atom_locs
	#print 'pdbqt atoms', pdbqt_atom_locs

	#Translate pdb serials to corresponding pdbqt serials
	#Serial number dictionary, { pdb : pdbqt ... }
	serial_dict ={}
			
	for serial, value in pdb_atom_locs.iteritems():
		found = False
		for k,v in pdbqt_atom_locs.iteritems():
			if value['xyz'] == v['xyz']:
				serial_dict[serial] = k
				found = True
		if found == False:
			serial_dict[serial] = ''

	print 'serial_dict',serial_dict
			
	pdbqt_conects = []
	for dic in pdb_conects:
		d = {'serial':serial_dict[dic['serial']], 
			'connections':[ serial_dict[ele] for ele in dic['connections'] ] }
		pdbqt_conects.append(d)

	print 'pdbqt_conects',pdbqt_conects
	
	f = open(pdbqtfile, 'a')
	f.write('\n')
	for dic in pdbqt_conects:
		if dic['serial'] != '':
		#Build string to write:
			line = 'CONECT'
			for num in dic['connections']:
				if len(str(num)) > 4: 
					raise Exception, 'Too many atoms found - conect serial too long'
				else: 
					spacing = ' '* (5 - len(str(num)))
					line += spacing + str(num)
			line += '\n'
			f.write(line)		
	f.close()
		
		
		
def main():
	home_dir = os.getcwd()
	working_dir = sys.argv[1]
	output_dir =os.path.split(os.path.abspath(working_dir))[0] +os.path.sep+'PDBQTs'
	if os.path.exists(output_dir) != True: os.makedirs(output_dir)
	
	try:
		os.chdir(working_dir)
	except IndexError:
		print 'No input directory given'
		raise
		
	#files = random.sample((glob.glob('*.pdb')),1)
	files = ['YYYF-SCCTconst.pdb']
	for f in files:
		f = os.path.abspath(f)
		pdbqtfile =run_pdbqt_script(f, home_dir,output_dir)
		#add_conects(f, pdbqtfile, pdb_parser)
	
	
if __name__ == '__main__': main()
