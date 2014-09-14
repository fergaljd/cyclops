#!/usr/bin/env python
'''Module to allow analysis of energies, hydrophobicity and associated plots, 
and ramachandran plotting. Uses openbabel obenergy programme to calculate 
energies, and BioPython data to estimate hydrophobicitys.
Ramachandran plot code has been written from first principles.

Recommended to import with:
from PepLibGen.Analysis import analysis as a
'''
import sys
import os
import glob
import operator
import random
import subprocess
import cPickle

import matplotlib.pyplot as plt
import numpy as np

from Bio.SeqUtils.ProtParamData import kd #BioPython dict of hydrophobicity

class ExternalProgrammeError(Exception): pass



#-############################# ENERGY #########################################
output_file = 'energy_output.csv'

findpatterns = {'name':'.pdb','bond_stretch_energy':\
	'TOTAL BOND STRETCHING ENERGY', \
	'bond_angle_energy':'TOTAL ANGLE BENDING ENERGY',\
	'bond_tors_energy':'TOTAL TORSIONAL ENERGY',\
	'bond_vdw_energy':'TOTAL VAN DER WAALS ENERGY','bond_energy':'TOTAL ENERGY'}

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

	
def get_energies(pdbfile):
	'''Run obenergy on a pdb file. Outputs a list of dictionaries - 
	each dictionary having filanames and energy types and values
	
	Will look for a file called energies.txt containing a pickled version
	of the energies saved previously'''
	
	energies = [] #Will be a list of dictionaries of molecular energies

	try:
		proc = subprocess.Popen(['obenergy',pdbfile], stderr=subprocess.PIPE,\
								stdout=subprocess.PIPE)
		returncode = proc.wait()
		output_msg = proc.communicate()
		if returncode != 1:
			print output_msg[1] #stderror
			print 'return code', returncode
			raise ExternalProgrammeError
		elif returncode == 1:
			raw_energy_data = output_msg[0]
			
	except OSError:
		print 'Could not run program - check current working directory and path'
		raise
	
	energy_dict={'name':os.path.split(pdbfile)[1]}
	#Loop to extract energy numbers from stdout of obenergy
	for line in raw_energy_data.split('\n'):
		line.strip()
		for key,value in findpatterns.items():
			if value == '.pdb':
				pass
			elif value in line:
			#Go through the line containing the energy value, and find the bit 
			#that's a float
				line = line.split(' ')
				for pos_energy in line:
					try:
						energy = float(pos_energy)
						energy_dict[key] = energy
					except ValueError: pass
					
	#Test if all energies were found:
	if len(energy_dict) != len(findpatterns):
		notfound = []
		for key in findpatterns.keys():
			if key not in energy_dict:
				notfound.append(key)
		raise Exception, 'Energy values not found: %s' % ', '.join(notfound)
						

	
	return energy_dict

def do_energy_plotting(energies):
	'''Produce plots with matplotlib of peptide energies supplied as a 
	dictionary'''
	
	#Sorts list of energies by total energy before plotting
	energies.sort(key = operator.itemgetter('bond_energy'))
	
	#Write energy data into lists for plotting
	fig = plt.figure()
	log_plot = fig.add_subplot(2,1,1)
	lin_plot = fig.add_subplot(2,1,2)
	
	for key, value in findpatterns.iteritems():
		if key == 'name': pass
		else:
			if key == 'bond_energy': plotmax = max([d[key] for d in energies])
			plot_data = [d[key] for d in energies]
			logplot_ref, = log_plot.semilogy(plot_data)
			plot_ref, = lin_plot.plot(plot_data)
			logplot_ref.set_label(key)
			plot_ref.set_label(key)
	
	for a in log_plot, lin_plot:
		a.grid(True)
		a.axis([0,len(energies), 0, plotmax])
		a.set_xlabel('Sorted bond energy data')
		a.set_ylabel('Bond energy (KJ/mol)')
		a.legend(loc = 'upper left', prop={'size':'small'})
		
	#plt.savefig('energies.png')

	plt.show()
		
		
def calc_lib_energies(pdbfolder, do_max=False):
	'''Given a library folder of peptide PDBs, will return a sorted list of 
	energies do_max can be given set to an integer if only a random sample of 
	energies is required, will sample <do_max> pdbs'''
	#Attempt to read pickled energy values - if all energies are wanted
	try:
		energies = get_saved_data('energies.obj', pdbfolder)
		#Test that the expected number of energy values exist:
		if len(energies) == len(glob.glob(os.path.join(pdbfolder,'*.pdb'))):
			print '%s energies read from %s' %(len(energies),'energies.obj')
			return energies
		else:
			print 'Number of energies read does not match number of pdbs'
		
	except IOError:
		print 'No existing energies found - will calculate them' 
		
	if do_max == False:
		input = glob.glob(os.path.join(pdbfolder, '*.pdb'))
	else:
		try:
			input = random.sample(glob.glob(os.path.join(pdbfolder, '*.pdb')), \
																		do_max)
		except ValueError:
			print "'do_max' must be either False or an integer"
			raise
			
	energies = []
	print 'Calculating %s energies...' % len(input)
	count = 0
	for structure in input:
		
		energies.append(get_energies(structure))
		count +=1
		if count % 1000 == 0: print '%s files processed' % count
		

	energies.sort(key = operator.itemgetter('bond_energy'))

	#Pickle to file - only if all energies are read.
	if do_max == False:
		pickler = write_object(energies, 'energies.obj', pdbfolder)
		if pickler == True:
			print 'Serialised energy values to file %s' % 'energies.obj'
	
	return energies
	

	
#-############################# HYDROPHOBICITY #################################
	

def hydrophobes(pepstring):
	'''Analyze peptide string for hydrophibicity. 
	Return result as dictionary'''
	
	hydrophob_values = []	
	for letter in pepstring:
		hydrophobicity = kd[letter]
		hydrophob_values.append(hydrophobicity)
	
	results = {}
	results['name'] = pepstring
	results['hydrophobicitys'] = hydrophob_values
	results['average_hydrophobicity'] =\
									sum(hydrophob_values)/len(hydrophob_values)
	results['num_hydrophobic_residues'] = \
						len([resi for resi in hydrophob_values if resi < 0])
	return results

def do_hydrophobicity_plotting(library_hydrophobes):
	'''Produce plots of average hydrophobicity, hydrophobic residues, and top 
	hydrophobes, given a list of dictionaries of hydrophobicity properties
	'''
	average_hydrophobicity = [] 
	number_hydrophobic = []
	for peptide in library_hydrophobes:
		average_hydrophobicity.append(peptide['average_hydrophobicity'])
		number_hydrophobic.append(peptide['num_hydrophobic_residues'])
	
	average_hydrophobicity.sort()	
	number_hydrophobic.sort()
	
	#Prepare plots of total hydrophobicity, and number of hydrophobicresidues
	x1 = \
		np.linspace(1, len(average_hydrophobicity), len(average_hydrophobicity))
	y1 = average_hydrophobicity
	y2 = number_hydrophobic
	
	plt.figure(1)
	plt.subplot(3,1,1)
	line1, = plt.plot(x1,y1, 'b-')
	line1.set_label('Average hydrophobicity')
	plt.legend(loc='lower right')
	plt.ylabel('Average hydrophobicity')
	plt.title('Hydrophibicity of Peptide Library', fontsize = 15)
	plt.grid(True)
	
	plt.subplot(3,1,2)
	line2, = plt.plot(x1,y2, 'r-')
	line2.set_label('Number of hydrophobic residues')
	plt.legend(loc='lower right')
	plt.ylabel('Number of hydrophobic residues')
	plt.grid(True)
	
	#Make a barchart showing 10 most and least hydrophobic peptides.
	bardata = sorted(library_hydrophobes, 
							key=operator.itemgetter('average_hydrophobicity'))
	top_hydrophobes = [d['average_hydrophobicity'] for d in bardata[:10]]
	bottom_hydrophobes = [d ['average_hydrophobicity'] for d in bardata[-10:] ]
	plotdata = bardata[:10] + bardata[-10:]
	
	plt.subplot(3,1,3)
	plt.bar(range(len(top_hydrophobes)), top_hydrophobes, color = 'r')
	plt.bar(range(len(top_hydrophobes), len(plotdata)), \
												bottom_hydrophobes, color='b')
	plt.xticks(range(len(plotdata)),[ d['name'] for d in plotdata],
																rotation = 45)
	plt.grid(True)
	plt.ylabel('10 most and least hydrophobic peptides')
	plt.savefig('hydrophobicity.png')
	
	
	plt.show()	
	
def library_hydrophobicity(pdbfolder):
	'''Take folder of pdb files as input, and process input list of files
	Also produces various plots of hydrophobic properties'''
	

	library_hydrophobes = []
	for peptide in glob.glob(os.path.join(pdbfolder,'*.pdb')):
		pepstring = os.path.split(peptide)[1].split('-')[0]
		residue_data = hydrophobes(pepstring)
		library_hydrophobes.append(residue_data)
	
	library_hydrophobes.sort(key =operator.itemgetter('average_hydrophobicity'))
	return library_hydrophobes
	
#-############################# RAMACHANDRAN PLOT ##############################

def calc_dihedral_angle(locs):
	'''Calculates the dihedral angle between two planes assuming that points 
	A, B, C are in plane 1, and B, C, D are in plane 2. Each point consists
	of tuples of (x,y,z) locations 
	For the peptide ramachandran plot, phi A,B,C,D values are the locations of 
	C'1, N, Ca, C'2 and for the psi plot they are N1, Ca, C', N2
	Returns answer in radians'''
	
	a,b,c,d = locs
	
	x1,y1,z1 = a[0],a[1],a[2]
	x2,y2,z2 = b[0],b[1],b[2]
	x3,y3,z3 = c[0],c[1],c[2]
	x4,y4,z4 = d[0],d[1],d[2]
	
	#Create 3 vectors = b1, b2, b3
	
	b1 = np.array( [x2-x1, y2-y1, z2-z1] )
	b2 = np.array( [x3-x2, y3-y2, z3-z2] )
	b3 = np.array( [x4-x3, y4-y3, z4-z3] )

	# angle = atan2(|b2|b1 . [b2 x b3], [b1 x b2].[b2 x b3])
	# from http://en.wikipedia.org/wiki/Torsion_angle
	
	#Calculate angle (in radians)
	angle = np.arctan2(np.dot(np.linalg.norm(b2)*b1,\
					np.cross(b2,b3)),np.dot(np.cross(b1,b2),np.cross(b2,b3)))
		
	angle = np.degrees(angle)
	return angle	

	
def follow_paths(conects,name,paths):
	'''Recursive function to follow paths, returning a list of complete paths
	Will handle loops in path'''
	
	newpaths = [] #New path list will be built from old paths list
	
	#Initialise, if first run.
	if paths == []:
		#Every atom an be reached from every other, 
		#doesn't matter where we startpos
		startpos = conects[0]['serial'] 						 
		for num in conects[0]['connections']:
			paths.append([startpos, num])
	try:
		for path in paths:
			if path[-1] != 'f' and path[-1] != 'l':
				#Follow the next connections in the path. Assumes atom 
				#serial numbers start at 1 and go up
				conect_rec = []
				for d in conects:
					if d['serial'] == path[-1]:
						conect_rec = list(d['connections']) #Create a copy
				#print 'Paths',paths
				conect_rec.remove(path[-2]) #prevent return from way we came
				
				if len(conect_rec) > 0: #If new conections exist
					
					for num in conect_rec:
						newpath = list(path) #Create copy of path to work with
						if num in newpath:
							#We've found a loop
							newpath.append(num)
							newpath.append('l')
						else:
							newpath.append(num)
						newpaths.append(newpath)
				elif len(conect_rec) == 0:
					#path is followed completely
					newpath = list(path) #Dont modify old path - build new one
					newpath.append('f')
					newpaths.append(newpath)
				else:
					raise Exception, 'Something very strange has happened'				
		
			else:
				#Keep old path unchanged
				newpath = list(path)
				newpaths.append(newpath)
	except ValueError:
		print 'ValueError, conects = %s path = %s ' % (conects, paths)
		print 'Filename', name
		raise
			
	#print newpaths
	#Check if all paths are completed, return paths, or run again.
	num_paths = len(newpaths)
	completed = [path for path in newpaths if path[-1] == 'f' or path[-1]== 'l']
	num_completed = len(completed)
	if num_paths == num_completed:
		paths = []
		return newpaths
	elif num_paths > num_completed:
		return follow_paths(conects,name, newpaths)
	else:
		raise Exception, 'Weird Error!'	
	
def pdb_ramachandran_parse(pdbfile):
	'''Takes a PDB file as input, attmepts to identify all Calphas, Ns and 
	C's required for ramachandran plot. Written to deal with peptides generated 
	by openbabel which have no, or broken residue information
	
	Returns (x,y,z) locations necessary for phi and psi angles.'''
	f = open(pdbfile)
	atoms = {}
	conects = []
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
	
	#Try and walk through all possible connection paths. The longest loop will 
	#be used
	#to pick out Ca, C' and Ns.
	paths = follow_paths(conects, pdbfile, [])
	#Leave out 'l' or 'f' with [:-1]
	#longest_loop = max([path for path in paths if path[-1] =='l'], key=len)[:-1]	
	#Improve longest loop code by only mesuring LOOP lenghth, 
	#not length of a whole path that contains a loop
	loops = []
	for path in paths:
		if path[-1] == 'l':
			#Atom that closes loop
			endloop = path[-2]
			#Atom that starts loop  = atom that closes loop
			startloop = path.index(endloop)
			loops.append(path[startloop:-1])
	longest_loop = max(loops, key=len)


	#Find out what's in the longest loop:
	atomloop = ''
	for num in longest_loop: 
		try: 
			atomloop+=atoms[num]['element']
		except KeyError:
			print 'Key error finding key %s in dictionary \n %s' % (num, atoms)
			raise
			
	#Get serial numbers of all atoms needed to plot phi and psi
	pattern = 'CNCCN' #Locate this in loop for psi and phi (Calpha is in middle)
	plot_serials = []
	pos = atomloop.find(pattern, pos)
	if pos == -1: 
			pass
	else: 
		#Add start of pattern before pos to end, to catch overlapping patterns
		atomloop += atomloop[:pos]
	
	while True:
		pattern_serials = []
		if pos == -1: break
		for num, letter in ennumerate(pattern):
			if pos+num >= len(longest_loop): 
				pattern_serials.append(longest_loop[pos+num-len(longest_loop)])
			else: 
				pattern_serials.append(longest_loop[pos+num])
		plot_serials.append(pattern_serials)
		pos = atomloop.find(pattern, pos+1)
		
	#Get x,y,z coordinates to plot from serial numbers	
	plot_coords = []
	for loc in plot_serials:
		i = []
		for serial in loc:
			i.append(atoms[serial]['xyz'])
		plot_coords.append(i)
	
	
	#Break up coords into phi and psi coords 
	#Phi coords need first CNCC of CNCCN
	#Psi coords need NCCN from CNCCN
	phi_coords = []
	psi_coords = []
	for li in plot_coords:
		phi_coords.append([ele for ele in li if li.index(ele) < len(pattern)-1])
		psi_coords.append([ele for ele in li if li.index(ele) >= \
													len(li)-(len(pattern)-1)])

	return phi_coords, psi_coords

def do_ram_plot(phidata, psidata):
	'''Use matplotlib to draw the ramachadran plot for the input peptides'''
	fig = plt.figure()
	ramachandran = fig.add_subplot(1,1,1)
	
	ref = ramachandran.plot(phidata, psidata, 'b.')
	ramachandran.axis([-180,180,-180,180])
	ramachandran.set_ylabel('Psi')
	ramachandran.set_xlabel('Psi')
	ramachandran.set_xticks([-180,0,180])
	ramachandran.set_yticks([-180,0,180])
	ramachandran.grid(True, linestyle='-',markevery=180)
	plt.show()
	
	
def library_ram_plot(pdbfolder, do_max=False, force_recalc=False):
	'''Iterates through folder of pdbs, calculates phi and psi angles for each
		structure using pdb_ramachandran_parse(), and draws plot using do_plot()
		function
		
		Can also write and read phi data to and from a file'''
	
	#Read from file
	try:
		phis, psis = get_saved_data('ramdata.obj',pdbfolder)
		#Test if data is consistent...
		if len(phis) == len(psis) and force_recalc == False:
			print 'Phis, psis read in from file...'
			do_ram_plot(phis, psis)
			return phis, psis
				
	except IOError:
		print 'Phi and psidata file not found - calculating...'
	
	
	if do_max == False:	
		input = glob.glob(os.path.join(pdbfolder, '*.pdb'))
	else:
		try:
			input = random.sample(glob.glob(os.path.join(pdbfolder, '*.pdb')), 
																		do_max)
		except ValueError:
			print "'do_max' must be either False or an integer"
			raise

	phis = []
	psis = []
	anglesperfile = []
	psiequalsphi = []
	
	count = 0
	for f in input:
	#Get phi and psi for peptides in input:
		pepname = os.path.split(f)[1].split('-')[0]
		try:
			count +=1
			phi_coords, psi_coords = pdb_ramachandran_parse(f)
			numangles = len(phi_coords)
			if numangles <= len(pepname) and len(phi_coords) == len(psi_coords):
				anglesperfile.append({'name':f,'num':numangles})
				for philocs, psilocs in zip(phi_coords, psi_coords):
					if philocs != psilocs:
						phis.append(calc_dihedral_angle(philocs))
						psis.append(calc_dihedral_angle(psilocs))
					else:
						psiequalsphi.append(f)
				
			if count % 1000 == 0: print '%s files processed' % count
		except:
			print 'Error processing file %s' % f
			raise
			
	print 'Total processed: %s structures' % count
	
	for i in range(300): #300 is arbitrary limit
		numangles = len([ele for ele in anglesperfile if ele['num'] ==i])
		print 'Number where %s angles found: %s' % (i, numangles)
		if numangles == 0: break
	
	print 'Psi equals phi for %s peptides' % len(psiequalsphi)
	
	#Write to file. If do_max is set - don't bother.
	if do_max == False:
		pickler = write_object((phis,psis),'ramdata.obj',pdbfolder)
		if pickler == True:
			print 'Serialised phi and psi data to file %s' % 'ramdata.obj'
	
	do_ram_plot(phis, psis)	
	return phis, psis
	