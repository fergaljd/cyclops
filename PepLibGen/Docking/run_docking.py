
#!/usr/bin/env python
# Script to carry out docking of a peptide library with Vina.
#python run_docking.py <ligand dir> <receptor dir>

import os
import multiprocessing
import glob
import subprocess
import sys
import time

class ReceptorFilesError(Exception): pass
class InputArgumentsError(Exception):pass

num_cpus = multiprocessing.cpu_count()
usage = 'Usage: python run_docking.py <ligand_dir> <receptor_dir>'
if len(sys.argv) != 3: 
	raise InputArgumentsError, usage

# #################            PROCESS INPUT FILES         #############################	
	
#Directory containing PDBQT files
if os.path.exists(sys.argv[1]): pdbqt_dir = sys.argv[1]
else: raise InputArgumentsError, '%s is not a valid directory' % sys.argv[1]


#Directory containing a directory for each receptor PDB, containing structure file, and 
#vina config file
if os.path.exists(sys.argv[2]): receptor_dir = sys.argv[2]
else: raise InputArgumentsError, '%s is not a valid directory' % sys.argv[2]



#Build list of absolute paths of PDBQT files...
pdbqt_files = [pdbqt for pdbqt in glob.glob(pdbqt_dir+os.path.sep+'*.pdbqt') ]
				


#Build list of absolute paths of receptor PDB and config (in a tuple for each receptor)
receptor_pdb_files = [pdb for pdb in glob.glob(receptor_dir+os.path.sep+'*.pdbqt')]

receptor_config_files=[]
for pdb in receptor_pdb_files:
	found_conf = False
	for conf_file in glob.glob(receptor_dir + os.path.sep + '*.conf'):
		if conf_file.split('.')[0] == pdb.split('.')[0]:
			receptor_config_files.append(conf_file)
			found_conf = True
	if found_conf == False:
		raise ReceptorFilesError, 'No config file found for %s', pdb


	
receptor_files = zip(receptor_pdb_files, receptor_config_files)
if len(receptor_files) == 0:
	raise InputArgumentsError,\
	'No receptor files found in input directory %s' % receptor_dir
elif len(pdbqt_files) == 0:
	raise InputArgumentsError,\
	'No ligand .pdbqt files found in input directory %s' % pdbqt_dir

	
	
out_dir = receptor_dir + 'docking_output'
if os.path.exists(out_dir) ==False: os.makedirs(out_dir)

print '\n\nBegan Docking at', time.ctime()

# ###############                     DOCKING          ################################

receptorcount = 0
for receptor in receptor_files:
	receptorcount+=1
	ligcount = 0
	for ligand in pdbqt_files:
		lig_name = os.path.split(ligand)[1].split('.')[0]
		rec_name = os.path.split(receptor[0])[1].split('.')[0]
		
		print 'Docking %s to %s' % (lig_name,rec_name)
		ligcount +=1
		print 'Docking ligand %s of %s to receptor %s of %s: %s' \
			% (ligcount, len(pdbqt_files), \
			receptorcount, len(receptor_files), time.ctime())
		
		
		receptorf = '--receptor '+ receptor[0]
		configf = '--config '+receptor[1]
		ligandf =  '--ligand ' +ligand
		outf = '--out ' + out_dir + os.path.sep + rec_name+'_'+ lig_name+'.pdbqt'
		logf = '--log ' + out_dir + os.path.sep + rec_name+'_'+ lig_name+'.log.txt'
				
		runstring ='vina'+' '+configf+' '+receptorf+' '+ligandf+' '+outf+' '+logf+' '
		try:
			print runstring
			docking =subprocess.Popen(runstring,
							stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
			prog_output = docking.communicate()
			docking.wait()
			if docking.returncode != 0: 
				print prog_output[1]
				
			else: print 'Docking complete'
		except OSError:
			print 'Could not run program - check current working directory and path'
			raise


#Attempt to tar up results... use command tar -czxf <outfile.tgz> <directory>
result_dir = os.path.split(out_dir)[0]
tar_runstring = 'tar -czf '+result_dir+os.path.sep+'docking_results.tgz ' + out_dir

subprocess.call(tar_runstring, shell=True)
		