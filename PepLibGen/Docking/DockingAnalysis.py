#!/usr/bin/env python
'''Docking Analysis - analyse cluster docking results

Recommended to import with:
from PepLibGen.Docking import DockingAnalysis as da'''

import os
import glob
import cPickle
import csv
import operator
import shutil
import math
import re
from os import path

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

from PepLibGen.tools import tools
from PepLibGen.Docking import LigplotTools as lig

class VinaLogParseError(Exception): pass
class ParseResultsFoldersError(Exception): pass
class PDBFormatException(Exception):pass


def get_affinity(vina_logfile, num_modes=1):
	'''Takes in a vina log file and extracts the top 
	<num_modes> binding affinity '''

	logfile = open(vina_logfile, 'r')
	allmodes = []
	for line in logfile:
		line = line.strip().split()

		#Results line should consist of all numbers...
		try:
			if line == []: 
				raise ValueError
			line = [float(num) for num in line]
			try:
				mode, affinity = int(line[0]), line[1]
			except IndexError:
				print line
				raise
			allmodes.append((mode, affinity))
		except ValueError:
			pass

	logfile.close()	
	
	if num_modes > 0 and len(allmodes) > 0:
		if num_modes == 1: return allmodes[0]
		else: return allmodes[:num_modes]
	else: 
		raise VinaLogParseError,'No affinity results found in %s' % vina_logfile
			
def idict_to_list(interaction):
	'''Takes a docking interaction of the form {'receptor':receptor,
	'ligand':ligand 'affinity':affinity} and returns a list
	[receptor,ligand,affinity]'''

	return [interaction[key] for key in interaction.keys()]	
		
def lig_receptor_ids(results):
	'''Returns 2 dictionarys, one for ligands, the other for receptors 
	Each will have the names of the ligands/receptors as keys, and a list of 
	indexes of the associated docking interactions in the results as values'''
	ligands = {}
	receptors = {}
	
	for index,interaction in enumerate(results):
		ligand = interaction['ligand']
		receptor = interaction['receptor']
		
		if ligand not in ligands: 
			ligands[ligand] = [index]
		elif ligand in ligands:
			ligands[ligand].append(index)
		
			
		if receptor not in receptors:
			receptors[receptor] = [index]
		elif receptor in receptors:
			receptors[receptor].append(index)

	return receptors, ligands
			
def analyse_results(results, ranked_by='affinity'):
	'''
	Takes in results list and does some analysis - giving each interaction
	a ranking by affinity, per receptor, and a specificity rating for each 
	interaction, low if the ligand does not bind well with other receptors
	
	Can pass other things to ranked_by, including norm_affinity
	'''
	
	receptor_interactions, ligand_interactions = lig_receptor_ids(results)
	rank_key = ranked_by+'_rank'
	rank_specifity_key = ranked_by+'_rank_specifity'
	specifity_key = ranked_by+'_specifity'
	ranked_results = []
		
	print 'Ranking results...'
	i = 1
	
	for rec, interactions in receptor_interactions.items():
	#Build a list of interactions specific to a receptor, sort, and add rank
		try:
			rec_results = sorted([results[inter] for inter in interactions], \
			key=operator.itemgetter(ranked_by))
		except KeyError:
			print len(results)
			for res in results:
				if not ranked_by in res:
					print '%s not in %s' % (ranked_by, str(res))
					raise
			raise
		rank = 1
		for interaction in rec_results:
			interaction[rank_key] = rank
			rank+=1
		ranked_results.extend(rec_results)
		print 'Ranked %s receptors' % i
		i+=1

	#results now include ranks
	results = tuple(ranked_results)
	
	results_withspecificity = []
	tot_num = 0
	for num_ligs,(lig, interactions) in enumerate(ligand_interactions.items()):
		for interaction_index in interactions:
		#Iterate through each interaction index associated with a ligand
		#For each interaction, calculate a specifity which is the smallest 
		#difference in binding rank between the bound receptor, and any other
		#receptor
		#Interactions with high specifity represent ligands that bind to
		#particular receptors without binding well to others
			rank_specificity = \
			min([abs(results[interaction_index][rank_key]-\
			results[index][rank_key]) for index in interactions \
			if index != interaction_index])
			
			affinity_specificity = \
			min([abs(results[interaction_index][ranked_by]-\
			results[index][ranked_by])\
			for index in interactions if index != interaction_index])
			
			interaction = results[interaction_index]
			interaction[rank_specifity_key]=rank_specificity
			interaction[specifity_key]=affinity_specificity
			results_withspecificity.append(interaction)
		if num_ligs % 1000 == 0 and num_ligs != 0: 
			print '%s ligand specificities calculated' % (num_ligs+1)
		tot_num = num_ligs
	print '%s ligand specificities calculated' % (tot_num+1)
		
	return results, receptor_interactions, ligand_interactions	

def seperate_results(top_dir, name, new_folder, test = True):
	'''
	Walks through the directory structure starting with top_dir, and moves 
	any subfolders containing <name> to <new_folder>
	
	Returns a list of moved folders
	'''
	movedirs = []
	for dir, subdirs, files in os.walk(top_dir):
		if name in path.basename(dir): 
			movedirs.append(path.basename(dir))
			if test != True: shutil.move(dir, new_folder)
		
	return movedirs		 

def walk_results_folder(results_directory, max_itercount=0):
	'''
	Walks through results directory, pulling affinity out of files
	Returns a list of {receptor:X,ligand:Y,affinity:Z} values
	'''
	itercount = 0
	results = []
	for dir, subdirs, files in os.walk(results_directory):
		#Find folder with the appropriate log file
		logs = glob.glob(path.join(dir, 'log.*'))
		if len(logs) > 1:
			raise ParseResultsFoldersError, 'More than 1 log file found in %s'\
			% dir
		elif len(logs) == 1: 
			log = logs[0]
			#Get receptor and ligand names	
			receptor, ligand = os.path.split(dir)[1].split('_',1)
			#Get affinity from log file
			mode, affinity = get_affinity(log,num_modes = 1)
			results.append({'receptor':receptor, \
			'ligand':ligand,'affinity':affinity})
			itercount +=1
			if itercount % 1000 == 0:
				print '%s results parsed' % itercount	
		elif len(logs) == 0:
			pass
		if itercount == max_itercount and max_itercount != 0: break
	tools.write_object(results, 'results.obj',results_directory)
	return results
	
def write_results(results, out_dir, out_name='all_docking_results.csv', 
top_out_name='best_docking_results.csv', sort_by='affinity'):
	'''Writes a list of dictionaries (results) to csv files
	
	Will write 2 lists: all results sorted by 'sort_by', and a list of 
	top 10 results for each receptor
	'''
	out_results = list(results) #take a copy
	out_results.sort(key=operator.itemgetter('receptor',sort_by))
	
	outfile = path.join(out_dir, out_name)
	top_outfile = path.join(out_dir,top_out_name)

	#Write all results to csv in out_dir
	filehandle = open(outfile, 'w')
	resultwriter = csv.writer( filehandle )
	resultwriter.writerow(out_results[0].keys())
	
	for interaction in out_results:
		resultwriter.writerow(idict_to_list(interaction))
	filehandle.close()
		
	#Break results into receptors, sort, and write
	filehandle = open(top_outfile, 'w')
	topresultwriter = csv.writer(filehandle)
	#First row will be names of keys, for ease of reading
	topresultwriter.writerow(out_results[0].keys())
	
	#Write top 10 results for each receptor
	count = 0
	receptor = ''
	for interaction in out_results:
	#relies on results being sorted by receptor...
		if receptor != interaction['receptor']:
			receptor = interaction['receptor']
			count = 0
		if count <= 10:
			topresultwriter.writerow(idict_to_list(interaction))
			count += 1
										
	filehandle.close()	
	
def results_histogram(results, filter_by=False, plot_score='combined_score', 
bin_size=50, show=True, norm=False):
	'''
	Uses matplotlib to draw histograms of results data. By default will
	try and plot all 'combined_score' values in the results, with a binsize of 
	50
	
	filter_by will give a pattern to match with results for plotting - e.g if
	SXXX is given, all 4-peptides starting with Serine will be plotted
	
	Will return the data used to construct the plot.
	Can turn of writing of plot by setting show to false
	'''
	data = [] #Holder for histogram data
	if filter_by:
		title = filter_by +' ' +'ligand ' + plot_score
	else:
		title = 'ligand ' + plot_score
	
	#convert filter_by pattern into re pattern
	#e.g. SXXX-SCSCconst would become ^S???-SCSC
	# CXXC-SS11 would become ^C??C-SS
	groups = ['SS','SCNT','SCSC','SCCT','HT']
	try:
		if filter_by == 'lig':
			filter_by = re.compile('.*lig$')
			title = 'Natural ligand combined scores'
		else:
			pattern, bond_def = filter_by.split('-')
			for group in groups:
				if group in bond_def: 
					bond_def = group
			filter_by =\
			re.compile('^' + pattern.replace('X','.') + '-' + bond_def)

	except (ValueError, AttributeError):
		try:
			filter_by = re.compile('^' + filter_by.replace('X','.'))
		except AttributeError: #If filter_by == False
			pass
		

	
	#Extract data from results
	for inter in results:
		if filter_by:
			if filter_by.match(inter['ligand']):
				data.append(inter[plot_score])
		else:
			data.append(inter[plot_score])
	
	if show == True:
		n, bins, patches = plt.hist(data, bin_size, label=title, normed=norm, \
		histtype='barstacked')
		
		# add a 'best fit' line
		mu = np.mean(data)
		sigma = np.std(data)
		y = mlab.normpdf(bins, mu, sigma)
		l = plt.plot(bins, y, 'r--', linewidth=1)
		
		plt.ylabel('Probability density of ligands (%s total)' % len(data))
		plt.xlabel(plot_score)
		#plt.title(title)
		plt.legend(loc='best')
		plt.grid(True)
		
		#plt.savefig('%s_hist.png' % filter_by)
		
		plt.show()
	
	return data
	
def combined_results_histogram(results, plot_by='combined_score', bins=100, 
groups='', title='', norm=False):
		'''
		Draw stacked histograms of results data for each type 
		'''
		plot_data = []
		data_min = -1
		data_max = 1
		if groups == '': 
			groups = ['FDA','lig','CXXC','SXXX-SCSC','SXXX-SCNT','SXXX-SCCT']
		for group in groups:
			plot_data.append(results_histogram(results, filter_by=group, \
			plot_score=plot_by, bin_size=bins, show=False))
			
		labels = [g + ' ' + plot_by for g in groups]
			
		#get ranges of data
		data_min = min([ p[plot_by] for p in results])
		data_max = max([ p[plot_by] for p in results])
		#Hack for nice affinity plots...
		#if plot_by == 'affinity': data_min, data_max = -20,10
			
		n, bins, patches = plt.hist(plot_data, bins, range=(data_min,data_max),\
		histtype='barstacked',label = labels, normed=norm)
				
		plt.ylabel('Number of ligands ')
		plt.xlabel(plot_by)
		plt.title(title)
		plt.legend(loc='best')
		plt.grid(True)
		plt.show()
		
		return plot_data
	
def distributions_bar_chart(results, groups, dist=[0.005,0.01, 0.025], 
score='combined_score', best='low'):
	'''
	Pass in results and will blot a bar chart, with one bar per group to
	show the distribution of results for each group.
	
	The bar will be split into len(dist)+1 segments, with each segment 
	proportional to the number of results in a particular results band, 
	specified by dist
	'''
	
	#Will be a list of lists corresponding to one list per dist segment
	data = [] 
	grouped_results = {}
	width = 0.35
	
	#Will hold the cut-off scores for determining the bars
	cut_scores = []
	
	N = len(groups)
	ind = np.arange(N) #x locations for groups
	
	#Use dist to put together score boundaries for bar chart...
	if best == 'low':
		results.sort(key=operator.itemgetter(score))
	elif best == 'high':
		results.sort(key=operator.itemgetter(score), reverse=True)
	else:
		raise Exception, '"best" must be specified as low or high'
	top = results[1][score]
	
	last_cut = ''
	for i, num in enumerate(dist):
		if i ==0:
			cut_scores.append( (top, results[int(len(results)*num)][score]) )
			last_cut = results[int(len(results)*num)][score] 
		else:
			cut_scores.append((last_cut, results[int(len(results)*num)][score]))
			last_cut = results[int(len(results)*num)][score] 
		
	#Put together results...
	for group in groups:
		filter = '' #Holder for re to match groups in results
		if group == 'lig':
			filter = re.compile('.*lig$')
		elif group == 'true_lig':
			filter = re.compile('.*lig$')
			grouped_results[group] = [res for res in results \
			if filter.match(res['ligand']) \
			and res['ligand'].split('_')[0] == res['receptor'] ]
		else:
			try:
				pattern, bond_def = group.split('-')
				filter = re.compile('^' + pattern.replace('X','.') + \
				'-' + bond_def)
			except ValueError: #No bond_def to split off...
				pattern = group
				filter =  re.compile('^' + pattern.replace('X','.'))	
		#Extract data from results
		if group != 'true_lig':
			grouped_results[group] = \
			[res for res in results if filter.match(res['ligand'])]
			#print len(grouped_results[group])

	for cut in cut_scores:
		#print 'Cut', cut
		fractions = []
		for group in groups:
			if best == 'low':
				tot = [ r for r in grouped_results[group]\
				if r[score] >= cut[0] and r[score] < cut[1] ]
				percent = \
				(float(len(tot)) / float(len(grouped_results[group]))) * 100
				#print 'Percent %s group %s' % (percent, group)
				fractions.append( percent )
			elif best == 'high':
				tot = [ r for r in grouped_results[group] \
				if r[score] <= cut[0] and r[score] > cut[1] ]
				percent = \
				(float(len(tot)) / float(len(grouped_results[group]))) * 100
				#print 'Percent %s group %s' % (percent, group)
				fractions.append( percent )
		data.append(fractions)
		
	plots = []	
	colours = ['b','g','r','c','m','y','k']
	labels = [ 'Top ' + str(s*100) + '% of ligands by '+score for s in dist ]
	
	bottoms = [0] * len(groups)
	for i, point in enumerate(data):
		#print bottoms, point
		p = plt.bar(ind, point, width=0.35, label=labels[i], color=colours[i], bottom=bottoms)
		#Add points to bottoms
		bottoms = [ sum(a) for a in zip(*[bottoms, point]) ]
	
	plt.ylabel('Percentage of ligands')
	plt.xlabel('Groups')
	plt.title('Distribution of top docking scores by ligand type')
	plt.xticks(ind+width/2., groups)
	plt.legend(loc='best')
	
	plt.show()
	
	return data
	

			
def get_results(top_dir, max_itercount=0, force_recalc=False, normalise=True,
analyse=True, write=True):
	'''Walks through all directories from the top directory looking for files 
    called log.*. It will then try to parse the file as a vina log file,
	read the receptor name and ligand name from the containing folder 
	Will attempt to write a python object containing the results data
	if sucessful, and will attempt to read any existing python object when 
	being run, unless the force_recalc method is not False
	
	Can optionally normalise results with the normalise_results function, and/or
	analyse results with analyse_results.
	
	Returns a list of tuples of (receptor, ligand, best_affinity)
	
	Will try to read the results to a python object file - can be overridden 
	with the force_recalc option. Will write results to an opject file called
	results.obj in the top_dir, or norm_results.obj if normalised'''
	
	objfile = 'results.obj'
	ligplot_obfile='ligplot_scores.obj'
	final_objfile = 'final_results.obj'
	weight = 0.5 #For combining scores...

	#Look for current results
	try:
		if max_itercount == 0 and force_recalc == False:
			results = tools.get_saved_data(objfile,top_dir)
			print 'Results loaded from %s' % path.join(top_dir,objfile)
		else:
			print 'Parsing from root directory %s' % top_dir
			force_recalc=True
	except IOError:
		print 'No previously written results found '
		print 'Parsing from root directory %s' % top_dir
		force_recalc=True
	if force_recalc:
		results = walk_results_folder(top_dir, max_itercount)
	#Sort results
	results = list(results)
	results.sort(key=operator.itemgetter('receptor','ligand'))
	#New normalisation = get normalised data, and add it to non-normalised
	#Take ligand and receptor names, and scores, pass to numpy  normaliser
	#Pass in  list of (receptor, ligand, affinity) to normaliser, get back
	#(receptor, ligand, normalised_affinity)
	
	print 'Len docking results', len(results)
	
	if normalise:
		norm_holder = []
		for interaction in results:
			norm_holder.append((interaction['receptor'],\
			interaction['ligand'],interaction['affinity']))
		print 'Normalising docking scores'
		norm_scores = tools.normaliser(norm_holder)
		for normed_inter, inter in zip(norm_scores, results):
			inter['norm_affinity'] = normed_inter[2]

	#Try and add ligplot scores (and normalise)...
	try: 
		#Try and get written Ligplot results...
		print 'Attempting to read from %s...' % \
		path.join(top_dir,ligplot_obfile)
		ligplot_scores = tools.get_saved_data(ligplot_obfile, top_dir)
		if len(ligplot_scores) != len(results):
			raise ParseResultsFoldersError, \
			'Inconsistent docking results vs. ligplot scores %s vs %s'\
			% (len(results), len(ligplot_scores))
		ligplot_scores.sort(key=operator.itemgetter('receptor','ligand'))
		for ligplot_inter, inter in zip(ligplot_scores, results):
			inter['ligplot_score'] = ligplot_inter['ligplot_score']
		if normalise:
		#Normalise ligplot scores too, before adding
			norm_holder_lig = []
			for interaction in ligplot_scores:
				norm_holder_lig.append((interaction['receptor'],\
				interaction['ligand'],interaction['ligplot_score']))
			print 'Normalising ligplot scores'
			norm_lig_scores = tools.normaliser(norm_holder_lig)
			for normed_ligplot, inter in zip(norm_lig_scores, results):
				inter['norm_ligplot_score'] = normed_ligplot[2]
				#Also add in combined normalised docking and ligplot scores
				inter['combined_score'] = \
				weight*inter['norm_affinity'] + \
				(1-weight)*(-1)*inter['norm_ligplot_score']
	except IOError:
		print 'Could not find ligplot results - ignoring'

	if analyse:
		for analysis in ['affinity','norm_affinity', 'norm_ligplot_score',\
		'combined_score', 'ligplot_score']:
			results, receptor_interactions, ligand_interactions \
			= analyse_results(results, ranked_by=analysis)
	
	if max_itercount == 0:
		written = tools.write_object(results, final_objfile, top_dir)
		if written == True:
			print 'Wrote %s results to file %s' % (len(results), final_objfile)
		
	print 'Writing results to %s' % top_dir
	write_results(results, top_dir, sort_by='combined_score')
	 				
	return results 


	
def combine_final_results(*args, **kwargs):
	'''
	Pass in final results lists with ligplot scores and affinity data and will
	renormalise, analyse, and write.
	'''
	#out_file='comb_final_results.obj', out_dir='', 
	out_results = []
	lig_norm_holder =[]
	dock_norm_holder = []
	out_file = 'comb_final_results.obj'
	out_dir = ''
	
	print kwargs
	
	arg_keys = kwargs.keys()
	if 'out_file' in arg_keys:
		out_file = kwargs['out_file']
	if 'out_dir' in arg_keys:
		out_dir = kwargs['out_dir']
	
	
	for res in args:
		out_results.extend(res)
		
	out_results.sort(key=operator.itemgetter('receptor','ligand'))
		
	#Normalise affinities and ligplot scores...
	for interaction in out_results:
		dock_norm_holder.append((interaction['receptor'],interaction['ligand'],\
		interaction['affinity']))
	print 'Normalising docking scores'
	dock_norm_scores = tools.normaliser(dock_norm_holder)
	for normed_inter, inter in zip(dock_norm_scores, out_results):
		inter['norm_affinity'] = normed_inter[2]
	
	for interaction in out_results:
		lig_norm_holder.append((interaction['receptor'],interaction['ligand'],\
		interaction['ligplot_score']))
	print 'Normalising ligplot scores'
	ligplot_norm_scores = tools.normaliser(lig_norm_holder)
	for normed_inter, inter in zip(ligplot_norm_scores, out_results):
		inter['norm_ligplot_score'] = normed_inter[2]
		
	tools.write_object(out_results, out_file, out_dir)
	return out_results
		
		
		
		
		