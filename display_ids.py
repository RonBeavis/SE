#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# displays the results of a job, either to the terminal or a file
#

#
# some handy lists of masses, in integer milliDaltons
#
#from __future__ import print_function
#from libcpp cimport bool as bool_t

import re
import json
import hashlib
import statistics
import scipy.stats
import math

#
# retrieves a list of modification masses and names for use in displays
# if 'common_mods_md.txt' is not available, a warning is thrown and a default
# list is used
#
def get_modifications():
	try:
		f = open('common_mods_md.txt','r')
	except:
		print('Warning: common_mods_md.txt is not available so default values used')
		return { 15995:'Oxidation',57021:'Carbamidomethyl',42011:'Acetyl',31990:'Dioxidation',
		28031:'Dimethyl',14016:'Methyl',984:'Deamidation',43006:'Carbamyl',79966:'Phosphoryl',
		-17027:'Ammonia-loss',-18011:'Water-loss',6020:'Label:+6 Da',
		10008:'Label:+10 Da',8014:'Label:+8 Da',4025:'Label:+4 Da'}
	mods = {}
	for l in f:
		l = l.strip()
		vs = l.split('\t')
		if len(vs) < 2:
			continue
		mods[int(vs[0])] = vs[1]
	f.close()
	return mods
#
# generates a simple terminal display for the results of a job
#

def simple_display(_ids,_scores,_spectra,_kernel,_job_stats,_params):
	ko = _params['kernel order'].copy()
	modifications = get_modifications()
	print('\n1. Identifications:')
	valid_ids = 0
	for j in _ids:
		if len(_ids[j]) == 0:
			print('     %i: no ids' % (j+1))
		else:
			valid_ids += 1
			for i in _ids[j]:
				kern = _kernel[i]
				print('     %i: %s' % (j+1,kern[ko['lb']]))
				print('         %i %s %i' % (kern[ko['beg']],kern[ko['seq']],kern[ko['end']]))
				line = '         '
				lseq = list(kern[ko['seq']])
				for k in kern['mods']:
					for c in k:
						if k[c] in modifications:
							pos = int(c)-int(kern[ko['beg']])
							line += '%s%s+%s;' % (lseq[pos],c,modifications[k[c]])
						else:
							pos = int(c)-int(kern[ko['beg']])
							line += '%s%s#%.3f;' % (lseq[pos],c,k[c]/1000)
				line = re.sub('\;$','',line)
				if re.search('[\+\=\#]',line):
					print(line)
				else:
					print('         unmodified')
				print('         ions = %i, Δm = %.3f' % (_scores[j],(_spectra[j]['pm']-(kern[ko['pm']]))/1000.0))
	print('    spectra = %i, ids = %i' % (len(_ids),valid_ids))			
	print('\n2. Input parameters:')
	for j in sorted(_params,reverse=True):
		if j == 'kernel order':
			continue
		print('     %s: %s' % (j,str(_params[j])))
	print('\n3. Job statistics:')
	for j in sorted(_job_stats,reverse=True):
		if j.find('time') == -1:
			print('    %s: %s' % (j,str(_job_stats[j])))
		else:
			print('    %s: %.3f s' % (j,_job_stats[j]))
#
# generates a TSV file for the results of a job
#
def display_parameters(_params):
	print('\nInput parameters:')
	for j in sorted(_params,reverse=True):
		if j == 'kernel order':
			continue
		print('     %s: %s' % (j,str(_params[j])))

def tsv_file(_ids,_scores,_spectra,_kernel,_job_stats,_params):
	ko = _params['kernel order'].copy()
	modifications = get_modifications()
	proton = 1.007276
	outfile = _params['output file']
	print('\n1. Job statistics:')
	for j in sorted(_job_stats,reverse=True):
		if j.find('time') == -1:
			print('    %s: %s' % (j,str(_job_stats[j])))
		else:
			print('    %s: %.3f s' % (j,_job_stats[j]))
	ofile = open(outfile,'w')
	if not ofile:
		print('Error: specified output file "%s" could not be opened\n       nothing written to file' % (outfile))
		return False
	valid_only = False
	if 'output valid only' in _params:
		valid_only = _params['output valid only']
	use_bcid = False
	if 'output bcid' in _params:
		use_bcid = _params['output bcid']
	valid_ids = 0
	line = 'PSM\tspectrum\tscan\trt\tm/z\tz\tprotein\tstart\tend\tpre\tsequence\tpost\tmodifications\tions\tscore\tdM\tppm'
	if use_bcid:
		line += '\tbcid'
	line += '\n'
	ofile.write(line)
	psm = 1
	z_list = {}
	ptm_list = {}
	ptm_aaa = {}
	parent_delta = []
	parent_delta_ppm = []
	parent_a = [0,0]
	pscore = 0.0
	pscore_min = 200.0
	vresults = 0
	res = _params['fragment mass tolerance']
	sfactor = 20
	sadjust = 1
	if res > 100:
		sfactor = 40
		sadjust = 0.5
	for j in _ids:
		pscore = 0.0
		rt = ''
		scan = ''
		line = '-----------------------------------------------------------------------'
		if 'rt' in _spectra[j]:
			rt = '%.1f' % _spectra[j]['rt']
		if 'sc' in _spectra[j]:
			scan = '%i' % _spectra[j]['sc']
		if len(_ids[j]) == 0:
			if valid_only:
				psm += 1
				continue
			line = '%i\t%i\t%s\t%s\t%.3f\t%i\t\t\t\t\t\t\t\t\t\n' % (psm,j+1,scan,rt,proton + (_spectra[j]['pm']/1000.0)/_spectra[j]['pz'],_spectra[j]['pz'])
			psm += 1
		else:
			sline = (json.dumps(_spectra[j])).encode()
			vresults = 0
			pscore = 0.0
			line = '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
			for i in _ids[j]:
				kern = _kernel[i]
				lseq = list(kern[ko['seq']])
				pmass = int(kern[ko['pm']]/1000)
				cells = int(pmass-200)
				total_ions = 2*(len(lseq) - 1)
				if total_ions > sfactor:
					total_ions = sfactor
				if total_ions < _scores[j]:
					total_ions = _scores[j] + 1
				rv = scipy.stats.hypergeom(cells,total_ions,len(_spectra[j]['sms'])/3)
				p = rv.pmf(_scores[j])
				pscore = -100.0*math.log10(p)*sadjust
				if pscore < pscore_min and valid_only:
					break
				valid_ids += 1
				z = _spectra[j]['pz']
				if z in z_list:
					z_list[z] += 1
				else:
					z_list[z] = 1
				mhash = hashlib.sha256()
				lb = kern[ko['lb']]
				if lb.find('-rev') != -1:
					lb = '-rev'
				line = '%i\t%i\t%s\t%s\t%.3f\t%i\t%s\t' % (psm,j+1,scan,rt,proton + (_spectra[j]['pm']/1000.0)/_spectra[j]['pz'],_spectra[j]['pz'],lb)
				psm += 1
				line += '%i\t%i\t%s\t%s\t%s\t' % (kern[ko['beg']],kern[ko['end']],kern[ko['pre']],kern[ko['seq']],kern[ko['post']])
				for k in kern[ko['mods']]:
					for c in k:
						if k[c] in modifications:
							aa = lseq[int(c)-int(kern[ko['beg']])]
							ptm = modifications[k[c]]
							if ptm in ptm_list:
								ptm_list[ptm] += 1
							else:
								ptm_list[ptm] = 1
							if ptm in ptm_aaa:
								if aa in ptm_aaa[ptm]:
									ptm_aaa[ptm][aa] += 1
								else:
									ptm_aaa[ptm].update({aa:1})
							else:
								ptm_aaa[ptm] = {aa:1}

							line += '%s%s+%s;' % (aa,c,modifications[k[c]])
						else:
							ptm = '%.3f' % float(k[c])/1000.0
							aa = lseq[int(c)-int(kern[ko['beg']])]
							if ptm in ptm_list:
								ptm_list[ptm] += 1
							else:
								ptm_list[ptm] = 1
							if ptm in ptm_aaa:
								if aa in ptm_aaa[ptm]:
									ptm_aaa[ptm][aa] += 1
								else:
									ptm_aaa[ptm].update({aa:1})
							else:
								ptm_aaa[ptm] = {aa:1}
							line += '%s%s#%.3f;' % (aa,c,float(k[c])/1000)
				delta = _spectra[j]['pm']-kern[ko['pm']]
				ppm = 1e6*delta/kern[ko['pm']]
				if delta/1000.0 > 0.9:
					parent_a[1] += 1
					ppm = 1.0e6*(delta-1003.0)/kern[ko['pm']]
					parent_delta_ppm.append(ppm)
				else:
					parent_a[0] += 1
					parent_delta.append(delta/1000.0)
					parent_delta_ppm.append(ppm)
				line += '\t%i\t%.0f\t%.3f\t%i' % (_scores[j],pscore,delta/1000,round(ppm,0))
				mhash.update(sline+(json.dumps(kern)).encode())
				if use_bcid:
					line += '\t%s' % (mhash.hexdigest())
				line += '\n'
				ofile.write(line)
	ofile.close()
	print('\n2. Output parameters:')
	print('    output file: %s' % (_params['output file']))
	print('    PSMs: %i' % (valid_ids))
	print('    charges:')
	for z in sorted(z_list):
		print('        %i: %i' % (z,z_list[z]))
	print('    modifications:')
	for ptm in sorted(ptm_list, key=lambda s: s.casefold()):
		aa_line = ''
		for aa in sorted(ptm_aaa[ptm]):
			aa_line += '%s[%i] ' % (aa,ptm_aaa[ptm][aa])
		print('        %s: %s= %i' % (ptm,aa_line,ptm_list[ptm]))
	print('    parent delta mean (Da): %.3f' % (statistics.mean(parent_delta)))
	print('    parent delta sd (Da): %.3f' % (statistics.stdev(parent_delta)))
	print('    parent delta mean (ppm): %.1f' % (statistics.mean(parent_delta_ppm)))
	print('    parent delta sd (ppm): %.1f' % (statistics.stdev(parent_delta_ppm)))
	print('    parent A: A0 = %i, A1 = %i' % (parent_a[0],parent_a[1]))

