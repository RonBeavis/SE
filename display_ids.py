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
from scipy.stats import hypergeom,tmean,tstd
import math
import copy

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
# generates a TSV file for the results of a job
#
def display_parameters(_params):
	print('\nInput parameters:')
	summary = ''
	for j in sorted(_params,reverse=True):
		summary += '     %s: %s\n' % (j,str(_params[j]))
	print(summary)
	f = open(_params['output file'] + '.summary','w')
	ss = [re.sub('^(.+?)\: +',r'\1	',l.lstrip()) for l in summary.split('\n')]
	f.write('\n'.join(ss))
	f.close()

def find_limits(_w,_ids,_spectra,_kernel,_st,_mins):
	bins = {}
	for j in _ids:
		if not _ids[j]:
			continue
		for i in _ids[j]:
			if (j,i) not in _st:
				continue
			if _st[(j,i)] < _mins:
				continue
			kern = _kernel[i]
			if kern['lb'].find('decoy-') != -1:
				continue
			delta = int(_spectra[j]['pm']-kern['pm'])
			if abs(delta) <= _w:
				delta = int(0.5 + 1.0e6*delta/_spectra[j]['pm'])
				if delta in bins:
					bins[delta] += 1
				else:
					bins[delta] = 1
			break
	max_bin = 0
	max_m = min(bins)
	first = -100
	last = 99
	for m in bins:
		if max_bin < bins[m]:
			max_bin = bins[m]
			max_m = m
	try:
		first = min(bins)
		last = max(bins)
	except:
		return (first,last+1,None)
	if max_bin < 200:
		return (first,last+1,None)
	min_bin = int(0.5 + float(max_bin)/100.0)
	low = None
	high = last
	a = first
	while a <= last:
		if a not in bins:
			bins[a] = 0
		a += 1
	m = max_m
	low = first
	while m >= first+1:
		if bins[m] <= min_bin and bins[m-1] <= min_bin:
			low = m
			break
		m -= 1
	low += 1
	high = last
	m = max_m
	while m <= last - 1:
		if bins[m] <= min_bin and bins[m+1] <= min_bin:
			high = m
			break
		m += 1
	high -= 1
	print(low,high,bins)
	return (low,high,bins)

def generate_scores(_ids,_scores,_spectra,_kernel,_params):
	res = _params['fragment mass tolerance']
	sfactor = 20
	sadjust = 1
	if res > 100:
		sfactor = 40
	sd = {}
	for j in _ids:
		p_score = 0.0
		if not _ids[j]:
			continue
		for i in _ids[j]:
			kern = _kernel[i]
			lseq = list(kern['seq'])
			pmass = int(kern['pm']/1000)
			cells = int(pmass-200)
			if cells > 1500:
				cells = 1500
			total_ions = 2*(len(lseq) - 1)
			if total_ions > sfactor:
				total_ions = sfactor
			if total_ions < _scores[j]:
				total_ions = _scores[j] + 1
			sc = len(_spectra[j]['sms'])/3
			if _scores[j] >= sc:
				sc = _scores[j] + 2
			rv = hypergeom(cells,total_ions,sc)
			p = rv.pmf(_scores[j])
			pscore = -100.0*math.log10(p)*sadjust
			sd[(j,i)] = pscore
	return sd

def create_header():
	return 	'PSM\tspectrum\tscan\trt\tm/z\tz\tprotein\tstart\tend\tpre\tsequence\tpost\tmodifications\tions\tscore\tsignal\tdM\tppm\tn\tsav\trs\tmaf'

def tsv_file(_ids,_stuples,_spectra,_kernel,_job_stats,_params):
	if len(_ids) == 0:
		ofile = open(_params['output file'],'w')
		if not ofile:
			print('Error: specified output file "%s" could not be opened\n       nothing written to file' % (_params['output file']))
			return False
		ofile.write(create_header() + '\n')
		print('\n2. Output parameters:')
		print('    output file: %s' % (_params['output file']))
		print('    PSMs: %i' % (0))
		ofile.close()
		return
	_scores = {}
	_intensities = {}
	for st in _stuples:
		_scores[st] = _stuples[st][0]
		_intensities[st] = _stuples[st][1]
	proteins = set([])
	pscore_min = 200.0
	print('     applying statistics')
	score_tuples = generate_scores(_ids,_scores,_spectra,_kernel,_params)
	(low,high,bins) = find_limits(int(_params['parent mass tolerance']),_ids,_spectra,_kernel,score_tuples,pscore_min)
	_params['output low ppm'] = low
	_params['output high ppm'] = high
	_params['output histogram ppm'] = bins
	outfile = _params['output file']
	print('     storing results in "%s"' % (outfile))
	modifications = get_modifications()
	proton = 1.007276
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
	line = create_header()
	if use_bcid:
		line += '\tbcid'
	line += '\n'
	ofile.write(line)
	psm = 1
	z_list = {}
	ptm_list = {}
	ptm_aaa = {}
	unique_psms = set([])
	parent_delta = []
	parent_delta_ppm = []
	parent_a = [0,0]
	pscore = 0.0
	vresults = 0
	res = _params['fragment mass tolerance']
	minimum_intensity = _params['minimum identified intensity']
	sfactor = 20
	sadjust = 1
	PSMs = 0
	SAVs = 0
	DECOYs = 0
	sav_mafs = {}
	if res > 100:
		sfactor = 40
		sadjust = 0.5
	for j in _ids:
		pscore = 0.0
		rt = ''
		scan = ''
		line = ''
		if 'rt' in _spectra[j]:
			rt = '%.1f' % _spectra[j]['rt']
		if 'sc' in _spectra[j]:
			scan = '%i' % _spectra[j]['sc']
		if len(_ids[j]) == 0:
			if valid_only:
				psm += 1
				continue
			line = '%i\t%i\t%s\t%s\t%.3f\t%i\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n' % (psm,j+1,scan,rt,proton + (_spectra[j]['pm']/1000.0)/_spectra[j]['pz'],_spectra[j]['pz'])
			psm += 1
		else:
			sline = (json.dumps(_spectra[j])).encode()
			vresults = 0
			pscore = 0.0
			line = ''
			x = 0
			for i in _ids[j]:
				kn = _kernel[i]
				kerns = [kn]
				if 'vlb' in kn:
					for a,lb in enumerate(kn['vlb']):
						nk = copy.deepcopy(kn)
						nk['lb'] = lb
						nk['pre'] = kn['vpre'][a]
						nk['post'] = kn['vpost'][a]
						nk['beg'] = kn['vbeg'][a]
						nk['end'] = kn['vend'][a]
						kerns.append(nk)
				for kern in kerns:
					line = ''
					lseq = list(kern['seq'])
					pmass = int(kern['pm']/1000)
					pscore = score_tuples[(j,i)]*_intensities[j]/100.0
					if pscore < pscore_min and valid_only:
						break
					delta = _spectra[j]['pm']-kern['pm']
					ppm = 1e6*delta/kern['pm']
					if delta/1000.0 > 0.9:
						ppm = 1.0e6*(delta-1003.0)/kern['pm']
					if ppm < low or ppm > high:
						continue
					if _intensities[j] < minimum_intensity:
						continue
					if x == 0 and delta/1000.0 > 0.9:
						parent_a[1] += 1
						parent_delta_ppm.append(ppm)
					elif x == 0:
						parent_a[0] += 1
						parent_delta.append(delta/1000.0)
						parent_delta_ppm.append(ppm)
					x += 1
					valid_ids += 1
					z = _spectra[j]['pz']
					if z in z_list:
						z_list[z] += 1
					else:
						z_list[z] = 1
					mhash = hashlib.sha256()
					lb = kern['lb']
					if lb.find('decoy-') == 0:
						DECOYs += 1
					unique_psms.add(scan)
					proteins.add(lb)
					line = '%i\t%i\t%s\t%s\t%.3f\t%i\t%s\t' % (psm,j+1,scan,rt,proton + (_spectra[j]['pm']/1000.0)/_spectra[j]['pz'],_spectra[j]['pz'],lb)
					psm += 1
					line += '%i\t%i\t%s\t%s\t%s\t' % (kern['beg'],kern['end'],kern['pre'],kern['seq'],kern['post'])
					for k in kern['mods']:
						for c in k:
							diff = kern['beg'] - kn['beg']
							if k[c] in modifications:
								aa = lseq[int(c+diff)-int(kern['beg'])]
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

								line += '%s%s+%s;' % (aa,c+diff,modifications[k[c]])
							else:
								ptm = '%.3f' % (float(k[c])/1000.0)
								aa = lseq[int(c+diff)-int(kern['beg'])]
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
								line += '%s%s#%.3f;' % (aa,c+diff,float(k[c])/1000)
					line = re.sub(';$','',line)
					line += '\t%i\t%.0f\t%.0f\t%.3f\t%i' % (_scores[j],pscore,_intensities[j],delta/1000,round(ppm,0))
					line += '\t%i' % (sum(kern['ns']))
					if 'sav' in kern:
						line += '\t%s%i%s\t%s\t%.2f' % (kern['res'],kern['pos'],kern['sav'],kern['rsn'],kern['maf'])
						sav_mafs[kern['rsn']] = kern['maf']
						SAVs += 1
					else:
						line += '\t\t\t'
					mhash.update(sline+(json.dumps(kern)).encode())
					if use_bcid:
						line += '\t%s' % (mhash.hexdigest())
					line += '\n'
					ofile.write(line)
					PSMs += 1
				
	ofile.close()
	hist = [0]*101
	for a in _intensities:
		v = int(_intensities[a])
		if v <= 100:
			hist[v] += 1
	hist[0] = 0
	v = 0
	total = sum(hist)
	int_hist = []
	for v in range(100):
			int_hist.append((v,hist[v],float(sum(hist[0:v])/total)))
			
	summary = ''
	if PSMs == 0:
		print('\n2. Output parameters:')
		summary += '    output file: %s\n' % (_params['output file'])
		summary += '    PSMs: %i\n' % (PSMs)
		print(summary)
		f = open(_params['output file'] + '.summary','a')
		ss = [l.lstrip() for l in summary.split('\n')]
		f.write('\n'.join(ss))
		f.close()
		return

	print('\n2. Output parameters:')
	summary += '    output file: %s\n' % (_params['output file'])
	summary += '    PSMs:\n'
	summary += '          total: %i\n' % (PSMs)
	summary += '          unique: %i\n' % (len(unique_psms))
	summary += '    proteins: %i\n' % len(proteins)
	summary += '    parent ppms: (%i,%i)\n' % (_params['output low ppm'],_params['output high ppm'])
	summary += '    charges:\n'
	for z in sorted(z_list):
		summary += '          %i: %i\n' % (z,z_list[z])
	summary += '    modifications:\n'
	for ptm in sorted(ptm_list, key=lambda s: s.casefold()):
		aa_line = ''
		for aa in sorted(ptm_aaa[ptm]):
			aa_line += '%s[%i] ' % (aa,ptm_aaa[ptm][aa])
		summary += '          %s: %s= %i\n' % (ptm,aa_line,ptm_list[ptm])
	if DECOYs > 0:
		summary += '    decoys:\n'
		summary += '       total: %i\n' % (DECOYs)

	summary += '    SAVs:\n'
	summary += '       total: %i\n' % (SAVs)
	if SAVs > 0:
		summary += '       unique: %i\n' % (len(sav_mafs))
		power = 1.0
		for maf in sav_mafs:
			if sav_mafs[maf] is not None and sav_mafs[maf] != 0.0:
				power *= sav_mafs[maf]
		summary += '       power: %.2e:1\n' % (1.0/power)
	if len(parent_delta) > 10:
		summary += '    parent delta mean (Da): %.3f\n' % (tmean(parent_delta))
		summary += '    parent delta sd (Da): %.3f\n' % (tstd(parent_delta))
		summary += '    parent delta mean (ppm): %.1f\n' % (tmean(parent_delta_ppm))
		summary += '    parent delta sd (ppm): %.1f\n' % (tstd(parent_delta_ppm))
	total = float(parent_a[0]+parent_a[1])
	if total > 0:
		summary += '    parent A: A0 = %i (%.1f%%), A1 = %i (%.1f%%)\n' % (parent_a[0],100*parent_a[0]/total,parent_a[1],100*parent_a[1]/total)
	else:
		summary += '    parent A: A0 = %i, A1 = %i' % (parent_a[0],parent_a[1])
	if len(_intensities) > 10:
		arr = [_intensities[a] for a in _intensities if _intensities[a] >= minimum_intensity]
		summary += '    signal mean (%%): %.1f\n' % (tmean(arr))
		summary += '    signal sd (%%): %.1f\n' % (tstd(arr))
	print(summary)
	f = open(_params['output file'] + '.summary','a')
	ss = [re.sub('^(.+?)\: +',r'\1	',l.lstrip()) for l in summary.split('\n')]
	f.write('\n'.join(ss))
	f.close()
	return

