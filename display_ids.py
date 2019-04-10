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

import re
import json
import hashlib

modifications =	{ 15995:'oxidation',57021:'carbamidomethyl',42011:'acetyl',31990:'dioxidation',
		 984:'deamidation',43006:'carbamyl',79966:'phosphoryl',79966:'phosphoryl',-17027:'ammonia-loss',
		 18011:'water-loss',6020:'Label:+6 Da',10008:'Label:+10 Da',8014:'Label:+8 Da',4025:'Label:+4 Da'}

#
# generates a simple terminal display for the results of a job
#

def simple_display(_ids,_scores,_spectra,_kernel,_job_stats,_params):
	print('\n1. Identifications:')
	valid_ids = 0
	for j in _ids:
		if len(_ids[j]) == 0:
			print('     %i: no ids' % (j+1))
		else:
			valid_ids += 1
			for i in _ids[j]:
				kern = _kernel[i]
				print('     %i: %s' % (j+1,kern['lb']))
				print('         %i %s %i' % (kern['beg'],kern['seq'],kern['end']))
				line = '         '
				lseq = list(kern['seq'])
				for k in kern['mods']:
					for c in k:
						if k[c] in modifications:
							line += '%s%s+%s;' % (lseq[int(c)-int(kern['beg'])],c,modifications[k[c]])
						else:
							line += '%s%s#%.3f;' % (lseq[int(c)-int(kern['beg'])],c,k[c]/1000)
				line = re.sub('\;$','',line)
				if re.search('[\+\=\#]',line):
					print(line)
				else:
					print('         unmodified')
				print('         ions = %i, Δm = %.3f' % (_scores[j],(_spectra[j]['pm']-(kern['pm']))/1000.0))
	print('    spectra = %i, ids = %i' % (len(_ids),valid_ids))			
	print('\n2. Input parameters:')
	for j in sorted(_params,reverse=True):
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

def tsv_file(_ids,_scores,_spectra,_kernel,_job_stats,_params):
	proton = 1.007276
	print('\n1. Input parameters:')
	for j in sorted(_params,reverse=True):
		print('     %s: %s' % (j,str(_params[j])))
	print('\n2. Job statistics:')
	for j in sorted(_job_stats,reverse=True):
		if j.find('time') == -1:
			print('    %s: %s' % (j,str(_job_stats[j])))
		else:
			print('    %s: %.3f s' % (j,_job_stats[j]))
	ofile = open(_params['output file'],'w')
	if not ofile:
		print('Error: specified output file "%s" could not be opened\n       nothing written to file' % (_of))
		return False
	valid_ids = 0
	line = 'PSM\tspectrum\tscan\trt\tm/z\tz\tprotein\tstart\tend\tsequence\tmodifications\tions\tscore\tdM\tbcid\n'
	ofile.write(line)
	psm = 1
	for j in _ids:
		rt = ''
		scan = ''
		if 'rt' in _spectra[j]:
			rt = '%.1f' % _spectra[j]['rt']
		if 'sc' in _spectra[j]:
			scan = '%i' % _spectra[j]['sc']
		if len(_ids[j]) == 0:
			line = '%i\t%i\t%s\t%s\t%.3f\t%i\t\t\t\t\t\t\t\t\t\n' % (psm,j+1,scan,rt,proton + (_spectra[j]['pm']/1000.0)/_spectra[j]['pz'],_spectra[j]['pz'])
			psm += 1
		else:
			sline = (json.dumps(_spectra[j])).encode()
			valid_ids += 1
			for i in _ids[j]:
				mhash = hashlib.sha256()
				kern = _kernel[i]
				line = '%i\t%i\t%s\t%s\t%.3f\t%i\t%s\t' % (psm,j+1,scan,rt,proton + (_spectra[j]['pm']/1000.0)/_spectra[j]['pz'],_spectra[j]['pz'],kern['lb'])
				psm += 1
				line += '%i\t%i\t%s\t' % (kern['beg'],kern['end'],kern['seq'])
				lseq = list(kern['seq'])
				for k in kern['mods']:
					for c in k:
						if k[c] in modifications:
							line += '%s%s+%s;' % (lseq[int(c)-int(kern['beg'])],c,modifications[k[c]])
						else:
							line += '%s%s#%.3f;' % (lseq[int(c)-int(kern['beg'])],c,k[c]/1000)
				line += '\t%i\t%i\t%.3f' % (_scores[j],100.0*_scores[j]/float(len(kern['seq'])),(_spectra[j]['pm']-(kern['pm']))/1000.0)
				mhash.update(sline+(json.dumps(kern)).encode())
				line += '\t%s\n' % (mhash.hexdigest())
		ofile.write(line)
	ofile.close()
	print('\n3. Output parameters:')
	print('    output file: %s' % (_params['output file']))
	print('    spectra: %i' % (len(_ids)))			
	print('    ids: %i' % (valid_ids))			

