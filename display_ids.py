#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
import re

modifications =	{ 15995:'oxidation',57021:'carbamidomethyl',42011:'acetyl',31990:'dioxidation',
		 984:'deamidation',43006:'carbamyl',79966:'phosphoryl',79966:'phosphoryl',-17027:'ammonia-loss',
		 18011:'water-loss',6020:'Label:+6 Da',10008:'Label:+10 Da',8014:'Label:+8 Da',4025:'Label:+4 Da'}

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
				print('         score = %i, Δm = %.3f' % (_scores[j],(_spectra[j]['pm']-(kern['pm']))/1000.0))
	print('    spectra = %i, ids = %i' % (len(_ids),valid_ids))			
	print('\n2. Input parameters:')
	for j in _params:
		print('     %s: %s' % (j,str(_params[j])))
	print('\n3. Job statistics:')
	for j in _job_stats:
		if j.find('time') == -1:
			print('    %s: %s' % (j,str(_job_stats[j])))
		else:
			print('    %s: %.3f s' % (j,_job_stats[j]))

