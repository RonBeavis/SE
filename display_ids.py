#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

def simple_display(_ids,_scores,_spectra,_kernel,_job_stats,_params):
	print('\n1. Identifications:')
	for j in _ids:
		if len(_ids[j]) == 0:
			print('     %i: no ids' % (j+1))
		else:
			for i in _ids[j]:
				t = _kernel[i]
				print('     %i: %s, %i %s %i, score = %i' % (j+1,t['lb'],t['beg'],t['seq'],t['end'],_scores[j]))
	print('\n2. Input parameters:')
	for j in _params:
		print('     %s: %s' % (j,str(_params[j])))
	print('\n3. Job statistics:')
	for j in _job_stats:
		if j.find('time') == -1:
			print('     %s: %s' % (j,str(_job_stats[j])))
		else:
			print('     %s: %.3f s' % (j,_job_stats[j]))

