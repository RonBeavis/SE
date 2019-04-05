#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

def simple_display(ids,spectra,kernel,job_stats,params):
	print('\n1. Identifications:')
	for j in ids:
		if len(ids[j]) == 0:
			print('     %i: no ids' % (j+1))
		else:
			for i in ids[j]:
				t = kernel[i]
				print('     %i: %s, %i %s %i' % (j+1,t['lb'],t['beg'],t['seq'],t['end']))
	print('\n2. Input parameters:')
	for j in params:
		print('%s : %s' % (j,str(params[j])))
	print('\n3. Job statistics:')
	for j in job_stats:
		if j.find('time') == -1:
			print('%s : %s' % (j,str(job_stats[j])))
		else:
			print('%s : %.3f s' % (j,job_stats[j]))

