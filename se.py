#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# runs a peptide identification mapping job, based on command line
# specified parameters
#

import time
import sys
import datetime
from load_spectra import load_spectra
from load_kernel import load_kernel
from load_params import load_params,load_defaults
from perform_ids import perform_ids
from display_ids import simple_display,tsv_file

def main():
	print('started ...\n')
	start = time.time()
	job_stats = {'Software': 'SE','Software version': '2019.04.10.1'}
	job_stats['Start'] = str(datetime.datetime.now())
#	
#	load command line parameters and test for errors
#
	print('Loading parameters')
	(params,ret) = load_params(sys.argv)
	if not ret:
		print('\n... exited')
		exit()
#
#	load spectra from files, using command line specified list
#
	print('Loading spectra')
	spectra = []
	sfs = params['spectra file'].split(',')
	for sf in sfs:
		spectra += load_spectra(sf,params)
	job_stats['Spectra'] = len(spectra)
#
#	load kernels from files, using command line specified list
#
	delta = time.time()-start
	job_stats['Load time spectra'] = delta
	print('   %.3f s' % (delta))
	start = time.time()

	print('\nLoading kernel')
	kfs = params['kernel file'].split(',')
	kernel = []
	kmass = []
	spectrum_list = {}
	k = 0
	qn = 0
	kernel_order = {}
	for kf in kfs:
		(kern,km,s_list,k1,qn,kernel_order) = load_kernel(kf,spectra,params,qn)
		if kernel_order is not None:
			params['kernel order'] = kernel_order
		kernel += kern
		kmass += km
		for s in s_list:
			if s in spectrum_list:
				spectrum_list[s] += s_list[s]
			else:
				spectrum_list[s] = s_list[s]
		k += k1
	job_stats['Kernels'] = k
	job_stats['KS-intersection'] = len(kernel)
	delta = time.time()-start
	job_stats['Load time kernel'] = delta
	print('   %.3f s' % (delta))
	start = time.time() 
#
#	generate identifications
#
	print('\nRunning ids')
	(ids,scores) = perform_ids(spectra,kmass,spectrum_list,params)
	delta = time.time()-start
	job_stats['Search time'] = time.time()-start
	print('   %.3f s' % (delta))
	start = time.time() 
	job_stats['Search time (/1000)'] = 1000*(time.time()-start)/len(spectra)
	job_stats['End'] = str(datetime.datetime.now())
#
#	output the results
#
	if 'output file' in params:
		tsv_file(ids,scores,spectra,kernel,job_stats,params)
	else:
		simple_display(ids,scores,spectra,kernel,job_stats,params)
	print('\n... done')

if __name__== "__main__":
	main()

