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
from display_ids import tsv_file,display_parameters

def main():
	start = time.time()
	job_stats = {'Software': 'SE','Software version': '2019.06.01.3'}
	job_stats['Start'] = str(datetime.datetime.now())
#	
#	load command line parameters and test for errors
#
	(params,ret,version) = load_params(sys.argv)
	if not ret:
		if version:
			print('%s: %s' % (job_stats['Software'],job_stats['Software version']))
		exit()
	display_parameters(params)
#
#	load spectra from files, using command line specified list
#
	print('\nLoading spectra')
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

	print('\nLoading kernels')
	kfs = params['kernel file'].split(',')
	kernel = []
	kmass = []
	spectrum_list = {}
	k = 0
	rcounts = 0
	(kernel,kmass,spectrum_list,k,rcounts) = load_kernel(kfs,spectra,params)

	job_stats['Kernels'] = k
	job_stats['KS-intersection'] = len(kernel)
	delta = time.time()-start
	job_stats['Load time kernel'] = delta
	print('   %.3f s' % (delta))
	job_stats['Redundant peptides'] = rcounts
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
	if len(spectra) > 0:
		job_stats['Search time (/1000)'] = 1000*(time.time()-start)/len(spectra)
	else:
		job_stats['Search time (/1000)'] = 0
	job_stats['End'] = str(datetime.datetime.now())
#
#	output the results
#
	print('\nCalulating statistics and storing results')

	if 'output file' in params:
		tsv_file(ids,scores,spectra,kernel,job_stats,params)
	else:
		tsv_file(ids,scores,spectra,kernel,job_stats,params)
	print('\n... done')

if __name__== "__main__":
	main()

