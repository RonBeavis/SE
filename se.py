#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import time
import sys
import datetime
from load_spectra import load_spectra
from load_kernel import read_kernel
from load_params import load_params
from perform_ids import perform_ids
from display_ids import simple_display

def main():
	start = time.time()
	job_stats = {}
	job_stats['Start'] = str(datetime.datetime.now())
	params = {'fragment mass tolerance': 400,
		'parent mass tolerance': 10,
#		'p mods':{'C':[57021,57021],'U':[57021,57021],'K':[0,6020],'R':[0,6020]},
		'p mods':{'C':[57021],'U':[57021]},
		'v mods':{'M':[15995]},
		'o mods':{'nt-ammonia':True,'nt-water':True}}
	print('Loading parameters')
	(params,ret) = load_params(sys.argv,params)
	if not ret:
		exit()
	print('\nLoading spectra')
	spectra = load_spectra(params['spectra file'])
	job_stats['S-dimension'] = len(spectra)
	print('\nLoading kernel')
	(kernel,spectrum_list,k) = read_kernel(params['kernel file'],spectra,params)
	job_stats['K-dimension'] = k
	job_stats['KS-intersection'] = len(kernel)

	job_stats['Load time'] = time.time()-start
	start = time.time()
	print('\nRunning ids')
	(ids,scores) = perform_ids(spectra,kernel,spectrum_list,params)
	job_stats['Search time'] = time.time()-start
	job_stats['Search time (/1000)'] = 1000*(time.time()-start)/len(spectra)
	job_stats['End'] = str(datetime.datetime.now())

	simple_display(ids,scores,spectra,kernel,job_stats,params)

if __name__== "__main__":
  main()

