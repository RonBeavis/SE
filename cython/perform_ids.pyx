#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# performs peptide-spectrum matches, based on the pre-loaded spectra
# and kernel
#
from __future__ import print_function
from libcpp cimport bool as bool_t

import sys

#
# runs through the list of spectra (_s) and matches each one
# to the best alternative present in the kernel list (_k)
# using a preprepared subset of the kernel(_list) and the
# specified job parameters (_param)
# the only externally called method
#

def perform_ids(_s,_k,_list,_param):
#
#	initialize local variables
#
	cdef dict ko = _param['kernel order'].copy()
	cdef dict ids = {}
	cdef dict scores = {}
	cdef long a = 0
	cdef float ires = float(_param['fragment mass tolerance'])
	cdef long score = 0
	cdef long b_score = 6
	cdef long best_score = 0
	cdef list ident = []
#
#	indicate progress to user
#
	print('.',end='')
	sys.stdout.flush()

#
#	iterate through spectra and perform PSM scoring
# 
	cdef dict s = {}
	for s in _s:
#
#		skip spectra with 0 kernels
#
		if a not in _list:
			a += 1
			continue
#
#		initialize spectrum variables
#
		ks = _list[a]
		best_score = b_score
		ident = []
		s_set = set(s['sms'])
#
#		iterate through kernels on the list
#		and track scoring
#
		for k in ks:
			score = score_id(s_set,_k[k])
			if score > best_score:
				best_score = score
				ident = []
				ident.append(k)
			elif score == best_score:
				ident.append(k)
#
#		record PSM results
#
		ids[a] = ident
		scores[a] = best_score
		a += 1
#
#		indicate progress to user
#
		if a % 10000 == 0:
			print('.',end='')
			sys.stdout.flush()
	return (ids,scores)

#
#	count the number of matches between a normalized kernel and spectrum pair
#

cdef long score_id(set _s,list _k):
	cdef long c = 0
	for k in _k:
		if k in _s:
			c += 1
	return c

