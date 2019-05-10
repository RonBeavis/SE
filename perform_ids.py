#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# performs peptide-spectrum matches, based on the pre-loaded spectra
# and kernel
#
#from __future__ import print_function
#from libcpp cimport bool as bool_t

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
	ids = {}
	scores = {}
	a = 0
	ires = float(_param['fragment mass tolerance'])
	score = 0
	b_score = 6
	if ires > 50:
		b_score = 7
	if ires > 100:
		b_score = 8
	best_score = 0
	best_intensity = 0.0
	ident = []
#
#	indicate progress to user
#
	print('.',end='')
	sys.stdout.flush()

#
#	iterate through spectra and perform PSM scoring
# 
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
		best_intensity = 0.0
		ident = []
#		s_set = set(s['sms'])
		s_set = dict(zip(s['sms'],s['ims']))
#
#		iterate through kernels on the list
#		and track scoring
#
		for k in ks:
			(score,intensity) = score_id(s_set,_k[k])
			f = 100.0*intensity/s['isum']
			if score > best_score:
				best_score = score
				ident = []
				best_intensity = f
				ident.append(k)
			elif score == best_score and f > best_intensity:				
				best_score = score
				ident = []
				best_intensity = f
				ident.append(k)
			elif score == best_score and f == best_intensity:
				best_intensity = f
				ident.append(k)
#
#		record PSM results
#
		ids[a] = ident
		scores[a] = (best_score,best_intensity)
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

def score_id(_s,_k):
	c = 0
	i = 0
	for k in _k:
		if k in _s:
			c += 1
			i += _s[k]
	return (c,i)

