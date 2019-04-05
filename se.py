#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import json
import re
import time
import sys
import gzip
from load_spectra import load_spectra
from load_kernel import read_kernel

def perform_ids(_s,_k,_list,_param):
	ids = {}
	a = 0
	res = _param['res']
	for s in _s:
		ks = _list[a]
		best_score = 5
		ident = []
		for k in ks:
			tk = _k[k]
			score = score_id(s['ms'],tk['ms'],res)
			if score > best_score:
				best_score = score
				ident = []
				ident.append(k)
			elif score == best_score:
				ident.append(k)
		ids[a] = ident
		a += 1
	return ids

def score_id(_s,_k,res):
	score = 0
	s = 0
	k = 0
	ls = len(_s)
	lk = len(_k)
	while s < ls and k < lk:
		if abs(_s[s] - _k[k]) < res:
			score += 1
			s += 1
			k += 1
		else:
			if _s[s] > _k[k]:
				k += 1
			else:
				s += 1
	return score

params = {'res': 400,'pres': 20}

start = time.time()

lib = sys.argv[1]
sfile = sys.argv[2]

spectra = load_spectra(sfile)
print('S-dimension = %i' % (len(spectra)))

(kernel,spectrum_list,k) = read_kernel(lib,spectra,params)
print('K-dimension = %i' % (k))
print('KS-dimension = %i' % len(kernel))

ids = perform_ids(spectra,kernel,spectrum_list,params)
print(ids)

print('time = %.3f' % (time.time()-start))

