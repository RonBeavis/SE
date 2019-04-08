#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
import sys

def perform_ids(_s,_k,_list,_param):
	ids = {}
	scores = {}
	a = 0
	res = _param['fragment mass tolerance']
	ires = 1.0*res
	score = 0
	best_score = 5
	okerns = []
	for k in _k:
		okerns.append([int(0.5+i/ires) for i in k['ms']])
	for s in _s:
		if a not in _list:
			a += 1
			continue
		ks = _list[a]
		best_score = 5
		ident = []
		sps = s['ms']
		tps = []
		for p in sps:
			val = int(0.5+p/ires)
			tps.append(val)
			tps.append(val-1)
			tps.append(val+1)
		s_set = set(tps)
		for k in ks:
			score = score_id(s_set,okerns[k],ires)
			if score > best_score:
				best_score = score
				ident = []
				ident.append(k)
			elif score == best_score:
				ident.append(k)
		ids[a] = ident
		scores[a] = best_score
		a += 1
		if a % 1000 == 0:
			print('.',end='')
			sys.stdout.flush()
	return (ids,scores)

def score_id(_s,_k,_r):
	score = 0
	for k in _k:
		if k in _s:
			score += 1
	return score

