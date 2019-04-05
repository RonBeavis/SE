#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

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

