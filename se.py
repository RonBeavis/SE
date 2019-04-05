import json
import re
import time
import sys
import gzip
from load_spectra import load_spectra

isotopes = {'p': 1.007276,'H': 1.007825,'C': 12.0,'N': 14.003074,'O': 15.994915}

def read_kernel(_f,_s,_param):
	if _f.find('.gz') == len(_f) - 3:
		f = gzip.open(_f,'rt', encoding='utf-8')
	else:
		f = open(_f,'r', encoding='utf-8')
	t = 0
	pms = []
	qs = []
	spectrum_list = {}
	a = 0
	res = _param['pres']
	for s in _s:
		pms.append(s['pm'])
		spectrum_list[a] = []
		a += 1
	qn = 0
	for l in f:
		m = re.search('"pm": (.+?)\,',l)
		if m is None:
			continue
		pm = int(m.group(1))
		m = re.search('"seq": "(.+?)"',l)
		if m is None:
			continue
		seq = m.group(1)
		eseq = enumerate(seq)
		mod = [x for x, v in eseq if v == 'M']
		mod = [x for x, v in eseq if v == 'C']
		a = 0
		ok = False
		for s in pms:
			if abs(s-pm) < res:
				spectrum_list[a].append(qn)
				ok = True
			a += 1
		if ok:
			js = json.loads(l)
			ms = js['bs']+js['ys']
			ms.sort()
			js['ms'] = ms
			qs.append(js)
			qn += 1
		t += 1
	return (qs,spectrum_list,t)


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

