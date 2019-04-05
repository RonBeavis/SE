#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import json
import re
import gzip

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
	res = _param['parent mass tolerance']
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


