#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import json
import re
import gzip

modifications =	{ 'oxidation':15995,'carbamidomethyl':57021,'acetyl':42011,'dioxidation':31990,
		 'deamidation':984, 'carbamyl':43006,'phosphoryl':79966 }
isotopes =	{ 'p':1.007276,'C13':13.003355,'H1':1.007825,'H2':2.014102,'C12':12.0,'N14':14.003074,
		 'O16':15.994915,'S32':31.972072 }
amino_acids =	{ 'A':71037,'R':156101,'N':114043,'D':115027,'C':103009,'E':129043,'Q':128059,'G':57021,
		 'H':137059,'I':113084,'L':113084,'K':128095,'M':131040,'F':147068,'P':97053,'S':87032,
		 'T':101048,'U':150954,'W':186079,'Y':163063,'V':99068 }

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
	ammonia = int(round(1000*(isotopes['N14'] + 3*isotopes['H1']),0))
	water = int(round(1000*(isotopes['O16'] + 2*isotopes['H1']),0))
	c13 = int(round(1000*(isotopes['C13'] - isotopes['C12']),0))
	acetyl = modifications['acetyl']
	pmods = {}
	vmods = {}
	if 'p mods' in _param:
		pmods = _param['p mods']
	if 'v mods' in _param:
		vmods = _param['v mods']
	p_pos = {}
	v_pos = {}
	p_total = 0
	for l in f:
		p_pos = {}
		v_pos = {}
		p_total = 0
		m = re.search('"pm": (.+?)\,.+"beg": (\d+).+"seq": "(.+?)"',l)
		if m is None:
			continue
		pm = int(m.group(1))
		beg = int(m.group(2))
		seq = m.group(3)
		eseq = enumerate(seq)
		keep = False
		for p in pmods:
			p_pos[p] = []
			if seq.find(p) != -1:
				p_pos[p] = [x for x, y in eseq if y == p]
				p_total += len(p_pos[p])*pmods[p][0]
				keep = True
		if not keep:
			p_pos = {}
		keep = False
		for v in vmods:
			v_pos[v] = []
			if seq.find(v) != -1:
				v_pos[v] = [x for x, y in eseq if y == v]
				keep = True
		if not keep:
			v_pos = {}
		a = 0
		ok = False
		delta = 0
		b_mods = []
		y_mods = []
		ammonia_loss = False
		if seq.find('Q') == 0 or seq.find('C') == 0:
			ammonia_loss = True
		for s in pms:
			delta = s-pm-p_total
			if abs(delta) < res:
				spectrum_list[a].append(qn)
				ok = True
			elif beg < 3 and abs(delta-acetyl) < res:
				spectrum_list[a].append(qn)
				if acetyl not in b_mods:
					b_mods.append(acetyl)
				ok = True
			elif ammonia_loss and abs(delta+ammonia) < res:
				spectrum_list[a].append(qn)
				if -1*ammonia not in b_mods:
					b_mods.append(-1*ammonia)
				ok = True
			if pm > 1500000:
				delta -= c13
				if abs(delta) < res:
					spectrum_list[a].append(qn)
					ok = True
				elif beg < 3 and abs(delta-acetyl) < res:
					spectrum_list[a].append(qn)
					if acetyl not in b_mods:
						b_mods.append(acetyl)
					ok = True
				elif ammonia_loss and abs(delta+ammonia) < res:
					spectrum_list[a].append(qn)
					if -1*ammonia not in b_mods:
						b_mods.append(-1*ammonia)
					ok = True
			a += 1
		if ok:
			js = json.loads(l)
			js['mods'] = []
			if len(p_pos) > 0:
				js = update_ions(js,pmods,p_pos)
			if len(b_mods):
				js = update_bions(js,b_mods)
			if len(y_mods):
				js = update_yions(js,y_mods)
			ms = js['bs']+js['ys']
			ms.sort()
			js['ms'] = ms
			qs.append(js)
			qn += 1
		t += 1
	return (qs,spectrum_list,t)

def update_ions(_js,_mods,_pos):
	js = _js
	t = len(js['bs'])
	mods = {}
	for m in _mods:
		a = 0
		delta = int(0)
		pmod = _mods[m][0]
		while a < t:
			if a in _pos[m]:
				mods[str(js['beg']+a)] = pmod
				js['mods'].append(mods)
				delta += pmod
			js['bs'][a] += delta
			a += 1
		a = 0
		delta = 0
		while a < t:
			if t-a in _pos[m]:
				delta += pmod
			js['ys'][a] += delta
			a += 1
	return js

def update_bions(_js,_bmods):
	js = _js
	t = len(js['bs'])
	mods = {}
	for b in _bmods:
		mods[str(js['beg'])] = b
		js['mods'].append(mods)
		a = 0
		while a < t:
			js['bs'][a] += b
			a += 1
	return js

def update_yions(_js,_ymods):
	js = _js
	t = len(js['ys'])
	mods = {}
	for b in _ymods:
		mods[str(js['end'])] = b
		js['mods'].append(mods)
		a = 0
		while a < t:
			js['ys'][a] += b
			a += 1
	return js

