#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import json
import re
import gzip
import sys
import bisect
import itertools

modifications =	{ 'oxidation':15995,'carbamidomethyl':57021,'acetyl':42011,'dioxidation':31990,
		 'deamidation':984, 'carbamyl':43006,'phosphoryl':79966 }
isotopes =	{ 'p':1.007276,'C13':13.003355,'H1':1.007825,'H2':2.014102,'C12':12.0,'N14':14.003074,
		 'O16':15.994915,'S32':31.972072 }
amino_acids =	{ 'A':71037,'R':156101,'N':114043,'D':115027,'C':103009,'E':129043,'Q':128059,'G':57021,
		 'H':137059,'I':113084,'L':113084,'K':128095,'M':131040,'F':147068,'P':97053,'S':87032,
		 'T':101048,'U':150954,'W':186079,'Y':163063,'V':99068 }

def normalize(_v):
	return int(round(1000*_v,0))

def read_kernel(_f,_s,_param):
	if _f.find('.gz') == len(_f) - 3:
		f = gzip.open(_f,'rt', encoding='utf-8')
	else:
		f = open(_f,'r', encoding='utf-8')
	t = 0
	pms = []
	qs = []
	spectrum_list = {}
	res = _param['parent mass tolerance']
	ires = float(res)
	qn = 0
	ammonia = normalize(isotopes['N14'] + 3*isotopes['H1'])
	water = normalize(isotopes['O16'] + 2*isotopes['H1'])
	c13 = normalize(isotopes['C13'] - isotopes['C12'])
	acetyl = modifications['acetyl']
	p_mods = {}
	v_mods = {}
	if 'p mods' in _param:
		p_mods = _param['p mods']
	if 'v mods' in _param:
		v_mods = _param['v mods']
	p_pos = {}
	v_pos = {}
	p_total = 0
	nt_ammonia = True
	if 'nt-ammonia' in _param['o mods']:
		nt_ammonia = _param['o mods']['nt-ammonia']
	nt_water = True
	if 'nt-water' in _param['o mods']:
		nt_water = _param['o mods']['nt-water']
	(s_index,s_masses) = create_index(_s,ires)
	for l in f:
		if t % 1000 == 0:
			print('.',end='')
			sys.stdout.flush()
		v_pos = {}
		m = re.search('"pm": (.+?)\,.+"beg": (\d+).+"seq": "(.+?)"',l)
		if m is None:
			continue
		pm = int(m.group(1))
		beg = int(m.group(2))
		seq = m.group(3)
		n_term = seq[:1]
		lp_pos = []
		lp_total = []
		lp_len = 0
		for p in p_mods:
			lp_len = len(p_mods[p])
			break
		lp = 0
		while lp < lp_len:
			p_pos = {}
			p_total = 0
			keep = False
			for p in p_mods:
				p_pos[p] = []
				if p == '[' and p_mods[p][lp] != 0:
					p_pos[p] = [0]
					p_total += p_mods[p][lp]
				elif p == ']' and p_mods[p][lp] != 0:
					p_pos[p] = [len(seq)-1]
					p_total += p_mods[p][lp]
				elif seq.find(p) != -1 and p_mods[p][lp] != 0:
					p_pos[p] = [x for x, y in enumerate(seq) if y == p]
					p_total += len(p_pos[p])*p_mods[p][lp]
					keep = True
			if not keep:
				p_pos = {}
			lp_pos.append(p_pos)
			lp_total.append(p_total)
			lp += 1
		keep = False
		for v in v_mods:
			v_pos[v] = []
			if v == '[' and v_mods[v] != 0:
				v_pos[v] = [0]
			elif v == ']' and v_mods[v] != 0:
				v_pos[v] = [len(seq)-1]
			if seq.find(v) != -1 and v_mods[v] != 0:
				v_pos[v] = [x for x, y in enumerate(seq) if y == v]
				keep = True
		if not keep:
			v_pos = {}
		for v in v_mods:
			if v in v_pos:
				tvs = list(itertools.combinations(v_pos[v],2))
		ok = False
		delta = 0
		b_mods = []
		y_mods = []
		ammonia_loss = False
		if nt_ammonia and (n_term == 'Q' or n_term == 'C'):
			ammonia_loss = True
		water_loss = False
		if nt_water and n_term == 'E':
			water_loss = True
		v_stack = []
		vs_pos = {}
		for v in v_mods:
			vs_pos[v] = []
		vs_item = [vs_pos,0]
		v_stack.append(vs_item)
		for v in v_mods:
			if v not in v_pos:
				continue
			l_pos = list(itertools.combinations(v_pos[v],1))
			for lt in l_pos:
				vs_pos = {}
				vs_pos[v] = list(lt)
				v_stack.append([vs_pos,len(vs_pos)*v_mods[v][0]])
			l_pos = list(itertools.combinations(v_pos[v],2))
			for lt in l_pos:
				vs_pos = {}
				vs_pos[v] = list(lt)
				v_stack.append([vs_pos,len(vs_pos)*v_mods[v][0]])
		for vp in v_stack:
			lp = 0
			vs_pos = vp[0]
			vs_total= vp[1]
			for lp in range(lp_len):
				p_pos = lp_pos[lp]
				p_total = lp_total[lp]+vs_total
				tmass = pm+p_total
				pms = get_spectra(s_index,tmass,ires)
				if tmass > 1500:
					pms += get_spectra(s_index,tmass+c13,ires)
				appended = False
				for s in pms:
					delta = s_masses[s]-tmass
					if abs(delta) < res or abs(delta-c13) < res:
						if s not in spectrum_list:
							spectrum_list[s] = []
						spectrum_list[s].append(qn)
						if not appended:
							qs.append(load_json(l,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods))
							appended = True
				if appended:
					qn += 1
				appended = False
				if beg < 3:
					tmass = pm+p_total+acetyl
					pms = get_spectra(s_index,tmass,ires)
					if tmass > 1500:
						pms += get_spectra(s_index,tmass+c13,ires)
					for s in pms:
						delta = s_masses[s]-pm-p_total
						if abs(delta-acetyl) < res or abs(delta-acetyl-c13) < res:
							if s not in spectrum_list:
								spectrum_list[s] = []
							spectrum_list[s].append(qn)
							if acetyl not in b_mods:
								b_mods.append(acetyl)
							if not appended:
								qs.append(load_json(l,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods))
								appended = True

				if appended:
					qn += 1
				appended = False
				if ammonia_loss or water_loss:
					dvalue = ammonia
					if water_loss:
						d_value = water
					tmass = pm+p_total-dvalue
					pms = get_spectra(s_index,tmass,ires)
					if pm+p_total-dvalue > 1500:
						pms += get_spectra(s_index,tmass+c13,ires)
					for s in pms:
						delta = s_masses[s]-pm-p_total
						if abs(delta+dvalue) < res or abs(delta+dvalue-c13) < res:
							if s not in spectrum_list:
								spectrum_list[s] = []
							spectrum_list[s].append(qn)
							if dvalue not in b_mods:
								b_mods.append(-1*dvalue)
							if not appended:
								qs.append(load_json(l,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods))
								appended = True

				if appended:
					qn += 1

		t += 1
	return (qs,spectrum_list,t)

def load_json(l,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods):
	js = json.loads(l)
	js['mods'] = []
	if len(p_pos) > 0:
		js = update_ions(js,p_mods,p_pos,lp)
	if len(vs_pos) > 0:
		js = update_ions(js,v_mods,vs_pos,0)
	if len(b_mods):
		js = update_bions(js,b_mods)
	if len(y_mods):
		js = update_yions(js,y_mods)
	ms = js['bs']+js['ys']
	ms.sort()
	js['ms'] = ms
	return js	

def get_spectra(_index,_mass,_r):
	iv = int(0.5+_mass/_r)
	pms = []
	if iv in _index:
		pms += _index.get(iv)
	return pms

def get_spectra_b(_index,_mass,_r):
	li = len(_index)
	pms = []
	last = _mass + _r
	s = bisect.bisect_left(_index,(_mass-_r,))
	while s < li and _index[s][0] <= last:
		pms.append(_index[s][1])
		s += 1
	return pms
 
def update_ions(_js,_mods,_pos,_lp):
	js = _js
	t = len(js['bs'])
	mods = {}
	tmod = 0.0
	for m in _mods:
		if len(_pos[m]) == 0:
			continue
		a = 0
		delta = 0
		pmod = _mods[m][_lp]
		while a < t:
			if a in _pos[m]:
				mods[str(js['beg']+a)] = pmod
				tmod += pmod
				delta += pmod
			if delta == 0:
				a += 1
				continue
			js['bs'][a] += delta
			a += 1
		if a in _pos[m]:
			mods[str(js['beg']+a)] = pmod
			tmod += pmod
		a = 0
		delta = 0
		while a < t:
			if t-a in _pos[m]:
				delta += pmod
			if delta == 0:
				a += 1
				continue
			js['ys'][a] += delta
			a += 1
	js['mods'].append(mods)
	js['pm'] += tmod
	return js

def update_bions(_js,_bmods):
	js = _js
	t = len(js['bs'])
	mods = {}
	for b in _bmods:
		js['pm'] += b
		mods[str(js['beg'])] = b
		a = 0
		while a < t:
			js['bs'][a] += b
			a += 1
	js['mods'].append(mods)
	return js

def update_yions(_js,_ymods):
	js = _js
	t = len(js['ys'])
	mods = {}
	for b in _ymods:
		js['pm'] += b
		mods[str(js['end'])] = b
		a = 0
		while a < t:
			js['ys'][a] += b
			a += 1
	js['mods'].append(mods)
	return js

def create_index(_sp,_r):
	index = {}
	masses = []
	a = 0
	pm = 0
	m = 0
	for s in _sp:
		m = s['pm']
		masses.append(m)
		pm = int(0.5 + m/_r)
		if pm in index:
			index[pm].append(a)
		else:
			index[pm] = [a]
		if pm-1 in index:
			index[pm-1].append(a)
		else:
			index[pm-1] = [a]
		if pm+1 in index:
			index[pm+1].append(a)
		else:
			index[pm+1] = [a]
		a += 1
	return (index,masses)

def create_index_b(_sp,_r):
	index = []
	masses = []
	a = 0
	pm = 0
	m = 0
	for s in _sp:
		m = s['pm']
		masses.append(m)
		index.append((m,a))
		a += 1
	index.sort()
	return (index,masses)

