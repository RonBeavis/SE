#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# contains methods for loading and editing sequence kernels
# Reads kernels in either plain text or gzip'd text
#

# uncomment bisect to use the alternate spectrum indexing methods
# create_index_b and get_spectra_b
#import bisect 

import json
import re
import gzip
import sys
import itertools

#
# some handy lists of masses, in integer milliDaltons
#

modifications =	{ 'oxidation':15995,'carbamidomethyl':57021,'acetyl':42011,'dioxidation':31990,
		 'deamidation':984, 'carbamyl':43006,'phosphoryl':79966 }
isotopes =	{ 'p':1.007276,'C13':13.003355,'H1':1.007825,'H2':2.014102,'C12':12.0,'N14':14.003074,
		 'O16':15.994915,'S32':31.972072 }
amino_acids =	{ 'A':71037,'R':156101,'N':114043,'D':115027,'C':103009,'E':129043,'Q':128059,'G':57021,
		 'H':137059,'I':113084,'L':113084,'K':128095,'M':131040,'F':147068,'P':97053,'S':87032,
		 'T':101048,'U':150954,'W':186079,'Y':163063,'V':99068 }
#
# method to generate a list of kernels that are potential matches for the input list of spectra (_s)
# the kernels are read from a file (_f) and a set of parameters (_param) govern the way that
# kernels should be tested and modified, based on the experiment that was performed
# this method is the only one called externally
#

def load_kernel(_f,_s,_param,_qi):
	if _f.find('.gz') == len(_f) - 3:
		f = gzip.open(_f,'rt', encoding='utf-8')
	else:
		f = open(_f,'r', encoding='utf-8')
#
#	set kernel offset when loading multiple kernel files
#
	qn = _qi

#
# 	retrieve information from the _param dictionary and
#	create faster local variables
#
	depth = 3
	if 'ptm depth' in _param:
		depth = _param['ptm depth']
	if depth > 10:
		depth = 10
	res = _param['parent mass tolerance']
	ires = float(res)
	nt_ammonia = True
	if 'nt-ammonia' in _param['mods o']:
		nt_ammonia = _param['mods o']['nt-ammonia']
	nt_water = True
	if 'nt-water' in _param['mods o']:
		nt_water = _param['mods o']['nt-water']
	use_c13 = True
	if 'c13' in _param:
		use_c13 = _param.get('c13')
	acetyl = modifications['acetyl']
	p_mods = {}
	v_mods = {}
	if 'mods p' in _param:
		p_mods = _param['mods p']
	if 'mods v' in _param:
		v_mods = _param['mods v']
#
# 	create local variables for specific masses
#
	ammonia = normalize(isotopes['N14'] + 3*isotopes['H1'])
	water = normalize(isotopes['O16'] + 2*isotopes['H1'])
	c13 = normalize(isotopes['C13'] - isotopes['C12'])
#
# 	initialize some variable outside of the main iteration
#
	t = 0
	pms = []
	qs = []
	spectrum_list = {}
	p_pos = {}
	v_pos = {}
	p_total = 0
	(s_index,s_masses) = create_index(_s,ires)
#
# 	show activity to the user
#
	print('.',end='')
	sys.stdout.flush()
#	lines = f.readlines()
	
	for l in f:
#
# 		show activity to the user
#
		if t % 10000 == 0:
			print('.',end='')
			sys.stdout.flush()
		v_pos = {}
#
# 		retrieve kernel mass, protein coordinate and peptide sequence
#		this is much faster than doing a jsons.loads at this point
#
		m = re.search('"pm": (.+?)\,.+"beg": (\d+).+"seq": "(.+?)"',l)
		if m is None:
			continue
		pm = int(m.group(1))
		beg = int(m.group(2))
		seq = m.group(3)
#
# 		generate fixed modification information
#
		(lp_pos,lp_total,lp_len) = generate_lpstack(p_mods,seq)
#
# 		generate variable modification information
#
		v_pos = generate_vd(v_mods,seq)
		v_stack = generate_vstack(v_mods,v_pos,depth)
		ok = False
		delta = 0
		b_mods = []
		y_mods = []
#
# 		check for special case peptide N-terminal cyclization at Q, C or E
#
		n_term = seq[:1]
		ammonia_loss = False
		if nt_ammonia and (n_term == 'Q' or n_term == 'C'):
			ammonia_loss = True
		water_loss = False
		if nt_water and n_term == 'E':
			water_loss = True
		for vp in v_stack:
			lp = 0
			vs_pos = vp[0]
			vs_total= vp[1]
			for lp in range(lp_len):
				p_pos = lp_pos[lp]
				p_total = lp_total[lp]+vs_total
				tmass = pm+p_total
				pms = get_spectra(s_index,tmass,ires)
				if use_c13 and tmass > 1500:
					pms += get_spectra(s_index,tmass+c13,ires)
				appended = False
				js_base = None
				for s in pms:
					delta = s_masses[s]-tmass
					if abs(delta) < res or abs(delta-c13) < res:
						if s not in spectrum_list:
							spectrum_list[s] = []
						spectrum_list[s].append(qn)
						if not appended:
							if js_base is None:
								js_base = load_base(l,p_pos,p_mods,lp)
							qs.append(load_json(js_base,b_mods,y_mods,lp,vs_pos,v_mods))
							appended = True
				if appended:
					qn += 1
				appended = False
				if beg < 3:
					tmass = pm+p_total+acetyl
					pms = get_spectra(s_index,tmass,ires)
					if use_c13 and tmass > 1500:
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
								if js_base is None:
									js_base = load_base(l,p_pos,p_mods,lp)
								qs.append(load_json(js_base,b_mods,y_mods,lp,vs_pos,v_mods))
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
					if use_c13 and tmass > 1500:
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
								if js_base is None:
									js_base = load_base(l,p_pos,p_mods,lp)
								qs.append(load_json(js_base,b_mods,y_mods,lp,vs_pos,v_mods))
								appended = True

				if appended:
					qn += 1

		t += 1
	return (qs,spectrum_list,t,qn)
#
# method to convert floating point masses in Daltons to integer masses in milliDaltons
#

def normalize(_v):
	return int(round(1000*_v,0))

#
# method to create an array with each possible variable modification state
# in a single element
#

def generate_vstack(_mods,_pos,_depth = 3):
	v_stack = []
	vs_pos = {}
	master_list = []
#
#	create an empty vs_pos dict (unmodified)
#	and generate a list of the possible modifications
#	where each element is a tuple of (residue,position)
#
	for v in _mods:
		vs_pos[v] = []
		if v not in _pos:
			continue
		else:
			for p in _pos[v]:
				master_list.append((v,p))
	vs_item = [vs_pos,0]
	v_stack.append(vs_item)
	d = 1
#
#	iterate to the specified depth of modification
#
	while d <= _depth:
#
#		generate a list of all possible combinations of "d" modifications
#
		m_list = list(itertools.combinations(master_list,d))
#
#		iterate through the modification combinations to create
#		the structures used to update the b and y ion lists
#		and supply the peptide mass change caused by the modifications
#
		for ml in m_list:
			dm = 0
			vs_pos = {}
			for v in _mods:
				if v not in _pos:
					continue
				vs_pos[v] = [x[1] for x in ml if x[0] == v]
				dm += _mods[v][0]*len(vs_pos[v])
			v_stack.append([vs_pos,dm])
		d += 1
	return v_stack

#
# method to locate possible variable modification sites in a sequence
#

def generate_vd(_mods,_seq):
	keep = False
	v_pos = {}
	for v in _mods:
		v_pos[v] = []
		if v == '[' and _mods[v] != 0:
			v_pos[v] = [0]
			keep = True
		elif v == ']' and _mods[v] != 0:
			v_pos[v] = [len(seq)-1]
			keep = True
		if _seq.find(v) != -1 and _mods[v] != 0:
			v_pos[v] = [x for x, y in enumerate(_seq) if y == v]
			keep = True
	if not keep:
		v_pos = {}
	return v_pos

#
# method to locate fixed modification sites in a sequence and create an array
# with this information for each set of modification states to be tested
#

def generate_lpstack(_mods,_seq):
	lp_len = 0
	lp_pos = []
	lp_total = []
	for p in _mods:
		lp_len = len(_mods[p])
		break
	lp = 0
	while lp < lp_len:
		p_pos = {}
		p_total = 0
		keep = False
		for p in _mods:
			p_pos[p] = []
			if p == '[' and _mods[p][lp] != 0:
				p_pos[p] = [0]
				p_total += _mods[p][lp]
			elif p == ']' and _mods[p][lp] != 0:
				p_pos[p] = [len(seq)-1]
				p_total += _mods[p][lp]
			elif _seq.find(p) != -1 and _mods[p][lp] != 0:
				p_pos[p] = [x for x, y in enumerate(_seq) if y == p]
				p_total += len(p_pos[p])*_mods[p][lp]
				keep = True
		if not keep:
			p_pos = {}
		lp_pos.append(p_pos)
		lp_total.append(p_total)
		lp += 1
	return (lp_pos,lp_total,lp_len)

#
# generate a JSON object from the kernel entry line (_l) and modify it
# using the information generated from the allowed modification lists
#

def load_base(_l,_p_pos,_p_mods,_lp):
	jin = json.loads(_l)
	jin['mods'] = []
	if len(_p_pos) > 0:
		jin = update_ions(jin,_p_mods,_p_pos,_lp)
	return jin	

def load_json(_js,_b_mods,_y_mods,_lp,_vs_pos,_v_mods):
	jin = _js.copy()
	if len(_vs_pos) > 0:
		jin = update_ions(jin,_v_mods,_vs_pos,0)
	if len(_b_mods):
		jin = update_bions(jin,_b_mods)
	if len(_y_mods):
		jin = update_yions(jin,_y_mods)
	ms = jin['bs']+jin['ys']
	ms.sort()
	jin['ms'] = ms
	if jin.pop('bs',None) is None:
		print('error removing bs')
	if jin.pop('ys',None) is None:
		print('error removing ys')
#	v = []
#	for j in jin:
#		v.append(jin[j])
	return jin	

#
# method to update ion series based on a set of modifications specified by:
# _pos - the array of locations (_pos), 
# _mod - the array of modifications, &
# _lp - the position in _pos
#
def update_ions(_js,_mods,_pos,_lp):
	jin = _js.copy()
	t = len(jin['bs'])
	mods = {}
	tmod = 0.0
	beg = jin['beg']
	for m in _mods:
		if len(_pos[m]) == 0:
			continue
		a = 0
		delta = 0
		pmod = _mods[m][_lp]
		while a < t:
			if a in _pos[m]:
				mods[beg+a] = pmod
				tmod += pmod
				delta += pmod
			if delta == 0:
				a += 1
				continue
			jin['bs'][a] += delta
			a += 1
		if a in _pos[m]:
			mods[beg+a] = pmod
			tmod += pmod
		a = 0
		delta = 0
		while a < t:
			if t-a in _pos[m]:
				delta += pmod
			if delta == 0:
				a += 1
				continue
			jin['ys'][a] += delta
			a += 1
	jin['mods'].append(mods)
	jin['pm'] += tmod
	return jin

#
# two ion series update functions that deal with modifications to either
# end of the sequence
# update_bions for N-terminal modifications
# update_yions for C-terminal modification
# Note: these could be handled by update_ions, but they are
#       broken out as separate methods for clarity and easier
#       debugging
#
def update_bions(_js,_bmods):
	js = _js.copy()
	t = len(js['bs'])
	mods = {}
	beg = js['beg']
	for b in _bmods:
		js['pm'] += b
		mods[beg] = b
		a = 0
		while a < t:
			js['bs'][a] += b
			a += 1
	js['mods'].append(mods)
	return js

def update_yions(_js,_ymods):
	js = _js.copy()
	t = len(js['ys'])
	mods = {}
	end = js['end']
	for b in _ymods:
		js['pm'] += b
		mods[end] = b
		a = 0
		while a < t:
			js['ys'][a] += b
			a += 1
	js['mods'].append(mods)
	return js

#
# Two methods for the spectrum mass indexing system
# that uses a dictionary and binning to perform the matches
#

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

def get_spectra(_index,_mass,_r):
	iv = int(0.5+_mass/_r)
	pms = []
	if iv in _index:
		pms += _index.get(iv)
	return pms

#
# Two methods for an alternate spectrum mass indexing system
# that uses "bisect" to find matches - more exact but slower
#
#def create_index_b(_sp,_r):
#	index = []
#	masses = []
#	a = 0
#	pm = 0
#	m = 0
#	for s in _sp:
#		m = s['pm']
#		masses.append(m)
#		index.append((m,a))
#		a += 1
#	index.sort()
#	return (index,masses)

#def get_spectra_b(_index,_mass,_r):
#	li = len(_index)
#	pms = []
#	last = _mass + _r
#	s = bisect.bisect_left(_index,(_mass-_r,))
#	while s < li and _index[s][0] <= last:
#		pms.append(_index[s][1])
#		s += 1
#	return pms

