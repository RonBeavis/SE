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
#from __future__ import print_function
#from libcpp cimport bool as bool_t

import ujson
#import json
import re
import gzip
import sys
import itertools
import copy
#
# method to import a list of isotopic masses necessary to do some of the calculations.
# if 'isotopes.txt' is not available, a warning is thrown and a default list is used.
#
def load_isotopes():
	try:
		f = open('isotopes.txt','r')
	except:
		print('Warning: isotopes.txt is not available so default values used')
		return { 'p':1.007276,'1H':1.007825,'2H':2.014102,'12C':12.0,'13C':13.003355,'14N':14.003074,
		 	'16O':15.994915,'32S':31.972072 }
	iso = {}
	for l in f:
		l = l.strip()
		vs = l.split('\t')
		if len(vs) < 2:
			continue
		iso[vs[0]] = float(vs[1])
	f.close()
	return iso
#
# method to generate a list of kernels that are potential matches for the input list of spectra (_s)
# the kernels are read from a file (_f) and a set of parameters (_param) govern the way that
# kernels should be tested and modified, based on the experiment that was performed
# this method is the only one called externally
#

def load_kernel(_f,_s,_param,_qi):
	isotopes = load_isotopes()
	if _f.find('.gz') == len(_f) - 3:
		f = gzip.open(_f,'rt', encoding='utf-8')
	else:
		f = open(_f,'r', encoding='utf-8')
#
#	set kernel offset when loading multiple kernel files
#
	qn = _qi
	sms_list = []
	for s in _s:
		sms_list.append(set(s['sms']))

#
# 	retrieve information from the _param dictionary and
#	create faster local variables
#
	depth = 3
	if 'ptm depth' in _param:
		depth = _param.get('ptm depth')
	if depth > 10:
		depth = 10
	c_limit = 5
	res = 50
	max_ppm = _param.get('parent mass tolerance')
	ires = float(res)
	fres = float(_param.get('fragment mass tolerance'))
	nt_ammonia = True
	if 'nt-ammonia' in _param['mods o']:
		nt_ammonia = _param['mods o']['nt-ammonia']
	nt_water = True
	if 'nt-water' in _param['mods o']:
		nt_water = _param['mods o']['nt-water']
	use_c13 = True
	if 'c13' in _param:
		use_c13 = _param.get('c13')
	acetyl = 42011
	p_mods = {}
	v_mods = {}
	if 'mods p' in _param:
		p_mods = _param.get('mods p')
	if 'mods v' in _param:
		v_mods = _param.get('mods v')
#
# 	create local variables for specific masses
#
	ammonia = normalize(isotopes.get('14N') + 3*isotopes.get('1H'))
	water = normalize(isotopes.get('16O') + 2*isotopes.get('1H'))
	c13 = normalize(isotopes.get('13C') - isotopes.get('12C'))
#
# 	initialize some variable outside of the main iteration
#
	t = 0
	pms = []
	qs = []
	qm = []
	spectrum_list = {}
	p_pos = {}
	v_pos = {}
	p_total = 0
	(s_index,s_masses) = create_index(_s,ires)
	kernel_order = None
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
		js_master = ujson.loads(l)
#		m = re.search('"pm": (.+?)\,.+"beg": (\d+).+"seq": "(.+?)"',l)
#		if m is None:
#			continue
		if 'pm' not in js_master:
			continue
		if kernel_order is None:
			kernel_order = {}
			x = 0
			js_master['mods'] = None
			for j in js_master:
				if j == 'bs' or j == 'ys':
					continue
				kernel_order[j] = x
				x += 1
			js_master.pop('mods')
		pm = js_master['pm']
		beg = js_master['beg']
		seq = js_master['seq']
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
		isIaa = []
		for v in p_mods:
			if v == 'C':
				for p in p_mods[v]:
					if p == 57021:
						isIaa.append(True)
					else:
						isIaa.append(False)
		n_term = seq[:1]
		q_ammonia_loss = False
		c_ammonia_loss = False
		if nt_ammonia and n_term == 'Q':
			q_ammonia_loss = True
		if nt_ammonia and n_term == 'C':
			c_ammonia_loss = True
		water_loss = False
		if nt_water and n_term == 'E':
			water_loss = True
#
#		Make copies of the arrays in js_master
#
		js_bs = list(js_master['bs'])
		js_ys = list(js_master['ys'])
		js_pm = js_master['pm']
#		for vp in v_stack:
#			lp = 0
#			vs_pos = vp[0]
#			vs_total= vp[1]
#			for lp in range(lp_len):
		for lp in range(lp_len):
			for vp in v_stack:
				vs_pos = vp[0]
				vs_total= vp[1]
				bIaa = False
				if 'C' in p_mods:
					bIaa = isIaa[lp]
				if c_ammonia_loss and not bIaa:
					c_ammonia_loss = False
				b_mods = []
				y_mods = []
				p_pos = lp_pos[lp]
				p_total = lp_total[lp]+vs_total
				tmass = pm+p_total
				slots = []
				ok_c13 = False
				if use_c13 and tmass > 1500000:
					slots = [c13]
					pms = get_spectra(s_index,tmass,ires,True,slots)
					ok_c13 = True
				else:
					pms = get_spectra(s_index,tmass,ires,False,slots)
				appended = False
				jv = None
				jc = None
				for s in pms:
					delta = s_masses[s]-tmass
					if(delta > 900):
						ppm = abs(1e6*(delta-c13)/tmass)
					else:
						ppm = abs(1e6*delta/tmass)
					if ppm < max_ppm:
#					if abs(delta) < res or (ok_c13 and abs(delta-c13) < res):
						if jv is None:
							(jv,jm) = load_json(js_master,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods,fres)
#
#							replace arrays in js_master
#
							js_master['bs'] = list(js_bs)
							js_master['ys'] = list(js_ys)
							js_master['pm'] = js_pm
							js_master['mods'] = []
						sms = sms_list[s]
						c = 0
						for k in jm:
							if k in sms:
								c += 1
						if c >= c_limit:
							if not appended:
								qs.append(jv)
								qm.append(jm)
								appended = True
							if s not in spectrum_list:
								spectrum_list[s] = []
							spectrum_list[s].append(qn)
				if appended:
					qn += 1
				appended = False
				if '[' in p_mods:
					if p_mods['['][lp] != 0:
						continue
				if beg < 4:
					b_mods = []
					y_mods = []
					tmass = pm+p_total+acetyl
					slots = []
					ok_c13 = False
					if use_c13 and tmass > 1500000:
						slots = [c13]
						pms = get_spectra(s_index,tmass,ires,True,slots)
						ok_c13 = True
					else:
						pms = get_spectra(s_index,tmass,ires,False,slots)
					jv = None
					jc = None
					for s in pms:
						delta = s_masses[s]-pm-p_total-acetyl
						if(delta > 900):
							ppm = abs(1e6*(delta-c13)/tmass)
						else:
							ppm = abs(1e6*delta/tmass)
						if ppm < max_ppm:
#						if abs(delta-acetyl) < res or (ok_c13 and abs(delta-acetyl-c13) < res):
							if acetyl not in b_mods:
								b_mods.append(acetyl)
							if jv is None:
								(jv,jm) = load_json(js_master,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods,fres)
#
#								replace arrays in js_master
#
								js_master['bs'] = list(js_bs)
								js_master['ys'] = list(js_ys)
								js_master['pm'] = js_pm
								js_master['mods'] = []
							sms = sms_list[s]
							c = 0
							for k in jm:
								if k in sms:
									c += 1
							if c >= c_limit:
								if not appended:
									qs.append(jv)
									qm.append(jm)
									appended = True
								if s not in spectrum_list:
									spectrum_list[s] = []
								spectrum_list[s].append(qn)

				if appended:
					qn += 1
				appended = False
				if q_ammonia_loss or c_ammonia_loss or water_loss:
					b_mods = []
					y_mods = []
					dvalue = ammonia
					if water_loss:
						dvalue = water
					tmass = pm+p_total-dvalue
					slots = []
					ok_c13 = False
					if use_c13 and tmass > 1500000:
						slots = [c13]
						pms = get_spectra(s_index,tmass,ires,True,slots)
						ok_c13= True
					else:
						pms = get_spectra(s_index,tmass,ires,False,slots)
					jv = None
					jc = None
					for s in pms:
						delta = s_masses[s]-pm-p_total+dvalue
						if(delta > 900):
							ppm = abs(1e6*(delta-c13)/tmass)
						else:
							ppm = abs(1e6*delta/tmass)
						if ppm < max_ppm:
#						if abs(delta+dvalue) < res or (ok_c13 and abs(delta+dvalue-c13) < res):
							if dvalue not in b_mods:
								b_mods.append(-1*dvalue)
							if jv is None:
								(jv,jm) = load_json(js_master,p_pos,p_mods,b_mods,y_mods,lp,vs_pos,v_mods,fres)
#
#								replace arrays in js_master
#
								js_master['bs'] = list(js_bs)
								js_master['ys'] = list(js_ys)
								js_master['pm'] = js_pm
								js_master['mods'] = []
							sms = sms_list[s]
							c = 0
							for k in jm:
								if k in sms:
									c += 1
							if c >= c_limit:
								if not appended:
									qs.append(jv)
									qm.append(jm)
									appended = True
								if s not in spectrum_list:
									spectrum_list[s] = []
								spectrum_list[s].append(qn)

				if appended:
					qn += 1

		t += 1
	return (qs,qm,spectrum_list,t,qn,kernel_order)
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
				dm += _mods[v][0]*len(vs_pos.get(v))
			v_stack.append([vs_pos,dm])
		d += 1
	return v_stack

#
# method to locate possible variable modification sites in a sequence
#

def generate_vd(_mods,_seq):
	keep = False
	v_pos = {}
	ls = len(_seq)
	for v in _mods:
		v_pos[v] = []
		if v == '[' and _mods.get(v) != 0:
			v_pos[v] = [0]
			keep = True
		elif v == ']' and _mods.get(v) != 0:
			v_pos[v] = [ls-1]
			keep = True
		if _seq.find(v) != -1 and _mods.get(v) != 0:
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
	ls = len(_seq)
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
			pms = _mods[p]
			if p == '[' and pms[lp] != 0:
				p_pos[p] = [0]
				p_total += pms[lp]
				keep = True
			elif p == ']' and pms[lp] != 0:
				p_pos[p] = [ls-1]
				p_total += pms[lp]
				keep = True
			elif _seq.find(p) != -1 and pms[lp] != 0:
				p_pos[p] = [x for x, y in enumerate(_seq) if y == p]
				p_total += len(p_pos[p])*pms[lp]
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

def load_json(_l,_p_pos,_p_mods,_b_mods,_y_mods,_lp,_vs_pos,_v_mods,_fres):
	jin = _l
	jin['mods'] = []
	if _p_pos:
		jin = update_ions(jin,_p_mods,_p_pos,_lp)
	if _vs_pos:
		jin = update_ions(jin,_v_mods,_vs_pos,0)
	if _b_mods:
		jin = update_bions(jin,_b_mods)
	if _y_mods:
		jin = update_yions(jin,_y_mods)
#	dm = jin['pm'] - 400000
	ms = jin['bs']+jin['ys']
#	ms += [m/2 for m in jin['bs'] if m > 800000 and m > dm]
#	ms += [m/2 for m in jin['ys'] if m > 800000 and m > dm]
#	ms.sort()
	ms = [int(0.5+i/_fres) for i in ms]
	v = []
	for j in jin:
		if j == 'bs' or j == 'ys':
			continue
		v.append(jin[j])
	return (v,ms)	

#
# method to update ion series based on a set of modifications specified by:
# _pos - the array of locations (_pos), 
# _mod - the array of modifications, &
# _lp - the position in _pos
#
def update_ions(_js,_mods,_pos,_lp):
	jin = _js
	t = len(jin['bs'])
	mod_tuples = []
	tmod = 0.0
	beg = jin['beg']
	for m in _mods:
		if not _pos[m]:
			continue
		a = 0
		delta = 0
		pmod = _mods[m][_lp]
		for a in range(t):
			if a in _pos[m]:
				mod_tuples.append((beg+a,pmod))
				tmod += pmod
				delta += pmod
			if delta == 0:
				continue
			jin['bs'][a] += delta
		if t in _pos[m]:
			mod_tuples.append((beg+t,pmod))
			tmod += pmod
		a = 0
		delta = 0
		for a in range(t):
			if t-a in _pos[m]:
				delta += pmod
			if delta == 0:
				continue
			jin['ys'][a] += delta
	for tup in mod_tuples:
		jin['mods'].append({tup[0]:tup[1]})
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
	js = _js
	t = len(js['bs'])
	mods = {}
	beg = js['beg']
	mod_tuples = []
	for b in _bmods:
		js['pm'] += b
		mod_tuples.append((beg,b))
		for a in range(t):
			js['bs'][a] += b
	for tup in mod_tuples:
		js['mods'].append({tup[0]:tup[1]})
	return js

def update_yions(_js,_ymods):
	js = _js
	t = len(js['ys'])
	mods = {}
	end = js['end']
	mod_tuples = []
	for b in _ymods:
		js['pm'] += b
		mod_tuples.append((end,b))
		for a in range(t):
			js['ys'][a] += b
	for tup in mod_tuples:
		js['mods'].append({tup[0]:tup[1]})
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

def get_spectra(_index,_mass,_r,_use_slots = False,_slots = []):
	iv = int(0.5+_mass/_r)
	pms = []
	if iv in _index:
		pms += _index.get(iv)
	if _use_slots:
		for s in _slots:
			iv = int(0.5+(_mass+s)/_r)
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

