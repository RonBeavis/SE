#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import gzip
import json
import re
import sys

def load_spectra(_in):
	if _in.find('.mgf') == len(_in)-4:
		return load_mgf(_in)
	elif _in.find('.mgf.gz') == len(_in)-7:
		return load_mgf(_in)
	elif _in.find('.jsms') == len(_in)-5:
		return load_jsms(_in)
	elif _in.find('.jsms.gz') == len(_in)-8:
		return load_jsms(_in)
	return load_mgf(_in)

def load_jsms(_in):
	sp = []
	if _in.find('.gz') == len(_in) - 3:
		f = gzip.open(_in,'rt',encoding = 'utf8')
	else:
		f = open(_in,'r',encoding = 'utf8')
	proton = 1.007276
	for l in f:
		js = json.loads(l)
		if 'lv' in js and 'pz' in js:
			js['pm'] = int(round(1000*(js['pm']*js['pz']-proton*js['pz']),0))
			ms = js['ms']
			vs = []
			for m in ms:
				vs.append(int(round(1000*(m-protein),0)))
			js['ms'] = vs
			sp.append(js)
			if len(sp) % 1000 == 0:
				print('.',end='')
				sys.stdout.flush()
	f.close()
	sp = clean_up(sp)
	return sp

def load_mgf(_in):
	if _in.find('.gz') == len(_in) - 3:
		ifile = gzip.open(_in,'rt')
	else:
		ifile = open(_in,'r')
	s = 0
	js = {}
	sp = []
	Ms = []
	Is = []
	amIn = 0
	proton = 1.007276
	for line in ifile:
		line = line.rstrip()
		if line.find('BEGIN IONS') == 0:
			js = {}
			Ms = []
			Is = []
			amIn = 1
		elif amIn == 0:
			continue
		elif line.find('END IONS') == 0:
			js['np'] = len(Ms)
			js['pm'] = int(round(1000*(js['pm']*js['pz']-proton*js['pz']),0))
			js['ms'] = Ms
			js['is'] = Is
			sp.append(js)
			if len(sp) % 1000 == 0:
				print('.',end='')
				sys.stdout.flush()
			amIn = 0
			s += 1
		elif line.find('PEPMASS') == 0:
			line = re.sub('^PEPMASS\=','',line)
			vs = line.split(' ')
			js['pm'] = float(vs[0])
			if len(vs) > 1:
				js['pi'] = float(vs[1])
		elif line.find('RTINSECONDS') == 0:
			line = re.sub('^RTINSECONDS\=','',line)
			js['rt'] = float('%.3f' % float(line))
		elif line.find('CHARGE=') == 0:
			line = re.sub('^CHARGE\=','',line)
			if line.find('+'):
				line = line.replace('+','')
				js['pz'] = int(line)
			else:
				line = line.replace('-','')
				js['pz'] = -1*int(line)
		elif line.find('TITLE=') == 0:
			line = re.sub('^TITLE=','',line)
			line = re.sub('\"','\'',line)
			js['ti'] = line
			if line.find('scan=') != -1:
				js['sc'] = int(re.sub('.+scan\=','',line))
			else:
				js['sc'] = s+1
		else:
			vs = line.split(' ')
			if len(vs) < 2 or len(vs) > 3:
				continue
			z = 0
			if float(vs[1]) <= 0.0:
				continue
			m = int(round(1000.0*(float(vs[0])-proton),0))
			if m <= 0:
				continue
			i = float(vs[1])
			if i <= 0:
				continue						
			Ms.append(m)
			Is.append(i)
	sp = clean_up(sp)
	return sp

def clean_up(_sp,l = 50):
	sp = _sp
	a = 0
	deleted = 0
	for s in sp:
		sMs = [x for _,x in sorted(zip(s['is'],s['ms']),reverse = True)]
		sIs = s['is']
		sIs.sort(reverse = True)
		sMs = sMs[:l]
		sIs = sIs[:l]
		minI = sIs[0]/100
		while sIs[-1] < minI:
			sMs.delete(-1)
			sIs.delete(-1)
			deleted += 1
		sIs = [x for _,x in sorted(zip(sMs,sIs))]
		sMs.sort()
		sp[a]['ms'] = sMs
		sp[a]['is'] = sIs
		a += 1
	return sp

