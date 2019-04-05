#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import gzip
import json
import re

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
		o = json.loads(l)
		if 'lv' in o and 'pz' in o:
			o['pm'] = int(round(1000*(o['pm'] - proton)*o['pz'],0))
			ms = o['ms']
			vs = []
			for m in ms:
				vs.append(int(round(1000*(m-protein),0)))
			vs.sort()
			o['ms'] = vs
			sp.append(o)
	f.close()
	return sp

def load_mgf(_in):
	if _in.find('.gz') == len(_in) - 3:
		ifile = gzip.open(_in,'rt')
	else:
		ifile = open(_in,'r')
	s = 0
	js = {}
	spectra = []
	Ms = []
	amIn = 0
	proton = 1.007276
	for line in ifile:
		line = line.rstrip()
		if line.find('BEGIN IONS') == 0:
			js = {}
			Ms = []
			amIn = 1
		elif amIn == 0:
			continue
		elif line.find('END IONS') == 0:
			js['np'] = len(Ms)
			js['pm'] = 1000.0*((js['pm']-proton)*js['pz'])
			js['pm'] = int(round(js['pm'],0))
			if len(Ms) > 0:
				Ms.sort()
				js['ms'] = Ms
			spectra.append(js)
			amIn = 0
			s += 1
		elif line.find('PEPMASS') == 0:
			line = re.sub('^PEPMASS\=','',line)
			vs = line.split(' ')
			js['pm'] = float('%.4f' % float(vs[0]))
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
	return spectra

