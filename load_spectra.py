#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004

#
# loads an array with spectra obtained from named data files
#
#from __future__ import print_function
#from libcpp cimport bool as bool_t

import gzip
import ujson
#import json
import re
import sys
import struct
import base64
import hashlib
import gzip
import zlib
import xml.sax

#
# method matches file name with parser and loads array
# only externally called method
#

def load_spectra(_in,_param):
	test = _in.lower()
	if test.find('.mgf') == len(test)-4:
		return load_mgf(_in,_param)
	elif test.find('.mgf.gz') == len(test)-7:
		return load_mgf(_in,_param)
	elif test.find('.jsms') == len(test)-5:
		return load_jsms(_in,_param)
	elif test.find('.jsms.gz') == len(test)-8:
		return load_jsms(_in,_param)
	elif test.find('.mzml') == len(test)-5:
		return load_mzml(_in,_param)
	elif test.find('.mzml.gz') == len(test)-8:
		return load_mzml(_in,_param)
	return load_mgf(_in,_param)

#
# JSMS parser
#

def load_jsms(_in,_param):
	sp = []
	if _in.find('.gz') == len(_in) - 3:
		f = gzip.open(_in,'rt',encoding = 'utf8')
	else:
		f = open(_in,'r',encoding = 'utf8')
	proton = 1.007276
	res = float(_param['fragment mass tolerance'])
	m = 0.0
	for l in f:
		js = ujson.loads(l)
		if 'lv' in js and 'pz' in js and js['pm']*js['pz'] > 600:
			js['pm'] = int(0.5+1000*(js['pm']*js['pz']-proton*js['pz']))
			ms = js['ms']
			vs = []
			for m in ms:
				vs.append(int(0.5 + 1000.0*(float(m)-proton)))
			js['ms'] = vs
			js = clean_one(js,50,res)
			sp.append(js)
			if len(sp) % 10000 == 0:
				print('.',end='')
				sys.stdout.flush()
	f.close()
	return sp
#
# MGF parser
#
def load_mgf(_in,_param):
	if _in.find('.gz') == len(_in) - 3:
		ifile = gzip.open(_in,'rt')
	else:
		ifile = open(_in,'r')
	s = 1
	js = {}
	jc = {}
	sp = []
	Ms = []
	Is = []
	MsPos = 0
	amIn = 0
	proton = 1.007276
	found = set([])
	print('.',end='')
	sys.stdout.flush()
	digit = set(['1','2','3','4','5','6','7','8','9'])
	res = float(_param['fragment mass tolerance'])
	m = 0.0
	i = 0.0
	for line in ifile:
		fl = line[:1]
		if amIn and fl in digit:
			vs = line.split(' ')
			if len(vs) < 2 or len(vs) > 3:
				continue
			if float(vs[1]) <= 0.0:
				continue
			m = int(0.5 + 1000.0*(float(vs[0])-proton))
			if m <= 0:
				continue
			i = float(vs[1])
			if i <= 0.0:
				continue
			if MsPos % 1000 == 0:
				ts = [0 for x in range(1000)]
				Ms.extend(ts)
				Is.extend(ts)
			Ms[MsPos] = m
			Is[MsPos] = i
			MsPos += 1
			continue
		gr = re.match('([\w ]+)',line)
		if not gr:
			continue
		tag = gr.group(1)
		if tag == 'BEGIN IONS':
			js = {}
			Ms = []
			Is = []
			MsPos = 0
			amIn = 1
		elif amIn == 0:
			continue
		elif tag == 'END IONS':
			if js['pm']*js['pz'] > 600:
				js['np'] = len(Ms)
				js['pm'] = int(0.5 + 1000*(js['pm']*js['pz']-proton*js['pz']))
				js['ms'] = Ms[:MsPos]
				js['is'] = Is[:MsPos]
				jc = clean_one(js,50,res)
				sp.append(jc)
			if len(sp) % 10000 == 0:
				print('.',end='')
				sys.stdout.flush()
			amIn = 0
			s += 1
		elif tag == 'PEPMASS':
			line = re.sub('^PEPMASS\=','',line)
			vs = line.split(' ')
			js['pm'] = float(vs[0])
			if len(vs) > 1:
				js['pi'] = float(vs[1])
		elif tag == 'RTINSECONDS':
			line = re.sub('^RTINSECONDS\=','',line)
			js['rt'] = float('%.3f' % float(line))
		elif  tag == 'CHARGE':
			line = re.sub('^CHARGE\=','',line)
			if line.find('+'):
				line = line.replace('+','')
				js['pz'] = int(line)
			else:
				line = line.replace('-','')
				js['pz'] = -1*int(line)
		elif  tag == 'TITLE':
			line = re.sub('^TITLE=','',line)
			line = re.sub('\"','\'',line)
			js['ti'] = line
			if line.find('scan=') != -1:
				js['sc'] = int(re.sub('.+scan\=','',line))
			else:
				js['sc'] = s+1
	return sp
#
#	mzML parser
#

#
#	mzMLHandler necessary to XML SAX handler
#
class mzMLHandler(xml.sax.ContentHandler):
	def __init__(self):
		self.cTag = ''
		self.isSpectrum = False
		self.isSelectedIon = False
		self.isMzArray = False
		self.isIntArray = False
		self.isBinaryDataArray = False
		self.isBinary = False
		self.isZlib = False
		self.isScan = False
		self.jsms = {}
		self.floatBytes = 8
		self.content = ''
		self.mhash = None
		self.spectra = []
		self.n = 0
		self.proton = 1.007276
		self.res = 10
	
	def getSpectra(self):
		return self.spectra	
	def setRes(self,_r):
		self.res = _r
		return self.res
	def cleanOne(self):
		s = self.jsms
		ires = self.res
		a = 0
		i_max = 0.0
		_l = 50
		pm = s['pm']
		pz = s['pz']
		m = 0
		for a,i in enumerate(s['is']):
			m = s['ms'][a]
			if m < 150000 or abs(pm-m) < 45000 or abs(pm/pz- m) < 2000 :
				continue
			if i > i_max:
				i_max = float(i)
		sMs = []
		sIs = []
		i_max =float(i_max)/100.0
		for a,m in enumerate(s['ms']):
			if m < 150000 or abs(pm-m) < 45000 or abs(pm/pz- m) < 2000 :
				continue
			if float(s['is'][a])/i_max > 1.0:
				i = s['is'][a]/i_max
				if sMs:
					if abs(sMs[-1] - m) < 500:
						if sIs[-1] < i:
							del sMs[-1]
							del sIs[-1]
							sMs.append(m)
							sIs.append(i)
					else:
						sMs.append(m)
						sIs.append(i)
				else:
					sMs.append(m)
					sIs.append(i)
		sMs = [x for _,x in sorted(zip(sIs,sMs), key=lambda pair: pair[0],reverse = True)]
		sIs.sort(reverse = True)
		max_l = 2 * int(0.5 + float(pm)/100000.0)
		if max_l > _l:
			max_l = _l
		sMs = sMs[:max_l]
		sIs = sIs[:max_l]
		sIs = [x for _,x in sorted(zip(sMs,sIs), key=lambda pair: pair[0])]
		sMs.sort()
		s.pop('ms')
		s.pop('is')
		tps = []
		ips = []
		val = 0
	#
	#		generate a normalized set of spectrum masses
	#
		for a,m in enumerate(sMs):
			val = int(0.5+float(m)/ires)
			tps.append(val)
			tps.append(val-1)
			tps.append(val+1)
			ips.extend([sIs[a],sIs[a],sIs[a]])
		s['sms'] = tps
		s['ims'] = ips
		s['isum'] = sum(sIs)
		return s

	def startElement(self, tag, attrs):
		self.cTag = tag
		if tag == 'spectrum':
			self.isSpectrum = True
			if 'scan' in attrs:
				self.jsms['sc'] = int(attrs['scan'])
			elif 'index' in attrs:
				self.jsms['sc'] = int(attrs['index'])
		if tag == 'precursor' and 'spectrumRef' in attrs:
			self.jsms['ti'] = attrs['spectrumRef']
			if self.jsms['ti'].find('scan=') != -1:
				grp = re.match(r'.+scan=(\d+)',self.jsms['ti'])
				self.jsms['sc'] = int(grp.group(1))
		if tag == 'scan':
			self.isScan = True
		if self.isScan and tag == 'cvParam':
			if attrs['name'] == 'filter string' and 'ti' not in self.jsms:
				self.jsms['ti'] = attrs['value']
			if attrs['name'] == 'scan start time':
				self.jsms['rt'] = float('%.3f' % (60.0*float(attrs['value'])))
		if self.isSpectrum and tag == 'cvParam':
			if attrs['name'] == 'ms level':
				self.jsms['lv'] = int(attrs['value'])
		if tag == 'selectedIon':
			self.isSelectedIon = True
		if self.isSelectedIon and tag == 'cvParam':
			if attrs['name'] == 'selected ion m/z':
				self.jsms['pm'] = float('%.4f' % float(attrs['value']))
			if attrs['name'] == 'charge state':
				self.jsms['pz'] = int(attrs['value'])
			if attrs['name'] == 'peak intensity':
				self.jsms['pi'] = float(attrs['value'])
		if tag == 'binaryDataArray':
			self.isBinaryDataArray = True
		if self.isBinaryDataArray and tag == 'cvParam':
			if attrs['name'] == 'm/z array':
				self.isMzArray = True
			if attrs['name'] == 'intensity array':
				self.isIntArray = True
			if attrs['name'] == '32-bit float':
				self.floatBytes = 4
			if attrs['name'] == '64-bit float':
				self.floatBytes = 8
			if attrs['name'] == 'zlib compression':
				self.isZlib = True
		if (self.isMzArray or self.isIntArray) and tag == 'binary':
			self.isBinary = True
	def characters(self, content):
		if self.isMzArray or self.isIntArray:
			self.content += content

	def endElement(self, tag):
		if tag == 'spectrum':
			self.isSpectrum = False
			if len(self.jsms) > 0 and self.jsms['lv'] >= 2:
				a = 0
				Ms = []
				Is = []
				Zs = []
				while a < self.jsms['np']:
					if self.jsms['is'][a] <= 0.0:
						a += 1
						continue
					m = int(0.5 + 1000.0*(float(self.jsms['ms'][a])-self.proton))
					Ms.append(m)
					Is.append(self.jsms['is'][a])
					if 'zs' in self.jsms:
						Zs.append(self.jsms['zs'][a])
					a += 1
				if self.jsms['pm']*self.jsms['pz'] > 600:
					self.jsms['ms'] = Ms
					self.jsms['is'] = Is
					self.jsms['pm'] = int(0.5 + 1000*(self.jsms['pm']*self.jsms['pz']-self.proton*self.jsms['pz']))
					if 'zs' in self.jsms:
						self.jsms['zs'] = Zs
					self.jsms['np'] = len(Ms)
					sp = self.cleanOne()
					self.spectra.append(sp)
				self.n += 1
				if self.n % 10000 == 0:
					print('.',end='',flush=True)

			self.jsms = {}
		if tag == 'scan':
			self.isScan = False
		if tag == 'selectedIon':
			self.isSelectedIon = False
		if tag == 'binaryDataArray':
			self.isBinaryDataArray = False
		if tag == 'binary':
			self.isBinary = False
			if self.isMzArray:
				str = self.content.strip()
				if self.isZlib:
					d = base64.standard_b64decode(str.encode())
					data = zlib.decompress(d)
					count = len(data)/int(self.floatBytes)
					if self.floatBytes == 4:
						result = struct.unpack('<%if' % (count),data)
					else:
						result = struct.unpack('<%id' % (count),data)
				else:
					data = base64.standard_b64decode(str.encode())
					count = len(data)/int(self.floatBytes)
					if self.floatBytes == 4:
						result = struct.unpack('<%if' % (count),data)
					else:
						result = struct.unpack('<%id' % (count),data)

				self.jsms['ms'] = result
				self.jsms['np'] = len(result)
				self.isMzArray = False
				self.isZlib = False
				self.content = ''
			if self.isIntArray:
				str = self.content.strip()
				if self.isZlib:
					d = base64.standard_b64decode(str.encode())
					data = zlib.decompress(d)
					count = len(data)/int(self.floatBytes)
					if self.floatBytes == 4:
						result = struct.unpack('<%if' % (count),data)
					else:
						result = struct.unpack('<%id' % (count),data)
				else:
					data = base64.standard_b64decode(str.encode())
					count = len(data)/int(self.floatBytes)
					if self.floatBytes == 4:
						result = struct.unpack('<%if' % (count),data)
					else:
						result = struct.unpack('<%id' % (count),data)
				self.jsms['is'] = result
				self.isIntArray = False
				self.content = ''

#
# Creates the xml.sax parser using the mzMLHandler class
# and parses the mzML file
#

def load_mzml(_in,_param):
	fpath = _in
	res = float(_param['fragment mass tolerance'])
	parser = xml.sax.make_parser()
	parser.setContentHandler(mzMLHandler())
	parser.getContentHandler().setRes(res)
	parser.parse(open(fpath,"r"))
	sp = parser.getContentHandler().getSpectra()
	return sp

def clean_one(_sp,_l,_ires):
	s = _sp.copy()
	a = 0
	i_max = 0.0
	pm = s['pm']
	pz = s['pz']
	m = 0
	for a,i in enumerate(s['is']):
		m = s['ms'][a]
		if m < 150000 or abs(pm-m) < 45000 or abs(pm/pz- m) < 2000 :
			continue
		if i > i_max:
			i_max = float(i)
	sMs = []
	sIs = []
	i_max =float(i_max)/100.0
	for a,m in enumerate(s['ms']):
		if m < 150000 or abs(pm-m) < 45000 or abs(pm/pz- m) < 2000 :
			continue
		if float(s['is'][a])/i_max > 1.0:
			i = s['is'][a]/i_max
			if sMs:
				if abs(sMs[-1] - m) < 500:
					if sIs[-1] < i:
						del sMs[-1]
						del sIs[-1]
						sMs.append(m)
						sIs.append(i)
				else:
					sMs.append(m)
					sIs.append(i)
			else:
				sMs.append(m)
				sIs.append(i)
	sMs = [x for _,x in sorted(zip(sIs,sMs), key=lambda pair: pair[0],reverse = True)]
	sIs.sort(reverse = True)
	max_l = 2 * int(0.5 + float(pm)/100000.0)
	if max_l > _l:
		max_l = _l
	sMs = sMs[:max_l]
	sIs = sIs[:max_l]
	sIs = [x for _,x in sorted(zip(sMs,sIs), key=lambda pair: pair[0])]
	sMs.sort()
#	sp[a]['ms'] = sMs
#	sp[a]['is'] = sIs
	s.pop('ms')
	s.pop('is')
	tps = []
	ips = []
	val = 0
#
#		generate a normalized set of spectrum masses
#
	for a,m in enumerate(sMs):
		val = int(0.5+float(m)/_ires)
		tps.append(val)
		tps.append(val-1)
		tps.append(val+1)
		ips.extend([sIs[a],sIs[a],sIs[a]])
	s['sms'] = tps
	s['ims'] = ips
	s['isum'] = sum(sIs)
	return s

