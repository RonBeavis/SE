#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004

#
# loads an array with spectra obtained from named data files
#

import gzip
import json
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

def load_spectra(_in):
	test = _in.lower()
	if test.find('.mgf') == len(test)-4:
		return load_mgf(_in)
	elif test.find('.mgf.gz') == len(test)-7:
		return load_mgf(_in)
	elif test.find('.jsms') == len(test)-5:
		return load_jsms(_in)
	elif test.find('.jsms.gz') == len(test)-8:
		return load_jsms(_in)
	elif test.find('.mzml') == len(test)-5:
		return load_mzml(_in)
	elif test.find('.mzml.gz') == len(test)-8:
		return load_mzml(_in)
	return load_mgf(_in)

#
# JSMS parser
#

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
				vs.append(int(round(1000*(m-proton),0)))
			js['ms'] = vs
			sp.append(js)
			if len(sp) % 10000 == 0:
				print('.',end='')
				sys.stdout.flush()
	f.close()
	sp = clean_up(sp)
	return sp
#
# MGF parser
#
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
	found = set([])
	print('.',end='')
	sys.stdout.flush()
	for line in ifile:
		line = line.rstrip()
		if 'BEGIN IONS' not in found and line.find('BEGIN IONS') == 0:
			js = {}
			Ms = []
			Is = []
			amIn = 1
			found.add('BEGIN IONS')
		elif amIn == 0:
			continue
		elif line.find('END IONS') == 0:
			js['np'] = len(Ms)
			js['pm'] = int(round(1000*(js['pm']*js['pz']-proton*js['pz']),0))
			js['ms'] = Ms
			js['is'] = Is
			sp.append(js)
			if len(sp) % 10000 == 0:
				print('.',end='')
				sys.stdout.flush()
			amIn = 0
			found = set([])
			s += 1
		elif 'PEPMASS' not in found and line.find('PEPMASS') == 0:
			line = re.sub('^PEPMASS\=','',line)
			vs = line.split(' ')
			js['pm'] = float(vs[0])
			if len(vs) > 1:
				js['pi'] = float(vs[1])
			found.add('PEPMASS')
		elif 'RTINSECONDS' not in found and line.find('RTINSECONDS') == 0:
			line = re.sub('^RTINSECONDS\=','',line)
			js['rt'] = float('%.3f' % float(line))
			found.add('RTINSECONDS')
		elif 'CHARGE=' not in found and line.find('CHARGE=') == 0:
			line = re.sub('^CHARGE\=','',line)
			if line.find('+'):
				line = line.replace('+','')
				js['pz'] = int(line)
			else:
				line = line.replace('-','')
				js['pz'] = -1*int(line)
			found.add('CHARGE=')
		elif 'TITLE=' not in found and line.find('TITLE=') == 0:
			line = re.sub('^TITLE=','',line)
			line = re.sub('\"','\'',line)
			js['ti'] = line
			if line.find('scan=') != -1:
				js['sc'] = int(re.sub('.+scan\=','',line))
			else:
				js['sc'] = s+1
			found.add('TITLE=')
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
	
	def getSpectra(self):
		return self.spectra	
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
					m = int(round(1000.0*(float(self.jsms['ms'][a])-self.proton),0))
					Ms.append(m)
					Is.append(self.jsms['is'][a])
					if 'zs' in self.jsms:
						Zs.append(self.jsms['zs'][a])
					a += 1
				self.jsms['ms'] = Ms
				self.jsms['is'] = Is
				self.jsms['pm'] = int(round(1000*(self.jsms['pm']*self.jsms['pz']-self.proton*self.jsms['pz']),0))
				if 'zs' in self.jsms:
					self.jsms['zs'] = Zs
				self.jsms['np'] = len(Ms)
				self.spectra.append(self.jsms)
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

def load_mzml(_in):
	fpath = _in
	parser = xml.sax.make_parser()
	parser.setContentHandler(mzMLHandler())
	parser.parse(open(fpath,"r"))
	sp = parser.getContentHandler().getSpectra()
	clean_up(sp)
	return sp

#
# cleans up spectra to conform to search engine requirements
#

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
		while sIs[-1] < minI and len(sMs) > 0:
			sMs.pop(-1)
			sIs.pop(-1)
			deleted += 1
		sIs = [x for _,x in sorted(zip(sMs,sIs))]
		sMs.sort()
		sp[a]['ms'] = sMs
		sp[a]['is'] = sIs
		a += 1
	return sp

