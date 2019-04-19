import json
import re
import datetime
import hashlib
import mysql.connector
import sys
import random

mhash = hashlib.sha256()
dhash = hashlib.sha256()

isotopes = {'p' : 1.007276,
'H' : 1.007825,
'C' : 12.0,
'N' : 14.003074,
'O': 15.994915 }

a_to_m = {'A':71.037114,
'R':156.101111,
'N':114.042927,
'D':115.026943,
'C':103.009185,
'E':129.042593,
'Q':128.058578,
'G':57.021464,
'H':137.058912,
'I':113.084064,
'L':113.084064,
'K':128.094963,
'M':131.040485,
'F':147.068414,
'P':97.052764,
'S':87.032028,
'T':101.047679,
'U':150.95363,
'W':186.079313,
'Y':163.06332,
'V':99.068414 }

def create_protein(_label,_p,_l,_ld,_c,_term):
	entries = 0
	rejects = 0
	if _p.find('*') != -1:
		return (0,0,0)
	if len(_p) < 20:
		return (0,0,0)
	_p = re.sub('X','',_p)
	_p = re.sub('B','N',_p)
	_p = re.sub('Z','Q',_p)
	rs = list(_p)
	ms = []
	for r in rs:
		ms.append(a_to_m[r])
	l_rs = len(rs)
	a = 0
	term = _term
	dB = 0.0
	dY = 2*isotopes['H'] + isotopes['O']
	water = 2*isotopes['H'] + isotopes['O']
	proton = isotopes['p']
	redundant = {}
	non_tryptic = 0
	while a < l_rs:
		i = a
		m = 0
		pep = ''
		mpep = []
		start = a + 1
		non_tryptic = False
		while m <= 2 and i < l_rs:
			pep += rs[i]
			mpep.append(ms[i])
			if term.find(rs[i]) != -1 or i+1 == l_rs:
				if m == 0:
					if a < 120:
						non_tryptic = True
						a += 1
					else:
						a = i + 1
				m += 1
				if i - start < 6:
					i += 1
					continue
				sql = 'select z1total,z2total,z3total,z4total from protein_omega_count where seq=%(seq)s'
				_c.execute(sql,{'seq':pep})
				vs = _c.fetchall()
				bail = True
				zs = [0]*4
				for v in vs:
					y = 0
					for x in v:
						if x >= 10:
							bail = False
						zs[y] += x
						y += 1
				if bail == True:
					i += 1
					if start > 1 and term.find(rs[start-2]) != -1:
						rejects += 1
					continue
				frags = []
				ion = dB
				bions = []
				parent = water
				js = {}
				decoy_pep = mpep.copy()
				for p in mpep:
					ion += p
					parent += p
					bions.append(int(round(1000.0*ion,0)))
				ion = dB
				random.shuffle(decoy_pep)
				decoy_bions = []
				decoy_js = {}
				decoy_frags = []
				for p in decoy_pep:
					ion += p
					decoy_bions.append(int(round(1000.0*ion,0)))
				yions = []
				ion = dY
				mpep.reverse()
				decoy_pep = mpep.copy()
				for p in mpep:
					ion += p
					yions.append(int(round(1000.0*ion,0)))
				ion = dY
				random.shuffle(decoy_pep)
				decoy_yions = []
				for p in decoy_pep:
					ion += p
					decoy_yions.append(int(round(1000.0*ion,0)))
				mpep.reverse()
				js['lv'] = 0
				js['pm'] = int(round(1000.0*parent,0))
				js['lb'] = _label
				if start - 2 < 0:
					js['pre'] = '['
				else:
					js['pre'] = rs[start-2]
				if i + 1 >= l_rs:
					js['post'] = ']'
				else:
					js['post'] = rs[i+1]
				js['beg'] = start
				js['end'] = i + 1
				js['seq'] = pep
				js['ns'] = zs
				ss = '%i %i' % (start,i+1)
				js['bs'] = bions[:-1]
				js['ys'] = yions[:-1]
				if ss in redundant:
					print('.',end='')
					continue
				redundant[ss] = True
				line = json.dumps(js)
				mhash.update(line.strip().encode())
				_l.write(line + '\n')
				decoy_js = js.copy()
				decoy_js['lb'] = 'decoy-' + js['lb']
				decoy_js['bs'] = decoy_bions[:-1]
				decoy_js['ys'] = decoy_yions[:-1]
				line = json.dumps(decoy_js)
				dhash.update(line.strip().encode())
				_ld.write(line + '\n')
				if start != 1 and term.find(js['pre']) != -1:
					entries += 1
				else:
					non_tryptic += 1
			i += 1
	return (entries,rejects,non_tryptic)

def generate_kernel(_f,_l,_ld,_c,_term):
	f = open(_f,'r')
	ofile = open(_l,'a')
	dfile = open(_ld,'a')
	label = ''
	protein = ''
	p = 0
	entries = 0
	rejects = 0
	non_tryptic = 0
	for l in f:
		l = l.strip()
		if l.find('>') == 0:
			vs = create_protein(label,protein,ofile,dfile,_c,_term)
			entries += vs[0]
			rejects += vs[1]
			non_tryptic += vs[2]
			if l.find('>sp|') == 0 or l.find('>tp|') == 0 or l.find('>gi|') == 0:
				m = re.search('\>(\w\w\|\w+\|)',l)
			else:
				m = re.search('\>(.+?)[\._\s]',l)
#			m = re.search('\>(.+?)$',l)
			if m is None:
				print(l)
			label = m.group(1)
			protein = ''
			p += 1
			if p % 500 == 0:
				print(p,entries,rejects,'%.1f' % (100*float(rejects)/float(entries+rejects)),non_tryptic)
		else:
			protein += l
	f.close()
	if len(protein) > 0 and len(label) > 0:
		vs = create_protein(label,protein,ofile,dfile,_c,_term)
		entries += vs[0]
		rejects += vs[1]
		non_tryptic += vs[2]
	ofile.close()
	return (entries,rejects,non_tryptic)

def add_format(_l,_v):
	ofile = open(_l,'w')
	info = {'format':'jsms 1.0','source':_l,'created':str(datetime.datetime.now())}
	line = json.dumps(info) + '\n'
	ofile.write(line)
	if _v:
		mhash.update(line.strip().encode())
	else:
		dhash.update(line.strip().encode())
	ofile.close()

def add_validation(_l,_v):
	ofile = open(_l,'a')
	if _v:
		info = {"validation" : "sha256", "value" : mhash.hexdigest()}
	else:
		info = {"validation" : "sha256", "value" : dhash.hexdigest()}
	line = json.dumps(info) + '\n'
	ofile.write(line)
	ofile.close()

try:
	conn = mysql.connector.connect(host='192.168.1.4',
					database='peakdb',
					user='thehome87team',
					password='hunt45andpeck') 
except:
	print('{"error":"Could not connect to database"}\n')
	exit()

curses = conn.cursor()

reverse = False
term = 'KR'
if len(sys.argv) > 2:
	fasta = sys.argv[1]
	a = 2
	while a < len(sys.argv):
		if len(sys.argv[a]) < 2:
			continue
		u = sys.argv[a][2:]
		if sys.argv[a].find('-t') == 0:
			term = u
		a += 1
else:
	print('kernelizer.py FASTA -tKR')
	exit()
print('Terminator: %s' % (term))
lib = re.sub('/fasta','/kernels',fasta)
lib = re.sub('\.fasta','.' + term + '.kernel',lib)
lib = re.sub('\.pep\.all','',lib)
lib_decoy = re.sub('\.kernel','.decoy.kernel',lib)
add_format(lib,True)
add_format(lib_decoy,False)
vs = generate_kernel(fasta,lib,lib_decoy,curses,term)
curses.close()
conn.close()
add_validation(lib,True)
add_validation(lib_decoy,False)
print('total entries = %i, rejects = %i (%.1f), non_tryptic = %i' % (vs[0],vs[1],100*float(vs[1])/float(vs[0]+vs[1]),vs[2]))
