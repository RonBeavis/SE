#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# loads and interprets command line parameters
# generates error messages and returns false if parameters fail simple tests
#

import os
import json

#
# loads default parameters from a JSON formatted file if the -d parameter
# is specified on the command line
# the parameters in this file are overriden by the parameters
# specified on the command line
#

def load_defaults(_param):
	if 'parameter file' not in _param:
		return _param
	if not os.path.isfile(_param['parameter file']):
		print('Parameter file "%s" not present' % (_param['parameter file']))
		return _param
	jfile = open(_param['parameter file'],'r')
	if not jfile:
		print('Parameter file "%s" could not be opened' % (_param['parameter file']))
		return _param
	js  = '';
	for j in jfile:
		js += j
	jfile.close()
	try:
		param = json.loads(js)
	except:
		print('Parameter file "%s" not a valid JSON file' % (_param['parameter file']))
		return _param
	for p in _param:
		param[p] = _param[p]
			
	return param

#
# loads the command line parameters from sys.argv
# tests the parameters and deals with loading
# default parameters from a specified file
#
			
def load_params(_argv):
	params = {}
	ret = True
	help = False
	additional_spectra = ''
	additional_kernels = ''
	for v in _argv:
		u = v[2:]
		if v.find('-c') == 0:
			params['c13'] = False
			if u == 'yes':
				params['c13'] = True
				
		if v.find('-k') == 0:
			params['kernel file'] = u
		if v.find('-s') == 0:
			params['spectra file'] = u
		if v.find('-K') == 0:
			additional_kernels = u
		if v.find('-S') == 0:
			additional_spectra = u
		if v.find('-o') == 0:
			params['output file'] = u
		if v.find('-d') == 0:
			params['parameter file'] = u
		if v.find('-p') == 0:
			try:
				params['parent mass tolerance'] = int(u)
			except:
				params['parent mass tolerance'] = None
		if v.find('-f') == 0:
			try:
				params['fragment mass tolerance'] = int(u)
			except:
				params['fragment mass tolerance'] = None
		if v.find('-h') != -1:
			ret = False
			help = True
	if len(_argv) == 1:
		help = True
		ret = False
	if help:
		print('''
	>python3 se.py -kKERNeL -sSPECTRA (-p20) (-f400) (-oFILE) (-dFILE) (-h) (-cV)
	   where:
		   -c: use C13 isotope-error testing (yes/no)
		   -d: default parameter file (JSON)
		   -f: fragment mass tolerance in mDa (400)
	           -h: show the help list
	           -k: proteome kernel file list
		   -K: additional proteome kernel file list
		   -o: output file (tsv)
		   -p: parent mass tolerance in mDa (20)
		   -s: spectrum file list (JSMS, MGF, mzML)
		   -S: additional spectrum file list (JSMS, MGF,mzML)
''')
		return (params,False)
#
#	load parameters file, if -d specified
#
	params = load_defaults(params)
#
#	test parameters for obvious problems
#
	pval = 'kernel file'
	if pval not in params:
		print(''' no %s (-k) specified on command line''' % (pval))
		ret = False
	if pval in params:
		kfs = params[pval].split(',')
		for k in kfs:
			if not os.path.isfile(k):
				print(''' %s (-k) "%s" does not exist''' % (pval,k))
				ret = False
	pval = 'spectra file'
	if pval not in params:
		print(''' no %s (-s) specified on command line''' % (pval))
		ret = False
	if pval in params:
		sfs = params[pval].split(',')
		for s in sfs:
			if not os.path.isfile(s):
				print(''' %s (-s) "%s" does not exist''' % (pval,s))
				ret = False
	pval = 'parent mass tolerance'
	if pval in params and (params[pval] is None or params[pval] < 1) :
		print(''' %s (-p) bad value\n   must be an integer > 0''' % (pval))
		ret = False
	pval = 'parent mass tolerance'
	if pval in params and (params[pval] is None or params[pval] < 1) :
		print(''' %s (-f) bad value\n   must be an integer > 0''' % (pval))
		ret = False
#
#	add in any default values that may have been missed
#
	para_min = {'fragment mass tolerance': 400,
		'parent mass tolerance': 10,
		'mods p': {'C':[57021],'U':[57021]},
		'mods v': {'M':[15995]},
		'mods o': {'nt-ammonia':True,'nt-water':True},
		'c13': False,
		'output valid only': False}
	for p in para_min:
		if p not in params:
			params[p] = para_min[p]
	if ret:
		print('.')
	if 'kernel file' in params:
		if len(additional_kernels) > 0:
			params['kernel file'] += ',%s' % additional_kernels
	elif len(additional_kernels) > 0:
		params['kernel file'] = '%s' % additional_kernels

	if 'spectra file' in params:
		if len(additional_spectra) > 0:
			params['spectra file'] += ',%s' % additional_spectra
	elif len(additional_spectra) > 0:
		params['spectra file'] += '%s' % additional_spectra
	return (params,ret)

