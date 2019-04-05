#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

def load_params(_argv,_d):
	params = _d
	ret = True
	for v in _argv:
		u = v[2:]
		if v.find('-k') == 0:
			params['kernel'] = u
		if v.find('-s') == 0:
			params['spectra'] = u
		if v.find('-p') == 0:
			params['pres'] = int(u)
		if v.find('-f') == 0:
			params['res'] = int(u)
		if v.find('-h') != -1:
			print('''
	>python3 se.py -kKERNeL -sSPECTRA (-p20) (-f400)
	   where:
	           -k: proteome kernel file
		   -s: spectrum file (JSMS, MGF)
		   -p: parent mass tolerance in mDa (20)
		   -f: fragment mass tolerance in mDa (400)''')
			ret = False
			break
	return (params,ret)

