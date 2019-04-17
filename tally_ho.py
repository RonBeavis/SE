import re
import sys

f = open(sys.argv[1],'r')
count_tally = {}
rev_count_tally = {}
score_tally = {}
rev_score_tally = {}
scores = {}
counts = {}
for l in f:
	if l.find('spectrum') != -1:
		continue
	l.rstrip()
	vs = l.split('\t')
	brev = False
	if vs[6].find('rev-') == 0:
		brev = True
	score = int(50*int(float(vs[14])/50))
	count = int(vs[13])
	scores[score] = 1
	counts[count] = 1
	if brev:
		if score in rev_score_tally:
			rev_score_tally[score] += 1
		else:
			rev_score_tally[score] = 1
		if count in rev_count_tally:
			rev_count_tally[count] += 1
		else:
			rev_count_tally[count] = 1
	else:
		if score in score_tally:
			score_tally[score] += 1
		else:
			score_tally[score] = 1
		if count in count_tally:
			count_tally[count] += 1
		else:
			count_tally[count] = 1
print('scores\treversed\tforward')
a_score_tally = {}
rev_a_score_tally = {}
for sc in sorted(scores):
	forward = score_tally.get(sc,0)
	rev = rev_score_tally.get(sc,0)
	a_score_tally[sc] = forward
	rev_a_score_tally[sc] = rev
	for s in sorted(scores):
		if s > sc:
			frw = score_tally.get(s,0)
			rv = rev_score_tally.get(s,0)
			a_score_tally[sc] += frw
			rev_a_score_tally[sc] += rv

a_count_tally = {}
rev_a_count_tally = {}
for sc in sorted(counts):
	forward = count_tally.get(sc,0)
	rev = rev_count_tally.get(sc,0)
	a_count_tally[sc] = forward
	rev_a_count_tally[sc] = rev
	for s in sorted(counts):
		if s > sc:
			frw = count_tally.get(s,0)
			rv = rev_count_tally.get(s,0)
			a_count_tally[sc] += frw
			rev_a_count_tally[sc] += rv
			
for sc in sorted(scores):
	line = '%i\t' % (sc)
	rev = 0
	forward = 0
	if sc in rev_score_tally:
		line += '%i\t' % (rev_score_tally[sc])
		rev = rev_score_tally[sc]
	else:
		line += '0\t'
	if sc in score_tally:
		line += '%i\t' % (score_tally[sc])
		forward = score_tally[sc]
	else:
		line += '0\t'
	if forward+rev > 0:
		line += '%.1f\t' % (100*float(rev)/float(rev+forward))
	else:
		line += '-\t'
	line += '%.3f'% (100*float(rev_a_score_tally[sc])/float(rev_a_score_tally[sc]+a_score_tally[sc]))
	print(line)
print('counts\treversed\tforward')
for sc in sorted(counts):
	rev = 0
	forward = 0
	line = '%i\t' % (sc)
	if sc in rev_count_tally:
		line += '%i\t' % (rev_count_tally[sc])
		rev = rev_count_tally[sc]
	else:
		line += '0\t'
	if sc in count_tally:
		line += '%i\t' % (count_tally[sc])
		forward = count_tally[sc]
	else:
		line += '0\t'
	if forward+rev > 0:
		line += '%.1f\t' % (100*float(rev)/float(rev+forward))
	else:
		line += '-\t'
	line += '%.3f'% (100*float(rev_a_count_tally[sc])/float(rev_a_count_tally[sc]+a_count_tally[sc]))
	print(line)

