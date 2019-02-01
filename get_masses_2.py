#!/usr/bin/env python

import numpy as np
from os.path import isfile

distances = {
'NGC0672': 7.6, 'NGC0891': 9.2, 'NGC0925': 9.1, 'NGC0949': 11.3,
'UGC2082': 14.4, 'NGC1003': 11.6, 'NGC2403': 3.2, 'UGC4278': 13.6,
'NGC2541': 12.0, 'NGC3198': 14.5, 'NGC4062': 16.9, 'NGC4244': 4.4,
'NGC4258': 7.6, 'NGC4274': 19.4, 'NGC4414': 17.8, 'NGC4448': 9.7,
'NGC4559': 7.9, 'NGC4565': 10.8, 'UGC7774': 24.4, 'NGC4631': 7.6,
'NGC5023': 6.6, 'NGC5055': 8.5, 'NGC5229': 5.1, 'NGC5585': 8.7 }

for galaxy in sorted(distances.keys()):
	distance = distances[galaxy]
	pfile = galaxy+'-gprof-lr.txt'
	if not isfile(pfile): continue
	vel, prof = np.loadtxt(pfile, usecols=(1,2), unpack=True)
	total = np.sum(prof)*np.abs(vel[1]-vel[0])
	mass = total * 2.36e5 * distance**2
	for line in open('globprof-log.txt'):
		if 'W50' in line and galaxy in line:
			sline = line.split()
			w50 = float(sline[3][:-1])
			w20 = float(sline[6])
	print '{id} & {d:.1f} & {s:.1f} & {m:.1f} & {w50:.1f} & {w20:.1f} \\\\'.format(id=galaxy, d=distance, s=total, m=mass/1.e8, w50=w50, w20=w20)

