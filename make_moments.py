#!/usr/bin/env python

from glob import glob
from astropy.io import fits
from os import getcwd
from subprocess import call
from tqdm import tqdm

# This sets the clip value for source finding (in sigmas)
clipval = 4.

#cubelist = sorted(glob('*cube.fits'))
cubelist = sorted(glob('*LR-cube.fits'))
pipeline_cmd = '/data/LIBRA_1/hea189/SoFiA/SoFiA-1.2.0/sofia_pipeline.py'

def getkernelxy(hdu, scale):
	header = hdu[0].header
	bmaj = header['BMAJ']
	bmin = header['BMIN']
	bpa = header['BPA']
	pixx = abs(header['CDELT1'])
	pixy = header['CDELT2']
	if abs(bpa) < 45.:
		kx = round(scale*bmin/pixx)
		ky = round(scale*bmaj/pixy)
	else:
		kx = round(scale*bmaj/pixx)
		ky = round(scale*bmin/pixy)
	return int(kx), int(ky)

sofia_template = 'sf_params.txt'
cwd = getcwd()

pbar = tqdm(cubelist)
for cube in pbar:
	name = cube[:10]
	pbar.set_description('Processing %s'%name)

	if 'LR' in name:
		scale = 1.
	else:
		scale = 1.
	
	hdu = fits.open(cube)
	kx, ky = getkernelxy(hdu, scale)
	hdu.close()

	sofia_file = 'MOMENTS/'+name+'-sfpars.txt'
	sofia_log = 'MOMENTS/'+name+'-sflog.txt'
	sfout = open(sofia_file,'w')
	
	for line in open(sofia_template):
		sline = line.split()
		if len(sline) == 0:
			# blank line in file
			print >>sfout, ""
			continue
		if sline[0] == 'SCfind.kernels':
			print >>sfout, "SCfind.kernels\t=\t[[0, 0, 0, 'b'], [0, 0, 1, 'b'], [0, 0, 2, 'b'], [{x:d}, {y:d}, 0, 'b'], [{x:d}, {y:d}, 1, 'b'], [{x:d}, {y:d}, 2, 'b']]".format(x=kx,y=ky)
		elif sline[0] == 'import.inFile':
			print >>sfout, "import.inFile\t=\t{cwd}/{cube}".format(cwd=cwd,cube=cube)
		elif sline[0] == 'SCfind.threshold':
			print >>sfout, "SCfind.threshold\t=\t{clip:.1f}".format(clip=clipval)
		elif sline[0] == 'writeCat.basename':
			print >>sfout, "writeCat.basename\t=\t{name}".format(name=name)
		else:
			print >>sfout, line.rstrip()
	sfout.close()
	logf = open(sofia_log,'w')
	call([pipeline_cmd,sofia_file],stdout=logf,stderr=logf)
	logf.close()

