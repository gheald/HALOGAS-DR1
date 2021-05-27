#!/usr/bin/env python3

import sys, os
import matplotlib
if 'DISPLAY' not in os.environ: matplotlib.use('agg')
from astropy.io import fits
import aplpy
import numpy as np
import glob
import matplotlib.pylab as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval

cubes = sorted(glob.glob('*-cube.fits'))
if len(sys.argv) > 1:
	galaxies = sys.argv[1:]
else:
	galaxies = []
	for cube in cubes:
		name = cube[:7]
		if name not in galaxies:
			galaxies.append(name)
#print galaxies
no_halpha = ['NGC4244', 'UGC2082', 'NGC0949', 'NGC4448', 'UGC7774']
no_zoom = []
double_zoom = ['NGC0672', 'NGC0891', 'NGC0925', 'NGC1003', 'NGC0949', 'NGC2541', 'NGC3198', 'NGC4062', 'NGC4258', 'NGC4274', 'NGC4414', 'NGC4448', 'NGC4559', 'NGC4565', 'NGC4631', 'NGC5023', 'NGC5055', 'NGC5229', 'NGC5585', 'UGC2082', 'UGC4278', 'UGC7774']

mom_dir = 'MOMENTS/'
plot_dir = 'PLOTS/'
optical_file = 'halpha-images.dat'
cont_file = 'continuum-images.dat'
pol_file = 'pol-images.dat'
vals_file = 'galaxy-values.dat'

opt = {}
zs = ZScaleInterval(contrast=1.)
for line in open(optical_file):
	sline = line.split()
	galname = sline[0]
	optpath = '../HALOGAS-OPTICAL/'+sline[1]
	hzs = fits.open(optpath)
	#lims = np.asarray(sline[2].split(','),dtype=float)
	lims = zs.get_limits(hzs[0].data)
	lims = (lims[0], -3.*lims[0])
	print((galname,lims))
	opt[galname] = [optpath, lims]

cont = {}
for line in open(cont_file):
	sline = line.split()
	galname = sline[0]
	noise = float(sline[2])
	cont[galname] = ['../HALOGAS-CONTINUUM/'+sline[1], noise]

pol = {}
polgals = []
polvecs = []
for line in open(pol_file):
	sline = line.split()
	galname = sline[0]
	polgals.append(galname)
	polpath = '../HALOGAS-CONTINUUM/'+sline[1]
	polclip = float(sline[2])
	for vline in open('../HALOGAS-CONTINUUM/'+sline[3]):
		if vline[0] == 'L':
			vsline = vline.split()
			x1 = float(vsline[2])
			y1 = float(vsline[3])
			x2 = float(vsline[5])
			y2 = float(vsline[6])
			tline = np.array([[x1,x2],[y1,y2]])
			polvecs.append(tline)
	pol[galname] = [polpath, polclip, polvecs]

gvals = {}
for line in open(vals_file):
	if line[0] == '#': continue
	sline = line.split()
	galname = sline[0]
	gvals[galname] = {}
	gvals[galname]['CENTER'] = SkyCoord(sline[1], sline[2], frame='fk5')
	gvals[galname]['PA'] = float(sline[3]) - 90. # so that PA=0 points north
	gvals[galname]['DIAMETER'] = float(sline[4])
	gvals[galname]['VSYS'] = float(sline[5])
	gvals[galname]['VROT'] = float(sline[6])
	gvals[galname]['INCL'] = float(sline[7])

#print opt
#print cont
#print gvals

for galaxy in galaxies:
	if 'not_specified' in cont[galaxy]: continue
	if galaxy in no_halpha: continue
	print(galaxy)
	if galaxy in no_zoom:
		zoomscale = 2.
		print("NO ZOOM ENABLED")
	elif galaxy in double_zoom:
		zoomscale = 0.5
		print("DOUBLE ZOOM ENABLED")
	else:
		zoomscale = 1.
		print("SINGLE ZOOM ENABLED")
	print("Establishing figure ...")
	fig = plt.figure(figsize=(14.,14.75))
	hdr = fits.open(opt[galaxy][0])[0].header
	hdr_c = fits.open(cont[galaxy][0])[0].header
	o = aplpy.FITSFigure(opt[galaxy][0],figure=fig,subplot=[0.1,0.1,0.9,0.9])
	#o.recenter(hdr['CRVAL1'],hdr['CRVAL2'],radius=0.25*zoomscale)
	o.recenter(gvals[galaxy]['CENTER'].ra.value,gvals[galaxy]['CENTER'].dec.value,radius=0.25*zoomscale)
	print("Plotting optical grayscale")
	#o.show_grayscale(vmin=opt[galaxy][1][0],vmax=opt[galaxy][1][1],invert=True,stretch='arcsinh',aspect='auto')
	o.show_colorscale(vmin=opt[galaxy][1][0],vmax=opt[galaxy][1][1],aspect='auto',cmap='Blues')
	#o.show_contour(mom_dir+galaxy+'-LR_coldens.fits',levels=1.e19*2.**np.arange(8),colors='k')
	o.show_contour(cont[galaxy][0],levels=3.*cont[galaxy][1]*2.**np.arange(8),colors='k')
	o.show_contour(cont[galaxy][0],levels=np.array([-3.*cont[galaxy][1]]),colors='k',linestyles='--')
	o.tick_labels.set_xformat('hh:mm:ss')
	o.tick_labels.set_yformat('dd:mm')
	o.tick_labels.set_font(size=24)
	o.axis_labels.set_font(size=24)
	o.add_grid()
	o.grid.set_color('black')
	o.grid.set_alpha(0.2)
	o.add_beam(major=hdr_c['BMAJ'], minor=hdr_c['BMIN'], angle=hdr_c['BPA'])
	o.beam.set_color('black')
	o.beam.set_frame(True)
	o.beam.set_corner('bottom right')
	o.add_label(0.05, 0.95, galaxy.replace('C0','C').replace('C','C '),
		relative=True, ha='left', va='center', weight='semibold',
		bbox=dict(boxstyle='round', ec='black', fc='white'),
		size=28, zorder=200)
	if galaxy in polgals:
		print('PLOTTING POLARIZATION VALUES')
		o.show_contour(data=pol[galaxy][0], levels=np.array([pol[galaxy][1]]), linewidths=1.5, colors='green')
		o.show_lines(pol[galaxy][2], color='red', linewidths=2)
	print("Saving figure ...")
	plt.savefig(plot_dir+galaxy+'-halphacont.pdf',bbox_inches='tight',dpi=80)
	#plt.show()
	print("Closing and cleaning up ...")
	plt.close(fig)
	print("Done with "+galaxy)

