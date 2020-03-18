#!/usr/bin/env python

import sys, os
import matplotlib
if 'DISPLAY' not in os.environ: matplotlib.use('agg')
import aplpy
import numpy as np
import glob
import matplotlib.pylab as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits

cubes = sorted(glob.glob('*-cube.fits'))
if len(sys.argv) > 1:
	galaxies = sys.argv[1:]
else:
	galaxies = []
	for cube in cubes:
		name = cube[:7]
		if name not in galaxies:
			galaxies.append(name)

plot_dir = 'PLOTS/'
vals_file = 'galaxy-values.dat'

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

for galaxy in galaxies:
    print galaxy
    # Figure out limits since axhline does not work
    h = fits.open(galaxy+'-PV.fits')
    hdr = h[0].header
    xmin = ((1.-hdr['CRPIX1'])*hdr['CDELT1'])+hdr['CRVAL1']
    xmax = ((hdr['NAXIS1']-hdr['CRPIX1'])*hdr['CDELT1'])+hdr['CRVAL1']
    ymin = ((1.-hdr['CRPIX2'])*hdr['CDELT2'])+hdr['CRVAL2']
    ymax = ((hdr['NAXIS2']-hdr['CRPIX2'])*hdr['CDELT2'])+hdr['CRVAL2']
    hline = np.array([[xmin, xmax], [gvals[galaxy]['VSYS'], gvals[galaxy]['VSYS']]])
    print hline
    vline = np.array([[0., 0.], [ymin, ymax]])
    print vline
    print "Establishing figure ..."
    fig = plt.figure(figsize=(14.,5.))
    f = aplpy.FITSFigure(galaxy+'-PV.fits', figure=fig)#, subplot=[0.1,0.05,0.8,0.9])
    f.show_grayscale(vmin=-0.0004,vmax=0.001,invert=True,aspect='auto')
    f.show_contour(galaxy+'-PV.fits',levels=8.*0.0002*2**np.arange(10.),colors='white',alpha=0.8)
    f.show_contour(galaxy+'-PV.fits',levels=2.*0.0002*2**np.arange(2.),colors='black',alpha=0.8)
    print 'SHOWING LINES'
    f.show_lines([hline, vline], color='gray', alpha=0.8)
    plt.savefig(plot_dir+galaxy+'-pv.pdf',bbox_inches='tight',dpi=200)
    plt.close(fig)
    print 'Done with',galaxy
