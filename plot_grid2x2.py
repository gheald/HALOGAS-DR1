#!/usr/bin/env python

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
no_zoom = ['NGC2403', 'NGC4274', 'UGC4278', 'NGC4631', 'NGC4258', 'NGC5055']
double_zoom = ['NGC0949', 'NGC4062', 'NGC4274', 'NGC4448', 'NGC5023', 'NGC5229', 'UGC2082', 'UGC7774']

mom_dir = 'MOMENTS/'
plot_dir = 'PLOTS/'
optical_file = 'optical-images.dat'
cont_file = 'continuum-images.dat'
vals_file = 'galaxy-values.dat'

opt = {}
for line in open(optical_file):
	sline = line.split()
	galname = sline[0]
	optpath = '../HALOGAS-OPTICAL/'+sline[1]
	lims = np.asarray(sline[2].split(','),dtype=float)
	opt[galname] = [optpath, lims]

cont = {}
for line in open(cont_file):
	sline = line.split()
	galname = sline[0]
	cont[galname] = '../HALOGAS-CONTINUUM/'+sline[1]

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

print opt
print cont
print gvals

for galaxy in galaxies:
	print galaxy
	if galaxy in no_zoom:
		zoomscale = 2.
		print "NO ZOOM ENABLED"
	elif galaxy in double_zoom:
		zoomscale = 0.5
		print "DOUBLE ZOOM ENABLED"
	else:
		zoomscale = 1.
		print "SINGLE ZOOM ENABLED"
	print "Opening HR column density map"
	h_hr = fits.open(mom_dir+galaxy+'-HR_coldens.fits')
	hdr_hr = h_hr[0].header
	print "Opening LR column density map"
	h_lr = fits.open(mom_dir+galaxy+'-LR_coldens.fits')
	hdr_lr = h_lr[0].header
	print "Establishing figure ..."
	fig = plt.figure(figsize=(14.,14.75))
	f = aplpy.FITSFigure(mom_dir+galaxy+'-HR_coldens_toplot.fits',figure=fig,subplot=[0.1,0.5,0.4,0.4])
	f.recenter(hdr_hr['CRVAL1'],hdr_hr['CRVAL2'],radius=0.25*zoomscale)
	print "Plotting column density grayscale"
	f.show_grayscale(vmin=0.01,vmax=2.5,invert=True)
	f.show_ellipses(gvals[galaxy]['CENTER'].ra.value, gvals[galaxy]['CENTER'].dec.value,
		gvals[galaxy]['DIAMETER']/3600., 1./3600., angle=gvals[galaxy]['PA'],
		edgecolor='black', zorder=100, alpha=0.8, facecolor='none')
	f.show_contour(mom_dir+galaxy+'-LR_coldens.fits',levels=1.e19*2.**np.arange(8),colors='k')
	f.tick_labels.set_xformat('hh:mm:ss')
	f.tick_labels.set_yformat('dd:mm')
	f.axis_labels.hide_x()
	f.tick_labels.hide_x()
	f.add_colorbar()
	f.colorbar.set_location('top')
	f.colorbar.set_axis_label_text('Column density ($\mathregular{10^{21}}$ atoms cm$\mathregular{^{-2}}$)')
	f.add_grid()
	f.grid.set_color('black')
	f.grid.set_alpha(0.2)
	f.add_beam()
	f.add_beam()
	f.beam[0].set_color('black')
	f.beam[0].set_alpha(0.5)
	f.beam[0].set_frame(True)
	f.beam[1].set_color('black')
	f.beam[1].set_alpha(0.5)
	f.beam[1].set_major(hdr_lr['BMAJ'])
	f.beam[1].set_minor(hdr_lr['BMIN'])
	f.beam[1].set_angle(hdr_lr['BPA'])
	o = aplpy.FITSFigure(opt[galaxy][0],figure=fig,subplot=[0.1,0.1,0.4,0.4])
	o.recenter(hdr_hr['CRVAL1'],hdr_hr['CRVAL2'],radius=0.25*zoomscale)
	print "Plotting optical grayscale"
	o.show_grayscale(vmin=opt[galaxy][1][0],vmax=opt[galaxy][1][1],invert=True,stretch='sqrt')
	o.show_contour(mom_dir+galaxy+'-LR_coldens.fits',levels=1.e19*2.**np.arange(8),colors='k')
	o.tick_labels.set_xformat('hh:mm:ss')
	o.tick_labels.set_yformat('dd:mm')
	o.add_grid()
	o.grid.set_color('black')
	o.grid.set_alpha(0.2)
	g = aplpy.FITSFigure(mom_dir+galaxy+'-LR_mom1m_toplot.fits',figure=fig,subplot=[0.5,0.5,0.4,0.4])
	g.recenter(hdr_hr['CRVAL1'],hdr_hr['CRVAL2'],radius=0.25*zoomscale)
	print "Plotting velocity field colormap"
	g.show_colorscale()
	g.show_ellipses(gvals[galaxy]['CENTER'].ra.value, gvals[galaxy]['CENTER'].dec.value,
		gvals[galaxy]['DIAMETER']/3600., 1./3600., angle=gvals[galaxy]['PA'],
		edgecolor='black', zorder=100, alpha=0.8, facecolor='none')
	g.show_contour(mom_dir+galaxy+'-LR_mom1m_toplot.fits', levels=gvals[galaxy]['VSYS'],
		colors='white', zorder=10, alpha=1.)
	max_dvel = gvals[galaxy]['VROT']*np.sin(gvals[galaxy]['INCL']*np.pi/180.)
	dlevels = np.arange(20.,max_dvel,20.)
	g.show_contour(mom_dir+galaxy+'-LR_mom1m_toplot.fits',
		levels=gvals[galaxy]['VSYS']+dlevels,
		colors='gray', zorder=10, alpha=0.8)
	g.show_contour(mom_dir+galaxy+'-LR_mom1m_toplot.fits',
		levels=gvals[galaxy]['VSYS']-np.flip(dlevels),
		colors='gray', zorder=10, alpha=0.8)
	g.tick_labels.set_xformat('hh:mm:ss')
	g.tick_labels.set_yformat('dd:mm')
	g.tick_labels.hide_y()
	g.axis_labels.hide_y()
	g.tick_labels.hide_x()
	g.axis_labels.hide_x()
	g.add_colorbar()
	g.colorbar.set_location('top')
	g.colorbar.set_axis_label_text('Velocity (km s$\mathregular{^{-1}}$)')
	g.add_grid()
	g.grid.set_color('black')
	g.grid.set_alpha(0.2)
	g.add_beam()
	g.beam.set_color('black')
	g.beam.set_frame(True)
	c = aplpy.FITSFigure(cont[galaxy],figure=fig,subplot=[0.5,0.1,0.4,0.4])
	c.recenter(hdr_hr['CRVAL1'],hdr_hr['CRVAL2'],radius=0.25*zoomscale)
	print "Plotting continuum grayscale"
	c.show_grayscale(vmin=-0.0001,vmax=0.001,invert=True)
	c.tick_labels.set_xformat('hh:mm:ss')
	c.tick_labels.set_yformat('dd:mm')
	c.axis_labels.hide_y()
	c.tick_labels.hide_y()
	c.add_grid()
	c.grid.set_color('black')
	c.grid.set_alpha(0.2)
	print "Saving figure ..."
	plt.savefig(plot_dir+galaxy+'-grid2x2.pdf',bbox_inches='tight',dpi=200)
	#plt.show()
	print "Closing and cleaning up ..."
	h_hr.close()
	h_lr.close()
	plt.close(fig)
	print "Done with", galaxy
