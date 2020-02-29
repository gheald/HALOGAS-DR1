#!/usr/bin/env python

from astropy.io import fits
import aplpy
import numpy as np
import glob
import matplotlib.pylab as plt

cubes = sorted(glob.glob('*-cube.fits'))
galaxies = []
for cube in cubes:
	name = cube[:7]
	if name not in galaxies:
		galaxies.append(name)
#print galaxies
no_zoom = ['NGC0891', 'NGC0925', 'NGC2403', 'NGC4274', 'UGC4278', 'NGC4631', 'NGC4258', 'NGC5055']

mom_dir = 'MOMENTS/'
plot_dir = 'PLOTS/'
optical_file = 'optical-images.dat'
cont_file = 'continuum-images.dat'

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

print opt
print cont

for galaxy in galaxies:
	if galaxy in no_zoom:
		zoomscale = 2.
	else:
		zoomscale = 1.
	h_hr = fits.open(mom_dir+galaxy+'-HR_coldens.fits')
	hdr_hr = h_hr[0].header
	h_lr = fits.open(mom_dir+galaxy+'-LR_coldens.fits')
	hdr_lr = h_lr[0].header
	fig = plt.figure(figsize=(14.,14.75))
	f = aplpy.FITSFigure(mom_dir+galaxy+'-HR_coldens_toplot.fits',figure=fig,subplot=[0.1,0.5,0.4,0.4])
	f.recenter(hdr_hr['CRVAL1'],hdr_hr['CRVAL2'],radius=0.25*zoomscale)
	f.show_grayscale(vmin=0.01,vmax=2.5,invert=True)
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
	o.show_grayscale(vmin=opt[galaxy][1][0],vmax=opt[galaxy][1][1],invert=True,stretch='sqrt')
	o.show_contour(mom_dir+galaxy+'-LR_coldens.fits',levels=1.e19*2.**np.arange(8),colors='k')
	o.tick_labels.set_xformat('hh:mm:ss')
	o.tick_labels.set_yformat('dd:mm')
	o.add_grid()
	o.grid.set_color('black')
        o.grid.set_alpha(0.2)
	g = aplpy.FITSFigure(mom_dir+galaxy+'-LR_mom1m_toplot.fits',figure=fig,subplot=[0.5,0.5,0.4,0.4])
	g.recenter(hdr_hr['CRVAL1'],hdr_hr['CRVAL2'],radius=0.25*zoomscale)
	g.show_colorscale()
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
	c.show_grayscale(vmin=-0.0001,vmax=0.001,invert=True)
	c.tick_labels.set_xformat('hh:mm:ss')
	c.tick_labels.set_yformat('dd:mm')
	c.axis_labels.hide_y()
	c.tick_labels.hide_y()
	c.add_grid()
	c.grid.set_color('black')
        c.grid.set_alpha(0.2)
	plt.savefig(plot_dir+galaxy+'-grid2x2.png',bbox_inches='tight',dpi=200)
	#plt.show()
	h_hr.close()
	h_lr.close()
	plt.close(fig)

