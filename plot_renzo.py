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
no_zoom = ['NGC2403', 'NGC4631', 'NGC5055']
double_zoom = ['NGC0949', 'NGC4062', 'NGC4274', 'NGC4448', 'NGC5023', 'NGC5229', 'UGC2082', 'UGC4278', 'UGC7774']

mom_dir = 'MOMENTS/'
plot_dir = 'PLOTS/'
optical_file = 'optical-images.dat'
cont_file = 'continuum-images.dat'
pol_file = 'pol-images.dat'
vals_file = 'galaxy-values.dat'
renzo_file = 'renzo-values.dat'

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

pol = {}
polgals = []
polvecs = []
plot_pol_products = False
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

renzo = {}
for line in open(renzo_file):
	if line[0] == '#': continue
	sline = line.split()
	galname = sline[0]
	renzo[galname] = {}
	renzo[galname]['W50'] = float(sline[1])
	renzo[galname]['W20'] = float(sline[2])
	renzo[galname]['CUBESTD'] = float(sline[3])*1.e-3 # Jy/beam

def getvelaxis(h):
	v = (np.arange(h['NAXIS3'])+1.-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
	return v/1000. # in km/s

def get_renzo_contours(gal, r, vsys, width):
	galcube = gal+'-LR-cube.fits'
	galh = fits.open(galcube)
	vax = getvelaxis(galh[0].header)
	galh.close()
	halfwidth = r[width]/2.
	vels = np.linspace(vsys-halfwidth, vsys+halfwidth, num=7)
	vp = np.zeros(vels.shape)
	for i in range(len(vels)):
		mv = np.abs(vax-vels[i])
		pix = np.where(mv==np.min(mv))
		vp[i] = pix[0]
	print(vels)
	print(vp)
	return vels, vp

for galaxy in galaxies:
	if 'not_specified' in cont[galaxy]: continue
	if gvals[galaxy]['INCL'] < 80.: continue
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
	print("Opening LR column density map")
	h_lr = fits.open(mom_dir+galaxy+'-LR_coldens.fits')
	hdr_lr = h_lr[0].header
	print("Establishing figure ...")
	fig = plt.figure(figsize=(14.,14.75))
	o = aplpy.FITSFigure(opt[galaxy][0],figure=fig,subplot=[0.1,0.1,0.4,0.4])
	o.recenter(hdr_lr['CRVAL1'],hdr_lr['CRVAL2'],radius=0.25*zoomscale)
	print("Plotting optical grayscale")
	o.show_grayscale(vmin=opt[galaxy][1][0],vmax=opt[galaxy][1][1],invert=True,stretch='arcsinh',aspect='auto')
	o.show_contour(mom_dir+galaxy+'-LR_coldens.fits',levels=np.array([1.e19]),colors='k')
	rclev = np.array([5.*renzo[galaxy]['CUBESTD']])
	vels, vpix = get_renzo_contours(galaxy, renzo[galaxy], gvals[galaxy]['VSYS'], 'W50')
	hcube = fits.open(galaxy+'-LR-cube.fits')
	cubedata = np.copy(hcube[0].data)
	hcube[0].data = cubedata[0,:,:]
	cubehdr = hcube[0].header
	del cubehdr['CRVAL3']
	del cubehdr['CDELT3']
	del cubehdr['CRPIX3']
	del cubehdr['CTYPE3']
	hcube[0].header = cubehdr
	for i in range(len(vels)):
		print('Plotting renzo contour at %f (velocity axis pixel %d)'%(vels[i],vpix[i]))
		# Because aplpy is broken and apparently won't be fixed:
		slicedata = cubedata[int(vpix[i]),:,:]
		hcube[0].data = slicedata
		hcube.writeto('aplpy_tmp.fits', overwrite=True)
		place = ((vels[i]-vels.min())/(vels.max()-vels.min()))*rclev[0]
		print(place, rclev[0]+place)
		o.show_contour('aplpy_tmp.fits', levels=rclev, vmin=place, vmax=rclev[0]+place, cmap='rainbow_r')
	o.tick_labels.set_xformat('hh:mm:ss')
	o.tick_labels.set_yformat('dd:mm')
	o.tick_labels.set_font(size='x-large')
	o.axis_labels.set_font(size='x-large')
	o.add_grid()
	o.grid.set_color('black')
	o.grid.set_alpha(0.2)
	o.add_beam(major=hdr_lr['BMAJ'], minor=hdr_lr['BMIN'], angle=hdr_lr['BPA'])
	o.beam.set_color('black')
	o.beam.set_frame(True)
	o.beam.set_corner('bottom right')
	o.add_label(0.05, 0.95, galaxy.replace('C0','C').replace('C','C '),
		relative=True, ha='left', va='center', weight='semibold',
		bbox=dict(boxstyle='round', ec='black', fc='white'),
		size='xx-large', zorder=200)
	print("Saving figure ...")
	plt.savefig(plot_dir+galaxy+'-renzo.pdf',bbox_inches='tight',dpi=80)
	#plt.show()
	print("Closing and cleaning up ...")
	h_lr.close()
	plt.close(fig)
	print("Done with "+galaxy)

