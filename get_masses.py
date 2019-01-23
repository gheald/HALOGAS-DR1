#!/usr/bin/env python

from subprocess import call
from astropy.io import fits
import numpy as np
import glob
import astropy.units as u

mom_dir = 'MOMENTS/'

def get_galaxies():
	cubelist = sorted(glob.glob('*-cube.fits'))
	galaxies = []
	for cube in cubelist:
		galaxy = cube[:7]
		if galaxy not in galaxies:
			galaxies.append(galaxy)
	return galaxies

galaxies = get_galaxies()

#cm_per_Mpc = 3.086e24
hydrogen_mass = 1.6737236e-24 * u.g

distances = {
'NGC0672': 7.6, 'NGC0891': 9.2, 'NGC0925': 9.1, 'NGC0949': 11.3,
'UGC2082': 14.4, 'NGC1003': 11.6, 'NGC2403': 3.2, 'UGC4278': 13.6,
'NGC2541': 12.0, 'NGC3198': 14.5, 'NGC4062': 16.9, 'NGC4244': 4.4,
'NGC4258': 7.6, 'NGC4274': 19.4, 'NGC4414': 17.8, 'NGC4448': 9.7,
'NGC4559': 7.9, 'NGC4565': 10.8, 'UGC7774': 24.4, 'NGC4631': 7.6,
'NGC5023': 6.6, 'NGC5055': 8.5, 'NGC5229': 5.1, 'NGC5585': 8.7 }

jout = open('junk.txt','w')

print 'Galaxy & Distance & MHI & W50 & W20 \\\\'

for galaxy in galaxies:
	distance = distances[galaxy] * u.Mpc
	cf = mom_dir+galaxy+'-LR_coldens.fits'
	call(['fits','op=xyin','in=%s'%cf,'out=coldens'],stdout=jout,stderr=jout)
	for line in open(mom_dir+galaxy+'.region'):
		reg = line.rstrip()
		call(['immask','in=coldens','region=%s'%reg,'logic=AND'],stdout=jout,stderr=jout)
	call(['imblr','in=coldens','out=coldensz'],stdout=jout,stderr=jout)
	call(['fits','op=xyout','in=coldensz','out=coldensz.fits'],stdout=jout,stderr=jout)
	hdu = fits.open('coldensz.fits')
	hdr = hdu[0].header
	pix_x = abs(hdr['CDELT1'])*np.pi/180.*u.m/u.m
	pix_y = hdr['CDELT2']*np.pi/180.*u.m/u.m
	data = hdu[0].data
	total_coldens = np.sum(data)/((u.cm)**2)
	total_atoms = total_coldens * pix_x * pix_y * distance**2
	total_hi_mass = total_atoms * hydrogen_mass
	w50 = w20 = 0.
	for line in open('globprof-log.txt'):
		if 'W50' in line and galaxy in line:
			sline = line.split()
			w50 = float(sline[3][:-1])
			w20 = float(sline[6])
	print galaxy, '&', '%.1f'%distance.value, '&', '%.1f'%(total_hi_mass.to('solMass').value/1.e8), '&', '%.1f'%w50, '&', '%.1f'%w20, '\\\\'
	call(['rm','-rf','coldens','coldensz','coldensz.fits'],stdout=jout,stderr=jout)

jout.close()

