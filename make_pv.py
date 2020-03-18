#!/usr/bin/env python

from pvextractor import PathFromCenter, extract_pv_slice
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob, sys
import numpy as np

cubes = sorted(glob.glob('*-cube.fits'))
if len(sys.argv) > 1:
	galaxies = sys.argv[1:]
else:
	galaxies = []
	for cube in cubes:
		name = cube[:7]
		if name not in galaxies:
			galaxies.append(name)


vals_file = 'galaxy-values.dat'

gvals = {}
for line in open(vals_file):
	if line[0] == '#': continue
	sline = line.split()
	galname = sline[0]
	gvals[galname] = {}
	gvals[galname]['CENTER'] = SkyCoord(sline[1], sline[2], frame='fk5')
	gvals[galname]['PA'] = float(sline[3])
	gvals[galname]['DIAMETER'] = float(sline[4])
	gvals[galname]['VSYS'] = float(sline[5])
	gvals[galname]['VROT'] = float(sline[6])
	gvals[galname]['INCL'] = float(sline[7])

velrange = {'NGC0891': (100.,900.),
            'NGC4274': (600.,1200.),
            'NGC4448': (-999.,600.),
            'NGC4559': (-999.,1000.),
            'NGC4631': (-999.,850.),
            'NGC5055': (200.,800.),
            'NGC5229': (-999.,470.),
            'UGC2082': (550.,850.),
            'UGC4278': (400.,700.),
            'UGC7774': (350.,680.)}

for galaxy in galaxies:
    pv_path = PathFromCenter(center=gvals[galaxy]['CENTER'],
        length=gvals[galaxy]['DIAMETER']*u.arcsec,
        angle=gvals[galaxy]['PA']*u.deg,
        width=1.*u.arcsec)

    pv_slice = extract_pv_slice(galaxy+'-LR-cube.fits', pv_path)

    hdr = pv_slice.header

    vaxis = (np.arange(hdr['NAXIS2'])+1.-hdr['CRPIX2'])*hdr['CDELT2']+hdr['CRVAL2']
    vaxis /= 1000.
    if galaxy in velrange.keys():
        print 'REDUCING VELOCITY RANGE'
        vmin, vmax = velrange[galaxy]
        if vmin < -998.: vmin = np.min(vaxis)
        gvp = np.where(np.logical_and(vaxis<=vmax,vaxis>=vmin))
        print np.min(gvp), np.max(gvp), len(vaxis)
        data = pv_slice.data
        newdata = data[gvp,:]
        pv_slice.data = newdata
        voffset = np.min(gvp)
    else:
        voffset = 0

    # Tweak header values to prepare for plotting in the next step
    hdr['CRPIX1'] = hdr['NAXIS1']/2
    hdr['CDELT1'] *= 60.
    hdr['CDELT2'] /= 1000.
    hdr['CUNIT1'] = 'arcmin'
    hdr['CTYPE2'] = 'VELOCITY'
    hdr['CUNIT2'] = 'km/s'
    hdr['CRVAL2'] /= 1000.
    hdr['CRPIX2'] -= voffset
    pv_slice.header = hdr

    pv_slice.writeto(galaxy+'-PV.fits', overwrite=True)
