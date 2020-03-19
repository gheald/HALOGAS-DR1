#!/usr/bin/env python

from astropy.io import fits
from subprocess import call
from datetime import datetime
from os.path import isfile
from os import remove
import argparse
import numpy as np

def logprint(s2p,lf):
	print >>lf, s2p
	print s2p

def removeaxis(hdu):
	header = hdu[0].header
	data = np.squeeze(hdu[0].data)
	hdu[0].data = data
	#del(header['NAXIS4'])
	del(header['CTYPE4'])
	del(header['CRVAL4'])
	del(header['CRPIX4'])
	del(header['CDELT4'])
	header['NAXIS']=3
	hdu[0].header=header

def main(args):
	logf = open(args.log_file,'w',1)
	logprint('Script run on %s'%str(datetime.now()),logf)
	for line in open(args.location_file):
		sline = line.split()
		cube_id = sline[0]
		if args.galaxy not in cube_id and args.galaxy != '__ALL__':
			continue
		cubename = '%s-cube.fits'%cube_id
		logprint('\n####\n',logf)
		if isfile(cubename):
			if not args.overwrite:
				logprint('Cube %s already exists, skipping'%cubename,logf)
				continue
			else:
				logprint('Cube %s already exists, replacing'%cubename,logf)
				remove(cubename)
		mir_loc = sline[1]
		logprint('Converting %s to fits format'%mir_loc,logf)
		logprint('Output FITS name is %s'%cubename,logf)
		logprint('\n####\n',logf)
		call(['fits','op=xyout','in=%s'%mir_loc,'out=%s'%cubename],stdout=logf,stderr=logf)
		# Check the number of axes
		hdu = fits.open(cubename,mode='update')
		header = hdu[0].header
		naxis = header['NAXIS']
		logprint('Cube %s has %d axes'%(cubename,naxis),logf)
		# If there are four, remove the Stokes axis
		if naxis == 4:
			logprint('Removing Stokes axis from %s'%cubename,logf)
			removeaxis(hdu)
		# Add header entries that say it's HALOGAS data and give ADS reference(s)
		if 'NGC0891' in cube_id:
			ads_ref = 'http://adsabs.harvard.edu/abs/2007AJ....134.1019O'
			header.insert('EPOCH',('COMMENT',ads_ref),after=True)
			header.insert('EPOCH',('COMMENT','Data for NGC 891 adopted by HALOGAS described by:'),after=True)
		elif 'NGC2403' in cube_id:
			ads_ref = 'http://adsabs.harvard.edu/abs/2002AJ....123.3124F'
			header.insert('EPOCH',('COMMENT',ads_ref),after=True)
			header.insert('EPOCH',('COMMENT','Data for NGC 2403 adopted by HALOGAS described by:'),after=True)
		ads_ref = 'http://adsabs.harvard.edu/abs/2011A%26A...526A.118H'
		header.insert('EPOCH',('COMMENT',ads_ref),after=True)
		header.insert('EPOCH',('COMMENT','Please include this reference in any publication:'),after=True)
		header.insert('EPOCH',('COMMENT','Data Release 1 (v2): https://doi.org/10.5281/zenodo.3715549'),after=True)
		header.insert('EPOCH',('COMMENT','Project website: http://www.astron.nl/halogas'),after=True)
		header.insert('EPOCH',('COMMENT','Data released as part of HALOGAS DR1'),after=True)
		# Flush and close the FITS file
		hdu[0].header = header
		hdu.flush()
		hdu.close()
	# Close the log file
	logf.close()


ap = argparse.ArgumentParser()
ap.add_argument('--galaxy','-g',help='Name of galaxy to process (which can either be one galaxy or all of them) [default all]',default='__ALL__')
ap.add_argument('--location_file','-l',help='Name of file containing cube locations [default cube-locations.dat]',default='cube-locations.dat')
ap.add_argument('--log_file','-f',help='Name of output log file [default log.txt]',default='log.txt')
ap.add_argument('--overwrite','-o',help='Overwrite any existing FITS files that already exist? [default False]',default=False,action='store_true')
args = ap.parse_args()
main(args)

