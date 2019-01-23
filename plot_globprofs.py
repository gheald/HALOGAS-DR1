#!/usr/bin/env python

from astropy.io import fits
import glob
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from tqdm import tqdm
import argparse
from os.path import isfile

def get_intercept(v, p, i, l):
	va = v[i-1:i+2]
	pa = p[i-1:i+2]
	pf = np.polyfit(va,pa,1)
	return (1./pf[0])*(l-pf[1])

def get_widths(v, p, fit=False):
	maxval = max(p)
	i50 = np.where(p>=0.5*maxval)[0]
	i20 = np.where(p>=0.2*maxval)[0]
	if fit:
		e50_min = get_intercept(v, p, min(i50), 0.5*maxval)
		e50_max = get_intercept(v, p, max(i50), 0.5*maxval)
		e20_min = get_intercept(v, p, min(i20), 0.2*maxval)
		e20_max = get_intercept(v, p, max(i20), 0.2*maxval)
		w50 = e50_max - e50_min
		w20 = e20_max - e20_min
	else:
		v50 = v[i50]
		v20 = v[i20]
		w50 = max(v50) - min(v50)
		w20 = max(v20) - min(v20)
	return w50, w20

def get_galaxies():
	cubelist = sorted(glob.glob('*-cube.fits'))
	galaxies = []
	for cube in cubelist:
		galaxy = cube[:7]
		if galaxy not in galaxies: 
			galaxies.append(galaxy)
	return galaxies

def main(args):

	mom_dir = 'MOMENTS/'
	
	if args.galaxies == '__ALL__':
		galaxies = get_galaxies()
	else:
		galaxies = args.galaxies.split(',')

	logf = open('globprof-log.txt','w')
	print >>logf, args

	for galaxy in tqdm(galaxies): 
		if not isfile(galaxy+'-gprof-lr.txt') or not args.skip:
			tqdm.write('%s: Opening LR cubes'%galaxy)
			hdu_cube_lr = fits.open(galaxy+'-LR-cube.fits')
			hdu_mask_lr = fits.open(mom_dir+galaxy+'-LR_mask.fits')
			data_lr = hdu_cube_lr[0].data
			mask_lr = hdu_mask_lr[0].data
			data_lr[mask_lr==0.] = 0.
			hdu_cube_lr[0].data = data_lr
			tqdm.write('%s: Writing masked LR cube'%galaxy)
			hdu_cube_lr.writeto('tmp-lr.fits')
			tqdm.write('%s: Converting masked LR cube to miriad'%galaxy)
			call(['fits','op=xyin','in=tmp-lr.fits','out=lr'], stderr=logf, stdout=logf)
			if isfile(mom_dir+galaxy+'.region'):
				tqdm.write('%s: Applying region mask'%galaxy)
				for line in open(mom_dir+galaxy+'.region'):
					reg = line.rstrip()
					call(['immask','in=lr','region=%s'%reg,'logic=AND'], stderr=logf, stdout=logf)
			tqdm.write('%s: Computing LR global profile'%galaxy)
			call(['imspec','in=lr','axes=ra,dec','device=/null','options=list,eformat,guaranteespaces,noheader','log=%s-gprof-lr.txt'%galaxy], stderr=logf, stdout=logf)
		if not isfile(galaxy+'-gprof-hr.txt') or not args.skip:
			tqdm.write('%s: Opening HR cubes'%galaxy)
			hdu_cube_hr = fits.open(galaxy+'-HR-cube.fits')
			hdu_mask_hr = fits.open(mom_dir+galaxy+'-HR_mask.fits')
			data_hr = hdu_cube_hr[0].data
			mask_hr = hdu_mask_hr[0].data
			data_hr[mask_hr==0.] = 0.
			hdu_cube_hr[0].data = data_hr
			tqdm.write('%s: Writing masked HR cube'%galaxy)
			hdu_cube_hr.writeto('tmp-hr.fits')
			tqdm.write('%s: Converting masked HR cube to miriad'%galaxy)
			call(['fits','op=xyin','in=tmp-hr.fits','out=hr'], stderr=logf, stdout=logf)
			if isfile(mom_dir+galaxy+'.region'):
				tqdm.write('%s: Applying region mask'%galaxy)
				for line in open(mom_dir+galaxy+'.region'):
					reg = line.rstrip()
					call(['immask','in=hr','region=%s'%reg,'logic=AND'], stderr=logf, stdout=logf)
			tqdm.write('%s: Computing HR global profile'%galaxy)
			call(['imspec','in=hr','axes=ra,dec','device=/null','options=list,eformat,guaranteespaces,noheader','log=%s-gprof-hr.txt'%galaxy], stderr=logf, stdout=logf)
		velaxis_lr, prof_lr = np.loadtxt('%s-gprof-lr.txt'%galaxy,usecols=(1,2),unpack=True)
		velaxis_hr, prof_hr = np.loadtxt('%s-gprof-hr.txt'%galaxy,usecols=(1,2),unpack=True)
		w50, w20 = get_widths(velaxis_lr, prof_lr, fit=args.fitwidth)
		tqdm.write('%s: W50 = %f, W20 = %f'%(galaxy, w50, w20))
		print >>logf, '%s: W50 = %f, W20 = %f'%(galaxy, w50, w20)
		if not args.noplot:
			tqdm.write('%s: Producing plot'%galaxy)
			fig = plt.figure(figsize=(8.,5.))
			plt.plot(velaxis_hr,prof_hr,'k--')
			plt.plot(velaxis_lr,prof_lr,'k-')
			plt.xlim((min(velaxis_lr),max(velaxis_lr)))
			plt.ylim(bottom=0)
			plt.xlabel('Velocity (km s$\mathregular{^{-1}}$)')
			plt.ylabel('Flux density (Jy)')
			plt.savefig(galaxy+'-profs.png',dpi=120,bbox_inches='tight')
			plt.close(fig)
		call(['rm','-rf','tmp-lr.fits','tmp-hr.fits','lr','hr'])

	logf.close()

ap = argparse.ArgumentParser()
ap.add_argument('--galaxies','-g',help='Comma separated list of galaxies to deal with [default all]',default='__ALL__')
ap.add_argument('--skip','-s',help='Skip initial steps if globprof text files are present? [default False]',default=False,action='store_true')
ap.add_argument('--noplot','-n',help='Produce no plots? [default False]',default=False,action='store_true')
ap.add_argument('--fitwidth','-f',help='Fit W50 and W20 instead of simple method? [default False]',default=False,action='store_true')
args = ap.parse_args()
main(args)

