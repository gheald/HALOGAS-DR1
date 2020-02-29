#!/bin/bash

# This value is the LR column density below which no velocities will be plotted 
# NB: this is column density BEFORE primary beam correction!
plot_clip=1e19

for f in `ls MOMENTS/*0891*nrch.fits | awk '{split($0,a,"_");print a[1]}'`; do
	echo $f
	\rm -rf ${f}_mom[01]?*.fits ${f}_colden*.fits
	cubename=`echo $f | awk '{split($0,a,"/");print a[2]}'`-cube.fits
	vellims=`./get_vel_limits.py $cubename`
	vmin=`echo $vellims | awk '{split($0,a," ");print a[1]}'`
	vmax=`echo $vellims | awk '{split($0,a," ");print a[2]}'`
	fits op=xyin in=${f}_mom0.fits out=mom0
	fits op=xyin in=${f}_mom1.fits out=mom1
	fits op=xyin in=${f}_nrch.fits out=nrch
	maths exp=mom0 mask=nrch.gt.2.and.mom1.gt.${vmin}.and.mom1.lt.${vmax} out=mom0msk
	maths exp=mom1 mask=nrch.gt.2.and.mom1.gt.${vmin}.and.mom1.lt.${vmax} out=mom1msk
	imblr in=mom0msk out=mom0mskz
	# No imblr for mom1 because those pixels should be blank
	linmos in=mom0mskz out=mom0mskzlm
	bmaj=`gethd in=mom0msk/bmaj format=arcsec`
	bmin=`gethd in=mom0msk/bmin format=arcsec`
	maths exp=mom0mskz*1.823e18*606*1e3/${bmaj}/${bmin} out=coldens
	maths exp=mom0mskzlm*1.823e18*606*1e3/${bmaj}/${bmin} out=coldenslm
	maths exp=coldenslm/1e21 out=coldens_toplot
	maths exp=mom1msk mask=coldens.gt.${plot_clip} out=mom1msk_toplot
	fits op=xyout in=mom0mskzlm out=${f}_mom0m.fits
	fits op=xyout in=coldenslm out=${f}_coldens.fits
	fits op=xyout in=coldens_toplot out=${f}_coldens_toplot.fits
	fits op=xyout in=mom1msk out=${f}_mom1m.fits
	fits op=xyout in=mom1msk_toplot out=${f}_mom1m_toplot.fits
	\rm -rf mom0 mom1 nrch mom0msk mom0mskz mom0mskzlm coldens coldenslm coldens_toplot mom1msk mom1msk_toplot
done

