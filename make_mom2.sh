#!/bin/bash

for f in `ls MOMENTS/*nrch.fits | awk '{split($0,a,"_");print a[1]}'`; do
	echo $f
	\rm -rf ${f}_mom2m.fits
	cubename=`echo $f | awk '{split($0,a,"/");print a[2]}'`-cube.fits
	maskname=${f}_mask.fits
	fits op=xyin in=${f}_mom1m.fits out=mom1m
	fits op=xyin in=${cubename} out=cube
	fits op=xyin in=${maskname} out=cubemask
	maths exp=cube mask=cubemask.gt.0 out=cubem
	imblr in=cubem out=cubemz
	moment mom=2 in=cubemz out=mom2 rngmsk=true 
	maths exp=mom2 mask=mom2.and.mom1m out=mom2m
	fits op=xyout in=mom2m out=${f}_mom2m.fits
	\rm -rf mom1m mom2 mom2m cube cubemask cubem cubemz
done

