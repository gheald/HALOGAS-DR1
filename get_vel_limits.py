#!/usr/bin/env python

from astropy.io import fits
import sys
from numpy import arange

h = fits.open(sys.argv[1])
hd = h[0].header
vaxis = (1.+arange(hd['NAXIS3'])-hd['CRPIX3'])*hd['CDELT3']+hd['CRVAL3']
vaxis /= 1000.
print min(vaxis),max(vaxis)

