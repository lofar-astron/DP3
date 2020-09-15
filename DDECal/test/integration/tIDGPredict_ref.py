#!/usr/bin/env python3
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs import WCS
import numpy as np
import os

os.system("rm dummy-image.fits dummy-dirty.fits")
os.system("wsclean -size 1024 1024 -scale 0.01 -name dummy tDDECal.MS")

directions = {
    "radec_off": ( 800, 128 ),
    "ra_off": ( 800, 512 ),
    "dec_off": ( 512, 128 ),
    "center": ( 512, 512 )
    }
brightness = {
    "radec_off": 10,
    "ra_off": 20,
    "dec_off": 20,
    "center": 10
    }

hdu = fits.open("dummy-image.fits")[0]

hdu.data *= 0

for direction in directions:
    x, y = directions[direction]
    hdu.data[0, 0, y, x] = brightness[direction]

os.system("rm -rf foursources.fits")
hdu.writeto("foursources.fits")

os.system("gzip foursources.fits")
os.system("rm dummy-image.fits dummy-dirty.fits")
