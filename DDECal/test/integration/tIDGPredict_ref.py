#!/usr/bin/env python3
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs import WCS
import numpy as np
import os

os.system("rm dummy-image.fits dummy-dirty.fits")
os.system("wsclean -size 512 512 -scale 0.01 -name dummy tDDECal.MS")

directions = {
    "radec_off": SkyCoord("01h33m40        +32d09m0       "),
    "ra_off": SkyCoord("01h29m40        +33d09m35.13240"),
    "dec_off": SkyCoord("01h37m41.299440 +32d09m0       "),
    "center": SkyCoord("01h37m41.299440 +33d09m35.13240")
    }
brightness = {
    "radec_off": 10,
    "ra_off": 20,
    "dec_off": 20,
    "center": 10
    }

hdu = fits.open("dummy-image.fits")[0]

mywcs = WCS(hdu.header)

hdu.data *= 0

for direction in directions:
    source_pixels = np.array(skycoord_to_pixel(directions[direction], mywcs))
    ra_pix, dec_pix = np.round(source_pixels)
    hdu.data[0, 0, int(dec_pix), int(ra_pix)] = brightness[direction]

hdu.writeto("foursources.fits")

os.system("gzip foursources.fits")
os.system("rm dummy-image.fits dummy-dirty.fits")
