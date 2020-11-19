#!/usr/bin/env python3

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

from astropy.io import fits
from astropy.wcs import WCS
import os

os.system("rm dummy-image.fits dummy-dirty.fits")
# -channel-range 0 1 ensures the reference frequency is from the first channel.
os.system("wsclean -size 512 512 -scale 0.01 -channel-range 0 1 -name dummy tDDECal.MS")

sources = {
    "radec": ( 400, 64 ),
    "ra": ( 400, 256 ),
    "dec": ( 256, 64 ),
    "center": ( 256, 256 )
    }
brightness = {
    "radec": 10,
    "ra": 20,
    "dec": 20,
    "center": 10
    }
term_brightness = { 0:10, 1:20000, 2:30000 }

fits_files=[]

hdu = fits.open("dummy-image.fits")[0]

def write_fits(name):
    filename = name + "-model.fits"
    os.system("rm -rf " + filename)
    hdu.writeto(filename)
    fits_files.append(filename)

# Generate foursources.fits, which has all four sources.
hdu.data *= 0

for source in sources:
    x, y = sources[source]
    hdu.data[0, 0, y, x] = brightness[source]

write_fits("foursources")

# Generate files with a single source only.
for source in sources:
    hdu.data *= 0
    x, y = sources[source]
    hdu.data[0, 0, y, x] = brightness[source]
    write_fits(source)

# Generate files for testing polynomial-terms.
# They have a single pixel at the center position, with different values.
for term in term_brightness:
    hdu.data *= 0
    hdu.data[0, 0, 256, 256] = term_brightness[term]
    hdu.header
    write_fits("term" + str(term))

os.system("tar cfj resources/idg-fits-sources.tbz2 " + " ".join(fits_files))
os.system("rm dummy-image.fits dummy-dirty.fits " + " ".join(fits_files))
