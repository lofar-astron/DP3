#!/usr/bin/env python3
from astropy.io import fits
import os

os.system("rm dummy-image.fits dummy-dirty.fits")
os.system("wsclean -size 1024 1024 -scale 0.01 -name dummy tDDECal.MS")

sources = {
    "radec": ( 800, 128 ),
    "ra": ( 800, 512 ),
    "dec": ( 512, 128 ),
    "center": ( 512, 512 )
    }
brightness = {
    "radec": 10,
    "ra": 20,
    "dec": 20,
    "center": 10
    }
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

os.system("tar cfj resources/idg-fits-sources.tbz2 " + " ".join(fits_files))
os.system("rm dummy-image.fits dummy-dirty.fits")
