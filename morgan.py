from __future__ import print_function

from glob import glob
import os

import numpy as np

from astrometry.util.fits import *

fns = glob('sdss-dr10d-phot/phot-????????.fits')
fns.sort()
print(len(fns), 'phot files')

cols = '''
objid ra dec cmodelflux psfflux cmodelflux_ivar psfflux_ivar
pointsource treated_as_pointsource
w1_nanomaggies w1_nanomaggies_ivar w2_nanomaggies w2_nanomaggies_ivar
w3_nanomaggies w3_nanomaggies_ivar w4_nanomaggies w4_nanomaggies_ivar
'''.split()

for i,fn in enumerate(fns):
    print('Reading', i+1, 'of', len(fns), ':', fn)
    T = fits_table(fn, columns=cols)

    psf = np.logical_or(T.pointsource, T.treated_as_pointsource)
    T.sdss_flux = np.zeros((len(T), 5), np.float32)
    T.sdss_flux_ivar = np.zeros((len(T), 5), np.float32)

    T.sdss_flux[psf,:] = T.psfflux[psf,:]
    T.sdss_flux_ivar[psf,:] = T.psfflux_ivar[psf,:]
    gal = np.logical_not(psf)
    T.sdss_flux[gal,:] = T.cmodelflux[gal,:]
    T.sdss_flux_ivar[gal,:] = T.cmodelflux_ivar[gal,:]

    T.pointsource = psf
    T.delete_column('treated_as_pointsource')
    T.delete_column('cmodelflux')
    T.delete_column('cmodelflux_ivar')
    T.delete_column('psfflux')
    T.delete_column('psfflux_ivar')

    outfn = fn.replace('sdss-dr10d-phot', 'data/sdss-dr10d-morgan')
    T.writeto(outfn)
    print('Wrote', outfn)

