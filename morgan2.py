from __future__ import print_function

from glob import glob
import os

import numpy as np

from astrometry.util.fits import *

# DR10
sdssbase = '/clusterfs/riemann/raid006/dr10/boss/photoObj/301'
wisebase = 'sdss-dr10d-pobj/301'

wisecols = '''
has_wise_phot ra dec
pointsource treated_as_pointsource
w1_nanomaggies w1_nanomaggies_ivar w2_nanomaggies w2_nanomaggies_ivar
w3_nanomaggies w3_nanomaggies_ivar w4_nanomaggies w4_nanomaggies_ivar
'''.split()

sdsscols = '''
cmodelflux psfflux cmodelflux_ivar psfflux_ivar
'''.split()

for path,dirnames,filenames in os.walk(wisebase):
    print('Path', path)
    dirnames.sort()
    filenames.sort()
    WW = []
    for fn in filenames:
        wfn = os.path.join(path, fn)
        print('Reading', wfn)
        W = fits_table(wfn, columns=wisecols)
        I = np.flatnonzero(W.has_wise_phot)
        W.cut(I)
        sfn = wfn.replace(wisebase, sdssbase).replace('photoWiseForced', 'photoObj')
        print('Reading', sfn)
        S = fits_table(sfn, columns=sdsscols, rows=I)

        psf = np.logical_or(W.pointsource, W.treated_as_pointsource)
        W.sdss_flux = np.zeros((len(W), 5), np.float32)
        W.sdss_flux_ivar = np.zeros((len(W), 5), np.float32)

        W.sdss_flux[psf,:] = S.psfflux[psf,:]
        W.sdss_flux_ivar[psf,:] = S.psfflux_ivar[psf,:]
        gal = np.logical_not(psf)
        W.sdss_flux[gal,:] = S.cmodelflux[gal,:]
        W.sdss_flux_ivar[gal,:] = S.cmodelflux_ivar[gal,:]

        W.pointsource = psf
        W.delete_column('treated_as_pointsource')
        W.delete_column('has_wise_phot')

        WW.append(W)

    if len(WW):
        W = merge_tables(WW)
        outfn = path.replace(wisebase, '')
        outfn = os.path.join('data/sdss-dr10d-morgan', 'phot' + outfn.replace('/','-') + '.fits')
        W.writeto(outfn)
        del W
        print('Wrote', outfn)
    
