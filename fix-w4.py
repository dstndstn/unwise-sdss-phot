#! /usr/bin/env python
import numpy as np
import os
from glob import glob
from astrometry.util.fits import *

def raiseitup(x):
    raise x

for inbase,outbase in [('sdss-dr10c-pobj', 'sdss-dr10d-pobj'),
                       ('sdss-dr10c-phot', 'sdss-dr10d-phot'),
                       ]:
    for dirpath, dirnames, fns in os.walk(inbase, followlinks=True, onerror=raiseitup):
        dirnames.sort()
        fns.sort()
        print
        print 'dirpath', dirpath
        outdir = dirpath.replace(inbase, outbase)
        print 'outdir ', outdir
        try:
            os.makedirs(outdir)
        except:
            pass
    
        for fn in fns:
            pth = os.path.join(dirpath, fn)
            print
            print '   ', pth
            T = fits_table(pth)
            T.w4_nanomaggies_ivar *= 0.25
            # mag error = np.abs((-2.5 / np.log(10.)) * dflux / flux)
            #   dflux = 1/sqrt(flux_invvar)
            T.w4_mag_err *= 2.
            
            outpath = os.path.join(outdir, fn)
            T.writeto(outpath)
            print '-> ', outpath
            #print 'Wrote', len(T)

