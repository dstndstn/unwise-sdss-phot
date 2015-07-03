import os

import numpy as np

from astrometry.util.fits import *

A = fits_table('allsky-atlas.fits')
print len(A), 'tiles'
#A = A[:100]

basedir = 'data/unwise/unwise-comp'

bands = [1,2]
for band in bands:
    A.set('n_frames_w%i' % band, np.zeros(len(A), np.int32))
    A.set('n_moon_w%i'   % band, np.zeros(len(A), np.int32))
    A.set('npix_w%i'     % band, np.zeros(len(A), np.int64))
    A.set('sky_w%i'      % band, np.zeros(len(A), np.float32))

for i,a in enumerate(A):

    tile = a.coadd_id
    for band in bands:
        fn = os.path.join(basedir, tile[:3], tile,
                          'unwise-%s-w%i-frames.fits' % (tile, band))
        T = fits_table(fn)
        print len(T), 'from', fn

        #print 'included:', np.unique(T.included)
        #print 'use:', np.unique(T.use)

        T.cut(T.included == 1)
        T.cut(T.use == 1)
        print 'Cut to', len(T)

        #print 'npix:', T.npixoverlap

        npix = np.sum(T.npixoverlap)

        A.get('n_frames_w%i' % band)[i] = len(T)
        A.get('n_moon_w%i'   % band)[i] = np.sum(T.moon_masked)
        A.get('sky_w%i'      % band)[i] = (np.sum(T.npixoverlap * T.sky1) / float(npix))
        A.get('npix_w%i'     % band)[i] = npix

A.writeto('tile-stats.fits')


