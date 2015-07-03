from astrometry.util.fits import *
import numpy as np
import os

# DR10
#resolvedir = '/clusterfs/riemann/raid006/dr10/boss/resolve/2010-05-23/'
# v5b / pre-DR13
resolvedir = '/clusterfs/riemann/raid006/bosswork/boss/resolve/2013-07-29'

B = fits_table(os.path.join(resolvedir, 'window_blist.fits'),
               columns=['iprimary'])
print 'Read', len(B), 'balkans'
I = np.unique(B.iprimary)
print len(I), 'fields are primary'

T = fits_table(os.path.join(resolvedir, 'window_flist.fits'),
               rows=I)
print 'Read', len(T), 'fields'
T.cut(T.rerun == '301')
print 'Cut to', len(T), 'rerun 301'
T.writeto('cut-window-flist.fits')


