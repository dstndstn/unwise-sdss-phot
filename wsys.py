
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

from scipy.ndimage.morphology import binary_dilation, binary_fill_holes

from astrometry.util.fits import *
from astrometry.util.plotutils import *
from astrometry.libkd.spherematch import *

ps = PlotSequence('wsys')

T = fits_table('wise_systematics_map.fits')
W = fits_table('tile-stats.fits')

# plt.clf()
# plt.plot(T.ra, T.dec, 'r.')
# plt.plot(W.ra, W.dec, 'bx')
# ps.savefig()

rr,dd = np.meshgrid(np.arange(360), np.arange(-90,91))
print 'rr,dd', rr.shape, dd.shape

TI,TJ,d = match_radec(rr.ravel(), dd.ravel(), T.ra, T.dec, 2., nearest=True)

bdy = np.zeros(rr.shape, bool)
bdy.flat[TI] = True
bdy = binary_fill_holes(bdy)
dilated = binary_dilation(bdy, np.ones((3,3)))
dilated &= np.logical_not(bdy)
bdy = dilated
del dilated

J,I,d = match_radec(rr.ravel(), dd.ravel(), W.ra, W.dec, 2., nearest=True)
print len(I), len(J)

plt.figure(figsize=(12,5))
plt.subplots_adjust(left=0.05, right=0.98, top=0.9)

for col in [#'n_frames_w1', 'n_frames_w2',
            'npix_w1', 'npix_w2',
            'n_moon_w1', 'n_moon_w2',
            'sky_w1', 'sky_w2']:
    x = W.get(col)

    if 'moon' in col:
        y = W.get(col.replace('moon', 'frames'))
        x = x.astype(np.float32) / y.astype(np.float32)
        col = col.replace('n_moon', 'frac_moon')

    xmap = np.zeros(rr.shape, x.dtype)
    xmap.flat[J] = x[I]

    if col.startswith('npix'):
        xmap /= (2048.**2)
        mx = 100
    elif col == 'sky_w1':
        mx = 100
    elif col == 'sky_w2':
        mx = 200
    else:
        mx = xmap.max()

    plt.clf()
    plt.imshow(xmap, interpolation='nearest', origin='lower',
               extent=[0,360,-90,90], vmin=0, vmax=mx)
    ax = plt.axis()
    plt.plot(rr[bdy], dd[bdy], 'r.')
    plt.axis(ax)
    plt.colorbar()
    plt.title(col)
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    ps.savefig()
    
