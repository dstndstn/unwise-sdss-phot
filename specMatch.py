import matplotlib
matplotlib.use('Agg')
from astrometry.util.fits import *
from astrometry.util.starutil_numpy import *
from astrometry.util.plotutils import *
from astrometry.libkd.spherematch import *
import pylab as plt
import numpy as np
import os

from tractor.cache import Cache

ps = PlotSequence('spec')

S = fits_table('specObj-dr10.fits')#, rows=np.arange(1000))
print 'Number of spectra:', len(S)
pdir = 'sdss-dr10d-phot'
T = fits_table('sdss-dr10d-atlas.fits')
print 'Tiles:', len(T)
S.xyz = radectoxyz(S.plug_ra, S.plug_dec)

r = 2.75/3600. * 2048 * np.sqrt(2.)/2. * 1.1
print 'radius:', r

Skd = tree_build_radec(S.plug_ra, S.plug_dec)

print 'Matching...'
#I,J,d = match_radec(S.plug_ra, S.plug_dec, T.ra, T.dec, r, nearest=True)
Tkd = tree_build_radec(T.ra, T.dec)
I,J,d = trees_match(Skd, Tkd, deg2dist(r), nearest=True)
tree_free(Tkd)
print 'Number of matches:', len(I)

K = np.unique(J)
print 'Reading', len(K), 'unique tiles...'
W = []
r = arcsec2dist(2.)
Nspec = 0

wroteSample = False
Imax = 0

for inum,i in enumerate(K):
    print 'Tile', inum+1, 'of', len(K)
    fn = os.path.join(pdir, 'phot-%s.fits' % T.coadd_id[i])
    if not os.path.exists(fn):
        print 'Does not exist:', fn
        continue
    print '  Reading', fn
    w = fits_table(fn, columns=['ra','dec'])
    if len(w) == 0:
        print 'No objects'
        continue
    #print '  building tree'
    wkd = tree_build_radec(w.ra, w.dec)
    #print '  matching'
    I,J,d = trees_match(Skd, wkd, r, nearest=True)
    print '  matched', len(J)
    tree_free(wkd)
    if len(J) == 0:
        continue

    w = fits_table(fn, rows=J, columns=[
        'objid', 'ra', 'dec', 'objc_type', 'run', 'camcol', 'field', 'id',
        'pointsource', 'treated_as_pointsource', 'x', 'y', 'coadd_id',
        'w1_nanomaggies', 'w1_nanomaggies_ivar', 'w1_mag', 'w1_mag_err',
        'w1_prochi2', 'w1_pronpix', 'w1_profracflux', 'w1_proflux',
        'w1_npix', 'w1_pronexp',
        'w2_nanomaggies', 'w2_nanomaggies_ivar', 'w2_mag', 'w2_mag_err',
        'w2_prochi2', 'w2_pronpix', 'w2_profracflux', 'w2_proflux',
        'w2_npix', 'w2_pronexp',
        'w3_nanomaggies', 'w3_nanomaggies_ivar', 'w3_mag', 'w3_mag_err',
        'w3_prochi2', 'w3_pronpix', 'w3_profracflux', 'w3_proflux',
        'w3_npix', 'w3_pronexp',
        'w4_nanomaggies', 'w4_nanomaggies_ivar', 'w4_mag', 'w4_mag_err',
        'w4_prochi2', 'w4_pronpix', 'w4_profracflux', 'w4_proflux',
        'w4_npix', 'w4_pronexp'])
    w.coadd_id = np.array([T.coadd_id[i]] * len(w))

    Scols = S.get_columns()
    for c in w.get_columns():
        x = w.get(c)
        sc = 'wise_' + c
        if not sc in Scols:
            X = np.zeros((len(S),)+x.shape[1:], x.dtype)
            S.set(sc, X)
        else:
            X = S.get(sc)
        X[I] = x
    Nspec += len(I)
    print '  Matched', len(I), 'total', Nspec

    Imax = max(Imax, I.max())
    if not wroteSample and Nspec > 100000:
        wroteSample = True
        outfn = 'data/specmatch-dr10-sample.fits'
        S[:Imax+1].writeto(outfn)
        print 'Wrote sample to', outfn

outfn = 'data/specmatch-dr10.fits'
S.writeto(outfn)
print 'Wrote', outfn
tree_free(Skd)
    
