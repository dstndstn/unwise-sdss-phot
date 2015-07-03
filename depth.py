if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import numpy as np
import pylab as plt

import sys
from glob import glob

from astrometry.util.fits import *
from astrometry.util.util import *
from astrometry.util.plotutils import *
from astrometry.libkd.spherematch import *

from unwise_coadd import tile_to_radec

def mag_for_sn():
    fns = glob('sdss3-phot/phot-????????.fits')
    fns.sort()
    #fns = list(reversed(fns))
    #fns = fns[:20]
    R = fits_table()
    R.tile = np.array([fn.replace('sdss3-phot/phot-','').replace('.fits','')
                       for fn in fns])
    rd = np.array([tile_to_radec(tile) for tile in R.tile])
    R.ra  = rd[:,0]
    R.dec = rd[:,1]

    nsig = np.arange(4.5, 6.51, 0.5)

    bands = [1,2,3,4]

    Warrs = {}
    for band in bands:
        nsigarrs = []
        nsigmedarrs = []
        nsigfmeanarrs = []
        nsigfmedarrs = []
        for nlo,nhi in zip(nsig, nsig[1:]):
            arr = np.zeros(len(R), np.float32)
            nsigarrs.append(arr)
            R.set('w%i_for_sn_%.1f_%.1f' % (band,nlo,nhi), arr)
    
            arr = np.zeros(len(R), np.float32)
            nsigmedarrs.append(arr)
            R.set('w%i_for_sn_%.1f_%.1f_med' % (band,nlo,nhi), arr)
    
            arr = np.zeros(len(R), np.float32)
            nsigfmeanarrs.append(arr)
            R.set('w%i_for_sn_%.1f_%.1f_fmean' % (band,nlo,nhi), arr)
    
            arr = np.zeros(len(R), np.float32)
            nsigfmedarrs.append(arr)
            R.set('w%i_for_sn_%.1f_%.1f_fmed' % (band,nlo,nhi), arr)
    
        pct = [0., 10., 25., 50., 75., 90., 100.]
        pctarrs = []
        for p in pct:
            arr = np.zeros(len(R), np.float32)
            pctarrs.append(arr)
            R.set('w%iiv_pct%i' % (band,p), arr)

        nn = [np.zeros(len(R), np.float32) for i in range(10)]
        R.set('w%i_nexp_min' % band, nn[0])
        R.set('w%i_nexp_min_in' % band, nn[1])
        R.set('w%i_nexp_1pct'% band, nn[2])
        R.set('w%i_nexp_2pct'% band, nn[3])
        R.set('w%i_nexp_3pct'% band, nn[4])
        R.set('w%i_nexp_4pct'% band, nn[5])
        R.set('w%i_nexp_5pct'% band, nn[6])
        R.set('w%i_nexp_mean'% band, nn[7])
        R.set('w%i_nexp_med' % band, nn[8])
        R.set('w%i_nexp_max' % band, nn[9])

        Warrs[band] = [(nsigarrs, nsigmedarrs, nsigfmeanarrs, nsigfmedarrs), pctarrs, nn]

    keepi = []
    for i,fn in enumerate(fns):
        print
        print (i+1), 'of', len(fns), ':', fn
        cols = ['x','y']
        for band in bands:
            cols.extend([s.replace('1', '%i'%band) for s in
                         ['w1_nanomaggies', 'w1_nanomaggies_ivar', 'w1_mag',
                          'w1_pronexp']])
        T = fits_table(fn, columns=cols)
        if T is None or len(T) == 0:
            continue
        keepi.append(i)

        inbounds = np.flatnonzero((T.x > -0.5) * (T.x < 2047.5) *
                                  (T.y > -0.5) * (T.y < 2047.5))

        for band in bands:
            #I = np.flatnonzero(T.w1_nanomaggies_ivar == 0)
            #print len(I), 'sources with zero invvar'
            wmag  = T.get('w%i_mag' % band)
            wflux = T.get('w%i_nanomaggies' % band)
            wiv   = T.get('w%i_nanomaggies_ivar' % band)
            wn    = T.get('w%i_pronexp' % band)
            wsn   = wflux * np.sqrt(wiv)

            (nsigarrs, pctarrs, nn) = Warrs[band]

            for nlo,nhi, arr,medarr,farr,fmedarr in zip(
                nsig, nsig[1:], *nsigarrs):
                I = np.flatnonzero((wsn >= nlo) * (wsn < nhi))
                arr[i] = np.mean(wmag[I])
                medarr[i] = np.median(wmag[I])
                farr[i] = -2.5 * (np.log10(np.mean(wflux[I])) - 9.)
                fmedarr[i] = -2.5 * (np.log10(np.median(wflux[I])) - 9.)

            for p,arr in zip(pct, pctarrs):
                arr[i] = np.percentile(wiv, p)

            nn[0][i] = np.min(wn)
            nn[1][i] = np.min(wn[inbounds])
            nn[2][i] = np.percentile(wn, 1)
            nn[3][i] = np.percentile(wn, 2)
            nn[4][i] = np.percentile(wn, 3)
            nn[5][i] = np.percentile(wn, 4)
            nn[6][i] = np.percentile(wn, 5)
            nn[7][i] = np.mean(wn)
            nn[8][i] = np.median(wn)
            nn[9][i] = np.max(wn)

    R.cut(np.array(keepi))
    return R

fn = 'spatial-err-stats-3.fits'
if not os.path.exists(fn):
    T = mag_for_sn()
    T.writeto(fn)
else:
    T = fits_table(fn)

ps = PlotSequence('errs')
plt.figure(figsize=(9,4))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

W,H = 1000,500
targetwcs = anwcs_create_hammer_aitoff(90., 0., 1.0, W, H, 0)
print 'Target WCS:', targetwcs

xx,yy = np.meshgrid(np.arange(W), np.arange(H))
ok,rr,dd = targetwcs.pixelxy2radec(xx.ravel(), yy.ravel())
#print 'ok:', np.unique(ok)
#rr = rr.reshape(xx.shape)
#dd = dd.reshape(xx.shape)
#ok = ok.reshape(xx.shape)
K = np.flatnonzero(ok == 0)
rr = rr.flat[K]
dd = dd.flat[K]

J,I,d = match_radec(rr, dd, T.ra, T.dec, 2., nearest=True)
J = K[J]

for col in [c for c in T.get_columns() if 'nexp' in c]:
    X = T.get(col)
    Y = np.zeros(xx.shape, X.dtype)
    Y.flat[J] = X[I]
    Y[np.logical_not(np.isfinite(Y))] = 0.

    ii = np.flatnonzero(Y)
    if len(ii) == 0:
        continue
    mn = Y.flat[ii].min()
    mx = np.percentile(Y.flat[ii], 99)

    plt.clf()
    plt.imshow(Y, interpolation='nearest', origin='lower',
               cmap='jet', vmin=0, vmax=mx)
    plt.colorbar()
    plt.title(col.replace('_',' ').replace('w','W').replace('nexp', 'number of exposures,'))
    plt.xticks([])
    plt.yticks([])
    ps.savefig()



for col in [c for c in T.get_columns() if c.startswith('w1')]:
    X = T.get(col)

    Y = np.zeros(xx.shape, X.dtype)
    Y.flat[J] = X[I]

    Y[np.logical_not(np.isfinite(Y))] = 0.

    ii = np.flatnonzero(Y)
    if len(ii) == 0:
        continue
    
    mn = Y[Y != 0].min()

    mx = np.percentile(Y[Y!= 0], 99)

    plt.clf()
    plt.imshow(Y, interpolation='nearest', origin='lower',
               cmap='jet', vmin=mn, vmax=mx)
    plt.colorbar()
    plt.title(col)
    plt.xticks([])
    plt.yticks([])
    ps.savefig()

    if not 'iv' in col:

        yi = Y[Y != 0]
        p5 = np.percentile(yi, 5)
        p10 = np.percentile(yi, 10)
        print ' 5th percentile of', col, ':', p5
        print '10th percentile of', col, ':', p10


        continue

    lo = np.percentile(Y[Y != 0], 10)
    print '10th percentile of', col, ':', lo
    
    # "What is the S/N of a W1=19.1 source?"

    targetflux = 10.**((19.1 - 22.5)/-2.5)
    err = 1. / np.sqrt(Y)
    sn = targetflux / err
    sn[Y == 0] = 0.

    mn = sn[sn > 0].min()
    mx = np.percentile(sn[sn > 0], 99)
    mn -= (mx-mn)*0.1
    
    plt.clf()
    plt.imshow(sn, interpolation='nearest', origin='lower',
               cmap='jet', vmin=mn, vmax=mx)
    plt.colorbar()
    plt.title('S/N of a W1=19.1 source (%s)' % col)
    plt.xticks([])
    plt.yticks([])
    ps.savefig()

    sn = sn[sn > 0]
    print
    print col
    for blo,bhi in [(0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0),
                    (2.0, 2.5), (2.5, 3.0), (3.0, 1e9)]:
        L = np.flatnonzero((sn > blo) * (sn <= bhi))
        print 'S/N in [%.1f, %.1f]: %.1f %%' % (blo,bhi, 100.*len(L)/float(len(sn)))

    
sys.exit(0)



T = fits_table('sdss2-phot/phot-1800p530.fits')
#T.about()
print len(T), 'WISE'

print 'W objid', T.objid.dtype

S = fits_table('sdss2-phot-temp/photoobjs-1800p530.fits')
print len(S), 'SDSS'
print 'S objid', S.objid.dtype

objidmap = dict(zip(S.objid, np.arange(len(S))))
I = np.array([objidmap[o] for o in T.objid])
S.cut(I)
print len(S), 'SDSS on objid'

flux = T.w1_nanomaggies
dflux = 1./np.sqrt(T.w1_nanomaggies_ivar)

flux2 = T.w2_nanomaggies
dflux2 = 1./np.sqrt(T.w2_nanomaggies_ivar)

ps = PlotSequence('depth')

# plt.clf()
# plt.hist(flux, 100, log=True, range=(-1e5, 1e6))
# ps.savefig()
# 
# plt.clf()
# plt.hist(flux, 100, range=(-1000, 1000))
# ps.savefig()
# 
# plt.clf()
# plt.hist(np.log10(np.maximum(1., flux)), 100)
# plt.xlabel('log Flux')
# ps.savefig()
# 
# plt.clf()
# plt.hist(dflux, 100, range=(0,25))
# plt.xlabel('delta-flux')
# ps.savefig()
# 
# plt.clf()
# loghist(flux, dflux, 100, range=((-1000,1000),(0, 25)))
# plt.xlabel('flux')
# plt.ylabel('delta-flux')
# ps.savefig()
# 
# plt.clf()
# loghist(flux, flux/dflux, 100, range=((-1000,1000),(-10,100)))
# plt.xlabel('flux')
# plt.ylabel('n sigma')
# ps.savefig()
# 
# plt.clf()
# loghist(flux, flux/dflux, 100, range=((-100,250),(-10,20)))
# plt.xlabel('flux')
# plt.ylabel('n sigma')
# ps.savefig()


for nsigma,magname in [(flux / dflux, 'W1'),
                       (flux2 / dflux2, 'W2'),
                       (np.hypot(flux/dflux, flux2/dflux2), 'hypot(W1,W2)'),
                       ((flux/dflux) + 0.5*(flux2/dflux2), 'W1+W2/2'),
                       ((flux/dflux) + (flux2/dflux2), 'W1+W2')]:
    plt.clf()
    lp,lt = [],[]
    for i,(nlo,nhi) in enumerate([(2.0, 2.5),
                                  (2.5, 3.0),
                                  (3.0, 3.5),
                                  (3.5, 4.0),
                                  (4.0, 4.5),
                                  (4.5, 5.0),
                                  (5.0, 5.5),]):
        I = np.flatnonzero((nsigma >= nlo) * (nsigma < nhi))
        if len(I) == 0:
            continue
        #cmap = matplotlib.cm.RdYlBu
        cmap = matplotlib.cm.jet

        plt.subplot(2,1,1)
        n,b,p = plt.hist(T.w1_mag[I], 80, range=(17.5, 20), histtype='step',
                         color=cmap(i/float(6)))
        plt.xlabel('W1')
        plt.xlim(17.5, 20)
        plt.subplot(2,1,2)
        n,b,p = plt.hist(T.w2_mag[I], 80, range=(15, 20), histtype='step',
                         color=cmap(i/float(6)))
        plt.xlim(15, 20)
        plt.xlabel('W2')

        print 'Median W1 mag for N sigma (in %s) in [%.1f, %.1f]: %.2f' % (magname, nlo, nhi, np.median(T.w1_mag[I]))
        print 'Median W2 mag for N sigma (in %s) in [%.1f, %.1f]: %.2f' % (magname, nlo, nhi, np.median(T.w2_mag[I]))

        lp.append(p[0])
        lt.append('Nsigma in [%.1f, %.1f]' % (nlo,nhi))

    plt.figlegend(reversed(lp), reversed(lt), loc='upper right')
    plt.suptitle('WISE mags by Nsigma: %s' % magname)
    ps.savefig()


wflux = (T.w1_nanomaggies + 0.5 * T.w2_nanomaggies) / 1.5

w1var = 1./T.w1_nanomaggies_ivar
w2var = 1./T.w2_nanomaggies_ivar
wflux_var = (w1var + 0.25 * w2var) / 1.5**2
wsn = wflux / np.sqrt(wflux_var)

rmag = -2.5 * (np.log10(S.modelflux[:,2]) - 9)
plt.clf()
loghist(rmag, wsn, 100, range=((17,24),(-3,20)))
plt.axhline(0, color='b')
#plt.axvline(22, color='b')
plt.xlabel('r mag')
plt.ylabel('WISE S/N')
plt.title('forced photometry: W1 + 0.5*W2 s/n vs r-band mag')
ps.savefig()

mags = np.arange(19, 24, 1.)

plt.clf()
plt.subplots_adjust(hspace=0.3)
lp,lt = [],[]
for i,(maglo,maghi) in enumerate(zip(mags, mags[1:])):
    I = np.flatnonzero((rmag > maglo) * (rmag <= maghi))
    if len(I) == 0:
        continue
    lo,hi = (-3,10)
    n,b,p = plt.hist(np.clip(wsn[I], lo, hi), 65, range=(lo,hi),
                     histtype='step',
                     color=cmap(i/float(3)))
    lp.append(p[0])
    lt.append('r mag in [%.1f, %.1f]' % (maglo, maghi))
    plt.xlabel('S/N in W1+0.5*W2')
    plt.xlim(lo,hi)
plt.ylim(0, 500)
#plt.figlegend(reversed(lp), reversed(lt), loc='upper right')
plt.figlegend(lp, lt, loc='upper right')
#plt.suptitle('W1+W2/2 mag in [%.1f, %.1f]' % (maglo, maghi))
plt.suptitle('WISE S/N wrt r mag')
ps.savefig()


wmag = -2.5*(np.log10(wflux) - 9)

mags = np.arange(17, 20.1, 0.5)

plt.clf()
plt.subplots_adjust(hspace=0.3)
lp,lt = [],[]
for i,(maglo,maghi) in enumerate(zip(mags, mags[1:])):
    I = np.flatnonzero((wmag > maglo) * (wmag <= maghi))
    if len(I) == 0:
        continue
    plt.subplot(2,1,1)
    lo,hi = (0,10)
    n,b,p = plt.hist(np.clip((flux/dflux)[I], lo, hi), 50, range=(lo,hi),
                     histtype='step',
                     color=cmap(i/float(6)))
    lp.append(p[0])
    lt.append('Wmag in [%.1f, %.1f]' % (maglo, maghi))
    plt.xlabel('S/N in W1')
    plt.xlim(lo,hi)
    plt.subplot(2,1,2)
    plt.hist(np.clip((flux2/dflux2)[I], lo, hi), 50, range=(lo, hi), 
             histtype='step',
             color=cmap(i/float(6)))
    plt.xlabel('S/N in W2')
    plt.xlim(lo,hi)
plt.subplot(2,1,1)
plt.ylim(0, 320)
plt.subplot(2,1,2)
plt.ylim(0, 600)
#plt.figlegend(reversed(lp), reversed(lt), loc='upper right')
plt.figlegend(lp, lt, loc='lower right')
#plt.suptitle('W1+W2/2 mag in [%.1f, %.1f]' % (maglo, maghi))
plt.suptitle('WISE S/N wrt W1+0.5*W2 mag')
ps.savefig()


