import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='computer modern roman')
matplotlib.rc('font', **{'sans-serif': 'computer modern sans serif'})
#matplotlib.rc('font', **{'sans-serif': 'computer modern roman'})
import numpy as np
import pylab as plt
import os
import sys

from collections import Counter

from astrometry.util.file import *
from astrometry.util.fits import *
from astrometry.util.multiproc import *
from astrometry.util.plotutils import *
from astrometry.util.miscutils import *
from astrometry.util.util import *
from astrometry.util.resample import *
from astrometry.libkd.spherematch import *
from astrometry.util.starutil_numpy import *
from astrometry.sdss.dr10 import *

from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.morphology import binary_dilation

import fitsio

from forcedphot import read_wise_sources

from tractor.sdss import *

from forcedphot import treat_as_pointsource

def sind(t):
    return np.sin(np.deg2rad(t))
def cosd(t):
    return np.cos(np.deg2rad(t))

from matplotlib.colors import LinearSegmentedColormap

lightgray = LinearSegmentedColormap('lightgray',
                                    {'red':   ((0., np.nan, 1.), (1., 0.5, np.nan)),
                                     'green': ((0., np.nan, 1.), (1., 0.5, np.nan)),
                                     'blue':  ((0., np.nan, 1.), (1., 0.5, np.nan))})




#http://data.sdss3.org/sas/bosswork/boss/qso/DR12Q/Final/DR12Q.fits
def sequels_qso_forcedphot():

    stages = SequelsStages()
    prereqs = { 'wise': 'sdss',
                'sdss': None,
                }
    from astrometry.util.stages import *

    plt.figure(1, figsize=(10,10))
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                        hspace=0.1, wspace=0.1)

    plt.figure(2, figsize=(3, 3))
    #plt.subplots_adjust(left=0.005, right=0.995, bottom=0.005, top=0.995,
    #                    hspace=0.1, wspace=0.1)
    plt.subplots_adjust(left=0., right=1., bottom=0., top=1.)

    plt.figure(1)
    
    runstage('wise', 'sequels-%s.pickle', stages, prereqs=prereqs, force=['wise'])
    
    # T = fits_table('/clusterfs/riemann/raid006/bosswork/boss/qso/SEQUELS/SEQUELS_v2_1.fits',
    #                column_map={'class':'clazz'})
    # print len(T), 'quasars'
    # T.about()
    # 
    # I = np.flatnonzero(T.zwarning == 0)
    # print len(I), 'with Zwarning = 0'
    # 
    # uu = np.unique(T.eboss_target0)
    # # print 'eboss_target0:', uu
    # bits = 0
    # for u in uu:
    #     bits |= u
    # print 'all bits: 0x%x' % bits
    # 
    # keepbits = (1<<64)-1
    # # 19 DR9_CALIB_TARGET
    # keepbits -= (1<<19)
    # 
    # targetbits = T.eboss_target0 & keepbits
    # 
    # print sum((targetbits & 1024) > 0), 'with QSO_CORE bit set'
    # I = np.flatnonzero((targetbits & 1024) > 0)
    # print 'Bit combos:'
    # for u in np.unique(targetbits[I]):
    #     for i in range(62):
    #         if ((1 << i) & int(u)):
    #             print i,
    #     print
    # 
    # print sum((targetbits  == 1024)), 'with ONLY QSO_CORE bit set'
    # 
    # #print 'Classes:', np.unique(T.clazz)
    # #qso = 'QSO   '
    # #I = np.flatnonzero((T.zwarning == 0) * (T.clazz == qso))
    # #print len(I), 'QSO with Zwarning = 0'

class SequelsStages(object):

    def __init__(self):
        pass

    def __call__(self, stage, **kwargs):
        print 'stage', stage
        func = getattr(self, 'stage_'+stage)
        print 'func', func
        return func(**kwargs)

    def stage_sdss(self, **kwargs):
        #T = fits_table('/home/schlegel/wise1ext/sdss/zans-plates7027-7032-zscan.fits', column_map={'class':'clazz'})
        T = fits_table('zans-plates7027-7032-zscan.fits', column_map={'class':'clazz'})
        print 'Read', len(T), 'zscan'
    
        I = np.flatnonzero(T.zwarning == 0)
        print len(I), 'with Zwarning = 0'
    
        print 'Classes:', np.unique(T.clazz)
    
        qso = 'QSO   '
        I = np.flatnonzero((T.zwarning == 0) * (T.clazz == qso))
        print len(I), 'QSO with Zwarning = 0'
        
        print 'Typescans:', np.unique(T.typescan)
        qsoscan = 'QSO    '
        I = np.flatnonzero((T.zwarning == 0) * (T.clazz == qso) * (T.typescan == qsoscan))
        print len(I), 'QSO, scan QSO, with Zwarning = 0'
        
        # for zcut in [2, 2.3, 2.5]:
        #     I = np.flatnonzero((T.zwarning == 0) * (T.clazz == qso) * (T.typescan == qsoscan) * (T.z > zcut))
        #     print len(I), 'QSO, scan QSO, with Zwarning = 0, z >', zcut
        #     
        #     I = np.flatnonzero((T.zwarning == 0) * (T.clazz == qso) * (T.typescan == qsoscan) * (T.zscan > zcut))
        #     print len(I), 'QSO, scan QSO, with Zwarning = 0, z and zscan >', zcut
            
        #zcut = 2.3
        #zcut = 2.0
        zcut = 1.0
        I = np.flatnonzero((T.zwarning == 0) * (T.clazz == qso) * (T.typescan == qsoscan) * (T.zscan > zcut))
        print len(I), 'QSO, scan QSO, with Zwarning = 0, z and zscan >', zcut
        T.cut(I)

        A = fits_table('ancil-QSO-eBOSS-W3-ADM-dr8.fits')
        print 'Read', len(A), 'targets'
        I,J,d = match_radec(T.plug_ra, T.plug_dec, A.ra, A.dec, 1./3600.)
        print len(I), 'matched'
        print len(np.unique(I)), 'unique fibers'
        print len(np.unique(J)), 'unique targets'
        
        A.cut(J)
        T.cut(I)

        I = np.flatnonzero(A.w3bitmask == 4)
        print 'Selected by WISE only:', len(I)
        A.cut(I)
        T.cut(I)
        T.ra = T.plug_ra
        T.dec = T.plug_dec

        I = np.argsort(-T.zscan)
        T.cut(I)
        A.cut(I)

        T.about()
        A.about()
    
        T.run = A.run
        T.camcol = A.camcol
        T.field = A.field
    
        sdss = DR10()
        sdss.useLocalTree()
        #sdss.setFitsioReadBZ2()
        sdss.saveUnzippedFiles('data/tmp')

        #rows,cols = 6,6
        rows,cols = 10,10
        pixscale = 0.396 / 3600.
        #sz = 50
        #sz = 45
        sz = 40
        stamps,Istamp,stampwcs,rawstamps = get_sdss_cutouts(
            sdss, T, pixscale, sz, rows*cols, get_rawvals=True)
    
        ps = PlotSequence('sequels')
    
        plt.clf()
        for i,stamp in enumerate(stamps):
            plt.subplot(rows, cols, i+1)
            dimshow(stamp)
            plt.xticks([]); plt.yticks([])
        ps.savefig()

        return dict(ps=ps, T=T, stamps=stamps, Istamp=Istamp, stampwcs=stampwcs,
                    rows=rows, cols=cols, sdss=sdss, rawstamps=rawstamps)

    def stage_wise(self, ps=None, T=None, Istamp=None, stampwcs=None,
                   rows=None, cols=None, stamps=None, rawstamps=None, **kwargs):
        atlas = fits_table('w3-atlas.fits')

        ps.suffixes = ['png','pdf']

        from unwise_coadd import get_coadd_tile_wcs
        from tractor import *

        rgbimgs = []
        rgbmods = []
        rgbumods = []

        sdss = DR10()
        sdss.useLocalTree()
        #sdss.setFitsioReadBZ2()
        sdss.saveUnzippedFiles('data/tmp')

        print 'stamps:', len(stamps)
        print 'rawstamps:', len(rawstamps)

        for si,(swcs,ra,dec,stamp,ti,rawstamp) in enumerate(
            zip(stampwcs, T.ra[Istamp], T.dec[Istamp], stamps, Istamp, rawstamps)):

            if si != 66:
                continue

            photdir = 'w3-phot'
            #photdir = 'sdss-dr10d-phot'

            sh,sw = int(swcs.get_height()), int(swcs.get_width())
            re1 = np.zeros((sh,sw), np.float32)
            re2 = np.zeros((sh,sw), np.float32)
            remod1 = np.zeros((sh,sw), np.float32)
            remod2 = np.zeros((sh,sw), np.float32)
            reumod1 = np.zeros((sh,sw), np.float32)
            reumod2 = np.zeros((sh,sw), np.float32)
            ren = np.zeros((sh,sw), np.uint8)

            k = np.argmin(degrees_between(atlas.ra, atlas.dec, ra, dec))
            print 'nearest coadd tile:', atlas.coadd_id[k]
            rc,dc = atlas.ra[k], atlas.dec[k]
            co = atlas.coadd_id[k]

            # run,camcol,field = T.run[ti], T.camcol[ti], T.field[ti]
            # po = sdss.readPhotoObj(run, camcol, field).getTable()
            # print 'photoObj:', po
            # print len(po), 'photoObjs'
            fn = os.path.join(photdir + '-temp', 'photoobjs-%s.fits' % co)
            po = fits_table(fn)
            print len(po), 'photoObjs'
            ok, x,y = swcs.radec2pixelxy(po.ra, po.dec)
            I = np.flatnonzero((x >= -1) * (y >= -1) * (x < sw) * (y < sw))
            po.x = x-1
            po.y = y-1
            po.cut(I)
            print len(po), 'within stamp'
            Ssrcs = po

            Wsrcs = []

            print 'Stamp WCS bounds', swcs.radec_bounds()

            cowcs = get_coadd_tile_wcs(rc, dc)
            if not cowcs.is_inside(ra, dec):
                print 'RA,Dec not inside tile', atlas.coadd_id[k]
                continue
    
            print 'source', si, 'RA,Dec', ra, dec, 'tile index', k, co
            # ok, x,y = cowcs.radec2pixelxy(ra, dec)
            # print 'x,y', x,y
            # x = int(np.round(x)-1)
            # y = int(np.round(y)-1)
            # S = 6
            # ix0 = max(0, x-S)
            # iy0 = max(0, y-S)
            # slc = (slice(max(0, y-S), min(H, y+S+1)),
            #        slice(max(0, x-S), min(W, x+S+1)))
            # imextent = [slc[1].start - 0.5, slc[1].stop - 0.5,
            #             slc[0].start - 0.5, slc[0].stop - 0.5]
            
            imgs = []
            for band in [1,2]:
                fn = 'unwise-coadds/%s/%s/unwise-%s-w%i-img-m.fits' % (co[:3], co, co, band)
                if not os.path.exists(fn):
                    print 'Does not exist:', fn
                    continue
                img = fitsio.read(fn)
                imgs.append(img)
            if len(imgs) != 2:
                continue

            fn = os.path.join(photdir, 'phot-%s.fits' % co)
            print 'Reading', fn
            phot = fits_table(fn)
            #phot.about()
            print len(phot), 'photometry results for', co
            ok,x,y = swcs.radec2pixelxy(phot.ra, phot.dec)
            I = np.flatnonzero((x >= -1) * (y >= -1) * (x < sw) * (y < sw))
            phot.x = x-1
            phot.y = y-1
            phot.cut(I)
            print len(phot), 'within stamp'

            phot.ptsrc = np.logical_or(phot.pointsource,
                                       phot.treated_as_pointsource)
            phot.wiseonly = np.zeros(len(phot), bool)

            fn = os.path.join(photdir, 'phot-wise-%s.fits' % co)
            print 'Reading', fn, 'WISE-only'
            if os.path.exists(fn):
                wphot = fits_table(fn)
                ok,x,y = swcs.radec2pixelxy(wphot.ra, wphot.dec)
                I = np.flatnonzero((x >= -1) * (y >= -1) * (x < sw) * (y < sw))
                wphot.cut(I)
                print len(wphot), 'within stamp'
                wphot.ptsrc = np.ones(len(wphot), bool)
                wphot.wiseonly = np.ones(len(wphot), bool)
                wphot.w1_nanomaggies = wphot.fit_flux_w1
                wphot.w2_nanomaggies = wphot.fit_flux_w2
    
                phot = merge_tables([phot,wphot], columns='fillzero')
            else:
                print 'DOES NOT EXIST:', fn

            phot.prob_psf = np.zeros((len(phot),5))
            phot.prob_psf[:,2] = phot.ptsrc
            phot.psfflux = np.zeros((len(phot),5))
            phot.cmodelflux = np.zeros((len(phot),5))
            phot.devflux = np.zeros((len(phot),5))
            phot.expflux = np.zeros((len(phot),5))

            mods = []
            umods = []
            for band in [1,2]:
                psffn = 'psf-allwise-con3.fits'
                P = fits_table(psffn, hdu=band)
                psf = GaussianMixturePSF(P.amp, P.mean, P.var)
                twcs = ConstantFitsWcs(cowcs)
                                            
                tim = Image(data=img, invvar=np.ones_like(img),
                            psf=psf, wcs=twcs,
                            photocal=LinearPhotoCal(1., band='w'))

                srcs = get_tractor_sources_dr9(None,None,None,'r',
                                               objs=phot, bands=[],
                                               extrabands=['w'],
                                               nanomaggies=True,
                                               fixedComposites=True)
                print 'srcs', srcs
                cat = Catalog(*srcs)
                cat.freezeAllRecursive()
                cat.thawPathsTo('w')
                print len(phot), 'fluxes'
                print 'Cat params:', cat.numberOfParams()
                cat.printThawedParams()
                cat.setParams(phot.get('w%i_nanomaggies' % band))
                print 'srcs', srcs

                tractor = Tractor([tim], srcs)
                mod = tractor.getModelImage(0)
                mods.append(mod)

                # find closest source
                ii = np.argmin(np.hypot(phot.ra - ra, phot.dec - dec))
                # render model of just that source
                tractor = Tractor([tim], [srcs[ii]])
                mod = tractor.getModelImage(0)
                umods.append(mod)
    
            try:
                Yo,Xo,Yi,Xi,nil = resample_with_wcs(swcs, cowcs, [], 3)
                re1[Yo,Xo] += imgs[0][Yi,Xi]
                re2[Yo,Xo] += imgs[1][Yi,Xi]
                remod1[Yo,Xo] += mods[0][Yi,Xi]
                remod2[Yo,Xo] += mods[1][Yi,Xi]
                reumod1[Yo,Xo] += umods[0][Yi,Xi]
                reumod2[Yo,Xo] += umods[1][Yi,Xi]
                ren[Yo,Xo] += 1
            except:
                import traceback
                traceback.print_exc()
                continue


            fn = os.path.join(photdir + '-temp', 'wise-sources-%s.fits' % co)
            W = fits_table(fn)
            print len(W), 'WISE sources in', co
            ok,x,y = swcs.radec2pixelxy(W.ra, W.dec)
            I = np.flatnonzero((x >= -1) * (y >= -1) * (x < sw) * (y < sw))
            W.x = x-1
            W.y = y-1
            W.cut(I)
            print len(W), 'within stamp'
            Wsrcs.append(W)
    
            re1 /= np.maximum(ren, 1)
            re2 /= np.maximum(ren, 1)
            remod1 /= np.maximum(ren, 1)
            remod2 /= np.maximum(ren, 1)
            reumod1 /= np.maximum(ren, 1)
            reumod2 /= np.maximum(ren, 1)
            
            rgb = np.zeros((sh,sw,3), np.float32)
            rgb[:,:,0] = (re2  + 10.) / 50.
            rgb[:,:,2] = (re1  + 10.) / 50.
            rgb[:,:,1] = (rgb[:,:,0] + rgb[:,:,2]) / 2.
            Wsrcs = merge_tables(Wsrcs)
            rgbimgs.append((rgb.copy(),Wsrcs))

            rgb[:,:,0] = (remod2  + 10.) / 50.
            rgb[:,:,2] = (remod1  + 10.) / 50.
            rgb[:,:,1] = (rgb[:,:,0] + rgb[:,:,2]) / 2.
            rgbmods.append(rgb.copy())

            rgb[:,:,0] = (reumod2  + 10.) / 50.
            rgb[:,:,2] = (reumod1  + 10.) / 50.
            rgb[:,:,1] = (rgb[:,:,0] + rgb[:,:,2]) / 2.
            rgbumods.append(rgb.copy())

            mima = dict(interpolation='nearest', origin='lower')

            plt.clf()
            plt.subplot(2,2,1)
            plt.imshow(stamp, **mima)
            ax = plt.axis()
            ssty = dict(mew=2, mfc='none', mec='r', ms=15)
            plt.plot(Ssrcs.x, Ssrcs.y, 'o', **ssty)
            plt.axis(ax)
            plt.subplot(2,2,2)
            plt.imshow(np.clip(rgbmods[-1],0,1), **mima)
            plt.subplot(2,2,3)
            rgb,Wsrcs = rgbimgs[-1]
            plt.imshow(np.clip(rgb,0,1), **mima)
            ax = plt.axis()
            plt.plot(Wsrcs.x, Wsrcs.y, 'm+', mew=2, ms=15)
            plt.plot(Ssrcs.x, Ssrcs.y, 'o', **ssty)
            plt.axis(ax)
            plt.subplot(2,2,4)
            plt.imshow(np.clip(rgbumods[-1],0,1), **mima)
            plt.suptitle('%i: z=%.2f, %.2f, %s, %.3f,%.3f' % (si, T.z[ti], T.zscan[ti], co, T.ra[ti], T.dec[ti]))
            ps.savefig()

            plt.figure(2)

            i,r,g = rawstamp
            i *= 1.0
            r *= 1.3
            g *= 2.5
            rgb = np.dstack([i,r,g])
            lo,hi = -0.1, 0.5
            rgb = (rgb - lo) / (hi - lo)

            plt.clf()
            dimshow(np.clip(rgb, 0., 1.), ticks=False)
            plt.axis('off')
            ps.savefig()

            # plt.clf()
            # dimshow(stamp, ticks=False)
            # plt.axis('off')
            # ps.savefig()

            plt.clf()
            rgb,Wsrcs = rgbimgs[-1]
            dimshow(np.clip(rgb,0,1), ticks=False)
            ax = plt.axis()
            plt.plot(Wsrcs.x, Wsrcs.y, 'rx', mew=2, ms=15)
            #plt.plot(Ssrcs.x, Ssrcs.y, 'o', **ssty)
            plt.axis(ax)
            plt.axis('off')
            ps.savefig()

            plt.clf()
            dimshow(np.clip(rgbmods[-1],0,1), ticks=False)
            ax = plt.axis()
            #plt.plot(Wsrcs.x, Wsrcs.y, 'rx', mew=2, ms=15)
            plt.plot(Ssrcs.x, Ssrcs.y, '+', **ssty)
            plt.axis(ax)
            plt.axis('off')
            ps.savefig()

            plt.clf()
            dimshow(np.clip(rgbumods[-1],0,1), ticks=False)
            plt.axis('off')
            ps.savefig()

            # plt.clf()
            # dimshow(np.clip(np.dstack(rawstamp), 0., 1.), ticks=False)
            # plt.axis('off')
            # ps.savefig()

            plt.figure(1)

            
        plt.clf()
        for i,(rgb,Wsrcs) in enumerate(rgbimgs):
            plt.subplot(rows,cols, 1+i)
            plt.imshow(np.clip(rgb, 0, 1), **mima)
            ax = plt.axis()
            plt.plot(Wsrcs.x, Wsrcs.y, 'r+', mew=2)
            plt.axis(ax)
        ps.savefig()

        
        plt.clf()
        for i,rgb in enumerate(rgbmods):
            plt.subplot(rows,cols, 1+i)
            plt.imshow(np.clip(rgb, 0, 1), **mima)
        ps.savefig()
        
        plt.clf()
        for i,rgb in enumerate(rgbumods):
            plt.subplot(rows,cols, 1+i)
            plt.imshow(np.clip(rgb, 0, 1), **mima)
        ps.savefig()




    

class AtlasTiles(object):
    def __init__(self):
        self.A = fits_table('allsky-atlas.fits')
        self.decrings = np.unique(self.A.dec)

    def cut_to_unique(self, tilera, tiledec, ra, dec):
        '''
        Given a tile RA,Dec center and arrays of ra,dec, returns an index array of the
        elements uniquely within the unique area of the tile.
        '''
        decabove = self.decrings[self.decrings > (tiledec + 0.1)]
        if len(decabove):
            #print 'Dec ring above:', decabove[0]
            dechi = (tiledec + decabove[0]) / 2.
        else:
            dechi = 90.
        decbelow = self.decrings[self.decrings < (tiledec - 0.1)]
        if len(decbelow):
            #print 'Dec ring below:', decbelow[-1]
            declo = (tiledec + decbelow[-1]) / 2.
        else:
            declo = -90.
        #print 'Dec bounds', declo, dechi

        raring = self.A.ra[self.A.dec == tiledec]
        
        dra = raring - tilera
        dra += ( 360. * (dra < -180.))
        dra += (-360. * (dra >  180.))
        dralo = np.max(dra[dra < 0]) / 2.
        drahi = np.min(dra[dra > 0]) / 2.

        I = np.flatnonzero((dec >= declo) * (dec < dechi))
        if len(I) == 0:
            return I
        dra = ra[I] - tilera
        dra += ( 360. * (dra < -180.))
        dra += (-360. * (dra >  180.))
        I = I[(dra >= dralo) * (dra < drahi)]
        return I


def bright(T, tdir, pdir, unwdir):
    ps = PlotSequence('bright')
    ps.suffixes = ['png', 'pdf']

    plt.figure(figsize=(5.5,5))
    #plt.subplots_adjust(left=0.15, right=0.99, bottom=0.1, top=0.95, hspace=0.05, wspace=0.05)
    #plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, hspace=0.05, wspace=0.05)
    #plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
    plt.subplots_adjust(left=0.005, right=0.95, bottom=0.005, top=0.95)

    # Bright-star bins
    brightbins = np.array([-10., -2.0, -1.5, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

    B = fits_table('data/wphot/allwise-bright-2.fits')
    #B.cut(B.w1mpro < 3.)
    B.cut(B.w1mpro < brightbins[-1])

    bin0 = np.floor(np.min(B.w1mpro))
    brightbins[0] = bin0

    II = []
    for blo,bhi in zip(brightbins, brightbins[1:]):
        I = np.flatnonzero((B.w1mpro >= blo) * (B.w1mpro < bhi))
        print len(I), 'in bin', blo, bhi
        # Cut:
        I = I[:20]
        II.append(I)
    B.cut(np.hstack(II))

    B.cut(np.argsort(B.w1mpro))
    print len(B), 'bright: W1', B.w1mpro

    AT = AtlasTiles()

    radius = 0.2

    statsOnly = True
    #statsOnly = False

    if statsOnly:
        ps = PlotSequence('brightstats')
        ps1 = PlotSequence('bright1')
        plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, hspace=0.05, wspace=0.05)

        if os.path.exists('bright-hists.pickle'):
            X = unpickle_from_file('bright-hists.pickle')
            hists = X['hists']

            meths = ['Forced Photometry', 'AllWISE']
            bands = ['W1', 'W2']
            brbins = np.arange(len(brightbins)-1)
            for meth in meths:
                for band in bands:
                    for brbin in brbins:
                        key = (meth, band, brbin)
                        if not key in hists:
                            continue
                        H,nstars,radbins,magbins,bright = hists[key]
                        # print 'radius bins', radbins
                        # print 'mag bins', magbins
                        # print 'N rad bins:', len(radbins)
                        # print 'N mag bins:', len(magbins)
                        # print 'Hist shape', H.shape

                        H /= nstars
    
                        weight = 1. / (np.pi * (radbins[1:]**2 - radbins[:-1]**2) / 3600.**2)
                        mlo,mhi = magbins.min(), magbins.max()
                        rlo,rhi = radbins.min(), radbins.max()

                        dmag = magbins[1] - magbins[0]

                        (blo,bhi) = bright
                        print 'N stars', nstars, 'in brightness bin', blo,bhi
    
                        # plt.clf()
                        # dimshow(H, cmap='hot', extent=[rlo,rhi, mlo,mhi])
                        # plt.gca().set_aspect('auto')
                        # plt.xlabel('Radius from bright star (arcsec)')
                        # plt.ylabel('Mag')
                        # plt.ylim(mhi, mlo)
                        # plt.title('%s %s: mag range %.1f to %.1f' % (meth, band, blo,bhi))
                        # ps.savefig()
    
                        plt.clf()
                        dens = H * weight[np.newaxis,:] / dmag
                        #mx = np.percentile(dens, 99)

                        mx = 12000 # stars / sq deg / mag

                        dimshow(dens, cmap='hot', vmax=mx, vmin=0, extent=[rlo,rhi, mlo,mhi])
                        plt.gca().set_aspect('auto')
                        plt.xlabel('Radius from bright star (arcsec)')
                        plt.ylabel('Mag')
                        plt.ylim(mhi, mlo)
                        plt.title('%s %s: mag range %.1f to %.1f: density' %
                                  (meth, band, blo,bhi))
                        plt.colorbar()
                        ps.savefig()
            

            sys.exit(0)
            
        hists = {}
        radius = 0.4

    else:
        # Cut
        B = B[:50]

    for i, (tile, ra, dec, w1mpro) in enumerate(zip(B.coadd_id, B.ra, B.dec, B.w1mpro)):
        print
        print 'Bright star', i, 'of', len(B), 'RA,Dec', ra, dec, 'w1mpro', w1mpro
        print

        cosdec = np.cos(np.deg2rad(dec))
        rlo,rhi = ra - radius/cosdec, ra + radius/cosdec
        dlo,dhi = dec - radius, dec + radius

        pixscale = 2.75
        sz = int(np.ceil(radius * 3600 / pixscale))
        coW = coH = 2*sz + 1
        targetwcs = Tan(ra, dec, 1. +coW/2, 1.+coH/2, -pixscale/3600.,0.,0.,pixscale/3600, coW,coH)

        A = fits_table('allsky-atlas.fits')
        A.cut(np.abs(A.dec - dec) < 0.9 + radius)
        A.cut(np.abs(A.ra  - ra ) < (0.9 + radius)/cosdec)
        print len(A), 'tiles nearby'

        WW,SS = [],[]
        for tile,tilera,tiledec in zip(A.coadd_id, A.ra, A.dec):
            allwise = True
            W,S = read_sources(allwise, tile, tdir, pdir, read_if_missing=False)
            if W is None or S is None:
                continue
            print len(W), 'AllWISE'
            print len(S), 'SDSS'
            I = AT.cut_to_unique(tilera, tiledec, W.ra, W.dec)
            if len(I) == 0:
                continue
            W.cut(I)
            I = AT.cut_to_unique(tilera, tiledec, S.ra, S.dec)
            if len(I) == 0:
                continue
            S.cut(I)
            W.cut((np.abs(W.ra - ra) < (radius / cosdec)) * (np.abs(W.dec - dec) < radius))
            S.cut((np.abs(S.ra - ra) < (radius / cosdec)) * (np.abs(S.dec - dec) < radius))
            WW.append(W)
            SS.append(S)
        W = merge_tables(WW)
        S = merge_tables(SS)
        del WW
        del SS
        # HACK -- arcsec_between, maybe, goes weird when N=1.
        if len(W) < 2 or len(S) < 2:
            continue

        if statsOnly:

            print 'len W:', len(W)
            print 'len S:', len(S)

            brbin = np.flatnonzero((w1mpro >= brightbins[:-1]) * (w1mpro < brightbins[1:]))
            if len(brbin) == 0:
                continue
            brbin = brbin[0]
            brightbin = (brightbins[brbin], brightbins[brbin+1])
            print 'Bright star bin:', brightbin

            minmag, maxmag = 10, 20
            magstep = 0.1
            minrad, maxrad = 0, radius * 3600.
            radstep = 10.

            magbins = np.arange(minmag, maxmag+magstep*0.01, magstep)
            radbins = np.arange(minrad, maxrad+radstep*0.01, radstep)
            #print 'magbins', magbins
            #print 'radbins', radbins

            #for band in [1,2]:
            for band in [1]:
                magname = 'W%i' % band

                r1 = arcsec_between(W.ra, W.dec, ra, dec)
                mag1 = W.get('w%impro' % band)
                meth1 = 'AllWISE'

                r2 = arcsec_between(S.ra, S.dec, ra, dec)
                mag2 = S.get('w%i_mag' % band)
                meth2 = 'Forced Photometry'

                for jj,(r,mag,meth) in enumerate([(r1,mag1,meth1), (r2,mag2,meth2)]):
                    print 'len r:', len(r)
                    print 'len mag:', len(mag)
                    print 'N radius bins:', len(radbins)
                    print 'N mag bins:', len(magbins)
                    print 'r', r.shape
                    print 'r', np.clip(mag, minmag, maxmag).shape
                    H,nil,nil = np.histogram2d(r, np.clip(mag, minmag, maxmag),
                                               bins=(radbins,magbins))
                    # print 'N radius bins:', len(radbins)
                    # print 'N mag bins:', len(magbins)
                    # print 'Raw H shape', H.shape
                    H = H.T
                    print 'H.T shape', H.shape
                    if len(H) == 2:
                        continue
                    key = (meth, magname, brbin)
                    if key in hists:
                        hh = hists[key]
                        hists[key] = (hh[0]+H, hh[1]+1) + hh[2:]
                    else:
                        hists[key] = H,1,radbins,magbins, brightbin

                    if jj == 1:
                        plt.clf()
                        weight = 1. / (np.pi * (radbins[1:]**2 - radbins[:-1]**2) / 3600.**2)
                        dmag = magbins[1] - magbins[0]
                        dens = H * weight[np.newaxis,:] / dmag
                        mx = 12000 # stars / sq deg / mag
                        mlo,mhi = magbins.min(), magbins.max()
                        rlo,rhi = radbins.min(), radbins.max()
                        dimshow(dens, cmap='hot', vmax=mx, vmin=0, extent=[rlo,rhi, mlo,mhi])
                        plt.gca().set_aspect('auto')
                        plt.xlabel('Radius from bright star (arcsec)')
                        plt.ylabel('Mag')
                        plt.ylim(mhi, mlo)
                        plt.title('%s %s: star (%.4f, %.4f)' %
                                  (meth, band, ra, dec))
                        plt.colorbar()
                        ps1.savefig()


            #plt.clf()
            #sn = W.w1flux / W.w1sigflux
            #r = arcsec_between(W.ra, W.dec, ra, dec)
            #plt.plot(r, np.clip(sn, minsn, maxsn), 'k.', alpha=0.2)
            #plt.xlabel('Radius (arcsec)')
            #plt.ylabel('S/N')
            #plt.title('AllWISE W1 S/N')
            #plt.ylim(minsn, maxsn)
            #ps.savefig()
            #
            #plt.clf()
            #sn = S.w1_nanomaggies * np.sqrt(S.w1_nanomaggies_ivar)
            #r = arcsec_between(S.ra, S.dec, ra, dec)
            #plt.plot(r, np.clip(sn, minsn, maxsn), 'k.', alpha=0.2)
            ##plt.plot(r, sn, 'k.')
            #plt.xlabel('Radius (arcsec)')
            #plt.ylabel('S/N')
            #plt.title('Tractor W1 S/N')
            #plt.ylim(minsn, maxsn)
            #ps.savefig()

            continue

        sdss = DR10()
        sdss.useLocalTree()
        sdss.saveUnzippedFiles('data/unzip')
    
        wfn = 'resolve-2010-05-23-cut/window_flist.fits'
        #wfn = 'photoResolve-new/window_flist.fits'
        RCF = radec_to_sdss_rcf(ra, dec, radius=np.hypot(np.sqrt(2.)*radius*60., np.hypot(14.,10.)/2.), tablefn=wfn)
        print len(RCF), 'SDSS fields nearby'
        RCF = [(run,camcol,field) for run,camcol,field,nil,nil in RCF
               if sdss.get_rerun(run, field) == '301']
        print len(RCF), 'in rerun 301'
        sco = np.zeros((coH,coW), np.float32)
        sn  = np.zeros((coH,coW), np.int32)
        band = 'r'
        for run,camcol,field in RCF:
            print 'RCF', run, camcol, field
            fn = sdss.retrieve('frame', run, camcol, field, band)
            frame = sdss.readFrame(run, camcol, field, band)
            h,w = frame.getImageShape()
            wcs = AsTransWrapper(frame.getAsTrans(), w, h)
            try:
                Yo,Xo,Yi,Xi,nil = resample_with_wcs(targetwcs, wcs, [], 3)
            except:
                continue
            if len(Yo) == 0:
                continue
            img = frame.getImage()
            # smooth the SDSS image a bit before resampling
            img = gaussian_filter(img, 4.)
            sco[Yo,Xo] += img[Yi,Xi]
            sn [Yo,Xo] += 1
        sco /= np.maximum(sn, 1)

        # Skip stars that have low SDSS coverage
        if (np.sum(sn == 0) > 0.1 * coH * coW):
            print 'Skipping: low SDSS coverage'
            continue

        w1co = np.zeros((coH,coW), np.float32)
        w1n  = np.zeros((coH,coW), np.int32)
        for a in A:
            tile = a.coadd_id
            print 'tile', tile
            fn = os.path.join(unwdir, tile[:3], tile, 'unwise-%s-w1-img-u.fits' % tile)
            wcs = Tan(fn)
            try:
                Yo,Xo,Yi,Xi,nil = resample_with_wcs(targetwcs, wcs, [], 3)
            except:
                continue
            if len(Yo) == 0:
                continue
            img = fitsio.read(fn)
            w1co[Yo,Xo] += img[Yi,Xi]
            w1n [Yo,Xo] += 1
        w1co /= np.maximum(w1n, 1)
    
        plt.clf()
        mn,mx = -8., 500.
        dimshow(w1co, vmin=mn, vmax=mx, ticks=False)
        ax = plt.axis()
        #plt.title('W1 mpro = %.3f' % w1mpro)
        ps.savefig()

        sca = dict(s=20, alpha=0.6, linewidths=0.5)

        def plotsn(ra, dec, sn):
            ok,x,y = targetwcs.radec2pixelxy(ra, dec)
            x -= 1
            y -= 1
            #I = (sn < 10)
            #p1 = plt.plot(x[I], y[I], 'bo', mew=0.5, alpha=0.5, ms=5)
            #p1 = plt.plot(x[I], y[I], 'bo', mew=0.5, alpha=0.5, ms=2)
            I = ((sn >= 10) * (sn < 30))
            #p2 = plt.plot(x[I], y[I], 'go', mew=0.5, alpha=0.5, ms=7)
            p2 = plt.plot(x[I], y[I], 'go', mew=0.5, alpha=0.5, ms=5)
            I = ((sn >= 30) * (sn < 100))
            #p3 = plt.plot(x[I], y[I], 'o', color=(1.,0.66,0.), mew=0.5, alpha=0.5, ms=10)
            p3 = plt.plot(x[I], y[I], 'o', color=(1.,0.66,0.), mew=0.5, alpha=0.5, ms=7)
            I = (sn > 100)
            #p4 = plt.plot(x[I], y[I], 'r*', mew=0.5, alpha=0.5, ms=15)
            p4 = plt.plot(x[I], y[I], 'r*', mew=0.5, alpha=0.5, ms=10)
            #plt.legend([p1[0],p2[0],p3[0],p4[0]], ('S/N < 10', 'S/N in [10,30)', 'S/N in [30,100)', 'S/N > 100'))
            return ([p2[0],p3[0],p4[0]], ('S/N 10 to 30', 'S/N 30 to 100', 'S/N $>$ 100'))

        def plotmag(ra, dec, mag):
            ok,x,y = targetwcs.radec2pixelxy(ra, dec)
            x -= 1
            y -= 1
            #I = (mag > 16)
            #p1 = plt.plot(x[I], y[I], 'bo', mew=0.5, alpha=0.5, ms=2)
            I = ((mag > 14) * (mag <= 16))
            p2 = plt.plot(x[I], y[I], 'go', mew=0.5, alpha=0.5, ms=5)
            I = ((mag > 12) * (mag <= 14))
            p3 = plt.plot(x[I], y[I], 'o', color=(1.,0.66,0.), mew=0.5, alpha=0.5, ms=7)
            I = (mag < 12)
            p4 = plt.plot(x[I], y[I], 'r*', mew=0.5, alpha=0.5, ms=10)
            # return ([p1[0],p2[0],p3[0],p4[0]],
            #         ('W1 $>$ 16', 'W1 14 to 16', 'W1 12 to 14', 'W1 $<$ 12'))
            return ([p2[0],p3[0],p4[0]],
                    ('W1 mag 14 to 16', 'W1 mag 12 to 14', 'W1 mag $<$ 12'))

        plt.clf()
        dimshow(w1co, vmin=mn, vmax=mx, cmap=lightgray, ticks=False)
        ax = plt.axis()
        #plt.title('W1 mpro = %.3f' % w1mpro)
        sn = W.w1flux / W.w1sigflux
        # ok,x,y = targetwcs.radec2pixelxy(W.ra, W.dec)
        # plt.scatter(x-1, y-1, c=np.round(W.w1mpro).astype(int), vmin=8, vmax=18, **sca)
        # plt.colorbar()
        #plotsn(W.ra, W.dec, sn)
        plotmag(W.ra, W.dec, W.w1mpro)
        plt.axis(ax)
        ps.savefig()
    
        plt.clf()
        mn,mx = 0., 0.25
        dimshow(sco, vmin=mn, vmax=mx, ticks=False)
        ax = plt.axis()
        ps.savefig()

        plt.clf()
        mn,mx = 0., 0.25
        dimshow(sco, vmin=mn, vmax=mx, cmap=lightgray, ticks=False)
        ax = plt.axis()
        sn = S.w1_nanomaggies * np.sqrt(S.w1_nanomaggies_ivar)
        #plt.scatter(x-1, y-1, c=np.round(sn).astype(int), vmin=0, vmax=20, **sca)
        # ok,x,y = targetwcs.radec2pixelxy(S.ra, S.dec)
        # plt.scatter(x-1, y-1, c=np.round(S.w1_mag).astype(int), vmin=8, vmax=18, **sca)
        # plt.colorbar()
        #lp,lt = plotsn(S.ra, S.dec, sn)
        print 'S w1_mag', S.w1_mag.min(), S.w1_mag.max()
        J = np.flatnonzero((S.w1_mag > -5) * (S.w1_mag < 30))
        S.cut(J)
        print 'S w1_mag', S.w1_mag.min(), S.w1_mag.max()
        lp,lt = plotmag(S.ra, S.dec, S.w1_mag)
        plt.figlegend(lp, lt, 'upper right', prop=dict(size=16))

        plt.axis(ax)
        ps.savefig()

    pickle_to_file(dict(hists=hists, brightbins=brightbins), 'bright-hists.pickle')



def density_hist_plots(ps, hists=None, whists=None, **kwargs):
    plt.figure(figsize=(5,4))
    plt.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95, hspace=0.05, wspace=0.05)

    small = 0.5
    for i in range(4):
        hi = 2*i
        plt.clf()
        hist,bins,name,dolog = hists[hi]
        xlo,xhi = min(bins), max(bins)
        xhi = 25
        p1 = plt.plot(np.repeat(bins, 2)[1:-1], np.repeat(np.maximum(hist, small), 2), 'b-')
        hist,bins,name,dolog = whists[hi]
        p2 = plt.plot(np.repeat(bins, 2)[1:-1], np.repeat(np.maximum(hist, small), 2), 'r-', lw=3, alpha=0.4)
        #plt.xlabel(name)
        plt.ylabel('Number of sources')
        if i == 0:
            plt.figlegend((p1[0], p2[0]), ('Forced photometry', 'AllWISE'), 'upper left')
        plt.xlabel('W%i mag' % (i+1))
        plt.xlim(xlo, xhi)
        plt.yscale('log', nonposy='clip')
        #plt.ylim(small, 2. * max(hist))
        plt.ylim(small, 1e8)
        ps.savefig()


    for i in range(4):
        hi = 2*i + 1
        plt.clf()
        hist,bins,name,dolog = hists[hi]
        p1 = plt.plot(np.repeat(bins, 2)[1:-1], np.repeat(np.maximum(hist, small), 2), 'b-')

        # HACK
        ylo,yhi = plt.ylim()

        hist,bins,name,dolog = whists[hi]
        p2 = plt.plot(np.repeat(bins, 2)[1:-1], np.repeat(np.maximum(hist, small), 2), 'r-', lw=3, alpha=0.4)

        plt.ylim(ylo,yhi)

        plt.ylabel('Number of sources')
        if i == 0:
            plt.figlegend((p1[0], p2[0]), ('Forced photometry', 'AllWISE'), 'upper right')
        plt.xlabel('W%i flux signal-to-noise' % (i+1))
        #plt.xlim(min(bins), max(bins))
        plt.xlim(-5, 15)
        #plt.ylim(0, 1.5e7)
        plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ps.savefig()
        
    return

    # small = 0.5
    # for hi in range(len(hists)):
    #     hist,bins,name,dolog = hists[hi]
    #     plt.clf()
    #     plt.plot(np.repeat(bins, 2)[1:-1], np.repeat(np.maximum(hist, small), 2), 'k-')
    #     plt.xlabel(name)
    #     plt.xlim(min(bins), max(bins))
    #     if dolog:
    #         plt.yscale('log', nonposy='clip')
    #         plt.ylim(small, 2. * max(hist))
    #         ps.savefig()
    #     else:
    #         ps.savefig()
    # 
    # for hi in range(len(whists)):
    #     hist,bins,name,dolog = whists[hi]
    #     plt.clf()
    #     plt.plot(np.repeat(bins, 2)[1:-1], np.repeat(np.maximum(hist, small), 2), 'k-')
    #     plt.xlabel(name)
    #     plt.xlim(min(bins), max(bins))
    #     if dolog:
    #         plt.yscale('log', nonposy='clip')
    #         plt.ylim(small, 2. * max(hist))
    #         ps.savefig()
    #     else:
    #         ps.savefig()


def density_map_plots(ps, mapnames=None, maps=None,
                      mapralo=None, maprahi=None, mapdeclo=None, mapdechi=None,
                      mapdra=None, mapddec=None, maxes={}, ticks={},
                      outline=None,
                      **kwargs):

    print 'Map RA', mapralo, maprahi
    print 'Map Dec', mapdeclo, mapdechi

    #plt.figure(figsize=(12,4))
    #plt.subplots_adjust(left=0.15, right=0.99, bottom=0.1, top=0.95)#, hspace=0.05, wspace=0.05)
    #plt.subplots_adjust(left=0.05, right=0.99, bottom=0.05, top=0.97)#, hspace=0.05, wspace=0.05)
    #plt.figure(figsize=(11,3))
    #plt.subplots_adjust(left=0.07, right=0.99, bottom=0.15, top=0.99)#, hspace=0.05, wspace=0.05)

    plt.figure(figsize=(9.5,3))
    #plt.subplots_adjust(left=0.07, right=0.99, bottom=0.16, top=0.99)
    plt.subplots_adjust(left=0.07, right=0.99, bottom=0.15, top=0.99)

    # for i,(name,themap) in enumerate(zip(mapnames, maps)):
    #     hdr = fitsio.FITSHDR()
    #     hdr.add_record(dict(name='MAPRALO', value=mapralo))
    #     hdr.add_record(dict(name='MAPRAHI', value=maprahi))
    #     hdr.add_record(dict(name='MAPDRA', value=mapdra))
    #     hdr.add_record(dict(name='MAPDECLO', value=mapdeclo))
    #     hdr.add_record(dict(name='MAPDECHI', value=mapdechi))
    #     hdr.add_record(dict(name='MAPDDEC', value=mapddec))
    #     hdr.add_record(dict(name='MAPNAME', value=name))
    #     fitsio.write('density-map-%i.fits' % i, themap, header=hdr, clobber=True)

    rabins  = np.arange(mapralo , maprahi +0.01*mapdra,  mapdra )
    decbins = np.arange(mapdeclo, mapdechi+0.01*mapddec, mapddec)
    print 'RA bins', len(rabins)
    print 'Dec bins', len(decbins)

    dec = (decbins[:-1] + decbins[1:]) / 2.
    cosdec = np.cos(np.deg2rad(dec))

    for name,themap in zip(mapnames, maps):
        print
        print 'Map:', name
        #print 'Map shape', themap.shape
        print 'Map sum:', themap.sum()
        print 'Map non-zero area:', ((themap > 0) * cosdec[:,np.newaxis] * mapdra * mapddec).sum(), 'sq deg'

        plt.clf()
        mapx = themap / cosdec[:,np.newaxis] / (mapdra * mapddec)

        mx = np.percentile(mapx, 99.5)
        if name in maxes:
            mx = maxes[name]

        # if outline is not None:
        #     rgb = antigray(np.clip((mapx - 0) / (mx - 0), 0., 1.))
        #     print 'rgb', rgb.shape, rgb.dtype
        #     rgb[:,:,0][outline] = 0
        #     rgb[:,:,1][outline] = 1.
        #     rgb[:,:,2][outline] = 0
        #     plim = dimshow(rgb, extent=[mapralo, maprahi, mapdeclo, mapdechi])
        # else:
        #     plim = dimshow(mapx, cmap=antigray,
        #                    extent=[mapralo, maprahi, mapdeclo, mapdechi], vmin=0, vmax=mx)
        #plt.title(name)
        plim = dimshow(mapx, cmap=antigray,
                       extent=[mapralo, maprahi, mapdeclo, mapdechi], vmin=0, vmax=mx)
        if outline is not None:
            H,W = mapx.shape
            rgba = np.zeros((H,W,4), np.uint8)
            rgba[:,:,1] = rgba[:,:,3] = outline * 255
            plt.imshow(rgba, extent=[mapralo,maprahi,mapdeclo,mapdechi],
                       interpolation='nearest', origin='lower')

        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        plt.xlim(maprahi, mapralo)
        tt = ticks.get(name, None)
        #plt.colorbar(fraction=0.07, pad=0.03, aspect=15, ticks=tt)

        # from dstn's thesis:plotstyle.py:colorbar_axes
        frac = 0.07
        pad  = 0.02
        aspect = 20
        parent = plt.gca()
        pb = parent.get_position(original=True).frozen()
        # new parent box, padding, child box
        (pbnew, padbox, cbox) = pb.splitx(1.0-(frac+pad), 1.0-frac)
        cbox = cbox.anchored('C', cbox)
        parent.set_position(pbnew)
        parent.set_anchor((1.0, 0.5))
        cax = parent.get_figure().add_axes(cbox)
        cax.set_aspect(aspect, anchor=((0.0, 0.5)), adjustable='box')
        parent.get_figure().sca(parent)

        c = plt.colorbar(plim, cax=cax, ticks=tt)

        ps.savefig()


def raiseitup(x):
    raise x

def density_pobjs(T, tdir, pdir, podir):
    mapdra,mapddec = 0.1, 0.1
    # Whole sky
    mapralo, maprahi = 0., 360.
    mapdeclo, mapdechi = -26., 86.

    rabins  = np.arange(mapralo , maprahi +0.01*mapdra,  mapdra )
    decbins = np.arange(mapdeclo, mapdechi+0.01*mapddec, mapddec)

    mapW = int(np.ceil((maprahi  - mapralo ) / mapdra ))
    mapH = int(np.ceil((mapdechi - mapdeclo) / mapddec))

    pomaps = {}
    pomapnames = {}

    for dirpath, dirnames, fns in os.walk(podir, followlinks=True, onerror=raiseitup):
        dirnames.sort()
        fns.sort()
        for fn in fns:
            pth = os.path.join(dirpath, fn)
            print pth
            T = fits_table(pth)

            cuts = [('photoObjs: has\_wise\_phot', T.has_wise_phot),
                    ]
            for i in range(1,5):
                npix = (T.get('w%i_npix' % i) > 0)
                cuts.append(('W%i npix $>$ 0' % i, npix))

            for i,(name,cut) in enumerate(cuts):
                if not i in pomaps:
                    pomaps[i] = np.zeros((mapH, mapW), np.int32)
                    pomapnames[i] = name
                r = T.ra[cut]
                if len(r) == 0:
                    continue
                an_hist2d(r, T.dec[cut], pomaps[i], mapralo, maprahi, mapdeclo, mapdechi)
                #H,nil,nil = np.histogram2d(r, T.dec[cut], bins=(rabins,decbins))
                #pomaps[i] += H.T
    # dict -> list
    pomaps      = [pomaps     [i] for i in range(len(pomaps     ))]
    pomapnames  = [pomapnames [i] for i in range(len(pomapnames ))]

    X = dict(pomaps=pomaps, pomapnames=pomapnames,
             mapralo=mapralo, maprahi=maprahi, mapdra=mapdra,
             mapdeclo=mapdeclo, mapdechi=mapdechi, mapddec=mapddec)
    pickle_to_file(X, 'data/wphot/density-pomaps-2.pickle')

            

def density(T, tdir, pdir):
    from unwise_coadd import tile_to_radec
    ps = PlotSequence('density')
    #ps.suffixes = ['png', 'pdf']
    plt.figure(figsize=(8,5))
    plt.subplots_adjust(left=0.15, right=0.99, bottom=0.1, top=0.95, hspace=0.05, wspace=0.05)

    fn = 'density-sdss-dr10.pickle'
    if False and os.path.exists(fn):
        X = unpickle_from_file(fn)
        #density_hist_plots(ps, **X)
        ps = PlotSequence('map')
        density_map_plots(ps, **X)
        return

    #mapdra,mapddec = 0.5, 0.5
    mapdra,mapddec = 0.1, 0.1

    # Whole sky
    mapralo, maprahi = 0., 360.
    mapdeclo, mapdechi = -26., 86.

    # 20 tiles
    #mapralo,  maprahi  = 0, 6.
    #mapdeclo, mapdechi = -10., -5.

    # ~200 tiles
    #mapralo,  maprahi  = 10., 30.
    #mapdeclo, mapdechi = 10., 30.

    # ~1000 tiles
    #mapralo,  maprahi  = 150., 200.
    #mapdeclo, mapdechi =  5.,   55.

    rabins  = np.arange(mapralo , maprahi +0.01*mapdra,  mapdra )
    decbins = np.arange(mapdeclo, mapdechi+0.01*mapddec, mapddec)

    T.cut((T.ra > mapralo-1) * (T.ra < maprahi+1) *
          (T.dec > mapdeclo-1) * (T.dec < mapdechi+1))

    mapW = int(np.ceil((maprahi  - mapralo ) / mapdra ))
    mapH = int(np.ceil((mapdechi - mapdeclo) / mapddec))
    print 'map size', mapW, mapH
    maps = {}
    mapnames = {}
    wmaps = {}
    wmapnames = {}

    assert(len(rabins)  == mapW+1)
    assert(len(decbins) == mapH+1)

    hists = {}
    whists = {}

    bright = []
    wbright = []

    A = AtlasTiles()

    for it,(tile,tilera,tiledec) in enumerate(zip(T.coadd_id, T.ra, T.dec)):
        print
        print 'tile', (it+1), 'of', len(T)
        allwise = True
        try:
            #W,S = read_sources(allwise, tile, tdir, pdir, match_to_sdss=False,
            #                   read_if_missing=False)

            sfn = os.path.join(pdir, 'phot-%s.fits' % tile)
            print 'Reading', sfn
            S = fits_table(sfn, columns=['ra','dec','objc_type','treated_as_pointsource', 'pointsource',
                                         'w1_nanomaggies', 'w1_nanomaggies_ivar', 'w1_npix', 'w1_mag',
                                         'w2_nanomaggies', 'w2_nanomaggies_ivar', 'w2_npix', 'w2_mag',
                                         'w3_nanomaggies', 'w3_nanomaggies_ivar', 'w3_npix', 'w3_mag',
                                         'w4_nanomaggies', 'w4_nanomaggies_ivar', 'w4_npix', 'w4_mag',
                                         ])
            S.star = np.logical_or(S.objc_type == 6, S.treated_as_pointsource == 1)
            S.gal = np.logical_not(S.star)
            # We updated objc_type based on treated_as_pointsource, but recorded the
            # original objc_type in .pointsource.
            #S.sdss_star = (S.objc_type == 6)
            S.sdss_star = (S.pointsource == 1)
            S.sdss_gal  = np.logical_not(S.pointsource)
            S.w1_snr = S.w1_nanomaggies * np.sqrt(S.w1_nanomaggies_ivar)
            S.w2_snr = S.w2_nanomaggies * np.sqrt(S.w2_nanomaggies_ivar)
            S.w3_snr = S.w3_nanomaggies * np.sqrt(S.w3_nanomaggies_ivar)
            S.w4_snr = S.w4_nanomaggies * np.sqrt(S.w4_nanomaggies_ivar)

            W = None

        except:
            print 'Failed to read_sources:'
            import traceback
            traceback.print_exc()
            continue

        if S is None:
            continue

        if S is not None:
            print 'Read', len(S), 'SDSS'
            I = A.cut_to_unique(tilera, tiledec, S.ra, S.dec)
            if len(I) == 0:
                continue
            S.cut(I)
    
            cuts = [
                ('All SDSS', np.ones((len(S)), bool)),
                ('SDSS stars', S.sdss_star),
                ('SDSS galaxies', S.sdss_gal),
                ('Treated as point sources', S.star),
                ('Treated as galaxies', S.gal),
                ('W1 npix > 0', S.w1_npix > 0),
                ('W2 npix > 0', S.w2_npix > 0),
                ('W3 npix > 0', S.w3_npix > 0),
                ('W4 npix > 0', S.w4_npix > 0),
                ('W1 invvar > 0', S.w1_nanomaggies_ivar > 0),
                ('W2 invvar > 0', S.w2_nanomaggies_ivar > 0),
                ('W3 invvar > 0', S.w3_nanomaggies_ivar > 0),
                ('W4 invvar > 0', S.w4_nanomaggies_ivar > 0),
                ]
            for i in range(1,5):
                #good = (S.get('w%i_npix' % i) > 0)
                good = (S.get('w%i_nanomaggies_ivar' % i) > 0)
                snr = S.get('w%i_snr' % i)
                cuts.append(('W%i SNR $\leq$ 0' % i, np.flatnonzero(good * (snr <= 0))))
                cuts.append(('W%i SNR in (0,1]' % i,
                             np.flatnonzero(good * (snr > 0) * (snr <= 1))))
                cuts.append(('W%i SNR in (1,2]' % i,
                             np.flatnonzero(good * (snr > 1) * (snr <= 2))))
                cuts.append(('W%i SNR in (2,3]' % i,
                             np.flatnonzero(good * (snr > 2) * (snr <= 3))))
                cuts.append(('W%i SNR in (3,4]' % i,
                             np.flatnonzero(good * (snr > 3) * (snr <= 4))))
                cuts.append(('W%i SNR in (4,5]' % i,
                             np.flatnonzero(good * (snr > 4) * (snr <= 5))))
                cuts.append(('W%i SNR in (5,10]' % i,
                             np.flatnonzero(good * (snr > 5) * (snr <= 10))))
                cuts.append(('W%i SNR $>$ 10' % i, np.flatnonzero(good * (snr > 10))))

            for i,(name,cut) in enumerate(cuts):
                # instantiate even if empty...
                if not i in maps:
                    maps[i] = np.zeros((mapH, mapW), np.int32)
                    mapnames[i] = name
                r = S.ra[cut]
                if len(r) == 0:
                    continue
                an_hist2d(r, S.dec[cut], maps[i], mapralo, maprahi, mapdeclo, mapdechi)
                #H,nil,nil = np.histogram2d(r, S.dec[cut], bins=(rabins,decbins))
                #H,nil,nil = np.histogram2d(r, S.dec[cut], range=((mapralo,maprahi),(mapdeclo,mapdechi)), bins=(mapW,mapH))
                #maps[i] += H.T

            # Grab bright stars
            I = np.flatnonzero(reduce(np.logical_or, [
                #S.w1_mag <= 12., S.w2_mag <= 12., S.w3_mag <= 11., S.w4_mag <= 7.]))
                S.w1_mag <= 9., S.w2_mag <= 9., S.w3_mag <= 8., S.w4_mag <= 4.]))
            if len(I):
                bb = S[I]
                bb.coadd_id = np.array([tile] * len(bb))
                bright.append(bb)

            hi = 0
            for i in range(1,5):
                good = np.flatnonzero(S.get('w%i_nanomaggies_ivar' % i) > 0)
                #I = np.flatnonzero(S.get('w%i_npix' % i) > 0)

                mag = S.get('w%i_mag' % i)[good]
                bins = np.arange(-5., 30.001, 0.1)
                hist,be = np.histogram(mag, bins=bins)
                if not hi in hists:
                    name = 'Forced Photometry W%i mag' % i
                    hists[hi] = (hist, bins, name, True)
                else:
                    h,nil,nil,nil = hists[hi]
                    h += hist
                    #print 'This hist', sum(hist), 'total', sum(h)
                hi += 1

                #flux = S.get('w%i_nanomaggies'      % i)[I]
                #ivar = S.get('w%i_nanomaggies_ivar' % i)[I]
                #snr = flux * np.sqrt(ivar)
                snr = S.get('w%i_snr' % i)[good]
                bins = np.arange(-5., 15.001, 0.1)
                hist,be = np.histogram(snr, bins=bins)
                if not hi in hists:
                    name = 'Forced Photometry W%i flux S/N' % i
                    hists[hi] = (hist, bins, name, False)
                else:
                    h,nil,nil,nil = hists[hi]
                    h += hist
                hi += 1

                
        if W is not None:
            print 'Read', len(W), 'AllWISE'
            I = A.cut_to_unique(tilera, tiledec, W.ra, W.dec)
            print len(I), 'within tile'
            if len(I) == 0:
                continue
            W.cut(I)
            W.bx = np.floor((W.ra  - mapralo ) / mapdra ).astype(int)
            W.by = np.floor((W.dec - mapdeclo) / mapddec).astype(int)
            W.cut((W.bx >= 0) * (W.bx < mapW) * (W.by >= 0) * (W.by < mapH))
            print 'Cut to', len(W), 'within map'
    
            wcuts = [
                ('AllWISE', np.ones((len(W)), bool)),
                ]
            for i in range(1,5):
                flux = W.get('w%iflux' % i)
                sigflux = W.get('w%isigflux' % i)
                snr = flux / sigflux
                wcuts.append(('AllWISE W%i SNR $\leq$ 0' % i, np.flatnonzero(snr <= 0)))
                wcuts.append(('AllWISE W%i SNR in (0,1]' % i,
                             np.flatnonzero((snr > 0) * (snr <= 1))))
                wcuts.append(('AllWISE W%i SNR in (1,2]' % i,
                             np.flatnonzero((snr > 1) * (snr <= 2))))
                wcuts.append(('AllWISE W%i SNR in (2,3]' % i,
                             np.flatnonzero((snr > 2) * (snr <= 3))))
                wcuts.append(('AllWISE W%i SNR in (3,4]' % i,
                             np.flatnonzero((snr > 3) * (snr <= 4))))
                wcuts.append(('AllWISE W%i SNR in (4,5]' % i,
                             np.flatnonzero((snr > 4) * (snr <= 5))))
                wcuts.append(('AllWISE W%i SNR in (5,10]' % i,
                             np.flatnonzero((snr > 5) * (snr <= 10))))
                wcuts.append(('AllWISE W%i SNR $>$ 10' % i, np.flatnonzero(snr > 10)))
            
            for i,(name,cut) in enumerate(wcuts):
                if not i in wmaps:
                    wmaps[i] = np.zeros((mapH, mapW), np.int32)
                    wmapnames[i] = name
                r = W.ra[cut]
                if len(r) == 0:
                    continue
                H,nil,nil = np.histogram2d(r, W.dec[cut], bins=(rabins,decbins))
                wmaps[i] += H.T
    
            I = np.flatnonzero(reduce(np.logical_or, [
                #W.w1mpro <= 12., W.w2mpro <= 12., W.w3mpro <= 11., W.w4mpro <= 7.]))
                W.w1mpro <= 9., W.w2mpro <= 9., W.w3mpro <= 8., W.w4mpro <= 4.]))
            if len(I):
                bb = W[I]
                bb.coadd_id = np.array([tile] * len(bb))
                wbright.append(bb)
    
            hi = 0
            for i in range(1,5):
                mag = W.get('w%impro' % i)
                bins = np.arange(-5, 30.001, 0.1)
                hist,be = np.histogram(mag, bins=bins)
                if not hi in whists:
                    name = 'AllWISE W%i mpro (mag)' % i
                    whists[hi] = (hist, bins, name, True)
                else:
                    h,nil,nil,nil = whists[hi]
                    h += hist
                hi += 1
    
                flux = W.get('w%iflux' % i)
                sigflux = W.get('w%isigflux' % i)
                snr = flux / sigflux
                bins = np.arange(-5, 15.001, 0.1)
                hist,be = np.histogram(snr, bins=bins)
                if not hi in whists:
                    name = 'AllWISE W%i flux S/N' % i
                    whists[hi] = (hist, bins, name, False)
                else:
                    h,nil,nil,nil = whists[hi]
                    h += hist
                hi += 1


    if len(bright):
        B = merge_tables(bright)
        print len(B), 'bright'
        B.writeto('data/wphot/bright-2.fits')
    if len(wbright):
        B = merge_tables(wbright)
        print len(B), 'AllWISE bright'
        for col in ['w1proflux', 'w1prodflux', 'w1prosnr', 'bx', 'by']:
            B.delete_column(col)
        B.writeto('data/wphot/allwise-bright-2.fits')

    # dict -> list
    if len(maps):
        maps      = [maps     [i].astype(np.int32) for i in range(len(maps     ))]
        mapnames  = [mapnames [i] for i in range(len(mapnames ))]
        X = dict(maps=maps, mapnames=mapnames,
                 mapralo=mapralo, maprahi=maprahi, mapdra=mapdra,
                 mapdeclo=mapdeclo, mapdechi=mapdechi, mapddec=mapddec)
        pickle_to_file(X, 'data/wphot/density-maps-2.pickle')

        X = dict(hists=hists)
        pickle_to_file(X, 'data/wphot/density-hists-2.pickle')


    if len(wmaps):
        wmaps     = [wmaps    [i] for i in range(len(wmaps    ))]
        wmapnames = [wmapnames[i] for i in range(len(wmapnames))]
        X = dict(wmaps=wmaps, wmapnames=wmapnames,
                 mapralo=mapralo, maprahi=maprahi, mapdra=mapdra,
                 mapdeclo=mapdeclo, mapdechi=mapdechi, mapddec=mapddec)
        pickle_to_file(X, 'data/wphot/density-wmaps-2.pickle')

        X = dict(whists=whists)
        pickle_to_file(X, 'data/wphot/density-whists-2.pickle')

    #mapnames = [name for name,cut in cuts + wcuts]
    #density_plots(ps, **X)
    # X = dict(maps=maps, mapnames=mapnames, hists=hists, whists=whists,
    #          mapralo=mapralo, maprahi=maprahi, mapdra=mapdra,
    #          mapdeclo=mapdeclo, mapdechi=mapdechi, mapddec=mapddec)
    # pickle_to_file(X, 'data/wphot/density.pickle')




def wdensity(T, tdir, pdir, tag, mask=None):
    #mapdra,mapddec = 0.5, 0.5
    mapdra,mapddec = 0.1, 0.1

    # Whole sky
    mapralo, maprahi = 0., 360.
    mapdeclo, mapdechi = -26., 86.

    # 20 tiles
    #mapralo,  maprahi  = 0, 6.
    #mapdeclo, mapdechi = -10., -5.

    # ~200 tiles
    #mapralo,  maprahi  = 10., 30.
    #mapdeclo, mapdechi = 10., 30.

    # ~1000 tiles
    #mapralo,  maprahi  = 150., 200.
    #mapdeclo, mapdechi =  5.,   55.

    rabins  = np.arange(mapralo , maprahi +0.01*mapdra,  mapdra )
    decbins = np.arange(mapdeclo, mapdechi+0.01*mapddec, mapddec)

    T.cut((T.ra > mapralo-1) * (T.ra < maprahi+1) *
          (T.dec > mapdeclo-1) * (T.dec < mapdechi+1))

    mapW = int(np.ceil((maprahi  - mapralo ) / mapdra ))
    mapH = int(np.ceil((mapdechi - mapdeclo) / mapddec))
    print 'map size', mapW, mapH
    assert(len(rabins)  == mapW+1)
    assert(len(decbins) == mapH+1)

    wmaps = {}
    wmapnames = {}
    whists = {}
    wbright = []

    A = AtlasTiles()

    for it,(tile,tilera,tiledec) in enumerate(zip(T.coadd_id, T.ra, T.dec)):
        print
        print 'tile', (it+1), 'of', len(T)
        fn = os.path.join(tdir, 'wise-sources-%s.fits' % tile)
        print fn
        if not os.path.exists(fn):
            from unwise_coadd import get_coadd_tile_wcs
            from wise.allwisecat import allwise_catalog_wcs
            wcs = get_coadd_tile_wcs(tilera, tiledec)
            W = allwise_catalog_wcs(wcs)
            print 'Read', len(W), 'AllWISE'
            W.writeto(fn)
            print 'Wrote', fn
        else:
            W = fits_table(fn, columns=[
                'ra', 'dec',
                'w1mpro', 'w2mpro', 'w3mpro', 'w4mpro',
                'w1sigmpro', 'w2sigmpro', 'w3sigmpro', 'w4sigmpro',
                'w1flux', 'w2flux', 'w3flux', 'w4flux',
                'w1sigflux', 'w2sigflux', 'w3sigflux', 'w4sigflux',])
                #'w1gmag', 'w1gerr', 'w2gmag', 'w2gerr', 'w3gmag',
                #'w3gerr', 'w4gmag', 'w4gerr', 'bx', 'by',
                #'coadd_id']
            print 'Read', len(W), 'AllWISE'

        I = A.cut_to_unique(tilera, tiledec, W.ra, W.dec)
        print len(I), 'within tile'
        if len(I) == 0:
            continue
        W.cut(I)
        W.bx = np.floor((W.ra  - mapralo ) / mapdra ).astype(int)
        W.by = np.floor((W.dec - mapdeclo) / mapddec).astype(int)
        W.cut((W.bx >= 0) * (W.bx < mapW) * (W.by >= 0) * (W.by < mapH))
        print 'Cut to', len(W), 'within map'

        if mask is not None:
            W.cut(mask[W.by, W.bx])
            print 'Cut to', len(W), 'within mask'

        wcuts = [
            ('AllWISE', np.ones((len(W)), bool)),
            ]
        for i in range(1,5):
            flux = W.get('w%iflux' % i)
            sigflux = W.get('w%isigflux' % i)
            snr = flux / sigflux
            wcuts.append(('AllWISE W%i SNR $\leq$ 0' % i, np.flatnonzero(snr <= 0)))
            wcuts.append(('AllWISE W%i SNR in (0,1]' % i,
                         np.flatnonzero((snr > 0) * (snr <= 1))))
            wcuts.append(('AllWISE W%i SNR in (1,2]' % i,
                         np.flatnonzero((snr > 1) * (snr <= 2))))
            wcuts.append(('AllWISE W%i SNR in (2,3]' % i,
                         np.flatnonzero((snr > 2) * (snr <= 3))))
            wcuts.append(('AllWISE W%i SNR in (3,4]' % i,
                         np.flatnonzero((snr > 3) * (snr <= 4))))
            wcuts.append(('AllWISE W%i SNR in (4,5]' % i,
                         np.flatnonzero((snr > 4) * (snr <= 5))))
            wcuts.append(('AllWISE W%i SNR in (5,10]' % i,
                         np.flatnonzero((snr > 5) * (snr <= 10))))
            wcuts.append(('AllWISE W%i SNR $>$ 10' % i, np.flatnonzero(snr > 10)))
        
        for i,(name,cut) in enumerate(wcuts):
            if not i in wmaps:
                wmaps[i] = np.zeros((mapH, mapW), np.int32)
                wmapnames[i] = name
            r = W.ra[cut]
            if len(r) == 0:
                continue
            an_hist2d(r, W.dec[cut], wmaps[i], mapralo, maprahi, mapdeclo, mapdechi)

        I = np.flatnonzero(reduce(np.logical_or, [
            #W.w1mpro <= 12., W.w2mpro <= 12., W.w3mpro <= 11., W.w4mpro <= 7.]))
            W.w1mpro <= 9., W.w2mpro <= 9., W.w3mpro <= 8., W.w4mpro <= 4.]))
        if len(I):
            bb = W[I]
            bb.coadd_id = np.array([tile] * len(bb))
            wbright.append(bb)

        hi = 0
        for i in range(1,5):
            mag = W.get('w%impro' % i)
            bins = np.arange(-5, 30.001, 0.1)
            hist,be = np.histogram(mag, bins=bins)
            if not hi in whists:
                name = 'AllWISE W%i mpro (mag)' % i
                whists[hi] = (hist, bins, name, True)
            else:
                h,nil,nil,nil = whists[hi]
                h += hist
            hi += 1

            flux = W.get('w%iflux' % i)
            sigflux = W.get('w%isigflux' % i)
            snr = flux / sigflux
            bins = np.arange(-5, 15.001, 0.1)
            hist,be = np.histogram(snr, bins=bins)
            if not hi in whists:
                name = 'AllWISE W%i flux S/N' % i
                whists[hi] = (hist, bins, name, False)
            else:
                h,nil,nil,nil = whists[hi]
                h += hist
            hi += 1

    # dict -> list
    if len(wmaps):
        wmaps     = [wmaps    [i] for i in range(len(wmaps    ))]
        wmapnames = [wmapnames[i] for i in range(len(wmapnames))]
        X = dict(wmaps=wmaps, wmapnames=wmapnames,
                 mapralo=mapralo, maprahi=maprahi, mapdra=mapdra,
                 mapdeclo=mapdeclo, mapdechi=mapdechi, mapddec=mapddec)
        pickle_to_file(X, 'data/wphot/density-wmaps%s.pickle' % tag)

        X = dict(whists=whists)
        pickle_to_file(X, 'data/wphot/density-whists%s.pickle' % tag)

    if len(wbright):
        B = merge_tables(wbright)
        print len(B), 'AllWISE bright'
        for col in ['bx', 'by']:
            B.delete_column(col)
        B.writeto('data/wphot/allwise-bright%s.fits' % tag)






def profiles(T, tdir, pdir, unwdir):
    ps = PlotSequence('profiles')
    ps.suffixes = ['png', 'pdf']

    plt.figure(figsize=(4,4))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=0.05, wspace=0.05)

    sdss = DR10()
    sdss.useLocalTree()
    #sdss.setFitsioReadBZ2()
    sdss.saveUnzippedFiles('data/unzip')

    for tile in T.coadd_id:
        fn = os.path.join(tdir, 'photoobjs-%s.fits' % tile)
        print 'Reading', fn
        cols = ['objid', 'ra', 'dec', 'fracdev', 'objc_type', 'modelflux',
                'theta_dev', 'theta_deverr', 'ab_dev', 'ab_deverr', 'phi_dev_deg',
                'theta_exp', 'theta_experr', 'ab_exp', 'ab_experr', 'phi_exp_deg',
                'resolve_status', 'nchild', 'flags', 'objc_flags',
                'run','camcol','field','id'
                ]
        # useful to have these in the outputs...
        cols += ['psfflux', 'psfflux_ivar', 'cmodelflux', 'cmodelflux_ivar',
                 'modelflux', 'modelflux_ivar']
        S = fits_table(fn, columns=cols)
        print 'Got', len(S)

        fn = os.path.join(unwdir, tile[:3], tile, 'unwise-%s-w1-img-u.fits' % tile)
        print 'unWISE', fn
        wiseimg = fitsio.read(fn)
        print 'Image', wiseimg.shape
        wisewcs = Tan(fn)

        fn = os.path.join(unwdir, tile[:3], tile, 'unwise-%s-w2-img-u.fits' % tile)
        print 'unWISE', fn
        wiseimg2 = fitsio.read(fn)

        # cut to in-bounds
        ok,wx,wy = wisewcs.radec2pixelxy(S.ra, S.dec)
        margin = 10
        h,w = wiseimg.shape
        S.cut((wx > margin) * (wy > margin) * (wx <= w-margin) * (wy <= h-margin))
    
        # r band
        b = bandnum = 2
        S.treat = treat_as_pointsource(S, bandnum)
        S.cut(S.objc_type == 3)
        print len(S), 'galaxies'

        S.dev = (S.fracdev[:,b] >= 0.5)
        S.exp = (S.fracdev[:,b] <  0.5)

        S.theta = np.zeros(len(S))
        S.theta[S.dev] = S.theta_dev[S.dev,b]
        S.theta[S.exp] = S.theta_exp[S.exp,b]

        Ibig = np.argsort(-S.theta)

        # skip the really big ones...
        Ibig = Ibig[100:]

        pixscale = 0.396 / 3600.
        #sz = 40
        sz = 60
        rows,cols = 5,5

        allwstamps = []
        allnwstamps = []
        allsstamps = []
        allrawsdss = []

        frameCache = {}
        for j in range(10):
            stamps,I,wcses,rawsdss = get_sdss_cutouts(
                sdss, S[Ibig], pixscale, sz, rows*cols,
                frameCache=frameCache, get_rawvals=True)
            plt.clf()
            for i,stamp in enumerate(stamps):
                plt.subplot(rows, cols, i+1)
                dimshow(stamp)
                s = S[Ibig[I[i]]]
                plt.title('%s %.3g' % (np.where(s.dev, 'dev', 'exp'), s.theta))
                plt.xticks([]); plt.yticks([])
            ps.savefig()

            allsstamps.append(stamps)
            allrawsdss.append(rawsdss)
            
            wstamps = []
            nwstamps = []
            for i,(targetwcs) in enumerate(wcses):
                Yo,Xo,Yi,Xi,[r1,r2] = resample_with_wcs(targetwcs, wisewcs,
                                                        [wiseimg,wiseimg2], 3)
                w,h = targetwcs.imagew, targetwcs.imageh
                rwise = np.zeros((h,w,3), np.float32)
                rwise[Yo,Xo,2] = r1
                rwise[Yo,Xo,0] = r2
                rwise[:,:,1] = (rwise[:,:,0] + rwise[:,:,2])/2.
                wstamps.append(rwise)
                rwise = np.zeros((h,w,3), np.float32)
                rwise[Yo,Xo,2] = wiseimg[Yi,Xi]
                rwise[Yo,Xo,0] = wiseimg2[Yi,Xi]
                rwise[:,:,1] = (rwise[:,:,0] + rwise[:,:,2])/2.
                nwstamps.append(rwise)
            allwstamps.append(wstamps)
            allnwstamps.append(nwstamps)

            # imw = dict(vmin=-5, vmax=500, cmap='gray')
            # plt.clf()
            # for i,wise in enumerate(wstamps):
            #     plt.subplot(rows, cols, i+1)
            #     #dimshow(wise, **imw)
            #     dimshow(np.clip(wise / 250, 0., 1.))
            #     print 'flux range', wise.min(), wise.max()
            #     plt.xticks([]); plt.yticks([])
            # ps.savefig()

            for Wstamps in [nwstamps,
                            wstamps]:
                plt.clf()
                for i,wise in enumerate(Wstamps):
                    Q = 20
                    m2 = 0.
                    alpha = 0.003
                    m = -5.
                    r,g,b = wise[:,:,0], wise[:,:,1], wise[:,:,2]
                    r = np.maximum(0, r - m)
                    g = np.maximum(0, g - m)
                    b = np.maximum(0, b - m)
                    Iim = (r+g+b)/3.
                    fI = np.arcsinh(alpha * Q * (Iim - m2)) / np.sqrt(Q)
                    Iim += (Iim == 0.) * 1e-6
                    R = fI * r / Iim
                    G = fI * g / Iim
                    B = fI * b / Iim
                    RGB = np.dstack((R,G,B))
                    print 'min,max', RGB.min(), RGB.max()
                    print 'median:', np.median(RGB.ravel())
                    #mn,mx = 1.0, 2.0
                    mn,mx = 0., 1.
                    RGB = np.clip((RGB - mn) / (mx - mn), 0., 1.)
                    plt.subplot(rows, cols, i+1)
                    dimshow(RGB)
                    #plt.hist(RGB.ravel(), 50, range=(0,5))
                    plt.xticks([]); plt.yticks([])
                ps.savefig()

            Ibig = Ibig[I[-1]+1:]

        pickle_to_file(dict(allwstamps=allwstamps, allnwstamps=allnwstamps,
                            allsstamps=allsstamps, allrawsdss=allrawsdss),
                       'wstamps.pickle')

        break

def get_sdss_cutout(sdss, ra, dec, run, camcol, field, pixscale, sz, frameCache={},
                    get_rawvals=False):
    rgbims = []
    targetwcs = Tan(ra, dec, sz+1, sz+1, -pixscale, 0., 0., pixscale,
                    2.*sz+1, 2.*sz+1)
    for bandnum in [3,2,1]:
        key = (run, camcol, field, bandnum)
        frame = frameCache.get(key, None)
        if frame is None:
            frame = sdss.readFrame(run, camcol, field, bandnum)
            print 'Got frame', frame
            frameCache[key] = frame
        h,w = frame.getImageShape()
        x,y = frame.astrans.radec_to_pixel(ra, dec)
        x,y = int(x), int(y)
        # add some margin for resampling
        sz2 = sz + 5
        if x < sz2 or y < sz2 or x >= (w-sz2) or y >= (h-sz2):
            break
        stamp = frame.getImageSlice((slice(y-sz2, y+sz2+1), slice(x-sz2, x+sz2+1)))
        sh,sw = stamp.shape
        wcs = AsTransWrapper(frame.astrans, sw, sh, x0=x-sz2, y0=y-sz2)
        Yo,Xo,Yi,Xi,[rim] = resample_with_wcs(targetwcs, wcs, [stamp], 3)
        targetim = np.zeros((2*sz+1, 2*sz+1), np.float32)
        targetim[Yo,Xo] = rim
        rgbims.append(targetim)
    if len(rgbims) < 3:
        return None

    if get_rawvals:
        rawvals = [x.copy() for x in rgbims]

    r,g,b = rgbims
    # i
    r *= 1.0
    # r
    #g *= 1.5
    g *= 1.3
    # g
    b *= 2.5
    m = -0.02
    r = np.maximum(0, r - m)
    g = np.maximum(0, g - m)
    b = np.maximum(0, b - m)
    I = (r+g+b)/3.
    alpha = 1.5
    Q = 20
    m2 = 0.
    fI = np.arcsinh(alpha * Q * (I - m2)) / np.sqrt(Q)
    I += (I == 0.) * 1e-6
    R = fI * r / I
    G = fI * g / I
    B = fI * b / I
    maxrgb = reduce(np.maximum, [R,G,B])
    J = (maxrgb > 1.)
    R[J] = R[J]/maxrgb[J]
    G[J] = G[J]/maxrgb[J]
    B[J] = B[J]/maxrgb[J]
    ss = 0.5
    RGBblur = np.clip(np.dstack([
        gaussian_filter(R, ss),
        gaussian_filter(G, ss),
        gaussian_filter(B, ss)]), 0., 1.)

    rtn = [RGBblur, targetwcs]
    if get_rawvals:
        rtn.append(rawvals)
    return rtn

def get_sdss_cutouts(sdss, S, pixscale, sz, N, frameCache={},
                     get_rawvals=False):
    stamps = []
    wcses = []
    ii = []
    rawvals = []
    for i,s in enumerate(S):
        X = get_sdss_cutout(sdss, s.ra, s.dec, s.run, s.camcol, s.field,
                            pixscale, sz, frameCache=frameCache,
                            get_rawvals=get_rawvals)
        if X is None:
            continue
        if get_rawvals:
            RGBblur, targetwcs, rv = X
            rawvals.append(rv)
        else:
            RGBblur, targetwcs = X
        stamps.append(RGBblur)
        ii.append(i)
        wcses.append(targetwcs)
        if len(stamps) >= N:
            break

    rtn = [stamps, np.array(ii), wcses]
    if get_rawvals:
        rtn.append(rawvals)
    return rtn


def treated_as_pointsource(T, tdir, pdir):
    ps = PlotSequence('ptsrc')
    ps.suffixes = ['png', 'pdf']

    plt.figure(figsize=(4,4))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=0.05, wspace=0.05)

    sdss = DR10()
    sdss.useLocalTree()
    sdss.setFitsioReadBZ2()
    
    allstamps = []

    #for tile in T.coadd_id:
    # More than enough in the first tile...
    tile = T.coadd_id[0]
    fn = os.path.join(tdir, 'photoobjs-%s.fits' % tile)
    print 'Reading', fn
    cols = ['objid', 'ra', 'dec', 'fracdev', 'objc_type', 'modelflux',
            'theta_dev', 'theta_deverr', 'ab_dev', 'ab_deverr', 'phi_dev_deg',
            'theta_exp', 'theta_experr', 'ab_exp', 'ab_experr', 'phi_exp_deg',
            'resolve_status', 'nchild', 'flags', 'objc_flags',
            'run','camcol','field','id'
            ]
    # useful to have these in the outputs...
    cols += ['psfflux', 'psfflux_ivar', 'cmodelflux', 'cmodelflux_ivar',
             'modelflux', 'modelflux_ivar']
    S = fits_table(fn, columns=cols)
    print 'Got', len(S)

    # r band
    b = 2
    gal = (S.objc_type == 3)
    S.cut(gal)
    gal = (S.objc_type == 3)
    dev = gal * (S.fracdev[:,b] >= 0.5)
    exp = gal * (S.fracdev[:,b] <  0.5)
    print sum(dev), 'deV,', sum(exp), 'exp'
    print 'Total', len(S), 'gals'

    thetasn = np.zeros(len(S))
    S.theta_deverr[dev,b] = np.maximum(1e-6, S.theta_deverr[dev,b])
    S.theta_experr[exp,b] = np.maximum(1e-5, S.theta_experr[exp,b])
    # theta_experr nonzero: 1.28507e-05
    # theta_deverr nonzero: 1.92913e-06
    thetasn[dev] = S.theta_dev[dev,b] / S.theta_deverr[dev,b]
    thetasn[exp] = S.theta_exp[exp,b] / S.theta_experr[exp,b]

    aberrzero = np.zeros(len(S), bool)
    aberrzero[dev] = (S.ab_deverr[dev,b] == 0.)
    aberrzero[exp] = (S.ab_experr[exp,b] == 0.)

    maxtheta = np.zeros(len(S), bool)
    maxtheta[dev] = (S.theta_dev[dev,b] >= 29.5)
    maxtheta[exp] = (S.theta_exp[exp,b] >= 59.0)

    # theta S/N > modelflux for dev, 10*modelflux for exp
    bigthetasn = (thetasn > (S.modelflux[:,b] * (1.*dev + 10.*exp)))

    print sum(gal * (thetasn < 3.)), 'have low S/N in theta'
    print sum(gal * (S.modelflux[:,b] > 1e4)), 'have big flux'
    print sum(aberrzero), 'have zero a/b error'
    print sum(maxtheta), 'have the maximum theta'
    print sum(bigthetasn), 'have large theta S/N vs modelflux'

    bit_thetasn = 1
    bit_bigflux = 2
    bit_aberr   = 4
    bit_maxtheta = 8
    bit_thetaflux = 16
    
    S.badbits = (bit_thetasn * (thetasn < 3) +
                 bit_bigflux * (S.modelflux[:,b] > 1e4) +
                 bit_aberr   * aberrzero  +
                 bit_maxtheta * maxtheta +
                 bit_thetaflux * bigthetasn)

    S.cut(S.badbits > 0)
    print 'Cut to', len(S), 'treated at point source'

    bits = [bit_thetasn, bit_bigflux, bit_aberr, bit_maxtheta, bit_thetaflux]
    for bit in bits:
        Sb = S[S.badbits == bit]
        print len(Sb), 'with only bit', bit, 'set'
        rows,cols = 5,5
        #Sb.about()
        sz = 25

        pixscale = 0.396 / 3600.
        stamps,nil,nil = get_sdss_cutouts(sdss, Sb, pixscale, sz, rows*cols)

        plt.clf()
        for i,stamp in enumerate(stamps):
            plt.subplot(rows, cols, i)
            dimshow(stamp)
            plt.xticks([]); plt.yticks([])
        ps.savefig()

        allstamps.append(stamps)


    plt.figure(figsize=(8,4))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=0.05, wspace=0.05)
    rows,cols = len(allstamps),10
    sp = 1
    for row,stamps in enumerate(allstamps):
        for col,stamp in enumerate(stamps[:cols]):
            plt.subplot(rows, cols, sp)
            sp += 1
            dimshow(stamp)
            plt.xticks([]); plt.yticks([])
    ps.savefig()
        


def model_plots(F, S, UW, unwdir):
    ps = PlotSequence('model')
    ps.suffixes = ['png', 'pdf']

    plt.figure(figsize=(3,3))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)

    sdss = DR10()
    sdss.useLocalTree()
    
    I = np.flatnonzero((F.w1_mag > 16) * (F.w1_mag < 17))
    print 'Mag 16-17 objs:'
    for i,(r,w1,ra,dec,run,cc,field,tile) in enumerate(zip(
        F.r_psf[I], F.w1_mag[I], F.ra[I], F.dec[I], F.run[I], F.camcol[I], F.field[I],
        F.tile[I])[:18]):
        print '  %12.3f,  %12.3f,  %12.4f %12.4f' % (r,w1, ra,dec), 'RCF', run,cc,field, 'tile', tile
        if i in [0,1,2,8,11,12,13,14]:
            continue
        print 'index', i
        
        frame = sdss.readFrame(run, cc, field, 'r')
        img = frame.getImage()
        print 'Image', img.shape
        h,w = img.shape
        ast = frame.getAsTrans()
        wcs = AsTransWrapper(ast, w, h)
        x,y = ast.radec_to_pixel(ra, dec)
        print 'Pix', x,y

        pixscale = 0.4 / 3600
        stampsz = 150
        targetwcs = Tan(ra, dec, stampsz/2 + 0.5, stampsz/2 + 0.5,
                        -pixscale, 0., 0., pixscale, stampsz, stampsz)
        r0,r1,d0,d1 = targetwcs.radec_bounds()
        print 'RA,Dec bounds:', r0,r1,d0,d1
        margin = 10. * 2.75 / 3600.
        cosdec = np.cos(np.deg2rad(dec))
        r0 -= margin/cosdec
        r1 += margin/cosdec
        d0 -= margin
        d1 += margin

        rsdss = np.zeros((stampsz,stampsz))
        Yo,Xo,Yi,Xi,nil = resample_with_wcs(targetwcs, wcs, [], 3)
        rsdss[Yo,Xo] = img[Yi,Xi]

        fn = os.path.join(unwdir, tile[:3], tile, 'unwise-%s-w1-img-u.fits' % tile)
        print 'unWISE', fn
        wiseimg = fitsio.read(fn)
        print 'Image', wiseimg.shape
        wisewcs = Tan(fn)
        fn = os.path.join(unwdir, tile[:3], tile, 'unwise-%s-w1-invvar-u.fits.gz' % tile)
        print 'unWISE', fn
        wiseiv = fitsio.read(fn)
        print 'Invvar', wiseiv.shape
        
        Yo,Xo,Yi,Xi,[limg] = resample_with_wcs(targetwcs, wisewcs, [wiseimg], 3)
        rwise = np.zeros((stampsz,stampsz))
        rwise[Yo,Xo] = limg
        rinverr = np.zeros((stampsz,stampsz))
        rinverr[Yo,Xo] = np.sqrt(wiseiv[Yi,Xi])

        plt.clf()
        #plt.subplot(2,3,4)
        plt.imshow(rsdss, interpolation='nearest', origin='lower',
                   vmin=-0.025, vmax=0.25, cmap='gray')
        #plt.title('SDSS image')
        plt.xticks([]); plt.yticks([])
        ax = plt.axis()
        J = np.flatnonzero((S.ra > r0) * (S.ra < r1) * (S.dec > d0) * (S.dec < d1))
        ok,sx,sy = targetwcs.radec2pixelxy(S.ra[J], S.dec[J])
        ptsrc = np.logical_or(S.pointsource[J] == 1, S.treated_as_pointsource[J] == 1)
        extsrc = np.logical_not(ptsrc)
        galcc = (1., 0.65, 0)
        plt.plot(sx[ptsrc]-1, sy[ptsrc]-1, 'o', mec='r', mfc='none', ms=12)
        plt.plot(sx[extsrc]-1, sy[extsrc]-1, 'o', mec=galcc, mfc='none', ms=16)
        plt.axis(ax)
        ps.savefig()

        SS = S[J]
        sband = 'r'
        wanyband = wband = 'w'
        SS.cmodelflux = np.ones((len(SS), 5))
        SS.psfflux = SS.cmodelflux
        SS.devflux = SS.cmodelflux
        SS.expflux = SS.cmodelflux
        cat = get_tractor_sources_dr9(None, None, None, bandname=sband, objs=SS,
                                      bands=[], nanomaggies=True, extrabands=[wband],
                                      fixedComposites=True, useObjcType=True,
                                      objCuts=False)
        print 'Created', len(cat), 'sources'
        psffn = 'psf-allwise-con3.fits'
        band = 1
        P = fits_table(psffn, hdu=band)
        psf = GaussianMixturePSF(P.amp, P.mean, P.var)

        J = np.flatnonzero((UW.ra > r0) * (UW.ra < r1) * (UW.dec > d0) * (UW.dec < d1))
        wcat = [PointSource(RaDecPos(UW.ra[i], UW.dec[i]), NanoMaggies(w=UW.flux[i]))
                for i in J]

        modscale = 2.75 / 3600
        modsz = int(np.ceil(stampsz * pixscale / modscale)) + 1
        modwcs = Tan(ra, dec, modsz/2 + 0.5, modsz/2 + 0.5,
                     -modscale, 0., 0., modscale, modsz, modsz)

        tim = Image(data=np.zeros((modsz,modsz)), invvar=np.ones((modsz,modsz)),
                    psf=psf, photocal=LinearPhotoCal(1., band=wband),
                    wcs=ConstantFitsWcs(modwcs))

        tractor = Tractor([tim], cat)
        tractor.freezeAllRecursive()
        tractor.thawPathsTo(wband)
        tractor.setParams(SS.w1_nanomaggies)
        tractor = Tractor([tim], cat + wcat)
        mod = tractor.getModelImage(0)

        Yo,Xo,Yi,Xi,[limg] = resample_with_wcs(targetwcs, modwcs, [mod], 3)
        rmod = np.zeros((stampsz,stampsz))
        rmod[Yo,Xo] = limg

        imw = dict(interpolation='nearest', origin='lower',
                   vmin=-5, vmax=50, cmap='gray')
        
        #plt.subplot(2,3,2)
        plt.clf()
        plt.imshow(rwise, **imw)
        plt.xticks([]); plt.yticks([])
        #plt.title('WISE image')
        ps.savefig()

        #plt.subplot(2,3,1)
        plt.clf()
        plt.imshow(rwise, **imw)
        plt.xticks([]); plt.yticks([])
        ax = plt.axis()
        J = np.flatnonzero((W.ra > r0) * (W.ra < r1) * (W.dec > d0) * (W.dec < d1))
        ok,x,y = targetwcs.radec2pixelxy(W.ra[J], W.dec[J])
        plt.plot(x-1, y-1, 'g+', ms=12, mew=2)
        plt.plot(sx[ptsrc]-1, sy[ptsrc]-1, 'o', mec='r', mfc='none', ms=6)
        plt.plot(sx[extsrc]-1, sy[extsrc]-1, 'o', mec=galcc, mfc='none', ms=8)
        plt.axis(ax)
        #plt.title('WISE image')
        ps.savefig()

        # imchi = dict(interpolation='nearest', origin='lower',
        #              vmin=-5, vmax=5, cmap='RdBu')
        # 
        # plt.subplot(2,2,3)
        # plt.imshow(-(rwise - rmod) * rinverr, **imchi)

        #plt.subplot(2,3,5)
        plt.clf()
        plt.imshow(rmod, **imw)
        plt.xticks([]); plt.yticks([])
        #plt.title('WISE model')
        ps.savefig()

        Yo,Xo,Yi,Xi,[limg] = resample_with_wcs(modwcs, wisewcs, [wiseimg], 3)
        mimg = np.zeros((modsz,modsz))
        mimg[Yo,Xo] = limg
        merr = np.zeros((modsz,modsz))
        merr[Yo,Xo] = 1. / np.sqrt(wiseiv[Yi,Xi])

        noise = np.random.normal(size=mod.shape) * merr
        #plt.subplot(2,3,6)
        plt.clf()
        plt.imshow(mod + noise, **imw)
        plt.xticks([]); plt.yticks([])
        #plt.title('WISE model + noise')
        ps.savefig()

        #plt.subplot(2,3,3)
        plt.clf()
        plt.imshow(mimg, **imw)
        plt.xticks([]); plt.yticks([])
        #plt.title('WISE image')
        ps.savefig()


def lrg_comparison(tdir, pdir, T, ps):
    # Matching vs forced photometry
    ffn = 'F.fits'
    mfn = 'M.fits'
    if os.path.exists(ffn) and os.path.exists(mfn):
        F = fits_table(ffn)
        M = fits_table(mfn)
    else:
        allwise = True
        MM,FF = [],[]
        for tile in T.coadd_id:
            W,S = read_sources(allwise, tile, tdir, pdir)
            # unique matches only
            I,J,d,counts = match_radec(W.ra, W.dec, S.ra, S.dec, 4./3600.,
                                       nearest=True, count=True)
            K = np.flatnonzero(counts == 1)
            I = I[K]
            J = J[K]
    
            M = W[I]
            SM = S[J]
            for c in SM.get_columns():
                if any([c.startswith('w%i_' % i) for i in [1,2,3,4]]):
                    SM.delete_column(c)
            M.add_columns_from(SM)
            MM.append(M)
    
            FF.append(S)
        M = merge_tables(MM)
        F = merge_tables(FF)
    
        print len(M), 'matched'
        print len(F), 'forced-phot'
        # Remove nearby pairs in both sets -- to de-dup margins
        I,J,d = match_radec(F.ra, F.dec, F.ra, F.dec, 1./3600, notself=True)
        K = np.ones(len(F), bool)
        K[np.minimum(I, J)] = False
        F.cut(K)
        print len(F), 'forced-phot after removing self-matches'
    
        I,J,d = match_radec(M.ra, M.dec, M.ra, M.dec, 1./3600, notself=True)
        K = np.ones(len(M), bool)
        K[np.minimum(I, J)] = False
        M.cut(K)
        print len(M), 'matched after removing self-matches'

        M.writeto(mfn)
        F.writeto(ffn)

    print len(M), 'matched'
    print len(F), 'forced-phot'

    # plt.clf()
    # plt.subplot(1,2,1)
    # plothist(M.ra, M.dec, 200, doclf=False)
    # plt.subplot(1,2,2)
    # plothist(F.ra, F.dec, 200, doclf=False)
    # ps.savefig()

    print 'Forced:'
    F.cut(F.dr_psf < 0.25)
    print 'Cut to', len(F), 'with r-mag err < 0.25'
    F.cut(F.di_psf < 0.25)
    print 'Cut to', len(F), 'with i-mag err < 0.25'

    print 'Matched:'
    M.cut(M.dr_psf < 0.25)
    print 'Cut to', len(M), 'with r-mag err < 0.25'
    M.cut(M.di_psf < 0.25)
    print 'Cut to', len(M), 'with i-mag err < 0.25'

    M.w1_mag = M.w1mpro

    plt.figure(1, figsize=(6,3.8))
    #plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.93)
    plt.subplots_adjust(left=0.08, right=0.99, bottom=0.12, top=0.96)
    plt.figure(2, figsize=(4,3.8))
    #plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.93)
    #plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.97)
    plt.subplots_adjust(left=0.1, right=0.96, bottom=0.12, top=0.96)

    # LRG
    mags = [13, #16,
            17, #17.5,
            18., 18.5, ]

    ### HACK
    #mags = [0.]
    #[(13,19)] + 

    fontsize = 16

    for mi,(mlo,mhi) in enumerate(zip(mags, mags[1:])):
        FI = np.flatnonzero((F.w1_mag > mlo) * (F.w1_mag < mhi))
        MI = np.flatnonzero((M.w1_mag > mlo) * (M.w1_mag < mhi))

        Mri = M.r_mag[MI] - M.i_mag[MI]
        Mrw = M.r_mag[MI] - M.w1_mag[MI]
        Fri = F.r_mag[FI] - F.i_mag[FI]
        Frw = F.r_mag[FI] - F.w1_mag[FI]

        plt.figure(1)
        plt.clf()
        #ha = dict(bins=100, range=((-0.25,2.25),(1,7)), doclf=False,
        ha = dict(bins=100, range=((-0.5,2.25),(1,7)), doclf=False,
                  imshowargs=dict(cmap=antigray), docolorbar=False, hot=False)
        mx = 0
        if len(Mri):
            plt.subplot(1,2,1)
            H,xe,ye = loghist(Mri, Mrw, **ha)
            mx = max(mx, H.max())
            plt.xticks([0,1,2])
            plt.xlabel('r - i (mag)', fontsize=fontsize)
            plt.ylabel('r - W1 (mag)', fontsize=fontsize)
            #plt.title('Matched')
            #plt.title('W1 in [%g,%g] (matched)' % (mlo,mhi))
        if len(Fri):
            plt.subplot(1,2,2)
            H,xe,ye = loghist(Fri, Frw, **ha)
            mx = max(mx, H.max())
            plt.xticks([0,1,2])
            plt.xlabel('r - i (mag)', fontsize=fontsize)
            #plt.title('Forced phot')
            #plt.title('W1 in [%g,%g] (forced phot)' % (mlo,mhi))
        #plt.suptitle('W1 in [%g,%g]' % (mlo,mhi))
        for sp in [1,2]:
            plt.subplot(1,2,sp)
            ii = plt.gci()
            if ii is None:
                continue
            ii.set_clim(0.3, np.log10(mx))

        ri0 = 1.
        rw0 = 4.
        angle = -np.arctan(2.)
        uumin, uumax = -0.5, 2.0

        if True:#mhi == 17:
            boxv = np.array([-0.45, -0.45, 1.65, 1.65, -0.45])
            boxu = np.array([uumin, uumax, uumax, uumin, uumin])
            slu = np.array([uumin, uumax])
            slv = np.array([-0.1, -0.1])
            boxri = np.cos(-angle) * boxu - np.sin(-angle) * boxv + ri0
            boxrw = np.sin(-angle) * boxu + np.cos(-angle) * boxv + rw0
            slri = np.cos(-angle) * slu - np.sin(-angle) * slv + ri0
            slrw = np.sin(-angle) * slu + np.cos(-angle) * slv + rw0
            ax = plt.axis()
            plt.plot(boxri, boxrw, 'r-')
            plt.plot(slri, slrw, 'r--')
            print 'Box:', boxri, boxrw
            plt.axis(ax)
        
        ps.savefig()


        dri = Mri - ri0
        drw = Mrw - rw0
        Muu = np.cos(angle) * dri - np.sin(angle) * drw
        Mvv = np.sin(angle) * dri + np.cos(angle) * drw
        I = np.flatnonzero((Muu > uumin) * (Muu < uumax))
        Muu = Muu[I]
        Mvv = Mvv[I]

        dri = Fri - ri0
        drw = Frw - rw0
        Fuu = np.cos(angle) * dri - np.sin(angle) * drw
        Fvv = np.sin(angle) * dri + np.cos(angle) * drw
        I = np.flatnonzero((Fuu > uumin) * (Fuu < uumax))
        Fuu = Fuu[I]
        Fvv = Fvv[I]

        plt.figure(2)
        plt.clf()
        ll,lt = [],[]
        ha = dict(bins=50, range=(-0.3, 1.8), histtype='step')
        if len(Mvv):
            n,b,p = plt.hist(Mvv+0.1, lw=3, color=(1,0.5,0.5), **ha)
            ll.append(p[0])
            lt.append('Matched')
        if len(Fvv):
            n,b,p = plt.hist(Fvv+0.1, color='b', **ha)
            ll.append(p[0])
            lt.append('Forced')
        plt.xlim(-0.3, 1.9)
        plt.xlabel('Distance from stellar locus (mag)', fontsize=fontsize)
        #plt.figlegend(ll, lt, 'upper right')
        #plt.legend(reversed(ll), reversed(lt), 'upper right')
        if mi == 0:
            plt.figlegend(reversed(ll), reversed(lt), 'upper right',
                          prop=dict(size=fontsize))
        #plt.title('W1 in [%g,%g]' % (mlo,mhi))

        
        ps.savefig()

    return

    
    ### QSO-like selection

    ps = PlotSequence('qso')
    ps.suffixes = ['png', 'pdf']

    imags = [13, 20, 21, 22, 23]
    wmags = [13, 16, 17, 18., 18.5, 19]#, 19.5, 20]

    for mlo,mhi,cutband in (#[(13,22,'i')] + zip(imags, imags[1:], ['i']*len(imags)) +
                            [(13,19,'w1')] + zip(wmags, wmags[1:], ['w1']*len(wmags))):
        #FI = np.flatnonzero((F.r_mag > mlo) * (F.r_mag < mhi) * F.star)
        #MI = np.flatnonzero((M.r_mag > mlo) * (M.r_mag < mhi) * M.star)

        #FI = np.flatnonzero((F.i_mag > mlo) * (F.i_mag < mhi) * F.star)
        #MI = np.flatnonzero((M.i_mag > mlo) * (M.i_mag < mhi) * M.star)
        mag = F.get('%s_mag' % cutband)
        FI = np.flatnonzero((mag > mlo) * (mag < mhi) * F.star)
        mag = M.get('%s_mag' % cutband)
        MI = np.flatnonzero((mag > mlo) * (mag < mhi) * M.star)

        Mgi = M.g_mag[MI] - M.i_mag[MI]
        Miw = M.i_mag[MI] - M.w1_mag[MI]
        Fgi = F.g_mag[FI] - F.i_mag[FI]
        Fiw = F.i_mag[FI] - F.w1_mag[FI]

        plt.figure(1)
        plt.clf()
        ha1 = dict(bins=100, range=((0,4),(0,6)), doclf=False,
                  imshowargs=dict(cmap=antigray), docolorbar=False, hot=False)
        mx = 0
        if len(Mri):
            plt.subplot(1,2,1)
            H,xe,ye = loghist(Mgi, Miw, **ha1)
            mx = max(mx, H.max())
            #plt.xticks([0,1,2])
            plt.ylabel('i - W1')
            plt.xlabel('g - i')
            plt.title('Matched')
        if len(Fri):
            plt.subplot(1,2,2)
            H,xe,ye = loghist(Fgi, Fiw, **ha1)
            mx = max(mx, H.max())
            #plt.xticks([0,1,2])
            plt.xlabel('g - i')
            plt.title('Forced phot')
        plt.suptitle('%s in [%g,%g]' % (cutband,mlo,mhi))
        for sp in [1,2]:
            plt.subplot(1,2,sp)
            ii = plt.gci()
            if ii is None:
                continue
            ii.set_clim(0.3, np.log10(mx))

        x0 = 1.0
        y0 = 1.8
        slope = 1.3 / 1.8
        angle = -np.arctan(slope)
        uumin, uumax = 0., 4.0
        vvmin, vvmax = -0.5, 3.0

        boxv = np.array([vvmin, vvmin, vvmax, vvmax, vvmin])
        boxu = np.array([uumin, uumax, uumax, uumin, uumin])
        slu = np.array([uumin, uumax])
        slv = np.array([0., 0.])
        boxri = np.cos(-angle) * boxu - np.sin(-angle) * boxv + x0
        boxrw = np.sin(-angle) * boxu + np.cos(-angle) * boxv + y0
        slri = np.cos(-angle) * slu - np.sin(-angle) * slv + x0
        slrw = np.sin(-angle) * slu + np.cos(-angle) * slv + y0
        ax = plt.axis()
        plt.plot(boxri, boxrw, 'r-')
        plt.plot(slri, slrw, 'r--')
        plt.axis(ax)

        # for x in [1,1.1,1.2,1.3,1.4,1.5, 2.5,2.6,2.7,2.8,2.9,3.0]:
        #     plt.axvline(x, color='r')
        # for y in [1.5,1.6,1.7,1.8,1.9,2.0, 3.0,3.1,3.2,3.3,3.4]:
        #     plt.axhline(y, color='r')

        # # Plot a few error bars
        # x = np.arange(5)
        # y = np.arange(1,6)
        # xx,yy = np.meshgrid(x,y)
        # for x,y in zip(xx.ravel(), yy.ravel()):
        #     #J = np.argmin(np.hypot(Fgi - x, Fiw - y) + 1e12*np.logical_not(np.logical_and(np.isfinite(Fgi), np.isfinite(Fiw))))
        #     #J = np.argmin(np.where(np.logical_and(np.isfinite(Fgi), np.isfinite(Fiw)),
        #     #                       np.hypot(Fgi - x, Fiw - y), 1e12))
        #     J = np.flatnonzero((np.abs(Fgi - x) < 0.5) * (np.abs(Fiw - y) < 0.5))
        #     if len(J) == 0:
        #         continue
        #     dg = np.median(F.dg_mag[FI][J])
        #     di = np.median(F.di_mag[FI][J])
        #     dw = np.median(F.w1_mag_err[FI][J])
        #     #print 'Plotting error bar:', Fgi[J], Fiw[J], 'dg', dg, 'di', di, 'dw', dw
        #     # plt.errorbar([Fgi[J]], [Fiw[J]], xerr=dg, fmt='o', color='g')
        #     # plt.errorbar([Fgi[J]], [Fiw[J]], xerr=di, fmt='o', color='r')
        #     # plt.errorbar([Fgi[J]], [Fiw[J]], yerr=dw, fmt='o', color='m')
        #     plt.errorbar([x], [y], xerr=dg, fmt='o', color='g')
        #     plt.errorbar([x], [y], xerr=di, fmt='o', color='r')
        #     plt.errorbar([x], [y], yerr=dw, fmt='o', color='m')

        ps.savefig()

         
        dx = Mgi - x0
        dy = Miw - y0
        Muu = np.cos(angle) * dx - np.sin(angle) * dy
        Mvv = np.sin(angle) * dx + np.cos(angle) * dy
        I = np.flatnonzero((Muu > uumin) * (Muu < uumax))
        Muu = Muu[I]
        Mvv = Mvv[I]
        
        dx = Fgi - x0
        dy = Fiw - y0
        Fuu = np.cos(angle) * dx - np.sin(angle) * dy
        Fvv = np.sin(angle) * dx + np.cos(angle) * dy
        I = np.flatnonzero((Fuu > uumin) * (Fuu < uumax))
        Fuu = Fuu[I]
        Fvv = Fvv[I]
        
        plt.figure(2)
        plt.clf()
        ll,lt = [],[]
        ha = dict(bins=50, range=(vvmin, vvmax), histtype='step')
        if len(Mvv):
            n,b,p = plt.hist(Mvv+0.1, lw=3, color=(1,0.5,0.5), **ha)
            ll.append(p[0])
            lt.append('Matched')
        if len(Fvv):
            n,b,p = plt.hist(Fvv+0.1, color='b', **ha)
            ll.append(p[0])
            lt.append('Forced')
        plt.xlim(vvmin, vvmax)
        plt.xlabel('Distance from stellar locus (mag)')
        plt.legend(reversed(ll), reversed(lt), 'upper right')
        ps.savefig()


        #Fsnr = F.w1_snr[FI]
        #Msnr = M.w1prosnr[MI]

        # predicted stellar-locus W1 given g,i
        Fw1pred = F.i_mag[FI] - (y0 + (Fgi - x0) * slope)
        Mw1pred = M.i_mag[MI] - (y0 + (Mgi - x0) * slope)

        Fchi = (F.w1_mag[FI] - Fw1pred) / F.w1_mag_err[FI]
        Mchi = (M.w1mpro[MI] - Mw1pred) / M.w1sigmpro[MI]

        Mchi = Mchi[np.isfinite(Mchi)]
        Fchi = Fchi[np.isfinite(Fchi)]
        Mchi *= -1
        Fchi *= -1

        plt.figure(2)
        plt.clf()
        ll,lt = [],[]
        #ha = dict(bins=50, range=(-3,10), histtype='step')
        ha = dict(bins=50, histtype='step')

        xx = np.append(Mchi, Fchi)
        lo,hi = [np.percentile(xx,p) for p in [2,98]]
        ha.update(range=(lo,hi))

        if len(Mchi):
            n,b,p = plt.hist(Mchi, lw=3, color=(1,0.5,0.5), **ha)
            ll.append(p[0])
            lt.append('Matched')
        if len(Fchi):
            n,b,p = plt.hist(Fchi, color='b', **ha)
            ll.append(p[0])
            lt.append('Forced')
        #plt.xlim(vvmin, vvmax)
        plt.xlim(lo,hi)

        plt.xlabel('Distance from stellar locus (sigmas)')
        plt.legend(reversed(ll), reversed(lt), 'upper right')
        ps.savefig()
        
def photoz(T, tdir, pdir):
    ps = PlotSequence('photoz')
    # z in [0.2, 0.5]
    #M = fits_table('zgals.fits')
    # Just RA,Dec cut, no z cut
    M = fits_table('gals2.fits')
    print len(M), 'galaxies'

    plt.clf()
    plt.plot(M.ra, M.dec, 'r.')
    plt.plot(T.ra, T.dec, 'bo')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    ps.savefig()

    allwise = True
    FF = []
    for tile in T.coadd_id:
        W,F = read_sources(allwise, tile, tdir, pdir)
        FF.append(F)
    F = merge_tables(FF)
    print len(F), 'forced-phot'
    F.cut(np.isfinite(F.r_mag))

    '''
    select ra,dec, run,camcol,field,obj, objid,
    psfMag_u, psfMag_g, psfMag_r, psfMag_i, psfMag_z,
    psfMagErr_u,psfMagErr_g,psfMagErr_r,psfMagErr_i,psfMagErr_z,
    modelmag_u, modelmag_g, modelmag_r, modelmag_i, modelmag_z,
    modelmagerr_u, modelmagerr_g, modelmagerr_r, modelmagerr_i, modelmagerr_z,
    fracdev_r, exprad_r, devrad_r
    into mydb.MyTable_7 from PhotoPrimary
    where ra between 172 and 188
    and dec between 40 and 44
    and (((fracdev_r > 0.5) and (devrad_r between 0.5 and 1.5))
    or ((fracdev_r <= 0.5) and (exprad_r between 0.5 and 1.5)))

    # from Galaxy -> biggals2.fits
    '''
    B = fits_table('biggals2.fits')

    I,J,d = match_radec(B.ra, B.dec, F.ra, F.dec, 1./3600.)
    print len(I), 'matches to medium-sized galaxies'

    BB = B[I]
    FF = F[J]
    K = np.isfinite(FF.w1_mag)
    print sum(K), 'with good W1'
    BB.cut(K)
    FF.cut(K)

    K = np.flatnonzero((FF.w1_mag_err < 1) * (FF.dg_mag < 1) * (FF.di_mag < 1))
    FF.cut(K)
    BB.cut(K)
    c1 = FF.i_mag - FF.w1_mag
    c2 = FF.g_mag - FF.i_mag

    ax = [1, 3.2, 2.5, 5.0]

    plt.clf()
    plt.scatter(c2, c1, c=BB.fracdev_r, edgecolors='none', alpha=0.5)
    plt.xlabel('g - i')
    plt.ylabel('i - W1')
    plt.colorbar()
    plt.axis(ax)
    plt.title('SDSS galaxies $r_e$ in $[0.5,1.5]$" (color: fracdev\_r)')
    ps.savefig()

    K = np.flatnonzero(np.logical_or(BB.fracdev_r < 0.1, BB.fracdev_r > 0.9))
    plt.clf()
    plt.scatter(c2[K], c1[K], c=BB.fracdev_r[K], edgecolors='none', alpha=0.5)
    plt.xlabel('g - i')
    plt.ylabel('i - W1')
    plt.colorbar()
    plt.axis(ax)
    plt.title('SDSS galaxies $r_e$ in $[0.5,1.5]$", extreme fracdev (color: fracdev\_r)')
    ps.savefig()

    plt.clf()
    plt.hist(BB.fracdev_r, 100, range=(0,1))
    plt.xlabel('fracdev')
    plt.title('SDSS galaxies $r_e$ in $[0.5,1.5]$"')
    ps.savefig()

    ff = np.linspace(0, 1, 100)
    KK = []
    for flo,fhi in zip(ff, ff[1:]):
        K = np.flatnonzero((BB.fracdev_r >= flo) * (BB.fracdev_r < fhi))
        print len(K), 'in', flo,fhi
        KK.append(np.random.permutation(K)[:100])
    K = np.hstack(KK)
    
    plt.clf()
    plt.scatter(c2[K], c1[K], c=BB.fracdev_r[K], edgecolors='none', alpha=0.5)
    plt.xlabel('g - i')
    plt.ylabel('i - W1')
    plt.colorbar()
    plt.axis(ax)
    plt.title('SDSS galaxies $r_e$ in $[0.5,1.5]$", (color: fracdev\_r)')
    ps.savefig()

    K = np.sort(K)
    plt.clf()
    plt.scatter(c2[K], c1[K], c=BB.fracdev_r[K], edgecolors='none', alpha=0.5)
    plt.xlabel('g - i')
    plt.ylabel('i - W1')
    plt.colorbar()
    plt.axis(ax)
    plt.title('SDSS galaxies $r_e$ in $[0.5,1.5]$", (color: fracdev\_r; sorted)')
    ps.savefig()


    for K,ss in [
        (np.flatnonzero(BB.fracdev_r > 0.9), 'fracdev $> 0.9$'),
        (np.flatnonzero(BB.fracdev_r < 0.1), 'fracdev $< 0.1$'),
        #(np.flatnonzero(BB.fracdev_r > 0.95), 'fracdev $> 0.95$'),
        #(np.flatnonzero(BB.fracdev_r < 0.05), 'fracdev $< 0.05$'),
        ]:
        loghist(c2[K], c1[K], range=((0,3.2),(2,5)), bins=100, hot=False,
                imshowargs=dict(cmap='gray'))
        plt.xlabel('g - i')
        plt.ylabel('i - W1')
        plt.title('SDSS galaxies $r_e$ in $[0.5,1.5]$", %s' % ss)
        ps.savefig()
        # Match to spectroscopic galaxies...
        I,J,d = match_radec(M.ra, M.dec, FF.ra[K], FF.dec[K], 1./3600.)
        print ss, 'has', len(I), 'spectroscopic matches'
        ax = plt.axis()
        plt.scatter(c2[K][J], c1[K][J], c=M.z[I], edgecolors='none', alpha=0.5)
        plt.axis(ax)
        plt.colorbar()
        ps.savefig()

    I,J,d = match_radec(M.ra, M.dec, F.ra, F.dec, 1./3600.)
    print len(I), 'matces'

    MM = M[I]
    FF = F[J]

    K = np.isfinite(FF.w1_mag)
    print sum(K), 'with good W1'
    MM.cut(K)
    FF.cut(K)

    K = np.flatnonzero((FF.w1_mag_err < 1) * (FF.dg_mag < 1) * (FF.di_mag < 1))
    FF.cut(K)
    MM.cut(K)
    c1 = FF.i_mag - FF.w1_mag
    c2 = FF.g_mag - FF.i_mag

    ax = [1, 3.2, 2.5, 5.0]

    plt.clf()
    plt.scatter(FF.g_mag - FF.r_mag, FF.i_mag - FF.z_mag, c=MM.z, edgecolors='none', alpha=0.5)
    plt.xlabel('g - r')
    plt.ylabel('i - z')
    plt.colorbar()
    plt.axis([0.7, 2.5, 0.1, 0.7])
    plt.title('SDSS galaxies')
    ps.savefig()

    plt.clf()
    plt.scatter(c2, c1, c=MM.z, edgecolors='none', alpha=0.5)
    plt.xlabel('g - i')
    plt.ylabel('i - W1')
    plt.colorbar()
    plt.axis(ax)
    plt.title('SDSS galaxies in SDSS/WISE colors')
    ps.savefig()


    #I = np.flatnonzero((MM.z > 0.3) * (MM.z < 0.4))
    #J = np.flatnonzero((MM.z > 0.4) * (MM.z < 0.5))
    I = np.flatnonzero((MM.z > 0.325) * (MM.z < 0.375))
    J = np.flatnonzero((MM.z > 0.425) * (MM.z < 0.475))
    plt.clf()
    #plt.errorbar(c2[I], c1[I], fmt='b.', xerr=np.hypot(FF.dg_mag[I], FF.di_mag[I]),
    #             yerr=np.hypot(FF.di_mag[I], FF.w1_mag_err[I]), alpha=0.5)
    #p1 = plt.plot(c2[I], c1[I], 'bo', mec='none', alpha=0.5)
    #p2 = plt.plot(c2[J], c1[J], 'ro', mec='none', alpha=0.5)

    IJ = np.append(I,J)
    IJ = np.random.permutation(IJ)
    zz = MM.z[IJ]
    zz = np.array(['b' if z < 0.4 else 'r' for z in zz])

    plt.scatter(c2[IJ], c1[IJ], c=zz, edgecolors='none', alpha=0.5)

    p1 = plt.plot([],[], 'bo', mec='none', alpha=0.5)
    p2 = plt.plot([],[], 'ro', mec='none', alpha=0.5)
    plt.legend((p1[0], p2[0]), ('$0.325 < z < 0.375$', '$0.425 < z < 0.475$'), loc='upper left')

    plt.xlabel('g - i')
    plt.ylabel('i - W1')
    plt.title('SDSS galaxies')
    plt.axis(ax)
    ps.savefig()


    for zlo,zhi in [(0.25, 0.3), (0.3, 0.35), (0.35,0.4), (0.4, 0.45), (0.45,0.5)]:
        I = np.flatnonzero((MM.z > zlo) * (MM.z < zhi))
        plt.clf()
        plt.errorbar(c2[I], c1[I], fmt='b.', xerr=np.hypot(FF.dg_mag[I], FF.di_mag[I]),
                     yerr=np.hypot(FF.di_mag[I], FF.w1_mag_err[I]), alpha=0.5)
        plt.xlabel('g - i')
        plt.ylabel('i - W1')
        plt.title('Redshift $%.2f < z < %.2f$' % (zlo,zhi))
        plt.axis(ax)
        ps.savefig()



    for band,(ylo,yhi) in zip([1,2,3,4], [(3,8.5),(3,8.5),(5,11),(8,14)]):

        w = FF.get('w%i_mag' % band)
        werr = FF.get('w%i_mag_err' % band)
        I = np.flatnonzero(werr < 1)

        MI = MM[I]
        FI = FF[I]
        w = w[I]
        werr = werr[I]

        plt.clf()
        plt.errorbar(MI.z, FI.g_mag - w, yerr=np.hypot(werr, FI.dg_mag), fmt='b.', color='c')
        plt.errorbar(MI.z, FI.g_mag - w, yerr=werr, fmt='b.', color='b')

        subs = np.unique(MI.subclazz)
        ll,lt = [],[]
        for sub,cc in zip(subs, 'ygrkcmb'):
            if len(sub.strip()) == 0:
                continue
            J = (MI.subclazz == sub)
            p = plt.errorbar(MI.z[J], FI.g_mag[J] - w[J], yerr=werr[J], fmt='o'+cc, color=cc)
            ll.append(p[0])
            lt.append(sub.strip())
        plt.legend(ll,lt, loc='lower right', prop=dict(size='x-small'))
        plt.xlabel('Redshift z')
        plt.ylabel('g - W%i' % band)
        plt.title('SDSS galaxies in WISE: W%i' % band)
        plt.ylim(ylo,yhi)
        ps.savefig()

        # plt.clf()
        # plt.errorbar(MM.z[I], (FF.r_mag - w)[I], yerr=werr[I], fmt='b.', color='b')
        # plt.xlabel('Redshift z')
        # plt.ylabel('r - W%i' % band)
        # ps.savefig()
        # 
        # plt.clf()
        # plt.errorbar(MM.z[I], (FF.i_mag - w)[I], yerr=werr[I], fmt='b.', color='b')
        # plt.xlabel('Redshift z')
        # plt.ylabel('i - W%i' % band)
        # ps.savefig()


def mstars(T, tdir, pdir):
    ps = PlotSequence('stars')
    M = fits_table('mstars.fits')
    print len(M), 'M stars'

    print 'RA', M.ra.min(), M.ra.max()
    print 'Dec', M.dec.min(), M.dec.max()

    M.giant = np.array(['III' in s for s in M.subclazz])
    M.dwarf = np.array(['V' in s for s in M.subclazz])

    allwise = True
    FF = []
    for tile in T.coadd_id:
        W,F = read_sources(allwise, tile, tdir, pdir)
        FF.append(F)
    F = merge_tables(FF)
    print len(F), 'forced-phot'

    F.cut(np.isfinite(F.g_mag) * np.isfinite(F.i_mag))

    I,J,d = match_radec(M.ra, M.dec, F.ra, F.dec, 1./3600.)
    print len(I), 'matches'

    MM = M[I]
    FF = F[J]

    K = np.isfinite(FF.w1_mag)
    print sum(K), 'with good W1'
    MK = MM[K]
    FK = FF[K]

    for band,(ylo,yhi) in zip([1,2,3,4], [(2,9), (1,8), (None,None),(None,None)]):
        plt.clf()
        w = FK.get('w%i_mag' % band)
        werr = FK.get('w%i_mag_err' % band)
        p1 = plt.plot((FK.g_mag - FK.i_mag), (FK.i_mag - w), 'b.')
        ax = plt.axis()
        p1 = plt.errorbar((FK.g_mag - FK.i_mag), (FK.i_mag - w),
                          yerr=werr, fmt='b.', color='b')
        #p2 = plt.plot((FK.g_mag - FK.i_mag)[MK.giant], (FK.i_mag - w)[MK.giant], 'ro')
        #p3 = plt.plot((FK.g_mag - FK.i_mag)[MK.dwarf], (FK.i_mag - w)[MK.dwarf], 'go')
        p2 = plt.errorbar((FK.g_mag - FK.i_mag)[MK.giant], (FK.i_mag - w)[MK.giant],
                          yerr=werr[MK.giant], fmt='ro', color='r')
        p3 = plt.errorbar((FK.g_mag - FK.i_mag)[MK.dwarf], (FK.i_mag - w)[MK.dwarf],
                          yerr=werr[MK.dwarf], fmt='go', color='g')
        plt.axis(ax)
        if ylo is not None:
            plt.ylim(ylo,yhi)
        plt.xlim(1,5)
        plt.xlabel('g - i')
        plt.ylabel('i - W%i' % band)
        plt.legend([p2[0],p3[0],p1[0]], ('Giant', 'Dwarf', 'Unknown'))
        plt.title('M stars')
        ps.savefig()


def psfcutoff():
    ps = PlotSequence('cut')

    tile = '1800p000'
    
    plt.clf()
    #plt.loglog(T1.w1_nanomaggies, T2.w1_nanomaggies, 'b.', alpha=0.1)
    lo,hi = 0.1,10.
    loghist(np.log10(T1.w1_nanomaggies),
            np.clip(T2.w1_nanomaggies / T1.w1_nanomaggies, lo, hi),
            200, range=((-3,8),(lo,hi)))
    ps.savefig()
    
    T1 = fits_table('vanilla/phot-%s.fits' % tile)
    T2 = fits_table('psfcutoff-c/phot-%s.fits' % tile)
    print len(T1), len(T2), 'sources'

    plt.clf()
    #plt.loglog(T1.w1_nanomaggies, T2.w1_nanomaggies, 'b.', alpha=0.1)
    lo,hi = 0.8,1.2
    loghist(np.log10(T1.w1_nanomaggies),
            np.clip(T2.w1_nanomaggies / T1.w1_nanomaggies, lo, hi),
            200, range=((-3,8),(lo,hi)))
    ax = plt.axis()
    mag = T1.w1_mag
    I = np.flatnonzero(np.isfinite(mag))
    mags = np.arange(np.ceil(mag[I].min()), np.floor(mag[I].max())+1)
    fluxes = []
    stds = []
    for mlo,mhi in zip(mags, mags[1:]):
        I = np.flatnonzero((mag > mlo) * (mag < mhi))
        fluxes.append(np.mean(T1.w1_nanomaggies[I]))
        stds.append(1./np.sqrt(np.mean(T1.w1_nanomaggies_ivar[I])))
    fluxes = np.array(fluxes)
    plt.errorbar(np.log10(fluxes), np.ones_like(fluxes), yerr=stds/fluxes)
    plt.axis(ax)
    
    ps.savefig()

    I = np.logical_or(T1.pointsource, T1.treated_as_pointsource)
    for cut in [I, np.logical_not(I)]:
        plt.clf()
        loghist(np.log10(T1.w1_nanomaggies[cut]),
                np.clip(T2.w1_nanomaggies[cut] / T1.w1_nanomaggies[cut], lo, hi),
                200, range=((-3,8),(lo,hi)))
        ps.savefig()



    T1 = fits_table('psfcutoff-c/phot-%s.fits' % tile)
    T2 = fits_table('psfcutoff-5c/phot-%s.fits' % tile)
    print len(T1), len(T2), 'sources'

    plt.clf()
    plt.loglog(T1.w1_nanomaggies, T2.w1_nanomaggies, 'b.', alpha=0.1)
    ps.savefig()
    
    plt.clf()
    #plt.loglog(T1.w1_nanomaggies, T2.w1_nanomaggies, 'b.', alpha=0.1)
    lo,hi = 0.8,1.2
    loghist(np.log10(T1.w1_nanomaggies),
            np.clip(T2.w1_nanomaggies / T1.w1_nanomaggies, lo, hi),
            200, range=((-3,8),(lo,hi)))
    ax = plt.axis()
    plt.errorbar(np.log10(fluxes), np.ones_like(fluxes), yerr=stds/fluxes)
    plt.axis(ax)
    ps.savefig()


    T1 = fits_table('psfcutoff-5c/phot-%s.fits' % tile)
    T2 = fits_table('psfcutoff-10c/phot-%s.fits' % tile)
    print len(T1), len(T2), 'sources'

    plt.clf()
    lo,hi = 0.8,1.2
    loghist(np.log10(T1.w1_nanomaggies),
            np.clip(T2.w1_nanomaggies / T1.w1_nanomaggies, lo, hi),
            200, range=((-3,8),(lo,hi)))
    ax = plt.axis()
    plt.errorbar(np.log10(fluxes), np.ones_like(fluxes), yerr=stds/fluxes)
    plt.axis(ax)
    ps.savefig()

    T1 = fits_table('psfcutoff-10c/phot-%s.fits' % tile)
    T2 = fits_table('psfcutoff-10d/phot-%s.fits' % tile)
    print len(T1), len(T2), 'sources'

    plt.clf()
    lo,hi = 0.8,1.2
    loghist(np.log10(T1.w1_nanomaggies),
            np.clip(T2.w1_nanomaggies / T1.w1_nanomaggies, lo, hi),
            200, range=((-3,8),(lo,hi)))
    ax = plt.axis()
    plt.errorbar(np.log10(fluxes), np.ones_like(fluxes), yerr=stds/fluxes)
    plt.axis(ax)
    ps.savefig()
    

def coadd_demo():
    from unwise_coadd import zeropointToScale
    from astrometry.util.resample import *
    
    T = fits_table('unwise-coadds/180/1800p000/unwise-1800p000-w1-frames.fits')
    print len(T), 'frames'
    T.cut(T.use == 1)
    print len(T), 'used'

    cowcs = Tan('unwise-coadds/180/1800p000/unwise-1800p000-w1-img-m.fits', 0)
    
    ps = PlotSequence('codemo')

    plt.figure(figsize=(4,4))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    
    mx = (T.coextent[:,0] + T.coextent[:,1]) / 2
    T.cut((mx > 500) * (mx < 1500))
    my = (T.coextent[:,2] + T.coextent[:,3]) / 2
    T.cut((my > 500) * (my < 1500))
    print len(T), 'central'

    #cowcs = cowcs.get_subimage(1000, 1000, 100, 100) #25, 25)#50, 50)
    cowcs = cowcs.get_subimage(1033, 1000, 33, 33)
    coH,coW = cowcs.get_height(), cowcs.get_width()
    coadd = np.zeros((coH, coW), np.float32)
    con = np.zeros((coH, coW), np.float32)
    
    for t in T:
        fn = ('wise-frames/%s/%s/%03i/%s%03i-w1-int-1b.fits' %
              (t.scan_id[-2:], t.scan_id, t.frame_num, t.scan_id,
               t.frame_num))

        try:
            wcs = Sip(fn, 0)
        except RuntimeError:
            import traceback
            traceback.print_exc()
            print 'carrying on'
            continue
        try:
            Yo,Xo,Yi,Xi,nil = resample_with_wcs(cowcs, wcs, [], 3)
        except:
            continue
        print len(Yo), 'pixels overlap'
        if len(Yo) < (coH*coW/2):
            continue
        
        I = fitsio.read(fn)

        #print 'Image median', np.median(I), 'vs sky1', t.sky1
        #I -= t.sky1

        I -= np.median(I)
        
        zp = t.zeropoint
        zpscale = 1. / zeropointToScale(zp)
        I *= zpscale
        iv = t.weight
        sig1 = np.sqrt(1./iv)

        maskfn = fn.replace('-int-', '-msk-')
        maskfn = maskfn + '.gz'
        mask = fitsio.read(maskfn)
        badbits = [0,1,2,3,4,5,6,7, 9, 
                   10,11,12,13,14,15,16,17,18,
                   21,26,27,28]
        if t.phase == 3:
            # 3-band cryo phase:
            ## 19 pixel is "hard-saturated"
            ## 23 for W3 only: static-split droop residual present
            badbits.append(19)
        maskbits = sum([1<<bit for bit in badbits])
        goodmask = ((mask & maskbits) == 0)
        goodmask[np.logical_not(np.isfinite(I))] = False

        I[np.logical_not(goodmask)] = 0.

        # diffraction-limited PSF
        convimg = gaussian_filter(I, 0.88)

        rimg = np.zeros((coH, coW), np.float32)
        rimg[Yo,Xo] += convimg[Yi,Xi]

        coadd[Yo,Xo] += convimg[Yi,Xi] * goodmask[Yi,Xi]
        con[Yo,Xo] += goodmask[Yi,Xi]

        r,d = cowcs.pixelxy2radec(np.array([0,0,coW,coW]),
                                  np.array([0,coH,coH,0]))
        print 'r,d', r,d
        ok,x,y = wcs.radec2pixelxy(r,d)
        print 'x,y', x,y
        x0 = x.min()
        x1 = x.max()
        y0 = y.min()
        y1 = y.max()

        keep = np.zeros(I.shape, bool)
        keep[Yi,Xi] = True
        keep = binary_dilation(keep)
        I[np.logical_not(keep)] = 0.
        convimg[np.logical_not(keep)] = 0.
        
        plt.clf()
        plt.imshow(I, interpolation='nearest', origin='lower',
                   vmin=-3*sig1, vmax=10*sig1)
        ax = [x0,x1,y0,y1]
        plt.axis(ax)
        ps.savefig()

        plt.clf()
        plt.imshow(convimg, interpolation='nearest', origin='lower',
                   vmin=-3*sig1, vmax=10*sig1)
        plt.axis(ax)
        ps.savefig()
        
        plt.clf()
        plt.imshow(rimg, interpolation='nearest', origin='lower',
                   vmin=-3*sig1, vmax=10*sig1)
        ps.savefig()

        plt.clf()
        plt.imshow(coadd / np.maximum(con,1),
                   interpolation='nearest', origin='lower',
                   vmin=-3*sig1, vmax=10*sig1)
        ps.savefig()

        
def main():

    ps = PlotSequence('l1b')

    tile = '1384p106'

    T1 = fits_table('vanilla2/phot-%s.fits' % tile)
    T2 = fits_table('l1b-sky/phot-%s.fits' % tile)
    print len(T1), len(T2), 'sources'

    mn = 1e-2
    plt.clf()
    plt.loglog(np.maximum(mn, T1.w1_nanomaggies),
               np.maximum(mn, T2.w1_nanomaggies), 'b.', alpha=0.1)
    plt.xlabel('Coadd flux (nmgy)')
    plt.ylabel('L1b flux (nmgy)')
    ps.savefig()

    plt.clf()
    plt.loglog(T1.w1_nanomaggies,
               T2.w1_nanomaggies / T1.w1_nanomaggies,
               'b.', alpha=0.1)
    plt.xlabel('Coadd flux (nmgy)')
    plt.ylabel('L1b flux / Coadd flux')
    plt.ylim(0.1, 10.)
    plt.xlim(mn, T1.w1_nanomaggies.max())
    ps.savefig()

    chi = ((T2.w1_nanomaggies - T1.w1_nanomaggies) /
           np.sqrt(1./T1.w1_nanomaggies_ivar + 1./T2.w1_nanomaggies_ivar))
    
    plt.clf()
    plt.semilogx(T1.w1_nanomaggies, np.clip(chi,-10,10),
                 'b.', alpha=0.1)
    plt.xlabel('Coadd flux (nmgy)')
    plt.ylabel('(L1b flux - Coadd flux) / Error')
    plt.ylim(-10, 10)
    ps.savefig()

    plt.clf()
    n,b,p = plt.hist(np.clip(chi, -10,10), range=(-10,10), bins=100, histtype='step', color='b')
    xx = np.linspace(-10,10,500)
    db = b[1]-b[0]
    plt.plot(xx, db * 1./np.sqrt(2.*np.pi) * np.exp(-0.5 * xx**2)*len(T1))
    plt.xlabel('(L1b flux - Coadd flux) / Error')
    ps.savefig()

    mags = np.arange(25, 11, -1)
    for mhi,mlo in zip(mags, mags[1:]):
        I = np.flatnonzero((T1.w1_mag >= mlo) * (T1.w1_mag < mhi))
        plt.clf()
        n,b,p = plt.hist(np.clip(chi[I], -10,10),
                         range=(-10,10), bins=100, histtype='step', color='b')
        xx = np.linspace(-10,10,500)
        db = b[1]-b[0]
        plt.plot(xx, db * 1./np.sqrt(2.*np.pi) * np.exp(-0.5 * xx**2)*sum(n))
        plt.xlabel('(L1b flux - Coadd flux) / Error')
        plt.title('W1 mag %g to %g' % (mlo,mhi))
        ps.savefig()
        

    
    sys.exit(0)

    
    dataset = 'sdss-dr10d'
    pdir = '%s-phot' % dataset
    tdir = '%s-phot-temp' % dataset
    T = fits_table('%s-atlas.fits' % dataset)
    T.cut((T.ra > 172.5) * (T.ra < 187.5) * (T.dec > 40) * (T.dec < 50))
    mstars(T, tdir, pdir)

    #sequels_qso_forcedphot()
    #psfcutoff()

    coadd_demo()


    # python -c "from astrometry.util.fits import *; T=fits_table('allsky-atlas.fits'); T.cut(T.dec > -27); T.cut(np.lexsort((T.ra, T.dec))); T.writeto('north-atlas.fits')"

    #dataset = 'sdss-dr10d'
    dataset = 'north'
    pdir = '%s-phot' % dataset
    tdir = '%s-phot-temp' % dataset
    T = fits_table('%s-atlas.fits' % dataset)
    
    # ~100 sq deg
    T.index = np.arange(len(T))
    #T.cut((T.ra > 172.5) * (T.ra < 187.5) * (T.dec > 40) * (T.dec < 50))
    #T.cut((T.ra > 175) * (T.ra < 185) * (T.dec > 40) * (T.dec < 45))
    if False:
        T.cut((T.ra > 172.5) * (T.ra < 187.5) * (T.dec > 40) * (T.dec < 50))
        T = T[:5]
        ps = PlotSequence('lrg')
        ps.suffixes = ['png', 'pdf']
        lrg_comparison(tdir, pdir, T, ps)
        sys.exit(0)

    print len(T), 'tiles'
    #print ' '.join('%i'%i for i in T.index)
    #T = T[:10]
    print 'Cut to', len(T), 'tiles'
    unwdir = 'data/unwise/unwise-comp'

    #profiles(T, tdir, pdir, unwdir)
    #treated_as_pointsource(T, tdir, pdir)
    bright(T, tdir, pdir, unwdir)
    #density(T, tdir, pdir)
    #wdensity(T, tdir, pdir, '-north')
    sys.exit(0)

    if False:
        S = unpickle_from_file('data/wphot/density-maps-dr10d.pickle')
        sdss = S['maps'][0]
        wmask = (sdss > 0)
        wdensity(T, tdir, pdir, '-masked', mask=wmask)
    #density_pobjs(T, tdir, pdir, '%s-pobj' % dataset)

    maxes = { 'AllWISE': 70000,
              'All SDSS': 70000,
              'SDSS stars': 70000,
              'SDSS galaxies': 70000,
              'Treated as point sources': 70000,
              #'Treated as galaxies':
              }

    tt = [0, 20000, 40000, 60000]
    ttg = [0, 2000, 4000, 6000]
    ticks = { 'AllWISE': tt,
              'All SDSS': tt,
              'SDSS stars': tt,
              'SDSS galaxies': tt,
              'Treated as point sources': tt,
              'Treated as galaxies': ttg,
              }
    if True:
        ps = PlotSequence('wmap')
        ps.suffixes = ['png', 'pdf']
        #X = unpickle_from_file('data/wphot/density-sdss-dr10.pickle')
        #print 'Maps:', X.keys()
        #wmaps = []
        #wmapnames = []
        #names = X['mapnames']
        # for nm,themap in zip(names, X['maps']):
        #     if not nm.startswith('AllWISE'):
        #         continue
        #     wmaps.append(themap)
        #     wmapnames.append(nm)
        #X['maps'] = wmaps
        #X['mapnames'] = wmapnames

        #X = unpickle_from_file('data/wphot/density-wmaps-dr10d.pickle')

        X = unpickle_from_file('data/wphot/density-wmaps-north.pickle')
        Nmaps = 1
        X['maps'] = X['wmaps'][:Nmaps]
        X['mapnames'] = X['wmapnames'][:Nmaps]

        density_map_plots(ps, maxes=maxes, ticks=ticks, **X)

        S = unpickle_from_file('data/wphot/density-maps-dr10d.pickle')
        forced = (S['maps'][0] > 0)
        f2 = binary_dilation(forced, structure=np.ones((3,3), bool),
                             iterations=5)
        outline = f2 - forced
        X.update(outline=outline)
        
        # wmap = X['maps'][0]
        # forced = S['maps'][0]
        # mapralo = X['mapralo']
        # maprahi = X['maprahi']
        # mapdra  = X['mapdra']
        # mapdeclo = X['mapdeclo']
        # mapdechi = X['mapdechi']
        # mapddec  = X['mapddec']
        # rabins  = np.arange(mapralo , maprahi +0.01*mapdra,  mapdra )
        # decbins = np.arange(mapdeclo, mapdechi+0.01*mapddec, mapddec)
        # dec = (decbins[:-1] + decbins[1:]) / 2.
        # cosdec = np.cos(np.deg2rad(dec))
        # print 'AllWISE Map sum:', wmap.sum()
        # print 'Map non-zero area:', ((wmap > 0) * cosdec[:,np.newaxis] * mapdra * mapddec).sum(), 'sq deg'
        # print 'SDSS Map sum:', forced.sum()
        # print 'Map non-zero area:', ((forced > 0) * cosdec[:,np.newaxis] * mapdra * mapddec).sum(), 'sq deg'
        # wmap[forced == 0] = 0
        # print 'Masked AllWISE Map sum:', wmap.sum()
        # print 'Map non-zero area:', ((wmap > 0) * cosdec[:,np.newaxis] * mapdra * mapddec).sum(), 'sq deg'
        # sys.exit(0)

        #X = unpickle_from_file('data/wphot/density-pomaps-dr10d.pickle')
        #X['maps'] = X['pomaps']
        #X['mapnames'] = X['pomapnames']

        density_map_plots(ps, maxes=maxes, ticks=ticks, **X)

    if False:
        ps = PlotSequence('density')
        ps.suffixes = ['png', 'pdf']
        #Y = unpickle_from_file('data/wphot/density-sdss-dr10.pickle')
        #Y = unpickle_from_file('data/wphot/density-whists-dr10d.pickle')
        Y = unpickle_from_file('data/wphot/density-whists-masked.pickle')
        X = {}
        X['whists'] = Y['whists']
        X.update(unpickle_from_file('data/wphot/density-hists-dr10d.pickle'))
        print 'Hists:', X.keys()
        density_hist_plots(ps, **X)

    if False:
        ps = PlotSequence('map')
        ps.suffixes = ['png', 'pdf']
        X = unpickle_from_file('data/wphot/density-maps-dr10d.pickle')
        print 'Maps:', X.keys()
        Nmaps = 5
        X['maps'] = X['maps'][:Nmaps]
        X['mapnames'] = X['mapnames'][:Nmaps]
        density_map_plots(ps, maxes=maxes, ticks=ticks, **X)

    sys.exit(0)
    
    ps = PlotSequence('comp')
    ps.suffixes = ['png', 'pdf']
    

    if True:
        Wallsky = None
        #for allwise in [False, True]:
        for allwise in [True]:
            Wallsky = compare_vs_wise_cat(allwise, T, tdir, pdir, ps, Wallsky)

    #photoz(T, tdir, pdir)
    #lrg_comparison(tdir, pdir, T, ps)
    sys.exit(0)

    #ps.skipto(36)

    # Match
    if False:
        allwise = True
        WW,SS = [],[]
        for tile in T.coadd_id:
            W,S = read_sources(allwise, tile, tdir, pdir)
            fn = os.path.join('sdss3-phot-temp', 'photoobjs-%s.fits' % tile)
            P = fits_table(fn, columns=['psfflux','ra','dec'])
            I,J,d = match_radec(P.ra, P.dec, S.ra, S.dec, 0.1/3600.)
            print len(S), 'phot objs'
            print len(P), 'photoObjs'
            print len(I), 'matches'
            S.psfflux = np.zeros((len(S),5))
            S.psfflux[J] = P.psfflux[I]
            # S.psfflux = P.psfflux
            I,J,d,counts = match_radec(W.ra, W.dec, S.ra, S.dec, 4./3600.,
                                       nearest=True, count=True)
            K = np.flatnonzero(counts == 1)
            I = I[K]
            J = J[K]
            
            print len(I), 'matches'
            W.cut(I)
            S.cut(J)
            WW.append(W)
            SS.append(S)
        W = merge_tables(WW)
        F = merge_tables(SS)
        stars = np.flatnonzero(np.logical_or(F.objc_type == 6, F.treated_as_pointsource == 1))
        print len(stars), 'stars'
        F.cut(stars)
        W.cut(stars)
    
        F.r_psf = -2.5 * (np.log10(F.psfflux[:,2]) - 9)
        F.i_psf = -2.5 * (np.log10(F.psfflux[:,3]) - 9)
        F.g_psf = -2.5 * (np.log10(F.psfflux[:,1]) - 9)
    
        dnm = 1./np.sqrt(F.psfflux_ivar[:,2])
        nm = F.psfflux[:,2]
        F.dr_psf = np.abs((-2.5 / np.log(10.)) * dnm / nm)
        F.cut(F.dr_psf < 0.5)
        print 'Cut to', len(F), 'with r-mag err < 0.5'
    
        plt.clf()
        loghist(F.g_psf - F.r_psf, F.r_psf - F.i_psf, 100,
                range=((-5,5),(-5,5)))
        plt.xlabel('g - r')
        plt.ylabel('r - i')
        ps.savefig()
    
        plt.clf()
        plothist(F.g_psf - F.r_psf, F.r_psf - F.i_psf, 100,
                range=((-5,5),(-5,5)))
        plt.xlabel('g - r')
        plt.ylabel('r - i')
        ps.savefig()
    
        mags = np.arange(14, 19+1)
        for mlo,mhi in [(13,25)] + zip(mags, mags[1:]):
            I = np.flatnonzero((F.w1_mag > mlo) * (F.w1_mag < mhi))
            plt.clf()
            #loghist(F.r_psf[I] - F.w1_mag[I], F.w1_mag[I] - F.w2_mag[I], 100,
            #        range=((-5,10),(-4,3)))
            #plt.xlabel('r - W1')
            #plt.ylabel('W1 - W2')
            loghist(F.r_psf[I] - F.i_psf[I], F.r_psf[I] - F.w1_mag[I], 100,
                    range=((-5,5),(-10,10)))
            plt.ylabel('r - W1')
            plt.xlabel('r - i')
            plt.title('Matched; W1 in [%g,%g]' % (mlo,mhi))
            ps.savefig()
    #ps.skipto(44)


    cfn = 'F.fits'
    if os.path.exists(cfn):
        F = fits_table(cfn)
    else:
        # Find unmatched SDSS sources
        allwise = True
        FF = []
        # all WISE sources
        WW = []
        # unmatched
        UW = []
        # all SDSS sources
        SS = []
        for tile in T.coadd_id:
            W,S = read_sources(allwise, tile, tdir, pdir)
            WW.append(W)
            SS.append(S)
            fn = os.path.join('sdss3-phot-temp', 'photoobjs-%s.fits' % tile)
            # match to photoobj to get psfflux column
            P = fits_table(fn, columns=['psfflux', 'ra','dec'])
            I,J,d = match_radec(P.ra, P.dec, S.ra, S.dec, 0.1/3600.)
            print len(S), 'phot objs'
            print len(P), 'photoObjs'
            print len(I), 'matches'
            S.psfflux = np.zeros((len(S),5))
            S.psfflux[J] = P.psfflux[I]
    
            S.tile = np.array([tile] * len(S))
            I,J,d = match_radec(W.ra, W.dec, S.ra, S.dec, 4./3600.)
            # Keep only SDSS objects with no match
            print len(I), 'matches'
            K = np.ones(len(S), bool)
            K[J] = False
            print 'Keeping', sum(K), 'unmatched'
            FF.append(S[K])
            # Keep WISE objects with no match also
            K = np.ones(len(W), bool)
            K[I] = False
            UW.append(W[K])
        F = merge_tables(FF)
        W = merge_tables(WW)
        S = merge_tables(SS)
        UW = merge_tables(UW)
        UW.flux = 10.**((UW.w1mpro - 22.5) / -2.5)
        
        # Remove self-matches to avoid double-counting margin objects
        I,J,d = match_radec(F.ra, F.dec, F.ra, F.dec, 4./3600., notself=True)
        print len(I), 'self-matches'
        K = np.ones(len(F), bool)
        K[I] = False
        K[J] = False
        F.cut(K)
        print 'Keeping', len(F), 'unmatched'
    
        gals = np.flatnonzero((F.objc_type == 3) * (F.treated_as_pointsource == 0))
        print len(gals), 'galaxies'
        stars = np.flatnonzero(np.logical_or(F.objc_type == 6, F.treated_as_pointsource == 1))
        print len(stars), 'stars'
        #F.cut(stars)
    
        # Spatial distribution of unmatched
        # plt.clf()
        # plothist(F.ra, F.dec, 200)
        # ps.savefig()
    
        F.g_psf = -2.5 * (np.log10(F.psfflux[:,1]) - 9)
        F.r_psf = -2.5 * (np.log10(F.psfflux[:,2]) - 9)
        F.i_psf = -2.5 * (np.log10(F.psfflux[:,3]) - 9)
    
        F.g_gal = -2.5 * (np.log10(F.modelflux[:,1]) - 9)
        F.r_gal = -2.5 * (np.log10(F.modelflux[:,2]) - 9)
        F.i_gal = -2.5 * (np.log10(F.modelflux[:,3]) - 9)
    
        F.g_mag = np.zeros(len(F))
        F.r_mag = np.zeros(len(F))
        F.i_mag = np.zeros(len(F))
        F.g_mag[stars] = F.g_psf[stars]
        F.g_mag[gals ] = F.g_gal[gals ]
        F.r_mag[stars] = F.r_psf[stars]
        F.r_mag[gals ] = F.r_gal[gals ]
        F.i_mag[stars] = F.i_psf[stars]
        F.i_mag[gals ] = F.i_gal[gals ]
    
        # Mag distributions of unmatched
        # plt.clf()
        # plt.hist(F.w1_mag[np.isfinite(F.w1_mag)], 100, histtype='step', color='b', range=(13,25))
        # plt.hist(F.w2_mag[np.isfinite(F.w2_mag)], 100, histtype='step', color='g', range=(13,25))
        # plt.hist(F.r_psf[np.isfinite(F.r_psf)], 100, histtype='step', color='r', range=(13,25))
        # plt.hist(F.i_psf[np.isfinite(F.i_psf)], 100, histtype='step', color='m', range=(13,25))
        # plt.hist(F.g_psf[np.isfinite(F.g_psf)], 100, histtype='step', color='c', range=(13,25))
        # plt.xlabel('Mag')
        # ps.savefig()
    
        F.w1_snr = F.w1_nanomaggies * np.sqrt(F.w1_nanomaggies_ivar)

        F.writeto(cfn)

    if False:
        model_plots(F, S, UW, unwdir)

    # What are the objects with large S/N but no WISE catalog detection?
    # Edges of bright stars, and blends

    plt.clf()
    I = np.isfinite(F.w1_mag)
    loghist(F.w1_mag[I], F.w1_snr[I], 100, range=((16,24),(-1,10)))
    plt.axhline(5., color='b')
    plt.axhline(4., color='b')
    plt.axhline(0., color='b')
    plt.axvline(17., color='b')
    plt.axvline(17.5, color='b')
    plt.axvline(18., color='b')
    plt.axvline(18.5, color='b')
    plt.xlabel('W1 mag')
    plt.ylabel('W1 S/N')
    ps.savefig()
    
    dnm = 1./np.sqrt(F.psfflux_ivar[:,2])
    nm = F.psfflux[:,2]
    F.dr_psf = np.abs((-2.5 / np.log(10.)) * dnm / nm)

    dnm = 1./np.sqrt(F.psfflux_ivar[:,3])
    nm = F.psfflux[:,3]
    F.di_psf = np.abs((-2.5 / np.log(10.)) * dnm / nm)

    F.cut(F.dr_psf < 0.5)
    print 'Cut to', len(F), 'with r-mag err < 0.5'
    F.cut(F.dr_psf < 0.25)
    print 'Cut to', len(F), 'with r-mag err < 0.25'
    F.cut(F.di_psf < 0.5)
    print 'Cut to', len(F), 'with i-mag err < 0.5'
    F.cut(F.di_psf < 0.25)
    print 'Cut to', len(F), 'with i-mag err < 0.25'

    # plt.clf()
    # loghist(F.g_psf - F.r_psf, F.r_psf - F.i_psf, 100,
    #         range=((-5,5),(-5,5)))
    # plt.xlabel('g - r')
    # plt.ylabel('r - i')
    # ps.savefig()

    plt.clf()
    plothist(F.g_psf - F.r_psf, F.r_psf - F.i_psf, 100,
             range=((-1,2),(-1,2)))
    #range=((-5,5),(-5,5)))
    plt.xlabel('g - r')
    plt.ylabel('r - i')
    ps.savefig()

    #mags = np.arange(14, 20+1)
    mags = [13, 16, 17, 17.5, 18., 18.5, 19., 20.]
    for mlo,mhi in [(13,25)] + zip(mags, mags[1:]):
        I = np.flatnonzero((F.w1_mag > mlo) * (F.w1_mag < mhi))
        plt.clf()
        loghist(F.r_mag[I] - F.i_mag[I], F.r_mag[I] - F.w1_mag[I], 100,
                range=((-1,3),(-3,8)))
                #range=((-5,5),(-10,10)))
        ax = plt.axis()
        plt.plot([0.75, 2], [3, 5.75], 'b-')
        plt.plot([0.5, 2], [3, 5.75], 'b-')
        plt.plot([0.5, 2], [3, 6], 'g-')
        plt.axis(ax)

        plt.ylabel('r - W1')
        plt.xlabel('r - i')
        plt.title('W1 in [%g,%g]' % (mlo,mhi))
        ps.savefig()

        plt.clf()

        ri = F.r_mag[I] - F.i_mag[I]
        rw = F.r_mag[I] - F.w1_mag[I]
        
        ri0 = 1.
        rw0 = 4.
        dri = ri - ri0
        drw = rw - rw0
        angle = -np.arctan(2.)
        uu = np.cos(angle) * dri - np.sin(angle) * drw
        vv = np.sin(angle) * dri + np.cos(angle) * drw

        loghist(uu, vv, 100)
        ps.savefig()

        I = np.flatnonzero((uu > -0.5) * (uu < 1.0))
        plt.clf()
        loghist(uu[I], vv[I], 100)
        ps.savefig()

        plt.clf()
        plt.hist(vv[I], 50, range=(-0.5, 1.5), histtype='step')
        ps.savefig()

    mags = np.arange(14, 20+1)
    for mlo,mhi in [(13,25)] + zip(mags, mags[1:]):
        I = np.flatnonzero((F.w1_mag > mlo) * (F.w1_mag < mhi))
        snr = F.w1_snr[I]
        print 'Mag range', mlo,mhi, 'SNR', snr.min(), np.median(snr), snr.max()
        plt.clf()
        loghist(F.g_mag[I] - F.i_mag[I], F.i_mag[I] - F.w1_mag[I], 100,
                range=((-1,5),(-4,7)))
        plt.ylabel('i - W1')
        plt.xlabel('g - i')
        plt.title('W1 in [%g,%g]' % (mlo,mhi))
        ps.savefig()

    # Quasar cut -- cut to "star" too
    mags = np.arange(14, 20+1)
    for mlo,mhi in [(13,25)] + zip(mags, mags[1:]):
        I = np.flatnonzero((F.w1_mag > mlo) * (F.w1_mag < mhi) *
                           np.logical_or(F.objc_type == 6,
                                         F.treated_as_pointsource == 1))
        plt.clf()
        loghist(F.g_mag[I] - F.i_mag[I], F.i_mag[I] - F.w1_mag[I], 100,
                range=((-1,5),(-4,7)))
        plt.ylabel('i - W1')
        plt.xlabel('g - i')
        plt.title('W1 in [%g,%g], point sources' % (mlo,mhi))
        ps.savefig()

    sns = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 100.]
    for snlo,snhi in zip(sns, sns[1:]):
        I = np.flatnonzero((F.w1_snr > snlo) * (F.w1_snr < snhi))
        mag = F.w1_mag[I]
        print 'S/N range', snlo,snhi, 'mag', mag.min(), np.median(mag), mag.max()
        plt.clf()
        loghist(F.r_mag[I] - F.i_mag[I], F.r_mag[I] - F.w1_mag[I], 100,
                range=((-1,3),(-3,8)))
        plt.ylabel('r - W1')
        plt.xlabel('r - i')
        plt.title('W1 SNR in [%g,%g]' % (snlo, snhi))
        ps.savefig()


        


def read_sources(allwise, tile, tdir, pdir, read_if_missing=True,
                 match_to_sdss=True):
    from unwise_coadd import get_coadd_tile_wcs, tile_to_radec
    r,d = tile_to_radec(tile)
    wcs = get_coadd_tile_wcs(r, d)
    r0,r1,d0,d1 = wcs.radec_bounds()
    #print 'RA,Dec bounds:', r0,r1,d0,d1

    if allwise:
        # sdss3+ phot used AllWISE catalogs.
        fn = os.path.join(tdir, 'wise-sources-%s.fits' % tile)
    else:
        fn = os.path.join(tdir, 'wise-allsky-sources-%s.fits' % tile)

    if os.path.exists(fn):
        W = fits_table(fn)
        if not 'w1sigmpro' in W.get_columns():
            cmd = 'mv %s .' % fn
            print 'File does not have w1sigmpro; running:', cmd
            os.system(cmd)
    if not os.path.exists(fn):
        if not read_if_missing:
            return None,None
        assert(0)
        # this API changed
        W = read_wise_sources(fn, r0,r1,d0,d1,
                              extracols=['w1sigmpro','w2sigmpro',
                                         'w3sigmpro','w4sigmpro',
                                         'w1gmag','w1gerr', 'w2gmag','w2gerr',
                                         'w3gmag','w3gerr', 'w4gmag','w4gerr',
                                         ],
                              allwise=allwise)
        W.writeto(fn)
    print len(W), 'from', fn

    W.w1proflux = 10.**((W.w1mpro - 22.5) / -2.5)
    W.w1prodflux = np.abs(W.w1sigmpro * W.w1mpro * (np.log(10.) / -2.5))
    W.w1prosnr = W.w1proflux / W.w1prodflux

    sfn = os.path.join(pdir, 'phot-%s.fits' % tile)
    if not os.path.exists(sfn):
        print 'No such file:', sfn
        return None,None
    print 'Reading', sfn
    S = fits_table(sfn)

    S.star = np.logical_or(S.objc_type == 6, S.treated_as_pointsource == 1)
    S.gal = np.logical_not(S.star)
    S.sdss_star = (S.objc_type == 6)
    S.sdss_gal  = (S.objc_type == 3)

    if not match_to_sdss:
        return W,S

    # , columns=['ra','dec','objc_type','treated_as_pointsource',
    #            'pointsource', 'coadd_id',
    #            'w1_mag', 'w2_mag', 
    #            'w3_mag', 'w4_mag',
    #            'w1_mag_err', 'w2_mag_err',
    #            'w3_mag_err', 'w4_mag_err',
    #            'w1_nanomaggies', 'w2_nanomaggies', 
    #            'w3_nanomaggies', 'w4_nanomaggies',
    #            'w1_nanomaggies_ivar', 'w2_nanomaggies_ivar', 
    #            'w3_nanomaggies_ivar', 'w4_nanomaggies_ivar',
    #            'run', 'camcol', 'field',
    #            'psfflux_ivar'])

    fn = os.path.join(tdir, 'photoobjs-%s.fits' % tile)
    print 'Reading', fn
    P = fits_table(fn, columns=['psfflux','ra','dec'])

    if len(P) == 0 or len(S) == 0:
        return W,None

    I,J,d = match_radec(P.ra, P.dec, S.ra, S.dec, 0.1/3600.)
    S.psfflux = np.zeros((len(S),5))
    S.psfflux[J] = P.psfflux[I]

    for iband,band in [(1,'g'),(2,'r'),(3,'i'),(4,'z')]:
        psfflux = S.psfflux[:,iband]
        galflux = S.modelflux[:,iband]
        psfmag = -2.5 * (np.log10(psfflux) - 9)
        galmag = -2.5 * (np.log10(galflux) - 9)
        S.set('%s_psf' % band, psfmag)
        S.set('%s_gal' % band, galmag)
        # combo mag
        mag = np.zeros(len(S))
        mag[S.star] = psfmag[S.star]
        mag[S.gal ] = galmag[S.gal ]
        S.set('%s_mag' % band, mag)

        dnm = 1./np.sqrt(S.psfflux_ivar[:,iband])
        nm = psfflux
        S.set('d%s_psf' % band, np.abs((-2.5 / np.log(10.)) * dnm / nm))
        dnm = 1./np.sqrt(S.modelflux_ivar[:,iband])
        nm = S.modelflux[:,iband]
        S.set('d%s_gal' % band, np.abs((-2.5 / np.log(10.)) * dnm / nm))
        dmag = np.zeros(len(S))
        dmag[S.star] = S.get('d%s_psf' % band)[S.star]
        dmag[S.gal ] = S.get('d%s_gal' % band)[S.gal ]
        S.set('d%s_mag' % band, dmag)

    S.w1_snr = S.w1_nanomaggies * np.sqrt(S.w1_nanomaggies_ivar)
    S.w2_snr = S.w2_nanomaggies * np.sqrt(S.w2_nanomaggies_ivar)
    S.w3_snr = S.w3_nanomaggies * np.sqrt(S.w3_nanomaggies_ivar)
    S.w4_snr = S.w4_nanomaggies * np.sqrt(S.w4_nanomaggies_ivar)

    print len(S), 'from', fn
    return W,S

def compare_vs_wise_cat(allwise, T, tdir, pdir, ps, Wallsky):
    WW = []
    SS = []

    for tile in T.coadd_id:
        W,S = read_sources(allwise, tile, tdir, pdir)
        I,J,d,counts = match_radec(W.ra, W.dec, S.ra, S.dec, 4./3600.,
                                   nearest=True, count=True)
        print 'Counts:', np.unique(counts)
        print len(I), 'I', len(counts), 'counts'
        K = np.flatnonzero(counts == 1)
        I = I[K]
        J = J[K]
        
        print len(I), 'matches'
        W.cut(I)
        S.cut(J)
        WW.append(W)
        SS.append(S)
    
    W = merge_tables(WW)
    F = merge_tables(SS)

    if allwise is False:
        Wallsky = W
    
    print 'Forced:'
    F.about()
    print 'WISE:'
    W.about()

    gals = np.flatnonzero((F.objc_type == 3) * (F.treated_as_pointsource == 0))
    print len(gals), 'galaxies'
    
    stars = np.flatnonzero(np.logical_or(F.objc_type == 6, F.treated_as_pointsource == 1))
    print len(stars), 'stars'
    
    # Hack
    # I = []
    # if len(I):
    #     tile = T.coadd_id[0]
    #     fn = os.path.join(unwdir, tile[:3], tile, 'unwise-%s-w%i-img-m.fits' % (tile, band))
    #     print 'reading', fn
    #     img = fitsio.read(fn)
    #     imh,imw = img.shape
    #     wcs = Tan(fn)
    #     for i in I[:10]:
    #         r,d = W.ra[i], W.dec[i]
    #         ok,x,y = wcs.radec2pixelxy(r, d)
    #         print 'x,y', x,y
    #         ix = int(np.round(x)) - 1
    #         iy = int(np.round(y)) - 1
    #         #S = 50
    #         S = 25
    #         if ix < S or iy < S or ix > imw-S or iy > imh-S:
    #             continue
    #         plt.clf()
    #         ylo,yhi = iy-S, iy+S
    #         xlo,xhi = ix-S, ix+S
    #         plt.imshow(img[ylo:yhi+1, xlo:xhi+1], interpolation='nearest', origin='lower',
    #                    cmap='gray', extent=[xlo-0.5,xhi+0.5,ylo-0.5,yhi+0.5],
    #                    vmin=-10, vmax=50)
    #         ax = plt.axis()
    #     
    #         J = np.flatnonzero((np.abs(W.ra - r) < 0.05) * (np.abs(W.dec - d) < 0.05))
    #         ok,x,y = wcs.radec2pixelxy(W.ra[J], W.dec[J])
    #         plt.plot(x-1, y-1, 'ro', mec='r', mfc='none', ms=20)
    #     
    #         J = np.flatnonzero((np.abs(F.ra - r) < 0.05) * (np.abs(F.dec - d) < 0.05))
    #         ok,x,y = wcs.radec2pixelxy(F.ra[J], F.dec[J])
    #         plt.plot(x-1, y-1, 'gx', ms=20)
    #     
    #         plt.axis(ax)
    #         ps.savefig()
    
    plt.figure(figsize=(4,4))

    # leave room for plt.title
    spa = dict(bottom=0.1, top=0.93, right=0.99)
    #spa = dict(bottom=0.1, top=0.98, right=0.99)
    
    for band in [1, 2, 3, 4]:
        wisempro = W.get('w%impro' % band)
        tracmag = F.get('w%i_mag' % band)

        lo, cathi = 10.2, 17.8
        if band == 3:
            lo,cathi = 8, 13
        elif band == 4:
            lo,cathi = 4.5, 9.5

        ha = dict(bins=200,
                  imshowargs=dict(cmap=antigray), hot=False,
                  docolorbar=False)

        # Compare 2MASS-aperture mags
        gmag = W.get('w%igmag' % band)
        has_gmag = np.flatnonzero(np.isfinite(gmag))

        # NOTE, we do a nasty little thing of using the local
        # variables in this loop in the code following the loop, so
        # ensure that the point sources are always the last element!
        
        for cut,txt,wisemag,wisemagname in [
                (has_gmag, ': 2MASS galaxies', gmag, 'gmag'),
                (None,'', None, None),
                (gals,': galaxies', None, None),
                (stars, ': point sources', None, None),
                                ]:
            if wisemag is None:
                wisemag = wisempro
    
            if cut is not None:
                x = wisemag[cut]
                y = tracmag[cut]
                print 'Cut', txt, ':', len(wisemag), '->', len(x)
            else:
                x = wisemag
                y = tracmag

            plt.subplots_adjust(left=0.12, **spa)

            plt.clf()
            loghist(x, y, range=((lo,cathi),(lo,cathi)), **ha)
            wmagname = 'mag'
            if wisemagname is not None:
                wmagname = wisemagname
            plt.xlabel('WISE W%i %s' % (band, wmagname))
            plt.ylabel('Tractor W%i mag' % band)

            # if allwise:
            #     plt.title('AllWISE')
            # else:
            #     plt.title('All-sky')
            #plt.title('WISE catalog vs Tractor forced photometry' + txt)
            plt.title('Tractor vs WISE' + txt)

            if cut is stars:
                magbins = np.arange(lo, cathi+2)
                xx, d = [],[]
                md = []
                
                if 'w%isigmpro' in W.get_columns():
                    we = W.get('w%isigmpro' % band)
                else:
                    we = F.get('w%i_mag_err' % band)
                    print 'Using Tractor errors'
                we = we[cut]

                for mlo,mhi in zip(magbins, magbins[1:]):
                    if band == 2 and mlo >= 17:
                        continue
                    I = np.flatnonzero((y > mlo) * (y <= mhi))
                    err = np.median(we[I])
                    xx.append((mlo+mhi)/2.)
                    d.append(err)
                    J = np.flatnonzero((y > mlo) * (y <= mhi) * (np.abs(x-y)<0.5))
                    md.append(np.median((y-x)[J]))

                plt.errorbar(xx, xx, yerr=d, fmt=None, ecolor='r', mew=1.5)
    
            plt.axis([cathi,lo,cathi,lo])

            ps.savefig()

        plt.subplots_adjust(left=0.16, **spa)
    
        #I = np.flatnonzero((x < 15) * (x > 10) * (np.abs(x - y) < 0.25))
        #med = -np.median((x-y)[I])
        #print 'Median', med

        ylo,yhi = -0.5,0.5
        #if band == 4:
        #    ylo,yhi = -0.5,1.0

        plt.clf()
        loghist(x, y - x, range=((lo,cathi),(ylo, yhi)), **ha)
        plt.xlabel('WISE W%i mag' % band)
        plt.ylabel('Tractor - WISE W%i mag' % band)
        #plt.errorbar(xx, np.zeros_like(xx), yerr=d, fmt=None, ecolor='r')
        #plt.axhline(med, color='r', alpha=0.5)
        plt.errorbar(xx, md, yerr=d, fmt='r.', ecolor='r', mew=1.5)
        plt.axhline(0., color='r', alpha=0.5)
        plt.axis([cathi,lo,ylo,yhi])

        # if allwise:
        #     plt.title('AllWISE')
        # else:
        #     plt.title('All-sky')

        ps.savefig()


        if allwise and Wallsky:
            I,J,d = match_radec(W.ra, W.dec, Wallsky.ra, Wallsky.dec, 4./3600.)
            v = Wallsky.get('w%impro' % band)
            plt.clf()
            loghist(wisempro[I], wisempro[I]-v[J], range=((lo,cathi),(ylo,yhi)), **ha)
            plt.xlabel('AllWISE W%i mag' % band)
            plt.ylabel('All-Sky - AllWISE W%i mag' % band)
            plt.axhline(0., color='r', alpha=0.5)
            plt.axis([cathi,lo,ylo,yhi])
            ps.savefig()



    return Wallsky

# 
#     # Tractor CMD
#     t2 = T.get('w%i_mag' % band)
#     ok2 = np.isfinite(t2)
#     t2 = t2[ok2]
#     rmag2 = T.modelflux[:,2]
#     rmag2 = rmag2[ok2]
#     rmag2 = -2.5 * (np.log10(rmag2) - 9.)
# 
#     plt.clf()
#     H,xe,ye = loghist(rmag2 - t2, rmag2, range=((-5,10),(12,25)), **ha)
#     plt.xlabel('r - W%i (mag)' % band)
#     plt.ylabel('r (mag)')
#     #plt.title('SDSS/WISE Tractor forced-photometry')
#     plt.axis([-5,10,25,12])
#     ps.savefig()
# 
#     # Catalog-match CMD
#     rmag = T.modelflux[:,2]
#     rmag = rmag[J][ok]
#     rmag = -2.5 * (np.log10(rmag) - 9.)
#     
#     plt.clf()
#     loghist(rmag - w, rmag, range=((-5,10),(12,25)),
#             imshowargs=dict(vmax=np.log10(np.max(H))), **ha)
#     plt.xlabel('r - W%i (mag)' % band)
#     plt.ylabel('r (mag)')
#     #plt.title('SDSS/WISE catalog matches')
#     plt.axis([-5,10,25,12])
#     ps.savefig()


    





if __name__ == '__main__':
    main()

sys.exit(0)

for band in []:#[1]: #,2,3,4]:
    coadd_id = '1384p454'
    for phase in ['a','d']:
        fn1 = 'c/phot-%s-%i%s.fits' % (coadd_id, band, phase)
        T1 = fits_table(fn1)
    
        W = fits_table('sequels-phot-temp/wise-sources-%s.fits' % coadd_id)
        print len(W), 'WISE'
    
        P = T1
        plt.clf()
        ok = np.isfinite(P.w1_mag)
        lo,hi = 10,25
        cathi = 18
        ha = dict(bins=100, histtype='step', range=(lo,hi), log=True)
        tsty = dict(color=(0.8,0.8,1.0), lw=3)
        csty = dict(color='b')
        a = ha.copy()
        a.update(tsty)
        n,b,p1 = plt.hist(P.w1_mag[ok], **a)
        a = ha.copy()
        a.update(csty)
        n,b,p2 = plt.hist(W.w1mpro, **a)
        # legend only
        p1 = plt.plot([1,1],[1,1], **tsty)
        p2 = plt.plot([1,1],[1,1], **csty)
        plt.xlabel('W1 mag (Vega)')
        plt.ylabel('Number of sources')
        plt.title('WISE catalog vs Tractor forced photometry depths')
        plt.legend((p1[0],p2[0]), ('W1 (Tractor)', 'W1 (WISE catalog)'), loc='lower right')
        plt.ylim(1., 2e4)
        plt.xlim(lo,hi)
        ps.savefig()
    
        I,J,d = match_radec(P.ra, P.dec, W.ra, W.dec, 4./3600.)
        print len(I), 'matches'
    
        ha = dict(bins=200,
                  )#imshowargs=dict(cmap=antigray), hot=False)
    
        plt.clf()
        loghist(W.w1mpro[J], P.w1_mag[I], range=((lo,cathi),(lo,cathi)), **ha)
        plt.xlabel('WISE W1 mag')
        plt.ylabel('Tractor W1 mag')
        plt.title('WISE catalog vs Tractor forced photometry')
        plt.axis([cathi,lo,cathi,lo])
        ps.savefig()
    
        plt.clf()
        P.r_mag = -2.5 * (np.log10(P.modelflux[:,2]) - 9.)
        loghist(P.r_mag - P.w1_mag, P.r_mag, range=((-5,10),(12,25)), **ha)
        plt.xlabel('r - W1 (mag)')
        plt.ylabel('r (mag)')
        plt.title('Tractor forced-photometered SDSS/WISE')
        ps.savefig()



for band in [1]:#,2,3,4]:
    coadd_id = '1384p454'

    pairs = [('a','e'), ('a','f'), ('a','g'), ('f','e')]

    for c1,c2 in pairs:
        fn1 = 'c/phot-%s-%i%s.fits' % (coadd_id, band, c1)
        T1 = fits_table(fn1)

        fn2 = 'c/phot-%s-%i%s.fits' % (coadd_id, band, c2)
        if not os.path.exists(fn1) or not os.path.exists(fn2):
            print 'not found:', fn1, fn2
            continue
        T1 = fits_table(fn1)
        T2 = fits_table(fn2)

        plt.clf()
        plt.plot(T1.get('w%i_nanomaggies' % band), T2.get('w%i_nanomaggies' % band), 'b.',
                 alpha=0.2)
        plt.xscale('symlog')
        plt.yscale('symlog')
        plt.title('%s - %s' % (fn1,fn2))
        ps.savefig()

sys.exit(0)

for band in [1,2,3,4]:
    fn1 = 'c/phot-1384p454-%ib.fits' % band
    fn2 = 'c/phot-1384p454-%ic.fits' % band
    if not os.path.exists(fn1) or not os.path.exists(fn2):
        print 'not found:', fn1, fn2
        continue
    T1 = fits_table(fn1)
    T2 = fits_table(fn2)

    plt.clf()
    plt.plot(T1.get('w%i_nanomaggies' % band), T2.get('w%i_nanomaggies' % band), 'b.',
             alpha=0.2)
    plt.xscale('symlog')
    plt.yscale('symlog')
    plt.title('%s - %s' % (fn1,fn2))
    ps.savefig()
    
