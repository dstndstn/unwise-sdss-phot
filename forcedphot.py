#! /usr/bin/env python

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os
import sys
from glob import glob
import fitsio


'''
Reviewer's questions:

python forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d vanilla
 -> 1800p000
python forcedphot.py --dataset sdss-dr10d 7536 --band 1 -d vanilla
 -> 2400p605

L1b: (on riemann)

python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d l1b --l1b --tiledir unwise-coadds-x/ > log 2>&1 &



PSF cutoff:

python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-5 --minrad 5 > log5 2>&1 &
python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-10 --minrad 10 > log10 2>&1 &
python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-10b --minrad 10 --minsig1 0.02 > log10b 2>&1 &

after rev 25401 (minRadius):

python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-c > logc 2>&1 &
python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-5c --minrad 5 > log5c 2>&1 &
python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-10c --minrad 10 > log10c 2>&1 &
python -u forcedphot.py --dataset sdss-dr10d 1578 --band 1 -d psfcutoff-10d --minrad 10 --minsig1 0.02 > log10d 2>&1 &


save_fits?

'''






# Notes to self:
#  -sort *-atlas.fits files by Dec and then RA, for better photoObj caching.
#  -W4 uncertainties

'''
See the file "NOTES" for details on the versions used for SEQUELS and
eBOSS targeting.

See the file "README" about the output products.
'''

'''
Relevant files/directories are:

DATASET-atlas.fits
  the coadd tiles to process

(--tiledir)
tiledir ("unwise-coadds")/xxx/xxxx[pm]xxx/unwise-xxxx[pm]xxx-wW-img-m.fits
                             and {invvar,std,n}-m.fits
  WISE coadd tiles

(--photoobjsdir) photoobjdir: photoObjs/301/R/C/...
  SDSS photoObj files

(--resolvedir) ("photoResolve")/window_flist.fits
  SDSS list of  files

(--tempdir) tempoutdir: ("DATASET-phot-temp")/photoobjs-TILE.fits
                                             /wise-sources-TILE.fits
  SDSS photoObjs files, and WISE catalog sources

(-d)
outdir ("DATASET-phot")/phot-TILE.fits
                       /phot-unsplit-TILE.fits
                       /phot-wise-TILE.fits
  phot-TILE.fits: WISE forced photometry for photoobjs-TILE.fits objects
  phot-unsplit-TILE.fits: objects that straddle tiles
  phot-wise-TILE.fits: photometry for WISE-only sources too

(--pobj)
pobjoutdir ("DATASET-pobj")/RERUN/RUN/CAMCOL/photoWiseForced-R-C-F.fits
  phot-TILE.fits entries split back into SDSS fields

photoobj-lengths.sqlite3
  Database holding the length of PhotoObj files


About --split:

With --finish --split:
   --will overwrite existing "pobj" output files!
   --will write OUTDIR/phot-unsplit-*.fits files


If you use --split when running tiles, it will still write out the
normal OUTDIR/phot-TILE.fits files, as well as
POBJ/RERUN/RUN/CAMCOL/photoWiseForced-R-C-F.fits files, and
OUTDIR/phot-unsplit-TILE.fits.  After all tiles have finished, one
must still run --finish on the phot-unsplit-TILE.fits files.


About version sdssv4:

cp sdss3-atlas.fits sdssv4-atlas.fits
# AllWISE, and photoObjs v5b
mv sdss3-phot-temp sdssv4-phot-temp
mkdir data/sdssv4-phot data/sdssv4-pobj data/sdssv4-logs
ln -s data/sdssv4-phot/ .
ln -s data/sdssv4-pobj/ .
ln -s data/sdssv4-logs/ .

module load astrometry_net/0.50
module load tractor/0.1

unwise-coadds -> data/unwise/unwise-comp
photoObjs-new -> /clusterfs/riemann/raid007/ebosswork/eboss/photoObj.v5b
photoResolve-new -> /clusterfs/riemann/raid006/bosswork/boss/resolve/2013-07-29
window_flist.fits -> /clusterfs/riemann/raid006/bosswork/boss/resolve/2013-07-29/window_flist.fits
unwise_coadd.py -> /home/dstn/unwise/unwise-coadd.py

seq 0 8342 | qdo load wphot-sdssv4 -
qdo launch wphot-sdssv4 1 --batchopts "-q batch -o sdssv4-logs -j oe -t 0-9 -l mem=6GB" --script="scripts/qdo-sdssv4.sh"

About version sdss-dr13:

# at rev 25233
python cut-window-flist.py
mkdir resolve-2013-07-29
cp cut-window-flist.fits resolve-2013-07-29/window_flist.fits
mv cut-window-flist.fits window_flist-cut.fits
python ~/unwise/sdss-footprint-wise.py
mv sdss-atlas.fits sdss-dr13-atlas.fits

mkdir data/sdss-dr13-phot data/sdss-dr13-pobj data/sdss-dr13-logs data/sdss-dr13-phot-temp
ln -s data/sdss-dr13-phot/ .
ln -s data/sdss-dr13-pobj/ .
ln -s data/sdss-dr13-logs/ .
ln -s data/sdss-dr13-phot-temp/ .
cd sdss-dr13-phot-temp
for x in ~/wphot/sdss-dr10d-phot-temp/wise-sources-*; do ln -s $x .; done
# rev 25245: grab photoObjs...
python -u forcedphot.py -B 8 --dataset sdss-dr13 --split -v -b 1234 --photoobjsdir /clusterfs/riemann/raid007/ebosswork/eboss/photoObj.v5b --resolvedir resolve-2013-07-29 --photoobjs-only 0-8038

Running tagged version of wphot: wphot-sdss-dr13
Tractor: v1.0-457-gbe2f2b9
Astrometry: 0.50-97-g64f30f5

# rev 25239
seq 0 8038 | qdo load wphot -
qdo launch wphot 1 --batchopts "-q batch -o sdss-dr13-logs -j oe -t 0-99 -l mem=4GB" --script="scripts/qdo-sdss-dr13.sh"

# When that finished, bump up memory request to 15 GB,
qdo retry wphot
qdo launch wphot 1 --batchopts "-q batch -o sdss-dr13-logs -j oe -t 0-49 -l mem=15GB" --script="scripts/qdo-sdss-dr13.sh"

python -u forcedphot.py --dataset sdss-dr13 -v -b 1234 --finish \
  --photoobjsdir /clusterfs/riemann/raid007/ebosswork/eboss/photoObj.v5b \
  --resolvedir resolve-2013-07-29 \
  sdss-dr13-phot/phot-unsplit-*.fits > fin-dr13.log 2>&1 &


About version sdss-dr10:

cp sdssv4-atlas.fits sdss-dr10-atlas.fits
mkdir data/sdss-dr10-logs data/sdss-dr10-phot data/sdss-dr10-phot-temp data/sdss-dr10-pobj
ln -s data/sdss-dr10-logs .
ln -s data/sdss-dr10-phot .
ln -s data/sdss-dr10-phot-temp .
ln -s data/sdss-dr10-pobj .

# sdssv3 and sdssv4 used the AllWISE catalogs; symlink
cd data/sdss-dr10-phot-temp; for x in ../sdss3-phot-temp/wise-sources-*.fits; do ln -s $x .; done

BOSS_PHOTOOBJ=/clusterfs/riemann/raid006/dr10/boss/photoObj
PHOTO_RESOLVE=/clusterfs/riemann/raid006/dr10/boss/resolve/2010-05-23

Start with:
module switch tree/dr10
python -u forcedphot.py --photoobjsdir $BOSS_PHOTOOBJ --resolvedir $PHOTO_RESOLVE --dataset sdss-dr10 --photoobjs-only 1-8343

seq 0 8342 | qdo load wphot-sdss-dr10 -
qdo launch wphot-sdss-dr10 1 --batchopts "-q batch -o sdss-dr10-logs -j oe -t 0-99 -l mem=4GB" --script="scripts/qdo-sdss-dr10.sh"

Running tagged version of wphot: wphot-sdss-dr10
Tractor        version: commit aa0dc8382435b6a1cdcccf8f8234e5c8c06d34bb Date:   Tue Sep 30 03:29:18 2014 -0700
Astrometry.net version: commit dcd92cc9e589dfb29c28d08f4cea8dea32e5c2be Date:   Mon Sep 29 11:16:49 2014 -0700

Midway through (bogged down by i/o), switched to:
python cut-window-flist.py
mkdir resolve-2010-05-23-cut
mv cut-window-flist.fits resolve-2010-05-23-cut/window_flist.fits

After that run (4 GB memory) finished, ran --todo, started a new qdo queue wphot-b,
bumped memory up to 16 GB,
qdo launch wphot-b 1 --batchopts "-q batch -o sdss-dr10-logs -j oe -t 0-49 -l mem=16GB" --script="scripts/qdo-sdss-dr10.sh"

That seems to have finished them all off!

python -u forcedphot.py --dataset sdss-dr10 -v -b 1234 --finish \
  --photoobjsdir /clusterfs/riemann/raid006/dr10/boss/photoObj \
  --resolvedir resolve-2010-05-23-cut \
  sdss-dr10-phot/phot-unsplit-*.fits > fin.log 2>&1 &

Then realized slight error in RA,Dec region for WISE-only catalog object search;
rm sdss-dr10-phot-temp/wise-sources-*
python -c "from astrometry.util.fits import *; import numpy as np; T=fits_table('sdss-dr10-atlas.fits');
  I=np.lexsort((T.ra, T.dec)); T[I].writeto('sdss-dr10b-atlas.fits')"
ln -s sdss-dr10-phot-temp/ sdss-dr10b-phot-temp
python -u forcedphot.py --dataset sdss-dr10b -v -b 1234 --wise-only > w.log 2>&1 &


About version sdss-dr10c:

-using new sdss-dr10c-atlas.fits (from unwise/sdss-footprint-wise.py)

mkdir data/sdss-dr10c-{logs,phot,pobj}
ln -s data/sdss-dr10c-logs .
ln -s data/sdss-dr10c-phot .
ln -s data/sdss-dr10c-pobj .

ln -s data/sdss-dr10-phot-temp sdss-dr10c-phot-temp

Running tagged version of wphot: wphot-sdss-dr10c
Tractor        version: commit 3b17cdb2bcd668b52de157bf06b6b247827c6299
Date:   Tue Oct 7 03:48:20 2014 -0700
Astrometry.net version: commit 45a088f0f5b6fdf3379078bb9f8920cb2768d361
Date:   Sun Oct 5 08:36:56 2014 -0400

seq 0 8038 | qdo load wphot -
qdo launch wphot 1 --batchopts "-q batch -o sdss-dr10c-logs -j oe -t 0-99 -l mem=4GB" --script="scripts/qdo-sdss-dr10c.sh"
(at svn rev 25161)

Succeeded  : 7522
Failed     : 517

qdo retry wphot
qdo launch wphot 1 --batchopts "-q batch -o sdss-dr10c-logs -j oe -t 0-99 -l mem=15GB" --script="scripts/qdo-sdss-dr10c.sh"
(at svn rev 25165)

# on big-mem node (actually, copied to bbq):
python -u forcedphot.py --dataset sdss-dr10c -v -b 1234 --finish \
  --photoobjsdir /clusterfs/riemann/raid006/dr10/boss/photoObj \
  --resolvedir resolve-2010-05-23-cut \
  sdss-dr10c-phot/phot-unsplit-*.fits > fin.log 2>&1 &

And finally, upon finding that the W4 errors were underestimated (darn
correlated errors from upsampling), I created "sdss-dr10d":

python -u fix-w4.py > fix.log 2>&1


'''

if __name__ == '__main__':
    d = os.environ.get('PBS_O_WORKDIR')
    if d is not None:
        os.chdir(d)
        sys.path.append(os.getcwd())

from astrometry.util.file import *
from astrometry.util.fits import *
from astrometry.util.multiproc import *
from astrometry.util.plotutils import *
from astrometry.util.miscutils import *
from astrometry.util.util import *
from astrometry.util.resample import *
from astrometry.libkd.spherematch import *
from astrometry.util.starutil_numpy import *
from astrometry.util.ttime import *
from astrometry.sdss import *

from tractor import *
from tractor.sdss import *

from wise.unwise import *
from wise.wisecat import wise_catalog_radecbox
from wise.allwisecat import allwise_catalog_radecbox

import logging

photoobjdir = None
resolvedir  = None

if __name__ == '__main__':
    tiledir = 'unwise-coadds'

    outdir = '%s-phot'
    tempoutdir = '%s-phot-temp'
    pobjoutdir = '%s-pobj'

    Time.add_measurement(MemMeas)

def read_wise_sources(wfn, wcs, extracols=[], allwise=True, margin=20):

    if allwise:
        catalog_radecbox = allwise_catalog_radecbox
    else:
        catalog_radecbox = wise_catalog_radecbox


    print 'looking for', wfn
    if os.path.exists(wfn):
        WISE = fits_table(wfn)
        print 'Read', len(WISE), 'WISE sources nearby'
    else:
        cols = ['ra','dec',
                'w1mpro','w2mpro', 'w3mpro','w4mpro',
                'w1sigmpro','w2sigmpro', 'w3sigmpro','w4sigmpro',
                'w1flux','w2flux', 'w3flux','w4flux',
                'w1sigflux','w2sigflux', 'w3sigflux','w4sigflux',
                'w1gmag','w1gerr', 'w2gmag','w2gerr',
                'w3gmag','w3gerr', 'w4gmag','w4gerr',
                ]
        cols += extracols

        r0,r1,d0,d1 = wcs.radec_bounds()
        print 'wcs bounds:', r0,r1,d0,d1
        margin = margin * wcs.pixel_scale() / 3600.
        cosdec = min([np.cos(np.deg2rad(x)) for x in [d0-margin, d1+margin]])
        rmargin = margin / cosdec
        print 'margin', margin
        print 'cosdec', cosdec
        # tan.radec_bounds on wrap-around returns the smallest < 360 and largest > 0,
        # and r0 > r1.
        d0 -= margin
        d1 += margin
        if d0-margin < -90 or d1 + margin > 90.:
            r0,r1 = 0, 360
        if r0 > r1:
            # wrap-around; glue together 0-r1 and r0-360
            Wa = catalog_radecbox(0., r1+rmargin, d0, d1, cols=cols)
            Wb = catalog_radecbox(r0-rmargin, 360., d0, d1, cols=cols)
            WW = []
            if Wa is not None:
                print 'Got', len(Wa), 'in RA 0 to', r1
                WW.append(Wa)
            if Wb is not None:
                print 'Got', len(Wb), 'in RA', r0, 'to 360'
                WW.append(Wb)
            WISE = merge_tables(WW)
        else:
            WISE = catalog_radecbox(r0-rmargin, r1+rmargin, d0, d1, cols=cols)
        WISE.writeto(wfn)
        print 'Found', len(WISE), 'WISE sources nearby; wrote to', wfn
    return WISE

photoobj_length_db = 'photoobj-lengths.sqlite3'
def get_photoobj_length(rerun, run, camcol, field, save=True, all=False):
    import sqlite3
    timeout = 60.
    create = not os.path.exists(photoobj_length_db)

    conn = sqlite3.connect(photoobj_length_db, timeout)
    c = conn.cursor()
    if create:
        print 'Creating db table in', photoobj_length_db
        c.execute('create table photoObjs (rerun text, run integer, ' +
                  'camcol integer, field integer, N integer)')
        conn.commit()

    if all:
        print 'Fetching all photoObj lengths'
        pobjs = {}
        for row in c.execute('select rerun, run, camcol, field, N from photoObjs'):
            rerun, run, camcol, field, N = row
            pobjs[(int(rerun), int(run), int(camcol), int(field))] = int(N)
        return pobjs

    print 'Getting photoObj length for', rerun, run, camcol, field
    #print type(rerun), type(run), type(camcol), type(field)

    # Ensure types (when read from FITS tables, they can be np.int16, eg)
    run = int(run)
    camcol = int(camcol)
    field = int(field)
    rerun = str(rerun)

    c.execute('select N from photoObjs where rerun=? and run=? and camcol=? '
              + 'and field=?', (rerun, run, camcol, field))
    row = c.fetchone()
    if row is None:
        # This photoObj is unknown
        sdss = DR9(basedir=photoobjdir)
        sdss.useLocalTree()
        pofn = sdss.retrieve('photoObj', run, camcol, field, rerun=rerun)
        #pofn = get_photoobj_filename(photoobjdir, rerun, run,camcol,field)
        F = fitsio.FITS(pofn)
        N = F[1].get_nrows()

        if save:
            c.execute('insert into photoObjs values (?,?,?,?,?)',
                      (rerun, run, camcol, field, N))
            conn.commit()
    else:
        print 'Row:', row
        N = row[0]
    conn.close()
    return N


def read_photoobjs(wcs, margin, cols=None):
    if cols is None:
        cols = ['objid', 'ra', 'dec', 'fracdev', 'objc_type', 'modelflux',
                'theta_dev', 'theta_deverr', 'ab_dev', 'ab_deverr', 'phi_dev_deg',
                'theta_exp', 'theta_experr', 'ab_exp', 'ab_experr', 'phi_exp_deg',
                'resolve_status', 'nchild', 'flags', 'objc_flags',
                'run','camcol','field','id'
                ]
        # useful to have these in the outputs...
        cols += ['psfflux', 'psfflux_ivar', 'cmodelflux', 'cmodelflux_ivar',
                 'modelflux', 'modelflux_ivar', 'raerr', 'decerr']
    wfn = os.path.join(resolvedir, 'window_flist.fits')
    sdss = DR9(basedir=photoobjdir)
    sdss.useLocalTree(photoObjs=photoobjdir, resolve=resolvedir)

    # in astrometry.sdss.fields
    return read_photoobjs_in_wcs(wcs, margin, cols=cols, wfn=wfn, sdss=sdss)

class BrightPointSource(PointSource):
    '''
    A class to use a pre-computed (constant) model, if available.
    '''
    def __init__(self, *args):
        super(BrightPointSource, self).__init__(*args)
        self.pixmodel = None
    def getUnitFluxModelPatch(self, *args, **kwargs):
        if self.pixmodel is not None:
            return self.pixmodel
        return super(BrightPointSource, self).getUnitFluxModelPatch(*args, **kwargs)

def set_bright_psf_mods(cat, WISE, T, brightcut, band, tile, wcs, sourcerad):
    mag = WISE.get('w%impro' % band)
    I = np.flatnonzero(mag < brightcut)
    if len(I) == 0:
        return
    BW = WISE[I]
    BW.nm = NanoMaggies.magToNanomaggies(mag[I])
    print len(I), 'catalog sources brighter than mag', brightcut
    I,J,d = match_radec(BW.ra, BW.dec, T.ra, T.dec, 4./3600., nearest=True)
    print 'Matched to', len(I), 'catalog sources (nearest)'
    if len(I) == 0:
        return

    fn = 'wise-psf-avg-pix-bright.fits'
    psfimg = fitsio.read(fn, ext=band-1).astype(np.float32)
    psfimg = np.maximum(0, psfimg)
    psfimg /= psfimg.sum()
    print 'PSF image', psfimg.shape
    print 'PSF image range:', psfimg.min(), psfimg.max()
    ph,pw = psfimg.shape
    pcx,pcy = ph/2, pw/2
    assert(ph == pw)
    phalf = ph/2

    ## HACK -- read an L1b frame to get the field rotation...
    thisdir = get_unwise_tile_dir(tiledir, tile.coadd_id)
    framesfn = os.path.join(thisdir, 'unwise-%s-w%i-frames.fits' % (tile.coadd_id, band))
    F = fits_table(framesfn)
    print 'intfn', F.intfn[0]
    #fwcs = fits_table(F.intfn[
    wisedir = 'wise-frames'
    scanid,frame = F.scan_id[0], F.frame_num[0]
    scangrp = scanid[-2:]
    fn = os.path.join(wisedir, scangrp, scanid, '%03i' % frame, 
                      '%s%03i-w%i-int-1b.fits' % (scanid, frame, band))
    fwcs = Tan(fn)
    # Keep CD matrix, set CRVAL/CRPIX to star position
    fwcs.set_crpix(pcx+1, pcy+1)
    fwcs.set_imagesize(float(pw), float(ph))

    for i,j in zip(I, J):
        if not isinstance(cat[j], BrightPointSource):
            print 'Bright source matched non-point source', cat[j]
            continue

        fwcs.set_crval(BW.ra[i], BW.dec[i])
        L=3
        Yo,Xo,Yi,Xi,rims = resample_with_wcs(wcs, fwcs, [psfimg], L)
        x0,x1 = int(Xo.min()), int(Xo.max())
        y0,y1 = int(Yo.min()), int(Yo.max())
        mod = np.zeros((1+y1-y0, 1+x1-x0), np.float32)
        mod[Yo-y0, Xo-x0] += rims[0]

        pat = Patch(x0, y0, mod)
        cat[j].pixmodel = pat

        cat[j].fixedRadius = phalf
        sourcerad[j] = max(sourcerad[j], phalf)

def treat_as_pointsource(T, bandnum, setObjcType=True):
    # Cut galaxies based on signal-to-noise of theta (effective
    # radius) measurement.
    b = bandnum
    gal = (T.objc_type == 3)
    dev = gal * (T.fracdev[:,b] >= 0.5)
    exp = gal * (T.fracdev[:,b] <  0.5)
    stars = (T.objc_type == 6)
    print sum(dev), 'deV,', sum(exp), 'exp, and', sum(stars), 'stars'
    print 'Total', len(T), 'sources'

    thetasn = np.zeros(len(T))
    T.theta_deverr[dev,b] = np.maximum(1e-6, T.theta_deverr[dev,b])
    T.theta_experr[exp,b] = np.maximum(1e-5, T.theta_experr[exp,b])
    # theta_experr nonzero: 1.28507e-05
    # theta_deverr nonzero: 1.92913e-06
    thetasn[dev] = T.theta_dev[dev,b] / T.theta_deverr[dev,b]
    thetasn[exp] = T.theta_exp[exp,b] / T.theta_experr[exp,b]

    aberrzero = np.zeros(len(T), bool)
    aberrzero[dev] = (T.ab_deverr[dev,b] == 0.)
    aberrzero[exp] = (T.ab_experr[exp,b] == 0.)

    maxtheta = np.zeros(len(T), bool)
    maxtheta[dev] = (T.theta_dev[dev,b] >= 29.5)
    maxtheta[exp] = (T.theta_exp[exp,b] >= 59.0)

    # theta S/N > modelflux for dev, 10*modelflux for exp
    bigthetasn = (thetasn > (T.modelflux[:,b] * (1.*dev + 10.*exp)))

    print sum(gal * (thetasn < 3.)), 'have low S/N in theta'
    print sum(gal * (T.modelflux[:,b] > 1e4)), 'have big flux'
    print sum(aberrzero), 'have zero a/b error'
    print sum(maxtheta), 'have the maximum theta'
    print sum(bigthetasn), 'have large theta S/N vs modelflux'
    
    badgals = gal * reduce(np.logical_or,
                           [thetasn < 3.,
                            T.modelflux[:,b] > 1e4,
                            aberrzero,
                            maxtheta,
                            bigthetasn,
                            ])
    print 'Found', sum(badgals), 'bad galaxies'
    if setObjcType:
        T.objc_type[badgals] = 6
    return badgals

def _get_photoobjs(tile, wcs, bandnum, existOnly):
    objfn = os.path.join(tempoutdir, 'photoobjs-%s.fits' % tile.coadd_id)
    T = None
    if os.path.exists(objfn):
        if existOnly:
            print 'Exists:', objfn
            return
        print 'Reading', objfn
        try:
            T = fits_table(objfn)
        except:
            print 'Failed to read', objfn, '-- will re-create'
            import traceback
            traceback.print_exc()
            T = None
    if T is None:
        print 'Did not find', objfn, '-- reading photoObjs'
        T = read_photoobjs(wcs, 1./60.)
        if T is None:
            return None
        T.writeto(objfn)
        print 'Wrote', objfn
        if existOnly:
            return

    # record "pointsource" before computing treat_as_pointsource; it modifies objc_type
    T.pointsource = (T.objc_type == 6)
    badgals = treat_as_pointsource(T, bandnum)
    T.treated_as_pointsource = badgals
    return T


def _unwise_l1b_tractor_images(unwdir, l1bdir, coadd_id, band, bandname,
                               psffn, l1bsky):
    from unwise_coadd import walk_wcs_boundary, zeropointToScale

    # All this is required to get the L1b image extents
    # (thanks to a typo in the unWISE coadd code)
    coimfn = os.path.join(unwdir, 'unwise-%s-w%i-img-m.fits' %
                          (coadd_id, band))
    cowcs = Tan(coimfn)
    # Intermediate world coordinates (IWC) polygon
    W = cowcs.get_width()
    r,d = walk_wcs_boundary(cowcs, step=W, margin=0)
    ok,u,v = cowcs.radec2iwc(r,d)
    copoly = np.array(list(reversed(zip(u,v))))

    sky = 0.
    if l1bsky:
        hdr = fitsio.read_header(coimfn)
        sky = hdr['UNW_SKY']
        print 'Read sky value of', sky, 'from coadd header'

    print 'Reading PSF from', psffn
    P = fits_table(psffn, hdu=band)
    psf = GaussianMixturePSF(P.amp, P.mean, P.var)
    
    framefn = os.path.join(unwdir, 'unwise-%s-w%i-frames.fits' %
                           (coadd_id, band))
    F = fits_table(framefn)
    print len(F), 'frames'
    F.cut(F.included == 1)
    print 'Cut to', len(F)
    tims = []
    for f in F:
        intfn = os.path.join(l1bdir, f.intfn)
        uncfn = intfn.replace('-int-', '-unc-')
        uncfn = uncfn + '.gz'
        maskfn = intfn.replace('-int-', '-msk-')
        maskfn = maskfn + '.gz'
        print
        print 'Image', intfn
        if not all([os.path.exists(fn) for fn in [intfn, uncfn, maskfn]]):
            print 'Does not exist:', intfn
            cmd = ('(wget -r -N -nH -np -nv --cut-dirs=4 -A "*w%i*" "http://irsa.ipac.caltech.edu/ibe/data/wise/merge/merge_p1bm_frm/%s")' %
                   (band, os.path.dirname(f.intfn.replace('wise-frames/', ''))+'/'))
            print cmd
            os.system(cmd)


        # Image extent...
        wcs = Sip(intfn)
        h,w = wcs.get_height(), wcs.get_width()
        r,d = walk_wcs_boundary(wcs, step=2.*w, margin=0)
        ok,u,v = cowcs.radec2iwc(r, d)
        poly = np.array(list(reversed(zip(u,v))))
        intersects = polygons_intersect(copoly, poly)
        if not intersects:
            print 'Image does not intersect target'
            print 'Coadd poly:', copoly
            print 'Image poly:', poly
            print 'Clipped:', clip_polygon(copoly, poly)
            continue
        cpoly = np.array(clip_polygon(copoly, poly))
        if len(cpoly) == 0:
            print 'No overlap between coadd and image polygons'
            continue
        rd = np.array([cowcs.iwc2radec(u,v) for u,v in cpoly])
        ok,x,y = np.array(wcs.radec2pixelxy(rd[:,0], rd[:,1]))
        x -= 1
        y -= 1
        x0,y0 = [np.floor(v.min(axis=0)).astype(int) for v in [x,y]]
        x1,y1 = [np.ceil (v.max(axis=0)).astype(int) for v in [x,y]]
        imextent = [np.clip(x0, 0, w-1),
                    np.clip(x1, 0, w-1),
                    np.clip(y0, 0, h-1),
                    np.clip(y1, 0, h-1)]
        #print 'imextent', imextent

        # inclusive
        x0,x1,y0,y1 = [int(x) for x in imextent]
        x1 += 1
        y1 += 1
        slc = (slice(y0, y1), slice(x0, x1))

        img = fitsio.FITS(intfn)[0][slc]
        #unc = fitsio.FITS(uncfn)[0][slc]
        mask = fitsio.FITS(maskfn)[0][slc]

        print 'Read', img.shape, img.dtype, 'image'
        print 'Image median', np.median(img), 'vs sky', f.sky1
        
        iv = f.weight
        zp = f.zeropoint
        zpscale = 1. / zeropointToScale(zp)
        #print 'Zeropoint:', zp, '-> scale', zpscale

        badbits = [0,1,2,3,4,5,6,7, 9, 
                   10,11,12,13,14,15,16,17,18,
                   21,26,27,28]
        if f.phase == 3:
            # 3-band cryo phase:
            ## 19 pixel is "hard-saturated"
            ## 23 for W3 only: static-split droop residual present
            badbits.append(19)
            if band == 3:
                badbits.append(23)

        maskbits = sum([1<<bit for bit in badbits])
        goodmask = ((mask & maskbits) == 0)
        #goodmask[unc == 0] = False
        goodmask[np.logical_not(np.isfinite(img))] = False
        #goodmask[np.logical_not(np.isfinite(unc))] = False

        invvar = np.zeros_like(img)
        invvar[goodmask] = iv

        # also mask out regions outside the coadd frame.
        #print 'Clipped polygon:', cpoly
        xx,yy = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1))
        r,d = wcs.pixelxy2radec(xx, yy)
        ok,u,v = cowcs.radec2iwc(r, d)
        inside = point_in_poly(u, v, cpoly)
        print np.sum(inside), 'of', xx.size, 'pixels are inside coadd bounds'
        #print 'inside:', inside.shape
        invvar[np.logical_not(inside)] = 0.

        # Add unWISE mask
        fn = os.path.join(unwdir, 
                          'unwise-%s-w%i-mask' % (coadd_id, band),
                          'unwise-mask-%s-%s%03i-w%i-1b.fits.gz' %
                          (coadd_id, f.scan_id, f.frame_num, band))
        print 'Looking for unWISE mask', fn
        if os.path.exists(fn):
            umask = fitsio.FITS(fn)[0][slc]
        else:
            tgzfn = os.path.join(unwdir, 
                                 'unwise-%s-w%i-mask.tgz' % (coadd_id, band))
            print 'Looking for tgz', tgzfn
            maskdir = 'unwise-%s-w%i-mask' % (coadd_id, band)
            maskfn = maskdir + ('/unwise-mask-%s-%s%03i-w%i-1b.fits.gz' %
                                (coadd_id, f.scan_id, f.frame_num, band))
            # extract in temp dir
            import tempfile
            tempdir = tempfile.mkdtemp()
            cmd = 'tar xf %s -C %s %s' % (tgzfn, tempdir, maskfn)
            print cmd
            from astrometry.util.run_command import run_command
            rtn,txt,err = run_command(cmd)
            if rtn:
                raise RuntimeError('Failed to untar mask file')
            print txt
            umask = fitsio.FITS(os.path.join(tempdir, maskfn))[0][slc]

            os.unlink(os.path.join(tempdir, maskfn))
            os.rmdir(os.path.join(tempdir, maskdir))
            os.rmdir(tempdir)

        n0 = np.sum(invvar > 0)
        invvar[umask > 0] = 0.
        print 'Masked', n0 - np.sum(invvar > 0), 'pixels from unWISE'

        img[invvar == 0] = 0
        assert(np.all(np.isfinite(img)))
        assert(np.all(np.isfinite(invvar)))
        
        img -= (f.sky1 + f.sky2)
        # -> nanomaggies
        img *= zpscale
        # the invvar is *already* in nanomaggies units.

        if sky != 0:
            print 'Subtracting off coadd sky level of', sky
            img -= sky

        # WCS for the sub-image...
        h,w = img.shape
        wcs = wcs.get_subimage(x0, y0, w, h)
        
        twcs = ConstantFitsWcs(wcs)
        sky = 0.
        tsky = ConstantSky(sky)

        tim = Image(data=img, invvar=invvar, psf=psf, wcs=twcs,
                    sky=tsky, photocal=LinearPhotoCal(1., band=bandname),
                    name='WISE L1b %s-%i W%i' % (f.scan_id, f.frame_num, band))
        tim.sig1 = np.sqrt(1./iv)
        tims.append(tim)
    return tims
        
def _unwise_tractor_image(thisdir, coadd_id, band, bandname, psffn):
    '''
    bandname: what to call the band in the PhotoCal object; must match name used
        in the source objects.
    '''
    imfn = os.path.join(thisdir, 'unwise-%s-w%i-img-m.fits'       %
                        (coadd_id, band))
    ivfn = os.path.join(thisdir, 'unwise-%s-w%i-invvar-m.fits.gz' %
                        (coadd_id, band))
    #ppfn = os.path.join(thisdir, 'unwise-%s-w%i-std-m.fits.gz'    %
    #                    (coadd_id, band))
    nifn = os.path.join(thisdir, 'unwise-%s-w%i-n-m.fits.gz'      %
                        (coadd_id, band))

    print 'Reading', imfn
    wcs = Tan(imfn)
    #ra,dec = wcs.radec_center()
    img = fitsio.read(imfn)
    print 'Reading', ivfn
    invvar = fitsio.read(ivfn)
    #print 'Reading', ppfn
    #pp = fitsio.read(ppfn)
    print 'Reading', nifn
    nims = fitsio.read(nifn)
    print 'Median # ims:', np.median(nims)

    assert(np.all(np.isfinite(img)))
    assert(np.all(np.isfinite(invvar)))
    #assert(np.all(np.isfinite(pp)))

    good = (nims > 0)
    invvar[np.logical_not(good)] = 0.

    sig1 = 1./np.sqrt(np.median(invvar[good]))
    assert(sig1 > 0)

    twcs = ConstantFitsWcs(wcs)
    sky = 0.
    tsky = ConstantSky(sky)

    # Load the average PSF model (generated by wise_psf.py)
    print 'Reading PSF from', psffn
    P = fits_table(psffn, hdu=band)
    psf = GaussianMixturePSF(P.amp, P.mean, P.var)
    
    tim = Image(data=img, invvar=invvar, psf=psf, wcs=twcs,
                sky=tsky, photocal=LinearPhotoCal(1., band=bandname),
                name='Coadd %s W%i' % (coadd_id, band))
    tim.nims = nims
    tim.sig1 = sig1
    
    return tim

def one_tile(tile, opt, savepickle, ps, tiles, tiledir, tempoutdir,
             T=None, hdr=None):

    bands = opt.bands
    outfn = opt.output % (tile.coadd_id)
    savewise_outfn = opt.save_wise_output % (tile.coadd_id)

    sband = 'r'
    bandnum = 'ugriz'.index(sband)

    tt0 = Time()
    print
    print 'Coadd tile', tile.coadd_id

    thisdir = get_unwise_tile_dir(tiledir, tile.coadd_id)
    fn = os.path.join(thisdir, 'unwise-%s-w%i-img-m.fits' % (tile.coadd_id, bands[0]))
    print 'Reading', fn
    wcs = Tan(fn)
    racenter,deccenter = wcs.radec_center()
    H,W = wcs.get_height(), wcs.get_width()

    if T is None:
        T = _get_photoobjs(tile, wcs, bandnum, opt.photoObjsOnly)
        if T is None:
            print 'Empty tile'
            return
        if opt.photoObjsOnly:
            return
    print len(T), 'objects'
    if len(T) == 0:
        return

    defaultflux = 1.

    # hack
    T.psfflux    = np.zeros((len(T),5), np.float32) + defaultflux
    T.cmodelflux = T.psfflux
    T.devflux    = T.psfflux
    T.expflux    = T.psfflux

    ok,T.x,T.y = wcs.radec2pixelxy(T.ra, T.dec)
    T.x = (T.x - 1.).astype(np.float32)
    T.y = (T.y - 1.).astype(np.float32)
    margin = 20.
    I = np.flatnonzero((T.x >= -margin) * (T.x < W+margin) *
                       (T.y >= -margin) * (T.y < H+margin))
    T.cut(I)
    print 'Cut to margins: N objects:', len(T)
    if len(T) == 0:
        return

    # Use pixelized PSF models for bright sources?
    bright_mods = ((1 in bands) and (opt.bright1 is not None))

    wanyband = wband = 'w'

    classmap = {}
    if bright_mods:
        classmap = {PointSource: BrightPointSource}

    print 'Creating tractor sources...'
    cat = get_tractor_sources_dr9(None, None, None, bandname=sband, objs=T,
                                  bands=[], nanomaggies=True, extrabands=[wband],
                                  fixedComposites=True, useObjcType=True,
                                  classmap=classmap)
    print 'Created', len(T), 'sources'
    assert(len(cat) == len(T))
    T.delete_column('index')

    pixscale = wcs.pixel_scale()
    # crude intrinsic source radii, in pixels
    sourcerad = np.zeros(len(cat))
    for i in range(len(cat)):
        src = cat[i]
        if isinstance(src, PointSource):
            continue
        elif isinstance(src, HoggGalaxy):
            sourcerad[i] = (src.nre * src.shape.re / pixscale)
        elif isinstance(src, FixedCompositeGalaxy):
            sourcerad[i] = max(src.shapeExp.re * ExpGalaxy.nre,
                               src.shapeDev.re * DevGalaxy.nre) / pixscale
    print 'sourcerad range:', min(sourcerad), max(sourcerad)

    # Find WISE-only catalog sources
    wfn = os.path.join(tempoutdir, 'wise-sources-%s.fits' % (tile.coadd_id))
    # expand search region by margin
    WISE = read_wise_sources(wfn, wcs)

    for band in bands:
        mag = WISE.get('w%impro' % band)
        nm = NanoMaggies.magToNanomaggies(mag)
        WISE.set('w%inm' % band, nm)
        print 'Band', band, 'max WISE catalog flux:', max(nm)
        print '  (min mag:', mag.min(), ')'

    unmatched = np.ones(len(WISE), bool)
    I,J,d = match_radec(WISE.ra, WISE.dec, T.ra, T.dec, 4./3600.)
    unmatched[I] = False
    UW = WISE[unmatched]
    print 'Got', len(UW), 'unmatched WISE sources'
    del unmatched, d
    # (I and J are used in this next block...)
    
    # Record WISE fluxes for catalog matches.
    # (this provides decent initialization for 'minsb' approx.)
    wiseflux = {}
    for band in bands:
        wiseflux[band] = np.zeros(len(T))
        if len(I) == 0:
            continue
        # X[I] += Y[J] with duplicate I doesn't work.
        #wiseflux[band][J] += WISE.get('w%inm' % band)[I]
        lhs = wiseflux[band]
        rhs = WISE.get('w%inm' % band)[I]
        print 'Band', band, 'max matched WISE flux:', max(rhs)
        for j,f in zip(J, rhs):
            lhs[j] += f

    ok,UW.x,UW.y = wcs.radec2pixelxy(UW.ra, UW.dec)
    UW.x -= 1.
    UW.y -= 1.

    if opt.savewise:
        fitwiseflux = {}
        for band in bands:
            fitwiseflux[band] = np.zeros(len(UW))
    
    T.coadd_id = np.array([tile.coadd_id] * len(T))

    inbounds = np.flatnonzero((T.x >= -0.5) * (T.x < W-0.5) *
                              (T.y >= -0.5) * (T.y < H-0.5))

    print 'Before looping over bands:', Time()-tt0
   
    for band in bands:
        tb0 = Time()
        print
        print 'Coadd tile', tile.coadd_id
        print 'Band', band
        wband = 'w%i' % band

        if opt.l1b:
            tims = _unwise_l1b_tractor_images(thisdir, '.', tile.coadd_id,
                                              band, wanyband, opt.psffn, opt.l1b_sky)
            tim = tims[0]
        else:
            tim = _unwise_tractor_image(thisdir, tile.coadd_id, band, wanyband,
                                        opt.psffn)
            tims = [tim]
            
        # Surface-brightness approximation
        minsig = getattr(opt, 'minsig%i' % band)
        sig1 = tim.sig1
        minsb = sig1 * minsig
        print 'Sigma1:', sig1, 'minsig', minsig, 'minsb', minsb

        # Render the PSF profile for figuring out source radii for
        # approximation purposes.
        R = 100
        psf = tim.psf
        psf.radius = R
        pat = psf.getPointSourcePatch(0., 0.)
        assert(pat.x0 == pat.y0)
        assert(pat.x0 == -R)
        psfprofile = pat.patch[R, R:]
        #print 'PSF profile:', psfprofile

        # Reset default flux based on min radius
        defaultflux = minsb / psfprofile[opt.minradius]
        print 'Setting default flux', defaultflux

        # Set WISE-only source radii based on flux
        UW.rad = np.zeros(len(UW), int)
        wnm = UW.get('w%inm' % band)
        for r,pro in enumerate(psfprofile):
            if pro == 0:
                continue
            flux = minsb / pro
            UW.rad[wnm > flux] = r
        UW.rad = np.maximum(UW.rad + 1, 3)

        # Set fluxes of SDSS objects based on WISE catalog matches.
        wf = wiseflux[band]
        I = np.flatnonzero(wf > defaultflux)
        wfi = wf[I]
        print 'Initializing', len(I), 'fluxes based on catalog matches'
        for i,flux in zip(I, wf[I]):
            assert(np.isfinite(flux))
            cat[i].getBrightness().setBand(wanyband, flux)

        # Set radii of SDSS objects based on WISE flux
        rad = np.zeros(len(I), int)
        for r,pro in enumerate(psfprofile):
            if pro == 0:
                continue
            flux = minsb / pro
            rad[wfi > flux] = r
        srad2 = np.zeros(len(cat), int)
        srad2[I] = rad
        del rad
        del psfprofile

        # Set radii
        for i in range(len(cat)):
            src = cat[i]
            # set fluxes
            b = src.getBrightness()
            if b.getBand(wanyband) <= defaultflux:
                b.setBand(wanyband, defaultflux)
                
            R = max([opt.minradius, sourcerad[i], srad2[i]])
            # "sourcerad" is used to select which sources are in-range
            sourcerad[i] = R
            if isinstance(src, PointSource):
                src.fixedRadius = R
                #src.minradius = opt.minradius
                src.minRadius = opt.minradius
                
            elif (isinstance(src, HoggGalaxy) or
                  isinstance(src, FixedCompositeGalaxy)):
                src.halfsize = R
                
        # Use pixelized PSF models for bright sources?
        bright_mods = ((band == 1) and (opt.bright1 is not None))

        if bright_mods:
            set_bright_psf_mods(cat, WISE, T, opt.bright1, band, tile, wcs, sourcerad)

        flux_invvars = np.zeros(len(cat))
        fskeys = ['prochi2', 'pronpix', 'profracflux', 'proflux', 'npix']
        if not opt.l1b:
            fskeys.append('pronexp')
        fitstats = dict([(k, np.zeros(len(cat))) for k in fskeys])

        if ps:
            tag = '%s W%i' % (tile.coadd_id, band)
            img = tim.getImage()
            sky = tim.getSky().getValue()
            plt.clf()
            n,b,p = plt.hist(img.ravel(), bins=100,
                             range=(-10*sig1, 20*sig1), log=True,
                             histtype='step', color='b')
            mx = max(n)
            plt.ylim(0.1, mx)
            plt.xlim(-10*sig1, 20*sig1)
            plt.axvline(sky, color='r')
            plt.title('%s: Pixel histogram' % tag)
            ps.savefig()

            if bright_mods:
                mod = np.zeros_like(img)
                for src in cat:
                    if src.pixmodel:
                        src.pixmodel.add(mod, scale=src.getBrightness().getBand(wanyband))

                plt.clf()
                plt.imshow(mod, interpolation='nearest', origin='lower',
                           cmap='gray',
                           vmin=-3*sig1, vmax=10*sig1)
                plt.colorbar()
                plt.title('%s: bright star models' % tag)
                ps.savefig()
    
                plt.clf()
                plt.imshow(img - mod, interpolation='nearest', origin='lower',
                           cmap='gray',
                           vmin=-3*sig1, vmax=10*sig1)
                plt.colorbar()
                plt.title('%s: data - bright star models' % tag)
                ps.savefig()


        if savepickle:
            mods = []
            cats = []

        # SDSS and WISE source margins beyond the image margins (+ source radii )
        smargin = 1
        wmargin = 1

        # Relevant SDSS sources:
        m = smargin + sourcerad
        I = np.flatnonzero(((T.x+m) >= -0.5) * ((T.x-m) < (W-0.5)) *
                           ((T.y+m) >= -0.5) * ((T.y-m) < (H-0.5)))
        inbox = ((T.x[I] >= -0.5) * (T.x[I] < (W-0.5)) *
                 (T.y[I] >= -0.5) * (T.y[I] < (H-0.5)))
        # Source indices inside the box
        srci = I[inbox]
        # Source indices in the margins
        margi = I[np.logical_not(inbox)]

        # sources in the box
        subcat = [cat[i] for i in srci]

        # include *copies* of sources in the margins
        # (that way we automatically don't save the results)
        subcat.extend([cat[i].copy() for i in margi])
        assert(len(subcat) == len(I))

        # add WISE-only sources in the margins
        m = wmargin + UW.rad
        J = np.flatnonzero(((UW.x+m) >= -0.5) * ((UW.x-m) < (W-0.5)) *
                           ((UW.y+m) >= -0.5) * ((UW.y-m) < (H-0.5)))

        if opt.savewise:
            jinbox = ((UW.x[J] >= -0.5) * (UW.x[J] < (W-0.5)) *
                      (UW.y[J] >= -0.5) * (UW.y[J] < (H-0.5)))
            uwcat = []
        wnm = UW.get('w%inm' % band)
        nomag = 0
        for ji,j in enumerate(J):
            if not np.isfinite(wnm[j]):
                nomag += 1
                continue
            ptsrc = PointSource(RaDecPos(UW.ra[j], UW.dec[j]),
                                      NanoMaggies(**{wanyband: wnm[j]}))
            ptsrc.radius = UW.rad[j]
            subcat.append(ptsrc)
            if opt.savewise:
                if jinbox[ji]:
                    uwcat.append((j, ptsrc))
                
        print 'WISE-only:', nomag, 'of', len(J), 'had invalid mags'
        print 'Sources:', len(srci), 'in the box,', len(I)-len(srci), 'in the margins, and', len(J), 'WISE-only'
        print 'Creating a Tractor with images', [t.shape for t in tims], 'and', len(subcat), 'sources'
        tractor = Tractor(tims, subcat)

        #### TEST
        if ps and opt.l1b:
            # Coadd data and models.
            codat = np.zeros((H,W), np.float32)
            coiv  = np.zeros((H,W), np.float32)
            for tim in tims:
                dat = tim.getImage()
                ie  = tim.getInvError()
                try:
                    Yo,Xo,Yi,Xi,nil = resample_with_wcs(wcs, tim.wcs.wcs, [],3)
                except:
                    continue
                if len(Yo) == 0:
                    continue
                iv = ie[Yi,Xi]**2
                codat[Yo,Xo] += dat[Yi,Xi] * iv
                coiv [Yo,Xo] += iv
            codat /= coiv
            codat[coiv == 0] = 0.
            dat = codat

            plt.clf()
            plt.imshow(dat, interpolation='nearest', origin='lower',
                       cmap='gray', vmin=-3*sig1, vmax=10*sig1)
            plt.colorbar()
            plt.title('Coadded data')
            ps.savefig()




        print 'Running forced photometry...'
        t0 = Time()
        tractor.freezeParamsRecursive('*')

        if opt.sky:
            tractor.thawPathsTo('sky')
            print 'Initial sky values:'
            for tim in tractor.getImages():
                print tim.getSky()

        tractor.thawPathsTo(wanyband)

        wantims = (savepickle or (ps is not None) or opt.save_fits)

        kwa = {}
        if not opt.l1b:
            kwa.update(fitstat_extras=[('pronexp', [tim.nims])])

        if opt.ceres:
            from tractor.ceres_optimizer import CeresOptimizer
            tractor.optimizer = CeresOptimizer(BW=opt.ceresblock,
                                               BH=opt.ceresblock)

        R = tractor.optimize_forced_photometry(
            minsb=minsb, mindlnp=1., sky=opt.sky, minFlux=None,
            fitstats=True, 
            variance=True, shared_params=False,
            wantims=wantims, negfluxval=0.1*sig1, **kwa)
        print 'That took', Time()-t0

        if opt.ceres:
            #print 'Ceres status:', R.ceres_status
            #'termination': 2
            term = R.ceres_status['termination']
            print 'Ceres termination status:', term
            # Running out of memory can cause failure to converge and term status = 2.
            # Fail completely in this case.
            if term != 0:
                raise RuntimeError('Ceres terminated with status %i' % term)

        if wantims:
            ims0 = R.ims0
            ims1 = R.ims1
        IV,fs = R.IV, R.fitstats

        if opt.sky:
            print 'Fit sky values:'
            for tim in tractor.getImages():
                print tim.getSky()

        if opt.savewise:
            for (j,src) in uwcat:
                fitwiseflux[band][j] = src.getBrightness().getBand(wanyband)

        if opt.save_fits:
            (dat,mod,ie,chi,roi) = ims1[0]

            tag = 'fit-%s-w%i' % (tile.coadd_id, band)
            pat = os.path.join(opt.outdir, '%s-%%s.fits' % tag)
            for name,data in [('data',dat),('mod',mod),('chi',chi)]:
                fn = pat % name
                fitsio.write(fn, data, clobber=True)
                print 'Wrote', fn

        if ps:
            tag = '%s W%i' % (tile.coadd_id, band)

        if ps and not opt.l1b:
            (dat,mod,ie,chi,roi) = ims1[0]
        if ps and opt.l1b:
            # Coadd data and models.
            codat = np.zeros((H,W), np.float32)
            coiv  = np.zeros((H,W), np.float32)
            comod = np.zeros((H,W), np.float32)

            assert(len(tims) == len(ims1))

            for tim,(dat,mod,ie,chi,roi) in zip(tims, ims1):
                try:
                    Yo,Xo,Yi,Xi,nil = resample_with_wcs(wcs, tim.wcs.wcs, [],3)
                except:
                    continue
                if len(Yo) == 0:
                    continue
                iv = ie[Yi,Xi]**2
                codat[Yo,Xo] += dat[Yi,Xi] * iv
                coiv [Yo,Xo] += iv
                comod[Yo,Xo] += mod[Yi,Xi] * iv

            codat /= coiv
            comod /= coiv
            codat[coiv == 0] = 0.
            comod[coiv == 0] = 0.
            dat = codat
            mod = comod
            chi = (codat - comod) * np.sqrt(coiv)

        if ps:
            plt.clf()
            plt.imshow(dat, interpolation='nearest', origin='lower',
                       cmap='gray', vmin=-3*sig1, vmax=10*sig1)
            plt.colorbar()
            plt.title('%s: data' % tag)
            ps.savefig()

            # plt.clf()
            # plt.imshow(1./ie, interpolation='nearest', origin='lower',
            #            cmap='gray', vmin=0, vmax=10*sig1)
            # plt.colorbar()
            # plt.title('%s: sigma' % tag)
            # ps.savefig()

            plt.clf()
            plt.imshow(mod, interpolation='nearest', origin='lower',
                       cmap='gray', vmin=-3*sig1, vmax=10*sig1)
            plt.colorbar()
            plt.title('%s: model' % tag)
            ps.savefig()

            plt.clf()
            plt.imshow(chi, interpolation='nearest', origin='lower',
                       cmap='gray', vmin=-5, vmax=+5)
            plt.colorbar()
            plt.title('%s: chi' % tag)
            ps.savefig()

            # plt.clf()
            # plt.imshow(np.round(chi), interpolation='nearest', origin='lower',
            #            cmap='jet', vmin=-5, vmax=+5)
            # plt.colorbar()
            # plt.title('Chi')
            # ps.savefig()

            plt.clf()
            plt.imshow(chi, interpolation='nearest', origin='lower',
                       cmap='gray', vmin=-20, vmax=+20)
            plt.colorbar()
            plt.title('%s: chi 2' % tag)
            ps.savefig()

            plt.clf()
            n,b,p = plt.hist(chi.ravel(), bins=100,
                             range=(-10, 10), log=True,
                             histtype='step', color='b')
            mx = max(n)
            plt.ylim(0.1, mx)
            plt.axvline(0, color='r')
            plt.title('%s: chi' % tag)
            ps.savefig()

            # fn = ps.basefn + '-chi.fits'
            # fitsio.write(fn, chi, clobber=True)
            # print 'Wrote', fn

        if savepickle:
            if ims1 is None:
                mod = None
            else:
                im,mod,ie,chi,roi = ims1[0]
            mods.append(mod)
            cats.append((
                srci, margi, UW.x[J], UW.y[J],
                T.x[srci], T.y[srci], T.x[margi], T.y[margi],
                [src.copy() for src in cat],
                [src.copy() for src in subcat]))

        if len(srci):
            # Save fit stats
            flux_invvars[srci] = IV[:len(srci)]
            for k in fskeys:
                x = getattr(fs, k)
                fitstats[k][srci] = np.array(x)

        if bright_mods:
            # Reset pixelized models
            for src in cat:
                if isinstance(src, BrightPointSource):
                    src.pixmodel = None

        nm = np.array([src.getBrightness().getBand(wanyband) for src in cat])
        nm_ivar = flux_invvars
        T.set(wband + '_nanomaggies', nm.astype(np.float32))
        T.set(wband + '_nanomaggies_ivar', nm_ivar.astype(np.float32))
        dnm = np.zeros(len(nm_ivar), np.float32)
        okiv = (nm_ivar > 0)
        dnm[okiv] = (1./np.sqrt(nm_ivar[okiv])).astype(np.float32)
        okflux = (nm > 0)
        mag = np.zeros(len(nm), np.float32)
        mag[okflux] = (NanoMaggies.nanomaggiesToMag(nm[okflux])).astype(np.float32)
        dmag = np.zeros(len(nm), np.float32)
        ok = (okiv * okflux)
        dmag[ok] = (np.abs((-2.5 / np.log(10.)) * dnm[ok] / nm[ok])).astype(np.float32)

        mag[np.logical_not(okflux)] = np.nan
        dmag[np.logical_not(ok)] = np.nan
        
        T.set(wband + '_mag', mag)
        T.set(wband + '_mag_err', dmag)
        for k in fskeys:
            T.set(wband + '_' + k, fitstats[k].astype(np.float32))

        if ps:
            I,J,d = match_radec(WISE.ra, WISE.dec, T.ra, T.dec, 4./3600.)

            plt.clf()
            lo,cathi = 10,18
            if band == 3:
                lo,cathi = 8, 13
            elif band == 4:
                #lo,cathi = 4.5, 10.5
                lo,cathi = 4.5, 12
            loghist(WISE.get('w%impro'%band)[I], T.get(wband+'_mag')[J],
                    range=((lo,cathi),(lo,cathi)), bins=200)
            plt.xlabel('WISE W%i mag' % band)
            plt.ylabel('Tractor W%i mag' % band)
            plt.title('WISE catalog vs Tractor forced photometry')
            plt.axis([cathi,lo,cathi,lo])
            ps.savefig()

        print 'Tile', tile.coadd_id, 'band', wband, 'took', Time()-tb0

    T.cut(inbounds)

    if 4 in bands:
        print 'Fixing W4'
        # Post-facto fix W4 error estimates due to upsampling/covariance
        T.w4_nanomaggies_ivar *= 0.25
        # mag error = np.abs((-2.5 / np.log(10.)) * dflux / flux)
        #   dflux = 1/sqrt(flux_invvar)
        T.w4_mag_err *= 2.

    # ??
    T.delete_column('psfflux')
    T.delete_column('cmodelflux')
    T.delete_column('devflux')
    T.delete_column('expflux')
    T.treated_as_pointsource = T.treated_as_pointsource.astype(np.uint8)
    T.pointsource = T.pointsource.astype(np.uint8)

    T.writeto(outfn, header=hdr)
    print 'Wrote', outfn

    if opt.savewise:
        for band in bands:
            UW.set('fit_flux_w%i' % band, fitwiseflux[band])
        UW.writeto(savewise_outfn)
        print 'Wrote', savewise_outfn

    if savepickle:
        fn = opt.output % (tile.coadd_id)
        fn = fn.replace('.fits','.pickle')
        pickle_to_file((mods, cats, T, sourcerad), fn)
        print 'Pickled', fn

    if opt.splitrcf:
        unsplitoutfn = opt.unsplitoutput % (tile.coadd_id)
        cols,dropcols = _get_output_column_names(opt.bands)
        for c in T.get_columns():
            if not c in cols:
                T.delete_column(c)
        splitrcf(tile, tiles, wcs, T, unsplitoutfn, cols, dropcols, hdr)

    print 'Tile', tile.coadd_id, 'took', Time()-tt0

def _bounce_split((tile, T, wcs, outfn, unsplitoutfn, cols, dropcols, hdr)):
    try:
        print 'Reading', outfn
        objs = fits_table(outfn, columns=cols)
        print 'Read', len(objs), 'from', outfn
        print 'Writing unsplit to', unsplitoutfn
        splitrcf(tile, T, wcs, objs, unsplitoutfn, cols, dropcols, hdr)
    except:
        import traceback
        print 'Exception processing', tile.coadd_id
        traceback.print_exc()
        #raise

def _write_output(T, fn, cols, dropcols, hdr):
    cols = ['has_wise_phot'] + [c for c in cols if not c in ['id']+dropcols]
    T.writeto(fn, columns=cols, header=hdr)

def splitrcf(tile, tiles, wcs, T, unsplitoutfn, cols, dropcols, hdr):
    # -Write out any run,camcol,field that is totally contained
    # within this coadd tile, and does not touch any other coadd
    # tile.
    from unwise_coadd import get_coadd_tile_wcs, walk_wcs_boundary

    print 'Splitting tile', tile.coadd_id
    
    # Find nearby tiles
    I,J,d = match_radec(tiles.ra, tiles.dec,
                        np.array([tile.ra]), np.array([tile.dec]), 2.5)
    neartiles = tiles[I]
    neartiles.cut(neartiles.coadd_id != tile.coadd_id)
    print len(neartiles), 'tiles nearby:', neartiles.coadd_id

    # Decribe all boundaries in Intermediate World Coords with respect
    # to this tile's WCS.
    
    rr,dd = walk_wcs_boundary(wcs)
    ok,uu,vv = wcs.radec2iwc(rr, dd)
    mybounds = np.array(zip(uu,vv))
        
    bounds = []
    for t in neartiles:
        w = get_coadd_tile_wcs(t.ra, t.dec)
        rr,dd = walk_wcs_boundary(w)
        ok,uu,vv = wcs.radec2iwc(rr, dd)
        bounds.append(np.array(zip(uu,vv)))

    RCF = np.unique(zip(T.run, T.camcol, T.field))
    print 'Unique run/camcol/field:', RCF

    wfn = os.path.join(resolvedir, 'window_flist.fits')
    # For --split: figure out which fields are completely within the tile.
    W = fits_table(wfn, columns=['node', 'incl', 'mu_start', 'mu_end',
                                 'nu_start', 'nu_end', 'ra', 'dec',
                                 'rerun', 'run', 'camcol', 'field'])
    straddle = np.zeros(len(T), bool)

    # HACK -- rerun
    rr = '301'

    for run,camcol,field in RCF:
            
        I = np.flatnonzero((W.run == run) * (W.camcol == camcol) *
                           (W.field == field) * (W.rerun == rr))
        assert(len(I) == 1)
        wi = W[I[0]]
        rd = np.array([munu_to_radec_deg(mu, nu, wi.node, wi.incl)
                       for mu,nu in [(wi.mu_start, wi.nu_start),
                                     (wi.mu_start, wi.nu_end),
                                     (wi.mu_end,   wi.nu_end),
                                     (wi.mu_end,   wi.nu_start)]])
        ok,uu,vv = wcs.radec2iwc(rd[:,0], rd[:,1])
        poly = np.array(zip(uu,vv))
        #assert(polygons_intersect(poly, mybounds))

        J = np.flatnonzero((T.run == run) * (T.camcol == camcol) *
                           (T.field == field))

        strads = False
        for b in bounds:
            if polygons_intersect(poly, b):
                print 'Field', run, camcol, field, 'straddles tiles'
                straddle[J] = True
                strads = True
                break
        if strads:
            continue

        print 'Field', run, camcol, field, 'totally within tile', tile.coadd_id

        myoutdir = os.path.join(pobjoutdir, rr, '%i'%run, '%i'%camcol)
        if not os.path.exists(myoutdir):
            try:
                os.makedirs(myoutdir)
            except:
                import traceback
                print 'Exception creating dir', myoutdir
                traceback.print_exc()
                
        outfn = os.path.join(myoutdir, 'photoWiseForced-%06i-%i-%04i.fits' %
                             (run, camcol, field))

        N = get_photoobj_length(rr, run, camcol, field)
        print 'PhotoObj has', N, 'rows'

        Ti = T[J]
        
        P = fits_table()
        P.has_wise_phot = np.zeros(N, bool)
        I = Ti.id - 1
        P.has_wise_phot[I] = True
        for col in Ti.get_columns():
            if col in dropcols:
                continue
            tval = Ti.get(col)
            X = np.zeros(N, tval.dtype)
            X[I] = tval
            P.set(col, X)
        
        _write_output(P, outfn, cols, dropcols, hdr)
        print 'Wrote', outfn

    # -Write out the remaining objects to be --finish'd later.
    S = T[straddle]
    print len(S), 'objects straddle tiles'
    if len(S):
        rcf = np.unique(zip(S.run, S.camcol, S.field))
        print 'Unique rcf of straddled sources:', rcf
        S.writeto(unsplitoutfn)
        print 'Wrote to', unsplitoutfn
    print 'Done'
    

def todo(A, opt, ps):
    need = []
    for i in range(len(A)):
        outfn = opt.output % (A.coadd_id[i])
        #outfn = opt.unsplitoutput % (A.coadd_id[i])
        print 'Looking for', outfn
        if not os.path.exists(outfn):
            need.append(i)
    print ' '.join('%i' %i for i in need)
            
    # Collapse contiguous ranges
    strings = []
    if len(need):
        start = need.pop(0)
        end = start
        while len(need):
            x = need.pop(0)
            if x == end + 1:
                # extend this run
                end = x
            else:
                # run finished; output and start new one.
                if start == end:
                    strings.append('%i' % start)
                else:
                    strings.append('%i-%i' % (start, end))
                start = end = x
        # done; output
        if start == end:
            strings.append('%i' % start)
        else:
            strings.append('%i-%i' % (start, end))
        print ','.join(strings)
    else:
        print 'Done (party now)'

        
def summary(A, opt, ps):
    plt.clf()
    missing = []
    for i in range(len(A)):
        r,d = A.ra[i], A.dec[i]
        dd = 1024 * 2.75 / 3600.
        dr = dd / np.cos(np.deg2rad(d))
        outfn = opt.output % (A.coadd_id[i])
        rr,dd = [r-dr,r-dr,r+dr,r+dr,r-dr], [d-dd,d+dd,d+dd,d-dd,d-dd]
        print 'Looking for', outfn
        if not os.path.exists(outfn):
            missing.append((i,rr,dd,r,d))
        plt.plot(rr, dd, 'k-')
    for i,rr,dd,r,d in missing:
        plt.plot(rr, dd, 'r-')
        plt.text(r, d, '%i' % i, rotation=90, color='b', va='center', ha='center')
    plt.title('missing tiles')
    plt.axis([118, 212, 44,61])
    ps.savefig()

    print 'Missing tiles:', [m[0] for m in missing]

    rdfn = 'rd.fits'
    if not os.path.exists(rdfn):
        fns = glob(os.path.join(tempoutdir, 'photoobjs-*.fits'))
        fns.sort()
        TT = []
        for fn in fns:
            T = fits_table(fn, columns=['ra','dec'])
            print len(T), 'from', fn
            TT.append(T)
        T = merge_tables(TT)
        print 'Total of', len(T)
        T.writeto(rdfn)
    else:
        T = fits_table(rdfn)
        print 'Got', len(T), 'from', rdfn
    
    plt.clf()
    loghist(T.ra, T.dec, 500, range=((118,212),(44,61)))
    plt.xlabel('RA')
    plt.ylabel('Dec')
    ps.savefig()

    ax = plt.axis()
    for i in range(len(A)):
        r,d = A.ra[i], A.dec[i]
        dd = 1024 * 2.75 / 3600.
        dr = dd / np.cos(np.deg2rad(d))
        plt.plot([r-dr,r-dr,r+dr,r+dr,r-dr], [d-dd,d+dd,d+dd,d-dd,d-dd], 'r-')
    plt.axis(ax)
    ps.savefig()

def _get_output_column_names(bands):
    cols = ['ra','dec', 'raerr', 'decerr', 'objid', 'x','y', 
            'treated_as_pointsource', 'pointsource', 'coadd_id', 'modelflux']
    for band in bands:
        for k in ['nanomaggies', 'nanomaggies_ivar', 'mag', 'mag_err',
                  'prochi2', 'pronpix', 'profracflux', 'proflux', 'npix',
                  'pronexp']:
            cols.append('w%i_%s' % (band, k))
    cols.extend(['run','camcol','field','id'])
    # columns to drop from the photoObj-parallels
    dropcols = ['run', 'camcol', 'field', 'modelflux']
    return cols, dropcols

def finish(A, opt, args, ps, hdr):
    '''
    A: atlas
    '''
    # Find all *-phot.fits outputs
    # Determine which photoObj files are involved
    # Collate and resolve objs measured in multiple tiles
    # Expand into photoObj-parallel files
    if len(args):
        fns = args
    else:
        fns = glob(os.path.join(outdir, 'phot-????????.fits'))
        fns.sort()
        print 'Found', len(fns), 'photometry output files'

    cols,dropcols = _get_output_column_names(opt.bands)

    if opt.splitrcf:
        from unwise_coadd import get_coadd_tile_wcs
        args = []
        for i in range(len(A)):
            #outfn = opt.output % A.coadd_id[i]
            #if not outfn in fns:
            #    continue
            outfn = None
            coadd_id = A.coadd_id[i]
            for fn in fns:
                if coadd_id in fn:
                    outfn = fn
                    break
            if outfn is None:
                print 'Did not find input file for', coadd_id
                continue
            unsplitoutfn = opt.unsplitoutput % A.coadd_id[i]
            if os.path.exists(unsplitoutfn):
                print 'Exists:', unsplitoutfn
                continue
            print 'File', outfn
            tile = A[i]
            wcs = get_coadd_tile_wcs(tile.ra, tile.dec)
            args.append((tile, A, wcs, outfn, unsplitoutfn, cols, dropcols, hdr))

        if opt.no_threads:
            map(_bounce_split, args)
        else:
            mp = multiproc(1)
            mp.map(_bounce_split, args)
        return

    flats = []
    fieldmap = {}
    for ifn,fn in enumerate(fns):
        print
        print 'Reading', (ifn+1), 'of', len(fns), fn
        T = fits_table(fn, columns=cols + opt.ftag)
        print 'Read', len(T), 'entries'

        if opt.flat is not None:
            flats.append(T)
            continue

        rcf = np.unique(zip(T.run, T.camcol, T.field))
        #!!!!
        #rcf = np.unique(T.run * 100000 + T.camcol * 10000 + T.field)
        #rcf = zip(rcf / 100000, (rcf % 100000) / 10000, (rcf % 10000))

        for run,camcol,field in rcf:
            if not (run,camcol,field) in fieldmap:
                fieldmap[(run,camcol,field)] = []
            Tsub = T[(T.run == run) * (T.camcol == camcol) * (T.field == field)]
            for col in dropcols:
                Tsub.delete_column(col)
            print '  ', len(Tsub), 'in', run,camcol,field, ', joining', [t.coadd_id[0] for t in fieldmap[(run,camcol,field)]]
            fieldmap[(run,camcol,field)].append(Tsub)

    # WISE coadd tile CRPIX-1 (x,y in the phot-*.fits files are 0-indexed)
    # (and x,y are based on the first-band (W1 usually) WCS)
    cx,cy = 1023.5, 1023.5

    if opt.flat is not None:
        F = merge_tables(flats)
        print 'Total of', len(F), 'measurements'
        r2 = (F.x - cx)**2 + (F.y - cy)**2
        I,J,d = match_radec(F.ra, F.dec, F.ra, F.dec, 1e-6, notself=True)
        print 'Matched', len(I), 'duplicates'
        keep = np.ones(len(F), bool)
        keep[np.where(r2[I] > r2[J], I, J)] = False
        F.cut(keep)
        print 'Cut to', len(F)
        F.delete_column('x')
        F.delete_column('y')
        F.writeto(opt.flat)
        return

    keys = fieldmap.keys()
    keys.sort()

    # HACK
    rr = '301'

    lengths = get_photoobj_length(None, None, None, None, all=True)
    print 'Got', len(lengths), 'photoObj lengths'
    
    args = []
    for i,(run,camcol,field) in enumerate(keys):
        TT = fieldmap.get((run,camcol,field))
        N = lengths.get((int(rr), int(run), int(camcol), int(field)), None)
        if N is None:
            N = get_photoobj_length(rr, run, camcol, field)
        myoutdir = os.path.join(pobjoutdir, rr, '%i'%run, '%i'%camcol)
        if not os.path.exists(myoutdir):
            os.makedirs(myoutdir)
        outfn = os.path.join(myoutdir, 'photoWiseForced-%06i-%i-%04i.fits' % (run, camcol, field))
        args.append((i, len(fieldmap), TT, N, rr, run, camcol, field, outfn,
                     cols, dropcols, hdr, cx, cy))

    if opt.no_threads:
        map(_finish_one, args)
    else:
        mp = multiproc(8)
        mp.map(_finish_one, args)
    print 'Done'

def _finish_one((i, Ntotal, TT, N, rr, run, camcol, field, outfn, cols, dropcols,
                 hdr, cx, cy)):
    print
    print (i+1), 'of', Ntotal, ': R,C,F', (run,camcol,field)
    print len(TT), 'tiles for', (run,camcol,field)

    resolve = (len(TT) > 1)

    if os.path.exists(outfn):
        print 'Output file already exists.  Updating.'
        P = fits_table(outfn)
        assert(N == len(P))
        print 'Read', N, 'from', outfn

        P.R2 = np.empty(N, np.float32)
        P.R2[:] = 1e9
        I = np.flatnonzero((P.x != 0) * (P.y != 0))
        P.R2[I] = (P.x[I] - cx)**2 + (P.y[I] - cy)**2
        resolve = True

    else:
        P = fits_table()
        P.has_wise_phot = np.zeros(N, bool)
        if resolve:
            # Resolve duplicate measurements (in multiple tiles)
            # based on || (x,y) - center ||^2
            P.R2 = np.empty(N, np.float32)
            P.R2[:] = 1e9

    for T in TT:
        coadd = T.coadd_id[0]
        if resolve:
            #print 'Coadd', coadd, ': T:'
            #T.about()
            I = T.id - 1
            R2 = (T.x - cx)**2 + (T.y - cy)**2
            J = (R2 < P.R2[I])
            I = I[J]
            P.R2[I] = R2[J].astype(np.float32)
            print '  tile', coadd, ':', len(I), 'are closest'
            T.cut(J)
            # Note here that we just choose which indices will be
            # *overwritten* by this "T" -- only those closer than any
            # existing 'T'; and those rows may be in turn overwritten
            # by a later one.
        print '  ', len(T), 'from', coadd
        if len(T) == 0:
            continue
        #print 'Coadd', coadd, ': T:'
        #T.about()
        I = T.id - 1
        P.has_wise_phot[I] = True
        pcols = P.get_columns()
        for col in T.get_columns():
            if col in pcols:
                pval = P.get(col)
                #print '  ', col, pval.dtype
                pval[I] = (T.get(col)).astype(pval.dtype)
            else:
                tval = T.get(col)
                X = np.zeros(N, tval.dtype)
                X[I] = tval
                P.set(col, X)

    if resolve:
        P.delete_column('R2')
        ##
        coadds = np.unique(P.coadd_id)
        print 'Coadds:', coadds
        for c in coadds:
            I = np.flatnonzero((P.coadd_id == c))
            print '  ', len(I), 'from', c
    _write_output(P, outfn, cols, dropcols, hdr)
    print 'Wrote', outfn

def check(T, opt, args, ps):
    print 'Walking', pobjoutdir
    decstep = 0.1
    rastep = 0.1
    ralo,rahi = 0., 360.
    declo,dechi = -30., 90.
    Nra = int((rahi - ralo) / rastep)
    Ndec = int((dechi - declo) / decstep)
    hist = np.zeros((Ndec, Nra), np.int16)

    def binimg(img, b):
        hh,ww = img.shape
        hh = int(hh / b) * b
        ww = int(ww / b) * b
        binx = reduce(np.add, [img[:, i:hh:b] for i in range(b)])
        return reduce(np.add, [img[i:ww:b, :] for i in range(b)])

    plt.figure(figsize=(12,5))
    def _plotit():
        fitsio.write('sdss-phot-density.fits', hist, clobber=True)
        plt.clf()
        plt.imshow(binimg(hist, 8), interpolation='nearest', origin='lower',
                   extent=(ralo,rahi,declo,dechi))
        plt.colorbar()
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        ps.savefig()

    n = 0
    for dirpath,dirnames,fns in os.walk(pobjoutdir, followlinks=True):
        dirnames.sort()
        fns.sort()
        for fn in fns:
            pth = os.path.join(dirpath, fn)
            print 'Reading', (n+1), pth
            T = fits_table(pth, columns=['ra','dec','has_wise_phot'])
            T.has_wise_phot = (T.has_wise_phot.astype(np.uint8) == ord('T'))
            #print 'has_wise_phot:', np.unique(T.has_wise_phot)
            #print 'has_wise_phot:', np.unique(T.has_wise_phot.astype(np.uint8))
            print '    Got', len(T), ',', sum(T.has_wise_phot), 'with photometry'
            T.cut(T.has_wise_phot)
            if len(T) == 0:
                continue
            for ira,idec in zip(((T.ra  - ralo )/ rastep).astype(int),
                                ((T.dec - declo)/decstep).astype(int)):
                hist[idec,ira] += 1
            n += 1
            if n and n % 1000 == 0:
                _plotit()
    _plotit()

def main():
    import optparse

    global outdir
    global tempoutdir
    global pobjoutdir
    global tiledir
    global photoobjdir
    global resolvedir

    import datetime
    import sys
    import os
    print 'forcedphot.py starting at', datetime.datetime.now().isoformat()
    print 'command-line args:', sys.argv
    print 'git describe:'
    os.system('git describe')

    parser = optparse.OptionParser('%prog [options]')
    parser.add_option('--minsig1', dest='minsig1', default=0.1, type=float)
    parser.add_option('--minsig2', dest='minsig2', default=0.1, type=float)
    parser.add_option('--minsig3', dest='minsig3', default=0.1, type=float)
    parser.add_option('--minsig4', dest='minsig4', default=0.1, type=float)
    parser.add_option('-d', dest='outdir', default=None,
                      help='Output directory')
    parser.add_option('--tempdir', default=None,
                      help='"Temp"-file output directory')
    parser.add_option('--pobj', dest='pobjdir', default=None,
                      help='Output directory for photoObj-parallels')
    parser.add_option('-o', dest='output', default=None, help='Output filename pattern')
    parser.add_option('-b', '--band', dest='bands', action='append', type=int, default=[],
                      help='Add WISE band (default: 1,2)')

    parser.add_option('--tiledir', type=str, help='Set input unWISE coadds dir; default %s' % tiledir)

    parser.add_option('--photoobjs-only', dest='photoObjsOnly',
                      action='store_true', default=False,
                      help='Ensure photoobjs file exists and then quit?')
    parser.add_option('--wise-only', dest='wiseOnly',
                      action='store_true', default=False,
                      help='Ensure WISE file exists and then quit?')

    parser.add_option('--photoobjsdir', help='Set photoObj input directory',
                      default='photoObjs')
    parser.add_option('--resolvedir', help='Set resolve input directory',
                      default='photoResolve')

    parser.add_option('-p', dest='pickle', default=False, action='store_true',
                      help='Save .pickle file for debugging purposes')
    parser.add_option('--plots', dest='plots', default=False, action='store_true')

    parser.add_option('--save-fits', dest='save_fits', default=False, action='store_true')

    parser.add_option('--plotbase', dest='plotbase', help='Base filename for plots')

    parser.add_option('--finish', dest='finish', default=False, action='store_true')

    parser.add_option('--check', dest='check', default=False, action='store_true')

    parser.add_option('--ftag', dest='ftag', action='append', default=[],
                      help='Tag-along extra columns in --finish phase')

    parser.add_option('--flat', dest='flat', type='str', default=None,
                      help='Just write a flat-file of (deduplicated) results, not photoObj-parallels')

    parser.add_option('--summary', dest='summary', default=False, action='store_true')
    parser.add_option('--todo', dest='todo', default=False, action='store_true')

    parser.add_option('--no-ceres', dest='ceres', action='store_false', default=True,
                       help='Use scipy lsqr rather than Ceres Solver?')

    parser.add_option('--ceres-block', '-B', dest='ceresblock', type=int, default=10,
                      help='Ceres image block size (default: %default)')

    parser.add_option('--minrad', dest='minradius', type=int, default=2,
                      help='Minimum radius, in pixels, for evaluating source models; default %default')

    parser.add_option('--sky', dest='sky', action='store_true', default=False,
                      help='Fit sky level also?')

    parser.add_option('--save-wise', dest='savewise', action='store_true', default=False,
                      help='Save WISE catalog source fits also?')
    parser.add_option('--save-wise-out', dest='save_wise_output', default=None)

    parser.add_option('--dataset', dest='dataset', default='sequels',
                      help='Dataset (region of sky) to work on')

    parser.add_option('--errfrac', dest='errfrac', type=float,
                      help='Add this fraction of flux to the error model.')

    parser.add_option('--bright1', dest='bright1', type=float, default=None,
                      help='Subtract WISE model PSF for stars brighter than this in W1')

    parser.add_option('--split', dest='splitrcf', action='store_true', default=False,
                      help='Split outputs into run/camcol/field right away?')

    parser.add_option('--tile', dest='tile', action='append', default=[],
                      help='Run a single tile')

    parser.add_option('--psf', dest='psffn', default='psf-allwise-con3.fits')

    parser.add_option('-v', dest='verbose', default=False, action='store_true')

    parser.add_option('--no-threads', action='store_true')

    parser.add_option('--l1b', action='store_true', help='Use individual exposures (L1b images), not coadds')
    parser.add_option('--l1b-sky', action='store_true', default=False,
                      help='Subtract coadd sky level from l1b frames?')
    
    opt,args = parser.parse_args()

    opt.unsplitoutput = None
    
    if opt.tiledir:
        tiledir = opt.tiledir

    photoobjdir = opt.photoobjsdir
    print 'Set photoObjs directory:', photoobjdir
    resolvedir = opt.resolvedir
    print 'Set resolve directory:', resolvedir
    

    if len(opt.bands) == 0:
        opt.bands = [1,2]

    # Allow specifying bands like "123"
    bb = []
    for band in opt.bands:
        for s in str(band):
            bb.append(int(s))
    opt.bands = bb
    print 'Bands', opt.bands

    lvl = logging.INFO
    if opt.verbose:
        lvl = logging.DEBUG
    logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

    # sequels-atlas.fits: written by wise-coadd.py
    fn = '%s-atlas.fits' % opt.dataset
    print 'Reading', fn
    T = fits_table(fn)

    version = {}
    from astrometry.util.run_command import run_command
    for key,cmd in [('Revision', 'git describe'),
                    ('URL', 'git config --get remote.origin.url'),
                    ]:
        rtn,txt,err = run_command(cmd)
        if rtn:
            raise RuntimeError('Failed to get version string (%s):' % cmd + txt + err)
        txt = txt.strip()
        version[key] = txt

    hdr = fitsio.FITSHDR()
    hdr.add_record(dict(name='WPHO_VER', value=version['Revision'],
                        comment='Git describe'))
    hdr.add_record(dict(name='WPHO_URL', value=version['URL'], comment='Git URL'))
    hdr.add_record(dict(name='WPHO_DATE', value=datetime.datetime.now().isoformat(),
                        comment='forced phot run time'))
    hdr.add_record(dict(name='WPHO_SKY', value=opt.sky, comment='fit sky?'))
    for band in opt.bands:
        minsig = getattr(opt, 'minsig%i' % band)
        hdr.add_record(dict(name='WPHO_MNS%i' % band, value=minsig,
                            comment='min surf brightness in sig, band %i' % band))
    hdr.add_record(dict(name='WPHO_CERE', value=opt.ceres, comment='use Ceres?'))
    hdr.add_record(dict(name='WPHO_ERRF', value=opt.errfrac, comment='error flux fraction'))
    if opt.ceres:
        hdr.add_record(dict(name='WPHO_CEBL', value=opt.ceresblock,
                        comment='Ceres blocksize'))
    
    if opt.plotbase is None:
        opt.plotbase = opt.dataset + '-phot'
    ps = PlotSequence(opt.plotbase)

    outdir     = outdir     % opt.dataset
    tempoutdir = tempoutdir % opt.dataset
    pobjoutdir = pobjoutdir % opt.dataset

    if opt.pobjdir is not None:
        pobjoutdir = opt.pobjdir

    if opt.outdir is not None:
        outdir = opt.outdir
    else:
        # default
        opt.outdir = outdir

    if opt.tempdir is not None:
        tempoutdir = opt.tempdir
    else:
        # default
        opt.tempdir = tempoutdir
        
    if opt.output is None:
        opt.output = os.path.join(outdir, 'phot-%s.fits')
    if opt.unsplitoutput is None:
        opt.unsplitoutput = os.path.join(outdir, 'phot-unsplit-%s.fits')
    if opt.save_wise_output is None:
        opt.save_wise_output = opt.output.replace('phot-', 'phot-wise-')

    if opt.summary:
        summary(T, opt, ps)
        sys.exit(0)

    if opt.todo:
        todo(T, opt, ps)
        sys.exit(0)
        
    if opt.finish:
        finish(T, opt, args, ps, hdr)
        sys.exit(0)

    if opt.check:
        check(T, opt, args, ps)
        sys.exit(0)
        
    for dirnm in [outdir, tempoutdir, pobjoutdir]:
        if not os.path.exists(dirnm):
            try:
                os.makedirs(dirnm)
            except:
                pass

    disable_galaxy_cache()

    tiles = []
    arr = os.environ.get('PBS_ARRAYID')
    if arr is not None:
        arr = int(arr)
        tiles.append(arr)

    if len(opt.tile):
        for t in opt.tile:
            I = np.flatnonzero(T.coadd_id == t)
            if len(I) == 0:
                print 'Failed to find tile id', t, 'in dataset', opt.dataset
                return -1
            assert(len(I) == 1)
            tiles.append(I[0])

    for a in args:
        if '-' in a:
            aa = a.split('-')
            if len(aa) != 2:
                print 'With arg containing a dash, expect two parts'
                print aa
                sys.exit(-1)
            start = int(aa[0])
            end = int(aa[1])
            for i in range(start, end+1):
                tiles.append(i)
        else:
            tiles.append(int(a))

    if len(tiles) == 0:
        tiles.append(0)

    for i in tiles:
        if opt.plots:
            plot = ps
        else:
            plot = None
        print
        print 'Tile index', i, 'coadd', T.coadd_id[i]

        if opt.wiseOnly:
            # Find WISE-only catalog sources
            tile = T[i]
            thisdir = get_unwise_tile_dir(tiledir, tile.coadd_id)
            bands = opt.bands
            fn = os.path.join(thisdir, 'unwise-%s-w%i-img-m.fits' %
                              (tile.coadd_id, bands[0]))
            print 'Reading', fn
            wcs = Tan(fn)
            wfn = os.path.join(tempoutdir, 'wise-sources-%s.fits' % (tile.coadd_id))
            WISE = read_wise_sources(wfn, wcs)
            continue
        
        one_tile(T[i], opt, opt.pickle, plot, T, tiledir, tempoutdir, hdr=hdr)

if __name__ == '__main__':
    main()

