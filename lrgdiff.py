import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
import fitsio
from astrometry.libkd.spherematch import *
from astrometry.util.util import *
from astrometry.util.fits import *
from astrometry.util.plotutils import *
from astrometry.util.miscutils import *
from astrometry.util.resample import *

from astrometry.sdss import *

from tractor import *
from tractor.sdss import *

ps = PlotSequence('lrg')


cheat = True

if False:
    T=fits_table('tractor_lrgs.fits')
    W=fits_table('allwise_lrgs.fits')
    plt.clf()
    plt.plot(T.ra, T.dec, 'rx')
    plt.plot(W.ra, W.dec, 'b.')
    plt.savefig('/tmp/1.png')
    I,J,d = match_radec(T.ra, T.dec, W.ra, W.dec, 4./3600.)
    matched = np.zeros(len(T), bool)
    matched[I] = True
    k=np.flatnonzero(matched == False)


T = fits_table('bad_spec_data.fits') # from Abhi Mar 14/2014
print 'Got', len(T), 'bad spectra'

T.forced_w1mag = np.zeros(len(T))
T.forced_w1magerr = np.zeros(len(T))
T.forced_w2mag = np.zeros(len(T))
T.forced_w2magerr = np.zeros(len(T))

A = fits_table('sdss2-atlas.fits')
I,J,d = match_radec(T.ra, T.dec, A.ra, A.dec, 1.6, nearest=True)
#matched = np.zeros(len(A), bool)
#matched[J] = True
#A.cut(matched)
A.cut(J)
print 'Cut to', len(A), 'unWISE tiles'

F = fits_table('window_flist.fits')
print 'Read', len(F), 'fields'
margin = 0.2
F.cut((F.ra  < (T.ra.max () + margin)) *
      (F.ra  > (T.ra.min () - margin)) *
      (F.dec < (T.dec.max() + margin)) *
      (F.dec > (T.dec.min() - margin)))
print 'Cut to', len(F), 'within range'

F.rdcorners = np.empty((len(F), 4, 2))
F.rdcorners[:,0,:] = np.vstack(munu_to_radec_deg(F.mu_start, F.nu_start, F.node, F.incl)).T
F.rdcorners[:,1,:] = np.vstack(munu_to_radec_deg(F.mu_end  , F.nu_start, F.node, F.incl)).T
F.rdcorners[:,2,:] = np.vstack(munu_to_radec_deg(F.mu_end  , F.nu_end  , F.node, F.incl)).T
F.rdcorners[:,3,:] = np.vstack(munu_to_radec_deg(F.mu_start, F.nu_end  , F.node, F.incl)).T
#print 'rdcorners', F.rdcorners.shape

sdss = DR9()

plt.figure(figsize=(2,2))
plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

# fn = 'unwise-1497p015-w1-img-u.fits'
# I = fitsio.read(fn)
# H,W = I.shape
# wcs = Tan(fn)

print >>sys.stderr, '<table><tr>' + ''.join(['<th>%s</th>' % x for x in ['index','RA,Dec','SDSS','W1','W2','u','g','r','i','z']]) + '</tr>'
for kk in range(len(T)):

    tile = A.coadd_id[kk]
    r,d = T.ra[kk], T.dec[kk]

    Wimgs = []
    skip = False
    for band in [1,2]:
        pth = os.path.join(tile[:3], tile,
                           'unwise-%s-w%i-img-u.fits' % (tile,band))
        fn = os.path.join('data','unwise','unwise-comp', pth)
        print 'WISE file', fn
        if not os.path.exists(fn):
            cmd = 'rsync -Rvz carver:unwise-comp/./%s data/unwise/unwise-comp 2>&1' % pth
            print cmd
            os.system(cmd)
        # Wimg = fitsio.read(fn)
        # Wimgs.append(Wimg)
        # H,W = Wimg.shape
        wcs = Tan(fn)

        if not wcs.is_inside(r,d):
            skip = True
            break

        ok,x,y = wcs.radec2pixelxy(r,d)
        x -= 1.
        y -= 1.
        S = 7
        x0 = max(int(x-S),0)
        y0 = max(int(y-S),0)

        # isub = Wimg[y0:,x0:]
        # hs,ws = isub.shape
        # if hs < 2*S+1 or ws < 2*S+1:
        #     continue
        # isub = isub[:2*S+1, :2*S+1]

        Wimg = fitsio.FITS(fn)[0][y0:y0+2*S+1, x0:x0+2*S+1]
        Wimgs.append(Wimg)
        
    if skip:
        continue
        
    wfn = 'wise-sources-%s.fits' % tile
    if not os.path.exists(wfn):
        print 'rsync -arvz sdss2-phot-temp/%s broiler:/tmp' % wfn
        continue
    W = fits_table(wfn)
    W.cut((np.abs(W.ra - r) < 0.01) *
          (np.abs(W.dec - d) < 0.01))
    print 'cut to', len(W), 'WISE sources nearby'
    
    fns = []

    for Wimg in Wimgs:
        circargs = dict(ms=20, mec='r', mfc='none')

        plt.clf()
        plt.imshow(Wimg, interpolation='nearest',
                   origin='lower', cmap='gray', vmin=-20, vmax=100)
        #plt.colorbar()
        ax = plt.axis()
        plt.plot(x-x0, y-y0, 'o', **circargs)
        # WISE catalog objects
        if len(W):
            ok,wx,wy = wcs.radec2pixelxy(W.ra, W.dec)
            plt.plot(wx-x0-1, wy-y0-1, 'x', ms=20, mec='g', mfc='none')
        plt.axis(ax)
        fn = ps.getnext()
        fns.append(fn)
        plt.savefig(fn)

    #http://data.sdss3.org/sas/ebosswork/eboss/wise-sdss-forced-phot/1000/1/
    # WISE forced-photometry catalog.
    run,camcol,field,oid = T.run[kk], T.camcol[kk], T.field[kk], T.id[kk]
    fn = 'photoWiseForced-%06i-%i-%04i.fits' % (run, camcol, field)
    if not os.path.exists(fn):
        print 'File missing:', fn
        uname = os.environ['SDSS3_DATA_USER']
        passwd = os.environ['SDSS3_DATA_PASS']
        cmd = ('wget --user %s --password %s http://data.sdss3.org/sas/ebosswork/eboss/wise-sdss-forced-phot/%i/%i/%s 2>&1' % 
               (uname, passwd, run, camcol, fn))
        print cmd
        os.system(cmd)

    FF = fits_table(fn)
    #FF.oid = np.bitwise_and(FF.objid.astype(np.int64), 0xffff)
    #FF.oid = np.array([o & 0xffff for o in FF.objid])
    #print 'forced-phot obj ids:', FF.oid.min(), FF.oid.max()
    #print 'looking for id', oid
    #I = np.flatnonzero(FF.oid == oid)
    #print 'Found', len(I), 'forced-phot matches'
    #assert(len(I) == 1)
    #ff = FF[I[0]]
    ff = FF[oid-1]
    T.forced_w1mag[kk] = ff.w1_mag
    T.forced_w2mag[kk] = ff.w2_mag
    T.forced_w1magerr[kk] = ff.w1_mag_err
    T.forced_w2magerr[kk] = ff.w2_mag_err
        

    Fk = F[(F.ra  < (r + margin)) *
           (F.ra  > (r - margin)) *
           (F.dec < (d + margin)) *
           (F.dec > (d - margin))]
        
    ii = []
    for i in range(len(Fk)):
        if not point_in_poly(r, d, Fk.rdcorners[i,:,:])[0]:
            continue
        ii.append(i)
    iii = np.argmax(Fk.score[np.array(ii)])
    f = Fk[ii[iii]]
    #print 'Best-scored field containing:', f.run, f.camcol, f.field

    for b in 'ugriz':

        if cheat:
            plotfn = ps.getnext()
            fns.append(plotfn)
            continue
        
        fn = sdss.retrieve('frame', f.run, f.camcol, f.field, band=b)
        frame = sdss.readFrame(f.run, f.camcol, f.field, b, filename=fn)
        img = frame.getImage()
        astrans = frame.getAsTrans()

        H,W = img.shape
        sdsswcs = AsTransWrapper(astrans, W, H)

        sz = 200
        pixscale = 0.2 / 3600.
        targetwcs = Tan(r, d, sz/2. + 0.5, sz/2. + 0.5,
                        -pixscale, 0., 0., pixscale, sz, sz)

        targetim = np.zeros((sz,sz))
        try:
            (Yo,Xo, Yi,Xi, rims) = resample_with_wcs(targetwcs, sdsswcs, [], 3)
        except OverlapError:
            continue
        except:
            import traceback
            traceback.print_exc()
            raise
        targetim[Yo,Xo] = img[Yi,Xi]

        plt.clf()
        plt.imshow(targetim, vmin=-0.2, vmax=1., cmap='gray',
                   interpolation='nearest', origin='lower')
        ax = plt.axis()
        plt.plot(sz/2.-0.5, sz/2.-0.5, 'o', **circargs)
        plt.axis(ax)
        plt.savefig(plotfn)
        

    print >>sys.stderr, '<tr><td>%i</td><td>%g,%g</td><td><img src="http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra=%f&dec=%f&WIDTH=200&HEIGHT=200&scale=0.2"></td>' % (kk, r, d, r, d)
    for fn in fns:
        print >>sys.stderr, '<td><img src="%s"></td>' % fn
    print >>sys.stderr, '</tr>'
    print >>sys.stderr, '<tr><td></td><td></td><td></td>'
    #<td></td><td></td>'
    print >>sys.stderr, ''.join(['<td align="center">%.2f &plusmn; %.2f</td>' % (m,e) for m,e in
                                 [(T.forced_w1mag[kk], T.forced_w1magerr[kk]),
                                  (T.forced_w2mag[kk], T.forced_w2magerr[kk])]])
    print >>sys.stderr, ''.join(['<td align="center">%g</td>' % m for m in T.mag[kk,:]]) + '</tr>'
    
print >>sys.stderr, '</table>'
