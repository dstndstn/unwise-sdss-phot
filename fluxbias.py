import sys
import fitsio
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

from tractor import *
from astrometry.util.plotutils import *
from astrometry.util.fits import *

def expected_flux(F, s):
    lam = F / 2.

    mn = lam
    st = np.sqrt(lam)
    lo = int(np.floor(max(0, mn - 10.*st)))
    hi = int(np.ceil (       mn + 10.*st))
    print 'lam', lam, 'lo,hi', lo,hi

    f1,f2 = np.meshgrid(np.arange(lo, hi+1), np.arange(lo, hi+1))
    #print 'f1', f1.dtype

    poiss = []
    p = 1.
    p *= np.exp(-lam)
    for f in range(1, lo):
        p *= lam / f
    for f in range(lo, hi+1):
        if f > 0:
            p *= lam / f
        poiss.append(p)
    poiss = np.array(poiss)
    #print 'Poisson factors:', poiss
    #np.exp(-2.*lam) * 
    fexp = 2. * np.sum(
        (2.*f1*f2 + s*(f1+f2)) / (f1+f2+2.*s) *
        poiss[:,np.newaxis] * poiss[np.newaxis,:])
    print '<f>', fexp
    return fexp
    
def analytic(ps):
    for sky in [10, 30, 100, 300, 1000]:
        fexp = []
        ftrue = []
        for F in np.logspace(0., 3., 15):
            print 'Flux', F
            f = expected_flux(F, sky)
            fexp.append(f)
            ftrue.append(F)
        fexp = np.array(fexp)
        ftrue = np.array(ftrue)
        plt.clf()
        plt.plot(ftrue, fexp / ftrue, 'b.')
        plt.xlabel('True flux')
        plt.ylabel('Measured flux / True flux')
        plt.title('Sky=%g' % sky)
        plt.xscale('log')
        plt.xlim(min(ftrue)*0.9, max(ftrue)*1.1)
        ps.savefig()

    
if __name__ == '__main__':
    ps = PlotSequence('fluxbias')

    analytic(ps)
    sys.exit(0)
    
    
    #psf = NCircularGaussianPSF([4.], [1.])
    psf = NCircularGaussianPSF([0.88], [1.])
    #W,H = 50,50
    W,H = 100,100
    cx,cy = W/2., H/2.
    sig1 = 1.
    img = np.zeros((H,W), np.float32)
    tim = Image(data=img, inverr=np.ones_like(img)/sig1, psf=psf,
                photocal=LinearPhotoCal(1.))

    #src = PointSource(PixPos(cx, cy), Flux(1.))
    # right between pixels in x,y
    src = PointSource(PixPos(int(cx)+0.5, int(cy)+0.5), Flux(1.))
    #src = PointSource(PixPos(int(cx)+0.5, int(cy)), Flux(1.))

    tractor = Tractor([tim], [src])

    trueimg = tractor.getModelImage(0)

    plt.clf()
    dimshow(trueimg)
    plt.colorbar()
    ps.savefig()

    # For a range of fluxes,
    # Draw multiple Poisson realizations
    # Add Gaussian "read noise"
    # Estimate pre-pixel error ~ flux
    # Forced-photometer
    
    #tractor.printThawedParams()
    
    tractor.freezeAllRecursive()
    tractor.thawPathsTo('brightness')

    from tractor.ceres_optimizer import CeresOptimizer
    B = 10
    tractor.optimizer = CeresOptimizer(BW=B, BH=B)

    print 'Thawed params:'
    tractor.printThawedParams()

    p0 = tractor.getParams()
    print 'P0:', p0

    # WISE W1: gain 3.1 e-/DN, sky 14 DN/pix, 24 exposures
    truesky = 3.1 * 14. * 24
    #truesky = 30.

    readnoise = 3.1

    #for measure_sky in [True, False]:
    for measure_sky in [False]:

        if measure_sky:
            tt = 'with measured sky'
        else:
            tt = 'with true sky'
            
        trueflux = []
        measflux = []
        fluxerr = []
        skyest = []
        
        #for flux in [10., 30., 100., 300., 1000., 3000.]:
        for flux in np.logspace(1., 5., 10):
        #for flux in np.logspace(1., 4., 13):

            for isample in range(100):

                # pimg = np.random.poisson(flux * trueimg)
                # samplesky = np.random.poisson(truesky, size=pimg.shape)
                # pimg += samplesky
                # # Gain?
                # img = pimg.astype(np.float32)
                # img += readnoise * np.random.normal(size=pimg.shape)

                #img = flux * trueimg + truesky
                pimg = np.random.poisson(flux * trueimg)
                #print 'pimg flux: sum', np.sum(pimg), 'vals', np.unique(pimg)
                img = pimg.astype(np.float32) + truesky
                print 'Image range:', img.min(), img.max()
                
                # Poisson process: variance = mean
                iv = np.maximum(0, 1./np.maximum(img, 1.))
                #print 'Invvar range:', iv.min(), iv.max()
                
                if measure_sky:
                    # Typical dumb sky estimate
                    sky = np.median(img)
                    print 'Estimated sky:', sky
                else:
                    sky = truesky
                img -= sky

                #print 'Sky:', truesky, 'vs sampled', np.mean(samplesky.astype(np.float32))
                
                tractor.setParams(p0)
    
                tim.setImage(img)
                tim.setInvvar(iv)
                
                R = tractor.optimize_forced_photometry(
                    variance=True, shared_params=False)
    
                mod = tractor.getModelImage(0)
    
                trueflux.append(flux)
                measflux.append(tractor.getParams()[0])
                skyest.append(sky)
                fluxerr.append(1./np.sqrt(R.IV[0]))
                
                if isample == 0 and False: # and measure_sky:
                    plt.clf()
                    plt.subplot(2,2,1)
                    #dimshow(np.log10(np.maximum(1e-3, flux * trueimg)))
                    tru = flux * trueimg
                    #mn,mx = np.percentile(tru, [25,99])
                    #mx = max(mx, truesky)
                    #ima = dict(vmin=mn, vmax=mx)
                    ima = dict(vmin=-100, vmax=1000)
                    dimshow(tru, **ima)
                    plt.title('true image')
                    plt.subplot(2,2,2)
                    #dimshow(np.log10(np.maximum(1e-3, tim.getImage())))
                    dimshow(tim.getImage(), **ima)
                    plt.colorbar()
                    plt.title('noisy image')
                    plt.subplot(2,2,3)
                    #dimshow(np.log10(np.clip(1./tim.getInvError(), 1e-3, 1e3)))
                    dimshow(1./tim.getInvError())
                    plt.colorbar()
                    plt.title('image sigma')
                    plt.subplot(2,2,4)
                    chi = (img - mod) * tim.getInvError()
                    dimshow(chi, cmap='RdBu', vmin=-5, vmax=5)
                    plt.title('chi')
                    plt.suptitle('true flux = %g' % flux)
                    ps.savefig()
    
        trueflux = np.array(trueflux)
        measflux = np.array(measflux)
        fluxerr  = np.array(fluxerr)
        
        plt.clf()
        plt.plot(trueflux, measflux, 'b.')
        plt.xlabel('True flux')
        plt.ylabel('Measured flux')
        plt.yscale('symlog')
        plt.xscale('log')
        ax = plt.axis()
        lo,hi = plt.xlim()
        plt.plot([lo,hi],[lo,hi], 'k-', alpha=0.25)
        plt.title(tt)
        ps.savefig()
    
        plt.clf()
        plt.semilogx(trueflux, (measflux - trueflux) / trueflux,
                     'b.', alpha=0.01)
        plt.xlabel('True flux')
        plt.ylabel('Measured flux fractional error')
        # average over n samples
        u,I = np.unique(trueflux, return_inverse=True)
        for i in range(len(u)):
            mf = measflux[i == I]
            y = (mf - u[i]) / u[i]
            #plt.plot(u[i], (np.mean(mf) - u[i]) / u[i], 'bx')
            plt.errorbar(u[i], np.mean(y), yerr=np.std(y), ecolor='b', fmt='bx')
        plt.xlim(min(u)*0.9, max(u)*1.1)
        #plt.ylim(-0.5, 0.5)
        plt.ylim(-0.05, 0.05)
        plt.axhline(0, color='k', alpha=0.25)
        plt.title(tt)
        ps.savefig()
    
        plt.clf()
        plt.semilogx(trueflux, (measflux - trueflux) / fluxerr, 'b.')
        plt.xlabel('True flux')
        plt.ylabel('Measured flux bias (sigma)')
        #plt.ylim(-0.5, 0.5)
        plt.axhline(0, color='k', alpha=0.25)
        plt.title(tt)
        ps.savefig()
        
        if measure_sky:
            plt.clf()
            plt.semilogx(trueflux, skyest, 'b.')
            plt.xlabel('True flux')
            plt.ylabel('Estimated sky')
            ps.savefig()
        
