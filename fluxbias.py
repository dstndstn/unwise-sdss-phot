import fitsio
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

from tractor import *
from astrometry.util.plotutils import *
from astrometry.util.fits import *

if __name__ == '__main__':
    ps = PlotSequence('fluxbias')

    psf = NCircularGaussianPSF([4.], [1.])
    #W,H = 50,50
    W,H = 100,100
    cx,cy = W/2., H/2.
    sig1 = 1.
    img = np.zeros((H,W), np.float32)
    tim = Image(data=img, inverr=np.ones_like(img)/sig1, psf=psf,
                photocal=LinearPhotoCal(1.))

    src = PointSource(PixPos(cx, cy), Flux(1.))

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
    
    readnoise = 1.

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

    truesky = 10.

    for measure_sky in [True, False]:

        if measure_sky:
            tt = 'with measured sky'
        else:
            tt = 'with true sky'
            
        trueflux = []
        measflux = []
        fluxerr = []
        skyest = []
        
        #for flux in [10., 30., 100., 300., 1000., 3000.]:
        for flux in np.logspace(2., 6., 11):
            for isample in range(10):
                pimg = np.random.poisson(flux * trueimg)
                pimg += np.random.poisson(truesky, size=pimg.shape)
                # Gain?
                img = pimg.astype(np.float32)
                img += readnoise * np.random.normal(size=pimg.shape)
    
                # Poisson process: variance = mean
                iv = np.maximum(0, 1./np.maximum(img, 1.))
                print 'Invvar range:', iv.min(), iv.max()
                
                if measure_sky:
                    # Typical dumb sky estimate
                    sky = np.median(img)
                    print 'Estimated sky:', sky
                else:
                    sky = truesky
                img -= sky
                
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
                
                if isample == 0 and measure_sky:
                    plt.clf()
                    plt.subplot(2,2,1)
                    #dimshow(np.log10(np.maximum(1e-3, flux * trueimg)))
                    tru = flux * trueimg
                    mn,mx = np.percentile(tru, [25,99])
                    mx = max(mx, truesky)
                    ima = dict(vmin=mn, vmax=mx)
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
        plt.semilogx(trueflux, (measflux - trueflux) / trueflux, 'b.')
        plt.xlabel('True flux')
        plt.ylabel('Measured flux fractional error')
        plt.ylim(-0.5, 0.5)
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
        
