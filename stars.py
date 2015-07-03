import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', serif='computer modern roman')
matplotlib.rc('font', **{'sans-serif': 'computer modern sans serif'})
import numpy as np
import pylab as plt
import os
import sys

from astrometry.util.file import *
from astrometry.util.fits import *
from astrometry.util.multiproc import *
from astrometry.util.plotutils import *
from astrometry.util.miscutils import *
from astrometry.util.util import *
from astrometry.util.resample import *
from astrometry.libkd.spherematch import *
from astrometry.util.starutil_numpy import *

import fitsio

from forcedphot import read_wise_sources

def gals():
    '''
    select s.class as clazz, s.subclass as subclazz, s.ra, s.dec, s.survey, s.programname,
    s.run, s.camcol, s.field, s.obj,
    s.psfMag_u, s.psfMag_g, s.psfMag_r, s.psfMag_i, s.psfMag_z,
    s.psfMagErr_u, s.psfMagErr_g, s.psfMagErr_r, s.psfMagErr_i, s.psfMagErr_z,
    s.modelmag_u, s.modelmag_g, s.modelmag_r, s.modelmag_i, s.modelmag_z,
    s.modelmagerr_u, s.modelmagerr_g, s.modelmagerr_r, s.modelmagerr_i, s.modelmagerr_z,
    s.z, s.zerr,
    p.fracdev_r, p.exprad_r, p.devrad_r
    into mydb.MyTable_6 from SpecPhoto as s join PhotoPrimary as p on s.objid = p.objid
    where s.class = 'GALAXY' and s.zwarning = 0

    modhead gals.fits+1 TDIM1 '(6)'
    modhead gals.fits+1 TDIM2 '(21)'
    modhead gals.fits+1 TDIM5 '(6)'
    modhead gals.fits+1 TDIM6 '(23)'
    '''
    G = fits_table('gals.fits')
    ps = PlotSequence('gals')
    print len(G), 'gals'

    I = np.flatnonzero((G.ra > 172) * (G.ra < 188) * (G.dec > 40) * (G.dec < 44))
    G[I].writeto('gals2.fits')

    
    plt.clf()
    plt.hist(G.z, 50)
    ps.savefig()

    plt.clf()
    ha = dict(bins=50, histtype='step', range=(0,10))
    I = np.flatnonzero(G.fracdev_r > 0.5)
    plt.hist(G.devrad_r[I], color='b', **ha)
    I = np.flatnonzero(G.fracdev_r <= 0.5)
    plt.hist(G.exprad_r[I], color='r', **ha)
    plt.xlabel('radius (arcsec)')
    ps.savefig()
    
    G.cut((G.z > 0.2) * (G.z < 0.5))
    print len(G), 'in redshift range'

    plt.clf()
    I = np.flatnonzero(G.fracdev_r > 0.5)
    plt.hist(G.devrad_r[I], color='b', **ha)
    I = np.flatnonzero(G.fracdev_r <= 0.5)
    plt.hist(G.exprad_r[I], color='r', **ha)
    plt.xlabel('radius (arcsec)')
    plt.title('Redshift [0.2,0.5]')
    ps.savefig()
    
    zrange = (0.2,0.5)

    plt.clf()
    loghist(G.z, G.modelmag_g - G.modelmag_r, 100, range=(zrange,(0.5,2.5)),
            imshowargs=dict(cmap=antigray), hot=False)
    plt.xlabel('Redshift z')
    plt.ylabel('g - r')
    plt.title('SDSS galaxies')
    ps.savefig()

    plt.clf()
    loghist(G.z, G.modelmag_r - G.modelmag_i, 100, range=(zrange,(-0.5,1.5)),
            imshowargs=dict(cmap=antigray), hot=False)
    plt.xlabel('Redshift z')
    plt.ylabel('r - i')
    plt.title('SDSS galaxies')
    ps.savefig()

    plt.clf()
    loghist(G.z, G.modelmag_i - G.modelmag_z, 100, range=(zrange,(-0.5,1.5)),
            imshowargs=dict(cmap=antigray), hot=False)
    plt.xlabel('Redshift z')
    plt.ylabel('i - z')
    plt.title('SDSS galaxies')
    ps.savefig()

    Gx = G[(G.ra > 172.5) * (G.ra < 187.5) * (G.dec > 40) * (G.dec < 50)]

    print 'Subclasses', np.unique(Gx.subclazz)
    print 'Surveys', np.unique(Gx.survey)
    print 'Programs', np.unique(Gx.programname)

    Gx.writeto('zgals.fits')
    

def stars():
    '''
    select class as clazz,subclass as subclazz,ra,dec, survey,programname,
    run,camcol,field,obj, psfMag_u, psfMag_g, psfMag_r, psfMag_i, psfMag_z, 
    psfMagErr_u,psfMagErr_g,psfMagErr_r,psfMagErr_i,psfMagErr_z
    into mydb.MyTable_4 from SpecPhoto
    where class = 'STAR' and zwarning = 0

    # CasJobs produces WRONG TDIM cards... fitsio chokes...
    modhead stars.fits+1 TDIM1 '(4)'
    modhead stars.fits+1 TDIM2 '(19)'
    modhead stars.fits+1 TDIM5 '(6)'
    modhead stars.fits+1 TDIM6 '(23)'
    '''
    S = fits_table('stars.fits')
    ps = PlotSequence('stars')
    print len(S), 'stars'
    
    S.cut((S.psfmag_g > -9999) * (S.psfmag_r > -9999) * (S.psfmag_i > -9999))
    print len(S), 'with good mags'
    
    ha = dict(bins=100, range=((-1,4),(-1,4)))
    plt.clf()
    loghist(S.psfmag_g - S.psfmag_r, S.psfmag_r - S.psfmag_i, **ha)
    ps.savefig()

    print 'Subclass:', S.subclazz.dtype, S.subclazz.shape

    S.subclazz = np.array([s.strip() for s in S.subclazz])
    print 'Subclass:', S.subclazz.dtype, S.subclazz.shape

    #K = S[S.subclazz[:,0] == 'K']
    print 'subclass:', S.subclazz[0]
    #isk = [s.strip().startswith('K') for s in S.subclazz]
    #print 'isk', isk
    K = S[np.array([s[0] == 'K' for s in S.subclazz])]
    print len(K), 'K stars'
    plt.clf()
    loghist(K.psfmag_g - K.psfmag_r, K.psfmag_r - K.psfmag_i, **ha)
    ps.savefig()
    K3 = K[np.array(['III' in s for s in K.subclazz])]
    print len(K3), 'K3 stars'
    plt.clf()
    loghist(K3.psfmag_g - K3.psfmag_r, K3.psfmag_r - K3.psfmag_i, **ha)
    ps.savefig()

    M = S[np.array([s[0] == 'M' for s in S.subclazz])]
    print len(M), 'M stars'
    print 'Subclasses:', np.unique(M.subclazz)
    plt.clf()
    loghist(M.psfmag_g - M.psfmag_r, M.psfmag_r - M.psfmag_i, **ha)
    ps.savefig()

    M3 = M[np.array(['III' in s for s in M.subclazz])]
    print len(M3), 'M3 stars'
    plt.clf()
    loghist(M3.psfmag_g - M3.psfmag_r, M3.psfmag_r - M3.psfmag_i, **ha)
    ps.savefig()

    plt.clf()
    plt.plot(K.ra, K.dec, 'b.')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.title('K')
    ps.savefig()

    Kx = K[(K.ra > 172.5) * (K.ra < 187.5) * (K.dec > 40) * (K.dec < 50)]
    print len(Kx), 'in RA,Dec region'
    Kx.writeto('kstars.fits')

    Mx = M[(M.ra > 172.5) * (M.ra < 187.5) * (M.dec > 40) * (M.dec < 50)]
    print len(Mx), 'in RA,Dec region'
    Mx.writeto('mstars.fits')
    
    
    # M5 = M[np.array(['V' in s for s in M.subclazz])]
    # print len(M5), 'M5 stars'
    # plt.clf()
    # loghist(M5.psfmag_g - M5.psfmag_r, M5.psfmag_r - M5.psfmag_i, **ha)
    # ps.savefig()
    
    



#stars()
gals()

