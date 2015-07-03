from __future__ import print_function

import os

import numpy as np

from astrometry.util.fits import *

W = fits_table('window_flist-cut.fits', columns=['run','camcol','field'])
print('Read', len(W), 'fields')
#W.cut(W.rerun == '301')
#print('Cut to', len(W), 'rerun 301')

#istart = 0
#k = 0
# Resume at i=235,001 after writing k=14 nic-014.fits
#k=15
#istart = 235001
#k=23
#istart = 323001
#k=33
#istart = 415455
k=45
istart = 493799

ilast = 0

TT = []
for i,(run,camcol,field) in enumerate(zip(W.run, W.camcol, W.field)):

    if i < istart:
        continue

    print()
    print(i+1, 'of', len(W), ': Run,camcol,field', run,camcol,field)
    fn = os.path.join(os.environ['BOSS_PHOTOOBJ'], '301', '%i'%run,
                      '%i'%camcol, 'photoObj-%06i-%i-%04i.fits' % (run,camcol,field))
    #print('Reading', fn)
    T = fits_table(fn, columns=['psfflux', 'psfflux_ivar', 'resolve_status'])
    #'run','camcol','field','objid', 'id'])
    if T is None:
        continue
    T.index = np.arange(len(T))
    #print('Read', len(T))
    T.cut((T.resolve_status & 256) > 0)
    #print(len(T), 'PRIMARY')
    i_flux = T.psfflux[:,3]
    z_flux = T.psfflux[:,4]
    i_ivar = T.psfflux_ivar[:,3]
    z_ivar = T.psfflux_ivar[:,4]
    I = np.flatnonzero((z_flux / i_flux > 2.5**2) *
                       (i_flux > 0) *
                       (i_flux**2 * i_ivar > 25.) *
                       (z_flux**2 * z_ivar > 25.))
    print(len(I), 'pass cut')
    if len(I) == 0:
        continue

    # Re-read SDSS photoObjs and grab ALL columns.
    inds = T.index[I]
    T = fits_table(fn, rows=inds)
    #T.cut(I)

    unwdir = '/project/projectdirs/cosmo/data/unwise/unwise-phot/sdss-collab/sdss-dr13-pobj'
    fn = os.path.join(unwdir, '%i'%run, '%i'%camcol, 'photoWiseForced-%06i-%i-%04i.fits' % (run, camcol, field))
    #print('Reading', fn)
    U = fits_table(fn, rows=inds)

    # *almost* all the time, has_wise_phot=True; exception is
    # sdss-dr13-pobj/4135/4/photoWiseForced-004135-4-0169.fits
    # but I haven't investigated.
    T.cut(U.has_wise_phot)
    U.cut(U.has_wise_phot)

    assert(np.all(T.objid == U.objid))

    U.delete_column('ra')
    U.delete_column('dec')
    U.delete_column('raerr')
    U.delete_column('decerr')
    U.delete_column('objid')
    U.delete_column('has_wise_phot')
    
    T.add_columns_from(U)

    # if not np.all(T.objid == U.objid):
    #     print('T objid:', T.objid)
    #     print('U objid:', U.objid)
    #     print('Has wise phot:', U.has_wise_phot)
    #assert(np.all(T.has_wise_phot))

    TT.append(T.copy())
    print('Total of', sum([len(x) for x in TT]), 'sources')
    if i - ilast >= 1000:
        ilast = i
        T = merge_tables(TT)
        print('Total of', len(T), 'pass cut')
        N = 10000
        if len(T) > N:
            fn = 'nic-%03i.fits' % k
            k += 1
            T.writeto(fn)
            print('Wrote', len(T), 'to', fn)
            TT = []
        else:
            T.writeto('nic-x.fits')
            TT = [T]

T = merge_tables(TT)

fn = 'nic-%03i.fits' % k
k += 1
T.writeto(fn)
print('Wrote', len(T), 'to', fn)

#T.writeto('nic.fits')


