from __future__ import print_function
import os
import numpy as np
from astrometry.util.fits import *
from astrometry.sdss.common import cas_flags, photo_flags1_map, photo_flags2_map

tag = 'nunez'

W = fits_table('window_flist-cut.fits', columns=['run','camcol','field'])
print('Read', len(W), 'fields')

istart = 0
k = 0
ilast = 0

TT = []
for i,(run,camcol,field) in enumerate(zip(W.run, W.camcol, W.field)):
    if i < istart:
        continue
    print()
    print(i+1, 'of', len(W), ': Run,camcol,field', run,camcol,field)
    fn = os.path.join(os.environ['BOSS_PHOTOOBJ'], '301', '%i'%run,
                      '%i'%camcol, 'photoObj-%06i-%i-%04i.fits' % (run,camcol,field))

    T = fits_table(fn, columns=['objc_type', 'flags', 'flags2',
                                'psfmagerr', 'resolve_status',])
    if T is None:
        continue
    T.index = np.arange(len(T))
    T.flags_r = T.flags[:,2]
    T.flags_i = T.flags[:,3]
    T.flags2_r = T.flags2[:,2]
    T.flags2_i = T.flags2[:,3]
    T.cut((T.resolve_status & 256) > 0)
    T.cut(T.objc_type == 3)

    def casToPhoto(casval):
        flag1val = 0
        flag2val = 0
        casvals = 0
        for k,v in cas_flags.items():
            if v & casval:
                print('CAS flag', k)
                casvals |= v
                if k in photo_flags1_map:
                    flag1val |= photo_flags1_map[k]
                elif k in photo_flags2_map:
                    flag2val |= photo_flags2_map[k]
                else:
                    print('Flag not found:', k)
                    assert(False)
        assert(casvals == casval)
        print('Flag values 0x%x 0x%x' % (flag1val,flag2val))
        return flag1val,flag2val

    f1,f2 = casToPhoto(0x10000000)
    assert(f2 == 0)
    T.cut((T.flags_r & f1) > 0)
    T.cut((T.flags_i & f1) > 0)

    f1,f2 = casToPhoto(0x800a0)
    assert(f2 == 0)
    T.cut((T.flags_r & f1) == 0)
    T.cut((T.flags_i & f1) == 0)

    f1,f2 = casToPhoto(0x400000000000)
    assert(f1 == 0)
    T.cut(np.logical_or((T.flags2_r & f2) == 0, T.psfmagerr[:,2] <= 0.2))
    T.cut(np.logical_or((T.flags2_i & f2) == 0, T.psfmagerr[:,3] <= 0.2))

    print(len(T), 'pass cut')
    if len(T) == 0:
        continue

    # Re-read SDSS photoObjs and grab ALL columns.
    inds = T.index
    T = fits_table(fn, rows=inds, columns=['ra','dec','raerr','decerr',
                                        'cmodelmag', 'cmodelmagerr',
                                        'psfmag','psfmagerr', 'flags', 'flags2', 'objid'])

    unwdir = '/project/projectdirs/cosmo/data/unwise/unwise-phot/sdss-collab/sdss-dr13-pobj'
    fn = os.path.join(unwdir, '%i'%run, '%i'%camcol, 'photoWiseForced-%06i-%i-%04i.fits' % (run, camcol, field))
    U = fits_table(fn, rows=inds, columns=['treated_as_pointsource', 'pointsource', 'w1_mag', 'w1_mag_err',
                                        'has_wise_phot', 'objid'])

    print('has_wise_phot:', np.unique(U.has_wise_phot))

    # *almost* all the time, has_wise_phot=True; exception is
    # sdss-dr13-pobj/4135/4/photoWiseForced-004135-4-0169.fits
    # but I haven't investigated.
    T.cut(U.has_wise_phot)
    U.cut(U.has_wise_phot)

    assert(np.all(T.objid == U.objid))

    U.delete_column('objid')
    T.add_columns_from(U)

    TT.append(T.copy())
    print('Total of', sum([len(x) for x in TT]), 'sources')
    if i - ilast >= 1000:
        ilast = i
        T = merge_tables(TT)
        print('Total of', len(T), 'pass cut')
        N = 10000
        if len(T) > N:
            fn = '%s-%03i.fits' % (tag,k)
            k += 1
            T.writeto(fn)
            print('Wrote', len(T), 'to', fn)
            TT = []
        else:
            T.writeto('%s-x.fits' % tag)
            TT = [T]

T = merge_tables(TT)
fn = '%s-%03i.fits' % (tag,k)
k += 1
T.writeto(fn)
print('Wrote', len(T), 'to', fn)

