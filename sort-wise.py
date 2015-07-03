import os
from astrometry.util.fits import *
import fitsio

indir = 'allwise-cats'
outdir = 'data/allwise-sorted'

rdpat = os.path.join(indir, 'wise-allwise-cat-part%02i-radec.fits')
catpat = os.path.join(indir, 'wise-allwise-cat-part%02i.fits')
outrdpat = os.path.join(outdir, 'wise-allwise-cat-part%02i-radec.fits')
outpat = os.path.join(outdir, 'wise-allwise-cat-part%02i.fits')

for part in range(18, 49):
    fn = rdpat % part
    T = fits_table(fn)
    print 'Read', len(T), 'from', fn
    I = np.argsort(T.ra)
    Nchunk = 100000
    #Nchunk = 10
    outrdfn = outrdpat % part
    outfn = outpat % part
    #if os.path.exists(outrdfn):
    #    os.unlike(outrdfn)
    #if os.path.exists(outfn):
    #    os.unlike(outfn)
    #I = I[:95]

    T[I].writeto(outrdfn)
    del T

    catfn = catpat % part

    first = True
    row0 = 0
    while row0 < len(I):
        rows = I[row0:row0+Nchunk]
        print 'Rows', row0, 'to', row0+len(rows)-1
        #print 'rows:', rows
        #T[rows].writeto(outrdfn, append=not(first), append_to_hdu=-1)
        print 'Reading', catfn
        C = fits_table(catfn, rows=rows)
        print 'Writing', outfn
        C.writeto(outfn, append=not(first), append_to_hdu=-1)
        row0 += len(rows)
        #I = I[Nchunk:]
        first = False

    break
