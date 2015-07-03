# unwise-sdss-phot
Forced photometry of the unWISE images given SDSS objects

Results described in 
http://arxiv.org/abs/1410.7397
and available at http://unwise.me




Forced photometry of the WISE images given SDSS detections
----------------------------------------------------------

The data described here are available at:

http://unwise.me/

or for SDSS folks,
riemann:/clusterfs/riemann/raid000/dstn/wise-sdss-forced-phot


Processing summary:
-------------------

I took all SDSS "v5b" survey_primary objects, and found all unWISE
coadd tiles touching them.

The "unWISE coadds" are new coadds of the WISE imaging, described in
[Lang 2014, AJ, 147, 108].

The unWISE tiles have the same centers as the official WISE Atlas
Images; they are on equal-Dec rings, roughly 1.5 x 1.5 degrees in
size, with a bit of overlap.

For each tile, I gathered all the SDSS objects within and near the
tile, and then cross-matched to the WISE catalog ("AllWISE release")
to find WISE-only sources.  I treat the WISE-only sources as point
sources, since no profile-fit measurements in the form I need are
available in the WISE catalog.  For SDSS sources, I apply some tests
for sources that do not have reliably measured sizes, and treat those
as point sources in the forced photometry; there is a flag for that in
the output files.

With that list of sources, I apply the WISE PSF model to predict the
appearance of each object in the WISE image, and fit for a linear
combination of those bases that best match the WISE image.  The
best-fit linear weights are the measured fluxes.  The unWISE coadds
already have a scalar sky background level subtracted, so I do not
allow the sky level to vary during fitting.

After fitting each tile, I need to split the sources back into files
that match the SDSS field-by-field file layout.  SDSS fields that are
entirely contained within one WISE tile are written in one phase, and
then the remaining ones that span multiple tiles are resolved (only
the measurement with the closest tile center is kept) and written out
in a second phase.  SDSS fields that contain no primary objects will
not have a corresponding photoWiseForced file.


Result files:
-------------

The results are row-by-row parallel to the SDSS "v5b" "photoObj"
files; the files

/clusterfs/riemann/raid000/dstn/wise-sdss-forced-phot/1000/1/photoWiseForced-001000-1-0100.fits

and
 
/clusterfs/riemann/raid007/ebosswork/eboss/photoObj.v5b/301/1000/1/photoObj-001000-1-0100.fits

each have 368 rows, and describe row-by-row the same objects.

The "photoWiseForced" files have the following columns:

* has_wise_phot

True for SDSS sources that were photometered in WISE.  The object must
be PRIMARY in SDSS for this to be set.  When this column is False, all
other columns have value zero.

* ra, dec
* objid

Copied straight from the SDSS photoObj file.

* x,y

Zero-indexed pixel position on the unWISE 2048x2048 image tile.

* treated_as_pointsource

The SDSS source is a galaxy (objc_type == 3) but was treated as a
point source for the purposes of forced photometry.  If you want an
optical/WISE color, it would be best to use the SDSS PSF mags, not the
model mags, for these objects.

* pointsource

The SDSS source is a point source (objc_type == 6).

* coadd_id

The unWISE coadd tile name.

And then for each WISE band:

* w1_nanomaggies
* w1_nanomaggies_ivar

The WISE flux and formal error as an inverse-variance.  Take the
formal errors with a grain of salt, of course.  NOTE that these are in
the native WISE photometric system: Vega, NOT AB, so some would say
these should not be called "nanomaggies".  A source with magnitude
22.5 in the VEGA system would have a w1_nanomaggies flux of 1.

The WISE team's suggested conversions to AB are available here:
http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab

And are (or were, when this was written),

mag_AB = mag_Vega + dm,

Band  dm
W1    2.699
W2    3.339
W3    5.174
W4    6.620

* w1_mag
* w1_mag_err

Vega magnitude and formal error in the forced photometry.  These are
simple conversions from the "nanomaggies" columns.

*  w1_prochi2
*  w1_pronpix

Profile-weighted chi-squared and number-of-pixels values.
"Profile-weighted" means that these are weighted according to the
brightness of the source in the WISE images; they are supposed to mean
"how bad is the fit, where I am?"  The "pronpix" is probably not
useful in this dataset -- it effectively counts what fraction of the
source was inside the image, and should be nearly 1. everywhere.

* w1_proflux
* w1_profracflux

"proflux": profile-weighted, how much flux comes from other sources in
the model?  "profracflux" = proflux divided by the flux of this
source.  This is an indication of blending or confusion; a source with
a large profracflux has a significant contribution from a nearby
source.

* w1_npix

How many pixels were included in the fit.


