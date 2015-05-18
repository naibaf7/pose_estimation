=====================================================
| README file for                                   |
|   Cornell 3D Location Recognition Datasets        |
|   http://www.cs.cornell.edu/projects/3Dlocation/  |
|    v1.0, April 2012                               |
=====================================================

This document describes the contents and format of the Cornell 3D
Location Recognition Datasets.  This dataset is designed especially
for algorithms for 2D-to-3D matching, and as such includes a 3D SfM
model for each dataset.

If you use this dataset in a publication, please cite the paper:

  Yunpeng Li, Noah Snavely, and Daniel P. Huttenlocher.  "Location
  Recognition using Prioritized Feature Matching."  Proceedings of
  ECCV 2010.

Please read the following carefully before using the dataset, as there
are several important notes.  *See especially the special notes on the
Dubrovnik6K dataset at the end of this document.*

-== I. Contents of each dataset ==-

Each dataset consists of a database -- including SIFT keys a 3D
structure-from-motion (SfM) model -- and a set of query images,
consisting of a set of SIFT keys.  Out of respect for Flickr users, we
are not distributing the original images with this archive, but we can
send the images separately if you require access to them.  (Many of
the images are available on Flickr as well.)

In particular, the provided files are as follows:

  list.db.txt -- the list of database image files
  list.query.txt -- the list of query image files

  db/ -- directory containing *squeezed* database SIFT key files (see below)
  query/ -- directory containing query SIFT key files

  bundle/ -- directory with 3D SfM models
    bundle/list.orig.txt -- full list of images used to build the SfM model
    bundle/bundle.orig.out -- original SfM model, before removing test images
    bundle/bundle.db.out -- database SfM model, after removing test images
                            (use this model for testing if you would like
                            to use our SfM results)


-== II. Description of file formats ==-

This section describes file formats for each of the data types above,
with special notes as appropriate.

  A. List files (list.db.txt, list.query.txt).
  
     List files specify filenames to images in jpg format, one per
     line (keep in mind that the actual jpg files are not distributed
     unless requested).  In addition, if the focal length of the image
     has been estimated from Exif tags, then that is also included.

     Images without known focal length information are specified with
     a line with a single field, the image name.  Example:

       query/10970812@N05_2553027508.jpg

     Images with known focal length information are specified with a
     line with three fields: the image name, a zero, and the Exif
     focal length.  (The second field is always zero but may change in
     future datasets.)  Example:

       query/11289373@N03_2733280477.jpg 0 1280.00000

  B. Key files

     We are providing (gzipped) SIFT key files for the database images
     and query images.  One important note is that the database keys
     are *not* the original set of SIFT keys, but the set of SIFT keys
     after removing features that did not end up getting used for
     reconstruction (i.e., SIFT features that do not correspond to a
     point in the 3D model).  We will make the *original* database key
     files available as a separate download (see the project webpage).

     The format of the key files is the same as that produced by David
     Lowe's SIFT extractor, described here:

       http://www.cs.ubc.ca/~lowe/keypoints/

  C. 3D SfM models

     We provide two 3D SfM models, computed using Bundler.  The
     formats of these files are described here:

     http://phototour.cs.washington.edu/bundler/bundler-v0.4-manual.html#S6

     Description of files:

       - bundle.orig.out: Original bundle file, containing camera
         parameters for all database and query images.  This can be
         used to compare recovered poses for query images to a "gold
         standard".

       - bundle.db.out: Database bundle file with query images removed
         (as well as any links between points and query images).  This
         model is used during testing, to localize new images.
         Because it is registered with bundle.orig.out, camera
         parameters for images localized to bundle.out can be directly
         compared with those in bundle.orig.out.  However, such
         comparisons should be take with a grain of salt, as the scale
         factor for the bundle files may not be meaningful.


-== III. Accessing the original photos ==-  

While we are not distributing the original images with this archive,
there are two ways to access the originals, if you need them.  If you
only need a few images (e.g., to show in a publication), then you can
access photos directly on Flickr as follows: each photo filename has
the form <flickruser>_<flickrphotoID>.jpg.  You can access Flickr
photo <flickrphotoID> through the Web with a URL of the form:

   http://www.flickr.com/photo_zoom.gne?id=<flickrphotoID>

     e.g.:

   http://www.flickr.com/photo_zoom.gne?id=2000000000

If you need all of the photos (for instance, because you'd like to
re-extract features), then please send a request to the email
addresses below.


-== IV. Special notes for the Dubrovnik6K dataset ==-

This dataset consists of a single large component for the Old Town of
Dubrovnik.  This reconstruction is not in any particular georegistered
coordinate system, but is approximately metric---i.e. 1 unit in
'bundle file coordinates' is approximately 1 meter.  This scaling was
performed by approximately aligning the SfM reconstruction to a map.



Questions?  Email snavely@cs.cornell.edu and yuli@cs.cornell.edu.
