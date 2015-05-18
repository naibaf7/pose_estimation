README for the 3D Project: Robust large scale localization

Use the Makefile to compile the project. The program has command line options and should be self-explanatory.

-== LIBRARIES: ==-
OpenMP
FLANN, Version 1.8.4
OpenCV (Core, Highgui, Imgproc), Version 2.4.9


-== SOURCE CODE: ==-
The code is structured with the following files (excluding code we did not change or code ourselves, and those not critically relevant):
-pose_estimation.cpp: Wraps all the high level functions and the command line interface.
-benchmark.cpp: Contains the code for benchmarking the two approaches.
-import_export.cpp: Code for result visualization and import/export. The output is Meshlab compatible.
-P3p.cpp and P4pf.cpp: Code for the pose estimation calculations given 3 points with known focal length or 4 points with unknown focal length. Reprogrammed by us but with given starting code.
-pose_utils.cpp: Small helper source for random functions and string modification operations.
query_loader.cpp: Contains the code to load a query image with all parameters (image size, focal length, SIFT descriptors). It also contains the query class that stores final pose solutions. Additionally, the query loader can automatically grab the original images from flickr.com.
-query_processor.cpp: Common class for the basic and advanced approach, including the reprojection/inlier test.
-query_processor_basic.cpp: Basic approach implementation of RANSAC.
-query_processor_advanced.cpp: Advanced approach, including co-occurence and backmatching. This file in which the algorithm by Y. Li et al. is implemented.


-== PROJECT DATA: ==-
The project data is distributed separately. It is best to use the Dubrovnik datasets, which are known to work.
See here for the bundler output:
http://www.cs.cornell.edu/projects/p2f/
Additionally, an INFO file (binary) of the dataset without the query poses is required. We cannot provide this file, due to its large size.

It is explained in the "tmp", "out" and "data" folders README how to extract and place the data in order to get it working.
