"""
pyGeoNet

This is the python implementation of GeoNet originally written in MATLAB
and maintained at https://sites.google.com/site/geonethome/ 

This version of GeoNet follows the same algorithm as GeoNet with a few differences enumerated below.

1. The geotiff files i.e. lidar datasets are read using gdal libraries.

2. The internal and external hull created within GeoNet are not created in
pyGeoNet. This step has no effect on the final results obtained.

3. The slope and curvature are computed using gradient function in numpy
library and this might have some differences.

4. The biggest difference between GeoNet[MATLAB] and pyGeoNet is the
flow accumulation. GRASS GIS is used in the background to read the
filtered DEM and then compute flow accumulation using the r.watershed
functionality in GRASS GIS. This might produce different flow accumulation
values than GeoNet[MATLAB].

5. The fast marching algorithm used for computing geodesics is also different.
pyGeoNet uses scikit-skfmm library to compute the cost function.

6. The skeleton threshold used to remove rouge skeleton pixels is also not
the same as in GeoNet. This is more of a bug in pyGeoNet rather than a
difference.


Author: Harish Sangireddy
email: harish2rb@gmail.com


"""
