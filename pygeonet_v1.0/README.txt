pyGeoNet
Author: Harish Sangireddy
email: harish2rb@gmail.com


This is the python version of GeoNet (pyGeoNet) originally written in MATLAB
and maintained at https://sites.google.com/site/geonethome/ 

This version of GeoNet follows the same algorithm as GeoNet with a few 
differences enumerated below.

1. The geotiff files i.e. lidar datasets are read using gdal libraries.

2. The slope and curvature are computed using gradient function in numpy
library and this might have some differences.

3. The biggest difference between GeoNet[MATLAB] and pyGeoNet is the
flow accumulation. GRASS GIS is used in the background to read the
filtered DEM and then compute flow accumulation using the r.watershed
functionality in GRASS GIS. This might produce different flow accumulation
values than GeoNet[MATLAB].

5. The fast marching algorithm used for computing geodesics is also different.
pyGeoNet uses scikit-skfmm library to compute the cost function.

***********************************************
Known limitations and open issues.
***********************************************
1. Thresholdqxx = 1: This is identified within GeoNet as the deviation point
from a qqplot and should be done the same in python.. Couldn't find a easy
way to do it in python.. but hoping I will solve in beta 2.0 release.

2. The skeleton threshold used to remove rouge skeleton pixels is also not
the same as in GeoNet. 


****************************************************************************
Source Code: https://github.com/harish2rb/pyGeoNet

Requirements: Numpy and a C/C++ compiler (gcc, MinGW, MSVC), Scikit-fmm,

Bugs, questions, patches, feature requests, discussion:
  Email list: https://sites.google.com/site/geonethome/
  Send an email to harish2rb@gmail.com, paola@austin.utexas.edu  for more information.

Installing:
 $ see the installation doc.. this is not a normal python library


****************************************************************************
VERSION HISTORY:
0.0.1: September 1 2014
  Initial Release

0.0.2: September 26th 2014
  Minor changes to code structure.
 
0.0.3 -- 9
	many small bug fixes
	establishing the pygeonet_defaults.py,pygeonet_prepare.py

1.0.0: March 2014
  vectorized the Perona Malik filtering
  compiled and tested for Numpy 1.9
  compiled and tested for Scikit-fmm 0.0.6
  Fixed bugs for GRASS GIS 7.0 stable release
  GRASS GIS functions used:
    g.proj
	g.mapsets
	g.mapset
	gsetup.init
	r.in.gdal
	r.out.gdal
	r.watershed 
		(for datasets > 4000 rows/cols : with disk swap memory option)
	r.maplac
	r.to.vect
	r.stream.basins
	r.drain 
		(for compute discrete geodesics on geodesic distance array)
	r.thin
	v.in.ogr
	v.out.ogr

