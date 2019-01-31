pyGeoNet
========

This is a python version of the MATLAB based feature extraction toolbox called GeoNet
https://sites.google.com/site/geonethome/

The algorithm remains the same in this version, but the python implementation draws and 
builds code from a variety of open sources libraries.

This is work in progress and will be ready for release soon
More updates on installation and usage can be found here:
http://harish2rb.github.io/pyGeoNet/


Running pyGeoNet with Grass GIS on Mac OS steps
* Download GRASS GIS 7.4 (current stable) from https://grass.osgeo.org/download/software/mac-osx/

Make a conda environment with python 2.7 
* `conda create -n grass_geonet27 python=2.7  gdal numpy pandas  requests statsmodels matplotlib`
* `source activate grass_geonet27`
* `pip install scikit-fmm`

Make sure you have your xcode updated
`xcode-select --install`

Open Grass GIS and add `r.stream.basins` to the add ons in GRASS GIS

Now we need to make GRASS GIS python environment available to pyGeonet

`cd pygeonet_V2.1`
`nano $HOME/.grass7/bashrc`

Add the following to bashrc
```
# GRASS GIS Dependencies
export GISBASE="/Applications/GRASS-7.4.1.app/Contents/Resources"
export PATH="$PATH:$GISBASE/bin:$GISBASE/scripts"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GISBASE/lib"
# for parallel session management, we use process ID (PID) as lock file number:
export GIS_LOCK=$$
export GRASS_ADDON_PATH="$GISBASE/addons"
# path to GRASS settings file
export GISRC="$HOME/.grassrc7"
export PYTHONPATH="$PYTHONPATH:$GISBASE/etc/python"
# GRASS WTF variables
export GISDBASE="/Applications/GRASS-7.4.1.app/Contents/Resources"
export LOCATION_NAME=demolocation
export MAPSET=PERMANENT
export GRASS_DB_ENCODING=utf-8
export DEBUG=0
export GUI=text
export GRASS_VERSION='7.4.1'
export LANG=C
```
Source the bashrc
`source $HOME/.grass7/bashrc`

Assuming you are in `grass_geonet27` conda env, and have added all the above steps
Made the necessary changes in the `prepare_pygeonet_inputs.py`
You can run the pygeonet module by typing `python run_pygeonet_processing.py`
