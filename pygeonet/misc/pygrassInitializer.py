import os, sys
import subprocess
from osgeo import gdal

if "C:\\GRASSGIS7svn\\etc\\python" not in sys.path:
    print "Adding path"
    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python","grass"))
    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python","grass","script"))
    sys.path.append(os.path.join(os.environ['GISBASE'], "lib"))
    #sys.path.append(os.path.join("C:\\Users\\Harish\\Documents\\grassdata", "lib"))
    os.environ['GISDBASE'] = "C:\\Users\\Harish\\Documents\\grassdata"
else:
    print "Path present"

#PATH="$PATH:$GISBASE/bin:$GISBASE/script:$GISBASE/lib"
#PYTHONPATH="${PYTHONPATH}:$GISBASE/etc/python/"
#PYTHONPATH="${PYTHONPATH}:$GISBASE/etc/python/grass"
#PYTHONPATH="${PYTHONPATH}:$GISBASE/etc/python/grass/script"
#LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GISBASE/lib"
#GIS_LOCK= "$$"
#GISRC="$HOME/.grassrc6"

grass7bin = 'C:\\GRASSGIS7svn\\grass70svn.bat'
startcmd = grass7bin + ' --config path'
p = subprocess.Popen(startcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
if p.returncode != 0:
    print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
    sys.exit(-1)
gisbase = out.strip('\n')
print 'gisbase :',gisbase

GISBASE="C:\GRASSGIS7svn"
os.environ['GISBASE'] = GISBASE

print sys.path

##if not os.environ.has_key('HOME'):
##         os.environ['HOME'] = os.path.join(os.environ['HOMEDRIVE'], \
##                                           os.environ['HOMEPATH'])
##gisdbase = os.path.join(os.environ['HOME'])


print sys.path

import grass.script as g
import grass.script.setup as gsetup

gisbase = "C:\\GRASSGIS7svn"
gisdb = "C:\\Users\\Harish\\Documents\\grassdata"
location = "australia"
mapset   = "missioncreek"
print gsetup.init(gisbase,gisdb, location, mapset)
print g.gisenv()

##location = 'USA'
##mapset = 'missioncreek'
##
##print 'gsetup'
##print gsetup.write_gisrc(gisbase,location,mapset)
###desc = g.parse_command('db.describe', flags='c', table="bxltot")
###dict.keys(desc)
##
##geotiff = 'c:\\Mystuff\\grassgisdatabase\\dem_2012_mission_v1.tif'
##print geotiff
##ds = gdal.Open(geotiff, gdal.GA_ReadOnly)
##geotransform = ds.GetGeoTransform()
##DEM_arr = ds.GetRasterBand(1).ReadAsArray()
###****************************region adjustment***********************************
### We create a temporary region that is only valid in this python session
##g.use_temp_region()
##rows = DEM_arr.shape[0]
##cols = DEM_arr.shape[1]
##resolution = 1
##print rows,cols
##n = 4928050 #some arbitrary value
##s = n - resolution*rows
##e = 609000  #some arbitrary value
##w = e - resolution*cols
##print "g.region command"
##print g.run_command('g.region', flags = 'd', \
##              n = n ,s = s, e = e, w = w,\
##              res = resolution, rows = rows ,cols = cols)
##
##stop
###Input the DEM ascii file into grass
##iteration = 1
##pathname = 'C:\\Users\\Harish\\Dropbox\\pyGeoNet1.0'
##fullpath = os.path.abspath(pathname)
##output_dir = 'test_results'
##print g.run_command('r.in.gdal', overwrite = True, flags='f', input = geotiff ,\
##              output='test_DEM')
##print g.run_command('r.out.png', input='test_DEM@user1', \
##              output = fullpath + '\\'+output_dir +'\\'+'DEM'+str(iteration))
###Flow computation for massive grids (float version) 
##g.run_command('r.terraflow', overwrite = True, elevation = 'test_DEM@user1', filled = 'flooded_DEM',\
##          direction = 'DEM_flow_direction',swatershed = 'DEM_sink_watershed', \
##              accumulation = 'DEM_flow_accum', tci = 'DEM_tci')
##g.run_command('r.out.png', input='flooded_DEM@user1', \
##              output = fullpath + '\\'+output_dir +'\\'+'flooded_DEM'+str(iteration))
##g.run_command('r.out.png', input='DEM_flow_direction@user1', \
##              output = fullpath + '\\'+output_dir +'\\'+'flow_direction'+str(iteration))
##g.run_command('r.out.png', input='DEM_sink_watershed@user1', \
##              output = fullpath + '\\'+output_dir +'\\'+'watershed'+str(iteration))
##g.run_command('r.out.png', input='DEM_flow_accum@user1', \
##              output = fullpath +'\\'+ output_dir +'\\'+'flow_accumulation'+str(iteration))
##g.run_command('r.out.png', input='DEM_tci@user1', \
##              output = fullpath + '\\'+output_dir +'\\'+'tci'+str(iteration))
##g.run_command('r.out.ascii',flags='h',input='DEM_tci@user1',\
##              output=fullpath +'\\data'+ '\\DEM_flow_accum',null='0')
##f = open(fullpath +'\\data'+ '\\DEM_flow_accum', 'r')
##Flow_accum_arr = numpy.loadtxt(f,unpack = True)
##f.close()
##
#curr_raster='myraster'
#print g.run_command('r.in.gdal',overwrite = True, flags='f',input =geotiff, output='dem_2012_mission_v1')
#print g.run_command('r.out.png', input='myraster@user1', output = 'C:\\Mystuff\\grassgisdatabase\\myraster')
# define output directory for files resulting from GRASS calculation:
#fullpath = 'C:\\Mystuff\\grassgisdatabase\\gisoutput\\'
#gisdb = 'C:\\Mystuff\\grassgisdatabase\\gisoutput'
#try:
#    os.stat(gisdb)
#except:
#    os.mkdir(gisdb)

#print gisdb
# define where to store the resulting files and format
#print g.run_command('g.region', rast='dem_2012_mission_v1@missioncreek')
#g.run_command('r.external.out',directory='C:\\Mystuff\\grassgisdatabase\\gisoutput\\', format='GTiff')

# do the GRASS GIS job
#Flow computation for massive grids (float version) 
#print g.run_command('r.watershed', overwrite = True, elevation = 'dem_2012_mission_v1@Missioncreek',\
#               threshold='2000',accumulation='facv15',drainage='fdrv15',memory='1024')

#g.run_command('r.out.tiff', input='flooded_DEM', output = fullpath +'flooded')
#g.run_command('r.out.tiff', input='DEM_flow_direction', output = fullpath +'direction')
#g.run_command('r.out.tiff', input='DEM_sink_watershed', output = fullpath +'watershed')
#g.run_command('r.out.tiff', input='DEM_flow_accum', output = fullpath +'accum')
#g.run_command('r.out.tiff', input='DEM_tci', output = fullpath +'tci')
