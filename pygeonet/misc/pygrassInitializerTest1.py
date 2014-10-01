import os, sys

GISBASE="C:\GRASSGIS7svn"
os.environ['GISBASE'] = GISBASE

gisbase = os.environ['GISBASE'] = GISBASE
gisdbase = 'C:\\Users\\Harish\\Documents\\grassdata\\'#os.path.join(os.environ['HOME'])
location = "australia"
mapset   = "PERMANENT"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

##if "C:\\GRASSGIS7svn\\etc\\python" not in sys.path:
##    print "Adding path"
##    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
##    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python","grass"))
##    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python","grass","script"))
##    sys.path.append(os.path.join(os.environ['GISBASE'], "lib"))
##    #sys.path.append(os.path.join("C:\\Users\\Harish\\Documents\\grassdata", "lib"))
##    #os.environ['GISDBASE'] = "C:\\Users\\Harish\\Documents\\grassdata"
##else:
##    print "Path present"

import grass.script as g
from grass.script import setup as gsetup

#print gsetup.write_gisrc(gisdbase, location, mapset)
#print gsetup.init(gisbase, gisdbase, location, mapset)
print g.gisenv()
print g.read_command('v.in.ogr', dsn="C:\predictAg\Rx\farmshapefile.shp", layer="a", output="a")
print g.parse_command('v.hull', input="a", output="b")
print g.parse_command('v.out.ogr', input="b", dsn="C:\\Python_Test\\farmGrass.shp", format="ESRI_Shapefile")
#grass.run_command("v.hull", input="C:\\predictAg\\Rx\farmshapefile.shp", output="C:\\Python_Test\5.shp")
