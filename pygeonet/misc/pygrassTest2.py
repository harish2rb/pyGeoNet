import sys, os, numpy, argparse

gisbase = os.environ['GISBASE'] = "C:\GRASSGIS7svn"
if not os.environ.has_key('HOME'):
         os.environ['HOME'] = os.path.join(os.environ['HOMEDRIVE'], \
                                           os.environ['HOMEPATH'])
gisdbase = os.path.join(os.environ['HOME'])
sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))


import grass.script as g
import grass.script.setup as gsetup

class grasswrapper():
	def __init__(self, dbBase="", location="", mapset="PERMANENT"):
		'''
		Wrapper of "python.setup.init" defined in GRASS python.
		Initialize system variables to run scripts without starting GRASS explicitly.

		@param dbBase: path to GRASS database (default: '').
		@param location: location name (default: '').
		@param mapset: mapset within given location (default: 'PERMANENT')

		@return: Path to gisrc file.
		'''
		self.gisbase = os.environ['GISBASE']
		self.gisdb = dbBase
		self.loc = location
		self.mapset = mapset
		gsetup.init(self.gisbase, self.gisdb, self.loc, self.mapset)


print g.run_command("g.list", _type="vect")
