import shutil
import subprocess
import sys
import os
import numpy as np
from matplotlib import cm
import prepare_pygeonet_defaults as defaults
import prepare_pygeonet_inputs as parameters
import pygeonet_rasterio as pyg_rio
import pygeonet_plot as pyg_plt
import grass.script as g
import grass.script.setup as gsetup


# Flow accumulation is computed by calling GRASS GIS functions.
def flowaccumulation(filteredDemArray):
    """
    Runs the flow accumulation section of the Geonet
    :param filteredDemArray: The filtered dem array
    :return: The flow accumulation results
    """
    print(sys.platform)
    if sys.platform.startswith('win'):
        # MS Windows
        grass7bin = r'C:\Program Files\GRASS GIS 7.2.1\grass72.bat'
        # uncomment when using standalone WinGRASS installer
        # grass7bin = r'C:\Program Files (x86)\GRASS GIS 7.2.0\grass72.bat'
        # this can be avoided if GRASS executable is added to PATH
    elif sys.platform == 'darwin':
        # Mac OS X
        # TODO: this have to be checked, maybe unix way is good enough
        grass7bin = '/Applications/GRASS-7.4.1.app/Contents/MacOS/Grass'

    mswin = sys.platform.startswith('win')
    if mswin:
        gisdbdir = os.path.join(os.path.expanduser("~"), "Documents\grassdata")
    else:
        gisdbdir = os.path.join(os.path.expanduser("~"), "grassdata")
    locationGeonet = 'geonet'
    mapsetGeonet = 'geonetuser'
    if sys.platform.startswith("win"):
        import ctypes
        SEM_NOGPFAULTERRORBOX = 0x0002  # From MSDN
        ctypes.windll.kernel32.SetErrorMode(SEM_NOGPFAULTERRORBOX);
        CREATE_NO_WINDOW = 0x08000000  # From Windows API
        subprocess_flags = CREATE_NO_WINDOW
    else:
        subprocess_flags = 0
    # -----------------------------------------------------------------
    startcmd = [grass7bin, '--config path']
    p = subprocess.Popen(startcmd, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >> sys.stderr, "ERROR: Cannot find GRASS GIS 7 " \
                             "start script (%s)" % startcmd
        sys.exit(-1)

    # say hello
    g.message('--- GRASS GIS 7: Current GRASS GIS 7 environment:')
    print g.gisenv()

    gisbase = out.strip('\n\r')
    gisdb = os.path.join(os.path.expanduser("~"), "grassdata")
    mswin = sys.platform.startswith('win')
    if mswin:
        gisdbdir = os.path.join(os.path.expanduser("~"), "Documents\grassdata")
    else:
        gisdbdir = os.path.join(os.path.expanduser("~"), "grassdata")
    originalGeotiff = os.path.join(parameters.demDataFilePath, parameters.demFileName)
    geotiff = parameters.pmGrassGISfileName
    locationGeonet = 'geonet'
    grassGISlocation = os.path.join(gisdbdir, locationGeonet)
    if os.path.exists(grassGISlocation):
        print("Cleaning existing Grass location")
        shutil.rmtree(grassGISlocation)
    gsetup.init(gisbase, gisdbdir, locationGeonet, 'PERMANENT')
    mapsetGeonet = 'geonetuser'
    print('Making the geonet location')
    g.run_command('g.proj', georef=geotiff, location=locationGeonet)
    print('Existing Mapsets after making locations:')
    g.read_command('g.mapsets', flags='l')
    print('Setting GRASSGIS environ')
    gsetup.init(gisbase, gisdbdir, locationGeonet, 'PERMANENT')

    print('Making mapset now')
    g.run_command('g.mapset', flags='c', mapset=mapsetGeonet, \
                  location=locationGeonet, dbase=gisdbdir)
    # gsetup initialization 
    gsetup.init(gisbase, gisdbdir, locationGeonet, mapsetGeonet)

    # Read the filtered DEM
    print('Import filtered DEM into GRASSGIS and '
          'name the new layer with the DEM name')
    tmpfile = parameters.demFileName  # this reads something like skunk.tif
    geotiffmapraster = tmpfile.split('.')[0]
    print('GRASSGIS layer name: {}'.format(geotiffmapraster))
    g.run_command('r.in.gdal', input=geotiff,
                  output=geotiffmapraster, overwrite=True)
    gtf = parameters.geotransform

    # Flow computation for massive grids (float version)
    print("Calling the r.watershed command from GRASS GIS")
    subbasinThreshold = defaults.thresholdAreaSubBasinIndexing
    if (not hasattr(parameters, 'xDemSize')) or (not hasattr(parameters, 'yDemSize')):
        parameters.yDemSize = np.size(filteredDemArray, 0)
        parameters.xDemSize = np.size(filteredDemArray, 1)
    if parameters.xDemSize > 4000 or parameters.yDemSize > 4000:
        print('using swap memory option for large size DEM')
        g.run_command('r.watershed', flags='am', overwrite=True, \
                      elevation=geotiffmapraster, \
                      threshold=subbasinThreshold, \
                      drainage='dra1v23')
        g.run_command('r.watershed', flags='am', overwrite=True, \
                      elevation=geotiffmapraster, \
                      threshold=subbasinThreshold, \
                      accumulation='acc1v23')
    else:
        g.run_command('r.watershed', flags='a', overwrite=True, \
                      elevation=geotiffmapraster, \
                      threshold=subbasinThreshold, \
                      accumulation='acc1v23', \
                      drainage='dra1v23')

    print('Identify outlets by negative flow direction')
    g.run_command('r.mapcalc', overwrite=True, \
                  expression='outletmap = if(dra1v23 >= 0,null(),1)')
    print('Convert outlet raster to vector')
    g.run_command('r.to.vect', overwrite=True, \
                  input='outletmap', output='outletsmapvec', \
                  type='point')

    print('Delineate basins according to outlets')
    g.run_command('r.stream.basins', overwrite=True, \
                  direction='dra1v23', points='outletsmapvec', \
                  basins='outletbains')
    # Save the outputs as TIFs
    outlet_filename = geotiffmapraster + '_outlets.tif'
    g.run_command('r.out.gdal', overwrite=True, \
                  input='outletmap', type='Float32', \
                  output=os.path.join(parameters.geonetResultsDir,
                                      outlet_filename), \
                  format='GTiff')
    outputFAC_filename = geotiffmapraster + '_fac.tif'
    g.run_command('r.out.gdal', overwrite=True, \
                  input='acc1v23', type='Float64', \
                  output=os.path.join(parameters.geonetResultsDir,
                                      outputFAC_filename), \
                  format='GTiff')
    outputFDR_filename = geotiffmapraster + '_fdr.tif'
    g.run_command('r.out.gdal', overwrite=True, \
                  input="dra1v23", type='Float64', \
                  output=os.path.join(parameters.geonetResultsDir,
                                      outputFDR_filename), \
                  format='GTiff')
    outputBAS_filename = geotiffmapraster + '_basins.tif'
    g.run_command('r.out.gdal', overwrite=True, \
                  input="outletbains", type='Int16', \
                  output=os.path.join(parameters.geonetResultsDir,
                                      outputBAS_filename), \
                  format='GTiff')

    # this reads something like skunk.tif
    tmpfile = parameters.demFileName
    geotiffmapraster = tmpfile.split('.')[0]
    print('GRASSGIS layer name: {}'.format(geotiffmapraster))
    gtf = parameters.geotransform
    # Save the outputs as TIFs
    outlet_filename = geotiffmapraster + '_outlets.tif'
    outputFAC_filename = geotiffmapraster + '_fac.tif'
    outputFDR_filename = geotiffmapraster + '_fdr.tif'
    outputBAS_filename = geotiffmapraster + '_basins.tif'
    # plot the flow directions
    nanDemArrayfdr = pyg_rio.read_geotif_generic(parameters.geonetResultsDir, outputFDR_filename)
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(nanDemArrayfdr, 'Flow directions DEM')
    """
    Output drainage raster map contains drainage direction.
    Provides the "aspect" for each cell measured CCW from East.
    Multiplying positive values by 45 will give the direction
    in degrees that the surface runoff will travel from that cell.
    The value 0 (zero) indicates that the cell is a depression area
    (defined by the depression input map).
    
    Negative values indicate that surface runoff is leaving the boundaries
    of the current geographic region. The absolute value of these
    negative cells indicates the direction of flow.
    """
    outlets = np.where(nanDemArrayfdr < 0)
    print("Number of outlets : {}".format(str(len(outlets[0]))))
    print ([[outlets[0][i], outlets[1][i]] for i in range(len(outlets[0]))])
    # plot the flow accumulation
    nanDemArrayfac = pyg_rio.read_geotif_generic(parameters.geonetResultsDir, outputFAC_filename)
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(nanDemArrayfac, 'Flow accumulations DEM')
    # getting the bigbasins from the r.streams.basins modules
    nanDemArraybasins = pyg_rio.read_geotif_generic(parameters.geonetResultsDir,
                                                    outputBAS_filename)
    nanDemArraybasins[np.isnan(filteredDemArray)] = 0
    # write outlet info into a csv file
    outlet_tablename = geotiffmapraster + '_outlets.csv'
    outlet_tablelocation = os.path.join(parameters.geonetResultsDir, outlet_tablename)
    if os.path.exists(outlet_tablelocation):
        os.remove(outlet_tablelocation)
    with open(outlet_tablelocation, 'a') as f:
        f.write('BasinID,YIndex,XIndex\n')
        for i in range(len(outlets[0])):
            f.write(str(nanDemArraybasins[outlets[0][i], outlets[1][i]]) + \
                    ',' + str(outlets[0][i]) + ',' + str(outlets[1][i]) + '\n')
    # outlets locations in projection of the input dataset
    outletsxx = outlets[1]
    outletsxxfloat = [float(x) + 0.5 for x in outletsxx]
    outletsyy = outlets[0]
    outletsyyfloat = [float(x) + 0.5 for x in outletsyy]
    """
    outletsxxProj = np.array(outletsxxfloat) * parameters.demPixelScale + \
                    parameters.xLowerLeftCoord + float(0.0164)
    outletsyyProj = parameters.yLowerLeftCoord - np.array(outletsyyfloat) * \
                    parameters.demPixelScale + \
                    parameters.yDemSize * parameters.demPixelScale + float(0.0155)
    
    # The extra decimal digits is essentially a hack into
    # Grass GIS r.water.outlet routine, which only, works
    # with atleast 4 significant digits
    """
    outletsxxProj = float(gtf[0]) + float(gtf[1]) * np.array(outletsxxfloat)
    outletsyyProj = float(gtf[3]) + float(gtf[5]) * np.array(outletsyyfloat)
    # plotting log10 flow accumulation with outlets
    drainageMeasure = np.log10(nanDemArrayfac)
    if defaults.doPlot == 1:
        pyg_plt.raster_point_plot(drainageMeasure, outlets, 'flowArray with outlets')
    # plotting subbasins with outlets
    if defaults.doPlot == 1:
        pyg_plt.raster_point_plot(nanDemArraybasins, outlets, 'basinIndexArray with outlets',
                                  cm.Dark2)
    return {'outlets': outlets,
            'fac': nanDemArrayfac,
            'fdr': nanDemArrayfdr,
            'outletsxxProj': outletsxxProj,
            'outletsyyProj': outletsyyProj,
            'bigbasins': nanDemArraybasins}
    # end of flow accumulation
