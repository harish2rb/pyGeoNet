#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Harish
#
# Created:     23/08/2014
# Copyright:   (c) Harish 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("pygeonet_diffusion", ["pygeonet_diffusion.pyx"])]

setup( name = 'Hello world app', \
cmdclass = {'build_ext': build_ext}, \
ext_modules = ext_modules)
