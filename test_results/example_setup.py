from distutils.core import setup, Extension

modules1 = Extension("_example", \
                     sources =["example_interface_wrap.c",\
                                      "example.c"])

setup( name = 'Hello world example app',\
       version = '1.0',
       description = 'simple example from swig',
       ext_modules = [modules1])
