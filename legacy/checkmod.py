#! /usr/bin/env python
# This small script will inform you about the installed packages on
# your local machine.
# Please make sure you have all the required libraries installed
# pyGeoNet
# Author: Harish Sangireddy
# Austin, Texas

import pip

installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version)
     for i in installed_packages])
print(installed_packages_list)
