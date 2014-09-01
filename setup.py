from distutils.core import setup

setup(
    name='pyGeoNet 1.0',
    version='1.0.0',
    author='Harish Sangireddy',
    author_email='hsangireddy@utexas.edu',
    packages=['pygeonet', 'pygeonet.test'],
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    #url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE.txt',
    description='Geomorphic feature extraction toolbox',
    long_description=open('C:\\Users\\Harish\\Dropbox\\pyGeoNet1.0\\README.txt').read(),
    install_requires=[
        "Django >= 1.1.1",
        "caldav == 0.1.4",
    ],
)
