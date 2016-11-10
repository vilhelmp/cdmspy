#!/usr/bin/env python
"""Python setup.py file for CDMSpy
"""

from distutils.core import setup

setup(
    name='cdmspy',
    version='0.1',
    author='Magnus Persson',
    author_email='magnusp@vilhelm.nu',
    packages=['cdmspy'],
    license='BSD',
    description='Functions to query the CDMS database and read CDMS ascii files.',
#    install_requires=['astropy'],
)
