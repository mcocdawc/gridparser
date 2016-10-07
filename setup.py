#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Setup file for the chemcoord package.
"""

from __future__ import with_statement
from __future__ import absolute_import
from setuptools import setup, find_packages
from io import open


MAIN_PACKAGE = 'gridparser'
DESCRIPTION = "Python module for dealing with chemical coordinates."
VERSION = '0.0.0'
LICENSE = 'LGPLv3'
AUTHOR = 'Oskar Weser'
EMAIL = 'oskar.weser@gmail.com'
URL = 'https://github.com/mcocdawc/chemcoord'
REQUIRES = [
    'numpy', 'pandas', 'copy', 'math', 'collections', 'os', 'sys',
    'io']
KEYWORDS = ['molcas']

CLASSIFIERS = [
    'Development Status :: 2 - Pre-Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    (
    'License :: OSI Approved :: '
    'GNU Lesser General Public License v3 (LGPLv3)'),
    'Natural Language :: English',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Physics']


def readme():
    '''Return the contents of the README.md file.'''
    with open('README.md') as freadme:
        return freadme.read()


def setup_package():
    setup(
        name=MAIN_PACKAGE,
        version=VERSION,
        url=URL,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email=EMAIL,
        include_package_data=True,
        keywords=KEYWORDS,
        license=LICENSE,
        long_description=readme(),
        classifiers=CLASSIFIERS,
        packages=find_packages(),
        requires=REQUIRES
    )


if __name__ == "__main__":
    setup_package()
