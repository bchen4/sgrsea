"""Description

Setup script for crispr

Copyrigh (c) Beibei Chen <beibei.chen@utsouthwestern.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (See the file COPYING included with
the distribution)

@status: alpha
@version: $Development$
@author: Beibei Chen
@contact: beibei.chen@utsouthwestern.edu

"""
import os
import sys
from setuptools import setup, Extension, find_packages

command_classes = {}

#try:
#  from numpy import get_include as numpy_get_include
#  numpy_include_dir = [numpy_get_include()]
#except:
#  numpy_include_dir = []
#  sys.stderr.write("CRITICAL: Numpy must be installed!\n")
#  sys.exit(1)
#
#try:
#  from pandas import get_include as pandas_get_include
#  pandas_include_dir = [pandas_get_include()]
#except:
#  pandas_include_dir = []
#  sys.stderr.write("CRITICAL: Pandas must be installed!\n")
#  sys.exit(1)
#
#try:
#  from statsmodels import get_include as statsmodels_get_include
#  statsmodels_include_dir = [statsmodels_get_include()]
#except:
#  statsmodels_include_dir = []
#  sys.stderr.write("CRITICAL: Statsmodels must be installed!\n")
#  sys.exit(1)

def main():
  #if not 2.7 < float(sys.version[:3]) < 2.8:
  #  sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
  #  sys.exit(1)

  #ext_modules = 
  setup(name='sgrsea',
      version='0.1',
      description='Identify enriched genes in CRISPR-Cas9 experiment',
      #url='https://github.com/QBRC/Crispr',
      author='Beibei Chen',
      author_email='beibei.chen@utsouthwestern.edu',
      license='BSD',
      package_dir={'sgrsea':'sgrsea'},#,'bin':'bin','test':'test'},
      packages=['sgrsea'],
      scripts=['bin/sgrsea'],
      classifiers=[
        'Development Status :: 1 - Dev',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python',
        ],
      install_requires=[
        'numpy >= 1.7.0',
        'pandas >= 0.15.0',
        'statsmodels',
        ],
      cmdclass = command_classes,
      #ext_modules = ext_modules
      )

if __name__ == '__main__':
  main()



