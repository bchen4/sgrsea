"""Description

Setup script for sgRSEA

Copyrigh (c) Beibei Chen <beibei.chen@utsouthwestern.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License (See the file COPYING included with
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

try:
  import numpy
except:
  sys.stderr.write("CRITICAL: Numpy must be installed!\n")
  sys.exit(1)

try:
  import pandas
except:
  sys.stderr.write("CRITICAL: Pandas must be installed!\n")
  sys.exit(1)

try:
  import statsmodels
except:
  sys.stderr.write("CRITICAL: Statsmodels must be installed!\n")
  sys.exit(1)

def main():
  #if not 2.7 < float(sys.version[:3]) < 2.8:
  #  sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
  #  sys.exit(1)

  #ext_modules = 
  setup(name='sgrsea',
      version='v0.0.1',
      description='Enrichment Analysis of CRISPR/Cas9 Knockout Screen Data',
      url='http://bchen4.github.io/sgrsea',
      author='Beibei Chen',
      author_email='beibei.chen@utsouthwestern.edu',
      license='MIT',
      package_dir={'sgrsea':'sgrsea'},#,'test':'test'},
      packages=['sgrsea'],
      scripts=['bin/sgrsea'],
      classifiers=[
        'Development Status :: 1 - Dev',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
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
      )

if __name__ == '__main__':
  main()



