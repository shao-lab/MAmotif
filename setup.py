import os
from setuptools import setup, find_packages
from mamotif import __version__

INSTALL_REQUIRES = ['MAnorm>=1.0', 'motifscan>=1.1']

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Environment :: Console',
               'Intended Audience :: Education',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Unix',
               'Programming Language :: Python',
               'Programming Language :: Python :: 2',
               'Programming Language :: Python :: 2.7',
               'Topic :: Scientific/Engineering :: Bio-Informatics',
               ]

with open('README.md', 'r') as fin:
    long_description = fin.read()

setup(name="MAmotif",
      version=__version__,
      packages=find_packages(),
      scripts=['bin/mamotif','bin/mamotif-integrate'],
      author="Hayden Sun",
      author_email="sunhongduo@picb.ac.cn",
      description="",
      long_description=long_description,
      url="",
      license='',
      install_requires=INSTALL_REQUIRES,
      classifiers=CLASSIFIERS,
      zip_safe=False
      )
