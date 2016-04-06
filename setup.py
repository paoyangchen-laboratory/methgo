#!/usr/bin/env python
from setuptools import setup

def readme():
    with open("README.rst") as infile:
        return infile.read()

setup(name="MethGo",
      version="0.5.0",
      description="A comprehensive tool for analyzing whole-genome bisulfite sequencing data",
      long_description=readme(),
      classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      keywords="whole genome bisulfite sequencing",
      url="http://paoyangchen-laboratory.github.io/methgo",
      author="Wen-Wei Liao",
      author_email="gattacaliao@gmail.com",
      license="MIT",
      pakages=["methgo"],
      install_requires=[
        "cython",
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "pysam",
        "biopython",
        "pyfasta",
        "pybedtools",
      ],
      entry_points={
        'console_scripts': ["methgo=methgo.command_line:main"],
      },
      include_package_data=True,
      zip_safe=False)
