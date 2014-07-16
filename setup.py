from setuptools import setup, find_packages
import sys
import os.path

import numpy as np

# Must be one line or PyPI will cut it off
DESC = ("Compute perceptual similarity between sRGB colors according to the "
        "CAM02-UCS formula given by Luo et al (2006)")

LONG_DESC = open("README.rst").read()

# defines __version__
exec(open("pycam02ucs/version.py").read())

setup(
    name="pycam02ucs",
    version=__version__,
    description=DESC,
    long_description=LONG_DESC,
    author="Nathaniel J. Smith",
    author_email="njs@pobox.com",
    url="https://github.com/njsmith/pycam02ucs",
    license="MIT",
    classifiers =
      [ "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        ],
    packages=find_packages(),
    install_requires=["numpy"],
)
