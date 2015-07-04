from setuptools import setup, find_packages
import sys
import os.path

import numpy as np

# Must be one line or PyPI will cut it off
DESC = ("A powerful, accurate, and easy-to-use Python library for "
        "doing colorspace conversions")

LONG_DESC = open("README.rst").read()

# defines __version__
exec(open("colorspacious/version.py").read())

setup(
    name="colorspacious",
    version=__version__,
    description=DESC,
    long_description=LONG_DESC,
    author="Nathaniel J. Smith",
    author_email="njs@pobox.com",
    url="https://github.com/njsmith/colorspacious",
    license="MIT",
    classifiers =
      [ "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        ],
    packages=find_packages(),
    install_requires=["numpy"],
)
