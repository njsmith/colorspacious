from setuptools import setup, Extension, find_packages
import sys
import os.path

if os.path.exists(".this_is_a_checkout"):
    USE_CYTHON = True
else:
    # Don't depend on Cython in builds-from-sdist
    USE_CYTHON = False

# Must be one line or PyPI will cut it off
DESC = ("Compute perceptual similarity between sRGB colors according to the "
        "CAM02-UCS formula given by Luo et al (2006)")

LONG_DESC = open("README.rst").read()

if USE_CYTHON:
    cython_ext = "pyx"
else:
    cython_ext = "c"
ext_modules = [
    Extension("pycam02ucs._ciecam02", ["pycam02ucs/_ciecam02.%s" % (cython_ext,)])
]
if USE_CYTHON:
    from Cython.Build import cythonize
    #import pdb; pdb.set_trace()
    ext_modules = cythonize(ext_modules)

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
    ext_modules=ext_modules,
)
