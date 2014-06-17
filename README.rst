pycam02ucs
==========

Compute perceptual similarity between sRGB colors according to the
CAM02-UCS formula given by:

Luo, M. R., Cui, G., & Li, C. (2006). Uniform colour spaces based on
CIECAM02 colour appearance model. Color Research & Application, 31(4),
320–330. doi:10.1002/col.20227

We can also calculate the "large color difference" (LCD) and "short
color difference" (SCD) metrics estimated by the same publication. The
UCS ("uniform color space") metric is defined as a compromise between
these two.

CIECAM02 calculations are performed using code released by Billy Biggs
under the MIT license, and downloaded from: http://scanline.ca/ciecam02/

.. image:: https://travis-ci.org/njsmith/pycam02ucs.png?branch=master
   :target: https://travis-ci.org/njsmith/pycam02ucs
.. image:: https://coveralls.io/repos/njsmith/pycam02ucs/badge.png?branch=master
   :target: https://coveralls.io/r/njsmith/pycam02ucs?branch=master

Documentation:
  http://pycam02ucs.readthedocs.org/ (TODO)

Installation:
  ``pip install pycam02ucs`` (requires a C compiler)

Downloads:
  http://pypi.python.org/pypi/pycam02ucs/ (TODO)

Code and bug tracker:
  https://github.com/njsmith/zs

Contact:
  Nathaniel J. Smith <nathaniel.smith@ed.ac.uk>

Developer dependencies (only needed for hacking on source):
  * Cython: needed to build from checkout
  * nose: needed to run tests

License:
  MIT, see LICENSE.txt for details.
