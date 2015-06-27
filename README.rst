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

We also have a complete CIECAM02 implementation. Maybe that should get
higher billing.

.. image:: https://travis-ci.org/njsmith/pycam02ucs.png?branch=master
   :target: https://travis-ci.org/njsmith/pycam02ucs
.. image:: https://coveralls.io/repos/njsmith/pycam02ucs/badge.png?branch=master
   :target: https://coveralls.io/r/njsmith/pycam02ucs?branch=master

Documentation:
  TODO

Installation:
  ``python setup.py install`` (requires a C compiler)

Downloads:
  TODO

Code and bug tracker:
  https://github.com/njsmith/pycam02ucs

Contact:
  Nathaniel J. Smith <njs@pobox.com>

Developer dependencies (only needed for hacking on source):
  * nose: needed to run tests

License:
  MIT, see LICENSE.txt for details.

Other Python packages with similar functionality that you might also
like to consider:
  * ``colour``: http://colour-science.org/
  * ``colormath``: http://python-colormath.readthedocs.org/
  * ``ciecam02``: https://pypi.python.org/pypi/ciecam02/
  * ``ColorPy``: http://markkness.net/colorpy/ColorPy.html
