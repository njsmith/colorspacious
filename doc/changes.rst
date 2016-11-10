Changes
=======

v1.1.0 (not yet released)
-------------------------

* **BUG AFFECTING CALCULATIONS:** In previous versions, it turns out
  that the CAM02-LCD and CAM02-SCD spaces were accidentally swapped –
  so if you asked for CAM02-LCD you got SCD, and vice-versa. This has
  now been corrected. (Thanks to Github user TFiFiE for catching
  this!)

* Fixed setup.py to be compatible with both python 2 and python 3.

* Miscellaneous documentation improvements.


v1.0.0
------

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.33086.svg
   :target: http://dx.doi.org/10.5281/zenodo.33086

Notable changes since v0.1.0 include:

* **BUG AFFECTING CALCULATIONS:** the sRGB viewing conditions
  (``colorspacious.CIECAM02Space.sRGB``), which are used by default in
  all calculations involving CIECAM02 or CAM02-UCS, were previously
  incorrect -- the :math:`L_A` parameter was supposed to be :math:`(64
  / \pi) / 5`, but instead was incorrectly calculated as :math:`(64 /
  \pi) * 5`. The effect of this was to assume much brighter ambient
  lighting than actually specified by the sRGB standard (i.e., the
  sRGB standard assumes that you are looking at your monitor in a dim
  environment, like a movie theatre; we were calculating as if you
  were looking at your monitor in an environment that was 125 times
  lighter -- something like, outside on an overcast day). This bug is
  corrected in this release.

  Fortunately this turns out to have had a negligible effect on
  viridis and the other matplotlib colormaps that were computed using
  the buggy code. Once the bug is corrected, the old colormaps'
  perceptual uniformity is no long analytically exactly perfect, but
  the deviations are numerically negligible, so there's no need to
  regenerate the colormaps. (Indeed, the buggy viewing conditions,
  while different from those specified in IEC 61966-2-1:1999, are
  probably still within the range of realistic viewing conditions
  where these colormaps will be used.)

  If it is necessary to reproduce results using the old code, then
  this can be accomplished by instantiating a custom
  :class:`CIECAM02Space` object::

      from colorspacious import CIECAM02Space
      # almost, but not quite, the sRGB viewing conditions:
      buggy_space = CIECAM02Space(
          XYZ100_w="D65",
          Y_b=20,
          # bug: should be (64 / np.pi) / 5
          L_A=(64 / np.pi) * 5)

  This can be used directly, or to create custom colorspace
  specifications to use with :func:`cspace_convert`. E.g., to convert
  from sRGB1 to JCh using the buggy viewing conditions::

      cspace_convert(..., "sRGB1",
                     {"name": "JCh", "ciecam02_space": buggy_space})

  Or to convert from XYZ100 to CAM02-UCS using the buggy viewing
  conditions::

      cspace_convert(..., "XYZ100",
                     {"name": "CAM02-UCS", "ciecam02_space": buggy_space})

  Similar code has been added to `viscm
  <https://github.com/matplotlib/viscm>`_ to allow reproduction and
  editing of viridis and related colormaps that were designed using
  the old code.

* :func:`colorspacious.deltaE` is now available as a convenience
  function for computing the perceptual distance between colors.

* Substantially improved docs (i.e. there is now actually a
  comprehensive manual).

* Better test coverage (currently at 100% statement and branch
  coverage).

* Miscellaneous bug fixes.


v0.1.0
------

Initial release.
