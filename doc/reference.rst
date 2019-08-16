Reference
=========

.. currentmodule:: colorspacious

Conversion functions
--------------------

.. autofunction:: cspace_convert
.. autofunction:: cspace_converter

.. _supported-colorspaces:

Specifying colorspaces
----------------------

Colorspacious knows about a wide variety of colorspaces, some of which
take additional parameters, and it can convert freely between any of
them. Here's an image showing all the known spaces, and the conversion
paths used. (This graph is generated directly from the source code:
when you request a conversion between two spaces,
:func:`cspace_convert` automatically traverses this graph to find the
best conversion path. This makes it very easy to add support for new
colorspaces.)

.. ipython:: python
   :suppress:

   import colorspacious
   with open("_static/colorspacious-graph.dot", "w") as f:
       colorspacious.conversion.GRAPH.dump_dot(f)
   import subprocess
   subprocess.check_call(["dot", "-Tsvg", "_static/colorspacious-graph.dot",
                          "-o", "_static/colorspacious-graph.svg"])

.. image:: /_static/colorspacious-graph.svg

The most general and primitive way to specify a colorspace is via a
dict, e.g., all the following are valid arguments that can be passed
to :func:`cspace_convert`::

   {"name": "XYZ100"}
   {"name": "CIELab", "XYZ100_w": "D65"}
   {"name": "CIELab", "XYZ100_w": [95.047, 100, 108.883]}

These dictionaries always have a ``"name"`` key specifying the
colorspace. Every bold-faced string in the above image is a recognized
colorspace name. Some spaces take additional parameters beyond the
name, such as the CIELab whitepoint above. These additional parameters
are indicated by the italicized strings in the image above.

There are also several shorthands accepted, to let you avoid writing
out long dicts in most cases. In particular:

* Any :class:`CIECAM02Space` object ``myspace`` is expanded to::

    {"name": "CIECAM02",
     "ciecam02_space": myspace}

* Any :class:`LuoEtAl2006UniformSpace` object ``myspace`` is expanded
  to::

    {"name": "J'a'b'",
     "ciecam02_space": CIECAM02.sRGB,
     "luoetal2006_space": myspace}

* The string ``"CIELab"`` expands to: ``{"name": "CIELab", "XYZ100_w": "D65"}``
* The string ``"CIELCh"`` expands to: ``{"name": "CIELCh", "XYZ100_w": "D65"}``
* the string ``"CIECAM02"`` expands to ``CIECAM02Space.sRGB``, which in turn
  expands to ``{"name": "CIECAM02", "ciecam02_space":
  CIECAM02Space.sRGB}``.
* The strings ``"CAM02-UCS"``, ``"CAM02-SCD"``, ``"CAM02-LCD"`` expand
  to the global instance objects :data:`CAM02UCS`, :data:`CAM02SCD`,
  :data:`CAM02LCD`, which in turn expand to ``"J'a'b'"`` dicts as
  described above.
* Any string consisting only of characters from the set "JChQMsH" is
  expanded to::

    {"name": "CIECAM02-subset",
     "axes": <the string provided>
     "ciecam02_space": CIECAM02.sRGB}

  This allows you to directly use common shorthands like ``"JCh"`` or
  ``"JMh"`` as first-class colorspaces.

Any other string ``"foo"`` expands to ``{"name": "foo"}``. So for any
space that doesn't take parameters, you can simply say ``"sRGB1"`` or
``"XYZ100"`` or whatever and ignore all these complications.

And, as one final trick, any alias can also be used as the ``"name"``
field in a colorspace dict, in which case its normal expansion is
used to provide overrideable defaults for parameters. For example::

  # You write:
  {"name": "CAM02-UCS",
   "ciecam02_space": my_ciecam02_space}

  # Colorspacious expands this to:
  {"name": "J'a'b'",
   "ciecam02_space": my_ciecam02_space,
   "luoetal2006_space": CAM02UCS}

Or::

  # You write:
  {"name": "JCh",
   "ciecam02_space": my_ciecam02_space}

  # Colorspacious expands this to:
  {"name": "CIECAM02-subset",
   "axes": "JCh",
   "ciecam02_space": my_ciecam02_space}


Well-known colorspaces
......................

**sRGB1**, **sRGB100**: The standard `sRGB colorspace
<https://en.wikipedia.org/wiki/SRGB>`_. If you have generic "RGB"
values with no further information specified, then usually the right
thing to do is to assume that they are in the sRGB space; the sRGB
space was originally designed to match the behavior of common consumer
monitors, and these days common consumer monitors are designed to
match sRGB. Use ``sRGB1`` if you have or want values that are
normalized to fall between 0 and 1, and use ``sRGB255`` if you have or
want values that are normalized to fall between 0 and 255.

**XYZ100**, **XYZ1**: The standard `CIE 1931 XYZ color space
<https://en.wikipedia.org/wiki/CIE_1931_color_space>`_. Use ``XYZ100``
if you have or want values that are normalized to fall between 0 and
100 (roughly speaking -- values greater than 100 are valid in certain
cases). Use ``XYZ1`` if you have or want values that are normalized to
fall between 0 and 1 (roughly). This is a space which is
"linear-light", i.e. related by a linear transformation to the photon
counts in a spectral power distribution. In particular, this means
that linear interpolation in this space is a valid way to simulate
physical mixing of lights.

**sRGB1-linear**: A linear-light version of **sRGB1**, i.e., it has
had gamma correction applied, but is still represented in terms of the
standard sRGB primaries.

**xyY100**, **xyY1**: The standard `CIE 1931 xyY color space
<https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_xy_chromaticity_diagram_and_the_CIE_xyY_color_space>`_. *The
x and y values are always normalized to fall between 0 and 1.* Use
``xyY100`` if you have or want a Y value that falls between 0 and 100,
and use ``xyY1`` if you have or want a Y value that falls between 0
and 1.

**CIELab**: The standard `CIE 1976 L*a*b* color space
<https://en.wikipedia.org/wiki/Lab_color_space>`_. L* is scaled to
vary from 0 to 100; This color space is unbounded. 
For colors inside the sRGB gamut, the values of a* and b* will be roughly in the range [-100, 100].
This space takes a parameter, *XYZ100_w*, which sets
the reference white point, and may be specified either directly as a
tristimulus value or as a string naming one of the well-known standard
illuminants like ``"D65"``.

**CIELCh**: Cylindrical version of **CIELab**. Accepts the same
parameters. h* is in degrees.


Simulation of color vision deficiency
.....................................

We provide simulation of common (and not so common) forms of color
vision deficiency (also known as "colorblindness"), using the model
described by :cite:`Machado-CVD`.

This is generally done by specifying a colorspace like::

  {"name": "sRGB1+CVD",
   "cvd_type": <type>,
   "severity": <severity>}

where ``<type>`` is one of the following strings:

* ``"protanomaly"``: A common form of red-green colorblindness;
  affects ~2% of white men to some degree (less common among other
  ethnicities, much less common among women, see Tables 1.5 and 1.6 in
  :cite:`Sharpe-CVD`).
* ``"deuteranomaly"``: The most common form of red-green
  colorblindness; affects ~6% of white men to some degree (less common
  among other ethnicities, much less common among women, see Tables
  1.5 and 1.6 in :cite:`Sharpe-CVD`).
* ``"tritanomaly"``: A very rare form of colorblindness affecting
  blue/yellow discrimination -- so rare that its detailed effects and
  even rate of occurrence are not well understood. Affects <0.1% of
  people, possibly much less (:cite:`Sharpe-CVD`, page 47). Also, the
  name we use here is somewhat misleading because only full
  trit\ **anopia** has been documented, and partial trit\ **anomaly**
  likely does not exist (:cite:`Sharpe-CVD`, page 45). What this means
  is that while Colorspacious will happily allow any severity value to
  be passed, probably only severity = 100 corresponds to any real
  people.

And ``<severity>`` is any number between 0 (indicating regular vision)
and 100 (indicating complete dichromacy).

.. warning:: If you have an image, e.g. a photo, and you want to
   "convert it to simulate colorblindness", then this is done with an
   incantation like::

     cspace_convert(img, some_cvd_space, "sRGB1")

   Notice that these arguments are given in the *opposite order* from
   what you might naively expect. See :ref:`tutorial-cvd` for
   explanation and worked examples.


CIECAM02
........

`CIECAM02 <https://en.wikipedia.org/wiki/CIECAM02>`_ is a
standardized, rather complex, state-of-the-art color appearance model,
i.e., it's not useful for describing the voltage that should be
applied to a phosphorescent element in your monitor (like RGB was
originally designed to do), and it's not useful for modelling physical
properties of light (like XYZ), but it is very useful to tell you what
a color will look like subjectively to a human observer, under a
certain set of viewing conditions. Unfortunately this makes it rather
complicated, because human vision is rather complicated.

If you just want a better replacement for traditional ad hoc spaces
like "Hue/Saturation/Value", then use the string ``"JCh"`` for your
colorspace (see :ref:`tutorial-perception` for a tutorial) and be
happy.

If you want the full power of CIECAM02, or just to understand what
*exactly* is happening when you type ``"JCh"``, then read on.

First, you need to specify your viewing conditions. For many purposes,
you can use the default :attr:`CIECAM02Space.sRGB`
object. Or if you want to specify different viewing conditions, you
can instantiate your own :class:`CIECAM02Space` object:

.. autoclass:: CIECAM02Space

   .. attribute:: sRGB

      A class-level constant representing the viewing conditions
      specified in the sRGB standard. (The sRGB standard defines two
      things: how a standard monitor should respond to different RGB
      values, and a standard set of viewing conditions in which you
      are supposed to look at such a monitor, and that attempt to
      approximate the average conditions in which people actually do
      look at such monitors. This object encodes the latter.)

   The CIECAM02Space object has some low-level methods you can use
   directly if you want, though usually it'll be easier to just use
   :func:`cspace_convert`:

   .. automethod:: XYZ100_to_CIECAM02
   .. automethod:: CIECAM02_to_XYZ100

.. autoclass:: CIECAM02Surround

   A namedtuple holding the CIECAM02 surround parameters, :math:`F`,
   :math:`c`, and :math:`N_c`.

   The CIECAM02 standard surrounds are available as constants defined
   on this class; for most purposes you'll just want to use one of
   them:

   * :data:`CIECAM02Surround.AVERAGE`
   * :data:`CIECAM02Surround.DIM`
   * :data:`CIECAM02Surround.DARK`

.. autoclass:: NegativeAError

Now that you have a :class:`CIECAM02Space` object, what can you do
with it?

First, you can pass it directly to :func:`cspace_convert` as an input
or output space (which is a shorthand for using a space like
``{"name": "CIECAM02", "ciecam02_space": <whatever>}``).

The plain vanilla ``"CIECAM02"`` space is weird and special: unlike
all the other spaces supported by colorspacious, it does not represent
values with ordinary NumPy arrays. This is because there are just too
many perceptual correlates, and trying to keep track of whether M is
at index 4 or 5 would be way too obnoxious. Instead, it returns an
object of class :class:`JChQMsH`:

.. autoclass:: JChQMsH

   A namedtuple with a mnemonic name: it has attributes ``J``, ``C``,
   ``h``, ``Q``, ``M``, ``s``, and ``H``, each of which holds a scalar
   or NumPy array representing lightness, chroma, hue angle,
   brightness, colorfulness, saturation, and hue composition,
   respectively.

Alternatively, because you usually only want a subset of these, you
can take advantage of the ``"CIECAM02-subset"`` space, which takes the
perceptual correlates you want as a parameter. So for example if you
just want JCh, you can write::

  {"name": "CIECAM02-subset",
   "axes": "JCh",
   "ciecam02_space": CIECAM02.sRGB}

When using ``"CIECAM02-subset"``, you don't have to worry about
:class:`JChQMsH` -- it just takes and returns regular NumPy arrays,
like all the other colorspaces.

And as a convenience, all strings composed of the character JChQMsH
are automatically treated as specifying CIECAM02-subset spaces, so you
can write::

  "JCh"

and it expands to::

  {"name": "CIECAM02-subset",
   "axes": "JCh",
   "ciecam02_space": CIECAM02.sRGB}

or you can write::

  {"name": "JCh",
   "ciecam02_space": my_space}

and it expands to::

  {"name": "CIECAM02-subset",
   "axes": "JCh",
   "ciecam02_space": my_space}


Perceptually uniform colorspaces based on CIECAM02
..................................................

The :math:`J'a'b'` spaces proposed by :cite:`CAM02-UCS` are
high-quality, approximately perceptually uniform spaces based on
CIECAM02. They propose three variants: CAM02-LCD optimized for "large
color differences" (e.g., estimating the similarity between blue and
green), CAM02-SCD optimized for "small color differences" (e.g.,
estimating the similarity between light blue with a faint greenish
cast and light blue with a faint purpleish cast), and CAM02-UCS which
attempts to provide a single "uniform color space" that is less
optimized for either case but provides acceptable performance in
general.

Colorspacious represents these spaces as instances of
:class:`LuoEtAl2006UniformSpace`:

.. autoclass:: LuoEtAl2006UniformSpace

Because these spaces are defined as transformations from CIECAM02, to
have a fully specified color space you must also provide some
particular CIECAM02 viewing conditions, e.g.::

  {"name": "J'a'b'",
   "ciecam02_space": CIECAM02.sRGB,
   "luoetal2006_space": CAM02UCS}

As usual, you can also pass any instance of
:class:`LuoEtAl2006UniformSpace` and it will be expanded into a dict
like the above, or for the three common variants you can pass the
strings ``"CAM02-UCS"``, ``"CAM02-LCD"``, or ``"CAM02-SCD"``.

.. versionchanged:: 1.1.0

   In v1.0.0 and earlier, colorspacious's definitions of the
   ``CAM02-LCD`` and ``CAM02-SCD`` spaces were swapped compared to
   what they should have been based on the :cite:`CAM02-UCS` – i.e.,
   if you asked for LCD, you got SCD, and vice-versa. (``CAM02-UCS``
   was correct, though). Starting in 1.1.0, all three spaces are now
   correct.

Color difference computation
----------------------------

.. autofunction:: deltaE

For examples, see :ref:`tutorial-deltaE` in the tutorial.

Utilities
---------

You probably won't need these, but just in case they're useful:

.. autofunction:: standard_illuminant_XYZ100
.. autofunction:: as_XYZ100_w

.. autofunction:: machado_et_al_2009_matrix
