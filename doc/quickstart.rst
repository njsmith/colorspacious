Quickstart
==========

.. currentmodule:: colorspacious

Here's some cut-and-pasteable examples to give an idea what
:mod:`colorspacious` is good for.

First we need to import it. The main function is
:func:`cspace_convert`:

.. ipython:: python

   from colorspacious import cspace_convert

This allows us to convert between many color spaces. For example,
suppose we want to know how the color with coordinates (128, 128, 128)
in `sRGB <https://en.wikipedia.org/wiki/SRGB>`_ space (represented
with values between 0 and 255) maps to `XYZ
<https://en.wikipedia.org/wiki/CIE_1931_color_space>`_ space
(represented with values between 0 and 100):

.. ipython:: python

   cspace_convert([128, 128, 128], "sRGB255", "XYZ100")

We can also conveniently work on whole images. Let's load one up as an
example.

.. ipython:: python

   import matplotlib.pyplot as plt
   from matplotlib.cbook import get_sample_data
   hopper_sRGB = plt.imread(get_sample_data("grace_hopper.png"))

This image appears to have been loaded as a 3-dimensional NumPy array,
where the last dimension contains the R, G, and B values.

.. ipython:: python

   hopper_sRGB.shape
   hopper_sRGB[:2, :2, :]
   @savefig hopper_sRGB.png width=4in
   plt.imshow(hopper_sRGB)

We can pass such an array directly to :func:`cspace_convert`, using
the ``"sRGB1"`` space now because the values appear to be encoded on a
scale ranging from 0-1. Let's calculate the corresponding lightness,
chroma, and hue, using the `CIECAM02
<https://en.wikipedia.org/wiki/CIECAM02>`_ model, partially desaturate
the image by reducing the chroma, and then convert back to sRGB to
look at the result. (Note that the CIECAM02 model in general requires
the specification of a number of viewing condition parameters; here we
accept the default, which happens to match the viewing conditions
specified in the sRGB standard):

.. ipython:: python

   hopper_desat_JCh = cspace_convert(hopper_sRGB, "sRGB1", "JCh")
   hopper_desat_JCh[..., 1] /= 2
   @savefig hopper_desaturated.png width=4in
   plt.imshow(cspace_convert(hopper_desat_JCh, "JCh", "sRGB1"))

Notice how the colors are more muted, but not entirely gone. Of course
we could also go all the way, for a highly accurate greyscale
conversion:

.. ipython:: python

   hopper_greyscale_JCh = cspace_convert(hopper_sRGB, "sRGB1", "JCh")
   hopper_greyscale_JCh[..., 1] = 0
   hopper_greyscale_sRGB = cspace_convert(hopper_greyscale_JCh, "JCh", "sRGB1")
   @savefig hopper_greyscale_unclipped.png width=4in
   plt.imshow(hopper_greyscale_sRGB)

But notice the small cyan patches on her collar and hat, though --
this occurs due to floating point rounding error creating a few points
with sRGB values that are greater than 1, which causes matplotlib to
render the points in a strange way:

.. ipython:: python

   hopper_greyscale_sRGB[np.any(hopper_greyscale_sRGB > 1, axis=-1), :]

Colorspacious doesn't do anything to clip such values, since they can
sometimes be useful for further processing -- e.g. when chaining
multiple conversions together, you don't want to clip between
intermediate steps, because this might introduce errors. But you can
easily clip them yourself for proper display:

.. ipython:: python

   @savefig hopper_greyscale_clipped.png width=4in
   plt.imshow(np.clip(hopper_greyscale_sRGB, 0, 1))

No more cyan splotches!

We can also simulate various sorts of `color vision deficiency,
a.k.a. "colorblindness"
<https://en.wikipedia.org/wiki/Color_blindness>`_. For example, ~5% of
white men have some degree of deuteranomaly. Here's a simulation of
what this image looks like to someone with a moderate degree of this
condition. Notice the use of the extended syntax for describing color
spaces that require extra parameters beyond just the name:

.. ipython:: python

   cvd_space = {"name": "sRGB1+CVD",
                "cvd_type": "deuteranomaly",
                "severity": 50}
   hopper_deuteranomaly_sRGB = cspace_convert(hopper_sRGB, cvd_space, "sRGB1")

   @savefig hopper_deuteranomaly.png width=4in
   plt.imshow(np.clip(hopper_deuteranomaly_sRGB, 0, 1))

Notice that contrary to what you might expect, we simulate CVD by
asking :func:`cspace_convert` to convert *from* a special CVD space
*to* the standard sRGB space. The way to think about this is, we want
to know what these specific RGB values -- as shown on an sRGB monitor
and viewed by someone with CVD -- would look like if they were instead
shown on an sRGB monitor and viewed by someone with normal color
vision. This may make more sense if you consider the way it makes it
easy to perform other operations, like query the lightness/chroma/hue
of a particular set of RGB coordinates as viewed by someone with CVD
as opposed to someone with normal color vision:

.. ipython:: python

   cspace_convert([1, 0, 0], cvd_space, "JCh")
   cspace_convert([1, 0, 0], "sRGB1", "JCh")

The model of CVD we use allows a "severity" scaling factor, specified
as a number between 0 and 100. A severity of 100 corresponds to
complete dichromacy:

.. ipython:: python

   cvd_space = {"name": "sRGB1+CVD",
                "cvd_type": "deuteranomaly",
                "severity": 100}
   hopper_deuteranopia_sRGB = cspace_convert(hopper_sRGB, cvd_space, "sRGB1")

   @savefig hopper_deuteranopia.png width=4in
   plt.imshow(np.clip(hopper_deuteranopia_sRGB, 0, 1))
