Tutorial
========

.. currentmodule:: colorspacious

Colorspacious is a Python library that lets you easily convert between
colorspaces like sRGB, XYZ, CIEL*a*b*, CIECAM02, CAM02-UCS, etc. If
you have no idea what these are or what each is good for, and reading
this list makes you feel like you're swimming in alphabet soup, then
this video provides a basic orientation and some examples. (The
overview of color theory starts at ~3:35.)

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/xAoljeRJ3lU" frameborder="0" allowfullscreen></iframe>

Now let's see some cut-and-pasteable examples of what
:mod:`colorspacious` is good for. We'll start by loading up some
utility modules for numerics and plotting that we'll use later:

.. ipython:: python

   import numpy as np
   import matplotlib
   import matplotlib.pyplot as plt

Now we need to import :mod:`colorspacious`. The main function we'll
use is :func:`cspace_convert`:

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

Colorspacious knows about a :ref:`wide variety of colorspaces
<supported-colorspaces>`, and you can convert between any of them by
naming them in a call to :func:`cspace_convert`.

We can also conveniently work on whole images. Let's load one up as an
example.

.. ipython:: python

   # if you want this file, try:
   #    hopper_sRGB = plt.imread(matplotlib.cbook.get_sample_data("grace_hopper.png"))
   hopper_sRGB = plt.imread("grace_hopper.png")

What have we got here?

.. ipython:: python

   hopper_sRGB.shape
   hopper_sRGB[:2, :2, :]
   @savefig hopper_sRGB.png width=4in
   plt.imshow(hopper_sRGB)

It looks like this image has been loaded as a 3-dimensional NumPy
array, where the last dimension contains the R, G, and B values (in
that order).

We can pass such an array directly to :func:`cspace_convert`. For
example, we can convert the whole image to XYZ space. This time we'll
specify that our input space is ``"sRGB1"`` instead of ``"sRGB255"``,
because the values appear to be encoded on a scale ranging from 0-1:

.. ipython:: python

   hopper_XYZ = cspace_convert(hopper_sRGB, "sRGB1", "XYZ100")
   hopper_XYZ.shape
   hopper_XYZ[:2, :2, :]


.. _tutorial-perception:

Perceptual transformations
--------------------------

RGB space is a useful way to store and transmit images, but because
the RGB values are basically just a raw record of what voltages should
be applied to some phosphors in a monitor, it's often difficult to
predict how a given change in RGB values will affect what an image
looks like to a person.

Suppose we want to desaturate an image -- that is, we want to replace
each color by a new color that has the same lightness (so white stays
white, black stays black, etc.), and the same hue (so each shade of
blue stays the same shade of blue, rather than turning into purple or
red), but the "chroma" is reduced (so colors are more muted). This is
very difficult to do when working in RGB space. So let's take our
colors and re-represent them in terms of lightness, chroma, and hue,
using the state-of-the-art `CIECAM02
<https://en.wikipedia.org/wiki/CIECAM02>`_ model.

The three axes in this space are conventionally called "J" (for
lightness), "C" (for chroma), and "h" (for hue). (The CIECAM-02
standard also defines a whole set of other axes with subtly different
meanings -- `see Wikipedia for details
<https://en.wikipedia.org/wiki/CIECAM02>`_ -- but for now we'll stick
to these three.) To desaturate our image, we're going to switch from
sRGB space to JCh space, reduce all the "C" values by a factor of 2,
and then convert back to sRGB to look at the result. (Note that the
CIECAM02 model in general requires the specification of a number of
viewing condition parameters; here we accept the default, which happen
to match the viewing conditions specified in the sRGB standard). All
this takes more words to describe than it does to implement:

.. ipython:: python

   hopper_desat_JCh = cspace_convert(hopper_sRGB, "sRGB1", "JCh")
   # This is in "JCh" space, and we want to modify the "C" channel, so
   # that's channel 1.
   hopper_desat_JCh[..., 1] /= 2
   hopper_desat_sRGB = cspace_convert(hopper_desat_JCh, "JCh", "sRGB1")

Let's see what this looks like. First we'll define a little utility
function to plot several images together:

.. ipython:: python

   def compare_hoppers(*new):
       image_width = 2.0  # inches
       total_width = (1 + len(new)) * image_width
       height = image_width / hopper_sRGB.shape[1] * hopper_sRGB.shape[0]
       fig = plt.figure(figsize=(total_width, height))
       ax = fig.add_axes((0, 0, 1, 1))
       ax.imshow(np.column_stack((hopper_sRGB,) + new))

And now we'll use it to look at the desaturated image we computed above:

.. ipython:: python

   @savefig hopper_desaturated.png width=6in
   compare_hoppers(hopper_desat_sRGB)

The original version is on the left, with our modified version on the
right. Notice how in the version with reduced chroma, the colors are
more muted, but not entirely gone. Of course we could also reduce the
chroma all the way to zero, for a highly accurate greyscale
conversion:

.. ipython:: python

   hopper_greyscale_JCh = cspace_convert(hopper_sRGB, "sRGB1", "JCh")
   hopper_greyscale_JCh[..., 1] = 0
   hopper_greyscale_sRGB = cspace_convert(hopper_greyscale_JCh, "JCh", "sRGB1")
   @savefig hopper_greyscale_unclipped.png width=6in
   compare_hoppers(hopper_greyscale_sRGB)

But notice the small cyan patches on her collar and hat --
this occurs due to floating point rounding error creating a few points
with sRGB values that are greater than 1, which causes matplotlib to
render the points in a strange way:

.. ipython:: python

   hopper_greyscale_sRGB[np.any(hopper_greyscale_sRGB > 1, axis=-1), :]

Colorspacious doesn't do anything to clip such values, since they can
sometimes be useful for further processing -- e.g. when chaining
multiple conversions together, you don't want to clip between
intermediate steps, because this might introduce errors. And
potentially you might want to handle them in some clever way
(e.g. rescaling your whole image). But in this case, where the values
are only just barely over 1, then simply clipping them to 1 is
probably the best approach, and you can easily do this yourself:

.. ipython:: python

   @savefig hopper_greyscale_clipped.png width=6in
   compare_hoppers(np.clip(hopper_greyscale_sRGB, 0, 1))

No more cyan splotches!

To explore, try applying other transformations. E.g., you could darken
the image by rescaling the lightness channel "J" by a factor of 2
(``image_JCh[..., 0] /= 2``), or try replacing each hue by its
complement (``image_JCh[..., 2] *= -1``).


.. _tutorial-cvd:

Simulating colorblindness
-------------------------

Another useful thing we can do by converting colorspaces is to
simulate various sorts of `color vision deficiency,
a.k.a. "colorblindness"
<https://en.wikipedia.org/wiki/Color_blindness>`_. For example,
deuteranomaly is the name for the most common form of red-green
colorblindness, and affects ~5% of white men to varying
amounts. Here's a simulation of what this image looks like to someone
with a moderate degree of this condition. Notice the use of the
extended syntax for describing color spaces that require extra
parameters beyond just the name:

.. ipython:: python

   cvd_space = {"name": "sRGB1+CVD",
                "cvd_type": "deuteranomaly",
                "severity": 50}
   hopper_deuteranomaly_sRGB = cspace_convert(hopper_sRGB, cvd_space, "sRGB1")

   @savefig hopper_deuteranomaly.png width=6in
   compare_hoppers(np.clip(hopper_deuteranomaly_sRGB, 0, 1))

Notice that contrary to what you might expect, we simulate CVD by
asking :func:`cspace_convert` to convert *from* a special CVD space
*to* the standard sRGB space. The way to think about this is that we
have a set of RGB values that will be viewed under certain conditions,
i.e. displayed on an sRGB monitor and viewed by someone with CVD. And
we want to find a new set of RGB values that will look the same under
a different set of viewing conditions, i.e., displayed on an sRGB
monitor and viewed by someone with normal color vision. So we are
starting in the ``sRGB1+CVD`` space, and converting to the normal
``sRGB1`` space.

This way of doing things is especially handy when you want to
perform other operations. For example, we might want to use the JCh
space described above to ask "what (approximate) lightness/chroma/hue
would someone with this form of CVD perceive when looking at a monitor
displaying a certain RGB value?". For example, taking a "pure red" color:

.. ipython:: python

   cspace_convert([1, 0, 0], cvd_space, "JCh")

If we compare this to someone with normal color vision, we see that
the person with CVD will perceive about the same lightness, but
desaturated and with a shifted hue:

.. ipython:: python

   cspace_convert([1, 0, 0], "sRGB1", "JCh")

The model of CVD we use allows a "severity" scaling factor, specified
as a number between 0 and 100. A severity of 100 corresponds to
complete dichromacy:

.. ipython:: python

   cvd_space = {"name": "sRGB1+CVD",
                "cvd_type": "deuteranomaly",
                "severity": 100}
   hopper_deuteranopia_sRGB = cspace_convert(hopper_sRGB, cvd_space, "sRGB1")

   @savefig hopper_deuteranopia.png width=8in
   compare_hoppers(np.clip(hopper_deuteranomaly_sRGB, 0, 1),
                   np.clip(hopper_deuteranopia_sRGB, 0, 1))

Here the leftmost and center images are repeats of ones we've seen
before: the leftmost image is the original, and the center image is
the moderate deuteranomaly simulation that we computed above. The
image on the right is the new image illustrating the more severe
degree of red-green colorblindness -- notice how the red in the flag
and her medals is muted in the middle image, but in the image on the
right it's disappeared completely.

You can also set the ``"cvd_type"`` to ``"protanomaly"`` to simulate
the other common form of red-green colorblindness, or to
``"tritanomaly"`` to simulate an extremely rare form of blue-yellow
colorblindness. Here's what moderate and severe protanomaly look like
when simulated by colorspacious:

.. ipython:: python

   cvd_space = {"name": "sRGB1+CVD",
                "cvd_type": "protanomaly",
                "severity": 50}
   hopper_protanomaly_sRGB = cspace_convert(hopper_sRGB, cvd_space, "sRGB1")

   cvd_space = {"name": "sRGB1+CVD",
                "cvd_type": "protanomaly",
                "severity": 100}
   hopper_protanopia_sRGB = cspace_convert(hopper_sRGB, cvd_space, "sRGB1")

   @savefig hopper_protanopia.png width=8in
   compare_hoppers(np.clip(hopper_protanomaly_sRGB, 0, 1),
                   np.clip(hopper_protanopia_sRGB, 0, 1))

Because deuteranomaly and protanomaly are both types of red-green
colorblindness, this is similar (but not quite identical) to the image
we saw above.


.. _tutorial-deltaE:

Color similarity
----------------

Suppose we have two colors, and we want to know how different they
will look to a person -- often known as computing the `"delta E"
<http://www.colorwiki.com/wiki/Delta_E:_The_Color_Difference>`_
between them. One way to do this is to map both colors into a
"perceptually uniform" colorspace, and then compute the Euclidean
distance. Colorspacious provides a convenience function to do just
this:

.. ipython:: python

   from colorspacious import deltaE

   deltaE([1, 0.5, 0.5], [0.5, 1, 0.5])

   deltaE([255, 127, 127], [127, 255, 127], input_space="sRGB255")

By default, these computations are done using the CAM02-UCS
perceptually uniform space (see :cite:`CAM02-UCS` for details), but if
you want to use the (generally inferior) `CIEL*a*b*
<https://en.wikipedia.org/wiki/Lab_color_space>`_, then just say the
word:

.. ipython:: python

   deltaE([1, 0.5, 0.5], [0.5, 1, 0.5], uniform_space="CIELab")
