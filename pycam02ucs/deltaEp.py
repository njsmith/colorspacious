# This file is part of pycam02ucs
# Copyright (C) 2014-2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

from .conversion import convert_cspace

def deltaEp(color1, color2,
            source_space="sRGB", uniform_space="CAM02-UCS"):
    """Computes the :math:`\delta E` distance between pairs of colors.

    This function is vectorized, i.e., color1, color2 may be arrays with shape
    (..., 3), in which case we compute the distance between corresponding
    pairs of colors.

    :param source_space: The space the colors start out in. Can be anything
       recognized by :func:`convert_cspace`.
    :param uniform_space: Which space to perform the distance measurement
       in. This should be a uniform space like CAM02-UCS where
       Euclidean distance approximates similarity judgements, because
       otherwise the results of this function won't be very meaningful, but in
       fact any color space known to :func:`convert_cspace` will be accepted.
    """

    JpKapbp1 = convert(color1, source_space, uniform_space)

    JpKapbp2 = convert(color2, source_space, uniform_space)

    return np.sqrt(np.sum((JpKapbp1 - JpKapbp2) ** 2, axis=-1))
