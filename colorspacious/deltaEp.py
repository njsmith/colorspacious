# This file is part of colorspacious
# Copyright (C) 2014-2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

from .conversion import cspace_convert

def deltaEp(color1, color2,
            source_space="sRGB1", uniform_space="CAM02-UCS"):
    """Computes the :math:`\delta E` distance between pairs of colors.

    This function is vectorized, i.e., color1, color2 may be arrays with shape
    (..., 3), in which case we compute the distance between corresponding
    pairs of colors.

    :param source_space: The space the colors start out in. Can be anything
       recognized by :func:`cspace_convert`. Default: "sRGB1"
    :param uniform_space: Which space to perform the distance measurement
       in. This should be a uniform space like CAM02-UCS where
       Euclidean distance approximates similarity judgements, because
       otherwise the results of this function won't be very meaningful, but in
       fact any color space known to :func:`cspace_convert` will be accepted.
    """

    uniform1 = cspace_convert(color1, source_space, uniform_space)

    uniform2 = cspace_convert(color2, source_space, uniform_space)

    return np.sqrt(np.sum((uniform1 - uniform2) ** 2, axis=-1))
