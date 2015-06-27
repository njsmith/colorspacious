# This file is part of pycam02ucs
# Copyright (C) 2014-2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

from .convert import convert
from .cam02ucs import CAM02_UCS
from .ciecam02 import ViewingConditions

def deltaEp(color1, color2,
            source_space="sRGB", uniform_space=CAM02_UCS,
            viewing_conditions=ViewingConditions.sRGB):
    """Computes the :math:`\delta E'` distance between pairs of colors.

    :math:`\delta E'` is color difference metric defined by Eq. (4) of Luo et
    al (2006); they show that it provides a good match to human color
    similarity judgements.

    Valid modes are "LCD" (which Luo et al tuned on their "large color
    difference" judgement data sets), "SCD" (which Luo et all tuned on their
    "small color difference" data sets), and "UCS" (which attempts define a
    single generic "uniform color space" which performs well across all their
    data sets). You can also pass in any LuoUniformSpace object as the mode.

    This function is vectorized, i.e., color1, color2 may be arrays with shape
    (..., 3), in which case we compute the distance between corresponding
    pairs of colors.

    :param source_space: The space the colors start out in. Can be anything
       recognized by :func:`convert`.
    :param cam02: Which :class:`CAM02` space to use. Default: CAM02.UCS.
    :param viewing_conditions: The :class:`ViewingConditions` to assume for
       CIECAM02 transformations.
    """

    JpKapbp1 = convert(color1, source_space, "J'/Ka'b'",
                       viewing_conditions=viewing_conditions,
                       cam02=cam02)

    JpKapbp2 = convert(color2, source_space, "J'/Ka'b'",
                       viewing_conditions=viewing_conditions,
                       cam02=cam02)

    return np.sqrt(np.sum((JpKapbp1 - JpKapbp2) ** 2, axis=-1))
