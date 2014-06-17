# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np

from .srgb import sRGB_to_XYZ
from ._ciecam02 import XYZ_to_JMh

# Equations (4) from Luo et al (2006)
def JMh_to_Jpapbp(J, M, h, c1, c2):
    Jp = (1 + 100 * c1) * J / (1 + c1 * J)
    Mp = (1. / c2) * np.log(1 + c2 * M)
    h_rad = np.deg2rad(h)
    ap = Mp * np.cos(h_rad)
    bp = Mp * np.sin(h_rad)
    return Jp, ap, bp

def deltaEp_JMh(J1, M1, h1, J2, M2, h2, KL, c1, c2):
    Jp1, ap1, bp1 = JMh_to_Jpapbp(J1, M1, h1, c1, c2)
    Jp2, ap2, bp2 = JMh_to_Jpapbp(J2, M2, h2, c1, c2)

    return np.sqrt(
        ((Jp1 - Jp2) / KL) ** 2
        + (ap1 - ap2) ** 2
        + (bp1 - bp2) ** 2
        )

def sRGB_to_JMh(R, G, B):
    return XYZ_to_JMh(*sRGB_to_XYZ(R, G, B))

def deltaEp_sRGB(R1, G1, B1, R2, G2, B2, mode="UCS"):
    """Computes the :math:`\delta E'` distance between pairs of sRGB colors.

    :math:`\delta E'` is color difference metric defined by Eq. (4) of Luo et
    al (2006); they show that it provides a good match to human color
    similarity judgements.

    Valid modes are "LCD" (which Luo et al tuned on their "large color
    difference" judgement data sets), "SCD" (which Luo et all tuned on their
    "small color difference" data sets), and "UCS" (which attempts define a
    single generic "uniform color space" which performs well across all their
    data sets).

    This function is vectorized, i.e., R1, G1, B1, R2, G2, B2 may be vectors,
    in which case we compute the distance between corresponding pairs of
    colors.
    """

    # Table II in Luo et al (2006)
    c1 = 0.007
    if mode == "LCD":
        KL = 0.77
        c2 = 0.0053
    elif mode == "SCD":
        KL = 1.24
        c2 = 0.0363
    elif mode == "UCS":
        KL = 1.00
        c2 = 0.0228
    else:
        raise ValueError("mode must be 'UCS', 'LCD', or 'SCD'")

    J1, M1, h1 = sRGB_to_JMh(R1, G1, B1)
    J2, M2, h2 = sRGB_to_JMh(R2, G2, B2)
    return deltaE_JMh(J1, M1, h1, J2, M2, h2, KL, c1, c2)
