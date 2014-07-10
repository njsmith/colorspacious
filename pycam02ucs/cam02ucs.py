# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np

from .srgb import sRGB_to_XYZ, XYZ_to_sRGB
from ._ciecam02 import XYZ_to_JMh

KL_c1_c2 = {'LCD': (0.77, 0.007, 0.0053),
            'SCD': (1.24, 0.007, 0.0363),
            'UCS': (1.00, 0.007, 0.0228),
             }

def get_KL_c1_c2(mode):
    if mode in KL_c1_c2:
        return KL_c1_c2[mode]
    else:
        raise ValueError("mode must be one of: %s" % " ".join(KL_c1_c2))

def JMh_to_XYZ(J, M, h):
    raise NotImplementedError    


###########  sRGB <=> JMh       ##################
# via XYZ

def sRGB_to_JMh(R, G, B):
    return XYZ_to_JMh(*sRGB_to_XYZ(R, G, B))

def JMh_to_sRGB(J, M, h):
    return sRGB_to_XYZ(*JMh_to_XYZ(J, M, h))


#############  JMh <=> J/K a' b' #################
# Equations (4) from Luo et al (2006)
# J/K a' b' is a Euclidean space corresponding to Luo et al (2006)'s perceptual
# distance metrics. Which metric depends on the values of KL, c1, and c2 used
# for conversion. 

def JMh_to_JKapbp(J, M, h, KL, c1, c2):
    Jp = (1 + 100 * c1) * J / (1 + c1 * J)
    JK = Jp / K
    Mp = (1. / c2) * np.log(1 + c2 * M)
    h_rad = np.deg2rad(h)
    ap = Mp * np.cos(h_rad)
    bp = Mp * np.sin(h_rad)
    return JK, ap, bp

def JKapbp_to_JMh(JK, ap, bp, KL, c1, c2):
    Jp = JK * KL
    J = Jp / (c1 * Jp - 100 * c1 - 1)
    # ap = Mp * cos(h)
    # bp = Mp * sin(h)
    # Mp = bp/sin(h)
    # ap = (bp/sin(h)) * cos(h)
    # ap = bp * sin(h)/cos(h)
    # cos(h) = bp/ap * sin(h)
    # solve numerically for h in the interval of possible hue degrees
    h = TODO
    Mp = bp/sin(h)
    assert np.allclose(Mp, ap/cos(h))
    M = (np.exp(c2*Mp) - 1)/c2
    return J, M, h


###########  sRGB <=> J/K a' b' ##################
# Full pipeline from sRGB to the uniform space J/K a' b'.
# sRGB -> XYZ -> JMh -> J/K a' b'

def sRGB_to_JKapbp(R, G, B, mode='UCS'):
    KL, c1, c2 = get_KL_c1_c2(mode)
    X, Y, Z = sRGB_to_XYZ(R, G, B)
    J, M, h = XYZ_to_JMh(X, Y, Z)
    return JMh_to_JKapbp(J, M, h, KL, c1, c2)

def JKapbp_to_sRGB(JK, ap, bp, mode='UCS'):
    KL, c1, c2 = get_KL_c1_c2(mode)
    J, M, h = JKapbp_to_JMh(JK, ap, bp, KL, c1, c2)
    X, Y, Z = JMh_to_XYZ(J, M, h)
    return XYZ_to_sRGB(X, Y, Z)


########## Similarity functions   ################
# These similarity functions are equivalent to converting the color
# to J/K a' b' and taking Euclidean distance.
    
def deltaEp_JMh(J1, M1, h1, J2, M2, h2, KL, c1, c2):
    JK1, ap1, bp1 = JMh_to_JKapbp(J1, M1, h1, KL, c1, c2)
    JK2, ap2, bp2 = JMh_to_JKapbp(J2, M2, h2, KL, c1, c2)

    return np.sqrt(
        (JK - JK2) ** 2
        + (ap1 - ap2) ** 2
        + (bp1 - bp2) ** 2
        )

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
    KL, c1, c2 = get_KL_c1_c2(mode)
    J1, M1, h1 = sRGB_to_JMh(R1, G1, B1)
    J2, M2, h2 = sRGB_to_JMh(R2, G2, B2)
    return deltaEp_JMh(J1, M1, h1, J2, M2, h2, KL, c1, c2)

