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


#############  JMh <=> J' a' b' #################
# Equations (4) from Luo et al (2006)

def JMh_to_Jpapbp(J, M, h, c1, c2):
    Jp = (1 + 100 * c1) * J / (1 + c1 * J)
    Mp = (1. / c2) * np.log(1 + c2 * M)
    h_rad = np.deg2rad(h)
    ap = Mp * np.cos(h_rad)
    bp = Mp * np.sin(h_rad)
    return Jp, ap, bp

def Jpapbp_to_JMh(Jp, ap, bp, c1, c2):
    raise NotImplementedError


###########  J' a' b' <=> J/K a' b' ##################
# J/K a' b' is a Euclidean space corresponding to Luo et al (2006)'s perceptual
# distance metrics. Which metric depends on the values of KL, c1, and c2 used
# for conversion.

def Jpapbp_to_JKapbp(Jp, ap, bp, KL):
    return Jp/KL, ap, bp

def JKapbp_to_Jpapbp(JK, ap, bp, KL):
    return Jp*KL, ap, bp
    

###########  sRGB <=> J/K a' b' ##################
# Full pipeline from sRGB to the uniform space J/K a' b'.
# sRGB -> XYZ -> JMh -> J' a' b' -> J/K a' b'

def sRGB_to_JKapbp(R, G, B, mode='UCS'):
    KL, c1, c2 = get_KL_c1_c2(mode)
    X, Y, Z = sRGB_to_XYZ(R, G, B)
    J, M, h = XYZ_to_JMh(X, Y, Z)
    Jp, ap, bp = JMh_to_Jpapbp(J, M, h, c1, c2)
    return Jpapbp_to_JKapbp(Jp, ap, bp, KL)

def JKapbp_to_sRGB(JK, ap, bp, mode='UCS'):
    KL, c1, c2 = get_KL_c1_c2(mode)
    Jp, ap, bp = JKapbp_to_Jpapbp(JK, ap, bp, KL)
    J, M, h = Jpapbp_to_JMh(Jp, ap, bp, c1, c2)
    X, Y, Z = JMh_to_XYZ(J, M, h)
    return XYZ_to_sRGB(X, Y, Z)


########## Similarity functions   ################
# These similarity functions are equivalent to converting the color
# to J/K a' b' and taking Euclidean distance.
    
def deltaEp_JMh(J1, M1, h1, J2, M2, h2, KL, c1, c2):
    Jp1, ap1, bp1 = JMh_to_Jpapbp(J1, M1, h1, c1, c2)
    Jp2, ap2, bp2 = JMh_to_Jpapbp(J2, M2, h2, c1, c2)

    return np.sqrt(
        ((Jp1 - Jp2) / KL) ** 2
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

def test_JKapbp():
    R1, G1, B1 = [1], [1], [1] # white
    R2, G2, B2 = [0], [0], [0] # black
    distance1 = deltaEp_sRGB(R1, G1, B1, R2, G2, B2, mode='UCS')

    JK1, ap1, bp1 = sRGB_to_JKapbp(R1, G1, B1, mode='UCS')
    JK2, ap2, bp2 = sRGB_to_JKapbp(R2, G2, B2, mode='UCS')
    distance2 = np.sqrt((JK1 - JK2)**2
                        + (ap1-ap2)**2
                        + (bp1-bp2)**2)

    assert abs(distance1-distance2) < .001

    distance1L = deltaEp_sRGB(R1, G1, B1, R2, G2, B2, mode='LCD')
    JK1, ap1, bp1 = sRGB_to_JKapbp(R1, G1, B1, mode='LCD')
    JK2, ap2, bp2 = sRGB_to_JKapbp(R2, G2, B2, mode='LCD')
    distance2L = np.sqrt((JK1 - JK2)**2
                         + (ap1-ap2)**2
                         + (bp1-bp2)**2)
    
    assert abs(distance1L-distance2L) < .001

    distance1S = deltaEp_sRGB(R1, G1, B1, R2, G2, B2, mode='SCD')
    JK1, ap1, bp1 = sRGB_to_JKapbp(R1, G1, B1, mode='SCD')
    JK2, ap2, bp2 = sRGB_to_JKapbp(R2, G2, B2, mode='SCD')
    distance2S = np.sqrt((JK1 - JK2)**2
                         + (ap1-ap2)**2
                         + (bp1-bp2)**2)
    
    assert abs(distance1S-distance2S) < .001    
