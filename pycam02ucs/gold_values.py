# This file is part of pycam02ucs
# Copyright (C) 2014-2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

# Test vectors for different conversions, sourced from various places

from collections import namedtuple

from .ciecam02 import CIECAM02Space, CIECAM02Surround, JChQMsH

################################################################
# CIECAM02
################################################################

CIECAM02TestVec = namedtuple("CIECAM02TestVec", ["XYZ100", "vc", "expected"])

_S = CIECAM02Surround
XYZ100_CIECAM02_gold = [
    # Gold values from
    #   https://github.com/igd-geo/pcolor/blob/master/de.fhg.igd.pcolor.test/src/de/fhg/igd/pcolor/test/CAMWorkedExample.java
    # apparently taken from CIE 159:2004 Section 9
    CIECAM02TestVec(XYZ100=[19.31, 23.93, 10.14],
                    vc=CIECAM02Space(XYZ100_w=[98.88, 90, 32.03],
                                     L_A=200,
                                     Y_b=18,
                                     surround=_S(F=1.0,
                                                 c=0.69,
                                                 N_c=1.0)),
                    expected=JChQMsH(h=191.0452, J=48.0314, Q=183.1240,
                                     s=46.0177, C=38.7789, M=38.7789,
                                     H=240.8885)),
    CIECAM02TestVec(XYZ100=[19.31, 23.93, 10.14],
                    vc=CIECAM02Space(XYZ100_w=[98.88, 90, 32.03],
                                     L_A=20, # <- different from above
                                     Y_b=18,
                                     surround=_S(F=1.0,
                                                 c=0.69,
                                                 N_c=1.0)),
                    expected=JChQMsH(h=185.3445, J=47.6856, Q=113.8401,
                                     s=51.1275, C=36.0527, M=29.7580,
                                     H=232.6630)),
    # gold values from Mark Fairchild's spreadsheet at
    #   http://rit-mcsl.org/fairchild//files/AppModEx.xls
    CIECAM02TestVec(XYZ100=[19.01, 20.00, 21.78],
                    vc=CIECAM02Space(XYZ100_w=[95.05, 100.0, 108.88],
                                     Y_b=20.0,
                                     L_A=318.30988618379,
                                     surround=_S(F=1.0,
                                                 c=0.69,
                                                 N_c=1.0)),
                    expected=JChQMsH(h=219.04841, J=41.73109, Q=195.37131,
                                     s=2.36031, C=0.10471, M=0.10884,
                                     H=278.06070)),
    CIECAM02TestVec(XYZ100=[57.06, 43.06, 31.96],
                    vc=CIECAM02Space(XYZ100_w=[95.05, 100.0, 108.88],
                                     L_A=31.830988618379,
                                     Y_b=20.0,
                                     surround=_S(F=1.0,
                                                 c=0.69,
                                                 N_c=1.0)),
                    # The H value here based on the corrected version of the
                    # spreadsheet that I sent Mark Fairchild on
                    # 2014-07-15... the original spreadsheet had it wrong, so
                    # if comparing be careful about which version you have!
                    expected=JChQMsH(h=19.55739, J=65.95523, Q=152.67220,
                                     s=52.24549, C=48.57050, M=41.67327,
                                     H=399.38837)),
]

################################################################
# Elementary conversions
################################################################

# Test values from http://davengrace.com/dave/cspace/
# In their notation, "sRGB'" is regular sRGB, and "sRGB" is sRGB_linear
# AFAICT these are accurate
sRGB_sRGB_linear_gold = [
    ([0.1, 0.2, 0.3],
     [0.010022825574869, 0.0331047665708851, 0.0732389558784054]),
    ([0.9, 0.8, 0.7],
     [0.787412289395617, 0.603827338855338, 0.447988412441883]),
    # make sure we have a test point with values less than the C linearity
    # cutoff.
    ([0.04, 0.02, 0.01],
     [0.00309597523219814, 0.00154798761609907, 0.000773993808049536]),
    ]

# Test values from http://davengrace.com/dave/cspace/
# This uses a rounded-off sRGB->XYZ100 matrix so has some error
sRGB_linear_XYZ100_gold = [
    ([0.00650735931, 0.00789021442, 0.114259116060], # sRGB_linear
     [2.61219, 1.52732, 10.96471]),                  # XYZ100
    ([0.03836396959, 0.01531740787, 0.014587362033], # sRGB_linear
     [2.39318, 2.01643, 1.64315]),                   # XYZ100
    ]

# The calculator at davengrace.com
#    http://davengrace.com/dave/cspace/
# has accurate XYZ100->sRGB and reasonable sRGB->XYZ100, but errorful
# XYZ100<->CIELAB.
#
# The calculator at
#    http://www.brucelindbloom.com/index.html?ColorCalculator.html
# has poor sRGB matrices, but accurate XYZ100<->CIELAB.
# For matching conventions to ours, make sure to check "Scale XYZ100" and
# "Scale Y", but not "Scale RGB".

# Test values from http://www.brucelindbloom.com/index.html?ColorCalculator.html
XYZ100_CIELAB_gold_D65 = [
    ([10, 20, 30],
     [51.8372, -56.3591, -13.1812]),
    ([80, 90, 10],
     [95.9968, -10.6593, 102.8625]),
    # make sure we have a test point with values below the linearity point
    ([0.5, 0.6, 0.4],
     [5.4198, -2.8790, 3.6230]),
    ]

# Test values from http://www.brucelindbloom.com/index.html?ColorCalculator.html
XYZ100_CIELAB_gold_D50 = [
    ([2.61219, 1.52732, 10.96471],  # XYZ100
     [12.7806, 26.1147, -52.4348]), # CIELAB
    ([2.39318, 2.01643, 1.64315],   # XYZ100
     [15.5732, 9.7574, 0.2281]),    # CIELAB
    # make sure we have a test point with values below the linearity point
    ([0.5, 0.6, 0.4],               # XYZ100
     [5.4198, -3.1711, 1.7953]),    # CIELAB
]

################################################################
# Multi-step conversions
################################################################

# http://www.brucelindbloom.com/index.html?ColorCalculator.html
# not very accurate
sRGB_xyY100_gold = [([1.012114, 0.554529, 0.567375],
                     [0.432011, 0.326015, 43.0600]),
                    ([0.2, 0.4, 0.6],
                     [0.210775, 0.222162, 12.5053])]

# http://davengrace.com/dave/cspace/
# (More accurate than the Lindbloom calculator at least for this...)
sRGB_CIELAB_gold_D65 = [
    ([0.2, 0.4, 0.6],
     [42.0099857768665, -0.147373964354935, -32.8445986139017]),
    ([0.1, 1.1, -0.1],
     [95.5971892939889, -92.1527304606657, 91.2351818053272])]

################################################################
# CAM02-UCS and friends
################################################################

# Unfortunately, we do not *have* independent gold test vectors for these,
# because there don't seem to be any other public implementations of these
# transforms. So these results are generated by our own code!
