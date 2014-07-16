# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.
import numpy as np
import scipy.optimize

from .srgb import sRGB_to_XYZ, XYZ_to_sRGB
from .ciecam02 import ViewingConditions

class LuoUniformSpace(object):
    def __init__(self, KL, c1, c2):
        self.KL = KL
        self.c1 = c1
        self.c2 = c2
        
    def JMh_to_JKapbp(self, (J, M, h)):
        J, M, h = np.asarray(J), np.asarray(M), np.asarray(h)
        Jp = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        JK = Jp / self.KL
        Mp = (1. / self.c2) * np.log(1 + self.c2 * M)
        h_rad = np.deg2rad(h)
        ap = Mp * np.cos(h_rad)
        bp = Mp * np.sin(h_rad)
        return JK, ap, bp

    def JKapbp_to_JMh(self, (JK, ap, bp)):
        JK, ap, bp = np.asarray(JK), np.asarray(ap), np.asarray(bp)
        Jp = JK * self.KL
        J = - Jp / (self.c1 * Jp - 100 * self.c1 - 1)
        # a' = M' * cos(h)
        # b' = M' * sin(h)
        # M' = b'/sin(h)
        # a' = (b'/sin(h)) * cos(h)
        # a' = b' * cos(h) / sin(h)
        # sin(h) = b'/a' * cos(h), 0 <= h <= 2pi and 0 <= M <= 100
        # Thanks Mathematica!
        h_rad = np.zeros(len(bp))
        h_rad_bp_negative = 2*np.pi + 2*np.arctan(-np.sqrt(1+ap**2/bp**2) - ap/bp)
        h_rad_bp_positive = 2*np.arctan(np.sqrt(1+ap**2/bp**2) - ap/bp)
        h_rad[np.where(bp < 0)] = h_rad_bp_negative[np.where(bp<0)]
        h_rad[np.where(bp >= 0)] = h_rad_bp_positive[np.where(bp>=0)]
        Mp = bp/np.sin(h_rad)
        assert np.allclose(Mp, ap/np.cos(h_rad))
        h = np.rad2deg(h_rad)
        M = (np.exp(self.c2*Mp) - 1) / self.c2
        return J, M, h
            
    def deltaEp_JMh(self, (J1, M1, h1), (J2, M2, h2)):
        JK1, ap1, bp1 = self.JMh_to_JKapbp((J1, M1, h1))
        JK2, ap2, bp2 = self.JMh_to_JKapbp((J2, M2, h2))

        return np.sqrt(
            (JK1 - JK2) ** 2
            + (ap1 - ap2) ** 2
            + (bp1 - bp2) ** 2
            )

ucs_space = LuoUniformSpace(1.00, 0.007, 0.0228)
lcd_space = LuoUniformSpace(1.24, 0.007, 0.0363)
scd_space = LuoUniformSpace(0.77, 0.007, 0.0053)    


########## Similarity functions   ################
    
def deltaEp_sRGB(R1, G1, B1, R2, G2, B2, mode='UCS'):
    """Computes the :math:`\delta E'` distance between pairs of sRGB colors.

    :math:`\delta E'` is color difference metric defined by Eq. (4) of Luo et
    al (2006); they show that it provides a good match to human color
    similarity judgements.

    Valid modes are "LCD" (which Luo et al tuned on their "large color
    difference" judgement data sets), "SCD" (which Luo et all tuned on their
    "small color difference" data sets), and "UCS" (which attempts define a
    single generic "uniform color space" which performs well across all their
    data sets). You can also pass in any LuoUniformSpace object as the mode.

    This function is vectorized, i.e., R1, G1, B1, R2, G2, B2 may be vectors,
    in which case we compute the distance between corresponding pairs of
    colors.
    """
    if mode == 'UCS':
        space = ucs_space
    elif mode == 'LCD':
        space = lcd_space
    elif mode == 'SCD':
        space = scd_space
    elif isinstance(mode, LuoUniformSpace):
        space = mode
    else:
        raise ValueError("Invalid mode passed to deltaEp_sRGB")
        
    X1, Y1, Z1 = srgb.sRGB_to_XYZ(R1, G1, B1)
    X2, Y2, Z2 = srgb.sRGB_to_XYZ(R2, G2, B2)
    J1, M1, h1 = _XYZ_to_JMh((X1, Y1, Z1))
    J2, M2, h2 = _XYZ_to_JMh((X2, Y2, Z2))
    return space.deltaEp_JMh((J1, M1, h1), (J2, M2, h2))

def _XYZ_to_JMh(XYZ):
    XYZ = np.asarray(XYZ).T
    vc = ViewingConditions.sRGB
    ciecam02_color = vc.XYZ_to_CIECAM02(XYZ)
    return ciecam02_color.J, ciecam02_color.M, ciecam02_color.h

def test_inversion_JMh_JKapbp(verbose=False):
    r = np.random.RandomState(0)
    for _ in xrange(100):
        R, G, B = r.rand(3, 100) # start with RGB to ensure real colors
        X, Y, Z = sRGB_to_XYZ(R, G, B)
        J, M, h = _XYZ_to_JMh((X, Y, Z))
        if verbose:
            print("JMh:", J, M, h)

        JK, ap, bp = ucs_space.JMh_to_JKapbp((J, M, h))
        if verbose:
            print("JK a' b':", JK, ap, bp)
        
        J_new, M_new, h_new = ucs_space.JKapbp_to_JMh((JK, ap, bp))
        if verbose:
            print("J'M'h':", J_new, M_new, h_new)
        
        assert np.allclose(J, J_new)
        assert np.allclose(M, M_new)
        assert np.allclose(h, h_new)

if __name__ == '__main__':
    import nose
    nose.runmodule()           