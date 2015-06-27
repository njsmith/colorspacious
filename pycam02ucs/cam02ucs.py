# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.
import numpy as np

from .ciecam02 import ViewingConditions

class LuoEtAl2006UniformSpace(object):
    def __init__(self, KL, c1, c2):
        self.KL = KL
        self.c1 = c1
        self.c2 = c2

    def JMh_to_JKapbp(self, JMh):
        JMh = np.asarray(JMh, dtype=float)
        J = JMh[..., 0]
        M = JMh[..., 1]
        h = JMh[..., 2]
        Jp = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        JK = Jp / self.KL
        Mp = (1. / self.c2) * np.log(1 + self.c2 * M)
        h_rad = np.deg2rad(h)
        ap = Mp * np.cos(h_rad)
        bp = Mp * np.sin(h_rad)
        return np.array([JK, ap, bp]).T

    def JKapbp_to_JMh(self, JKapbp):
        JKapbp = np.asarray(JKapbp)
        JK = JKapbp[..., 0]
        ap = JKapbp[..., 1]
        bp = JKapbp[..., 2]
        Jp = JK * self.KL
        J = - Jp / (self.c1 * Jp - 100 * self.c1 - 1)
        # a' = M' * cos(h)
        # b' = M' * sin(h)
        # M' = b'/sin(h)
        # a' = (b'/sin(h)) * cos(h)
        # a' = b' * cos(h) / sin(h)
        # sin(h) = b'/a' * cos(h), 0 <= h <= 2pi and 0 <= M <= 100
        # -> tan(h) = b'/a'
        # -> h = arctan2(b', a')
        h_rad = np.arctan2(bp, ap)
        Mp = bp/np.sin(h_rad)
        Mpp = ap/np.cos(h_rad)
        np.testing.assert_allclose(Mp, Mpp)
        h = np.rad2deg(h_rad) % 360
        M = (np.exp(self.c2*Mp) - 1) / self.c2
        return np.array([J, M, h]).T

CAM02_UCS = LuoEtAl2006UniformSpace(1.00, 0.007, 0.0228)
CAM02_LCD = LuoEtAl2006UniformSpace(1.24, 0.007, 0.0363)
CAM02_SCD = LuoEtAl2006UniformSpace(0.77, 0.007, 0.0053)


########## Similarity functions   ################

# def _XYZ_to_JMh(XYZ):
#     XYZ = np.asarray(XYZ)
#     vc = ViewingConditions.sRGB
#     ciecam02_color = vc.XYZ_to_CIECAM02(XYZ)
#     return np.array([ciecam02_color.J, ciecam02_color.M, ciecam02_color.h]).T

# def test_inversion_JMh_JKapbp(verbose=False):
#     r = np.random.RandomState(0)

#     def test(space, num_dims):
#         for _ in range(100):
#             RGB = r.rand(*(10,) * (num_dims - 1) + (3,))
#             XYZ = np.asarray(sRGB_to_XYZ(RGB))
#             JMh = np.array(_XYZ_to_JMh(XYZ))
#             if JMh.ndim == 1:
#                 JMh[0] = np.nan # test nans
#             elif JMh.ndim == 2:
#                 JMh[0, 0] = np.nan

#             if verbose:
#                 print("JMh:", JMh)

#             JKapbp = UCS_space.JMh_to_JKapbp(JMh)
#             if verbose:
#                 print("JK a' b':", JKapbp)

#             JMh_new = UCS_space.JKapbp_to_JMh(JKapbp)
#             np.testing.assert_allclose(JMh, JMh_new)

#     test(UCS_space, 1)
#     test(UCS_space, 2)

#     test(LCD_space, 1)
#     test(LCD_space, 2)

#     test(SCD_space, 1)
#     test(SCD_space, 2)
