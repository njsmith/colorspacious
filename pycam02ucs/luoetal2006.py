# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np

from .util import stacklast

class LuoEtAl2006UniformSpace(object):
    def __init__(self, KL, c1, c2):
        self.KL = KL
        self.c1 = c1
        self.c2 = c2

    def JMh_to_Jpapbp(self, JMh):
        JMh = np.asarray(JMh, dtype=float)
        J = JMh[..., 0]
        M = JMh[..., 1]
        h = JMh[..., 2]
        Jp = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        Jp = Jp / self.KL
        Mp = (1. / self.c2) * np.log(1 + self.c2 * M)
        h_rad = np.deg2rad(h)
        ap = Mp * np.cos(h_rad)
        bp = Mp * np.sin(h_rad)
        return stacklast(Jp, ap, bp)

    def Jpapbp_to_JMh(self, Jpapbp):
        Jpapbp = np.asarray(Jpapbp)
        Jp = Jpapbp[..., 0]
        ap = Jpapbp[..., 1]
        bp = Jpapbp[..., 2]
        Jp = Jp * self.KL
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
        return stacklast(J, M, h)

CAM02UCS = LuoEtAl2006UniformSpace(1.00, 0.007, 0.0228)
CAM02LCD = LuoEtAl2006UniformSpace(1.24, 0.007, 0.0363)
CAM02SCD = LuoEtAl2006UniformSpace(0.77, 0.007, 0.0053)
