# This file is part of colorspacious
# Copyright (C) 2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np

def stacklast(*arrs):
    arrs = [np.asarray(arr)[..., np.newaxis] for arr in arrs]
    return np.concatenate(arrs, axis=-1)

# Using color conventions: degrees 0-360
def color_cart2polar(a, b):
    h_rad = np.arctan2(b, a)
    h = np.rad2deg(h_rad) % 360
    r = np.hypot(a, b)
    return (r, h)

def color_polar2cart(r, h):
    h_rad = np.deg2rad(h)
    return (r * np.cos(h_rad), r * np.sin(h_rad))
