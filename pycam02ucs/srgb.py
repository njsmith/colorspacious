# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np

# https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation

def C_linear(C_srgb):
    out = np.empty(C_srgb.shape, dtype=float)
    linear_portion = (C_srgb < 0.04045)
    out[linear_portion] = C_srgb[linear_portion] / 12.92
    nonlinear_portion = ~linear_portion
    a = 0.055
    out[nonlinear_portion] = ((C_srgb[nonlinear_portion] + a)
                              / (1 + a)) ** 2.4
    return out

# http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
magic_matrix = np.array([
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041],
    ])

# RGB in 0-to-1 range; XYZ with its well-defined range (roughly 0-1)
def sRGB_to_XYZ(R, G, B):
    R = np.asarray(R, dtype=float)
    G = np.asarray(G, dtype=float)
    B = np.asarray(B, dtype=float)

    for arr in R, G, B:
        if np.any(arr < 0) or np.any(arr > 1):
            raise ValueError("RGB values must be in between 0 and 1")

    R_linear = C_linear(R)
    G_linear = C_linear(G)
    B_linear = C_linear(B)

    XYZ = np.dot(magic_matrix,
                 np.row_stack([R_linear, G_linear, B_linear]))
    return XYZ[0, :], XYZ[1, :], XYZ[2, :]
