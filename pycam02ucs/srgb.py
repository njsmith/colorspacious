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

def C_srgb(C_linear):
    out = np.empty(C_linear.shape, dtype=float)
    linear_portion = (C_linear <= 0.0031308)
    out[linear_portion] = C_linear[linear_portion] * 12.92
    nonlinear_portion = ~linear_portion
    a = 0.055
    out[nonlinear_portion] = (1+a) * C_linear ** (1/2.4) - a
    return out

# http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
sRGB_to_XYZ_matrix = np.array([
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041],
    ])

XYZ_to_sRGB_matrix = np.array([
    [ 3.2404542, -1.5371385, -0.4985314],
    [-0.9692660,  1.8760108,  0.0415560],
    [ 0.0556434, -0.2040259,  1.0572252],
    ])

def XYZ_to_sRGB(X, Y, Z):
    """ Convert XYZ to sRGB, where XYZ are normalized so that reference white
    D65 is X=.9505, Y=1, Z=1.0890 """
    X = np.asarray(X, dtype=float) 
    Y = np.asarray(Y, dtype=float) 
    Z = np.asarray(Z, dtype=float)

    RGB_linear = np.dot(XYZ_to_sRGB_matrix, np.row_stack([X, Y, Z]))

    R = C_srgb(RGB_linear[0, :])
    G = C_srgb(RGB_linear[1, :])
    B = C_srgb(RGB_linear[2, :])
    return R, G, B
    
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

    XYZ = np.dot(sRGB_to_XYZ_matrix,
                 np.row_stack([R_linear, G_linear, B_linear]))
    return XYZ[0, :], XYZ[1, :], XYZ[2, :]

# Test values calculated from http://davengrace.com/cgi-bin/cspace.pl """
# ((gold_RGB, gold_XYZ), ...)
_test_values = ((([18.99/255], [21.75/255], [94.93/255]),
                 ([0.0261219], [0.0152732], [0.1096471])),
                (([55.12/255], [33.14/255], [32.19/255]),
                 ([0.0239318], [0.0201643], [0.0164315])),
                 )
        
def test_sRGB_to_XYZ():
    (one_RGB, one_XYZ), (two_RGB, two_XYZ) = _test_values
     
    gold_X, gold_Y, gold_Z = one_XYZ
    conv_X, conv_Y, conv_Z = sRGB_to_XYZ(*one_RGB)
    assert np.allclose(gold_X, conv_X, rtol=.0001)
    assert np.allclose(gold_Y, conv_Y, rtol=.0001)
    assert np.allclose(gold_Z, conv_Z, rtol=.0001)

    gold_X, gold_Y, gold_Z = two_XYZ
    conv_X, conv_Y, conv_Z = sRGB_to_XYZ(*two_RGB)
    assert np.allclose(gold_X, conv_X, rtol=.0001)
    assert np.allclose(gold_Y, conv_Y, rtol=.0001)
    assert np.allclose(gold_Z, conv_Z, rtol=.0001)

def test_XYZ_to_sRGB():
    (one_RGB, one_XYZ), (two_RGB, two_XYZ) = _test_values

    gold_R, gold_G, gold_B = one_RGB
    conv_R, conv_G, conv_B = XYZ_to_sRGB(*one_XYZ)
    print gold_R, conv_R
    assert np.allclose(gold_R, conv_R, rtol=.0001)
    assert np.allclose(gold_G, conv_G, rtol=.0001)
    assert np.allclose(gold_B, conv_B, rtol=.0001)

    gold_R, gold_G, gold_B = two_RGB
    conv_R, conv_G, conv_B = XYZ_to_sRGB(*two_XYZ)
    assert np.allclose(gold_R, conv_R, rtol=.0001)
    assert np.allclose(gold_G, conv_G, rtol=.0001)
    assert np.allclose(gold_B, conv_B, rtol=.0001)

def test_inversion():
    R, G, B = [18.99/255], [21.75/255], [94.93/255]
    X, Y, Z = sRGB_to_XYZ(R, G, B)
    Rp, Gp, Bp = XYZ_to_sRGB(X, Y, Z)
    assert np.allclose(R, Rp)
    assert np.allclose(G, Gp)
    assert np.allclose(B, Bp)

    Xp, Yp, Zp = sRGB_to_XYZ(Rp, Gp, Bp)
    assert np.allclose(X, Xp)
    assert np.allclose(Y, Yp)
    assert np.allclose(Z, Zp)

if __name__ == '__main__':
    import nose
    nose.runmodule()
    