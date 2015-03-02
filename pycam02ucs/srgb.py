# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np

# https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation

def C_linear(C_srgb):
    linear_portion = (C_srgb < 0.04045)
    a = 0.055
    out = np.select([linear_portion, ~linear_portion],
                    [C_srgb / 12.92,
                     ((C_srgb + a) / (a + 1)) ** 2.4])
    return out

def C_srgb(C_linear):
    out = np.empty(C_linear.shape, dtype=float)
    linear_portion = (C_linear <= 0.0031308)
    a = 0.055
    out = np.select([linear_portion, ~linear_portion],
                    [C_linear * 12.92,
                     (1+a) * C_linear ** (1/2.4) - a])
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

def XYZ_to_sRGB(XYZ):
    """ Convert XYZ to sRGB, where XYZ are normalized so that reference white
    D65 is X=95.05, Y=100, Z=108.90 """
    XYZ = np.asarray(XYZ, dtype=float)
    # this is broadcasting matrix * array-of-vectors, where the vector is the
    # last dim
    RGB_linear = np.einsum("ij,...j->...i", XYZ_to_sRGB_matrix, XYZ / 100)
    RGB = C_srgb(RGB_linear)
    return RGB
    
# RGB in 0-to-1 range; XYZ with its well-defined range (roughly 0-100)
def sRGB_to_XYZ(RGB):
    RGB = np.asarray(RGB)
    if np.any(RGB < 0) or np.any(RGB > 1):
        raise ValueError("RGB values must be in between 0 and 1")
    RGB_linear = C_linear(RGB)
    # this is broadcasting matrix * array-of-vectors, where the vector is the
    # last dim
    XYZ = np.einsum("ij,...j->...i", sRGB_to_XYZ_matrix, RGB_linear)
    XYZ *= 100
    return XYZ

# Test values calculated from http://davengrace.com/cgi-bin/cspace.pl """
# ((gold_RGB, gold_XYZ), ...)
_test_values = (([18.99/255, 21.75/255, 94.93/255],  # RGB
                 [2.61219, 1.52732, 10.96471]), # XYZ
                ([55.12/255, 33.14/255, 32.19/255],  # RGB
                 [2.39318, 2.01643, 1.64315]), # XYZ
                 )
        
def test_sRGB_to_XYZ():
    (one_RGB, one_XYZ), (two_RGB, two_XYZ) = _test_values
     
    conv1_XYZ = sRGB_to_XYZ(one_RGB)
    assert np.allclose(one_XYZ, conv1_XYZ, rtol=.0001)

    conv2_XYZ = sRGB_to_XYZ(two_RGB)
    assert np.allclose(two_XYZ, conv2_XYZ, rtol=.0001)

    both_RGB = np.asarray([one_RGB, two_RGB])
    both_XYZ = np.asarray([one_XYZ, two_XYZ])
    both_conv_XYZ = sRGB_to_XYZ(both_RGB)
    assert np.allclose(both_XYZ, both_conv_XYZ, rtol=.0001)

    # make sure we can handle higher dimensional arrays
    # last dimension is always 3, other ones just get vectorized over
    RGB_3d = both_RGB[np.newaxis, ...]
    XYZ_3d = both_XYZ[np.newaxis, ...]
    conv_XYZ_3d = sRGB_to_XYZ(RGB_3d)
    assert np.allclose(XYZ_3d, conv_XYZ_3d, rtol=0.0001)

def test_XYZ_to_sRGB():
    (one_RGB, one_XYZ), (two_RGB, two_XYZ) = _test_values

    conv1_RGB = XYZ_to_sRGB(one_XYZ)
    assert np.allclose(one_RGB, conv1_RGB, rtol=.0001)

    conv2_RGB = XYZ_to_sRGB(two_XYZ)
    assert np.allclose(two_RGB, conv2_RGB, rtol=.0001)

    both_RGB = np.asarray([one_RGB, two_RGB])
    both_XYZ = np.asarray([one_XYZ, two_XYZ])
    both_conv_RGB = XYZ_to_sRGB(both_XYZ)
    assert np.allclose(both_RGB, both_conv_RGB, rtol=.0001)

    # make sure we can handle higher dimensional arrays
    # last dimension is always 3, other ones just get vectorized over
    RGB_3d = both_RGB[np.newaxis, ...]
    XYZ_3d = both_XYZ[np.newaxis, ...]
    conv_RGB_3d = XYZ_to_sRGB(XYZ_3d)
    assert np.allclose(RGB_3d, conv_RGB_3d, rtol=0.0001)

def test_inversion():
    RGB = [18.99/255, 21.75/255, 94.93/255]
    XYZ = sRGB_to_XYZ(RGB)
    RGBnew = XYZ_to_sRGB(XYZ)
    assert np.allclose(RGB, RGBnew)

    XYZnew = sRGB_to_XYZ(RGBnew)
    assert np.allclose(XYZ, XYZnew)

