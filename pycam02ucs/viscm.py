# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

# Simple script using CIECAM02 and CAM02-UCS to visualize properties of a
# matplotlib colormap

import numpy as np

from pycam02ucs import ViewingConditions
from pycam02ucs.cam02ucs import deltaEp_sRGB, UCS_space
from pycam02ucs.srgb import sRGB_to_XYZ

def _sRGB_to_CIECAM02(RGB):
    XYZ = sRGB_to_XYZ(RGB)
    return ViewingConditions.sRGB.XYZ_to_CIECAM02(XYZ)

def _CIECAM02_to_JKapbp(ciecam02):
    JMh = np.column_stack((ciecam02.J, ciecam02.M, ciecam02.h))
    return UCS_space.JMh_to_JKapbp(JMh)

def _show_cmap(ax, rgb):
    ax.imshow(rgb[np.newaxis, ...],
              # (left, right, bottom, top)
              # left and right are real
              # bottom and top are used to set the aspect ratio
              extent=(0, 1, 0, 0.2),
              )

# Matrices for simulating anomalous color vision from:
#   Machado, Oliveira, & Fernandes, A Physiologically-based Model for
#   Simulation of Color Vision Deficiency. doi: 10.1109/TVCG.2009.113
#
#   http://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html
# Python code with pre-typed matrices:
#   http://registry.gimp.org/files/color_vision_deficiency.py_1.txt
#
# The 05 variables are "0.5" strength (moderate anomaly), the 10 variables are
# "1.0" strength, i.e., true dichromacy.
#
# Most people with anomalous color vision (~5% of all men) fall somewhere on
# the deuteranomaly spectrum. A minority (~1% of all men) are either fully
# deuteranopic or fall on the protanomaly spectrum. A much smaller number fall
# on the tritanomaly spectrum (<0.1% of people) or have other more exotic
# anomalies.

PROTANOMALY_05 = [[0.458064, 0.679578, -0.137642],
                  [0.092785, 0.846313, 0.060902],
                  [-0.007494, -0.016807, 1.024301]]

DEUTERANOMALY_05 = [[0.547494, 0.607765, -0.155259],
                    [0.181692, 0.781742, 0.036566],
                    [-0.010410, 0.027275, 0.983136]]

TRITANOMALY_05 = [[1.017277, 0.027029, -0.044306],
                  [-0.006113, 0.958479, 0.047634],
                  [0.006379, 0.248708, 0.744913]]

PROTANOMALY_10 = [[0.152286, 1.052583, -0.204868],
                  [0.114503, 0.786281, 0.099216],
                  [-0.003882, -0.048116, 1.051998]]

DEUTERANOMALY_10 = [[0.367322, 0.860646, -0.227968],
                    [0.280085, 0.672501, 0.047413],
                    [-0.011820, 0.042940, 0.968881]]

TRITANOMALY_10 = [[1.255528, -0.076749, -0.178779],
                  [-0.078411, 0.930809, 0.147602],
                  [0.004733, 0.691367, 0.303900]]

def _apply_rgb_mat(mat, rgb):
    return np.clip(np.dot(mat, rgb.T).T, 0, 1)

# N=256 matches the default quantization for LinearSegmentedColormap, which
# reduces quantization/aliasing artifacts (esp. in the perceptual deltas
# plot).
def viscm(cm, name=None, N=256, N_dots=50, show_gamut=False):
    import mpl_toolkits.mplot3d
    import matplotlib.pyplot as plt

    if isinstance(cm, str):
        cm = plt.get_cmap(cm)
    if name is None:
        name = cm.name

    x = np.linspace(0, 1, N)
    x_dots = np.linspace(0, 1, N_dots)
    RGB = cm(x)[:, :3]
    RGB_dots = cm(x_dots)[:, :3]

    fig = plt.figure()
    fig.subplots_adjust(top=0.9, bottom=0.01, left=0.05, right=0.95,
                        wspace=0.3)

    fig.suptitle("Colormap evaluation: %s" % (name,), fontsize=24)

    ax = fig.add_subplot(8, 2, 1)
    _show_cmap(ax, RGB)
    ax.set_title("The colormap in its glory")
    ax.get_yaxis().set_visible(False)

    ax = fig.add_subplot(8, 2, 2)
    local_deltas = N * deltaEp_sRGB(RGB[:-1, :], RGB[1:, :])
    ax.plot(x[1:], local_deltas)
    arclength = np.sum(local_deltas) / N
    ax.set_title("Perceptual deltas (total: %0.2f)" % (arclength,))
    ax.set_ylim(0, ax.get_ylim()[1])
    # ax.text(0.05, 0.9, "Total length: %0.2f" % (arclength,),
    #         horizontalalignment="left",
    #         verticalalignment="top",
    #         transform=ax.transAxes)

    def anom(i, mat, name):
        ax = fig.add_subplot(8, 2, 3 + i)
        _show_cmap(ax, _apply_rgb_mat(mat, RGB))
        ax.text(0.95, 0.05, name,
                horizontalalignment="right",
                verticalalignment="bottom",
                transform=ax.transAxes)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    anom(0, DEUTERANOMALY_05, "Moderate deuteranomaly")
    anom(1, DEUTERANOMALY_10, "Complete deuteranopia")

    anom(2, PROTANOMALY_05, "Moderate protanomaly")
    anom(3, PROTANOMALY_10, "Complete protanopia")

    ciecam02 = _sRGB_to_CIECAM02(RGB)
    ax = fig.add_subplot(8, 3, 10)
    ax.plot(x, ciecam02.J, label="Lightness (J)")
    ax.set_title("Lightness (J)")
    ax.set_ylim(0, 105)
    ax = fig.add_subplot(8, 3, 11)
    ax.plot(x, ciecam02.M, label="Colourfulness (M)")
    ax.set_title("Colourfulness (M)")
    ax = fig.add_subplot(8, 3, 12)
    ax.plot(x, ciecam02.h, label="Hue angle (h)")
    ax.set_title("Hue angle (h)")

    JKapbp = _CIECAM02_to_JKapbp(ciecam02)
    ax = fig.add_subplot(2, 1, 2, projection="3d")
    ax.plot(JKapbp[:, 1], JKapbp[:, 2], JKapbp[:, 0])
    JKapbp_dots = _CIECAM02_to_JKapbp(_sRGB_to_CIECAM02(RGB_dots))
    ax.scatter(JKapbp_dots[:, 1],
               JKapbp_dots[:, 2],
               JKapbp_dots[:, 0],
               c=RGB_dots[:, :],
               s=80,
           )
    ax.set_xlabel("a' (green -> red)")
    ax.set_ylabel("b' (blue -> yellow)")
    ax.set_zlabel("J'/K (white -> black)")

    # Draw a wireframe indicating the sRGB gamut
    if show_gamut:
        GAMUT_POINTS = 20
        step = 1.0 / GAMUT_POINTS
        sRGB_quads = []
        # each entry in 'quads' is a 4x3 array where each row contains the
        # coordinates of a corner point
        for fixed in 0, 1:
            for i in range(GAMUT_POINTS):
                for j in range(GAMUT_POINTS):
                    sRGB_quads.append([[fixed, i * step, j * step],
                                       [fixed, (i+1) * step, j * step],
                                       [fixed, i * step, (j+1) * step],
                                       [fixed, (i+1) *step, (j+1) * step]])
                    sRGB_quads.append([[i * step, fixed, j * step],
                                       [(i+1) * step, fixed, j * step],
                                       [i * step, fixed, (j+1) * step],
                                       [(i+1) *step, fixed, (j+1) * step]])
                    sRGB_quads.append([[i * step, j * step, fixed],
                                       [(i+1) * step, j * step, fixed],
                                       [i * step, (j+1) * step, fixed],
                                       [(i+1) *step, (j+1) * step, fixed]])
        sRGB_quads = np.asarray(sRGB_quads)
        # work around colorspace transform bugginess in handling high-dim
        # arrays
        sRGB_quads_2d = sRGB_quads.reshape((-1, 3))
        CIECAM02_quads_2d = _sRGB_to_CIECAM02(sRGB_quads_2d)
        JKapbp_quads_2d = _CIECAM02_to_JKapbp(CIECAM02_quads_2d)
        JKapbp_quads = JKapbp_quads_2d.reshape((-1, 4, 3))
        gamut_patch = mpl_toolkits.mplot3d.art3d.Poly3DCollection(
            JKapbp_quads[:, :, [1, 2, 0]])
        gamut_patch.set_facecolor([0.5, 0.5, 0.5, 0.1])
        gamut_patch.set_edgecolor([0.2, 0.2, 0.2, 0.1])
        ax.add_collection3d(gamut_patch)

    # sRGB corners: a' goes from -37.4 to 45
    ax.set_xlim(-38, 46)
    # b' goes from -46.5 to 42
    ax.set_ylim(-47, 43)
    # J'/K goes from 0 to 100
    ax.set_zlim(-1, 101)
