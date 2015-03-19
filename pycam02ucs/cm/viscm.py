# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

# Simple script using CIECAM02 and CAM02-UCS to visualize properties of a
# matplotlib colormap

import numpy as np

# Most of this file doesn't actually need matpotlib, and I'm too lazy ATM to
# get matplotlib installed on travis. So this lets the travis build go
# through.
try:
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d
    from matplotlib.gridspec import GridSpec
except ImportError:
    pass

from pycam02ucs import ViewingConditions
from pycam02ucs.cam02ucs import deltaEp_sRGB, UCS_space
from pycam02ucs.srgb import sRGB_to_XYZ, XYZ_to_sRGB

def _sRGB_to_CIECAM02(RGB):
    XYZ = sRGB_to_XYZ(RGB)
    return ViewingConditions.sRGB.XYZ_to_CIECAM02(XYZ)

def _CIECAM02_to_JKapbp(ciecam02):
    JMh = np.column_stack((ciecam02.J, ciecam02.M, ciecam02.h))
    return UCS_space.JMh_to_JKapbp(JMh)

def _JKapbp_to_JMh(JKapbp):
    return UCS_space.JKapbp_to_JMh(JKapbp)

def _JMh_to_sRGB(JMh):
    XYZ = ViewingConditions.sRGB.CIECAM02_to_XYZ(J=JMh[..., 0],
                                                 M=JMh[..., 1],
                                                 h=JMh[..., 2])
    return XYZ_to_sRGB(XYZ)

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

# sRGB corners: a' goes from -37.4 to 45
AP_LIM = (-38, 46)
# b' goes from -46.5 to 42
BP_LIM = (-47, 43)
# J'/K goes from 0 to 100
JK_LIM = (-1, 101)

def _setup_JKapbp_axis(ax):
    ax.set_xlabel("a' (green -> red)")
    ax.set_ylabel("b' (blue -> yellow)")
    ax.set_zlabel("J'/K (white -> black)")
    ax.set_xlim(*AP_LIM)
    ax.set_ylim(*BP_LIM)
    ax.set_zlim(*JK_LIM)


def _vis_axes(editor=False):
    grid = GridSpec(5, 7,
                    width_ratios=[1, 1, 1, 1, 1, 1, 6],
                    height_ratios=[1, 1, 1, 1, 2])
    axes = {'cmap': grid[0, :3],
            'deltas': grid[0, 3:6],
            'deuteranomaly': grid[1, :3],
            'deuteranopia': grid[1, 3:6],
            'protanomaly': grid[2, :3],
            'protanopia': grid[2, 3:6],
            'lightness': grid[3, :2],
            'colourfulness': grid[3, 2:4],
            'hue': grid[3, 4:6]}

    if editor:
        axes['editor'] = grid[:, 6]

    axes = {key: plt.subplot(value) for (key, value) in axes.items()}
    axes['gamut'] = plt.subplot(grid[4, :6], projection='3d')

    return axes

# N=256 matches the default quantization for LinearSegmentedColormap, which
# reduces quantization/aliasing artifacts (esp. in the perceptual deltas
# plot).
def viscm(cm, name=None, N=256, N_dots=50, show_gamut=False,
          axes=None, editor=False):
    if isinstance(cm, str):
        cm = plt.get_cmap(cm)
    if name is None:
        name = cm.name

    if axes is None:
        fig = plt.figure()
        fig.suptitle("Colormap evaluation: %s" % (name,), fontsize=24)
        axes = _vis_axes(editor=editor)

    x = np.linspace(0, 1, N)
    x_dots = np.linspace(0, 1, N_dots)
    RGB = cm(x)[:, :3]
    RGB_dots = cm(x_dots)[:, :3]

    ax = axes['cmap']
    _show_cmap(ax, RGB)
    ax.set_title("The colormap in its glory")
    ax.get_yaxis().set_visible(False)

    ax = axes['deltas']
    local_deltas = N * deltaEp_sRGB(RGB[:-1, :], RGB[1:, :])
    ax.plot(x[1:], local_deltas)
    arclength = np.sum(local_deltas) / N
    ax.set_title("Perceptual deltas (total: %0.2f)" % (arclength,))
    ax.set_ylim(0, ax.get_ylim()[1])
    # ax.text(0.05, 0.9, "Total length: %0.2f" % (arclength,),
    #         horizontalalignment="left",
    #         verticalalignment="top",
    #         transform=ax.transAxes)

    def anom(ax, mat, name):
        _show_cmap(ax, _apply_rgb_mat(mat, RGB))
        ax.text(0.95, 0.05, name,
                horizontalalignment="right",
                verticalalignment="bottom",
                transform=ax.transAxes)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    anom(axes['deuteranomaly'], DEUTERANOMALY_05, "Moderate deuteranomaly")
    anom(axes['deuteranopia'], DEUTERANOMALY_10, "Complete deuteranopia")

    anom(axes['protanomaly'], PROTANOMALY_05, "Moderate protanomaly")
    anom(axes['protanopia'], PROTANOMALY_10, "Complete protanopia")

    ciecam02 = _sRGB_to_CIECAM02(RGB)
    ax = axes['lightness']
    ax.plot(x, ciecam02.J, label="Lightness (J)")
    ax.set_title("Lightness (J)")
    ax.set_ylim(0, 105)

    ax = axes['colourfulness']
    ax.plot(x, ciecam02.M, label="Colourfulness (M)")
    ax.set_title("Colourfulness (M)")

    ax = axes['hue']
    ax.plot(x, ciecam02.h, label="Hue angle (h)")
    ax.set_title("Hue angle (h)")

    JKapbp = _CIECAM02_to_JKapbp(ciecam02)
    ax = axes['gamut']
    ax.plot(JKapbp[:, 1], JKapbp[:, 2], JKapbp[:, 0])
    JKapbp_dots = _CIECAM02_to_JKapbp(_sRGB_to_CIECAM02(RGB_dots))
    ax.scatter(JKapbp_dots[:, 1],
               JKapbp_dots[:, 2],
               JKapbp_dots[:, 0],
               c=RGB_dots[:, :],
               s=80,
           )

    # Draw a wireframe indicating the sRGB gamut
    if show_gamut:
        gamut_patch = sRGB_gamut_patch()
        # That function returns a patch where each face is colored to match
        # the represented colors. For present purposes we want something
        # less... colorful.
        gamut_patch.set_facecolor([0.5, 0.5, 0.5, 0.1])
        gamut_patch.set_edgecolor([0.2, 0.2, 0.2, 0.1])
        ax.add_collection3d(gamut_patch)

    _setup_JKapbp_axis(ax)

    if editor:
        from .bezierbuilder import BezierBuilder
        from matplotlib.lines import Line2D

        ax = axes['editor']
        line, = ax.plot([-4, 40, -9.6], [-34, 4.6, 41], ls='--', c='#666666',
                       marker='x', mew=2, mec='#204a87')

        draw_sRGB_gamut_JK_slice(ax, 50)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)

        bezier = BezierBuilder(line)
        print(bezier.bezier_curve.get_data())

def sRGB_gamut_patch(resolution=20):
    step = 1.0 / resolution
    sRGB_quads = []
    sRGB_values = []
    # each entry in 'quads' is a 4x3 array where each row contains the
    # coordinates of a corner point
    for fixed in 0, 1:
        for i in range(resolution):
            for j in range(resolution):
                # R quad
                sRGB_quads.append([[fixed, i * step, j * step],
                                   [fixed, (i+1) * step, j * step],
                                   [fixed, (i+1) *step, (j+1) * step],
                                   [fixed, i * step, (j+1) * step],
                               ])
                sRGB_values.append((fixed, (i + 0.5) * step, (j + 0.5) * step,
                                    1))
                # G quad
                sRGB_quads.append([[i * step, fixed, j * step],
                                   [(i+1) * step, fixed, j * step],
                                   [(i+1) *step, fixed, (j+1) * step],
                                   [i * step, fixed, (j+1) * step],
                               ])
                sRGB_values.append(((i + 0.5) * step, fixed, (j + 0.5) * step,
                                    1))
                # B quad
                sRGB_quads.append([[i * step, j * step, fixed],
                                   [(i+1) * step, j * step, fixed],
                                   [(i+1) *step, (j+1) * step, fixed],
                                   [i * step, (j+1) * step, fixed],
                               ])
                sRGB_values.append(((i + 0.5) * step, (j + 0.5) * step, fixed,
                                    1))
    sRGB_quads = np.asarray(sRGB_quads)
    # work around colorspace transform bugginess in handling high-dim
    # arrays
    sRGB_quads_2d = sRGB_quads.reshape((-1, 3))
    CIECAM02_quads_2d = _sRGB_to_CIECAM02(sRGB_quads_2d)
    JKapbp_quads_2d = _CIECAM02_to_JKapbp(CIECAM02_quads_2d)
    JKapbp_quads = JKapbp_quads_2d.reshape((-1, 4, 3))
    gamut_patch = mpl_toolkits.mplot3d.art3d.Poly3DCollection(
        JKapbp_quads[:, :, [1, 2, 0]])
    gamut_patch.set_facecolor(sRGB_values)
    gamut_patch.set_edgecolor(sRGB_values)
    return gamut_patch

def sRGB_gamut_JK_slice(JK,
                        ap_lim=(-50, 50), bp_lim=(-50, 50), resolution=200):
    ap_grid, bp_grid = np.mgrid[ap_lim[0] : ap_lim[1] : resolution * 1j,
                                bp_lim[0] : bp_lim[1] : resolution * 1j]
    JK_grid = JK * np.ones((resolution, resolution))
    JKapbp = np.concatenate((JK_grid[:, :, np.newaxis],
                             ap_grid[:, :, np.newaxis],
                             bp_grid[:, :, np.newaxis]),
                            axis=2)
    JMh = _JKapbp_to_JMh(JKapbp)
    sRGB = _JMh_to_sRGB(JMh)
    sRGB[np.any((sRGB < 0) | (sRGB > 1), axis=-1)] = np.nan
    return sRGB

def draw_sRGB_gamut_JK_slice(ax, JK, ap_lim=(-50, 50), bp_lim=(-50, 50),
                             **kwargs):
    sRGB = sRGB_gamut_JK_slice(JK, ap_lim=ap_lim, bp_lim=bp_lim, **kwargs)
    im = ax.imshow(sRGB, aspect="equal",
                   extent=ap_lim + bp_lim, origin="lower")
    # Pure hue angles from CIECAM-02
    for color, angle in [("r", 20.14),
                         ("y", 90.00),
                         ("g", 164.25),
                         ("b", 237.53),
                     ]:
        x = np.cos(np.deg2rad(angle))
        y = np.sin(np.deg2rad(angle))
        ax.plot([0, x * 1000], [0, y * 1000], color + "--")
    ax.set_xlim(ap_lim)
    ax.set_ylim(bp_lim)
    return im

# def sRGB_gamut_J_slice(J,
#                        ap_lim=(-50, 50), bp_lim=(-50, 50), resolution=200):
#     a_grid, b_grid = np.mgrid[ap_lim[0] : ap_lim[1] : resolution * 1j,
#                               bp_lim[0] : bp_lim[1] : resolution * 1j]
#     J_grid = J * np.ones((resolution, resolution))
#     h = np.rad2deg(np.arctan2(b_grid, a_grid))
#     M = np.hypot(a_grid, b_grid)
#     XYZ = ViewingConditions.sRGB.CIECAM02_to_XYZ(J=J_grid, M=M, h=h)
#     sRGB = XYZ_to_sRGB(XYZ)
#     sRGB[np.any((sRGB < 0) | (sRGB > 1), axis=-1)] = np.nan
#     return sRGB

class viscm_editor(object):
    def __init__(self):
        fig, (ax0, ax1) = plt.subplots(2, 1)

        from .bezierbuilder import BezierBuilder
        from matplotlib.lines import Line2D

        line, = ax0.plot([-4, 40, -9.6], [-34, 4.6, 41], ls='--', c='#666666',
                         marker='x', mew=2, mec='#204a87')

        draw_sRGB_gamut_JK_slice(ax0, 50)
        ax0.set_xlim(-100, 100)
        ax0.set_ylim(-100, 100)

        self.bezier = BezierBuilder(line, update_callback=self.update)
        self.bezier._update_bezier()

    def update(self):
        print(np.sum(self.bezier.bezier_curve.get_data()))
