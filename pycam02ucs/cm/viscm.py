# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

# Simple script using CIECAM02 and CAM02-UCS to visualize properties of a
# matplotlib colormap

import os.path

import numpy as np

# Most of this file doesn't actually need matpotlib, and I'm too lazy ATM to
# get matplotlib installed on travis. So this lets the travis build go
# through.
try:
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d
    from matplotlib.gridspec import GridSpec
    from matplotlib.widgets import Button, Slider
    import matplotlib.colors
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    print("\nWarning! could not import matplotlib\n")
    pass

from pycam02ucs import ViewingConditions
from pycam02ucs.cam02ucs import UCS_space, SCD_space, LCD_space
from pycam02ucs.srgb import sRGB_to_XYZ, XYZ_to_sRGB
from pycam02ucs.cm.minimvc import Trigger

# Our preferred space (mostly here so we can easily tweak it when curious)
UNIFORM_SPACE = UCS_space

def _sRGB_to_CIECAM02(RGB):
    XYZ = sRGB_to_XYZ(RGB)
    return ViewingConditions.sRGB.XYZ_to_CIECAM02(XYZ)

def _CIECAM02_to_JKapbp(ciecam02):
    JMh = np.column_stack((ciecam02.J, ciecam02.M, ciecam02.h))
    return UNIFORM_SPACE.JMh_to_JKapbp(JMh)


def _JKapbp_to_JMh(JKapbp):
    return UNIFORM_SPACE.JKapbp_to_JMh(JKapbp)


def _JMh_to_sRGB(JMh):
    XYZ = ViewingConditions.sRGB.CIECAM02_to_XYZ(J=JMh[..., 0],
                                                 M=JMh[..., 1],
                                                 h=JMh[..., 2])
    return XYZ_to_sRGB(XYZ)


def _show_cmap(ax, rgb):
    ax.imshow(rgb[np.newaxis, ...], aspect="auto")

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
    return np.clip(np.einsum("...ij,...j->...i", mat, rgb), 0, 1)

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


# Adapt a matplotlib colormap to a linearly transformed version -- useful for
# visualizing how colormaps look given color deficiency.
# Kinda a hack, b/c we inherit from Colormap (this is required), but then
# ignore its implementation entirely.
class TransformedCMap(matplotlib.colors.Colormap):
    def __init__(self, transform, base_cmap):
        self.transform = np.asarray(transform)
        assert self.transform.shape[0] == self.transform.shape[1]
        # Convert RGB transformation matrix to RGBA transformation matrix
        if self.transform.shape == (3, 3):
            t = np.eye(4)
            t[:3, :3] = self.transform
            self.transform = t
        self.base_cmap = base_cmap

    def __call__(self, *args, **kwargs):
        fx = self.base_cmap(*args, **kwargs)
        tfx = _apply_rgb_mat(self.transform, fx)
        return tfx

    def set_bad(self, *args, **kwargs):
        self.base_cmap.set_bad(*args, **kwargs)

    def set_under(self, *args, **kwargs):
        self.base_cmap.set_under(*args, **kwargs)

    def set_over(self, *args, **kwargs):
        self.base_cmap.set_over(*args, **kwargs)

    def is_gray(self):
        return False

def _vis_axes():
    grid = GridSpec(10, 4,
                    left=0.02,
                    right=0.98,
                    bottom=0.02,
                    width_ratios=[1] * 4,
                    height_ratios=[1] * 10)
    axes = {'cmap': grid[0, 0],
            'deltas': grid[1:4, 0],

            'cmap-greyscale': grid[0, 1],
            'lightness-deltas': grid[1:4, 1],

            'deuteranomaly': grid[4, 0],
            'deuteranopia': grid[5, 0],
            'protanomaly': grid[4, 1],
            'protanopia': grid[5, 1],

            # 'lightness': grid[4:6, 1],
            # 'colourfulness': grid[4:6, 2],
            # 'hue': grid[4:6, 3],

            'image0': grid[0:3, 2],
            'image0-cb': grid[0:3, 3],
            'image1': grid[3:7, 2],
            'image1-cb': grid[3:7, 3],
            'image2': grid[7:, 2],
            'image2-cb': grid[7:, 3],
    }

    axes = {key: plt.subplot(value) for (key, value) in axes.items()}
    axes['gamut'] = plt.subplot(grid[6:, :2], projection='3d')
    axes['gamut-toggle'] = plt.axes([0.01, 0.01, 0.08, 0.025])

    return axes


# N=256 matches the default quantization for LinearSegmentedColormap, which
# reduces quantization/aliasing artifacts (esp. in the perceptual deltas
# plot).
class viscm(object):
    def __init__(self, cm, name=None, N=256, N_dots=50, show_gamut=False):
        if isinstance(cm, str):
            cm = plt.get_cmap(cm)
        if name is None:
            name = cm.name

        fig = plt.figure()
        fig.suptitle("Colormap evaluation: %s" % (name,), fontsize=24)
        axes = _vis_axes()

        x = np.linspace(0, 1, N)
        x_dots = np.linspace(0, 1, N_dots)
        RGB = cm(x)[:, :3]
        RGB_dots = cm(x_dots)[:, :3]

        ax = axes['cmap']
        _show_cmap(ax, RGB)
        ax.set_title("The colormap in its glory")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        def label(ax, s):
            ax.text(0.95, 0.05, s,
                    horizontalalignment="right",
                    verticalalignment="bottom",
                    transform=ax.transAxes)

        ax = axes['deltas']
        local_deltas = N * UNIFORM_SPACE.deltaEp_sRGB(RGB[:-1, :], RGB[1:, :])
        ax.plot(x[1:], local_deltas)
        arclength = np.sum(local_deltas) / N
        label(ax, "Perceptual deltas (total: %0.2f)" % (arclength,))
        ax.set_ylim(0, ax.get_ylim()[1])
        ax.get_xaxis().set_visible(False)

        ciecam02 = _sRGB_to_CIECAM02(RGB)
        JKapbp = _CIECAM02_to_JKapbp(ciecam02)

        ax = axes['cmap-greyscale']
        grey_RGB = _JMh_to_sRGB(np.column_stack((ciecam02.J,
                                                 np.zeros_like(ciecam02.M),
                                                 ciecam02.h)))
        _show_cmap(ax, grey_RGB)
        ax.set_title("Black-and-white printed")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        ax = axes['lightness-deltas']
        ax.axhline(0, linestyle="--", color="grey")
        lightness_deltas = N * np.diff(JKapbp[:, 0])
        ax.plot(x[1:], lightness_deltas)
        label(ax,
              "Perceptual lightness deltas (total: %0.2f)"
              % (np.sum(np.abs(lightness_deltas)) / N,))
        #ax.set_ylim(0, ax.get_ylim()[1])
        ax.get_xaxis().set_visible(False)

        # ax = axes['lightness']
        # ax.plot(x, ciecam02.J)
        # label(ax, "Lightness (J)")
        # ax.set_ylim(0, 105)

        # ax = axes['colourfulness']
        # ax.plot(x, ciecam02.M)
        # label(ax, "Colourfulness (M)")

        # ax = axes['hue']
        # ax.plot(x, ciecam02.h)
        # label(ax, "Hue angle (h)")
        # ax.set_ylim(0, 360)

        def anom(ax, mat, name):
            _show_cmap(ax, _apply_rgb_mat(mat, RGB))
            label(ax, name)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        anom(axes['deuteranomaly'], DEUTERANOMALY_05, "Moderate deuteranomaly")
        anom(axes['deuteranopia'], DEUTERANOMALY_10, "Complete deuteranopia")

        anom(axes['protanomaly'], PROTANOMALY_05, "Moderate protanomaly")
        anom(axes['protanopia'], PROTANOMALY_10, "Complete protanopia")

        ax = axes['gamut']
        ax.plot(JKapbp[:, 1], JKapbp[:, 2], JKapbp[:, 0])
        JKapbp_dots = _CIECAM02_to_JKapbp(_sRGB_to_CIECAM02(RGB_dots))
        ax.scatter(JKapbp_dots[:, 1],
                   JKapbp_dots[:, 2],
                   JKapbp_dots[:, 0],
                   c=RGB_dots[:, :],
                   s=80)

        # Draw a wireframe indicating the sRGB gamut
        self.gamut_patch = sRGB_gamut_patch()
        # That function returns a patch where each face is colored to match
        # the represented colors. For present purposes we want something
        # less... colorful.
        self.gamut_patch.set_facecolor([0.5, 0.5, 0.5, 0.1])
        self.gamut_patch.set_edgecolor([0.2, 0.2, 0.2, 0.1])
        ax.add_collection3d(self.gamut_patch)
        self.gamut_patch.set_visible(show_gamut)

        self.gamut_patch_toggle = Button(axes['gamut-toggle'], "Toggle gamut")
        def toggle(*args):
            self.gamut_patch.set_visible(not self.gamut_patch.get_visible())
            plt.draw()
        self.gamut_patch_toggle.on_clicked(toggle)

        _setup_JKapbp_axis(ax)

        images = []
        image_args = []
        example_dir = os.path.dirname(__file__) + "/examples/"

        images.append(np.loadtxt(example_dir + "hist2d.txt"))
        image_args.append({"aspect": "equal",
                           "origin": "lower",
                           "interpolation": "nearest",
                           "vmin": 0})

        images.append(np.loadtxt(example_dir
                                 + "st-helens_before-modified.txt.gz").T)
        image_args.append({})

        # Adapted from http://matplotlib.org/mpl_examples/images_contours_and_fields/pcolormesh_levels.py
        dx = dy = 0.05
        y, x = np.mgrid[-5 : 5 + dy : dy, -5 : 10 + dx : dx]
        z = np.sin(x) ** 10 + np.cos(10 + y * x) + np.cos(x) + 0.2 * y + 0.1 * x
        images.append(z)
        image_args.append({})

        deuter_cm = TransformedCMap(DEUTERANOMALY_05, cm)
        for i, (image, args) in enumerate(zip(images, image_args)):
            ax = axes['image%i' % (i,)]
            ax.imshow(image, cmap=cm, **args)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

            ax_cb = axes['image%i-cb' % (i,)]
            ax_cb.imshow(image, cmap=deuter_cm, **args)
            ax_cb.get_xaxis().set_visible(False)
            ax_cb.get_yaxis().set_visible(False)

        axes['image0'].set_title("Sample images")
        axes['image0-cb'].set_title("Moderate deuter.")

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
                                   [fixed, (i+1) * step, (j+1) * step],
                                   [fixed, i * step, (j+1) * step]])
                sRGB_values.append((fixed, (i + 0.5) * step, (j + 0.5) * step,
                                    1))
                # G quad
                sRGB_quads.append([[i * step, fixed, j * step],
                                   [(i+1) * step, fixed, j * step],
                                   [(i+1) * step, fixed, (j+1) * step],
                                   [i * step, fixed, (j+1) * step]])
                sRGB_values.append(((i + 0.5) * step, fixed, (j + 0.5) * step,
                                    1))
                # B quad
                sRGB_quads.append([[i * step, j * step, fixed],
                                   [(i+1) * step, j * step, fixed],
                                   [(i+1) * step, (j+1) * step, fixed],
                                   [i * step, (j+1) * step, fixed]])
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
    sRGBA = np.concatenate((sRGB, np.ones(sRGB.shape[:2] + (1,))),
                           axis=2)
    sRGBA[np.any((sRGB < 0) | (sRGB > 1), axis=-1)] = [0, 0, 0, 0]
    return sRGBA


def draw_pure_hue_angles(ax):
    # Pure hue angles from CIECAM-02
    for color, angle in [("r", 20.14),
                         ("y", 90.00),
                         ("g", 164.25),
                         ("b", 237.53)]:
        x = np.cos(np.deg2rad(angle))
        y = np.sin(np.deg2rad(angle))
        ax.plot([0, x * 1000], [0, y * 1000], color + "--")


def draw_sRGB_gamut_JK_slice(ax, JK, ap_lim=(-50, 50), bp_lim=(-50, 50),
                             **kwargs):
    sRGB = sRGB_gamut_JK_slice(JK, ap_lim=ap_lim, bp_lim=bp_lim, **kwargs)
    im = ax.imshow(sRGB, aspect="equal",
                   extent=ap_lim + bp_lim, origin="lower")
    draw_pure_hue_angles(ax)
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


def _viscm_editor_axes():
    grid = GridSpec(1, 2,
                    width_ratios=[5, 1],
                    height_ratios=[6, 1])
    axes = {'bezier': grid[0, 0],
            'cm': grid[0, 1]}

    axes = {key: plt.subplot(value) for (key, value) in axes.items()}
    return axes


class viscm_editor(object):
    def __init__(self, min_JK=15, max_JK=95, xp=None, yp=None):
        from pycam02ucs.cm.bezierbuilder import BezierModel, BezierBuilder

        axes = _viscm_editor_axes()

        ax_btn_wireframe = plt.axes([0.7, 0.15, 0.1, 0.025])
        self.btn_wireframe = Button(ax_btn_wireframe, 'Show 3D gamut')
        self.btn_wireframe.on_clicked(self.plot_3d_gamut)

        ax_btn_wireframe = plt.axes([0.81, 0.15, 0.1, 0.025])
        self.btn_save = Button(ax_btn_wireframe, 'Save colormap')
        self.btn_save.on_clicked(self.save_colormap)

        ax_btn_props = plt.axes([0.81, 0.1, 0.1, 0.025])
        self.btn_props = Button(ax_btn_props, 'Properties')
        self.btn_props.on_clicked(self.show_viscm)
        self.prop_windows = []

        axcolor = 'None'
        ax_jk_min = plt.axes([0.1, 0.1, 0.5, 0.03], axisbg=axcolor)
        ax_jk_min.imshow(np.linspace(0, 100, 101).reshape(1, -1), cmap='gray')
        ax_jk_min.set_xlim(0, 100)

        ax_jk_max = plt.axes([0.1, 0.15, 0.5, 0.03], axisbg=axcolor)
        ax_jk_max.imshow(np.linspace(0, 100, 101).reshape(1, -1), cmap='gray')

        self.jk_min_slider = Slider(ax_jk_min, r"$J/K_\mathrm{min}$", 0, 100, valinit=min_JK)
        self.jk_max_slider = Slider(ax_jk_max, r"$J/K_\mathrm{max}$", 0, 100, valinit=max_JK)

        self.jk_min_slider.on_changed(self._jk_update)
        self.jk_max_slider.on_changed(self._jk_update)

        # This is my favorite set of control points so far (just from playing
        # around with things):
        #   min_JK = 15
        #   max_JK = 95
        #   xp =
        #     [-4, 27.041103603603631, 84.311067635550557, 12.567076579094476, -9.6]
        #   yp =
        #     [-34, -41.447876447876524, 36.28563443264386, 25.357741755170423, 41]
        # -- njs, 2015-04-05

        if xp is None:
            xp = [-4, 38.289146128951984, 52.1923711457504,
                  39.050944362271053, 18.60872492130315, -9.6]

        if yp is None:
            yp = [-34, -34.34528254916614, -21.594701710471412,
                  31.701084689194829, 29.510846891948262, 41]

        self.bezier_model = BezierModel(xp, yp)
        self.cmap_model = BezierCMapModel(self.bezier_model,
                                          self.jk_min_slider.val,
                                          self.jk_max_slider.val)
        self.highlight_point_model = HighlightPointModel(self.cmap_model, 0.5)

        self.bezier_builder = BezierBuilder(axes['bezier'], self.bezier_model)
        self.bezier_gamut_viewer = GamutViewer2D(axes['bezier'],
                                                 self.highlight_point_model)
        tmp = HighlightPoint2DView(axes['bezier'],
                                   self.highlight_point_model)
        self.bezier_highlight_point_view = tmp

        draw_pure_hue_angles(axes['bezier'])
        axes['bezier'].set_xlim(-100, 100)
        axes['bezier'].set_ylim(-100, 100)

        self.cmap_view = CMapView(axes['cm'], self.cmap_model)
        self.cmap_highlighter = HighlightPointBuilder(
            axes['cm'],
            self.highlight_point_model)

    def plot_3d_gamut(self, event):
        fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
        self.wireframe_view = WireframeView(ax,
                                            self.cmap_model,
                                            self.highlight_point_model)
        plt.show()

    def save_colormap(self, event):
        import textwrap

        template = textwrap.dedent('''
        from matplotlib.colors import LinearSegmentedColormap
        from numpy import nan, inf

        # Used to reconstruct the colormap in pycam02ucs.cm.viscm
        parameters = {{'xp': {xp},
                      'yp': {yp},
                      'min_JK': {min_JK},
                      'max_JK': {max_JK}}}

        cm_data = {array_list}

        test_cm = LinearSegmentedColormap.from_list(__file__, cm_data)


        if __name__ == "__main__":
            import matplotlib.pyplot as plt
            import numpy as np

            try:
                from pycam02ucs.cm.viscm import viscm
                viscm(test_cm)
            except ImportError:
                print("pycam02ucs not found, falling back on simple display")
                plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',
                           cmap=test_cm)
            plt.show()
        ''')

        rgb, _ = self.cmap_model.get_sRGB(num=256)
        with open('/tmp/new_cm.py', 'w') as f:
            array_list = np.array_repr(rgb, max_line_width=78)
            array_list = array_list.replace('array(', '')[:-1]

            xp, yp = self.cmap_model.bezier_model.get_control_points()

            data = dict(array_list=array_list,
                        xp=xp,
                        yp=yp,
                        min_JK=self.cmap_model.min_JK,
                        max_JK=self.cmap_model.max_JK)

            f.write(template.format(**data))

            print("*" * 50)
            print("Saved colormap to /tmp/new_cm.py")
            print("*" * 50)

    def show_viscm(self, event):
        cm = LinearSegmentedColormap.from_list(
            'test_cm',
            self.cmap_model.get_sRGB(num=64)[0])
        self.prop_windows.append(viscm(cm, name='test_cm'))
        plt.show()

    def _jk_update(self, val):
        jk_min = self.jk_min_slider.val
        jk_max = self.jk_max_slider.val

        smallest, largest = min(jk_min, jk_max), max(jk_min, jk_max)
        if (jk_min > smallest) or (jk_max < largest):
            self.jk_min_slider.set_val(smallest)
            self.jk_max_slider.set_val(largest)

        self.cmap_model.set_JK_minmax(smallest, largest)

class BezierCMapModel(object):
    def __init__(self, bezier_model, min_JK, max_JK):
        self.bezier_model = bezier_model
        self.min_JK = min_JK
        self.max_JK = max_JK
        self.trigger = Trigger()

        self.bezier_model.trigger.add_callback(self.trigger.fire)

    def set_JK_minmax(self, min_JK, max_JK):
        self.min_JK = min_JK
        self.max_JK = max_JK
        self.trigger.fire()

    def get_JKapbp_at(self, at):
        ap, bp = self.bezier_model.get_bezier_points_at(at)
        JK = (self.max_JK - self.min_JK) * at + self.min_JK
        return JK, ap, bp

    def get_JKapbp(self, num=200):
        # ap, bp = self.bezier_model.get_bezier_points(num)
        # assert ap.ndim == bp.ndim == 1
        # JK = np.linspace(self.min_JK, self.max_JK, num=ap.shape[0])
        # return JK, ap, bp
        return self.get_JKapbp_at(np.linspace(0, 1, num))

    def get_sRGB(self, num=200):
        # Return sRGB and out-of-gamut mask
        JK, ap, bp = self.get_JKapbp(num=num)
        JMh = _JKapbp_to_JMh(np.column_stack((JK, ap, bp)))
        sRGB = _JMh_to_sRGB(JMh)
        oog = np.any((sRGB > 1) | (sRGB < 0), axis=-1)
        sRGB[oog, :] = np.nan
        return sRGB, oog


class CMapView(object):
    def __init__(self, ax, cmap_model):
        self.ax = ax
        self.cmap_model = cmap_model

        rgb_display, oog_display = self._drawable_arrays()
        self.image = self.ax.imshow(rgb_display, extent=(0, 0.2, 0, 1),
                                    origin="lower")
        self.gamut_alert_image = self.ax.imshow(oog_display,
                                                extent=(0.05, 0.15, 0, 1),
                                                origin="lower")
        self.ax.set_xlim(0, 0.2)
        self.ax.set_ylim(0, 1)
        self.ax.get_xaxis().set_visible(False)

        self.cmap_model.trigger.add_callback(self._refresh)

    def _drawable_arrays(self):
        rgb, oog = self.cmap_model.get_sRGB()
        rgb_display = rgb[:, np.newaxis, :]
        oog_display = np.empty((rgb.shape[0], 1, 4))
        oog_display[...] = [0, 0, 0, 0]
        oog_display[oog, :, :] = [0, 1, 1, 1]
        return rgb_display, oog_display

    def _refresh(self):
        rgb_display, oog_display = self._drawable_arrays()
        self.image.set_data(rgb_display)
        self.gamut_alert_image.set_data(oog_display)


class HighlightPointModel(object):
    def __init__(self, cmap_model, point):
        self._cmap_model = cmap_model
        self._point = point
        self.trigger = Trigger()

        self._cmap_model.trigger.add_callback(self.trigger.fire)

    def get_point(self):
        return self._point

    def set_point(self, point):
        self._point = point
        self.trigger.fire()

    def get_JKapbp(self):
        return self._cmap_model.get_JKapbp_at(self._point)


class HighlightPointBuilder(object):
    def __init__(self, ax, highlight_point_model):
        self.ax = ax
        self.highlight_point_model = highlight_point_model

        self.canvas = self.ax.figure.canvas
        self._in_drag = False
        self.canvas.mpl_connect("button_press_event", self._on_button_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("button_release_event",
                                self._on_button_release)

        self.marker_line = self.ax.axhline(highlight_point_model.get_point(),
                                           linewidth=3, color="r")

        self.highlight_point_model.trigger.add_callback(self._refresh)

    def _on_button_press(self, event):
        if event.inaxes != self.ax:
            return
        if event.button != 1:
            return
        self._in_drag = True
        self.highlight_point_model.set_point(event.ydata)

    def _on_motion(self, event):
        if self._in_drag and event.ydata is not None:
            self.highlight_point_model.set_point(event.ydata)

    def _on_button_release(self, event):
        if event.button != 1:
            return
        self._in_drag = False

    def _refresh(self):
        point = self.highlight_point_model.get_point()
        self.marker_line.set_data([0, 1], [point, point])
        self.canvas.draw()


class GamutViewer2D(object):
    def __init__(self, ax, highlight_point_model,
                 ap_lim=(-50, 50), bp_lim=(-50, 50)):
        self.ax = ax
        self.highlight_point_model = highlight_point_model
        self.ap_lim = ap_lim
        self.bp_lim = bp_lim

        self.bgcolors = {"light": (0.9, 0.9, 0.9),
                         "dark": (0.1, 0.1, 0.1)}
        # We want some hysteresis, so that there's no point where wiggling the
        # line back and forth causes background flickering.
        self.bgcolor_ranges = {"light": (0, 60), "dark": (40, 100)}
        self.bg_opposites = {"light": "dark", "dark": "light"}
        self.bg = "light"
        self.ax.set_axis_bgcolor(self.bgcolors[self.bg])

        self.image = self.ax.imshow([[[0, 0, 0]]], aspect="equal",
                                    extent=ap_lim + bp_lim,
                                    origin="lower")

        self.highlight_point_model.trigger.add_callback(self._refresh)

    def _refresh(self):
        JK, _, _ = self.highlight_point_model.get_JKapbp()
        low, high = self.bgcolor_ranges[self.bg]
        if not (low <= JK <= high):
            self.bg = self.bg_opposites[self.bg]
            self.ax.set_axis_bgcolor(self.bgcolors[self.bg])
        sRGB = sRGB_gamut_JK_slice(JK, self.ap_lim, self.bp_lim)
        self.image.set_data(sRGB)


class HighlightPoint2DView(object):
    def __init__(self, ax, highlight_point_model):
        self.ax = ax
        self.highlight_point_model = highlight_point_model

        _, ap, bp = self.highlight_point_model.get_JKapbp()
        self.marker = self.ax.plot([ap], [bp], "y.", mew=3)[0]

        self.highlight_point_model.trigger.add_callback(self._refresh)

    def _refresh(self):
        _, ap, bp = self.highlight_point_model.get_JKapbp()
        self.marker.set_data([ap], [bp])
        self.ax.figure.canvas.draw()


class WireframeView(object):
    def __init__(self, ax, cmap_model, highlight_point_model):
        self.ax = ax
        self.cmap_model = cmap_model
        self.highlight_point_model = highlight_point_model

        JK, ap, bp = self.cmap_model.get_JKapbp()
        self.line = self.ax.plot([0, 10], [0, 10])[0]
        #self.line = self.ax.plot(JK, ap, bp)[0]

        JK, ap, bp = self.highlight_point_model.get_JKapbp()
        self.marker = self.ax.plot([JK], [ap], [bp], "y.", mew=3)[0]

        gamut_patch = sRGB_gamut_patch()
        # That function returns a patch where each face is colored to match
        # the represented colors. For present purposes we want something
        # less... colorful.
        gamut_patch.set_facecolor([0.5, 0.5, 0.5, 0.1])
        gamut_patch.set_edgecolor([0.2, 0.2, 0.2, 0.1])
        self.ax.add_collection3d(gamut_patch)

        _setup_JKapbp_axis(self.ax)

        #self.cmap_model.trigger.add_callback(self._refresh_line)
        #self.highlight_point_model.trigger.add_callback(self._refresh_point)
        self._refresh_line()
        self._refresh_point()

    def _refresh_line(self):
        JK, ap, bp = self.cmap_model.get_JKapbp()
        self.line.set_data(ap, bp)
        self.line.set_3d_properties(zs=JK)
        self.ax.figure.canvas.draw()

    def _refresh_point(self):
        JK, ap, bp = self.highlight_point_model.get_JKapbp()
        self.marker.set_data([ap], [bp])
        self.marker.set_3d_properties(zs=[JK])
        self.ax.figure.canvas.draw()


if __name__ == "__main__":
    import sys
    import os

    ns = {'__name__': ''}

    if len(sys.argv) > 1:
        cmap_params = sys.argv[1]
        if os.path.isfile(cmap_params):
            with open(cmap_params) as f:
                code = compile(f.read(), cmap_params, 'exec')
                exec(code, globals(), ns)

    params = ns.get('parameters', {})
    viscm_editor(**params)
    plt.show()
