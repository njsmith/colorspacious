# BézierBuilder
#
# Copyright (c) 2013, Juan Luis Cano Rodríguez <juanlu001@gmail.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""BézierBuilder, an interactive Bézier curve explorer.

Just run it with

$ python bezier_builder.py

"""

import matplotlib

# WebAgg canvas requires matplotlib >= 1.3 and tornado
#matplotlib.use('webagg')

import numpy as np
from scipy.special import binom

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D


class BezierBuilder(object):
    """Bézier curve interactive builder.

    """
    def __init__(self, control_polygon, ax_bernstein=None, update_callback=None):
        """Constructor.

        Receives the initial control polygon of the curve.

        """
        self.control_polygon = control_polygon
        self.xp = list(control_polygon.get_xdata())
        self.yp = list(control_polygon.get_ydata())
        self.canvas = control_polygon.figure.canvas
        self.ax_main = control_polygon.axes
        self.ax_bernstein = ax_bernstein
        self.update_callback = update_callback

        # Event handler for mouse clicking
        self.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.canvas.mpl_connect('button_release_event', self.on_button_release)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.canvas.mpl_connect('key_release_event', self.on_key_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion_notify)

        # Create Bézier curve
        line_bezier = Line2D([], [],
                             c=control_polygon.get_markeredgecolor())
        self.bezier_curve = self.ax_main.add_line(line_bezier)

        self._shift_is_held = False
        self._ctrl_is_held = False
        self._index = None  # Active vertex

    def on_button_press(self, event):
        # Ignore clicks outside axes
        if event.inaxes != self.ax_main: return
        if self._shift_is_held or self._ctrl_is_held:
            res, ind = self.control_polygon.contains(event)
            if res:
                self._index = ind['ind'][0]
                if self._ctrl_is_held:
                    self._remove_point(event)
            else:
                return
        else:
            self._add_point(event)

    def on_button_release(self, event):
        if event.button != 1: return
        self._index = None

    def on_key_press(self, event):
        if event.key == 'shift':
            self._shift_is_held = True
        elif event.key == 'control':
            self._ctrl_is_held = True

    def on_key_release(self, event):
        if event.key == 'shift':
            self._shift_is_held = False
        elif event.key == 'control':
            self._ctrl_is_held = False

    def on_motion_notify(self, event):
        if event.inaxes != self.ax_main: return
        if self._index is None: return
        x, y = event.xdata, event.ydata

        self.xp[self._index] = x
        self.yp[self._index] = y
        self.control_polygon.set_data(self.xp, self.yp)

        self._update_bezier()

    def _add_point(self, event):
        self.xp.append(event.xdata)
        self.yp.append(event.ydata)
        self.control_polygon.set_data(self.xp, self.yp)

        # Rebuild Bézier curve and update canvas
        self._update_bernstein()
        self._update_bezier()

    def _remove_point(self, event):
        del self.xp[self._index]
        del self.yp[self._index]
        self.control_polygon.set_data(self.xp, self.yp)

        # Rebuild Bézier curve and update canvas
        self._update_bernstein()
        self._update_bezier()

    def _build_bezier(self):
        x, y = Bezier(list(zip(self.xp, self.yp))).T
        return x, y

    def _update_bezier(self):
        self.bezier_curve.set_data(*self._build_bezier())

        if self.update_callback:
            self.update_callback()

        self.canvas.draw()

    def _update_bernstein(self):
        N = len(self.xp) - 1
        t = np.linspace(0, 1, num=200)
        ax = self.ax_bernstein

        if not ax:
            return

        ax.clear()
        for kk in range(N + 1):
            ax.plot(t, Bernstein(N, kk)(t))
        if N > 0:
            ax.set_title("Bernstein basis, N = {}".format(N))
        else:
            ax.set_title("Bernstein basis")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)


def Bernstein(n, k):
    """Bernstein polynomial.

    """
    coeff = binom(n, k)

    def _bpoly(x):
        return coeff * x ** k * (1 - x) ** (n - k)

    return _bpoly


def Bezier(points, num=200):
    """Build Bézier curve from points.

    """
    N = len(points)
    t = np.linspace(0, 1, num=num)
    curve = np.zeros((num, 2))
    for ii in range(N):
        curve += np.outer(Bernstein(N - 1, ii)(t), points[ii])
    return curve


if __name__ == '__main__':
    # Initial setup
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Empty line
    line = Line2D([], [], ls='--', c='#666666',
                  marker='x', mew=2, mec='#204a87')
    ax1.add_line(line)

    # Canvas limits
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_title("Bézier curve")

    # Bernstein plot
    ax2.set_title("Bernstein basis")

    # Create BezierBuilder
    bezier_builder = BezierBuilder(line, ax2)

    fig.suptitle("BézierBuilder", fontsize=24)
    fig.text(0.052, 0.07, "Click to add points, Shift + Click & Drag to move them, "
             "Ctrl + Click to remove them.", color="#333333")
    fig.text(0.052, 0.03, "(c) 2013, Juan Luis Cano Rodríguez. Code available at "
             "https://github.com/Pybonacci/bezierbuilder/", color="#666666")
    fig.subplots_adjust(left=0.08, top=0.8, right=0.92, bottom=0.2)
    plt.show()
