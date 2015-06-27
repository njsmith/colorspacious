# This file is part of pycam02ucs
# Copyright (C) 2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np
from collections import defaultdict

from .testing import check_conversion
from .basics import (sRGB_to_sRGB_linear, sRGB_linear_to_sRGB,
                     sRGB_linear_to_XYZ, XYZ_to_sRGB_linear,
                     XYZ_to_xyY, xyY_to_XYZ,
                     XYZ_to_CIELAB, CIELAB_to_XYZ)

from .ciecam02 import ViewingConditions
from .cam02ucs import (LuoEtAl2006UniformSpace,
                       CAM02_UCS, CAM02_LCD, CAM02_SCD)

################################################################

class _PathKeeper(object):
    def __init__(self):
        # maps (source, target) : data
        self._edges = {}
        # maps (source, target) : [source, intermediate node 1, ..., target]
        self._paths = {}
        # {source: all nodes reachable from source}
        self._forward_reachable = defaultdict(set)
        # {target: all nodes that can reach target}
        self._backward_reachable = defaultdict(set)

    def add_connection(self, source, target, data):
        assert (source, target) not in self._edges
        self._edges[(source, target)] = data
        #print("Adding %s -> %s: %s" % (source, target, data))
        queue = [(source, target, [source, target])]
        while queue:
            #print("  To process: %s" % (queue,))
            source, target, path = queue.pop()
            old_path = self._paths.get((source, target))
            if old_path is None or len(path) < len(old_path):
                #print("    %s -> %s: %s beats %s"
                #      % (source, target, path, old_path))
                self._paths[(source, target)] = path
                self._forward_reachable[source].add(target)
                self._backward_reachable[target].add(source)
                for new_target in self._forward_reachable[target]:
                    queue.append(
                        (source, new_target,
                         path[:-1] + self._paths[(target, new_target)])
                    )
                for new_source in self._backward_reachable[source]:
                    queue.append(
                        (new_source, target,
                         self._paths[(new_source, source)][:-1] + path)
                    )
        #print("--- Done ---\n")

    def get_path(self, source, target):
        "Returns (nodes, edges)"
        nodes = self._paths.get((source, target))
        if nodes is None:
            raise ValueError("No path found from %r -> %r" % (source, target))
        edges = []
        for i in range(len(nodes) - 1):
            edges.append(self._edges[(nodes[i], nodes[i + 1])])
        return nodes, edges

    def _dump_dot(self, f):  # pragma: no cover
        f.write("digraph {\n")
        # Hack: technically this class isn't supposed to know about colors,
        # but the layout looks much nicer if we mention sRGB first so it goes
        # on top...
        if "sRGB" in self._forward_reachable:
            f.write("sRGB\n")
        for (source, target), data in self._edges.items():
            f.write("\"%s\" -> \"%s\"\n" % (source, target))
        f.write("}\n")

def test__PathKeeper():
    from nose.tools import assert_raises

    pk = _PathKeeper()
    assert_raises(ValueError, pk.get_path, "a", "b")

    pk.add_connection("a", "b", 1)
    assert pk.get_path("a", "b") == (["a", "b"], [1])

    pk.add_connection("b", "c", 2)
    assert pk.get_path("b", "c") == (["b", "c"], [2])
    assert pk.get_path("a", "c") == (["a", "b", "c"], [1, 2])

    pk.add_connection("y", "z", 25)
    pk.add_connection("z", "a", 26)
    assert pk.get_path("y", "c") == (["y", "z", "a", "b", "c"],
                                     [25, 26, 1, 2])

    pk.add_connection("y", "a", "shortcut")
    assert pk.get_path("y", "a") == (["y", "a"], ["shortcut"])
    assert pk.get_path("y", "c") == (["y", "a", "b", "c"], ["shortcut", 1, 2])

    # No automatic self-connections
    assert_raises(ValueError, pk.get_path, "a", "a")
    pk.add_connection("a", "a", "null")
    assert pk.get_path("a", "a") == (["a", "a"], ["null"])
    # We don't end up with ["null", 1] or anything
    assert pk.get_path("a", "b") == (["a", "b"], [1])
    # Even for paths added afterwards
    pk.add_connection("q", "a", "new")
    assert pk.get_path("q", "a") == (["q", "a"], ["new"])
    assert pk.get_path("q", "b") == (["q", "a", "b"], ["new", 1])

################################################################

_CONVERT_PATHS = _PathKeeper()

def _identity3d(x):
    x = np.asarray(x, dtype=float)
    if x.shape[-1] != 3:
        raise ValueError("Expected array with shape (..., 3)")
    return x

_CONVERT_PATHS.add_connection("XYZ", "XYZ", _identity3d)

_CONVERT_PATHS.add_connection("sRGB", "sRGB", _identity3d)
_CONVERT_PATHS.add_connection("sRGB", "sRGB-linear", sRGB_to_sRGB_linear)
_CONVERT_PATHS.add_connection("sRGB-linear", "sRGB", sRGB_linear_to_sRGB)

_CONVERT_PATHS.add_connection("sRGB-linear", "sRGB-linear", _identity3d)
_CONVERT_PATHS.add_connection("sRGB-linear", "XYZ", sRGB_linear_to_XYZ)
_CONVERT_PATHS.add_connection("XYZ", "sRGB-linear", XYZ_to_sRGB_linear)

_CONVERT_PATHS.add_connection("xyY", "xyY", _identity3d)
_CONVERT_PATHS.add_connection("XYZ", "xyY", XYZ_to_xyY)
_CONVERT_PATHS.add_connection("xyY", "XYZ", xyY_to_XYZ)

_CONVERT_PATHS.add_connection("CIELAB", "CIELAB", _identity3d)
_CONVERT_PATHS.add_connection("XYZ", "CIELAB", XYZ_to_CIELAB)
_CONVERT_PATHS.add_connection("CIELAB", "XYZ", CIELAB_to_XYZ)

# XX: CIELCh
# and J'/K M' h'

def _CIECAM02_to_XYZ(CIECAM02, viewing_conditions):
    return viewing_conditions.CIECAM02_to_XYZ(J=CIECAM02.J,
                                              C=CIECAM02.C,
                                              h=CIECAM02.h)

def _XYZ_to_CIECAM02(XYZ, viewing_conditions):
    return viewing_conditions.XYZ_to_CIECAM02(XYZ)

_CONVERT_PATHS.add_connection("CIECAM02", "CIECAM02", lambda x: x)
_CONVERT_PATHS.add_connection("XYZ", "CIECAM02", _XYZ_to_CIECAM02)
_CONVERT_PATHS.add_connection("CIECAM02", "XYZ", _CIECAM02_to_XYZ)

_CIECAM02_axes = set("JChQMsH")

class _Tag(object):
    def __init__(self, name):
        self._name = name
    def __repr__(self):
        return self._name
_CIECAM02_partial_tag = _Tag("(CIECAM02 subsets)")

def _CIECAM02_to_CIECAM02_partial(CIECAM02, axes):
    pieces = []
    for axis in axes:
        pieces.append(getattr(CIECAM02, axis)[..., np.newaxis])
    return np.concatenate(pieces, axis=-1)

def _CIECAM02_partial_to_XYZ(partial, viewing_conditions, axes):
    partial = np.asarray(partial, dtype=float)
    kwargs = {}
    if partial.shape[-1] != len(axes):
        raise ValueError("shape mismatch: last dimension of color array is "
                         "%s, but need %s for %r"
                         % partial.shape[-1], len(axes), axes)
    for i, coord in enumerate(axes):
        kwargs[coord] = partial[..., i]
    return viewing_conditions.CIECAM02_to_XYZ(**kwargs)

# We do *not* provide any CIECAM02-partial <-> CIECAM02-partial converter
# This will be implicitly created by going
#   CIECAM02-partial -> XYZ -> CIECAM02 -> CIECAM02-partial
# which is the correct way to do it.
_CONVERT_PATHS.add_connection("CIECAM02", _CIECAM02_partial_tag,
                              _CIECAM02_to_CIECAM02_partial)
_CONVERT_PATHS.add_connection(_CIECAM02_partial_tag, "XYZ",
                              _CIECAM02_partial_to_XYZ)

# Special case to give CAM02-UCS and friends a route
def _JMh_to_XYZ(JMh, viewing_conditions):
    return viewing_conditions.CIECAM02_to_XYZ(J=JMh[..., 0],
                                              M=JMh[..., 1],
                                              h=JMh[..., 2])

def _CIECAM02_to_JMh(CIECAM02):
    return _CIECAM02_to_CIECAM02_partial(CIECAM02, "JMh")

_CONVERT_PATHS.add_connection("JMh", "JMh", _identity3d)
_CONVERT_PATHS.add_connection("JMh", "XYZ", _JMh_to_XYZ)
_CONVERT_PATHS.add_connection("CIECAM02", "JMh", _CIECAM02_to_JMh)

def _LuoEtAl2006_to_JMh(JKapbp, uniform_space):
    return uniform_space.JKapbp_to_JMh(CAM02)

def _JMh_to_LuoEtAl2006(JMh, uniform_space):
    return uniform_space.JMh_to_JKapbp(JMh)

_LuoEtAl2006_tag = _Tag("(CAM02-UCS-like spaces)")

_CONVERT_PATHS.add_connection("JMh", _LuoEtAl2006_tag, _JMh_to_LuoEtAl2006)
_CONVERT_PATHS.add_connection(_LuoEtAl2006_tag, "JMh", _LuoEtAl2006_to_JMh)

def _tag_for_space(space):
    # We could use _CIECAM02_partial_tag for JMh as well, but it would make
    # explicit convert_cspace(..., "JMh", "JKapbp") pointlessly inefficient --
    # we'd end up going
    #   JMh (as partial) -> XYZ -> CIECAM02 -> JMh (as builtin hack) -> JKapbp
    if _CIECAM02_axes.issuperset(space) and space != "JMh":
        return _CIECAM02_partial_tag
    elif isinstance(space, LuoEtAlUniformSpace):
        return _LuoEtAl2006_tag
    else:
        return space

_ALIASES = {
    "CAM02-UCS": CAM02_UCS,
    "CAM02-LCD": CAM02_LCD,
    "CAM02-SCD": CAM02_SCD,
}

def convert_cspace(arr, start, end,
                   viewing_conditions=ViewingConditions.sRGB):
    start = _ALIASES.get(start, start)
    end = _ALIASES.get(end, end)

    start_tag = _tag_for_space(start)
    end_tag = _tag_for_space(end)

    nodes, converters = _CONVERT_PATHS.get_path(start_tag, end_tag)
    current = arr
    for i, converter in enumerate(converters):
        if converter in (XYZ_to_CIELAB, CIELAB_to_XYZ):
            current = converter(current, XYZ_w=viewing_conditions.XYZ_w)
        elif converter in (_CIECAM02_to_XYZ, _XYZ_to_CIECAM02, _JMh_to_XYZ):
            current = converter(current,
                                viewing_conditions=viewing_conditions)
        elif converter is _CIECAM02_to_CIECAM02_partial:
            # Our conversion graph is set up so that the CIECAM02 partial
            # conversions are a dead-end -- you never pass through them on the
            # way to anywhere else, only if you are starting and/or ending
            # there. This lets us recover the actual requested axes.
            assert i == len(converters) - 1
            current = converter(current, axes=end)
        elif converter is _CIECAM02_partial_to_XYZ:
            assert i == 0
            axes = start
            current = converter(current,
                                viewing_conditions=viewing_conditions,
                                axes=axes)
        elif converter is _LuoEtAl2006_to_JMh:
            current = converter(current, start)
        elif converter is _JMh_to_LuoEtAl2006:
            current = converter(current, end)
        else:
            current = converter(current)

    return current

def test_convert_cspace_long_paths():
    from .gold_values import sRGB_xyY_gold
    check_conversion(lambda x: convert_cspace(x, "sRGB", "xyY"),
                     lambda y: convert_cspace(y, "xyY", "sRGB"),
                     sRGB_xyY_gold,
                     a_min=0, a_max=1,
                     b_min=0, b_max=[1, 1, 100])

    from .gold_values import sRGB_CIELAB_gold_D65
    check_conversion(lambda x: convert_cspace(x, "sRGB", "CIELAB"),
                     lambda y: convert_cspace(y, "CIELAB", "sRGB"),
                     sRGB_CIELAB_gold_D65,
                     a_min=0, a_max=1,
                     b_min=[10, -30, 30], b_max=[90, 30, 30],
                     # Ridiculously low precision, but both Lindbloom and
                     # Grace's calculators have rounding errors in both the
                     # CIELAB coefficients and the sRGB matrices.
                     gold_rtol=1.5e-2)

    # Makes sure that CIELAB conversions are sensitive to whitepoint
    from .gold_values import XYZ_CIELAB_gold_D50
    from .ciecam02 import Illuminant
    vc = ViewingConditions(XYZ_w=Illuminant.D50,
                           Y_b=ViewingConditions.sRGB.Y_b,
                           L_A=ViewingConditions.sRGB.L_A,
                           surround=ViewingConditions.sRGB.surround)
    check_conversion(lambda x:
                     convert_cspace(x, "XYZ", "CIELAB",
                                    viewing_conditions=vc),
                     lambda y:
                     convert_cspace(y, "CIELAB", "XYZ",
                                    viewing_conditions=vc),
                     XYZ_CIELAB_gold_D50,
                     b_min=[10, -30, 30], b_max=[90, 30, 30])

    from .gold_values import XYZ_CIECAM02_gold
    for t in XYZ_CIECAM02_gold:
        # Check full-fledged CIECAM02 conversions
        xyY = convert_cspace(t.XYZ, "XYZ", "xyY")
        CIECAM02_got = convert_cspace(xyY, "xyY", "CIECAM02",
                                      viewing_conditions=t.vc)
        for i in range(len(CIECAM02_got)):
            assert np.allclose(CIECAM02_got[i], t.expected[i], atol=1e-5)
        xyY_got = convert_cspace(CIECAM02_got, "CIECAM02", "xyY",
                                 viewing_conditions=t.vc)
        assert np.allclose(xyY_got, xyY)

        # Check partial CIECAM02 conversions
        def stacklast(*arrs):
            arrs = [np.asarray(arr)[..., np.newaxis] for arr in arrs]
            return np.concatenate(arrs, axis=-1)
        JCh = stacklast(t.expected.J, t.expected.C, t.expected.h)
        xyY_got2 = convert_cspace(JCh, "JCh", "xyY", viewing_conditions=t.vc)
        assert np.allclose(xyY_got2, xyY)

        JCh_got = convert_cspace(xyY, "xyY", "JCh", viewing_conditions=t.vc)
        assert np.allclose(JCh_got, JCh, rtol=1e-4)

        # Check partial->partial CIECAM02
        QMH = stacklast(t.expected.Q, t.expected.M, t.expected.H)
        JCh_got2 = convert_cspace(QMH, "QMH", "JCh", viewing_conditions=t.vc)
        assert np.allclose(JCh_got2, JCh, rtol=1e-4)

        # And check the JMh special case too
        JMh = stacklast(t.expected.J, t.expected.M, t.expected.h)
        assert np.allclose(convert_cspace(JMh, "JMh", "JCh",
                                          viewing_conditions=t.vc),
                           JCh, rtol=1e-4)

        assert np.allclose(convert_cspace(JCh, "JCh", "JMh",
                                          viewing_conditions=t.vc),
                           JMh, rtol=1e-4)

        # TODO: test JKapbp
