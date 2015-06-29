# This file is part of pycam02ucs
# Copyright (C) 2015 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np
from collections import defaultdict

from .util import stacklast
from .testing import check_conversion
from .basics import (sRGB_to_sRGB_linear, sRGB_linear_to_sRGB,
                     sRGB_linear_to_XYZ100, XYZ100_to_sRGB_linear,
                     XYZ100_to_xyY100, xyY100_to_XYZ100,
                     XYZ100_to_CIELAB, CIELAB_to_XYZ100)

from .ciecam02 import CIECAM02Space
from .luoetal2006 import (LuoEtAl2006UniformSpace,
                          CAM02UCS, CAM02LCD, CAM02SCD)
from .cvd import machado_et_al_2009_matrix

from .transform_graph import Edge, MATCH, ANY, TransformGraph

__all__ = ["cspace_converter", "convert_cspace"]

################################################################

EDGES = []

def pair(a, b, a2b, b2a):
    if isinstance(a, str):
        a = {"name": a}
    if isinstance(b, str):
        b = {"name": b}
    return [Edge(a, b, a2b), Edge(b, a, b2a)]

EDGES += pair("sRGB", "sRGB255",
              lambda sRGB: np.asarray(sRGB) * 255.0,
              lambda sRGB255: np.asarray(sRGB255) / 255.0)

def _apply_rgb_mat(mat, rgb):
    return np.einsum("...ij,...j->...i", mat, rgb)

def _CVD_forward(sRGB, cvd_type, severity):
    mat = machado_et_al_2009_matrix(cvd_type, severity)
    return _apply_rgb_mat(mat, sRGB)

def _CVD_inverse(sRGB, cvd_type, severity):
    forward = machado_et_al_2009_matrix(cvd_type, severity)
    return _apply_rgb_mat(np.linalg.inv(forward), sRGB)

EDGES += pair({"name": "sRGB+CVD", "cvd_type": MATCH, "severity": MATCH},
              {"name": "sRGB-linear+CVD", "cvd_type": MATCH, "severity": MATCH},
              sRGB_to_sRGB_linear, sRGB_linear_to_sRGB)

EDGES += pair({"name": "sRGB-linear+CVD", "cvd_type": ANY, "severity": ANY},
              "sRGB-linear",
              _CVD_forward, _CVD_inverse)

EDGES += pair("sRGB", "sRGB-linear", sRGB_to_sRGB_linear, sRGB_linear_to_sRGB)

EDGES += pair("sRGB-linear", "XYZ100", sRGB_linear_to_XYZ100, XYZ100_to_sRGB_linear)

EDGES += pair("XYZ100", "XYZ1",
              lambda XYZ100: np.asarray(XYZ100) / 100.,
              lambda XYZ1: np.asarray(XYZ1) * 100.0)

EDGES += pair("XYZ100", "xyY100", XYZ100_to_xyY100, xyY100_to_XYZ100)

EDGES += pair("xyY100", "xyY1",
              lambda xyY100: np.asarray(xyY100) / [1.0, 1.0, 100.0],
              lambda xyY1: np.asarray(xyY1) * [1.0, 1.0, 100.0])

EDGES += pair("XYZ100", {"name": "CIELAB", "XYZ100_w": ANY},
              XYZ100_to_CIELAB, CIELAB_to_XYZ100)

def _XYZ100_to_CIECAM02(XYZ100, ciecam02_space):
    return ciecam02_space.XYZ100_to_CIECAM02(XYZ100)

def _CIECAM02_to_XYZ100(CIECAM02, ciecam02_space):
    return ciecam02_space.CIECAM02_to_XYZ100(J=CIECAM02.J,
                                              C=CIECAM02.C,
                                              h=CIECAM02.h)

EDGES += pair("XYZ100", {"name": "CIECAM02", "ciecam02_space": ANY},
              _XYZ100_to_CIECAM02, _CIECAM02_to_XYZ100)

_CIECAM02_axes = set("JChQMsH")

def _CIECAM02_to_CIECAM02_subset(CIECAM02, ciecam02_space, axes):
    pieces = []
    for axis in axes:
        pieces.append(getattr(CIECAM02, axis))
    return stacklast(*pieces)

def _CIECAM02_subset_to_XYZ100(subset, ciecam02_space, axes):
    subset = np.asarray(subset, dtype=float)
    kwargs = {}
    if subset.shape[-1] != len(axes):
        raise ValueError("shape mismatch: last dimension of color array is "
                         "%s, but need %s for %r"
                         % subset.shape[-1], len(axes), axes)
    for i, coord in enumerate(axes):
        kwargs[coord] = subset[..., i]
    return ciecam02_space.CIECAM02_to_XYZ100(**kwargs)

# We do *not* provide any CIECAM02-subset <-> CIECAM02-subset converter
# This will be implicitly created by going
#   CIECAM02-subset -> XYZ100 -> CIECAM02 -> CIECAM02-subset
# which is the correct way to do it.
EDGES += [
    Edge({"name": "CIECAM02",
          "ciecam02_space": MATCH},
         {"name": "CIECAM02-subset",
          "ciecam02_space": MATCH, "axes": ANY},
         _CIECAM02_to_CIECAM02_subset),
    Edge({"name": "CIECAM02-subset",
          "ciecam02_space": ANY, "axes": ANY},
         {"name": "XYZ100"},
         _CIECAM02_subset_to_XYZ100),
    ]

def _JMh_to_LuoEtAl2006(JMh, ciecam02_space, luoetal2006_space, axes):
    return luoetal2006_space.JMh_to_Jpapbp(JMh)

def _LuoEtAl2006_to_JMh(Jpapbp, ciecam02_space, luoetal2006_space, axes):
    return luoetal2006_space.Jpapbp_to_JMh(Jpapbp)

EDGES += pair({"name": "CIECAM02-subset",
                 "ciecam02_space": MATCH,
                 "axes": "JMh"},
              {"name": "J'a'b'",
                 "ciecam02_space": MATCH,
                 "luoetal2006_space": ANY},
              _JMh_to_LuoEtAl2006, _LuoEtAl2006_to_JMh)

GRAPH = TransformGraph(EDGES,
                       [["sRGB", "sRGB+CVD"],
                        ["sRGB-linear", "sRGB-linear+CVD"]])

ALIASES = {
    "CAM02-UCS": CAM02UCS,
    "CAM02-LCD": CAM02LCD,
    "CAM02-SCD": CAM02SCD,
    "CIECAM02": CIECAM02Space.sRGB,
    "CIELAB": {"name": "CIELAB", "XYZ100_w": CIECAM02Space.sRGB.XYZ100_w},
}

def norm_cspace_id(cspace):
    try:
        cspace = ALIASES[cspace]
    except (KeyError, TypeError):
        pass
    if isinstance(cspace, str):
        if _CIECAM02_axes.issuperset(cspace):
            return {"name": "CIECAM02-subset",
                    "ciecam02_space": CIECAM02Space.sRGB,
                    "axes": cspace}
        else:
            return {"name": cspace}
    elif isinstance(cspace, CIECAM02Space):
        return {"name": "CIECAM02",
                "ciecam02_space": cspace}
    elif isinstance(cspace, LuoEtAl2006UniformSpace):
        return {"name": "J'a'b'",
                "ciecam02_space": CIECAM02Space.sRGB,
                "luoetal2006_space": cspace}
    elif isinstance(cspace, dict):
        if cspace["name"] in ALIASES:
            base = ALIASES[cspace["name"]]
            if isinstance(base, dict) and base["name"] == cspace["name"]:
                # avoid infinite recursion
                return cspace
            else:
                base = norm_cspace_id(base)
                base = dict(base)
                del cspace["name"]
                base.update(cspace)
                return base
        return cspace
    else:
        raise ValueError("unrecognized color space %r" % (cspace,))

def cspace_converter(start, end):
    start = norm_cspace_id(start)
    end = norm_cspace_id(end)
    return GRAPH.get_transform(start, end)

def convert_cspace(arr, start, end):
    converter = cspace_converter(start, end)
    return converter(arr)

def check_convert_cspace(source_cspace, target_cspace, gold, **kwargs):
    def forward(source_values):
        return convert_cspace(source_values, source_cspace, target_cspace)
    def reverse(target_values):
        return convert_cspace(target_values, target_cspace, source_cspace)
    check_conversion(forward, reverse, gold, **kwargs)

def test_convert_cspace_long_paths():
    from .gold_values import sRGB_xyY100_gold
    check_convert_cspace("sRGB", "xyY100", sRGB_xyY100_gold,
                         a_min=0, a_max=1,
                         b_min=0, b_max=[1, 1, 100])

    from .gold_values import sRGB_xyY1_gold
    check_convert_cspace("sRGB", "xyY1", sRGB_xyY1_gold,
                         a_min=0, a_max=1,
                         b_min=0, b_max=1)

    from .gold_values import XYZ1_sRGB255_gold
    check_convert_cspace("XYZ1", "sRGB255", XYZ1_sRGB255_gold,
                         a_min=0, a_max=1,
                         b_min=0, b_max=255)

    from .gold_values import sRGB_CIELAB_gold_D65
    check_convert_cspace("sRGB", "CIELAB", sRGB_CIELAB_gold_D65,
                         a_min=0, a_max=1,
                         b_min=[10, -30, 30], b_max=[90, 30, 30],
                         # Ridiculously low precision, but both Lindbloom and
                         # Grace's calculators have rounding errors in both the
                         # CIELAB coefficients and the sRGB matrices.
                         gold_rtol=1.5e-2)

    # Makes sure that CIELAB conversions are sensitive to whitepoint
    from .gold_values import XYZ100_CIELAB_gold_D50
    check_convert_cspace("XYZ100", {"name": "CIELAB", "XYZ100_w": "D50"},
                         XYZ100_CIELAB_gold_D50,
                         b_min=[10, -30, 30], b_max=[90, 30, 30])

    from .gold_values import XYZ100_CIECAM02_gold
    for t in XYZ100_CIECAM02_gold:
        # Check full-fledged CIECAM02 conversions
        xyY100 = convert_cspace(t.XYZ100, "XYZ100", "xyY100")
        CIECAM02_got = convert_cspace(xyY100, "xyY100", t.vc)
        for i in range(len(CIECAM02_got)):
            assert np.allclose(CIECAM02_got[i], t.expected[i], atol=1e-5)
        xyY100_got = convert_cspace(CIECAM02_got, t.vc, "xyY100")
        assert np.allclose(xyY100_got, xyY100)

        # Check subset CIECAM02 conversions
        def subset(axes):
            return {"name": "CIECAM02-subset",
                    "axes": axes, "ciecam02_space": t.vc}
        JCh = stacklast(t.expected.J, t.expected.C, t.expected.h)
        xyY100_got2 = convert_cspace(JCh, subset("JCh"), "xyY100")
        assert np.allclose(xyY100_got2, xyY100)

        JCh_got = convert_cspace(xyY100, "xyY100", subset("JCh"))
        assert np.allclose(JCh_got, JCh, rtol=1e-4)

        # Check subset->subset CIECAM02
        # This needs only plain arrays so we can use check_convert_cspace
        QMH = stacklast(t.expected.Q, t.expected.M, t.expected.H)
        check_convert_cspace(subset("JCh"), subset("QMH"),
                             [(JCh, QMH)],
                             a_max=[100, 100, 360],
                             b_max=[100, 100, 400])

    # Check that directly transforming between two different viewing
    # conditions works.
    # This exploits the fact that the first two entries in our gold vector
    # have the same XYZ100.
    t1 = XYZ100_CIECAM02_gold[0]
    t2 = XYZ100_CIECAM02_gold[1]

    # "If we have a color that looks like t1.expected under viewing conditions
    # t1, then what does it look like under viewing conditions t2?"
    got2 = convert_cspace(t1.expected, t1.vc, t2.vc)
    for i in range(len(got2)):
        assert np.allclose(got2[i], t2.expected[i], atol=1e-5)

    JCh1 = stacklast(t1.expected.J, t1.expected.C, t1.expected.h)
    JCh2 = stacklast(t2.expected.J, t2.expected.C, t2.expected.h)
    JCh_space1 = {"name": "CIECAM02-subset", "axes": "JCh",
                  "ciecam02_space": t1.vc}
    JCh_space2 = {"name": "CIECAM02-subset", "axes": "JCh",
                  "ciecam02_space": t2.vc}
    check_convert_cspace(JCh_space1, JCh_space2,
                         [(JCh1, JCh2)],
                         a_max=[100, 100, 360],
                         b_max=[100, 100, 360])

    # J'a'b'
    from .gold_values import JMh_to_CAM02UCS_silver
    check_convert_cspace("JMh", "CAM02-UCS", JMh_to_CAM02UCS_silver,
                         a_max=[100, 100, 360],
                         b_min=[0, -30, -30], b_max=[100, 30, 30])
    from .gold_values import JMh_to_CAM02LCD_silver
    check_convert_cspace("JMh", "CAM02-LCD", JMh_to_CAM02LCD_silver,
                         a_max=[100, 100, 360],
                         b_min=[0, -30, -30], b_max=[100, 30, 30])
    from .gold_values import JMh_to_CAM02SCD_silver
    check_convert_cspace("JMh", "CAM02-SCD", JMh_to_CAM02SCD_silver,
                         a_max=[100, 100, 360],
                         b_min=[0, -30, -30], b_max=[100, 30, 30])
