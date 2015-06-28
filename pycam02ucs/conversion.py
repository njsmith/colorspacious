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
from .luoetal2006 import (LuoEtAl2006UniformSpace,
                          CAM02UCS, CAM02LCD, CAM02SCD)

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

EDGES += pair("sRGB", "sRGB-linear", sRGB_to_sRGB_linear, sRGB_linear_to_sRGB)

EDGES += pair("sRGB-linear", "XYZ", sRGB_linear_to_XYZ, XYZ_to_sRGB_linear)

EDGES += pair("XYZ", "xyY", XYZ_to_xyY, xyY_to_XYZ)

EDGES += pair("XYZ", {"name": "CIELAB", "XYZ_w": ANY},
              XYZ_to_CIELAB, CIELAB_to_XYZ)

# XX: CIELCh
# and J'/K M' h'

def _XYZ_to_CIECAM02(XYZ, viewing_conditions):
    return viewing_conditions.XYZ_to_CIECAM02(XYZ)

def _CIECAM02_to_XYZ(CIECAM02, viewing_conditions):
    return viewing_conditions.CIECAM02_to_XYZ(J=CIECAM02.J,
                                              C=CIECAM02.C,
                                              h=CIECAM02.h)

EDGES += pair("XYZ", {"name": "CIECAM02", "viewing_conditions": ANY},
              _XYZ_to_CIECAM02, _CIECAM02_to_XYZ)

_CIECAM02_axes = set("JChQMsH")

def _CIECAM02_to_CIECAM02_subset(CIECAM02, viewing_conditions, axes):
    pieces = []
    for axis in axes:
        pieces.append(getattr(CIECAM02, axis)[..., np.newaxis])
    return np.concatenate(pieces, axis=-1)

def _CIECAM02_subset_to_XYZ(subset, viewing_conditions, axes):
    subset = np.asarray(subset, dtype=float)
    kwargs = {}
    if subset.shape[-1] != len(axes):
        raise ValueError("shape mismatch: last dimension of color array is "
                         "%s, but need %s for %r"
                         % subset.shape[-1], len(axes), axes)
    for i, coord in enumerate(axes):
        kwargs[coord] = subset[..., i]
    return viewing_conditions.CIECAM02_to_XYZ(**kwargs)

# We do *not* provide any CIECAM02-subset <-> CIECAM02-subset converter
# This will be implicitly created by going
#   CIECAM02-subset -> XYZ -> CIECAM02 -> CIECAM02-subset
# which is the correct way to do it.
EDGES += [
    Edge({"name": "CIECAM02",
          "viewing_conditions": MATCH},
         {"name": "CIECAM02-subset",
          "viewing_conditions": MATCH, "axes": ANY},
         _CIECAM02_to_CIECAM02_subset),
    Edge({"name": "CIECAM02-subset",
          "viewing_conditions": ANY, "axes": ANY},
         {"name": "XYZ"},
         _CIECAM02_subset_to_XYZ),
    ]

def _JMh_to_LuoEtAl2006(JMh, viewing_conditions, luoetal2006_space, axes):
    return luoetal2006_space.JMh_to_JKapbp(JMh)

def _LuoEtAl2006_to_JMh(JKapbp, viewing_conditions, luoetal2006_space, axes):
    return luoetal2006_space.JKapbp_to_JMh(CAM02)

EDGES += pair({"name": "CIECAM02-subset",
                 "viewing_conditions": MATCH,
                 "axes": "JMh"},
              {"name": "J'a'b'",
                 "viewing_conditions": MATCH,
                 "luoetal2006_space": ANY},
              _JMh_to_LuoEtAl2006, _LuoEtAl2006_to_JMh)

GRAPH = TransformGraph(EDGES)

ALIASES = {
    "CAM02-UCS": CAM02UCS,
    "CAM02-LCD": CAM02LCD,
    "CAM02-SCD": CAM02SCD,
    "CIECAM02": ViewingConditions.sRGB,
    "CIELAB": {"name": "CIELAB", "XYZ_w": ViewingConditions.sRGB.XYZ_w},
}

def norm_cspace_id(cspace):
    try:
        cspace = ALIASES[cspace]
    except (KeyError, TypeError):
        pass
    if isinstance(cspace, str):
        if _CIECAM02_axes.issuperset(cspace):
            return {"name": "CIECAM02-subset",
                    "viewing_conditions": ViewingConditions.sRGB,
                    "axes": cspace}
        else:
            return {"name": cspace}
    elif isinstance(cspace, ViewingConditions):
        return {"name": "CIECAM02",
                "viewing_conditions": cspace}
    elif isinstance(cspace, LuoEtAl2006UniformSpace):
        return {"name": "J'a'b'",
                "viewing_conditions": ViewingConditions.sRGB,
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
    CIELAB_D50 = {"name": "CIELAB", "XYZ_w": Illuminant.D50}
    check_conversion(lambda x:
                     convert_cspace(x, "XYZ", CIELAB_D50),
                     lambda y:
                     convert_cspace(y, CIELAB_D50, "XYZ"),
                     XYZ_CIELAB_gold_D50,
                     b_min=[10, -30, 30], b_max=[90, 30, 30])

    from .gold_values import XYZ_CIECAM02_gold
    for t in XYZ_CIECAM02_gold:
        # Check full-fledged CIECAM02 conversions
        xyY = convert_cspace(t.XYZ, "XYZ", "xyY")
        CIECAM02_got = convert_cspace(xyY, "xyY", t.vc)
        for i in range(len(CIECAM02_got)):
            assert np.allclose(CIECAM02_got[i], t.expected[i], atol=1e-5)
        xyY_got = convert_cspace(CIECAM02_got, t.vc, "xyY")
        assert np.allclose(xyY_got, xyY)

        # Check subset CIECAM02 conversions
        def stacklast(*arrs):
            arrs = [np.asarray(arr)[..., np.newaxis] for arr in arrs]
            return np.concatenate(arrs, axis=-1)
        def subset(axes):
            return {"name": "CIECAM02-subset",
                    "axes": axes, "viewing_conditions": t.vc}
        JCh = stacklast(t.expected.J, t.expected.C, t.expected.h)
        xyY_got2 = convert_cspace(JCh, subset("JCh"), "xyY")
        assert np.allclose(xyY_got2, xyY)

        JCh_got = convert_cspace(xyY, "xyY", subset("JCh"))
        assert np.allclose(JCh_got, JCh, rtol=1e-4)

        # Check subset->subset CIECAM02
        QMH = stacklast(t.expected.Q, t.expected.M, t.expected.H)
        JCh_got2 = convert_cspace(QMH, subset("QMH"), subset("JCh"))
        assert np.allclose(JCh_got2, JCh, rtol=1e-4)

        # And check the JMh special case too
        JMh = stacklast(t.expected.J, t.expected.M, t.expected.h)
        assert np.allclose(convert_cspace(JMh, subset("JMh"), subset("JCh")),
                           JCh, rtol=1e-4)

        assert np.allclose(convert_cspace(JCh, subset("JCh"), subset("JMh")),
                           JMh, rtol=1e-4)

        # TODO: test JKapbp
        # TODO: check transforms between different viewing conditions
        #   the first two tests in our gold vector have the same XYZ
