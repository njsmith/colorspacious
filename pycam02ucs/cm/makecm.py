# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

# A simple approach to generating "optimal" colormaps
#
# Our optimality criteria:
# given a desired set of J and h values, choose M values such that:
# - we stay within the sRGB gamut
#     penalize for distance outside gamut, then clip and continue
# - we maximize arc length
# - we minimize some sort of curvature metric (sum of 3rd order differences?)
# - we stay as close as possible to perceptually uniform (all steps should be
#   the same size in CAM02-UCS space. Or maybe CAM02-SCD, dunno, should try
#   both.)

# Useful to know:
#   focal red: h = 20.14
#   focal yellow: h = 90.00
#   focal green: h = 164.25
#   focal blue: h = 237.53
