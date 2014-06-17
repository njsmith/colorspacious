# Stub code for potentially switching to using LCMS for CIECAM02 calculations;
# not working or used at present.

cdef extern from "lcms.h":
   #...

def sRGB_to_XYZ(R, G, B):
    sRGB_prof = cmsCreate_sRGBProfile()
    XYZ_prof = cmsCreateXYZProfile()
    xform = cmsCreateTransform(sRGB_prof, TYPE_RGB_DBL,
                               XYZ_prof, TYPE_XYZ_DBL,
                               INTENT_PERCEPTUAL, cmsFLAGS_NOTPRECALC)
    double RGB[3]
    cmsCIEXYZ XYZ
    cmsDoTransform(xform, RGB, &XYZ, 1)

def XYZ_to_CIECAM02(X, Y, Z):
    cmsViewingConditions vc
    vc.whitePoint = (95.05, 100.00, 108.88)
    vc.Yb = 20.0
    vc.La = 4.0
    vc.surround = AVG_SURROUND
    # full adaptation
    vc.D_value = 1.0

    # need to calculate M from J as
