# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import numpy as np
cimport numpy as np
np.import_array()
np.import_ufunc()

# We use the v0.3 code, because the v0.4 code totally changed the API in the
# .c file, but did not change the .h file, so it just... doesn't work at all.
cdef extern from "lib/ciecam02.h":
    cdef void xyz2jchqms_ciecam02(double *J, double *C, double *h,
                                  double *Q, double *M, double *s,
                                  double x, double y, double z,
                                  double xw, double yw, double zw,
                                  double yb, double la,
                                  double f, double c, double nc)

    cdef void jch2xyz_ciecam02(double *x, double *y, double *z,
                               double J, double C, double h,
                               double xw, double yw, double zw,
                               double yb, double la,
                               double f, double c, double nc)


# easier than trying to link it...
cdef extern from "lib/ciecam02.c":
    pass

cdef void XYZ_to_JMh_inner(char ** args, np.npy_intp * dimensions,
                           np.npy_intp * steps, void * data):
    cdef np.npy_intp i
    cdef char * X_arr = args[0]
    cdef np.npy_intp X_steps = steps[0]
    cdef char * Y_arr = args[1]
    cdef np.npy_intp Y_steps = steps[1]
    cdef char * Z_arr = args[2]
    cdef np.npy_intp Z_steps = steps[2]
    cdef char * J_arr = args[3]
    cdef np.npy_intp J_steps = steps[3]
    cdef char * M_arr = args[4]
    cdef np.npy_intp M_steps = steps[4]
    cdef char * h_arr = args[5]
    cdef np.npy_intp h_steps = steps[5]
    for i in range(dimensions[0]):
        xyz2jchqms_ciecam02(<double *> (J_arr + (i * J_steps)),
                            NULL,
                            <double *> (h_arr + (i * h_steps)),
                            NULL,
                            <double *> (M_arr + (i * M_steps)),
                            NULL,
                            (<double *> (X_arr + (i * X_steps)))[0] * 100,
                            (<double *> (Y_arr + (i * Y_steps)))[0] * 100,
                            (<double *> (Z_arr + (i * Z_steps)))[0] * 100,
                            # magic values copied from comments in ciecam02.h,
                            # allegedly appropriate for sRGB conversions
                            95.05, 100.00, 108.88,
                            20.0, 4.0,
                            1.0, 0.690, 1.0)

cdef np.PyUFuncGenericFunction * XYZ_to_JMh_loops = [&XYZ_to_JMh_inner]
cdef char * XYZ_to_JMh_types = [np.NPY_DOUBLE, np.NPY_DOUBLE, np.NPY_DOUBLE,
                                np.NPY_DOUBLE, np.NPY_DOUBLE, np.NPY_DOUBLE]
cdef void ** XYZ_to_JMh_data = [NULL]

cdef char * XYZ_to_JMh_docstring = (
"""Calculate CIECAM02 J, M, h values for given XYZ color.

J = lightness, M = colorfulness, h = hue angle (0-360).

Assumes the CIECAM02 "average" surround, with D65 whitepoint, 20% gray, and
"La value of 4 cd/m^2 to correspond to an ambient illumation of 64 lux", which
match the illuminant assumptions for sRGB.""")

XYZ_to_JMh = np.PyUFunc_FromFuncAndData(XYZ_to_JMh_loops,
                                        XYZ_to_JMh_data,
                                        XYZ_to_JMh_types,
                                        1, # number of loops
                                        3, # nin
                                        3, # nout
                                        np.PyUFunc_None, # reduction identity
                                        "XYZ_to_JMh",
                                        XYZ_to_JMh_docstring,
                                        # meaningless argument because numpy:
                                        0,
                                    )

################

cdef void JMh_to_XYZ_inner(char ** args, np.npy_intp * dimensions,
                           np.npy_intp * steps, void * data):
    cdef np.npy_intp i
    cdef char * J_arr = args[0]
    cdef np.npy_intp J_steps = steps[0]
    cdef char * M_arr = args[1]
    cdef np.npy_intp M_steps = steps[1]
    cdef char * h_arr = args[2]
    cdef np.npy_intp h_steps = steps[2]
    cdef char * X_arr = args[3]
    cdef np.npy_intp X_steps = steps[3]
    cdef char * Y_arr = args[4]
    cdef np.npy_intp Y_steps = steps[4]
    cdef char * Z_arr = args[5]
    cdef np.npy_intp Z_steps = steps[5]
    for i in range(dimensions[0]):
        jch2xyz_ciecam02(<double *> (X_arr + (i * X_steps)),
                         <double *> (Y_arr + (i * Y_steps)),
                         <double *> (Z_arr + (i * Z_steps)),
                         (<double *> (J_arr + (i * J_steps)))[0],
                         (<double *> (h_arr + (i * h_steps)))[0],
                         (<double *> (M_arr + (i * M_steps)))[0],
                         # magic values copied from comments in ciecam02.h,
                         # allegedly appropriate for sRGB conversions
                         95.05, 100.00, 108.88,
                         20.0, 4.0,
                         1.0, 0.690, 1.0)
        (<double *> (X_arr + (i * X_steps)))[0] /= 100
        (<double *> (Y_arr + (i * Y_steps)))[0] /= 100
        (<double *> (Z_arr + (i * Z_steps)))[0] /= 100

cdef np.PyUFuncGenericFunction * JMh_to_XYZ_loops = [&JMh_to_XYZ_inner]
cdef char * JMh_to_XYZ_types = [np.NPY_DOUBLE, np.NPY_DOUBLE, np.NPY_DOUBLE,
                                np.NPY_DOUBLE, np.NPY_DOUBLE, np.NPY_DOUBLE]
cdef void ** JMh_to_XYZ_data = [NULL]

cdef char * JMh_to_XYZ_docstring = (
"""Calculate XYZ values given CIECAM02 J, M, h values.

J = lightness, M = colorfulness, h = hue angle (0-360).

Assumes the CIECAM02 "average" surround, with D65 whitepoint, 20% gray, and
"La value of 4 cd/m^2 to correspond to an ambient illumation of 64 lux", which
match the illuminant assumptions for sRGB.""")

JMh_to_XYZ = np.PyUFunc_FromFuncAndData(JMh_to_XYZ_loops,
                                        JMh_to_XYZ_data,
                                        JMh_to_XYZ_types,
                                        1, # number of loops
                                        3, # nin
                                        3, # nout
                                        np.PyUFunc_None, # reduction identity
                                        "JMh_to_XYZ",
                                        JMh_to_XYZ_docstring,
                                        # meaningless argument because numpy:
                                        0,
                                    )
