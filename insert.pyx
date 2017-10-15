import numpy as np 
cimport numpy as np 
from libc.math cimport sin, cos

def insert(np.ndarray[np.float64_t, ndim = 2] slices, np.ndarray[np.float64_t, ndim = 1] thetas):

    thetas = thetas / 180 * np.pi
    cdef np.ndarray[np.float64_t, ndim = 1] slc
    cdef int s0, s1, sc, c0, c1, x, y
    cdef double xx, yy, fx, fy, ifx, ify, w, th
    cdef Py_ssize_t i, t


    s0 = slices.shape[0]
    s1 = slices.shape[1]
    assert(s0 == len(thetas))
    n0 = s1
    n1 = s1

    sc = s1 // 2
    c0 = n0 // 2
    c1 = n1 // 2

    cdef np.ndarray[np.float64_t, ndim = 2] model = np.zeros((n0, n1), dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 2] weight = np.zeros((n0, n1), dtype = np.float64)

    for i, th in enumerate(thetas):
        slc = slices[i]
        for t, v in enumerate(slc):
            x = 0
            y = t - sc

            xx = cos(th)*x - sin(th) * y
            yy = sin(th)*x + cos(th) * y

            xx += c0
            yy += c1

            if xx <= 0 or xx >= n0 - 1: continue
            if yy <= 0 or yy >= n1 - 1: continue

            xxx = <int>xx
            yyy = <int>yy

            fx = xx - xxx
            fy = yy - yyy
            ifx = 1. - fx
            ify = 1. - fy

            w = ifx * ify
            model[xxx, yyy] += w*v 
            weight[xxx, yyy] += w

            w = ifx * fy
            model[xxx, yyy + 1]  += w*v
            weight[xxx, yyy + 1] += w

            w = fx * fy
            model[xxx + 1, yyy + 1] += w*v
            weight[xxx + 1, yyy + 1] += w

            w = fx * ify
            model[xxx + 1, yyy] += w*v
            weight[xxx + 1, yyy] += w

    weight =  np.where(weight == 0., 1., weight)
    return model / weight



