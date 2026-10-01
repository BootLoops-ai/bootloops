#!/usr/bin/env python3
"""io_util.py — (de)serialization helpers for the cli (mp matrices <-> json)."""
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST


def mat2s(M, digits=40):
    return [[(mp.nstr(M[i, j].real, digits), mp.nstr(mp.im(M[i, j]), digits))
             for j in range(M.cols)] for i in range(M.rows)]


def s2mat(L):
    r, c = len(L), len(L[0])
    M = mp.matrix(r, c)
    for i in range(r):
        for j in range(c):
            M[i, j] = mp.mpc(mp.mpf(L[i][j][0]), mp.mpf(L[i][j][1]))
    return M


def vec2s(v, digits=40):
    return [(mp.nstr(mp.mpc(v[i]).real, digits), mp.nstr(mp.im(mp.mpc(v[i])), digits))
            for i in range(len(v))]


def coeffs2s(Dc, digits=60):
    return [(mp.nstr(mp.mpc(c).real, digits), mp.nstr(mp.im(mp.mpc(c)), digits))
            for c in Dc]


def s2coeffs(L):
    return [mp.mpf(re) if mp.mpf(im) == 0 else mp.mpc(mp.mpf(re), mp.mpf(im))
            for re, im in L]
