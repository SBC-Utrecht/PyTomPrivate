from numba import jit
import numpy as xp

@jit(nopython=False)
def vol2sf(data, res, r, b, m_x, m_y, m_z):
    """volTOsf(vol vol, double const r, int const b, float const m_x, float const m_y, float const m_z) -> vol"""
    from pytom.agnostic.interpolation import splineInterpolation

    for jj in range(0, 2*b):
        for kk in range(0, 2*b):
            the = xp.pi * (2 * jj + 1) / (4 * b)
            phi = xp.pi * kk / b
            x = r * xp.cos(phi) * xp.sin(the)
            y = r * xp.sin(phi) * xp.sin(the)
            z = r * xp.cos(the)
            res[jj * 2 * b + kk,0,0] = splineInterpolation(data, x + m_x, y + m_y, z + m_z)

    return res

@jit(nopython=False)
def fvol2sf(data, res, r, b):
    """fvolTOsf(vol_comp vol, double const r, int const b) -> vol_comp"""
    from pytom.agnostic.interpolation import splineInterpolation

    m_x, m_y, m_z = data.sizeX()//2, data.sizeY()//2, data.sizeZ()//2
    dr = data.real
    di = data.imag

    for i in range(2 * b * b):
        jj = i // (2 * b)
        kk = i % (2 * b)
        the = xp.pi * (2 * jj + 1) / (4 * b)
        phi = xp.pi * kk / b
        x = r * xp.cos(phi) * xp.sin(the)
        y = r * xp.sin(phi) * xp.sin(the)
        z = r * xp.cos(the)

        cr = splineInterpolation(dr, x + m_x, y + m_y, z + m_z)
        ci = splineInterpolation(di, x + m_x, y + m_y, z + m_z)

        res[jj * 2 * b + kk, 0, 0] = complex(cr,ci)
        res[(2 * b - jj - 1) * 2 * b + (kk+b if (kk<b) else kk-b), 0, 0] = complex(cr,-ci)

    return res

@jit(nopython=False)
def apply_shift_sf(data, res, b, r, M_PI, sizex, sizey, sizez, shiftX, shiftY, shiftZ):
    for i in range(2 * b * b):
        jj = i // (2 * b)
        kk = i % (2 * b)
        the = M_PI * (2 * jj + 1) / (4 * b)
        phi = M_PI * kk / b
        x = r * xp.cos(phi) * xp.sin(the)
        y = r * xp.sin(phi) * xp.sin(the)
        z = r * xp.cos(the)

        shift = complex(xp.cos(2 * M_PI / sizex * x * shiftX), -xp.sin(2 * M_PI / sizex * x * shiftX)) * \
                complex(xp.cos(2 * M_PI / sizey * y * shiftY), -xp.sin(2 * M_PI / sizey * y * shiftY)) * \
                complex(xp.cos(2 * M_PI / sizez * z * shiftZ), -xp.sin(2 * M_PI / sizez * z * shiftZ))

        c = data[jj * 2 * b + kk, 0, 0] * shift

        res[jj * 2 * b + kk, 0, 0] = c
        res[(2 * b - jj - 1) * 2 * b + (kk+b if (kk<b) else kk-b), 0, 0] = xp.conj(c)

    #return res

