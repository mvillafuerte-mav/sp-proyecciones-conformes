"""Mercator (Paso 1) y Cónica Conforme de Lambert (Paso 2), modelo esférico.

Las funciones reciben lat, lon en grados (números o arreglos) y devuelven
(x, y) en metros. Por dentro todo se hace con la coordenada compleja
w = λ + i q(φ), igual que en la wiki.
"""
import numpy as np

from .constantes import (R_M, PHI1, PHI2, LAT0_DEG, LON0_DEG, X0_M, Y0_M,
                         latitud_isometrica, coordenada_w)


def _w(lat, lon, lon0):
    # diferencia de longitudes llevada a [-180, 180)
    dlon = (np.asarray(lon) - lon0 + 180.0) % 360.0 - 180.0
    return coordenada_w(np.radians(lat), np.radians(dlon))


def mercator(lat, lon, lon0=0.0, R=R_M):
    """f(w) = R w  ->  x = R λ, y = R q(φ)."""
    z = R * _w(lat, lon, lon0)
    return z.real, z.imag


def lambert_constantes(phi1=PHI1, phi2=PHI2):
    """(n, F) para que σ(φ1) = σ(φ2) = 1. Paralelos en radianes.

    Si φ1 = φ2 (cono tangente) se usa el límite n = sen φ1.
    """
    q1, q2 = latitud_isometrica(phi1), latitud_isometrica(phi2)
    if np.isclose(phi1, phi2):
        n = np.sin(phi1)
    else:
        n = np.log(np.cos(phi1) / np.cos(phi2)) / (q2 - q1)
    F = np.cos(phi1) * np.exp(n * q1) / n
    return n, F


def lambert(lat, lon, phi1=PHI1, phi2=PHI2, lat0=LAT0_DEG, lon0=LON0_DEG,
            x0=X0_M, y0=Y0_M, R=R_M):
    """Cónica conforme de Lambert, f(w) = R F exp(i n w).

    Por defecto usa los parámetros de la Lambert de INEGI (EPSG:6372), pero
    sobre la esfera. La salida se rota -90° y se traslada para que el norte
    quede hacia arriba y el origen coincida con el de INEGI; eso no cambia
    la conformidad ni el factor de escala.
    """
    n, F = lambert_constantes(phi1, phi2)
    z = R * F * np.exp(1j * n * _w(lat, lon, lon0))
    rho0 = R * F * np.exp(-n * latitud_isometrica(np.radians(lat0)))
    z = -1j * z + 1j * rho0 + (x0 + 1j * y0)
    return z.real, z.imag
