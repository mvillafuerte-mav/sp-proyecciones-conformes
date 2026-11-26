"""Mercator (Paso 1) y Cónica Conforme de Lambert (Paso 2)."""
from .constantes import R_M, PHI1, PHI2, latitud_isometrica


def mercator(lat, lon, R=R_M):
    """f(w) = R w  ->  x = R λ, y = R q(φ). Devuelve (x, y) en metros."""
    # TODO Persona 1/4
    raise NotImplementedError


def lambert_constantes(phi1=PHI1, phi2=PHI2):
    """Devuelve (n, F) para que σ(φ1) = σ(φ2) = 1."""
    # TODO Persona 1
    raise NotImplementedError


def lambert(lat, lon, phi1=PHI1, phi2=PHI2, R=R_M):
    """f(w) = R F exp(i n w). Devuelve (x, y) en metros."""
    # TODO Persona 1/4
    raise NotImplementedError
