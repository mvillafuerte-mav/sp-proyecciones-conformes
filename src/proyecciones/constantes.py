"""Constantes que usamos todos. Si cambias algo aquí, avisa en el grupo."""
import numpy as np

# Radio medio de la Tierra (modelo esférico)
R_KM = 6371.0088
R_M = R_KM * 1000.0

# Lambert de INEGI (EPSG:6361). Revisar origen y falso este en epsg.io
PHI1_DEG = 17.5
PHI2_DEG = 29.5
PHI1 = np.radians(PHI1_DEG)
PHI2 = np.radians(PHI2_DEG)

# Extremos aproximados de latitud de México continental (verificar con el contorno)
LAT_MIN_DEG = 14.53
LAT_MAX_DEG = 32.72


def latitud_isometrica(phi):
    """q(φ) = ln tan(π/4 + φ/2), con φ en radianes."""
    return np.log(np.tan(np.pi / 4 + phi / 2))


def coordenada_w(lat, lon):
    """w = λ + i q(φ). Recibe lat, lon en radianes."""
    return lon + 1j * latitud_isometrica(lat)
