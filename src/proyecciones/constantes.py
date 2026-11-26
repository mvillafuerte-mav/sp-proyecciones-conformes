"""Constantes que usamos todos. Si cambias algo aquí, avisa en el grupo."""
import numpy as np

# Radio medio de la Tierra (modelo esférico)
R_KM = 6371.0088
R_M = R_KM * 1000.0

# Lambert de INEGI: EPSG:6372 (Mexico ITRF2008 / LCC), elipsoide GRS80.
# Ojo: la consigna dice EPSG:6361, pero ese código no existe en el registro
# EPSG; los parámetros que da corresponden a 6372 (y a 6362, ITRF92).
EPSG_LAMBERT_INEGI = 6372
PHI1_DEG = 17.5
PHI2_DEG = 29.5
PHI1 = np.radians(PHI1_DEG)
PHI2 = np.radians(PHI2_DEG)
LAT0_DEG = 12.0          # latitud de origen
LON0_DEG = -102.0        # meridiano central
X0_M = 2_500_000.0       # falso este
Y0_M = 0.0               # falso norte

# Extremos de latitud de México continental (verificados con el contorno del INEGI)
LAT_MIN_DEG = 14.532
LAT_MAX_DEG = 32.719


def latitud_isometrica(phi):
    """q(φ) = ln tan(π/4 + φ/2), con φ en radianes."""
    return np.log(np.tan(np.pi / 4 + phi / 2))


def coordenada_w(lat, lon):
    """w = λ + i q(φ). Recibe lat, lon en radianes."""
    return lon + 1j * latitud_isometrica(lat)
