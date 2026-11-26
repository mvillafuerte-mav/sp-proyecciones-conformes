"""Pruebas: comparar con pyproj (esfera R = 6371008.8) y caso armónico conocido."""
import numpy as np
from proyecciones.constantes import latitud_isometrica


def test_latitud_isometrica_ecuador():
    assert latitud_isometrica(0.0) == 0.0


def test_derivada_q_es_sec():
    phi = np.radians(20.0)
    h = 1e-6
    dq = (latitud_isometrica(phi + h) - latitud_isometrica(phi - h)) / (2 * h)
    assert np.isclose(dq, 1 / np.cos(phi))


# TODO Persona 4: mercator y lambert vs pyproj
# TODO Persona 2: Dirichlet con h = Re(e^w) en la frontera
