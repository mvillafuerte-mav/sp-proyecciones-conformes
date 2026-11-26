"""Pruebas: comparar con pyproj (esfera R = 6371008.8) y caso armónico conocido."""
import numpy as np
from proyecciones.constantes import latitud_isometrica


def test_latitud_isometrica_ecuador():
    assert np.isclose(latitud_isometrica(0.0), 0.0)


def test_derivada_q_es_sec():
    phi = np.radians(20.0)
    h = 1e-6
    dq = (latitud_isometrica(phi + h) - latitud_isometrica(phi - h)) / (2 * h)
    assert np.isclose(dq, 1 / np.cos(phi))


# TODO Persona 2: Dirichlet con h = Re(e^w) en la frontera


def test_haversine_un_grado_de_latitud():
    from proyecciones.datos import haversine
    from proyecciones.constantes import R_M
    d = haversine(20.0, -100.0, 21.0, -100.0)
    assert np.isclose(d, R_M * np.pi / 180)


def test_haversine_simetrica():
    from proyecciones.datos import haversine
    a = haversine(19.43, -99.13, 25.69, -100.32)
    b = haversine(25.69, -100.32, 19.43, -99.13)
    assert np.isclose(a, b)


def _puntos_prueba():
    rng = np.random.default_rng(1)
    return rng.uniform(14.5, 32.7, 500), rng.uniform(-117, -86.7, 500)


def test_mercator_igual_a_proj():
    from pyproj import Proj
    from proyecciones.clasicas import mercator
    lat, lon = _puntos_prueba()
    x, y = mercator(lat, lon)
    X, Y = Proj("+proj=merc +R=6371008.8 +lon_0=0")(lon, lat)
    assert np.max(np.hypot(x - X, y - Y)) < 1e-6


def test_lambert_igual_a_proj():
    from pyproj import Proj
    from proyecciones.clasicas import lambert
    lat, lon = _puntos_prueba()
    x, y = lambert(lat, lon)
    p = Proj("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 "
             "+x_0=2500000 +y_0=0 +R=6371008.8")
    X, Y = p(lon, lat)
    assert np.max(np.hypot(x - X, y - Y)) < 1e-6


def test_lambert_escala_uno_en_paralelos():
    from proyecciones.distorsion import sigma_lambert
    assert np.allclose(sigma_lambert(np.array([17.5, 29.5])), 1.0)


def test_sigma_igual_a_factores_de_proj():
    from pyproj import Proj
    from proyecciones.distorsion import sigma_lambert
    lat, lon = _puntos_prueba()
    p = Proj("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 "
             "+x_0=2500000 +y_0=0 +R=6371008.8")
    k = p.get_factors(lon, lat).parallel_scale
    assert np.allclose(sigma_lambert(lat), k, rtol=1e-9)
