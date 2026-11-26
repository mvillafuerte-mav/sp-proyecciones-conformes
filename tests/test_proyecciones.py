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


def _region_prueba():
    # cuadrilátero con una muesca, en lat/lon (grados)
    lon = np.array([-110, -90, -90, -100, -100, -110, -110.0])
    lat = np.array([15, 15, 30, 30, 24, 24, 15.0])
    return lat[:-1], lon[:-1]


def test_dirichlet_funcion_conocida():
    from proyecciones.dirichlet import Malla, w_de
    lat, lon = _region_prueba()
    ws = w_de(35.0, -120.0)
    exacta = lambda w: np.log(np.abs(w - ws))
    m = Malla(w_de(lat, lon), 0.004)
    u = m.resolver(exacta)
    assert np.max(np.abs(u - exacta(m.w))) < 1e-4


def test_optima_mejor_que_lambert_y_conforme():
    from proyecciones.dirichlet import ProyeccionOptima, verificar_escala
    from proyecciones.distorsion import distorsion, sigma_lambert, malla_interior
    lat, lon = _region_prueba()
    p = ProyeccionOptima(lat, lon, h=0.004)
    la, lo = malla_interior(np.r_[lat, lat[0]], np.r_[lon, lon[0]], paso=0.5)
    assert p.delta < distorsion(sigma_lambert(la))
    e_sigma, e_ne, e_ang = verificar_escala(p, np.array([20.0, 27.0]),
                                            np.array([-95.0, -95.0]), 1000.0)
    assert e_sigma < 1e-3 and e_ne < 1e-3 and e_ang < 0.05


def test_kmeans_lloyd_grupos_separados():
    from proyecciones.efecto import kmeans_lloyd
    rng = np.random.default_rng(0)
    X = np.vstack([rng.normal(0, 1, (50, 2)), rng.normal(100, 1, (50, 2))])
    et = kmeans_lloyd(X, X[[0, 50]])
    assert (et[:50] == 0).all() and (et[50:] == 1).all()


def test_conversion_con_modelos_guardados():
    from proyecciones.conversion import latlon_a_xy
    from proyecciones.datos import haversine
    lat, lon = np.array([19.4326, 25.6866]), np.array([-99.1332, -100.3161])
    for p in ("recomendada", "optima", "lambert", "mercator"):
        x, y = latlon_a_xy(lat, lon, p)
        assert x.shape == (2,)
    x, y = latlon_a_xy(lat, lon)
    d = np.hypot(x[1] - x[0], y[1] - y[0])
    assert abs(d / haversine(lat[0], lon[0], lat[1], lon[1]) - 1) < 0.002
    for p in ("recomendada", "lambert", "mercator"):
        x0, y0 = latlon_a_xy(19.4326, -99.1332, p)
        assert type(x0) is float and type(y0) is float


def test_conversion_con_nan_y_formas():
    from proyecciones.conversion import latlon_a_xy
    for p in ("recomendada", "lambert"):
        x, y = latlon_a_xy([19.43, np.nan], [-99.13, -99.0], p)
        assert np.isfinite(x[0]) and np.isnan(x[1]) and np.isnan(y[1])
    x, y = latlon_a_xy(np.array([[19.4, 20.0]]), np.array([[-99.1, -100.0]]))
    assert x.shape == (1, 2)
