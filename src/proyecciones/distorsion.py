"""Factor de escala σ y distorsión δ(Ω) = max ln σ − min ln σ (Paso 3).

lat en grados. Para Mercator y Lambert σ depende solo de la latitud.
"""
import numpy as np

from .constantes import latitud_isometrica
from .clasicas import lambert_constantes


def sigma_mercator(lat):
    """σ = sec φ."""
    return 1.0 / np.cos(np.radians(lat))


def sigma_lambert(lat, n=None, F=None):
    """σ = n F exp(−n q(φ)) / cos φ. Sin n, F usa los paralelos de INEGI."""
    if n is None or F is None:
        n, F = lambert_constantes()
    phi = np.radians(lat)
    return n * F * np.exp(-n * latitud_isometrica(phi)) / np.cos(phi)


def distorsion(sigma):
    """δ = max ln σ − min ln σ sobre los valores dados."""
    ls = np.log(np.asarray(sigma))
    return float(ls.max() - ls.min())


def error_max_escala(delta):
    """Error máximo de escala en %, con la escala centrada: 100 (e^{δ/2} − 1).

    Supone que el mapa se reescala para que σ quede entre e^{-δ/2} y e^{δ/2}.
    """
    return 100.0 * (np.exp(delta / 2) - 1.0)


def malla_interior(lat_contorno, lon_contorno, paso=0.05):
    """Puntos de una malla regular (en grados) dentro del contorno, más el
    propio contorno. Devuelve (lat, lon)."""
    from shapely import contains_xy
    from shapely.geometry import Polygon

    pol = Polygon(np.column_stack([lon_contorno, lat_contorno]))
    lo = np.arange(lon_contorno.min(), lon_contorno.max() + paso, paso)
    la = np.arange(lat_contorno.min(), lat_contorno.max() + paso, paso)
    LO, LA = np.meshgrid(lo, la)
    dentro = contains_xy(pol, LO.ravel(), LA.ravel())
    lat = np.concatenate([LA.ravel()[dentro], lat_contorno])
    lon = np.concatenate([LO.ravel()[dentro], lon_contorno])
    return lat, lon


def resumen():
    """δ y errores de escala de Mercator y Lambert sobre México continental."""
    import pandas as pd
    from .datos import cargar_contorno

    lat_c, lon_c = cargar_contorno()
    lat, _ = malla_interior(lat_c, lon_c)
    filas = []
    for nombre, s in [("Mercator", sigma_mercator(lat)),
                      ("Lambert INEGI", sigma_lambert(lat))]:
        d = distorsion(s)
        filas.append({
            "proyeccion": nombre,
            "delta": d,
            "sigma_min": s.min(),
            "sigma_max": s.max(),
            "error_max_escala_centrada_%": error_max_escala(d),
            "error_max_escala_tal_cual_%": 100 * np.max(np.abs(s - 1)),
        })
    return pd.DataFrame(filas)


if __name__ == "__main__":
    from .datos import RAIZ

    df = resumen()
    print(df.to_string(index=False))
    out = RAIZ / "resultados" / "tablas" / "distorsion_clasicas.csv"
    df.to_csv(out, index=False)
    print(f"guardado en {out}")
