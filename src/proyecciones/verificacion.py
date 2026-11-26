"""Verificación contra PROJ y efecto del elipsoide (sección 5 de la wiki).

    python -m proyecciones.verificacion

Guarda los números en resultados/tablas/verificacion_proj.csv.
"""
import numpy as np
import pandas as pd
from pyproj import CRS, Geod, Proj, Transformer

from .clasicas import lambert, mercator
from .datos import RAIZ, cargar_contorno, cargar_puntos, haversine
from .distorsion import distorsion, malla_interior, sigma_lambert, sigma_mercator

MERC_ESFERA = "+proj=merc +R=6371008.8 +lon_0=0"
LCC_ESFERA = ("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 "
              "+x_0=2500000 +y_0=0 +R=6371008.8")


def verificar(n_pares=20000, semilla=0):
    df = cargar_puntos()
    lat, lon = df["lat"].to_numpy(), df["lon"].to_numpy()
    res = {}

    # 1) contra PROJ en la esfera
    for nombre, f, p, s in [("mercator", mercator, Proj(MERC_ESFERA), sigma_mercator),
                            ("lambert", lambert, Proj(LCC_ESFERA), sigma_lambert)]:
        x, y = f(lat, lon)
        X, Y = p(lon, lat)
        res[f"{nombre}_dif_max_posicion_m"] = np.max(np.hypot(x - X, y - Y))
        k = p.get_factors(lon, lat).parallel_scale
        res[f"{nombre}_dif_max_sigma"] = np.max(np.abs(s(lat) - k))

    # 2) Lambert esférica contra EPSG:6372 sobre GRS80
    t = Transformer.from_crs(4326, 6372, always_xy=True)
    Xe, Ye = t.transform(lon, lat)
    x, y = lambert(lat, lon)
    d = np.hypot(x - Xe, y - Ye)
    res["esfera_vs_grs80_desplaz_max_km"] = d.max() / 1000
    res["esfera_vs_grs80_desplaz_mediana_km"] = np.median(d) / 1000

    la, lo = malla_interior(*cargar_contorno())
    k_elip = Proj(CRS.from_epsg(6372)).get_factors(lo, la).meridional_scale
    res["delta_lambert_grs80"] = distorsion(k_elip)
    res["delta_lambert_esfera"] = distorsion(sigma_lambert(la))

    # 3) haversine contra geodésica sobre GRS80
    rng = np.random.default_rng(semilla)
    i = rng.integers(0, len(lat), n_pares)
    j = rng.integers(0, len(lat), n_pares)
    m = i != j
    i, j = i[m], j[m]
    dh = haversine(lat[i], lon[i], lat[j], lon[j])
    _, _, dg = Geod(ellps="GRS80").inv(lon[i], lat[i], lon[j], lat[j])
    r = 100 * (dh - dg) / dg
    res["pares"] = len(i)
    res["haversine_vs_grs80_min_%"] = r.min()
    res["haversine_vs_grs80_mediana_%"] = np.median(r)
    res["haversine_vs_grs80_max_%"] = r.max()
    return res


if __name__ == "__main__":
    res = verificar()
    for k, v in res.items():
        print(f"{k:40s} {v:.6g}")
    out = RAIZ / "resultados" / "tablas" / "verificacion_proj.csv"
    pd.Series(res).to_csv(out, header=["valor"])
    print(f"guardado en {out}")
