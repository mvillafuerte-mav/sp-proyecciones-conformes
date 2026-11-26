"""Herramienta para el cliente: convierte (lat, lon) a (x, y) en metros.

Uso desde Python:

    from proyecciones.conversion import latlon_a_xy
    x, y = latlon_a_xy([19.4326, 25.6866, 20.9674], [-99.1332, -100.3161, -89.5926])

Uso desde la terminal, con un CSV que tenga columnas lat y lon:

    python -m proyecciones.conversion entrada.csv salida.csv
    python -m proyecciones.conversion entrada.csv salida.csv --proyeccion lambert

Proyecciones disponibles:
  "recomendada"        = "optima_envolvente"
  "optima_envolvente"  óptima de Chebyshev sobre la envolvente convexa de México
                       (error de escala ≤ 0.19 %, también para pares que cruzan
                       los golfos)
  "optima"             óptima sobre México continental (≤ 0.092 % dentro del
                       territorio, pero peor para pares que cruzan el mar)
  "lambert"            Lambert de INEGI (EPSG:6372) sobre la esfera
  "mercator"           Mercator esférica

Las óptimas usan los modelos ya ajustados de resultados/modelos/, así que no
hace falta descargar datos ni resolver nada.
"""
import argparse
import warnings
from pathlib import Path

import numpy as np

from .clasicas import lambert, mercator

MODELOS = Path(__file__).resolve().parents[2] / "resultados" / "modelos"
ARCHIVOS = {"optima": "optima_mexico.npz", "optima_envolvente": "optima_envolvente.npz"}
_cache = {}

# un punto a más de ~0.01 en w (≈ 60 km) de la malla está fuera de la región
DISTANCIA_AVISO = 0.01


def _modelo(nombre):
    if nombre not in _cache:
        from .dirichlet import ProyeccionOptima
        _cache[nombre] = ProyeccionOptima.cargar(MODELOS / ARCHIVOS[nombre])
    return _cache[nombre]


def latlon_a_xy(lat, lon, proyeccion="recomendada"):
    """lat, lon en grados (números o arreglos). Devuelve (x, y) en metros.

    Para las óptimas el origen está en (23.5° N, 102° O); para Lambert, en el
    origen de INEGI (falso este 2 500 000 m).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if proyeccion == "recomendada":
        proyeccion = "optima_envolvente"
    if proyeccion == "mercator":
        return mercator(lat, lon)
    if proyeccion == "lambert":
        return lambert(lat, lon)
    if proyeccion not in ARCHIVOS:
        raise ValueError(f"proyección desconocida: {proyeccion!r}")
    m = _modelo(proyeccion)
    lejos = m.distancia_a_malla(lat, lon) > DISTANCIA_AVISO
    if np.any(lejos):
        warnings.warn(f"{int(lejos.sum())} punto(s) quedan fuera de la región "
                      "donde se ajustó la proyección; ahí no hay garantía de escala.")
    x, y = m.proyectar(lat, lon)
    if lat.ndim == 0:
        return float(x[0]), float(y[0])
    return x.reshape(lat.shape), y.reshape(lat.shape)


def convertir_csv(entrada, salida, proyeccion="recomendada", col_lat="lat", col_lon="lon"):
    """Lee un CSV, agrega columnas x, y (metros) y lo guarda."""
    import pandas as pd

    df = pd.read_csv(entrada)
    df["x"], df["y"] = latlon_a_xy(df[col_lat].to_numpy(), df[col_lon].to_numpy(), proyeccion)
    df.to_csv(salida, index=False)
    return df


EJEMPLO = {
    "Ciudad de México": (19.4326, -99.1332),
    "Monterrey": (25.6866, -100.3161),
    "Mérida": (20.9674, -89.5926),
}


def ejemplo():
    """Las tres ciudades del ejemplo y su distancia en el plano contra haversine."""
    from .datos import haversine

    nombres = list(EJEMPLO)
    lat = np.array([EJEMPLO[n][0] for n in nombres])
    lon = np.array([EJEMPLO[n][1] for n in nombres])
    x, y = latlon_a_xy(lat, lon)
    filas = []
    for n, xi, yi in zip(nombres, x, y):
        filas.append(f"{n:18s} x = {xi:12.1f} m   y = {yi:12.1f} m")
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        d = np.hypot(x[i] - x[j], y[i] - y[j])
        dh = haversine(lat[i], lon[i], lat[j], lon[j])
        filas.append(f"{nombres[i]} - {nombres[j]}: plano {d / 1000:.3f} km, "
                     f"esfera {dh / 1000:.3f} km, error {100 * (d / dh - 1):+.3f} %")
    return "\n".join(filas)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Convierte lat/lon a x/y en metros.")
    ap.add_argument("entrada", nargs="?", help="CSV con columnas lat y lon")
    ap.add_argument("salida", nargs="?", help="CSV de salida")
    ap.add_argument("--proyeccion", default="recomendada",
                    choices=["recomendada", "optima_envolvente", "optima", "lambert", "mercator"])
    ap.add_argument("--lat", default="lat")
    ap.add_argument("--lon", default="lon")
    args = ap.parse_args()
    if args.entrada is None:
        print(ejemplo())
    else:
        convertir_csv(args.entrada, args.salida or args.entrada.replace(".csv", "_xy.csv"),
                      args.proyeccion, args.lat, args.lon)
