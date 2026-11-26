"""Contorno de México y cabeceras municipales del INEGI (Paso 5).

Fuente: servicio web del Catálogo Único de Claves Geoestadísticas del INEGI
(https://gaia.inegi.org.mx/wscatgeo/v2/). Los archivos descargados van en
data/raw/ y los procesados en data/processed/ (ninguno se sube a git).

Para descargar y procesar todo:
    python -m proyecciones.datos
(desde la carpeta src/, o con src/ en el PYTHONPATH)
"""
import json
import time
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from .constantes import R_M

RAIZ = Path(__file__).resolve().parents[2]
RAW = RAIZ / "data" / "raw"
PROCESSED = RAIZ / "data" / "processed"

API = "https://gaia.inegi.org.mx/wscatgeo/v2"
ARCHIVO_ESTADOS = RAW / "inegi_mgee.geojson"
ARCHIVO_CABECERAS_RAW = RAW / "inegi_cabeceras.csv"
ARCHIVO_CONTORNO = PROCESSED / "contorno_mexico.csv"
ARCHIVO_PUNTOS = PROCESSED / "cabeceras.csv"

# Regiones por entidad (clave INEGI). Norte, Centro y Sur.
REGIONES = {
    "Norte": ["02", "03", "05", "08", "10", "19", "25", "26", "28", "32"],
    "Centro": ["01", "06", "09", "11", "13", "14", "15", "16", "17", "18",
               "21", "22", "24", "29"],
    "Sur": ["04", "07", "12", "20", "23", "27", "30", "31"],
}
REGION_DE_ENTIDAD = {e: r for r, ents in REGIONES.items() for e in ents}


# ---------------------------------------------------------------- descarga

def _get_json(url, intentos=4):
    for k in range(intentos):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception:
            if k == intentos - 1:
                raise
            time.sleep(2 ** k)


def descargar_estados():
    """GeoJSON de las 32 entidades (Marco Geoestadístico)."""
    RAW.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(f"{API}/geo/mgee/", ARCHIVO_ESTADOS)


def _entero(x, defecto=0):
    """La API a veces trae '-' o nada en lugar de un número."""
    try:
        return int(x)
    except (TypeError, ValueError):
        return defecto


def _fila_cabecera(ent, m):
    """Busca la cabecera del municipio m.

    Las alcaldías de la CDMX no tienen cabecera (cve_cab = '----'); en ese
    caso se usa la localidad 0001, que es la sede de la alcaldía, y si no
    existe, la más poblada.
    """
    url = f"{API}/localidades/{ent}/{m['cve_mun']}"
    locs = [l for l in _get_json(url).get("datos", []) if l.get("latitud")]
    por_clave = {l["cve_loc"]: l for l in locs}
    cab = por_clave.get(m.get("cve_cab")) or por_clave.get("0001")
    if cab is None and locs:
        cab = max(locs, key=lambda l: _entero(l.get("pob_total")))
    if cab is None:
        print(f"sin cabecera: {m['cvegeo']} {m['nomgeo']}")
        return None
    return {
        "cvegeo": m["cvegeo"],
        "cve_ent": ent,
        "municipio": m["nomgeo"],
        "cabecera": cab["nomgeo"],
        "lat": float(cab["latitud"]),
        "lon": float(cab["longitud"]),
        "pob_municipio": _entero(m.get("pob_total"), None),
    }


def descargar_cabeceras(pausa=0.05):
    """Coordenadas de la cabecera de cada municipio.

    Para cada municipio se pide su lista de localidades y se toma la que
    tiene la clave de la cabecera (cve_cab).
    """
    filas = []
    for ent in range(1, 33):
        ent = f"{ent:02d}"
        for m in _get_json(f"{API}/mgem/{ent}").get("datos", []):
            fila = _fila_cabecera(ent, m)
            if fila:
                filas.append(fila)
            time.sleep(pausa)
        print(f"entidad {ent}: {len(filas)} cabeceras acumuladas")
    RAW.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(filas)
    df.to_csv(ARCHIVO_CABECERAS_RAW, index=False)
    (RAW / "LEEME_descarga.txt").write_text(
        f"Descargado el {date.today().isoformat()} de {API}\n"
        f"Cabeceras: {len(df)}\n", encoding="utf-8")
    return df


# ---------------------------------------------------------------- contorno

def _poligono_continental():
    import geopandas as gpd

    estados = gpd.read_file(ARCHIVO_ESTADOS)
    pais = estados.geometry.union_all()
    if pais.geom_type == "MultiPolygon":
        # nos quedamos con el polígono más grande (sin islas)
        pais = max(pais.geoms, key=lambda g: g.area)
    return pais


def procesar_contorno(tolerancia_grados=0.01):
    """Contorno de México continental, simplificado.

    tolerancia_grados: tolerancia de Douglas-Peucker en grados
    (0.01° ≈ 1 km). Se quitan los huecos interiores.
    """
    from shapely.geometry import Polygon

    pol = Polygon(_poligono_continental().exterior)
    if tolerancia_grados:
        pol = pol.simplify(tolerancia_grados, preserve_topology=True)
    lon, lat = np.asarray(pol.exterior.coords).T
    PROCESSED.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"lat": lat, "lon": lon}).to_csv(ARCHIVO_CONTORNO, index=False)
    return lat, lon


def cargar_contorno():
    """Contorno de México continental: (lat, lon) en grados, curva cerrada."""
    df = pd.read_csv(ARCHIVO_CONTORNO)
    return df["lat"].to_numpy(), df["lon"].to_numpy()


# ---------------------------------------------------------------- puntos

def procesar_puntos():
    """Agrega región y si el punto cae dentro del contorno continental."""
    from shapely import contains_xy
    from shapely.geometry import Polygon

    df = pd.read_csv(ARCHIVO_CABECERAS_RAW, dtype={"cvegeo": str, "cve_ent": str})
    # polígono a resolución completa con ~1 km de margen, para no perder
    # cabeceras que están en la costa o pegadas a la frontera
    pol = Polygon(_poligono_continental().exterior).buffer(0.01)
    df["region"] = df["cve_ent"].map(REGION_DE_ENTIDAD)
    df["continental"] = contains_xy(pol, df["lon"], df["lat"])
    df.to_csv(ARCHIVO_PUNTOS, index=False)
    return df


def cargar_puntos(solo_continentales=True):
    """Cabeceras municipales: DataFrame con lat, lon (grados), región, etc."""
    df = pd.read_csv(ARCHIVO_PUNTOS, dtype={"cvegeo": str, "cve_ent": str})
    if solo_continentales:
        df = df[df["continental"]].reset_index(drop=True)
    return df


# ---------------------------------------------------------------- distancias

def haversine(lat1, lon1, lat2, lon2, R=R_M):
    """Distancia sobre la esfera en metros. Ángulos en grados (acepta arreglos)."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


if __name__ == "__main__":
    if not ARCHIVO_ESTADOS.exists():
        descargar_estados()
    if not ARCHIVO_CABECERAS_RAW.exists():
        descargar_cabeceras()
    lat, lon = procesar_contorno()
    df = procesar_puntos()
    print(f"contorno: {len(lat)} vértices, lat {lat.min():.3f}..{lat.max():.3f}")
    print(f"cabeceras: {len(df)} ({df['continental'].sum()} continentales)")
    print(df.groupby("region").size())
