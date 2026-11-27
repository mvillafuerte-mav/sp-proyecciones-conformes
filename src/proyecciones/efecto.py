"""Efecto de la proyección en los análisis del cliente (Paso 5).

    python -m proyecciones.efecto

Compara Mercator, Lambert de INEGI, la óptima y la óptima sobre la envolvente
convexa contra la referencia sobre la esfera (haversine), con las cabeceras:
  1. error relativo en distancias, por región;
  2. k-means: asignaciones que cambian respecto a k-means sobre la esfera;
  3. k vecinos más cercanos: puntos cuyo conjunto de vecinos cambia;
  4. áreas: factor σ² en cada región.
Las tablas quedan en resultados/tablas/efecto_*.csv y los errores de cada par
en data/processed/errores_pares.csv (para la Figura 3).
"""
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import BallTree, NearestNeighbors

from .clasicas import lambert, mercator
from .datos import PROCESSED, RAIZ, cargar_puntos, haversine
from .dirichlet import proyectar_optima, proyectar_optima_envolvente
from .distorsion import sigma_lambert, sigma_mercator

PROYECCIONES = {
    "Mercator": mercator,
    "Lambert INEGI": lambert,
    "Óptima": proyectar_optima,
    "Óptima envolvente": proyectar_optima_envolvente,
}
REGIONES = ["Norte", "Centro", "Sur"]


def coordenadas_planas(df):
    """{nombre: arreglo (n, 2) en metros} para cada proyección."""
    lat, lon = df["lat"].to_numpy(), df["lon"].to_numpy()
    return {nombre: np.column_stack(f(lat, lon)) for nombre, f in PROYECCIONES.items()}


# ------------------------------------------------------------ distancias

def pares_aleatorios(df, n_por_grupo=20000, semilla=0):
    """Pares (i, j) al azar: dentro de cada región y en todo el país."""
    rng = np.random.default_rng(semilla)
    filas = []
    grupos = [(r, np.where(df["region"] == r)[0]) for r in REGIONES]
    grupos.append(("Nacional", np.arange(len(df))))
    for nombre, idx in grupos:
        i = rng.choice(idx, n_por_grupo)
        j = rng.choice(idx, n_por_grupo)
        ok = i != j
        filas.append(pd.DataFrame({"grupo": nombre, "i": i[ok], "j": j[ok]}))
    return pd.concat(filas, ignore_index=True)


def geodesica_dentro(df, pares, n=60):
    """True si el arco de círculo máximo de cada par queda dentro de México."""
    from shapely import contains, linestrings
    from shapely.geometry import Polygon
    from .datos import cargar_contorno

    lat_c, lon_c = cargar_contorno()
    pol = Polygon(np.column_stack([lon_c, lat_c]))
    p, l = np.radians(df["lat"].to_numpy()), np.radians(df["lon"].to_numpy())
    U = np.column_stack([np.cos(p) * np.cos(l), np.cos(p) * np.sin(l), np.sin(p)])
    a, b = U[pares["i"].to_numpy()], U[pares["j"].to_numpy()]
    t = np.linspace(0, 1, n)[None, :, None]
    v = a[:, None, :] * (1 - t) + b[:, None, :] * t
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    la = np.degrees(np.arcsin(v[..., 2]))
    lo = np.degrees(np.arctan2(v[..., 1], v[..., 0]))
    return contains(pol, linestrings(np.stack([lo, la], axis=2)))


def errores_distancia(df, xy, pares):
    """Error relativo de cada par en cada proyección, crudo y con la mejor
    escala constante (la mediana nacional del cociente)."""
    lat, lon = df["lat"].to_numpy(), df["lon"].to_numpy()
    i, j = pares["i"].to_numpy(), pares["j"].to_numpy()
    d_ref = haversine(lat[i], lon[i], lat[j], lon[j])
    out = pares.copy()
    out["d_km"] = d_ref / 1000
    out["dentro"] = geodesica_dentro(df, pares)
    nacional = (pares["grupo"] == "Nacional").to_numpy()
    for nombre, P in xy.items():
        razon = np.hypot(*(P[i] - P[j]).T) / d_ref
        k = np.median(razon[nacional])
        out[f"{nombre}|crudo"] = 100 * (razon - 1)
        out[f"{nombre}|ajustado"] = 100 * (razon / k - 1)
    return out


def resumen_distancias(err):
    filas = []
    cols = [c for c in err.columns if "|" in c]
    for grupo in REGIONES + ["Nacional"]:
        e = err[err["grupo"] == grupo]
        for c in cols:
            nombre, tipo = c.split("|")
            a = np.abs(e[c].to_numpy())
            a_dentro = np.abs(e.loc[e["dentro"], c].to_numpy())
            filas.append({"grupo": grupo, "proyeccion": nombre, "escala": tipo,
                          "mediana_%": np.median(a), "p95_%": np.percentile(a, 95),
                          "max_%": a.max(), "max_dentro_%": a_dentro.max(),
                          "sesgo_mediano_%": np.median(e[c])})
    return pd.DataFrame(filas)


# ------------------------------------------------------------ k-means

def _unitarios(df):
    p, l = np.radians(df["lat"].to_numpy()), np.radians(df["lon"].to_numpy())
    return np.column_stack([np.cos(p) * np.cos(l), np.cos(p) * np.sin(l), np.sin(p)])


def kmeans_lloyd(X, centros, esferico=False, iters=300):
    """Lloyd a partir de centros dados. Si esferico, X son vectores unitarios,
    se asigna por máximo producto punto (= mínima distancia geodésica) y los
    centros se renormalizan."""
    C = centros.copy()
    etiquetas = None
    for _ in range(iters):
        if esferico:
            nuevas = np.argmax(X @ C.T, axis=1)
        else:
            d2 = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
            nuevas = np.argmin(d2, axis=1)
        if etiquetas is not None and np.array_equal(nuevas, etiquetas):
            break
        etiquetas = nuevas
        for k in range(len(C)):
            miembros = X[etiquetas == k]
            if len(miembros):
                C[k] = miembros.mean(axis=0)
                if esferico:
                    C[k] /= np.linalg.norm(C[k])
    return etiquetas


def _kmeanspp(X, k, rng):
    """Centros iniciales k-means++ con distancia geodésica (X unitarios)."""
    idx = [rng.integers(len(X))]
    d = np.full(len(X), np.inf)
    for _ in range(k - 1):
        d = np.minimum(d, np.arccos(np.clip(X @ X[idx[-1]], -1, 1)) ** 2)
        idx.append(rng.choice(len(X), p=d / d.sum()))
    return np.array(idx)


def cambios_kmeans(df, xy, ks=(10, 20, 50), semillas=range(10)):
    """% de puntos cuyo grupo cambia respecto a k-means sobre la esfera.

    Para aislar el efecto de la proyección, cada corrida plana arranca de los
    mismos puntos iniciales que la esférica; las etiquetas se emparejan con
    el método húngaro antes de comparar.
    """
    U = _unitarios(df)
    filas = []
    for k in ks:
        for s in semillas:
            ini = _kmeanspp(U, k, np.random.default_rng(s))
            ref = kmeans_lloyd(U, U[ini], esferico=True)
            for nombre, P in xy.items():
                et = kmeans_lloyd(P, P[ini])
                M = np.zeros((k, k), dtype=int)
                np.add.at(M, (ref, et), 1)
                r, c = linear_sum_assignment(-M)
                iguales = M[r, c].sum()
                filas.append({"k": k, "semilla": s, "proyeccion": nombre,
                              "cambian_%": 100 * (1 - iguales / len(U))})
    return pd.DataFrame(filas)


# ------------------------------------------------------------ vecinos

def cambios_vecinos(df, xy, k=10):
    """% de puntos cuyo conjunto de k vecinos cambia respecto a haversine."""
    rad = np.radians(df[["lat", "lon"]].to_numpy())
    _, ref = BallTree(rad, metric="haversine").query(rad, k=k + 1)
    ref = [set(r[1:]) for r in ref]
    filas = []
    for nombre, P in xy.items():
        _, vec = NearestNeighbors(n_neighbors=k + 1).fit(P).kneighbors(P)
        cambia = np.array([set(v[1:]) != r for v, r in zip(vec, ref)])
        filas.append({"proyeccion": nombre, "k": k, "cambian_%": 100 * cambia.mean()})
        for reg in REGIONES:
            m = (df["region"] == reg).to_numpy()
            filas.append({"proyeccion": nombre, "k": k, "region": reg,
                          "cambian_%": 100 * cambia[m].mean()})
    return pd.DataFrame(filas)


# ------------------------------------------------------------ áreas

def factores_area(df):
    """σ² (cuánto se agranda un área) en las cabeceras, por región."""
    from .dirichlet import obtener_modelo
    lat, lon = df["lat"].to_numpy(), df["lon"].to_numpy()
    s2 = {"Mercator": sigma_mercator(lat) ** 2,
          "Lambert INEGI": sigma_lambert(lat) ** 2,
          "Óptima": np.exp(2 * obtener_modelo("optima").ln_sigma(lat, lon)),
          "Óptima envolvente": np.exp(2 * obtener_modelo("optima_envolvente").ln_sigma(lat, lon))}
    filas = []
    for nombre, v in s2.items():
        for reg in REGIONES + ["Nacional"]:
            m = np.ones(len(df), bool) if reg == "Nacional" else (df["region"] == reg).to_numpy()
            filas.append({"proyeccion": nombre, "grupo": reg,
                          "sigma2_min": v[m].min(), "sigma2_max": v[m].max(),
                          "sigma2_mediana": np.median(v[m])})
    return pd.DataFrame(filas)


# ------------------------------------------------------------ elipsoide

def contra_elipsoide(df, xy, pares, error_posicion_m=50.0):
    """Error contra la geodésica real sobre GRS80 (no contra la esfera), para
    nuestras proyecciones y para la Lambert oficial EPSG:6372 (elipsoidal).

    También da el umbral de distancia a partir del cual el error de la
    proyección supera un error de posición dado: d = error_posicion / error.
    """
    from pyproj import Geod, Transformer

    lat, lon = df["lat"].to_numpy(), df["lon"].to_numpy()
    i, j = pares["i"].to_numpy(), pares["j"].to_numpy()
    _, _, dg = Geod(ellps="GRS80").inv(lon[i], lat[i], lon[j], lat[j])
    todas = dict(xy)
    X, Y = Transformer.from_crs(4326, 6372, always_xy=True).transform(lon, lat)
    todas["EPSG:6372 (elipsoide)"] = np.column_stack([X, Y])
    filas = []
    for nombre, P in todas.items():
        r = 100 * np.abs(np.hypot(*(P[i] - P[j]).T) / dg - 1)
        med, mx = np.median(r), r.max()
        filas.append({"proyeccion": nombre, "mediana_%": med,
                      "p95_%": np.percentile(r, 95), "max_%": mx,
                      "umbral_tipico_km": error_posicion_m / (med / 100) / 1000,
                      "umbral_peor_km": error_posicion_m / (mx / 100) / 1000})
    return pd.DataFrame(filas)


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    tablas = RAIZ / "resultados" / "tablas"
    df = cargar_puntos()
    xy = coordenadas_planas(df)

    pares = pares_aleatorios(df)
    err = errores_distancia(df, xy, pares)
    err.to_csv(PROCESSED / "errores_pares.csv", index=False)
    rd = resumen_distancias(err)
    rd.to_csv(tablas / "efecto_distancias.csv", index=False)
    pd.Series({"pares": len(err), "pares_dentro_%": 100 * err["dentro"].mean()}).to_csv(
        tablas / "efecto_pares.csv", header=["valor"])
    print(f"pares: {len(err)}, con geodésica dentro de México: {100 * err['dentro'].mean():.1f} %")
    print("1) error en distancias (|%|)")
    print(rd.round(4).to_string(index=False))

    km = cambios_kmeans(df, xy)
    km.to_csv(tablas / "efecto_kmeans_corridas.csv", index=False)
    kmr = km.groupby(["k", "proyeccion"])["cambian_%"].agg(
        mediana="median", media="mean", maximo="max",
        corridas_mas_5pct=lambda x: int((x > 5).sum())).reset_index()
    kmr.to_csv(tablas / "efecto_kmeans.csv", index=False)
    print("\n2) k-means: % de asignaciones que cambian (10 semillas)")
    print(kmr.round(3).to_string(index=False))

    kn = cambios_vecinos(df, xy)
    kn.to_csv(tablas / "efecto_vecinos.csv", index=False)
    print("\n3) 10 vecinos: % de puntos cuyo conjunto cambia")
    print(kn.round(3).to_string(index=False))

    ar = factores_area(df)
    ar.to_csv(tablas / "efecto_areas.csv", index=False)
    print("\n4) factor de área σ²")
    print(ar.round(5).to_string(index=False))

    el = contra_elipsoide(df, xy, pares)
    el.to_csv(tablas / "efecto_elipsoide.csv", index=False)
    print("\n5) contra la geodésica sobre GRS80, y umbral para 50 m de error de posición")
    print(el.round(3).to_string(index=False))
