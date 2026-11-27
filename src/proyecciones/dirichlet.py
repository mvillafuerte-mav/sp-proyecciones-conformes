"""Proyección óptima de Chebyshev (Paso 4).

Pasos:
  1. Llevar el contorno de México al plano w = λ + i q(φ).
  2. Resolver Δh = 0 en Ω con h = ln cos φ en ∂Ω  (con R = 1; la constante c
     se agrega al final para centrar la escala).
  3. Conjugada armónica h~, integrando las ecuaciones de Cauchy-Riemann.
  4. f' = R exp(h + c + i h~) e integrar para obtener f.

Método: diferencias finitas de 5 puntos en una malla cuadrada del plano w.
Cerca de la frontera se usa la fórmula de Shortley-Weller (las distancias
reales al contorno en lugar de h), que mantiene el error en O(h²).

Antes se probó ajustar h como parte real de un polinomio complejo
(ver `ajuste_polinomial`), pero no converge: Baja California y el Golfo
forman un canal angosto fuera de Ω y el residuo se estanca en ~1e-3.

    python -m proyecciones.dirichlet
"""
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import breadth_first_order
from scipy.spatial import cKDTree

from .constantes import R_M, latitud_isometrica

LAT_CENTRO = 23.5      # punto de referencia: ahí f = 0 y el norte apunta hacia arriba
LON_CENTRO = -102.0

# vecinos: este, oeste, norte, sur  (di en λ, dj en q)
DIRS = [(1, 0), (-1, 0), (0, 1), (0, -1)]


# ------------------------------------------------------------ coordenadas

def w_de(lat, lon, lon0=LON_CENTRO):
    """w = λ + i q(φ), con λ medida desde lon0. Grados de entrada, radianes de salida."""
    return np.radians(np.asarray(lon) - lon0) + 1j * latitud_isometrica(np.radians(lat))


def latlon_de(w, lon0=LON_CENTRO):
    """Inversa de w_de."""
    lat = np.degrees(2 * np.arctan(np.exp(np.imag(w))) - np.pi / 2)
    return lat, np.degrees(np.real(w)) + lon0


def ln_cos_lat(w):
    """ln cos φ escrito con q = Im w:  cos φ = sech q."""
    return -np.log(np.cosh(np.imag(w)))


# ------------------------------------------------------------ malla

class Malla:
    """Malla cuadrada de paso h en el plano w, con los nodos dentro de Ω y,
    para cada nodo y cada dirección, la distancia al vecino o a la frontera."""

    def __init__(self, w_contorno, h):
        from shapely import contains_xy, get_coordinates, intersection, linestrings
        from shapely.geometry import LineString, Polygon

        self.h = h
        pts = np.column_stack([w_contorno.real, w_contorno.imag])
        self.poligono = Polygon(pts)
        borde = LineString(np.vstack([pts, pts[:1]]))

        x0, y0 = pts.min(axis=0) - h
        x1, y1 = pts.max(axis=0) + h
        nx = int(np.ceil((x1 - x0) / h)) + 1
        ny = int(np.ceil((y1 - y0) / h)) + 1
        I, J = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
        X, Y = x0 + I * h, y0 + J * h
        dentro = contains_xy(self.poligono, X, Y)

        self.ij = np.column_stack([I[dentro], J[dentro]])
        self.w = X[dentro] + 1j * Y[dentro]
        n = len(self.w)
        num = -np.ones((nx, ny), dtype=int)
        num[dentro] = np.arange(n)

        # para cada dirección: índice del vecino (-1 si es frontera) y distancia
        self.vecino = -np.ones((n, 4), dtype=int)
        self.dist = np.full((n, 4), h)
        self.w_borde = np.zeros((n, 4), dtype=complex)
        for k, (di, dj) in enumerate(DIRS):
            ii, jj = self.ij[:, 0] + di, self.ij[:, 1] + dj
            nb = num[ii, jj]
            self.vecino[:, k] = nb
            afuera = np.where(nb < 0)[0]
            if len(afuera) == 0:
                continue
            p = self.w[afuera]
            q = p + h * (di + 1j * dj)
            segs = linestrings(np.stack([np.column_stack([p.real, p.imag]),
                                         np.column_stack([q.real, q.imag])], axis=1))
            cortes = intersection(segs, borde)
            xy, idx = get_coordinates(cortes, return_index=True)
            d = np.abs((xy[:, 0] + 1j * xy[:, 1]) - p[idx])
            dmin = np.full(len(afuera), h)
            np.minimum.at(dmin, idx, d)
            dmin = np.clip(dmin, 1e-3 * h, h)
            self.dist[afuera, k] = dmin
            self.w_borde[afuera, k] = p + dmin * (di + 1j * dj)

    def resolver(self, datos_frontera):
        """Δh = 0 en los nodos, h = datos_frontera(w) en ∂Ω (Shortley-Weller)."""
        n = len(self.w)
        filas, cols, vals = [], [], []
        rhs = np.zeros(n)
        diag = np.zeros(n)
        for eje in (0, 1):                         # 0: dirección λ, 1: dirección q
            kp, km = 2 * eje, 2 * eje + 1          # +, -
            dp, dm = self.dist[:, kp], self.dist[:, km]
            for k, d in ((kp, dp), (km, dm)):
                a = 2.0 / (d * (dp + dm))
                diag -= a
                nb = self.vecino[:, k]
                es_nodo = nb >= 0
                filas.append(np.where(es_nodo)[0])
                cols.append(nb[es_nodo])
                vals.append(a[es_nodo])
                fr = ~es_nodo
                rhs[fr] -= a[fr] * datos_frontera(self.w_borde[fr, k])
        filas.append(np.arange(n)); cols.append(np.arange(n)); vals.append(diag)
        A = sp.csr_matrix((np.concatenate(vals),
                           (np.concatenate(filas), np.concatenate(cols))), shape=(n, n))
        u = spla.spsolve(A.tocsc(), rhs)
        self.u = u
        self._datos = datos_frontera
        return u

    def _valores_vecinos(self, k):
        nb = self.vecino[:, k]
        v = np.empty(len(self.w))
        v[nb >= 0] = self.u[nb[nb >= 0]]
        fr = nb < 0
        v[fr] = self._datos(self.w_borde[fr, k])
        return v

    def gradiente(self):
        """(h_λ, h_q) en cada nodo, diferencia centrada no uniforme, O(h²)."""
        g = []
        for eje in (0, 1):
            kp, km = 2 * eje, 2 * eje + 1
            dp, dm = self.dist[:, kp], self.dist[:, km]
            up, um = self._valores_vecinos(kp), self._valores_vecinos(km)
            g.append((dm**2 * (up - self.u) + dp**2 * (self.u - um)) / (dp * dm * (dp + dm)))
        return g[0], g[1]


# ------------------------------------------------------------ proyección

class ProyeccionOptima:
    """Proyección conforme de Chebyshev sobre la región del contorno dado.

    Atributos útiles: delta (distorsión δ en la malla), c (constante que centra
    la escala), y los valores en los nodos: w, h, hconj, fp (= f') y f.
    """

    def __init__(self, lat_contorno, lon_contorno, h=0.001, R=R_M):
        t0 = time.perf_counter()
        self.R = R
        self.paso = h
        wc = w_de(lat_contorno, lon_contorno)
        self.contorno = wc
        self.malla = m = Malla(wc, h)
        self.w = m.w
        self.h = m.resolver(ln_cos_lat)
        hl, hq = m.gradiente()
        self.gprima = hl - 1j * hq                 # g' = h_λ + i h~_λ = h_λ − i h_q

        # ln σ = h − ln cos φ (sin la constante c); en ∂Ω vale 0
        ls = self.h - ln_cos_lat(self.w)
        mx, mn = max(ls.max(), 0.0), min(ls.min(), 0.0)
        self.delta = mx - mn
        self.c = -(mx + mn) / 2                    # centra: ln σ en [−δ/2, δ/2]

        self._integrar(hl, hq)
        self._nodos = np.where(self.conectados)[0]
        wn = self.w[self._nodos]
        self._arbol = cKDTree(np.column_stack([wn.real, wn.imag]))
        self.segundos = time.perf_counter() - t0

    def _integrar(self, hl, hq):
        """h~ y f sobre los nodos.

        Cada arista de la malla da un incremento (Cauchy-Riemann para h~,
        ∫ f' dw para f). Como los incrementos traen error numérico, no se
        suman por un solo camino: se busca la función nodal que mejor los
        reproduce en todas las aristas a la vez (mínimos cuadrados), y así
        el error queda repartido de forma suave.
        """
        m, n = self.malla, len(self.w)
        ps, cs = [], []
        for k in (0, 2):                            # aristas hacia el este y el norte
            ok = m.vecino[:, k] >= 0
            ps.append(np.where(ok)[0]); cs.append(m.vecino[ok, k])
        p_, c_ = np.concatenate(ps), np.concatenate(cs)
        ne = len(p_)
        D = sp.csr_matrix((np.r_[-np.ones(ne), np.ones(ne)],
                           (np.r_[np.arange(ne), np.arange(ne)], np.r_[p_, c_])),
                          shape=(ne, n))
        # componente conexa del centro
        w0 = w_de(LAT_CENTRO, LON_CENTRO)
        raiz = int(np.argmin(np.abs(self.w - w0)))
        G = (D.T @ D).tocsr()
        orden = breadth_first_order(G, raiz, directed=False, return_predecessors=False)
        self.conectados = np.zeros(n, dtype=bool)
        self.conectados[orden] = True
        keep = self.conectados[p_] & self.conectados[c_]
        D = D[keep][:, orden]
        p_, c_ = p_[keep], c_[keep]
        dw = self.w[c_] - self.w[p_]
        L = (D.T @ D).tolil()
        L[0, :] = 0
        L[0, 0] = 1                                 # fija el valor en la raíz
        L = L.tocsc()
        resolver = spla.factorized(L)

        def integrar(incremento, valor_raiz):
            b = D.T @ incremento
            b[0] = valor_raiz
            x = np.full(n, np.nan, dtype=incremento.dtype)
            if np.iscomplexobj(incremento):
                x[orden] = resolver(b.real) + 1j * resolver(b.imag)
            else:
                x[orden] = resolver(b)
            return x

        # dh~ = −h_q dλ + h_λ dq   (trapecio)
        inc_h = (-(hq[p_] + hq[c_]) * dw.real + (hl[p_] + hl[c_]) * dw.imag) / 2
        self.hconj = integrar(inc_h, 0.0)
        self.fp = self.R * np.exp(self.h + self.c + 1j * self.hconj)
        fpp = self.fp * self.gprima
        # trapecio con corrección de Hermite: error local O(h^5)
        inc_f = (self.fp[p_] + self.fp[c_]) * dw / 2 - (fpp[c_] - fpp[p_]) * dw**2 / 12
        f = integrar(inc_f.astype(complex), 0.0)
        # f(w0) = 0 en el punto de referencia (no en el nodo más cercano)
        self.f = f - self.fp[raiz] * (w0 - self.w[raiz])

    # ---- guardar / cargar (para usar la proyección sin recalcularla)

    def guardar(self, ruta):
        """Guarda solo lo necesario para proyectar (nodos conectados)."""
        i = self._nodos
        np.savez_compressed(ruta, w=self.w[i], f=self.f[i], fp=self.fp[i],
                            gprima=self.gprima[i], h=self.h[i], c=self.c,
                            delta=self.delta, R=self.R, paso=self.paso,
                            contorno=self.contorno)

    @classmethod
    def cargar(cls, ruta):
        d = np.load(ruta)
        p = cls.__new__(cls)
        p.w, p.f, p.fp, p.gprima, p.h = d["w"], d["f"], d["fp"], d["gprima"], d["h"]
        p.c, p.delta, p.R, p.paso = float(d["c"]), float(d["delta"]), float(d["R"]), float(d["paso"])
        p.contorno = d["contorno"] if "contorno" in d else None
        p.conectados = np.ones(len(p.w), dtype=bool)
        p._nodos = np.arange(len(p.w))
        p._arbol = cKDTree(np.column_stack([p.w.real, p.w.imag]))
        return p

    def contiene(self, lat, lon):
        """Pertenencia al polígono de ajuste, incluyendo su frontera."""
        from shapely import intersects_xy
        from shapely.geometry import Polygon
        if self.contorno is None:
            raise ValueError("El modelo antiguo no contiene su dominio; regenérelo.")
        w = np.atleast_1d(w_de(lat, lon))
        pol = Polygon(np.column_stack([self.contorno.real, self.contorno.imag]))
        return intersects_xy(pol, w.real, w.imag)

    def distancia_a_malla(self, lat, lon):
        """Distancia (en unidades de w) al nodo más cercano; sirve para avisar
        si un punto cae fuera de la región donde se ajustó la proyección."""
        w = np.atleast_1d(w_de(lat, lon))
        d, _ = self._arbol.query(np.column_stack([w.real, w.imag]))
        return d

    def _nodo_cercano(self, w):
        w = np.atleast_1d(w)
        _, k = self._arbol.query(np.column_stack([w.real, w.imag]))
        return self._nodos[k]

    def ln_sigma(self, lat, lon):
        """ln σ en puntos arbitrarios (interpolado del nodo más cercano)."""
        w = w_de(lat, lon)
        i = self._nodo_cercano(w)
        hw = self.h[i] + (self.gprima[i] * (w - self.w[i])).real
        return hw + self.c - ln_cos_lat(w)

    def proyectar(self, lat, lon):
        """(lat, lon) en grados -> (x, y) en metros. Origen en (23.5° N, 102° O),
        norte hacia arriba en ese punto.

        Para puntos fuera de Ω (islas, o costa recortada por la
        simplificación) el valor se extrapola desde el nodo más cercano.
        La reconstrucción es aproximada y puede tener saltos milimétricos entre
        celdas. La cota del problema continuo no certifica esas interfaces.
        """
        w = np.atleast_1d(w_de(lat, lon))
        i = self._nodo_cercano(w)
        dw = w - self.w[i]
        gp = self.gprima[i]
        # Aproximación local con g' constante; no es el Taylor cúbico exacto.
        z = self.f[i] + self.fp[i] * (dw + gp * dw**2 / 2 + gp**2 * dw**3 / 6)
        return z.real, z.imag


# ------------------------------------------------------------ verificación

def validar(lat_c, lon_c, pasos=(0.004, 0.002, 0.001, 0.0005)):
    """Error del solver con una función armónica conocida, h = ln|w − w*|,
    con w* fuera de Ω (en el Pacífico frente a Baja California)."""
    import pandas as pd

    wc = w_de(lat_c, lon_c)
    w_estrella = w_de(30.0, -118.0)
    exacta = lambda w: np.log(np.abs(w - w_estrella))
    filas = []
    for h in pasos:
        t = time.perf_counter()
        m = Malla(wc, h)
        u = m.resolver(exacta)
        filas.append({"h": h, "nodos": len(u),
                      "error_max": np.max(np.abs(u - exacta(m.w))),
                      "segundos": time.perf_counter() - t})
    return pd.DataFrame(filas)


def verificar_escala(proy, lat, lon, paso_m=50.0):
    """Mide σ de verdad: proyecta dos puntos a ~paso_m al norte y al este y
    divide la distancia en el mapa entre la distancia sobre la esfera.
    Devuelve (error relativo máximo contra exp(ln σ), diferencia máxima
    entre σ norte y σ este, desviación máxima del ángulo recto en grados)."""
    from .datos import haversine

    dlat = np.degrees(paso_m / proy.R)
    dlon = dlat / np.cos(np.radians(lat))
    x, y = proy.proyectar(lat, lon)
    xn, yn = proy.proyectar(lat + dlat, lon)
    xe, ye = proy.proyectar(lat, lon + dlon)
    sn = np.hypot(xn - x, yn - y) / haversine(lat, lon, lat + dlat, lon)
    se = np.hypot(xe - x, ye - y) / haversine(lat, lon, lat, lon + dlon)
    s = np.exp(proy.ln_sigma(lat, lon))
    ang = np.degrees(np.angle((xe - x + 1j * (ye - y)) / (xn - x + 1j * (yn - y))))
    return (np.max(np.abs(sn / s - 1)), np.max(np.abs(sn / se - 1)),
            np.max(np.abs(np.abs(ang) - 90)))


def muestrear_frontera(w_vertices, m):
    """m puntos repartidos por longitud de arco sobre el polígono."""
    w = np.append(w_vertices, w_vertices[0])
    s = np.concatenate([[0], np.cumsum(np.abs(np.diff(w)))])
    t = np.linspace(0, s[-1], m, endpoint=False)
    return np.interp(t, s, w.real) + 1j * np.interp(t, s, w.imag)


def ajuste_polinomial(w_frontera, valores, grado):
    """Residuo máximo al ajustar h = Re p(w) con p polinomio de grado dado
    (mínimos cuadrados, base de Vandermonde con Arnoldi). Solo se usa para
    documentar que este enfoque no converge para México."""
    c = w_frontera.mean()
    Z = (w_frontera - c) / np.max(np.abs(w_frontera - c))
    M = len(Z)
    Q = np.zeros((M, grado + 1), dtype=complex)
    Q[:, 0] = 1
    for k in range(1, grado + 1):
        v = Z * Q[:, k - 1]
        for j in range(k):
            v = v - (Q[:, j].conj() @ v / M) * Q[:, j]
        Q[:, k] = v / (np.linalg.norm(v) / np.sqrt(M))
    A = np.hstack([Q.real, -Q.imag])
    x, *_ = np.linalg.lstsq(A, valores, rcond=None)
    return np.max(np.abs(A @ x - valores))


def contorno_envolvente(lat_c, lon_c, m=2000):
    """Envolvente convexa del contorno (incluye los golfos), densificada."""
    from shapely.geometry import Polygon
    lo, la = np.asarray(Polygon(np.column_stack([lon_c, lat_c])).convex_hull.exterior.coords).T
    z = muestrear_frontera(lo[:-1] + 1j * la[:-1], m)
    return z.imag, z.real


_MODELOS = {}


def obtener_modelo(nombre="optima"):
    """Una única carga de los modelos persistentes para conversión y análisis."""
    from pathlib import Path
    archivos = {"optima": "optima_mexico.npz", "optima_envolvente": "optima_envolvente.npz"}
    if nombre not in archivos:
        raise ValueError(f"Modelo desconocido: {nombre}")
    if nombre not in _MODELOS:
        ruta = Path(__file__).resolve().parents[2] / "resultados/modelos" / archivos[nombre]
        _MODELOS[nombre] = ProyeccionOptima.cargar(ruta)
    return _MODELOS[nombre]


def proyectar_optima(lat, lon):
    """(lat, lon) -> (x, y) con la proyección óptima para México continental."""
    return obtener_modelo("optima").proyectar(lat, lon)


def proyectar_optima_envolvente(lat, lon):
    """Óptima sobre la envolvente convexa de México. La cota es de escala
    local sobre la esfera; las distancias entre pares se verifican aparte."""
    return obtener_modelo("optima_envolvente").proyectar(lat, lon)


if __name__ == "__main__":
    import pandas as pd

    from .datos import RAIZ, cargar_contorno, cargar_puntos
    from .distorsion import error_max_escala

    tablas = RAIZ / "resultados" / "tablas"
    lat_c, lon_c = cargar_contorno()

    print("1) validación con h = ln|w − w*| conocida")
    val = validar(lat_c, lon_c)
    print(val.to_string(index=False))
    val.to_csv(tablas / "dirichlet_validacion.csv", index=False)

    print("\n2) México: δ según el paso de la malla")
    filas = []
    for h in (0.002, 0.001, 0.0005):
        p = ProyeccionOptima(lat_c, lon_c, h=h)
        filas.append({"h": h, "nodos": len(p.w), "delta": p.delta,
                      "segundos": p.segundos})
    conv = pd.DataFrame(filas)
    print(conv.to_string(index=False))
    conv.to_csv(tablas / "optima_convergencia.csv", index=False)

    print("\n3) verificación de punta a punta en las cabeceras (h = 0.001)")
    p = ProyeccionOptima(lat_c, lon_c, h=0.001)
    df = cargar_puntos()
    e_sigma, e_ne, e_ang = verificar_escala(p, df["lat"].to_numpy(),
                                            df["lon"].to_numpy(), paso_m=1000.0)
    ls_cab = p.ln_sigma(df["lat"].to_numpy(), df["lon"].to_numpy())

    print("\n4) ajuste polinomial (para documentar que no converge)")
    wf = muestrear_frontera(w_de(lat_c, lon_c), 8000)
    pol = {g: ajuste_polinomial(wf, ln_cos_lat(wf), g) for g in (40, 80, 160)}
    for g, r in pol.items():
        print(f"   grado {g}: residuo máximo en la frontera {r:.2e}")

    print("\n5) sensibilidad a la simplificación del contorno")
    from shapely.geometry import Polygon
    from .datos import _poligono_continental
    completo = Polygon(_poligono_continental().exterior)
    sens = []
    for tol in (0.005, 0.01, 0.02):
        lo_s, la_s = np.asarray(completo.simplify(tol, preserve_topology=True).exterior.coords).T
        sens.append({"tolerancia_grados": tol, "vertices": len(la_s),
                     "delta": ProyeccionOptima(la_s, lo_s, h=0.001).delta})
    sens = pd.DataFrame(sens)
    print(sens.to_string(index=False))
    sens.to_csv(tablas / "optima_sensibilidad_contorno.csv", index=False)

    p_env = ProyeccionOptima(*contorno_envolvente(lat_c, lon_c), h=0.001)

    # modelos ajustados para la herramienta del cliente y la Figura 2
    modelos = RAIZ / "resultados" / "modelos"
    modelos.mkdir(parents=True, exist_ok=True)
    p.guardar(modelos / "optima_mexico.npz")
    p_env.guardar(modelos / "optima_envolvente.npz")

    i_min = int(np.argmin(p.h - ln_cos_lat(p.w)))
    lat_min, lon_min = latlon_de(p.w[i_min])
    res = {
        "delta_optima": p.delta,
        "lat_min_sigma": float(lat_min),
        "lon_min_sigma": float(lon_min),
        "error_max_escala_%": error_max_escala(p.delta),
        "sigma_min": np.exp(-p.delta / 2),
        "sigma_max": np.exp(p.delta / 2),
        "ln_sigma_max_interior_sin_c": float(np.max(p.h - ln_cos_lat(p.w))),
        "delta_en_cabeceras": float(ls_cab.max() - ls_cab.min()),
        "verif_err_sigma": e_sigma,
        "verif_sigma_norte_vs_este": e_ne,
        "verif_angulo_grados": e_ang,
        "nodos": len(p.w),
        "nodos_conectados_%": 100 * p.conectados.mean(),
        "segundos": p.segundos,
        "delta_envolvente": p_env.delta,
        "error_max_escala_envolvente_%": error_max_escala(p_env.delta),
        "polinomio_residuo_g40": pol[40],
        "polinomio_residuo_g80": pol[80],
        "polinomio_residuo_g160": pol[160],
    }
    for k, v in res.items():
        print(f"{k:32s} {v:.6g}")
    pd.Series(res).to_csv(tablas / "optima_resumen.csv", header=["valor"])
