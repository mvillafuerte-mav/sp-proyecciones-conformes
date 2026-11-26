"""Tabla comparativa y figuras de la wiki.

    python -m proyecciones.figuras

Requiere haber corrido antes distorsion, dirichlet y efecto (lee sus tablas).
Escribe:
  wiki/tablas/tabla_comparativa.tex
  wiki/figuras/fig1_reticula.pdf, fig2_sigma.pdf, fig3_errores.pdf
  (y copias .png en resultados/figuras/)
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from .clasicas import lambert, mercator
from .datos import PROCESSED, RAIZ, cargar_contorno
from .distorsion import malla_interior, sigma_lambert, sigma_mercator

TABLAS = RAIZ / "resultados" / "tablas"
FIG_WIKI = RAIZ / "wiki" / "figuras"
FIG_RES = RAIZ / "resultados" / "figuras"
TEX_TABLAS = RAIZ / "wiki" / "tablas"

ORDEN = ["Mercator", "Lambert INEGI", "Óptima", "Óptima envolvente"]
# paleta categórica validada (orden fijo; pasa CVD en pares adyacentes)
COLOR = dict(zip(ORDEN, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
TINTA, TINTA_2, REJILLA = "#0b0b0b", "#52514e", "#d9d8d4"
# divergente azul <-> rojo con punto medio gris
DIVERGENTE = LinearSegmentedColormap.from_list(
    "div", ["#184f95", "#6da7ec", "#f0efec", "#ef8f8e", "#a8201f"])

plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": TINTA_2, "axes.labelcolor": TINTA, "text.color": TINTA,
    "xtick.color": TINTA_2, "ytick.color": TINTA_2, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "font.family": "serif", "mathtext.fontset": "cm", "savefig.dpi": 300,
})
ANCHO = 6.5   # pulgadas, ancho de texto de la wiki


def _guardar(fig, nombre):
    FIG_WIKI.mkdir(parents=True, exist_ok=True)
    FIG_RES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_WIKI / f"{nombre}.pdf", bbox_inches="tight")
    fig.savefig(FIG_RES / f"{nombre}.png", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------ tabla

def tabla_comparativa():
    clas = pd.read_csv(TABLAS / "distorsion_clasicas.csv").set_index("proyeccion")
    opt = pd.read_csv(TABLAS / "optima_resumen.csv", index_col=0)["valor"]
    dist = pd.read_csv(TABLAS / "efecto_distancias.csv")
    dist = dist[(dist["grupo"] == "Nacional") & (dist["escala"] == "crudo")].set_index("proyeccion")
    km = pd.read_csv(TABLAS / "efecto_kmeans.csv")
    km = km[km["k"] == 20].set_index("proyeccion")["mediana"]
    kn = pd.read_csv(TABLAS / "efecto_vecinos.csv")
    kn = kn[kn["region"].isna()].set_index("proyeccion")["cambian_%"]

    delta = {"Mercator": clas.loc["Mercator", "delta"],
             "Lambert INEGI": clas.loc["Lambert INEGI", "delta"],
             "Óptima": opt["delta_optima"],
             "Óptima envolvente": opt["delta_envolvente"]}
    filas = []
    for p in ORDEN:
        d = delta[p]
        filas.append({"Proyección": p, "delta": d,
                      "esc": 100 * (np.exp(d / 2) - 1),
                      "med": dist.loc[p, "mediana_%"], "p95": dist.loc[p, "p95_%"],
                      "max": dist.loc[p, "max_%"], "km": km[p], "kn": kn[p]})
    t = pd.DataFrame(filas)
    t.to_csv(TABLAS / "tabla_comparativa.csv", index=False)

    nombres = {"Mercator": "Mercator", "Lambert INEGI": "Lambert INEGI",
               "Óptima": "Óptima (México)", "Óptima envolvente": "Óptima (envolvente)"}
    lineas = [
        r"\begin{tabular}{lrrrrrrr}",
        r"\toprule",
        r" & & Error máx. & \multicolumn{3}{c}{Error en distancias (\%)} & k-means & 10 vecinos \\",
        r"\cmidrule(lr){4-6}",
        r"Proyección & $\delta$ & escala (\%) & mediana & P95 & máx. & cambios (\%) & cambian (\%) \\",
        r"\midrule",
    ]
    for _, r in t.iterrows():
        lineas.append(f"{nombres[r['Proyección']]} & {r['delta']:.5f} & {r['esc']:.2f} & "
                      f"{r['med']:.2f} & {r['p95']:.2f} & {r['max']:.2f} & "
                      f"{r['km']:.2f} & {r['kn']:.2f} \\\\")
    lineas += [r"\bottomrule", r"\end{tabular}"]
    TEX_TABLAS.mkdir(parents=True, exist_ok=True)
    (TEX_TABLAS / "tabla_comparativa.tex").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    return t


# ------------------------------------------------------------ figura 1

def fig_reticula():
    lat_c, lon_c = cargar_contorno()
    fig, axs = plt.subplots(1, 2, figsize=(ANCHO, 2.7))
    paneles = [("Mercator", lambda la, lo: mercator(la, lo, lon0=-102.0)),
               ("Lambert INEGI", lambert)]
    for ax, (nombre, f) in zip(axs, paneles):
        x0, y0 = f(23.5, -102.0)
        conv = lambda la, lo: tuple((v - c) / 1000 for v, c in zip(f(la, lo), (x0, y0)))
        lat_l = np.linspace(13, 34, 200)
        for lo in range(-120, -84, 4):
            ax.plot(*conv(lat_l, np.full_like(lat_l, lo)), color=REJILLA, lw=0.7, zorder=1)
        lon_l = np.linspace(-120, -84, 200)
        for la in range(14, 35, 4):
            ax.plot(*conv(np.full_like(lon_l, la), lon_l), color=REJILLA, lw=0.7, zorder=1)
        ax.plot(*conv(lat_c, lon_c), color=COLOR[nombre], lw=1.0, zorder=2)
        ax.set_title(nombre, color=TINTA)
        ax.set_aspect("equal")
        ax.set_xlabel("x (km)")
        ax.tick_params(labelsize=7.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axs[0].set_ylabel("y (km)")
    fig.tight_layout()
    _guardar(fig, "fig1_reticula")


# ------------------------------------------------------------ figura 2

def fig_sigma():
    from .dirichlet import ProyeccionOptima

    lat_c, lon_c = cargar_contorno()
    la, lo = malla_interior(lat_c, lon_c, paso=0.1)
    modelos = RAIZ / "resultados" / "modelos"
    ls = {"Mercator": np.log(sigma_mercator(la)),
          "Lambert INEGI": np.log(sigma_lambert(la)),
          "Óptima": ProyeccionOptima.cargar(modelos / "optima_mexico.npz").ln_sigma(la, lo),
          "Óptima envolvente": ProyeccionOptima.cargar(
              modelos / "optima_envolvente.npz").ln_sigma(la, lo)}
    titulos = {"Mercator": "Mercator", "Lambert INEGI": "Lambert INEGI",
               "Óptima": "Óptima (México)", "Óptima envolvente": "Óptima (envolvente)"}
    fig, axs = plt.subplots(2, 2, figsize=(ANCHO, 4.6))
    for ax, nombre in zip(axs.ravel(), ORDEN):
        v = ls[nombre]
        v = 100 * (np.exp(v - (v.max() + v.min()) / 2) - 1)   # escala centrada, en %
        lim = np.abs(v).max()
        sc = ax.scatter(lo, la, c=v, cmap=DIVERGENTE, vmin=-lim, vmax=lim,
                        s=1.6, marker="s", linewidths=0, rasterized=True)
        ax.plot(lon_c, lat_c, color=TINTA_2, lw=0.5)
        ax.set_aspect(1 / np.cos(np.radians(23.5)))
        ax.set_title(titulos[nombre], color=TINTA)
        ax.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
        cb.set_label(r"$\sigma$ relativa (%)", fontsize=7.5)
        cb.ax.tick_params(labelsize=7)
        cb.outline.set_linewidth(0.4)
    for ax in axs[1]:
        ax.set_xlabel("longitud (°)")
    for ax in axs[:, 0]:
        ax.set_ylabel("latitud (°)")
    fig.tight_layout()
    _guardar(fig, "fig2_sigma")


# ------------------------------------------------------------ figura 3

def fig_errores():
    err = pd.read_csv(PROCESSED / "errores_pares.csv")
    grupos = ["Norte", "Centro", "Sur", "Nacional"]
    fig, ax = plt.subplots(figsize=(ANCHO, 2.9))
    ancho = 0.18
    for g_i, g in enumerate(grupos):
        e = err[err["grupo"] == g]
        for p_i, p in enumerate(ORDEN):
            a = np.abs(e[f"{p}|crudo"].to_numpy())
            a = np.maximum(a, 1e-5)
            q5, q25, q50, q75, q95 = np.percentile(a, [5, 25, 50, 75, 95])
            x = g_i + (p_i - 1.5) * (ancho + 0.02)
            ax.add_patch(plt.Rectangle((x - ancho / 2, q25), ancho, q75 - q25,
                                       facecolor=COLOR[p], edgecolor="none", zorder=3))
            ax.plot([x, x], [q5, q25], color=COLOR[p], lw=1.0, zorder=2)
            ax.plot([x, x], [q75, q95], color=COLOR[p], lw=1.0, zorder=2)
            ax.plot([x - ancho / 2, x + ancho / 2], [q50, q50], color="#fcfcfb", lw=1.2, zorder=4)
            ax.plot(x, a.max(), marker="_", ms=6, mew=1.0, color=COLOR[p], zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 40)
    ax.set_xlim(-0.6, len(grupos) - 0.4)
    ax.set_xticks(range(len(grupos)), grupos)
    ax.set_ylabel("|error| en distancias (%)")
    ax.grid(axis="y", color=REJILLA, lw=0.5, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    nombres = {"Mercator": "Mercator", "Lambert INEGI": "Lambert INEGI",
               "Óptima": "Óptima (México)", "Óptima envolvente": "Óptima (envolvente)"}
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[p]) for p in ORDEN]
    ax.legend(handles, [nombres[p] for p in ORDEN], ncol=4, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, 1.15), fontsize=7.5,
              handlelength=1.2, columnspacing=1.2)
    fig.tight_layout()
    _guardar(fig, "fig3_errores")


if __name__ == "__main__":
    print(tabla_comparativa().round(4).to_string(index=False))
    fig_reticula()
    fig_sigma()
    fig_errores()
    print("figuras listas")
