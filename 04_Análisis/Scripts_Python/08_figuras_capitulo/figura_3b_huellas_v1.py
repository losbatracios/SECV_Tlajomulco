#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figura_3b_huellas_v1.py — variante de la figura 3 con las huellas nominales

La figura 3 muestra las líneas de vuelo de 1992 a 1994 por sus centros de
toma. Con ese criterio, el municipio solo registra dos fotogramas. Esta
variante añade lo que el centro de toma no muestra: la huella nominal de los
fotogramas que alcanzan el área geoestadística municipal 14 097 aunque su
centro quede fuera.

No sustituye a la figura 3 ni la modifica: se guarda aparte, para que [2] y
[1.2] decidan si entra en el capítulo o si el dato va solo en una nota.

Escala de grises, JPG a 300 dpi, 16 cm de ancho (mínimo del libro: 12 cm).
Usa el mismo AGEM y el mismo corpus que linea_base_v8.py, y las huellas de
huella_v1.py, recalculadas aquí con la misma función.

Uso:  python3 figura_3b_huellas_v1.py
"""

import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from shapely.geometry import Point
from shapely.ops import unary_union

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]
sys.path.insert(0, str(RAIZ / "04_Análisis/Scripts_Python/06_linea_base"))
from poligonos_v12 import cargar_oficial, SEIS, NOMBRES          # noqa: E402
from linea_base_v8 import cargar, CASO                           # noqa: E402
from huella_v1 import huella, rumbos_por_linea                   # noqa: E402

SALIDA = RAIZ / "04_Análisis/Visualizaciones/06_figuras_capitulo"
CM = 1 / 2.54
ANCHO = 16 * CM
HALO = [pe.withStroke(linewidth=2.5, foreground="white")]
plt.rcParams.update({"font.family": "serif", "font.size": 8,
                     "axes.edgecolor": "0.6", "savefig.dpi": 300})


def contorno(ax, geom, **kw):
    geoms = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for g in geoms:
        x, y = g.exterior.xy
        ax.plot(x, y, **kw)


def relleno(ax, geom, **kw):
    geoms = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for g in geoms:
        x, y = g.exterior.xy
        ax.fill(x, y, **kw)


def barra_escala(ax, metros, etiqueta):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x = x0 + (x1 - x0) * 0.04
    y = y0 + (y1 - y0) * 0.06
    ax.plot([x, x + metros], [y, y], color="black", lw=2.5,
            solid_capstyle="butt", zorder=20)
    ax.text(x, y + (y1 - y0) * 0.015, etiqueta, fontsize=7.5, zorder=20)


def main():
    agem = cargar_oficial()
    df, _ = cargar(agem)
    caso = agem[CASO]
    km2 = caso.area / 1e6

    d = df[df["Año"].isin([1992, 1993, 1994, 2003])].copy()
    d = d[d["den"].notna()]
    rumbos = rumbos_por_linea(d)

    tocan, cubierto = [], []
    for _, r in d.iterrows():
        h = huella(r["E"], r["N"], r["den"])
        inter = h.intersection(caso)
        if not inter.is_empty:
            tocan.append((h, inter, r))
            cubierto.append(inter)
    union = unary_union(cubierto)
    print(f"  fotogramas cuya huella toca el AGEM: {len(tocan)}")
    print(f"  superficie cubierta: {union.area / 1e6:.2f} km² = "
          f"{union.area / 1e6 / km2 * 100:.2f} % del AGEM")

    fig, ax = plt.subplots(figsize=(ANCHO, ANCHO * 0.92))

    for c in SEIS:
        es = c == CASO
        contorno(ax, agem[c], color="0.15" if es else "0.68",
                 lw=1.4 if es else 0.7, zorder=4 if es else 2)
    relleno(ax, caso, color="0.94", zorder=1)

    # Líneas de vuelo de 1992–1994, para situar las huellas
    dv = df[df["Año"].between(1992, 1994)]
    for _, g in dv.groupby(["Vuelo", "Bote", "Línea", "Escala"]):
        g = g.sort_values("Número de foto")
        ax.plot(g["E"], g["N"], "-", color="0.78", lw=0.4, zorder=3)

    # Huellas que tocan el AGEM y la parte de municipio que registran
    for h, inter, r in tocan:
        contorno(ax, h, color="0.25", lw=0.7, zorder=6)
        relleno(ax, inter, color="0.45", zorder=5)

    centros_dentro = [(r["E"], r["N"]) for _, _, r in tocan
                      if caso.contains(Point(r["E"], r["N"]))]
    if centros_dentro:
        ax.plot([p[0] for p in centros_dentro], [p[1] for p in centros_dentro],
                "o", mfc="white", mec="black", mew=1.0, ms=4, zorder=8)

    for c in SEIS:
        p = agem[c].representative_point()
        ax.text(p.x, p.y, NOMBRES[c], fontsize=7.5, ha="center",
                weight="bold" if c == CASO else "normal",
                color="0.1" if c == CASO else "0.45",
                zorder=9, path_effects=HALO)

    ax.set_aspect("equal")
    ax.axis("off")
    barra_escala(ax, 10000, "10 km")

    ley = [Patch(facecolor="0.45", edgecolor="0.25",
                 label="Superficie del municipio registrada por esas huellas"),
           Line2D([], [], color="0.25", lw=0.7,
                  label="Huella nominal que alcanza el municipio"),
           Line2D([], [], color="0.78", lw=0.8,
                  label="Línea de vuelo de 1992 a 1994"),
           Line2D([], [], color="black", marker="o", mfc="white", ms=4, lw=0,
                  label="Centro de toma dentro del municipio")]
    ax.legend(handles=ley, loc="upper center", bbox_to_anchor=(0.5, 0.02),
              ncol=2, frameon=False, fontsize=7.5)

    # El PNG intermedio va a una carpeta temporal: en la carpeta del proyecto
    # el borrado está deshabilitado y quedaría un archivo suelto.
    temporal = Path(tempfile.mkdtemp())
    png = temporal / "figura_3b.png"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    im = Image.open(png).convert("L")
    im.save(SALIDA / "figura_3b_huellas_1992_2003.jpg", quality=95, dpi=(300, 300))
    shutil.rmtree(temporal, ignore_errors=True)
    print(f"  figura_3b_huellas_1992_2003.jpg  {im.size[0]} × {im.size[1]} px "
          f"= {im.size[0] / 300 * 2.54:.1f} × {im.size[1] / 300 * 2.54:.1f} cm")


if __name__ == "__main__":
    main()
