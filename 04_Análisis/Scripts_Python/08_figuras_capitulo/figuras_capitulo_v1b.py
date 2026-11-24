"""
# v1b (16-11-2026): igual que v1; solo cambia el rótulo de la figura 1 a "Presa La Playa" (convención de nombre del capítulo).
figuras_capitulo_v1.py  [2. TJ Escritura], 2-11-2026

Figuras 1 a 5 del capítulo "Mirada extractiva" (figura 4 sin panel actual; figura 5 solo con el panel de 1971).
Escala de grises, JPG a 300 dpi, 16 cm de ancho (mínimo del libro: 12 cm).

Reutiliza sin cambios el procedimiento de [4]:
  poligonos_v12.cargar_oficial   AGEM del Marco Geoestadístico 2025, EPSG:32613
  linea_base_v8.cargar           corpus SIIV 07-07, municipio por centro de toma
Solo cambia la presentación. Las densidades se cotejan contra la tabla v8.

Insumos ráster:
  figura 4: _insumos/PE_1971_comun.png y PE_1993_comun.png, de recorte_par_estrella_fig4.py
            (misma carpeta), que reproduce los leave-one-out validados (11.48 m y 3.87 m)
  figura 5: 05_humedal_georef/_figuras_capitulo/HU1_1971_recorte_fig5.jpg, recorte en gris
            de HU1_1971_georef.tif reducido 4 veces (2.27 m por píxel)
Requiere geopandas, matplotlib, pyshp y pyproj (se ejecutó en el contenedor, no en la Mac).
Uso: python3 figuras_capitulo_v1.py <raiz_del_proyecto> <carpeta_de_salida> [1 2 3 4 5]
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from PIL import Image
import geopandas as gpd
from pyproj import Transformer
from shapely.geometry import Point

RAIZ = Path(sys.argv[1]); SAL = Path(sys.argv[2]); SAL.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(RAIZ / "04_Análisis/Scripts_Python/06_linea_base"))
from poligonos_v12 import cargar_oficial, SEIS, NOMBRES          # noqa: E402
import linea_base_v8 as lb                                         # noqa: E402

CM = 1 / 2.54
ANCHO = 16 * CM
plt.rcParams.update({"font.family": "serif", "font.size": 8,
                     "axes.edgecolor": "0.6", "savefig.dpi": 300})
HALO = [pe.withStroke(linewidth=2.5, foreground="white")]
CASO = "14097"
A_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32613", always_xy=True)

agem = cargar_oficial()
df, _ = lb.cargar(agem)


def guardar(fig, nombre):
    ruta_png = SAL / f"{nombre}.png"
    fig.savefig(ruta_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    im = Image.open(ruta_png).convert("L")
    im.save(SAL / f"{nombre}.jpg", quality=95, dpi=(300, 300))
    ruta_png.unlink()
    print(f"  {nombre}.jpg  {im.size[0]} x {im.size[1]} px")


def dibujar(ax, geom, **kw):
    gs = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for g in gs:
        x, y = g.exterior.xy
        ax.fill(x, y, **kw)


def escala_y_norte(ax, largo_m, texto, x0f=0.05, y0f=0.06):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    bx = x0 + (x1 - x0) * x0f; by = y0 + (y1 - y0) * y0f
    ax.plot([bx, bx + largo_m], [by, by], color="0.1", lw=2.2, solid_capstyle="butt", zorder=20)
    ax.plot([bx, bx + largo_m / 2], [by, by], color="white", lw=1.0, solid_capstyle="butt", zorder=21)
    ax.text(bx + largo_m / 2, by + (y1 - y0) * 0.018, texto, ha="center", va="bottom",
            fontsize=7, zorder=22, path_effects=HALO)
    nx = x0 + (x1 - x0) * 0.94; ny = y0 + (y1 - y0) * 0.86
    ax.annotate("N", xy=(nx, ny + (y1 - y0) * 0.07), xytext=(nx, ny),
                ha="center", va="center", fontsize=8, weight="bold",
                arrowprops=dict(arrowstyle="-|>", color="0.1", lw=1.0), zorder=22)


def limpiar(ax):
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def encuadre(fig, ax, margen=0.03, ancho_cm=16):
    ax.autoscale_view()
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    dx, dy = x1 - x0, y1 - y0
    ax.set_xlim(x0 - dx * margen, x1 + dx * margen); ax.set_ylim(y0 - dy * margen, y1 + dy * margen)
    ax.set_position([0, 0, 1, 1])
    fig.set_size_inches(ancho_cm * CM, ancho_cm * CM * dy / dx)


ETIQ = {"14039": ("Guadalajara", 0, 0), "14070": ("El Salto", 0, 0),
        "14097": ("Tlajomulco de Zúñiga", 0, 0), "14098": ("San Pedro\nTlaquepaque", 3500, -500),
        "14101": ("Tonalá", 0, 0), "14120": ("Zapopan", 0, 0)}


def rotular(ax, z=9):
    for c in SEIS:
        p = agem[c].representative_point(); t, dx, dy = ETIQ[c]
        ax.text(p.x + dx, p.y + dy, t, ha="center", va="center", linespacing=1.0,
                fontsize=8.5 if c == CASO else 7.5, weight="bold" if c == CASO else "normal",
                color="0.1" if c == CASO else "0.3", zorder=z, path_effects=HALO)


# ---------------------------------------------------------------------------
# Figura 1. Localización
# ---------------------------------------------------------------------------
def figura_1():
    vaso = gpd.read_file(RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef/vaso_1971_afin_v1.gpkg")
    planta = gpd.read_file(RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef/pisa_denue.geojson").to_crs(32613)
    cx, cy = A_UTM.transform(-103.29, 20.53)       # corredor del aeropuerto (par estrella)

    fig, ax = plt.subplots(figsize=(ANCHO, ANCHO * 1.05))
    for c in SEIS:
        es = c == CASO
        dibujar(ax, agem[c], facecolor="0.88" if es else "white",
                edgecolor="0.15" if es else "0.55", lw=1.4 if es else 0.7, zorder=3 if es else 2)
    rotular(ax, z=6)
    for g in vaso.geometry:
        dibujar(ax, g, facecolor="0.25", edgecolor="0.1", lw=0.6, zorder=7)
    vx, vy = vaso.geometry.union_all().centroid.coords[0]
    ax.annotate("Presa La Playa", xy=(vx, vy), xytext=(vx + 6000, vy - 5500),
                fontsize=7.5, ha="center", path_effects=HALO, zorder=9,
                arrowprops=dict(arrowstyle="-", color="0.2", lw=0.7))
    px_, py_ = planta.geometry.iloc[0].coords[0]
    ax.plot(px_, py_, marker="^", ms=5, color="0.05", zorder=8)
    ax.annotate("Planta registrada\nen el DENUE", xy=(px_, py_), xytext=(px_ - 11000, py_ + 4500),
                fontsize=7, ha="center", path_effects=HALO, zorder=9,
                arrowprops=dict(arrowstyle="-", color="0.2", lw=0.7))
    ax.plot(cx, cy, marker="s", ms=7, mfc="white", mec="0.05", mew=1.3, zorder=8)
    ax.annotate("Corredor del\naeropuerto", xy=(cx, cy), xytext=(cx + 5000, cy - 7500),
                fontsize=7.5, ha="center", path_effects=HALO, zorder=9,
                arrowprops=dict(arrowstyle="-", color="0.2", lw=0.7))
    limpiar(ax); encuadre(fig, ax)
    escala_y_norte(ax, 10000, "10 km")
    guardar(fig, "figura_1_localizacion")


# ---------------------------------------------------------------------------
# Figura 2. Densidad de fotogramas por municipio y periodo
# ---------------------------------------------------------------------------
V8 = {  # tabla validada por [4] (línea base v8): fotogramas por km²
    "14039": (0.689, 1.166), "14070": (0.248, 0.216), "14097": (0.334, 0.003),
    "14098": (0.445, 0.617), "14101": (0.386, 0.424), "14120": (0.321, 0.162)}


def figura_2():
    dens = {}
    for c in SEIS:
        km2 = agem[c].area / 1e6
        a = ((df["mun"] == c) & df["Año"].between(1969, 1983)).sum() / km2
        b = ((df["mun"] == c) & df["Año"].between(1992, 1994)).sum() / km2
        dens[c] = (round(a, 3), round(b, 3))
        assert dens[c] == V8[c], (c, dens[c], V8[c])
    print("  densidades replicadas: coinciden con la tabla v8")
    orden = sorted(SEIS, key=lambda c: dens[c][0])
    fig, ax = plt.subplots(figsize=(ANCHO, ANCHO * 0.62))
    y = np.arange(len(orden)); h = 0.36
    b1 = ax.barh(y + h / 2 + 0.02, [dens[c][0] for c in orden], h, color="0.72", label="1969 a 1983")
    b2 = ax.barh(y - h / 2 - 0.02, [dens[c][1] for c in orden], h, color="0.18", label="1992 a 1994")
    for bars in (b1, b2):
        for r in bars:
            v = r.get_width()
            ax.text(v + 0.015, r.get_y() + r.get_height() / 2, f"{v:.3f}", va="center", fontsize=7, color="0.15")
    ax.set_yticks(y)
    ax.set_yticklabels([NOMBRES[c] for c in orden])
    for lab, c in zip(ax.get_yticklabels(), orden):
        if c == CASO:
            lab.set_fontweight("bold")
    ax.set_xlabel("Fotogramas por kilómetro cuadrado (centro de toma dentro del área municipal)", fontsize=7.5)
    ax.set_xlim(0, 1.32)
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=7, colors="0.35")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.xaxis.grid(True, color="0.9", lw=0.5); ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)
    guardar(fig, "figura_2_densidad")


# ---------------------------------------------------------------------------
# Figura 3. Líneas de vuelo de 1992 a 1994
# ---------------------------------------------------------------------------
def figura_3():
    d = df[df["Año"].between(1992, 1994)].copy()
    fig, ax = plt.subplots(figsize=(ANCHO, ANCHO * 1.1))
    for c in SEIS:
        es = c == CASO
        dibujar(ax, agem[c], facecolor="0.9" if es else "white",
                edgecolor="0.15" if es else "0.6", lw=1.4 if es else 0.7, zorder=3 if es else 2)
    destacada = None
    for (vuelo, bote, linea, esc), g in d.groupby(["Vuelo", "Bote", "Línea", "Escala"]):
        g = g.sort_values("Número de foto")
        es713 = "713" in str(bote) and str(linea).strip() in ("1", "1.0") and int(g["Año"].iloc[0]) == 1993
        if es713:
            destacada = g
            continue
        ax.plot(g["E"], g["N"], "-", color="0.35", lw=0.6, zorder=5)
        ax.plot(g["E"], g["N"], ".", color="0.35", ms=1.6, zorder=6)
    assert destacada is not None, "no se encontró la línea 1 del bote 713VE"
    ax.plot(destacada["E"], destacada["N"], "-", color="0.0", lw=1.8, zorder=7)
    ax.plot(destacada["E"], destacada["N"], "o", color="0.0", ms=2.6, zorder=8)
    rotular(ax, z=9)
    from matplotlib.lines import Line2D
    ley = [Line2D([], [], color="0.35", lw=0.6, marker=".", ms=3, label="Línea de vuelo y centros de toma"),
           Line2D([], [], color="0.0", lw=1.8, marker="o", ms=3, label="Línea 1 del bote 713VE (1993)")]
    limpiar(ax); encuadre(fig, ax)
    ax.legend(handles=ley, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2,
              frameon=False, fontsize=7.5)
    escala_y_norte(ax, 10000, "10 km")
    guardar(fig, "figura_3_lineas_1992_1994")
    print(f"  líneas dibujadas: {d.groupby(['Vuelo','Bote','Línea','Escala']).ngroups}")


# ---------------------------------------------------------------------------
# Figura 5, panel de 1971 (recorte en gris del fotograma georreferenciado)
# ---------------------------------------------------------------------------
def figura_5_1971():
    import json
    base = RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef"
    ext = json.load(open(base / "_figuras_capitulo/HU1_1971_recorte_fig5.json"))
    img = np.asarray(Image.open(base / "_figuras_capitulo/HU1_1971_recorte_fig5.jpg"))
    vaso = gpd.read_file(base / "vaso_1971_afin_v1.gpkg")
    planta = gpd.read_file(base / "pisa_denue.geojson").to_crs(32613)
    fig, ax = plt.subplots(figsize=(ANCHO, ANCHO * img.shape[0] / img.shape[1]))
    ax.imshow(img, cmap="gray", vmin=0, vmax=255,
              extent=(ext["xmin"], ext["xmax"], ext["ymin"], ext["ymax"]), zorder=1)
    for g in vaso.geometry:
        gs = g.geoms if g.geom_type == "MultiPolygon" else [g]
        for p in gs:
            x, y = p.exterior.xy
            ax.plot(x, y, color="white", lw=2.6, zorder=3)
            ax.plot(x, y, color="black", lw=1.0, ls=(0, (4, 2)), zorder=4)
    px_, py_ = planta.geometry.iloc[0].coords[0]
    ax.plot(px_, py_, marker="^", ms=7, mfc="white", mec="black", mew=1.2, zorder=5)
    ax.annotate("Planta registrada en el DENUE", xy=(px_, py_), xytext=(px_ + 170, py_), va="center",
                fontsize=7.5, zorder=6, path_effects=HALO)
    ax.text(ext["xmin"] + 150, ext["ymin"] + 450, "Febrero de 1971", fontsize=8.5, weight="bold",
            zorder=6, path_effects=HALO)
    ax.set_xlim(ext["xmin"], ext["xmax"]); ax.set_ylim(ext["ymin"], ext["ymax"])
    limpiar(ax); ax.set_position([0, 0, 1, 1])
    escala_y_norte(ax, 1000, "1 km", x0f=0.40, y0f=0.06)
    guardar(fig, "figura_5_vaso_1971_panel")


# ---------------------------------------------------------------------------
# Figura 4. Par estrella: 1971 y 1993 en la misma retícula
# (insumo de recorte_par_estrella_fig4.py: PE1 similitud con _v2.points, PE2 afín)
# ---------------------------------------------------------------------------
def figura_4():
    import json
    from shapely.geometry import box
    ins = RAIZ / "04_Análisis/Visualizaciones/06_figuras_capitulo/_insumos"
    fig, axs = plt.subplots(1, 2, figsize=(ANCHO, ANCHO * 0.5))
    for ax, anio in zip(axs, ("1971", "1993")):
        e = json.load(open(ins / f"PE_{anio}_comun.json"))
        img = np.asarray(Image.open(ins / f"PE_{anio}_comun.png")).astype(float)
        # estiramiento lineal de contraste (percentiles 1 y 99, sin el borde negro): solo presentación
        lo, hi = np.percentile(img[img > 5], (1, 99))
        img = np.clip((img - lo) / (hi - lo) * 255, 0, 255)
        ext = (e["xmin"], e["xmax"], e["ymin"], e["ymax"])
        L = ext[1] - ext[0]
        for t, fx, fy in (("Tlajomulco de Zúñiga", 0.04, 0.42), ("El Salto", 0.62, 0.80)):
            ax.text(ext[0] + fx * L, ext[2] + fy * L, t, ha="left", va="center", fontsize=7.5,
                    style="italic", zorder=6, path_effects=HALO)
        ax.imshow(img, cmap="gray", vmin=0, vmax=255, extent=ext, zorder=1)
        marco = box(ext[0], ext[2], ext[1], ext[3])
        lim = agem[CASO].boundary.intersection(marco)
        for g in getattr(lim, "geoms", [lim]):
            if g.is_empty:
                continue
            x, y = g.xy
            ax.plot(x, y, color="white", lw=2.4, zorder=3)
            ax.plot(x, y, color="black", lw=0.9, ls=(0, (4, 2)), zorder=4)
        ax.text(ext[0] + 60, ext[3] - 60, f"Febrero de {anio}", ha="left", va="top", fontsize=8.5,
                weight="bold", zorder=6, path_effects=HALO)
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
        limpiar(ax)
    escala_y_norte(axs[1], 500, "500 m", x0f=0.06, y0f=0.06)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0.03)
    guardar(fig, "figura_4_corredor_1971_1993")


if __name__ == "__main__":
    import sys as _s
    pedidas = _s.argv[3:] or ["1", "2", "3", "4", "5"]
    for k in pedidas:
        {"1": figura_1, "2": figura_2, "3": figura_3, "4": figura_4, "5": figura_5_1971}[k]()
