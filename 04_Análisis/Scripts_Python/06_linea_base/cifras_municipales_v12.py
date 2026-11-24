#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cifras municipales recalculadas sobre el área geoestadística oficial
=====================================================================

PROPÓSITO
---------
Volver a calcular, sobre el área geoestadística municipal (AGEM) 14 097 del
Marco Geoestadístico 2025 del INEGI, las tres familias de cifras municipales
que hasta ahora se habían obtenido con una capa provisional:

  1. el recuento del corpus dentro del municipio, con su reparto por año,
     el peso de 1971 y los botes que aportan material;
  2. la parte del cuadro de la fotografía de 1971 del par comparativo que
     cae dentro del AGEM, por franjas de norte a sur;
  3. el conjunto A, es decir, las fotografías de 1993 que registran
     territorio municipal, con la media huella medida del fotograma.

Cada cifra se imprime junto a la que daba la capa provisional, de modo que
la comparación quede a la vista sin tener que buscarla en otra parte.

CRITERIOS
---------
Pertenencia al municipio: el centro de toma declarado en el índice de vuelos
cae dentro del AGEM (point-in-polygon). Es el mismo criterio de siempre.

Conjunto A: una fotografía de 1993 registra territorio municipal si su centro
está dentro del AGEM o a no más de una media huella de su borde. La media
huella es la mitad del lado del área de imagen medida en el fotograma,
multiplicada por la resolución en el terreno del ajuste: 1,209 m. No depende
del formato del negativo.

Distancias y áreas en EPSG:32613 (WGS 84 / UTM 13N).

Entrada : índice de vuelos (consulta 07-07) · AGEM oficial ·
          puntos de control del par comparativo
Requiere: numpy, pandas, pyproj, shapely, pyshp
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(RAIZ / "04_Análisis/Scripts_Python/05_regresion_par_estrella"))

from poligonos_v12 import cargar_oficial, cargar_provisional, A_UTM, A_GEO
from validacion_par_estrella_v2 import (
    leer_points, ajustar_similitud, RUTA_1971, FICHA,
)

CSV_CORPUS = (RAIZ / "02_Archivo_INEGI/Metadatos_SINFA"
                   / "07-07_consulta_ÍndicesdeVuelos260107080134.csv")

MEDIA_HUELLA_1993 = 1209.0        # medida en el fotograma (v5b, 14-10)
UMBRAL_DECLARADO_A = 1200.0       # umbral con que se armó el conjunto A (08-07)
INDICES_1971 = [0, 1, 2, 3, 5]    # ajuste definitivo, sin el punto del río

# Cifras publicadas en la bitácora con la capa provisional, para contrastar.
PROVISIONAL = dict(
    en_municipio=226, por_anio={1970: 3, 1971: 155, 1978: 43, 1981: 23, 1993: 2},
    peso_1971=68.6, botes=13, dentro_1993=2,
    conjunto_a=12, area_cuadro=44.72, dentro_cuadro=20.87, pct_cuadro=46.7,
    bandas=[0.0, 5.8, 48.1, 84.7, 94.7],
)


def cargar_corpus():
    df = pd.read_csv(CSV_CORPUS, encoding="latin-1")
    e, n = A_UTM.transform(df["Longitud"].values, df["Latitud"].values)
    df["E"], df["N"] = e, n
    return df


def marca(valor, referencia, unidad="", decimales=0):
    """Imprime el valor nuevo junto al anterior y señala si cambió."""
    igual = round(valor, decimales) == round(referencia, decimales)
    txt = "sin cambio" if igual else f"CAMBIA (antes {referencia:,.{decimales}f}{unidad})"
    return txt


# ---------------------------------------------------------------------------
# 1. El corpus dentro del municipio
# ---------------------------------------------------------------------------

def corpus_en_municipio(df, agem):
    dentro = np.array([agem.contains(Point(x, y)) for x, y in zip(df["E"], df["N"])])
    d = df[dentro]

    print(f"\n  fotogramas con centro dentro del AGEM 14 097   "
          f"{len(d):>6}   {marca(len(d), PROVISIONAL['en_municipio'])}")

    print(f"\n  reparto por año")
    print(f"    {'año':>6}{'fotogramas':>12}{'con la capa provisional':>28}")
    print("    " + "─" * 46)
    cuenta = d["Año"].value_counts().sort_index()
    for anio, n in cuenta.items():
        antes = PROVISIONAL["por_anio"].get(int(anio))
        txt = "—" if antes is None else (f"{antes}" + ("" if antes == n else "  ←"))
        print(f"    {int(anio):>6}{n:>12}{txt:>28}")

    peso = 100 * cuenta.get(1971, 0) / len(d)
    print(f"\n  peso de 1971 en el municipio                   "
          f"{peso:>5.1f} %   {marca(peso, PROVISIONAL['peso_1971'], ' %', 1)}")

    botes = d["Bote"].nunique()
    print(f"  botes con fotogramas en el municipio           "
          f"{botes:>6}   {marca(botes, PROVISIONAL['botes'])}")

    n93 = int((d["Año"] == 1993).sum())
    print(f"  fotogramas de 1993 con centro dentro           "
          f"{n93:>6}   {marca(n93, PROVISIONAL['dentro_1993'])}")
    return d


# ---------------------------------------------------------------------------
# 2. El cuadro de 1971 frente al AGEM
# ---------------------------------------------------------------------------

def huella_1971(agem):
    mapa, fuente, _ = leer_points(RUTA_1971)
    H = ajustar_similitud(fuente[INDICES_1971], mapa[INDICES_1971])
    f = FICHA["1971"]
    W, Hh = f["ancho_px"], f["alto_px"]
    esquinas_px = np.array([[0, 0], [W, 0], [W, -Hh], [0, -Hh]], dtype=float)
    q = np.column_stack([esquinas_px, np.ones(4)]) @ H.T
    esq = q[:, :2] / q[:, 2:3]
    huella = Polygon(esq)

    inter = huella.intersection(agem)
    area, dentro = huella.area / 1e6, inter.area / 1e6
    pct = 100 * dentro / area
    print(f"\n  área del cuadro                     {area:>8.2f} km²   "
          f"{marca(area, PROVISIONAL['area_cuadro'], ' km²', 2)}")
    print(f"  área dentro del AGEM                {dentro:>8.2f} km²   "
          f"{marca(dentro, PROVISIONAL['dentro_cuadro'], ' km²', 2)}")
    print(f"  proporción dentro                   {pct:>8.1f} %     "
          f"{marca(pct, PROVISIONAL['pct_cuadro'], ' %', 1)}")

    ys = esq[:, 1]
    y0, y1 = ys.min(), ys.max()
    x0, x1 = esq[:, 0].min() - 1e4, esq[:, 0].max() + 1e4
    print(f"\n  reparto por franjas de norte a sur (5 bandas de igual altura)")
    print(f"    {'banda':>6}{'latitud media':>16}{'dentro':>10}"
          f"{'con la capa provisional':>28}")
    print("    " + "─" * 60)
    etiquetas = ["1 (N)", "2", "3", "4", "5 (S)"]
    for i in range(5):
        a = y1 - (i + 1) * (y1 - y0) / 5
        b = y1 - i * (y1 - y0) / 5
        banda = huella.intersection(Polygon([(x0, a), (x1, a), (x1, b), (x0, b)]))
        if banda.is_empty or banda.area == 0:
            continue
        p = 100 * banda.intersection(agem).area / banda.area
        _, lat = A_GEO.transform(esq[:, 0].mean(), (a + b) / 2)
        antes = PROVISIONAL["bandas"][i]
        txt = f"{antes:.1f} %" + ("" if abs(antes - p) < 0.05 else "  ←")
        print(f"    {etiquetas[i]:>6}{lat:>16.5f}{p:>9.1f} %{txt:>28}")
    return huella


# ---------------------------------------------------------------------------
# 3. El conjunto A
# ---------------------------------------------------------------------------

def conjunto_a(df, agem):
    d = df[df["Año"] == 1993].copy()
    pts = [Point(x, y) for x, y in zip(d["E"], d["N"])]
    d["dentro"] = [agem.contains(p) for p in pts]
    d["dist"] = [0.0 if agem.contains(p) else p.distance(agem) for p in pts]

    print(f"\n  fotografías de 1993 en el corpus               {len(d):>6}")
    print(f"    con centro dentro del AGEM                   {int(d['dentro'].sum()):>6}")

    print(f"\n    {'umbral':<46}{'metros':>9}{'fotografías':>14}")
    print("    " + "─" * 69)
    for etiqueta, u in [("media huella nominal 1:10,000", 1150.0),
                        ("umbral declarado del conjunto A (08-07)", UMBRAL_DECLARADO_A),
                        ("media huella medida del fotograma", MEDIA_HUELLA_1993)]:
        print(f"    {etiqueta:<46}{u:>9,.0f}{int((d['dist'] <= u).sum()):>14}")

    franja = d[(d["dist"] > 1050) & (d["dist"] <= 1400)].sort_values("dist")
    print(f"\n  fotografías entre 1,050 y 1,400 m del AGEM ({len(franja)})")
    print(f"    {'Bote':>8}{'Línea':>7}{'Foto':>6}{'dist (m)':>10}"
          f"{'latitud':>11}{'longitud':>12}{'¿en el conjunto A?':>21}")
    print("    " + "─" * 75)
    for _, r in franja.iterrows():
        est = "sí" if r["dist"] <= MEDIA_HUELLA_1993 else "no"
        print(f"    {str(r['Bote']):>8}{str(r['Línea']):>7}"
              f"{int(r['Número de foto']):>6}{r['dist']:>10.0f}"
              f"{r['Latitud']:>11.4f}{r['Longitud']:>12.4f}{est:>21}")

    n = int((d["dist"] <= MEDIA_HUELLA_1993).sum())
    print(f"\n  conjunto A con la media huella medida          "
          f"{n:>6}   {marca(n, PROVISIONAL['conjunto_a'])}")
    return d


if __name__ == "__main__":
    print("=" * 92)
    print("CIFRAS MUNICIPALES SOBRE EL AGEM 14 097 DEL MARCO GEOESTADÍSTICO 2025")
    print("=" * 92)

    agem = cargar_oficial(["14097"])["14097"]
    prov = cargar_provisional()
    print(f"\n  AGEM oficial       {agem.area/1e6:>10,.2f} km²")
    print(f"  capa provisional   {prov.area/1e6:>10,.2f} km²"
          f"   (diferencia simétrica "
          f"{agem.symmetric_difference(prov).area/1e6:.3f} km²)")

    df = cargar_corpus()
    print(f"  corpus 07-07       {len(df):>10,} registros")

    print("\n" + "─" * 92)
    print("1. EL CORPUS DENTRO DEL MUNICIPIO")
    print("─" * 92)
    corpus_en_municipio(df, agem)

    print("\n" + "─" * 92)
    print("2. EL CUADRO DE 1971 DEL PAR COMPARATIVO FRENTE AL AGEM")
    print("─" * 92)
    huella_1971(agem)

    print("\n" + "─" * 92)
    print("3. CONJUNTO A — FOTOGRAFÍAS DE 1993 QUE REGISTRAN TERRITORIO MUNICIPAL")
    print("─" * 92)
    conjunto_a(df, agem)

    print("\n" + "=" * 92)
    print("FIN DEL SCRIPT")
    print("=" * 92)
