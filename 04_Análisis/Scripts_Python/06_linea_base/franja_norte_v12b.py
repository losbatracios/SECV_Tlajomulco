#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Franja norte del cuadro de 1971 sobre el AGEM oficial — versión corregida
==========================================================================

POR QUÉ EXISTE ESTE SCRIPT
--------------------------
`cifras_municipales_v12.py` calculaba la parte del cuadro de 1971 que cae
dentro del área geoestadística municipal reajustando la transformación con
un conjunto de puntos de control fijado a mano por su posición en el archivo
`.points`. Esa identificación resultó estar equivocada: el punto que se creía
sobre una carretera era el del cauce y viceversa, de modo que el ajuste
retiraba un punto válido y conservaba el que debía quedar fuera.

Aquí la transformación ya no se reajusta: se lee la matriz que produjo la
validación del par comparativo, donde los puntos se identifican por sus
coordenadas y no por su orden. Así hay una sola fuente para el ajuste y este
script solo hace lo suyo, que es cruzarlo con el polígono.

Lo que se recalcula es únicamente la relación del cuadro con el polígono. Las
otras cifras municipales —el recuento del corpus, el reparto por año y el
conjunto A— no dependen de ninguna transformación y no se ven afectadas: se
vuelven a imprimir aquí para dejarlo comprobado.

ESCENARIOS
----------
La validación dejó dos ajustes para 1971, según qué puntos de control entran:

  (i)   5 puntos, el conjunto reducido
  (ii)  7 puntos, los activos del archivo vigente

Se calculan los dos, porque la diferencia entre ellos es la medida de cuánto
depende esta cifra del ajuste. El modelo es el que fijó la regla de elección
para 1971: la similitud (Helmert).

Distancias y áreas en EPSG:32613 (WGS 84 / UTM 13N).

Requiere: numpy, pandas, pyproj, shapely, pyshp
"""

import json
import sys
from pathlib import Path

import numpy as np
from shapely.geometry import Polygon

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]
sys.path.insert(0, str(AQUI))
from poligonos_v12 import cargar_oficial, cargar_provisional, A_GEO

MATRICES = (RAIZ / "04_Análisis/Visualizaciones/03_regresion_par_estrella"
                 / "_validacion_v3/matrices_pixel_a_UTM32613.json")

# Dimensiones del ráster de 1971, de la inspección de calidad del corpus
# (04_Análisis/Visualizaciones/01_inspeccion_calidad/calidad_corpus.csv).
ANCHO_PX, ALTO_PX = 14038, 13617

MODELO = "similitud (Helmert)"      # el que fija la regla de elección para 1971
ESCENARIOS = {"(i) 5 puntos de control": "pe1_i",
              "(ii) 7 puntos de control": "pe1_ii"}

# Cifras previas, para contrastar.
RONDA3 = dict(area=44.72, dentro=20.87, pct=46.7,
              bandas=[0.0, 5.8, 48.1, 84.7, 94.7])     # conjunto desmentido
RONDA4 = dict(pct=46.7, bandas=[0.0, 5.9, 48.2, 84.7, 94.7])


def huella(H):
    esquinas_px = np.array([[0, 0], [ANCHO_PX, 0],
                            [ANCHO_PX, -ALTO_PX], [0, -ALTO_PX]], dtype=float)
    q = np.column_stack([esquinas_px, np.ones(4)]) @ np.asarray(H).T
    return Polygon(q[:, :2] / q[:, 2:3])


def bandas(h, poli, n=5):
    ys = np.array(h.exterior.coords)[:, 1]
    xs = np.array(h.exterior.coords)[:, 0]
    y0, y1 = ys.min(), ys.max()
    x0, x1 = xs.min() - 1e4, xs.max() + 1e4
    out = []
    for i in range(n):
        a = y1 - (i + 1) * (y1 - y0) / n
        b = y1 - i * (y1 - y0) / n
        franja = h.intersection(Polygon([(x0, a), (x1, a), (x1, b), (x0, b)]))
        if franja.is_empty or franja.area == 0:
            out.append((None, None))
            continue
        _, lat = A_GEO.transform(xs.mean(), (a + b) / 2)
        out.append((lat, 100 * franja.intersection(poli).area / franja.area))
    return out


if __name__ == "__main__":
    print("=" * 92)
    print("FRANJA NORTE DEL CUADRO DE 1971 SOBRE EL AGEM 14 097 (VERSIÓN CORREGIDA)")
    print("=" * 92)

    agem = cargar_oficial(["14097"])["14097"]
    prov = cargar_provisional()
    M = json.load(open(MATRICES, encoding="utf-8"))
    print(f"\n  matrices leídas de {MATRICES.relative_to(RAIZ)}")
    print(f"  modelo: {MODELO}")
    print(f"  ráster de 1971: {ANCHO_PX:,} × {ALTO_PX:,} px")

    for etiqueta, clave in ESCENARIOS.items():
        # Las matrices vienen en un diccionario plano, con la clave formada
        # por el escenario y el modelo separados por una barra vertical.
        H = M[f"{clave}|{MODELO}"]
        h = huella(H)
        inter = h.intersection(agem)
        print(f"\n{'─' * 92}")
        print(f"ESCENARIO {etiqueta}")
        print("─" * 92)
        print(f"  área del cuadro                {h.area/1e6:>8.2f} km²")
        print(f"  área dentro del AGEM oficial   {inter.area/1e6:>8.2f} km²"
              f"   ({100*inter.area/h.area:>5.1f} %)")
        pct_prov = 100 * h.intersection(prov).area / h.area
        print(f"  la misma cifra con la capa provisional    "
              f"{pct_prov:>5.1f} %   "
              f"({'igual' if abs(pct_prov - 100*inter.area/h.area) < 0.05 else 'difiere'})")

        print(f"\n  {'banda':>7}{'latitud media':>16}{'dentro':>10}"
              f"{'ronda 3':>11}{'ronda 4':>10}")
        print("  " + "─" * 54)
        etiquetas = ["1 (N)", "2", "3", "4", "5 (S)"]
        for i, (lat, p) in enumerate(bandas(h, agem)):
            if p is None:
                continue
            print(f"  {etiquetas[i]:>7}{lat:>16.5f}{p:>9.1f} %"
                  f"{RONDA3['bandas'][i]:>10.1f} %{RONDA4['bandas'][i]:>9.1f} %")

    print(f"\n{'=' * 92}")
    print("Las cifras de la ronda 3 salieron del conjunto de puntos que el 20-10")
    print("quedó desmentido. Las de este script sustituyen a aquellas.")
    print("=" * 92)
