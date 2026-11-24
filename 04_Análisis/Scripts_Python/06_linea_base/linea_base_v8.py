#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Línea base de cobertura: seis municipios, cuatro periodos
==========================================================

PROPÓSITO
---------
Dar una línea base contra la cual leer la cobertura fotográfica del municipio
del caso. Dos preguntas, separadas:

  1. ¿Cuánto material hay por municipio y por periodo, y a qué escala? Sin
     una referencia vecina, decir que un municipio tiene 226 fotogramas no
     dice nada: puede ser mucho o poco.

  2. En el periodo 1992–1994, ¿dónde empiezan y dónde terminan las líneas de
     vuelo? La pregunta es si el material se interrumpe en el borde del área
     geoestadística municipal de Tlajomulco o si esa coincidencia no existe.

Lo que estas cifras pueden sostener es una relación de registro: el archivo
registra o no registra un territorio en un periodo, y una interrupción
coincide o no coincide con un borde. No sostienen ninguna intención ni
ninguna causa. Los verbos admitidos son registra, no registra y coincide con.

PERIODOS
--------
Se agrupan por tramos del propio corpus, no por una periodización externa:

  1969–1983   el grueso del material, de CETENAL y del primer INEGI
  1984–1991   el hueco: el corpus no registra ningún fotograma
  1992–1994   el material posterior a la reforma del artículo 27
  2003        la última cobertura del corpus

CRITERIOS
---------
Pertenencia al municipio: el centro de toma declarado en el índice de vuelos
cae dentro del área geoestadística municipal (point-in-polygon). Un fotograma
cuyo centro cae fuera registra territorio del municipio si lo alcanza su
huella, pero esa corrección solo está medida para el par comparativo, así
que aquí se cuenta por centro y se dice así.

Identidad del fotograma: el identificador de archivo (bote, rollo, línea y
número de foto) no siempre distingue un fotograma — hay al menos un caso de
dos registros con la misma llave y distinta longitud —, de modo que la
identidad se toma como la llave más las coordenadas.

Línea de vuelo: la cuaterna vuelo, bote, línea y escala. La escala hace falta:
hay botes que reutilizan la numeración de líneas en dos bloques de distinta
altura de vuelo, y sin ella dos series separadas se fundirían en una sola
traza. Sus fotogramas se ordenan por número de foto, y sus extremos son el
primero y el último de esa secuencia.

Distancias y áreas en EPSG:32613 (WGS 84 / UTM 13N).

SALIDAS
-------
  linea_base_municipio_periodo.csv   la tabla del punto 1
  lineas_vuelo_1992_1994.csv         una fila por línea de vuelo
  mapa_lineas_1992_1994.png          mapa en escala de grises

Requiere: numpy, pandas, pyproj, shapely, pyshp, matplotlib
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Point, LineString

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]
sys.path.insert(0, str(AQUI))
from poligonos_v12 import cargar_oficial, SEIS, NOMBRES, A_UTM

CSV_CORPUS = (RAIZ / "02_Archivo_INEGI/Metadatos_SINFA"
                   / "07-07_consulta_ÍndicesdeVuelos260107080134.csv")
SALIDA = RAIZ / "04_Análisis/Visualizaciones/04_linea_base"

PERIODOS = [("1969–1983", 1969, 1983),
            ("1984–1991", 1984, 1991),
            ("1992–1994", 1992, 1994),
            ("2003", 2003, 2003)]

CASO = "14097"            # el municipio del caso
LLAVE = ["Bote", "Rollo", "Línea", "Número de foto", "Año"]


def denominador(x):
    """Denominador de una escala escrita como 1:N."""
    m = re.search(r"(\d[\d,\.]*)\s*$", str(x).replace(" ", ""))
    return int(m.group(1).replace(",", "").replace(".", "")) if m else None


def cargar(agem):
    df = pd.read_csv(CSV_CORPUS, encoding="latin-1")
    e, n = A_UTM.transform(df["Longitud"].values, df["Latitud"].values)
    df["E"], df["N"] = e, n
    df["den"] = df["Escala"].map(denominador)

    # Identidad: llave más coordenadas. Se informa cuántas llaves se repiten.
    rep = df.duplicated(subset=LLAVE, keep=False)
    if rep.any():
        g = df[rep].groupby(LLAVE).size()
        print(f"  llaves repetidas en el corpus: {len(g)}"
              f" ({int(rep.sum())} registros); se distinguen por coordenadas")
        for k, n_ in g.items():
            print(f"    {k}  ×{n_}")
    df["id"] = (df[LLAVE].astype(str).agg("|".join, axis=1)
                + "|" + df["Latitud"].round(6).astype(str)
                + "|" + df["Longitud"].round(6).astype(str))
    if df["id"].duplicated().any():
        print(f"  AVISO: {int(df['id'].duplicated().sum())} registros siguen"
              f" sin distinguirse ni con coordenadas")

    # Municipio de cada fotograma, por centro de toma.
    puntos = [Point(x, y) for x, y in zip(df["E"], df["N"])]
    df["mun"] = ""
    for c in SEIS:
        dentro = np.array([agem[c].contains(p) for p in puntos])
        df.loc[dentro, "mun"] = c
    return df, puntos


# ---------------------------------------------------------------------------
# 1. Tabla por municipio y periodo
# ---------------------------------------------------------------------------

def tabla(df, agem):
    filas = []
    for c in SEIS:
        km2 = agem[c].area / 1e6
        for etiqueta, a, b in PERIODOS:
            d = df[(df["mun"] == c) & (df["Año"].between(a, b))]
            fina = int(d["den"].min()) if len(d) and d["den"].notna().any() else None
            filas.append(dict(cvegeo=c, municipio=NOMBRES[c], km2=round(km2, 2),
                              periodo=etiqueta, fotogramas=len(d),
                              por_km2=round(len(d) / km2, 4),
                              escala_mas_fina=(f"1:{fina:,}" if fina else "—")))
    t = pd.DataFrame(filas)

    print(f"\n  {'municipio':24}{'km²':>10}", end="")
    for etiqueta, _, _ in PERIODOS:
        print(f"{etiqueta:>22}", end="")
    print()
    print("  " + "─" * (34 + 22 * len(PERIODOS)))
    for c in SEIS:
        sub = t[t.cvegeo == c]
        marca = " ←" if c == CASO else "  "
        print(f"  {NOMBRES[c]:22}{marca}{sub.km2.iloc[0]:>10,.0f}", end="")
        for etiqueta, _, _ in PERIODOS:
            r = sub[sub.periodo == etiqueta].iloc[0]
            celda = f"{r.fotogramas:>5}  {r.por_km2:>6.3f}/km²"
            print(f"{celda:>22}", end="")
        print()

    print(f"\n  escala más fina registrada por municipio y periodo")
    print(f"  {'municipio':24}", end="")
    for etiqueta, _, _ in PERIODOS:
        print(f"{etiqueta:>14}", end="")
    print()
    print("  " + "─" * (26 + 14 * len(PERIODOS)))
    for c in SEIS:
        sub = t[t.cvegeo == c]
        print(f"  {NOMBRES[c]:24}", end="")
        for etiqueta, _, _ in PERIODOS:
            print(f"{sub[sub.periodo == etiqueta].escala_mas_fina.iloc[0]:>14}", end="")
        print()

    SALIDA.mkdir(parents=True, exist_ok=True)
    t.to_csv(SALIDA / "linea_base_municipio_periodo.csv", index=False,
             encoding="utf-8")
    fuera = int((df["mun"] == "").sum())
    print(f"\n  fotogramas del corpus fuera de los seis municipios: {fuera:,}"
          f" de {len(df):,}")
    return t


# ---------------------------------------------------------------------------
# 2. Geometría de las líneas de vuelo de 1992–1994
# ---------------------------------------------------------------------------

def lineas_de_vuelo(df, agem):
    d = df[df["Año"].between(1992, 1994)].copy()
    print(f"\n  fotogramas de 1992–1994 en el corpus: {len(d):,}")

    caso = agem[CASO]
    otros = {c: agem[c] for c in SEIS if c != CASO}

    def donde(p):
        """Municipio que contiene el punto, o el más próximo de los seis."""
        for c in SEIS:
            if agem[c].contains(p):
                return NOMBRES[c], 0.0
        cerca = min(SEIS, key=lambda c: p.distance(agem[c]))
        return f"cerca de {NOMBRES[cerca]}", p.distance(agem[cerca])

    filas = []
    for (vuelo, bote, linea, escala), g in d.groupby(
            ["Vuelo", "Bote", "Línea", "Escala"]):
        g = g.sort_values("Número de foto")
        pts = [Point(x, y) for x, y in zip(g["E"], g["N"])]
        ini, fin = pts[0], pts[-1]
        traza = LineString([(p.x, p.y) for p in pts]) if len(pts) > 1 else None

        d_ini = 0.0 if caso.contains(ini) else ini.distance(caso)
        d_fin = 0.0 if caso.contains(fin) else fin.distance(caso)
        cruza = bool(traza is not None and traza.intersects(caso))
        dentro = int(sum(caso.contains(p) for p in pts))
        d_min = min(0.0 if caso.contains(p) else p.distance(caso) for p in pts)

        filas.append(dict(
            vuelo=vuelo, bote=bote, linea=linea, anio=int(g["Año"].iloc[0]),
            escala=escala, n=len(g),
            foto_ini=int(g["Número de foto"].iloc[0]),
            foto_fin=int(g["Número de foto"].iloc[-1]),
            lat_ini=round(g["Latitud"].iloc[0], 4),
            lon_ini=round(g["Longitud"].iloc[0], 4),
            lat_fin=round(g["Latitud"].iloc[-1], 4),
            lon_fin=round(g["Longitud"].iloc[-1], 4),
            largo_km=round(traza.length / 1000, 2) if traza else 0.0,
            dist_ini_m=round(d_ini), dist_fin_m=round(d_fin),
            dist_min_m=round(d_min), fotos_dentro=dentro,
            traza_cruza_el_agem=cruza,
            extremo_inicial_en=donde(ini)[0], extremo_final_en=donde(fin)[0],
        ))
    t = pd.DataFrame(filas).sort_values(["anio", "vuelo", "bote", "linea"])
    t.to_csv(SALIDA / "lineas_vuelo_1992_1994.csv", index=False, encoding="utf-8")

    print(f"  líneas de vuelo distintas (vuelo · bote · línea · escala): {len(t)}")
    print(f"\n  {'vuelo':22}{'bote':>7}{'lín':>5}{'año':>6}{'esc.':>10}{'n':>4}"
          f"{'largo':>8}{'d.ini':>8}{'d.fin':>8}{'d.mín':>8}{'dentro':>7}{'cruza':>7}")
    print("  " + "─" * 100)
    for _, r in t.iterrows():
        print(f"  {str(r.vuelo)[:21]:22}{r.bote:>7}{r.linea:>5}{r.anio:>6}"
              f"{str(r.escala):>10}{r.n:>4}{r.largo_km:>7.1f}k"
              f"{r.dist_ini_m:>8,}{r.dist_fin_m:>8,}{r.dist_min_m:>8,}"
              f"{r.fotos_dentro:>7}{('sí' if r.traza_cruza_el_agem else 'no'):>7}")

    print(f"\n  d.ini, d.fin y d.mín: distancia en metros al AGEM 14 097"
          f" (0 = dentro).")

    # ¿Las líneas terminan en el borde del AGEM del caso?
    print(f"\n  Relación de cada línea con el AGEM 14 097")
    print(f"    {'situación':58}{'líneas':>8}")
    print("    " + "─" * 66)
    cruzan = t[t.traza_cruza_el_agem]
    tocan = t[(~t.traza_cruza_el_agem) & (t.dist_min_m <= 2000)]
    lejos = t[(~t.traza_cruza_el_agem) & (t.dist_min_m > 2000)]
    print(f"    {'la traza cruza el AGEM':58}{len(cruzan):>8}")
    print(f"    {'no lo cruza, pero pasa a 2 km o menos':58}{len(tocan):>8}")
    print(f"    {'no se acerca a menos de 2 km':58}{len(lejos):>8}")

    if len(cruzan):
        print(f"\n  De las que cruzan el AGEM, dónde está su último fotograma")
        print(f"    {'vuelo · bote · línea':34}{'fotos dentro':>14}"
              f"{'extremo final':>30}")
        print("    " + "─" * 78)
        for _, r in cruzan.iterrows():
            print(f"    {str(r.vuelo)[:14]} · {r.bote} · {r.linea:<10}"
                  f"{r.fotos_dentro:>14}{r.extremo_final_en:>30}")

    # ¿Hay un bloque de vuelo que explique el corte? Un bloque es un bote a
    # una altura de vuelo: la unidad con que se programa una cobertura.
    print(f"\n  Bloques de vuelo en 1992–1994 (vuelo · bote · escala)")
    print(f"    {'vuelo · bote':26}{'escala':>9}{'año':>6}{'líneas':>8}"
          f"{'fotos':>7}{'en 14 097':>11}{'d.mín al AGEM':>15}{'centro del bloque':>22}")
    print("    " + "─" * 104)
    for (vuelo, bote, escala), g in t.groupby(["vuelo", "bote", "escala"]):
        lat = (g.lat_ini.mean() + g.lat_fin.mean()) / 2
        lon = (g.lon_ini.mean() + g.lon_fin.mean()) / 2
        print(f"    {str(vuelo)[:12] + ' · ' + str(bote):26}{str(escala):>9}"
              f"{int(g.anio.iloc[0]):>6}{len(g):>8}{int(g.n.sum()):>7}"
              f"{int(g.fotos_dentro.sum()):>11}{int(g.dist_min_m.min()):>13,} m"
              f"{lat:>11.3f},{lon:>9.3f}")
    return t, d


# ---------------------------------------------------------------------------
# 3. Mapa
# ---------------------------------------------------------------------------

def mapa(d, agem, t):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 9))

    def dibujar(geom, **kw):
        gs = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for g in gs:
            x, y = g.exterior.xy
            ax.plot(x, y, **kw)

    for c in SEIS:
        es_caso = c == CASO
        dibujar(agem[c], color="0.25" if es_caso else "0.65",
                lw=2.0 if es_caso else 0.9,
                zorder=3 if es_caso else 2)
        p = agem[c].representative_point()
        ax.annotate(NOMBRES[c], (p.x, p.y), ha="center", fontsize=7.5,
                    color="0.2" if es_caso else "0.5",
                    weight="bold" if es_caso else "normal", zorder=4)

    for _, g in d.groupby(["Vuelo", "Bote", "Línea", "Escala"]):
        g = g.sort_values("Número de foto")
        ax.plot(g["E"], g["N"], "-", color="0.1", lw=0.7, alpha=0.85, zorder=5)
        ax.plot(g["E"], g["N"], ".", color="0.1", ms=2.0, zorder=6)

    ax.set_aspect("equal")
    ax.set_title("Líneas de vuelo de 1992–1994 sobre las áreas geoestadísticas\n"
                 "municipales de la ZMG (Marco Geoestadístico 2025, INEGI)",
                 fontsize=10.5, color="0.15")
    ax.set_xlabel("Este, EPSG:32613 (m)", fontsize=8, color="0.35")
    ax.set_ylabel("Norte, EPSG:32613 (m)", fontsize=8, color="0.35")
    ax.tick_params(labelsize=7, colors="0.45")
    for s in ax.spines.values():
        s.set_color("0.8")
    ax.grid(True, color="0.92", lw=0.5)

    ruta = SALIDA / "mapa_lineas_1992_1994.png"
    fig.savefig(ruta, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n  mapa guardado: {ruta.relative_to(RAIZ)}")


if __name__ == "__main__":
    print("=" * 100)
    print("LÍNEA BASE DE COBERTURA — SEIS MUNICIPIOS DE LA ZMG, CUATRO PERIODOS")
    print("=" * 100)

    agem = cargar_oficial()
    print(f"\n  AGEM del Marco Geoestadístico 2025 (INEGI), EPSG:32613")
    df, _ = cargar(agem)

    print("\n" + "─" * 100)
    print("1. FOTOGRAMAS POR MUNICIPIO Y PERIODO (centro de toma dentro del AGEM)")
    print("─" * 100)
    tabla(df, agem)

    print("\n" + "─" * 100)
    print("2. GEOMETRÍA DE LAS LÍNEAS DE VUELO DE 1992–1994")
    print("─" * 100)
    t, d = lineas_de_vuelo(df, agem)

    print("\n" + "─" * 100)
    print("3. MAPA")
    print("─" * 100)
    mapa(d, agem, t)

    print("\n" + "=" * 100)
    print("FIN DEL SCRIPT")
    print("=" * 100)
