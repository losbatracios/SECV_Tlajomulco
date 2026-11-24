#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Qué registra el vaso de La Playa: el archivo aéreo y la serie satelital
========================================================================

PROPÓSITO
---------
Dos preguntas distintas sobre el mismo lugar, las dos en términos de registro:

  A. ¿Qué fotogramas del archivo cubren el vaso, en qué años y a qué escala?
     Es la operación sobre la mirada: qué miró el Estado y cuándo.
  B. ¿En qué meses registra agua superficial abierta la serie satelital?
     Es cronología sobre el territorio.

Ninguna de las dos atribuye causa. Los verbos son registra, no registra,
coincide con, precede a y es posterior a.

A. COBERTURA DEL ARCHIVO
------------------------
Para cada fotograma del índice de vuelos se construye su **huella nominal**:
un cuadrado centrado en el centro de toma declarado, de medio lado igual a
0.115 m (la mitad del formato de 23 cm) por el denominador de la escala
nominal. Un fotograma a 1:35,000 cubre así algo más de 8 km de lado.

Es una huella nominal, no medida: supone la escala declarada, no corrige la
inclinación de la toma ni el relieve, y se orienta con los ejes de la
proyección en vez de con el rumbo del vuelo. Sirve para decir qué fotogramas
*pueden* contener el vaso, no para medir cuánto de él sale en cada uno. Donde
el resultado depende de un borde, se dice.

B. SERIE SATELITAL
------------------
Se replican, desde el CSV crudo exportado del archivo de agua superficial,
los meses con observación suficiente y el agua registrada en cada uno. Un mes
mal observado no es un mes sin agua: se cuenta aparte. La serie detecta agua
abierta a 30 m, de modo que un humedal con vegetación emergente, agua somera
o canales más estrechos puede no quedar registrado. Por eso se escribe "no
registra agua superficial abierta", nunca "no había humedal".

Entrada : índice de vuelos (consulta 07-07) · `vaso_1971.gpkg` y los
          contornos por modelo · CSV mensual del archivo de agua superficial
CRS     : EPSG:32613 (WGS 84 / UTM 13N)
Salidas : `05_humedal_georef/_validacion/`
Requiere: numpy, pandas, pyproj, shapely, pyogrio
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
from pyogrio.raw import read as ogr_read
from shapely import wkb
from shapely.geometry import Point, Polygon, box

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]

CSV_CORPUS = (RAIZ / "02_Archivo_INEGI/Metadatos_SINFA"
                   / "07-07_consulta_ÍndicesdeVuelos260107080134.csv")
VIS = RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef"
SALIDA = VIS / "_validacion"
CSV_JRC = RAIZ / "05_Datos_abiertos/JRC_GSW/jrc_gsw_monthly_bbox_v1.csv"

LADO_NEGATIVO_M = 0.23
A_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32613", always_xy=True)

# Caja del vaso estimada el 25-08, que es el área sobre la que se extrajo la
# serie satelital y la que usó la tabla de la solicitud para el artículo.
CAJA = dict(sur=20.4647, norte=20.4796, oeste=-103.5442, este=-103.5100)

# Los 7 fotogramas de la tabla e5 de la Hoja de ruta (solicitud para el
# artículo). Se recalcula su cobertura con la huella nominal.
E5 = [
    ("1978", "96VE", 284, 7, 15, 30000, 20.4578, -103.522, "1"),
    ("1981", "1121VE", 714, 3, 18, 40000, 20.4627, -103.526, "1"),
    ("1978", "96VE", 284, 7, 14, 30000, 20.4595, -103.501, "2"),
    ("1978", "96VE", 284, 7, 16, 30000, 20.4579, -103.546, "2"),
    ("1981", "1121VE", 714, 3, 17, 40000, 20.4392, -103.526, "2"),
    ("1981", "1121VE", 714, 3, 19, 40000, 20.4861, -103.525, "2"),
    ("1978", "96VE", 284, 6, 3, 30000, 20.5027, -103.530, "3"),
]

UMBRAL_BIEN_OBSERVADO = 0.90
UMBRAL_PARCIAL = 0.50


def denominador(x):
    m = re.search(r"(\d[\d,\.]*)\s*$", str(x).replace(" ", ""))
    return int(m.group(1).replace(",", "").replace(".", "")) if m else None


def huella(e, n, den):
    """Cuadrado de medio lado 0.115 × denominador, centrado en el punto."""
    h = LADO_NEGATIVO_M / 2 * den
    return box(e - h, n - h, e + h, n + h)


def caja_proyectada(n=200):
    lon = np.linspace(CAJA["oeste"], CAJA["este"], n)
    lat = np.linspace(CAJA["sur"], CAJA["norte"], n)
    borde = ([(x, CAJA["sur"]) for x in lon] + [(CAJA["este"], y) for y in lat]
             + [(x, CAJA["norte"]) for x in lon[::-1]]
             + [(CAJA["oeste"], y) for y in lat[::-1]])
    e, nn = A_UTM.transform([p[0] for p in borde], [p[1] for p in borde])
    return Polygon(zip(e, nn))


def leer_poligono(ruta):
    _, _, g, _ = ogr_read(str(ruta))
    return wkb.loads(g[0])


# ---------------------------------------------------------------------------
# A. Cobertura del archivo
# ---------------------------------------------------------------------------

def cobertura(vaso, caja):
    centroide = vaso.centroid
    lon_c, lat_c = Transformer.from_crs("EPSG:32613", "EPSG:4326",
                                        always_xy=True).transform(centroide.x,
                                                                  centroide.y)
    print(f"\n  centroide del vaso: E {centroide.x:,.0f} N {centroide.y:,.0f}"
          f"  ({lat_c:.4f}, {lon_c:.4f})")

    df = pd.read_csv(CSV_CORPUS, encoding="latin-1")
    e, n = A_UTM.transform(df["Longitud"].values, df["Latitud"].values)
    df["E"], df["N"] = e, n
    df["den"] = df["Escala"].map(denominador)
    df = df[df["den"].notna()].copy()
    df["medio_lado"] = LADO_NEGATIVO_M / 2 * df["den"]
    dx = (df["E"] - centroide.x).abs()
    dy = (df["N"] - centroide.y).abs()
    df["cubre"] = (dx <= df["medio_lado"]) & (dy <= df["medio_lado"])
    df["dist_centro_m"] = np.hypot(dx, dy)

    cubren = df[df["cubre"]].copy()
    print(f"\n  fotogramas del corpus: {len(df):,}")
    print(f"  fotogramas cuya huella nominal contiene el centroide del vaso:"
          f" {len(cubren)}")

    print(f"\n  por año y escala")
    print(f"    {'año':>6}{'escala':>12}{'fotogramas':>12}")
    print("    " + "─" * 30)
    for (a, d), g in cubren.groupby(["Año", "den"]):
        print(f"    {int(a):>6}{'1:' + format(int(d), ','):>12}{len(g):>12}")
    anios = sorted(cubren["Año"].unique())
    print(f"\n  años con cobertura: {', '.join(str(int(a)) for a in anios)}")
    print(f"  el último año en que el archivo registra el vaso es {int(max(anios))}")

    print(f"\n  detalle de los fotogramas que cubren el centroide")
    print(f"    {'bote':>9}{'rollo':>8}{'lín':>5}{'foto':>6}{'año':>6}{'escala':>10}"
          f"{'dist. al centroide':>20}{'% del vaso':>12}{'% de la caja':>14}")
    print("    " + "─" * 92)
    filas = []
    for _, r in cubren.sort_values(["Año", "Bote", "Línea",
                                    "Número de foto"]).iterrows():
        h = huella(r["E"], r["N"], r["den"])
        pv = 100 * h.intersection(vaso).area / vaso.area
        pc = 100 * h.intersection(caja).area / caja.area
        print(f"    {str(r['Bote']):>9}{str(r['Rollo']):>8}{str(r['Línea']):>5}"
              f"{int(r['Número de foto']):>6}{int(r['Año']):>6}"
              f"{'1:' + format(int(r['den']), ','):>10}"
              f"{r['dist_centro_m']:>17,.0f} m{pv:>11.0f} %{pc:>13.0f} %")
        filas.append(dict(bote=r["Bote"], rollo=r["Rollo"], linea=r["Línea"],
                          foto=int(r["Número de foto"]), anio=int(r["Año"]),
                          escala=f"1:{int(r['den']):,}",
                          dist_centroide_m=round(r["dist_centro_m"]),
                          pct_vaso=round(pv, 1), pct_caja=round(pc, 1)))
    pd.DataFrame(filas).to_csv(SALIDA / "cobertura_vaso_corpus.csv", index=False,
                               encoding="utf-8")

    # Años en que el archivo registra algo del vaso aunque no su centroide.
    df["toca_vaso"] = [huella(x, y, d).intersects(vaso)
                       for x, y, d in zip(df["E"], df["N"], df["den"])]
    tocan = df[df["toca_vaso"]]
    print(f"\n  fotogramas cuya huella nominal toca alguna parte del vaso:"
          f" {len(tocan)}, en {', '.join(str(int(a)) for a in sorted(tocan['Año'].unique()))}")
    return cubren, df


def cobertura_e5(vaso, caja):
    print(f"\n  los 7 fotogramas de la tabla de la solicitud para el artículo")
    print(f"    {'prio':>5}{'año':>6}{'bote':>9}{'lín':>5}{'foto':>6}{'escala':>10}"
          f"{'¿centroide?':>13}{'% del vaso':>12}{'% de la caja':>14}{'e5 decía':>11}")
    print("    " + "─" * 92)
    e5_dice = {(1978, 15): "100 %", (1981, 18): "100 %", (1978, 14): "~70 %",
               (1978, 16): "~90 %", (1981, 17): "100 %", (1981, 19): "100 %",
               (1978, 3): "~55 %"}
    filas = []
    for anio, bote, rollo, linea, foto, den, lat, lon, prio in E5:
        x, y = A_UTM.transform(lon, lat)
        h = huella(x, y, den)
        pv = 100 * h.intersection(vaso).area / vaso.area
        pc = 100 * h.intersection(caja).area / caja.area
        dentro = h.contains(vaso.centroid)
        print(f"    {prio:>5}{anio:>6}{bote:>9}{linea:>5}{foto:>6}"
              f"{'1:' + format(den, ','):>10}{'sí' if dentro else 'no':>13}"
              f"{pv:>11.0f} %{pc:>13.0f} %{e5_dice.get((int(anio), foto), '—'):>11}")
        filas.append(dict(prioridad=prio, anio=anio, bote=bote, rollo=rollo,
                          linea=linea, foto=foto, escala=f"1:{den:,}",
                          contiene_centroide=dentro, pct_vaso=round(pv, 1),
                          pct_caja=round(pc, 1)))
    pd.DataFrame(filas).to_csv(SALIDA / "cobertura_vaso_e5.csv", index=False,
                               encoding="utf-8")
    return filas


# ---------------------------------------------------------------------------
# B. Serie satelital
# ---------------------------------------------------------------------------

def serie_satelital():
    df = pd.read_csv(CSV_JRC)
    df["total_px"] = df["px_water"] + df["px_no_water"] + df["px_no_obs"]
    total = df["total_px"].max()
    ha_px = 30 * 30 / 1e4
    orden = df.sort_values(["year", "month"])
    pri, ult = orden.iloc[0], orden.iloc[-1]
    print(f"\n  meses en el archivo: {len(df)}  "
          f"({int(pri['year'])}-{int(pri['month']):02d}"
          f" a {int(ult['year'])}-{int(ult['month']):02d})")
    print(f"  píxeles de la caja: {total:,.1f} de 30 m  →  "
          f"área de la caja {total * ha_px:,.0f} ha")
    print(f"  El README decía ~568 ha: la cifra correcta es"
          f" {total * ha_px:,.0f} ha (corregir).")

    df["obs"] = (df["px_water"] + df["px_no_water"]) / total
    bien = df[df["obs"] >= UMBRAL_BIEN_OBSERVADO]
    parcial = df[df["obs"] >= UMBRAL_PARCIAL]
    print(f"\n  meses con al menos el {UMBRAL_BIEN_OBSERVADO:.0%} de la caja"
          f" observada: {len(bien)} de {len(df)}")
    print(f"  meses con al menos el {UMBRAL_PARCIAL:.0%}: {len(parcial)} de {len(df)}")

    con_agua = df[df["px_water"] > 0]
    print(f"\n  meses con algún píxel de agua: {len(con_agua)}")
    if len(con_agua):
        print(f"    {'mes':>9}{'píxeles':>10}{'área':>12}{'% de la caja':>14}"
              f"{'observado':>12}")
        print("    " + "─" * 57)
        for _, r in con_agua.sort_values(["year", "month"]).iterrows():
            print(f"    {int(r['year'])}-{int(r['month']):02d}"
                  f"{r['px_water']:>10,.0f}{r['px_water'] * ha_px:>9,.1f} ha"
                  f"{100 * r['px_water'] / total:>13.1f} %{r['obs']:>11.0%}")
    antes = con_agua[con_agua["year"] < 2021]
    print(f"\n  meses con agua antes de 2021: {len(antes)}"
          f"  →  {'confirmado: ninguno' if len(antes) == 0 else 'ATENCIÓN: los hay'}")
    ep = con_agua[con_agua["year"] == 2021]
    if len(ep):
        meses = sorted(int(m) for m in ep["month"])
        print(f"  episodio de 2021: meses {meses[0]} a {meses[-1]}"
              f" ({len(meses)} meses), entre {ep['px_water'].min() * ha_px:,.1f}"
              f" y {ep['px_water'].max() * ha_px:,.1f} ha"
              f" (máximo {100 * ep['px_water'].max() / total:.1f} % de la caja)")

    print(f"\n  meses bien observados antes de 1993 (cuando la serie aún es rala)")
    temprano = bien[bien["year"] < 1993].sort_values(["year", "month"])
    print(f"    {', '.join(f'{int(r.year)}-{int(r.month):02d}' for _, r in temprano.iterrows())}")
    feb = bien[bien["month"] == 2]
    print(f"\n  febreros bien observados (el mismo mes que la fotografía de 1971):"
          f" {len(feb)}; agua registrada en {int((feb['px_water'] > 0).sum())}")

    primer = bien.sort_values(["year", "month"]).iloc[0]
    print(f"\n  Lectura, en términos de registro: la serie no registra agua")
    print(f"  superficial abierta en la caja del vaso en ninguno de los"
          f" {len(bien)} meses bien")
    print(f"  observados entre {int(primer['year'])}-{int(primer['month']):02d}"
          f" y 2020. La desaparición de la lámina que registra la fotografía de")
    print(f"  febrero de 1971 precede al primer mes observado: la serie no puede")
    print(f"  fecharla, solo acotarla entre febrero de 1971 y"
          f" {int(primer['year'])}-{int(primer['month']):02d}.")

    df[["year", "month", "px_water", "px_no_water", "px_no_obs", "obs"]].to_csv(
        SALIDA / "jrc_replica_mensual.csv", index=False, encoding="utf-8")
    return df, bien, con_agua


if __name__ == "__main__":
    SALIDA.mkdir(parents=True, exist_ok=True)
    print("=" * 94)
    print("QUÉ REGISTRA EL VASO DE LA PLAYA: EL ARCHIVO AÉREO Y LA SERIE SATELITAL")
    print("=" * 94)

    # Los contornos por modelo se reconstruyen aquí desde el contorno original
    # y las matrices que dejó la validación, en vez de leer los GPKG que
    # aquella escribió: así hay una sola fuente para la transformación.
    import json
    M = json.load(open(SALIDA / "matrices_humedal_pixel_a_UTM32613.json",
                       encoding="utf-8"))

    def reproyectar(poli, origen, destino):
        v = np.array(poli.exterior.coords)
        q = np.column_stack([v, np.ones(len(v))]) @ np.linalg.inv(
            np.array(origen)).T
        px = q[:, :2] / q[:, 2:3]
        q2 = np.column_stack([px, np.ones(len(px))]) @ np.array(destino).T
        return Polygon(q2[:, :2] / q2[:, 2:3])

    original = leer_poligono(VIS / "vaso_1971.gpkg")
    vaso = reproyectar(original, M["afín"], M["afín"])
    vaso_h = reproyectar(original, M["afín"], M["similitud (Helmert)"])
    caja = caja_proyectada()
    print(f"\n  contorno del vaso (modelo afín, el elegido): "
          f"{vaso.area / 1e6:.4f} km²")
    print(f"  caja estimada el 25-08, usada por la serie satelital: "
          f"{caja.area / 1e6:.4f} km²")
    print(f"  distancia entre los dos centroides: "
          f"{vaso.centroid.distance(caja.centroid):,.0f} m")

    print(f"\n{'─' * 94}\nA. QUÉ FOTOGRAMAS DEL ARCHIVO CUBREN EL VASO\n{'─' * 94}")
    cubren, df = cobertura(vaso, caja)
    cobertura_e5(vaso, caja)
    d_cent = vaso.centroid.distance(vaso_h.centroid)
    print(f"\n  El centroide del contorno se mueve {d_cent:.1f} m al cambiar de")
    print(f"  modelo, frente a huellas de 7 a 9 km de lado: la cobertura no depende")
    print(f"  del modelo.")

    print(f"\n{'─' * 94}\nB. LA SERIE SATELITAL SOBRE LA CAJA DEL VASO\n{'─' * 94}")
    serie_satelital()

    print("\n" + "=" * 94)
    print("FIN DEL SCRIPT")
    print("=" * 94)
