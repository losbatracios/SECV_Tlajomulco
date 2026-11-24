#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Huella nominal de los fotogramas sobre el AGEM 14 097
=====================================================

PROPÓSITO
---------
La línea base v8 cuenta un fotograma en un municipio cuando su centro de toma
cae dentro del área geoestadística municipal. Un fotograma cuyo centro cae
fuera puede registrar territorio del municipio si lo alcanza su huella. Este
script mide eso para los periodos posteriores a la reforma del artículo 27
(1992–1994 y 2003) y dice si la conclusión del conteo por centro cambia.

DEFINICIÓN DE HUELLA
--------------------
Huella nominal: cuadrado de 0.23 m × denominador de la escala, centrado en el
centro de toma declarado en el índice de vuelos. Es la huella que corresponde
al formato del negativo y a la escala nominal; no se mide sobre la imagen.

Se calcula en dos variantes, para ver si la orientación cambia el resultado:

  (a) alineada a los ejes de EPSG:32613;
  (b) girada al rumbo de su línea de vuelo, estimado con los centros de toma
      consecutivos de esa línea.

Ninguna de las dos es la huella real: el fotograma no está ortorrectificado,
la toma tiene inclinación y la escala nominal no es la efectiva.

CRITERIOS
---------
Un fotograma "toca" el AGEM si su huella interseca el polígono. El área
registrada de un fotograma es el área de esa intersección. El área registrada
de un conjunto es el área de la unión de sus huellas intersecada con el AGEM:
las huellas se solapan mucho (traslape nominal de vuelo), así que sumar áreas
por fotograma contaría dos veces el mismo terreno.

Una línea de vuelo se clasifica por su huella en tres clases: la huella toca
el AGEM; no lo toca pero queda a 2 km o menos; o queda a más de 2 km. Se
compara con la clasificación por traza de centros de la v8.

Áreas y distancias en EPSG:32613 (WGS 84 / UTM 13N).

SALIDAS  (04_Análisis/Visualizaciones/04_linea_base/)
  huella_fotogramas_post1992.csv    una fila por fotograma de 1992, 1993,
                                    1994 y 2003, con su huella y su área
  huella_resumen_anio_escala.csv    por año, escala y vuelo
  huella_lineas_1992_1994.csv       las líneas de vuelo reclasificadas

Requiere: numpy, pandas, pyproj, shapely, pyshp
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.affinity import rotate
from shapely.geometry import Point, box
from shapely.ops import unary_union

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
from poligonos_v12 import cargar_oficial, NOMBRES                # noqa: E402
from linea_base_v8 import cargar, SALIDA, CASO                   # noqa: E402

LADO_NEGATIVO_M = 0.23
ANIOS = [1992, 1993, 1994, 2003]
UMBRAL_CERCA_M = 2000.0


def huella(e, n, den, rumbo=None):
    """Cuadrado de 0.23 m × denominador centrado en (e, n), girado al rumbo."""
    semi = LADO_NEGATIVO_M / 2 * den
    cuadro = box(e - semi, n - semi, e + semi, n + semi)
    if rumbo is None:
        return cuadro
    return rotate(cuadro, rumbo, origin=(e, n), use_radians=False)


def rumbos_por_linea(df):
    """Rumbo de cada línea de vuelo, en grados, desde sus centros de toma."""
    r = {}
    for llave, g in df.groupby(["Vuelo", "Bote", "Línea", "Escala"]):
        g = g.sort_values("Número de foto")
        if len(g) < 2:
            r[llave] = 0.0
            continue
        de = g["E"].iloc[-1] - g["E"].iloc[0]
        dn = g["N"].iloc[-1] - g["N"].iloc[0]
        r[llave] = float(np.degrees(np.arctan2(dn, de)))
    return r


# ---------------------------------------------------------------------------
# 1. Huella por fotograma
# ---------------------------------------------------------------------------

def por_fotograma(df, agem):
    caso = agem[CASO]
    km2_agem = caso.area / 1e6
    print(f"\n  AGEM 14 097: {km2_agem:,.2f} km²")

    d = df[df["Año"].isin(ANIOS)].copy()
    print(f"  fotogramas del corpus en {ANIOS}: {len(d):,}")
    sin_escala = int(d["den"].isna().sum())
    if sin_escala:
        print(f"  AVISO: {sin_escala} sin denominador de escala legible; se omiten")
        d = d[d["den"].notna()]

    rumbos = rumbos_por_linea(d)

    filas = []
    for _, r in d.iterrows():
        llave = (r["Vuelo"], r["Bote"], r["Línea"], r["Escala"])
        h_ejes = huella(r["E"], r["N"], r["den"])
        h_giro = huella(r["E"], r["N"], r["den"], rumbos[llave])
        inter_e = h_ejes.intersection(caso)
        inter_g = h_giro.intersection(caso)
        p = Point(r["E"], r["N"])
        filas.append(dict(
            anio=int(r["Año"]), vuelo=r["Vuelo"], bote=r["Bote"],
            linea=r["Línea"], foto=r["Número de foto"], escala=r["Escala"],
            den=int(r["den"]), lado_huella_m=round(LADO_NEGATIVO_M * r["den"]),
            lat=round(r["Latitud"], 5), lon=round(r["Longitud"], 5),
            mun_centro=NOMBRES.get(r["mun"], "fuera de los seis"),
            centro_en_el_agem=bool(caso.contains(p)),
            dist_centro_agem_m=round(0.0 if caso.contains(p) else p.distance(caso)),
            rumbo_linea_grados=round(rumbos[llave], 1),
            toca_ejes=not inter_e.is_empty,
            area_ejes_km2=round(inter_e.area / 1e6, 2),
            toca_giro=not inter_g.is_empty,
            area_giro_km2=round(inter_g.area / 1e6, 2),
            dist_huella_agem_m=round(0.0 if not inter_e.is_empty
                                     else h_ejes.distance(caso)),
        ))
    t = pd.DataFrame(filas).sort_values(["anio", "vuelo", "bote", "linea", "foto"])
    t.to_csv(SALIDA / "huella_fotogramas_post1992.csv", index=False,
             encoding="utf-8")

    print(f"\n  fotogramas cuya huella alineada a los ejes toca el AGEM: "
          f"{int(t['toca_ejes'].sum())}")
    print(f"  fotogramas cuya huella girada al rumbo toca el AGEM:      "
          f"{int(t['toca_giro'].sum())}")
    print(f"  fotogramas con el centro dentro del AGEM:                 "
          f"{int(t['centro_en_el_agem'].sum())}")
    desacuerdo = t[t["toca_ejes"] != t["toca_giro"]]
    print(f"  fotogramas en que las dos variantes discrepan: {len(desacuerdo)}")
    if len(desacuerdo):
        print(desacuerdo[["anio", "vuelo", "bote", "linea", "foto", "escala",
                          "toca_ejes", "toca_giro"]].to_string(index=False))
    return t, d, caso, km2_agem, rumbos


# ---------------------------------------------------------------------------
# 2. Resumen por año, escala y vuelo
# ---------------------------------------------------------------------------

def area_union(sub, d, caso, rumbos, girada):
    """Área del AGEM cubierta por la unión de las huellas de sub."""
    polis = []
    for _, r in sub.iterrows():
        llave = (r["vuelo"], r["bote"], r["linea"], r["escala"])
        fila = d[(d["Vuelo"] == r["vuelo"]) & (d["Bote"] == r["bote"])
                 & (d["Línea"] == r["linea"]) & (d["Escala"] == r["escala"])
                 & (d["Número de foto"] == r["foto"])]
        e, n = fila["E"].iloc[0], fila["N"].iloc[0]
        polis.append(huella(e, n, r["den"], rumbos[llave] if girada else None))
    if not polis:
        return 0.0
    return unary_union(polis).intersection(caso).area / 1e6


def resumen(t, d, caso, km2_agem, rumbos):
    print(f"\n  área registrada = unión de las huellas que tocan el AGEM, "
          f"intersecada con el AGEM")

    filas = []
    for (anio, escala, vuelo), g in t.groupby(["anio", "escala", "vuelo"]):
        tocan = g[g["toca_ejes"]]
        a_ejes = area_union(tocan, d, caso, rumbos, False)
        a_giro = area_union(g[g["toca_giro"]], d, caso, rumbos, True)
        filas.append(dict(
            anio=anio, escala=escala, vuelo=vuelo, fotogramas=len(g),
            centro_dentro=int(g["centro_en_el_agem"].sum()),
            huella_toca=len(tocan),
            area_union_km2=round(a_ejes, 2),
            pct_agem=round(a_ejes / km2_agem * 100, 2),
            area_union_girada_km2=round(a_giro, 2),
            max_area_fotograma_km2=round(g["area_ejes_km2"].max(), 2),
            dist_min_huella_m=int(g["dist_huella_agem_m"].min()),
        ))
    r = pd.DataFrame(filas).sort_values(["anio", "escala", "vuelo"])
    r.to_csv(SALIDA / "huella_resumen_anio_escala.csv", index=False,
             encoding="utf-8")

    print()
    print(r.to_string(index=False))

    print("\n  por año:")
    for anio in ANIOS:
        g = t[t["anio"] == anio]
        if not len(g):
            print(f"    {anio}: sin fotogramas en el corpus")
            continue
        a = area_union(g[g["toca_ejes"]], d, caso, rumbos, False)
        print(f"    {anio}: {len(g):3d} fotogramas, "
              f"{int(g['toca_ejes'].sum()):3d} con huella que toca el AGEM, "
              f"{a:6.2f} km² = {a / km2_agem * 100:5.2f} % del AGEM")

    a_total = area_union(t[t["toca_ejes"]], d, caso, rumbos, False)
    a_total_g = area_union(t[t["toca_giro"]], d, caso, rumbos, True)
    print(f"\n  unión de todas las huellas de {ANIOS[0]} en adelante: "
          f"{a_total:.2f} km² = {a_total / km2_agem * 100:.2f} % del AGEM")
    print(f"  la misma unión con las huellas giradas al rumbo: "
          f"{a_total_g:.2f} km² = {a_total_g / km2_agem * 100:.2f} %")
    return r, a_total, a_total_g


# ---------------------------------------------------------------------------
# 3. ¿Cambia la clasificación de las líneas de vuelo?
# ---------------------------------------------------------------------------

def lineas(t):
    v8 = pd.read_csv(SALIDA / "lineas_vuelo_1992_1994.csv")
    print(f"\n  líneas de vuelo de 1992–1994 según la v8: {len(v8)}")
    print(f"    la traza de centros cruza el AGEM: "
          f"{int(v8['traza_cruza_el_agem'].sum())}")
    print(f"    no cruza y queda a {UMBRAL_CERCA_M / 1000:.0f} km o menos: "
          f"{int(((~v8['traza_cruza_el_agem']) & (v8['dist_min_m'] <= UMBRAL_CERCA_M)).sum())}")
    print(f"    queda a más de {UMBRAL_CERCA_M / 1000:.0f} km: "
          f"{int((v8['dist_min_m'] > UMBRAL_CERCA_M).sum())}")

    tl = t[t["anio"].between(1992, 1994)]
    filas = []
    for (vuelo, bote, linea, escala), g in tl.groupby(
            ["vuelo", "bote", "linea", "escala"]):
        toca = bool(g["toca_ejes"].any())
        d_min = int(g["dist_huella_agem_m"].min())
        clase = ("la huella toca el AGEM" if toca
                 else ("a 2 km o menos" if d_min <= UMBRAL_CERCA_M
                       else "a más de 2 km"))
        filas.append(dict(vuelo=vuelo, bote=bote, linea=linea, escala=escala,
                          anio_huella=int(g["anio"].iloc[0]),
                          n_huella=len(g),
                          fotos_que_tocan=int(g["toca_ejes"].sum()),
                          dist_min_huella_m=d_min,
                          clase_por_huella=clase))
    h = pd.DataFrame(filas)

    # Las dos tablas vienen de rutas distintas y los tipos de bote y línea no
    # coinciden; se comparan como texto para que el cruce sea uno a uno.
    for col in ("vuelo", "bote", "linea", "escala"):
        v8[col] = v8[col].astype(str)
        h[col] = h[col].astype(str)

    fusion = v8.merge(h, on=["vuelo", "bote", "linea", "escala"],
                      how="outer", indicator=True)
    if (fusion["_merge"] != "both").any():
        print("  AVISO: hay líneas que no cruzan entre las dos tablas")
        print(fusion.loc[fusion["_merge"] != "both",
                         ["vuelo", "bote", "linea", "escala", "_merge"]]
              .to_string(index=False))
    fusion = fusion.drop(columns="_merge")
    fusion["clase_por_centros"] = np.where(
        fusion["traza_cruza_el_agem"], "la traza cruza el AGEM",
        np.where(fusion["dist_min_m"] <= UMBRAL_CERCA_M,
                 "a 2 km o menos", "a más de 2 km"))
    fusion.to_csv(SALIDA / "huella_lineas_1992_1994.csv", index=False,
                  encoding="utf-8")

    print(f"\n  líneas de vuelo clasificadas por huella: {len(h)}")
    for k, v in h["clase_por_huella"].value_counts().items():
        print(f"    {k:26} {v:3d}")

    print("\n  comparación de las dos clasificaciones:")
    print(pd.crosstab(fusion["clase_por_centros"],
                      fusion["clase_por_huella"]).to_string())
    cambian = fusion[fusion["clase_por_centros"] != fusion["clase_por_huella"]]
    print(f"\n  líneas que cambian de clase: {len(cambian)} de {len(fusion)}")
    if len(cambian):
        print(cambian[["anio", "vuelo", "bote", "linea", "escala", "n",
                       "fotos_que_tocan", "dist_min_m", "dist_min_huella_m",
                       "clase_por_centros", "clase_por_huella"]]
              .to_string(index=False))
    return fusion


if __name__ == "__main__":
    print("=" * 100)
    print("HUELLA NOMINAL SOBRE EL AGEM 14 097 — 1992, 1993, 1994 Y 2003")
    print("=" * 100)

    agem = cargar_oficial()
    df, _ = cargar(agem)

    print("\n" + "─" * 100)
    print("1. HUELLA POR FOTOGRAMA")
    print("─" * 100)
    t, d, caso, km2, rumbos = por_fotograma(df, agem)

    print("\n" + "─" * 100)
    print("2. ÁREA REGISTRADA POR AÑO, ESCALA Y VUELO")
    print("─" * 100)
    resumen(t, d, caso, km2, rumbos)

    print("\n" + "─" * 100)
    print("3. ¿CAMBIA LA CLASIFICACIÓN DE LAS LÍNEAS DE VUELO?")
    print("─" * 100)
    lineas(t)

    print("\n" + "=" * 100)
    print("FIN DEL SCRIPT")
    print("=" * 100)
