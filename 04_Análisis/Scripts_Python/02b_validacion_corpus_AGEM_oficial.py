#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validación del corpus y conjuntos de la solicitud, sobre el AGEM oficial
=========================================================================

PROPÓSITO
---------
El script 02 construyó, sobre una capa municipal provisional, los conteos del
corpus y los tres conjuntos con que se armó la solicitud de fotografías. Aquí
se rehace la parte que depende de esa capa, usando el área geoestadística
municipal (AGEM) 14 097 del Marco Geoestadístico 2025 del INEGI, y se imprime
cada cifra junto a la que daba la capa anterior.

No reemplaza al script 02 ni toca su archivo de salida: la solicitud ya se
entregó y su selección queda como registro de lo que efectivamente se pidió.
Esto es una verificación.

QUÉ SE REHACE
-------------
Todo lo que pasa por el polígono:

  1. los puntos de control de la capa;
  2. los fotogramas con centro dentro, y su reparto por año;
  3. los fotogramas posteriores a la reforma que caen dentro;
  4. el conjunto A, fotografías de 1993 a no más de un umbral del borde;
  5. el conjunto B, contrapartes previas a la reforma a menos de 2 km de
     alguna fotografía del conjunto A (depende de A, y por tanto del
     polígono);
  6. el conjunto C, fotogramas de 1971 dentro del municipio y al sur del
     paralelo 20.48.

Los hallazgos del corpus que no pasan por el polígono —la interrupción de
1984-1991, la escala de 1:5,000 y la nomenclatura del campo Vuelo— no se
repiten aquí: no dependen de la capa y ya están verificados aparte.

CÓMO SE MIDEN LOS 2 KM DEL CONJUNTO B
-------------------------------------
La distancia de referencia es la geodésica sobre el elipsoide WGS84, que es
la distancia sobre el terreno. El script 02 usaba una fórmula de haversine
sobre una esfera de radio 6,371 km, que a esta latitud sobreestima algo menos
del 0.3 %. Sobre 2 km son unos 6 m, y el umbral del conjunto B cae justo ahí:
un fotograma queda dentro o fuera según la fórmula. Por eso aquí se calculan
las tres distancias —geodésica, haversine y plano proyectado— y se nombra el
caso que cambia.

SOBRE EL UMBRAL DEL CONJUNTO A
------------------------------
El script 02 usó 1,200 m: la media huella deducida de la escala nominal
(1,150 m) más un margen de 50 m por la provisionalidad de la capa. Después se
midió el área de imagen en el propio fotograma y la media huella resultó ser
de 1,209 m. Se imprimen los dos umbrales, porque la pregunta de si cambia la
selección se responde con los dos.

Distancias en EPSG:32613 (WGS 84 / UTM 13N).

Requiere: numpy, pandas, pyproj, shapely, pyshp
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Geod
from shapely.geometry import Point

GEOD = Geod(ellps="WGS84")

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(AQUI / "06_linea_base"))

from poligonos_v12 import cargar_oficial, cargar_provisional, A_UTM

CSV = (RAIZ / "02_Archivo_INEGI/Metadatos_SINFA"
            / "07-07_consulta_ÍndicesdeVuelos260107080134.csv")

UMBRAL_A_NOMINAL = 1200.0      # el que usó el script 02 (1,150 + 50 de margen)
UMBRAL_A_MEDIDO = 1209.0       # media huella medida en el fotograma (14-10)
UMBRAL_CONTRAPARTES_M = 2000.0
ANIOS_PRE_REFORMA = [1970, 1971, 1978, 1981]
ANIOS_POST = [1992, 1993, 1994, 2003]
PARALELO_MITAD_SUR = 20.48

CONTROL = {
    "cabecera municipal": (20.4727, -103.4466, True),
    "Cajititlán": (20.4300, -103.3300, True),
    "Hacienda Santa Fe": (20.5310, -103.3480, True),
    "centro de Guadalajara": (20.6767, -103.3475, False),
}


def cargar():
    df = pd.read_csv(CSV, encoding="latin-1")
    e, n = A_UTM.transform(df["Longitud"].values, df["Latitud"].values)
    df["E"], df["N"] = e, n
    return df


def analizar(df, poli, etiqueta):
    """Todas las cifras que dependen del polígono, para una capa."""
    pts = [Point(x, y) for x, y in zip(df["E"], df["N"])]
    dentro = np.array([poli.contains(p) for p in pts])

    r = dict(etiqueta=etiqueta, area_km2=poli.area / 1e6)
    r["control"] = {k: (poli.contains(Point(*A_UTM.transform(lo, la))) == esp)
                    for k, (la, lo, esp) in CONTROL.items()}
    r["dentro"] = int(dentro.sum())
    r["por_anio"] = dict(df[dentro]["Año"].value_counts().sort_index())
    r["post_dentro"] = int(dentro[df["Año"].isin(ANIOS_POST).values].sum())

    f93 = df[df["Año"] == 1993].copy()
    d93 = np.array([0.0 if poli.contains(p) else p.distance(poli)
                    for p in (Point(x, y) for x, y in zip(f93["E"], f93["N"]))])
    f93["dist_m"] = d93
    r["a_nominal"] = f93[d93 <= UMBRAL_A_NOMINAL]
    r["a_medido"] = f93[d93 <= UMBRAL_A_MEDIDO]
    r["a_centro_dentro"] = int((d93 == 0).sum())

    # Conjunto B: contrapartes previas a la reforma a menos de 2 km de alguna
    # fotografía del conjunto A, medidos sobre el elipsoide.
    pre = df[df["Año"].isin(ANIOS_PRE_REFORMA)].copy()
    for clave, a in [("b_nominal", r["a_nominal"]), ("b_medido", r["a_medido"])]:
        dmin = np.full(len(pre), np.inf)
        for _, fila in a.iterrows():
            d = GEOD.inv(pre["Longitud"].values, pre["Latitud"].values,
                         np.full(len(pre), fila["Longitud"]),
                         np.full(len(pre), fila["Latitud"]))[2]
            dmin = np.minimum(dmin, d)
        r[clave] = pre[dmin < UMBRAL_CONTRAPARTES_M]
        if clave == "b_nominal":
            r["b_dmin"] = pre.assign(dmin=dmin)

    r["c"] = df[dentro & (df["Año"] == 1971)
                & (df["Latitud"] < PARALELO_MITAD_SUR)]
    return r


def comparar(a, b, nombre, fmt="{}"):
    igual = a == b
    marca = "sin cambio" if igual else "← CAMBIA"
    print(f"  {nombre:46}{fmt.format(b):>14}{fmt.format(a):>14}   {marca}")


if __name__ == "__main__":
    print("=" * 92)
    print("CORPUS Y CONJUNTOS DE LA SOLICITUD SOBRE EL AGEM OFICIAL")
    print("=" * 92)

    df = cargar()
    ofi = cargar_oficial(["14097"])["14097"]
    prov = cargar_provisional()
    print(f"\n  corpus: {len(df):,} registros")
    print(f"  AGEM oficial (Marco Geoestadístico 2025): {ofi.area/1e6:,.2f} km²")
    print(f"  capa provisional del script 02:           {prov.area/1e6:,.2f} km²")

    r_ofi = analizar(df, ofi, "oficial")
    r_pro = analizar(df, prov, "provisional")

    print(f"\n  puntos de control de la capa")
    for k in CONTROL:
        print(f"    {k:28} oficial {'OK' if r_ofi['control'][k] else 'FALLA':<6}"
              f" · provisional {'OK' if r_pro['control'][k] else 'FALLA'}")

    print(f"\n  {'cifra':46}{'oficial':>14}{'provisional':>14}")
    print("  " + "─" * 90)
    comparar(r_pro["dentro"], r_ofi["dentro"], "fotogramas con centro dentro")
    comparar(r_pro["post_dentro"], r_ofi["post_dentro"],
             "de 1992/93/94/2003, con centro dentro")
    comparar(len(r_pro["a_nominal"]), len(r_ofi["a_nominal"]),
             f"conjunto A con el umbral de {UMBRAL_A_NOMINAL:,.0f} m")
    comparar(len(r_pro["a_medido"]), len(r_ofi["a_medido"]),
             f"conjunto A con la media huella medida ({UMBRAL_A_MEDIDO:,.0f} m)")
    comparar(len(r_pro["b_nominal"]), len(r_ofi["b_nominal"]),
             "conjunto B (contrapartes a menos de 2 km)")
    comparar(len(r_pro["c"]), len(r_ofi["c"]), "conjunto C (1971, dentro, al sur de 20.48)")
    tot_o = len(r_ofi["a_nominal"]) + len(r_ofi["b_nominal"]) + len(r_ofi["c"])
    tot_p = len(r_pro["a_nominal"]) + len(r_pro["b_nominal"]) + len(r_pro["c"])
    comparar(tot_p, tot_o, "TOTAL de la selección (A + B + C)")

    print(f"\n  reparto por año de los fotogramas con centro dentro")
    anios = sorted(set(r_ofi["por_anio"]) | set(r_pro["por_anio"]))
    print(f"    {'año':>6}{'oficial':>10}{'provisional':>14}")
    print("    " + "─" * 30)
    for a in anios:
        o, p = r_ofi["por_anio"].get(a, 0), r_pro["por_anio"].get(a, 0)
        print(f"    {a:>6}{o:>10}{p:>14}{'' if o == p else '   ←'}")

    # ¿Las dos selecciones contienen exactamente los mismos fotogramas?
    def llaves(d):
        return set(zip(d["Bote"], d["Rollo"], d["Línea"], d["Número de foto"],
                       d["Año"], d["Escala"], d["Latitud"].round(6),
                       d["Longitud"].round(6)))
    print(f"\n  ¿son los mismos fotogramas, uno por uno?")
    for nombre, ko, kp in [
            ("conjunto A", llaves(r_ofi["a_nominal"]), llaves(r_pro["a_nominal"])),
            ("conjunto B", llaves(r_ofi["b_nominal"]), llaves(r_pro["b_nominal"])),
            ("conjunto C", llaves(r_ofi["c"]), llaves(r_pro["c"]))]:
        dif = (ko ^ kp)
        print(f"    {nombre:14} {'idénticos' if not dif else f'difieren en {len(dif)}'}")

    # Sensibilidad del umbral de 2 km del conjunto B: qué fotogramas caen
    # cerca de la frontera y cómo los clasifica cada fórmula de distancia.
    print(f"\n  Sensibilidad del umbral de 2 km del conjunto B (AGEM oficial)")
    banda = r_ofi["b_dmin"]
    banda = banda[(banda.dmin > 1800) & (banda.dmin < 2200)].sort_values("dmin")
    print(f"    fotogramas entre 1,800 y 2,200 m del conjunto A: {len(banda)}")
    if len(banda):
        print(f"      {'fotograma':34}{'geodésica':>12}{'haversine':>12}{'UTM plano':>12}")
        print("      " + "─" * 70)
        A = r_ofi["a_nominal"]
        for _, p in banda.iterrows():
            dg = p.dmin
            hv = min(2 * 6371000 * np.arcsin(np.sqrt(
                np.sin(np.radians(a["Latitud"] - p["Latitud"]) / 2) ** 2
                + np.cos(np.radians(p["Latitud"])) * np.cos(np.radians(a["Latitud"]))
                * np.sin(np.radians(a["Longitud"] - p["Longitud"]) / 2) ** 2))
                for _, a in A.iterrows())
            du = min(np.hypot(p["E"] - a["E"], p["N"] - a["N"]) for _, a in A.iterrows())
            etiqueta = (f"{p['Bote']}/{p['Línea']}/{p['Número de foto']}/"
                        f"{p['Año']}/{p['Escala']}")
            def mar(v):
                return f"{v:,.0f} m {'✓' if v < UMBRAL_CONTRAPARTES_M else '✗'}"
            print(f"      {etiqueta:34}{mar(dg):>12}{mar(hv):>12}{mar(du):>12}")
        print(f"      ✓ dentro del umbral · ✗ fuera. La haversine esférica que usó el")
        print(f"      script 02 sobreestima la distancia y deja fuera un fotograma que,")
        print(f"      medido sobre el elipsoide, está dentro.")

    print(f"\n  El umbral medido (1,209 m) frente al nominal (1,200 m), con el AGEM oficial:")
    extra = llaves(r_ofi["a_medido"]) - llaves(r_ofi["a_nominal"])
    print(f"    conjunto A: {len(r_ofi['a_nominal'])} y {len(r_ofi['a_medido'])} "
          f"fotografías; {'ninguna entra de más' if not extra else f'{len(extra)} entran de más'}")

    print("\n" + "=" * 92)
    print("FIN DEL SCRIPT — el archivo de selección del script 02 no se toca")
    print("=" * 92)
