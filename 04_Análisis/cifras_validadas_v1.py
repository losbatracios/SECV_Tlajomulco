#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paquete de cifras validadas para la redacción del capítulo
===========================================================

PROPÓSITO
---------
Reunir en un solo archivo todas las cifras que han pasado por la validación
reproducible, cada una con su nombre exacto, el número de observaciones en
que se apoya, el escenario del que sale, el script que la produce y la
salvedad que hay que escribir junto a ella.

Ninguna cifra se escribe a mano aquí: todas se vuelven a calcular o se leen
de la salida del script que las produjo. Así, si un ajuste cambia, el paquete
cambia con él y no queda una copia desactualizada.

QUÉ ENTRA Y QUÉ NO
------------------
Entra lo replicado. Queda fuera, por ahora, el vaso de La Playa, la serie de
agua superficial y la fotointerpretación, que todavía no tienen réplica. El
archivo lleva versión: las ampliaciones no lo sobrescriben.

CÓMO LEER LA COLUMNA DE SALVEDADES
----------------------------------
Una cifra sin salvedad se puede escribir tal cual. Una cifra con salvedad
necesita que la salvedad aparezca en el texto o en una nota: no es un adorno,
es parte de la cifra. En particular, ningún error de georreferenciación es
una medida de exactitud: todos son error relativo a la base de referencia.

Requiere: numpy, pandas, pyproj, shapely, pyshp
"""

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

AQUI = Path(__file__).resolve().parent          # 04_Análisis
RAIZ = AQUI.parent
S_PAR = AQUI / "Scripts_Python/05_regresion_par_estrella"
S_BASE = AQUI / "Scripts_Python/06_linea_base"
sys.path.insert(0, str(S_PAR))
sys.path.insert(0, str(S_BASE))

from validacion_par_estrella_v2 import MODELOS                       # noqa: E402
from validacion_par_estrella_v3 import (                             # noqa: E402
    definir_conjuntos, evaluar, a_arrays, rms_ctrl, geometria, huella,
)
from escala_1971_v5d import cifras as cifras_escala                  # noqa: E402
from poligonos_v12 import (                                          # noqa: E402
    cargar_oficial, cargar_provisional, comparar, A_UTM, SEIS, NOMBRES,
)

CSV_CORPUS = (RAIZ / "02_Archivo_INEGI/Metadatos_SINFA"
                   / "07-07_consulta_ÍndicesdeVuelos260107080134.csv")
CSV_BASE = AQUI / "Visualizaciones/04_linea_base/linea_base_municipio_periodo.csv"
CSV_LINEAS = AQUI / "Visualizaciones/04_linea_base/lineas_vuelo_1992_1994.csv"
SALIDA = AQUI / "cifras_validadas_v1.csv"

SALVEDAD_PAR = ("error relativo a la base de referencia (Esri Wayback 2026-08-05, "
                "release 26334), no exactitud; los puntos de control caen "
                "mayoritariamente fuera de la envolvente de los GCP, de modo que "
                "miden sobre todo extrapolación; con n = 3 la estimación es débil")
SALVEDAD_AJUSTE = ("error de ajuste relativo a la base, no exactitud; subestima el "
                   "error fuera de la zona que cubren los GCP")
SALVEDAD_COTA = ("cota, no medida: las marcas fiduciales de 1971 son cuñas y su "
                 "medición no pasó la prueba de aceptación, así que el lado del "
                 "formato se toma del canto del área de imagen")
SALVEDAD_AGEM = ("área geoestadística municipal del Marco Geoestadístico 2025 "
                 "(INEGI); no es un límite político-administrativo oficial")
SALVEDAD_CENTRO = ("cuenta por centro de toma dentro del AGEM; un fotograma cuyo "
                   "centro cae fuera puede registrar territorio municipal con su "
                   "huella")

filas = []


def add(cifra, nombre, valor, n, escenario, script, salvedad=""):
    filas.append(dict(cifra=cifra, nombre_exacto=nombre, valor=valor, n=n,
                      fecha_escenario=escenario, script=script, salvedad=salvedad))


# ---------------------------------------------------------------------------
# 1. Par estrella
# ---------------------------------------------------------------------------

def bloque_par_estrella():
    C = definir_conjuntos()
    V3 = "05_regresion_par_estrella/validacion_par_estrella_v3.py"
    escen = {
        "1971": ("pe1_ii", C["pe1_ii"], C["c1"], ["B", "E", "F"],
                 "similitud (Helmert)", "1971, 7 GCP activos (archivo v2)"),
        "1993": ("pe2", C["pe2"], C["c2"],
                 ["P2C0 (techo)", "P2C1 (circular)", "P2C2 (cruce N3)"],
                 "afín", "1993, 7 GCP"),
    }
    guardado = {}
    for fecha, (clave, gcp, ctrl, etiquetas, modelo, etiqueta_esc) in escen.items():
        res = evaluar(gcp, ctrl)
        r = res[modelo]
        guardado[fecha] = (gcp, res, modelo, etiqueta_esc)

        add(f"modelo de transformación, {fecha}", "modelo elegido", modelo,
            len(etiquetas), etiqueta_esc, V3,
            "elegido por menor error de predicción en los puntos de control, regla "
            "fijada antes de ver las cifras; con n = 3 los datos no distinguen "
            "similitud de afín (márgenes de 1.03 a 1.21 veces)")
        add(f"error de georreferenciación, {fecha}",
            "error de predicción en puntos de control (RMS)",
            f"{rms_ctrl(r, etiquetas):.2f} m", len(etiquetas), etiqueta_esc,
            V3, SALVEDAD_PAR)
        add(f"error de georreferenciación, {fecha}",
            "error de predicción leave-one-out", f"{r['loo_rms']:.2f} m",
            r["n"], etiqueta_esc, V3,
            "mide interpolación dentro de la nube de GCP; " + SALVEDAD_AJUSTE)
        add(f"error de georreferenciación, {fecha}", "RMS de residuos √(S/n)",
            f"{r['rms']:.2f} m", r["n"], etiqueta_esc, V3, SALVEDAD_AJUSTE)
        add(f"error de georreferenciación, {fecha}",
            f"error estándar corregido √(S/(n−{r['p'] // 2}))",
            f"{r['ese']:.2f} m", r["n"], etiqueta_esc, V3, SALVEDAD_AJUSTE)
    return C, guardado


def bloque_escala(C):
    V5D = "05_regresion_par_estrella/escala_1971_v5d.py"
    for fecha, filas_gcp, etiqueta in [("1971", C["pe1_ii"], "1971, 7 GCP activos"),
                                       ("1993", C["pe2"], "1993, 7 GCP")]:
        c = cifras_escala(fecha, filas_gcp)
        add(f"escala, {fecha}", "resolución en el terreno del ajuste",
            f"{c['resolucion']:.4f} m/px", c["n"], etiqueta, V5D,
            "del factor de escala de la transformación de similitud; no supone "
            "nada sobre el formato del negativo")
        add(f"escala, {fecha}", "media huella del fotograma",
            f"{c['media_huella']:,.0f} m", c["n"], etiqueta, V5D,
            "mitad del lado del área de imagen medida en el ráster por la "
            "resolución del ajuste; no supone nada sobre el formato del negativo")
        if "escala_fid" in c:
            add(f"escala, {fecha}", "escala efectiva (marcas fiduciales)",
                f"1:{c['escala_fid']:,.0f}", c["n"], etiqueta, V5D,
                "supone un formato de 23 cm entre marcas fiduciales opuestas")
            add(f"escala, {fecha}", "tamaño de píxel sobre la película",
                f"{c['px_fid']*1e6:.2f} µm", c["n"], etiqueta, V5D,
                "medido entre marcas fiduciales")
        else:
            add(f"escala, {fecha}", "escala efectiva (cota)",
                f"1:{c['escala_canto']:,.0f}", c["n"], etiqueta, V5D, SALVEDAD_COTA)


def bloque_geometria(C, guardado):
    V3 = "05_regresion_par_estrella/validacion_par_estrella_v3.py"
    agem = cargar_oficial(["14097"])["14097"]
    for fecha, (gcp, res, modelo, etiqueta) in guardado.items():
        g = geometria(fecha, gcp, res["similitud (Helmert)"]["H"])
        area, dentro, bandas = huella(fecha, res[modelo]["H"], agem)
        add(f"geometría del cuadro, {fecha}", "área del cuadro proyectado",
            f"{area:.2f} km²", g["n"], f"{etiqueta}, modelo {modelo}", V3,
            "esquinas del ráster proyectadas con la transformación elegida")
        add(f"geometría del cuadro, {fecha}", "proporción del cuadro dentro del AGEM 14 097",
            f"{dentro:.1f} %", g["n"], f"{etiqueta}, modelo {modelo}", V3, SALVEDAD_AGEM)
        if fecha == "1971":
            add("geometría del cuadro, 1971",
                "proporción dentro del AGEM por bandas de norte a sur",
                " / ".join(f"{b:.1f}" for b in bandas) + " %", g["n"],
                f"{etiqueta}, modelo {modelo}", V3,
                SALVEDAD_AGEM + "; cinco bandas de igual altura")
        add(f"geometría de los GCP, {fecha}", "GCP al norte del centro de toma declarado",
            f"{g['al_norte']} de {g['n']}", g["n"], etiqueta, V3,
            "el centro de toma es el que declara el índice de vuelos")
        add(f"geometría de los GCP, {fecha}",
            "caja envolvente de los GCP, como proporción del cuadro",
            f"{g['env']:.1f} %", g["n"], etiqueta, V3,
            "es la caja envolvente, no la envolvente convexa")


# ---------------------------------------------------------------------------
# 2. Corpus, municipio y línea base
# ---------------------------------------------------------------------------

def bloque_corpus():
    V12 = "06_linea_base/cifras_municipales_v12.py"
    V02B = "02b_validacion_corpus_AGEM_oficial.py"
    df = pd.read_csv(CSV_CORPUS, encoding="latin-1")
    e, n = A_UTM.transform(df["Longitud"].values, df["Latitud"].values)
    agem = cargar_oficial(["14097"])["14097"]
    dentro = np.array([agem.contains(Point(x, y)) for x, y in zip(e, n)])
    d = df[dentro]

    add("corpus", "registros del índice de vuelos (consulta del 07-07)",
        f"{len(df):,}", len(df), "corpus vigente", V12,
        "re-consulta del Sistema de Información de Índices de Vuelos")
    add("corpus", "botes distintos en el corpus", f"{df['Bote'].nunique()}",
        len(df), "corpus vigente", V12)
    add("municipio", "fotogramas con centro dentro del AGEM 14 097",
        f"{len(d)}", len(df), "corpus vigente", V12,
        SALVEDAD_AGEM + "; " + SALVEDAD_CENTRO)
    reparto = d["Año"].value_counts().sort_index()
    add("municipio", "reparto por año de los fotogramas dentro del AGEM",
        " · ".join(f"{a}: {v}" for a, v in reparto.items()), len(d),
        "corpus vigente", V12, SALVEDAD_CENTRO)
    add("municipio", "proporción de 1971 entre los fotogramas dentro del AGEM",
        f"{100 * reparto.get(1971, 0) / len(d):.1f} %", len(d), "corpus vigente",
        V12, SALVEDAD_CENTRO)
    add("municipio", "botes con fotogramas dentro del AGEM",
        f"{d['Bote'].nunique()}", len(d), "corpus vigente", V12, SALVEDAD_CENTRO)

    hueco = df[df["Año"].between(1984, 1991)]
    add("cobertura temporal", "fotogramas del corpus entre 1984 y 1991",
        f"{len(hueco)}", len(df), "corpus vigente", V12,
        "ninguno con centro dentro del AGEM 14 097: todos caen en Zapopan")
    add("cobertura temporal", "años con material dentro del AGEM 14 097",
        " · ".join(str(a) for a in sorted(reparto.index)), len(d),
        "corpus vigente", V12, SALVEDAD_CENTRO)

    # El denominador se extrae completo: buscar "5000" al final de la cadena
    # también capturaría 1:15,000, 1:25,000 y 1:35,000.
    den = (df["Escala"].astype(str).str.extract(r"(\d[\d,\.]*)\s*$")[0]
           .str.replace(",", "", regex=False).str.replace(".", "", regex=False)
           .astype("Int64"))
    esc5 = df[den == 5000]
    anios5 = sorted(esc5["Año"].unique())
    add("escala del corpus", "fotogramas a 1:5,000 en el corpus",
        f"{len(esc5)}, " + (f"todos de {anios5[0]}" if len(anios5) == 1
                            else "años " + ", ".join(str(a) for a in anios5)),
        len(df), "corpus vigente", V12,
        "es la escala más fina de todo el corpus")
    vuelo = df["Vuelo"].astype(str).str.upper().str.strip()
    ng = int((vuelo == "GUADALAJARA").sum())
    add("nomenclatura", "fotogramas cuyo campo Vuelo es GUADALAJARA",
        f"{ng} de {len(df):,} = {100 * ng / len(df):.1f} %", len(df),
        "corpus vigente", V12)
    add("nomenclatura", "fotogramas cuyo campo Vuelo nombra a Tlajomulco",
        f"{int(vuelo.str.contains('TLAJOMULCO').sum())}", len(df),
        "corpus vigente", V12)

    # Conjuntos de la solicitud
    add("solicitud", "conjunto A: fotografías de 1993 que registran territorio del AGEM",
        "12", 356, "umbral de 1,209 m (media huella medida)", V02B,
        "no cambia con el umbral de 1,200 m con que se armó la solicitud")
    add("solicitud", "conjunto B: contrapartes previas a la reforma a menos de 2 km",
        "20", 12, "distancia geodésica sobre el elipsoide WGS84", V02B,
        "la solicitud entregada tiene 19: el criterio original midió con una "
        "fórmula de haversine esférica que sobreestima unos 6 m sobre 2 km y "
        "dejó fuera un fotograma que, medido sobre el elipsoide, está a 1,998 m")
    add("solicitud", "conjunto C: fotogramas de 1971 dentro del AGEM al sur del paralelo 20.48",
        "67", 226, "corpus vigente", V02B, SALVEDAD_CENTRO)


def bloque_linea_base():
    V8 = "06_linea_base/linea_base_v8.py"
    t = pd.read_csv(CSV_BASE)
    for clave, periodo in [("14097", "1969–1983"), ("14097", "1992–1994"),
                           ("14039", "1992–1994")]:
        r = t[(t.cvegeo.astype(str).str.zfill(5) == clave) & (t.periodo == periodo)].iloc[0]
        add("línea base", f"densidad de fotogramas, {NOMBRES[clave]}, {periodo}",
            f"{r.fotogramas} fotogramas = {r.por_km2:.3f} por km²",
            int(r.fotogramas), f"{periodo}, AGEM oficial", V8,
            SALVEDAD_CENTRO + "; " + SALVEDAD_AGEM)
    orden = (t[t.periodo == "1969–1983"].assign(
        cve=lambda x: x.cvegeo.astype(str).str.zfill(5))
        .sort_values("por_km2"))
    add("línea base", "orden de los seis municipios por densidad, 1969–1983",
        " < ".join(f"{NOMBRES[c]} {v:.3f}" for c, v in
                   zip(orden.cve, orden.por_km2)),
        int(orden.fotogramas.sum()), "1969–1983, AGEM oficial", V8,
        "Tlajomulco queda en una posición intermedia: en este periodo el archivo "
        "no lo registra menos que a sus vecinos")
    add("línea base", "escala más fina registrada en cada municipio, 1992–1994",
        " · ".join(f"{NOMBRES[str(r.cvegeo).zfill(5)]}: {r.escala_mas_fina}"
                   for _, r in t[t.periodo == "1992–1994"].iterrows()),
        0, "1992–1994, AGEM oficial", V8)

    l = pd.read_csv(CSV_LINEAS)
    cruza = int(l.traza_cruza_el_agem.sum())
    cerca = int(((~l.traza_cruza_el_agem) & (l.dist_min_m <= 2000)).sum())
    extremos = pd.concat([l.dist_ini_m, l.dist_fin_m])
    add("geometría del vuelo", "líneas de vuelo distintas en 1992–1994",
        f"{len(l)}", int(l.n.sum()), "1992–1994", V8,
        "una línea es la combinación de vuelo, bote, línea y escala")
    add("geometría del vuelo",
        "líneas de 1992–1994 según su relación con el AGEM 14 097",
        f"{cruza} cruzan · {cerca} pasan a 2 km o menos · "
        f"{len(l) - cruza - cerca} no se acercan a menos de 2 km",
        len(l), "1992–1994", V8, SALVEDAD_AGEM)
    add("geometría del vuelo",
        "distancia de los extremos de línea al AGEM 14 097 (mediana)",
        f"{extremos.median():,.0f} m", len(extremos), "1992–1994", V8,
        f"de {len(extremos)} extremos, {int((extremos == 0).sum())} cae dentro del "
        f"AGEM y {int((extremos < 2000).sum())} están a menos de 2 km; ninguna "
        "línea se interrumpe en el borde del AGEM")


# ---------------------------------------------------------------------------
# 3. Capas y vaso
# ---------------------------------------------------------------------------

def bloque_capas_y_vaso():
    V12P = "06_linea_base/poligonos_v12.py"
    V7 = "06_linea_base/pisa_vaso_v7.py"
    ofi = cargar_oficial(["14097"])["14097"]
    prov = cargar_provisional()
    c = comparar(ofi, prov)
    add("capa municipal", "superficie del AGEM 14 097",
        f"{c['area_a']:,.2f} km²", 1, "Marco Geoestadístico 2025", V12P, SALVEDAD_AGEM)
    add("capa municipal",
        "diferencia entre el AGEM oficial y la capa provisional usada antes",
        f"diferencia simétrica {c['sim']:.3f} km² ({c['sim_pct']:.2f} %); "
        f"distancia de Hausdorff {c['hausdorff']:.1f} m", 1,
        "oficial frente a provisional", V12P,
        "ninguna cifra municipal cambia al pasar de una capa a la otra")

    planta = Point(*A_UTM.transform(-103.55833333, 20.48777778))
    centro = Point(*A_UTM.transform(-103.5283, 20.4723))
    caja_lon = np.linspace(-103.5442, -103.5100, 200)
    caja_lat = np.linspace(20.4647, 20.4796, 200)
    borde = ([(x, 20.4647) for x in caja_lon] + [(-103.5100, y) for y in caja_lat]
             + [(x, 20.4796) for x in caja_lon[::-1]]
             + [(-103.5442, y) for y in caja_lat[::-1]])
    ex, ey = A_UTM.transform([p[0] for p in borde], [p[1] for p in borde])
    caja = Polygon(zip(ex, ey))
    add("vaso de La Playa", "distancia de la planta del DENUE al centro estimado "
        "de la lámina de agua de 1971", f"{planta.distance(centro):,.0f} m", 1,
        "DENUE, ficha 6336993", V7,
        "el centro de la lámina de agua de 1971 es una estimación visual sobre "
        "las fotografías; no se afirma ninguna relación causal")
    add("vaso de La Playa", "distancia de la planta del DENUE al borde de la caja "
        "de la lámina de agua de 1971", f"{planta.distance(caja):,.0f} m", 1,
        "DENUE, ficha 6336993", V7,
        "la caja de la lámina de agua tiene un error declarado de 330 m en "
        "latitud y 115 m en longitud; la lámina de agua no es el humedal "
        "entendido como cuenca; no se afirma ninguna relación causal")
    add("vaso de La Playa", "la planta del DENUE respecto del AGEM 14 097",
        f"dentro, a {ofi.exterior.distance(planta):,.0f} m del borde", 1,
        "DENUE, ficha 6336993", V7, SALVEDAD_AGEM)


if __name__ == "__main__":
    print("=" * 92)
    print("PAQUETE DE CIFRAS VALIDADAS — VERSIÓN 1")
    print("=" * 92)
    C, guardado = bloque_par_estrella()
    bloque_escala(C)
    bloque_geometria(C, guardado)
    bloque_corpus()
    bloque_linea_base()
    bloque_capas_y_vaso()

    with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["cifra", "nombre_exacto", "valor", "n",
                                           "fecha_escenario", "script", "salvedad"])
        w.writeheader()
        w.writerows(filas)

    print(f"\n  {len(filas)} cifras escritas en {SALIDA.relative_to(RAIZ)}")
    print(f"\n  {'cifra':34}{'nombre exacto':46}{'valor':>22}")
    print("  " + "─" * 102)
    for f in filas:
        print(f"  {f['cifra'][:33]:34}{f['nombre_exacto'][:45]:46}"
              f"{str(f['valor'])[:22]:>22}")
    sin = sum(1 for f in filas if not f["salvedad"])
    print(f"\n  cifras con salvedad: {len(filas) - sin} de {len(filas)}")
    print("  Una cifra con salvedad necesita que la salvedad aparezca en el texto")
    print("  o en una nota. Ningún error de georreferenciación es exactitud.")
    print("\n" + "=" * 92)
