#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paquete de cifras validadas — mantenimiento de la versión 2
===========================================================

PROPÓSITO
---------
Regenera cifras_validadas_v2.csv con tres cambios y sin perder ninguna fila:

  1. Añade una columna "estado". Las cifras que el borrador ya no usa quedan
     marcadas como retiradas, con el motivo y la cifra que las sustituye en la
     columna "salvedad". Ninguna fila se borra: el CSV es el registro de lo que
     se ha replicado, incluida la cifra que dejó de usarse.

  2. Añade las cifras de la línea base v8 que el borrador cita y que no
     estaban en el CSV.

  3. Añade las cifras de la huella nominal sobre el AGEM y las del inventario
     Landsat MSS.

Las 77 filas anteriores no se escriben a mano: se regeneran llamando a los
mismos bloques de cifras_validadas_v1.py y cifras_validadas_v2.py, de modo que
siguen saliendo de los scripts que las produjeron. Las cifras nuevas se
calculan aquí sobre el índice de vuelos, el AGEM y las salidas de los scripts
de esta ronda; ninguna se transcribe de un reporte.

El script lleva versión propia (v2b) y el producto conserva su nombre
(cifras_validadas_v2.csv), que es el que cita el capítulo.

Requiere: numpy, pandas, pyproj, shapely, pyshp, pyogrio
"""

import csv
import sys
from pathlib import Path

import pandas as pd
from shapely.geometry import Point

AQUI = Path(__file__).resolve().parent                  # 04_Análisis
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI / "Scripts_Python/06_linea_base"))

import cifras_validadas_v1 as v1                        # noqa: E402
import cifras_validadas_v2 as v2                        # noqa: E402
from poligonos_v12 import cargar_oficial, NOMBRES       # noqa: E402
from linea_base_v8 import cargar, CASO                  # noqa: E402

SALIDA = AQUI / "cifras_validadas_v2.csv"
BASE = AQUI / "Visualizaciones/04_linea_base"
MSS = AQUI / "Visualizaciones/06_mss"

V_V8 = "06_linea_base/linea_base_v8.py"
V_HUELLA = "06_linea_base/huella_v1.py"
V_MSS = "08_mss/mss_contraste_v1.py"

COLUMNAS = ["cifra", "nombre_exacto", "valor", "n", "fecha_escenario",
            "script", "salvedad", "estado"]

SALVEDAD_CENTRO = ("pertenencia por centro de toma declarado en el índice; un "
                   "fotograma con el centro fuera puede registrar territorio "
                   "del municipio con su huella")
SALVEDAD_HUELLA = ("huella nominal: cuadrado de medio lado 0.115 m × "
                   "denominador de la escala declarada, centrado en el centro "
                   "de toma; no corrige inclinación ni relieve; alineada a los "
                   "ejes de la proyección, con la variante girada al rumbo "
                   "reportada aparte")
SALVEDAD_MSS = ("el material disponible son 12 mosaicos anuales, no un ráster "
                "por escena: la serie tiene 12 observaciones y no 132; los "
                "chips vienen en 8 bits ya estirados, censurados en los dos "
                "extremos")
SALVEDAD_MSS_OSCURO = ("un área oscura en el infrarrojo cercano a 60 m no "
                       "equivale a una superficie de agua abierta: también la "
                       "producen suelo húmedo, vegetación densa y sombra")

# Cifras que el borrador ya no usa. La clave es el nombre exacto.
RETIRADAS = {
    "distancia de la planta del DENUE al borde de la caja de la lámina de agua de 1971":
        ("medida contra la caja envolvente de agosto, de 5.88 km², y no contra "
         "el contorno trazado, de 4.16 km²; la sustituye la distancia de 1,182 m "
         "al borde del vaso"),
    "distancia de la planta del DENUE al centro estimado de la lámina de agua de 1971":
        ("medida contra el centro de la caja envolvente de agosto; la sustituye "
         "la distancia de 3,536 m al centroide del vaso"),
    "ventana en que se puede acotar la desaparición de la lámina":
        ("la ventana sale del capítulo por decisión del 2-11; el inventario "
         "Landsat MSS no la acota, porque su material son 12 mosaicos anuales "
         "y no permite fechar nada"),
}


def add(cifra, nombre, valor, n, escenario, script, salvedad=""):
    v1.filas.append(dict(cifra=cifra, nombre_exacto=nombre, valor=valor, n=n,
                         fecha_escenario=escenario, script=script,
                         salvedad=salvedad))


def denominador_texto(den):
    return f"1:{int(den):,}"


# ---------------------------------------------------------------------------
# Cifras de la línea base v8 que el borrador cita y no estaban en el CSV
# ---------------------------------------------------------------------------

def bloque_v8_faltantes(df, agem):
    caso = agem[CASO]

    # 2003
    d03 = df[df["Año"] == 2003]
    reparto = d03["mun"].map(lambda c: NOMBRES.get(c, "fuera de los seis")) \
                        .value_counts()
    add("la mirada", "fotogramas del corpus de 2003 y municipio de su centro",
        f"{len(d03)}, " + " · ".join(f"{k}: {v}" for k, v in reparto.items()),
        len(d03), "2003", V_V8, SALVEDAD_CENTRO)

    # Escala más fina en Tlajomulco, 1969–1983
    dtj = df[(df["mun"] == CASO) & (df["Año"].between(1969, 1983))
             & df["den"].notna()]
    fina = dtj["den"].min()
    sub = dtj[dtj["den"] == fina]
    vuelos = sorted(set(sub["Vuelo"].astype(str)))
    anios = sorted(set(sub["Año"].astype(int)))
    add("la mirada",
        "escala más fina registrada dentro del AGEM 14 097 en 1969–1983",
        f"{denominador_texto(fina)}, vuelo {' / '.join(vuelos)} de "
        f"{' y '.join(str(a) for a in anios)} ({len(sub)} fotogramas)",
        len(dtj), "1969–1983", V_V8,
        "escala nominal declarada en el índice, no escala efectiva medida; "
        + SALVEDAD_CENTRO)

    # Bloques de 1992 y 1993. Un bloque es la pareja bote y escala: el bote
    # 691VE vuela dos escalas y sus dos series están a distancias distintas,
    # así que agrupar solo por bote escondería la mayor de las dos.
    g = df[(df["Bote"].astype(str).isin(["676VE", "691VE", "712VE"]))
           & (df["Año"].between(1992, 1993))].copy()
    g["d"] = [Point(x, y).distance(caso) for x, y in zip(g["E"], g["N"])]
    tramos = (g.groupby([g["Bote"].astype(str), "Escala"])["d"]
                .agg(["min", "size"]).reset_index())
    tramos = tramos.sort_values("min")
    add("la mirada",
        "distancia mínima al AGEM 14 097 del centro de toma más próximo de "
        "cada bloque de 1992 y 1993, por bote y escala",
        f"de {tramos['min'].min() / 1000:.1f} a {tramos['min'].max() / 1000:.1f} km · "
        + " · ".join(f"{r.Bote} {r.Escala} {r['min'] / 1000:.1f} km "
                     f"({int(r['size'])} fotogramas)"
                     for _, r in tramos.iterrows()),
        int(tramos["size"].sum()), "1992–1993", V_V8,
        "no incluye el bote 713VE de 1993, cuya línea 1 sí entra en el AGEM; "
        + SALVEDAD_CENTRO)

    # Barrido TEQUILA de 1994
    g = df[(df["Vuelo"].astype(str) == "TEQUILA") & (df["Año"] == 1994)]
    dmin = min(Point(x, y).distance(caso) for x, y in zip(g["E"], g["N"]))
    add("la mirada",
        "distancia mínima al AGEM 14 097 del centro de toma del barrido TEQUILA",
        f"{dmin:,.0f} m", len(g), "1994", V_V8, SALVEDAD_CENTRO)

    # Años con material en el corpus
    anios = sorted(set(df["Año"].astype(int)))
    add("cobertura temporal", "años distintos con material en el corpus",
        f"{len(anios)}: " + " · ".join(str(a) for a in anios),
        len(df), "corpus del 07-07", V_V8,
        "el corpus es la consulta por cuadrante, no el acervo completo del SINFA")

    # Línea 713VE/1
    g = df[(df["Bote"].astype(str) == "713VE") & (df["Línea"].astype(str) == "1")
           & (df["Año"] == 1993)].sort_values("Número de foto")
    fotos = list(g["Número de foto"].astype(int))
    dentro = [int(f) for f, x, y in zip(g["Número de foto"], g["E"], g["N"])
              if caso.contains(Point(x, y))]
    posiciones = [fotos.index(f) + 1 for f in dentro]
    add("geometría del vuelo",
        "fotogramas de la línea 713VE/1 con el centro de toma dentro del "
        "AGEM 14 097, y su posición en la línea",
        f"{len(dentro)} de {len(fotos)}: fotos {' y '.join(str(f) for f in dentro)}"
        f", en las posiciones {' y '.join(str(p) for p in posiciones)} de "
        f"{len(fotos)}; la línea termina dentro del AGEM",
        len(fotos), "1993", V_V8, SALVEDAD_CENTRO)
    return len(dentro)


# ---------------------------------------------------------------------------
# Huella nominal sobre el AGEM
# ---------------------------------------------------------------------------

def bloque_huella():
    t = pd.read_csv(BASE / "huella_fotogramas_post1992.csv")
    r = pd.read_csv(BASE / "huella_resumen_anio_escala.csv")
    li = pd.read_csv(BASE / "huella_lineas_1992_1994.csv")

    tocan = t[t["toca_ejes"]]
    por_anio = []
    for anio in sorted(t["anio"].unique()):
        g = t[t["anio"] == anio]
        por_anio.append(f"{anio}: {int(g['toca_ejes'].sum())} de {len(g)}")
    add("huella sobre el municipio",
        "fotogramas posteriores a la reforma cuya huella nominal toca el "
        "AGEM 14 097",
        f"{len(tocan)} de {len(t)} · " + " · ".join(por_anio),
        len(t), "1992, 1993, 1994 y 2003", V_HUELLA, SALVEDAD_HUELLA)

    g93 = r[r["anio"] == 1993]
    add("huella sobre el municipio",
        "reparto por escala de los fotogramas de 1993 cuya huella toca el "
        "AGEM 14 097",
        " · ".join(f"{row.escala} ({row.vuelo}): {int(row.huella_toca)} de "
                   f"{int(row.fotogramas)}" for row in g93.itertuples()),
        int(g93["fotogramas"].sum()), "1993", V_HUELLA, SALVEDAD_HUELLA)

    add("huella sobre el municipio",
        "mayor superficie del AGEM 14 097 registrada por un solo fotograma",
        f"{t['area_ejes_km2'].max():.2f} km²",
        len(tocan), "1993, 1:75,000, bote 720VE", V_HUELLA, SALVEDAD_HUELLA)

    add("huella sobre el municipio",
        "distancia mínima al AGEM 14 097 de la huella nominal, por año sin "
        "contacto",
        " · ".join(f"{anio}: {int(t[t['anio'] == anio]['dist_huella_agem_m'].min()):,} m"
                   for anio in (1992, 1994, 2003)),
        int(len(t[t["anio"].isin([1992, 1994, 2003])])),
        "1992, 1994 y 2003", V_HUELLA, SALVEDAD_HUELLA)

    cuenta = li["clase_por_huella"].value_counts()
    add("geometría del vuelo",
        "líneas de 1992–1994 según su relación con el AGEM 14 097, medidas "
        "por huella",
        " · ".join(f"{int(v)} {k}" for k, v in cuenta.items()),
        int(len(li)), "1992–1994", V_HUELLA,
        "la clasificación por traza de centros de la v8 da 1 que cruza, 2 a "
        "2 km o menos y 45 a más de 2 km; cambian de clase 4 líneas; "
        + SALVEDAD_HUELLA)


def bloque_huella_union():
    """La unión total se lee del resumen por año, escala y vuelo."""
    import importlib
    sys.path.insert(0, str(AQUI / "Scripts_Python/06_linea_base"))
    h = importlib.import_module("huella_v1")
    agem = cargar_oficial()
    df, _ = cargar(agem)
    t, d, caso, km2_agem, rumbos = h.por_fotograma(df, agem)
    total = h.area_union(t[t["toca_ejes"]], d, caso, rumbos, False)
    girada = h.area_union(t[t["toca_giro"]], d, caso, rumbos, True)
    add("huella sobre el municipio",
        "superficie del AGEM 14 097 cubierta por la unión de las huellas "
        "de 1992 en adelante",
        f"{total:.2f} km² = {total / km2_agem * 100:.2f} % del AGEM",
        int(t["toca_ejes"].sum()), "1992 en adelante", V_HUELLA,
        f"con las huellas giradas al rumbo de su línea, {girada:.2f} km² = "
        f"{girada / km2_agem * 100:.2f} %; " + SALVEDAD_HUELLA)
    return total, km2_agem


# ---------------------------------------------------------------------------
# Landsat MSS
# ---------------------------------------------------------------------------

def bloque_mss():
    c = pd.read_csv(MSS / "mss_contraste_anual.csv")
    e = pd.read_csv(MSS / "mss_cruce_escenas_fechas.csv")
    j = pd.read_csv(MSS / "mss_jrc_meses_comunes.csv")

    anios = sorted(c["anio"].astype(int))
    faltan = sorted(set(range(1972, 1986)) - set(anios))
    add("Landsat MSS",
        "escenas del inventario y material con valores por píxel",
        f"{len(e)} escenas en {len(anios)} años, pero {len(c)} imágenes "
        f"distintas (un mosaico anual por año)",
        len(e), f"{e['date'].min()} a {e['date'].max()}", V_MSS, SALVEDAD_MSS)

    add("Landsat MSS", "años sin cobertura en el periodo 1972–1985",
        " y ".join(str(a) for a in faltan), len(anios),
        "1972–1985", V_MSS,
        "la ausencia es de cobertura satelital, no del vaso")

    q = c["prob_contorno_mas_oscuro"].astype(float)
    add("Landsat MSS",
        "probabilidad de que un píxel del contorno del vaso sea más oscuro en "
        "el infrarrojo cercano que un píxel de su anillo de 500 m",
        f"de {q.min():.3f} a {q.max():.3f}, mediana {q.median():.3f}; "
        f"por encima de 0.5 en los {len(c)} años",
        len(c), "1972–1985, un mosaico por año", V_MSS,
        "medida con rangos, invariante a la censura del estiramiento; "
        + SALVEDAD_MSS_OSCURO)

    cumplen = c[c["cumple_criterio_1km2"]]
    no = c[~c["cumple_criterio_1km2"]]
    filas_si = int(e["year"].isin(cumplen["anio"]).sum())
    add("Landsat MSS",
        "años en que el área en el extremo oscuro dentro de la caja alcanza el "
        "parche de 1 km² del protocolo",
        f"{len(cumplen)} de {len(c)}; no lo alcanzan "
        + " · ".join(str(int(a)) for a in no["anio"])
        + f" ({filas_si} filas de la plantilla frente a {len(e) - filas_si})",
        len(c), "1972–1985", V_MSS,
        "el umbral del extremo oscuro lo fija el estiramiento del producto, "
        "declarado antes de ver las imágenes; " + SALVEDAD_MSS)

    obs = j[j["fraccion_observada"] > 0]
    add("Landsat MSS",
        "meses con observación de la serie JRC en la ventana común con el "
        "inventario MSS",
        f"{len(obs)} de {len(j)}: "
        + " · ".join(f"{int(r.year)}-{int(r.month):02d} "
                     f"({r.fraccion_observada:.0%} de la caja, "
                     f"{r.px_agua:.0f} píxeles de agua)" for r in obs.itertuples()),
        len(j), "1984-03 a 1985-12", V_MSS,
        "la serie JRC clasifica agua abierta a 30 m y en esa ventana solo "
        "observa dos meses; " + SALVEDAD_MSS_OSCURO)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 94)
    print("PAQUETE DE CIFRAS VALIDADAS — MANTENIMIENTO DE LA VERSIÓN 2")
    print("=" * 94)

    # 1. Regenerar las 77 filas anteriores con sus propios scripts
    C, guardado = v1.bloque_par_estrella()
    v1.bloque_escala(C)
    v1.bloque_geometria(C, guardado)
    v1.bloque_corpus()
    v1.bloque_linea_base()
    v1.bloque_capas_y_vaso()
    v2.bloque_vaso()
    v2.bloque_cobertura()
    v2.bloque_serie()
    heredadas = len(v1.filas)
    print(f"\n  filas regeneradas de la versión 2: {heredadas}")

    # 2. Cifras nuevas
    print("\n  cifras nuevas de la línea base v8")
    agem = cargar_oficial()
    df, _ = cargar(agem)
    bloque_v8_faltantes(df, agem)
    print("\n  cifras nuevas de la huella nominal")
    bloque_huella()
    bloque_huella_union()
    print("\n  cifras nuevas de Landsat MSS")
    bloque_mss()

    # 3. Estado de cada fila
    for f in v1.filas:
        f["estado"] = "vigente"
    retiradas = 0
    for f in v1.filas:
        motivo = RETIRADAS.get(f["nombre_exacto"])
        if motivo:
            f["estado"] = "retirada"
            f["salvedad"] = (f"RETIRADA: {motivo}"
                             + (f" · salvedad anterior: {f['salvedad']}"
                                if f["salvedad"] else ""))
            retiradas += 1
    if retiradas != len(RETIRADAS):
        print(f"\n  AVISO: se esperaban {len(RETIRADAS)} filas retiradas y se "
              f"marcaron {retiradas}; revisar los nombres exactos")

    with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNAS)
        w.writeheader()
        w.writerows(v1.filas)

    nuevas = len(v1.filas) - heredadas
    print()
    print("=" * 94)
    print(f"  {heredadas} filas de la versión anterior + {nuevas} nuevas "
          f"= {len(v1.filas)} en {SALIDA.name}")
    print(f"  filas marcadas como retiradas: {retiradas}")
    for f in v1.filas:
        if f["estado"] == "retirada":
            print(f"    · {f['nombre_exacto']} = {f['valor']}")
    print(f"\n  cifras nuevas")
    print(f"  {'cifra':28}{'nombre exacto':62}{'valor':>24}")
    print("  " + "─" * 114)
    for f in v1.filas[heredadas:]:
        print(f"  {f['cifra'][:27]:28}{f['nombre_exacto'][:61]:62}"
              f"{str(f['valor'])[:24]:>24}")
    sin = sum(1 for f in v1.filas if not f["salvedad"])
    print(f"\n  cifras con salvedad: {len(v1.filas) - sin} de {len(v1.filas)}")
    print("=" * 94)
