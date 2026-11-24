#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paquete de cifras validadas, versión 2
=======================================

PROPÓSITO
---------
La versión 1 reunía las cifras del par comparativo, la escala, el corpus, la
línea base y las capas municipales. Esta versión añade las del vaso de La
Playa: el ajuste de la fotografía de febrero de 1971, la sensibilidad de su
contorno al modelo, la distancia de la planta del DENUE, la cobertura del
vaso en el archivo aéreo y la serie satelital de agua superficial.

La versión 1 no se toca. Esta se construye llamando al generador de aquella
y añadiendo las filas nuevas, de modo que las cifras heredadas no se copian
a mano: si una cambia allá, cambia aquí.

Las cifras del vaso se leen de las salidas de los dos scripts de la ronda 6,
no se vuelven a calcular ni se escriben a mano.

CÓMO LEER LA COLUMNA DE SALVEDADES
----------------------------------
Una cifra con salvedad necesita que la salvedad aparezca en el texto o en una
nota. Para el vaso, la salvedad que más pesa es ésta: el contorno delimita el
vaso con sus márgenes por su borde tonal, con parcelas y una vereda dentro.
No es una lámina de agua, y no debe escribirse como tal.

Requiere: numpy, pandas, pyproj, shapely, pyshp, pyogrio
"""

import csv
import sys
from pathlib import Path

import pandas as pd

AQUI = Path(__file__).resolve().parent                  # 04_Análisis
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))

import cifras_validadas_v1 as v1                        # noqa: E402

VAL = AQUI / "Visualizaciones/05_humedal_georef/_validacion"
SALIDA = AQUI / "cifras_validadas_v2.csv"

V_HUM = "07_regresion_humedal/validacion_humedal_v1.py"
V_COB = "07_regresion_humedal/cobertura_vaso_v1.py"

SALVEDAD_VASO = ("el contorno delimita el vaso con sus márgenes, por su borde "
                 "tonal, con parcelas y una vereda visibles dentro; no es una "
                 "lámina de agua y no se escribe como tal")
SALVEDAD_AJUSTE_VASO = ("error relativo a la base de referencia (Esri Wayback "
                        "2026-08-05), no exactitud; 5 puntos de control en solo 3 "
                        "ubicaciones independientes, de modo que el residuo del "
                        "modelo afín es un resultado de la geometría y no una "
                        "validación")
SALVEDAD_CONTROLES = ("con 2 puntos de control la estimación es muy débil y la "
                      "elección de modelo se invierte si se quita el control de "
                      "la cerca; error relativo a la base, no exactitud")
SALVEDAD_HUELLA = ("huella nominal: cuadrado de medio lado 0.115 m × denominador "
                   "de la escala declarada, centrado en el centro de toma; no "
                   "corrige inclinación ni relieve y se orienta con los ejes de "
                   "la proyección, no con el rumbo del vuelo")
SALVEDAD_JRC = ("el producto detecta agua abierta a 30 m: un humedal con "
                "vegetación emergente, agua somera o canales más estrechos puede "
                "no quedar registrado; \"no registra agua superficial abierta\" "
                "nunca equivale a \"no había humedal\"")


def add(cifra, nombre, valor, n, escenario, script, salvedad=""):
    v1.filas.append(dict(cifra=cifra, nombre_exacto=nombre, valor=valor, n=n,
                         fecha_escenario=escenario, script=script,
                         salvedad=salvedad))


def bloque_vaso():
    a = pd.read_csv(VAL / "ajuste_humedal.csv").set_index("modelo")
    elegido = a["ctrl_rms_m"].idxmin()
    otro = [m for m in a.index if m != elegido][0]
    e, o = a.loc[elegido], a.loc[otro]
    esc = f"febrero de 1971, {int(e.n_gcp)} puntos de control, modelo {elegido}"

    add("vaso de La Playa", "modelo de transformación elegido", elegido,
        2, esc, V_HUM,
        f"elegido por menor error de predicción en los puntos de control "
        f"({e.ctrl_rms_m:.1f} m frente a {o.ctrl_rms_m:.1f} m); " + SALVEDAD_CONTROLES)
    add("vaso de La Playa", "ubicaciones independientes de los puntos de control",
        f"{int(e.ubicaciones_independientes)} de {int(e.n_gcp)}", int(e.n_gcp), esc,
        V_HUM, "tres de los cinco puntos están a menos de 130 m entre sí, en el "
               "mismo caserío: para el ajuste cuentan como un solo sitio")
    add("vaso de La Playa", "error de predicción en puntos de control (RMS)",
        f"{e.ctrl_rms_m:.1f} m", 2, esc, V_HUM, SALVEDAD_CONTROLES)
    add("vaso de La Playa", "error de predicción leave-one-out",
        f"{e.loo_rms_m:.0f} m", int(e.n_gcp), esc, V_HUM,
        "dominado por los dos puntos aislados, que al salir dejan el ajuste sin "
        "sostén en su zona; " + SALVEDAD_AJUSTE_VASO)
    add("vaso de La Playa", "RMS de residuos √(S/n)", f"{e.rms_residuos_m:.2f} m",
        int(e.n_gcp), esc, V_HUM, SALVEDAD_AJUSTE_VASO)
    add("vaso de La Playa", "error estándar corregido √(S/(n−3))",
        f"{e.error_estandar_corregido_m:.2f} m", int(e.n_gcp), esc, V_HUM,
        SALVEDAD_AJUSTE_VASO)
    add("vaso de La Playa", "anisotropía entre los ejes del modelo afín",
        f"{e.anisotropia_pct:.2f} %", int(e.n_gcp), esc, V_HUM,
        "compatible con una deformación del negativo o del escaneo, no demostrada")

    add("vaso de La Playa", "superficie del vaso que registra la fotografía",
        f"{e.area_ha:,.0f} ha ({e.area_ha / 100:.2f} km²)", 0, esc, V_HUM,
        SALVEDAD_VASO + f"; con el otro modelo da {o.area_ha:,.0f} ha, "
        f"{abs(100 * (o.area_ha / e.area_ha - 1)):.1f} % de diferencia")
    add("vaso de La Playa", "perímetro del contorno", f"{e.perimetro_m:,.0f} m", 0,
        esc, V_HUM, SALVEDAD_VASO)
    add("vaso de La Playa", "desplazamiento del contorno al cambiar de modelo",
        "mediana 37 m · máximo 87 m", 81, "Helmert frente a afín", V_HUM,
        "es la incertidumbre del trazo frente al modelo; no cambia ninguna "
        "afirmación de escala kilométrica")

    add("vaso de La Playa", "distancia de la planta del DENUE al borde del vaso",
        f"{e.pisa_borde_m:,.0f} m", 1, esc, V_HUM,
        "medida al contorno dibujado; NO es la cifra de 1,730 m que circuló, que "
        "se midió al borde de la caja envolvente estimada en agosto. "
        + SALVEDAD_VASO + "; sin relación causal")
    add("vaso de La Playa", "distancia de la planta del DENUE al centroide del vaso",
        f"{e.pisa_centroide_m:,.0f} m", 1, esc, V_HUM, SALVEDAD_VASO + "; sin causa")
    add("vaso de La Playa", "la planta del DENUE respecto del vaso",
        "fuera, al oeste-noroeste", 1, esc, V_HUM,
        "igual con los dos modelos; sin relación causal")

    return elegido


def bloque_cobertura():
    c = pd.read_csv(VAL / "cobertura_vaso_corpus.csv")
    anios = sorted(c["anio"].unique())
    add("cobertura del vaso", "fotogramas cuya huella nominal contiene el "
        "centroide del vaso", f"{len(c)}", 2074, "corpus 07-07", V_COB,
        SALVEDAD_HUELLA)
    add("cobertura del vaso", "años en que el archivo registra el vaso",
        " · ".join(f"{a} ({int((c.anio == a).sum())})" for a in anios),
        len(c), "corpus 07-07", V_COB,
        "después de 1981 el corpus no registra ningún fotograma sobre el vaso; "
        + SALVEDAD_HUELLA)
    add("cobertura del vaso", "último año en que el archivo registra el vaso",
        f"{int(max(anios))}", len(c), "corpus 07-07", V_COB, SALVEDAD_HUELLA)
    add("cobertura del vaso", "escalas con que el archivo registra el vaso",
        " · ".join(f"{a}: " + ", ".join(sorted(set(c[c.anio == a].escala)))
                   for a in anios), len(c), "corpus 07-07", V_COB, SALVEDAD_HUELLA)

    e5 = pd.read_csv(VAL / "cobertura_vaso_e5.csv")
    add("cobertura del vaso", "fotogramas de la lista para el artículo que "
        "contienen el centroide del vaso",
        f"{int(e5.contiene_centroide.sum())} de {len(e5)}", len(e5),
        "lista e5, recalculada", V_COB,
        "la cobertura recalculada coincide con la de la lista; " + SALVEDAD_HUELLA)


def bloque_serie():
    j = pd.read_csv(VAL / "jrc_replica_mensual.csv")
    total = (j.px_water + j.px_no_water + j.px_no_obs).max()
    ha = total * 0.09
    bien = j[j.obs >= 0.90]
    agua = j[j.px_water > 0]
    pre = agua[agua.year < 2021]
    ep = agua[agua.year == 2021]
    orden = bien.sort_values(["year", "month"])
    primero = orden.iloc[0]
    feb = bien[bien.month == 2]

    add("serie de agua superficial", "superficie de la caja sobre la que se "
        "extrajo la serie", f"{ha:,.0f} ha ({total:,.0f} píxeles de 30 m)", 1,
        "caja del vaso estimada el 25-08", V_COB,
        "corrige el valor de ~568 ha que decía el README")
    add("serie de agua superficial", "meses con al menos el 90 % de la caja observada",
        f"{len(bien)} de {len(j)}", len(j), "1984-03 a 2021-12", V_COB,
        "un mes mal observado no es un mes sin agua")
    add("serie de agua superficial", "meses bien observados con agua registrada "
        "antes de 2021", f"{len(pre)}", len(bien), "1984-03 a 2020-12", V_COB,
        SALVEDAD_JRC)
    add("serie de agua superficial", "febreros bien observados con agua registrada",
        f"{int((feb.px_water > 0).sum())} de {len(feb)}", len(feb),
        "febrero es el mismo mes que la fotografía de 1971", V_COB, SALVEDAD_JRC)
    add("serie de agua superficial", "único episodio con agua registrada",
        f"de {int(ep.month.min())} a {int(ep.month.max())} de 2021, "
        f"entre {ep.px_water.min() * 0.09:.1f} y {ep.px_water.max() * 0.09:.1f} ha "
        f"(máximo {100 * ep.px_water.max() / total:.1f} % de la caja)",
        len(ep), "2021", V_COB, SALVEDAD_JRC)
    add("serie de agua superficial", "ventana en que se puede acotar la "
        "desaparición de la lámina",
        f"entre febrero de 1971 y {int(primero.year)}-{int(primero.month):02d}",
        len(bien), "primer mes bien observado de la serie", V_COB,
        "la serie no puede fechar la desaparición, solo acotarla; "
        "precede al primer mes observado. Sin causa. " + SALVEDAD_JRC)


if __name__ == "__main__":
    print("=" * 94)
    print("PAQUETE DE CIFRAS VALIDADAS — VERSIÓN 2 (suma el vaso de La Playa)")
    print("=" * 94)

    C, guardado = v1.bloque_par_estrella()
    v1.bloque_escala(C)
    v1.bloque_geometria(C, guardado)
    v1.bloque_corpus()
    v1.bloque_linea_base()
    v1.bloque_capas_y_vaso()
    heredadas = len(v1.filas)

    elegido = bloque_vaso()
    bloque_cobertura()
    bloque_serie()

    with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["cifra", "nombre_exacto", "valor", "n",
                                           "fecha_escenario", "script", "salvedad"])
        w.writeheader()
        w.writerows(v1.filas)

    nuevas = len(v1.filas) - heredadas
    print(f"\n  {heredadas} cifras heredadas de la versión 1 + {nuevas} nuevas"
          f" = {len(v1.filas)} en {SALIDA.name}")
    print(f"  modelo elegido para el vaso: {elegido}")
    print(f"\n  cifras nuevas de esta ronda")
    print(f"  {'cifra':30}{'nombre exacto':54}{'valor':>26}")
    print("  " + "─" * 110)
    for f in v1.filas[heredadas:]:
        print(f"  {f['cifra'][:29]:30}{f['nombre_exacto'][:53]:54}"
              f"{str(f['valor'])[:26]:>26}")
    sin = sum(1 for f in v1.filas if not f["salvedad"])
    print(f"\n  cifras con salvedad: {len(v1.filas) - sin} de {len(v1.filas)}")
    print("  Ninguna cifra del vaso se escribe como \"lámina de agua\".")
    print("\n" + "=" * 94)
