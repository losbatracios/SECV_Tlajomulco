#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Escala y media huella de 1971 con el ajuste vigente
====================================================

PROPÓSITO
---------
Las cifras de escala de la toma de 1971 —resolución en el terreno, media
huella del fotograma y cota de escala efectiva— se calcularon el 14-10 a
partir de un ajuste cuyo conjunto de puntos resultó estar mal identificado.
Aquí se rehacen con los ajustes vigentes, sin volver a ajustar nada: los
conjuntos de puntos se toman de la validación v3, donde se identifican por
coordenadas y no por su posición en una lista.

La toma de 1993 no cambia y se imprime solo para dejarlo comprobado.

CÓMO SE OBTIENE CADA CIFRA
--------------------------
  resolución en el terreno   factor de escala de la transformación de
                             similitud ajustada sobre los puntos de control.
                             No supone nada sobre el negativo.

  media huella               mitad del lado del área de imagen medida en el
                             propio ráster, por la resolución en el terreno.
                             Tampoco supone nada sobre el negativo: es lo que
                             fija hasta dónde registra terreno la toma.

  escala efectiva            resolución en el terreno dividida entre el tamaño
                             de píxel sobre la película. Ese tamaño sí supone
                             un formato de 23 cm, y además hace falta saber
                             qué distancia del fotograma corresponde a esos
                             23 cm. En 1993 se midieron las cuatro marcas
                             fiduciales y esa es la distancia buena. En 1971
                             las marcas son cuñas que no pasaron la prueba de
                             aceptación, de modo que se usa el canto del área
                             de imagen y la cifra queda como COTA, no como
                             medida.

El área de imagen y las marcas fiduciales se midieron el 14-10 sobre los
rásteres originales (ver el script del borde del fotograma y el de fiduciales);
son mediciones del material, independientes de cualquier ajuste.

ESCENARIOS DE 1971
------------------
  28-09 (superado)  los 5 puntos del conjunto mal identificado
  (i)               5 puntos, sin el del cauce
  (ii)              7 puntos activos del archivo vigente  ← el que se cita

Requiere: numpy, pyproj, shapely, pyshp
"""

import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent / "06_linea_base"))

from validacion_par_estrella_v2 import MODELOS, parametros_similitud, FICHA
from validacion_par_estrella_v3 import definir_conjuntos, a_arrays, LADO_IMAGEN_PX

LADO_NEGATIVO_M = 0.23

# Distancia entre marcas fiduciales opuestas, medida el 14-10. En 1971 la
# medición no pasó su prueba de aceptación y no se usa.
FIDUCIALES_PX = {"1993": 13997.0, "1971": None}

# Cifras publicadas el 14-10 con el ajuste superado, para contrastar.
ANTES = dict(resolucion=0.4837, media_huella=3071, escala=26705)


def cifras(clave_fecha, filas):
    mapa, fuente = a_arrays(filas)
    H = MODELOS["similitud (Helmert)"]["ajustar"](fuente, mapa)
    resolucion, _ = parametros_similitud(H)

    lx, ly = LADO_IMAGEN_PX[clave_fecha]
    lado_canto = (lx + ly) / 2
    media_huella = lado_canto / 2 * resolucion

    px_canto = LADO_NEGATIVO_M / lado_canto
    salida = dict(n=len(mapa), resolucion=resolucion, lado_canto=lado_canto,
                  media_huella=media_huella,
                  px_canto=px_canto, escala_canto=resolucion / px_canto)

    fid = FIDUCIALES_PX[clave_fecha]
    if fid:
        px_fid = LADO_NEGATIVO_M / fid
        salida.update(px_fid=px_fid, escala_fid=resolucion / px_fid)
    return salida


def imprimir(titulo, c, nominal, citable=False):
    marca = "   ← la que se cita" if citable else ""
    print(f"\n  {titulo}{marca}")
    print(f"    puntos de control en el ajuste      {c['n']:>10}")
    print(f"    resolución en el terreno            {c['resolucion']:>10.4f} m/px")
    print(f"    media huella del fotograma          {c['media_huella']:>10,.0f} m")
    if "escala_fid" in c:
        print(f"    escala efectiva (marcas fiduciales) {'1:' + format(c['escala_fid'], ',.0f'):>10}"
              f"   ({100*(c['escala_fid']/nominal - 1):+.1f} % frente a la nominal 1:{nominal:,})")
        print(f"    tamaño de píxel sobre la película   {c['px_fid']*1e6:>10.2f} µm")
    else:
        print(f"    escala efectiva (COTA, del canto)   {'1:' + format(c['escala_canto'], ',.0f'):>10}"
              f"   ({100*(c['escala_canto']/nominal - 1):+.1f} % frente a la nominal 1:{nominal:,})")
        print(f"    tamaño de píxel sobre la película   {c['px_canto']*1e6:>10.2f} µm   (cota)")


if __name__ == "__main__":
    print("=" * 92)
    print("ESCALA Y MEDIA HUELLA CON LOS AJUSTES VIGENTES")
    print("=" * 92)
    C = definir_conjuntos()

    print("\n" + "─" * 92)
    print("1971 — 237BZ línea 56 foto 41 (nominal 1:30,000)")
    print("─" * 92)
    print(f"  área de imagen medida en el ráster: "
          f"{LADO_IMAGEN_PX['1971'][0]:,} × {LADO_IMAGEN_PX['1971'][1]:,} px")
    print(f"  marcas fiduciales: cuñas; la medición no pasó la prueba de")
    print(f"  aceptación del 14-10, así que la escala queda como cota.")

    c_err = cifras("1971", C["pe1_err"])
    c_i = cifras("1971", C["pe1_i"])
    c_ii = cifras("1971", C["pe1_ii"])
    imprimir("conjunto del 28-09 (mal identificado) — SUPERADO", c_err, 30000)
    imprimir("(i) 5 puntos, sin el del cauce", c_i, 30000)
    imprimir("(ii) 7 puntos activos del archivo vigente", c_ii, 30000, citable=True)

    print(f"\n  Qué cambia respecto de lo publicado el 14-10")
    print(f"    {'cifra':36}{'14-10':>12}{'vigente (ii)':>14}{'diferencia':>14}")
    print("    " + "─" * 76)
    for nombre, antes, ahora, fmt in [
            ("resolución en el terreno (m/px)", ANTES["resolucion"], c_ii["resolucion"], "{:.4f}"),
            ("media huella (m)", ANTES["media_huella"], c_ii["media_huella"], "{:,.0f}"),
            ("escala efectiva, cota (1:N)", ANTES["escala"], c_ii["escala_canto"], "{:,.0f}")]:
        print(f"    {nombre:36}{fmt.format(antes):>12}{fmt.format(ahora):>14}"
              f"{100*(ahora/antes - 1):>13.2f} %")

    print("\n" + "─" * 92)
    print("1993 — 713VE línea 1 foto 29 (nominal 1:10,000) — sin cambio")
    print("─" * 92)
    c93 = cifras("1993", C["pe2"])
    print(f"  área de imagen medida: "
          f"{LADO_IMAGEN_PX['1993'][0]:,} × {LADO_IMAGEN_PX['1993'][1]:,} px · "
          f"entre marcas fiduciales: {FIDUCIALES_PX['1993']:,.0f} px")
    imprimir("7 puntos de control (sin cambio desde el 14-10)", c93, 10000, citable=True)

    # El modelo elegido para 1993 es el afín, que no tiene un único factor de
    # escala: conviene saber cuánto se separan sus dos ejes del valor de la
    # similitud, del que salen la media huella y la escala efectiva.
    mapa, fuente = a_arrays(C["pe2"])
    Ha = MODELOS["afín"]["ajustar"](fuente, mapa)
    ejes = np.linalg.svd(np.asarray(Ha)[:2, :2], compute_uv=False)
    print(f"\n  Nota sobre el modelo elegido para 1993. La media huella y la escala")
    print(f"  efectiva salen del factor de escala de la similitud ({c93['resolucion']:.4f} m/px).")
    print(f"  El modelo elegido es el afín, que tiene dos factores: "
          f"{ejes[0]:.4f} y {ejes[1]:.4f} m/px")
    print(f"  ({100*(ejes[0]/ejes[1] - 1):.2f} % de diferencia entre ejes). La media huella")
    print(f"  cambiaría en menos de {abs(100*((ejes.mean())/c93['resolucion'] - 1)):.2f} %, "
          f"así que se conserva la cifra de la similitud.")

    print("\n" + "=" * 92)
    print("FIN DEL SCRIPT")
    print("=" * 92)
