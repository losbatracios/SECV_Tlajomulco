#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Marcas fiduciales y lado del formato
====================================

PROPÓSITO
---------
Separar dos medidas que suelen confundirse en un fotograma aéreo escaneado:

  · el CANTO del área expuesta, que es lo que determina cuánto terreno
    registra la toma, y
  · la distancia entre MARCAS FIDUCIALES opuestas, que es lo que define el
    formato nominal de 23 cm de la cámara métrica.

No coinciden. La emulsión se impresiona más allá de las marcas en unos
fotogramas y queda por dentro en otros, según el modelo de cámara. De esa
diferencia depende el tamaño de píxel sobre la película, y con él el
denominador de escala; no depende de ella, en cambio, la media huella sobre
el terreno, que se obtiene del canto medido y de la resolución del ajuste.

  tamaño de píxel  =  0.23 m / lado del formato (px)        ← usa las fiduciales
  escala efectiva  =  resolución del terreno / tamaño de píxel
  media huella     =  (canto en px / 2) × resolución del terreno   ← sin supuesto

Los dos fotogramas del par llevan sistemas fiduciales distintos:

  1993 (713VE)  cuatro cruces inscritas en un disco, en el punto medio de cada
                lado, por dentro del canto;
  1971 (237BZ)  cuatro cuñas de emulsión que sobresalen del canto hacia la
                banda marginal, sin cruz.

CRITERIOS DE DETECCIÓN
----------------------
Declarados antes de calcular, y no se ajusta ninguna marca a mano.

  cruces (1993)  dentro de una ventana centrada en el punto medio de cada
                 lado, la marca es el objeto claro (luminancia > 200) rodeado
                 de fondo oscuro (media local en 101 px < 70) de mayor área;
                 su posición es el centroide de ese objeto.

  cuñas  (1971)  en una banda de ±400 px alrededor del punto medio de cada
                 lado, el vértice es el punto en que el área expuesta alcanza
                 su posición más exterior, con el umbral del perfil del lado.

El criterio de las cuñas trae además su propia prueba de aceptación, fijada
de antemano: las distancias horizontal y vertical entre vértices opuestos
deben coincidir dentro del 0.3 %, y ningún vértice puede caer a menos de 50 px
del borde del archivo. Si no se cumple, la medida se descarta.

SALIDAS
-------
Posición de cada marca, distancias entre marcas opuestas, tamaño de píxel,
escala efectiva y media huella, por fotograma.

Requiere: numpy, Pillow
"""

import numpy as np
from PIL import Image

from escala_v5b_borde_fotograma import medir, escala_terreno, RASTER
from validacion_par_estrella_v2 import FICHA

Image.MAX_IMAGE_PIXELS = None

LADO_NEGATIVO_M = 0.23
PX_DECLARADO_XML = 18e-6

VENTANA = 900       # lado de la ventana de búsqueda de las cruces
BANDA = 400         # semiancho de la banda de búsqueda de las cuñas
TOL_CUADRATURA = 0.003      # 0.3 %, prueba de aceptación de las cuñas
MARGEN_MINIMO = 50          # px al borde del archivo, prueba de aceptación


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def media_local(a, k=101):
    """Media en una ventana cuadrada de lado k, por imagen integral."""
    s = np.cumsum(np.cumsum(a.astype(np.float64), axis=0), axis=1)
    s = np.pad(s, ((1, 0), (1, 0)))
    h, (A, L) = k // 2, a.shape
    y0 = np.clip(np.arange(A) - h, 0, A); y1 = np.clip(np.arange(A) + h + 1, 0, A)
    x0 = np.clip(np.arange(L) - h, 0, L); x1 = np.clip(np.arange(L) + h + 1, 0, L)
    S = (s[np.ix_(y1, x1)] - s[np.ix_(y0, x1)]
         - s[np.ix_(y1, x0)] + s[np.ix_(y0, x0)])
    return S / ((y1 - y0)[:, None] * (x1 - x0)[None, :])


def componentes(m):
    """Etiquetado de componentes conexas en 4-vecindad, por unión-búsqueda."""
    A, L = m.shape
    et = np.zeros((A, L), dtype=np.int32)
    padre = [0]

    def raiz(i):
        while padre[i] != i:
            padre[i] = padre[padre[i]]; i = padre[i]
        return i

    sig = 1
    for y in range(A):
        for x in np.where(m[y])[0]:
            arr = et[y - 1, x] if y else 0
            izq = et[y, x - 1] if x else 0
            if arr and izq:
                et[y, x] = min(arr, izq)
                ra, rb = raiz(arr), raiz(izq)
                if ra != rb:
                    padre[max(ra, rb)] = min(ra, rb)
            elif arr or izq:
                et[y, x] = arr or izq
            else:
                et[y, x] = sig; padre.append(sig); sig += 1
    for y in range(A):
        for x in np.where(et[y])[0]:
            et[y, x] = raiz(et[y, x])
    return et


# ---------------------------------------------------------------------------
# Detección de cruces (1993)
# ---------------------------------------------------------------------------

def cruces(a, bordes):
    L, R, T, B = bordes
    cx, cy, h = (L + R) // 2, (T + B) // 2, VENTANA // 2
    alto, ancho = a.shape
    salida = {}
    for nombre, (x, y) in {"izquierda": (L, cy), "derecha": (R, cy),
                           "arriba": (cx, T), "abajo": (cx, B)}.items():
        x0, y0 = max(0, x - h), max(0, y - h)
        sub = a[y0:min(alto, y + h), x0:min(ancho, x + h)]
        m = (sub > 200) & (media_local(sub) < 70)
        if not m.any():
            salida[nombre] = None
            continue
        et = componentes(m)
        areas = np.bincount(et.ravel()); areas[0] = 0
        ys, xs = np.where(et == int(areas.argmax()))
        salida[nombre] = (xs.mean() + x0, ys.mean() + y0, int(areas.max()))
    return salida


# ---------------------------------------------------------------------------
# Detección de cuñas (1971)
# ---------------------------------------------------------------------------

def cunas(a, bordes, umbral, holgura=500):
    L, R, T, B = bordes
    cx, cy = (L + R) // 2, (T + B) // 2
    salida = {}
    for nombre in ("izquierda", "derecha", "arriba", "abajo"):
        mejor = None
        if nombre in ("izquierda", "derecha"):
            for y in range(cy - BANDA, cy + BANDA):
                ini = max(0, L - holgura)
                seg = a[y, ini:R + holgura] > umbral
                if not seg.any():
                    continue
                i = np.where(seg)[0]
                x = ini + (i[0] if nombre == "izquierda" else i[-1])
                if mejor is None or (x < mejor[0] if nombre == "izquierda"
                                     else x > mejor[0]):
                    mejor = (x, y)
        else:
            for x in range(cx - BANDA, cx + BANDA):
                ini = max(0, T - holgura)
                seg = a[ini:B + holgura, x] > umbral
                if not seg.any():
                    continue
                i = np.where(seg)[0]
                y = ini + (i[0] if nombre == "arriba" else i[-1])
                if mejor is None or (y < mejor[1] if nombre == "arriba"
                                     else y > mejor[1]):
                    mejor = (x, y)
        salida[nombre] = mejor
    return salida


def acepta_cunas(p, forma):
    """Prueba de aceptación de las cuñas, fijada antes de calcular."""
    alto, ancho = forma
    dh = p["derecha"][0] - p["izquierda"][0]
    dv = p["abajo"][1] - p["arriba"][1]
    fallos = []
    if abs(dh / dv - 1) > TOL_CUADRATURA:
        fallos.append(f"horizontal y vertical difieren {100*abs(dh/dv-1):.2f} % "
                      f"(tolerancia {100*TOL_CUADRATURA:.1f} %)")
    for nombre, (x, y) in p.items():
        if min(x, y, ancho - x, alto - y) < MARGEN_MINIMO:
            fallos.append(f"el vértice {nombre} cae a menos de "
                          f"{MARGEN_MINIMO} px del borde del archivo")
    return (dh + dv) / 2, fallos


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------

def informe(fecha, tipo):
    f = FICHA[fecha]
    terreno = escala_terreno(fecha)
    m = medir(fecha, escala_lectura=1)
    bordes = (int(m["izquierda"]), int(m["derecha"]),
              int(m["arriba"]), int(m["abajo"]))
    canto = (m["ancho_img"] + m["alto_img"]) / 2
    media_huella = canto / 2 * terreno

    print(f"\n{'─' * 86}")
    print(f"{fecha} — {f['bote']} línea {f['linea']} foto {f['foto']}"
          f"   (nominal 1:{f['escala_nominal']:,};"
          f" fiduciales de {tipo})")
    print(f"{'─' * 86}")
    print(f"    canto del área de imagen      {canto:,.0f} px")
    print(f"    resolución del terreno        {terreno:.4f} m/px   (del ajuste)")
    print(f"    MEDIA HUELLA                  {media_huella:,.0f} m"
          f"   (no depende del formato del negativo)")

    a = np.asarray(Image.open(RASTER[fecha]).convert("L"), dtype=np.uint8)

    if tipo == "cruz":
        p = cruces(a, bordes)
        print(f"\n    marca          posición en el archivo      área (px)")
        print("    " + "─" * 56)
        for k, v in p.items():
            print(f"    {k:11}    ({v[0]:9,.1f}, {v[1]:9,.1f})  {v[2]:>10,}")
        dh = p["derecha"][0] - p["izquierda"][0]
        dv = p["abajo"][1] - p["arriba"][1]
        lado = (dh + dv) / 2
        print(f"\n    entre marcas opuestas: horizontal {dh:,.0f} px"
              f" · vertical {dv:,.0f} px  (relación {dh/dv:.4f})")
        aceptado = True
    else:
        p = cunas(a, bordes, umbral=int(
            (m["umbral_filas"] + m["umbral_cols"]) / 2))
        print(f"\n    vértice        posición en el archivo")
        print("    " + "─" * 44)
        for k, v in p.items():
            print(f"    {k:11}    ({v[0]:9,}, {v[1]:9,})")
        lado, fallos = acepta_cunas(p, a.shape)
        dh = p["derecha"][0] - p["izquierda"][0]
        dv = p["abajo"][1] - p["arriba"][1]
        print(f"\n    entre vértices opuestos: horizontal {dh:,} px"
              f" · vertical {dv:,} px  (relación {dh/dv:.4f})")
        aceptado = not fallos
        if fallos:
            print(f"\n    LA MEDIDA NO PASA LA PRUEBA DE ACEPTACIÓN:")
            for x in fallos:
                print(f"      · {x}")
            print(f"    La banda marginal de este fotograma tiene rayas y")
            print(f"    emulsión velada que el criterio confunde con el")
            print(f"    vértice. Se descarta y se conserva la medida del canto.")

    print(f"\n    {'medida':38}{'lado (px)':>12}{'píxel (µm)':>13}"
          f"{'escala efectiva':>19}")
    print("    " + "─" * 82)
    filas = [("canto del área de imagen", canto, True)]
    filas.append((f"entre marcas fiduciales ({tipo})", lado, aceptado))
    for etiqueta, lado_px, usar in filas:
        px = LADO_NEGATIVO_M / lado_px
        marca = "" if usar else "   ← descartada"
        print(f"    {etiqueta:38}{lado_px:>12,.0f}{px*1e6:>13.2f}"
              f"{'1:' + format(terreno/px, ',.0f'):>19}{marca}")
    print(f"    {'declarado en el XML del acervo':38}{'—':>12}"
          f"{PX_DECLARADO_XML*1e6:>13.2f}"
          f"{'1:' + format(terreno/PX_DECLARADO_XML, ',.0f'):>19}")
    return dict(canto=canto, lado=lado, aceptado=aceptado,
                terreno=terreno, media_huella=media_huella)


if __name__ == "__main__":
    print("=" * 86)
    print("MARCAS FIDUCIALES, LADO DEL FORMATO Y ESCALA EFECTIVA")
    print("=" * 86)
    r = {"1993": informe("1993", "cruz"), "1971": informe("1971", "cuña")}

    print("\n" + "=" * 86)
    print("RESUMEN")
    print("=" * 86)
    print("""
  Sin supuesto sobre el formato del negativo
    resolución del terreno y media huella salen del ajuste y del canto medido.

  Con el supuesto de un formato de 23 cm
    el tamaño de píxel y el denominador de escala dependen de qué distancia
    se tome como lado del formato. Donde las fiduciales se miden, esa es la
    distancia buena; donde no, el canto deja el denominador como cota.
""")
    print(f"  {'':10}{'media huella':>14}{'resolución':>13}"
          f"{'escala efectiva':>19}")
    print("  " + "─" * 58)
    for fecha in ("1971", "1993"):
        d = r[fecha]
        lado = d["lado"] if d["aceptado"] else d["canto"]
        den = d["terreno"] / (LADO_NEGATIVO_M / lado)
        nota = "" if d["aceptado"] else "  (cota, del canto)"
        print(f"  {fecha:10}{d['media_huella']:>12,.0f} m"
              f"{d['terreno']:>11.4f} m/px"
              f"{'1:' + format(den, ',.0f'):>19}{nota}")

    print("\n" + "=" * 86)
    print("FIN DEL SCRIPT")
    print("=" * 86)
