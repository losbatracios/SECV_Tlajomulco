#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_inspeccion_calidad_corpus.py
Proyecto Tlajomulco — capítulo "Mirada extractiva"

Inspección de calidad del corpus de 98 fotogramas del SINFA-INEGI
(entrega de agosto de 2026) y validación de la tabla maestra de
correspondencia contra los nombres de carpeta entregados por el INEGI.

El script hace tres cosas, todas reproducibles y sin intervención manual:

  1. VALIDACIÓN DE CORRESPONDENCIA
     Reempareja cada fila de la tabla maestra con su carpeta en disco
     usando la clave explícita (BOTE, LÍNEA, FOTO, AÑO, ESCALA) leída del
     nombre de carpeta, en lugar de la posición en la lista. Reporta toda
     fila cuya carpeta asignada no corresponda a su fotograma.

  2. BARRIDO DE CALIDAD
     Para cada fotograma calcula dimensiones, tamaño de píxel en la
     película, resolución en el terreno (GSD), estadísticos radiométricos
     y un proxy de textura, y genera un panorama reducido para revisión
     visual. Escribe un CSV con una fila por fotograma.

  3. DIAGNÓSTICO DE CANALES
     Verifica si los archivos codificados en RGB contienen información
     cromática real o solo una dominante uniforme (envejecimiento de la
     película / balance del escáner).

Convención de nombre de carpeta del INEGI:
    BOTE_LÍNEA_FOTO_ROLLO_VUELO-CARTA_AÑO_SEC_ESCALA_RUMBO
    ej. 0713VE_0001_029_0116_GUADAL-JAL_1993_02_0010K_123
    El campo ESCALA viene en miles: 0010K = 1:10,000.

Uso:
    python3 03_inspeccion_calidad_corpus.py <ruta_a_Imagenes> <ruta_salida>

Salidas en <ruta_salida>:
    calidad_corpus.csv          una fila por fotograma
    correspondencia_tabla.csv   resultado de la validación fila a fila
    panoramas/                  panorama reducido por fotograma

Requiere: Pillow, openpyxl.

Nota sobre el proxy de textura (lapvar): mide la varianza del laplaciano
en el centro del panorama reducido. Es útil para ORDENAR el corpus, no
para descartar. Un valor bajo puede significar una toma poco nítida, pero
también una superficie homogénea (un cuerpo de agua que ocupa medio
fotograma). Ninguna imagen debe descartarse por esta métrica sin revisión
visual.
"""

import csv
import os
import sys

from PIL import Image

Image.MAX_IMAGE_PIXELS = None      # los fotogramas superan el límite por defecto

LADO_NEGATIVO_UM = 230_000          # negativo métrico de 23 x 23 cm, en micras
LADO_PANORAMA = 1500                # px del lado mayor del panorama reducido
LADO_VENTANA_TEXTURA = 500          # px de la ventana central del proxy de textura


# ---------------------------------------------------------------- utilidades

def parsear_nombre(nombre):
    """Extrae los campos identificatorios del nombre de carpeta del INEGI."""
    p = nombre.split("_")
    return {
        "bote": p[0].lstrip("0"),
        "linea": int(p[1]),
        "foto": int(p[2]),
        "rollo_sec": int(p[3]),
        "vuelo_carta": p[4],
        "anio": int(p[5]),
        "escala": int(p[7].rstrip("K")) * 1000,
        "rumbo": int(p[8]),
    }


def leer_tabla_maestra(ruta_xlsx, hoja="Tabla_maestra"):
    import openpyxl
    ws = openpyxl.load_workbook(ruta_xlsx, read_only=True)[hoja]
    filas = list(ws.values)
    cab = list(filas[0])
    return [dict(zip(cab, f)) for f in filas[1:] if f[0] is not None]


# ------------------------------------------------- 1. validación de correspondencia

def validar_correspondencia(tabla, carpetas):
    """
    Reempareja cada fila con su carpeta por clave explícita.
    Devuelve (resultados, incidencias).
    """
    indice = {}
    for c in carpetas:
        d = parsear_nombre(c)
        indice.setdefault(
            (d["bote"], d["linea"], d["foto"], d["anio"], d["escala"]), []
        ).append(c)

    duplicadas = {k: v for k, v in indice.items() if len(v) > 1}
    if duplicadas:
        raise RuntimeError(f"Claves de carpeta duplicadas: {duplicadas}")

    resultados, incidencias = [], []
    for r in tabla:
        clave = (
            str(r["BOTE"]).lstrip("0"),
            int(r["LINEA"]),
            int(r["FOTO"]),
            int(r["AÑO"]),
            int(str(r["ESCALA"]).split(":")[1]),
        )
        correcta = indice.get(clave, [None])[0]
        asignada = r.get("nombre_carpeta_disco")
        resultados.append(
            {
                "llave": r["llave"],
                "conjunto": r["conjunto"],
                "carpeta_asignada": asignada,
                "carpeta_por_clave": correcta,
                "coincide": correcta == asignada,
            }
        )
        if correcta is None:
            incidencias.append((r["llave"], "sin carpeta para la clave"))
        elif correcta != asignada:
            incidencias.append((r["llave"], f"asignada {asignada} -> correcta {correcta}"))

    asignadas = [x["carpeta_por_clave"] for x in resultados]
    if sorted(a for a in asignadas if a) != sorted(carpetas):
        incidencias.append(("GLOBAL", "la correspondencia por clave no es biyectiva"))

    return resultados, incidencias


# ------------------------------------------------------------ 2. barrido de calidad

def perfil_fotograma(ruta_jpg, nombre, dir_panoramas, conjunto):
    """Estadísticos de calidad de un fotograma + panorama reducido."""
    meta = parsear_nombre(nombre)
    im = Image.open(ruta_jpg)
    ancho, alto = im.size
    modo = im.mode

    um_por_px = LADO_NEGATIVO_UM / max(ancho, alto)
    gsd_cm = um_por_px * 1e-6 * meta["escala"] * 100      # cm de terreno por píxel
    huella_m = 0.23 * meta["escala"]                       # lado del fotograma en el terreno

    im.draft("L", (LADO_PANORAMA, LADO_PANORAMA))          # decodifica a escala reducida
    im = im.convert("L")
    im.thumbnail((LADO_PANORAMA, LADO_PANORAMA), Image.LANCZOS)
    os.makedirs(dir_panoramas, exist_ok=True)
    im.save(os.path.join(dir_panoramas, f"{conjunto}_{nombre}.jpg"), quality=80)

    h = im.histogram()
    total = sum(h)
    media = sum(v * c for v, c in enumerate(h)) / total
    varianza = sum((v - media) ** 2 * c for v, c in enumerate(h)) / total

    acumulado = p1 = p99 = None
    acumulado = 0
    for v, c in enumerate(h):
        acumulado += c
        if p1 is None and acumulado >= 0.01 * total:
            p1 = v
        if p99 is None and acumulado >= 0.99 * total:
            p99 = v

    w, hh = im.size
    s = LADO_VENTANA_TEXTURA
    vent = im.crop((w // 2 - s // 2, hh // 2 - s // 2, w // 2 + s // 2, hh // 2 + s // 2))
    px = vent.load()
    lap = [
        4 * px[x, y] - px[x - 1, y] - px[x + 1, y] - px[x, y - 1] - px[x, y + 1]
        for y in range(1, s - 1, 2)
        for x in range(1, s - 1, 2)
    ]
    m = sum(lap) / len(lap)
    lapvar = sum((v - m) ** 2 for v in lap) / len(lap)

    return {
        "conjunto": conjunto,
        "carpeta": nombre,
        "anio": meta["anio"],
        "escala_denom": meta["escala"],
        "rumbo": meta["rumbo"],
        "ancho_px": ancho,
        "alto_px": alto,
        "mpx": round(ancho * alto / 1e6, 1),
        "modo": modo,
        "bytes": os.path.getsize(ruta_jpg),
        "um_por_px": round(um_por_px, 2),
        "gsd_cm": round(gsd_cm, 1),
        "huella_m": round(huella_m),
        "media_gris": round(media, 1),
        "desv_gris": round(varianza ** 0.5, 1),
        "p1": p1,
        "p99": p99,
        "rango_p1_p99": p99 - p1,
        "pct_claro_ge250": round(100 * sum(h[250:]) / total, 3),
        "pct_oscuro_le5": round(100 * sum(h[:6]) / total, 3),
        "lapvar_centro": round(lapvar, 1),
    }


# --------------------------------------------------------- 3. diagnóstico de canales

def diagnostico_canales(ruta_jpg, umbral_luminancia=30):
    """
    Distingue color real de dominante uniforme en archivos RGB.
    Usa diferencias medianas entre canales, no extremos: un solo píxel de
    una raya o una mota de polvo dispara el extremo sin que exista color.
    Excluye el marco negro del fotograma por luminancia.
    """
    im = Image.open(ruta_jpg)
    im.draft("RGB", (600, 600))
    im = im.convert("RGB")
    im.thumbnail((600, 600), Image.LANCZOS)
    px = [p for p in im.getdata() if sum(p) / 3 > umbral_luminancia]
    if not px:
        return None
    n = len(px)
    d_rg = sorted(abs(p[0] - p[1]) for p in px)
    d_gb = sorted(abs(p[1] - p[2]) for p in px)
    return {
        "media_R": round(sum(p[0] for p in px) / n, 1),
        "media_G": round(sum(p[1] for p in px) / n, 1),
        "media_B": round(sum(p[2] for p in px) / n, 1),
        "dif_RG_mediana": d_rg[n // 2],
        "dif_GB_mediana": d_gb[n // 2],
        "cromatico": d_rg[n // 2] > 15 and d_gb[n // 2] > 15,
    }


# ------------------------------------------------------------------------ principal

def main(dir_imagenes, dir_salida):
    os.makedirs(dir_salida, exist_ok=True)
    ruta_tabla = os.path.join(dir_imagenes, "tabla_maestra_98_imagenes.xlsx")
    tabla = leer_tabla_maestra(ruta_tabla)
    carpetas = [
        d for d in os.listdir(dir_imagenes)
        if os.path.isdir(os.path.join(dir_imagenes, d)) and not d.startswith("_")
    ]

    print(f"Filas en la tabla maestra: {len(tabla)}")
    print(f"Carpetas en disco:         {len(carpetas)}")

    resultados, incidencias = validar_correspondencia(tabla, carpetas)
    with open(os.path.join(dir_salida, "correspondencia_tabla.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, list(resultados[0].keys()))
        w.writeheader()
        w.writerows(resultados)
    print(f"Filas con carpeta mal asignada: {sum(1 for r in resultados if not r['coincide'])}")
    for llave, detalle in incidencias:
        print(f"  {llave}: {detalle}")

    conjunto_de = {r["carpeta_por_clave"]: r["conjunto"] for r in resultados}
    perfiles = []
    for nombre in sorted(carpetas):
        jpg = os.path.join(dir_imagenes, nombre, nombre + ".jpg")
        perfiles.append(
            perfil_fotograma(jpg, nombre, os.path.join(dir_salida, "panoramas"),
                             conjunto_de.get(nombre, "?"))
        )
    with open(os.path.join(dir_salida, "calidad_corpus.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, list(perfiles[0].keys()))
        w.writeheader()
        w.writerows(perfiles)
    print(f"Fotogramas perfilados: {len(perfiles)}")

    rgb = [p for p in perfiles if p["modo"] == "RGB"]
    print(f"Archivos codificados en RGB: {len(rgb)}")
    for p in rgb[:10]:
        d = diagnostico_canales(os.path.join(dir_imagenes, p["carpeta"], p["carpeta"] + ".jpg"))
        print(f"  {p['carpeta'][:38]:40s} RGB={d['media_R']}/{d['media_G']}/{d['media_B']}"
              f"  dif medianas RG={d['dif_RG_mediana']} GB={d['dif_GB_mediana']}"
              f"  cromático={d['cromatico']}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
