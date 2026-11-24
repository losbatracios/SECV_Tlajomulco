#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tamaño de píxel sobre la película y escala efectiva, medidos en el fotograma
=============================================================================

PROPÓSITO
---------
Determinar el tamaño de píxel sobre la película de dos fotografías aéreas
verticales escaneadas, midiendo en el propio ráster el área de imagen expuesta.

El negativo de una cámara aérea métrica tiene un formato de 23 × 23 cm. El
archivo entregado por el acervo no coincide con ese formato: incluye, además
del área expuesta, una banda marginal con las marcas fiduciales y la
rotulación del vuelo, de anchura distinta en cada lado. De ahí que el ancho
del archivo en píxeles no sirva para deducir el tamaño de píxel: hay que
medir el área de imagen, de borde a borde del fotograma.

Con esa medida se obtienen dos cifras:

  tamaño de píxel sobre la película  =  0.23 m / lado del área de imagen (px)
  escala efectiva (denominador)      =  resolución en el terreno (m/px)
                                        ───────────────────────────────────
                                        tamaño de píxel sobre la película

La resolución en el terreno proviene del ajuste de la georreferenciación
(factor de escala de la transformación de similitud estimada sobre los puntos
de control), no de la escala nominal declarada en el índice de vuelos.

La escala efectiva fija la mitad de la huella del fotograma sobre el terreno
(0.115 m × denominador), que es el radio dentro del cual una toma registra
territorio del municipio aunque su centro caiga fuera.

MÉTODO DE MEDICIÓN
------------------
La banda marginal es oscura y uniforme; el área expuesta es clara y texturada.
Se recorre un perfil de luminancia media por fila (sobre la franja central de
columnas) y por columna (sobre la franja central de filas); la franja central
evita las esquinas, las marcas fiduciales y la rotulación. El umbral se fija a
la mitad entre el nivel del margen (mediana del 1 % exterior del perfil) y el
nivel de la imagen (mediana del 50 % central), y el área de imagen es el tramo
contiguo más largo por encima de ese umbral.

El criterio se declara antes de ver las cifras: el área de imagen es el tramo
claro contiguo mayor, y no se ajusta ninguno de sus cuatro bordes a mano.

SALIDAS
-------
  - Cifras por fotograma: bordes medidos, lado del área de imagen, tamaño de
    píxel, escala efectiva y mitad de la huella.
  - Una imagen de verificación por fotograma: vista general con las cuatro
    líneas medidas y cuatro recortes a resolución completa, uno por borde,
    para confirmación visual.

Requiere: numpy, pandas, pyproj, shapely, Pillow
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from pyproj import Transformer

Image.MAX_IMAGE_PIXELS = None       # los fotogramas superan el límite por defecto

from validacion_par_estrella_v2 import (
    leer_points, ajustar_similitud, parametros_similitud,
    RUTA_1971, RUTA_1993, FICHA, RAIZ, AQUI,
)

LADO_NEGATIVO_M = 0.23              # formato del negativo, 23 × 23 cm
PX_DECLARADO_XML = 18e-6            # "escaneo fotogramétrico a 18 micras",
                                    # texto idéntico en el XML de ambas tomas

CSV_CORPUS = RAIZ / "02_Archivo_INEGI/Metadatos_SINFA/07-07_consulta_ÍndicesdeVuelos260107080134.csv"
POLIGONO = RAIZ / "05_Datos_abiertos/Capas_cartograficas/poligono_tlajomulco_14097.geojson"
IMAGENES = RAIZ / "02_Archivo_INEGI/Imágenes/_pares"
SALIDA = RAIZ / "04_Análisis/Visualizaciones/03_regresion_par_estrella/_borde_fotograma"

# Rásteres originales del par (solo lectura).
RASTER = {
    "1971": IMAGENES / "B_contrapartes_prereforma"
            / "0237BZ_0056_041_0200_F1312-Z15B_1971_02_0030K_270"
            / "0237BZ_0056_041_0200_F1312-Z15B_1971_02_0030K_270.jpg",
    "1993": IMAGENES / "A_1993_municipio"
            / "0713VE_0001_029_0116_GUADAL-JAL_1993_02_0010K_123"
            / "0713VE_0001_029_0116_GUADAL-JAL_1993_02_0010K_123.jpg",
}
PUNTOS = {"1971": RUTA_1971, "1993": RUTA_1993}

# Los 5 puntos de control definitivos de 1971 (el sexto, sobre el río, se
# retiró por una salvedad registrada antes de calcular el ajuste).
INDICES_1971 = [0, 1, 2, 3, 5]

A_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32613", always_xy=True)


# ---------------------------------------------------------------------------
# 1. Medición del área de imagen
# ---------------------------------------------------------------------------

def perfiles(a, fraccion_central=0.60):
    """Luminancia media por fila y por columna, sobre la franja central opuesta."""
    alto, ancho = a.shape
    m = (1.0 - fraccion_central) / 2.0
    c0, c1 = int(ancho * m), int(ancho * (1 - m))
    f0, f1 = int(alto * m), int(alto * (1 - m))
    return a[:, c0:c1].mean(axis=1), a[f0:f1, :].mean(axis=0)


def tramo_claro(perfil):
    """Primer y último índice del tramo claro contiguo más largo del perfil.

    Devuelve además el umbral y los dos niveles usados, para poder auditar la
    decisión sin volver a abrir la imagen.
    """
    n = len(perfil)
    borde = max(1, n // 100)
    nivel_margen = float(np.median(np.concatenate([perfil[:borde], perfil[-borde:]])))
    nivel_imagen = float(np.median(perfil[n // 4: 3 * n // 4]))
    umbral = (nivel_margen + nivel_imagen) / 2.0

    claro = perfil > umbral
    if not claro.any():
        raise ValueError("ningún tramo del perfil supera el umbral")
    idx = np.where(claro)[0]
    tramos = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    mayor = max(tramos, key=len)
    return int(mayor[0]), int(mayor[-1]), umbral, nivel_margen, nivel_imagen


def medir(fecha, escala_lectura=1):
    """Mide el área de imagen de un fotograma.

    `escala_lectura` es el divisor de la lectura (1 = resolución completa).
    Se usa el modo `draft` de JPEG, que decodifica directamente a 1/2, 1/4 u
    1/8 sin pasar por la imagen completa.
    """
    im = Image.open(RASTER[fecha])
    ancho0, alto0 = im.size
    if escala_lectura > 1:
        im.draft("L", (ancho0 // escala_lectura, alto0 // escala_lectura))
    a = np.asarray(im.convert("L"), dtype=np.uint8)
    alto, ancho = a.shape
    factor = ancho0 / ancho

    pf, pc = perfiles(a)
    f_ini, f_fin, u_f, mg_f, mi_f = tramo_claro(pf)
    c_ini, c_fin, u_c, mg_c, mi_c = tramo_claro(pc)

    return dict(
        fecha=fecha, archivo_px=(ancho0, alto0), leido_px=(ancho, alto), factor=factor,
        # bordes en píxeles del archivo original
        arriba=f_ini * factor, abajo=f_fin * factor,
        izquierda=c_ini * factor, derecha=c_fin * factor,
        alto_img=(f_fin - f_ini + 1) * factor,
        ancho_img=(c_fin - c_ini + 1) * factor,
        umbral_filas=u_f, umbral_cols=u_c,
        nivel_margen=(mg_f + mg_c) / 2, nivel_imagen=(mi_f + mi_c) / 2,
        bordes_leidos=(c_ini, c_fin, f_ini, f_fin),
    )


def escala_terreno(fecha):
    """Resolución en el terreno (m/px) del ajuste de similitud de esa toma."""
    mapa, fuente, _ = leer_points(PUNTOS[fecha])
    if fecha == "1971":
        mapa, fuente = mapa[INDICES_1971], fuente[INDICES_1971]
    M = ajustar_similitud(fuente, mapa)
    return parametros_similitud(M)[0]


# ---------------------------------------------------------------------------
# 2. Imagen de verificación
# ---------------------------------------------------------------------------

def imagen_verificacion(med, ancho_vista=1700, lado_recorte=420):
    """Vista general con las líneas medidas + un recorte por borde.

    Los recortes salen a resolución completa, centrados en el punto medio de
    cada borde, con la línea medida superpuesta: permiten confirmar a ojo que
    la línea cae en el canto del fotograma y no en la banda marginal.
    """
    SALIDA.mkdir(parents=True, exist_ok=True)
    fecha = med["fecha"]
    im = Image.open(RASTER[fecha])
    ancho0, alto0 = im.size

    # --- vista general
    # `draft` debe pedirse antes de que el archivo se decodifique: decodifica a
    # 1/2, 1/4 u 1/8 directamente y evita cargar los 200 Mpx en memoria.
    vista = Image.open(RASTER[fecha])
    vista.draft("L", (ancho_vista, int(ancho_vista * alto0 / ancho0)))
    vista = vista.convert("RGB")
    vista = vista.resize((ancho_vista, int(ancho_vista * alto0 / ancho0)),
                         Image.LANCZOS)
    k = vista.size[0] / ancho0
    d = ImageDraw.Draw(vista)
    for x in (med["izquierda"], med["derecha"]):
        d.line([(x * k, 0), (x * k, vista.size[1])], fill=(255, 40, 40), width=3)
    for y in (med["arriba"], med["abajo"]):
        d.line([(0, y * k), (vista.size[0], y * k)], fill=(255, 40, 40), width=3)
    ruta_vista = SALIDA / f"borde_{fecha}_vista.jpg"
    vista.save(ruta_vista, quality=88)

    # --- recortes a resolución completa, uno por borde
    # Se toman al 30 % del recorrido de cada borde, no en su punto medio: ahí
    # está la marca fiducial lateral y taparía la transición que hay que ver.
    im_full = Image.open(RASTER[fecha]).convert("L")
    h = lado_recorte // 2
    x0, x1 = med["izquierda"], med["derecha"]
    y0, y1 = med["arriba"], med["abajo"]
    cx = int(x0 + 0.30 * (x1 - x0))
    cy = int(y0 + 0.30 * (y1 - y0))
    bordes = [
        ("izquierda", int(x0), cy, "v"),
        ("derecha", int(x1), cy, "v"),
        ("arriba", cx, int(y0), "h"),
        ("abajo", cx, int(y1), "h"),
    ]
    panel = Image.new("RGB", (lado_recorte * 2, lado_recorte * 2), (20, 20, 20))
    for i, (nombre, x, y, orient) in enumerate(bordes):
        caja = (max(0, x - h), max(0, y - h),
                min(ancho0, x - h + lado_recorte), min(alto0, y - h + lado_recorte))
        rec = im_full.crop(caja).convert("RGB")
        dr = ImageDraw.Draw(rec)
        if orient == "v":
            px = x - caja[0]
            dr.line([(px, 0), (px, rec.size[1])], fill=(255, 40, 40), width=2)
        else:
            py = y - caja[1]
            dr.line([(0, py), (rec.size[0], py)], fill=(255, 40, 40), width=2)
        dr.text((8, 8), f"{nombre}  ({x}, {y})", fill=(255, 220, 60))
        panel.paste(rec, ((i % 2) * lado_recorte, (i // 2) * lado_recorte))
    ruta_panel = SALIDA / f"borde_{fecha}_recortes.jpg"
    panel.save(ruta_panel, quality=92)
    del im_full
    return ruta_vista, ruta_panel


# ---------------------------------------------------------------------------
# 3. Umbral del conjunto A con la escala medida
# ---------------------------------------------------------------------------

def recuento_conjunto_a(media_huella_m, umbrales_extra):
    """Fotografías de 1993 que registran territorio municipal, por umbral."""
    from shapely.geometry import shape, Point
    from shapely.ops import unary_union
    import json

    gj = json.load(open(POLIGONO, encoding="utf-8"))
    geoms = [shape(f["geometry"]) for f in gj["features"]]
    poli_geo = unary_union(geoms)

    df = pd.read_csv(CSV_CORPUS, encoding="latin-1")
    col_lat = [c for c in df.columns if "atitud" in c][0]
    col_lon = [c for c in df.columns if "ongitud" in c][0]
    col_ano = [c for c in df.columns if "ño" in c or "Ano" in c][0]
    d93 = df[df[col_ano].astype(str).str.strip() == "1993"].copy()

    # Distancia en metros, en EPSG:32613, del centro de toma al límite municipal.
    from shapely.ops import transform as shp_transform
    poli_utm = shp_transform(lambda x, y: A_UTM.transform(x, y), poli_geo)

    filas = []
    for _, r in d93.iterrows():
        lat, lon = float(r[col_lat]), float(r[col_lon])
        p = Point(*A_UTM.transform(lon, lat))
        dist = 0.0 if poli_utm.contains(p) else p.distance(poli_utm)
        filas.append(dict(lat=lat, lon=lon, dist=dist,
                          bote=str(r.get("Bote", "")).strip(),
                          linea=str(r.get("Línea", r.get("Linea", ""))).strip(),
                          foto=str(r.get("Número de foto", "")).strip()))
    t = pd.DataFrame(filas)

    print("\n  Fotografías de 1993 que registran territorio municipal")
    print(f"    {'umbral':<52}{'metros':>9}{'fotografías':>14}")
    print("    " + "─" * 75)
    for nombre, m in umbrales_extra:
        print(f"    {nombre:<52}{m:>9,.0f}{(t.dist <= m).sum():>14}")

    franja = t[(t.dist > 1050) & (t.dist <= 1400)].sort_values("dist")
    print(f"\n  Fotografías en la franja 1,050–1,400 m del límite ({len(franja)})")
    print(f"    {'Bote':>8}{'Línea':>7}{'Foto':>6}{'dist(m)':>9}"
          f"{'latitud':>11}{'longitud':>12}   ¿dentro de la media huella medida?")
    print("    " + "─" * 88)
    for _, r in franja.iterrows():
        dentro = "sí" if r.dist <= media_huella_m else "no"
        print(f"    {r.bote:>8}{r.linea:>7}{r.foto:>6}{r.dist:>9.0f}"
              f"{r.lat:>11.4f}{r.lon:>12.4f}{dentro:>38}")
    return t


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------

def main():
    print("=" * 88)
    print("ÁREA DE IMAGEN, TAMAÑO DE PÍXEL SOBRE LA PELÍCULA Y ESCALA EFECTIVA")
    print("=" * 88)

    resultados = {}
    for fecha in ("1971", "1993"):
        f = FICHA[fecha]
        terreno = escala_terreno(fecha)

        # Se mide dos veces, a resolución completa y a 1/4, para comprobar que
        # la cifra no depende del nivel de lectura.
        m = medir(fecha, escala_lectura=1)
        m4 = medir(fecha, escala_lectura=4)

        lado_px = (m["ancho_img"] + m["alto_img"]) / 2

        # La media huella sobre el terreno se obtiene de dos cantidades medidas
        # —el lado del área de imagen en píxeles y la resolución en el terreno
        # del ajuste— sin pasar por el formato del negativo. Es, por tanto,
        # independiente del supuesto de que el fotograma mida 23 cm.
        media_huella = (lado_px / 2) * terreno

        # El tamaño de píxel sobre la película y la escala efectiva sí dependen
        # de ese supuesto: son la misma medida expresada en unidades de película.
        px_pelicula = LADO_NEGATIVO_M / lado_px
        escala = terreno / px_pelicula

        resultados[fecha] = dict(
            m=m, terreno=terreno, lado_px=lado_px,
            px_pelicula=px_pelicula, escala=escala, media_huella=media_huella,
        )

        print(f"\n{'─' * 88}")
        print(f"{fecha} — {f['bote']} línea {f['linea']} foto {f['foto']}"
              f"   (escala nominal declarada 1:{f['escala_nominal']:,})")
        print(f"{'─' * 88}")
        print(f"    archivo entregado                {m['archivo_px'][0]:,} × {m['archivo_px'][1]:,} px")
        print(f"    niveles del perfil               margen {m['nivel_margen']:.0f}"
              f" · imagen {m['nivel_imagen']:.0f}"
              f" · umbral {(m['umbral_filas']+m['umbral_cols'])/2:.0f}")
        print(f"    bordes medidos (px del archivo)  columnas {m['izquierda']:.0f}–{m['derecha']:.0f}"
              f" · filas {m['arriba']:.0f}–{m['abajo']:.0f}")
        print(f"    área de imagen                   {m['ancho_img']:,.0f} × {m['alto_img']:,.0f} px"
              f"   (relación {m['ancho_img']/m['alto_img']:.4f})")
        print(f"    misma medida leyendo a 1/4       {m4['ancho_img']:,.0f} × {m4['alto_img']:,.0f} px"
              f"   (diferencia {abs(m4['ancho_img']-m['ancho_img']):.0f}"
              f" × {abs(m4['alto_img']-m['alto_img']):.0f} px)")
        print(f"\n    lado del área de imagen          {lado_px:,.0f} px")
        print(f"    resolución en el terreno         {terreno:.4f} m/px  (del ajuste)")
        print(f"    MEDIA HUELLA DEL FOTOGRAMA       {media_huella:,.0f} m"
              f"   (lado/2 × resolución; no depende del formato del negativo)")
        print(f"\n    — expresado en unidades de película, suponiendo un formato de"
              f" {LADO_NEGATIVO_M*100:.0f} × {LADO_NEGATIVO_M*100:.0f} cm:")
        print(f"    tamaño de píxel en la película   {px_pelicula*1e6:.2f} µm"
              f"   (declarado en el XML: {PX_DECLARADO_XML*1e6:.0f} µm,"
              f" diferencia {100*(px_pelicula-PX_DECLARADO_XML)/PX_DECLARADO_XML:+.1f} %)")
        print(f"    escala efectiva                  1:{escala:,.0f}"
              f"   ({100*(escala-f['escala_nominal'])/f['escala_nominal']:+.1f} % frente a la nominal)")

    # --- imágenes de verificación
    print(f"\n{'=' * 88}")
    print("IMÁGENES DE VERIFICACIÓN")
    print("=" * 88)
    for fecha in ("1971", "1993"):
        v, p = imagen_verificacion(resultados[fecha]["m"])
        print(f"  {fecha}:  {v.name}  ·  {p.name}")
    print(f"  carpeta: {SALIDA.relative_to(RAIZ)}")

    # --- conjunto A
    print(f"\n{'=' * 88}")
    print("UMBRAL DEL CONJUNTO A CON LA ESCALA MEDIDA")
    print("=" * 88)
    mh93 = resultados["1993"]["media_huella"]
    recuento_conjunto_a(mh93, [
        ("media huella nominal 1:10,000", 1150.0),
        ("umbral declarado del conjunto A (08-07)", 1200.0),
        ("media huella con la escala MEDIDA", mh93),
        ("media huella si el píxel fuera de 18 µm", (LADO_NEGATIVO_M/2) * (resultados['1993']['terreno']/PX_DECLARADO_XML)),
        ("media huella si el lado de 23 cm fuera el ancho del archivo",
         (LADO_NEGATIVO_M/2) * (resultados['1993']['terreno'] / (LADO_NEGATIVO_M/FICHA['1993']['ancho_px']))),
    ])

    print("\n" + "=" * 88)
    print("FIN DEL SCRIPT — área de imagen medida en los dos fotogramas")
    print("=" * 88)


if __name__ == "__main__":
    main()
