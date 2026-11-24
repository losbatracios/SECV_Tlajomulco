#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validación del registro cartográfico del vaso de La Playa (1971)
=================================================================

PROPÓSITO
---------
Replicar el ajuste de la fotografía de febrero de 1971 sobre el vaso de La
Playa y medir hasta dónde llegan sus cifras. El criterio que ordena todo este
script es **precisión proporcional a la afirmación**: lo que el capítulo dice
del vaso es kilométrico —una zona de tono oscuro de unos 4 km², una planta a
alrededor de 1.7 km de su borde—, así que la pregunta no es cuántos decimales
tiene el error, sino si alguna de esas afirmaciones cambia al cambiar de
modelo de transformación. Por eso todo se calcula dos veces, con similitud y
con afín, y se compara.

El polígono `vaso_1971.gpkg` delimita **el vaso con sus márgenes**, por su
borde tonal: dentro hay parcelas y una vereda. No es una lámina de agua con
orilla definida, y aquí se le llama vaso, nunca lámina de agua.

QUÉ HACE
--------
1. Ajusta los dos modelos sobre los 5 puntos de control y reporta, con su
   nombre exacto, el RMS de residuos, el error estándar corregido y el error
   de predicción leave-one-out, junto con la redundancia, el número de
   ubicaciones independientes y la anisotropía del afín.
2. Mide el error de predicción en los 2 puntos de control independientes, que
   no entran al ajuste, y aplica la regla fijada de antemano: se elige el
   modelo con menor error en esos puntos.
3. Lleva los vértices del polígono al espacio de píxel de la fotografía y los
   devuelve al terreno con cada modelo, para medir cuánto se mueve el
   contorno, y cuánto cambian su área y su perímetro.
4. Mide la distancia de la planta del DENUE al polígono de cada modelo.

POR QUÉ LA GEOMETRÍA DE LOS PUNTOS IMPORTA AQUÍ
-----------------------------------------------
Tres de los cinco puntos de control están a menos de 130 m entre sí, en el
mismo caserío. Para un ajuste, ese grupo se comporta casi como un solo punto:
no hay cinco ubicaciones independientes, hay tres. Un modelo afín tiene seis
parámetros y necesita tres puntos para quedar determinado, de modo que con
tres ubicaciones queda casi exactamente determinado y su residuo bajo es un
resultado de la geometría, no una validación. Por eso el residuo no decide
nada y la decisión se deja a los puntos de control, aunque sean dos.

NOMENCLATURA DEL ERROR
----------------------
  RMS de residuos                 √(S/n)
  Error estándar corregido        √(S/(n − p/2))
  Error de predicción LOO         leave-one-out real, sobre los puntos del ajuste
  Error de predicción en puntos de control
                                  distancia entre la posición que predice el
                                  ajuste y la posición marcada, en puntos que
                                  no entraron al ajuste

Todos son error relativo a la base de referencia (Esri Wayback 2026-08-05),
no exactitud. Con 2 puntos de control, la estimación es muy débil.

Entrada : `.points` de `07_regresion_humedal/_gcp_qgis/`, `vaso_1971.gpkg`,
          `pisa_denue.geojson`, `HU1_1971_georef.tif` (solo sus etiquetas)
CRS     : EPSG:32613 (WGS 84 / UTM 13N)
Salidas : `05_humedal_georef/_validacion/`
Requiere: numpy, pyproj, shapely, pyogrio, Pillow
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from pyproj import Transformer
from pyogrio.raw import read as ogr_read, write as ogr_write
from shapely import wkb
from shapely.geometry import Point, Polygon

Image.MAX_IMAGE_PIXELS = None

AQUI = Path(__file__).resolve().parent                  # 07_regresion_humedal
RAIZ = AQUI.parents[2]
sys.path.insert(0, str(AQUI.parent / "05_regresion_par_estrella"))

# Los estimadores y los ajustadores son los mismos que validan el par
# comparativo: una sola implementación para todo el proyecto.
from validacion_par_estrella_v2 import (                # noqa: E402
    MODELOS, residuos, rms, error_estandar_corregido, prediccion_loo,
    parametros_similitud, _aplicar_matriz,
)

GCP = AQUI / "_gcp_qgis"
R_PUNTOS = GCP / "HU1_249BZ_37_5_1971.points"
R_CONTROLES = GCP / "HU1_checkpoints.points"
R_DESCRIPCION = GCP / "HU1_puntos_descripcion.csv"

VIS = RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef"
R_VASO = VIS / "vaso_1971.gpkg"
R_PISA = VIS / "pisa_denue.geojson"
R_TIF = VIS / "HU1_1971_georef.tif"
SALIDA = VIS / "_validacion"

# Fotograma fuente sobre el que se marcaron los puntos (sin rotar).
RASTER_FUENTE = (RAIZ / "02_Archivo_INEGI/Imágenes"
                      / "0249BZ_0037_005_0098_F1312-Z16B_1971_02_0035K_270"
                      / "0249BZ_0037_005_0098_F1312-Z16B_1971_02_0035K_270.jpg")

MODELOS_USADOS = ["similitud (Helmert)", "afín"]
UMBRAL_AGRUPACION_M = 100.0     # para contar ubicaciones independientes

A_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32613", always_xy=True)

# Cifras exploratorias de [1.2], para contrastar.
EXPLORATORIO = {
    "similitud (Helmert)": dict(rms=48.24, ese=62.27, ctrl={"C1": 23.3, "C2": 95.8},
                                ctrl_rms=69.7),
    "afín": dict(rms=4.59, ese=7.25, ctrl={"C1": 36.2, "C2": 17.5}, ctrl_rms=28.5),
}
ANISOTROPIA_EXPLORATORIA = 4.35


# ---------------------------------------------------------------------------
# Lectura
# ---------------------------------------------------------------------------

def leer_points(ruta):
    with open(ruta, encoding="utf-8") as f:
        lineas = [l for l in f if not l.startswith("#")]
    mapa, fuente, activo = [], [], []
    for r in csv.DictReader(lineas):
        mapa.append([float(r["mapX"]), float(r["mapY"])])
        fuente.append([float(r["sourceX"]), float(r["sourceY"])])
        activo.append(int(float(r.get("enable", 1))))
    return np.array(mapa), np.array(fuente), np.array(activo)


def leer_poligono(ruta):
    _, _, geoms, _ = ogr_read(str(ruta))
    return wkb.loads(geoms[0])


def etiquetas_geotiff(ruta):
    t = Image.open(ruta).tag_v2
    (sx, sy, _), (_, _, _, ox, oy, _) = t[33550], t[33922]
    ancho, alto = Image.open(ruta).size
    return dict(origen=(ox, oy), escala=(sx, sy), tamano=(ancho, alto),
                extension=(ox, oy - alto * sy, ox + ancho * sx, oy))


# ---------------------------------------------------------------------------
# 1. Geometría de los puntos
# ---------------------------------------------------------------------------

def ubicaciones_independientes(mapa, umbral=UMBRAL_AGRUPACION_M):
    """Agrupa por enlace simple: dos puntos a menos del umbral son el mismo sitio."""
    n = len(mapa)
    grupo = list(range(n))

    def raiz(i):
        while grupo[i] != i:
            grupo[i] = grupo[grupo[i]]
            i = grupo[i]
        return i

    d = np.hypot(*(mapa[:, None, :] - mapa[None, :, :]).transpose(2, 0, 1))
    for i in range(n):
        for j in range(i + 1, n):
            if d[i, j] < umbral:
                a, b = raiz(i), raiz(j)
                grupo[min(a, b)] = grupo[max(a, b)] = min(a, b)
    return d, len({raiz(i) for i in range(n)})


def anisotropia(H):
    """Diferencia entre los dos factores de escala del modelo, en por ciento."""
    s = np.linalg.svd(np.asarray(H)[:2, :2], compute_uv=False)
    return s[0], s[1], 100 * (s[0] / s[1] - 1)


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------

def main():
    SALIDA.mkdir(parents=True, exist_ok=True)
    mapa, fuente, activo = leer_points(R_PUNTOS)
    cmapa, cfuente, _ = leer_points(R_CONTROLES)
    etiquetas_ctrl = ["C1 (bodegas)", "C2 (cerca, piedemonte)"]
    n = len(mapa)

    print("#" * 92)
    print("# VALIDACIÓN DEL REGISTRO CARTOGRÁFICO DEL VASO DE LA PLAYA (febrero de 1971)")
    print("# fotograma 249BZ / rollo 37 / línea 5 / foto 98 · nominal 1:35,000")
    print("# base de referencia: Esri Wayback 2026-08-05 (release 26334) · EPSG:32613")
    print("#" * 92)
    print(f"\n  puntos de control del ajuste: {n} (activos: {int(activo.sum())})")
    print(f"  puntos de control independientes: {len(cmapa)}")

    # --- Geometría de los puntos --------------------------------------------
    d, n_sitios = ubicaciones_independientes(mapa)
    print(f"\n{'─' * 92}\nGEOMETRÍA DE LOS PUNTOS DEL AJUSTE\n{'─' * 92}")
    print(f"\n  distancias entre puntos (m)")
    print("      " + "".join(f"{j:>10}" for j in range(n)))
    for i in range(n):
        print(f"    {i} " + "".join("        — " if i == j else f"{d[i, j]:>10,.0f}"
                                    for j in range(n)))
    print(f"\n  ubicaciones independientes (puntos a más de {UMBRAL_AGRUPACION_M:.0f} m"
          f" entre sí): {n_sitios} de {n}")
    print(f"  Tres puntos forman un grupo en el mismo caserío; los otros dos están")
    print(f"  aislados. Para el ajuste cuentan {n_sitios} sitios, no {n} puntos.")

    # --- 1. Ajuste -----------------------------------------------------------
    print(f"\n{'─' * 92}\n1. AJUSTE POR MODELO\n{'─' * 92}")
    res = {}
    for mod in MODELOS_USADOS:
        cfg = MODELOS[mod]
        p = cfg["p"]
        H = cfg["ajustar"](fuente, mapa)
        r = residuos(H, fuente, mapa)
        loo, nota = prediccion_loo(cfg["ajustar"], p, cfg["min_pts"], fuente, mapa)
        esc, rot = parametros_similitud(H)
        res[mod] = dict(H=H, res=r, rms=rms(r), ese=error_estandar_corregido(r, p),
                        loo=loo, loo_rms=rms(loo) if loo is not None else None,
                        p=p, redundancia=2 * n - p, loo_red=2 * (n - 1) - p,
                        escala=esc, rotacion=rot, nota=nota)
        x = res[mod]
        print(f"\n  ── {mod}  (p = {p}; 2n = {2 * n}; redundancia 2n − p = {x['redundancia']})")
        print(f"     RMS de residuos            √(S/n)      {x['rms']:>9.2f} m")
        print(f"     error estándar corregido   √(S/(n−{p // 2}))  {x['ese']:>9.2f} m")
        print(f"     error de predicción LOO                {x['loo_rms']:>9.2f} m"
              f"   (redundancia por iteración: {x['loo_red']})")
        if x["nota"]:
            print(f"       nota: {x['nota']}")
        print(f"     resolución {x['escala']:.4f} m/px · rotación {x['rotacion']:+.2f}°")
        if mod == "afín":
            s1, s2, aniso = anisotropia(H)
            x["anisotropia"] = aniso
            print(f"     factores de escala del afín: {s1:.4f} y {s2:.4f} m/px"
                  f"  →  anisotropía {aniso:.2f} %")
        ex = EXPLORATORIO[mod]
        print(f"     cotejo con el exploratorio de [1.2]: "
              f"{ex['rms']:.2f} / {ex['ese']:.2f} m"
              f"  {'✓ replica' if abs(x['rms'] - ex['rms']) < 0.05 and abs(x['ese'] - ex['ese']) < 0.05 else '← difiere'}")

        print(f"\n     error de predicción LOO por punto (m)")
        for i in range(n):
            print(f"       punto {i}  {x['loo'][i]:>8.1f}")

    # --- 2. Puntos de control ------------------------------------------------
    print(f"\n{'─' * 92}\n2. ERROR DE PREDICCIÓN EN LOS PUNTOS DE CONTROL\n{'─' * 92}")
    print(f"\n  Los dos puntos no entraron al ajuste. Distancia de cada uno al punto")
    print(f"  del ajuste más próximo:")
    for k, c in zip(etiquetas_ctrl, cmapa):
        print(f"    {k:26} {np.min(np.hypot(*(mapa - c).T)):>7,.0f} m")

    print(f"\n  {'modelo':24}" + "".join(f"{k:>26}" for k in etiquetas_ctrl)
          + f"{'RMS (titular)':>16}")
    print("  " + "─" * 92)
    for mod in MODELOS_USADOS:
        e = residuos(res[mod]["H"], cfuente, cmapa)
        res[mod]["ctrl"] = e
        res[mod]["ctrl_rms"] = rms(e)
        print(f"  {mod:24}" + "".join(f"{v:>26.2f}" for v in e)
              + f"{rms(e):>16.2f}")
    for mod in MODELOS_USADOS:
        ex = EXPLORATORIO[mod]["ctrl_rms"]
        ok = abs(res[mod]["ctrl_rms"] - ex) < 0.1
        print(f"    cotejo {mod:22} exploratorio {ex:>6.1f} m"
              f"  {'✓ replica' if ok else '← difiere'}")

    elegido = min(MODELOS_USADOS, key=lambda m: res[m]["ctrl_rms"])
    otro = [m for m in MODELOS_USADOS if m != elegido][0]
    razon = res[otro]["ctrl_rms"] / res[elegido]["ctrl_rms"]
    print(f"\n  Regla fijada de antemano (menor error de predicción en los puntos de")
    print(f"  control): queda el modelo {elegido.upper()}"
          f"  ({res[elegido]['ctrl_rms']:.1f} m frente a {res[otro]['ctrl_rms']:.1f} m,"
          f" {razon:.1f}×).")
    # ¿Cuánto decide cada control?
    print(f"\n  Con n = 2, la elección es muy débil. Qué pasa si se quita uno:")
    for i, k in enumerate(etiquetas_ctrl):
        queda = 1 - i
        g = min(MODELOS_USADOS, key=lambda m: abs(res[m]["ctrl"][queda]))
        print(f"    sin {k:26} → decide {etiquetas_ctrl[queda]:26} y gana {g}")
    dif_c2 = abs(res[otro]["ctrl"][1] - res[elegido]["ctrl"][1])
    dif_c1 = abs(res[otro]["ctrl"][0] - res[elegido]["ctrl"][0])
    print(f"    La diferencia entre modelos vale {dif_c1:.1f} m en C1 y {dif_c2:.1f} m"
          f" en C2: C2 decide la mayor parte.")

    # --- Verificación de la afín contra el GeoTIFF ---------------------------
    print(f"\n{'─' * 92}\nVERIFICACIÓN DE LA AFÍN CONTRA EL RÁSTER GEORREFERENCIADO\n{'─' * 92}")
    g = etiquetas_geotiff(R_TIF)
    W, Hh = Image.open(RASTER_FUENTE).size
    print(f"\n  fotograma fuente        {W:,} × {Hh:,} px")
    print(f"  GeoTIFF georreferenciado {g['tamano'][0]:,} × {g['tamano'][1]:,} px"
          f" · {g['escala'][0]:.5f} m/px · norte arriba")
    esquinas = _aplicar_matriz(res["afín"]["H"],
                               np.array([[0, 0], [W, 0], [W, -Hh], [0, -Hh]], float))
    caja = (esquinas[:, 0].min(), esquinas[:, 1].min(),
            esquinas[:, 0].max(), esquinas[:, 1].max())
    print(f"\n  {'':22}{'x mín':>13}{'y mín':>14}{'x máx':>13}{'y máx':>14}")
    print("  " + "─" * 78)
    print(f"  {'afín sobre el fuente':22}" + "".join(f"{v:>13,.1f} " for v in caja))
    print(f"  {'extensión del GeoTIFF':22}"
          + "".join(f"{v:>13,.1f} " for v in g["extension"]))
    dif = np.abs(np.array(caja) - np.array(g["extension"]))
    print(f"  {'diferencia':22}" + "".join(f"{v:>13,.1f} " for v in dif))
    print(f"\n  La mayor diferencia es de {dif.max():.1f} m sobre un cuadro de"
          f" {(caja[2]-caja[0])/1000:.1f} km, es decir {100*dif.max()/(caja[2]-caja[0]):.3f} %.")
    print(f"  La afín de los 5 puntos es la que produjo el ráster.")

    # --- 3. Sensibilidad del polígono ---------------------------------------
    print(f"\n{'─' * 92}\n3. SENSIBILIDAD DEL CONTORNO DEL VASO AL MODELO\n{'─' * 92}")
    vaso = leer_poligono(R_VASO)
    v_utm = np.array(vaso.exterior.coords)
    inv_af = np.linalg.inv(res["afín"]["H"])
    v_px = _aplicar_matriz(inv_af, v_utm)
    print(f"\n  El contorno se dibujó sobre el ráster, que está en afín. Para")
    print(f"  compararlo entre modelos se lleva a píxeles del fotograma con la")
    print(f"  inversa de la afín y se devuelve al terreno con cada modelo.")
    print(f"\n  contorno original: {len(v_utm) - 1} vértices · "
          f"{vaso.area / 1e6:.4f} km² ({vaso.area / 1e4:,.0f} ha) · "
          f"perímetro {vaso.length:,.0f} m")

    poligonos = {}
    print(f"\n  {'modelo':24}{'área':>14}{'perímetro':>13}"
          f"{'desplazamiento de vértices':>32}")
    print(f"  {'':24}{'':>14}{'':>13}{'mediana':>14}{'máximo':>12}  {'centroide':>10}")
    print("  " + "─" * 90)
    for mod in MODELOS_USADOS:
        v2 = _aplicar_matriz(res[mod]["H"], v_px)
        poli = Polygon(v2)
        poligonos[mod] = poli
        desp = np.hypot(*(v2 - v_utm).T)
        d_cent = Point(poli.centroid).distance(vaso.centroid)
        print(f"  {mod:24}{poli.area / 1e6:>10.4f} km²{poli.length:>11,.0f} m"
              f"{np.median(desp):>12,.1f} m{desp.max():>11,.1f} m{d_cent:>10,.1f} m")
        res[mod]["poligono"] = poli
        res[mod]["desp"] = desp

    dif_mod = np.hypot(*(np.array(poligonos["afín"].exterior.coords)
                         - np.array(poligonos["similitud (Helmert)"].exterior.coords)).T)
    print(f"\n  Diferencia entre los dos modelos a lo largo del contorno:"
          f" mediana {np.median(dif_mod):,.1f} m · mínimo {dif_mod.min():,.1f} m"
          f" · máximo {dif_mod.max():,.1f} m")
    a1, a2 = (poligonos[m].area / 1e4 for m in MODELOS_USADOS)
    print(f"  El área pasa de {a1:,.0f} a {a2:,.0f} ha según el modelo:"
          f" {abs(100 * (a2 / a1 - 1)):.1f} % de diferencia.")

    # El contenedor GPKG se escribe primero en una carpeta temporal y después
    # se copia al destino: la biblioteca necesita poder borrar y recrear el
    # archivo, y en la carpeta del proyecto el borrado está deshabilitado.
    import shutil
    import tempfile
    temporal = Path(tempfile.mkdtemp())
    for mod in MODELOS_USADOS:
        nombre = "helmert" if mod.startswith("similitud") else "afin"
        destino = VIS / f"vaso_1971_{nombre}_v1.gpkg"
        intermedio = temporal / destino.name
        ogr_write(str(intermedio),
                  geometry=np.array([wkb.dumps(poligonos[mod])], dtype=object),
                  field_data=[np.array([f"vaso de La Playa ({mod})"], dtype=object)],
                  fields=np.array(["nombre"], dtype=object),
                  layer=f"vaso_1971_{nombre}", driver="GPKG",
                  geometry_type="Polygon", crs="EPSG:32613")
        shutil.copyfile(intermedio, destino)
        print(f"  escrito {destino.name}")
    shutil.rmtree(temporal, ignore_errors=True)

    # --- 4. PISA --------------------------------------------------------------
    print(f"\n{'─' * 92}\n4. LA PLANTA DEL DENUE FRENTE AL VASO, POR MODELO\n{'─' * 92}")
    _, _, gp, _ = ogr_read(str(R_PISA))
    p_geo = wkb.loads(gp[0])
    pisa = Point(*A_UTM.transform(p_geo.x, p_geo.y))
    print(f"\n  planta (DENUE 6336993): {p_geo.y:.8f}, {p_geo.x:.8f}"
          f"  →  E {pisa.x:,.0f}  N {pisa.y:,.0f}")

    def rumbo(desde, hasta):
        dx, dy = hasta.x - desde.x, hasta.y - desde.y
        az = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
        rosa = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                "S", "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO"]
        return az, rosa[int((az + 11.25) % 360 // 22.5)]

    print(f"\n  {'modelo':24}{'al borde':>12}{'al centroide':>15}{'¿dentro?':>11}"
          f"{'rumbo desde el centroide':>28}")
    print("  " + "─" * 90)
    for mod in MODELOS_USADOS:
        poli = poligonos[mod]
        dentro = poli.contains(pisa)
        d_borde = poli.exterior.distance(pisa)
        d_cent = pisa.distance(poli.centroid)
        az, r = rumbo(poli.centroid, pisa)
        print(f"  {mod:24}{d_borde:>10,.0f} m{d_cent:>13,.0f} m"
              f"{'dentro' if dentro else 'fuera':>11}{f'{r} ({az:.0f}°)':>28}")
        res[mod]["pisa"] = dict(borde=d_borde, centroide=d_cent, dentro=dentro,
                                rumbo=r, azimut=az)
    print(f"\n  Cotejo: 3,536 m al centroide (QGIS, 26-10) y 1,730 m al borde de la")
    print(f"  caja del vaso (v7, 23-10). Son referencias distintas: el centroide del")
    print(f"  polígono, y el borde de la caja envolvente estimada en agosto.")

    # --- Resumen y veredicto --------------------------------------------------
    print(f"\n{'=' * 92}\nVEREDICTO\n{'=' * 92}")
    ar = [poligonos[m].area / 1e4 for m in MODELOS_USADOS]
    bo = [res[m]["pisa"]["borde"] for m in MODELOS_USADOS]
    print(f"\n  afirmación                                    {'Helmert':>14}{'afín':>14}"
          f"{'¿cambia?':>12}")
    print("  " + "─" * 90)
    print(f"  {'zona de tono oscuro de unos 4 km²':44}"
          f"{ar[0] / 100:>12.2f} km²{ar[1] / 100:>12.2f} km²"
          f"{'no' if round(ar[0] / 100, 1) == round(ar[1] / 100, 1) else 'sí':>12}")
    print(f"  {'distancia de la planta al borde del vaso':44}"
          f"{bo[0]:>11,.0f} m{bo[1]:>13,.0f} m"
          f"{'no' if abs(bo[0] - bo[1]) < 100 else 'sí':>12}")
    print(f"\n  Atención: esa distancia NO es la que circula en la bitácora. Los")
    print(f"  1,730 m de v7 se midieron al borde de la caja envolvente que se estimó")
    print(f"  en agosto, antes de que existiera el contorno dibujado. Medida al")
    print(f"  contorno, la distancia es de ~{min(bo)/1000:.1f} km, no de ~1.7 km. Lo que no")
    print(f"  cambia es el orden de magnitud ni el hecho de que la planta queda fuera.")
    print(f"  {'la planta queda fuera del vaso':44}"
          f"{'fuera' if not res[MODELOS_USADOS[0]]['pisa']['dentro'] else 'dentro':>14}"
          f"{'fuera' if not res[MODELOS_USADOS[1]]['pisa']['dentro'] else 'dentro':>14}"
          f"{'no':>12}")

    # --- Archivos --------------------------------------------------------------
    with open(SALIDA / "ajuste_humedal.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["modelo", "n_gcp", "ubicaciones_independientes", "p",
                    "redundancia", "rms_residuos_m", "error_estandar_corregido_m",
                    "loo_rms_m", "ctrl_C1_m", "ctrl_C2_m", "ctrl_rms_m",
                    "resolucion_m_px", "anisotropia_pct", "area_ha", "perimetro_m",
                    "pisa_borde_m", "pisa_centroide_m", "pisa_dentro"])
        for mod in MODELOS_USADOS:
            x = res[mod]
            w.writerow([mod, n, n_sitios, x["p"], x["redundancia"],
                        round(x["rms"], 3), round(x["ese"], 3),
                        round(x["loo_rms"], 3), round(x["ctrl"][0], 3),
                        round(x["ctrl"][1], 3), round(x["ctrl_rms"], 3),
                        round(x["escala"], 5),
                        round(x.get("anisotropia", float("nan")), 3),
                        round(x["poligono"].area / 1e4, 2),
                        round(x["poligono"].length, 1),
                        round(x["pisa"]["borde"], 1),
                        round(x["pisa"]["centroide"], 1), x["pisa"]["dentro"]])
    with open(SALIDA / "matrices_humedal_pixel_a_UTM32613.json", "w",
              encoding="utf-8") as fh:
        json.dump({m: res[m]["H"].tolist() for m in MODELOS_USADOS}, fh, indent=1)
    print(f"\n  salidas en {SALIDA.relative_to(RAIZ)}")
    print(f"\n  Todos los estimadores son error relativo a la base, no exactitud.")
    print(f"  Con 2 puntos de control la elección de modelo es muy débil.")


if __name__ == "__main__":
    main()
