#!/usr/bin/env python3
"""
mss_contraste_v1.py — Landsat MSS sobre el vaso de la Presa La Playa.
Proyecto: Tlajomulco · análisis reproducible · creado 2026-11-10

QUÉ HACE
  A. Cruza mss_clasificacion_plantilla.csv con mss_inventory.csv (revela las
     fechas) y tabula las escenas por año y por mes.
  B. Audita la estructura del material: cuántas imágenes distintas existen
     frente a cuántas filas tiene la plantilla.
  C. Mide el contraste del infrarrojo cercano en cada chip disponible, dentro
     del contorno vaso_1971.gpkg frente a un anillo exterior, con el mismo
     procedimiento en todos.
  D. Reporta la distribución de esa medida.
  E. Compara con la serie JRC en los meses en que ambos productos coinciden.

QUÉ MATERIAL HAY, Y QUÉ LIMITA LA MEDIDA
  No existe un ráster por escena. El único material con valores por píxel son
  12 chips anuales: cada uno es la mediana de las escenas de su año y está
  escrito en 8 bits ya estirados entre el percentil 2 y el 98 de la mediana de
  toda la serie (mss_vaso_v1.js, bloque 4). El estiramiento es el mismo para
  los doce, así que los valores son comparables entre años, pero está censurado
  en los dos extremos: todo lo más oscuro que el percentil 2 se escribe como 0
  y todo lo más claro que el 98 como 255. Por eso una mediana de nivel digital
  no es una medida utilizable: en varios años la mediana del contorno o la del
  anillo cae en un extremo censurado. Este script la calcula y la reporta, pero
  marcada como censurada, y la medida principal es otra.

MEDIDA PRINCIPAL  (definida por el estiramiento, no elegida aquí)
  Dos medidas, una absoluta y una relativa, porque la censura afecta a cada una
  de distinta manera y ninguna basta sola.

  (1) Absoluta — fracción de píxeles en el extremo oscuro, es decir con nivel
      digital 0, o sea con mediana anual del infrarrojo cercano por debajo del
      percentil 2 de toda la serie. Se calcula dentro del contorno y en el
      anillo, y el contraste es la diferencia entre las dos fracciones, en
      puntos porcentuales. El umbral lo fija el estiramiento del producto,
      declarado antes de ver cualquier imagen; no se ajusta aquí. Su límite:
      en un año de mosaico claro, el vaso puede ser más oscuro que su entorno
      sin llegar a ese extremo, y la medida lo registra como cero.

  (2) Relativa — probabilidad de que un píxel del contorno, tomado al azar, sea
      más oscuro que un píxel del anillo tomado al azar, con los empates
      contados a la mitad. Se calcula con los rangos, así que no cambia si a
      todos los píxeles de la imagen se les aplica la misma transformación
      monótona, recorte de extremos incluido: por eso sobrevive a la censura.
      Vale 0.5 cuando el contorno y el anillo son indistinguibles y 1 cuando
      todo el contorno es más oscuro que todo el anillo. Su límite: dice que
      hay diferencia, no cuánta, y no distingue agua de suelo húmedo o de
      vegetación densa.

REGLAS FIJADAS ANTES DE CALCULAR
  - El anillo de comparación es la corona entre el contorno y su dilatación de
    500 m, calculada en EPSG:32613.
  - El valor 0 también marca los píxeles sin dato de las esquinas del recorte,
    así que el script verifica primero que el contorno dilatado cae por
    completo dentro de la región exportada.
  - El área oscura se compara con el criterio visual declarado en el README
    (parche de ~1 km² dentro de la caja del vaso).

USO
    cd /ruta/a/Tlajomulco
    python3 04_Análisis/Scripts_Python/08_mss/mss_contraste_v1.py

SALIDAS  (04_Análisis/Visualizaciones/06_mss/)
    mss_cruce_escenas_fechas.csv
    mss_escenas_anio_mes.csv
    mss_contraste_anual.csv
    mss_jrc_meses_comunes.csv
"""

import pathlib
import sys

import numpy as np
import pandas as pd
from PIL import Image
from matplotlib.path import Path as MplPath
from pyogrio.raw import read as ogr_read
from pyproj import Transformer
from shapely import wkb as shapely_wkb
from shapely.geometry import box

Image.MAX_IMAGE_PIXELS = None

# ─── Rutas ────────────────────────────────────────────────────────────────────
RAIZ = pathlib.Path(__file__).resolve().parents[3]

INVENTARIO = RAIZ / "05_Datos_abiertos/Landsat_MSS/mss_inventory.csv"
PLANTILLA = RAIZ / "05_Datos_abiertos/Landsat_MSS/mss_clasificacion_plantilla.csv"
CHIPS = RAIZ / "05_Datos_abiertos/Landsat_MSS/chips"
VASO = RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef/vaso_1971.gpkg"
JRC = RAIZ / "04_Análisis/Visualizaciones/05_humedal_georef/_validacion/jrc_replica_mensual.csv"
SALIDA = RAIZ / "04_Análisis/Visualizaciones/06_mss"

# ─── Constantes del procedimiento ─────────────────────────────────────────────
CRS_VASO = "EPSG:32613"
CRS_CHIP = "EPSG:4326"
ANILLO_M = 500.0
ND_EXTREMO_OSCURO = 0                # nivel digital del extremo oscuro censurado
ND_EXTREMO_CLARO = 255               # nivel digital del extremo claro censurado
BUFFER_GEE_M = 3000.0                # buffer del recorte en mss_vaso_v1.js
LADO_PX_M = 60.0                     # resolución nominal del producto
CAJA_README = dict(sur=20.4647, norte=20.4796, oeste=-103.5442, este=-103.5100)
CRITERIO_PARCHE_KM2 = 1.0            # criterio visual declarado en el README

# Meses en que los dos productos pueden coincidir: la serie mensual JRC
# empieza en marzo de 1984 y el inventario MSS termina en diciembre de 1985.
VENTANA_COMUN = ((1984, 3), (1985, 12))


# ══════════════════════════════════════════════════════════════════════════════
# A. Cruce de la plantilla con el inventario
# ══════════════════════════════════════════════════════════════════════════════
def bloque_cruce():
    print("=" * 78)
    print("A. CRUCE DE LA CLASIFICACIÓN CIEGA CON EL INVENTARIO")
    print("=" * 78)

    inv = pd.read_csv(INVENTARIO)
    pla = pd.read_csv(PLANTILLA)

    # scene_code = SATELITE_TIER_YYYYMMDD, tal como lo construye mss_contacto_v1.py
    inv["scene_code"] = (inv["satellite"] + "_" + inv["tier"] + "_"
                         + inv["date"].str.replace("-", "", regex=False))

    print(f"  inventario: {len(inv)} escenas")
    print(f"  plantilla:  {len(pla)} filas")
    print(f"  scene_code duplicados en el inventario: {int(inv['scene_code'].duplicated().sum())}")

    d = pla.merge(
        inv[["scene_code", "date", "year", "month", "satellite", "tier",
             "cloud_cover_pct", "nir_band"]],
        on="scene_code", how="outer", indicator=True)
    conteo = d["_merge"].value_counts().to_dict()
    print(f"  cruce: {conteo}")
    if conteo.get("both", 0) != len(pla):
        print("  ATENCIÓN: el cruce no es completo.")
    d = d.drop(columns="_merge")
    d["clasificacion"] = d["clasificacion"].fillna("(sin clasificar)")
    d["nota"] = d["nota"].fillna("")
    d = d.sort_values("date").reset_index(drop=True)

    print("\n  clasificación revelada:")
    for k, v in d["clasificacion"].value_counts().items():
        print(f"    {k:22s} {v:3d}")
    print("  notas:")
    for k, v in d["nota"].replace("", "(vacía)").value_counts().items():
        print(f"    {v:3d}  {k}")

    # Tabla año × mes
    tabla = (pd.crosstab(d["year"], d["month"])
               .reindex(columns=range(1, 13), fill_value=0))
    tabla.columns = [f"m{c:02d}" for c in tabla.columns]
    tabla["total"] = tabla.sum(axis=1)
    print("\n  escenas por año y mes:")
    print(tabla.to_string())
    print(f"\n  años con escena: {sorted(d['year'].unique())}")
    print(f"  años sin cobertura en 1972–1985: "
          f"{sorted(set(range(1972, 1986)) - set(d['year'].unique()))}")
    print(f"  rango de fechas: {d['date'].min()} a {d['date'].max()}")

    return d, tabla


# ══════════════════════════════════════════════════════════════════════════════
# B. Auditoría de la estructura del material
# ══════════════════════════════════════════════════════════════════════════════
def bloque_estructura(d):
    print()
    print("=" * 78)
    print("B. ESTRUCTURA DEL MATERIAL CLASIFICADO")
    print("=" * 78)

    chips = {int(p.stem.split("_")[-1]): p for p in sorted(CHIPS.glob("mss_chip_*.tif"))}
    print(f"  chips en disco: {len(chips)} → años {sorted(chips)}")
    print("  cada chip es la mediana de todas las escenas de su año (mss_vaso_v1.js,")
    print("  bloque 4), de modo que a cada escena de un año le corresponde la misma imagen.")

    d = d.copy()
    d["chip_year"] = d["year"]
    d["chip_existe"] = d["chip_year"].isin(chips)
    print(f"\n  filas de la plantilla: {len(d)}")
    print(f"  imágenes distintas disponibles: {d.loc[d['chip_existe'], 'chip_year'].nunique()}")
    print(f"  filas sin chip propio: {int((~d['chip_existe']).sum())}")

    print("\n  repeticiones de la misma imagen en la plantilla:")
    rep = d.groupby("chip_year").size()
    for y, n in rep.items():
        print(f"    {y}: {n:2d} filas con la misma imagen")

    print("\n  coherencia de la nota dentro de cada año (la imagen es una sola):")
    incoherentes = []
    for y, g in d.groupby("chip_year"):
        u = g["nota"].replace("", "(vacía)").value_counts()
        marca = "" if len(u) == 1 else "  ← la misma imagen con notas distintas"
        if len(u) > 1:
            incoherentes.append(y)
        print(f"    {y}  n={len(g):2d}  notas distintas={len(u)}{marca}")
    print(f"\n  años con notas contradictorias sobre la misma imagen: {incoherentes}")

    return d, chips


# ══════════════════════════════════════════════════════════════════════════════
# C. Contraste del infrarrojo cercano por chip
# ══════════════════════════════════════════════════════════════════════════════
def leer_chip(path):
    """Devuelve (banda, lon_centros, lat_centros) del chip GeoTIFF."""
    im = Image.open(path)
    a = np.array(im)
    if a.ndim == 3:
        # visualize() escribe tres bandas idénticas en escala de grises
        if not (np.array_equal(a[:, :, 0], a[:, :, 1])
                and np.array_equal(a[:, :, 0], a[:, :, 2])):
            raise ValueError(f"{path.name}: las tres bandas no son idénticas")
        banda = a[:, :, 0]
    else:
        banda = a
    tp = im.tag_v2[33922]        # ModelTiepoint
    ps = im.tag_v2[33550]        # ModelPixelScale
    lon0, lat0 = float(tp[3]), float(tp[4])
    dlon, dlat = float(ps[0]), float(ps[1])
    alto, ancho = banda.shape
    lon = lon0 + (np.arange(ancho) + 0.5) * dlon
    lat = lat0 - (np.arange(alto) + 0.5) * dlat
    return banda.astype(float), lon, lat


def prob_mas_oscuro(v_dentro, v_anillo):
    """P(píxel del contorno más oscuro que píxel del anillo), empates a la mitad.

    Se calcula con rangos promediados, de modo que es invariante a cualquier
    transformación monótona común a las dos muestras, incluido el recorte de
    extremos del estiramiento de 8 bits.
    """
    v = np.concatenate([v_dentro, v_anillo])
    orden = np.argsort(v, kind="mergesort")
    v_ord = v[orden]
    rangos_ord = np.empty(len(v), dtype=float)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v_ord[j + 1] == v_ord[i]:
            j += 1
        rangos_ord[i:j + 1] = (i + j) / 2 + 1     # rango promedio, base 1
        i = j + 1
    rangos = np.empty(len(v), dtype=float)
    rangos[orden] = rangos_ord
    n1, n2 = len(v_dentro), len(v_anillo)
    r1 = rangos[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2        # pares en que el contorno es MÁS CLARO
    return 1 - u1 / (n1 * n2)          # se devuelve el complemento: más oscuro


def comprobar_direccion():
    """Verifica el sentido de prob_mas_oscuro con dos muestras conocidas."""
    oscuro = np.zeros(10)
    claro = np.full(10, 200.0)
    p1 = prob_mas_oscuro(oscuro, claro)
    p2 = prob_mas_oscuro(claro, oscuro)
    p3 = prob_mas_oscuro(claro, claro)
    assert abs(p1 - 1.0) < 1e-9, p1
    assert abs(p2 - 0.0) < 1e-9, p2
    assert abs(p3 - 0.5) < 1e-9, p3
    return p1, p2, p3


def area_px_ha(lon, lat):
    """Área real del píxel en hectáreas, medida sobre la retícula proyectada."""
    hacia_utm = Transformer.from_crs(CRS_CHIP, CRS_VASO, always_xy=True)
    i, j = len(lon) // 2, len(lat) // 2
    e0, n0 = hacia_utm.transform(lon[i], lat[j])
    e1, _ = hacia_utm.transform(lon[i + 1], lat[j])
    _, n1 = hacia_utm.transform(lon[i], lat[j + 1])
    return abs(e1 - e0) * abs(n1 - n0) / 1e4, abs(e1 - e0), abs(n1 - n0)


def region_recorte(caja):
    """Región exportada por mss_vaso_v1.js: la caja dilatada 3 km, en UTM."""
    hacia_utm = Transformer.from_crs(CRS_CHIP, CRS_VASO, always_xy=True)
    e0, n0 = hacia_utm.transform(caja["oeste"], caja["sur"])
    e1, n1 = hacia_utm.transform(caja["este"], caja["norte"])
    return box(min(e0, e1), min(n0, n1), max(e0, e1), max(n0, n1)).buffer(BUFFER_GEE_M)


def mascaras(lon, lat, poli, caja):
    """Máscaras booleanas (contorno, anillo, caja) sobre la retícula del chip."""
    LON, LAT = np.meshgrid(lon, lat)
    hacia_utm = Transformer.from_crs(CRS_CHIP, CRS_VASO, always_xy=True)
    E, N = hacia_utm.transform(LON.ravel(), LAT.ravel())
    pts = np.column_stack([E, N])

    dentro = MplPath(np.asarray(poli.exterior.coords)).contains_points(pts)
    dilatado = poli.buffer(ANILLO_M)
    en_dilatado = MplPath(np.asarray(dilatado.exterior.coords)).contains_points(pts)
    anillo = en_dilatado & ~dentro

    en_caja = ((LON.ravel() >= caja["oeste"]) & (LON.ravel() <= caja["este"])
               & (LAT.ravel() >= caja["sur"]) & (LAT.ravel() <= caja["norte"]))

    forma = LON.shape
    return (dentro.reshape(forma), anillo.reshape(forma), en_caja.reshape(forma))


def bloque_contraste(chips):
    print()
    print("=" * 78)
    print("C. CONTRASTE DEL INFRARROJO CERCANO POR CHIP")
    print("=" * 78)

    p1, p2, p3 = comprobar_direccion()
    print(f"  comprobación de la medida relativa: todo oscuro frente a todo claro "
          f"= {p1:.1f}; al revés = {p2:.1f}; iguales = {p3:.1f}")

    meta, _, geom, _ = ogr_read(str(VASO))
    if meta["crs"] != CRS_VASO:
        raise ValueError(f"vaso_1971.gpkg en {meta['crs']}, se esperaba {CRS_VASO}")
    poli = shapely_wkb.loads(geom[0])
    dilatado = poli.buffer(ANILLO_M)
    print(f"  contorno: {poli.geom_type}, {poli.area / 1e4:.1f} ha, CRS {meta['crs']}")
    print(f"  anillo: corona entre el contorno y su dilatación de {ANILLO_M:.0f} m")

    # Verificación previa: el anillo no puede tocar las esquinas sin dato
    region = region_recorte(CAJA_README)
    contenido = region.contains(dilatado)
    print(f"  el contorno dilatado cae dentro de la región exportada: {contenido}")
    print(f"    holgura al borde de la región: "
          f"{dilatado.distance(region.exterior):.0f} m")
    if not contenido:
        print("    ATENCIÓN: parte del anillo podría caer en píxeles sin dato.")

    filas = []
    primera = True
    for anio in sorted(chips):
        banda, lon, lat = leer_chip(chips[anio])
        if primera:
            px_ha, px_e, px_n = area_px_ha(lon, lat)
            print(f"\n  retícula del chip: {banda.shape[1]} × {banda.shape[0]} px, "
                  f"paso {abs(lon[1] - lon[0]):.6f}°")
            print(f"  píxel proyectado: {px_e:.1f} × {px_n:.1f} m = {px_ha:.4f} ha "
                  f"(nominal del producto: {LADO_PX_M:.0f} m)")
            m_dentro, m_anillo, m_caja = mascaras(lon, lat, poli, CAJA_README)
            print(f"  píxeles: contorno {int(m_dentro.sum())}, anillo "
                  f"{int(m_anillo.sum())}, caja del README {int(m_caja.sum())}")
            print(f"  área del contorno en píxeles: "
                  f"{int(m_dentro.sum()) * px_ha:.1f} ha "
                  f"(contorno vectorial: {poli.area / 1e4:.1f} ha)")
            primera = False

        v_in, v_an, v_caja = banda[m_dentro], banda[m_anillo], banda[m_caja]

        n0_in = int((v_in == ND_EXTREMO_OSCURO).sum())
        n0_an = int((v_an == ND_EXTREMO_OSCURO).sum())
        n0_caja = int((v_caja == ND_EXTREMO_OSCURO).sum())
        f0_in, f0_an = n0_in / len(v_in), n0_an / len(v_an)

        n255_in = int((v_in == ND_EXTREMO_CLARO).sum())
        n255_an = int((v_an == ND_EXTREMO_CLARO).sum())

        p_oscuro = prob_mas_oscuro(v_in, v_an)
        med_in, med_an = float(np.median(v_in)), float(np.median(v_an))
        censurada = (med_in in (ND_EXTREMO_OSCURO, ND_EXTREMO_CLARO)
                     or med_an in (ND_EXTREMO_OSCURO, ND_EXTREMO_CLARO))

        area_oscura_caja = n0_caja * px_ha / 100

        filas.append(dict(
            anio=anio,
            px_contorno=int(m_dentro.sum()),
            px_anillo=int(m_anillo.sum()),
            frac_extremo_oscuro_contorno=round(f0_in * 100, 1),
            frac_extremo_oscuro_anillo=round(f0_an * 100, 1),
            contraste_pp=round((f0_in - f0_an) * 100, 1),
            prob_contorno_mas_oscuro=round(p_oscuro, 3),
            area_oscura_contorno_ha=round(n0_in * px_ha, 1),
            area_oscura_caja_km2=round(area_oscura_caja, 2),
            cumple_criterio_1km2=bool(area_oscura_caja >= CRITERIO_PARCHE_KM2),
            frac_extremo_claro_contorno=round(n255_in / len(v_in) * 100, 1),
            frac_extremo_claro_anillo=round(n255_an / len(v_an) * 100, 1),
            nd_mediana_contorno=med_in,
            nd_mediana_anillo=med_an,
            mediana_censurada=censurada,
            chip=chips[anio].name,
        ))

    df = pd.DataFrame(filas)
    print()
    cols = ["anio", "frac_extremo_oscuro_contorno", "frac_extremo_oscuro_anillo",
            "contraste_pp", "area_oscura_contorno_ha", "area_oscura_caja_km2",
            "cumple_criterio_1km2", "prob_contorno_mas_oscuro"]
    print("  medidas por chip — absoluta (extremo oscuro) y relativa (rangos):")
    print(df[cols].to_string(index=False))
    print()
    print("  censura del estiramiento (porcentaje de píxeles en cada extremo):")
    print(df[["anio", "frac_extremo_oscuro_contorno", "frac_extremo_claro_contorno",
              "frac_extremo_oscuro_anillo", "frac_extremo_claro_anillo",
              "nd_mediana_contorno", "nd_mediana_anillo",
              "mediana_censurada"]].to_string(index=False))
    print(f"\n  chips con alguna mediana en un extremo censurado: "
          f"{int(df['mediana_censurada'].sum())} de {len(df)} "
          f"→ la mediana de nivel digital no se usa como cifra")
    return df


# ══════════════════════════════════════════════════════════════════════════════
# D. Distribución de la medida
# ══════════════════════════════════════════════════════════════════════════════
def bloque_distribucion(df, d):
    print()
    print("=" * 78)
    print("D. DISTRIBUCIÓN DE LA MEDIDA")
    print("=" * 78)

    c = df["contraste_pp"].astype(float)
    print(f"  n = {len(c)} chips anuales (no 132: hay una imagen por año)")
    print(f"  contraste en puntos porcentuales de píxeles en el extremo oscuro:")
    print(f"    mínimo    {c.min():6.1f}  (año {int(df.loc[c.idxmin(), 'anio'])})")
    print(f"    cuartil 1 {c.quantile(.25):6.1f}")
    print(f"    mediana   {c.median():6.1f}")
    print(f"    cuartil 3 {c.quantile(.75):6.1f}")
    print(f"    máximo    {c.max():6.1f}  (año {int(df.loc[c.idxmax(), 'anio'])})")
    print(f"  chips con contraste ≤ 0 (el vaso no es más oscuro que su anillo): "
          f"{sorted(df.loc[c <= 0, 'anio'].tolist())}")

    q = df["prob_contorno_mas_oscuro"].astype(float)
    print(f"\n  probabilidad de que el contorno sea más oscuro que el anillo:")
    print(f"    mínimo    {q.min():6.3f}  (año {int(df.loc[q.idxmin(), 'anio'])})")
    print(f"    mediana   {q.median():6.3f}")
    print(f"    máximo    {q.max():6.3f}  (año {int(df.loc[q.idxmax(), 'anio'])})")
    print(f"  chips con probabilidad ≤ 0.5 (contorno no más oscuro): "
          f"{sorted(df.loc[q <= 0.5, 'anio'].tolist())}")
    print(f"  chips con probabilidad ≥ 0.9: "
          f"{sorted(df.loc[q >= 0.9, 'anio'].tolist())}")

    a = df["area_oscura_caja_km2"].astype(float)
    print(f"\n  área en el extremo oscuro dentro de la caja del README (km²): "
          f"mín {a.min():.2f}, mediana {a.median():.2f}, máx {a.max():.2f}")
    print(f"  chips que no alcanzan el parche de {CRITERIO_PARCHE_KM2:.0f} km²: "
          f"{sorted(df.loc[~df['cumple_criterio_1km2'], 'anio'].tolist())}")

    pasan = df.loc[df["cumple_criterio_1km2"], "anio"].tolist()
    fallan = df.loc[~df["cumple_criterio_1km2"], "anio"].tolist()
    n_pasan = int(d["year"].isin(pasan).sum())
    n_fallan = int(d["year"].isin(fallan).sum())
    print(f"\n  si el criterio de ~1 km² se aplica con esta medida sobre la imagen "
          f"de cada año:")
    print(f"    años que lo alcanzan: {len(pasan)} de {len(df)} → {n_pasan} filas "
          f"de la plantilla")
    print(f"    años que no lo alcanzan: {len(fallan)} de {len(df)} → {n_fallan} "
          f"filas de la plantilla")
    print(f"    frente a las {len(d)} filas clasificadas como agua visible")

    print("\n  orden de los años por contraste (de menor a mayor):")
    for _, r in df.sort_values("contraste_pp").iterrows():
        n_esc = int((d["year"] == r["anio"]).sum())
        print(f"    {int(r['anio'])}  contraste {r['contraste_pp']:6.1f} pp  "
              f"extremo oscuro: contorno {r['frac_extremo_oscuro_contorno']:5.1f} % / "
              f"anillo {r['frac_extremo_oscuro_anillo']:5.1f} %  "
              f"área {r['area_oscura_caja_km2']:5.2f} km²  "
              f"P(más oscuro) {r['prob_contorno_mas_oscuro']:.3f}  "
              f"({n_esc} filas de ese año en la plantilla)")


# ══════════════════════════════════════════════════════════════════════════════
# E. Comparación con la serie JRC
# ══════════════════════════════════════════════════════════════════════════════
def bloque_jrc(d, df):
    print()
    print("=" * 78)
    print("E. COMPARACIÓN CON LA SERIE JRC EN LOS MESES COMUNES")
    print("=" * 78)

    jrc = pd.read_csv(JRC).drop_duplicates(subset=["year", "month"])
    (a0, m0), (a1, m1) = VENTANA_COMUN
    clave = jrc["year"] * 100 + jrc["month"]
    jrc = jrc[(clave >= a0 * 100 + m0) & (clave <= a1 * 100 + m1)].copy()
    print(f"  ventana común: {a0}-{m0:02d} a {a1}-{m1:02d}  ({len(jrc)} meses JRC)")

    mss = (d.groupby(["year", "month"])
             .agg(escenas_mss=("scene_code", "size"),
                  fechas=("date", lambda s: " ".join(sorted(s))),
                  clasificacion=("clasificacion", lambda s: "/".join(sorted(set(s)))))
             .reset_index())

    t = jrc.merge(mss, on=["year", "month"], how="left")
    t["escenas_mss"] = t["escenas_mss"].fillna(0).astype(int)
    t["px_agua"] = t["px_water"].round(1)
    t["fraccion_observada"] = t["obs"].round(3)
    t = t[["year", "month", "escenas_mss", "clasificacion", "px_agua",
           "fraccion_observada", "fechas"]]
    print()
    print(t.drop(columns="fechas").to_string(index=False))

    obs = t[t["fraccion_observada"] > 0]
    print(f"\n  meses con observación JRC en la ventana: {len(obs)} de {len(t)}")
    for _, r in obs.iterrows():
        print(f"    {int(r['year'])}-{int(r['month']):02d}: JRC observa "
              f"{r['fraccion_observada']:.1%} de la caja, agua abierta "
              f"{r['px_agua']:.0f} px; MSS aporta {int(r['escenas_mss'])} escena(s) "
              f"{r['fechas'] if isinstance(r['fechas'], str) else ''}")
    print(f"  meses con escena MSS y observación JRC a la vez: "
          f"{int(((t['escenas_mss'] > 0) & (t['fraccion_observada'] > 0)).sum())}")
    print(f"  píxeles de agua abierta JRC en toda la ventana: {t['px_agua'].sum():.0f}")

    for anio in (1984, 1985):
        fila = df[df["anio"] == anio]
        if len(fila):
            r = fila.iloc[0]
            print(f"  chip MSS {anio}: contraste {r['contraste_pp']:.1f} pp, "
                  f"P(más oscuro) {r['prob_contorno_mas_oscuro']:.3f}, "
                  f"área en el extremo oscuro dentro de la caja "
                  f"{r['area_oscura_caja_km2']:.2f} km²")
    return t


# ══════════════════════════════════════════════════════════════════════════════
def main():
    for p in (INVENTARIO, PLANTILLA, VASO, JRC):
        if not p.exists():
            sys.exit(f"ERROR: no se encuentra {p}")
    SALIDA.mkdir(parents=True, exist_ok=True)

    d, tabla = bloque_cruce()
    d, chips = bloque_estructura(d)
    df = bloque_contraste(chips)
    bloque_distribucion(df, d)
    t = bloque_jrc(d, df)

    d.merge(df[["anio", "contraste_pp", "prob_contorno_mas_oscuro",
                "area_oscura_caja_km2"]],
            left_on="year", right_on="anio", how="left").drop(columns="anio") \
     .to_csv(SALIDA / "mss_cruce_escenas_fechas.csv", index=False)
    tabla.to_csv(SALIDA / "mss_escenas_anio_mes.csv")
    df.to_csv(SALIDA / "mss_contraste_anual.csv", index=False)
    t.to_csv(SALIDA / "mss_jrc_meses_comunes.csv", index=False)

    print()
    print("=" * 78)
    print("SALIDAS")
    print("=" * 78)
    for n in ("mss_cruce_escenas_fechas.csv", "mss_escenas_anio_mes.csv",
              "mss_contraste_anual.csv", "mss_jrc_meses_comunes.csv"):
        print(f"  {(SALIDA / n).relative_to(RAIZ)}")


if __name__ == "__main__":
    main()
