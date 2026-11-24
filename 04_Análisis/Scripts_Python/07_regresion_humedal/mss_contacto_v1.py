#!/usr/bin/env python3
"""
mss_contacto_v1.py — Hojas de contacto Landsat MSS + plantilla de clasificación ciega
Proyecto: Tlajomulco · [3. TJ Datos] · despacho b1c.7
Fecha de creación: 2026-10-31

REQUISITOS (instalar si no están):
    pip install rasterio geopandas matplotlib numpy pandas pillow

PREREQUISITOS:
    1. mss_inventory.csv descargado de GEE a:
       05_Datos_abiertos/Landsat_MSS/mss_inventory.csv
    2. Chips TIF anuales descargados de GEE a:
       05_Datos_abiertos/Landsat_MSS/chips/mss_chip_YYYY.tif
    3. vaso_1971.gpkg disponible en:
       04_Análisis/Visualizaciones/05_humedal_georef/vaso_1971.gpkg

USO:
    cd /ruta/a/Tlajomulco
    python 04_Análisis/Scripts_Python/07_regresion_humedal/mss_contacto_v1.py

SALIDA:
    05_Datos_abiertos/Landsat_MSS/hojas_contacto/hoja_contacto_MSS.png
    05_Datos_abiertos/Landsat_MSS/mss_clasificacion_plantilla.csv
    (la plantilla también se regenera aquí para mantener coherencia con los chips)

PROTOCOLO DE CLASIFICACIÓN CIEGA:
    - Los chips se ordenan con semilla aleatoria fija (seed=42)
    - La fecha NO aparece en la hoja de contacto (solo el código de escena)
    - La clasificación (agua visible / no visible / no legible) se anota en la plantilla
    - Fecha revelada SOLO después de completar la clasificación
"""

import pathlib
import sys
import random
import csv
import json
from datetime import datetime

import numpy as np
import pandas as pd

# ─── Verificar dependencias opcionales ─────────────────────────────────────────
try:
    import rasterio
    from rasterio.plot import show as rio_show
    from rasterio.warp import reproject, Resampling
    from rasterio.crs import CRS
    import rasterio.mask
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    print("ADVERTENCIA: rasterio no disponible. Instalar con: pip install rasterio")

try:
    import geopandas as gpd
    from shapely.geometry import box
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    print("ADVERTENCIA: geopandas no disponible. Instalar con: pip install geopandas")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("ADVERTENCIA: matplotlib no disponible. Instalar con: pip install matplotlib")

# ─── Rutas ─────────────────────────────────────────────────────────────────────
ROOT = pathlib.Path(__file__).resolve().parents[3]   # raíz Tlajomulco/

INVENTORY_CSV   = ROOT / '05_Datos_abiertos/Landsat_MSS/mss_inventory.csv'
CHIPS_DIR       = ROOT / '05_Datos_abiertos/Landsat_MSS/chips'
VASO_GPKG       = ROOT / '04_Análisis/Visualizaciones/05_humedal_georef/vaso_1971.gpkg'
OUTPUT_DIR      = ROOT / '05_Datos_abiertos/Landsat_MSS/hojas_contacto'
PLANTILLA_CSV   = ROOT / '05_Datos_abiertos/Landsat_MSS/mss_clasificacion_plantilla.csv'

# Bbox del vaso (WGS84)
BBOX_WGS84 = (-103.5442, 20.4647, -103.5100, 20.4796)  # (xmin, ymin, xmax, ymax)

# Semilla fija para orden aleatorio en clasificación ciega
SEED_CLASIFICACION = 42

# Criterio de agua visible (documentado aquí y en README)
CRITERIO_AGUA = "parche oscuro en NIR de al menos ~1 km² dentro del bbox del vaso"


def verificar_prerequisitos():
    """Verifica que los archivos necesarios existan."""
    errores = []
    if not INVENTORY_CSV.exists():
        errores.append(f"  FALTA: {INVENTORY_CSV}")
        errores.append("  → Exportar desde GEE con mss_vaso_v1.js y copiar aquí")
    if not VASO_GPKG.exists():
        errores.append(f"  FALTA: {VASO_GPKG}")

    if errores:
        print("ERROR — archivos requeridos no encontrados:")
        for e in errores:
            print(e)
        sys.exit(1)

    if not CHIPS_DIR.exists():
        CHIPS_DIR.mkdir(parents=True, exist_ok=True)
        print(f"  Creado directorio de chips: {CHIPS_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"  Directorio de salida: {OUTPUT_DIR}")


def cargar_inventario():
    """Carga el CSV de inventario y devuelve DataFrame."""
    df = pd.read_csv(INVENTORY_CSV)
    df['year']  = pd.to_numeric(df['year'],  errors='coerce').astype(int)
    df['month'] = pd.to_numeric(df['month'], errors='coerce').astype(int)
    df['cloud_cover_pct'] = pd.to_numeric(df['cloud_cover_pct'], errors='coerce')
    df = df.sort_values('date').reset_index(drop=True)
    print(f"\n  Inventario: {len(df)} escenas encontradas")
    print(f"  Periodo: {df['date'].min()} — {df['date'].max()}")
    print(f"  Satélites: {df['satellite'].unique().tolist()}")
    print(f"  Tiers: {df['tier'].unique().tolist()}")
    return df


def verificar_umbral(df):
    """
    Aplica el umbral de viabilidad del análisis:
    Si hay <3 años distintos con escenas en 1972–1984 → b1c.7 queda fuera.
    """
    df_rango = df[(df['year'] >= 1972) & (df['year'] <= 1984)]
    anios_con_escenas = sorted(df_rango['year'].unique().tolist())
    n_anios = len(anios_con_escenas)

    print(f"\n  ──────────────────────────────────────────")
    print(f"  UMBRAL DE VIABILIDAD b1c.7")
    print(f"  Años con escenas en 1972–1984: {n_anios}")
    print(f"  Años: {anios_con_escenas}")
    if n_anios < 3:
        print(f"\n  ⚠ RESULTADO: b1c.7 QUEDA FUERA DE ALCANCE")
        print(f"  Menos de 3 años con escenas en 1972–1984.")
        print(f"  Reportar como límite metodológico en el capítulo.")
        print(f"  Formulación sugerida:")
        print(f"    'Landsat MSS registra solo {n_anios} año(s) con cobertura")
        print(f"     disponible sobre el vaso en 1972–1984 (colecciones LM01–LM05,")
        print(f"     C02/T1+T2); el análisis visual no es viable con esta densidad")
        print(f"     temporal y se declara fuera de alcance.'")
    else:
        print(f"\n  ✓ UMBRAL SUPERADO: análisis viable ({n_anios} años con escenas)")
    print(f"  ──────────────────────────────────────────")

    return n_anios, anios_con_escenas


def generar_plantilla_clasificacion(df, seed=SEED_CLASIFICACION):
    """
    Genera la plantilla CSV de clasificación ciega.
    - Ordena las escenas aleatoriamente con semilla fija
    - Oculta la fecha (incluye solo scene_code)
    - scene_code = SAT_TIER_YYYYMMDD (único y reproducible)
    """
    # Construir scene_code
    def make_code(row):
        d = str(row['date']).replace('-', '')
        return f"{row['satellite']}_{row['tier']}_{d}"

    df = df.copy()
    df['scene_code'] = df.apply(make_code, axis=1)

    # Orden aleatorio con semilla fija
    rng = random.Random(seed)
    indices = list(range(len(df)))
    rng.shuffle(indices)
    df_shuffled = df.iloc[indices].reset_index(drop=True)

    # Plantilla: solo scene_code, sin fecha
    plantilla = pd.DataFrame({
        'orden_clasificacion': range(1, len(df_shuffled) + 1),
        'scene_code':   df_shuffled['scene_code'],
        'clasificacion': '',      # COMPLETAR: agua_visible / no_visible / no_legible
        'nota':          ''       # COMPLETAR: observación libre
    })

    plantilla.to_csv(PLANTILLA_CSV, index=False, encoding='utf-8')
    print(f"\n  Plantilla de clasificación escrita: {PLANTILLA_CSV}")
    print(f"  Semilla: seed={seed}")
    print(f"  Criterio de 'agua_visible': {CRITERIO_AGUA}")
    print(f"  Valores permitidos en columna 'clasificacion':")
    print(f"    agua_visible | no_visible | no_legible")
    print(f"  IMPORTANTE: completar la clasificación SIN ver las fechas.")
    print(f"  Las fechas están en mss_inventory.csv; no consultarlo hasta terminar.")

    return df_shuffled, df['scene_code'].tolist()


def armar_hoja_contacto(df, df_shuffled):
    """
    Ensambla una hoja de contacto PNG con un chip por fila.
    - Busca chips en CHIPS_DIR (archivos mss_chip_YYYY.tif)
    - Overlay: contorno del bbox y del vaso_1971.gpkg
    - Escala de grises uniforme (percentil 2–98 global)
    - Muestra scene_code, NO la fecha
    - Orden: el del df_shuffled (aleatorio con semilla fija)
    """
    if not (HAS_RASTERIO and HAS_GEOPANDAS and HAS_MPL):
        print("\n  HOJA DE CONTACTO OMITIDA: faltan dependencias.")
        print("  Instalar: pip install rasterio geopandas matplotlib")
        return

    # Cargar vaso_1971.gpkg
    try:
        vaso_gdf = gpd.read_file(VASO_GPKG)
        print(f"  Polígono vaso_1971 cargado: {len(vaso_gdf)} geometría(s), CRS={vaso_gdf.crs}")
    except Exception as e:
        print(f"  ADVERTENCIA: no se pudo cargar vaso_1971.gpkg: {e}")
        vaso_gdf = None

    # Buscar chips disponibles (archivos por año: mss_chip_YYYY.tif)
    chips_disponibles = {}
    for f in sorted(CHIPS_DIR.glob('mss_chip_*.tif')):
        yr_str = f.stem.replace('mss_chip_', '')
        try:
            yr = int(yr_str)
            chips_disponibles[yr] = f
        except ValueError:
            pass

    if not chips_disponibles:
        print(f"  SIN CHIPS DISPONIBLES en {CHIPS_DIR}")
        print(f"  → Exportar chips con GEE (mss_vaso_v1.js → Tasks → mss_chip_YYYY)")
        print(f"  → Copiar TIF a {CHIPS_DIR}")
        print(f"  → Ejecutar este script de nuevo")
        return

    print(f"  Chips disponibles: años {sorted(chips_disponibles.keys())}")

    # Calcular estiramiento global sobre todos los chips disponibles
    all_vals = []
    for yr, path in chips_disponibles.items():
        try:
            with rasterio.open(path) as src:
                data = src.read(1).astype(float)
                if src.nodata is not None:
                    data[data == src.nodata] = np.nan
                valid = data[~np.isnan(data)]
                if len(valid) > 0:
                    all_vals.extend(valid[::10].tolist())  # muestra cada 10px
        except Exception as e:
            print(f"  ADVERTENCIA chip {yr}: {e}")

    if not all_vals:
        print("  ERROR: no se pudo leer ningún chip para calcular estiramiento.")
        return

    vmin = float(np.percentile(all_vals, 2))
    vmax = float(np.percentile(all_vals, 98))
    print(f"  Estiramiento global NIR: p2={vmin:.1f}, p98={vmax:.1f}")

    # Determinar años que aparecen en df_shuffled
    def year_from_code(code):
        # code = SAT_TIER_YYYYMMDD → año = caracteres 8–11 del final
        parts = code.split('_')
        try:
            return int(parts[2][:4])
        except Exception:
            return None

    df_shuffled_copy = df_shuffled.copy()
    df_shuffled_copy['chip_year'] = df_shuffled_copy['scene_code'].apply(year_from_code)
    df_shuffled_copy['chip_path'] = df_shuffled_copy['chip_year'].map(chips_disponibles)
    df_plot = df_shuffled_copy[df_shuffled_copy['chip_path'].notna()].reset_index(drop=True)

    if len(df_plot) == 0:
        print("  No hay chips TIF descargados que coincidan con las escenas del inventario.")
        print("  → Descargar chips de Drive y colocar en:", CHIPS_DIR)
        return

    n_chips = len(df_plot)
    n_cols  = min(4, n_chips)
    n_rows  = (n_chips + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(n_cols * 4, n_rows * 4),
                              constrained_layout=True)
    if n_rows * n_cols == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]

    # Bbox como geometría para overlay
    bbox_geom = box(*BBOX_WGS84)
    bbox_gdf  = gpd.GeoDataFrame(geometry=[bbox_geom], crs='EPSG:4326')

    for idx, row in df_plot.iterrows():
        r, c = divmod(idx, n_cols)
        ax = axes[r][c]

        try:
            with rasterio.open(row['chip_path']) as src:
                # Reproyectar vaso y bbox a CRS del chip si es necesario
                chip_crs = src.crs

                data = src.read(1).astype(float)
                if src.nodata is not None:
                    data[data == src.nodata] = np.nan

                extent = [src.bounds.left, src.bounds.right,
                          src.bounds.bottom, src.bounds.top]

                ax.imshow(data,
                          cmap='gray',
                          vmin=vmin, vmax=vmax,
                          extent=extent,
                          origin='upper',
                          interpolation='nearest')

                # Overlay bbox (en WGS84 o reproyectado)
                try:
                    bbox_repr = bbox_gdf.to_crs(chip_crs)
                    for geom in bbox_repr.geometry:
                        x, y = geom.exterior.xy
                        ax.plot(x, y, 'r-', linewidth=1.5, label='bbox' if idx == 0 else '')
                except Exception:
                    pass

                # Overlay vaso_1971
                if vaso_gdf is not None:
                    try:
                        vaso_repr = vaso_gdf.to_crs(chip_crs)
                        for geom in vaso_repr.geometry:
                            if geom.geom_type == 'Polygon':
                                x, y = geom.exterior.xy
                                ax.plot(x, y, 'y-', linewidth=1.5,
                                        label='vaso_1971' if idx == 0 else '')
                            elif geom.geom_type == 'MultiPolygon':
                                for poly in geom.geoms:
                                    x, y = poly.exterior.xy
                                    ax.plot(x, y, 'y-', linewidth=1.5)
                    except Exception:
                        pass

        except Exception as e:
            ax.text(0.5, 0.5, f"Error:\n{e}", transform=ax.transAxes,
                    ha='center', va='center', fontsize=6, color='red')

        # Título: scene_code SIN fecha
        ax.set_title(row['scene_code'], fontsize=7, pad=3)
        ax.set_xticks([])
        ax.set_yticks([])

    # Ocultar ejes vacíos
    total_cells = n_rows * n_cols
    for blank_idx in range(len(df_plot), total_cells):
        r, c = divmod(blank_idx, n_cols)
        axes[r][c].set_visible(False)

    # Leyenda global
    patches = [
        mpatches.Patch(color='red',    label='Bbox vaso (WGS84)'),
        mpatches.Patch(color='yellow', label='Vaso 1971 (~416 ha)')
    ]
    fig.legend(handles=patches, loc='lower center', ncol=2,
               fontsize=8, framealpha=0.8)

    fig.suptitle(
        f'Landsat MSS NIR — Vaso de La Playa · 1972–1985\n'
        f'Estiramiento uniforme p2–p98 · Orden aleatorio seed={SEED_CLASIFICACION}\n'
        f'CLASIFICACIÓN CIEGA: fechas ocultas · ver mss_clasificacion_plantilla.csv',
        fontsize=9
    )

    out_path = OUTPUT_DIR / 'hoja_contacto_MSS.png'
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  ✓ Hoja de contacto guardada: {out_path}")


def main():
    print("="*60)
    print("mss_contacto_v1.py — Hojas de contacto MSS + plantilla")
    print(f"Ejecutado: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)

    verificar_prerequisitos()
    df = cargar_inventario()
    n_anios, anios = verificar_umbral(df)

    # Generar plantilla independientemente del umbral
    df_shuffled, codes = generar_plantilla_clasificacion(df)

    # Intentar armar hoja de contacto
    armar_hoja_contacto(df, df_shuffled)

    print("\n" + "="*60)
    print("RESUMEN")
    print(f"  Escenas en inventario: {len(df)}")
    print(f"  Años con escenas (1972–1984): {n_anios} {anios}")
    if n_anios < 3:
        print(f"  → b1c.7 FUERA DE ALCANCE (< 3 años)")
    else:
        print(f"  → b1c.7 VIABLE")
    print(f"  Plantilla CSV: {PLANTILLA_CSV}")
    print(f"  Hojas de contacto: {OUTPUT_DIR}")
    print("="*60)
    print("\nPRÓXIMOS PASOS:")
    print("  1. Revisar la hoja de contacto (hoja_contacto_MSS.png)")
    print("  2. Completar mss_clasificacion_plantilla.csv con:")
    print("     agua_visible | no_visible | no_legible")
    print("  3. Registrar clasificación ANTES de ver las fechas")
    print("  4. Luego cruzar con mss_inventory.csv para revelar fechas")
    print("  5. Enviar resultados a [4. TJ Python] (ronda 6) para validación")
    print("\nNOTA: Todas las cifras son exploratorias hasta validación en [4].")


if __name__ == '__main__':
    main()
