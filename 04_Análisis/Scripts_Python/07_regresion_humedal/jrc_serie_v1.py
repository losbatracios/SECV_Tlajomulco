"""
jrc_serie_v1.py
---------------
Procesa el CSV exportado por jrc_serie_v1.js (GEE) y produce:
  1. CSV limpio con estadísticas mensuales.
  2. CSV de resumen anual (max / mean área con agua; meses con obs).
  3. CSV de la serie de febrero únicamente.
  4. Gráfica en escala de grises: área con agua por año (febrero marcado).

Requisitos
----------
  pip install pandas matplotlib --break-system-packages

Uso
---
  Desde la raíz del proyecto:
    python 04_Análisis/Scripts_Python/07_regresion_humedal/jrc_serie_v1.py

Entradas (INPUT_DIR = 05_Datos_abiertos/JRC_GSW/)
--------------------------------------------------
  jrc_gsw_monthly_bbox_v1.csv
    Columnas esperadas del CSV de GEE:
      year, month, px_water, px_no_water, px_no_obs, area_ha_water, pct_bbox_area
    Nota: GEE puede exportar columnas adicionales (.geo, system:index); se ignoran.

Salidas (mismo INPUT_DIR)
-------------------------
  jrc_monthly_clean.csv        — datos mensuales limpios
  jrc_annual_summary.csv       — resumen anual
  jrc_february_series.csv      — solo febrero (mes de la foto de 1971)
  jrc_agua_por_anio.png        — gráfica en escala de grises

Notas metodológicas
-------------------
  • Periodo: 1984-03 a 2021-12 (JRC/GSW1_4/MonthlyHistory, v1.4).
  • px_no_obs > 0 = cobertura parcial o nula: se declara como 'sin observación',
    NUNCA como 'sin agua'.
  • Área de referencia: CAJA del vaso de 1971. Actualizar con polígono validado.
  • Cifras exploratorias hasta que [4. TJ Python] las replique.
  • Verbos admitidos: registra, no registra, coincide con, precede a,
    es posterior a. Sin inferencia de causa.

Fecha: 2026-10-26 · Proyecto: Tlajomulco — capítulo SECV-UNAM
"""

import sys
import pathlib
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ─── Rutas ────────────────────────────────────────────────────────────────────
ROOT       = pathlib.Path(__file__).resolve().parents[4]
INPUT_DIR  = ROOT / '05_Datos_abiertos' / 'JRC_GSW'
INPUT_CSV  = INPUT_DIR / 'jrc_gsw_monthly_bbox_v1.csv'
OUTPUT_DIR = INPUT_DIR

# ─── Verificación de entrada ──────────────────────────────────────────────────
if not INPUT_CSV.exists():
    print(
        f"\n[ERROR] No se encontró el archivo de entrada:\n  {INPUT_CSV}\n"
        "  Descarga el CSV desde Google Drive (carpeta JRC_GSW_Tlajomulco)\n"
        "  y colócalo en 05_Datos_abiertos/JRC_GSW/ antes de ejecutar.\n"
    )
    sys.exit(1)

# ─── 1. Leer y limpiar ────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)

cols_needed = ['year', 'month', 'px_water', 'px_no_water', 'px_no_obs',
               'area_ha_water', 'pct_bbox_area']
df = df[[c for c in cols_needed if c in df.columns]].copy()

for col in ['year', 'month', 'px_water', 'px_no_water', 'px_no_obs']:
    df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
df['area_ha_water'] = pd.to_numeric(df['area_ha_water'], errors='coerce')
df['pct_bbox_area'] = pd.to_numeric(df['pct_bbox_area'], errors='coerce')

df = df.sort_values(['year', 'month']).reset_index(drop=True)

total_px = (df['px_water'] + df['px_no_water'] + df['px_no_obs']).median()
df['obs_disponible'] = df['px_no_obs'] < (total_px * 0.9)

out_monthly = OUTPUT_DIR / 'jrc_monthly_clean.csv'
df.to_csv(out_monthly, index=False)
print(f"[OK] CSV mensual guardado: {out_monthly}")
print(f"     Registros totales    : {len(df)}")
print(f"     Periodo              : {df['year'].min()}-{df['month'].min():02d}"
      f" a {df['year'].max()}-{df['month'].max():02d}")

# ─── 2. Resumen anual ─────────────────────────────────────────────────────────
df_obs = df[df['obs_disponible']].copy()

annual = df_obs.groupby('year').agg(
    meses_con_obs  = ('month', 'count'),
    area_max_ha    = ('area_ha_water', 'max'),
    area_media_ha  = ('area_ha_water', 'mean'),
    pct_max        = ('pct_bbox_area', 'max'),
    pct_media      = ('pct_bbox_area', 'mean'),
    meses_con_agua = ('px_water', lambda x: (x > 0).sum()),
).reset_index()

total_months_per_year = df.groupby('year').size().rename('meses_totales')
annual = annual.merge(total_months_per_year, on='year')
annual['meses_sin_obs'] = annual['meses_totales'] - annual['meses_con_obs']

out_annual = OUTPUT_DIR / 'jrc_annual_summary.csv'
annual.to_csv(out_annual, index=False)
print(f"[OK] Resumen anual guardado: {out_annual}")

# ─── 3. Serie de febrero ──────────────────────────────────────────────────────
df_feb = df[df['month'] == 2].copy()
out_feb = OUTPUT_DIR / 'jrc_february_series.csv'
df_feb.to_csv(out_feb, index=False)
print(f"[OK] Serie de febrero guardada: {out_feb}")
print(f"     Registros de febrero      : {len(df_feb)}")
print(f"     Febrero con observación   : {df_feb['obs_disponible'].sum()}")
print(f"     Febrero con agua (>0 ha)  : {(df_feb['area_ha_water'] > 0).sum()}")

# ─── 4. Gráfica en escala de grises ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4.5))
ax.set_facecolor('#ffffff')
fig.patch.set_facecolor('#ffffff')

ax.bar(annual['year'], annual['area_max_ha'],
       color='#888888', alpha=0.65, width=0.7,
       label='Área máx. anual con agua (ha)')
ax.plot(annual['year'], annual['area_media_ha'],
        color='#222222', linewidth=1.2, linestyle='--',
        label='Área media anual (ha)')

df_feb_obs   = df_feb[df_feb['obs_disponible']]
df_feb_noobs = df_feb[~df_feb['obs_disponible']]

ax.scatter(df_feb_obs['year'], df_feb_obs['area_ha_water'],
           color='#000000', zorder=5, s=30, marker='o',
           label='Febrero (con observación)')
if len(df_feb_noobs) > 0:
    ax.scatter(df_feb_noobs['year'], [0] * len(df_feb_noobs),
               color='#aaaaaa', zorder=5, s=30, marker='x',
               label='Febrero (sin observación)')

for _, row in df_feb_obs[df_feb_obs['area_ha_water'] > 0].iterrows():
    ax.annotate(f"feb\n{int(row['year'])}",
                xy=(row['year'], row['area_ha_water']),
                xytext=(0, 6), textcoords='offset points',
                fontsize=6.5, ha='center', color='#111111')

ax.set_xlabel('Año', fontsize=10)
ax.set_ylabel('Área con agua (ha)', fontsize=10)
ax.set_title(
    'Vaso de La Playa — agua registrada por año (JRC GSW v1.4, 1984–2021)\n'
    'Caja del vaso de 1971 · Exploratorio · Cifras no citables hasta replicación de [4]',
    fontsize=9, color='#333333'
)
ax.xaxis.set_major_locator(mticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(mticker.MultipleLocator(1))
ax.tick_params(axis='both', labelsize=9)
ax.legend(fontsize=8, loc='upper right', framealpha=0.8)

fig.text(
    0.01, 0.01,
    'Producto: JRC/GSW1_4/MonthlyHistory (Pekel et al. 2016). '
    'Área de referencia: caja 20.4647–20.4796 N, −103.5442–−103.5100 O. '
    'Los meses sin observación no son "sin agua".',
    fontsize=6.5, color='#555555', va='bottom'
)

plt.tight_layout(rect=[0, 0.04, 1, 1])
out_png = OUTPUT_DIR / 'jrc_agua_por_anio.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"[OK] Gráfica guardada: {out_png}")

# ─── 5. Resumen en pantalla ───────────────────────────────────────────────────
print("\n─── Resumen rápido ───────────────────────────────────────────────────")
print(f"Registros totales                : {len(df)}")
print(f"Meses con observación disponible : {df['obs_disponible'].sum()}")
print(f"Meses sin observación            : {(~df['obs_disponible']).sum()}")
print(f"Meses con agua > 0 ha            : {(df['area_ha_water'] > 0).sum()}")
print(f"Área máx. registrada (ha)        : {df['area_ha_water'].max():.1f}")
print(f"Área media (meses con obs, ha)   : {df_obs['area_ha_water'].mean():.1f}")
print(f"Febrero con observación          : {df_feb_obs.shape[0]} de {len(df_feb)}")
print(f"Febrero con agua > 0 ha          : {(df_feb_obs['area_ha_water'] > 0).sum()}")
print("\nNOTA: Cifras exploratorias. No citar hasta replicación de [4. TJ Python].")
print("      Los meses con obs_disponible = False son sin observación, no sin agua.")
